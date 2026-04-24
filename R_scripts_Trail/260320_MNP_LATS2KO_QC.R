library(tidyverse)
library(edgeR)
library(limma)
library(biomaRt)

## ── Setup ──────────────────────────────────────────────────────────────────────
BASE_DIR <- '/scratch/rmlab/rmlab_shared3/alferreiro/BRB_Seq'
ANA_DIR  <- file.path(BASE_DIR, '260320_MNP_WT_LATS2KO_Rot_TRULI_Timecourse_Analyses')
DATA_DIR <- file.path(ANA_DIR, 'data')
FIG_DIR  <- file.path(ANA_DIR, 'qc_figures')
LOG_DIR  <- file.path(ANA_DIR, 'logs')

for (d in c(DATA_DIR, FIG_DIR, LOG_DIR)) dir.create(d, showWarnings = FALSE, recursive = TRUE)
setwd(BASE_DIR)

log_file <- file.path(LOG_DIR, '260320_MNP_LATS2KO_QC.log')
cat('', file = log_file)
log_msg <- function(...) {
  msg <- paste0('[', format(Sys.time(), '%Y-%m-%d %H:%M:%S'), '] ', ...)
  message(msg); cat(msg, '\n', file = log_file, append = TRUE)
}
log_msg('Script started')

## ── Load counts ────────────────────────────────────────────────────────────────
COUNTS_FILE <- '/lts/rmlab/rmlab_shared3/alferreiro/BRB_Seq/NRP_FullCounts/MD912_MNP_LATS2KO_TRULI_Dedup_Counts.txt'
META_FILE   <- '/lts/rmlab/rmlab_shared3/alferreiro/BRB_Seq/NRP_FullCounts/MD912_MNP_LATS2KO_TRULI_metadata.csv'

sep_char <- if (grepl('\t', readLines(COUNTS_FILE, n = 1))) '\t' else ' '
counts   <- read.table(COUNTS_FILE, header = TRUE, row.names = 1,
                       sep = sep_char, check.names = FALSE)
meta     <- read.csv(META_FILE)
log_msg('Counts loaded: ', nrow(counts), ' genes x ', ncol(counts), ' samples')
log_msg('Metadata loaded: ', nrow(meta), ' rows')

## ── Strip Ensembl version numbers (GRCh38) ────────────────────────────────────
rownames(counts) <- sub('\\.\\d+$', '', rownames(counts))

## ── Align metadata to counts ──────────────────────────────────────────────────
if (all(meta$OriginalName %in% colnames(counts))) {
  counts           <- counts[, meta$OriginalName, drop = FALSE]
  colnames(counts) <- meta$SampleName
  log_msg('Columns matched via OriginalName → renamed to SampleName')
} else if (all(meta$SampleName %in% colnames(counts))) {
  counts <- counts[, meta$SampleName, drop = FALSE]
  log_msg('Columns matched via SampleName')
} else {
  stop('Count columns match neither OriginalName nor SampleName')
}
rownames(meta) <- meta$SampleName
meta <- meta[colnames(counts), ]
stopifnot(all(rownames(meta) == colnames(counts)))
log_msg('Metadata aligned: ', ncol(counts), ' samples matched')

## ── Factor levels ──────────────────────────────────────────────────────────────
meta <- meta %>% mutate(
  Genotype  = factor(Genotype,  levels = c('WT', 'LATS2KO')),
  Treatment = factor(Treatment, levels = c('Vehicle', 'TRULI', 'Rotenone', 'TRULI + Rotenone')),
  Timepoint = factor(as.character(Timepoint), levels = c('6', '16', '24', '48', '96')),
  Run       = factor(Run)
)

## ── Color palettes ─────────────────────────────────────────────────────────────
paletteGenotype  <- c('WT' = '#555555', 'LATS2KO' = '#3182BD')
paletteTreatment <- c('Vehicle' = '#969696', 'TRULI' = '#74C476',
                      'Rotenone' = '#FC8D59', 'TRULI + Rotenone' = '#31A354')
paletteTimepoint <- c('6' = '#FEE391', '16' = '#FEB24C', '24' = '#FC4E2A', '48' = '#800026', '96' = '#4D004B')

## ── Gene symbol lookup (biomaRt, GRCh38) ──────────────────────────────────────
log_msg('Querying biomaRt for gene symbols (GRCh38)...')
ensembl_ids <- rownames(counts)
mart <- tryCatch(
  useEnsembl(biomart = 'genes', dataset = 'hsapiens_gene_ensembl'),
  error = function(e) {
    log_msg('useEnsembl failed (', conditionMessage(e), '), trying useMart...')
    useMart('ensembl', dataset = 'hsapiens_gene_ensembl')
  }
)
bm_result <- getBM(
  attributes = c('ensembl_gene_id', 'external_gene_name'),
  filters    = 'ensembl_gene_id',
  values     = unique(ensembl_ids),
  mart       = mart
)
bm_map      <- bm_result %>%
  filter(external_gene_name != '') %>%
  distinct(ensembl_gene_id, .keep_all = TRUE) %>%
  { setNames(.$external_gene_name, .$ensembl_gene_id) }
gene_symbol <- bm_map[ensembl_ids]
names(gene_symbol) <- ensembl_ids
no_symbol <- is.na(gene_symbol) | gene_symbol == ''
gene_symbol[no_symbol] <- ensembl_ids[no_symbol]
log_msg('Gene symbols: ', sum(!no_symbol), ' mapped, ', sum(no_symbol), ' fell back to ID')

gene_data <- data.frame(
  ensembl_gene_id = ensembl_ids,
  gene_symbol     = unname(gene_symbol),
  row.names       = ensembl_ids
)

## ── DGEList, filter, normalise, save ──────────────────────────────────────────
dge      <- DGEList(counts = counts, samples = meta, genes = gene_data)
n_before <- nrow(dge)
keep     <- filterByExpr(dge, group = dge$samples$Group)
dge      <- dge[keep, , keep.lib.sizes = FALSE]
log_msg('filterByExpr: retained ', nrow(dge), ' / ', n_before, ' genes')
dge <- calcNormFactors(dge, method = 'TMM')
saveRDS(dge, file.path(DATA_DIR, '260320_MNP_LATS2KO_DGElist.RDS'))
log_msg('Saved DGEList: data/260320_MNP_LATS2KO_DGElist.RDS')

## ── QC: Overview ──────────────────────────────────────────────────────────────
threshold <- 2e6
p_overview <- ggplot(meta, aes(x = Timepoint, y = NumUnqDedupedReads)) +
  geom_hline(yintercept = threshold, linetype = 'dashed', colour = 'darkred', linewidth = 0.8) +
  geom_boxplot(outlier.shape = NA, alpha = 0.4, fill = 'gray95', colour = 'grey40') +
  geom_jitter(aes(colour = Treatment,
                  shape  = NumUnqDedupedReads < threshold),
              width = 0.2, size = 2.5, alpha = 0.85) +
  facet_wrap(~ Genotype, nrow = 1) +
  scale_colour_manual(values = paletteTreatment, name = 'Treatment') +
  scale_shape_manual(values = c('TRUE' = 4, 'FALSE' = 16),
                     name = 'Below 2M', labels = c('TRUE' = 'Yes', 'FALSE' = 'No')) +
  scale_y_continuous(labels = scales::unit_format(unit = 'M', scale = 1e-6)) +
  labs(title    = 'Sequencing Effort — MNP LATS2KO TRULI Timecourse',
       subtitle = 'Dashed line: 2M unique deduped reads; × = below threshold',
       x = 'Timepoint (hr)', y = 'Unique Deduped Reads') +
  theme_bw() + theme(legend.position = 'bottom')
ggsave(file.path(FIG_DIR, '260320_MNP_LATS2KO_UMIDeDuped_ReadCounts_Overview.png'),
       p_overview, width = 10, height = 5)
log_msg('Saved: qc_figures/260320_MNP_LATS2KO_UMIDeDuped_ReadCounts_Overview.png')

## ── QC: By group ──────────────────────────────────────────────────────────────
p_group <- ggplot(meta, aes(x = Treatment, y = NumUnqDedupedReads)) +
  geom_hline(yintercept = threshold, linetype = 'dashed', colour = 'darkred', linewidth = 0.8) +
  geom_boxplot(outlier.shape = NA, alpha = 0.4, fill = 'gray95') +
  geom_jitter(aes(colour = NumUnqDedupedReads < threshold), width = 0.2, size = 2, alpha = 0.8) +
  facet_grid(Genotype ~ Timepoint) +
  scale_colour_manual(values = c('TRUE' = 'red', 'FALSE' = 'black'),
                      name = 'Below Threshold', labels = c('TRUE' = '< 2M', 'FALSE' = '>= 2M')) +
  scale_y_continuous(labels = scales::unit_format(unit = 'M', scale = 1e-6)) +
  labs(title = 'Sequencing Effort by Group — MNP LATS2KO TRULI Timecourse',
       x = 'Treatment', y = 'Unique Deduped Reads') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        strip.text  = element_text(size = 8), legend.position = 'bottom')
ggsave(file.path(FIG_DIR, '260320_MNP_LATS2KO_UMIDeDuped_ReadCounts_ByGroup.png'),
       p_group, width = 12, height = 6)
log_msg('Saved: qc_figures/260320_MNP_LATS2KO_UMIDeDuped_ReadCounts_ByGroup.png')

log_msg('Script complete')
