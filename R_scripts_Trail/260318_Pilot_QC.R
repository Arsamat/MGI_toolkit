library(tidyverse)
library(edgeR)
library(limma)
library(biomaRt)

## ── Setup ──────────────────────────────────────────────────────────────────────
BASE_DIR  <- '/scratch/rmlab/rmlab_shared3/alferreiro/BRB_Seq'
ANA_DIR   <- file.path(BASE_DIR, '260318_Pilot_Analyses')
DATA_DIR  <- file.path(ANA_DIR, 'data')
FIG_DIR   <- file.path(ANA_DIR, 'qc_figures')
LOG_DIR   <- file.path(ANA_DIR, 'logs')

dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR,  showWarnings = FALSE, recursive = TRUE)
dir.create(LOG_DIR,  showWarnings = FALSE, recursive = TRUE)

setwd(BASE_DIR)

log_file <- file.path(LOG_DIR, '260318_Pilot_QC.log')
cat('', file = log_file)
log_msg <- function(...) {
  msg <- paste0('[', format(Sys.time(), '%Y-%m-%d %H:%M:%S'), '] ', ...)
  message(msg)
  cat(msg, '\n', file = log_file, append = TRUE)
}
log_msg('Script started')
log_msg('ANA_DIR  : ', ANA_DIR)
log_msg('DATA_DIR : ', DATA_DIR)
log_msg('FIG_DIR  : ', FIG_DIR)

## ── Load ───────────────────────────────────────────────────────────────────────
COUNTS_FILE <- '/lts/rmlab/rmlab_shared3/alferreiro/BRB_Seq/NRP_FullCounts/Pilot_12S_Counts_Dedup.txt'
META_FILE   <- '/lts/rmlab/rmlab_shared3/alferreiro/BRB_Seq/NRP_FullCounts/Pilot_12S_metadata.csv'

pilot.counts <- read.table(
  COUNTS_FILE,
  header = TRUE, row.names = 1, sep = ' ', check.names = FALSE
)
pilot.meta <- read.csv(META_FILE)

log_msg('Counts: ', nrow(pilot.counts), ' genes x ', ncol(pilot.counts), ' samples')
log_msg('Metadata: ', nrow(pilot.meta), ' samples')

## ── Align metadata to count column order ──────────────────────────────────────
# Pilot counts columns match SampleName (not OriginalName)
if (!all(colnames(pilot.counts) %in% pilot.meta$SampleName)) {
  missing <- setdiff(colnames(pilot.counts), pilot.meta$SampleName)
  stop('Count columns not found in pilot.meta$SampleName: ',
       paste(missing, collapse = ', '))
}
row.names(pilot.meta) <- pilot.meta$SampleName
pilot.meta <- pilot.meta[colnames(pilot.counts), ]
stopifnot(all(rownames(pilot.meta) == colnames(pilot.counts)))
log_msg('Metadata aligned: ', ncol(pilot.counts), ' samples matched')

## ── Factor levels ──────────────────────────────────────────────────────────────
pilot.meta <- pilot.meta %>%
  mutate(
    Genotype = factor(Genotype, levels = c(
      'WT', 'TDP43'
    )),
    Treatment = factor(Treatment, levels = c(
      'NoTreatment', 'Vehicle', 'Rotenone', 'Thapsigargin'
    )),
    Dose = factor(Dose, levels = c(
      'NoTreatment', 'Vehicle', 'EC10', 'EC50', 'EC90'
    )),
    Treatment2 = factor(Treatment2, levels = c(
      'NoTreatment', 'Vehicle',
      'Rotenone_EC10', 'Rotenone_EC50', 'Rotenone_EC90',
      'Thapsigargin_EC10', 'Thapsigargin_EC50', 'Thapsigargin_EC90'
    )),
    Timepoint = factor(Timepoint, levels = c('0', '6', '24', '72', '96')),
    Run       = factor(Run, levels = c(
      'Pilot_Spinoff', 'Pilot_Spinoff2', 'Pilot_Batch1', 'Pilot_Batch2'
    ))
  )

## ── Color palettes ─────────────────────────────────────────────────────────────
paletteGenotype <- c(
  'WT'    = '#555555',
  'TDP43' = '#E41A1C'
)

paletteTreatment <- c(
  'NoTreatment'       = '#D9D9D9',
  'Vehicle'           = '#969696',
  'Rotenone_EC10'     = '#FDCC8A',
  'Rotenone_EC50'     = '#FC8D59',
  'Rotenone_EC90'     = '#D7301F',
  'Thapsigargin_EC10' = '#BDD7E7',
  'Thapsigargin_EC50' = '#6BAED6',
  'Thapsigargin_EC90' = '#08519C'
)

paletteTimepoint <- c(
  '0'  = '#BDBDBD',
  '6'  = '#FEE391',
  '24' = '#FE9929',
  '72' = '#CC4C02',
  '96' = '#7F2704'
)

paletteRun <- c(
  'Pilot_Spinoff'  = '#1B9E77',
  'Pilot_Spinoff2' = '#66C2A5',
  'Pilot_Batch1'   = '#D95F02',
  'Pilot_Batch2'   = '#FDAE6B'
)

shapesGenotype <- c(
  'WT'    = 22,
  'TDP43' = 23
)

## ── Gene symbol lookup (biomaRt) ───────────────────────────────────────────────
# Pilot was aligned to GRCh37; IDs are unversioned (e.g. ENSG00000001617)
log_msg('Querying biomaRt for gene symbols (GRCh37)...')

ensembl_ids <- rownames(pilot.counts)

mart <- tryCatch(
  useEnsembl(biomart = 'genes', dataset = 'hsapiens_gene_ensembl', GRCh = 37),
  error = function(e) {
    log_msg('GRCh37 useEnsembl failed (', conditionMessage(e),
            '), trying current Ensembl...')
    useMart('ensembl', dataset = 'hsapiens_gene_ensembl')
  }
)

bm_result <- getBM(
  attributes = c('ensembl_gene_id', 'external_gene_name'),
  filters    = 'ensembl_gene_id',
  values     = unique(ensembl_ids),
  mart       = mart
)

bm_map <- bm_result %>%
  filter(external_gene_name != '') %>%
  distinct(ensembl_gene_id, .keep_all = TRUE) %>%
  { setNames(.$external_gene_name, .$ensembl_gene_id) }

gene_symbol        <- bm_map[ensembl_ids]
names(gene_symbol) <- ensembl_ids
no_symbol          <- is.na(gene_symbol) | gene_symbol == ''
gene_symbol[no_symbol] <- ensembl_ids[no_symbol]

log_msg('Gene symbol mapping: ', sum(!no_symbol), ' mapped, ', sum(no_symbol),
        ' fell back to Ensembl ID (total ', length(gene_symbol), ')')

gene_data <- data.frame(
  ensembl_gene_id = ensembl_ids,
  gene_symbol     = unname(gene_symbol),
  row.names       = ensembl_ids
)

## ── DGEList, filtering, normalisation ─────────────────────────────────────────
dge      <- DGEList(counts = pilot.counts, samples = pilot.meta, genes = gene_data)
n_before <- nrow(dge)
keep     <- filterByExpr(dge, group = dge$samples$Group)
dge      <- dge[keep, , keep.lib.sizes = FALSE]
n_after  <- nrow(dge)
log_msg('filterByExpr: retained ', n_after, ' / ', n_before,
        ' genes (dropped ', n_before - n_after, ')')

dge <- calcNormFactors(dge, method = 'TMM')
saveRDS(dge, file.path(DATA_DIR, '260318_Pilot_DGElist.RDS'))
log_msg('Saved DGEList RDS: data/260318_Pilot_DGElist.RDS')

## ── QC: Sequencing effort — all samples overview ──────────────────────────────
# All samples are Experiment type, so no separate controls plot.
# Instead, show an overview by Run coloured by Genotype.
threshold <- 2000000

p_overview <- ggplot(pilot.meta,
                     aes(x = Run, y = NumUnqDedupedReads, colour = Genotype)) +
  geom_hline(yintercept = threshold, linetype = 'dashed',
             colour = 'darkred', linewidth = 0.8) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, fill = 'gray95',
               colour = 'grey40') +
  geom_jitter(aes(colour = Genotype,
                  shape  = NumUnqDedupedReads < threshold),
              width = 0.2, size = 2, alpha = 0.8) +
  scale_colour_manual(values = paletteGenotype, name = 'Genotype') +
  scale_shape_manual(values = c('TRUE' = 4, 'FALSE' = 16),
                     name   = 'Below 2M threshold',
                     labels = c('TRUE' = 'Yes', 'FALSE' = 'No')) +
  scale_y_continuous(labels = scales::unit_format(unit = 'M', scale = 1e-6)) +
  labs(title    = 'Sequencing Effort — All Samples by Run',
       subtitle = 'Dashed line: 2 million unique deduped reads; × = below threshold',
       x = 'Run', y = 'Unique Deduped Reads (Millions)') +
  theme_bw() +
  theme(axis.text.x   = element_text(angle = 45, hjust = 1),
        legend.position = 'bottom')

ggsave(file.path(FIG_DIR, '260318_UMIDeDuped_ReadCounts_Overview.png'),
       p_overview, width = 9, height = 5)
log_msg('Saved: qc_figures/260318_UMIDeDuped_ReadCounts_Overview.png')

## ── QC: Sequencing effort — per Genotype × Treatment combo ───────────────────
meta_exp <- pilot.meta %>% filter(SampleType == 'Experiment')

# Combos: non-WT genotypes paired with active treatments
exp_combos <- meta_exp %>%
  filter(Genotype != 'WT', Treatment %in% c('Rotenone', 'Thapsigargin')) %>%
  distinct(Genotype, Treatment) %>%
  arrange(Genotype, Treatment)

for (i in seq_len(nrow(exp_combos))) {
  geno <- as.character(exp_combos$Genotype[i])
  trt  <- as.character(exp_combos$Treatment[i])

  # Include: both genotypes × (Vehicle + this treatment) at all timepoints
  subset_meta <- meta_exp %>%
    filter(Treatment %in% c('Vehicle', 'NoTreatment', trt)) %>%
    filter(Genotype %in% c('WT', geno))

  p_exp <- ggplot(subset_meta,
                  aes(x = Treatment2, y = NumUnqDedupedReads)) +
    geom_hline(yintercept = threshold, linetype = 'dashed',
               colour = 'darkred', linewidth = 0.8) +
    geom_boxplot(outlier.shape = NA, alpha = 0.5, fill = 'gray95') +
    geom_jitter(aes(colour = NumUnqDedupedReads < threshold),
                width = 0.2, size = 2, alpha = 0.8) +
    facet_grid(Genotype ~ Timepoint, scales = 'free_x', space = 'free_x') +
    scale_colour_manual(values = c('TRUE' = 'red', 'FALSE' = 'black'),
                        name   = 'Below Threshold',
                        labels = c('TRUE' = '< 2M', 'FALSE' = '>= 2M')) +
    scale_y_continuous(labels = scales::unit_format(unit = 'M', scale = 1e-6)) +
    labs(title    = paste0('Sequencing Effort — ', geno, ' + ', trt),
         subtitle = 'Dashed line and red points: < 2 million unique deduped reads',
         x = 'Treatment', y = 'Unique Deduped Reads (Millions)') +
    theme_bw() +
    theme(axis.text.x  = element_text(angle = 45, hjust = 1, size = 8),
          strip.text.y = element_text(angle = 0, size = 7, face = 'bold'),
          legend.position = 'bottom')

  fname <- paste0('260318_UMIDeDuped_ReadCounts_', geno, '_', trt, '.png')
  ggsave(file.path(FIG_DIR, fname), p_exp, width = 12, height = 6)
  log_msg('Saved: qc_figures/', fname)
}

log_msg('Script complete')
