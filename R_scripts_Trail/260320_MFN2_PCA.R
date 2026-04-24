library(tidyverse)
library(edgeR)
library(RColorBrewer)

## ── Setup ──────────────────────────────────────────────────────────────────────
BASE_DIR <- '/scratch/rmlab/rmlab_shared3/alferreiro/BRB_Seq'
ANA_DIR  <- file.path(BASE_DIR, '260320_MFN2_MNP_Rotenone_Analyses')
DATA_DIR <- file.path(ANA_DIR, 'data')
FIG_DIR  <- file.path(ANA_DIR, 'pca_figures')
LOG_DIR  <- file.path(ANA_DIR, 'logs')

for (d in c(FIG_DIR, LOG_DIR)) dir.create(d, showWarnings = FALSE, recursive = TRUE)
setwd(BASE_DIR)

log_file <- file.path(LOG_DIR, '260320_MFN2_PCA.log')
cat('', file = log_file)
log_msg <- function(...) {
  msg <- paste0('[', format(Sys.time(), '%Y-%m-%d %H:%M:%S'), '] ', ...)
  message(msg); cat(msg, '\n', file = log_file, append = TRUE)
}
log_msg('Script started')

## ── Load DGEList ───────────────────────────────────────────────────────────────
dge <- readRDS(file.path(DATA_DIR, '260320_MFN2_DGElist.RDS'))
log_msg('Loaded DGEList: ', nrow(dge), ' genes x ', ncol(dge), ' samples')

## ── Factor levels ──────────────────────────────────────────────────────────────
dge$samples <- dge$samples %>% mutate(
  Genotype   = factor(Genotype,   levels = c('WT', 'MFN2-Het', 'MFN2')),
  Treatment2 = factor(Treatment2, levels = c('Vehicle', 'Rotenone_0.4uM', 'Rotenone_2uM'))
)

## ── Palettes and shapes ────────────────────────────────────────────────────────
paletteGenotype   <- c('WT' = '#555555', 'MFN2-Het' = '#FC8D59', 'MFN2' = '#D7301F')
paletteTreatment2 <- c('Vehicle' = '#969696', 'Rotenone_0.4uM' = '#FDCC8A', 'Rotenone_2uM' = '#D7301F')
shapesGenotype    <- c('WT' = 22, 'MFN2-Het' = 23, 'MFN2' = 24)
shapesTreatment2  <- c('Vehicle' = 22, 'Rotenone_0.4uM' = 23, 'Rotenone_2uM' = 24)

## ── Helpers ────────────────────────────────────────────────────────────────────
pca_axis_label <- function(pca_result, axis) {
  pct <- round(100 * pca_result$sdev[axis]^2 / sum(pca_result$sdev^2), 1)
  paste0('PC', axis, ' (', pct, '% variance)')
}

pca_theme <- function() {
  theme_bw() %+replace% theme(
    panel.grid.major = element_line(colour = 'grey92'),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 12), axis.text = element_text(size = 10),
    legend.title = element_text(size = 11, face = 'bold'),
    legend.text = element_text(size = 9),
    plot.title = element_text(size = 14, face = 'bold'),
    plot.subtitle = element_text(size = 11)
  )
}

calculate_PCA <- function(dge_obj, n = 2000) {
  mat       <- edgeR::cpm(dge_obj, log = TRUE, prior.count = 1, normalized.lib.sizes = TRUE)
  top_genes <- head(order(apply(mat, 1, var), decreasing = TRUE), n)
  prcomp(t(mat[top_genes, ]), scale = TRUE)
}

## ── Run PCA ────────────────────────────────────────────────────────────────────
top_number <- 2000
log_msg('Calculating PCA on ', top_number, ' most variable genes...')
pca_res <- calculate_PCA(dge, top_number)

pca_dat <- data.frame(
  PC1                = pca_res$x[, 1],
  PC2                = pca_res$x[, 2],
  Genotype           = as.character(dge$samples$Genotype),
  Treatment2         = as.character(dge$samples$Treatment2),
  NumUnqDedupedReads = dge$samples$NumUnqDedupedReads,
  Lib_size           = dge$samples$lib.size,
  Replicate          = dge$samples$Replicate,
  Sample             = colnames(dge),
  stringsAsFactors   = FALSE
)
pca_dat <- pca_dat[pca_dat$Lib_size > 500000, ]

base_layers <- list(
  geom_hline(yintercept = 0, linetype = 'dashed', colour = 'black', linewidth = 0.75),
  geom_vline(xintercept = 0, linetype = 'dashed', colour = 'black', linewidth = 0.75),
  xlab(pca_axis_label(pca_res, 1)), ylab(pca_axis_label(pca_res, 2)),
  ggtitle('RNA-Seq', subtitle = paste0(top_number, ' most variable genes')),
  pca_theme(), coord_fixed()
)

## Plot 1: Genotype_by_Treatment2
p1 <- ggplot(pca_dat, aes(PC1, PC2)) +
  geom_point(aes(fill = Treatment2, colour = Treatment2, shape = Genotype), size = 4, alpha = 0.8) +
  scale_fill_manual(values = paletteTreatment2) +
  scale_colour_manual(values = paletteTreatment2) +
  scale_shape_manual(values = shapesGenotype) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) + base_layers
pdf(file.path(FIG_DIR, '260320_MFN2_PCA_Genotype_by_Treatment2.pdf'), width = 8, height = 8, useDingbats = FALSE)
print(p1); dev.off()
log_msg('Saved: pca_figures/260320_MFN2_PCA_Genotype_by_Treatment2.pdf')

## Plot 2: Treatment2_by_Genotype
p2 <- ggplot(pca_dat, aes(PC1, PC2)) +
  geom_point(aes(fill = Genotype, colour = Genotype, shape = Treatment2), size = 4, alpha = 0.8) +
  scale_fill_manual(values = paletteGenotype) +
  scale_colour_manual(values = paletteGenotype) +
  scale_shape_manual(values = shapesTreatment2) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) + base_layers
pdf(file.path(FIG_DIR, '260320_MFN2_PCA_Treatment2_by_Genotype.pdf'), width = 8, height = 8, useDingbats = FALSE)
print(p2); dev.off()
log_msg('Saved: pca_figures/260320_MFN2_PCA_Treatment2_by_Genotype.pdf')

## Plot 3: NumUnqDedupedReads
pca_cov <- pca_dat %>%
  arrange(desc(NumUnqDedupedReads)) %>%
  mutate(LowCoverage = NumUnqDedupedReads < 2e6)
p3 <- ggplot(pca_cov, aes(PC1, PC2)) +
  geom_point(aes(fill = NumUnqDedupedReads, colour = LowCoverage),
             shape = 21, size = 4, alpha = 0.8, stroke = 1.2) +
  scale_fill_gradient(low = '#fee8c8', high = '#7f0000', name = 'Unique Deduped\nReads',
                      labels = scales::unit_format(unit = 'M', scale = 1e-6)) +
  scale_colour_manual(values = c('TRUE' = 'red', 'FALSE' = 'grey30'),
                      name = '< 2M reads', labels = c('TRUE' = 'Yes', 'FALSE' = 'No')) +
  base_layers
pdf(file.path(FIG_DIR, '260320_MFN2_PCA_NumUnqDedupedReads.pdf'), width = 8, height = 8, useDingbats = FALSE)
print(p3); dev.off()
log_msg('Saved: pca_figures/260320_MFN2_PCA_NumUnqDedupedReads.pdf')

log_msg('Script complete')
