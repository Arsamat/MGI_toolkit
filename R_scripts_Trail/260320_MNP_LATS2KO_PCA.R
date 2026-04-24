library(tidyverse)
library(edgeR)
library(RColorBrewer)

## ── Setup ──────────────────────────────────────────────────────────────────────
BASE_DIR <- '/scratch/rmlab/rmlab_shared3/alferreiro/BRB_Seq'
ANA_DIR  <- file.path(BASE_DIR, '260320_MNP_WT_LATS2KO_Rot_TRULI_Timecourse_Analyses')
DATA_DIR <- file.path(ANA_DIR, 'data')
FIG_DIR  <- file.path(ANA_DIR, 'pca_figures')
LOG_DIR  <- file.path(ANA_DIR, 'logs')

for (d in c(FIG_DIR, LOG_DIR)) dir.create(d, showWarnings = FALSE, recursive = TRUE)
setwd(BASE_DIR)

log_file <- file.path(LOG_DIR, '260320_MNP_LATS2KO_PCA.log')
cat('', file = log_file)
log_msg <- function(...) {
  msg <- paste0('[', format(Sys.time(), '%Y-%m-%d %H:%M:%S'), '] ', ...)
  message(msg); cat(msg, '\n', file = log_file, append = TRUE)
}
log_msg('Script started')

dge <- readRDS(file.path(DATA_DIR, '260320_MNP_LATS2KO_DGElist.RDS'))
log_msg('Loaded DGEList: ', nrow(dge), ' genes x ', ncol(dge), ' samples')

dge$samples <- dge$samples %>% mutate(
  Genotype  = factor(Genotype,  levels = c('WT', 'LATS2KO')),
  Treatment = factor(Treatment, levels = c('Vehicle', 'TRULI', 'Rotenone', 'TRULI + Rotenone')),
  Timepoint = factor(as.character(Timepoint), levels = c('6', '16', '24', '48', '96'))
)

paletteGenotype  <- c('WT' = '#555555', 'LATS2KO' = '#3182BD')
paletteTreatment <- c('Vehicle' = '#969696', 'TRULI' = '#74C476',
                      'Rotenone' = '#FC8D59', 'TRULI + Rotenone' = '#31A354')
paletteTimepoint <- c('6' = '#FEE391', '16' = '#FEB24C', '24' = '#FC4E2A', '48' = '#800026', '96' = '#4D004B')
shapesGenotype   <- c('WT' = 22, 'LATS2KO' = 23)
shapesTreatment  <- c('Vehicle' = 22, 'TRULI' = 23, 'Rotenone' = 24, 'TRULI + Rotenone' = 25)
shapesTimepoint  <- c('6' = 21, '16' = 22, '24' = 23, '48' = 24, '96' = 25)

pca_axis_label <- function(pca_result, axis) {
  pct <- round(100 * pca_result$sdev[axis]^2 / sum(pca_result$sdev^2), 1)
  paste0('PC', axis, ' (', pct, '% variance)')
}
pca_theme <- function() {
  theme_bw() %+replace% theme(
    panel.grid.major = element_line(colour = 'grey92'), panel.grid.minor = element_blank(),
    axis.title = element_text(size = 12), axis.text = element_text(size = 10),
    legend.title = element_text(size = 11, face = 'bold'), legend.text = element_text(size = 9),
    plot.title = element_text(size = 14, face = 'bold'), plot.subtitle = element_text(size = 11)
  )
}
calculate_PCA <- function(dge_obj, n = 2000) {
  mat <- edgeR::cpm(dge_obj, log = TRUE, prior.count = 1, normalized.lib.sizes = TRUE)
  prcomp(t(mat[head(order(apply(mat, 1, var), decreasing = TRUE), n), ]), scale = TRUE)
}

top_number <- 2000
log_msg('Calculating PCA...')
pca_res <- calculate_PCA(dge, top_number)
pca_dat <- data.frame(
  PC1 = pca_res$x[,1], PC2 = pca_res$x[,2],
  Genotype           = as.character(dge$samples$Genotype),
  Treatment          = as.character(dge$samples$Treatment),
  Timepoint          = dge$samples$Timepoint,
  NumUnqDedupedReads = dge$samples$NumUnqDedupedReads,
  Lib_size           = dge$samples$lib.size,
  Replicate          = dge$samples$Replicate,
  Sample             = colnames(dge), stringsAsFactors = FALSE
)
pca_dat <- pca_dat[pca_dat$Lib_size > 500000, ]

base_layers <- list(
  geom_hline(yintercept = 0, linetype = 'dashed', colour = 'black', linewidth = 0.75),
  geom_vline(xintercept = 0, linetype = 'dashed', colour = 'black', linewidth = 0.75),
  xlab(pca_axis_label(pca_res, 1)), ylab(pca_axis_label(pca_res, 2)),
  ggtitle('RNA-Seq', subtitle = paste0(top_number, ' most variable genes')),
  pca_theme(), coord_fixed()
)

save_pdf <- function(p, name) {
  fpath <- file.path(FIG_DIR, paste0('260320_MNP_LATS2KO_PCA_', name, '.pdf'))
  pdf(fpath, width = 8, height = 8, useDingbats = FALSE); print(p); dev.off()
  log_msg('Saved: pca_figures/', basename(fpath))
}

p1 <- ggplot(pca_dat, aes(PC1, PC2)) +
  geom_point(aes(fill = Treatment, colour = Treatment, shape = Genotype), size = 4, alpha = 0.8) +
  scale_fill_manual(values = paletteTreatment) + scale_colour_manual(values = paletteTreatment) +
  scale_shape_manual(values = shapesGenotype) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) + base_layers
save_pdf(p1, 'Genotype_by_Treatment')

p2 <- ggplot(pca_dat, aes(PC1, PC2)) +
  geom_point(aes(fill = Genotype, colour = Genotype, shape = Treatment), size = 4, alpha = 0.8) +
  scale_fill_manual(values = paletteGenotype) + scale_colour_manual(values = paletteGenotype) +
  scale_shape_manual(values = shapesTreatment) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) + base_layers
save_pdf(p2, 'Treatment_by_Genotype')

p3 <- ggplot(pca_dat, aes(PC1, PC2)) +
  geom_point(aes(fill = Timepoint, colour = Timepoint, shape = Genotype), size = 4, alpha = 0.8) +
  scale_fill_manual(values = paletteTimepoint) + scale_colour_manual(values = paletteTimepoint) +
  scale_shape_manual(values = shapesGenotype) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) + base_layers
save_pdf(p3, 'Genotype_by_Timepoint')

p4 <- ggplot(pca_dat %>% arrange(desc(NumUnqDedupedReads)) %>% mutate(Low = NumUnqDedupedReads < 2e6),
             aes(PC1, PC2)) +
  geom_point(aes(fill = NumUnqDedupedReads, colour = Low), shape = 21, size = 4, alpha = 0.8, stroke = 1.2) +
  scale_fill_gradient(low = '#fee8c8', high = '#7f0000', name = 'Unique Deduped\nReads',
                      labels = scales::unit_format(unit = 'M', scale = 1e-6)) +
  scale_colour_manual(values = c('TRUE' = 'red', 'FALSE' = 'grey30'),
                      name = '< 2M reads', labels = c('TRUE' = 'Yes', 'FALSE' = 'No')) +
  base_layers
save_pdf(p4, 'NumUnqDedupedReads')

log_msg('Script complete')
