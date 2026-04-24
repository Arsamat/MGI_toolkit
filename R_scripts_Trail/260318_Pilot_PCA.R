library(tidyverse)
library(edgeR)
library(limma)
library(RColorBrewer)

## ── Setup ──────────────────────────────────────────────────────────────────────
BASE_DIR <- '/scratch/rmlab/rmlab_shared3/alferreiro/BRB_Seq'
ANA_DIR  <- file.path(BASE_DIR, '260318_Pilot_Analyses')
DATA_DIR <- file.path(ANA_DIR, 'data')
FIG_DIR  <- file.path(ANA_DIR, 'pca_figures')
LOG_DIR  <- file.path(ANA_DIR, 'logs')

dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR,  showWarnings = FALSE, recursive = TRUE)
dir.create(LOG_DIR,  showWarnings = FALSE, recursive = TRUE)

setwd(BASE_DIR)

log_file <- file.path(LOG_DIR, '260318_Pilot_PCA.log')
cat('', file = log_file)
log_msg <- function(...) {
  msg <- paste0('[', format(Sys.time(), '%Y-%m-%d %H:%M:%S'), '] ', ...)
  message(msg)
  cat(msg, '\n', file = log_file, append = TRUE)
}
log_msg('Script started')
log_msg('ANA_DIR : ', ANA_DIR)
log_msg('FIG_DIR : ', FIG_DIR)

## ── Load DGEList ───────────────────────────────────────────────────────────────
dge <- readRDS(file.path(DATA_DIR, '260318_Pilot_DGElist.RDS'))
log_msg('Loaded DGEList: ', nrow(dge), ' genes x ', ncol(dge), ' samples')

## ── Factor levels ──────────────────────────────────────────────────────────────
dge$samples <- dge$samples %>%
  mutate(
    Genotype = factor(Genotype, levels = c('WT', 'TDP43')),
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
    Run = factor(Run, levels = c(
      'Pilot_Spinoff', 'Pilot_Spinoff2', 'Pilot_Batch1', 'Pilot_Batch2'
    ))
  )

## ── Color palettes and shapes ──────────────────────────────────────────────────
paletteGenotype <- c(
  'WT'    = '#555555',
  'TDP43' = '#E41A1C'
)

paletteTreatment2 <- c(
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

shapesTreatment2 <- c(
  'NoTreatment'       = 21,
  'Vehicle'           = 22,
  'Rotenone_EC10'     = 23,
  'Rotenone_EC50'     = 24,
  'Rotenone_EC90'     = 25,
  'Thapsigargin_EC10' = 3,
  'Thapsigargin_EC50' = 4,
  'Thapsigargin_EC90' = 8
)

shapesTimepoint <- c(
  '0'  = 21,
  '6'  = 22,
  '24' = 23,
  '72' = 24,
  '96' = 25
)

## ── Helper functions ───────────────────────────────────────────────────────────
pca_axis_label <- function(pca_result, pca_axis) {
  pct <- round(100 * pca_result$sdev[pca_axis]^2 / sum(pca_result$sdev^2), 1)
  paste0('PC', pca_axis, ' (', pct, '% variance)')
}

pca_personal_theme <- function() {
  theme_bw() %+replace% theme(
    panel.grid.major = element_line(colour = 'grey92'),
    panel.grid.minor = element_blank(),
    axis.title       = element_text(size = 12),
    axis.text        = element_text(size = 10),
    legend.title     = element_text(size = 11, face = 'bold'),
    legend.text      = element_text(size = 9),
    plot.title       = element_text(size = 14, face = 'bold'),
    plot.subtitle    = element_text(size = 11)
  )
}

calculate_PCA <- function(myDGEList, numberTopGenes = 2000) {
  mat       <- edgeR::cpm(myDGEList, log = TRUE, prior.count = 1, normalized.lib.sizes = TRUE)
  top_genes <- head(order(apply(mat, 1, var), decreasing = TRUE), n = numberTopGenes)
  prcomp(t(mat[top_genes, ]), scale = TRUE)
}

calculate_PCA_postcorrect <- function(mat, numberTopGenes = 2000) {
  top_genes <- head(order(apply(mat, 1, var), decreasing = TRUE), n = numberTopGenes)
  prcomp(t(mat[top_genes, ]), scale = TRUE)
}

make_pca_data <- function(pca_res, dge_obj) {
  df <- data.frame(
    PC1                = pca_res$x[, 1],
    PC2                = pca_res$x[, 2],
    PC3                = pca_res$x[, 3],
    Lib_size           = dge_obj$samples$lib.size,
    NumUnqDedupedReads = dge_obj$samples$NumUnqDedupedReads,
    Genotype           = as.character(dge_obj$samples$Genotype),
    Treatment2         = as.character(dge_obj$samples$Treatment2),
    Timepoint          = dge_obj$samples$Timepoint,
    Run                = as.character(dge_obj$samples$Run),
    Replicate          = dge_obj$samples$Replicate,
    Sample             = colnames(dge_obj),
    stringsAsFactors   = FALSE
  )
  df[df$Lib_size > 500000, ]
}

save_pca_plots <- function(pca_res, pca_dat, suffix, top_number) {
  subtitle <- paste0(top_number, ' most variable genes')

  base_layers <- list(
    geom_hline(yintercept = 0, linetype = 'dashed', colour = 'black', linewidth = 0.75),
    geom_vline(xintercept = 0, linetype = 'dashed', colour = 'black', linewidth = 0.75),
    xlab(pca_axis_label(pca_res, 1)),
    ylab(pca_axis_label(pca_res, 2)),
    ggtitle('RNA-Seq', subtitle = subtitle),
    pca_personal_theme(),
    coord_fixed()
  )

  # 1) Shape by Genotype, color/fill by Run
  p1 <- ggplot(pca_dat, aes(PC1, PC2, text = Sample)) +
    geom_point(aes(fill = Run, color = Run, shape = Genotype), size = 4, alpha = 0.8) +
    scale_shape_manual(values = shapesGenotype) +
    scale_color_manual(values = paletteRun) +
    scale_fill_manual(values = paletteRun) +
    guides(fill = guide_legend(override.aes = list(shape = 21))) +
    base_layers

  # 2) Shape by Genotype, color/fill by Treatment2
  p2 <- ggplot(pca_dat, aes(PC1, PC2, text = Sample)) +
    geom_point(aes(fill = Treatment2, color = Treatment2, shape = Genotype), size = 4, alpha = 0.8) +
    scale_fill_manual(values = paletteTreatment2) +
    scale_color_manual(values = paletteTreatment2) +
    scale_shape_manual(values = shapesGenotype) +
    guides(fill = guide_legend(override.aes = list(shape = 21))) +
    base_layers

  # 3) Shape by Timepoint, color/fill by Treatment2
  p3 <- ggplot(pca_dat, aes(PC1, PC2, text = Sample)) +
    geom_point(aes(fill = Treatment2, color = Treatment2, shape = Timepoint), size = 4, alpha = 0.8) +
    scale_fill_manual(values = paletteTreatment2) +
    scale_color_manual(values = paletteTreatment2) +
    scale_shape_manual(values = shapesTimepoint) +
    guides(fill = guide_legend(override.aes = list(shape = 21))) +
    base_layers

  # 4) Shape by Genotype, color/fill by Timepoint
  p4 <- ggplot(pca_dat, aes(PC1, PC2, text = Sample)) +
    geom_point(aes(fill = Timepoint, color = Timepoint, shape = Genotype), size = 4, alpha = 0.8) +
    scale_fill_manual(values = paletteTimepoint) +
    scale_color_manual(values = paletteTimepoint) +
    scale_shape_manual(values = shapesGenotype) +
    guides(fill = guide_legend(override.aes = list(shape = 21))) +
    base_layers

  # 5) Shape by Treatment2, color/fill by Genotype
  p5 <- ggplot(pca_dat, aes(PC1, PC2, text = Sample)) +
    geom_point(aes(fill = Genotype, color = Genotype, shape = Treatment2), size = 4, alpha = 0.8) +
    scale_fill_manual(values = paletteGenotype) +
    scale_color_manual(values = paletteGenotype) +
    scale_shape_manual(values = shapesTreatment2) +
    guides(fill = guide_legend(override.aes = list(shape = 21))) +
    base_layers

  # 6) WT Vehicle samples highlighted by Run; all others in grey
  pca_rh <- pca_dat %>%
    mutate(RunHighlight = if_else(
      Genotype == 'WT' & Treatment2 == 'Vehicle',
      Run, 'OtherSamples'
    ))
  paletteRunHighlight <- c('OtherSamples' = '#D3D3D3', paletteRun)
  pca_rh <- bind_rows(
    pca_rh %>% filter(RunHighlight == 'OtherSamples'),
    pca_rh %>% filter(RunHighlight != 'OtherSamples')
  )
  p6 <- ggplot(pca_rh, aes(PC1, PC2, text = Sample)) +
    geom_point(aes(fill = RunHighlight, color = RunHighlight, shape = Genotype),
               size = 4, alpha = 0.8) +
    scale_fill_manual(values = paletteRunHighlight) +
    scale_color_manual(values = paletteRunHighlight) +
    scale_shape_manual(values = shapesGenotype) +
    guides(fill = guide_legend(override.aes = list(shape = 21))) +
    labs(title    = 'RNA-Seq',
         subtitle = paste0(top_number, ' most variable genes — WT Vehicle highlighted by Run')) +
    base_layers

  # 7) Fixed shape (21), fill by NumUnqDedupedReads (continuous gradient)
  pca_cov <- bind_rows(
    pca_dat %>% filter(NumUnqDedupedReads >= 2e6),
    pca_dat %>% filter(NumUnqDedupedReads <  2e6)
  ) %>% mutate(LowCoverage = NumUnqDedupedReads < 2e6)
  p7 <- ggplot(pca_cov, aes(PC1, PC2, text = Sample)) +
    geom_point(aes(fill = NumUnqDedupedReads, color = LowCoverage),
               shape = 21, size = 4, alpha = 0.8, stroke = 1.2) +
    scale_fill_gradient(low = '#fee8c8', high = '#7f0000',
                        name   = 'Unique Deduped\nReads',
                        labels = scales::unit_format(unit = 'M', scale = 1e-6)) +
    scale_color_manual(values = c('TRUE' = 'red', 'FALSE' = 'grey30'),
                       name   = '< 2M reads',
                       labels = c('TRUE' = 'Yes', 'FALSE' = 'No')) +
    base_layers

  plots <- list(
    Genotype_by_Run        = p1,
    Genotype_by_Treatment2 = p2,
    Timepoint_by_Treatment2 = p3,
    Genotype_by_Timepoint  = p4,
    Treatment2_by_Genotype = p5,
    WTVehicle_by_Run       = p6,
    NumUnqDedupedReads     = p7
  )

  for (nm in names(plots)) {
    fpath <- file.path(FIG_DIR, paste0('260318_PCA_', nm, suffix, '.pdf'))
    pdf(fpath, width = 8, height = 8, onefile = TRUE, useDingbats = FALSE)
    print(plots[[nm]])
    invisible(dev.off())
    log_msg('Saved: pca_figures/', basename(fpath))
  }
}

## ── run_pca_suite: pre + post-correction PCA + WT Vehicle Run distance ─────────
run_pca_suite <- function(dge_obj, lcpm_rbe_full, file_suffix, top_number = 2000) {
  log_msg('=== PCA suite [', file_suffix, '] — ', ncol(dge_obj), ' samples ===')

  ## ── Pre-correction PCA ──────────────────────────────────────────────────────
  log_msg('Calculating pre-correction PCA...')
  pca_pre <- calculate_PCA(dge_obj, numberTopGenes = top_number)
  pca_dat <- make_pca_data(pca_pre, dge_obj)
  save_pca_plots(pca_pre, pca_dat,
                 suffix = file_suffix, top_number = top_number)

  ## ── Post-correction PCA ─────────────────────────────────────────────────────
  log_msg('Calculating batch-corrected PCA...')
  lcpm_rbe <- lcpm_rbe_full[, colnames(dge_obj), drop = FALSE]
  pca_bc   <- calculate_PCA_postcorrect(lcpm_rbe, numberTopGenes = top_number)
  pca_dat2 <- data.frame(
    PC1                = pca_bc$x[, 1],
    PC2                = pca_bc$x[, 2],
    PC3                = pca_bc$x[, 3],
    Lib_size           = dge_obj$samples$lib.size,
    NumUnqDedupedReads = dge_obj$samples$NumUnqDedupedReads,
    Genotype           = as.character(dge_obj$samples$Genotype),
    Treatment2         = as.character(dge_obj$samples$Treatment2),
    Timepoint          = dge_obj$samples$Timepoint,
    Run                = as.character(dge_obj$samples$Run),
    Replicate          = dge_obj$samples$Replicate,
    Sample             = colnames(dge_obj),
    stringsAsFactors   = FALSE
  )
  pca_dat2 <- pca_dat2[pca_dat2$Lib_size > 500000, ]

  save_pca_plots(pca_bc, pca_dat2,
                 suffix = paste0('_BATCHCORRECTED', file_suffix),
                 top_number = top_number)

  ## ── WT Vehicle pairwise distance analysis (pre-correction) ──────────────────
  log_msg('Running WT Vehicle pairwise distance analysis by Run...')
  wt_veh <- pca_dat %>%
    filter(Genotype == 'WT', Treatment2 == 'Vehicle')
  log_msg('WT Vehicle subset: ', nrow(wt_veh), ' samples')

  if (nrow(wt_veh) < 2) {
    log_msg('Skipping distance analysis: fewer than 2 WT Vehicle samples')
    return(invisible(NULL))
  }

  n_pc <- min(5, ncol(pca_pre$x))
  pc5  <- pca_pre$x[wt_veh$Sample, seq_len(n_pc), drop = FALSE]
  D    <- as.matrix(dist(pc5, method = 'euclidean'))
  idx  <- which(upper.tri(D), arr.ind = TRUE)

  same_run <- wt_veh$Run[idx[, 1]] == wt_veh$Run[idx[, 2]]
  dvals    <- D[idx]

  dist_long <- data.frame(
    Distance = dvals,
    Group    = factor(
      ifelse(same_run, 'Within Run', 'Between Run'),
      levels = c('Within Run', 'Between Run')
    )
  )
  log_msg('Pairwise distances: ', nrow(idx), ' pairs')

  p_dist <- ggplot(dist_long, aes(x = Group, y = Distance)) +
    geom_boxplot(outlier.shape = NA, fill = 'gray90', width = 0.55) +
    geom_jitter(width = 0.2, size = 0.3, alpha = 0.08, colour = 'steelblue') +
    labs(title    = 'Pairwise Euclidean Distances (top 5 PCs, pre-correction)',
         subtitle = paste0('WT Vehicle samples only; n = ',
                           nrow(wt_veh), ' samples, ', nrow(idx), ' pairs'),
         x = NULL, y = 'Euclidean Distance') +
    theme_bw() +
    theme(axis.text.x = element_text(size = 11))

  png_name <- paste0('260318_WTVehicle_Run_PairwiseDistances', file_suffix, '.png')
  ggsave(file.path(FIG_DIR, png_name), p_dist, width = 5, height = 5)
  log_msg('Saved: pca_figures/', png_name)

  t_result <- t.test(Distance ~ Group, data = dist_long)
  txt_name <- paste0('260318_WTVehicle_Run_DistanceTtest', file_suffix, '.txt')
  sink(file.path(DATA_DIR, txt_name))
  cat('=== t-test: pairwise Euclidean distances (top 5 PCs) ~ Run group ===\n')
  cat('Samples: WT Vehicle; n =', nrow(wt_veh), '\n')
  cat('Note: pairs share samples so are not fully independent.\n\n')
  print(t_result)
  sink()
  log_msg('Saved: data/', txt_name)

  invisible(NULL)
}

## ── Batch correction ─────────────────────────────────────────────────────────
top_number <- 2000
log_msg('Running batch correction on full dataset...')
lcpm_full     <- edgeR::cpm(dge, log = TRUE, prior.count = 1, normalized.lib.sizes = TRUE)
mod_full      <- model.matrix(~ Genotype + Treatment2 + Timepoint,
                              data = droplevels(dge$samples))
lcpm_rbe_full <- limma::removeBatchEffect(lcpm_full, batch = dge$samples$Run,
                                          design = mod_full)

write.csv(lcpm_rbe_full,
          file.path(DATA_DIR, '260318_Pilot_Dedup_Counts_BatchCorrected.csv'))
log_msg('Saved batch-corrected logCPM: data/260318_Pilot_Dedup_Counts_BatchCorrected.csv')

## ── Run PCA suite ─────────────────────────────────────────────────────────────
run_pca_suite(dge, lcpm_rbe_full, file_suffix = '', top_number = top_number)

log_msg('Script complete')
