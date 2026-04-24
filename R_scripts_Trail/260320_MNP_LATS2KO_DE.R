library(tidyverse)
library(edgeR)
library(limma)
library(ggrepel)

BASE_DIR <- '/scratch/rmlab/rmlab_shared3/alferreiro/BRB_Seq'
ANA_DIR  <- file.path(BASE_DIR, '260320_MNP_WT_LATS2KO_Rot_TRULI_Timecourse_Analyses')
DATA_DIR <- file.path(ANA_DIR, 'data')
LOG_DIR  <- file.path(ANA_DIR, 'logs')
OUT_DE   <- file.path(ANA_DIR, 'DE_analysis')
OUT_G0   <- file.path(OUT_DE,  'Group0')
OUT_G1   <- file.path(OUT_DE,  'Group1')

for (d in c(LOG_DIR, OUT_G0, OUT_G1)) dir.create(d, showWarnings = FALSE, recursive = TRUE)
setwd(BASE_DIR)

PADJ_THRESH <- 0.05; LFC_THRESH <- 0.6

log_file <- file.path(LOG_DIR, '260320_MNP_LATS2KO_DE_Analysis.log')
cat('', file = log_file)
log_msg <- function(...) {
  msg <- paste0('[', format(Sys.time(), '%Y-%m-%d %H:%M:%S'), '] ', ...)
  message(msg); cat(msg, '\n', file = log_file, append = TRUE)
}
log_msg('Script started')

dge <- readRDS(file.path(DATA_DIR, '260320_MNP_LATS2KO_DGElist.RDS'))
log_msg('Loaded DGEList: ', nrow(dge), ' genes x ', ncol(dge), ' samples')

## ── Shared helpers ─────────────────────────────────────────────────────────────
safe_str  <- function(x) gsub('[^A-Za-z0-9_]', '_', x)
count_sig <- function(res) sum(res$adj.P.Val < PADJ_THRESH & abs(res$logFC) > LFC_THRESH, na.rm = TRUE)

make_volcano <- function(res, coef_name, title_str, fname, out_dir) {
  res <- res %>% mutate(sig = adj.P.Val < PADJ_THRESH & abs(logFC) > LFC_THRESH, logP = -log10(P.Value))
  top <- res %>% filter(sig) %>% slice_min(adj.P.Val, n = 50)
  p <- ggplot(res, aes(logFC, logP)) +
    geom_point(aes(colour = sig), alpha = 0.3, size = 1) +
    geom_text_repel(data = top, aes(label = gene_symbol), size = 2.5, max.overlaps = 20, segment.alpha = 0.5) +
    scale_colour_manual(values = c('TRUE' = 'firebrick3', 'FALSE' = 'grey70')) +
    geom_vline(xintercept = c(-LFC_THRESH, LFC_THRESH), linetype = 'dotted') +
    geom_hline(yintercept = -log10(PADJ_THRESH), linetype = 'dotted') +
    labs(title = title_str, x = paste0('logFC (', coef_name, ')'), y = '-log10(p)') +
    theme_bw() + theme(legend.position = 'none')
  ggsave(file.path(out_dir, fname), p, width = 8, height = 6)
  log_msg('    Saved: ', fname)
}

## ── Shared QC filter ───────────────────────────────────────────────────────────
dge_qc <- dge[, dge$samples$NumUnqDedupedReads >= 2e6]
log_msg('After QC filter: ', ncol(dge_qc), ' samples')

## ── GROUP0: LATS2KO vs WT in Vehicle-treated samples ──────────────────────────
log_msg('=== GROUP0: ~ Genotype (Vehicle only) ===')

dge_g0 <- dge_qc[, dge_qc$samples$Treatment == 'Vehicle']
dge_g0$samples <- droplevels(dge_g0$samples)
dge_g0$samples$Genotype <- factor(dge_g0$samples$Genotype, levels = c('WT', 'LATS2KO'))
dge_g0$samples$Genotype <- relevel(dge_g0$samples$Genotype, ref = 'WT')
log_msg('  Samples: ', ncol(dge_g0))
log_msg('  Genotypes: ', paste(levels(dge_g0$samples$Genotype), collapse = ', '))
log_msg('  Timepoints: ', paste(sort(unique(dge_g0$samples$Timepoint)), collapse = ', '))

sample_rows <- list()
sample_rows[['GROUP0_Vehicle']] <- data.frame(
  group        = 'GROUP0',
  context      = 'Treatment=Vehicle',
  model_design = '~ Genotype',
  n_samples    = ncol(dge_g0),
  sample_names = paste(dge_g0$samples$SampleName, collapse = ', '),
  stringsAsFactors = FALSE
)

design_g0 <- model.matrix(~ Genotype, data = dge_g0$samples)
log_msg('  Design matrix: ', nrow(design_g0), ' x ', ncol(design_g0))

fit_g0 <- tryCatch({
  v <- voom(dge_g0, design_g0, plot = FALSE); eBayes(lmFit(v, design_g0))
}, error = function(e) stop('GROUP0 voom/lmFit failed: ', conditionMessage(e)))

geno_coefs_g0 <- grep('^Genotype', colnames(fit_g0$coefficients), value = TRUE)
log_msg('  Genotype terms: ', length(geno_coefs_g0))

summary_rows_g0 <- list()
nsig_geno_g0    <- list()

for (coef in geno_coefs_g0) {
  if (all(is.na(fit_g0$coefficients[, coef]))) { log_msg('  SKIP (all NA): ', coef); next }
  log_msg('  [Genotype] ', coef)
  res <- tryCatch(topTable(fit_g0, coef = coef, n = Inf, sort.by = 'P'),
                  error = function(e) { log_msg('  SKIP topTable: ', conditionMessage(e)); NULL })
  if (is.null(res)) next
  geno_val <- sub('^Genotype', '', coef)
  n_sig    <- count_sig(res)
  log_msg('    Sig: ', n_sig)
  res$model_design <- '~ Genotype'
  res$term         <- 'Genotype'
  res$genotype     <- geno_val
  safe_coef        <- safe_str(coef)
  write.csv(res, file.path(OUT_G0, paste0('260320_MNP_LATS2KO_DE_GROUP0_Vehicle_', safe_coef, '.csv')),
            row.names = FALSE)
  make_volcano(res, coef,
               paste0('GROUP0 | Vehicle | Genotype | ', coef, '\n(n sig = ', n_sig, ')'),
               paste0('260320_MNP_LATS2KO_Volcano_GROUP0_Vehicle_', safe_coef, '.png'),
               OUT_G0)
  nsig_geno_g0[[geno_val]] <- n_sig
}

for (geno_val in names(nsig_geno_g0)) {
  summary_rows_g0[[geno_val]] <- data.frame(
    group            = 'GROUP0',
    context          = 'Treatment=Vehicle',
    model_design     = '~ Genotype',
    genotype_value   = geno_val,
    genotype_n_sig   = nsig_geno_g0[[geno_val]],
    stringsAsFactors = FALSE
  )
}

## ── GROUP1: Treatment × Timepoint per Genotype (Timepoint as factor) ──────────
log_msg('=== GROUP1: ~ Treatment + Timepoint + Treatment:Timepoint (Timepoint factor) ===')

## Set Treatment and Timepoint factor levels
dge_qc$samples$Treatment <- factor(dge_qc$samples$Treatment,
                                   levels = c('Vehicle', 'TRULI', 'Rotenone', 'TRULI + Rotenone'))
dge_qc$samples$Treatment <- relevel(dge_qc$samples$Treatment, ref = 'Vehicle')
dge_qc$samples$Timepoint <- factor(as.character(dge_qc$samples$Timepoint),
                                   levels = c('6', '16', '24', '48', '96'))
dge_qc$samples$Timepoint <- relevel(dge_qc$samples$Timepoint, ref = '6')

GENOTYPES    <- c('WT', 'LATS2KO')
MODEL_DESIGN <- '~ Treatment + Timepoint + Treatment:Timepoint'

## extract_factors for Treatment (factor) × Timepoint (factor)
## Interaction coefs: TreatmentX:TimepointY (both have level suffixes)
extract_factors <- function(coef) {
  if (grepl(':', coef)) {
    list(factor1_value = sub('^Treatment([^:]+):.+$', '\\1', coef),
         factor2_value = sub('^.+:Timepoint(.+)$',   '\\1', coef))
  } else if (grepl('^Treatment', coef)) {
    list(factor1_value = sub('^Treatment', '', coef), factor2_value = NA_character_)
  } else if (grepl('^Timepoint', coef)) {
    list(factor1_value = NA_character_, factor2_value = sub('^Timepoint', '', coef))
  } else {
    list(factor1_value = NA_character_, factor2_value = NA_character_)
  }
}

summary_rows_g1 <- list()

for (geno in GENOTYPES) {
  log_msg('  --- Genotype: ', geno, ' ---')

  dge_sub <- dge_qc[, dge_qc$samples$Genotype == geno]
  dge_sub$samples <- droplevels(dge_sub$samples)
  dge_sub$samples$Treatment <- relevel(dge_sub$samples$Treatment, ref = 'Vehicle')
  dge_sub$samples$Timepoint <- relevel(dge_sub$samples$Timepoint, ref = '6')
  log_msg('  Samples: ', ncol(dge_sub))
  log_msg('  Treatments: ', paste(levels(dge_sub$samples$Treatment), collapse = ', '))
  log_msg('  Timepoints: ', paste(levels(dge_sub$samples$Timepoint), collapse = ', '))

  sample_rows[[paste0('GROUP1_', geno)]] <- data.frame(
    group        = 'GROUP1',
    context      = paste0('Genotype=', geno),
    model_design = MODEL_DESIGN,
    n_samples    = ncol(dge_sub),
    sample_names = paste(dge_sub$samples$SampleName, collapse = ', '),
    stringsAsFactors = FALSE
  )

  design <- model.matrix(~ Treatment + Timepoint + Treatment:Timepoint, data = dge_sub$samples)
  log_msg('  Design matrix: ', nrow(design), ' x ', ncol(design))

  fit <- tryCatch({
    v <- voom(dge_sub, design, plot = FALSE); eBayes(lmFit(v, design))
  }, error = function(e) stop('GROUP1 voom/lmFit failed for ', geno, ': ', conditionMessage(e)))

  all_coefs   <- colnames(fit$coefficients)
  treat_coefs <- grep('^Treatment[^:]+$',   all_coefs, value = TRUE)
  time_coefs  <- grep('^Timepoint[^:]+$',   all_coefs, value = TRUE)
  inter_coefs <- grep('Treatment[^:]+:Timepoint|Timepoint[^:]+:Treatment', all_coefs, value = TRUE)
  log_msg('  Treatment terms: ', length(treat_coefs),
          ' | Timepoint terms: ', length(time_coefs),
          ' | Interaction terms: ', length(inter_coefs))

  file_prefix <- paste0('260320_MNP_LATS2KO_DE_GROUP1_', geno, '_')
  nsig_treat  <- list(); nsig_time <- list(); nsig_inter <- list()
  term_map    <- list('Treatment' = treat_coefs, 'Timepoint' = time_coefs,
                      'Treatment:Timepoint' = inter_coefs)

  for (term_label in names(term_map)) {
    for (coef in term_map[[term_label]]) {
      if (all(is.na(fit$coefficients[, coef]))) { log_msg('  SKIP (all NA): ', coef); next }
      log_msg('  [', term_label, '] ', coef)
      res <- tryCatch(topTable(fit, coef = coef, n = Inf, sort.by = 'P'),
                      error = function(e) { log_msg('  SKIP topTable: ', conditionMessage(e)); NULL })
      if (is.null(res)) next
      fv    <- extract_factors(coef)
      n_sig <- count_sig(res)
      log_msg('    Sig: ', n_sig)
      res$model_design  <- MODEL_DESIGN
      res$term          <- term_label
      res$factor1_value <- fv$factor1_value   # Treatment value (or NA)
      res$factor2_value <- fv$factor2_value   # Timepoint level (or NA)
      safe_coef <- safe_str(coef)
      write.csv(res, file.path(OUT_G1, paste0(file_prefix, safe_coef, '.csv')), row.names = FALSE)
      make_volcano(res, coef,
                   paste0('GROUP1 | ', geno, ' | ', term_label, ' | ', coef, '\n(n sig = ', n_sig, ')'),
                   paste0('260320_MNP_LATS2KO_Volcano_GROUP1_', geno, '_', safe_coef, '.png'),
                   OUT_G1)
      if (term_label == 'Treatment')         nsig_treat[[fv$factor1_value]] <- n_sig
      else if (term_label == 'Timepoint')   nsig_time[[fv$factor2_value]]  <- n_sig
      else nsig_inter[[paste0(fv$factor1_value, ':', fv$factor2_value)]]   <- n_sig
    }
  }

  for (key in names(nsig_inter)) {
    parts <- strsplit(key, ':')[[1]]
    summary_rows_g1[[paste0(geno, '_', key)]] <- data.frame(
      group             = 'GROUP1',
      Genotype          = geno,
      model_design      = MODEL_DESIGN,
      treatment_value   = parts[1],
      treatment_n_sig   = if (!is.null(nsig_treat[[parts[1]]])) nsig_treat[[parts[1]]] else NA_integer_,
      timepoint_value   = parts[2],
      timepoint_n_sig   = if (!is.null(nsig_time[[parts[2]]])) nsig_time[[parts[2]]] else NA_integer_,
      interaction_value = key,
      interaction_n_sig = nsig_inter[[key]],
      stringsAsFactors  = FALSE
    )
  }
}

## ── Write summaries ────────────────────────────────────────────────────────────
summary_g0 <- do.call(rbind, summary_rows_g0)
summary_g1 <- do.call(rbind, summary_rows_g1)

if (!is.null(summary_g0) && nrow(summary_g0) > 0) {
  write.csv(summary_g0,
            file.path(OUT_DE, '260320_MNP_LATS2KO_Summary_GROUP0_DEG_Counts.csv'), row.names = FALSE)
  log_msg('Saved: 260320_MNP_LATS2KO_Summary_GROUP0_DEG_Counts.csv')
}
if (!is.null(summary_g1) && nrow(summary_g1) > 0) {
  write.csv(summary_g1,
            file.path(OUT_DE, '260320_MNP_LATS2KO_Summary_GROUP1_DEG_Counts.csv'), row.names = FALSE)
  log_msg('Saved: 260320_MNP_LATS2KO_Summary_GROUP1_DEG_Counts.csv')
}
if (length(sample_rows) > 0) {
  write.csv(do.call(rbind, sample_rows),
            file.path(OUT_DE, '260320_MNP_LATS2KO_Summary_SampleInclusion.csv'), row.names = FALSE)
  log_msg('Saved: 260320_MNP_LATS2KO_Summary_SampleInclusion.csv')
}
log_msg('Script complete')
