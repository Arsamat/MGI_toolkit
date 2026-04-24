library(tidyverse)
library(edgeR)
library(limma)
library(ggrepel)

## ── Setup ──────────────────────────────────────────────────────────────────────
BASE_DIR <- '/scratch/rmlab/rmlab_shared3/alferreiro/BRB_Seq'
ANA_DIR  <- file.path(BASE_DIR, '260318_Pilot_Analyses')
DATA_DIR <- file.path(ANA_DIR, 'data')
LOG_DIR  <- file.path(ANA_DIR, 'logs')
OUT_DE   <- file.path(ANA_DIR, 'DE_analysis_wEC90')
OUT_G1   <- file.path(OUT_DE,  'Group1')

for (d in c(DATA_DIR, LOG_DIR, OUT_G1)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}
setwd(BASE_DIR)

PADJ_THRESH <- 0.05
LFC_THRESH  <- 0.6

log_file <- file.path(LOG_DIR, '260318_Pilot_DE_Analysis_wEC90.log')
cat('', file = log_file)
log_msg <- function(...) {
  msg <- paste0('[', format(Sys.time(), '%Y-%m-%d %H:%M:%S'), '] ', ...)
  message(msg)
  cat(msg, '\n', file = log_file, append = TRUE)
}
log_msg('Script started')
log_msg('Significance threshold: adj.P.Val < ', PADJ_THRESH, ' & |logFC| > ', LFC_THRESH)

## ── Load and subset ────────────────────────────────────────────────────────────
dge <- readRDS(file.path(DATA_DIR, '260318_Pilot_DGElist.RDS'))
log_msg('Loaded DGEList: ', nrow(dge), ' genes x ', ncol(dge), ' samples')

# Exclude: 96hr timepoint, NoTreatment samples, low-coverage samples (EC90 retained)
dge_sub <- dge[, dge$samples$SampleType        == 'Experiment' &
                  dge$samples$NumUnqDedupedReads >= 2e6         &
                  dge$samples$Timepoint          != '96'        &
                  dge$samples$Treatment          != 'NoTreatment']
dge_sub$samples <- droplevels(dge_sub$samples)
log_msg('After subsetting (Experiment, >=2M reads, excl. 96hr & NoTreatment (EC90 retained)): ',
        ncol(dge_sub), ' samples')
log_msg('Timepoints remaining: ',
        paste(levels(dge_sub$samples$Timepoint), collapse = ', '))
log_msg('Doses remaining: ',
        paste(levels(dge_sub$samples$Dose), collapse = ', '))

## ── Helper functions ───────────────────────────────────────────────────────────
safe_str <- function(x) gsub('[^A-Za-z0-9_]', '_', x)

count_sig <- function(res) {
  sum(res$adj.P.Val < PADJ_THRESH & abs(res$logFC) > LFC_THRESH, na.rm = TRUE)
}

make_volcano <- function(res, coef_name, title_str, fname, out_dir) {
  res <- res %>%
    mutate(
      sig_status = adj.P.Val < PADJ_THRESH & abs(logFC) > LFC_THRESH,
      logP       = -log10(P.Value)
    )

  top_labels <- res %>%
    filter(sig_status) %>%
    slice_min(order_by = adj.P.Val, n = 50)

  p <- ggplot(res, aes(x = logFC, y = logP)) +
    geom_point(aes(color = sig_status), alpha = 0.3, size = 1) +
    geom_text_repel(data          = top_labels,
                    aes(label     = gene_symbol),
                    size          = 2.5,
                    max.overlaps  = 20,
                    segment.alpha = 0.5) +
    scale_color_manual(values = c('TRUE' = 'firebrick3', 'FALSE' = 'grey70')) +
    geom_vline(xintercept = c(-LFC_THRESH, LFC_THRESH), linetype = 'dotted') +
    geom_hline(yintercept = -log10(PADJ_THRESH), linetype = 'dotted') +
    labs(title = title_str,
         x     = paste0('logFC (', coef_name, ')'),
         y     = '-log10(p-value)') +
    theme_bw() +
    theme(legend.position = 'none')

  ggsave(file.path(out_dir, fname), p, width = 8, height = 6)
  log_msg('    Saved: ', fname)
}

extract_factors <- function(coef) {
  if (grepl(':', coef)) {
    list(factor1_value = sub('^Genotype([^:]+):.+$',  '\\1', coef),
         factor2_value = sub('^.+:Treatment(.+)$', '\\1', coef))
  } else if (grepl('^Genotype', coef)) {
    list(factor1_value = sub('^Genotype', '', coef),
         factor2_value = NA_character_)
  } else if (grepl('^Treatment', coef)) {
    list(factor1_value = NA_character_,
         factor2_value = sub('^Treatment', '', coef))
  } else {
    list(factor1_value = NA_character_, factor2_value = NA_character_)
  }
}

## ── Accumulators ───────────────────────────────────────────────────────────────
summary_g1_rows <- list()
sample_rows     <- list()

## ── GROUP1: Genotype × Treatment interaction (iMN only) ───────────────────────
log_msg('=== GROUP1: Genotype x Treatment (iMN) ===')

g1_design_str <- '~ Genotype + Treatment + Genotype:Treatment + Run'

dge_g1 <- dge_sub[, dge_sub$samples$CellType == 'iMN']
dge_g1$samples <- droplevels(dge_g1$samples)
dge_g1$samples$Genotype  <- relevel(dge_g1$samples$Genotype,  ref = 'WT')
dge_g1$samples$Treatment <- relevel(dge_g1$samples$Treatment, ref = 'Vehicle')
dge_g1$samples$Run       <- as.factor(dge_g1$samples$Run)

log_msg('  Samples: ', ncol(dge_g1),
        ' | Genotypes: ', paste(levels(dge_g1$samples$Genotype), collapse = ', '),
        ' | Treatments: ', paste(levels(dge_g1$samples$Treatment), collapse = ', '),
        ' | Timepoints: ', paste(levels(dge_g1$samples$Timepoint), collapse = ', '),
        ' | Runs: ', paste(levels(dge_g1$samples$Run), collapse = ', '))

sample_rows[['GROUP1_iMN']] <- data.frame(
  group        = 'GROUP1',
  context      = 'CellType=iMN',
  model_design = g1_design_str,
  n_samples    = ncol(dge_g1),
  sample_names = paste(dge_g1$samples$SampleName, collapse = ', '),
  stringsAsFactors = FALSE
)

if (nlevels(dge_g1$samples$Genotype) < 2 || nlevels(dge_g1$samples$Treatment) < 2) {
  stop('Fewer than 2 genotypes or treatments — cannot fit GROUP1 model')
}

design1 <- tryCatch(
  model.matrix(~ Genotype + Treatment + Genotype:Treatment + Run,
               data = dge_g1$samples),
  error = function(e) stop('model.matrix failed: ', conditionMessage(e))
)

fit1 <- tryCatch({
  v1 <- voom(dge_g1, design1, plot = FALSE)
  f  <- lmFit(v1, design1)
  eBayes(f)
}, error = function(e) stop('voom/lmFit failed: ', conditionMessage(e)))

all_coefs   <- colnames(fit1$coefficients)
geno_coefs  <- grep('^Genotype[^:]+$',  all_coefs, value = TRUE)
treat_coefs <- grep('^Treatment[^:]+$', all_coefs, value = TRUE)
inter_coefs <- grep('Genotype[^:]+:Treatment|Treatment[^:]+:Genotype',
                    all_coefs, value = TRUE)

log_msg('  Genotype terms: ',    length(geno_coefs),
        ' | Treatment terms: ',  length(treat_coefs),
        ' | Interaction terms: ', length(inter_coefs))

g1_nsig_geno  <- list()
g1_nsig_treat <- list()
g1_nsig_inter <- list()

term_coef_map <- list(
  'Genotype'           = geno_coefs,
  'Treatment'          = treat_coefs,
  'Genotype:Treatment' = inter_coefs
)

for (term_label in names(term_coef_map)) {
  for (coef in term_coef_map[[term_label]]) {
    if (all(is.na(fit1$coefficients[, coef]))) {
      log_msg('  SKIP coef (all NA): ', coef); next
    }
    log_msg('  [', term_label, '] ', coef)

    res <- tryCatch(
      topTable(fit1, coef = coef, n = Inf, sort.by = 'P'),
      error = function(e) { log_msg('  SKIP topTable: ', conditionMessage(e)); NULL }
    )
    if (is.null(res)) next

    fv    <- extract_factors(coef)
    n_sig <- count_sig(res)
    log_msg('    Sig: ', n_sig)

    res$model_design  <- g1_design_str
    res$term          <- term_label
    res$factor1_value <- fv$factor1_value
    res$factor2_value <- fv$factor2_value

    safe_coef     <- safe_str(coef)
    csv_fname     <- paste0('260318_wEC90_DE_GROUP1_iMN_', safe_coef, '.csv')
    volcano_fname <- paste0('260318_wEC90_Volcano_GROUP1_iMN_', safe_coef, '.png')

    write.csv(res, file.path(OUT_G1, csv_fname), row.names = FALSE)
    log_msg('    Saved: ', csv_fname)

    make_volcano(res, coef,
                 title_str = paste0('GROUP1 | iMN | ', term_label, ' | ', coef,
                                    '\n(n sig = ', n_sig, ')'),
                 fname     = volcano_fname,
                 out_dir   = OUT_G1)

    if (term_label == 'Genotype') {
      g1_nsig_geno[[fv$factor1_value]] <- n_sig
    } else if (term_label == 'Treatment') {
      g1_nsig_treat[[fv$factor2_value]] <- n_sig
    } else {
      inter_key <- paste0(fv$factor1_value, ':', fv$factor2_value)
      g1_nsig_inter[[inter_key]] <- n_sig
    }
  }
}

for (inter_key in names(g1_nsig_inter)) {
  parts     <- strsplit(inter_key, ':')[[1]]
  geno_lev  <- parts[1]
  treat_lev <- parts[2]
  summary_g1_rows[[length(summary_g1_rows) + 1]] <- data.frame(
    group             = 'GROUP1',
    CellType          = 'iMN',
    model_design      = g1_design_str,
    factor1_value     = geno_lev,
    factor1_n_sig     = if (!is.null(g1_nsig_geno[[geno_lev]]))   g1_nsig_geno[[geno_lev]]  else NA_integer_,
    factor2_value     = treat_lev,
    factor2_n_sig     = if (!is.null(g1_nsig_treat[[treat_lev]])) g1_nsig_treat[[treat_lev]] else NA_integer_,
    interaction_value = inter_key,
    interaction_n_sig = g1_nsig_inter[[inter_key]],
    stringsAsFactors  = FALSE
  )
}

## ── Write summary tables ───────────────────────────────────────────────────────
if (length(summary_g1_rows) > 0) {
  df_sum_g1 <- do.call(rbind, summary_g1_rows)
  write.csv(df_sum_g1,
            file.path(OUT_DE, '260318_wEC90_Summary_GROUP1_DEG_Counts.csv'),
            row.names = FALSE)
  log_msg('Saved: DE_analysis/260318_wEC90_Summary_GROUP1_DEG_Counts.csv (',
          nrow(df_sum_g1), ' rows)')
}

if (length(sample_rows) > 0) {
  df_samples <- do.call(rbind, sample_rows)
  write.csv(df_samples,
            file.path(OUT_DE, '260318_wEC90_Summary_SampleInclusion.csv'),
            row.names = FALSE)
  log_msg('Saved: DE_analysis/260318_wEC90_Summary_SampleInclusion.csv (',
          nrow(df_samples), ' rows)')
}

log_msg('Script complete')
