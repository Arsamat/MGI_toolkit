library(tidyverse)
library(edgeR)
library(limma)
library(ggrepel)

BASE_DIR <- '/scratch/rmlab/rmlab_shared3/alferreiro/BRB_Seq'
ANA_DIR  <- file.path(BASE_DIR, '260320_WT_MNP_Rotenone_Timecourse')
DATA_DIR <- file.path(ANA_DIR, 'data')
LOG_DIR  <- file.path(ANA_DIR, 'logs')
OUT_DE   <- file.path(ANA_DIR, 'DE_analysis')
OUT_G1   <- file.path(OUT_DE,  'Group1')

for (d in c(LOG_DIR, OUT_G1)) dir.create(d, showWarnings = FALSE, recursive = TRUE)
setwd(BASE_DIR)

PADJ_THRESH <- 0.05
LFC_THRESH  <- 0.6

log_file <- file.path(LOG_DIR, '260320_WT_MNP_DE_Analysis.log')
cat('', file = log_file)
log_msg <- function(...) {
  msg <- paste0('[', format(Sys.time(), '%Y-%m-%d %H:%M:%S'), '] ', ...)
  message(msg); cat(msg, '\n', file = log_file, append = TRUE)
}
log_msg('Script started')

dge <- readRDS(file.path(DATA_DIR, '260320_WT_MNP_DGElist.RDS'))
log_msg('Loaded DGEList: ', nrow(dge), ' genes x ', ncol(dge), ' samples')

dge_sub <- dge[, dge$samples$SampleType == 'Experiment' &
                  dge$samples$NumUnqDedupedReads >= 2e6]
dge_sub$samples <- droplevels(dge_sub$samples)
dge_sub$samples$Timepoint <- relevel(dge_sub$samples$Timepoint, ref = '6')
dge_sub$samples$Treatment <- relevel(dge_sub$samples$Treatment, ref = 'Vehicle')
log_msg('After subsetting: ', ncol(dge_sub), ' samples')
log_msg('Treatments: ', paste(levels(dge_sub$samples$Treatment), collapse = ', '))
log_msg('Timepoints: ', paste(levels(dge_sub$samples$Timepoint), collapse = ', '))

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

extract_factors <- function(coef) {
  if (grepl(':', coef)) {
    list(factor1_value = sub('^Timepoint([^:]+):.+$', '\\1', coef),
         factor2_value = sub('^.+:Treatment(.+)$',    '\\1', coef))
  } else if (grepl('^Timepoint', coef)) {
    list(factor1_value = sub('^Timepoint', '', coef), factor2_value = NA_character_)
  } else if (grepl('^Treatment', coef)) {
    list(factor1_value = NA_character_, factor2_value = sub('^Treatment', '', coef))
  } else {
    list(factor1_value = NA_character_, factor2_value = NA_character_)
  }
}

summary_rows <- list()
sample_rows  <- list()

design <- model.matrix(~ Timepoint + Treatment + Timepoint:Treatment, data = dge_sub$samples)
log_msg('Design matrix: ', nrow(design), ' x ', ncol(design))

fit <- tryCatch({
  v <- voom(dge_sub, design, plot = FALSE)
  eBayes(lmFit(v, design))
}, error = function(e) stop('voom/lmFit failed: ', conditionMessage(e)))

all_coefs   <- colnames(fit$coefficients)
time_coefs  <- grep('^Timepoint[^:]+$', all_coefs, value = TRUE)
treat_coefs <- grep('^Treatment[^:]+$', all_coefs, value = TRUE)
inter_coefs <- grep('Timepoint[^:]+:Treatment|Treatment[^:]+:Timepoint', all_coefs, value = TRUE)
log_msg('Timepoint terms: ', length(time_coefs), ' | Treatment terms: ', length(treat_coefs),
        ' | Interaction terms: ', length(inter_coefs))

sample_rows[['GROUP1_MNP']] <- data.frame(
  group = 'GROUP1', context = 'CellType=MNP',
  model_design = '~ Timepoint + Treatment + Timepoint:Treatment',
  n_samples = ncol(dge_sub), sample_names = paste(dge_sub$samples$SampleName, collapse = ', '),
  stringsAsFactors = FALSE
)

nsig_time <- list(); nsig_treat <- list(); nsig_inter <- list()
term_map  <- list('Timepoint' = time_coefs, 'Treatment' = treat_coefs, 'Timepoint:Treatment' = inter_coefs)

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
    res$model_design <- '~ Timepoint + Treatment + Timepoint:Treatment'
    res$term <- term_label; res$factor1_value <- fv$factor1_value; res$factor2_value <- fv$factor2_value
    safe_coef <- safe_str(coef)
    write.csv(res, file.path(OUT_G1, paste0('260320_WT_MNP_DE_GROUP1_MNP_', safe_coef, '.csv')), row.names = FALSE)
    make_volcano(res, coef, paste0('GROUP1 | MNP | ', term_label, ' | ', coef, '\n(n sig = ', n_sig, ')'),
                 paste0('260320_WT_MNP_Volcano_GROUP1_MNP_', safe_coef, '.png'), OUT_G1)
    if (term_label == 'Timepoint')      nsig_time[[fv$factor1_value]]  <- n_sig
    else if (term_label == 'Treatment') nsig_treat[[fv$factor2_value]] <- n_sig
    else nsig_inter[[paste0(fv$factor1_value, ':', fv$factor2_value)]] <- n_sig
  }
}

for (key in names(nsig_inter)) {
  parts <- strsplit(key, ':')[[1]]
  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    group = 'GROUP1', CellType = 'MNP',
    model_design = '~ Timepoint + Treatment + Timepoint:Treatment',
    factor1_value = parts[1],
    factor1_n_sig = if (!is.null(nsig_time[[parts[1]]])) nsig_time[[parts[1]]] else NA_integer_,
    factor2_value = parts[2],
    factor2_n_sig = if (!is.null(nsig_treat[[parts[2]]])) nsig_treat[[parts[2]]] else NA_integer_,
    interaction_value = key, interaction_n_sig = nsig_inter[[key]],
    stringsAsFactors = FALSE
  )
}

if (length(summary_rows) > 0) {
  write.csv(do.call(rbind, summary_rows),
            file.path(OUT_DE, '260320_WT_MNP_Summary_GROUP1_DEG_Counts.csv'), row.names = FALSE)
  log_msg('Saved: 260320_WT_MNP_Summary_GROUP1_DEG_Counts.csv')
}
if (length(sample_rows) > 0) {
  write.csv(do.call(rbind, sample_rows),
            file.path(OUT_DE, '260320_WT_MNP_Summary_SampleInclusion.csv'), row.names = FALSE)
  log_msg('Saved: 260320_WT_MNP_Summary_SampleInclusion.csv')
}

log_msg('Script complete')
