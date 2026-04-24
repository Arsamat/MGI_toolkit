library(tidyverse)
library(fgsea)
library(org.Hs.eg.db)
library(msigdbr)
library(patchwork)

## ── Setup ──────────────────────────────────────────────────────────────────────
BASE_DIR <- '/scratch/rmlab/rmlab_shared3/alferreiro/BRB_Seq'
setwd(BASE_DIR)
message('[', Sys.time(), '] Script started')

## ── Custom gene sets ──────────────────────────────────────────────────────────
pLATS_UP <- c(
  "TARDBP","LATS1","LATS2","STK3","STK4","SAV1","MOB1A","MOB1B","NF2","FRMD6",
  "WWC1","WWC2","AMOT","AMOTL1","AMOTL2","TAOK1","TAOK2","TAOK3",
  "MAP4K1","MAP4K2","MAP4K3","MAP4K4","MAP4K5","MAP4K6","MAP4K7",
  "RASSF1","RASSF5","PTPN14","STRN","STRIP1","RHOA","ROCK1","ROCK2",
  "YES1","WWTR1","YAP1",
  "YWHAB","YWHAE","YWHAG","YWHAH","YWHAQ","YWHAZ","SFN",
  "STMN2","UNC13A"
)
YAP_TAZ_TARGETS <- c(
  "CCN2","CCN1","ANKRD1","AREG","BIRC5","SERPINE1","SPP1",
  "THBS1","FSTL1","ITGB1","ITGA5","ITGB3","COL1A1","COL3A1",
  "COL5A1","MYC","CCND1","FOSL1","FN1","FOXM1","AXL","AMOTL2",
  "PDGFB","FGF2","HK2","SDC4","CTNNB1","CLU","EPHA2","LGALS1",
  "ITGB6","EGFR"
)
YAP_TAZ_REPRESSED <- c(
  "CDKN1A","CDKN2A","CDKN2B","GDF15","IGFBP3","IGFBP5","IGFBP7",
  "BTG2","DUSP5","PPP1R15A","TP53INP1","IER3","GADD45A",
  "GADD45B","TNFRSF10B","FOSL2","MMP1","MMP3","IL6","IL1A",
  "IL1B","CXCL8","CXCL5","ICAM1","SOCS3","PLAU","PLAUR",
  "CXCL1","CDKN2C","SAT1","SOD2"
)

## ── Build pathway list ────────────────────────────────────────────────────────
message('[', Sys.time(), '] Building pathway list...')
h_df         <- msigdbr(species = "Homo sapiens", category = "H")
hallmark_list <- split(as.character(h_df$entrez_gene), h_df$gs_name)

custom_sets_sym <- list(
  CUSTOM_HIPPO_pLATS_UP    = pLATS_UP,
  CUSTOM_YAP_TAZ_TARGETS   = YAP_TAZ_TARGETS,
  CUSTOM_YAP_TAZ_REPRESSED = YAP_TAZ_REPRESSED
)
all_custom_syms   <- unique(unlist(custom_sets_sym))
entrez_custom_map <- AnnotationDbi::mapIds(
  org.Hs.eg.db, keys = all_custom_syms,
  column = "ENTREZID", keytype = "SYMBOL", multiVals = "first"
)
custom_list <- lapply(custom_sets_sym, function(syms) {
  ez <- entrez_custom_map[syms]
  as.character(ez[!is.na(ez)])
})
all_pathways <- c(hallmark_list, custom_list)
message('[', Sys.time(), '] Pathways: ', length(all_pathways), ' gene sets')

## ── Helpers ───────────────────────────────────────────────────────────────────
clean_desc <- function(x) {
  x <- sub("^HALLMARK_", "", x)
  x <- sub("^CUSTOM_",   "", x)
  gsub("_", " ", x)
}

## Strip leading YYMMDD_ prefix from a filename stem
strip_date <- function(x) sub("^\\d{6}_", "", x)

run_gsea_for_file <- function(df, sym2entrez, label) {
  df$entrez <- sym2entrez[df$gene_symbol]
  df_m <- df %>%
    filter(!is.na(entrez)) %>%
    group_by(entrez) %>%
    slice_max(abs(logFC), n = 1, with_ties = FALSE) %>%
    ungroup()
  ranked <- sort(setNames(df_m$logFC, df_m$entrez), decreasing = TRUE)

  if (length(ranked) < 10) {
    message('  [', label, '] Too few mapped genes (', length(ranked), '), skipping')
    return(NULL)
  }

  set.seed(42)
  gsea_res <- tryCatch(
    fgsea::fgseaMultilevel(all_pathways, stats = ranked,
                           minSize = 10, maxSize = 500, eps = 0),
    error = function(e) { message('  Error: ', e$message); NULL }
  )
  if (is.null(gsea_res) || nrow(gsea_res) == 0) return(NULL)

  res_df <- as.data.frame(gsea_res) %>%
    dplyr::rename(ID = pathway, p.adjust = padj) %>%
    filter(!is.na(p.adjust) & p.adjust < 0.05) %>%
    arrange(p.adjust) %>%
    mutate(
      Description = clean_desc(ID),
      leadingEdge = sapply(leadingEdge, paste, collapse = ';')
    )

  if (nrow(res_df) == 0) NULL else res_df
}

save_barplot <- function(sig_df, label, out_path) {
  plot_df <- sig_df %>%
    arrange(p.adjust) %>%
    head(10) %>%
    mutate(
      neg_log10_p = -log10(p.adjust),
      Direction   = ifelse(NES > 0, "Up (NES > 0)", "Down (NES < 0)"),
      Description = factor(Description, levels = rev(Description))
    )
  p <- ggplot(plot_df, aes(x = neg_log10_p, y = Description, fill = Direction)) +
    geom_col() +
    geom_vline(xintercept = -log10(0.05), linetype = "dashed", colour = "grey40") +
    scale_fill_manual(values = c("Up (NES > 0)" = "#CB181D", "Down (NES < 0)" = "#2171B5")) +
    labs(title    = paste0("Top pathways: ", label),
         subtitle = "GSEA Hallmark + custom gene sets; FDR < 0.05",
         x        = expression(-log[10](adj.p)), y = NULL, fill = "Direction") +
    theme_bw() +
    theme(plot.title    = element_text(size = 11, face = "bold"),
          plot.subtitle = element_text(size = 8),
          axis.text.y   = element_text(size = 8),
          legend.position = "bottom")
  ggsave(out_path, p, width = 9, height = 5, dpi = 150)
  message('  Saved barplot: ', basename(out_path))
}

save_dotplot <- function(gsea_sig_list, dir_label, out_path) {
  all_sig <- bind_rows(lapply(names(gsea_sig_list), function(lb) {
    df <- gsea_sig_list[[lb]]
    if (is.null(df) || nrow(df) == 0) return(NULL)
    df %>% dplyr::select(ID, Description, NES, p.adjust) %>% mutate(term = lb)
  }))

  if (is.null(all_sig) || nrow(all_sig) == 0) {
    message('  No significant pathways for dotplot in: ', dir_label)
    return(invisible(NULL))
  }

  top_ids <- all_sig %>%
    group_by(ID, Description) %>%
    summarise(n_terms = n(), mean_abs_NES = mean(abs(NES)), .groups = "drop") %>%
    arrange(desc(n_terms), desc(mean_abs_NES)) %>%
    head(20) %>%
    pull(ID)

  all_terms <- unique(all_sig$term)
  nes_mat <- outer(
    top_ids, all_terms,
    Vectorize(function(pid, trm) {
      v <- all_sig$NES[all_sig$ID == pid & all_sig$term == trm]
      if (length(v) == 0) 0 else v[1]
    })
  )
  rownames(nes_mat) <- top_ids; colnames(nes_mat) <- all_terms

  term_order <- if (length(all_terms) > 1) {
    all_terms[hclust(dist(t(nes_mat)), method = "average")$order]
  } else all_terms

  path_order_ids <- top_ids[hclust(dist(nes_mat), method = "average")$order]
  path_desc_order <- all_sig %>%
    dplyr::distinct(ID, Description) %>%
    dplyr::slice(match(path_order_ids, ID)) %>%
    pull(Description)

  dot_df <- all_sig %>%
    filter(ID %in% top_ids) %>%
    mutate(
      neg_log10_p = -log10(p.adjust),
      term        = factor(term, levels = term_order),
      Description = factor(Description, levels = path_desc_order)
    )

  p <- ggplot(dot_df, aes(x = term, y = Description, size = neg_log10_p, colour = NES)) +
    geom_point() +
    scale_colour_gradient2(low = "#2171B5", mid = "white", high = "#CB181D",
                           midpoint = 0, name = "NES") +
    scale_size_continuous(name = expression(-log[10](adj.p)), range = c(1, 6)) +
    labs(title    = paste0("Top 20 pathways: ", dir_label),
         subtitle = "GSEA Hallmark + custom; FDR < 0.05; size = -log10(adj.p), colour = NES",
         x = NULL, y = NULL) +
    theme_bw() +
    theme(axis.text.x   = element_text(angle = 45, hjust = 1, size = 7),
          axis.text.y   = element_text(size = 8),
          plot.title    = element_text(size = 12, face = "bold"),
          plot.subtitle = element_text(size = 9),
          legend.position = "right")

  n_terms <- length(unique(dot_df$term))
  ggsave(out_path, p,
         width  = max(8, 2 + 0.6 * n_terms),
         height = 8, dpi = 150)
  message('  Saved dotplot: ', basename(out_path))
}

## ── Target Group directories ──────────────────────────────────────────────────
## Each entry: path to a Group0 or Group1 subdirectory
group_dirs <- c(
  file.path(BASE_DIR, '260318_Pilot_Analyses/DE_analysis_wEC90/Group1'),
  file.path(BASE_DIR, '260320_iPSC_WT_LATS2KO_Rot_Tun_Timecourse_Analyses/DE_analysis/Group1'),
  file.path(BASE_DIR, '260320_MFN2_MNP_Rotenone_Analyses/DE_analysis/Group1'),
  file.path(BASE_DIR, '260320_MNP_WT_LATS2KO_Rot_TRULI_Timecourse_Analyses/DE_analysis/Group0'),
  file.path(BASE_DIR, '260320_MNP_WT_LATS2KO_Rot_TRULI_Timecourse_Analyses/DE_analysis/Group1'),
  file.path(BASE_DIR, '260320_WT_MNP_Rotenone_Timecourse/DE_analysis/Group1')
)

## ── Process each Group directory ─────────────────────────────────────────────
for (group_dir in group_dirs) {
  if (!dir.exists(group_dir)) {
    message('[', Sys.time(), '] SKIP (does not exist): ', group_dir)
    next
  }

  ## Short label for plot titles: strip leading date from analysis dir + Group subdir
  dir_label <- paste0(
    strip_date(basename(dirname(dirname(group_dir)))), '/',
    basename(group_dir)
  )

  message('[', Sys.time(), '] === Processing: ', dir_label, ' ===')

  ## Find DE result CSVs (exclude Summary files)
  csv_files <- list.files(group_dir, pattern = '\\.csv$', full.names = TRUE)
  csv_files <- csv_files[!grepl('Summary', basename(csv_files))]

  if (length(csv_files) == 0) {
    message('  No DE CSVs found, skipping')
    next
  }
  message('  Found ', length(csv_files), ' DE CSV(s)')

  gsea_out_dir <- file.path(group_dir, 'GSEA')
  dir.create(gsea_out_dir, showWarnings = FALSE, recursive = TRUE)

  ## Build sym2entrez once for all files in this directory
  all_dfs  <- lapply(csv_files, read.csv, stringsAsFactors = FALSE)
  all_syms <- unique(unlist(lapply(all_dfs, function(d) d$gene_symbol)))
  sym2entrez <- AnnotationDbi::mapIds(
    org.Hs.eg.db, keys = all_syms,
    column = "ENTREZID", keytype = "SYMBOL", multiVals = "first"
  )
  message('  Mapped ', sum(!is.na(sym2entrez)), '/', length(sym2entrez), ' symbols to Entrez')

  gsea_sig_dir <- list()   # label -> sig df, for combined dotplot

  for (i in seq_along(csv_files)) {
    f     <- csv_files[[i]]
    stem  <- tools::file_path_sans_ext(basename(f))
    label <- strip_date(stem)
    df    <- all_dfs[[i]]
    message('[', Sys.time(), '] GSEA: ', label, ' (', nrow(df), ' genes)')

    sig_df <- run_gsea_for_file(df, sym2entrez, label)

    if (is.null(sig_df)) {
      message('  No significant pathways')
      gsea_sig_dir[[label]] <- data.frame(
        ID = character(), Description = character(),
        NES = numeric(), p.adjust = numeric(), leadingEdge = character()
      )
      next
    }

    message('  Sig pathways: ', nrow(sig_df))

    ## Save full GSEA results CSV
    csv_out <- file.path(gsea_out_dir, paste0('260323_GSEA_', label, '.csv'))
    write.csv(sig_df, csv_out, row.names = FALSE)
    message('  Saved: ', basename(csv_out))

    ## Save barplot (top 10)
    bar_out <- file.path(gsea_out_dir, paste0('260323_Barplot_', label, '.png'))
    save_barplot(sig_df, label, bar_out)

    gsea_sig_dir[[label]] <- sig_df %>%
      dplyr::select(ID, Description, NES, p.adjust)
  }

  ## Combined dotplot for this Group directory
  dot_out <- file.path(gsea_out_dir,
    paste0('260323_Dotplot_Top20_', gsub('/', '_', dir_label), '.png'))
  save_dotplot(gsea_sig_dir, dir_label, dot_out)
}

message('[', Sys.time(), '] Script complete')
