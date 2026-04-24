#!/usr/bin/env Rscript

#########
# Usage
#########

# Replacement for /ref/bmlab/software/umi_dup_shared/Mtx_to_Counts.R
#
# Fixes: original script reads mapping file with sep="\t" but map files are
# space-separated with Windows line endings (\r\n), causing map$Barcode to be
# NULL and all column names to be set to NA.
#
# Three arguments required:
#   1) Directory from STARsolo with mtx files (Gene/raw/)
#   2) Sample mapping file (SampleName and Barcode columns, any whitespace separator)
#   3) Output directory for counts files
#
# Example:
#   Rscript Mtx_to_Counts_Fixed.R \
#     NRP_Production_1b/MD915_Pool1/Star/MD915_Pool1_Solo.out/Gene/raw/ \
#     NRP_Production_1b/MD915_Pool1_map.txt \
#     NRP_Production_1b/MD915_Pool1/Counts_Files/

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
    stop(paste("",
        "Three arguments required:",
        "  1) <Directory from STARsolo with mtx files>",
        "  2) SampleMapping.txt  (SampleName and Barcode columns)",
        "  3) Output directory for counts files",
        "",
        "Example:",
        "  Rscript Mtx_to_Counts_Fixed.R Solo.out/Gene/raw/ Map.txt Counts/",
        "", sep = "\n"))
}

library(data.table)
library(Matrix)

matrix_dir <- gsub("/$", "", args[1])
map_file   <- args[2]
output_dir <- gsub("/$", "", args[3])

##########
# Read mapping file
# - sep="" splits on any whitespace (handles both space- and tab-separated files)
# - strip.white removes leading/trailing whitespace from fields
# - fileEncoding handles Windows \r\n line endings
##########

map <- read.table(map_file, header = TRUE, sep = "", stringsAsFactors = FALSE,
                  strip.white = TRUE, fill = TRUE, fileEncoding = "latin1")

# Strip any residual carriage returns from all character columns
map <- as.data.frame(lapply(map, function(x) if (is.character(x)) gsub("\r", "", x) else x),
                     stringsAsFactors = FALSE)

if (!all(c("SampleName", "Barcode") %in% colnames(map))) {
    stop(paste("Mapping file must have 'SampleName' and 'Barcode' columns.",
               "Columns found:", paste(colnames(map), collapse = ", ")))
}

cat("Mapping file loaded:", nrow(map), "samples\n")

##########
# Read STARsolo barcodes and features
##########

barcode.names  <- fread(file.path(matrix_dir, "barcodes.tsv"),
                        header = FALSE, stringsAsFactors = FALSE, data.table = FALSE)
feature.names  <- fread(file.path(matrix_dir, "features.tsv"),
                        header = FALSE, stringsAsFactors = FALSE, data.table = FALSE)

# Match STAR barcode order to SampleNames via the mapping file
ordered_names <- map$SampleName[match(barcode.names$V1, map$Barcode)]

n_matched <- sum(!is.na(ordered_names))
n_na      <- sum(is.na(ordered_names))

cat("Barcodes matched:", n_matched, "of", nrow(barcode.names), "\n")

if (n_na > 0) {
    warning(paste(n_na, "barcodes in barcodes.tsv had no match in the mapping file."))
    cat("Unmatched barcodes:\n")
    cat(paste(barcode.names$V1[is.na(ordered_names)], collapse = "\n"), "\n")
}

##########
# Read and write Raw Counts (NoDedup)
##########

noDedup <- as.data.frame(as.matrix(readMM(file.path(matrix_dir, "umiDedup-NoDedup.mtx"))))
colnames(noDedup) <- ordered_names
rownames(noDedup) <- feature.names$V1

write.table(noDedup, file = file.path(output_dir, "Raw_Counts.txt"),
            sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)

totalCounts <- sum(colSums(noDedup))
cat("\n")
cat(paste0(c("Total_Mapped_Counts_StarSolo:", totalCounts), collapse = "\t"))
cat("\n")

##########
# Read and write Deduplicated Counts (1MM_All)
##########

dedup <- as.data.frame(as.matrix(readMM(file.path(matrix_dir, "umiDedup-1MM_All.mtx"))))
colnames(dedup) <- ordered_names
rownames(dedup) <- feature.names$V1

write.table(dedup, file = file.path(output_dir, "Dedup_Counts.txt"),
            sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)

totalDedupCounts <- sum(colSums(dedup))
cat(paste0(c("Total_Deduplicated_Counts_StarSolo:", totalDedupCounts), collapse = "\t"))
cat("\n")

percentUnique <- sprintf("%.2f", 100 * totalDedupCounts / totalCounts)
cat(paste0(c("Percent_Unique_Counts_StarSolo:", percentUnique), collapse = "\t"))
cat("\n")
