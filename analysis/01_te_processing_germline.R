#!/usr/bin/env Rscript

# TE Data Processing Script
# This script processes transposable element (TE) data and saves various data tables as R objects # nolint

#### SETUP ####
r_dir <- "/Users/briannelaverty/Documents/R_Malkin/te/data/R_obj/germline/"
plot_dir <- "/Users/briannelaverty/Library/Mobile Documents/com~apple~CloudDocs/Malkin_lab/figures/te/graphs/germline/"
r_dir_files <- "/Users/briannelaverty/Library/Mobile Documents/com~apple~CloudDocs/Malkin_lab/figures/te/graphs/germline/files/"
hpc_dir <- "/Users/briannelaverty/Documents/R_Malkin/te/hpc_data/processed/"

# Create files directory if it doesn't exist
if (!dir.exists(r_dir_files)) {
  dir.create(r_dir_files, recursive = TRUE)
}

# Configuration flags
INCLUDE_ALU <- TRUE
INCLUDE_SVA <- TRUE

# Processing options
PROCESS_COMMON_TES <- FALSE  # Set to TRUE to process common TEs (slower)
GENERATE_FREQUENCY_PLOTS <- FALSE  # Set to TRUE to generate TE frequency distribution and HostSeq filtering plots (slower)

# Pipeline tracking options
ENABLE_DETAILED_TRACKING <- TRUE # Set to TRUE to enable detailed pipeline tracking (slower)

#### LIBRARIES ####
suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
  library(ggsignif)
  library(tidyr)
  library(umap)
  library(pheatmap)
  library(purrr)
  library(RColorBrewer)
  library(caret)
  library("org.Hs.eg.db")
  library(clusterProfiler)
  library(AnnotationDbi)
  library(stringr)
  library(scales)
  library(glmnet)
  library(reshape2)
  library(ggrepel)
  library(GO.db)
  library(tibble)
  library(VennDiagram)
  library(eulerr)
  library(enrichplot)
  library(ComplexHeatmap)
  library(circlize)
  library(dplyr)
})

# Source functions
source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/functions_te.R")

#### LOAD DATA ####
cat("Loading raw germline data...\n")
te_raw <- fread("/Users/briannelaverty/Documents/R_Malkin/te/data/final/germline_annotSV_full.tsv", sep="\t", header=TRUE)

# Store initial counts for pipeline tracking (before any filtering)
if (ENABLE_DETAILED_TRACKING) {
  step0_raw_te_count <- nrow(te_raw)
  step0_sample_count <- length(unique(te_raw$Samples_ID))
}

te_split <- fread("/Users/briannelaverty/Documents/R_Malkin/te/data/final/germline_annotSV_split.tsv", sep="\t", header=TRUE)
chr_length <- read.csv("/Users/briannelaverty/Documents/R_Malkin/te/data/raw/chromosome_length.csv", header=TRUE)
clinical <- read.delim("/Users/briannelaverty/Documents/R_Malkin/clinical/te_clinical.csv", sep=",", header=TRUE)
metrics <- read.delim("/Users/briannelaverty/Documents/R_Malkin/te/data/final/combined_metrics.txt", sep="\t", header=TRUE)
hg37_genes <- read.delim("/Users/briannelaverty/Documents/R_Malkin/te/data/final/hg37_genes.tsv", sep="\t", header=TRUE)
nonproband <- read.delim("/Users/briannelaverty/Documents/R_Malkin/te/data/final/non_proband_normals", sep="\t", header=TRUE, colClasses = c("character"))
noconsent <- read.delim("/Users/briannelaverty/Documents/R_Malkin/te/data/final/kics_samples_exclude", sep="\t", header=TRUE, colClasses = c("character"))
ancestry <- read.delim("/Users/briannelaverty/Documents/R_Malkin/te/data/final/ancestry.txt", sep="\t", header=TRUE)
unique_l1_tp <- as.data.table(readr::read_csv("/Users/briannelaverty/Documents/R_Malkin/TE/data/final/unique_l1_tp.csv", col_names=c("sample", "ID", "caller"), locale=readr::locale(encoding="UTF-8"), show_col_types=FALSE))
unique_l1_fp <- as.data.table(readr::read_csv("/Users/briannelaverty/Documents/R_Malkin/TE/data/final/unique_l1_fp.csv", locale=readr::locale(encoding="UTF-8"), show_col_types=FALSE))
# Keep only first 3 columns and rename to match expected structure
if (ncol(unique_l1_fp) >= 3) {
  unique_l1_fp <- unique_l1_fp[, 1:3]
  colnames(unique_l1_fp) <- c("sample", "ID", "caller")
}
unique_l1_master <- as.data.table(readr::read_csv("/Users/briannelaverty/Documents/R_Malkin/TE/data/final/unique_l1_master.csv", locale=readr::locale(encoding="UTF-8"), show_col_types=FALSE))
l1_merge_fp <- fread("/Users/briannelaverty/Documents/R_Malkin/te/IGV_tracking/L1_merge_sample_id.csv")

#### PREP DATA ####
cat("Preparing data...\n")

# Metrics
metrics <- prep_metrics(metrics)

# Clinical
clinical <- prep_clinical(clinical) 

# Host seq with cancer
hostseq_cancer <- prep_hostseq(clinical)

# Chromosome length
chr_length$chr <- factor(chr_length$chr, levels=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"X","Y")) # nolint
chr_lengths <- setNames(chr_length$length, chr_length$chr)

# Gene size
gene_size <- add_gene_size(hg37_genes)

# Ancestry
ancestry <- prep_ancestry(ancestry)

# Save all prepared data objects
cat("Saving prepared data objects...\n")
save(ancestry, file = paste0(r_dir, "ancestry.RData"))
save(gene_size, file = paste0(r_dir, "gene_size.RData"))
save(chr_lengths, file = paste0(r_dir, "chr_lengths.RData"))
save(hostseq_cancer, file = paste0(r_dir, "hostseq_cancer.RData"))
save(clinical, file = paste0(r_dir, "clinical.RData"))
save(metrics, file = paste0(r_dir, "metrics.RData"))

# Merge clinical and ancestry data and save as CSV
cat("Merging clinical and ancestry data...\n")
clinical_ancestry <- merge_dfs(clinical, ancestry, include_all_x = TRUE, print_info = TRUE, dataset_name = "ancestry")
processed_dir <- "/Users/briannelaverty/Documents/R_Malkin/te/hpc_data/processed/"
if (!dir.exists(processed_dir)) {
  dir.create(processed_dir, recursive = TRUE)
  cat("Created directory:", processed_dir, "\n")
}
clinical_ancestry_file <- paste0(processed_dir, "clinical_ancestry.csv")
write.csv(clinical_ancestry, clinical_ancestry_file, row.names = FALSE)
cat("✓ Saved merged clinical and ancestry data to:", clinical_ancestry_file, "\n")

# Filter TE types based on configuration
cat("\n========================================\n")
cat("STEP 1: Filtering TE types\n")
cat("========================================\n")
te_raw <- filter_te_types(te_raw, include_alu = INCLUDE_ALU, include_sva = INCLUDE_SVA)
if (ENABLE_DETAILED_TRACKING) {
  step1 <- init_step_tracking(te_raw, "Samples_ID")
  step1_te_count <- step1$te_count
  step1_sample_count <- step1$sample_count
  step1_te_loss <- step0_raw_te_count - step1_te_count
}

# Germline has no IGV review - skip false positive filtering
cat("\n========================================\n")
cat("STEP 2: No false positive filtering (germline has no IGV review)\n")
cat("========================================\n")
cat("Skipping false positive filtering step for germline data.\n")

if (ENABLE_DETAILED_TRACKING) {
  step2_te_count <- step1_te_count
  step2_sample_count <- step1_sample_count
  step2_te_loss <- 0
  step2_sample_loss <- 0
}

# Prepare te_raw
cat("\n========================================\n")
cat("STEP 3: Quality filtering (prep_te)\n")
cat("========================================\n")

# Always capture prep_te output for processed_output.txt
prep_te_stdout <- capture.output({
  te_raw_prepped <- prep_te(te_raw, nonproband, noconsent, hostseq_cancer, metrics,
                           c("mean_cov", "avg_quality", "pct_chimeras"), c(20, 25, 2), type="N",
                           export_filtered = TRUE, output_dir = r_dir_files)
})
write_stdout_to_file(prep_te_stdout, paste0(r_dir_files, "/processed_output.txt"),
                     "STDOUT for prep_te", append = FALSE)  # First write: create new file
# Print the output to console as well
cat(paste(prep_te_stdout, collapse = "\n"), "\n")

if (ENABLE_DETAILED_TRACKING) {
  step3 <- init_step_tracking(te_raw_prepped, "sample")
  step3_te_count <- step3$te_count
  step3_sample_count <- step3$sample_count
  step3_te_loss <- step2_te_count - step3_te_count
  step3_sample_loss <- step2_sample_count - step3_sample_count
}
save(te_raw_prepped, file = paste0(r_dir, "te_raw_prepped.RData"))

#### SPLIT HOSTSEQ SAMPLES ####
cat("\n========================================\n")
cat("STEP 3.5: Splitting HostSeq samples\n")
cat("========================================\n")
cat("Splitting HostSeq samples into filtering (66%) and analysis (33%) groups...\n")

# Merge ancestry data for stratified splitting
cat("Merging ancestry data for stratified split...\n")
te_raw_prepped_with_ancestry <- merge_dfs(te_raw_prepped, ancestry, include_all_x = TRUE, print_info = FALSE, dataset_name = "ancestry")

# Split HostSeq samples by adding labels (stratified by ancestry)
hostseq_split <- split_hostseq_samples(
  te_data = te_raw_prepped_with_ancestry,
  filter_pct = 66,
  seed = 123
)

# Save split information
save(hostseq_split, file = paste0(r_dir, "hostseq_split.RData"))

# Use the labeled data for processing
te_raw_labeled <- hostseq_split$te_data

#### PLOT HOSTSEQ ANCESTRY DISTRIBUTIONS (BEFORE ANY PROCESSING) ####
cat("\n========================================\n")
cat("Plotting HostSeq ancestry distributions\n")
cat("========================================\n")
cat("Generating ancestry pie charts for: All HostSeq, Filter group (66%), and Analysis group (34%)\n")
cat("This must be done BEFORE any processing/filtering to capture both groups\n\n")

# Plot HostSeq ancestry pie charts using te_raw_labeled which has both filter and analysis groups
# Create plots for both predicted_ancestry_thres and mapped_label
tryCatch({
  hostseq_ancestry_pies <- plot_hostseq_ancestry_pies(
    te_data = te_raw_labeled,
    output_dir = plot_dir,
    plot_prefix = "germline",
    ancestry_col = "predicted_ancestry_thres"
  )
  cat("✓ HostSeq ancestry pie charts (predicted_ancestry_thres) saved to:", plot_dir, "\n\n")
}, error = function(e) {
  cat("Warning: Could not create HostSeq ancestry pie charts (predicted_ancestry_thres):", e$message, "\n\n")
})

tryCatch({
  hostseq_ancestry_pies_mapped <- plot_hostseq_ancestry_pies(
    te_data = te_raw_labeled,
    output_dir = plot_dir,
    plot_prefix = "germline",
    ancestry_col = "mapped_label"
  )
  cat("✓ HostSeq ancestry pie charts (mapped_label) saved to:", plot_dir, "\n\n")
}, error = function(e) {
  cat("Warning: Could not create HostSeq ancestry pie charts (mapped_label):", e$message, "\n\n")
})

#### REMOVE ANCESTRY COLUMNS BEFORE PROCESSING ####
cat("\n========================================\n")
cat("Removing ancestry columns from labeled data\n")
cat("========================================\n")
cat("Ancestry was needed for stratified HostSeq split and pie charts.\n")
cat("Removing now so process_te_data_germline() can do proper merge.\n\n")

# Remove ancestry columns (they will be merged again in process_te_data_germline)
ancestry_cols <- c("predicted_ancestry_thres", "mapped_label", "base_sample")
te_raw_labeled <- te_raw_labeled %>%
  select(-any_of(ancestry_cols))

cat("✓ Ancestry columns removed\n\n")

#### PLOT TE FREQUENCY DISTRIBUTIONS (OPTIONAL) ####
if (GENERATE_FREQUENCY_PLOTS) {
  cat("\n========================================\n")
  cat("STEP 3.6: Plotting TE frequency distributions\n")
  cat("========================================\n")
  cat("Creating histograms of TE frequency in gnomAD before filtering...\n")

  freq_plots <- plot_te_frequency_distributions(
    te_data = te_raw_prepped,
    output_dir = plot_dir,
    output_prefix = "germline"
  )

  #### PLOT HOSTSEQ TE FREQUENCY ####
  cat("\n========================================\n")
  cat("STEP 3.7: Plotting HostSeq TE frequency distribution (filter group 66%)\n")
  cat("========================================\n")
  cat("Creating histogram of HostSeq TEs and their frequency in HostSeq filter group (66%)...\n")

  # Filter to HostSeq filter samples only (using labeled data)
  te_hostseq_only <- te_raw_labeled %>%
    filter(hostseq_group == "filter")

  # Convert to expanded format for this plot
  te_hostseq_expand <- te_hostseq_only %>%
    distinct(sample, SV_chrom, SV_start, ALT, .keep_all = TRUE)

  hostseq_freq_plot <- plot_te_hostseq_frequency(
    te_expand = te_hostseq_expand,
    output_dir = plot_dir,
    output_prefix = "hostseq"
  )

  #### PLOT GERMLINE TE FREQUENCY IN HOSTSEQ ####
  cat("\n========================================\n")
  cat("STEP 3.8: Plotting germline TEs and their HostSeq frequency (filter group 66%)\n")
  cat("========================================\n")
  cat("Creating histogram of germline TEs and their frequency in HostSeq filter group (66%)...\n")

  # Convert to expanded format for this plot
  te_filter_expand <- te_germline_filter %>%
    distinct(sample, SV_chrom, SV_start, ALT, .keep_all = TRUE)

  germline_hostseq_freq_plot <- plot_te_hostseq_frequency(
    te_expand = te_filter_expand,
    output_dir = plot_dir,
    output_prefix = "germline"
  )

  #### HOSTSEQ FILTER SENSITIVITY ANALYSIS ####
  cat("\n========================================\n")
  cat("STEP 3.9: HostSeq filter sensitivity analysis\n")
  cat("========================================\n")
  cat("Analyzing sensitivity of common TE filtering to HostSeq sample size...\n")

  sensitivity_plot <- analyze_hostseq_filter_sensitivity(
    te_data = te_raw_prepped,
    rare_gnomad = 3,
    rare_hostseq = 3,
    sample_sizes = NULL,  # Uses default: 10%, 20%, ..., 100%
    output_dir = plot_dir,
    output_prefix = "germline",
    seed = 123
  )
} else {
  cat("\n========================================\n")
  cat("STEP 3.6-3.9: Skipping frequency plots (GENERATE_FREQUENCY_PLOTS = FALSE)\n")
  cat("========================================\n")
}

##### PROCESS RARE TE DATA ####
cat("\n========================================\n")
cat("STEP 4: Processing rare TE data\n")
cat("========================================\n")

# Rare TE: filter common, create count matrix, merge with clinical
# Always capture process_te_data_germline output for processed_output.txt
process_te_stdout <- capture.output({
  final_te_count <- process_te_data_germline(te_raw_labeled, clinical, metrics, ancestry,
                                            apply_filter_common = TRUE, rare_gnomad = 3, rare_hostseq = 3,
                                            split_by_gene = FALSE, apply_process_combinations = TRUE)
})
write_stdout_to_file(process_te_stdout, paste0(r_dir_files, "/processed_output.txt"),
                     "STDOUT for process_te_data_germline (final_te_count)", append = TRUE)
# Print the output to console as well
cat(paste(process_te_stdout, collapse = "\n"), "\n")
save(final_te_count, file = paste0(r_dir, "final_te_count_rare.RData"))

# Rare TE: expanded format (one row per full TE, no count matrix)
invisible(capture.output({
  suppressMessages({
    final_te_count_expand <- process_te_data_germline(te_raw_labeled, clinical, metrics, ancestry,
                                                    apply_filter_common = TRUE, rare_gnomad = 3, rare_hostseq = 3,
                                                    split_by_gene = FALSE, apply_process_combinations = FALSE)
  })
}))
save(final_te_count_expand, file = paste0(r_dir, "final_te_count_expand_rare.RData"))

# Rare TE: split by gene format
# COMMENTED OUT: Using split data from annotSV for better gene annotations
# suppressMessages({
#   final_te_count_split <- process_te_data_germline(te_germline_filter, clinical, metrics, ancestry,
#                                                  apply_filter_common = TRUE, rare_gnomad = 3, rare_hostseq = 3,
#                                                  split_by_gene = TRUE, apply_process_combinations = FALSE)
# })
# save(final_te_count_split, file = paste0(r_dir, "final_te_count_split_rare.RData"))

#### PROCESS COMMON TE DATA ####
if (PROCESS_COMMON_TES) {
  cat("\n========================================\n")
  cat("STEP 5: Processing common TE data\n")
  cat("========================================\n")
  cat("Filtering common TEs (no frequency filters)...\n")

  # Common TE: keep common, create count matrix
  invisible(capture.output({
    suppressMessages({
      final_te_count_common <- process_te_data_germline(te_raw_labeled, clinical, metrics, ancestry,
                                                      apply_filter_common = FALSE, split_by_gene = FALSE,
                                                      apply_process_combinations = TRUE)
    })
  }))
  save(final_te_count_common, file = paste0(r_dir, "final_te_count_common.RData"))

  # Common TE: expanded format (one row per full TE, no count matrix)
  invisible(capture.output({
    suppressMessages({
      final_te_count_expand_common <- process_te_data_germline(te_raw_labeled, clinical, metrics, ancestry,
                                                             apply_filter_common = FALSE, split_by_gene = FALSE,
                                                             apply_process_combinations = FALSE)
    })
  }))
  save(final_te_count_expand_common, file = paste0(r_dir, "final_te_count_expand_common.RData"))

  # Common TE: split by gene format
  invisible(capture.output({
    suppressMessages({
      final_te_count_split_common <- process_te_data_germline(te_raw_labeled, clinical, metrics, ancestry,
                                                            apply_filter_common = FALSE, split_by_gene = TRUE,
                                                            apply_process_combinations = FALSE)
    })
  }))
  save(final_te_count_split_common, file = paste0(r_dir, "final_te_count_split_common.RData"))
} else {
  cat("\n========================================\n")
  cat("STEP 5: Skipping common TE processing (PROCESS_COMMON_TES = FALSE)\n")
  cat("========================================\n")
}

#### PROCESS SPLIT TE DATA ####
cat("\n========================================\n")
cat("STEP 6: Processing split TE data for gene analysis\n")
cat("========================================\n")
cat("Processing gene-split TE data...\n")

# Filter TE types for split data (suppress output)
suppressMessages({
  te_split <- filter_te_types(te_split, include_alu = INCLUDE_ALU, include_sva = INCLUDE_SVA)
})

# Prepare te_split (suppress output)
suppressMessages({
  te_split_prepped <- prep_te(te_split, nonproband, noconsent, hostseq_cancer, metrics,
                             c("mean_cov", "avg_quality", "pct_chimeras"), c(20, 25, 2), type="N")
})

# Add HostSeq labels to split data using the same split
te_split_prepped$hostseq_group <- case_when(
  te_split_prepped$sample %in% hostseq_split$filter_samples ~ "filter",
  te_split_prepped$sample %in% hostseq_split$analysis_samples ~ "analysis",
  TRUE ~ NA_character_
)

cat("Split data prepared with", nrow(te_split_prepped), "TEs (with labeled HostSeq groups)\n")

# Split TE data for gene analysis - filtering will handle removing filter group automatically
invisible(capture.output({
  suppressMessages({
    final_te_count_expand_split <- process_te_data_germline(te_split_prepped, clinical, metrics, ancestry,
                                                          apply_filter_common = TRUE, rare_gnomad = 3, rare_hostseq = 3,
                                                          apply_process_combinations = FALSE)
  })
}))
save(final_te_count_expand_split, file = paste0(r_dir, "final_te_count_expand_split.RData"))

# Write pipeline summary tables to processed output file
pipeline_summary <- capture.output({
  cat("\n" , rep("=", 80), "\n")
  cat("DETAILED PIPELINE TRACKING SUMMARY - GERMLINE\n")
  cat(rep("=", 80), "\n")

  # TE Table with detailed reasons
  cat("\nTE PROCESSING PIPELINE:\n")
  cat(sprintf("%-6s %-10s %-10s %-8s %s\n", "Step", "Count", "Lost", "Lost%", "Detailed Reason"))
  cat(rep("-", 100), "\n")

  if (ENABLE_DETAILED_TRACKING && exists("step0_raw_te_count") && exists("step1_te_count") && exists("step2_te_count") && exists("step3_te_count")) {
    cat(sprintf("%-6s %-10d %-10s %-8s %s\n", "0", step0_raw_te_count, "-", "-", "Raw TE insertions from annotSV pipeline"))
    cat(sprintf("%-6s %-10d %-10d %-8.1f %s\n", "1", step1_te_count, step1_te_loss, (step1_te_loss/step0_raw_te_count)*100,
        "TE type filtering (germline data includes all TE types)"))
    cat(sprintf("%-6s %-10d %-10d %-8.1f %s\n", "2", step2_te_count, step2_te_loss, (step2_te_loss/step1_te_count)*100,
        "No false positive filtering (germline has no IGV review)"))
    cat(sprintf("%-6s %-10d %-10d %-8.1f %s\n", "3", step3_te_count, step3_te_loss, (step3_te_loss/step2_te_count)*100,
        "Quality control: Removed TEs from samples failing metrics"))
    cat(sprintf("    %s\n", "(mean_cov <20, avg_quality <30, pct_chimeras >2), consent, or"))
    cat(sprintf("    %s\n", "inappropriate sample type filters"))
  } else {
    cat("Detailed step tracking disabled (set ENABLE_DETAILED_TRACKING = TRUE to enable)\n")
  }

  # Sample Table with detailed reasons
  cat("\nSAMPLE PROCESSING PIPELINE:\n")
  cat(sprintf("%-6s %-10s %-10s %-8s %s\n", "Step", "Count", "Lost", "Lost%", "Detailed Reason"))
  cat(rep("-", 100), "\n")

  if (ENABLE_DETAILED_TRACKING && exists("step1_sample_count") && exists("step2_sample_count") && exists("step3_sample_count")) {
    cat(sprintf("%-6s %-10d %-10s %-8s %s\n", "1", step1_sample_count, "-", "-",
        "Unique samples with at least one TE insertion called"))
    cat(sprintf("%-6s %-10d %-10d %-8.1f %s\n", "2", step2_sample_count, step2_sample_loss, (step2_sample_loss/step1_sample_count)*100,
        "No false positive filtering (germline)"))
    cat(sprintf("%-6s %-10d %-10d %-8.1f %s\n", "3", step3_sample_count, step3_sample_loss, (step3_sample_loss/step2_sample_count)*100,
        "Quality control: Excluded samples with poor sequencing"))
    cat(sprintf("    %s\n", "metrics, no consent, or inappropriate sample type"))
  } else {
    cat("Detailed step tracking disabled (set ENABLE_DETAILED_TRACKING = TRUE to enable)\n")
  }

  cat("\n", rep("=", 80), "\n\n")
})

# Write to file and console
write_stdout_to_file(pipeline_summary, paste0(r_dir_files, "/processed_output.txt"),
                     "DETAILED PIPELINE TRACKING SUMMARY - GERMLINE", append = TRUE)
cat(pipeline_summary, sep = "\n")

#### EXPORT CSV FILES ####
cat("\n========================================\n")
cat("Exporting CSV files\n")
cat("========================================\n")
# Export split dataframes to CSV (germline has no _selected version)
write.csv(final_te_count_expand_split$te_aff, paste0(r_dir_files, "te_aff_split.csv"), row.names = FALSE, quote = FALSE)
cat("Exported te_aff_split.csv to files directory\n")

write.csv(final_te_count_expand_split$te_aff_unaff, paste0(r_dir_files, "te_aff_unaff_split.csv"), row.names = FALSE, quote = FALSE)
cat("Exported te_aff_unaff_split.csv to files directory\n")

write.csv(final_te_count_expand_split$te_lfs, paste0(r_dir_files, "te_lfs_split.csv"), row.names = FALSE, quote = FALSE)
cat("Exported te_lfs_split.csv to files directory\n")

write.csv(final_te_count_expand_split$te_kics_hostseq, paste0(r_dir_files, "te_kics_hostseq_split.csv"), row.names = FALSE, quote = FALSE)
cat("Exported te_kics_hostseq_split.csv to files directory\n")

# Export to HPC directory
cat("\n========================================\n")
cat("Exporting to HPC directory...\n")
cat("========================================\n")

# Export split format (one row per TE-gene overlap)
cat("\nExporting split format files...\n")
write.csv(final_te_count_expand_split$te_aff, paste0(hpc_dir, "te_aff_split.csv"), row.names = FALSE, quote = FALSE)
write.csv(final_te_count_expand_split$te_aff_unaff, paste0(hpc_dir, "te_aff_unaff_split.csv"), row.names = FALSE, quote = FALSE)
write.csv(final_te_count_expand_split$te_lfs, paste0(hpc_dir, "te_lfs_split.csv"), row.names = FALSE, quote = FALSE)
write.csv(final_te_count_expand_split$te_kics_hostseq, paste0(hpc_dir, "te_kics_hostseq_split.csv"), row.names = FALSE, quote = FALSE)
cat("✓ Exported split format files\n")

# Export count matrix format (one row per sample with TE counts)
cat("\nExporting count matrix format files...\n")
write.csv(final_te_count$te_all, paste0(hpc_dir, "te_all.csv"), row.names = FALSE, quote = FALSE)
write.csv(final_te_count$te_aff, paste0(hpc_dir, "te_aff.csv"), row.names = FALSE, quote = FALSE)
write.csv(final_te_count$te_lfs, paste0(hpc_dir, "te_lfs.csv"), row.names = FALSE, quote = FALSE)
write.csv(final_te_count$te_kics, paste0(hpc_dir, "te_kics.csv"), row.names = FALSE, quote = FALSE)
write.csv(final_te_count$te_taylor, paste0(hpc_dir, "te_taylor.csv"), row.names = FALSE, quote = FALSE)
write.csv(final_te_count$te_hostseq, paste0(hpc_dir, "te_hostseq.csv"), row.names = FALSE, quote = FALSE)
cat("✓ Exported count matrix format files\n")

# Export expand format (one row per complete TE insertion)
cat("\nExporting expand format files...\n")
write.csv(final_te_count_expand$te_all, paste0(hpc_dir, "te_all_expand.csv"), row.names = FALSE, quote = FALSE)
write.csv(final_te_count_expand$te_aff, paste0(hpc_dir, "te_aff_expand.csv"), row.names = FALSE, quote = FALSE)
write.csv(final_te_count_expand$te_lfs, paste0(hpc_dir, "te_lfs_expand.csv"), row.names = FALSE, quote = FALSE)
write.csv(final_te_count_expand$te_kics, paste0(hpc_dir, "te_kics_expand.csv"), row.names = FALSE, quote = FALSE)
write.csv(final_te_count_expand$te_taylor, paste0(hpc_dir, "te_taylor_expand.csv"), row.names = FALSE, quote = FALSE)
write.csv(final_te_count_expand$te_hostseq, paste0(hpc_dir, "te_hostseq_expand.csv"), row.names = FALSE, quote = FALSE)

# Create and export derived expand datasets
te_aff_unaff_expand <- final_te_count_expand$te_all %>%
  filter(cohort != "Taylor" & cohort != "HostSeq")
write.csv(te_aff_unaff_expand, paste0(hpc_dir, "te_aff_unaff_expand.csv"), row.names = FALSE, quote = FALSE)

cat("\nCreating KICS + HostSeq combined dataset:\n")
cat("  KICS samples:", nrow(final_te_count_expand$te_kics), "rows\n")
cat("  HostSeq samples:", nrow(final_te_count_expand$te_hostseq), "rows\n")
if (nrow(final_te_count_expand$te_kics) > 0) {
  cat("  KICS unique samples:", length(unique(final_te_count_expand$te_kics$sample)), "\n")
}
if (nrow(final_te_count_expand$te_hostseq) > 0) {
  cat("  HostSeq unique samples:", length(unique(final_te_count_expand$te_hostseq$sample)), "\n")
}
te_kics_hostseq_expand <- rbind(
  final_te_count_expand$te_kics,
  final_te_count_expand$te_hostseq
)
cat("  Combined total:", nrow(te_kics_hostseq_expand), "rows\n")
write.csv(te_kics_hostseq_expand, paste0(hpc_dir, "te_kics_hostseq_expand.csv"), row.names = FALSE, quote = FALSE)
cat("✓ Exported expand format files\n")

# Extract nohits from final_te_count (samples with total == 0)
cat("\nExporting nohits samples...\n")
nohits_aff <- final_te_count$te_aff[final_te_count$te_aff$total == 0, "sample", drop = FALSE]
write.table(nohits_aff, paste0(hpc_dir, "nohits_te_aff.csv"), row.names = FALSE, col.names = FALSE, sep = ",", quote = FALSE)
cat("✓ Exported nohits samples\n")

cat("\n✓ All files exported to HPC directory:", hpc_dir, "\n")

cat("\n========================================\n")
cat("PROCESSING COMPLETE\n")
cat("========================================\n")
cat("All R objects saved to:", r_dir, "\n\n")

cat("Main output files:\n")
cat("  • Preprocessed data: te_raw_prepped.RData\n")
cat("  • Rare TEs (count matrix): final_te_count_rare.RData\n")
cat("  • Rare TEs (expanded): final_te_count_expand_rare.RData\n")
cat("  • Rare TEs (gene-split): final_te_count_split_rare.RData\n")
cat("  • Common TEs (count matrix): final_te_count_common.RData\n")
cat("  • Common TEs (expanded): final_te_count_expand_common.RData\n")
cat("  • Common TEs (gene-split): final_te_count_split_common.RData\n")
cat("  • Gene-split expanded: final_te_count_expand_split.RData\n")
cat("  • Clinical/metrics/ancestry: clinical.RData, metrics.RData, ancestry.RData\n")
cat("  • CSV exports: te_aff.csv, te_all_all.csv\n\n")

if (ENABLE_DETAILED_TRACKING) {
  cat("Detailed processing logs written to:\n")
  cat("  • ", r_dir_files, "processed_output.txt\n\n")
}
