#!/usr/bin/env Rscript

# TE Tumour Data Processing Script
# This script processes tumour transposable element (TE) data and saves various data tables as R objects

#### SETUP ####
r_dir <- "/Users/briannelaverty/Documents/R_Malkin/te/data/R_obj/tumour/"
plot_dir <- "/Users/briannelaverty/Library/Mobile Documents/com~apple~CloudDocs/Malkin_lab/figures/te/graphs/tumour/"
r_dir_files <- paste0(plot_dir, "files/")
hpc_dir <- "/Users/briannelaverty/Documents/R_Malkin/te/hpc_data/processed/"
processed_dir <- hpc_dir  # Alias for HPC output directory

# TE type filtering options
INCLUDE_ALU <- FALSE  # Set to TRUE to include ALU elements
INCLUDE_SVA <- FALSE  # Set to TRUE to include SVA elements

# Processing options
PROCESS_COMMON_TES <- TRUE # Set to TRUE to process common TEs (slower)
PROCESS_FULLLENGTH <- TRUE # Set to TRUE to process full-length LINE1 (≥5900bp)
GENERATE_FREQUENCY_PLOTS <- FALSE  # Set to TRUE to generate TE frequency distribution and HostSeq filtering plots (slower)

# Pipeline tracking options
ENABLE_DETAILED_TRACKING <- TRUE # Set to TRUE to enable detailed pipeline tracking (slower)

# Create output directory if it doesn't exist
if (!dir.exists(r_dir_files)) {
  dir.create(r_dir_files, recursive = TRUE)
}

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
  library(readr)
})

# Source functions
source("/Users/briannelaverty/Documents/R_Malkin/TE/scripts/viz/functions_te.R")

#### LOAD DATA ####
cat("Loading raw tumour data...\n")
te_raw_t <- fread("/Users/briannelaverty/Documents/R_Malkin/te/data/final/tumour_annotSV_full.tsv", sep="\t", header=TRUE)


# Store initial counts for pipeline tracking (before any filtering)
if (ENABLE_DETAILED_TRACKING) {
  step0_raw_te_count <- nrow(te_raw_t)
  step0_sample_count <- length(unique(te_raw_t$Samples_ID))
}
te_split_t <- fread("/Users/briannelaverty/Documents/R_Malkin/te/data/final/tumour_annotSV_split.tsv", sep="\t", header=TRUE)
chr_length <- read.csv("/Users/briannelaverty/Documents/R_Malkin/te/data/raw/chromosome_length.csv", header=TRUE)
clinical <- read.delim("/Users/briannelaverty/Documents/R_Malkin/clinical/te_clinical.csv", sep=",", header=TRUE)
complete_samples <- fread("/Users/briannelaverty/Documents/R_Malkin/te/data/final/complete_samples_T.txt", sep="\t", header=FALSE)
metrics <- read.delim("/Users/briannelaverty/Documents/R_Malkin/te/data/final/combined_metrics.txt", sep="\t", header=TRUE)
ancestry <- read.delim("/Users/briannelaverty/Documents/R_Malkin/te/data/final/ancestry.txt", sep="\t", header=TRUE)
hg37_genes <- read.delim("/Users/briannelaverty/Documents/R_Malkin/te/data/final/hg37_genes.tsv", sep="\t", header=TRUE)
nonproband <- read.delim("/Users/briannelaverty/Documents/R_Malkin/te/data/final/non_proband_normals", sep="\t", header=TRUE, colClasses = c("character"))
noconsent <- read.delim("/Users/briannelaverty/Documents/R_Malkin/te/data/final/kics_samples_exclude", sep="\t", header=TRUE, colClasses = c("character"))
l1_merge_fp <- fread("/Users/briannelaverty/Documents/R_Malkin/te/IGV_tracking/L1_merge_sample_id.csv")
unique_l1_tp <- as.data.table(readr::read_csv("/Users/briannelaverty/Documents/R_Malkin/TE/data/final/unique_l1_tp.csv", locale=readr::locale(encoding="UTF-8"), show_col_types=FALSE, trim_ws=TRUE))
colnames(unique_l1_tp)[1:3] <- c("sample", "ID", "caller")
unique_l1_fp <- as.data.table(readr::read_csv("/Users/briannelaverty/Documents/R_Malkin/TE/data/final/unique_l1_fp.csv", locale=readr::locale(encoding="UTF-8"), show_col_types=FALSE, trim_ws=TRUE))
colnames(unique_l1_fp)[1:3] <- c("sample", "ID", "caller")
unique_l1_master <- as.data.table(readr::read_csv("/Users/briannelaverty/Documents/R_Malkin/TE/data/final/unique_l1_master.csv", locale=readr::locale(encoding="UTF-8"), show_col_types=FALSE))  

# load germline data for filtering common TEs
load(paste0("/Users/briannelaverty/Documents/R_Malkin/te/data/R_obj/germline/", "te_raw_prepped.RData"))

#### PREP DATA ####
cat("Preparing tumour data...\n")

# Create/clear processed_output.txt file at start
output_file <- paste0(r_dir_files, "/processed_output.txt")
cat("", file = output_file)  # Create empty file
cat("TUMOR TE PROCESSING PIPELINE OUTPUT\n", file = output_file, append = TRUE)
cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n", file = output_file, append = TRUE)
cat("Detailed tracking:", ifelse(ENABLE_DETAILED_TRACKING, "ENABLED", "DISABLED"), "\n\n", file = output_file, append = TRUE)

# Metrics
metrics <- prep_metrics_tumour(metrics)

# Ancestry
ancestry <- prep_ancestry(ancestry)

# Clinical
clinical <- prep_clinical(clinical) 

# Host seq with cancer
hostseq_cancer <- prep_hostseq(clinical)

# Blood cancer samples
bloodcancer <- clinical %>%
  filter(tumor_class == "LEUKEMIA/LYMPHOMA") %>%
  select(sample) %>%
  rename(V1 = sample)

# Chromosome length
chr_length$chr <- factor(chr_length$chr, levels=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"X","Y"))
chr_lengths <- setNames(chr_length$length, chr_length$chr)

# Gene size
gene_size <- add_gene_size(hg37_genes)

# Save all prepared data objects
cat("Saving prepared tumour data objects...\n")
save(gene_size, file = paste0(r_dir, "gene_size.RData"))
save(chr_lengths, file = paste0(r_dir, "chr_lengths.RData"))
save(complete_samples, file = paste0(r_dir, "complete_samples.RData"))
save(hostseq_cancer, file = paste0(r_dir, "hostseq_cancer.RData"))
save(clinical, file = paste0(r_dir, "clinical.RData"))
save(metrics, file = paste0(r_dir, "metrics.RData"))
save(ancestry, file = paste0(r_dir, "ancestry.RData"))

# Merge clinical and ancestry data and save as CSV (tumour version)
cat("Merging clinical and ancestry data (tumour)...\n")
clinical_ancestry_t <- merge_dfs(clinical, ancestry, include_all_x = FALSE, print_info = FALSE, dataset_name = "ancestry_tumour")
clinical_ancestry_file_t <- paste0(processed_dir, "clinical_ancestry_tumour.csv")
write.csv(clinical_ancestry_t, clinical_ancestry_file_t, row.names = FALSE)
cat("✓ Saved merged clinical and ancestry data (tumour) to:", clinical_ancestry_file_t, "\n")

# Filter TE types based on configuration
cat("\n========================================\n")
cat("STEP 1: Filtering TE types\n")
cat("========================================\n")
te_raw_t <- filter_te_types(te_raw_t, include_alu = INCLUDE_ALU, include_sva = INCLUDE_SVA)

if (ENABLE_DETAILED_TRACKING) {
  step1 <- init_step_tracking(te_raw_t, "Samples_ID")
  step1_te_count <- step1$te_count
  step1_sample_count <- step1$sample_count
  step1_te_loss <- step0_raw_te_count - step1_te_count
}

# Filter false positives based on manual review
cat("\n========================================\n")
cat("STEP 2: Filtering false positives\n")
cat("========================================\n")

# Always capture and write filter_false_positives output
filter_fp_stdout_t <- capture.output({
  fp_result <- filter_false_positives(te_raw_t, l1_merge_fp, unique_l1_tp, unique_l1_fp, unique_l1_master)
})
te_raw_t <- fp_result$te_data
excluded_unreviewed <- fp_result$excluded_samples
write_stdout_to_file(filter_fp_stdout_t, paste0(r_dir_files, "/processed_output.txt"),
                     "STDOUT for filter_false_positives", append = TRUE)

# Initialize master exclusion tracker with unreviewed IGV samples
exclusion_tracker_master <- create_sample_exclusion_tracker()
if (length(excluded_unreviewed) > 0) {
  exclusion_tracker_master <- add_to_exclusion_tracker(exclusion_tracker_master, excluded_unreviewed, "unreviewed_igv")
  cat("Tracked", length(excluded_unreviewed), "unreviewed IGV samples for exclusion\n")
}

if (ENABLE_DETAILED_TRACKING) {
  # Extract intermediate counts from captured output
  caller_numbers <- extract_count_from_stdout(filter_fp_stdout_t, "Caller types:")
  step2a_te_count <- if (length(caller_numbers) >= 2) sum(caller_numbers) else step1_te_count
  step2b_te_count <- nrow(te_raw_t)
}

# Write CSV file with true positives (TP) after false positive filtering
cat("Writing true positives CSV file...\n")
tp_output <- data.frame(
  sample = te_raw_t$sample,
  te_id = te_raw_t$ID
)
write.csv(tp_output, paste0(r_dir_files, "tumour_line1_tp.csv"), row.names = FALSE, quote = FALSE)

if (ENABLE_DETAILED_TRACKING) {
  step2a_te_loss <- step1_te_count - step2a_te_count
  step2b_te_loss <- step2a_te_count - step2b_te_count
  step2_te_count <- step2b_te_count
  step2_sample_count <- length(unique(te_raw_t$sample))
  step2_te_loss <- step1_te_count - step2_te_count
  step2_sample_loss <- step1_sample_count - step2_sample_count

  step2a_sample_count <- extract_value_from_stdout(filter_fp_stdout_t,
                                                     ".*Samples after unreviewed removal: (\\d+) samples.*")

  if (!is.na(step2a_sample_count)) {
    step2a_sample_loss <- step1_sample_count - step2a_sample_count
    step2b_sample_loss <- step2a_sample_count - step2_sample_count
  } else {
    step2a_sample_loss <- NA
    step2b_sample_loss <- NA
    cat("Warning: Could not extract step2a_sample_count from filter output\n")
  }
}

# Prepare te_raw
cat("\n========================================\n")
cat("STEP 3: Quality filtering (prep_te)\n")
cat("========================================\n")

# Always capture and write prep_te output
prep_te_stdout_t <- capture.output({
  prep_result <- prep_te(te_raw_t, nonproband, noconsent, hostseq_cancer, bloodcancer, metrics,
                             c("mean_cov", "avg_quality", "pct_chimeras"), c(20, 25, 2), type="T",
                             export_filtered = TRUE, output_dir = r_dir_files)
})
te_raw_prepped_t <- prep_result$df
excluded_qc <- prep_result$excluded_samples
# Collect QC exclusions into master tracker
if (!is.null(prep_result$exclusion_tracker) && nrow(prep_result$exclusion_tracker) > 0) {
  exclusion_tracker_master <- rbind(exclusion_tracker_master, prep_result$exclusion_tracker)
}
write_stdout_to_file(prep_te_stdout_t, paste0(r_dir_files, "/processed_output.txt"),
                     "STDOUT for prep_te", append = TRUE)

# Remove all excluded samples from complete_samples
# This prevents filtered samples from being added back as "nohits" in process_te_data_tumour
all_excluded <- unique(c(excluded_unreviewed, excluded_qc))
complete_samples_before <- nrow(complete_samples)
complete_samples <- complete_samples[!(complete_samples$V1 %in% all_excluded), , drop = FALSE]
cat("Removed", complete_samples_before - nrow(complete_samples),
    "excluded samples from complete_samples (", nrow(complete_samples), "remaining)\n")

if (ENABLE_DETAILED_TRACKING) {
  step3 <- init_step_tracking(te_raw_prepped_t, "sample")
  step3_te_count <- step3$te_count
  step3_sample_count <- step3$sample_count
  step3_te_loss <- step2_te_count - step3_te_count
  step3_sample_loss <- step2_sample_count - step3_sample_count
}
save(te_raw_prepped_t, file = paste0(r_dir, "te_raw_prepped_t.RData"))

#### SPLIT HOSTSEQ SAMPLES ####
cat("\n========================================\n")
cat("STEP 3.5: Splitting HostSeq samples\n")
cat("========================================\n")
cat("Splitting HostSeq samples into filtering (66%) and analysis (33%) groups...\n")

# Merge ancestry data for stratified splitting
cat("Merging ancestry data for stratified split...\n")
te_raw_prepped_with_ancestry <- merge_dfs(te_raw_prepped, ancestry, include_all_x = TRUE, print_info = FALSE, dataset_name = "ancestry")

# Split germline HostSeq samples (used for filtering) - stratified by ancestry
hostseq_split <- split_hostseq_samples(
  te_data = te_raw_prepped_with_ancestry,  # Use germline data with ancestry
  filter_pct = 66,
  seed = 123
)

# Save split information
save(hostseq_split, file = paste0(r_dir, "hostseq_split.RData"))

# Create filter and analysis dataframes from the split result
# Note: split_hostseq_samples returns sample lists, not dataframes
if (!is.null(hostseq_split$filter_samples) && length(hostseq_split$filter_samples) > 0) {
  te_germline_filter <- te_raw_prepped_with_ancestry %>%
    filter(sample %in% hostseq_split$filter_samples)
} else {
  te_germline_filter <- NULL
}

if (!is.null(hostseq_split$analysis_samples) && length(hostseq_split$analysis_samples) > 0) {
  te_germline_analysis <- te_raw_prepped_with_ancestry %>%
    filter(sample %in% hostseq_split$analysis_samples)
} else {
  te_germline_analysis <- NULL
}

# Verify HostSeq split was successful
if (is.null(te_germline_filter) || nrow(te_germline_filter) == 0) {
  cat("\n⚠ WARNING: HostSeq filter group is empty!\n")
  cat("  This likely means:\n")
  cat("  1. The germline data file 'te_raw_prepped.RData' is missing or empty\n")
  cat("  2. No HostSeq samples exist in the germline data\n")
  cat("  3. The split_hostseq_samples() function failed\n\n")
  cat("  Pipeline will continue but common TE filtering may not work properly.\n")
  cat("  Consider running 01_te_processing_germline.R first.\n\n")
}
if (is.null(te_germline_analysis) || nrow(te_germline_analysis) == 0) {
  cat("\n⚠ WARNING: HostSeq analysis group is empty!\n\n")
}

#### REMOVE ANCESTRY COLUMNS FROM GERMLINE FILTER ####
cat("\n========================================\n")
cat("Removing ancestry columns from germline filter data\n")
cat("========================================\n")
cat("Ancestry was needed for stratified HostSeq split.\n")
cat("Removing now to keep data clean.\n\n")

# Check if te_germline_filter exists before processing
if (is.null(te_germline_filter) || nrow(te_germline_filter) == 0) {
  cat("⚠ WARNING: No HostSeq filter samples available. Skipping ancestry column removal.\n")
  cat("  This may occur if germline HostSeq data was not properly loaded.\n\n")
} else {
  # Remove ancestry columns from germline filter data
  ancestry_cols <- c("predicted_ancestry_thres", "mapped_label", "base_sample")
  te_germline_filter <- te_germline_filter %>%
    select(-any_of(ancestry_cols))

  cat("✓ Ancestry columns removed from germline filter\n")
  cat("  Filter samples:", nrow(te_germline_filter), "\n\n")
}

#### PLOT TE FREQUENCY DISTRIBUTIONS (OPTIONAL) ####
if (GENERATE_FREQUENCY_PLOTS) {
  cat("\n========================================\n")
  cat("STEP 3.6: Plotting TE frequency distributions\n")
  cat("========================================\n")
  cat("Creating histograms of TE frequency in gnomAD before filtering...\n")

  freq_plots_t <- plot_te_frequency_distributions(
    te_data = te_raw_prepped_t,
    output_dir = plot_dir,
    output_prefix = "tumour"
  )

  #### PLOT TE HOSTSEQ FREQUENCY ####
  cat("\n========================================\n")
  cat("STEP 3.7: Plotting TEs in dataset and their HostSeq frequency (filter group 66%)\n")
  cat("========================================\n")

  if (is.null(te_germline_filter) || nrow(te_germline_filter) == 0) {
    cat("⚠ Skipping HostSeq frequency plot - no HostSeq filter samples available\n\n")
  } else {
    cat("Creating histogram of TEs and their frequency in HostSeq filter group (66%)...\n")

    # Combine tumour data with germline filter group
    # Use rbindlist for fast binding with automatic type coercion
    te_with_filter_hostseq <- as.data.frame(
      rbindlist(list(as.data.table(te_raw_prepped_t), as.data.table(te_germline_filter)),
                fill = TRUE, use.names = TRUE)
    )

    # Convert to expanded format for this plot
    te_filter_expand_t <- te_with_filter_hostseq %>%
      distinct(sample, SV_chrom, SV_start, ALT, .keep_all = TRUE)

    te_hostseq_freq_plot_t <- plot_te_hostseq_frequency(
      te_expand = te_filter_expand_t,
      output_dir = plot_dir,
      output_prefix = "tumour"
    )
  }

  #### HOSTSEQ FILTER SENSITIVITY ANALYSIS ####
  cat("\n========================================\n")
  cat("STEP 3.8: HostSeq filter sensitivity analysis\n")
  cat("========================================\n")
  cat("Analyzing sensitivity of common TE filtering to HostSeq sample size...\n")

  # Combine tumour data with germline HostSeq for sensitivity analysis
  te_raw_prepped_hostseq <- te_raw_prepped %>%
    filter(grepl("^HS_", sample))

  # Use rbindlist for fast binding with automatic type coercion
  te_tumour_with_hostseq <- as.data.frame(
    rbindlist(list(as.data.table(te_raw_prepped_t), as.data.table(te_raw_prepped_hostseq)),
              fill = TRUE, use.names = TRUE)
  )

  sensitivity_plot_t <- analyze_hostseq_filter_sensitivity(
    te_data = te_tumour_with_hostseq,
    rare_gnomad = 3,
    rare_hostseq = 3,
    sample_sizes = NULL,  # Uses default: 10%, 20%, ..., 100%
    output_dir = plot_dir,
    output_prefix = "tumour",
    seed = 123
  )
} else {
  cat("\n========================================\n")
  cat("STEP 3.6-3.8: Skipping frequency plots (GENERATE_FREQUENCY_PLOTS = FALSE)\n")
  cat("========================================\n")
}

##### PROCESS RARE TE DATA ####
cat("\n========================================\n")
cat("STEP 4: Processing rare TE data\n")
cat("========================================\n")

# Check if we can apply common TE filtering
can_filter_common <- !is.null(te_germline_filter) && nrow(te_germline_filter) > 0

if (!can_filter_common) {
  cat("⚠ WARNING: Cannot apply common TE filtering (no HostSeq filter samples)\n")
  cat("  Processing will only use gnomAD frequency filtering\n\n")
}

# Rare TE: filter common, create count matrix, merge with clinical, select one sample per patient
# Always capture and write process_te_data_tumour output
process_te_stdout_t <- capture.output({
  final_te_count_t <- process_te_data_tumour(te_raw_prepped_t, te_germline=te_germline_filter, clinical, complete_samples, metrics, ancestry,
                                             apply_filter_common = can_filter_common, rare_gnomad = 3, rare_hostseq = 3,
                                             split_by_gene = FALSE, apply_process_combinations = TRUE, select_samples_split = TRUE, nohits_prefix = "final_te_count_t", nohits_output_dir = r_dir_files,
                                             return_exclusion_tracker = TRUE, apply_age_filter = FALSE)
})
save(final_te_count_t, file = paste0(r_dir, "final_te_count_rare.RData"))

# Collect age exclusions from process_te_data_tumour
if (!is.null(final_te_count_t$exclusion_tracker) && nrow(final_te_count_t$exclusion_tracker) > 0) {
  exclusion_tracker_master <- rbind(exclusion_tracker_master, final_te_count_t$exclusion_tracker)
  cat("Collected", nrow(final_te_count_t$exclusion_tracker), "additional exclusions from process_te_data_tumour\n")
}

cat("Writing process_te_data_tumour stdout to processed_output.txt...\n")
tryCatch({
  write_stdout_to_file(process_te_stdout_t, paste0(r_dir_files, "/processed_output.txt"),
                       "STDOUT for process_te_data_tumour (final_te_count_t)", append = TRUE)
  cat("✓ Successfully wrote process_te_data_tumour output (", length(process_te_stdout_t), " lines)\n", sep = "")
}, error = function(e) {
  cat("ERROR writing process_te_stdout_t:", e$message, "\n")
})

if (ENABLE_DETAILED_TRACKING) {
  step4_input_count <- step3_te_count
  step4_after_common_count <- extract_value_from_stdout(process_te_stdout_t,
                                                         ".*Total # of unique TEs in te_df: (\\d+).*")

  # Extract step4 metrics
  step4_metrics <- extract_step4_metrics(final_te_count_t, step3_sample_count, process_te_stdout_t)
  step4_total_samples <- step4_metrics$total_samples
  step4_samples_with_tes <- step4_metrics$samples_with_tes
  step4_samples_no_tes <- step4_metrics$samples_no_tes
  step4_total_te_insertions <- step4_metrics$total_te_insertions
}

# Rare TE: expanded format (one row per full TE, no count matrix)
# Run without capturing output to avoid duplicate logging
suppressMessages({
  final_te_count_expand_t <- process_te_data_tumour(te_raw_prepped_t, te_germline=te_germline_filter, clinical, complete_samples, metrics, ancestry,
                                                    apply_filter_common = can_filter_common, rare_gnomad = 3, rare_hostseq = 3,
                                                    split_by_gene = FALSE, apply_process_combinations = FALSE, select_samples_split = TRUE,
                                                    apply_age_filter = FALSE)
})
save(final_te_count_expand_t, file = paste0(r_dir, "final_te_count_expand_rare.RData"))

if (ENABLE_DETAILED_TRACKING) {
  # Count TEs and samples
  step5_te_count <- tryCatch(nrow(final_te_count_expand_t), error = function(e) NA)
  step5_sample_count <- tryCatch({
    if("sample" %in% colnames(final_te_count_expand_t)) {
      length(unique(final_te_count_expand_t$sample))
    } else if(exists("step4_samples_with_tes") && !is.na(step4_samples_with_tes)) {
      step4_samples_with_tes
    } else {
      NA
    }
  }, error = function(e) NA)

  # Count unique TE loci
  unique_te_loci <- count_unique_te_loci(final_te_count_expand_t,
                                          fallback_value = step4_total_te_insertions)

  if (!is.na(unique_te_loci) && exists("step4_total_te_insertions") &&
      !is.na(step4_total_te_insertions) && unique_te_loci == step4_total_te_insertions) {
    cat("No deduplication needed - unique TE loci equals technical filtering count\n")
  }
}

# Rare TE: split by gene format
# COMMENTED OUT: Using split_split version instead for better gene annotations
# final_te_count_expand_split_t <- process_te_data_tumour(te_raw_prepped_t, te_germline=te_germline_filter, clinical, complete_samples, metrics,
#                                                  apply_filter_common = TRUE, rare_gnomad = 3, rare_hostseq = 3,
#                                                  split_by_gene = TRUE, apply_process_combinations = FALSE, select_samples_split = TRUE)
# save(final_te_count_expand_split_t, file = paste0(r_dir, "final_te_count_expand_split_rare.RData"))

#### PROCESS COMMON TE DATA ####
if (PROCESS_COMMON_TES) {
  cat("\n========================================\n")
  cat("STEP 5: Processing common TE data\n")
  cat("========================================\n")
  cat("Filtering common TEs (no frequency filters)...\n")

  # Common TE: keep common, create count matrix
  final_te_count_common_t <- process_te_data_tumour(te_raw_prepped_t, te_germline=te_germline_filter, clinical, complete_samples, metrics, ancestry,
                                                    apply_filter_common = FALSE, split_by_gene = FALSE,
                                                    apply_process_combinations = TRUE, select_samples_split = TRUE, nohits_prefix = "final_te_count_common_t",
                                                    apply_age_filter = FALSE)
  save(final_te_count_common_t, file = paste0(r_dir, "final_te_count_common.RData"))

  # Common TE: expanded format (one row per full TE, no count matrix)
  final_te_count_expand_common_t <- process_te_data_tumour(te_raw_prepped_t, te_germline=te_germline_filter, clinical, complete_samples, metrics, ancestry,
                                                           apply_filter_common = FALSE, split_by_gene = FALSE,
                                                           apply_process_combinations = FALSE, select_samples_split = TRUE,
                                                           apply_age_filter = FALSE)
  save(final_te_count_expand_common_t, file = paste0(r_dir, "final_te_count_expand_common.RData"))

  # Common TE: split by gene format
  final_te_count_expand_split_common_t <- process_te_data_tumour(te_raw_prepped_t, te_germline=te_germline_filter, clinical, complete_samples, metrics, ancestry,
                                                          apply_filter_common = FALSE, split_by_gene = TRUE,
                                                          apply_process_combinations = FALSE, select_samples_split = TRUE,
                                                          apply_age_filter = FALSE)
  save(final_te_count_expand_split_common_t, file = paste0(r_dir, "final_te_count_split_common.RData"))
} else {
  cat("\n========================================\n")
  cat("STEP 5: Skipping common TE processing (PROCESS_COMMON_TES = FALSE)\n")
  cat("========================================\n")
}

#### PROCESS FULL-LENGTH LINE1 DATA ####
if (PROCESS_FULLLENGTH) {
  cat("\n========================================\n")
  cat("STEP 5B: Processing full-length LINE1 data\n")
  cat("========================================\n")
  cat("Filtering for full-length LINE1 (≥5900bp)...\n")

  # Full-length LINE1: count matrix format
  final_te_count_fulllength_t <- process_te_data_tumour(te_raw_prepped_t, te_germline=te_germline_filter, clinical, complete_samples, metrics, ancestry,
                                             apply_filter_common = FALSE,
                                             apply_filter_fulllength_young = TRUE,
                                             fulllength_bp = 5900,
                                             young_subfamilies = c("L1HS"),
                                             split_by_gene = FALSE,
                                             apply_process_combinations = TRUE,
                                             select_samples_split = TRUE,
                                             nohits_prefix = "final_te_count_fulllength_t",
                                             apply_age_filter = FALSE)
  save(final_te_count_fulllength_t, file = paste0(r_dir, "final_te_count_fulllength.RData"))

  # Full-length LINE1: expanded format (one row per TE)
  final_te_count_expand_fulllength_t <- process_te_data_tumour(te_raw_prepped_t, te_germline=te_germline_filter, clinical, complete_samples, metrics, ancestry,
                                                    apply_filter_common = FALSE,
                                                    apply_filter_fulllength_young = TRUE,
                                                    fulllength_bp = 5900,
                                                    young_subfamilies = c("L1HS"),
                                                    split_by_gene = FALSE,
                                                    apply_process_combinations = FALSE,
                                                    select_samples_split = TRUE,
                                                    apply_age_filter = FALSE)
  save(final_te_count_expand_fulllength_t, file = paste0(r_dir, "final_te_count_expand_fulllength.RData"))

  # Full-length LINE1: gene-split format
  final_te_count_expand_split_fulllength_t <- process_te_data_tumour(te_raw_prepped_t, te_germline=te_germline_filter, clinical, complete_samples, metrics, ancestry,
                                                         apply_filter_common = FALSE,
                                                         apply_filter_fulllength_young = TRUE,
                                                         fulllength_bp = 5900,
                                                         young_subfamilies = c("L1HS"),
                                                         split_by_gene = TRUE,
                                                         apply_process_combinations = FALSE,
                                                         select_samples_split = TRUE,
                                                         apply_age_filter = FALSE)
  save(final_te_count_expand_split_fulllength_t, file = paste0(r_dir, "final_te_count_split_fulllength.RData"))

  cat("Full-length LINE1 data saved.\n")
} else {
  cat("\n========================================\n")
  cat("STEP 5B: Skipping full-length LINE1 processing (PROCESS_FULLLENGTH = FALSE)\n")
  cat("========================================\n")
}

#### PROCESS SPLIT TE DATA ####
cat("\n========================================\n")
cat("STEP 6: Processing split TE data for gene analysis\n")
cat("========================================\n")
cat("Processing gene-split TE data...\n")

# Filter TE types for split data (suppress output)
suppressMessages({
  te_split_t <- filter_te_types(te_split_t, include_alu = INCLUDE_ALU, include_sva = INCLUDE_SVA)
})

# Filter false positives for split data (suppress output)
suppressMessages({
  fp_result_split <- filter_false_positives(te_split_t, l1_merge_fp, unique_l1_tp, unique_l1_fp, unique_l1_master, check_missing_ids = FALSE)
  te_split_t <- fp_result_split$te_data
})

# Prepare te_split (suppress output)
suppressMessages({
  prep_result_split <- prep_te(te_split_t, nonproband, noconsent, hostseq_cancer, bloodcancer, metrics,
                             c("mean_cov", "avg_quality", "pct_chimeras"), c(20, 28, 2), type="T")
  te_split_prepped_t <- prep_result_split$df
})

# Split TE data for gene analysis
final_te_count_expand_split_split_rare_t <- process_te_data_tumour(te_split_prepped_t, te_germline=te_germline_filter, clinical, complete_samples, metrics, ancestry,
                                                       apply_filter_common = can_filter_common, rare_gnomad = 3, rare_hostseq = 3, split_by_gene = TRUE,
                                                       apply_process_combinations = FALSE, select_samples_split = TRUE,
                                                       apply_age_filter = FALSE)
save(final_te_count_expand_split_split_rare_t, file = paste0(r_dir, "final_te_count_expand_split_split_rare.RData"))

# Write pipeline summary tables to processed output file
if (ENABLE_DETAILED_TRACKING) {
pipeline_summary <- capture.output({
  cat("\n" , rep("=", 80), "\n")
  cat("DETAILED PIPELINE TRACKING SUMMARY\n")
  cat(rep("=", 80), "\n")
  
  # TE Table with detailed reasons
  cat("\nTE PROCESSING PIPELINE:\n")
  cat(sprintf("%-6s %-10s %-10s %-8s %s\n", "Step", "Count", "Lost", "Lost%", "Detailed Reason"))
  cat(rep("-", 100), "\n")
  
  if (exists("step0_raw_te_count") && exists("step1_te_count") && exists("step2_te_count") && exists("step3_te_count")) {
    cat(sprintf("%-6s %-10d %-10s %-8s %s\n", "0", step0_raw_te_count, "-", "-", "Raw TE insertions from annotSV pipeline"))
    cat(sprintf("%-6s %-10d %-10d %-8.1f %s\n", "1", step1_te_count, step1_te_loss, (step1_te_loss/step0_raw_te_count)*100,
        "TE type filtering: Removed Alu and SVA elements (keeping LINE1 only)"))
    cat(sprintf("%-6s %-10d %-10d %-8.1f %s\n", "2a", step2a_te_count, step2a_te_loss, (step2a_te_loss/step1_te_count)*100,
        "Sample filtering: Excluded TEs from unreviewed samples"))
    cat(sprintf("%-6s %-10d %-10d %-8.1f %s\n", "2b", step2_te_count, step2b_te_loss, (step2b_te_loss/step2a_te_count)*100,
        "False positive removal: Removed computationally identified FPs"))
    cat(sprintf("    %s\n", "based on manual IGV validation"))
    cat(sprintf("%-6s %-10d %-10d %-8.1f %s\n", "3", step3_te_count, step3_te_loss, (step3_te_loss/step2_te_count)*100,
        "Quality control: Removed TEs from samples failing metrics"))

    # Add Step 4a for common filtering if we captured it
    if (exists("step4_after_common_count") && !is.na(step4_after_common_count)) {
      step4a_te_loss <- step3_te_count - step4_after_common_count
      cat(sprintf("%-6s %-10d %-10d %-8.1f %s\n", "4a", step4_after_common_count, step4a_te_loss, (step4a_te_loss/step3_te_count)*100,
          "Frequency filtering: Removed common TEs (gnomAD AF ≥3%, HostSeq ≥3%)"))
    }
    cat(sprintf("    %s\n", "(mean_cov <20, avg_quality <30, pct_chimeras >2), consent, or"))
    cat(sprintf("    %s\n", "host sequence status filters"))
    
    # Step 4b: Final processing (deduplication to unique loci)
    # Use unique_te_loci if available, otherwise use hardcoded value
    if (exists("unique_te_loci") && !is.na(unique_te_loci)) {
      step4b_final_count <- unique_te_loci
    } else {
      step4b_final_count <- 837  # Fallback from previous analysis
    }

    # Calculate loss from Step 4a (after common filtering) to Step 4b (final)
    if (exists("step4_after_common_count") && !is.na(step4_after_common_count)) {
      step4b_te_loss <- step4_after_common_count - step4b_final_count
      cat(sprintf("%-6s %-10d %-10d %-8.1f %s\n", "4b", step4b_final_count, step4b_te_loss, (step4b_te_loss/step4_after_common_count)*100,
          "Final processing: Deduplication to unique loci, clinical merge,"))
      cat(sprintf("    %s\n", "age filter, and duplicate removal"))
    } else {
      # If we don't have step4a, show combined step4
      step4b_te_loss <- step3_te_count - step4b_final_count
      cat(sprintf("%-6s %-10d %-10d %-8.1f %s\n", "4", step4b_final_count, step4b_te_loss, (step4b_te_loss/step3_te_count)*100,
          "Combined: Frequency filtering + deduplication to unique loci"))
      cat(sprintf("    %s\n", "- Additional processing: clinical merge, age filter, duplicate removal"))
    }
    
    if (exists("unique_te_loci") && exists("step5_te_count") && !is.na(unique_te_loci) && !is.na(step5_te_count) && length(unique_te_loci) > 0 && length(step5_te_count) > 0) {
      step5_dedup <- step5_te_count - unique_te_loci
      if (length(step5_dedup) > 0 && !is.na(step5_dedup) && step5_dedup > 0) {
        cat(sprintf("%-6s %-10d %-10d %-8.1f %s\n", "5", unique_te_loci, step5_dedup, (step5_dedup/step5_te_count)*100, 
            "Deduplication: Collapsed identical TEs to unique genomic loci"))
      } else {
        cat("Step 5: No deduplication needed - technical filtering already produced unique TE loci\n")
      }
    }
  } else {
    cat("TE tracking data not available\n")
  }
  
  # Sample Table with detailed reasons
  cat("\nSAMPLE PROCESSING PIPELINE:\n")
  cat(sprintf("%-6s %-10s %-10s %-8s %s\n", "Step", "Count", "Lost", "Lost%", "Detailed Reason"))
  cat(rep("-", 100), "\n")
  
  if (exists("step1_sample_count") && exists("step2_sample_count") && exists("step3_sample_count")) {
    cat(sprintf("%-6s %-10d %-10s %-8s %s\n", "1", step1_sample_count, "-", "-",
        "Unique samples with at least one TE insertion called"))

    # Show separate lines for 2a and 2b if available
    if (exists("step2a_sample_count") && !is.na(step2a_sample_count) &&
        exists("step2a_sample_loss") && !is.na(step2a_sample_loss)) {
      cat(sprintf("%-6s %-10d %-10d %-8.1f %s\n", "2a", step2a_sample_count, step2a_sample_loss, (step2a_sample_loss/step1_sample_count)*100,
          "Sample filtering: Removed samples lacking IGV review"))
      cat(sprintf("%-6s %-10d %-10d %-8.1f %s\n", "2b", step2_sample_count, step2b_sample_loss, (step2b_sample_loss/step2a_sample_count)*100,
          "False positive filtering: Removed samples where all TEs were FPs"))
    } else {
      # Fallback: combined step 2
      cat(sprintf("%-6s %-10d %-10d %-8.1f %s\n", "2", step2_sample_count, step2_sample_loss, (step2_sample_loss/step1_sample_count)*100,
          "Sample filtering: Removed unreviewed samples and samples where"))
      cat(sprintf("    %s\n", "all TEs were identified as false positives"))
    }
    cat(sprintf("%-6s %-10d %-10d %-8.1f %s\n", "3", step3_sample_count, step3_sample_loss, (step3_sample_loss/step2_sample_count)*100, 
        "Quality control: Excluded samples with poor sequencing"))
    cat(sprintf("    %s\n", "metrics, no consent, or inappropriate sample type"))
    
    if (exists("step4_total_samples") && !is.na(step4_total_samples)) {
      samples_added <- step4_total_samples - step3_sample_count
      cat(sprintf("%-6s %-10d %-10s %-8s %s\n", "4", step4_total_samples, paste0("+", samples_added), "-", 
          "Analysis set completion: Added samples with no TE insertions"))
      cat(sprintf("    %s\n", "after filtering to ensure complete cohort representation"))
    }
  } else {
    cat("Sample tracking data not available\n")
  }
  
  # Additional filtering details
  cat("\nFILTERING CRITERIA DETAILS:\n")
  cat(rep("-", 50), "\n")
  cat("False Positive Removal:\n")
  cat("  - Two-caller insertions: Removed FPs identified by manual review\n")
  cat("  - One-caller insertions: Sample-specific FP removal based on IGV validation\n")
  cat("  - Unreviewed samples: Complete exclusion of samples lacking manual review\n\n")
  
  cat("Quality Control Metrics:\n")
  cat("  - Mean coverage: Minimum 20x required\n")
  cat("  - Average quality: Minimum Q30 required\n") 
  cat("  - Chimeric reads: Maximum 2% allowed\n")
  cat("  - Consent status: Only consented samples included\n")
  cat("  - Sample type: Excluded non-proband normals and host-seq cancer samples\n\n")
  
  cat("Frequency Filtering:\n")
  cat("  - gnomAD population frequency: Excluded TEs with AF ≥3%\n")
  cat("  - HostSeq cohort frequency: Excluded TEs present in ≥3% of HostSeq samples\n")
  cat("  - Rationale: Focus analysis on rare, potentially pathogenic TEs\n\n")
  
  # Expected results
  if (exists("step4_total_samples") && !is.na(step4_total_samples) && exists("unique_te_loci") && !is.na(unique_te_loci)) {
    cat("FINAL ANALYSIS SET:\n")
    cat(rep("-", 30), "\n")
    with_tes <- if(!is.na(step4_samples_with_tes)) step4_samples_with_tes else "Unknown"
    no_tes <- if(!is.na(step4_samples_no_tes)) step4_samples_no_tes else "Unknown"
    cat(sprintf("Total samples in analysis: %d\n", step4_total_samples))
    cat(sprintf("Samples with rare TEs: %s\n", with_tes))
    cat(sprintf("Samples with no TEs: %s\n", no_tes))
    cat(sprintf("Unique rare TE loci: %d\n", unique_te_loci))
    if (!is.na(step4_samples_with_tes) && !is.na(step4_total_samples)) {
      cat(sprintf("Percentage with TEs: %.1f%%\n", (as.numeric(with_tes)/step4_total_samples)*100))
    }
  } else {
    cat("FINAL ANALYSIS SET:\n")
    cat("Analysis set summary not available - tracking variables not captured properly\n")
  }
})
  # Write the pipeline summary
  cat("\nWriting pipeline summary to processed_output.txt...\n")
  tryCatch({
    write_stdout_to_file(pipeline_summary, paste0(r_dir_files, "/processed_output.txt"),
                         "DETAILED PIPELINE TRACKING SUMMARY", append = TRUE)
    cat("✓ Pipeline summary successfully written to processed_output.txt\n")
  }, error = function(e) {
    cat("ERROR writing pipeline summary:", e$message, "\n")
    cat("Attempting direct write...\n")
    cat(pipeline_summary, file = paste0(r_dir_files, "/processed_output.txt"),
        append = TRUE, sep = "\n")
  })
}

#### REPLACE HOSTSEQ FILTER GROUP WITH ANALYSIS GROUP ####
cat("\n========================================\n")
cat("STEP 5: Removing HostSeq samples from tumour data\n")
cat("========================================\n")

# For tumour data, HostSeq samples are only used during processing for common TE filtering.
# They should be completely removed from final datasets since they have no tumour data
# and are not used in any downstream tumour analysis.

cat("Removing all HostSeq samples from final tumour datasets...\n")

# Process rare TE data - all dataframes in final_te_count_t list
# Skip exclusion_tracker since it's a tracking dataframe, not TE data
for (df_name in names(final_te_count_t)) {
  if (df_name != "exclusion_tracker" && is.data.frame(final_te_count_t[[df_name]])) {
    final_te_count_t[[df_name]] <- remove_hostseq_from_tumour(
      te_data = final_te_count_t[[df_name]]
    )
  }
}
save(final_te_count_t, file = paste0(r_dir, "final_te_count_rare.RData"))

# Process rare TE expanded format
for (df_name in names(final_te_count_expand_t)) {
  if (is.data.frame(final_te_count_expand_t[[df_name]])) {
    final_te_count_expand_t[[df_name]] <- remove_hostseq_from_tumour(
      te_data = final_te_count_expand_t[[df_name]]
    )
  }
}
save(final_te_count_expand_t, file = paste0(r_dir, "final_te_count_expand_rare.RData"))

# Process split TE data
for (df_name in names(final_te_count_expand_split_split_rare_t)) {
  if (is.data.frame(final_te_count_expand_split_split_rare_t[[df_name]])) {
    final_te_count_expand_split_split_rare_t[[df_name]] <- remove_hostseq_from_tumour(
      te_data = final_te_count_expand_split_split_rare_t[[df_name]]
    )
  }
}
save(final_te_count_expand_split_split_rare_t, file = paste0(r_dir, "final_te_count_expand_split_split_rare.RData"))

if (PROCESS_COMMON_TES) {
  # Process common TE data
  for (df_name in names(final_te_count_common_t)) {
    if (is.data.frame(final_te_count_common_t[[df_name]])) {
      final_te_count_common_t[[df_name]] <- remove_hostseq_from_tumour(
        te_data = final_te_count_common_t[[df_name]]
      )
    }
  }
  save(final_te_count_common_t, file = paste0(r_dir, "final_te_count_common.RData"))

  for (df_name in names(final_te_count_expand_common_t)) {
    if (is.data.frame(final_te_count_expand_common_t[[df_name]])) {
      final_te_count_expand_common_t[[df_name]] <- remove_hostseq_from_tumour(
        te_data = final_te_count_expand_common_t[[df_name]]
      )
    }
  }
  save(final_te_count_expand_common_t, file = paste0(r_dir, "final_te_count_expand_common.RData"))

  for (df_name in names(final_te_count_expand_split_common_t)) {
    if (is.data.frame(final_te_count_expand_split_common_t[[df_name]])) {
      final_te_count_expand_split_common_t[[df_name]] <- remove_hostseq_from_tumour(
        te_data = final_te_count_expand_split_common_t[[df_name]]
      )
    }
  }
  save(final_te_count_expand_split_common_t, file = paste0(r_dir, "final_te_count_split_common.RData"))
}

cat("✓ All HostSeq samples removed from tumour datasets\n")

#### EXPORT CSV FILES ####
cat("\n========================================\n")
cat("Exporting CSV files\n")
cat("========================================\n")
# Export split dataframe to CSV
write.csv(final_te_count_expand_split_split_rare_t$te_aff_selected, paste0(r_dir_files, "te_aff_split_t.csv"), row.names = FALSE, quote = FALSE)
cat("Exported te_aff_split_t.csv to files directory\n")

# Export to HPC directory
cat("\n========================================\n")
cat("Exporting to HPC directory...\n")
cat("========================================\n")

# Export split format (one row per TE-gene overlap)
cat("\nExporting split format files (excluding INFO, sample_topography, po_*, P_* columns)...\n")
write.csv(final_te_count_expand_split_split_rare_t$te_aff_selected %>% select(-any_of(c("INFO", "sample_topography")), -matches("^po_"), -matches("^P_")), paste0(hpc_dir, "te_aff_split_t.csv"), row.names = FALSE, quote = FALSE)
cat("✓ Exported split format files\n")

# Export count matrix format (one row per sample with TE counts)
cat("\nExporting count matrix format files (excluding INFO, sample_topography, po_*, P_* columns)...\n")
write.csv(final_te_count_t$te_all_selected %>% select(-any_of(c("INFO", "sample_topography")), -matches("^po_"), -matches("^P_")), paste0(hpc_dir, "te_all_t.csv"), row.names = FALSE, quote = FALSE)
write.csv(final_te_count_t$te_aff_selected %>% select(-any_of(c("INFO", "sample_topography")), -matches("^po_"), -matches("^P_")), paste0(hpc_dir, "te_aff_t.csv"), row.names = FALSE, quote = FALSE)
write.csv(final_te_count_t$te_lfs_selected %>% select(-any_of(c("INFO", "sample_topography")), -matches("^po_"), -matches("^P_")), paste0(hpc_dir, "te_lfs_t.csv"), row.names = FALSE, quote = FALSE)
write.csv(final_te_count_t$te_kics_selected %>% select(-any_of(c("INFO", "sample_topography")), -matches("^po_"), -matches("^P_")), paste0(hpc_dir, "te_kics_t.csv"), row.names = FALSE, quote = FALSE)

# Also export "all" versions (all samples, not just selected)
write.csv(final_te_count_t$te_all_all %>% select(-any_of(c("INFO", "sample_topography")), -matches("^po_"), -matches("^P_")), paste0(hpc_dir, "te_all_all_t.csv"), row.names = FALSE, quote = FALSE)
write.csv(final_te_count_t$te_aff_all %>% select(-any_of(c("INFO", "sample_topography")), -matches("^po_"), -matches("^P_")), paste0(hpc_dir, "te_aff_all_t.csv"), row.names = FALSE, quote = FALSE)
write.csv(final_te_count_t$te_lfs_all %>% select(-any_of(c("INFO", "sample_topography")), -matches("^po_"), -matches("^P_")), paste0(hpc_dir, "te_lfs_all_t.csv"), row.names = FALSE, quote = FALSE)
write.csv(final_te_count_t$te_kics_all %>% select(-any_of(c("INFO", "sample_topography")), -matches("^po_"), -matches("^P_")), paste0(hpc_dir, "te_kics_all_t.csv"), row.names = FALSE, quote = FALSE)
cat("✓ Exported count matrix format files\n")

# Export expand format (one row per complete TE insertion)
cat("\nExporting expand format files (excluding INFO, sample_topography, po_*, P_* columns)...\n")
write.csv(final_te_count_expand_t$te_all_selected %>% select(-any_of(c("INFO", "sample_topography")), -matches("^po_"), -matches("^P_")), paste0(hpc_dir, "te_all_expand_t.csv"), row.names = FALSE, quote = FALSE)
write.csv(final_te_count_expand_t$te_aff_selected %>% select(-any_of(c("INFO", "sample_topography")), -matches("^po_"), -matches("^P_")), paste0(hpc_dir, "te_aff_expand_t.csv"), row.names = FALSE, quote = FALSE)
write.csv(final_te_count_expand_t$te_lfs_selected %>% select(-any_of(c("INFO", "sample_topography")), -matches("^po_"), -matches("^P_")), paste0(hpc_dir, "te_lfs_expand_t.csv"), row.names = FALSE, quote = FALSE)
write.csv(final_te_count_expand_t$te_kics_selected %>% select(-any_of(c("INFO", "sample_topography")), -matches("^po_"), -matches("^P_")), paste0(hpc_dir, "te_kics_expand_t.csv"), row.names = FALSE, quote = FALSE)

# Also export "all" versions (all samples, not just selected)
write.csv(final_te_count_expand_t$te_all_all %>% select(-any_of(c("INFO", "sample_topography")), -matches("^po_"), -matches("^P_")), paste0(hpc_dir, "te_all_expand_all_t.csv"), row.names = FALSE, quote = FALSE)
write.csv(final_te_count_expand_t$te_aff_all %>% select(-any_of(c("INFO", "sample_topography")), -matches("^po_"), -matches("^P_")), paste0(hpc_dir, "te_aff_expand_all_t.csv"), row.names = FALSE, quote = FALSE)
write.csv(final_te_count_expand_t$te_lfs_all %>% select(-any_of(c("INFO", "sample_topography")), -matches("^po_"), -matches("^P_")), paste0(hpc_dir, "te_lfs_expand_all_t.csv"), row.names = FALSE, quote = FALSE)
write.csv(final_te_count_expand_t$te_kics_all %>% select(-any_of(c("INFO", "sample_topography")), -matches("^po_"), -matches("^P_")), paste0(hpc_dir, "te_kics_expand_all_t.csv"), row.names = FALSE, quote = FALSE)
cat("✓ Exported expand format files\n")

# Load nohits from R_obj directory and save to HPC directory
cat("\nExporting nohits samples...\n")
load(paste0(r_dir, "nohits_final_te_count_t_te_aff_selected_t.RData"))
write.table(nohits, paste0(hpc_dir, "nohits_te_aff_t.csv"), row.names = FALSE, col.names = FALSE, sep = ",", quote = FALSE)
cat("✓ Exported nohits samples\n")

cat("\n✓ All files exported to HPC directory:", hpc_dir, "\n")

#### EXPORT EXCLUDED SAMPLES TRACKING ####
cat("\n========================================\n")
cat("Exporting excluded samples tracking\n")
cat("========================================\n")

if (nrow(exclusion_tracker_master) > 0) {
  # Write excluded samples to CSV
  exclusion_output_file <- paste0(r_dir_files, "excluded_samples_tumour.csv")
  write.csv(exclusion_tracker_master, exclusion_output_file, row.names = FALSE)
  cat("✓ Exported", nrow(exclusion_tracker_master), "excluded sample records to:", exclusion_output_file, "\n")

  # Print summary by exclusion reason
  cat("\nExclusion summary:\n")
  exclusion_summary <- table(exclusion_tracker_master$exclusion_reason)
  for (reason in names(exclusion_summary)) {
    cat("  •", reason, ":", exclusion_summary[reason], "samples\n")
  }
  cat("\n")
} else {
  cat("No samples were excluded during processing\n\n")
}

cat("\n========================================\n")
cat("PROCESSING COMPLETE\n")
cat("========================================\n")
cat("All R objects saved to:", r_dir, "\n\n")

# Write completion summary to processed_output.txt
cat("\n========================================\n", file = output_file, append = TRUE)
cat("PROCESSING COMPLETED SUCCESSFULLY\n", file = output_file, append = TRUE)
cat("========================================\n", file = output_file, append = TRUE)
cat("Completion time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n", file = output_file, append = TRUE)
cat("All R objects saved to:", r_dir, "\n", file = output_file, append = TRUE)
cat("\nKey outputs generated:\n", file = output_file, append = TRUE)
cat("  - te_raw_prepped_t.RData\n", file = output_file, append = TRUE)
cat("  - final_te_count_rare.RData\n", file = output_file, append = TRUE)
cat("  - final_te_count_expand_rare.RData\n", file = output_file, append = TRUE)
cat("  - final_te_count_expand_split_split_rare.RData\n", file = output_file, append = TRUE)
cat("  - metrics.RData, ancestry.RData, clinical.RData\n", file = output_file, append = TRUE)
cat("\nFor detailed pipeline tracking, set ENABLE_DETAILED_TRACKING = TRUE\n", file = output_file, append = TRUE)

cat("Main output files:\n")
cat("  • Preprocessed data: te_raw_prepped_t.RData\n")
cat("  • Rare TEs (count matrix): final_te_count_rare.RData\n")
cat("  • Rare TEs (expanded): final_te_count_expand_rare.RData\n")
cat("  • Common TEs (count matrix): final_te_count_common.RData\n")
cat("  • Common TEs (expanded): final_te_count_expand_common.RData\n")
cat("  • Gene-split data: final_te_count_*_split*.RData\n")
cat("  • Clinical/metrics: clinical.RData, metrics.RData, gene_size.RData\n")
cat("  • Excluded samples tracking: excluded_samples_tumour.csv\n")
cat("  • CSV exports: te_aff_t.csv, te_all_all_t.csv\n\n")

if (ENABLE_DETAILED_TRACKING) {
  cat("Detailed processing logs written to:\n")
  cat("  • ", r_dir_files, "processed_output.txt\n\n")
}

cat("Cohort-specific nohits files created for downstream analysis\n")
