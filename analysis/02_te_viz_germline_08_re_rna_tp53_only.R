#!/usr/bin/env Rscript

# Germline TE Visualization - TP53-Stratified RE-RNA Analysis ONLY
# Runs just the TP53-stratified section from 02_te_viz_germline_08_re_rna.R
# Uses 3-level TP53 classification: Germline, Somatic, WT

# Source common setup and load data
source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/00_viz_common_setup.R")
REQUIRED_DATA <- c("count_matrix", "expand", "split", "rna", "clinical")
source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/00_viz_load_data_germline.R")

# Initialize module-specific text output
init_module_sink(paste0(plot_dir, "reg_element/"), "RE_RNA_TP53")

cat("Running 02_te_viz_germline_08_re_rna_tp53_only.R...\n")
cat("TP53-stratified RE-RNA analysis with 3-level classification\n\n")

#### CONFIGURATION ####

# Databases to run
DATABASES_TO_RUN <- DEFAULT_DATABASES  # All 5 databases

# GSEA parameter settings
Q_VALUES <- c(0.05, 0.1, 0.25)
MIN_SAMPLES_VALUES <- c(3, 5)

# Create output directories
re_rna_dir <- paste0(plot_dir, "reg_element/re_rna/")
re_rna_files_dir <- paste0(re_rna_dir, "files/")
dir.create(re_rna_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(re_rna_files_dir, showWarnings = FALSE, recursive = TRUE)

cat("Output directories:\n")
cat("  Plots:", re_rna_dir, "\n")
cat("  Files:", re_rna_files_dir, "\n\n")

#### LOAD RNA DATA ####
cat("===== LOADING RNA DATA =====\n")

# Load rna_filtered from RNA script (required)
rna_filtered_file <- paste0(plot_dir, "rna/files/rna_filtered_germline.rds")
if (!file.exists(rna_filtered_file)) {
  stop("Required file not found: ", rna_filtered_file, "\nRun 02_te_viz_germline_06_rna.R first.")
}

cat("Loading rna_filtered from RNA script...\n")
rna_filtered <- readRDS(rna_filtered_file)
cat("Successfully loaded rna_filtered:", nrow(rna_filtered), "genes,", ncol(rna_filtered) - 1, "samples\n\n")

#### LOAD RE DATA ####
cat("===== LOADING RE DATA =====\n")

re_report_path_germline <- "/Users/briannelaverty/Documents/R_Malkin/te/data/final/germline_annotSV_output.SV_RE_intersect.report"

# Preprocess RE data for affected cohort (TP53-stratified analysis)
cat("Preprocessing RE data for affected cohort (TP53 analysis)...\n")
re_aff_preprocessed <- preprocess_re_for_rna(
  re_report_path = re_report_path_germline,
  valid_samples = unique(te_aff$sample)
)
cat("Affected RE data:", nrow(re_aff_preprocessed$re_split), "rows,",
    length(unique(re_aff_preprocessed$re_split$sample_id)), "samples\n\n")

#### TP53-STRATIFIED RE-RNA ANALYSIS ####
cat("===== TP53-STRATIFIED RE-RNA ANALYSIS =====\n")

# Add 3-level TP53 classification (Germline/Somatic/WT)
te_aff <- add_tp53_3level(te_aff)
cat("TP53_3level distribution:\n")
print(table(te_aff$TP53_3level, useNA = "always"))
cat("\n")

# Get TP53 groups from te_aff (affected cohort only)
# Using 3-level classification: Germline, Somatic, WT
tp53_groups <- list(
  "TP53germline" = te_aff %>% filter(TP53_3level == "Germline") %>% pull(sample) %>% unique(),
  "TP53somatic" = te_aff %>% filter(TP53_3level == "Somatic") %>% pull(sample) %>% unique(),
  "TP53wt" = te_aff %>% filter(TP53_3level == "WT") %>% pull(sample) %>% unique()
)

cat("TP53 group sample counts:\n")
for (grp in names(tp53_groups)) {
  cat("  ", grp, ":", length(tp53_groups[[grp]]), "samples\n")
}
cat("\n")

# Store results for combined analysis
tp53_de_results <- list()
tp53_gsea_results <- list()

for (tp53_group in names(tp53_groups)) {
  samples_in_group <- tp53_groups[[tp53_group]]

  cat("\n--- Analyzing TP53 group:", tp53_group, "---\n")
  cat("Samples in group:", length(samples_in_group), "\n")

  if (length(samples_in_group) < 5) {
    cat("Skipping", tp53_group, "- too few samples\n")
    next
  }

  # Filter RE data (affected cohort) to these samples
  re_tp53_split <- re_aff_preprocessed$re_split %>%
    filter(sample_id %in% samples_in_group)

  if (nrow(re_tp53_split) == 0) {
    cat("No RE data for", tp53_group, ". Skipping.\n")
    next
  }

  cat("RE entries in group:", nrow(re_tp53_split), "\n")

  # Recalculate gene counts for this subset
  gene_counts_tp53 <- re_tp53_split %>%
    count(gene_reg, name = "n_occurrences") %>%
    filter(n_occurrences >= 3) %>%
    arrange(desc(n_occurrences))

  re_preprocessed_tp53 <- list(re_split = re_tp53_split, gene_counts = gene_counts_tp53)

  for (min_samples in MIN_SAMPLES_VALUES) {
    cat("\n  min_samples =", min_samples, "\n")

    # Run differential expression test
    re_rna_tp53_results <- test_rna_by_gene_re_status(
      re_report_path = re_report_path_germline,
      rna_data = rna_filtered,
      sample_type = "germline",
      min_samples_per_group = min_samples,
      valid_samples = samples_in_group,
      re_preprocessed = re_preprocessed_tp53
    )

    if (nrow(re_rna_tp53_results) > 0) {
      # Store for combined analysis
      tp53_de_results[[paste0(tp53_group, "_min", min_samples)]] <- re_rna_tp53_results

      # Save differential expression results
      output_file <- paste0(re_rna_files_dir, "re_rna_differential_", tp53_group, "_min", min_samples, ".csv")
      write.csv(re_rna_tp53_results, output_file, row.names = FALSE)
      cat("  Saved:", nrow(re_rna_tp53_results), "genes tested\n")

      # Calculate logFC for GSEA
      if (all(c("mean_with_re", "mean_without_re") %in% colnames(re_rna_tp53_results))) {
        re_rna_with_logfc <- re_rna_tp53_results %>%
          mutate(
            logFC = log2((mean_with_re + 1) / (mean_without_re + 1))
          ) %>%
          filter(!is.na(logFC) & is.finite(logFC))

        if (nrow(re_rna_with_logfc) >= 10) {
          gene_list <- re_rna_with_logfc$logFC
          names(gene_list) <- re_rna_with_logfc$gene
          gene_list <- sort(gene_list, decreasing = TRUE)

          cat("  Running GSEA with", length(gene_list), "genes\n")

          gsea_results <- run_gsea_parameter_sweep_multi_db(
            gene_logfc = gene_list,
            output_dir = re_rna_dir,
            output_prefix = paste0("re_rna_gsea_", tp53_group, "_min", min_samples),
            databases = DATABASES_TO_RUN,
            q_values = Q_VALUES
          )

          tp53_gsea_results[[paste0(tp53_group, "_min", min_samples)]] <- gsea_results
        }
      }
    } else {
      cat("  No differential expression results\n")
    }
  }
}

#### TP53 COMBINED ANALYSIS ####
cat("\n===== TP53 COMBINED ANALYSIS =====\n")

# Create combined TP53 x TE boxplots for each min_samples and q threshold
if (length(tp53_de_results) > 0) {
  cat("Creating combined TP53 x TE boxplots...\n")

  for (min_samples in MIN_SAMPLES_VALUES) {
    for (q_thresh in Q_VALUES) {
      combined_boxplot_path <- paste0(re_rna_dir, "re_rna_de_boxplot_TP53combined_min", min_samples, "_q", q_thresh, ".png")
      tryCatch({
        plot_de_genes_tp53_combined_boxplot(
          de_results_list = tp53_de_results,
          rna_data = rna_filtered,
          tp53_sample_groups = tp53_groups,
          output_path = combined_boxplot_path,
          p_cutoff = q_thresh,
          min_samples_key = paste0("min", min_samples),
          use_raw_p = TRUE,
          max_genes = 10
        )
      }, error = function(e) {
        cat("  Warning: Could not create combined boxplot for min", min_samples, "q", q_thresh, ":", e$message, "\n")
      })
    }
  }
}

# Create combined GSEA comparison plots across TP53 groups
if (length(tp53_gsea_results) > 0) {
  cat("\nCreating TP53 GSEA comparison plots...\n")

  for (min_samples in MIN_SAMPLES_VALUES) {
    for (q_thresh in Q_VALUES) {
      tryCatch({
        plot_gsea_tp53_comparison(
          gsea_results_list = tp53_gsea_results,
          output_dir = re_rna_dir,
          output_prefix = "re_rna_gsea_TP53comparison",
          databases = DATABASES_TO_RUN,
          q_thresh = q_thresh,
          min_samples_key = paste0("min", min_samples)
        )
      }, error = function(e) {
        cat("  Warning: Could not create GSEA comparison for min", min_samples, "q", q_thresh, ":", e$message, "\n")
      })
    }
  }
}

#### SUMMARY ####
cat("\n===== SUMMARY =====\n")

cat("TP53-stratified RE-RNA analysis completed:\n")
cat("  - 3-level TP53 classification: Germline, Somatic, WT\n")
cat("  - Differential expression for each TP53 group\n")
cat("  - GSEA pathway analysis for each TP53 group\n")
cat("  - Combined TP53 x TE status boxplots\n")
cat("  - TP53 GSEA comparison plots\n\n")

cat("Databases used:", paste(DATABASES_TO_RUN, collapse = ", "), "\n")
cat("Min samples:", paste(MIN_SAMPLES_VALUES, collapse = ", "), "\n")
cat("Q-value thresholds:", paste(Q_VALUES, collapse = ", "), "\n")

cat("\nOutput locations:\n")
cat("  Plots:", re_rna_dir, "\n")
cat("  CSV files:", re_rna_files_dir, "\n")

cat("\n Script completed successfully\n")

# Close module-specific sink
close_module_sink()
