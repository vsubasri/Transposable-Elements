#!/usr/bin/env Rscript

# Germline TE Visualization - Specific TEs
# Unique TEs in LFS, specific TEs by TP53/Cancer status

# Source common setup and load data
source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/00_viz_common_setup.R")
source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/00_viz_load_data_germline.R")

cat("Running 02_te_viz_germline_03_specific_tes.R...\n")

#### UNIQUE TE IN LFS ####
# Defensive check: only call if sig_te has rows and required columns exist in te_aff_expand
if (exists("sig_te") && nrow(sig_te) > 0 && all(c("ALT") %in% colnames(te_aff_expand))) {
  sig_te_samples <- find_sig_te_samples(sig_te, te_aff_expand)
} else {
  sig_te_samples <- NULL
  warning("sig_te is empty or required columns are missing in te_aff_expand; skipping find_sig_te_samples.")
}


#### SPECIFIC TEs BY TP53 STATUS ####
cat("\n===== SPECIFIC TEs BY TP53 STATUS =====\n")

# Affected only (te_aff_expand)
cat("\n--- Analyzing AFFECTED samples (te_aff_expand) ---\n")
tryCatch({
  # Test specific TE insertions for differential representation by TP53 status
  # Run with multiple thresholds
  for (min_samples in c(3, 5, 10)) {
    cat("\n--- Testing with min_samples =", min_samples, "---\n")
    te_specific_results_aff <- test_specific_tes_by_group(
      te_expand = te_aff_expand,
      te_count = te_aff,
      group_column = "TP53_status",
      min_samples_with = min_samples,
      min_samples_without = min_samples,
      output_dir = r_dir_files,
      output_prefix = paste0("specific_tes_germline_aff_min", min_samples)
    )
  }

  if (!is.null(te_specific_results_aff)) {
    write_output(quote(head(te_specific_results_aff$full_results, 20)), "Top 20 TEs by adjusted p-value (affected)")
  }

  # Test by gene (all genes)
  cat("\n--- Testing BY GENE (all genes) ---\n")
  for (min_samples in c(3, 5, 10)) {
    cat("\n--- Testing with min_samples =", min_samples, "---\n")
    te_specific_by_gene_aff <- test_specific_tes_by_gene(
      te_expand = te_aff_expand,
      te_count = te_aff,
      group_column = "TP53_status",
      min_samples_with = min_samples,
      min_samples_without = min_samples,
      gene_filter = NULL,
      output_dir = r_dir_files,
      output_prefix = paste0("specific_tes_germline_aff_by_gene_min", min_samples)
    )
  }

  # Test by gene (cancer genes only)
  cat("\n--- Testing BY GENE (cancer genes only) ---\n")
  for (min_samples in c(3, 5, 10)) {
    cat("\n--- Testing with min_samples =", min_samples, "---\n")
    te_specific_by_gene_cancer_aff <- test_specific_tes_by_gene(
      te_expand = te_aff_expand,
      te_count = te_aff,
      group_column = "TP53_status",
      min_samples_with = min_samples,
      min_samples_without = min_samples,
      gene_filter = genes,
      output_dir = r_dir_files,
      output_prefix = paste0("specific_tes_germline_aff_by_gene_cancer_min", min_samples)
    )
  }
}, error = function(e) {
  cat("Warning: Could not perform specific TE testing (affected):", e$message, "\n")
})

# Affected + Unaffected (te_aff_unaff_expand) - includes unaffected controls
cat("\n--- Analyzing AFFECTED + UNAFFECTED samples (te_aff_unaff_expand) ---\n")
tryCatch({
  # Create count matrix for te_aff_unaff if needed
  te_aff_unaff <- te_all %>% filter(cohort != "Taylor" & cohort != "HostSeq")

  # Test specific TE insertions for differential representation by TP53 status
  # Including unaffected controls provides more power for detection

  # Test by coordinates (original approach)
  for (min_samples in c(3, 5, 10)) {
    cat("\n--- Testing with min_samples =", min_samples, "---\n")
    te_specific_results_aff_unaff <- test_specific_tes_by_group(
      te_expand = te_aff_unaff_expand,
      te_count = te_aff_unaff,
      group_column = "TP53_status",
      min_samples_with = min_samples,
      min_samples_without = min_samples,
      output_dir = r_dir_files,
      output_prefix = paste0("specific_tes_germline_aff_unaff_min", min_samples)
    )
  }

  # Test by gene (grouped approach) - all genes
  for (min_samples in c(3, 5, 10)) {
    cat("\n--- Testing gene-grouped with min_samples =", min_samples, "---\n")
    te_specific_by_gene_aff_unaff <- test_specific_tes_by_gene(
      te_expand = te_aff_unaff_expand,
      te_count = te_aff_unaff,
      group_column = "TP53_status",
      min_samples_with = min_samples,
      min_samples_without = min_samples,
      gene_filter = NULL,
      output_dir = r_dir_files,
      output_prefix = paste0("specific_tes_germline_aff_unaff_by_gene_min", min_samples)
    )
  }

  # Test by gene (grouped approach) - cancer genes only
  for (min_samples in c(3, 5, 10)) {
    cat("\n--- Testing gene-grouped (cancer genes only) with min_samples =", min_samples, "---\n")
    te_specific_by_gene_cancer_aff_unaff <- test_specific_tes_by_gene(
      te_expand = te_aff_unaff_expand,
      te_count = te_aff_unaff,
      group_column = "TP53_status",
      min_samples_with = min_samples,
      min_samples_without = min_samples,
      gene_filter = genes,
      output_dir = r_dir_files,
      output_prefix = paste0("specific_tes_germline_aff_unaff_by_gene_cancer_min", min_samples)
    )
  }

  if (!is.null(te_specific_results_aff_unaff)) {
    write_output(quote(head(te_specific_results_aff_unaff$full_results, 20)), "Top 20 TEs by adjusted p-value (affected + unaffected)")
  }
}, error = function(e) {
  cat("Warning: Could not perform specific TE testing (aff + unaff):", e$message, "\n")
})


#### SPECIFIC TEs BY CANCER STATUS (LFS) ####
cat("\n===== SPECIFIC TEs BY CANCER STATUS (LFS) =====\n")
tryCatch({
  # Test for TEs specific to cancer status in LFS samples (affected vs unaffected)
  # Create Cancer_status column in te_aff_unaff_expand
  te_aff_unaff_expand_cancer <- te_aff_unaff_expand %>%
    mutate(Cancer_status = ifelse(tumor_type != "U", "Affected", "Unaffected"))

  te_aff_unaff_cancer <- te_all %>%
    filter(cohort != "Taylor" & cohort != "HostSeq") %>%
    mutate(Cancer_status = ifelse(tumor_type != "U", "Affected", "Unaffected"))

  # Test by coordinates (original approach)
  for (min_samples_cancer in c(3, 5, 10)) {
    cat("\n--- Testing with min_samples =", min_samples_cancer, "---\n")
    te_specific_results_cancer <- test_specific_tes_by_group(
      te_expand = te_aff_unaff_expand_cancer,
      te_count = te_aff_unaff_cancer,
      group_column = "Cancer_status",
      min_samples_with = min_samples_cancer,
      min_samples_without = min_samples_cancer,
      output_dir = r_dir_files,
      output_prefix = paste0("specific_tes_germline_cancer_min", min_samples_cancer)
    )
  }

  # Test by gene (grouped approach) - all genes
  for (min_samples_cancer in c(3, 5, 10)) {
    cat("\n--- Testing gene-grouped with min_samples =", min_samples_cancer, "---\n")
    te_specific_by_gene_cancer <- test_specific_tes_by_gene(
      te_expand = te_aff_unaff_expand_cancer,
      te_count = te_aff_unaff_cancer,
      group_column = "Cancer_status",
      min_samples_with = min_samples_cancer,
      min_samples_without = min_samples_cancer,
      gene_filter = NULL,
      output_dir = r_dir_files,
      output_prefix = paste0("specific_tes_germline_cancer_by_gene_min", min_samples_cancer)
    )
  }

  # Test by gene (grouped approach) - cancer genes only
  for (min_samples_cancer in c(3, 5, 10)) {
    cat("\n--- Testing gene-grouped (cancer genes only) with min_samples =", min_samples_cancer, "---\n")
    te_specific_by_gene_cancer_cancer <- test_specific_tes_by_gene(
      te_expand = te_aff_unaff_expand_cancer,
      te_count = te_aff_unaff_cancer,
      group_column = "Cancer_status",
      min_samples_with = min_samples_cancer,
      min_samples_without = min_samples_cancer,
      gene_filter = genes,
      output_dir = r_dir_files,
      output_prefix = paste0("specific_tes_germline_cancer_by_gene_cancer_min", min_samples_cancer)
    )
  }

  if (!is.null(te_specific_results_cancer)) {
    write_output(quote(head(te_specific_results_cancer$full_results, 20)), "Top 20 TEs by adjusted p-value (Cancer-specific)")
  }
}, error = function(e) {
  cat("Warning: Could not perform specific TE testing (Cancer):", e$message, "\n")
})



cat("✓ Script completed successfully\n")
