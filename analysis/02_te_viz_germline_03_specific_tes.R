#!/usr/bin/env Rscript

# Germline TE Visualization - Specific TEs
# Unique TEs in LFS, specific TEs by TP53/Cancer status

# Source common setup and load data
source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/00_viz_common_setup.R")
REQUIRED_DATA <- c("count_matrix", "expand", "split", "clinical", "genes")
source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/00_viz_load_data_germline.R")

# Initialize module-specific text output
init_module_sink(paste0(plot_dir, "cancer_genes/"), "SPECIFIC_TES")

cat("Running 02_te_viz_germline_03_specific_tes.R...\n")

# Load cancer predisposition genes
genes <- cpg$V1

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
      output_prefix = paste0("specific_tes_aff_TP53_fullins_min", min_samples)
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
      te_expand = te_aff_split,
      te_count = te_aff,
      group_column = "TP53_status",
      min_samples_with = min_samples,
      min_samples_without = min_samples,
      gene_filter = NULL,
      output_dir = r_dir_files,
      output_prefix = paste0("specific_tes_aff_TP53_gene_min", min_samples)
    )
  }

  # Test by gene (cancer genes only)
  cat("\n--- Testing BY GENE (cancer genes only) ---\n")
  for (min_samples in c(3, 5, 10)) {
    cat("\n--- Testing with min_samples =", min_samples, "---\n")
    te_specific_by_gene_cancer_aff <- test_specific_tes_by_gene(
      te_expand = te_aff_split,
      te_count = te_aff,
      group_column = "TP53_status",
      min_samples_with = min_samples,
      min_samples_without = min_samples,
      gene_filter = genes,
      output_dir = r_dir_files,
      output_prefix = paste0("specific_tes_aff_TP53_cancergene_min", min_samples)
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
      output_prefix = paste0("specific_tes_affunaff_TP53_fullins_min", min_samples)
    )
  }

  # Test by gene (grouped approach) - all genes
  for (min_samples in c(3, 5, 10)) {
    cat("\n--- Testing gene-grouped with min_samples =", min_samples, "---\n")
    te_specific_by_gene_aff_unaff <- test_specific_tes_by_gene(
      te_expand = te_aff_unaff_split,
      te_count = te_aff_unaff,
      group_column = "TP53_status",
      min_samples_with = min_samples,
      min_samples_without = min_samples,
      gene_filter = NULL,
      output_dir = r_dir_files,
      output_prefix = paste0("specific_tes_affunaff_TP53_gene_min", min_samples)
    )
  }

  # Test by gene (grouped approach) - cancer genes only
  for (min_samples in c(3, 5, 10)) {
    cat("\n--- Testing gene-grouped (cancer genes only) with min_samples =", min_samples, "---\n")
    te_specific_by_gene_cancer_aff_unaff <- test_specific_tes_by_gene(
      te_expand = te_aff_unaff_split,
      te_count = te_aff_unaff,
      group_column = "TP53_status",
      min_samples_with = min_samples,
      min_samples_without = min_samples,
      gene_filter = genes,
      output_dir = r_dir_files,
      output_prefix = paste0("specific_tes_affunaff_TP53_cancergene_min", min_samples)
    )
  }

  if (!is.null(te_specific_results_aff_unaff)) {
    write_output(quote(head(te_specific_results_aff_unaff$full_results, 20)), "Top 20 TEs by adjusted p-value (affected + unaffected)")
  }
}, error = function(e) {
  cat("Warning: Could not perform specific TE testing (aff + unaff):", e$message, "\n")
})


#### SPECIFIC TEs BY CANCER STATUS (KICS + LFS) ####
cat("\n===== SPECIFIC TEs BY CANCER STATUS (KICS + LFS) =====\n")
tryCatch({
  # Test for TEs specific to cancer status in KICS+LFS samples (affected vs unaffected)
  # Create Cancer_status column in te_aff_unaff_expand and te_aff_unaff_split
  te_aff_unaff_expand_cancer <- te_aff_unaff_expand %>%
    mutate(Cancer_status = ifelse(tumor_type != "U", "Affected", "Unaffected"))

  te_aff_unaff_split_cancer <- te_aff_unaff_split %>%
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
      output_prefix = paste0("specific_tes_affunaff_Cancer_fullins_min", min_samples_cancer)
    )
  }

  # Test by gene (grouped approach) - all genes
  for (min_samples_cancer in c(3, 5, 10)) {
    cat("\n--- Testing gene-grouped with min_samples =", min_samples_cancer, "---\n")
    te_specific_by_gene_cancer <- test_specific_tes_by_gene(
      te_expand = te_aff_unaff_split_cancer,
      te_count = te_aff_unaff_cancer,
      group_column = "Cancer_status",
      min_samples_with = min_samples_cancer,
      min_samples_without = min_samples_cancer,
      gene_filter = NULL,
      output_dir = r_dir_files,
      output_prefix = paste0("specific_tes_affunaff_Cancer_gene_min", min_samples_cancer)
    )
  }

  # Test by gene (grouped approach) - cancer genes only
  for (min_samples_cancer in c(3, 5, 10)) {
    cat("\n--- Testing gene-grouped (cancer genes only) with min_samples =", min_samples_cancer, "---\n")
    te_specific_by_gene_cancer_cancer <- test_specific_tes_by_gene(
      te_expand = te_aff_unaff_split_cancer,
      te_count = te_aff_unaff_cancer,
      group_column = "Cancer_status",
      min_samples_with = min_samples_cancer,
      min_samples_without = min_samples_cancer,
      gene_filter = genes,
      output_dir = r_dir_files,
      output_prefix = paste0("specific_tes_affunaff_Cancer_cancergene_min", min_samples_cancer)
    )
  }

  if (!is.null(te_specific_results_cancer)) {
    write_output(quote(head(te_specific_results_cancer$full_results, 20)), "Top 20 TEs by adjusted p-value (Cancer-specific)")
  }
}, error = function(e) {
  cat("Warning: Could not perform specific TE testing (Cancer - KICS+LFS):", e$message, "\n")
})


#### SPECIFIC TEs BY CANCER STATUS (LFS only) ####
cat("\n===== SPECIFIC TEs BY CANCER STATUS (LFS only) =====\n")
tryCatch({
  # Test for TEs specific to cancer status in LFS samples only (affected vs unaffected)
  # Create Cancer_status column in LFS datasets
  te_lfs_expand_cancer <- te_lfs_expand %>%
    mutate(Cancer_status = ifelse(tumor_type != "U", "Affected", "Unaffected"))

  te_lfs_split_cancer <- te_lfs_split %>%
    mutate(Cancer_status = ifelse(tumor_type != "U", "Affected", "Unaffected"))

  te_lfs_cancer <- te_lfs %>%
    mutate(Cancer_status = ifelse(tumor_type != "U", "Affected", "Unaffected"))

  # Test by coordinates (original approach)
  for (min_samples_cancer in c(3, 5, 10)) {
    cat("\n--- Testing with min_samples =", min_samples_cancer, "---\n")
    te_specific_results_lfs_cancer <- test_specific_tes_by_group(
      te_expand = te_lfs_expand_cancer,
      te_count = te_lfs_cancer,
      group_column = "Cancer_status",
      min_samples_with = min_samples_cancer,
      min_samples_without = min_samples_cancer,
      output_dir = r_dir_files,
      output_prefix = paste0("specific_tes_lfs_Cancer_fullins_min", min_samples_cancer)
    )
  }

  # Test by gene (grouped approach) - all genes
  for (min_samples_cancer in c(3, 5, 10)) {
    cat("\n--- Testing gene-grouped with min_samples =", min_samples_cancer, "---\n")
    te_specific_by_gene_lfs_cancer <- test_specific_tes_by_gene(
      te_expand = te_lfs_split_cancer,
      te_count = te_lfs_cancer,
      group_column = "Cancer_status",
      min_samples_with = min_samples_cancer,
      min_samples_without = min_samples_cancer,
      gene_filter = NULL,
      output_dir = r_dir_files,
      output_prefix = paste0("specific_tes_lfs_Cancer_gene_min", min_samples_cancer)
    )
  }

  # Test by gene (grouped approach) - cancer genes only
  for (min_samples_cancer in c(3, 5, 10)) {
    cat("\n--- Testing gene-grouped (cancer genes only) with min_samples =", min_samples_cancer, "---\n")
    te_specific_by_gene_cancer_lfs_cancer <- test_specific_tes_by_gene(
      te_expand = te_lfs_split_cancer,
      te_count = te_lfs_cancer,
      group_column = "Cancer_status",
      min_samples_with = min_samples_cancer,
      min_samples_without = min_samples_cancer,
      gene_filter = genes,
      output_dir = r_dir_files,
      output_prefix = paste0("specific_tes_lfs_Cancer_cancergene_min", min_samples_cancer)
    )
  }

  if (!is.null(te_specific_results_lfs_cancer)) {
    write_output(quote(head(te_specific_results_lfs_cancer$full_results, 20)), "Top 20 TEs by adjusted p-value (Cancer-specific - LFS only)")
  }
}, error = function(e) {
  cat("Warning: Could not perform specific TE testing (Cancer - LFS only):", e$message, "\n")
})


#### SPECIFIC TEs BY COHORT (KICS vs HostSeq) ####
cat("\n===== SPECIFIC TEs BY COHORT (KICS vs HostSeq) =====\n")
tryCatch({
  # Compare KICS (pediatric cancer patients) vs HostSeq (healthy controls)
  # Create count matrix for kics + hostseq by combining the two cohort matrices
  # (te_all doesn't contain HostSeq samples, so we need to rbind them manually)
  te_kics_hostseq_count <- rbind(te_kics, te_hostseq)

  # Test by coordinates (individual TE insertions)
  for (min_samples in c(3, 5, 10)) {
    cat("\n--- Testing with min_samples =", min_samples, "---\n")
    te_specific_results_kics_hostseq <- test_specific_tes_by_group(
      te_expand = te_kics_hostseq_expand,
      te_count = te_kics_hostseq_count,
      group_column = "cohort",
      min_samples_with = min_samples,
      min_samples_without = min_samples,
      output_dir = r_dir_files,
      output_prefix = paste0("specific_tes_kics_hostseq_fullins_min", min_samples)
    )
  }

  # Test by gene (all genes)
  for (min_samples in c(3, 5, 10)) {
    cat("\n--- Testing gene-grouped with min_samples =", min_samples, "---\n")
    te_specific_by_gene_kics_hostseq <- test_specific_tes_by_gene(
      te_expand = te_kics_hostseq,
      te_count = te_kics_hostseq_count,
      group_column = "cohort",
      min_samples_with = min_samples,
      min_samples_without = min_samples,
      gene_filter = NULL,
      output_dir = r_dir_files,
      output_prefix = paste0("specific_tes_kics_hostseq_gene_min", min_samples)
    )
  }

  # Test by gene (cancer genes only)
  for (min_samples in c(3, 5, 10)) {
    cat("\n--- Testing gene-grouped (cancer genes only) with min_samples =", min_samples, "---\n")
    te_specific_by_gene_cancer_kics_hostseq <- test_specific_tes_by_gene(
      te_expand = te_kics_hostseq,
      te_count = te_kics_hostseq_count,
      group_column = "cohort",
      min_samples_with = min_samples,
      min_samples_without = min_samples,
      gene_filter = genes,
      output_dir = r_dir_files,
      output_prefix = paste0("specific_tes_kics_hostseq_cancergene_min", min_samples)
    )
  }

  if (!is.null(te_specific_results_kics_hostseq)) {
    write_output(quote(head(te_specific_results_kics_hostseq$full_results, 20)), "Top 20 TEs by adjusted p-value (KICS vs HostSeq)")
  }
}, error = function(e) {
  cat("Warning: Could not perform specific TE testing (KICS vs HostSeq):", e$message, "\n")
})


#### SPECIFIC TEs BY TUMOR TYPE (One-vs-Rest) ####
cat("\n===== SPECIFIC TEs BY TUMOR TYPE (One-vs-Rest) =====\n")
cat("Note: Multiway Fisher test (global null) available as fisher_test_by_tumor_type() in functions_te.R\n")
tryCatch({
  # Test by insertion (fullins)
  cat("\n--- Testing BY INSERTION (fullins) ---\n")
  for (min_samples_te in c(3, 5)) {
    cat("\n--- Testing with min_samples_te =", min_samples_te, "---\n")
    te_specific_results_tt <- test_specific_tes_by_tumor_type(
      te_expand = te_aff_expand,
      te_count = te_aff,
      min_samples_tt = 10,
      min_samples_te = min_samples_te,
      output_dir = r_dir_files,
      output_prefix = paste0("specific_tes_aff_tumourtype_fullins_min", min_samples_te)
    )
  }

  # Test by gene (all genes)
  cat("\n--- Testing BY GENE (all genes) ---\n")
  for (min_samples_te in c(3, 5)) {
    cat("\n--- Testing with min_samples_te =", min_samples_te, "---\n")
    te_specific_by_gene_tt <- test_specific_tes_by_tumor_type_gene(
      te_expand = te_aff_split,
      te_count = te_aff,
      min_samples_tt = 10,
      min_samples_te = min_samples_te,
      gene_filter = NULL,
      output_dir = r_dir_files,
      output_prefix = paste0("specific_tes_aff_tumourtype_gene_min", min_samples_te)
    )
  }

  # Test by gene (cancer genes only)
  cat("\n--- Testing BY GENE (cancer genes only) ---\n")
  for (min_samples_te in c(3, 5)) {
    cat("\n--- Testing with min_samples_te =", min_samples_te, "---\n")
    te_specific_by_gene_cancer_tt <- test_specific_tes_by_tumor_type_gene(
      te_expand = te_aff_split,
      te_count = te_aff,
      min_samples_tt = 10,
      min_samples_te = min_samples_te,
      gene_filter = genes,
      output_dir = r_dir_files,
      output_prefix = paste0("specific_tes_aff_tumourtype_cancergene_min", min_samples_te)
    )
  }
}, error = function(e) {
  cat("Warning: Could not perform specific TE testing (tumor type):", e$message, "\n")
})


#### SPECIFIC TEs BY TAYLOR SUBTYPE ####
# Note: Taylor subtype analysis moved to 02_te_viz_germline_10_taylor.R
# (tumor_type_subclass has >2 groups which requires different handling)
cat("\n===== SPECIFIC TEs BY TAYLOR SUBTYPE =====\n")
cat("Skipping: Taylor subtype-specific TE analysis is performed in 02_te_viz_germline_10_taylor.R\n")


################################################################################
#### LOGISTIC REGRESSION ANALYSES (controlling for covariates) ####
################################################################################
cat("\n\n################################################################################\n")
cat("#### LOGISTIC REGRESSION ANALYSES (controlling for covariates) ####\n")
cat("################################################################################\n")
cat("Note: GLM analyses control for ancestry to account for population stratification\n\n")

# Define covariates for GLM analyses (same as covar_med but with mapped_label instead of predicted_ancestry_thres)
glm_covariates <- c("med_cov", "total_reads", "med_read_len", "pct_chimeras", "avg_quality", "tumor_type", "mapped_label")

#### GLM: SPECIFIC TEs BY TP53 STATUS ####
cat("\n===== GLM: SPECIFIC TEs BY TP53 STATUS =====\n")

# Affected only (te_aff_expand)
cat("\n--- GLM: Analyzing AFFECTED samples (te_aff_expand) ---\n")
tryCatch({
  # Test specific TE insertions controlling for ancestry
  for (min_samples in c(3, 5, 10)) {
    cat("\n--- GLM Testing with min_samples =", min_samples, "---\n")
    te_specific_results_aff_glm <- test_specific_tes_by_group_glm(
      te_expand = te_aff_expand,
      te_count = te_aff,
      group_column = "TP53_status",
      covariates = glm_covariates,
      min_samples_with = min_samples,
      min_samples_without = min_samples,
      output_dir = r_dir_files,
      output_prefix = paste0("specific_tes_aff_TP53_fullins_glm_min", min_samples)
    )
  }

  if (!is.null(te_specific_results_aff_glm)) {
    write_output(quote(head(te_specific_results_aff_glm$full_results, 20)), "Top 20 TEs by adjusted p-value (affected, GLM)")
  }

  # Test by gene (all genes) - GLM
  cat("\n--- GLM: Testing BY GENE (all genes) ---\n")
  for (min_samples in c(3, 5, 10)) {
    cat("\n--- GLM Testing with min_samples =", min_samples, "---\n")
    te_specific_by_gene_aff_glm <- test_specific_tes_by_gene_glm(
      te_expand = te_aff_split,
      te_count = te_aff,
      group_column = "TP53_status",
      covariates = glm_covariates,
      min_samples_with = min_samples,
      min_samples_without = min_samples,
      gene_filter = NULL,
      output_dir = r_dir_files,
      output_prefix = paste0("specific_tes_aff_TP53_gene_glm_min", min_samples)
    )
  }

  # Test by gene (cancer genes only) - GLM
  cat("\n--- GLM: Testing BY GENE (cancer genes only) ---\n")
  for (min_samples in c(3, 5, 10)) {
    cat("\n--- GLM Testing with min_samples =", min_samples, "---\n")
    te_specific_by_gene_cancer_aff_glm <- test_specific_tes_by_gene_glm(
      te_expand = te_aff_split,
      te_count = te_aff,
      group_column = "TP53_status",
      covariates = glm_covariates,
      min_samples_with = min_samples,
      min_samples_without = min_samples,
      gene_filter = genes,
      output_dir = r_dir_files,
      output_prefix = paste0("specific_tes_aff_TP53_cancergene_glm_min", min_samples)
    )
  }

  # Generate plots for significant results (q < 0.1) - use min_samples=5 results
  cat("\n--- Generating plots for significant results (q < 0.1) ---\n")

  # Plot for coordinate-level results
  if (!is.null(te_specific_results_aff_glm)) {
    p_fullins <- plot_specific_te_presence(te_specific_results_aff_glm, q_threshold = 0.1,
                                           id_column = "te_id",
                                           title = "TE Presence by TP53 Status (Affected, Coordinates)")
    if (!is.null(p_fullins)) {
      ggsave(paste0(plot_dir, "cancer_genes/specific_tes_aff_TP53_fullins_glm_presence.png"),
             plot = p_fullins, width = 10, height = 8)
    }
  }

  # Plot for gene-level results
  if (!is.null(te_specific_by_gene_aff_glm)) {
    p_gene <- plot_specific_te_presence(te_specific_by_gene_aff_glm, q_threshold = 0.1,
                                        id_column = "Gene_name",
                                        title = "TE Presence by TP53 Status (Affected, All Genes)")
    if (!is.null(p_gene)) {
      ggsave(paste0(plot_dir, "cancer_genes/specific_tes_aff_TP53_gene_glm_presence.png"),
             plot = p_gene, width = 10, height = 8)
    }
  }

  # Plot for cancer gene results
  if (!is.null(te_specific_by_gene_cancer_aff_glm)) {
    p_cancer <- plot_specific_te_presence(te_specific_by_gene_cancer_aff_glm, q_threshold = 0.1,
                                          id_column = "Gene_name",
                                          title = "TE Presence by TP53 Status (Affected, Cancer Genes)")
    if (!is.null(p_cancer)) {
      ggsave(paste0(plot_dir, "cancer_genes/specific_tes_aff_TP53_cancergene_glm_presence.png"),
             plot = p_cancer, width = 10, height = 8)
    }
  }
}, error = function(e) {
  cat("Warning: Could not perform GLM specific TE testing (affected):", e$message, "\n")
})

# Affected + Unaffected - GLM
cat("\n--- GLM: Analyzing AFFECTED + UNAFFECTED samples ---\n")
tryCatch({
  te_aff_unaff <- te_all %>% filter(cohort != "Taylor" & cohort != "HostSeq")

  # Test by coordinates - GLM
  for (min_samples in c(3, 5, 10)) {
    cat("\n--- GLM Testing with min_samples =", min_samples, "---\n")
    te_specific_results_aff_unaff_glm <- test_specific_tes_by_group_glm(
      te_expand = te_aff_unaff_expand,
      te_count = te_aff_unaff,
      group_column = "TP53_status",
      covariates = glm_covariates,
      min_samples_with = min_samples,
      min_samples_without = min_samples,
      output_dir = r_dir_files,
      output_prefix = paste0("specific_tes_affunaff_TP53_fullins_glm_min", min_samples)
    )
  }

  # Test by gene (all genes) - GLM
  for (min_samples in c(3, 5, 10)) {
    cat("\n--- GLM Testing gene-grouped with min_samples =", min_samples, "---\n")
    te_specific_by_gene_aff_unaff_glm <- test_specific_tes_by_gene_glm(
      te_expand = te_aff_unaff_split,
      te_count = te_aff_unaff,
      group_column = "TP53_status",
      covariates = glm_covariates,
      min_samples_with = min_samples,
      min_samples_without = min_samples,
      gene_filter = NULL,
      output_dir = r_dir_files,
      output_prefix = paste0("specific_tes_affunaff_TP53_gene_glm_min", min_samples)
    )
  }

  # Test by gene (cancer genes only) - GLM
  for (min_samples in c(3, 5, 10)) {
    cat("\n--- GLM Testing gene-grouped (cancer genes only) with min_samples =", min_samples, "---\n")
    te_specific_by_gene_cancer_aff_unaff_glm <- test_specific_tes_by_gene_glm(
      te_expand = te_aff_unaff_split,
      te_count = te_aff_unaff,
      group_column = "TP53_status",
      covariates = glm_covariates,
      min_samples_with = min_samples,
      min_samples_without = min_samples,
      gene_filter = genes,
      output_dir = r_dir_files,
      output_prefix = paste0("specific_tes_affunaff_TP53_cancergene_glm_min", min_samples)
    )
  }

  # Generate plots for significant results (q < 0.1)
  cat("\n--- Generating plots for significant results (q < 0.1) ---\n")
  if (!is.null(te_specific_results_aff_unaff_glm)) {
    p <- plot_specific_te_presence(te_specific_results_aff_unaff_glm, q_threshold = 0.1,
                                   id_column = "te_id",
                                   title = "TE Presence by TP53 Status (Aff+Unaff, Coordinates)")
    if (!is.null(p)) ggsave(paste0(plot_dir, "cancer_genes/specific_tes_affunaff_TP53_fullins_glm_presence.png"), plot = p, width = 10, height = 8)
  }
  if (!is.null(te_specific_by_gene_aff_unaff_glm)) {
    p <- plot_specific_te_presence(te_specific_by_gene_aff_unaff_glm, q_threshold = 0.1,
                                   id_column = "Gene_name",
                                   title = "TE Presence by TP53 Status (Aff+Unaff, All Genes)")
    if (!is.null(p)) ggsave(paste0(plot_dir, "cancer_genes/specific_tes_affunaff_TP53_gene_glm_presence.png"), plot = p, width = 10, height = 8)
  }
  if (!is.null(te_specific_by_gene_cancer_aff_unaff_glm)) {
    p <- plot_specific_te_presence(te_specific_by_gene_cancer_aff_unaff_glm, q_threshold = 0.1,
                                   id_column = "Gene_name",
                                   title = "TE Presence by TP53 Status (Aff+Unaff, Cancer Genes)")
    if (!is.null(p)) ggsave(paste0(plot_dir, "cancer_genes/specific_tes_affunaff_TP53_cancergene_glm_presence.png"), plot = p, width = 10, height = 8)
  }
}, error = function(e) {
  cat("Warning: Could not perform GLM specific TE testing (aff + unaff):", e$message, "\n")
})


#### GLM: SPECIFIC TEs BY CANCER STATUS (KICS + LFS) ####
cat("\n===== GLM: SPECIFIC TEs BY CANCER STATUS (KICS + LFS) =====\n")
tryCatch({
  te_aff_unaff_expand_cancer <- te_aff_unaff_expand %>%
    mutate(Cancer_status = ifelse(tumor_type != "U", "Affected", "Unaffected"))

  te_aff_unaff_split_cancer <- te_aff_unaff_split %>%
    mutate(Cancer_status = ifelse(tumor_type != "U", "Affected", "Unaffected"))

  te_aff_unaff_cancer <- te_all %>%
    filter(cohort != "Taylor" & cohort != "HostSeq") %>%
    mutate(Cancer_status = ifelse(tumor_type != "U", "Affected", "Unaffected"))

  # Test by coordinates - GLM
  for (min_samples_cancer in c(3, 5, 10)) {
    cat("\n--- GLM Testing with min_samples =", min_samples_cancer, "---\n")
    te_specific_results_cancer_glm <- test_specific_tes_by_group_glm(
      te_expand = te_aff_unaff_expand_cancer,
      te_count = te_aff_unaff_cancer,
      group_column = "Cancer_status",
      covariates = glm_covariates,
      min_samples_with = min_samples_cancer,
      min_samples_without = min_samples_cancer,
      output_dir = r_dir_files,
      output_prefix = paste0("specific_tes_affunaff_Cancer_fullins_glm_min", min_samples_cancer)
    )
  }

  # Test by gene (all genes) - GLM
  for (min_samples_cancer in c(3, 5, 10)) {
    cat("\n--- GLM Testing gene-grouped with min_samples =", min_samples_cancer, "---\n")
    te_specific_by_gene_cancer_glm <- test_specific_tes_by_gene_glm(
      te_expand = te_aff_unaff_split_cancer,
      te_count = te_aff_unaff_cancer,
      group_column = "Cancer_status",
      covariates = glm_covariates,
      min_samples_with = min_samples_cancer,
      min_samples_without = min_samples_cancer,
      gene_filter = NULL,
      output_dir = r_dir_files,
      output_prefix = paste0("specific_tes_affunaff_Cancer_gene_glm_min", min_samples_cancer)
    )
  }

  # Test by gene (cancer genes only) - GLM
  for (min_samples_cancer in c(3, 5, 10)) {
    cat("\n--- GLM Testing gene-grouped (cancer genes only) with min_samples =", min_samples_cancer, "---\n")
    te_specific_by_gene_cancer_cancer_glm <- test_specific_tes_by_gene_glm(
      te_expand = te_aff_unaff_split_cancer,
      te_count = te_aff_unaff_cancer,
      group_column = "Cancer_status",
      covariates = glm_covariates,
      min_samples_with = min_samples_cancer,
      min_samples_without = min_samples_cancer,
      gene_filter = genes,
      output_dir = r_dir_files,
      output_prefix = paste0("specific_tes_affunaff_Cancer_cancergene_glm_min", min_samples_cancer)
    )
  }

  # Generate plots for significant results (q < 0.1)
  cat("\n--- Generating plots for significant results (q < 0.1) ---\n")
  if (!is.null(te_specific_results_cancer_glm)) {
    p <- plot_specific_te_presence(te_specific_results_cancer_glm, q_threshold = 0.1,
                                   id_column = "te_id",
                                   title = "TE Presence by Cancer Status (KICS+LFS, Coordinates)")
    if (!is.null(p)) ggsave(paste0(plot_dir, "cancer_genes/specific_tes_affunaff_Cancer_fullins_glm_presence.png"), plot = p, width = 10, height = 8)
  }
  if (!is.null(te_specific_by_gene_cancer_glm)) {
    p <- plot_specific_te_presence(te_specific_by_gene_cancer_glm, q_threshold = 0.1,
                                   id_column = "Gene_name",
                                   title = "TE Presence by Cancer Status (KICS+LFS, All Genes)")
    if (!is.null(p)) ggsave(paste0(plot_dir, "cancer_genes/specific_tes_affunaff_Cancer_gene_glm_presence.png"), plot = p, width = 10, height = 8)
  }
  if (!is.null(te_specific_by_gene_cancer_cancer_glm)) {
    p <- plot_specific_te_presence(te_specific_by_gene_cancer_cancer_glm, q_threshold = 0.1,
                                   id_column = "Gene_name",
                                   title = "TE Presence by Cancer Status (KICS+LFS, Cancer Genes)")
    if (!is.null(p)) ggsave(paste0(plot_dir, "cancer_genes/specific_tes_affunaff_Cancer_cancergene_glm_presence.png"), plot = p, width = 10, height = 8)
  }
}, error = function(e) {
  cat("Warning: Could not perform GLM specific TE testing (Cancer - KICS+LFS):", e$message, "\n")
})


#### GLM: SPECIFIC TEs BY CANCER STATUS (LFS only) ####
cat("\n===== GLM: SPECIFIC TEs BY CANCER STATUS (LFS only) =====\n")
tryCatch({
  te_lfs_expand_cancer <- te_lfs_expand %>%
    mutate(Cancer_status = ifelse(tumor_type != "U", "Affected", "Unaffected"))

  te_lfs_split_cancer <- te_lfs_split %>%
    mutate(Cancer_status = ifelse(tumor_type != "U", "Affected", "Unaffected"))

  te_lfs_cancer <- te_lfs %>%
    mutate(Cancer_status = ifelse(tumor_type != "U", "Affected", "Unaffected"))

  # Test by coordinates - GLM
  for (min_samples_cancer in c(3, 5, 10)) {
    cat("\n--- GLM Testing with min_samples =", min_samples_cancer, "---\n")
    te_specific_results_lfs_cancer_glm <- test_specific_tes_by_group_glm(
      te_expand = te_lfs_expand_cancer,
      te_count = te_lfs_cancer,
      group_column = "Cancer_status",
      covariates = glm_covariates,
      min_samples_with = min_samples_cancer,
      min_samples_without = min_samples_cancer,
      output_dir = r_dir_files,
      output_prefix = paste0("specific_tes_lfs_Cancer_fullins_glm_min", min_samples_cancer)
    )
  }

  # Test by gene (all genes) - GLM
  for (min_samples_cancer in c(3, 5, 10)) {
    cat("\n--- GLM Testing gene-grouped with min_samples =", min_samples_cancer, "---\n")
    te_specific_by_gene_lfs_cancer_glm <- test_specific_tes_by_gene_glm(
      te_expand = te_lfs_split_cancer,
      te_count = te_lfs_cancer,
      group_column = "Cancer_status",
      covariates = glm_covariates,
      min_samples_with = min_samples_cancer,
      min_samples_without = min_samples_cancer,
      gene_filter = NULL,
      output_dir = r_dir_files,
      output_prefix = paste0("specific_tes_lfs_Cancer_gene_glm_min", min_samples_cancer)
    )
  }

  # Test by gene (cancer genes only) - GLM
  for (min_samples_cancer in c(3, 5, 10)) {
    cat("\n--- GLM Testing gene-grouped (cancer genes only) with min_samples =", min_samples_cancer, "---\n")
    te_specific_by_gene_cancer_lfs_cancer_glm <- test_specific_tes_by_gene_glm(
      te_expand = te_lfs_split_cancer,
      te_count = te_lfs_cancer,
      group_column = "Cancer_status",
      covariates = glm_covariates,
      min_samples_with = min_samples_cancer,
      min_samples_without = min_samples_cancer,
      gene_filter = genes,
      output_dir = r_dir_files,
      output_prefix = paste0("specific_tes_lfs_Cancer_cancergene_glm_min", min_samples_cancer)
    )
  }
  # Generate plots for significant results (q < 0.1) - LFS Cancer status
  cat("\n--- Generating plots for significant results (q < 0.1) - Cancer Status (LFS only) ---\n")
  if (exists("te_specific_results_lfs_cancer_glm") && !is.null(te_specific_results_lfs_cancer_glm)) {
    p_fullins <- plot_specific_te_presence(te_specific_results_lfs_cancer_glm, q_threshold = 0.1,
                                           id_column = "te_id",
                                           title = "TE Presence by Cancer Status (LFS, Coordinates)")
    if (!is.null(p_fullins)) {
      ggsave(paste0(plot_dir, "cancer_genes/specific_tes_lfs_Cancer_fullins_glm_presence.png"),
             plot = p_fullins, width = 10, height = 8)
      cat("Saved coordinate-level presence plot\n")
    }
  }
  if (exists("te_specific_by_gene_lfs_cancer_glm") && !is.null(te_specific_by_gene_lfs_cancer_glm)) {
    p_gene <- plot_specific_te_presence(te_specific_by_gene_lfs_cancer_glm, q_threshold = 0.1,
                                        id_column = "Gene_name",
                                        title = "TE Presence by Cancer Status (LFS, All Genes)")
    if (!is.null(p_gene)) {
      ggsave(paste0(plot_dir, "cancer_genes/specific_tes_lfs_Cancer_gene_glm_presence.png"),
             plot = p_gene, width = 10, height = 8)
      cat("Saved gene-level presence plot\n")
    }
  }
  if (exists("te_specific_by_gene_cancer_lfs_cancer_glm") && !is.null(te_specific_by_gene_cancer_lfs_cancer_glm)) {
    p_cancer_gene <- plot_specific_te_presence(te_specific_by_gene_cancer_lfs_cancer_glm, q_threshold = 0.1,
                                               id_column = "Gene_name",
                                               title = "TE Presence by Cancer Status (LFS, Cancer Genes)")
    if (!is.null(p_cancer_gene)) {
      ggsave(paste0(plot_dir, "cancer_genes/specific_tes_lfs_Cancer_cancergene_glm_presence.png"),
             plot = p_cancer_gene, width = 10, height = 8)
      cat("Saved cancer gene-level presence plot\n")
    }
  }
}, error = function(e) {
  cat("Warning: Could not perform GLM specific TE testing (Cancer - LFS only):", e$message, "\n")
})


#### GLM: SPECIFIC TEs BY COHORT (KICS vs HostSeq) ####
cat("\n===== GLM: SPECIFIC TEs BY COHORT (KICS vs HostSeq) =====\n")
tryCatch({
  # Create count matrix for kics + hostseq by combining the two cohort matrices
  # (te_all doesn't contain HostSeq samples, so we need to rbind them manually)
  te_kics_hostseq_count <- rbind(te_kics, te_hostseq)

  # Test by coordinates - GLM
  for (min_samples in c(3, 5, 10)) {
    cat("\n--- GLM Testing with min_samples =", min_samples, "---\n")
    te_specific_results_kics_hostseq_glm <- test_specific_tes_by_group_glm(
      te_expand = te_kics_hostseq_expand,
      te_count = te_kics_hostseq_count,
      group_column = "cohort",
      covariates = glm_covariates,
      min_samples_with = min_samples,
      min_samples_without = min_samples,
      output_dir = r_dir_files,
      output_prefix = paste0("specific_tes_kics_hostseq_fullins_glm_min", min_samples)
    )
  }

  # Test by gene (all genes) - GLM
  for (min_samples in c(3, 5, 10)) {
    cat("\n--- GLM Testing gene-grouped with min_samples =", min_samples, "---\n")
    te_specific_by_gene_kics_hostseq_glm <- test_specific_tes_by_gene_glm(
      te_expand = te_kics_hostseq,
      te_count = te_kics_hostseq_count,
      group_column = "cohort",
      covariates = glm_covariates,
      min_samples_with = min_samples,
      min_samples_without = min_samples,
      gene_filter = NULL,
      output_dir = r_dir_files,
      output_prefix = paste0("specific_tes_kics_hostseq_gene_glm_min", min_samples)
    )
  }

  # Test by gene (cancer genes only) - GLM
  for (min_samples in c(3, 5, 10)) {
    cat("\n--- GLM Testing gene-grouped (cancer genes only) with min_samples =", min_samples, "---\n")
    te_specific_by_gene_cancer_kics_hostseq_glm <- test_specific_tes_by_gene_glm(
      te_expand = te_kics_hostseq,
      te_count = te_kics_hostseq_count,
      group_column = "cohort",
      covariates = glm_covariates,
      min_samples_with = min_samples,
      min_samples_without = min_samples,
      gene_filter = genes,
      output_dir = r_dir_files,
      output_prefix = paste0("specific_tes_kics_hostseq_cancergene_glm_min", min_samples)
    )
  }

  # Generate plots for significant results (q < 0.1) - KICS vs HostSeq
  cat("\n--- Generating plots for significant results (q < 0.1) - KICS vs HostSeq ---\n")
  if (exists("te_specific_results_kics_hostseq_glm") && !is.null(te_specific_results_kics_hostseq_glm)) {
    p_fullins <- plot_specific_te_presence(te_specific_results_kics_hostseq_glm, q_threshold = 0.1,
                                           id_column = "te_id",
                                           title = "TE Presence by Cohort (KICS vs HostSeq, Coordinates)")
    if (!is.null(p_fullins)) {
      ggsave(paste0(plot_dir, "cancer_genes/specific_tes_kics_hostseq_fullins_glm_presence.png"),
             plot = p_fullins, width = 10, height = 8)
      cat("Saved coordinate-level presence plot\n")
    }
  }
  if (exists("te_specific_by_gene_kics_hostseq_glm") && !is.null(te_specific_by_gene_kics_hostseq_glm)) {
    p_gene <- plot_specific_te_presence(te_specific_by_gene_kics_hostseq_glm, q_threshold = 0.1,
                                        id_column = "Gene_name",
                                        title = "TE Presence by Cohort (KICS vs HostSeq, All Genes)")
    if (!is.null(p_gene)) {
      ggsave(paste0(plot_dir, "cancer_genes/specific_tes_kics_hostseq_gene_glm_presence.png"),
             plot = p_gene, width = 10, height = 8)
      cat("Saved gene-level presence plot\n")
    }
  }
  if (exists("te_specific_by_gene_cancer_kics_hostseq_glm") && !is.null(te_specific_by_gene_cancer_kics_hostseq_glm)) {
    p_cancer_gene <- plot_specific_te_presence(te_specific_by_gene_cancer_kics_hostseq_glm, q_threshold = 0.1,
                                               id_column = "Gene_name",
                                               title = "TE Presence by Cohort (KICS vs HostSeq, Cancer Genes)")
    if (!is.null(p_cancer_gene)) {
      ggsave(paste0(plot_dir, "cancer_genes/specific_tes_kics_hostseq_cancergene_glm_presence.png"),
             plot = p_cancer_gene, width = 10, height = 8)
      cat("Saved cancer gene-level presence plot\n")
    }
  }
}, error = function(e) {
  cat("Warning: Could not perform GLM specific TE testing (KICS vs HostSeq):", e$message, "\n")
})


cat("✓ Script completed successfully\n")

# Close module-specific sink
close_module_sink()
