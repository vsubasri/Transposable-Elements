#!/usr/bin/env Rscript

# Tumour TE Visualization - Specific TEs
# Specific TEs by TP53 status

# Source common setup and load data
source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/00_viz_common_setup.R")
REQUIRED_DATA <- c("count_matrix", "expand", "split", "clinical", "genes")
source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/00_viz_load_data_tumour.R")

# Initialize module-specific text output
init_module_sink(paste0(plot_dir, "cancer_genes/"), "SPECIFIC_TES")

cat("Running 02_te_viz_tumour_03_specific_tes.R...\n")

# Load cancer predisposition genes
cpg <- read.csv("/Users/briannelaverty/Documents/R_Malkin/te/data/raw/kics_cpg.txt", header=FALSE, sep="\t")
genes <- cpg$V1

#### PATHWAY ANALYSIS - TEs WITH EXPRESSION EFFECTS (p < 0.25) ####
cat("\n===== UNIQUE TE IN LFS ANALYSIS =====\n")
tryCatch({
  specific_te_fisher <- as.data.frame(fisher_test_unique_te(te_aff_expand_t, min=3))
  write_output(quote(head(specific_te_fisher)), "Head of specific_te_fisher (UNIQUE TE IN LFS)")
  sig_te <- specific_te_fisher %>% filter(fisher_p_value_BH < 0.1)
  write_output(quote(sig_te), "Significant TE (FDR < 0.1) in LFS")
  # samples with significant TE
  sig_te_samples <- find_sig_te_samples(sig_te, te_aff_expand_t)
  write_output(quote(lapply(sig_te_samples, print_summary_sig_te_samples)), "Summary of samples with significant TE (LFS)")
  # genes effected
  sig_te_genes <- find_sig_te_genes(sig_te, te_aff_split_genes_t)
  write_output(quote(sig_te_genes), "Genes affected by significant TE (LFS)")
  if (exists("plot_sig_te_genes")) {
    titled_print(plot_sig_te_genes(sig_te_genes), "Significant TE Genes in LFS")
    ggsave(paste0(plot_dir, "counts_clinical_lfs/te_sig_te_genes_lfs.png"), width = 8, height = 5)
  }
}, error = function(e) {
  cat("Warning: Could not perform unique TE in LFS analysis:", e$message, "\n")
})

cat("\n===== GERMLINE VARIANTS ANALYSIS =====\n")
tryCatch({
  kics_germline_variants <- fread("/Users/briannelaverty/Documents/R_Malkin/TE/data/kics_germline_variants_panel.csv", sep = ",", header = TRUE)
  kics_germline_variants$sample <- sprintf("%04d", as.numeric(kics_germline_variants$KiCS_ID))
  kics_germline_variants <- subset(kics_germline_variants, sample %in% te_kics_t$base_sample)
  kics_germline_variants <- subset(kics_germline_variants, interpretation %in% c("Pathogenic", "Likely Pathogenic"))
  kics_germline_variants <- kics_germline_variants %>% rename("gene" = "geneSymbol")
  write_output(quote(table(kics_germline_variants$gene)), "Number of samples with P/LP germline variants by gene")
  te_variant_df <- add_variant_columns(data = te_kics_t, variant_data = kics_germline_variants, variant_type = "germline")
  write_output(quote(filter_and_print_g_columns(te_variant_df)), "Summary of germline variant columns in TE data")
  high <- subset(te_variant_df, sample %in% te_aff_t[te_aff_t$LINE1>100, "sample"])
  write_output(quote(high[,c(1,144:162)]), "High LINE1 samples with germline variant columns")
  te_aff_variant_t <- prep_variant_df(te_kics_t, te_aff_t)
  generate_plots(plot_count_kruskal, te_aff_variant_t, chr=NA, type=types, group="TP53", x_lab="Cluster")
  for (i in seq_along(types)) {
    titled_print(plot_count_kruskal(te_aff_variant_t, group="TP53", x_order=c("None", "Somatic", "Germline"), x_lab="TP53 variant", y_lab="Total TE count", type=types[i], chr=NA, log_scale=TRUE), 
                 paste("TE count by TP53 variant status -", types[i]))
    ggsave(paste0(plot_dir, "counts_clinical_lfs/te_count_tp53_variant_", types[i], ".png"), width = 9, height = 5)
  }
}, error = function(e) {
  cat("Warning: Could not perform germline variants analysis:", e$message, "\n")
})

cat("\n===== TP53 VARIANTS COUNT ANALYSIS =====\n")
tryCatch({
  kics_tp53_somatic_variants <- fread("/Users/briannelaverty/Documents/R_Malkin/TE/data/kics_tp53_somatic_variants.csv", sep=",", header=TRUE)
  somatic_conversion <- read.delim("/Users/briannelaverty/Documents/R_Malkin/te/data/final/somatic_variant_conversion.txt", sep="\t", header=TRUE)
  kics_tp53_somatic_variants_processed <- prep_variants_tumour(kics_tp53_somatic_variants, somatic_conversion, te_kics_t)
  te_aff_variant_t <- add_variant_columns(data=te_aff_t, variant_data=kics_tp53_somatic_variants_processed, variant_type="somatic")
  out <- te_aff_variant_t[te_aff_variant_t$s_TP53==1, "sample"]
  var <- unique(kics_tp53_somatic_variants_processed$sample)
  write_output(quote(var[!var%in% out]), "Samples with variant and base sample match but no TE calls")
  te_aff_variant_t <- te_aff_variant_t %>%
    mutate(TP53 = case_when(
      TP53_status == "Mutant" ~ "Germline",
      s_TP53 == 1 ~ "Somatic",
      TRUE ~ "WT"
    ))

  # Generate plots for all TE types
  for (i in seq_along(types)) {
    write_output(quote(plot_count_kruskal(te_aff_variant_t, group="TP53", x_order=c("WT", "Somatic", "Germline"), x_lab="TP53 variant", y_lab="Repeat count", type=types[i], chr=NA, log_scale=TRUE)),
                 paste0("TP53 variant (Somatic/Germline/WT) count (type=", ifelse(is.na(types[i]), "all", types[i]), ")"))
    titled_print(plot_count_kruskal(te_aff_variant_t, group="TP53", x_order=c("WT", "Somatic", "Germline"), x_lab="TP53 variant", y_lab="Repeat count", type=types[i], chr=NA, log_scale=TRUE),
                 paste0("TP53 variant (Somatic/Germline/WT) count (type=", ifelse(is.na(types[i]), "all", types[i]), ")"))
    ggsave(paste0(plot_dir, "counts_clinical_lfs/te_count_tp53_somatic_germline_", ifelse(is.na(types[i]), "all", types[i]), ".png"),
           plot=plot_count_kruskal(te_aff_variant_t, group="TP53", x_order=c("WT", "Somatic", "Germline"), x_lab="TP53 variant", y_lab="Repeat count", type=types[i], chr=NA, log_scale=TRUE), width=9, height=5)
  }
}, error = function(e) {
  cat("Warning: Could not perform TP53 variants count analysis:", e$message, "\n")
})

cat("\n===== MERGED TUMOURS ANALYSIS =====\n")
tryCatch({
  clonality <- read.csv("/Users/briannelaverty/Documents/R_Malkin/clinical/lfs_clonality.csv", header=TRUE)
  plot_multisample_scatter_nick(te_all_all_t, te_count_col ="total")
  plot_multisample_scatter_clonality(te_all_all_t, clonal_df=clonality, clonal_column="ssm_prop_clonal", legend_lab = "SSM clonality", te_count_col ="total")
  titled_print(plot_multisample_scatter_clonality(te_all_all_t, clonal_df=clonality, clonal_column="ssm_prop_clonal", legend_lab = "SSM clonality", te_count_col ="total"), "Merged Tumours: SSM clonality")
  ggsave(paste0(plot_dir, "counts_clinical_lfs/te_tumour_lfs_clonality.png"), width = 9, height = 5)
}, error = function(e) {
  cat("Warning: Could not perform merged tumours analysis:", e$message, "\n")
})

write_output(quote(NULL), "MULTISAMPLE ANALYSIS")
tryCatch({
  write_output(
    quote(plot_multisample_scatter_kics(te_all_all_t, te_count_col ="total")),
    "KICS Multisample Scatter"
  )
  titled_print(plot_multisample_scatter_kics(te_all_all_t, te_count_col ="total"), "KICS Multisample Scatter")
  ggsave(paste0(plot_dir, "counts_clinical_kics/te_multisample_kics.png"), width = 9, height = 5)

  write_output(quote(plot_individual_patient_samples(te_all_all_t)), "Individual Patient Samples")
  write_output(quote(plot_individual_patient_samples_facet(te_all_all_t)), "Individual Patient Samples (Faceted)")
  write_output(quote(plot_specific_sample(te_all_all_t, "0074", x_nudge=-100, y_nudge=1, log_scale=FALSE)), "Specific Sample 0074")
  write_output(quote(plot_multisample_sametime(te_all_all_t)), "Multisample Same Time")
  write_output(quote(plot_te_change(te_all_all, log_scale=TRUE)), "TE Change by Sample (KICS)")
  titled_print(plot_te_change(te_all_all, log_scale=TRUE), "TE change by sample (KICS)")
  ggsave(paste0(plot_dir, "counts_clinical_kics/te_tumour_kics_change.png"), width = 9, height = 5)
}, error = function(e) {
  cat("Warning: Could not perform multisample analysis:", e$message, "\n")
})

# NOTE: TE source analysis moved to 02_te_viz_tumour_14_fulllength_young_source.R

#### SPECIFIC TEs BY TP53 STATUS ####
cat("\n===== SPECIFIC TEs BY TP53 STATUS =====\n")

# Ensure output directory exists
specific_tes_dir <- paste0(plot_dir, "specific_tes/")
if (!dir.exists(specific_tes_dir)) {
  dir.create(specific_tes_dir, recursive = TRUE)
  cat("Created directory:", specific_tes_dir, "\n")
}

tryCatch({
  # Test by insertion (fullins)
  cat("\n--- Testing BY INSERTION (fullins) ---\n")
  for (min_samples in c(3, 5, 10)) {
    cat("\n--- Testing with min_samples =", min_samples, "---\n")
    te_specific_results_tp53 <- test_specific_tes_by_group(
      te_expand = te_aff_expand_t,
      te_count = te_aff_t,
      group_column = "TP53_status",
      min_samples_with = min_samples,
      min_samples_without = min_samples,
      output_dir = paste0(plot_dir, "specific_tes/"),
      output_prefix = paste0("specific_tes_aff_TP53_fullins_min", min_samples)
    )
  }

  if (!is.null(te_specific_results_tp53)) {
    write_output(quote(head(te_specific_results_tp53$full_results, 20)), "Top 20 TEs by adjusted p-value (TP53)")
  }

  # Test by gene (all genes)
  cat("\n--- Testing BY GENE (all genes) ---\n")
  for (min_samples in c(3, 5, 10)) {
    cat("\n--- Testing with min_samples =", min_samples, "---\n")
    te_specific_by_gene_tp53 <- test_specific_tes_by_gene(
      te_expand = te_aff_split_t,
      te_count = te_aff_t,
      group_column = "TP53_status",
      min_samples_with = min_samples,
      min_samples_without = min_samples,
      gene_filter = NULL,
      output_dir = paste0(plot_dir, "specific_tes/"),
      output_prefix = paste0("specific_tes_aff_TP53_gene_min", min_samples)
    )
  }

  # Test by gene (cancer genes only)
  cat("\n--- Testing BY GENE (cancer genes only) ---\n")
  for (min_samples in c(3, 5, 10)) {
    cat("\n--- Testing with min_samples =", min_samples, "---\n")
    te_specific_by_gene_cancer_tp53 <- test_specific_tes_by_gene(
      te_expand = te_aff_split_t,
      te_count = te_aff_t,
      group_column = "TP53_status",
      min_samples_with = min_samples,
      min_samples_without = min_samples,
      gene_filter = genes,
      output_dir = paste0(plot_dir, "specific_tes/"),
      output_prefix = paste0("specific_tes_aff_TP53_cancergene_min", min_samples)
    )
  }
}, error = function(e) {
  cat("Warning: Could not perform specific TE testing (TP53):", e$message, "\n")
})


#### SPECIFIC TEs BY CANCER STATUS (LFS) ####
# Note: Skipped for tumor samples - all tumor samples have cancer by definition
# Cancer status comparison only makes sense in germline data
cat("\n===== SPECIFIC TEs BY CANCER STATUS (LFS) =====\n")
cat("Skipped: Not applicable for tumor samples (all have cancer)\n")


#### SPECIFIC TEs BY TUMOR TYPE (One-vs-Rest) ####
cat("\n===== SPECIFIC TEs BY TUMOR TYPE (One-vs-Rest) =====\n")

# Check if there are enough tumor types with sufficient samples
min_samples_tt <- 5
min_samples_te_thresholds <- c(3, 5)

tumor_type_check <- te_aff_t %>%
  filter(!is.na(tumor_type)) %>%
  count(tumor_type) %>%
  filter(n >= min_samples_tt)

if (nrow(tumor_type_check) < 2) {
  # Write explanation file
  reason_file <- paste0(plot_dir, "specific_tes/specific_tes_tumourtype_NOT_RUN.txt")
  all_tumor_types <- te_aff_t %>%
    filter(!is.na(tumor_type)) %>%
    count(tumor_type)

  writeLines(c(
    "SPECIFIC TEs BY TUMOR TYPE (One-vs-Rest) - NOT RUN",
    "",
    paste("Reason: Insufficient tumor types with >=", min_samples_tt, "samples"),
    "",
    "Sample counts by tumor type:",
    paste(capture.output(print(all_tumor_types)), collapse = "\n"),
    "",
    paste("Tumor types with >=", min_samples_tt, "samples:", nrow(tumor_type_check)),
    paste(capture.output(print(tumor_type_check)), collapse = "\n"),
    "",
    "Requirements for tumor type testing:",
    paste("  1. At least 2 tumor types with >=", min_samples_tt, "samples each (allows one-vs-rest comparison)"),
    paste("  2. Individual TEs/genes must have >=", paste(min_samples_te_thresholds, collapse = " or "), "samples in the tumor type of interest"),
    "",
    "Note: Even if tumor types have sufficient samples, specific tests may return no results",
    "if no TEs/genes meet the minimum sample threshold within those tumor types.",
    "",
    paste("Date:", Sys.time())
  ), reason_file)
  cat("✗ Skipped tumor type testing - insufficient tumor types\n")
  cat("  Explanation saved to:", basename(reason_file), "\n")
} else {
  tryCatch({
    # Track if any results were found
    any_results_found <- FALSE

    # Test by insertion (fullins)
    cat("\n--- Testing BY INSERTION (fullins) ---\n")
    for (min_samples_te in c(3, 5)) {
      cat("\n--- Testing with min_samples_te =", min_samples_te, "---\n")
      te_specific_results_tt <- test_specific_tes_by_tumor_type(
        te_expand = te_aff_expand_t,
        te_count = te_aff_t,
        min_samples_tt = 5,
        min_samples_te = min_samples_te,
        output_dir = paste0(plot_dir, "specific_tes/"),
        output_prefix = paste0("specific_tes_aff_tumourtype_fullins_min", min_samples_te)
      )
      if (!is.null(te_specific_results_tt) && is.data.frame(te_specific_results_tt) && nrow(te_specific_results_tt) > 0) {
        any_results_found <- TRUE
      }
    }

    # Test by gene (all genes)
    cat("\n--- Testing BY GENE (all genes) ---\n")
    for (min_samples_te in c(3, 5)) {
      cat("\n--- Testing with min_samples_te =", min_samples_te, "---\n")
      te_specific_by_gene_tt <- test_specific_tes_by_tumor_type_gene(
        te_expand = te_aff_split_t,
        te_count = te_aff_t,
        min_samples_tt = 5,
        min_samples_te = min_samples_te,
        gene_filter = NULL,
        output_dir = paste0(plot_dir, "specific_tes/"),
        output_prefix = paste0("specific_tes_aff_tumourtype_gene_min", min_samples_te)
      )
      if (!is.null(te_specific_by_gene_tt) && is.data.frame(te_specific_by_gene_tt) && nrow(te_specific_by_gene_tt) > 0) {
        any_results_found <- TRUE
      }
    }

    # Test by gene (cancer genes only)
    cat("\n--- Testing BY GENE (cancer genes only) ---\n")
    for (min_samples_te in c(3, 5)) {
      cat("\n--- Testing with min_samples_te =", min_samples_te, "---\n")
      te_specific_by_gene_cancer_tt <- test_specific_tes_by_tumor_type_gene(
        te_expand = te_aff_split_t,
        te_count = te_aff_t,
        min_samples_tt = 5,
        min_samples_te = min_samples_te,
        gene_filter = genes,
        output_dir = paste0(plot_dir, "specific_tes/"),
        output_prefix = paste0("specific_tes_aff_tumourtype_cancergene_min", min_samples_te)
      )
      if (!is.null(te_specific_by_gene_cancer_tt) && is.data.frame(te_specific_by_gene_cancer_tt) && nrow(te_specific_by_gene_cancer_tt) > 0) {
        any_results_found <- TRUE
      }
    }

    # If no results found across all tests, create explanation file
    cat("\nChecking if any results were found... any_results_found =", any_results_found, "\n")

    if (!any_results_found) {
      all_tumor_types <- te_aff_t %>%
        filter(!is.na(tumor_type)) %>%
        count(tumor_type)

      no_results_file <- paste0(plot_dir, "specific_tes/specific_tes_tumourtype_NO_RESULTS.txt")
      cat("Creating NO_RESULTS file:", no_results_file, "\n")

      writeLines(c(
        "SPECIFIC TEs BY TUMOR TYPE (One-vs-Rest) - NO RESULTS FOUND",
        "",
        "Tests were run successfully but returned no significant results.",
        "",
        paste("Reason: No TEs or genes had >=", paste(min_samples_te_thresholds, collapse = " or "), "samples in any tumor type"),
        "",
        "Sample counts by tumor type (eligible types with >= 5 samples):",
        paste(capture.output(print(tumor_type_check)), collapse = "\n"),
        "",
        "All tumor types:",
        paste(capture.output(print(all_tumor_types)), collapse = "\n"),
        "",
        "What was tested:",
        "  - Insertion-level tests (fullins) with min_samples_te = 3 and 5",
        "  - Gene-level tests (all genes) with min_samples_te = 3 and 5",
        "  - Gene-level tests (cancer genes) with min_samples_te = 3 and 5",
        "",
        "Note: While tumor types had sufficient samples for testing, no individual",
        "TEs or genes met the minimum sample threshold within those tumor types.",
        "This suggests TEs are too rare or distributed across tumor types.",
        "",
        paste("Date:", Sys.time())
      ), no_results_file)
      cat("✓ Tests completed but no results found\n")
      cat("  Explanation saved to:", basename(no_results_file), "\n")
    } else {
      cat("✓ Some results were found - no NO_RESULTS file needed\n")
    }

  }, error = function(e) {
    cat("Warning: Could not perform specific TE testing (tumor type):", e$message, "\n")
    # Write error explanation file
    all_tumor_types <- te_aff_t %>%
      filter(!is.na(tumor_type)) %>%
      count(tumor_type)
    error_file <- paste0(plot_dir, "specific_tes/specific_tes_tumourtype_ERROR.txt")
    writeLines(c(
      "SPECIFIC TEs BY TUMOR TYPE (One-vs-Rest) - ERROR",
      "",
      paste("Error:", e$message),
      "",
      "Sample counts by tumor type:",
      paste(capture.output(print(all_tumor_types)), collapse = "\n"),
      "",
      paste("Date:", Sys.time())
    ), error_file)
  })
}


cat("✓ Script completed successfully\n")

# Close module-specific sink
close_module_sink()
