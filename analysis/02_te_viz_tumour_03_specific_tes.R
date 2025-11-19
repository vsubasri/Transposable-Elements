#!/usr/bin/env Rscript

# Tumour TE Visualization - Specific TEs
# Specific TEs by TP53 status

# Source common setup and load data
source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/00_viz_common_setup.R")
source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/00_viz_load_data_tumour.R")

cat("Running 02_te_viz_tumour_03_specific_tes.R...\n")

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
    ggsave(paste0(plot_dir, "other/te_sig_te_genes_lfs.png"), width = 8, height = 5)
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
  ggsave(paste0(plot_dir, "other/te_tumour_lfs_clonality.png"), width = 9, height = 5)
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
  ggsave(paste0(plot_dir, "other/te_multisample_kics.png"), width = 9, height = 5)

  write_output(quote(plot_individual_patient_samples(te_all_all_t)), "Individual Patient Samples")
  write_output(quote(plot_individual_patient_samples_facet(te_all_all_t)), "Individual Patient Samples (Faceted)")
  write_output(quote(plot_specific_sample(te_all_all_t, "0074", x_nudge=-100, y_nudge=1, log_scale=FALSE)), "Specific Sample 0074")
  write_output(quote(plot_multisample_sametime(te_all_all_t)), "Multisample Same Time")
  write_output(quote(plot_te_change(te_all_all, log_scale=TRUE)), "TE Change by Sample (KICS)")
  titled_print(plot_te_change(te_all_all, log_scale=TRUE), "TE change by sample (KICS)")
  ggsave(paste0(plot_dir, "other/te_tumour_kics_change.png"), width = 9, height = 5)
}, error = function(e) {
  cat("Warning: Could not perform multisample analysis:", e$message, "\n")
})

cat("\n===== TE SOURCE ANALYSIS =====\n")
tryCatch({
  te_aff_expand_line_t <- te_aff_expand_t %>% filter(ALT=="LINE1")
  te_aff_expand_line_t <- extract_info_fields(te_aff_expand_line_t)

  cat("Number of unique LINE1 sources (affected):", nrow(table(te_aff_expand_line_t$source)), "\n")
  cat("Number of LINE1 insertions (affected):", nrow(te_aff_expand_line_t), "\n")

  sources_with_multiple <- table(te_aff_expand_line_t$source)[table(te_aff_expand_line_t$source) > 1]
  if (length(sources_with_multiple) > 0) {
    cat("\nLINE1 sources with >1 occurrence (affected):\n")
    print(sources_with_multiple)
  } else {
    cat("\nNo LINE1 sources with >1 occurrence\n")
  }

  sources_transduction <- te_aff_expand_line_t %>% filter(source != "not_transduction")
  if (nrow(sources_transduction) > 0) {
    samples_with_sources <- table(sources_transduction[, c("sample", "source")])
    cat("\nSamples with LINE1 sources (transductions):\n")
    print(head(samples_with_sources, 50))
  } else {
    cat("\nNo LINE1 transductions found\n")
  }

  sources <- te_aff_expand_line_t %>%
    filter(source != "not_transduction") %>%
    count(sample, source) %>%
    filter(n > 0) %>%
    arrange(desc(n))

  if (nrow(sources) > 0) {
    cat("\nSummary of LINE1 sources per sample (affected):\n")
    print(head(sources, 20))

    # Create summary by source with TP53 status breakdown
    sources_summary <- te_aff_expand_line_t %>%
      filter(source != "not_transduction") %>%
      group_by(source) %>%
      summarise(
        n = n(),
        n_TP53_mut = sum(TP53_status == "Mutant", na.rm = TRUE),
        n_TP53_wt = sum(TP53_status == "WT", na.rm = TRUE),
        samples_TP53_mut = paste(unique(sample[!is.na(TP53_status) & TP53_status == "Mutant"]), collapse = ";"),
        samples_TP53_wt = paste(unique(sample[!is.na(TP53_status) & TP53_status == "WT"]), collapse = ";"),
        all_samples = paste(unique(sample), collapse = ";"),
        .groups = "drop"
      ) %>%
      # Replace empty strings with NA for cleaner output
      mutate(
        samples_TP53_mut = ifelse(samples_TP53_mut == "", NA_character_, samples_TP53_mut),
        samples_TP53_wt = ifelse(samples_TP53_wt == "", NA_character_, samples_TP53_wt)
      ) %>%
      arrange(desc(n)) %>%
      select(source, n, n_TP53_mut, n_TP53_wt, samples_TP53_mut, samples_TP53_wt, all_samples)

    # Save source analysis to CSV
    write.csv(sources_summary, paste0(r_dir_files, "line1_source_analysis_tumour.csv"), row.names = FALSE)
    cat("✓ Source analysis saved to: line1_source_analysis_tumour.csv\n")
  } else {
    cat("\nNo LINE1 sources found\n")
  }

  if (exists("plot_te_source")) {
    titled_print(plot_te_source(sources), "TE Source Plot (affected)")
    ggsave(paste0(plot_dir, "other/te_source_affected.png"), width = 7, height = 5)
  }
}, error = function(e) {
  cat("Warning: Could not perform TE source analysis:", e$message, "\n")
})

cat("\n===== SPECIFIC TEs BY TP53 STATUS =====\n")
tryCatch({
  # Test specific TE insertions for differential representation by TP53 status
  # Run with multiple thresholds
  for (min_samples in c(3, 5, 10)) {
    cat("\n--- Testing with min_samples =", min_samples, "---\n")
    te_specific_results <- test_specific_tes_by_group(
      te_expand = te_aff_expand_t,
      te_count = te_aff_t,
      group_column = "TP53_status",
      min_samples_with = min_samples,
      min_samples_without = min_samples,
      output_dir = r_dir_files,
      output_prefix = paste0("specific_tes_tumour_min", min_samples)
    )
  }

  if (!is.null(te_specific_results)) {
    write_output(quote(head(te_specific_results$full_results, 20)), "Top 20 TEs by adjusted p-value")
  }
}, error = function(e) {
  cat("Warning: Could not perform specific TE testing:", e$message, "\n")
})



cat("✓ Script completed successfully\n")
