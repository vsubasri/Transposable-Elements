#!/usr/bin/env Rscript

# Tumour TE Visualization - PCA/UMAP
# Ancestry analysis, PCA/UMAP by various variables

# Source common setup and load data
source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/00_viz_common_setup.R")
source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/00_viz_load_data_tumour.R")

cat("Running 02_te_viz_tumour_09_pca_umap.R...\n")

##### ANCESTRY ANALYSIS ####
#cat("\n========================================\n")
#cat("ANCESTRY ANALYSIS\n")
#cat("========================================\n\n")
#
## Load ancestry data
#cat("Loading ancestry data...\n")
#load(paste0(r_dir, "ancestry.RData"))
#
## Load location windows data
#cat("Loading location windows data...\n")
#location_100kb_t <- read.csv("/hpf/largeprojects/davidm/blaverty/te/ml/output/tumour/location_csv/100kb_complete_filtered_t.csv", stringsAsFactors = FALSE)
#
## Merge location windows with ancestry using same approach as processing script
#location_ancestry <- write_output(
#  quote(merge_location_ancestry(location_100kb_t, ancestry)),
#  "Merging location windows with ancestry"
#)
#
## Identify clinical columns to exclude from PCA/UMAP
#clinical_cols <- c("sample", "predicted_ancestry_thres", "mapped_label", "age", "sex", "cohort", "tumor_type",
#                   "TP53_germline", "affected", "cancer", "age_diagnosis", "TP53_status")
#exclude_cols <- intersect(clinical_cols, colnames(location_ancestry))
#
## Perform PCA on location windows
#pca_results <- write_output(
#  quote(perform_location_pca(location_ancestry, exclude_cols = exclude_cols)),
#  "Performing PCA on location windows"
#)
#
## Perform UMAP on location windows
#umap_results <- write_output(
#  quote(perform_location_umap(location_ancestry, exclude_cols = exclude_cols)),
#  "Performing UMAP on location windows"
#)
#
## Plot PCA colored by predicted_ancestry_thres
#p_pca_ancestry <- write_output(
#  quote(plot_pca_ancestry(pca_results, color_by = "predicted_ancestry_thres",
#                         output_dir = plot_dir, plot_prefix = "tumour")),
#  "Plotting PCA colored by predicted ancestry"
#)
#titled_print(p_pca_ancestry, "PCA of TE Location Windows by Predicted Ancestry (Tumour)")
#
## Plot PCA colored by mapped_label
#p_pca_mapped <- write_output(
#  quote(plot_pca_ancestry(pca_results, color_by = "mapped_label",
#                         output_dir = plot_dir, plot_prefix = "tumour")),
#  "Plotting PCA colored by mapped label"
#)
#titled_print(p_pca_mapped, "PCA of TE Location Windows by Mapped Label (Tumour)")
#
## Plot UMAP colored by predicted_ancestry_thres
#p_umap_ancestry <- write_output(
#  quote(plot_umap_ancestry(umap_results, color_by = "predicted_ancestry_thres",
#                          output_dir = plot_dir, plot_prefix = "tumour")),
#  "Plotting UMAP colored by predicted ancestry"
#)
#titled_print(p_umap_ancestry, "UMAP of TE Location Windows by Predicted Ancestry (Tumour)")
#
## Plot UMAP colored by mapped_label
#p_umap_mapped <- write_output(
#  quote(plot_umap_ancestry(umap_results, color_by = "mapped_label",
#                          output_dir = plot_dir, plot_prefix = "tumour")),
#  "Plotting UMAP colored by mapped label"
#)
#titled_print(p_umap_mapped, "UMAP of TE Location Windows by Mapped Label (Tumour)")
#
# Plot count LM grouped by predicted_ancestry_thres
# Check if we have enough samples with ancestry data
#if ("predicted_ancestry_thres" %in% colnames(te_aff_t) &&
#    sum(!is.na(te_aff_t$predicted_ancestry_thres)) >= 10) {
#
#  ancestry_palette <- scales::hue_pal()(length(unique(te_aff_t$predicted_ancestry_thres[!is.na(te_aff_t$predicted_ancestry_thres)])))
#
#  p_lm_ancestry_pred <- tryCatch({
#    write_output(
#      quote(plot_count_lm(te_aff_t, min_samples=5, residuals=FALSE, group="predicted_ancestry_thres",
#                          log_scale=TRUE, covariates = covar_med, x_lab="Predicted Ancestry",
#                          y_lab="Total TE count", type=NA, chr=NA, fill_palette=ancestry_palette)),
#      "Linear model by predicted ancestry thres (te_aff_t, all types, log scale)"
#    )
#  }, error = function(e) {
#    cat("Skipping ancestry LM plot - insufficient data or error:", e$message, "\n")
#    NULL
#  })
#
#  if (!is.null(p_lm_ancestry_pred)) {
#    titled_print(p_lm_ancestry_pred, "Linear model by predicted ancestry thres (te_aff_t, all types, log scale)")
#    ggsave(paste0(plot_dir, "te_count_lm_predicted_ancestry_thres_aff_all_t.png"), plot=p_lm_ancestry_pred, width = 9, height = 5)
#  }
#} else {
#  cat("Skipping predicted_ancestry_thres LM plot - insufficient samples with ancestry data\n")
#  cat("Samples with ancestry:", sum(!is.na(te_aff_t$predicted_ancestry_thres)), "\n")
#}
#
## Plot count LM grouped by mapped_label
#if ("mapped_label" %in% colnames(te_aff_t) &&
#    sum(!is.na(te_aff_t$mapped_label)) >= 10) {
#
#  mapped_label_palette <- scales::hue_pal()(length(unique(te_aff_t$mapped_label[!is.na(te_aff_t$mapped_label)])))
#
#  p_lm_ancestry_mapped <- tryCatch({
#    write_output(
#      quote(plot_count_lm(te_aff_t, min_samples=5, residuals=FALSE, group="mapped_label",
#                          log_scale=TRUE, covariates = covar_med, x_lab="Mapped Label",
#                          y_lab="Total TE count", type=NA, chr=NA, fill_palette=mapped_label_palette)),
#      "Linear model by mapped label (te_aff_t, all types, log scale)"
#    )
#  }, error = function(e) {
#    cat("Skipping mapped_label LM plot - insufficient data or error:", e$message, "\n")
#    NULL
#  })
#
#  if (!is.null(p_lm_ancestry_mapped)) {
#    titled_print(p_lm_ancestry_mapped, "Linear model by mapped label (te_aff_t, all types, log scale)")
#    ggsave(paste0(plot_dir, "te_count_lm_mapped_label_aff_all_t.png"), plot=p_lm_ancestry_mapped, width = 9, height = 5)
#  }
#} else {
#  cat("Skipping mapped_label LM plot - insufficient samples with ancestry data\n")
#  cat("Samples with mapped_label:", sum(!is.na(te_aff_t$mapped_label)), "\n")
#}
#
## Plot pie chart for ancestry
#cat("Plotting: clinical_ancestry_tumour.pdf\n")
#p_ancestry_pie_dataset <- plot_ancestry_pie(te_all_t, output_dir = NULL, plot_prefix = "tumour")
#print(p_ancestry_pie_dataset)
#ggsave(paste0(plot_dir, "dataset/clinical_ancestry_tumour.pdf"), width = 6, height = 5)
#
##### PCA/UMAP BY OTHER VARIABLES ####
#cat("\n========================================\n")
#cat("PCA/UMAP COLORED BY OTHER VARIABLES\n")
#cat("========================================\n\n")
#
## Variables to color by (includes both ancestry variables and clinical variables)
#color_variables <- c("mapped_label", "tumor_type", "cohort", "TP53_status", "sex")
#
## Plot PCA by each variable
#for (var in color_variables) {
#  if (var %in% colnames(location_ancestry)) {
#    cat("Plotting PCA colored by", var, "...\n")
#    p_pca <- write_output(
#      quote(plot_pca_by_variable(pca_results, color_by = var,
#                                  output_dir = plot_dir, plot_prefix = "tumour")),
#      paste0("PCA colored by ", var)
#    )
#    if (!is.null(p_pca)) {
#      titled_print(p_pca, paste0("PCA of TE Location Windows by ", var, " (Tumour)"))
#    }
#  } else {
#    cat("Skipping PCA for", var, "- column not found in data\n")
#  }
#}
#
## Plot UMAP by each variable
#for (var in color_variables) {
#  if (var %in% colnames(location_ancestry)) {
#    cat("Plotting UMAP colored by", var, "...\n")
#    p_umap <- write_output(
#      quote(plot_umap_by_variable(umap_results, color_by = var,
#                                   output_dir = plot_dir, plot_prefix = "tumour")),
#      paste0("UMAP colored by ", var)
#    )
#    if (!is.null(p_umap)) {
#      titled_print(p_umap, paste0("UMAP of TE Location Windows by ", var, " (Tumour)"))
#    }
#  } else {
#    cat("Skipping UMAP for", var, "- column not found in data\n")
#  }
#}



cat("✓ Script completed successfully\n")
