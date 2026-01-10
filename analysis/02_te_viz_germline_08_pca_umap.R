#!/usr/bin/env Rscript

# Germline TE Visualization - PCA/UMAP
# PCA and UMAP analyses by various variables

# Source common setup and load data
source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/00_viz_common_setup.R")
REQUIRED_DATA <- c("ancestry", "location", "clinical")
source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/00_viz_load_data_germline.R")

# Initialize module-specific text output
init_module_sink(paste0(plot_dir, "pca_umap/"), "PCA_UMAP")

cat("Running 02_te_viz_germline_08_pca_umap.R...\n")

#### PCA/UMAP BY OTHER VARIABLES ####
cat("\n========================================\n")
cat("PCA/UMAP COLORED BY OTHER VARIABLES\n")
cat("========================================\n\n")

# Merge location windows with ancestry
location_ancestry <- write_output(
  quote(merge_location_ancestry(location_100kb_g, ancestry)),
  "Merging location windows with ancestry"
)

# Identify clinical columns to exclude from PCA/UMAP
clinical_cols <- c("sample", "predicted_ancestry_thres", "mapped_label", "age", "sex", "cohort", "tumor_type",
                   "TP53_germline", "affected", "cancer", "age_diagnosis", "TP53_status")
exclude_cols <- intersect(clinical_cols, colnames(location_ancestry))

# Perform PCA and UMAP on location windows
pca_results <- write_output(
  quote(perform_location_pca(location_ancestry, exclude_cols = exclude_cols)),
  "Performing PCA on location windows"
)

umap_results <- write_output(
  quote(perform_location_umap(location_ancestry, exclude_cols = exclude_cols)),
  "Performing UMAP on location windows"
)

# Variables to color by (includes both ancestry variables and clinical variables)
color_variables <- c("mapped_label", "tumor_type", "cohort", "TP53_status", "sex")

# Plot PCA by each variable
for (var in color_variables) {
  if (var %in% colnames(location_ancestry)) {
    cat("Plotting PCA colored by", var, "...\n")
    p_pca <- write_output(
      quote(plot_pca_by_variable(pca_results, color_by = var,
                                  output_dir = plot_dir, plot_prefix = "germline")),
      paste0("PCA colored by ", var)
    )
    if (!is.null(p_pca)) {
      titled_print(p_pca, paste0("PCA of TE Location Windows by ", var, " (Germline)"))
    }
  } else {
    cat("Skipping PCA by", var, "(column not found)\n")
  }
}

# Plot UMAP by each variable
for (var in color_variables) {
  if (var %in% colnames(location_ancestry)) {
    cat("Plotting UMAP colored by", var, "...\n")
    p_umap <- write_output(
      quote(plot_umap_by_variable(umap_results, color_by = var,
                                   output_dir = plot_dir, plot_prefix = "germline")),
      paste0("UMAP colored by ", var)
    )
    if (!is.null(p_umap)) {
      titled_print(p_umap, paste0("UMAP of TE Location Windows by ", var, " (Germline)"))
    }
  } else {
    cat("Skipping UMAP by", var, "(column not found)\n")
  }
}

# UMAP for HostSeq cohort only
if ("cohort" %in% colnames(location_ancestry)) {
  cat("\n===== UMAP ANALYSIS - HOSTSEQ COHORT ONLY =====\n")
  location_ancestry_hostseq <- location_ancestry %>% filter(cohort == "HostSeq")

  if (nrow(location_ancestry_hostseq) > 0) {
    # Perform UMAP on HostSeq samples only
    umap_results_hostseq <- write_output(
      quote(perform_location_umap(location_ancestry_hostseq, exclude_cols = exclude_cols)),
      "Performing UMAP on location windows (HostSeq only)"
    )

    # Plot UMAP by each variable for HostSeq cohort
    for (var in color_variables) {
      if (var %in% colnames(location_ancestry_hostseq)) {
        cat("Plotting HostSeq UMAP colored by", var, "...\n")
        p_umap_hostseq <- write_output(
          quote(plot_umap_by_variable(umap_results_hostseq, color_by = var,
                                       output_dir = plot_dir, plot_prefix = "germline_hostseq")),
          paste0("HostSeq UMAP colored by ", var)
        )
        if (!is.null(p_umap_hostseq)) {
          titled_print(p_umap_hostseq, paste0("UMAP of TE Location Windows by ", var, " (HostSeq Only)"))
        }
      }
    }
  } else {
    cat("No HostSeq samples found - skipping HostSeq-only UMAP\n")
  }
}



cat("✓ Script completed successfully\n")

# Close module-specific sink
close_module_sink()
