#!/usr/bin/env Rscript

# Tumour TE Visualization - Taylor Cohort
# Taylor cohort-specific analyses

# Source common setup and load data
source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/00_viz_common_setup.R")
source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/00_viz_load_data_tumour.R")

cat("Running 02_te_viz_tumour_12_taylor.R...\n")

#### TAYLOR COHORT ANALYSIS ####
# Check if Taylor cohort data exists (Taylor is a germline cohort, not available in tumour data)
if (!exists("te_taylor_t")) {
  cat("⚠ NOTE: Taylor cohort is a germline cohort and not available in tumour analysis.\n")
  cat("  Skipping tumour Taylor cohort analysis.\n")
  cat("  For Taylor cohort analysis, see germline script: 02_te_viz_germline_10_taylor.R\n\n")
  cat("✓ Script completed successfully (skipped)\n")
} else {
  # This code will only run if Taylor tumour data exists
  write_output(quote({
    cat("Taylor Cohort Summary:\n")
    cat("Total samples:", nrow(te_taylor_t), "\n")
    cat("\nTumor type distribution:\n")
    print(table(te_taylor_t$tumor_type))
    cat("\nTumor type subclass distribution:\n")
    print(table(te_taylor_t$tumor_type_subclass))
  }), "Taylor Cohort Summary Statistics")

  write_output(quote(plot_count_kruskal(df = te_taylor_t, chr = NA, type = "total", group = "tumor_type_subclass", x_lab = "Tumor Type Subclass", y_lab = "Total TE Count", log_scale = FALSE)),
               "Taylor: TE Count by Tumor Type Subclass (Kruskal-Wallis)")
  p_kruskal <- plot_count_kruskal(df = te_taylor_t, chr = NA, type = "total", group = "tumor_type_subclass", x_lab = "Tumor Type Subclass", y_lab = "Total TE Count", log_scale = FALSE)
  titled_print(p_kruskal, "Taylor: TE Count by Tumor Type Subclass (Kruskal-Wallis)")
  ggsave(paste0(plot_dir, "other/taylor_te_count_by_subclass_kruskal.png"), plot = p_kruskal, width = 10, height = 6)

  write_output(quote(plot_count_age(df = te_taylor_t, type = "total", chr = NA, y_lab = "Total TE Count")),
               "Taylor: TE Count by Age at Diagnosis")
  p_age <- plot_count_age(df = te_taylor_t, type = "total", chr = NA, y_lab = "Total TE Count")
  titled_print(p_age, "Taylor: TE Count by Age at Diagnosis")
  ggsave(paste0(plot_dir, "other/taylor_te_count_by_age.png"), plot = p_age, width = 8, height = 6)

  cat("\n===== SCRIPT COMPLETED SUCCESSFULLY =====\n")
  cat("All analysis sections have been processed.\n")
  cat("Generated plots saved to:", plot_dir, "\n")
}
