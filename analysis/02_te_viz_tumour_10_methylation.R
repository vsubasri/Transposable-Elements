#!/usr/bin/env Rscript

# Tumour TE Visualization - Methylation
# TE-methylation probe overlap analysis

# Source common setup and load data
source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/00_viz_common_setup.R")
REQUIRED_DATA <- c("expand", "clinical")
source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/00_viz_load_data_tumour.R")

# Initialize module-specific text output
init_module_sink(paste0(plot_dir, "methylation/"), "METHYLATION")

cat("Running 02_te_viz_tumour_10_methylation.R...\n")

#### TE-METHYLATION PROBE OVERLAP ANALYSIS ####
cat("\n===== TE-METHYLATION PROBE OVERLAP ANALYSIS =====\n")
tryCatch({
  # Find TEs overlapping with methylation probe regions
  te_probe_overlaps_t <- find_te_probe_overlaps(
    te_expand = te_aff_expand_t,
    probe_file = "/Users/briannelaverty/Documents/R_Malkin/te/data/lfs_probe_info.csv",
    min_overlap_pct = 0,
    output_file = paste0(plot_dir, "methylation/methylation_probe_overlaps_tumour.csv")
  )

  if (!is.null(te_probe_overlaps_t)) {
    cat("\nTop probe regions with TE overlaps:\n")
    print(head(te_probe_overlaps_t %>% filter(n_overlaps > 0) %>% arrange(desc(n_overlaps)), 10))
  }
}, error = function(e) {
  cat("Warning: Could not perform TE-probe overlap analysis:", e$message, "\n")
})

cat("\n===== STRUCTURAL VARIANT - TE OVERLAP ANALYSIS =====\n")


cat("✓ Script completed successfully\n")

# Close module-specific sink
close_module_sink()
