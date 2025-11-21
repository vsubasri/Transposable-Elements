#!/usr/bin/env Rscript

# Tumour TE Visualization - LOH Timing
# TP53 LOH timing analysis

# Source common setup and load data
source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/00_viz_common_setup.R")
source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/00_viz_load_data_tumour.R")

cat("Running 02_te_viz_tumour_04_loh.R...\n")

cat("\n===== TP53 LOH TIMING ANALYSIS =====\n")
tryCatch({
  loh_time_mod <- loh_time %>%
    select(sample, time) %>%
    mutate(sample = ifelse(grepl("^KiCS", sample),
                           paste0(sub("^KiCS", "", sample), "_T"),
                           sample))
  timing <- merge(te_aff_t, loh_time_mod, by="sample")
  write_output(quote(head(timing)), "Head of merged timing data (TP53 LOH)")
  write_output(quote(wilcox.test(time ~ TP53_status, data = timing)), "Wilcoxon test for LOH timing by TP53 status")
  write.table(timing, file = "/Users/briannelaverty/Documents/R_Malkin/te/data/final/te_loh.csv", quote = FALSE, sep = ',', row.names = FALSE, col.names = TRUE)
}, error = function(e) {
  cat("Warning: Could not perform TP53 LOH timing analysis:", e$message, "\n")
})



cat("✓ Script completed successfully\n")
