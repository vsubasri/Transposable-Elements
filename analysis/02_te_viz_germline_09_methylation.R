#!/usr/bin/env Rscript

# Germline TE Visualization - Methylation
# TE-methylation probe overlap analysis

# Source common setup and load data
source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/00_viz_common_setup.R")
source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/00_viz_load_data_germline.R")

cat("Running 02_te_viz_germline_09_methylation.R...\n")

#### TE-METHYLATION PROBE OVERLAP ANALYSIS ####
cat("\n===== TE-METHYLATION PROBE OVERLAP ANALYSIS =====\n")

## Affected only (te_aff_expand)
#cat("\nAnalyzing affected samples (te_aff_expand)...\n")
#tryCatch({
#  te_probe_overlaps_aff <- find_te_probe_overlaps(
#    te_expand = te_aff_expand,
#    probe_file = "/Users/briannelaverty/Documents/R_Malkin/te/data/lfs_probe_info.csv",
#    min_overlap_pct = 0,
#    output_file = paste0(r_dir_files, "te_methylation_probe_overlaps_aff.csv")
#  )
#
#  if (!is.null(te_probe_overlaps_aff)) {
#    cat("\nTop probe regions with TE overlaps (affected only):\n")
#    print(head(te_probe_overlaps_aff %>% filter(n_overlaps > 0) %>% arrange(desc(n_overlaps)), 10))
#  }
#}, error = function(e) {
#  cat("Warning: Could not perform TE-probe overlap analysis (affected):", e$message, "\n")
#})

# Affected + Unaffected, excluding Taylor and HostSeq (te_aff_unaff_expand)
cat("\nAnalyzing affected + unaffected samples, excluding Taylor/HostSeq (te_aff_unaff_expand)...\n")
tryCatch({
  te_probe_overlaps_aff_unaff <- find_te_probe_overlaps(
    te_expand = te_aff_unaff_expand,
    probe_file = "/Users/briannelaverty/Documents/R_Malkin/te/data/aff_unaff_probe_info.csv",
    min_overlap_pct = 0,
    output_file = paste0(r_dir_files, "te_methylation_probe_overlaps_aff_unaff.csv")
  )

  if (!is.null(te_probe_overlaps_aff_unaff)) {
    cat("\nTop probe regions with TE overlaps (affected + unaffected):\n")
    print(head(te_probe_overlaps_aff_unaff %>% filter(n_overlaps > 0) %>% arrange(desc(n_overlaps)), 10))
  }
}, error = function(e) {
  cat("Warning: Could not perform TE-probe overlap analysis (aff + unaff):", e$message, "\n")
})

# lfs
#cat(paste0("Plotting: te_count_wilcox_lfs.pdf\n"))
#plots_wilcox_lfs <- vector("list", length(types))
#for (i in seq_along(types)) {
#  cat(paste0("\n===== Plotting: te_count_wilcox_lfs_", ifelse(is.na(types[i]), "all", types[i]), ".pdf =====\n"))
#  plots_wilcox_lfs[[i]] <- write_output(
#    quote(plot_count_wilcox(te_lfs, chr=NA, type=types[i], group="Cancer", x_lab="Cancer status", y_lab="Repeat frequency")),
#    paste0("Wilcoxon plot output (LFS, type=", ifelse(is.na(types[i]), "all", types[i]), ")")
#  )
#  titled_print(plots_wilcox_lfs[[i]], paste0("Wilcoxon plot (LFS, type=", ifelse(is.na(types[i]), "all", types[i]), ")"))
#  ggsave(paste0(plot_dir, "te_count_wilcox_lfs_", ifelse(is.na(types[i]), "all", types[i]), ".pdf"), plot=plots_wilcox_lfs[[i]], width = 3, height = 5)
#}

# sequencing cohort
#plots_kruskal_lfs_seq <- vector("list", length(types))
#for (i in seq_along(types)) {
#  cat(paste0("\n===== Plotting: te_count_kruskal_lfs_seq_", ifelse(is.na(types[i]), "all", types[i]), ".png =====\n"))
#  plots_kruskal_lfs_seq[[i]] <- write_output(
#    quote(plot_count_kruskal(te_lfs_sequencing, type=types[i], log_scale=TRUE, group="lfs_center", x_lab="Cluster", y_lab="Repeat frequency", chr=NA)),
#    paste0("Kruskal plot output (LFS sequencing, type=", ifelse(is.na(types[i]), "all", types[i]), ")")
#  )
#  titled_print(plots_kruskal_lfs_seq[[i]], paste0("Kruskal plot (LFS sequencing, type=", ifelse(is.na(types[i]), "all", types[i]), ")"))
#  ggsave(paste0(plot_dir, "te_count_kruskal_lfs_seq_", ifelse(is.na(types[i]), "all", types[i]), ".png"), plot=plots_kruskal_lfs_seq[[i]], width = 3, height = 5)
#}




cat("✓ Script completed successfully\n")
