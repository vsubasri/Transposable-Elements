#!/usr/bin/env Rscript

# Tumour TE Visualization - Descriptive Statistics
# Dataset description, basic statistics, full-length LINE1, ancestry distributions

# Source common setup and load data
source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/00_viz_common_setup.R")
REQUIRED_DATA <- c("count_matrix", "expand", "clinical")
source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/00_viz_load_data_tumour.R")

# Initialize module-specific text output
init_module_sink(paste0(plot_dir, "dataset/"), "DESCRIPTIVE")

cat("Running 02_te_viz_tumour_01_descriptive.R...\n")

# Define x-axis label for TP53 status (used in multiple plots)
x_tp53 <- expression("Somatic " * italic("TP53") * " status")

#### GENERAL STATS ####
cat("\n===== GENERAL STATS: Generating general statistics and plots... =====\n")

cat("Calculating total TE counts for all samples...\n")
# Identify samples in the top 99% quantile of total TE count (but don't print)
quantile_99 <- quantile(te_all_t$total, 0.99, na.rm = TRUE)
top_99_samples <- te_all_t$sample[te_all_t$total >= quantile_99]

# Remove outlier samples
to_remove <- te_all_t$sample[te_all_t$total > 100000]

# Cohort Summary Table
write_output(quote({
  cohort_table <- create_cohort_summary(te_all_t, hostseq_cancer = NULL)
  cohort_table
}), "Cohort Summary Statistics")

# Count TE occurrences (doesn't include no hit samples for summary)
write_output(quote(count_TE_occurrences(te_expand=te_all_expand_t, te_all=te_all_t, te_all_all=te_all_all_t, nohits_file="nohits_final_te_count_t_te_all_all_t", type="T")), "Count TE Occurrences for All Samples")

write_output(quote(count_TE_occurrences(te_expand=te_aff_expand_t %>% filter(sample %in% te_aff_t$sample), te_all=te_aff_t, te_all_all=te_aff_t, nohits_file="nohits_final_te_count_t_te_aff_selected_t", type="T")), "Count TE Occurrences for Affected KICS Samples")

write_output(quote(count_TE_occurrences(te_expand=te_kics_expand_t %>% filter(sample %in% te_kics_t$sample), te_all=te_kics_t, te_all_all=te_kics_t, nohits_file="nohits_final_te_count_t_te_kics_selected_t", type="T")), "Count TE Occurrences for KICS Samples")

cat("\n===== COMPARE COUNT TO ADULT =====\n")
te_gt1 <- nrow(te_aff_t[te_aff_t$total>0,])
te_0 <- nrow(te_aff_t[te_aff_t$total==0,])
write_output(quote({
  cat("Tumour samples with >1 TE:", te_gt1, "\n")
  cat("Tumour samples with 0 TE:", te_0, "\n")
}), "Tumour sample TE count summary")

# Calculate KICS cohort counts
kics_with_te <- nrow(te_kics_t[te_kics_t$total > 0,])
kics_without_te <- nrow(te_kics_t[te_kics_t$total == 0,])

write_output(quote({
  cat("KICS samples with >0 TE:", kics_with_te, "\n")
  cat("KICS samples with 0 TE:", kics_without_te, "\n")
}), "KICS sample TE count summary")

# compared to greenbaum paper which used totalrecall and xtea for LINE1
write_output(quote(fisher_test_and_plot(group1_yes = 2176 , group1_no = 2494, group2_yes = kics_with_te, group2_no = kics_without_te)),
             "Fisher test: Tumour vs Adult (at least 1 TE)")
titled_print(fisher_test_and_plot(group1_yes = 2176, group1_no = 2494, group2_yes = kics_with_te, group2_no = kics_without_te),
              "Fisher test: Tumour vs Adult (at least 1 TE)")
ggsave(paste0(plot_dir, "general/te_tumour_atleast1.png"), width = 9, height = 5)

# Sample statistics: percent with at least 1 TE
write_output(quote({
  kics_with_te_pct <- nrow(te_kics_t[te_kics_t$total > 0,])/nrow(te_kics_t) * 100
  lfs_with_te_pct <- nrow(te_lfs_t[te_lfs_t$total > 0,])/nrow(te_lfs_t) * 100
  cat("Percent of KICS tumour samples with at least 1 TE:", round(kics_with_te_pct, 2), "%\n")
  cat("Percent of LFS tumour samples with at least 1 TE:", round(lfs_with_te_pct, 2), "%\n")
}), "Percent of Samples with at Least 1 TE")

# Sample statistics: percent with at least 1 LINE1
write_output(quote({
  kics_with_l1_pct <- nrow(te_kics_t[te_kics_t$LINE1 > 0,])/nrow(te_kics_t) * 100
  lfs_with_l1_pct <- nrow(te_lfs_t[te_lfs_t$LINE1 > 0,])/nrow(te_lfs_t) * 100
  cat("Percent of KICS tumour samples with at least 1 LINE1:", round(kics_with_l1_pct, 2), "%\n")
  cat("Percent of LFS tumour samples with at least 1 LINE1:", round(lfs_with_l1_pct, 2), "%\n")
}), "Percent of Samples with at Least 1 LINE1")

# Summary statistics for tumour samples
write_output(quote({
  cat("Summary statistics for KICS tumour samples:\n")
  print(summary(te_kics_t$total))
  cat("Summary statistics for LFS tumour samples:\n")
  print(summary(te_lfs_t$total))
}), "Summary Statistics for Tumour Samples")


# TE count summary plot - all common (only if COMMON_MODE is TRUE)
if (COMMON_MODE && exists("te_all_common_t")) {
  write_output(quote(plot_te_counts_summary(te_all_common_t, y_lab="Repeat count", log_scale=TRUE, breaks=c(10,100,1000))), "TE Count Summary Plot - All Common")
  titled_print(plot_te_counts_summary(te_all_common_t, y_lab="Repeat count", log_scale=TRUE, breaks=c(10,100,1000)), "TE Count Summary Plot - All Common")
  ggsave(paste0(plot_dir, "general/te_count_type_all_common.png"), width = 9, height = 5)
} else {
  cat("Skipping common TE plot (COMMON_MODE = FALSE)\n")
}

# TE count summary plot - all rare
write_output(quote(plot_te_counts_summary(te_all_t, y_lab="Repeat count", log_scale=TRUE, breaks=c(10,100,500))), "TE Count Summary Plot - All Rare")
titled_print(plot_te_counts_summary(te_all_t, y_lab="Repeat count", log_scale=TRUE, breaks=c(10,100,500)), "TE Count Summary Plot - All Rare")
ggsave(paste0(plot_dir, "general/te_count_type_all_rare.png"), width = 9, height = 5)

# Plots for TEs shared by samples
write_output(quote(plot_te_counts_unique(te_aff_expand_t)), "TEs Shared by Samples")
titled_print(plot_te_counts_unique(te_aff_expand_t), "TEs Shared by Samples")

# Stacked bar plot for common TEs (only if COMMON_MODE is TRUE)
if (COMMON_MODE && exists("te_aff_expand_common_t")) {
  write_output(quote(stacked_bar_plot_num_samples(te_aff_expand_common_t, c(1, 5, 20))), "Stacked Bar Plot: # TEs in Range of Samples")
  titled_print(stacked_bar_plot_num_samples(te_aff_expand_common_t, c(1, 5, 20)), "Stacked Bar Plot: # TEs in Range of Samples")
  ggsave(paste0(plot_dir, "general/te_type_sample_common.png"), width = 9, height = 5)
} else {
  cat("Skipping common TE stacked bar plot (COMMON_MODE = FALSE)\n")
}


# Total TE calls plot
cat("Creating total TE calls plot...\n")
write_output(quote(plot_te_sum(te_kics_t)), "Total TE Calls Plot (KICS)")
titled_print(plot_te_sum(te_kics_t), "Total TE Calls Plot (KICS)")
ggsave(paste0(plot_dir, "general/te_total_calls_kics.png"), width = 9, height = 5)

cat("General statistics and plots for tumour data completed!\n")


# NOTE: Full-length LINE1 analysis moved to 02_te_viz_tumour_14_fulllength_young_source.R


cat("\n===== DESCRIBE DATASET =====\n")
tryCatch({
  clin <- te_all_t %>% dplyr::select(sample, tumor_type, sex, age_at_diagnosis, TP53_status)
  clin <- clin %>% distinct()
  titled_print(plot_pie_chart(clin, "tumor_type", "Tumour type"), "Tumour type distribution")
  ggsave(paste0(plot_dir, "dataset/clinical_germline_tt.png"), width = 6, height = 5)
  titled_print(plot_pie_chart(clin, "TP53_status", x_tp53), "TP53 status distribution")
  ggsave(paste0(plot_dir, "dataset/clinical_germline_tp53.png"), width = 6, height = 5)
  titled_print(plot_pie_chart(clin, "sex", "Sex"), "Sex distribution")
  ggsave(paste0(plot_dir, "dataset/clinical_sex.png"), width = 6, height = 5)
  titled_print(plot_histogram_age(clin, bin=2), "Age at diagnosis distribution")
  ggsave(paste0(plot_dir, "dataset/clinical_germline_age.png"), width = 9, height = 5)
}, error = function(e) {
  cat("Warning: Could not perform dataset description analysis:", e$message, "\n")
})


cat("\n===== RNA GENE EXPRESSION ANALYSIS =====\n")



cat("✓ Script completed successfully\n")

# Close module-specific sink
close_module_sink()
