#!/usr/bin/env Rscript

# Germline TE Visualization - Descriptive Statistics
# Dataset description, basic statistics, TE source, ancestry distributions

# Source common setup and load data
source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/00_viz_common_setup.R")
REQUIRED_DATA <- c("count_matrix", "expand", "clinical", "genes")
source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/00_viz_load_data_germline.R")

# Initialize module-specific text output
init_module_sink(paste0(plot_dir, "dataset/"), "DESCRIPTIVE")

cat("Running 02_te_viz_germline_01_descriptive.R...\n")

#### DESCRIBE DATASET ####
clin <- te_all %>% dplyr::select(sample, tumor_type, sex, age_at_diagnosis, TP53_status, cohort)
clin <- clin %>% distinct()

# Plot pie chart for tumor_type
cat("Plotting: clinical_germline_tt.pdf\n")
plot_pie_chart(clin, "tumor_type", "Tumour type")
ggsave(paste0(plot_dir, "dataset/clinical_germline_tt.pdf"), width = 6, height = 5)

# Plot pie chart for TP53_status
cat("Plotting: clinical_germline_tp53.pdf\n")
tp53_title <- expression("Germline " * italic("TP53") * " status")
plot_pie_chart(clin, "TP53_status", tp53_title)
ggsave(paste0(plot_dir, "dataset/clinical_germline_tp53.pdf"), width = 6, height = 5)

# Plot pie chart for sex
cat("Plotting: clinical_sex.pdf\n")
plot_pie_chart(clin, "sex", "Sex")
ggsave(paste0(plot_dir, "dataset/clinical_sex.pdf"), width = 6, height = 5)

# Plot pie chart for ancestry
cat("Plotting: clinical_ancestry.pdf\n")
p_ancestry_pie_dataset <- plot_ancestry_pie(te_all, output_dir = NULL, plot_prefix = "germline")
print(p_ancestry_pie_dataset)
ggsave(paste0(plot_dir, "dataset/clinical_ancestry.pdf"), width = 6, height = 5)

# NOTE: HostSeq ancestry pie charts (All/Filter/Analysis groups) are now generated
# in 01_te_processing_germline.R BEFORE the filter group is removed from the data.
# This ensures both filter (66%) and analysis (34%) groups are properly represented.

# Plot histogram for age_at_diagnosis
cat("Plotting: clinical_germline_age.pdf\n")
plot_histogram_age(clin, bin=2)
ggsave(paste0(plot_dir, "dataset/clinical_germline_age.pdf"), width = 9, height = 5)

# cohort - use te_all samples only
clin_te_all <- clin %>%
  filter(sample %in% te_all$sample) %>%
  mutate(new_cohort = case_when(
    cohort == "LFS_mut" ~ "LFS",
    cohort == "nick" ~ "LFS",
    TRUE ~ cohort
  ))
cat("Plotting: clinical_dataset.pdf (using te_all samples)\n")
plot_pie_chart(clin_te_all, "new_cohort", "Dataset")
ggsave(paste0(plot_dir, "dataset/clinical_dataset.pdf"), width = 9, height = 5)


#### GENERAL STATS ####

# Cohort Summary Table
write_output(quote({
  cohort_table <- create_cohort_summary(te_all, hostseq_cancer = hostseq_cancer)
  cat("Cohort Summary (Sample Counts):\n\n")
  print(cohort_table, row.names = FALSE)
}), "Cohort Summary Statistics")

# Dataset Summary: te_lfs and te_lfs_mut_wt
write_output(quote({
  cat("\n=== te_lfs Dataset Summary ===\n")
  cat("Definition: TP53_status == 'Mutant'\n")
  cat("Total samples:", length(unique(te_lfs$sample)), "\n")
  lfs_cohorts <- te_lfs %>%
    group_by(cohort, TP53_status) %>%
    summarise(n_samples = n_distinct(sample), .groups = "drop") %>%
    arrange(cohort, TP53_status)
  print(as.data.frame(lfs_cohorts), row.names = FALSE)

  cat("\n=== te_lfs_mut_wt Dataset Summary ===\n")
  cat("Definition: TP53_status == 'Mutant' OR cohort == 'LFS_wt'\n")
  cat("Total samples:", length(unique(te_lfs_mut_wt$sample)), "\n")
  lfs_mut_wt_cohorts <- te_lfs_mut_wt %>%
    group_by(cohort, TP53_status) %>%
    summarise(n_samples = n_distinct(sample), .groups = "drop") %>%
    arrange(cohort, TP53_status)
  print(as.data.frame(lfs_mut_wt_cohorts), row.names = FALSE)
}), "LFS Dataset Summaries (te_lfs and te_lfs_mut_wt)")

# TE count summary plots
if (PROCESS_COMMON_TES) {
  write_output(quote(plot_te_counts_summary(te_all_common, y="Repeat frequency", log_scale=TRUE, breaks=c(10,100,1000))), "TE Count Summary - All Common")
  p1 <- plot_te_counts_summary(te_all_common, y="Repeat frequency", log_scale=TRUE, breaks=c(10,100,1000))
  titled_print(p1, "TE count summary plot - all common")
  ggsave(paste0(plot_dir, "counts_clinical_kics/te_count_type_all_common.pdf"), plot = p1, width = 9, height = 5)
}

write_output(quote(plot_te_counts_summary(te_all, y="Rare repeat frequency", log_scale=TRUE, breaks=c(10,100,500))), "TE Count Summary - All Rare")
p2 <- plot_te_counts_summary(te_all, y="Rare repeat frequency", log_scale=TRUE, breaks=c(10,100,500))
titled_print(p2, "TE count summary plot - all rare")
ggsave(paste0(plot_dir, "counts_clinical_kics/te_count_type_all_rare.pdf"), plot = p2, width = 9, height = 5)

write_output(quote(plot_te_counts_summary(te_aff, y="Rare repeat frequency", log_scale=TRUE, breaks=c(10,100,500))), "TE Count Summary - Affected Rare")
p3 <- plot_te_counts_summary(te_aff, y="Rare repeat frequency", log_scale=TRUE, breaks=c(10,100,500))
titled_print(p3, "TE count summary plot - affected rare")
ggsave(paste0(plot_dir, "counts_clinical_kics/te_count_type_aff_rare.pdf"), plot = p3, width = 9, height = 5)

# TEs shared by samples
write_output(quote(plot_te_counts_unique(te_aff_expand)), "TEs Shared by Samples")
p4 <- plot_te_counts_unique(te_aff_expand)
titled_print(p4, "TE counts unique (te_aff_expand)")
ggsave(paste0(plot_dir, "counts_clinical_kics/te_counts_unique_aff_expand.pdf"), plot = p4, width = 9, height = 5)

if (PROCESS_COMMON_TES) {
  write_output(quote(stacked_bar_plot_num_samples(te_aff_expand_common, c(1, 5, 20))), "Stacked Bar Plot: # TEs in Range of Samples")
  p5 <- stacked_bar_plot_num_samples(te_aff_expand_common, c(1, 5, 20))
  titled_print(p5, "Stacked bar plot num samples (te_aff_expand_common)")
  ggsave(paste0(plot_dir, "counts_clinical_kics/te_type_sample_common.pdf"), plot = p5, width = 9, height = 5)
}

# Total TE calls
write_output(quote(plot_te_sum(te_kics)), "Total TE Calls - KICS")
p6 <- plot_te_sum(te_kics)
titled_print(p6, "Total TE calls plot (te_kics)")
ggsave(paste0(plot_dir, "counts_clinical_kics/te_count_total_te_calls_kics.pdf"), plot = p6, width = 9, height = 5)

# Median number of TEs per sample
write_output(quote({
  cat("Summary statistics for total TEs per sample:\n\n")
  cat("KICS:\n")
  print(summary(te_kics$total))
  cat("\nLFS:\n")
  print(summary(te_lfs$total))
}), "Summary Statistics: Total TEs per Sample")

#### TE SOURCE ####
cat("Extracting TE source information...\n")
te_aff_expand <- extract_info_fields(te_aff_expand)
te_lfs_expand <- extract_info_fields(te_lfs_expand)

write_output(quote(nrow(table(te_aff_expand$source))), "Number of unique sources in te_aff_expand")
write_output(quote(table(te_aff_expand$source)[table(te_aff_expand$source) > 5]), "Source counts > 5 in te_aff_expand")
write_output(quote(table(te_lfs_expand$source)[table(te_lfs_expand$source) > 1]), "Source counts > 1 in te_lfs_expand")

cat("✓ Descriptive statistics completed successfully\n")

# Close module-specific sink
close_module_sink()
