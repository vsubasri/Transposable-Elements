#!/usr/bin/env Rscript

# Tumour TE Visualization - Clinical Associations
# TE counts by clinical variables, P53 mutation, ancestry comparisons

# Source common setup and load data
source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/00_viz_common_setup.R")
source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/00_viz_load_data_tumour.R")

cat("Running 02_te_viz_tumour_02_clinical.R...\n")

#### TE COUNT ####
# Build types list based on configuration
types <- c(NA, "LINE1")
if (INCLUDE_ALU) types <- c(types, "ALU")
if (INCLUDE_SVA) types <- c(types, "SVA")
cat("TE types for analysis:", paste(types, collapse=", "), "\n")

x_tp53 <- expression("Germline " * italic("TP53") * " status")

# affected
cat("\n===== TE COUNT of affected by TP53_status - wilcox test =====\n")
plots_wilcox_aff <- vector("list", length(types))
for (i in seq_along(types)) {
  plots_wilcox_aff[[i]] <- plot_count_wilcox(te_aff_t, type=types[i], group="TP53_status", y_lab="Repeat count", chr=NA, x_lab=x_tp53, log_scale=FALSE)
  write_output(quote(plot_count_wilcox(te_aff_t, type=types[i], group="TP53_status", y_lab="Repeat count", chr=NA, x_lab=x_tp53, log_scale=FALSE)),
               paste0("Wilcoxon plot (affected, type=", ifelse(is.na(types[i]), "all", types[i]), ")"))
  titled_print(plots_wilcox_aff[[i]], paste0("Wilcoxon plot (affected, type=", ifelse(is.na(types[i]), "all", types[i]), ")"))
  ggsave(paste0(plot_dir, "counts_cohort/te_count_wilcox_aff_", ifelse(is.na(types[i]), "all", types[i]), ".png"), plot=plots_wilcox_aff[[i]], width = 3, height = 5)
}

cat("\n===== TE count of affected by TP53 status - linear model =====\n")
plots_lm_aff <- lapply(seq_along(types), function(i) {
  y_label <- ifelse(is.na(types[i]), "Total TE count", paste0(types[i], " count"))
  p <- write_output(
    quote(plot_count_lm(te_aff_t, type=types[i], group="TP53_status", y_lab=y_label, breaks=c(10,100,1000,10000), covariates=covar_med, residuals=FALSE, log_scale=TRUE, x_lab=x_tp53, min_samples=5, chr=NA)),
    paste0("Linear model plot (affected, type=", ifelse(is.na(types[i]), "all", types[i]), ")")
  )
  titled_print(p, paste0("Linear model plot (affected, type=", ifelse(is.na(types[i]), "all", types[i]), ")"))
  ggsave(paste0(plot_dir, "counts_cohort/te_count_lm_aff_", ifelse(is.na(types[i]), "all", types[i]), ".png"), plot=p, width=5, height=5)
  p
})

cat("\n===== TE count of affected by TP53 status - bootstap test =====\n")
boot_aff <- bootstrap_test(te_aff_t, "TP53_status", "total", mean, n_bootstraps = 10000, step_size=100)
write_output(quote(bootstrap_test(te_aff_t, "TP53_status", "total", mean, n_bootstraps = 10000, step_size=100)),
             "Bootstrap test (affected, total)")
titled_print(boot_aff$dist_plot, "Bootstrap test (affected, total) - Distribution Plot")
ggsave(paste0(plot_dir, "counts_cohort/te_count_bootstrap_aff_dist_total.png"), plot=boot_aff$dist_plot, width=3, height=5)
titled_print(boot_aff$pval_plot, "Bootstrap test (affected, total) - P-value Convergence Plot")
ggsave(paste0(plot_dir, "counts_cohort/te_count_bootstrap_aff_pval_total.png"), plot=boot_aff$pval_plot, width=3, height=5)

# lfs
cat("\n===== TE COUNT of LFS group by cancer status =====\n")
plots_wilcox_lfs <- vector("list", length(types))
for (i in seq_along(types)) {
  plots_wilcox_lfs[[i]] <- plot_count_wilcox(te_lfs_t, type=types[i], chr=NA, group="Cancer", x_lab="Cancer", y_lab="Total TE count")
  write_output(quote(plot_count_wilcox(te_lfs_t, type=types[i], chr=NA, group="Cancer", x_lab="Cancer", y_lab="Total TE count")),
               paste0("Wilcoxon plot (LFS, type=", ifelse(is.na(types[i]), "all", types[i]), ")"))
  titled_print(plots_wilcox_lfs[[i]], paste0("Wilcoxon plot (LFS, type=", ifelse(is.na(types[i]), "all", types[i]), ")"))
  ggsave(paste0(plot_dir, "counts_clinical_lfs/te_count_wilcox_lfs_", ifelse(is.na(types[i]), "all", types[i]), ".png"), plot=plots_wilcox_lfs[[i]], width = 3, height = 5)
}

cat("\n===== TE COUNT of LFS by cluster =====\n")
plots_kruskal_cluster <- vector("list", length(types))
for (i in seq_along(types)) {
  write_output(quote(plots_kruskal_cluster[[i]] <- plot_count_kruskal(te_lfs_t, type=types[i], chr=NA, group="cluster", x_lab="Cluster", y_lab="Repeat count")),
               paste0("Kruskal plot (LFS cluster, type=", ifelse(is.na(types[i]), "all", types[i]), ")"))
  titled_print(plots_kruskal_cluster[[i]], paste0("Kruskal plot (LFS cluster, type=", ifelse(is.na(types[i]), "all", types[i]), ")"))
  ggsave(paste0(plot_dir, "counts_clinical_lfs/te_count_lfs_cluster_", ifelse(is.na(types[i]), "all", types[i]), ".png"), plot=plots_kruskal_cluster[[i]], width = 9, height = 5)
}

cat("\n===== TE COUNT: By inheritance =====\n")
te_lfs_t <- te_lfs_t[!grepl("DELETE", te_lfs_t$inheritance), ]
plots_kruskal_inheritance <- vector("list", length(types))
for (i in seq_along(types)) {
  write_output(quote(plots_kruskal_inheritance[[i]] <- plot_count_kruskal(te_lfs_t, type=types[i], chr=NA, group="inheritance", x_lab="Inheritance", y_lab="Repeat count")),
               paste0("Kruskal plot (LFS inheritance, type=", ifelse(is.na(types[i]), "all", types[i]), ")"))
  titled_print(plots_kruskal_inheritance[[i]], paste0("Kruskal plot (LFS inheritance, type=", ifelse(is.na(types[i]), "all", types[i]), ")"))
  ggsave(paste0(plot_dir, "counts_clinical_lfs/te_count_lfs_inheritance_", ifelse(is.na(types[i]), "all", types[i]), ".png"), plot=plots_kruskal_inheritance[[i]], width = 9, height = 5)
}

cat("\n===== TE COUNT of LFS by variant location =====\n")
plots_kruskal_varloc <- vector("list", length(types))
for (i in seq_along(types)) {
  write_output(quote(plots_kruskal_varloc[[i]] <- plot_count_kruskal(te_lfs_t, type=types[i], chr=NA, group="Variant_location", x_lab="TP53 variant location", y_lab="Total TE count")),
               paste0("Kruskal plot (LFS variant location, type=", ifelse(is.na(types[i]), "all", types[i]), ")"))
  titled_print(plots_kruskal_varloc[[i]], paste0("Kruskal plot (LFS variant location, type=", ifelse(is.na(types[i]), "all", types[i]), ")"))
  ggsave(paste0(plot_dir, "counts_clinical_lfs/te_count_lfs_varloc_", ifelse(is.na(types[i]), "all", types[i]), ".png"), plot=plots_kruskal_varloc[[i]], width = 9, height = 5)
}

cat("\n===== TE COUNT of LFS by variant classification =====\n")
plots_kruskal_varclass <- vector("list", length(types))
for (i in seq_along(types)) {
  write_output(quote(plots_kruskal_varclass[[i]] <- plot_count_kruskal(te_lfs_t, type=types[i], chr=NA, group="Variant_Classification", x_lab="TP53 variant classification", y_lab="Total TE count")),
               paste0("Kruskal plot (LFS variant classification, type=", ifelse(is.na(types[i]), "all", types[i]), ")"))
  titled_print(plots_kruskal_varclass[[i]], paste0("Kruskal plot (LFS variant classification, type=", ifelse(is.na(types[i]), "all", types[i]), ")"))
  ggsave(paste0(plot_dir, "counts_clinical_lfs/te_count_lfs_varclass_", ifelse(is.na(types[i]), "all", types[i]), ".png"), plot=plots_kruskal_varclass[[i]], width = 9, height = 5)
}


cat("\n===== COUNT PER CHROMOSOME =====\n")
for (i in seq_along(types)) {
  write_output(quote(plot_count_perchr_notest(te_kics_t, type=types[i], y_lab="Total TE count normalized by chromosome length", log_scale=FALSE)),
               paste0("Count per chromosome (KICS, type=", ifelse(is.na(types[i]), "all", types[i]), ")"))
  titled_print(plot_count_perchr_notest(te_kics_t, type=types[i], y_lab="Total TE count normalized by chromosome length", log_scale=FALSE),
               paste0("Count per chromosome (KICS, type=", ifelse(is.na(types[i]), "all", types[i]), ")"))
  ggsave(paste0(plot_dir, "other/te_count_perchr_kics_", ifelse(is.na(types[i]), "all", types[i]), ".png"), width = 9, height = 5)
}

cat("\n===== COUNT PER CHROMOSOME BY TP53 STATUS =====\n")
for (i in seq_along(types)) {
  write_output(quote(plot_count_perchr(te_aff_t, type=types[i], group="TP53_status", y_lab="TE count normalized by chromosome length")),
               paste0("Count per chromosome (affected, type=", ifelse(is.na(types[i]), "all", types[i]), ", by TP53_status)"))
  titled_print(plot_count_perchr(te_aff_t, type=types[i], group="TP53_status", y_lab="TE count normalized by chromosome length"),
               paste0("Count per chromosome (affected, type=", ifelse(is.na(types[i]), "all", types[i]), ", by TP53_status)"))
  ggsave(paste0(plot_dir, "other/te_count_perchr_aff_tp53_", ifelse(is.na(types[i]), "all", types[i]), ".png"),
         plot=plot_count_perchr(te_aff_t, type=types[i], group="TP53_status", y_lab="TE count normalized by chromosome length"), width=9, height=5)
}

# By Cancer (LFS)
for (i in seq_along(types)) {
  write_output(quote(plot_count_perchr(te_lfs_t, type=types[i], group="Cancer", y_lab="TE count normalized by chromosome length")),
               paste0("Count per chromosome (LFS, type=", ifelse(is.na(types[i]), "all", types[i]), ", by Cancer)"))
  titled_print(plot_count_perchr(te_lfs_t, type=types[i], group="Cancer", y_lab="TE count normalized by chromosome length"),
               paste0("Count per chromosome (LFS, type=", ifelse(is.na(types[i]), "all", types[i]), ", by Cancer)"))
  ggsave(paste0(plot_dir, "other/te_count_perchr_lfs_cancer_", ifelse(is.na(types[i]), "all", types[i]), ".png"),
         plot=plot_count_perchr(te_lfs_t, type=types[i], group="Cancer", y_lab="TE count normalized by chromosome length"), width=9, height=5)
}

cat("\n===== BY TUMOUR TYPE (AFF) =====\n")
for (i in seq_along(types)) {
  write_output(quote(plot_count_kruskal_nogroup(te_aff_t, column="tumor_type", x_lab="Tumour type", y_lab="Repeat count", type=types[i], chr=NA, min=3, log_scale = TRUE)),
               paste0("Kruskal (no group) by tumour type (affected, type=", ifelse(is.na(types[i]), "all", types[i]), ")"))
  titled_print(plot_count_kruskal_nogroup(te_aff_t, column="tumor_type", x_lab="Tumour type", y_lab = "Repeat count", type=types[i], chr=NA, min=3, log_scale = TRUE),
               paste0("Kruskal (no group) by tumour type (affected, type=", ifelse(is.na(types[i]), "all", types[i]), ")"))
  ggsave(paste0(plot_dir, "counts_tumour_type/te_count_tt_aff_", ifelse(is.na(types[i]), "all", types[i]), ".png"),
         plot=plot_count_kruskal_nogroup(te_aff_t, column="tumor_type", x_lab="Tumour type", y_lab="Repeat count", type=types[i], chr=NA, min=3, log_scale = TRUE), width=9, height=5)
}

cat("\n===== BY TUMOUR TYPE (KICS) =====\n")
for (i in seq_along(types)) {
  write_output(quote(plot_count_kruskal_nogroup(te_kics_t, column="tumor_type", x_lab="Tumour type", y_lab="Repeat count", type=types[i], chr=NA, min=3, log_scale=TRUE, breaks=c(10, 100, 1000, 10000, 100000))),
               paste0("Kruskal (no group) by tumour type (KICS, type=", ifelse(is.na(types[i]), "all", types[i]), ")"))
  titled_print(plot_count_kruskal_nogroup(te_kics_t, column="tumor_type", x_lab="Tumour type", y_lab="Repeat count", type=types[i], chr=NA, min=3, log_scale=TRUE, breaks=c(10, 100, 1000, 10000, 100000)),
               paste0("Kruskal (no group) by tumour type (KICS, type=", ifelse(is.na(types[i]), "all", types[i]), ")"))
  ggsave(paste0(plot_dir, "counts_tumour_type/te_count_tt_kics_", ifelse(is.na(types[i]), "all", types[i]), ".png"),
         plot=plot_count_kruskal_nogroup(te_kics_t, column="tumor_type", x_lab="Tumour type", y_lab="Repeat count", type=types[i], chr=NA, min=3, log_scale=TRUE, breaks=c(10, 100, 1000, 10000, 100000)), width=9, height=5)
}

cat("\n===== BY TUMOUR TYPE (LFS) =====\n")
for (i in seq_along(types)) {
  write_output(quote(plot_count_kruskal_nogroup(te_lfs_t, column="tumor_type", x_lab="Tumour type", y_lab="Repeat count", type=types[i], chr=NA, min=3, log_scale=TRUE)),
               paste0("Kruskal (no group) by tumour type (LFS, type=", ifelse(is.na(types[i]), "all", types[i]), ")"))
  titled_print(plot_count_kruskal_nogroup(te_lfs_t, column="tumor_type", x_lab="Tumour type", y_lab="Repeat count", type=types[i], chr=NA, min=3, log_scale=TRUE),
               paste0("Kruskal (no group) by tumour type (LFS, type=", ifelse(is.na(types[i]), "all", types[i]), ")"))
  ggsave(paste0(plot_dir, "counts_tumour_type/te_count_tt_lfs_", ifelse(is.na(types[i]), "all", types[i]), ".png"),
         plot=plot_count_kruskal_nogroup(te_lfs_t, column="tumor_type", x_lab="Tumour type", y_lab="Repeat count", type=types[i], chr=NA, min=3, log_scale=TRUE), width=9, height=5)
}

# Export top 10% samples by total insertions
cat("\n===== EXPORTING TOP 10% SAMPLES BY TOTAL INSERTIONS =====\n")
top10_threshold <- quantile(te_aff_t$total, 0.90, na.rm = TRUE)
te_aff_top10 <- te_aff_t %>% filter(total >= top10_threshold)
cat("Top 10% threshold:", top10_threshold, "\n")
cat("Number of samples in top 10%:", nrow(te_aff_top10), "\n")
write.csv(te_aff_top10, paste0(r_dir_files, "te_aff_top10_samples.csv"), row.names = FALSE, quote = FALSE)
cat("Exported te_aff_top10_samples.csv\n")

# Export outlier samples by tumor type (types with >=3 samples)
cat("\n===== EXPORTING OUTLIER SAMPLES BY TUMOR TYPE =====\n")
outlier_samples <- te_aff_t %>%
  group_by(tumor_type) %>%
  filter(n() >= 3) %>%
  mutate(
    Q1 = quantile(total, 0.25, na.rm = TRUE),
    Q3 = quantile(total, 0.75, na.rm = TRUE),
    IQR = Q3 - Q1,
    lower_bound = Q1 - 1.5 * IQR,
    upper_bound = Q3 + 1.5 * IQR,
    is_outlier = total < lower_bound | total > upper_bound
  ) %>%
  filter(is_outlier) %>%
  ungroup() %>%
  select(-Q1, -Q3, -IQR, -lower_bound, -upper_bound, -is_outlier)

cat("Number of outlier samples:", nrow(outlier_samples), "\n")
if (nrow(outlier_samples) > 0) {
  cat("Outliers by tumor type:\n")
  print(table(outlier_samples$tumor_type))
  write.csv(outlier_samples, paste0(r_dir_files, "te_aff_outliers_by_tumor_type.csv"), row.names = FALSE, quote = FALSE)
  cat("Exported te_aff_outliers_by_tumor_type.csv\n")
} else {
  cat("No outliers detected\n")
}

cat("\n===== COUNT BY TUMOUR TYPE AND TP53 STATUS =====\n")
for (i in seq_along(types)) {
  write_output(quote(plot_count_tt(te_aff_t, chr=NA, type=types[i], group="TP53_status", min=3, y_lab="Repeat count", legend_title=x_tp53, log_scale=TRUE)),
               paste0("Count by tumour type and TP53_status (affected, type=", ifelse(is.na(types[i]), "all", types[i]), ")"))
  titled_print(plot_count_tt(te_aff_t, chr=NA, type=types[i], group="TP53_status", min=3, y_lab="Repeat count", legend_title=x_tp53, log_scale=TRUE),
               paste0("Count by tumour type and TP53_status (affected, type=", ifelse(is.na(types[i]), "all", types[i]), ")"))
  ggsave(paste0(plot_dir, "counts_tumour_type/te_count_tt_tp53_", ifelse(is.na(types[i]), "all", types[i]), ".png"),
         plot=plot_count_tt(te_aff_t, chr=NA, type=types[i], group="TP53_status", min=3, y_lab="Repeat count", legend_title=x_tp53, log_scale=TRUE), width=9, height=5)
}

cat("\n===== LFS BY TP53 VARIANT LOCATION =====\n")
for (i in seq_along(types)) {
  write_output(quote(plot_count_kruskal_nogroup(te_lfs_t, column="Variant_location", x_lab="TP53 Variant Location", y_lab="Repeat count", type=types[i], chr=NA, min=3, log_scale=TRUE)),
               paste0("Kruskal by TP53 variant location (LFS, type=", ifelse(is.na(types[i]), "all", types[i]), ")"))
  titled_print(plot_count_kruskal_nogroup(te_lfs_t, column="Variant_location", x_lab="TP53 Variant Location", y_lab="Repeat count", type=types[i], chr=NA, min=3, log_scale=TRUE),
               paste0("Kruskal by TP53 variant location (LFS, type=", ifelse(is.na(types[i]), "all", types[i]), ")"))
  ggsave(paste0(plot_dir, "counts_clinical_lfs/te_count_lfs_variant_location_", ifelse(is.na(types[i]), "all", types[i]), ".png"),
         plot=plot_count_kruskal_nogroup(te_lfs_t, column="Variant_location", x_lab="TP53 Variant Location", y_lab="Repeat count", type=types[i], chr=NA, min=3, log_scale=TRUE), width=9, height=5)
}

cat("\n===== LFS BY TP53 CLUSTER =====\n")
for (i in seq_along(types)) {
  write_output(quote(plot_count_kruskal_nogroup(te_lfs_t, column="cluster", x_lab="TP53 Cluster", y_lab="Repeat count", type=types[i], chr=NA, min=3, log_scale=TRUE)),
               paste0("Kruskal by TP53 cluster (LFS, type=", ifelse(is.na(types[i]), "all", types[i]), ")"))
  titled_print(plot_count_kruskal_nogroup(te_lfs_t, column="cluster", x_lab="TP53 Cluster", y_lab="Repeat count", type=types[i], chr=NA, min=3, log_scale=TRUE),
               paste0("Kruskal by TP53 cluster (LFS, type=", ifelse(is.na(types[i]), "all", types[i]), ")"))
  ggsave(paste0(plot_dir, "counts_clinical_lfs/te_count_lfs_cluster_", ifelse(is.na(types[i]), "all", types[i]), ".png"),
         plot=plot_count_kruskal_nogroup(te_lfs_t, column="cluster", x_lab="TP53 Cluster", y_lab="Repeat count", type=types[i], chr=NA, min=3, log_scale=TRUE), width=9, height=5)
}

cat("\n===== SPECIFIC TEs BY TUMOUR TYPE (Fisher test) =====\n")
write_output(quote({
  specific_te_bytt_tt <- as.data.frame(fisher_test_by_tumor_type(te_aff_expand_t, min_samples_tt=5, min_samples_te=5))
}), "Specific TE by tumour type (Fisher test)")

cat("\n===== COUNT BY CLINICAL VARIABLES =====\n")
for (i in seq_along(types)) {
  write_output(quote(plot_count_age(te_kics_t, chr=NA, type=types[i], y_lab="Repeat count")),
               paste0("Count by age (KICS, type=", ifelse(is.na(types[i]), "all", types[i]), ")"))
  titled_print(plot_count_age(te_kics_t, chr=NA, type=types[i], y_lab="Repeat count"),
               paste0("Count by age (KICS, type=", ifelse(is.na(types[i]), "all", types[i]), ")"))
  ggsave(paste0(plot_dir, "counts_clinical_kics/te_count_age_kics_", ifelse(is.na(types[i]), "all", types[i]), ".png"),
         plot=plot_count_age(te_kics_t, chr=NA, type=types[i], y_lab="Repeat count"), width=9, height=5)
}

cat("\n===== WILCOXON TESTS BY TREATMENT =====\n")
# Standardize treatment values to fix case inconsistency
te_aff_t$treatment <- tolower(te_aff_t$treatment)
for (i in seq_along(types)) {
  plot_wilcox_treatment <- plot_count_wilcox(te_aff_t, chr=NA, type=types[i], group="treatment", x_lab="Treatment status", y_lab="Repeat count", log_scale=TRUE)
  write_output(quote(plot_count_wilcox(te_aff_t, chr=NA, type=types[i], group="treatment", x_lab="Treatment status", y_lab="Repeat count", log_scale=TRUE)),
               paste0("Wilcoxon by treatment (affected, type=", ifelse(is.na(types[i]), "all", types[i]), ")"))
  titled_print(plot_wilcox_treatment,
               paste0("Wilcoxon by treatment (affected, type=", ifelse(is.na(types[i]), "all", types[i]), ")"))
  ggsave(paste0(plot_dir, "counts_clinical_kics/te_count_wilcox_treatment_", ifelse(is.na(types[i]), "all", types[i]), ".png"),
         plot=plot_wilcox_treatment, width=9, height=5)
}

cat("\n===== LINEAR MODEL BY TREATMENT =====\n")
for (i in seq_along(types)) {
  y_label <- ifelse(is.na(types[i]), "Total TE count", paste0(types[i], " count"))
  p <- plot_count_lm(te_aff_t, min_samples=3, residuals=FALSE, log_scale=TRUE, group="treatment", covariates = covar_med, x_lab="Disease state", y_lab=y_label, type=types[i], chr=NA)
  write_output(quote(plot_count_lm(te_aff_t, min_samples=3, residuals=FALSE, log_scale=TRUE, group="treatment", covariates = covar_med, x_lab="Disease state", y_lab=y_label, type=types[i], chr=NA)),
               paste0("Linear model by treatment (affected, type=", ifelse(is.na(types[i]), "all", types[i]), ")"))
  titled_print(p, paste0("Linear model by treatment (affected, type=", ifelse(is.na(types[i]), "all", types[i]), ")"))
  ggsave(paste0(plot_dir, "counts_clinical_kics/te_count_lm_treatment_", ifelse(is.na(types[i]), "all", types[i]), ".png"), plot=p, width=9, height=5)
}

cat("\n===== KRUSKAL TEST BY LESION TYPE =====\n")
for (i in seq_along(types)) {
  p <- plot_count_kruskal(te_aff_t, chr=NA, type=types[i], group="lesion_type", x_lab="Lesion type", y_lab="Repeat count", log_scale=TRUE)
  write_output(quote(plot_count_kruskal(te_aff_t, chr=NA, type=types[i], group="lesion_type", x_lab="Lesion type", y_lab="Repeat count", log_scale=TRUE)),
               paste0("Kruskal by lesion type (affected, type=", ifelse(is.na(types[i]), "all", types[i]), ")"))
  titled_print(p, paste0("Kruskal by lesion type (affected, type=", ifelse(is.na(types[i]), "all", types[i]), ")"))
  ggsave(paste0(plot_dir, "counts_clinical_kics/te_count_kruskal_lesion_", ifelse(is.na(types[i]), "all", types[i]), ".png"), plot=p, width=9, height=5)
}

cat("\n===== LINEAR MODEL BY LESION TYPE =====\n")
for (i in seq_along(types)) {
  y_label <- ifelse(is.na(types[i]), "Total TE count", paste0(types[i], " count"))
  p <- plot_count_lm(te_aff_t, min_samples=5, residuals=FALSE, log_scale=TRUE, group="lesion_type", covariates = covar_med, x_lab="Disease state", y_lab=y_label, type=types[i], chr=NA)
  write_output(quote(plot_count_lm(te_aff_t, min_samples=5, residuals=FALSE, log_scale=TRUE, group="lesion_type", covariates = covar_med, x_lab="Disease state", y_lab=y_label, type=types[i], chr=NA)),
               paste0("Linear model by lesion type (affected, type=", ifelse(is.na(types[i]), "all", types[i]), ")"))
  titled_print(p, paste0("Linear model by lesion type (affected, type=", ifelse(is.na(types[i]), "all", types[i]), ")"))
  ggsave(paste0(plot_dir, "counts_clinical_kics/te_count_lm_lesion_", ifelse(is.na(types[i]), "all", types[i]), ".png"), plot=p, width=9, height=5)
}

cat("\n===== KRUSKAL TEST BY DISEASE STATE =====\n")
for (i in seq_along(types)) {
  p <- plot_count_kruskal(te_aff_t, chr=NA, type=types[i], group="disease_state", x_lab="Disease state", y_lab="Repeat count", log_scale=TRUE)
  write_output(quote(plot_count_kruskal(te_aff_t, chr=NA, type=types[i], group="disease_state", x_lab="Disease state", y_lab="Repeat count", log_scale=TRUE)),
               paste0("Kruskal by disease state (affected, type=", ifelse(is.na(types[i]), "all", types[i]), ")"))
  titled_print(p, paste0("Kruskal by disease state (affected, type=", ifelse(is.na(types[i]), "all", types[i]), ")"))
  ggsave(paste0(plot_dir, "counts_clinical_kics/te_count_kruskal_disease_", ifelse(is.na(types[i]), "all", types[i]), ".png"), plot=p, width=9, height=5)
}

cat("\n===== LINEAR MODEL BY DISEASE STATE =====\n")
for (i in seq_along(types)) {
  y_label <- ifelse(is.na(types[i]), "Total TE count", paste0(types[i], " count"))
  p <- plot_count_lm(te_aff_t, min_samples=5, residuals=FALSE, log_scale=TRUE, group="disease_state", covariates = covar_med, x_lab="Disease state", y_lab=y_label, type=types[i], chr=NA)
  write_output(quote(plot_count_lm(te_aff_t, min_samples=5, residuals=FALSE, log_scale=TRUE, group="disease_state", covariates = covar_med, x_lab="Disease state", y_lab=y_label, type=types[i], chr=NA)),
               paste0("Linear model by disease state (affected, type=", ifelse(is.na(types[i]), "all", types[i]), ")"))
  titled_print(p, paste0("Linear model by disease state (affected, type=", ifelse(is.na(types[i]), "all", types[i]), ")"))
  ggsave(paste0(plot_dir, "counts_clinical_kics/te_count_lm_disease_", ifelse(is.na(types[i]), "all", types[i]), ".png"), plot=p, width=9, height=5)
}

cat("\n===== KRUSKAL TEST BY SEX =====\n")
for (i in seq_along(types)) {
  p <- plot_count_kruskal(te_aff_t, chr=NA, type=types[i], group="sex", x_lab="Sex", y_lab="Repeat count", log_scale=TRUE)
  write_output(quote(plot_count_kruskal(te_aff_t, chr=NA, type=types[i], group="sex", x_lab="Sex", y_lab="Repeat count", log_scale=TRUE)),
               paste0("Kruskal by sex (affected, type=", ifelse(is.na(types[i]), "all", types[i]), ")"))
  titled_print(p, paste0("Kruskal by sex (affected, type=", ifelse(is.na(types[i]), "all", types[i]), ")"))
  ggsave(paste0(plot_dir, "counts_clinical_kics/te_count_kruskal_sex_", ifelse(is.na(types[i]), "all", types[i]), ".png"), plot=p, width=9, height=5)
}

cat("\n===== LINEAR MODEL BY SEX =====\n")
for (i in seq_along(types)) {
  y_label <- ifelse(is.na(types[i]), "Total TE count", paste0(types[i], " count"))
  p <- plot_count_lm(te_aff_t, min_samples=5, residuals=FALSE, group="sex", log_scale=TRUE, covariates = covar_med, x_lab="Sex", y_lab=y_label, type=types[i], chr=NA)
  write_output(quote(plot_count_lm(te_aff_t, min_samples=5, residuals=FALSE, group="sex", log_scale=TRUE, covariates = covar_med, x_lab="Sex", y_lab=y_label, type=types[i], chr=NA)),
               paste0("Linear model by sex (affected, type=", ifelse(is.na(types[i]), "all", types[i]), ")"))
  titled_print(p, paste0("Linear model by sex (affected, type=", ifelse(is.na(types[i]), "all", types[i]), ")"))
  ggsave(paste0(plot_dir, "counts_clinical_kics/te_count_lm_sex_", ifelse(is.na(types[i]), "all", types[i]), ".png"), plot=p, width=9, height=5)
}


cat("\n===== P53 MUTATION AND FITNESS =====\n")
cat("Plotting TE count by p53 mutation (Kruskal test)...\n")
for (i in seq_along(types)) {
  tryCatch({
    plot_obj <- plot_count_kruskal_nogroup(te_lfs_t, column="mutation", min=3, chr=NA, type=types[i], x_lab="p53 mutation", y_lab="Repeat count", log_scale = TRUE)
    titled_print(plot_obj, paste0("TE count by p53 mutation (", ifelse(is.na(types[i]), "all types", types[i]), ")"))
    ggsave(paste0(plot_dir, "counts_clinical_lfs/te_count_p53_mutation_", ifelse(is.na(types[i]), "all", types[i]), ".png"), plot=plot_obj, width = 9, height = 5)
  }, error = function(e) {
    cat("Skipping p53 mutation plot for type", ifelse(is.na(types[i]), "all", types[i]), "due to insufficient data variation:", e$message, "\n")
  })
}
cat("Plotting TE frequency by mutation along p53 gene...\n")
hotspot <- c(175, 213, 245, 248, 273, 282, 337)
for (i in seq_along(types)) {
  plot_title <- paste0("TE locations by p53 mutation (", ifelse(is.na(types[i]), "all types", types[i]), ")")
  p <- plot_te_locations(te_lfs_t, chr=NA, type=types[i], hotspot=hotspot, log_scale=TRUE)
  titled_print(p, plot_title)
  ggsave(
    paste0(plot_dir, "survival_burden/te_locations_lfs_", ifelse(is.na(types[i]), "all", types[i]), ".png"),
    plot = p, width = 9, height = 5
  )
}


cat("✓ Script completed successfully\n")
