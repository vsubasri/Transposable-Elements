#!/usr/bin/env Rscript

# LFS TE Count by Tumor Type - Colored by Cohort
# Standalone script for generating tumor type plots with cohort coloring
# Uses te_lfs_t data, tumor types with >=3 total samples, colored by cohort

#### SOURCE COMMON SETUP ####
source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/00_viz_common_setup.R")
source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/00_viz_load_data_tumour.R")

#### SETUP OUTPUT ####
output_dir <- paste0(plot_dir, "counts_clinical_lfs/")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Close the PDF device opened by data loading and open our own
dev.off()
pdf(file = paste0(output_dir, "lfs_cohort_plots.pdf"), width = 12, height = 8)

#### GENERATE PLOTS ####
cat("\n===== LFS TE COUNT BY TUMOR TYPE - COLORED BY COHORT =====\n")

# TE types to analyze
types <- c(NA, "LINE1")  # NA = all types combined

# Use te_lfs_t - filter to samples with valid cohort and tumor_type
# Exclude KiCS (not in LFS), keep only tumor types with at least 3 total samples
te_valid <- te_lfs_t %>%
  filter(!is.na(cohort), !is.na(tumor_type), cohort != "KiCS") %>%
  group_by(tumor_type) %>%
  filter(n() >= 3) %>%
  ungroup()

cat("Total samples:", nrow(te_valid), "\n")
cat("Cohort distribution:\n")
print(table(te_valid$cohort))
cat("\nTumor types (>=3 samples):", paste(unique(te_valid$tumor_type), collapse=", "), "\n")
cat("\nTumor type by cohort:\n")
print(table(te_valid$tumor_type, te_valid$cohort))

for (i in seq_along(types)) {
  type_label <- ifelse(is.na(types[i]), "all", types[i])
  y_label <- ifelse(is.na(types[i]), "Repeat count", paste0(types[i], " count"))
  count_col <- ifelse(is.na(types[i]), "total", types[i])

  cat("\n--- Type:", type_label, "---\n")

  tryCatch({
    # Custom plot: separate boxplots per cohort for each tumor type
    p <- ggplot(te_valid, aes(x = tumor_type, y = .data[[count_col]], fill = cohort)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(aes(color = cohort, shape = cohort),
                  position = position_jitterdodge(jitter.width = 0.2), size = 2.5, alpha = 0.8) +
      scale_y_log10() +
      scale_fill_manual(values = c("nick" = "#E41A1C", "SJ" = "#377EB8")) +
      scale_color_manual(values = c("nick" = "#E41A1C", "SJ" = "#377EB8")) +
      scale_shape_manual(values = c("nick" = 16, "SJ" = 17)) +
      labs(x = "Tumor type", y = y_label, fill = "Cohort", color = "Cohort", shape = "Cohort") +
      guides(color = "none", shape = "none") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    titled_print(p, paste0("LFS TE by tumor type & cohort (type=", type_label, ")"))
    ggsave(paste0(output_dir, "te_count_tt_lfs_cohort_", type_label, ".png"),
           plot=p, width=10, height=6)
    cat("Saved: te_count_tt_lfs_cohort_", type_label, ".png\n")
  }, error = function(e) {
    cat("Error:", e$message, "\n")
  })
}

#### CLEANUP ####
dev.off()
cat("\n✓ LFS cohort plots complete\n")
cat("Output directory:", output_dir, "\n")
