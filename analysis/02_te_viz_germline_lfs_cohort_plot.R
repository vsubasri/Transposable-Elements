#!/usr/bin/env Rscript

# Germline LFS TE Count by Tumor Type - Colored by Cohort
# Standalone script for generating tumor type plots with cohort coloring
# Uses te_lfs data, tumor types with >=3 total samples, colored by cohort

#### SOURCE COMMON SETUP ####
source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/00_viz_common_setup.R")

#### LOAD ONLY REQUIRED DATA ####
cat("Loading minimal data for LFS cohort plot...\n")

# Set data directory based on TEST_MODE
if (TEST_MODE) {
  r_dir <- "/Users/briannelaverty/Documents/R_Malkin/te/data/R_obj/test_data/germline/"
  cat("*** TEST MODE ACTIVE: Using test data from", r_dir, "***\n\n")
} else {
  r_dir <- "/Users/briannelaverty/Documents/R_Malkin/te/data/R_obj/germline/"
}

# Set output directory
plot_dir <- "/Users/briannelaverty/Library/Mobile Documents/com~apple~CloudDocs/Malkin_lab/figures/te/graphs/germline/"

# Load only the count matrix data (contains te_lfs and te_aff)
load(paste0(r_dir, "final_te_count_rare.RData"))

if (TEST_MODE) {
  te_lfs <- final_te_count_test$te_lfs
  te_aff <- final_te_count_test$te_aff
} else {
  te_lfs <- final_te_count$te_lfs
  te_aff <- final_te_count$te_aff
}

cat("✓ Loaded te_lfs:", nrow(te_lfs), "samples\n")
cat("✓ Loaded te_aff:", nrow(te_aff), "samples\n")

#### SETUP OUTPUT ####
output_dir <- paste0(plot_dir, "counts_clinical_lfs/")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

pdf(file = paste0(output_dir, "lfs_cohort_plots.pdf"), width = 12, height = 8)

#### GENERATE PLOTS ####
cat("\n===== GERMLINE LFS TE COUNT BY TUMOR TYPE - COLORED BY COHORT =====\n")

# TE types to analyze
types <- c(NA, "LINE1")  # NA = all types combined

# Use te_lfs - filter to samples with valid cohort and tumor_type
# Exclude HostSeq and Taylor (not LFS clinical), keep only tumor types with at least 3 total samples
te_lfs_valid <- te_lfs %>%
  filter(!is.na(cohort), !is.na(tumor_type), !cohort %in% c("HostSeq", "Taylor")) %>%
  group_by(tumor_type) %>%
  filter(n() >= 3) %>%
  ungroup()

cat("Total samples:", nrow(te_lfs_valid), "\n")
cat("Cohort distribution:\n")
print(table(te_lfs_valid$cohort))
cat("\nTumor types (>=3 samples):", paste(unique(te_lfs_valid$tumor_type), collapse=", "), "\n")
cat("\nTumor type by cohort:\n")
print(table(te_lfs_valid$tumor_type, te_lfs_valid$cohort))

# Get unique cohorts for color palette
cohorts <- unique(te_lfs_valid$cohort)
n_cohorts <- length(cohorts)
cat("\nNumber of cohorts:", n_cohorts, "\n")

# Create color palette based on number of cohorts
if (n_cohorts == 2) {
  cohort_colors <- c("#E41A1C", "#377EB8")
  cohort_shapes <- c(16, 17)
} else {
  cohort_colors <- scales::hue_pal()(n_cohorts)
  cohort_shapes <- rep(c(16, 17, 15, 18), length.out = n_cohorts)
}
names(cohort_colors) <- cohorts
names(cohort_shapes) <- cohorts

for (i in seq_along(types)) {
  type_label <- ifelse(is.na(types[i]), "all", types[i])
  y_label <- ifelse(is.na(types[i]), "Repeat count", paste0(types[i], " count"))
  count_col <- ifelse(is.na(types[i]), "total", types[i])

  cat("\n--- Type:", type_label, "---\n")

  tryCatch({
    # Custom plot: separate boxplots per cohort for each tumor type
    p <- ggplot(te_lfs_valid, aes(x = tumor_type, y = .data[[count_col]], fill = cohort)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(aes(color = cohort, shape = cohort),
                  position = position_jitterdodge(jitter.width = 0.2), size = 2.5, alpha = 0.8) +
      scale_y_log10() +
      scale_fill_manual(values = cohort_colors) +
      scale_color_manual(values = cohort_colors) +
      scale_shape_manual(values = cohort_shapes) +
      labs(x = "Tumor type", y = y_label, fill = "Cohort", color = "Cohort", shape = "Cohort") +
      guides(color = "none", shape = "none") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    titled_print(p, paste0("Germline LFS TE by tumor type & cohort (type=", type_label, ")"))
    ggsave(paste0(output_dir, "te_count_tt_lfs_cohort_", type_label, ".png"),
           plot=p, width=10, height=6)
    cat("Saved: te_count_tt_lfs_cohort_", type_label, ".png\n")
  }, error = function(e) {
    cat("Error:", e$message, "\n")
  })
}

#### ANCESTRY PLOTS ####
cat("\n===== GERMLINE LFS TE COUNT BY ANCESTRY - COLORED BY COHORT =====\n")

# Filter to samples with valid ancestry
te_lfs_ancestry <- te_lfs %>%
  filter(!is.na(cohort), !is.na(predicted_ancestry_thres), !cohort %in% c("HostSeq", "Taylor")) %>%
  group_by(predicted_ancestry_thres) %>%
  filter(n() >= 3) %>%
  ungroup()

cat("Total samples with ancestry:", nrow(te_lfs_ancestry), "\n")
cat("Cohort distribution:\n")
print(table(te_lfs_ancestry$cohort))
cat("\nAncestry groups (>=3 samples):", paste(unique(te_lfs_ancestry$predicted_ancestry_thres), collapse=", "), "\n")
cat("\nAncestry by cohort:\n")
print(table(te_lfs_ancestry$predicted_ancestry_thres, te_lfs_ancestry$cohort))

# Get unique cohorts for color palette (reuse from above if same)
cohorts_anc <- unique(te_lfs_ancestry$cohort)
n_cohorts_anc <- length(cohorts_anc)

if (n_cohorts_anc == 2) {
  cohort_colors_anc <- c("#E41A1C", "#377EB8")
  cohort_shapes_anc <- c(16, 17)
} else {
  cohort_colors_anc <- scales::hue_pal()(n_cohorts_anc)
  cohort_shapes_anc <- rep(c(16, 17, 15, 18), length.out = n_cohorts_anc)
}
names(cohort_colors_anc) <- cohorts_anc
names(cohort_shapes_anc) <- cohorts_anc

for (i in seq_along(types)) {
  type_label <- ifelse(is.na(types[i]), "all", types[i])
  y_label <- ifelse(is.na(types[i]), "Repeat count", paste0(types[i], " count"))
  count_col <- ifelse(is.na(types[i]), "total", types[i])

  cat("\n--- Type:", type_label, "---\n")

  tryCatch({
    # Custom plot: separate boxplots per cohort for each ancestry group
    p <- ggplot(te_lfs_ancestry, aes(x = predicted_ancestry_thres, y = .data[[count_col]], fill = cohort)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(aes(color = cohort, shape = cohort),
                  position = position_jitterdodge(jitter.width = 0.2), size = 2.5, alpha = 0.8) +
      scale_y_log10() +
      scale_fill_manual(values = cohort_colors_anc) +
      scale_color_manual(values = cohort_colors_anc) +
      scale_shape_manual(values = cohort_shapes_anc) +
      labs(x = "Ancestry", y = y_label, fill = "Cohort", color = "Cohort", shape = "Cohort") +
      guides(color = "none", shape = "none") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    titled_print(p, paste0("Germline LFS TE by ancestry & cohort (type=", type_label, ")"))
    ggsave(paste0(output_dir, "te_count_ancestry_lfs_cohort_", type_label, ".png"),
           plot=p, width=10, height=6)
    cat("Saved: te_count_ancestry_lfs_cohort_", type_label, ".png\n")
  }, error = function(e) {
    cat("Error:", e$message, "\n")
  })
}

#### ANCESTRY PLOTS (MAPPED LABEL) ####
cat("\n===== GERMLINE LFS TE COUNT BY ANCESTRY (MAPPED LABEL) - COLORED BY COHORT =====\n")

# Filter to samples with valid mapped_label
te_lfs_mapped <- te_lfs %>%
  filter(!is.na(cohort), !is.na(mapped_label), !cohort %in% c("HostSeq", "Taylor")) %>%
  group_by(mapped_label) %>%
  filter(n() >= 3) %>%
  ungroup()

cat("Total samples with mapped_label:", nrow(te_lfs_mapped), "\n")
cat("Cohort distribution:\n")
print(table(te_lfs_mapped$cohort))
cat("\nMapped label groups (>=3 samples):", paste(unique(te_lfs_mapped$mapped_label), collapse=", "), "\n")
cat("\nMapped label by cohort:\n")
print(table(te_lfs_mapped$mapped_label, te_lfs_mapped$cohort))

# Get unique cohorts for color palette
cohorts_map <- unique(te_lfs_mapped$cohort)
n_cohorts_map <- length(cohorts_map)

if (n_cohorts_map == 2) {
  cohort_colors_map <- c("#E41A1C", "#377EB8")
  cohort_shapes_map <- c(16, 17)
} else {
  cohort_colors_map <- scales::hue_pal()(n_cohorts_map)
  cohort_shapes_map <- rep(c(16, 17, 15, 18), length.out = n_cohorts_map)
}
names(cohort_colors_map) <- cohorts_map
names(cohort_shapes_map) <- cohorts_map

for (i in seq_along(types)) {
  type_label <- ifelse(is.na(types[i]), "all", types[i])
  y_label <- ifelse(is.na(types[i]), "Repeat count", paste0(types[i], " count"))
  count_col <- ifelse(is.na(types[i]), "total", types[i])

  cat("\n--- Type:", type_label, "---\n")

  tryCatch({
    # Custom plot: separate boxplots per cohort for each mapped_label group
    p <- ggplot(te_lfs_mapped, aes(x = mapped_label, y = .data[[count_col]], fill = cohort)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(aes(color = cohort, shape = cohort),
                  position = position_jitterdodge(jitter.width = 0.2), size = 2.5, alpha = 0.8) +
      scale_y_log10() +
      scale_fill_manual(values = cohort_colors_map) +
      scale_color_manual(values = cohort_colors_map) +
      scale_shape_manual(values = cohort_shapes_map) +
      labs(x = "Ancestry (Mapped Label)", y = y_label, fill = "Cohort", color = "Cohort", shape = "Cohort") +
      guides(color = "none", shape = "none") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    titled_print(p, paste0("Germline LFS TE by mapped_label & cohort (type=", type_label, ")"))
    ggsave(paste0(output_dir, "te_count_mapped_label_lfs_cohort_", type_label, ".png"),
           plot=p, width=10, height=6)
    cat("Saved: te_count_mapped_label_lfs_cohort_", type_label, ".png\n")
  }, error = function(e) {
    cat("Error:", e$message, "\n")
  })
}

#### TE_AFF ANCESTRY PLOTS (PREDICTED ANCESTRY) ####
cat("\n===== GERMLINE AFF TE COUNT BY ANCESTRY (PREDICTED) - COLORED BY COHORT =====\n")

# Filter to samples with valid predicted_ancestry_thres
te_aff_ancestry <- te_aff %>%
  filter(!is.na(cohort), !is.na(predicted_ancestry_thres)) %>%
  group_by(predicted_ancestry_thres) %>%
  filter(n() >= 3) %>%
  ungroup()

cat("Total samples with ancestry:", nrow(te_aff_ancestry), "\n")
cat("Cohort distribution:\n")
print(table(te_aff_ancestry$cohort))
cat("\nAncestry groups (>=3 samples):", paste(unique(te_aff_ancestry$predicted_ancestry_thres), collapse=", "), "\n")
cat("\nAncestry by cohort:\n")
print(table(te_aff_ancestry$predicted_ancestry_thres, te_aff_ancestry$cohort))

# Get unique cohorts for color palette
cohorts_aff_anc <- unique(te_aff_ancestry$cohort)
n_cohorts_aff_anc <- length(cohorts_aff_anc)

if (n_cohorts_aff_anc == 2) {
  cohort_colors_aff_anc <- c("#E41A1C", "#377EB8")
  cohort_shapes_aff_anc <- c(16, 17)
} else {
  cohort_colors_aff_anc <- scales::hue_pal()(n_cohorts_aff_anc)
  cohort_shapes_aff_anc <- rep(c(16, 17, 15, 18), length.out = n_cohorts_aff_anc)
}
names(cohort_colors_aff_anc) <- cohorts_aff_anc
names(cohort_shapes_aff_anc) <- cohorts_aff_anc

for (i in seq_along(types)) {
  type_label <- ifelse(is.na(types[i]), "all", types[i])
  y_label <- ifelse(is.na(types[i]), "Repeat count", paste0(types[i], " count"))
  count_col <- ifelse(is.na(types[i]), "total", types[i])

  cat("\n--- Type:", type_label, "---\n")

  tryCatch({
    p <- ggplot(te_aff_ancestry, aes(x = predicted_ancestry_thres, y = .data[[count_col]], fill = cohort)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(aes(color = cohort, shape = cohort),
                  position = position_jitterdodge(jitter.width = 0.2), size = 2.5, alpha = 0.8) +
      scale_y_log10() +
      scale_fill_manual(values = cohort_colors_aff_anc) +
      scale_color_manual(values = cohort_colors_aff_anc) +
      scale_shape_manual(values = cohort_shapes_aff_anc) +
      labs(x = "Ancestry", y = y_label, fill = "Cohort", color = "Cohort", shape = "Cohort") +
      guides(color = "none", shape = "none") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    titled_print(p, paste0("Germline AFF TE by ancestry & cohort (type=", type_label, ")"))
    ggsave(paste0(output_dir, "te_count_ancestry_aff_cohort_", type_label, ".png"),
           plot=p, width=10, height=6)
    cat("Saved: te_count_ancestry_aff_cohort_", type_label, ".png\n")
  }, error = function(e) {
    cat("Error:", e$message, "\n")
  })
}

#### TE_AFF ANCESTRY PLOTS (MAPPED LABEL) ####
cat("\n===== GERMLINE AFF TE COUNT BY ANCESTRY (MAPPED LABEL) - COLORED BY COHORT =====\n")

# Filter to samples with valid mapped_label
te_aff_mapped <- te_aff %>%
  filter(!is.na(cohort), !is.na(mapped_label)) %>%
  group_by(mapped_label) %>%
  filter(n() >= 3) %>%
  ungroup()

cat("Total samples with mapped_label:", nrow(te_aff_mapped), "\n")
cat("Cohort distribution:\n")
print(table(te_aff_mapped$cohort))
cat("\nMapped label groups (>=3 samples):", paste(unique(te_aff_mapped$mapped_label), collapse=", "), "\n")
cat("\nMapped label by cohort:\n")
print(table(te_aff_mapped$mapped_label, te_aff_mapped$cohort))

# Get unique cohorts for color palette
cohorts_aff_map <- unique(te_aff_mapped$cohort)
n_cohorts_aff_map <- length(cohorts_aff_map)

if (n_cohorts_aff_map == 2) {
  cohort_colors_aff_map <- c("#E41A1C", "#377EB8")
  cohort_shapes_aff_map <- c(16, 17)
} else {
  cohort_colors_aff_map <- scales::hue_pal()(n_cohorts_aff_map)
  cohort_shapes_aff_map <- rep(c(16, 17, 15, 18), length.out = n_cohorts_aff_map)
}
names(cohort_colors_aff_map) <- cohorts_aff_map
names(cohort_shapes_aff_map) <- cohorts_aff_map

for (i in seq_along(types)) {
  type_label <- ifelse(is.na(types[i]), "all", types[i])
  y_label <- ifelse(is.na(types[i]), "Repeat count", paste0(types[i], " count"))
  count_col <- ifelse(is.na(types[i]), "total", types[i])

  cat("\n--- Type:", type_label, "---\n")

  tryCatch({
    p <- ggplot(te_aff_mapped, aes(x = mapped_label, y = .data[[count_col]], fill = cohort)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(aes(color = cohort, shape = cohort),
                  position = position_jitterdodge(jitter.width = 0.2), size = 2.5, alpha = 0.8) +
      scale_y_log10() +
      scale_fill_manual(values = cohort_colors_aff_map) +
      scale_color_manual(values = cohort_colors_aff_map) +
      scale_shape_manual(values = cohort_shapes_aff_map) +
      labs(x = "Ancestry (Mapped Label)", y = y_label, fill = "Cohort", color = "Cohort", shape = "Cohort") +
      guides(color = "none", shape = "none") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    titled_print(p, paste0("Germline AFF TE by mapped_label & cohort (type=", type_label, ")"))
    ggsave(paste0(output_dir, "te_count_mapped_label_aff_cohort_", type_label, ".png"),
           plot=p, width=10, height=6)
    cat("Saved: te_count_mapped_label_aff_cohort_", type_label, ".png\n")
  }, error = function(e) {
    cat("Error:", e$message, "\n")
  })
}

#### CLEANUP ####
dev.off()
cat("\n✓ Germline LFS and AFF cohort plots complete\n")
cat("Output directory:", output_dir, "\n")
