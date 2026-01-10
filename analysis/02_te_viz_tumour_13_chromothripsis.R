#!/usr/bin/env Rscript

# Tumour TE Visualization - Chromothripsis Analysis
# Association between TE burden and chromothripsis / SV metrics
# Runs analysis for both ShatterSeek and GRIDSS SV sources

# Source common setup and load data
source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/00_viz_common_setup.R")
REQUIRED_DATA <- c("count_matrix", "clinical")
source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/00_viz_load_data_tumour.R")

cat("Running 02_te_viz_tumour_13_chromothripsis.R...\n")

# ============================================================================
# SETUP
# ============================================================================

# Data paths
chrom_path <- "/Users/briannelaverty/Documents/R_Malkin/te/data/final/master_chromSummary_fdr_chromothripsis.tsv"
gridss_path <- "/Users/briannelaverty/Documents/R_Malkin/te/data/sv/vis_sv_data_combined.tsv"

chromosomes <- c(as.character(1:22), "X")

# Helper function to format p-values
format_p <- function(p) {
  if (is.na(p)) return("p = NA")
  if (p < 0.001) return("p < 0.001")
  return(paste0("p = ", formatC(p, format = "f", digits = 3)))
}

# ============================================================================
# SECTION 1: LOAD CHROMOTHRIPSIS DATA (shared by both analyses)
# ============================================================================
cat("\n=== Loading Chromothripsis Data ===\n")

chrom_data_raw <- read.delim(chrom_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
cat("Loaded chromothripsis data:", nrow(chrom_data_raw), "rows,", ncol(chrom_data_raw), "columns\n")
cat("Unique samples:", length(unique(chrom_data_raw$sample)), "\n")

# Create chromothripsis status only (no SV metrics - those will come from source-specific data)
chromothripsis_status_only <- chrom_data_raw %>%
  group_by(sample) %>%
  summarise(
    chromothripsis_overall = as.integer(any(chromothripsis %in% c("high", "low"), na.rm = TRUE)),
    chromothripsis_count = sum(chromothripsis %in% c("high", "low"), na.rm = TRUE),
    .groups = "drop"
  )

cat("Samples with chromothripsis:", sum(chromothripsis_status_only$chromothripsis_overall == 1), "\n")
cat("Samples without chromothripsis:", sum(chromothripsis_status_only$chromothripsis_overall == 0), "\n")

# ============================================================================
# SECTION 2: LOAD AND PROCESS SHATTERSEEK SV DATA
# ============================================================================
cat("\n=== Processing ShatterSeek SV Data ===\n")

shatterseek_sv <- chrom_data_raw %>%
  group_by(sample) %>%
  summarise(
    DEL = sum(number_DEL, na.rm = TRUE),
    DUP = sum(number_DUP, na.rm = TRUE),
    h2hINV = sum(number_h2hINV, na.rm = TRUE),
    t2tINV = sum(number_t2tINV, na.rm = TRUE),
    TRA = sum(number_TRA, na.rm = TRUE),
    CNV_segments = sum(number_CNV_segments, na.rm = TRUE),
    SVs_sample = first(number_SVs_sample),
    .groups = "drop"
  )

shatterseek_metrics <- c("DEL", "DUP", "h2hINV", "t2tINV", "TRA", "CNV_segments", "SVs_sample")
cat("ShatterSeek samples:", nrow(shatterseek_sv), "\n")
cat("ShatterSeek SV metrics:", paste(shatterseek_metrics, collapse = ", "), "\n")

# ============================================================================
# SECTION 3: LOAD AND PROCESS GRIDSS SV DATA
# ============================================================================
cat("\n=== Processing GRIDSS SV Data ===\n")

gridss_raw <- read.delim(gridss_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
cat("Loaded GRIDSS data:", nrow(gridss_raw), "rows\n")

# Filter out SGL type
gridss_filtered <- gridss_raw %>% filter(Type != "SGL")
cat("After removing SGL:", nrow(gridss_filtered), "rows\n")

# Extract base sample ID (first part before underscore) to match TE data
gridss_filtered$sample <- sub("_.*", "", gridss_filtered$SampleId)

# Aggregate SV counts per sample
gridss_sv <- gridss_filtered %>%
  group_by(sample) %>%
  summarise(
    BND = sum(Type == "BND", na.rm = TRUE),
    DEL = sum(Type == "DEL", na.rm = TRUE),
    DUP = sum(Type == "DUP", na.rm = TRUE),
    INV = sum(Type == "INV", na.rm = TRUE),
    INS = sum(Type == "INS", na.rm = TRUE),
    INF = sum(Type %in% c("INF", "inf"), na.rm = TRUE),
    SVs_sample = n(),
    .groups = "drop"
  )

gridss_metrics <- c("BND", "DEL", "DUP", "INV", "INS", "INF", "SVs_sample")
cat("GRIDSS samples:", nrow(gridss_sv), "\n")
cat("GRIDSS SV metrics:", paste(gridss_metrics, collapse = ", "), "\n")

# ============================================================================
# SECTION 4: PREPARE TE DATA
# ============================================================================
cat("\n=== Preparing TE Data ===\n")

# Prepare TE counts from te_all_t
te_counts_sample <- te_all_t %>%
  select(sample_full = sample, sample = base_sample.x, total_TE = total,
         tumor_type = tumor_type) %>%
  distinct()

cat("TE samples:", nrow(te_counts_sample), "\n")
cat("Samples with TE = 0:", sum(te_counts_sample$total_TE == 0), "\n")
cat("Samples with TE > 0:", sum(te_counts_sample$total_TE > 0), "\n")

# Add LINE1/ALU/SVA counts from te_all_t (already has these columns)
if (all(c("LINE1", "ALU", "SVA") %in% colnames(te_all_t))) {
  te_counts_sample <- te_counts_sample %>%
    left_join(
      te_all_t %>%
        select(sample_full = sample, LINE1_count = LINE1, ALU_count = ALU, SVA_count = SVA) %>%
        distinct(),
      by = "sample_full"
    ) %>%
    mutate(
      LINE1_count = replace_na(LINE1_count, 0),
      ALU_count = replace_na(ALU_count, 0),
      SVA_count = replace_na(SVA_count, 0)
    )
}

# ============================================================================
# SECTION 5: PREPARE PER-CHROMOSOME DATA
# ============================================================================
cat("\n=== Preparing Per-Chromosome Data ===\n")

# Create per-chromosome TE counts from te_all_t (which has chr1, chr2, etc. columns)
chrom_cols <- paste0("chr", chromosomes)
available_chrom_cols <- intersect(chrom_cols, colnames(te_all_t))

if (length(available_chrom_cols) > 0) {
  te_counts_per_chrom <- te_all_t %>%
    select(sample = base_sample.x, all_of(available_chrom_cols)) %>%
    distinct() %>%
    pivot_longer(cols = all_of(available_chrom_cols), names_to = "chrom_col", values_to = "te_count") %>%
    mutate(chrom = gsub("chr", "", chrom_col)) %>%
    select(sample, chrom, te_count)

  cat("Per-chromosome TE counts created:", nrow(te_counts_per_chrom), "sample-chromosome combinations\n")
  cat("Unique samples:", length(unique(te_counts_per_chrom$sample)), "\n")
} else {
  te_counts_per_chrom <- NULL
  cat("Warning: No chromosome columns found in te_all_t\n")
}

# Create per-chromosome chromothripsis + SV data (keep per-chromosome level)
chrom_level_data <- chrom_data_raw %>%
  mutate(
    chrom = as.character(chrom),
    chrom_chromothripsis_status = case_when(
      chromothripsis == "high" ~ "Chromothripsis+",
      chromothripsis == "low" ~ "Chromothripsis+",
      TRUE ~ "Chromothripsis-"
    ),
    chrom_chromothripsis_level = case_when(
      chromothripsis == "high" ~ "High",
      chromothripsis == "low" ~ "Low",
      TRUE ~ "None"
    )
  ) %>%
  filter(chrom %in% chromosomes)

cat("Per-chromosome chromothripsis data:", nrow(chrom_level_data), "sample-chromosome combinations\n")
cat("Chromosomes with chromothripsis (high):", sum(chrom_level_data$chromothripsis == "high", na.rm = TRUE), "\n")
cat("Chromosomes with chromothripsis (low):", sum(chrom_level_data$chromothripsis == "low", na.rm = TRUE), "\n")

# ============================================================================
# MAIN ANALYSIS FUNCTION - Runs for each SV source
# ============================================================================

run_chromothripsis_analysis <- function(sv_source, sv_data, sv_metrics, output_subdir) {

  cat("\n")
  cat("============================================================\n")
  cat("  RUNNING ANALYSIS FOR:", toupper(sv_source), "\n")
  cat("============================================================\n")

  # Create output directories
  base_dir <- file.path(plot_dir, "chromothripsis", output_subdir)
  files_dir <- file.path(base_dir, "files")
  cancer_type_dir <- file.path(base_dir, "by_cancer_type")
  dir.create(base_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(files_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(cancer_type_dir, showWarnings = FALSE, recursive = TRUE)

  # Open PDF
  pdf(file.path(base_dir, "chromothripsis_graphs.pdf"), width = 10, height = 8, onefile = TRUE)

  # Initialize p-value collection
  all_pvalues <- data.frame(
    plot_name = character(), test_type = character(), comparison = character(),
    n = integer(), statistic = numeric(), statistic_name = character(),
    p_value = numeric(), stringsAsFactors = FALSE
  )

  # Join TE data with chromothripsis status and SV data
  te_chrom_sample <- chromothripsis_status_only %>%
    left_join(te_counts_sample, by = "sample") %>%
    left_join(sv_data, by = "sample") %>%
    filter(!is.na(total_TE))

  cat("\nSamples with TE + chromothripsis + SV data:", nrow(te_chrom_sample), "\n")
  cat("  - with TE = 0:", sum(te_chrom_sample$total_TE == 0), "\n")
  cat("  - with TE > 0:", sum(te_chrom_sample$total_TE > 0), "\n")

  # Create chromothripsis status factor
  te_chrom_sample$chromothripsis_status <- factor(
    ifelse(te_chrom_sample$chromothripsis_overall == 1, "Chromothripsis+", "Chromothripsis-"),
    levels = c("Chromothripsis-", "Chromothripsis+")
  )

  # Create TE-positive subset
  te_positive <- te_chrom_sample %>% filter(total_TE > 0)

  # Save matched samples info
  write.csv(
    data.frame(
      sv_source = sv_source,
      total_samples = nrow(te_chrom_sample),
      te_zero = sum(te_chrom_sample$total_TE == 0),
      te_positive = sum(te_chrom_sample$total_TE > 0),
      chromothripsis_pos = sum(te_chrom_sample$chromothripsis_overall == 1),
      chromothripsis_neg = sum(te_chrom_sample$chromothripsis_overall == 0)
    ),
    file.path(files_dir, "sample_summary.csv"), row.names = FALSE
  )

  # ==========================================================================
  # CHROMOTHRIPSIS STATUS ANALYSIS (same for both sources)
  # ==========================================================================
  cat("\n--- TE vs Chromothripsis Status ---\n")

  # Wilcoxon test - all samples
  wilcox_result <- wilcox.test(total_TE ~ chromothripsis_status, data = te_chrom_sample)
  y_max <- max(te_chrom_sample$total_TE, na.rm = TRUE)

  all_pvalues <- rbind(all_pvalues, data.frame(
    plot_name = "TE vs Chromothripsis (Wilcoxon)", test_type = "Wilcoxon",
    comparison = "Chromothripsis+ vs Chromothripsis-", n = nrow(te_chrom_sample),
    statistic = NA, statistic_name = NA, p_value = wilcox_result$p.value
  ))

  p1 <- ggplot(te_chrom_sample, aes(x = chromothripsis_status, y = total_TE, fill = chromothripsis_status)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.5, size = 1) +
    scale_fill_manual(values = c("Chromothripsis-" = colours[1], "Chromothripsis+" = colours[2])) +
    scale_y_continuous(trans = "log1p") +
    geom_signif(comparisons = list(c("Chromothripsis-", "Chromothripsis+")),
                annotations = format_p(wilcox_result$p.value),
                y_position = log1p(y_max) * 1.1, tip_length = 0.02, textsize = 4) +
    labs(x = "Chromothripsis Status", y = "Total TE Count (log scale)",
         title = paste0("TE Burden by Chromothripsis Status (", sv_source, ")")) +
    theme(legend.position = "none")

  ggsave(file.path(base_dir, "te_vs_chromothripsis_overall_wilcox.png"), plot = p1, width = 6, height = 6)
  titled_print(p1, "TE vs Chromothripsis (Wilcoxon)")

  # LM test - all samples
  lm_result <- lm(total_TE ~ chromothripsis_status, data = te_chrom_sample)
  lm_p <- summary(lm_result)$coefficients[2, 4]

  all_pvalues <- rbind(all_pvalues, data.frame(
    plot_name = "TE vs Chromothripsis (LM)", test_type = "LM",
    comparison = "Chromothripsis+ vs Chromothripsis-", n = nrow(te_chrom_sample),
    statistic = summary(lm_result)$r.squared, statistic_name = "R2", p_value = lm_p
  ))

  p2 <- ggplot(te_chrom_sample, aes(x = chromothripsis_status, y = total_TE, fill = chromothripsis_status)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.5, size = 1) +
    scale_fill_manual(values = c("Chromothripsis-" = colours[1], "Chromothripsis+" = colours[2])) +
    scale_y_continuous(trans = "log1p") +
    geom_signif(comparisons = list(c("Chromothripsis-", "Chromothripsis+")),
                annotations = format_p(lm_p),
                y_position = log1p(y_max) * 1.1, tip_length = 0.02, textsize = 4) +
    labs(x = "Chromothripsis Status", y = "Total TE Count (log scale)",
         title = paste0("TE Burden by Chromothripsis Status - LM (", sv_source, ")")) +
    theme(legend.position = "none")

  ggsave(file.path(base_dir, "te_vs_chromothripsis_overall_lm.png"), plot = p2, width = 6, height = 6)
  titled_print(p2, "TE vs Chromothripsis (LM)")

  # TE-positive only - Wilcoxon
  if (nrow(te_positive) >= 10) {
    wilcox_pos <- wilcox.test(total_TE ~ chromothripsis_status, data = te_positive)
    y_max_pos <- max(te_positive$total_TE, na.rm = TRUE)

    all_pvalues <- rbind(all_pvalues, data.frame(
      plot_name = "TE-Positive: TE vs Chromothripsis (Wilcoxon)", test_type = "Wilcoxon",
      comparison = "TE>0: Chromothripsis+ vs Chromothripsis-", n = nrow(te_positive),
      statistic = NA, statistic_name = NA, p_value = wilcox_pos$p.value
    ))

    p3 <- ggplot(te_positive, aes(x = chromothripsis_status, y = total_TE, fill = chromothripsis_status)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(width = 0.2, alpha = 0.5, size = 1) +
      scale_fill_manual(values = c("Chromothripsis-" = colours[1], "Chromothripsis+" = colours[2])) +
      scale_y_continuous(trans = "log1p") +
      geom_signif(comparisons = list(c("Chromothripsis-", "Chromothripsis+")),
                  annotations = format_p(wilcox_pos$p.value),
                  y_position = log1p(y_max_pos) * 1.1, tip_length = 0.02, textsize = 4) +
      labs(x = "Chromothripsis Status", y = "Total TE Count (log scale)",
           title = paste0("TE>0 Only (n=", nrow(te_positive), ") - Wilcoxon (", sv_source, ")")) +
      theme(legend.position = "none")

    ggsave(file.path(base_dir, "te_vs_chromothripsis_overall_wilcox_gt0.png"), plot = p3, width = 6, height = 6)
    titled_print(p3, "TE-Positive: Wilcoxon")

    # TE-positive only - LM
    lm_pos <- lm(total_TE ~ chromothripsis_status, data = te_positive)
    lm_pos_p <- summary(lm_pos)$coefficients[2, 4]

    all_pvalues <- rbind(all_pvalues, data.frame(
      plot_name = "TE-Positive: TE vs Chromothripsis (LM)", test_type = "LM",
      comparison = "TE>0: Chromothripsis+ vs Chromothripsis-", n = nrow(te_positive),
      statistic = summary(lm_pos)$r.squared, statistic_name = "R2", p_value = lm_pos_p
    ))

    p4 <- ggplot(te_positive, aes(x = chromothripsis_status, y = total_TE, fill = chromothripsis_status)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(width = 0.2, alpha = 0.5, size = 1) +
      scale_fill_manual(values = c("Chromothripsis-" = colours[1], "Chromothripsis+" = colours[2])) +
      scale_y_continuous(trans = "log1p") +
      geom_signif(comparisons = list(c("Chromothripsis-", "Chromothripsis+")),
                  annotations = format_p(lm_pos_p),
                  y_position = log1p(y_max_pos) * 1.1, tip_length = 0.02, textsize = 4) +
      labs(x = "Chromothripsis Status", y = "Total TE Count (log scale)",
           title = paste0("TE>0 Only (n=", nrow(te_positive), ") - LM (", sv_source, ")")) +
      theme(legend.position = "none")

    ggsave(file.path(base_dir, "te_vs_chromothripsis_overall_lm_gt0.png"), plot = p4, width = 6, height = 6)
    titled_print(p4, "TE-Positive: LM")
  }

  # ==========================================================================
  # SV METRIC SCATTER PLOTS - All samples AND TE>0
  # ==========================================================================
  cat("\n--- TE vs SV Metrics ---\n")

  for (metric in sv_metrics) {
    cat("  Processing:", metric, "\n")

    if (!metric %in% colnames(te_chrom_sample)) {
      cat("    Skipping - metric not found\n")
      next
    }

    # --- ALL SAMPLES scatter plot ---
    tryCatch({
      plot_data_all <- te_chrom_sample %>%
        filter(!is.na(.data[[metric]]))

      if (nrow(plot_data_all) >= 10) {
        cor_all <- cor.test(plot_data_all$total_TE, plot_data_all[[metric]],
                            method = "spearman", exact = FALSE)

        all_pvalues <- rbind(all_pvalues, data.frame(
          plot_name = paste("TE vs", metric, "(All)"), test_type = "Spearman",
          comparison = paste("All samples: TE vs", metric), n = nrow(plot_data_all),
          statistic = cor_all$estimate, statistic_name = "rho", p_value = cor_all$p.value
        ))

        x_max <- max(plot_data_all[[metric]], na.rm = TRUE)
        y_max <- max(plot_data_all$total_TE, na.rm = TRUE)

        p_all <- ggplot(plot_data_all, aes(x = .data[[metric]], y = total_TE)) +
          geom_point(alpha = 0.6, size = 2, color = colours[2]) +
          geom_smooth(method = "lm", se = TRUE, color = "red") +
          scale_x_continuous(trans = "log1p") +
          scale_y_continuous(trans = "log1p") +
          annotate("text", x = x_max * 0.5, y = y_max * 0.8,
                   label = paste0("rho = ", round(cor_all$estimate, 3), "\n",
                                  format_p(cor_all$p.value), "\nn = ", nrow(plot_data_all)),
                   hjust = 0, size = 3.5) +
          labs(x = paste(metric, "(log scale)"), y = "Total TE Count (log scale)",
               title = paste0("TE vs ", metric, " - All Samples (", sv_source, ")"))

        ggsave(file.path(base_dir, paste0("te_vs_", metric, "_scatter.png")),
               plot = p_all, width = 7, height = 6)
        titled_print(p_all, paste("TE vs", metric, "- All"))
      }
    }, error = function(e) cat("    Error (all):", e$message, "\n"))

    # --- TE>0 ONLY scatter plot ---
    tryCatch({
      plot_data_gt0 <- te_positive %>%
        filter(!is.na(.data[[metric]]) & .data[[metric]] > 0)

      if (nrow(plot_data_gt0) >= 10) {
        cor_gt0 <- cor.test(plot_data_gt0$total_TE, plot_data_gt0[[metric]],
                            method = "spearman", exact = FALSE)

        all_pvalues <- rbind(all_pvalues, data.frame(
          plot_name = paste("TE vs", metric, "(TE>0)"), test_type = "Spearman",
          comparison = paste("TE>0: TE vs", metric), n = nrow(plot_data_gt0),
          statistic = cor_gt0$estimate, statistic_name = "rho", p_value = cor_gt0$p.value
        ))

        x_max <- max(plot_data_gt0[[metric]], na.rm = TRUE)
        y_max <- max(plot_data_gt0$total_TE, na.rm = TRUE)

        p_gt0 <- ggplot(plot_data_gt0, aes(x = .data[[metric]], y = total_TE)) +
          geom_point(alpha = 0.6, size = 2, color = colours[2]) +
          geom_smooth(method = "lm", se = TRUE, color = "red") +
          scale_x_continuous(trans = "log1p") +
          scale_y_continuous(trans = "log1p") +
          annotate("text", x = x_max * 0.5, y = y_max * 0.8,
                   label = paste0("rho = ", round(cor_gt0$estimate, 3), "\n",
                                  format_p(cor_gt0$p.value), "\nn = ", nrow(plot_data_gt0)),
                   hjust = 0, size = 3.5) +
          labs(x = paste(metric, "(log scale)"), y = "Total TE Count (log scale)",
               title = paste0("TE vs ", metric, " - TE>0 Only (", sv_source, ")"))

        ggsave(file.path(base_dir, paste0("te_vs_", metric, "_scatter_gt0.png")),
               plot = p_gt0, width = 7, height = 6)
        titled_print(p_gt0, paste("TE vs", metric, "- TE>0"))
      }
    }, error = function(e) cat("    Error (gt0):", e$message, "\n"))
  }

  # ==========================================================================
  # TE TYPE ANALYSIS (LINE1, ALU, SVA) - Graphs only, not included in all_pvalues_summary
  # ==========================================================================
  cat("\n--- TE Type Analysis (graphs only, not in summary) ---\n")

  te_types <- c("LINE1", "ALU", "SVA")

  for (te_type in te_types) {
    count_col <- paste0(te_type, "_count")

    if (!count_col %in% colnames(te_chrom_sample)) next
    if (sum(te_chrom_sample[[count_col]], na.rm = TRUE) == 0) next

    cat("  Processing:", te_type, "\n")

    # TE type vs chromothripsis - graph only
    tryCatch({
      wilcox_type <- wilcox.test(te_chrom_sample[[count_col]] ~ te_chrom_sample$chromothripsis_status)
      y_max_type <- max(te_chrom_sample[[count_col]], na.rm = TRUE)

      p_type <- ggplot(te_chrom_sample,
                       aes(x = chromothripsis_status, y = .data[[count_col]], fill = chromothripsis_status)) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(width = 0.2, alpha = 0.5, size = 1) +
        scale_fill_manual(values = c("Chromothripsis-" = colours[1], "Chromothripsis+" = colours[2])) +
        scale_y_continuous(trans = "log1p") +
        geom_signif(comparisons = list(c("Chromothripsis-", "Chromothripsis+")),
                    annotations = format_p(wilcox_type$p.value),
                    y_position = log1p(y_max_type) * 1.1, tip_length = 0.02, textsize = 4) +
        labs(x = "Chromothripsis Status", y = paste(te_type, "Count (log scale)"),
             title = paste0(te_type, " by Chromothripsis Status (", sv_source, ")")) +
        theme(legend.position = "none")

      ggsave(file.path(base_dir, paste0(te_type, "_vs_chromothripsis_wilcox.png")),
             plot = p_type, width = 6, height = 6)
      titled_print(p_type, paste(te_type, "vs Chromothripsis"))

    }, error = function(e) cat("    Error:", e$message, "\n"))

    # TE type vs SV metrics - graph only
    for (metric in sv_metrics) {
      if (!metric %in% colnames(te_chrom_sample)) next

      tryCatch({
        plot_data <- te_chrom_sample %>%
          filter(!is.na(.data[[metric]]) & .data[[metric]] > 0 & .data[[count_col]] > 0)

        if (nrow(plot_data) >= 10) {
          cor_type <- cor.test(plot_data[[count_col]], plot_data[[metric]],
                               method = "spearman", exact = FALSE)

          x_max <- max(plot_data[[metric]], na.rm = TRUE)
          y_max <- max(plot_data[[count_col]], na.rm = TRUE)

          p_type_sv <- ggplot(plot_data, aes(x = .data[[metric]], y = .data[[count_col]])) +
            geom_point(alpha = 0.6, size = 2, color = colours[2]) +
            geom_smooth(method = "lm", se = TRUE, color = "red") +
            scale_x_continuous(trans = "log1p") +
            scale_y_continuous(trans = "log1p") +
            annotate("text", x = x_max * 0.5, y = y_max * 0.8,
                     label = paste0("rho = ", round(cor_type$estimate, 3), "\n",
                                    format_p(cor_type$p.value)),
                     hjust = 0, size = 3.5) +
            labs(x = paste(metric, "(log scale)"), y = paste(te_type, "Count (log scale)"),
                 title = paste0(te_type, " vs ", metric, " (", sv_source, ")"))

          ggsave(file.path(base_dir, paste0(te_type, "_vs_", metric, "_scatter.png")),
                 plot = p_type_sv, width = 7, height = 6)
          titled_print(p_type_sv, paste(te_type, "vs", metric))
        }
      }, error = function(e) NULL)
    }
  }

  # ==========================================================================
  # BY CANCER TYPE ANALYSIS
  # ==========================================================================
  cat("\n--- By Cancer Type Analysis ---\n")

  cancer_type_counts <- te_chrom_sample %>%
    filter(!is.na(tumor_type)) %>%
    group_by(tumor_type) %>%
    summarise(n = n(), .groups = "drop") %>%
    filter(n >= 5)

  cancer_type_correlations <- data.frame(
    cancer_type = character(), sv_metric = character(), n_samples = integer(),
    spearman_rho = numeric(), spearman_p = numeric(), stringsAsFactors = FALSE
  )

  # All cancers combined
  for (metric in sv_metrics) {
    if (!metric %in% colnames(te_chrom_sample)) next

    plot_data <- te_chrom_sample %>%
      filter(!is.na(.data[[metric]]) & .data[[metric]] > 0 & total_TE > 0)

    if (nrow(plot_data) >= 5) {
      cor_result <- cor.test(plot_data$total_TE, plot_data[[metric]],
                             method = "spearman", exact = FALSE)

      cancer_type_correlations <- rbind(cancer_type_correlations, data.frame(
        cancer_type = "All Cancers", sv_metric = metric, n_samples = nrow(plot_data),
        spearman_rho = cor_result$estimate, spearman_p = cor_result$p.value
      ))

      x_max <- max(plot_data[[metric]], na.rm = TRUE)
      y_max <- max(plot_data$total_TE, na.rm = TRUE)

      p <- ggplot(plot_data, aes(x = .data[[metric]], y = total_TE)) +
        geom_point(alpha = 0.6, size = 2, color = colours[2]) +
        geom_smooth(method = "lm", se = TRUE, color = "red") +
        scale_x_continuous(trans = "log1p") +
        scale_y_continuous(trans = "log1p") +
        annotate("text", x = x_max * 0.5, y = y_max * 0.8,
                 label = paste0("rho = ", round(cor_result$estimate, 3), "\n",
                                format_p(cor_result$p.value), "\nn = ", nrow(plot_data)),
                 hjust = 0, size = 3.5) +
        labs(x = paste(metric, "(log scale)"), y = "Total TE Count (log scale)",
             title = paste0("All Cancers: TE vs ", metric, " (", sv_source, ")"))

      ggsave(file.path(cancer_type_dir, paste0("all_cancers_te_vs_", metric, ".png")),
             plot = p, width = 7, height = 6)
      titled_print(p, paste("All Cancers: TE vs", metric))
    }
  }

  # By cancer type
  for (cancer in cancer_type_counts$tumor_type) {
    cancer_data <- te_chrom_sample %>% filter(tumor_type == cancer)
    cancer_clean <- gsub("[^a-zA-Z0-9]", "_", cancer)

    for (metric in sv_metrics) {
      if (!metric %in% colnames(cancer_data)) next

      plot_data <- cancer_data %>%
        filter(!is.na(.data[[metric]]) & .data[[metric]] > 0 & total_TE > 0)

      if (nrow(plot_data) >= 5) {
        cor_result <- cor.test(plot_data$total_TE, plot_data[[metric]],
                               method = "spearman", exact = FALSE)

        cancer_type_correlations <- rbind(cancer_type_correlations, data.frame(
          cancer_type = cancer, sv_metric = metric, n_samples = nrow(plot_data),
          spearman_rho = cor_result$estimate, spearman_p = cor_result$p.value
        ))

        x_max <- max(plot_data[[metric]], na.rm = TRUE)
        y_max <- max(plot_data$total_TE, na.rm = TRUE)

        p <- ggplot(plot_data, aes(x = .data[[metric]], y = total_TE)) +
          geom_point(alpha = 0.6, size = 2, color = colours[2]) +
          geom_smooth(method = "lm", se = TRUE, color = "red") +
          scale_x_continuous(trans = "log1p") +
          scale_y_continuous(trans = "log1p") +
          annotate("text", x = x_max * 0.5, y = y_max * 0.8,
                   label = paste0("rho = ", round(cor_result$estimate, 3), "\n",
                                  format_p(cor_result$p.value), "\nn = ", nrow(plot_data)),
                   hjust = 0, size = 3.5) +
          labs(x = paste(metric, "(log scale)"), y = "Total TE Count (log scale)",
               title = paste0(cancer, ": TE vs ", metric, " (", sv_source, ")"))

        ggsave(file.path(cancer_type_dir, paste0(cancer_clean, "_te_vs_", metric, ".png")),
               plot = p, width = 7, height = 6)
        titled_print(p, paste(cancer, ": TE vs", metric))
      }
    }
  }

  # Save cancer type correlations
  if (nrow(cancer_type_correlations) > 0) {
    cancer_type_correlations$p_adj <- p.adjust(cancer_type_correlations$spearman_p, method = "fdr")
    cancer_type_correlations <- cancer_type_correlations[order(cancer_type_correlations$spearman_p), ]
    write.csv(cancer_type_correlations, file.path(files_dir, "sv_te_correlations_by_cancer_type.csv"), row.names = FALSE)
  }

  # ==========================================================================
  # SAVE P-VALUE SUMMARY
  # ==========================================================================
  all_pvalues$p_adj <- p.adjust(all_pvalues$p_value, method = "fdr")
  all_pvalues <- all_pvalues[order(all_pvalues$p_value), ]
  write.csv(all_pvalues, file.path(files_dir, "all_pvalues_summary.csv"), row.names = FALSE)

  cat("\n=== P-Value Summary ===\n")
  cat("Total tests:", nrow(all_pvalues), "\n")
  cat("Significant at p < 0.05:", sum(all_pvalues$p_value < 0.05, na.rm = TRUE), "\n")
  cat("Significant at FDR < 0.05:", sum(all_pvalues$p_adj < 0.05, na.rm = TRUE), "\n")

  # ==========================================================================
  # CHROMOSOME-LEVEL ANALYSIS
  # ==========================================================================
  cat("\n--- Chromosome-Level Analysis ---\n")

  # Create per_chromosome directory
  per_chrom_dir <- file.path(base_dir, "per_chromosome")
  dir.create(per_chrom_dir, showWarnings = FALSE, recursive = TRUE)

  # Initialize chromosome-level p-value collection
  chrom_pvalues <- data.frame(
    chromosome = character(), analysis_type = character(), comparison = character(),
    test = character(), n = integer(), rho = numeric(), statistic = numeric(),
    p_value = numeric(), stringsAsFactors = FALSE
  )

  # Skip if per-chromosome TE data not available
  if (is.null(te_counts_per_chrom)) {
    cat("  Skipping chromosome-level analysis - per-chromosome TE data not available\n")
  } else {

    # Merge per-chromosome TE counts with chromothripsis data
    chrom_te_merged <- chrom_level_data %>%
      left_join(te_counts_per_chrom, by = c("sample", "chrom")) %>%
      mutate(te_count = replace_na(te_count, 0)) %>%
      left_join(chromothripsis_status_only %>% select(sample, chromothripsis_overall), by = "sample") %>%
      mutate(
        overall_status = factor(
          ifelse(chromothripsis_overall == 1, "Chromothripsis+", "Chromothripsis-"),
          levels = c("Chromothripsis-", "Chromothripsis+")
        )
      )

    cat("  Merged chromosome-level data:", nrow(chrom_te_merged), "rows\n")

    # --------------------------------------------------------------------------
    # A. TE on Chromothriptic vs Non-Chromothriptic Chromosomes
    # --------------------------------------------------------------------------
    cat("\n  --- A. TE on Chromothriptic vs Non-Chromothriptic Chromosomes ---\n")

    # Overall comparison: all chromosome-level observations
    chrom_te_for_test <- chrom_te_merged %>%
      mutate(chrom_status = factor(chrom_chromothripsis_status, levels = c("Chromothripsis-", "Chromothripsis+")))

    if (length(unique(chrom_te_for_test$chrom_status)) == 2) {
      wilcox_chrom <- wilcox.test(te_count ~ chrom_status, data = chrom_te_for_test)

      chrom_pvalues <- rbind(chrom_pvalues, data.frame(
        chromosome = "ALL", analysis_type = "TE_vs_chrom_chromothripsis",
        comparison = "TE on Chromothripsis+ vs Chromothripsis- chromosomes",
        test = "Wilcoxon", n = nrow(chrom_te_for_test), rho = NA,
        statistic = wilcox_chrom$statistic, p_value = wilcox_chrom$p.value
      ))

      # Box plot
      y_max_chrom <- max(chrom_te_for_test$te_count, na.rm = TRUE)
      p_chrom_box <- ggplot(chrom_te_for_test, aes(x = chrom_status, y = te_count, fill = chrom_status)) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(width = 0.2, alpha = 0.3, size = 0.5) +
        scale_fill_manual(values = c("Chromothripsis-" = colours[1], "Chromothripsis+" = colours[2])) +
        scale_y_continuous(trans = "log1p") +
        geom_signif(comparisons = list(c("Chromothripsis-", "Chromothripsis+")),
                    annotations = format_p(wilcox_chrom$p.value),
                    y_position = log1p(y_max_chrom) * 1.1, tip_length = 0.02, textsize = 4) +
        labs(x = "Chromosome Chromothripsis Status", y = "TE Count per Chromosome (log scale)",
             title = paste0("TE on Chromothriptic vs Non-Chromothriptic Chromosomes (", sv_source, ")"),
             subtitle = paste0("n = ", nrow(chrom_te_for_test), " chromosome-sample observations")) +
        theme(legend.position = "none")

      ggsave(file.path(base_dir, "te_vs_chrom_status_boxplot.png"), plot = p_chrom_box, width = 7, height = 6)
      titled_print(p_chrom_box, "TE on Chrom+ vs Chrom- chromosomes")

      cat("    Wilcoxon p =", format_p(wilcox_chrom$p.value), "\n")
    }

    # Stratify by overall sample chromothripsis status
    for (overall_stat in c("Chromothripsis+", "Chromothripsis-")) {
      subset_data <- chrom_te_for_test %>% filter(overall_status == overall_stat)

      if (length(unique(subset_data$chrom_status)) == 2 && nrow(subset_data) >= 20) {
        wilcox_strat <- wilcox.test(te_count ~ chrom_status, data = subset_data)

        chrom_pvalues <- rbind(chrom_pvalues, data.frame(
          chromosome = "ALL", analysis_type = "TE_vs_chrom_chromothripsis_stratified",
          comparison = paste0("TE on Chrom+ vs Chrom- chromosomes (", overall_stat, " samples)"),
          test = "Wilcoxon", n = nrow(subset_data), rho = NA,
          statistic = wilcox_strat$statistic, p_value = wilcox_strat$p.value
        ))

        cat("    Stratified (", overall_stat, " samples): p =", format_p(wilcox_strat$p.value), "\n")
      }
    }

    # --------------------------------------------------------------------------
    # B. Genome-Wide Heatmaps
    # --------------------------------------------------------------------------
    cat("\n  --- B. Genome-Wide Heatmaps ---\n")

    # Create sample x chromosome matrix for TE counts
    te_matrix <- chrom_te_merged %>%
      select(sample, chrom, te_count) %>%
      pivot_wider(names_from = chrom, values_from = te_count, values_fill = 0) %>%
      column_to_rownames("sample")

    # Reorder columns to chromosome order
    chrom_order <- c(as.character(1:22), "X")
    te_matrix <- te_matrix[, intersect(chrom_order, colnames(te_matrix))]

    # Create sample x chromosome matrix for chromothripsis status
    chrom_status_matrix <- chrom_te_merged %>%
      select(sample, chrom, chrom_chromothripsis_level) %>%
      pivot_wider(names_from = chrom, values_from = chrom_chromothripsis_level, values_fill = "None") %>%
      column_to_rownames("sample")
    chrom_status_matrix <- chrom_status_matrix[, intersect(chrom_order, colnames(chrom_status_matrix))]

    # Get sample annotations
    sample_annotations <- chrom_te_merged %>%
      select(sample, overall_status) %>%
      distinct() %>%
      left_join(te_counts_sample %>% select(sample, total_TE, tumor_type), by = "sample")

    if (nrow(te_matrix) >= 5) {
      tryCatch({
        # Convert chromothripsis status to numeric for clustering (None=0, Low=1, High=2)
        chrom_status_numeric <- chrom_status_matrix
        chrom_status_numeric[chrom_status_numeric == "None"] <- 0
        chrom_status_numeric[chrom_status_numeric == "Low"] <- 1
        chrom_status_numeric[chrom_status_numeric == "High"] <- 2
        chrom_status_numeric <- apply(chrom_status_numeric, 2, as.numeric)
        rownames(chrom_status_numeric) <- rownames(chrom_status_matrix)

        # Cluster based on chromothripsis pattern and get row order
        chrom_dist <- dist(chrom_status_numeric, method = "euclidean")
        chrom_clust <- hclust(chrom_dist, method = "ward.D2")
        row_order_by_chrom <- chrom_clust$order

        # Row annotation
        row_anno <- rowAnnotation(
          Overall_Status = sample_annotations$overall_status[match(rownames(te_matrix), sample_annotations$sample)],
          Total_TE = anno_barplot(sample_annotations$total_TE[match(rownames(te_matrix), sample_annotations$sample)],
                                  width = unit(2, "cm")),
          col = list(Overall_Status = c("Chromothripsis-" = colours[1], "Chromothripsis+" = colours[2]))
        )

        # TE count heatmap - use row order from chromothripsis clustering
        ht_te <- Heatmap(
          as.matrix(te_matrix),
          name = "TE Count",
          col = colorRamp2(c(0, 1, max(te_matrix, na.rm = TRUE)), c("white", "lightblue", "darkblue")),
          cluster_columns = FALSE,
          cluster_rows = FALSE,
          row_order = row_order_by_chrom,
          show_row_names = FALSE,
          column_title = paste0("TE Count per Chromosome (", sv_source, ")"),
          right_annotation = row_anno
        )

        png(file.path(base_dir, "genome_heatmap_te.png"), width = 12, height = 10, units = "in", res = 150)
        draw(ht_te)
        dev.off()

        cat("    TE count heatmap saved (ordered by chromothripsis pattern)\n")

        # Chromothripsis status heatmap - same row order
        status_colors <- c("None" = "white", "Low" = "orange", "High" = "red")

        ht_chrom <- Heatmap(
          as.matrix(chrom_status_matrix),
          name = "Chromothripsis",
          col = status_colors,
          cluster_columns = FALSE,
          cluster_rows = FALSE,
          row_order = row_order_by_chrom,
          show_row_names = FALSE,
          column_title = paste0("Chromothripsis Status per Chromosome (", sv_source, ")"),
          right_annotation = row_anno
        )

        png(file.path(base_dir, "genome_heatmap_chromothripsis.png"), width = 12, height = 10, units = "in", res = 150)
        draw(ht_chrom)
        dev.off()

        cat("    Chromothripsis status heatmap saved\n")

      }, error = function(e) cat("    Heatmap error:", e$message, "\n"))
    }

    # --------------------------------------------------------------------------
    # C. Per-Chromosome Correlations
    # --------------------------------------------------------------------------
    cat("\n  --- C. Per-Chromosome Correlations ---\n")

    # Per-chromosome SV metrics available in chrom_level_data
    per_chrom_sv_metrics <- c("number_DEL", "number_DUP", "number_h2hINV", "number_t2tINV",
                               "number_TRA", "number_CNV_segments")

    for (chr in chromosomes) {
      chr_dir <- file.path(per_chrom_dir, paste0("chr", chr))
      dir.create(chr_dir, showWarnings = FALSE, recursive = TRUE)

      chr_data <- chrom_te_merged %>% filter(chrom == chr)

      if (nrow(chr_data) < 10) {
        cat("    chr", chr, ": Skipping - insufficient samples (n=", nrow(chr_data), ")\n")
        next
      }

      chr_results <- data.frame()

      # --- TE vs per-chromosome SV metrics ---
      chr_scatter_plots <- list()

      for (metric in per_chrom_sv_metrics) {
        if (!metric %in% colnames(chr_data)) next

        # Option 2: Filter to samples with BOTH TE > 0 AND SV metric > 0
        # This ensures meaningful correlations by excluding zero-zero pairs
        metric_data <- chr_data %>%
          filter(!is.na(.data[[metric]]) & .data[[metric]] > 0 & te_count > 0)

        # Clean metric name for display (remove "number_" prefix)
        metric_clean <- gsub("number_", "", metric)

        # Only run correlation if we have at least 5 samples with both TE > 0 and SV > 0
        if (nrow(metric_data) >= 5 && var(metric_data$te_count, na.rm = TRUE) > 0 &&
            var(metric_data[[metric]], na.rm = TRUE) > 0) {
          cor_result <- cor.test(metric_data$te_count, metric_data[[metric]],
                                 method = "spearman", exact = FALSE)

          chrom_pvalues <- rbind(chrom_pvalues, data.frame(
            chromosome = chr, analysis_type = "TE_vs_SV_filtered",
            comparison = paste("TE vs", metric_clean, "(filtered)"), test = "Spearman",
            n = nrow(metric_data), rho = cor_result$estimate,
            statistic = cor_result$statistic, p_value = cor_result$p.value
          ))

          chr_results <- rbind(chr_results, data.frame(
            comparison = paste("TE vs", metric_clean, "(filtered)"), test = "Spearman",
            n = nrow(metric_data), rho = cor_result$estimate, p_value = cor_result$p.value
          ))

          # Create scatter plot for this metric
          x_max <- max(metric_data[[metric]], na.rm = TRUE)
          y_max <- max(metric_data$te_count, na.rm = TRUE)

          p_scatter <- ggplot(metric_data, aes(x = .data[[metric]], y = te_count)) +
            geom_point(alpha = 0.6, size = 2, color = colours[2]) +
            geom_smooth(method = "lm", se = TRUE, color = "red") +
            scale_x_continuous(trans = "log1p") +
            scale_y_continuous(trans = "log1p") +
            annotate("text", x = x_max * 0.5, y = y_max * 0.8,
                     label = paste0("rho = ", round(cor_result$estimate, 3), "\n",
                                    format_p(cor_result$p.value), "\nn = ", nrow(metric_data)),
                     hjust = 0, size = 3.5) +
            labs(x = paste(metric_clean, "(log scale)"), y = "TE Count (log scale)",
                 title = paste0("Chr", chr, ": TE vs ", metric_clean))

          chr_scatter_plots[[metric_clean]] <- p_scatter

          # Save individual scatter plot (filtered)
          ggsave(file.path(chr_dir, paste0("chr", chr, "_te_vs_", metric_clean, "_filtered.png")),
                 plot = p_scatter, width = 6, height = 5)
        }

        # Unfiltered Spearman correlation (all samples, no TE>0 or SV>0 filter)
        metric_data_unfiltered <- chr_data %>% filter(!is.na(.data[[metric]]))

        if (nrow(metric_data_unfiltered) >= 10 && var(metric_data_unfiltered$te_count, na.rm = TRUE) > 0 &&
            var(metric_data_unfiltered[[metric]], na.rm = TRUE) > 0) {
          cor_result_unfiltered <- cor.test(metric_data_unfiltered$te_count, metric_data_unfiltered[[metric]],
                                            method = "spearman", exact = FALSE)

          chrom_pvalues <- rbind(chrom_pvalues, data.frame(
            chromosome = chr, analysis_type = "TE_vs_SV_unfiltered",
            comparison = paste("TE vs", metric_clean, "(unfiltered)"), test = "Spearman",
            n = nrow(metric_data_unfiltered), rho = cor_result_unfiltered$estimate,
            statistic = cor_result_unfiltered$statistic, p_value = cor_result_unfiltered$p.value
          ))

          chr_results <- rbind(chr_results, data.frame(
            comparison = paste("TE vs", metric_clean, "(unfiltered)"), test = "Spearman",
            n = nrow(metric_data_unfiltered), rho = cor_result_unfiltered$estimate, p_value = cor_result_unfiltered$p.value
          ))

          # Create scatter plot for unfiltered data
          x_max_uf <- max(metric_data_unfiltered[[metric]], na.rm = TRUE)
          y_max_uf <- max(metric_data_unfiltered$te_count, na.rm = TRUE)

          p_scatter_unfiltered <- ggplot(metric_data_unfiltered, aes(x = .data[[metric]], y = te_count)) +
            geom_point(alpha = 0.6, size = 2, color = colours[3]) +
            geom_smooth(method = "lm", se = TRUE, color = "red") +
            scale_x_continuous(trans = "log1p") +
            scale_y_continuous(trans = "log1p") +
            annotate("text", x = x_max_uf * 0.5, y = y_max_uf * 0.8,
                     label = paste0("rho = ", round(cor_result_unfiltered$estimate, 3), "\n",
                                    format_p(cor_result_unfiltered$p.value), "\nn = ", nrow(metric_data_unfiltered)),
                     hjust = 0, size = 3.5) +
            labs(x = paste(metric_clean, "(log scale)"), y = "TE Count (log scale)",
                 title = paste0("Chr", chr, ": TE vs ", metric_clean, " (unfiltered)"))

          # Save unfiltered scatter plot
          ggsave(file.path(chr_dir, paste0("chr", chr, "_te_vs_", metric_clean, "_unfiltered.png")),
                 plot = p_scatter_unfiltered, width = 6, height = 5)
        }

        # Option 4: Fisher's exact test for binary analysis
        # Test: "Do samples with SV > 0 on this chromosome have higher rate of TE > 0?"
        # This works even with sparse data where many samples have zeros
        fisher_data <- chr_data %>%
          mutate(
            has_te = te_count > 0,
            has_sv = .data[[metric]] > 0
          ) %>%
          filter(!is.na(has_sv))

        if (nrow(fisher_data) >= 10) {
          # Create 2x2 contingency table: rows = has_sv, cols = has_te
          contingency <- table(fisher_data$has_sv, fisher_data$has_te)

          # Need at least some variation in both dimensions
          if (nrow(contingency) == 2 && ncol(contingency) == 2) {
            fisher_result <- fisher.test(contingency)

            # Calculate proportions for interpretation
            sv_pos_te_rate <- sum(fisher_data$has_sv & fisher_data$has_te) / sum(fisher_data$has_sv)
            sv_neg_te_rate <- sum(!fisher_data$has_sv & fisher_data$has_te) / sum(!fisher_data$has_sv)

            chrom_pvalues <- rbind(chrom_pvalues, data.frame(
              chromosome = chr, analysis_type = "TE_vs_SV_binary",
              comparison = paste("TE presence vs", metric_clean, "presence"),
              test = "Fisher",
              n = nrow(fisher_data),
              rho = fisher_result$estimate,  # odds ratio
              statistic = NA,
              p_value = fisher_result$p.value
            ))

            chr_results <- rbind(chr_results, data.frame(
              comparison = paste("TE vs", metric_clean, "(binary)"),
              test = "Fisher",
              n = nrow(fisher_data),
              rho = fisher_result$estimate,  # odds ratio
              p_value = fisher_result$p.value
            ))
          }
        }
      }

      # --- TE vs chromosome-specific chromothripsis status ---
      if (length(unique(chr_data$chrom_chromothripsis_status)) == 2) {
        wilcox_chr <- wilcox.test(te_count ~ chrom_chromothripsis_status, data = chr_data)

        chrom_pvalues <- rbind(chrom_pvalues, data.frame(
          chromosome = chr, analysis_type = "TE_vs_chrom_chromothripsis",
          comparison = paste0("TE by chr", chr, " chromothripsis status"),
          test = "Wilcoxon", n = nrow(chr_data), rho = NA,
          statistic = wilcox_chr$statistic, p_value = wilcox_chr$p.value
        ))

        chr_results <- rbind(chr_results, data.frame(
          comparison = paste0("TE by chr", chr, " status"), test = "Wilcoxon",
          n = nrow(chr_data), rho = NA, p_value = wilcox_chr$p.value
        ))

        # Box plot for chromosome-specific chromothripsis
        y_max_chr <- max(chr_data$te_count, na.rm = TRUE)
        p_box_chr <- ggplot(chr_data, aes(x = chrom_chromothripsis_status, y = te_count,
                                           fill = chrom_chromothripsis_status)) +
          geom_boxplot(outlier.shape = NA) +
          geom_jitter(width = 0.2, alpha = 0.5, size = 1) +
          scale_fill_manual(values = c("Chromothripsis-" = colours[1], "Chromothripsis+" = colours[2])) +
          scale_y_continuous(trans = "log1p") +
          geom_signif(comparisons = list(c("Chromothripsis-", "Chromothripsis+")),
                      annotations = format_p(wilcox_chr$p.value),
                      y_position = log1p(y_max_chr) * 1.1, tip_length = 0.02, textsize = 4) +
          labs(x = paste0("Chr", chr, " Chromothripsis Status"), y = "TE Count (log scale)",
               title = paste0("Chr", chr, ": TE by Chromosome Chromothripsis Status")) +
          theme(legend.position = "none")

        ggsave(file.path(chr_dir, paste0("chr", chr, "_te_vs_chrom_chromothripsis.png")),
               plot = p_box_chr, width = 6, height = 5)
      }

      # --- TE on this chromosome vs overall chromothripsis status ---
      if (length(unique(chr_data$overall_status)) == 2) {
        wilcox_overall <- wilcox.test(te_count ~ overall_status, data = chr_data)

        chrom_pvalues <- rbind(chrom_pvalues, data.frame(
          chromosome = chr, analysis_type = "TE_vs_overall_chromothripsis",
          comparison = paste0("TE on chr", chr, " by overall chromothripsis"),
          test = "Wilcoxon", n = nrow(chr_data), rho = NA,
          statistic = wilcox_overall$statistic, p_value = wilcox_overall$p.value
        ))

        chr_results <- rbind(chr_results, data.frame(
          comparison = paste0("TE on chr", chr, " by overall status"), test = "Wilcoxon",
          n = nrow(chr_data), rho = NA, p_value = wilcox_overall$p.value
        ))

        # Box plot for overall chromothripsis status
        y_max_overall <- max(chr_data$te_count, na.rm = TRUE)
        p_box_overall <- ggplot(chr_data, aes(x = overall_status, y = te_count, fill = overall_status)) +
          geom_boxplot(outlier.shape = NA) +
          geom_jitter(width = 0.2, alpha = 0.5, size = 1) +
          scale_fill_manual(values = c("Chromothripsis-" = colours[1], "Chromothripsis+" = colours[2])) +
          scale_y_continuous(trans = "log1p") +
          geom_signif(comparisons = list(c("Chromothripsis-", "Chromothripsis+")),
                      annotations = format_p(wilcox_overall$p.value),
                      y_position = log1p(y_max_overall) * 1.1, tip_length = 0.02, textsize = 4) +
          labs(x = "Overall Chromothripsis Status", y = "TE Count (log scale)",
               title = paste0("Chr", chr, ": TE by Overall Chromothripsis Status")) +
          theme(legend.position = "none")

        ggsave(file.path(chr_dir, paste0("chr", chr, "_te_vs_overall_chromothripsis.png")),
               plot = p_box_overall, width = 6, height = 5)
      }

      # Save per-chromosome results CSV
      if (nrow(chr_results) > 0) {
        write.csv(chr_results, file.path(chr_dir, paste0("chr", chr, "_te_vs_sv_correlations.csv")), row.names = FALSE)
      }

      cat("    chr", chr, ": ", nrow(chr_results), " tests completed\n")
    }

    # --------------------------------------------------------------------------
    # Summary Visualizations
    # --------------------------------------------------------------------------
    cat("\n  --- Summary Visualizations ---\n")

    if (nrow(chrom_pvalues) > 0) {
      # Apply BH correction
      chrom_pvalues$p_adj <- p.adjust(chrom_pvalues$p_value, method = "fdr")

      # Save master file
      write.csv(chrom_pvalues, file.path(files_dir, "chromosome_pvalues_correlations.csv"), row.names = FALSE)
      cat("    Saved chromosome_pvalues_correlations.csv with", nrow(chrom_pvalues), "tests\n")

      # Per-chromosome correlation heatmap (TE vs SV metrics)
      sv_correlations <- chrom_pvalues %>%
        filter(analysis_type == "TE_vs_SV") %>%
        select(chromosome, comparison, rho) %>%
        mutate(sv_metric = gsub("TE vs ", "", comparison)) %>%
        select(chromosome, sv_metric, rho)

      if (nrow(sv_correlations) > 0) {
        corr_matrix <- sv_correlations %>%
          pivot_wider(names_from = sv_metric, values_from = rho)

        if (nrow(corr_matrix) > 1) {
          corr_matrix_vals <- corr_matrix %>%
            column_to_rownames("chromosome") %>%
            as.matrix()

          # Reorder rows to chromosome order
          row_order <- intersect(c(as.character(1:22), "X"), rownames(corr_matrix_vals))
          corr_matrix_vals <- corr_matrix_vals[row_order, , drop = FALSE]

          tryCatch({
            ht_corr <- Heatmap(
              corr_matrix_vals,
              name = "Spearman rho",
              col = colorRamp2(c(-0.5, 0, 0.5), c("blue", "white", "red")),
              cluster_columns = TRUE,
              cluster_rows = FALSE,
              column_title = paste0("TE-SV Correlations by Chromosome (", sv_source, ")"),
              row_title = "Chromosome",
              cell_fun = function(j, i, x, y, width, height, fill) {
                val <- corr_matrix_vals[i, j]
                if (!is.na(val)) {
                  grid.text(sprintf("%.2f", val), x, y, gp = gpar(fontsize = 8))
                }
              }
            )

            png(file.path(base_dir, "per_chrom_correlation_heatmap.png"), width = 10, height = 12, units = "in", res = 150)
            draw(ht_corr)
            dev.off()

            cat("    Per-chromosome correlation heatmap saved\n")
          }, error = function(e) cat("    Correlation heatmap error:", e$message, "\n"))
        }
      }

      # Summary: significant results
      sig_results <- chrom_pvalues %>% filter(p_adj < 0.05)
      cat("\n    Significant results (FDR < 0.05):", nrow(sig_results), "\n")
      if (nrow(sig_results) > 0) {
        print(sig_results %>% arrange(p_adj) %>% head(10))
      }
    }
  }

  # Close PDF
  dev.off()

  cat("\nCompleted analysis for", sv_source, "\n")
  cat("Output directory:", base_dir, "\n")

  return(invisible(NULL))
}

# ============================================================================
# RUN ANALYSIS FOR BOTH SV SOURCES
# ============================================================================

# Run ShatterSeek analysis
run_chromothripsis_analysis(
  sv_source = "shatterseek",
  sv_data = shatterseek_sv,
  sv_metrics = shatterseek_metrics,
  output_subdir = "shatterseek"
)

# Run GRIDSS analysis
run_chromothripsis_analysis(
  sv_source = "gridss",
  sv_data = gridss_sv,
  sv_metrics = gridss_metrics,
  output_subdir = "gridss"
)

cat("\n")
cat("========================================\n")
cat("  CHROMOTHRIPSIS ANALYSIS COMPLETE\n")
cat("========================================\n")
cat("ShatterSeek output:", file.path(plot_dir, "chromothripsis/shatterseek"), "\n")
cat("GRIDSS output:", file.path(plot_dir, "chromothripsis/gridss"), "\n")

cat("\nCompleted 02_te_viz_tumour_13_chromothripsis.R\n")
