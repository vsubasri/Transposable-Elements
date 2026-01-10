#!/usr/bin/env Rscript

# Tumour TE Visualization - LOH Timing Analysis
# Comprehensive TP53 LOH timing analysis with descriptive stats, scatter plots,
# boxplots, and linear models

# Source common setup and load data
source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/00_viz_common_setup.R")
REQUIRED_DATA <- c("count_matrix", "loh", "clinical")
source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/00_viz_load_data_tumour.R")

# Initialize module-specific text output
init_module_sink(paste0(plot_dir, "loh/"), "LOH")

cat("Running 02_te_viz_tumour_04_loh.R...\n")

# =============================================================================
# DATA PREPARATION
# =============================================================================

cat("\n===== DATA PREPARATION =====\n")

# Define color schemes
colours <- c("#5FBFF9", "#235789")
colours_3 <- c("#DDD8C4", "#5FBFF9", "#235789")
colours_5 <- brewer.pal(5, "Set2")
colours_8 <- brewer.pal(8, "Set2")
lfs_colour <- "#5FBFF9"

# Prepare LOH time data - standardize sample names
tryCatch({
  # Function to transform KiCS sample names
  transform_kics_sample <- function(sample) {
    if (grepl("^KiCS", sample)) {
      # Remove KiCS prefix
      s <- sub("^KiCS", "", sample)
      # Split on first underscore to get ID and rest
      parts <- strsplit(s, "_", fixed = TRUE)[[1]]
      id <- parts[1]
      rest <- paste(parts[-1], collapse = "-")
      paste0(id, "_", rest, "_T")
    } else {
      sample
    }
  }

  loh_time_mod <- loh_time %>%
    select(sample, time) %>%
    mutate(sample = sapply(sample, transform_kics_sample, USE.NAMES = FALSE)) %>%
    # Clean near-zero time values (floating-point artifacts)
    mutate(time = ifelse(time < 1e-10, 0, time))

  cat("\n=== LOH sample name transformation ===\n")
  cat("Original:", paste(head(loh_time$sample, 10), collapse=", "), "\n")
  cat("Modified:", paste(head(loh_time_mod$sample, 10), collapse=", "), "\n")

  # Merge LOH time with TE data to create data_loh
  data_loh <- merge(te_all_all_t, loh_time_mod, by = "sample")

  # Use te_all_all_t as data_clinical
  data_clinical <- te_all_all_t

  cat("LOH data merged successfully\n")
  cat("  - data_loh samples:", nrow(data_loh), "\n")
  cat("  - data_clinical samples:", nrow(data_clinical), "\n")

  # Show available columns
  cat("\nAvailable columns in data_loh:\n")
  cat(paste(colnames(data_loh), collapse = ", "), "\n")

  write_output(quote(head(data_loh)), "Head of merged LOH data")

  # Export merged data
  write.table(data_loh, file = paste0(r_dir_files, "te_loh.csv"),
              quote = FALSE, sep = ',', row.names = FALSE, col.names = TRUE)

}, error = function(e) {
  cat("Warning: Could not prepare LOH data:", e$message, "\n")
})

# =============================================================================
# SECTION 1: (HISTOGRAMS REMOVED)
# =============================================================================

cat("\n===== SECTION 1: SKIPPED (histograms removed) =====\n")

# =============================================================================
# SECTION 2: SCATTER PLOTS (LOH TIME vs TE COUNTS)
# =============================================================================

cat("\n===== SECTION 2: SCATTER PLOTS =====\n")

# LOH time vs total TEs (with jitter for overlapping points)
tryCatch({
  if ("time" %in% colnames(data_loh) && "total" %in% colnames(data_loh)) {
    y_max <- max(data_loh$total, na.rm = TRUE)
    p <- ggplot(data_loh, aes(x = time, y = total)) +
      geom_jitter(width = 0.02, height = 0, size = 2) +
      geom_smooth(method = "lm", col = "#5FBFF9") +
      coord_cartesian(xlim = c(0, 1), ylim = c(0, y_max * 1.1)) +
      labs(x = "LOH Time", y = "Total TE Count")
    titled_print(p, "LOH Time vs Total TEs")
    ggsave(paste0(plot_dir, "loh/scatter_time_total.png"), plot = p, width = 6, height = 5)
    cat("  - LOH time vs total scatter created (n =", nrow(data_loh), ")\n")
  }
}, error = function(e) {
  cat("Warning: Could not create LOH time vs total scatter:", e$message, "\n")
})

# LOH time vs LINE1 (with jitter for overlapping points)
tryCatch({
  if ("time" %in% colnames(data_loh) && "LINE1" %in% colnames(data_loh)) {
    y_max <- max(data_loh$LINE1, na.rm = TRUE)
    p <- ggplot(data_loh, aes(x = time, y = LINE1)) +
      geom_jitter(width = 0.02, height = 0, size = 2) +
      geom_smooth(method = "lm", col = "#5FBFF9") +
      coord_cartesian(xlim = c(0, 1), ylim = c(0, max(50, y_max * 1.1))) +
      labs(x = "LOH Time", y = "LINE1 Count")
    titled_print(p, "LOH Time vs LINE1")
    ggsave(paste0(plot_dir, "loh/scatter_time_line1.png"), plot = p, width = 6, height = 5)
    cat("  - LOH time vs LINE1 scatter created (n =", nrow(data_loh), ")\n")
  }
}, error = function(e) {
  cat("Warning: Could not create LOH time vs LINE1 scatter:", e$message, "\n")
})

# Age at diagnosis vs LOH time (with jitter for overlapping points)
tryCatch({
  if ("time" %in% colnames(data_loh) && "age_at_diagnosis" %in% colnames(data_loh)) {
    age_max <- max(data_loh$age_at_diagnosis, na.rm = TRUE)
    p <- ggplot(data_loh, aes(x = age_at_diagnosis, y = time)) +
      geom_jitter(width = 0, height = 0.02, size = 2) +
      geom_smooth(method = "lm", col = "#5FBFF9") +
      coord_cartesian(xlim = c(0, age_max * 1.1), ylim = c(0, 1)) +
      labs(x = "Age at Diagnosis", y = "LOH Time")
    titled_print(p, "Age at Diagnosis vs LOH Time")
    ggsave(paste0(plot_dir, "loh/scatter_age_time.png"), plot = p, width = 6, height = 5)
    cat("  - Age vs LOH time scatter created (n =", nrow(data_loh), ")\n")
  }
}, error = function(e) {
  cat("Warning: Could not create age vs LOH time scatter:", e$message, "\n")
})

# Age at diagnosis vs total TEs (use data_loh)
tryCatch({
  if ("age_at_diagnosis" %in% colnames(data_loh) && "total" %in% colnames(data_loh)) {
    age_max <- max(data_loh$age_at_diagnosis, na.rm = TRUE)
    y_max <- max(data_loh$total, na.rm = TRUE)
    p <- scatter_template(data_loh, "age_at_diagnosis", "total",
                          c(0, age_max * 1.1), c(0, max(y_max * 1.1, 100)),
                          "Age at Diagnosis (days)", "Total TE Count")
    titled_print(p, "Age at Diagnosis vs Total TEs (LOH samples)")
    ggsave(paste0(plot_dir, "loh/scatter_age_total.png"), plot = p, width = 6, height = 5)
    cat("  - Age vs total TEs scatter created\n")
  }
}, error = function(e) {
  cat("Warning: Could not create age vs total scatter:", e$message, "\n")
})

# =============================================================================
# SECTION 3: FACETED SCATTER PLOTS
# =============================================================================

cat("\n===== SECTION 3: FACETED SCATTER PLOTS =====\n")

# Age at diagnosis vs LOH time by TP53 status (with jitter)
tryCatch({
  if (all(c("time", "age_at_diagnosis", "TP53_status") %in% colnames(data_loh))) {
    age_max <- max(data_loh$age_at_diagnosis, na.rm = TRUE)
    p <- ggplot(data_loh, aes(x = age_at_diagnosis, y = time)) +
      geom_jitter(width = 0, height = 0.02, size = 2) +
      geom_smooth(method = "lm", col = "#5FBFF9") +
      coord_cartesian(xlim = c(0, age_max * 1.1), ylim = c(0, 1)) +
      labs(x = "Age at Diagnosis", y = "LOH Time") +
      facet_wrap(~TP53_status)
    titled_print(p, "Age vs LOH Time by TP53 Status")
    ggsave(paste0(plot_dir, "loh/scatter_facet_age_time_tp53.png"), plot = p, width = 10, height = 5)
    cat("  - Faceted age vs LOH time scatter created (n =", nrow(data_loh), ")\n")
  }
}, error = function(e) {
  cat("Warning: Could not create faceted age vs LOH time scatter:", e$message, "\n")
})

# LOH time vs total by TP53 status (with jitter)
tryCatch({
  if (all(c("time", "total", "TP53_status") %in% colnames(data_loh))) {
    y_max <- max(data_loh$total, na.rm = TRUE)
    p <- ggplot(data_loh, aes(x = time, y = total)) +
      geom_jitter(width = 0.02, height = 0, size = 2) +
      geom_smooth(method = "lm", col = "#5FBFF9") +
      coord_cartesian(xlim = c(0, 1), ylim = c(0, y_max * 1.1)) +
      labs(x = "LOH Time", y = "Total TE Count") +
      facet_wrap(~TP53_status)
    titled_print(p, "LOH Time vs Total by TP53 Status")
    ggsave(paste0(plot_dir, "loh/scatter_facet_time_total_tp53.png"), plot = p, width = 10, height = 5)
    cat("  - Faceted LOH time vs total scatter created (n =", nrow(data_loh), ")\n")
  }
}, error = function(e) {
  cat("Warning: Could not create faceted LOH time vs total scatter:", e$message, "\n")
})

# =============================================================================
# SECTION 4: BOXPLOTS WITH STATISTICAL TESTS
# =============================================================================

cat("\n===== SECTION 4: BOXPLOTS WITH STATISTICAL TESTS =====\n")

# Variant location vs LOH time
tryCatch({
  if (all(c("Variant_location", "time") %in% colnames(data_loh))) {
    # Filter to groups with enough samples
    data_loh_filtered <- filter_by_min_samples(data_loh, "Variant_location", 3)
    n_levels <- length(unique(data_loh_filtered$Variant_location))

    if (n_levels >= 2) {
      p <- ggplot(data_loh_filtered, aes(x = Variant_location, y = time, fill = Variant_location)) +
        geom_boxplot(outlier.shape = NA) +
        scale_fill_manual(values = colours_3[1:min(n_levels, 3)]) +
        geom_jitter(position = position_jitterdodge(jitter.width = 0.2), color = "black", size = 1.5) +
        labs(x = "Variant Location", y = "LOH Time") +
        guides(fill = "none")
      titled_print(p, "Variant Location vs LOH Time")
      ggsave(paste0(plot_dir, "loh/boxplot_variant_location_time.png"), plot = p, width = 6, height = 5)
      cat("  - Variant location boxplot created\n")
    } else {
      cat("  - Skipping variant location boxplot (insufficient groups)\n")
    }
  }
}, error = function(e) {
  cat("Warning: Could not create variant location boxplot:", e$message, "\n")
})

# Sex vs LOH time (t-test)
tryCatch({
  if (all(c("sex", "time") %in% colnames(data_loh))) {
    data_loh_sex <- data_loh %>% filter(!is.na(sex) & !is.na(time))
    n_sex_levels <- length(unique(data_loh_sex$sex))

    if (n_sex_levels == 2 && nrow(data_loh_sex) >= 4) {
      # Perform t-test
      t_result <- t.test(time ~ sex, data = data_loh_sex)
      p_val <- t_result$p.value
      p_label <- ifelse(p_val < 0.001, "p < 0.001", paste0("p = ", round(p_val, 3)))

      p <- ggplot(data_loh_sex, aes(x = sex, y = time, fill = sex)) +
        geom_boxplot(outlier.shape = NA) +
        scale_fill_manual(values = colours) +
        geom_jitter(position = position_jitterdodge(jitter.width = 0.2), color = "black", size = 1.5) +
        labs(x = "Sex", y = "LOH Time", subtitle = p_label) +
        guides(fill = "none")
      titled_print(p, "Sex vs LOH Time (t-test)")
      ggsave(paste0(plot_dir, "loh/boxplot_ttest_sex_time.png"), plot = p, width = 5, height = 5)
      cat("  - Sex t-test boxplot created\n")
    } else {
      cat("  - Skipping sex t-test (need exactly 2 groups)\n")
    }
  }
}, error = function(e) {
  cat("Warning: Could not create sex t-test boxplot:", e$message, "\n")
})

# TP53 status vs LOH time (Wilcoxon)
tryCatch({
  if (all(c("TP53_status", "time") %in% colnames(data_loh))) {
    data_loh_tp53 <- data_loh %>% filter(!is.na(TP53_status) & !is.na(time))
    n_tp53_levels <- length(unique(data_loh_tp53$TP53_status))

    if (n_tp53_levels == 2 && nrow(data_loh_tp53) >= 4) {
      # Perform Wilcoxon test
      w_result <- wilcox.test(time ~ TP53_status, data = data_loh_tp53)
      p_val <- w_result$p.value
      p_label <- ifelse(p_val < 0.001, "p < 0.001", paste0("p = ", round(p_val, 3)))

      p <- ggplot(data_loh_tp53, aes(x = TP53_status, y = time, fill = TP53_status)) +
        geom_boxplot(outlier.shape = NA) +
        scale_fill_manual(values = colours) +
        geom_jitter(position = position_jitterdodge(jitter.width = 0.2), color = "black", size = 1.5) +
        labs(x = "TP53 Status", y = "LOH Time", subtitle = p_label) +
        guides(fill = "none")
      titled_print(p, "TP53 Status vs LOH Time (Wilcoxon)")
      ggsave(paste0(plot_dir, "loh/boxplot_wilcox_tp53_time.png"), plot = p, width = 5, height = 5)

      # Also output the Wilcoxon test results
      write_output(quote(wilcox.test(time ~ TP53_status, data = data_loh_tp53)),
                   "Wilcoxon test for LOH timing by TP53 status")
      cat("  - TP53 Wilcoxon boxplot created\n")
    } else {
      cat("  - Skipping TP53 Wilcoxon (need exactly 2 groups)\n")
    }
  }
}, error = function(e) {
  cat("Warning: Could not create TP53 Wilcoxon boxplot:", e$message, "\n")
})

# Tumor type vs LOH time (Kruskal-Wallis)
tryCatch({
  if (all(c("tumor_type", "time") %in% colnames(data_loh))) {
    # Filter to groups with at least 3 samples
    data_loh_tumor <- filter_by_min_samples(data_loh, "tumor_type", 3)
    data_loh_tumor <- data_loh_tumor %>% filter(!is.na(time))
    n_tumor_levels <- length(unique(data_loh_tumor$tumor_type))

    if (n_tumor_levels >= 2) {
      # Perform Kruskal-Wallis test
      kw_result <- kruskal.test(time ~ tumor_type, data = data_loh_tumor)
      p_val <- kw_result$p.value
      p_label <- ifelse(p_val < 0.001, "p < 0.001", paste0("p = ", round(p_val, 3)))

      tumor_colours <- colours_8[1:min(n_tumor_levels, 8)]
      p <- ggplot(data_loh_tumor, aes(x = tumor_type, y = time, fill = tumor_type)) +
        geom_boxplot(outlier.shape = NA) +
        scale_fill_manual(values = tumor_colours) +
        geom_jitter(position = position_jitterdodge(jitter.width = 0.2), color = "black", size = 1.5) +
        labs(x = "Tumor Type", y = "LOH Time", subtitle = p_label) +
        guides(fill = "none") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      titled_print(p, "Tumor Type vs LOH Time (Kruskal-Wallis)")
      ggsave(paste0(plot_dir, "loh/boxplot_kruskal_tumor_time.png"), plot = p, width = 8, height = 5)
      cat("  - Tumor type Kruskal-Wallis boxplot created\n")
    } else {
      cat("  - Skipping tumor type Kruskal-Wallis (insufficient groups after filtering)\n")
    }
  }
}, error = function(e) {
  cat("Warning: Could not create tumor type Kruskal-Wallis boxplot:", e$message, "\n")
})

# =============================================================================
# SECTION 5: LINEAR MODELS
# =============================================================================

cat("\n===== SECTION 5: LINEAR MODELS =====\n")

# Filter to TP53 mutants only for linear model analysis
tryCatch({
  if ("TP53_status" %in% colnames(data_loh)) {
    data_loh_mutant <- data_loh %>%
      filter(TP53_status == "Mutant" | TP53_status == "LFS")

    cat("Samples for linear model (TP53 mutants):", nrow(data_loh_mutant), "\n")

    if (nrow(data_loh_mutant) >= 5) {
      # Simple linear model: total TEs vs LOH time
      cat("\n--- Simple Linear Model: Total TEs vs LOH Time ---\n")
      if (all(c("total", "time") %in% colnames(data_loh_mutant))) {
        model_fit <- linear_model_one_variable(data_loh_mutant, "total", "time")
        cat("\nCorrelation:", model_fit$correlation, "\n")
        cat("R-squared:", model_fit$r_squared, "\n")
        cat("P-value:", model_fit$p_value, "\n")

        # Get axis limits
        y_max <- max(data_loh_mutant$total, na.rm = TRUE)

        # Plot with linear model stats (using jitter for overlapping points)
        p <- ggplot(data_loh_mutant, aes(x = time, y = total)) +
          geom_jitter(width = 0.02, height = 0, size = 2) +
          geom_smooth(method = "lm", col = "#5FBFF9") +
          coord_cartesian(xlim = c(0, 1), ylim = c(0, y_max * 1.1)) +
          labs(x = "LOH Time", y = "Total TE Count") +
          annotate("text", x = 0.7, y = y_max * 0.9,
                   label = paste("r =", format(model_fit$correlation, digits = 3),
                                 "\nR² =", format(model_fit$r_squared, digits = 3),
                                 "\np =", format(model_fit$p_value, digits = 3)),
                   size = 4, hjust = 0)
        titled_print(p, "LOH Time vs Total TEs (Linear Model)")
        ggsave(paste0(plot_dir, "loh/scatter_lm_time_total.png"), plot = p, width = 7, height = 5)
        cat("  - Simple LM scatter created (n =", nrow(data_loh_mutant), ")\n")

        # Plot with shape and color variables
        if (all(c("TP53_status", "tumor_type") %in% colnames(data_loh_mutant))) {
          p2 <- ggplot(data_loh_mutant, aes(x = time, y = total, shape = TP53_status, color = tumor_type)) +
            geom_jitter(width = 0.02, height = 0, size = 2) +
            geom_smooth(method = "lm", col = "#5FBFF9", aes(group = 1)) +
            coord_cartesian(xlim = c(0, 1), ylim = c(0, min(900, y_max * 1.1))) +
            labs(x = "LOH Time", y = "Total TE Count")
          titled_print(p2, "LOH Time vs Total TEs (by TP53 status and tumor type)")
          ggsave(paste0(plot_dir, "loh/scatter_lm_time_total_grouped.png"), plot = p2, width = 8, height = 5)
          cat("  - Grouped LM scatter created (n =", nrow(data_loh_mutant), ")\n")
        }
      }

      # Multiple variable linear model
      cat("\n--- Multiple Variable Linear Model ---\n")
      if (all(c("sex", "TP53_status", "total", "time") %in% colnames(data_loh_mutant))) {
        variables <- c("sex", "TP53_status")
        tryCatch({
          model_multi <- linear_model_multiple_variables(data_loh_mutant, "total", "time", variables, 4)
          cat("\nMultiple LM Results:\n")
          cat("Correlation:", model_multi$correlation, "\n")
          cat("R-squared:", model_multi$r_squared, "\n")
          cat("P-value:", model_multi$p_value, "\n")
        }, error = function(e) {
          cat("Warning: Multiple variable LM failed:", e$message, "\n")
        })
      }
    } else {
      cat("  - Insufficient TP53 mutant samples for linear model analysis\n")
    }
  }
}, error = function(e) {
  cat("Warning: Could not perform linear model analysis:", e$message, "\n")
})

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n===== SCRIPT COMPLETED =====\n")
cat("Plots saved to:", paste0(plot_dir, "loh/"), "\n")
cat("Data exported to:", paste0(r_dir_files, "te_loh.csv"), "\n")

# Close module-specific sink
close_module_sink()
