#!/usr/bin/env Rscript

# Tumour TE Visualization - Full-length LINE1 and Source Analysis
# Full-length LINE1 (>=5900bp) and Source/Transduction analysis

# Source common setup and load data
source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/00_viz_common_setup.R")
REQUIRED_DATA <- c("count_matrix", "expand", "clinical", "ancestry")
source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/00_viz_load_data_tumour.R")

# Create output directory
output_dir <- paste0(plot_dir, "fulllength_source/")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(paste0(output_dir, "files/"), showWarnings = FALSE, recursive = TRUE)

# Initialize module-specific text output
init_module_sink(output_dir, "FULLLENGTH_YOUNG_SOURCE")

cat("Running 02_te_viz_tumour_14_fulllength_young_source.R...\n")

# Define x-axis label for TP53 status (3-way: WT/Somatic/Germline)
x_tp53 <- expression(italic("TP53") * " status")

# Define covariates
covar_full <- c("predicted_ancestry_thres", "age_at_diagnosis", "med_cov", "total_reads", "med_read_len", "pct_chimeras", "avg_quality", "tumor_type_grouped")
covar_no_ancestry <- c("age_at_diagnosis", "med_cov", "total_reads", "med_read_len", "pct_chimeras", "avg_quality", "tumor_type_grouped")
covar_no_tumour_type <- c("predicted_ancestry_thres", "age_at_diagnosis", "med_cov", "total_reads", "med_read_len", "pct_chimeras", "avg_quality")

# Helper function to create 3-way TP53 grouping (WT/Somatic/Germline)
# Uses add_tp53_3level from functions_te.R which reads somatic TP53 variants from file
create_tp53_3way <- function(df) {
  # Use the existing add_tp53_3level function which handles somatic TP53 properly
  if (exists("add_tp53_3level")) {
    df <- add_tp53_3level(df)
    # Rename TP53_3level to TP53_3way for consistency
    if ("TP53_3level" %in% colnames(df)) {
      df <- df %>% rename(TP53_3way = TP53_3level)
    }
  } else {
    # Fallback: use simple binary TP53 status if function not available
    warning("add_tp53_3level function not found, using binary TP53 status")
    df <- df %>%
      mutate(TP53_3way = ifelse(TP53_status == "Mutant", "Germline", "WT"))
  }
  return(df)
}

# Helper function for repeated downsampling analysis comparing TP53 groups
# Downsamples larger groups to match the smallest group size, repeats n_iterations times
# Returns distribution of median differences and p-values
repeated_downsample_tp53 <- function(df, value_col = "n", group_col = "TP53_3way",
                                      n_iterations = 1000, seed = 42) {
  set.seed(seed)

  # Check if required columns exist
  if (!value_col %in% colnames(df)) {
    warning(sprintf("Column '%s' not found in data. Available columns: %s",
                    value_col, paste(head(colnames(df), 10), collapse = ", ")))
    return(NULL)
  }
  if (!group_col %in% colnames(df)) {
    warning(sprintf("Column '%s' not found in data", group_col))
    return(NULL)
  }

  # Remove NA values using base R to avoid dplyr issues
  df_clean <- df[!is.na(df[[group_col]]) & !is.na(df[[value_col]]), ]

  # Get group sizes
  group_sizes <- df_clean %>% group_by(.data[[group_col]]) %>% summarise(n = n()) %>% deframe()
  min_size <- min(group_sizes)
  groups <- names(group_sizes)

  if (length(groups) < 2) {
    return(NULL)
  }

  cat(sprintf("  Group sizes: %s\n", paste(names(group_sizes), group_sizes, sep="=", collapse=", ")))
  cat(sprintf("  Downsampling to n=%d per group, %d iterations\n", min_size, n_iterations))

  # Run repeated downsampling
  results <- lapply(1:n_iterations, function(i) {
    # Downsample each group to min_size
    df_downsampled <- df_clean %>%
      group_by(.data[[group_col]]) %>%
      slice_sample(n = min_size, replace = FALSE) %>%
      ungroup()

    # Calculate medians for each group
    medians <- df_downsampled %>%
      group_by(.data[[group_col]]) %>%
      summarise(median_val = median(.data[[value_col]], na.rm = TRUE)) %>%
      deframe()

    # Kruskal-Wallis test across all groups
    kw_test <- kruskal.test(as.formula(paste(value_col, "~", group_col)), data = df_downsampled)

    # Pairwise Wilcoxon tests
    pairwise_pvals <- list()
    for (j in 1:(length(groups)-1)) {
      for (k in (j+1):length(groups)) {
        g1 <- groups[j]
        g2 <- groups[k]
        pair_data <- df_downsampled %>% filter(.data[[group_col]] %in% c(g1, g2))
        wt <- wilcox.test(as.formula(paste(value_col, "~", group_col)), data = pair_data)
        pairwise_pvals[[paste(g1, "vs", g2)]] <- wt$p.value
      }
    }

    list(
      medians = medians,
      kw_pval = kw_test$p.value,
      pairwise = pairwise_pvals
    )
  })

  # Aggregate results
  kw_pvals <- sapply(results, function(x) x$kw_pval)

  # Get pairwise comparisons
  pairwise_names <- names(results[[1]]$pairwise)
  pairwise_pvals <- lapply(pairwise_names, function(name) {
    sapply(results, function(x) x$pairwise[[name]])
  })
  names(pairwise_pvals) <- pairwise_names

  # Get median distributions for each group
  median_distributions <- lapply(groups, function(g) {
    sapply(results, function(x) x$medians[g])
  })
  names(median_distributions) <- groups

  list(
    group_sizes = group_sizes,
    min_size = min_size,
    n_iterations = n_iterations,
    kw_pvals = kw_pvals,
    pairwise_pvals = pairwise_pvals,
    median_distributions = median_distributions
  )
}

# Helper function to plot downsampling results (boxplot + histogram)
plot_downsample_results <- function(ds_results, title = "Repeated Downsampling", y_lab = "TE count") {
  if (is.null(ds_results)) return(NULL)

  # Create median distribution plot
  median_df <- do.call(rbind, lapply(names(ds_results$median_distributions), function(g) {
    data.frame(
      group = g,
      median = ds_results$median_distributions[[g]]
    )
  }))

  # Order groups: WT, Somatic, Germline
  if (all(c("WT", "Somatic", "Germline") %in% unique(median_df$group))) {
    median_df$group <- factor(median_df$group, levels = c("WT", "Somatic", "Germline"))
  }

  # Calculate summary stats for annotation
  summary_stats <- median_df %>%
    group_by(group) %>%
    summarise(
      mean_median = mean(median),
      sd_median = sd(median),
      .groups = "drop"
    )

  # Calculate proportion of significant p-values
  prop_sig_kw <- mean(ds_results$kw_pvals < 0.05)

  pairwise_sig <- sapply(ds_results$pairwise_pvals, function(pvals) mean(pvals < 0.05))

  # Create boxplot
  p_box <- ggplot(median_df, aes(x = group, y = median, fill = group)) +
    geom_boxplot(alpha = 0.7, outlier.shape = 21, outlier.size = 1) +
    scale_fill_manual(values = c("WT" = "#4DAF4A", "Somatic" = "#FF7F00", "Germline" = "#E41A1C")) +
    labs(
      title = title,
      subtitle = sprintf("n=%d per group, %d iterations | KW p<0.05: %.1f%%",
                        ds_results$min_size, ds_results$n_iterations, prop_sig_kw * 100),
      x = expression(italic("TP53") * " status"),
      y = paste("Median", y_lab, "(across iterations)"),
      caption = paste(sapply(names(pairwise_sig), function(x)
        sprintf("%s p<0.05: %.1f%%", x, pairwise_sig[x] * 100)), collapse = " | ")
    ) +
    theme_bw() +
    theme(
      legend.position = "none",
      plot.caption = element_text(size = 8)
    )

  # Create histogram (faceted by group)
  p_hist <- ggplot(median_df, aes(x = median, fill = group)) +
    geom_histogram(bins = 30, alpha = 0.7, color = "black", linewidth = 0.2) +
    facet_wrap(~group, ncol = 1, scales = "free_y") +
    scale_fill_manual(values = c("WT" = "#4DAF4A", "Somatic" = "#FF7F00", "Germline" = "#E41A1C")) +
    labs(
      title = paste(title, "- Distribution"),
      subtitle = sprintf("n=%d per group, %d iterations", ds_results$min_size, ds_results$n_iterations),
      x = paste("Median", y_lab),
      y = "Frequency"
    ) +
    theme_bw() +
    theme(legend.position = "none")

  # Return both plots as a list
  return(list(boxplot = p_box, histogram = p_hist))
}

# Filter to LINE1 only
te_aff_expand_line_t <- te_aff_expand_t %>% filter(ALT == "LINE1")

################################################################################
#### FULL-LENGTH LINE1 ANALYSIS ####
################################################################################
cat("\n===== FULL-LENGTH LINE1 ANALYSIS =====\n")

# Summary of LINE1 lengths
write_output(quote(summary(te_aff_expand_line_t$SV_length)), "Summary of LINE1 SV_length (affected)")

# Filter to full-length (>=5900bp)
te_aff_expand_line_fulllength_t <- te_aff_expand_line_t %>% filter(SV_length >= 5900)
write_output(quote(nrow(te_aff_expand_line_fulllength_t)), "Number of full-length LINE1 (SV_length >= 5900)")

# Process combinations and add nohit samples
te_fl_processed <- process_all_combinations(te_aff_expand_line_fulllength_t)
load(paste0(r_dir, "nohits_final_te_count_t_te_aff_selected_t", ".RData"))
te_fl_processed <- as.data.frame(add_nohit_samples(te_fl_processed, nohits))
# Merge with clinical and metrics
te_fl_processed <- merge_dfs(te_fl_processed, clinical, include_all_x = FALSE)
te_fl_processed <- merge_dfs(te_fl_processed, metrics, include_all_x = TRUE)

# Create tumor_type_grouped with >=10 samples threshold
te_fl_processed <- te_fl_processed %>%
  mutate(tumor_type_grouped = ifelse(
    tumor_type %in% names(table(tumor_type))[table(tumor_type) >= 10],
    tumor_type, "Other"
  ))

# Create 3-way TP53 grouping
te_fl_processed <- create_tp53_3way(te_fl_processed)

cat("\n--- Full-length LINE1 by TP53 Status (3-way: WT/Somatic/Germline) ---\n")

# Kruskal test
cat("Plotting full-length LINE1 count by TP53 status (Kruskal test)...\n")
p_fl_tp53 <- plot_count_kruskal(te_fl_processed, type=NA, chr=NA, group="TP53_3way",
                                 log_scale=TRUE, x_lab=x_tp53, y_lab="Full-length LINE1 count",
                                 x_order=c("WT", "Somatic", "Germline"))
titled_print(p_fl_tp53, "Full-length LINE1 by TP53 status (3-way)")
ggsave(paste0(output_dir, "full_length_tp53_kruskal.png"), plot=p_fl_tp53, width = 5, height = 5)

# Linear model
cat("Plotting full-length LINE1 count by TP53 status (Linear model)...\n")
tryCatch({
  p_fl_tp53_lm <- plot_count_lm(te_fl_processed, type=NA, chr=NA, group="TP53_3way",
                                 y_lab="Full-length LINE1 count", covariates=covar_full,
                                 residuals=FALSE, log_scale=FALSE, x_lab=x_tp53, min_samples=3)
  titled_print(p_fl_tp53_lm, "Full-length LINE1 by TP53 status (LM)")
  ggsave(paste0(output_dir, "full_length_tp53_lm.png"), plot=p_fl_tp53_lm, width = 5, height = 5)
}, error = function(e) {
  cat("Warning: Could not create LM plot for TP53 status:", e$message, "\n")
})

# Repeated downsampling analysis
cat("Running repeated downsampling analysis for full-length LINE1 by TP53 status...\n")
tryCatch({
  ds_fl_tp53 <- repeated_downsample_tp53(te_fl_processed, value_col = "total", group_col = "TP53_3way",
                                          n_iterations = 1000, seed = 42)
  if (!is.null(ds_fl_tp53)) {
    p_fl_tp53_ds <- plot_downsample_results(ds_fl_tp53,
                                             title = "Full-length LINE1 by TP53 (Repeated Downsampling)",
                                             y_lab = "Full-length LINE1 count")
    # Print and save boxplot
    titled_print(p_fl_tp53_ds$boxplot, "Full-length LINE1 by TP53 status (Downsampling - Boxplot)")
    ggsave(paste0(output_dir, "full_length_tp53_downsample_boxplot.png"), plot=p_fl_tp53_ds$boxplot, width = 6, height = 5)
    # Print and save histogram
    titled_print(p_fl_tp53_ds$histogram, "Full-length LINE1 by TP53 status (Downsampling - Histogram)")
    ggsave(paste0(output_dir, "full_length_tp53_downsample_histogram.png"), plot=p_fl_tp53_ds$histogram, width = 6, height = 7)

    # Save results to file
    ds_summary <- data.frame(
      comparison = c("Kruskal-Wallis", names(ds_fl_tp53$pairwise_pvals)),
      prop_significant = c(
        mean(ds_fl_tp53$kw_pvals < 0.05),
        sapply(ds_fl_tp53$pairwise_pvals, function(x) mean(x < 0.05))
      ),
      median_pval = c(
        median(ds_fl_tp53$kw_pvals),
        sapply(ds_fl_tp53$pairwise_pvals, median)
      )
    )
    write.csv(ds_summary, paste0(output_dir, "files/full_length_tp53_downsample_results.csv"), row.names = FALSE)
    cat("Downsampling results:\n")
    print(ds_summary)
  }
}, error = function(e) {
  cat("Warning: Could not run downsampling analysis:", e$message, "\n")
})

cat("\n--- Full-length LINE1 by Ancestry ---\n")

if ("predicted_ancestry_thres" %in% colnames(te_fl_processed)) {
  te_fl_ancestry <- te_fl_processed %>%
    filter(!is.na(predicted_ancestry_thres) & predicted_ancestry_thres != "Unknown" & predicted_ancestry_thres != "Undefined")

  if (nrow(te_fl_ancestry) > 0 && length(unique(te_fl_ancestry$predicted_ancestry_thres)) > 1) {
    # Kruskal test
    cat("Plotting full-length LINE1 count by ancestry (Kruskal test)...\n")
    p_fl_ancestry <- plot_count_kruskal(te_fl_ancestry, type=NA, chr=NA, group="predicted_ancestry_thres",
                                         log_scale=TRUE, x_lab="Ancestry", y_lab="Full-length LINE1 count")
    titled_print(p_fl_ancestry, "Full-length LINE1 by ancestry")
    ggsave(paste0(output_dir, "full_length_ancestry_kruskal.png"), plot=p_fl_ancestry, width = 6, height = 5)

    # Linear model (without ancestry as covariate)
    cat("Plotting full-length LINE1 count by ancestry (Linear model)...\n")
    tryCatch({
      p_fl_ancestry_lm <- plot_count_lm(te_fl_ancestry, type=NA, chr=NA, group="predicted_ancestry_thres",
                                         y_lab="Full-length LINE1 count", covariates=covar_no_ancestry,
                                         residuals=FALSE, log_scale=FALSE, x_lab="Ancestry", min_samples=3,
                                         fill_palette=color_palette_7)
      titled_print(p_fl_ancestry_lm, "Full-length LINE1 by ancestry (LM)")
      ggsave(paste0(output_dir, "full_length_ancestry_lm.png"), plot=p_fl_ancestry_lm, width = 6, height = 5)
    }, error = function(e) {
      cat("Warning: Could not create LM plot for ancestry:", e$message, "\n")
    })
  }
}

cat("\n--- Full-length LINE1 by Tumour Type ---\n")

if ("tumor_type_grouped" %in% colnames(te_fl_processed)) {
  te_fl_tt <- te_fl_processed %>% filter(!is.na(tumor_type_grouped))

  if (nrow(te_fl_tt) > 0 && length(unique(te_fl_tt$tumor_type_grouped)) > 1) {
    # Kruskal test
    cat("Plotting full-length LINE1 count by tumour type (Kruskal test)...\n")
    tryCatch({
      p_fl_tt <- plot_count_kruskal(te_fl_tt, type=NA, chr=NA, group="tumor_type_grouped",
                                     log_scale=TRUE, x_lab="Tumour type", y_lab="Full-length LINE1 count")
      titled_print(p_fl_tt, "Full-length LINE1 by tumour type")
      ggsave(paste0(output_dir, "full_length_tumour_type_kruskal.png"), plot=p_fl_tt, width = 8, height = 5)
    }, error = function(e) {
      cat("Warning: Could not create Kruskal plot for tumour type:", e$message, "\n")
    })

    # Linear model (without tumour type as covariate)
    cat("Plotting full-length LINE1 count by tumour type (Linear model)...\n")
    tryCatch({
      p_fl_tt_lm <- plot_count_lm(te_fl_tt, type=NA, chr=NA, group="tumor_type_grouped",
                                   y_lab="Full-length LINE1 count", covariates=covar_no_tumour_type,
                                   residuals=FALSE, log_scale=FALSE, x_lab="Tumour type", min_samples=3)
      titled_print(p_fl_tt_lm, "Full-length LINE1 by tumour type (LM)")
      ggsave(paste0(output_dir, "full_length_tumour_type_lm.png"), plot=p_fl_tt_lm, width = 8, height = 5)
    }, error = function(e) {
      cat("Warning: Could not create LM plot for tumour type:", e$message, "\n")
    })
  }
}

################################################################################
#### SOURCE/TRANSDUCTION ANALYSIS ####
################################################################################
cat("\n===== SOURCE/TRANSDUCTION ANALYSIS =====\n")

tryCatch({
  te_aff_expand_line_t <- extract_info_fields(te_aff_expand_line_t)

  cat("Number of unique LINE1 sources (affected):", length(unique(te_aff_expand_line_t$source)), "\n")
  cat("Number of LINE1 insertions (affected):", nrow(te_aff_expand_line_t), "\n")

  # Sources with multiple occurrences
  sources_with_multiple <- table(te_aff_expand_line_t$source)[table(te_aff_expand_line_t$source) > 1]
  if (length(sources_with_multiple) > 0) {
    cat("\nLINE1 sources with >1 occurrence (affected):\n")
    print(sources_with_multiple)
  } else {
    cat("\nNo LINE1 sources with >1 occurrence\n")
  }

  # Filter to transductions only
  sources_transduction <- te_aff_expand_line_t %>% filter(source != "not_transduction")

  # Create tumor_type_grouped with >=10 samples threshold
  if (nrow(sources_transduction) > 0 && "tumor_type" %in% colnames(sources_transduction)) {
    sources_transduction <- sources_transduction %>%
      mutate(tumor_type_grouped = ifelse(
        tumor_type %in% names(table(tumor_type))[table(tumor_type) >= 10],
        tumor_type, "Other"
      ))
  }

  if (nrow(sources_transduction) > 0) {
    cat("\nNumber of transductions:", nrow(sources_transduction), "\n")

    samples_with_sources <- table(sources_transduction[, c("sample", "source")])
    cat("\nSamples with LINE1 sources (transductions):\n")
    print(head(samples_with_sources, 50))

    # Summary of sources per sample
    sources <- sources_transduction %>%
      count(sample, source) %>%
      filter(n > 0) %>%
      arrange(desc(n))

    if (nrow(sources) > 0) {
      cat("\nSummary of LINE1 sources per sample (affected):\n")
      print(head(sources, 20))

      # Create summary by source with TP53 status breakdown
      sources_summary <- sources_transduction %>%
        group_by(source) %>%
        summarise(
          n = n(),
          n_TP53_mut = sum(TP53_status == "Mutant", na.rm = TRUE),
          n_TP53_wt = sum(TP53_status == "WT", na.rm = TRUE),
          samples_TP53_mut = paste(unique(sample[!is.na(TP53_status) & TP53_status == "Mutant"]), collapse = ";"),
          samples_TP53_wt = paste(unique(sample[!is.na(TP53_status) & TP53_status == "WT"]), collapse = ";"),
          all_samples = paste(unique(sample), collapse = ";"),
          .groups = "drop"
        ) %>%
        mutate(
          samples_TP53_mut = ifelse(samples_TP53_mut == "", NA_character_, samples_TP53_mut),
          samples_TP53_wt = ifelse(samples_TP53_wt == "", NA_character_, samples_TP53_wt)
        ) %>%
        arrange(desc(n)) %>%
        select(source, n, n_TP53_mut, n_TP53_wt, samples_TP53_mut, samples_TP53_wt, all_samples)

      # Save source analysis to CSV
      write.csv(sources_summary, paste0(output_dir, "files/line1_source_analysis.csv"), row.names = FALSE)
      cat("Source analysis saved to: files/line1_source_analysis.csv\n")

      # Plot source distribution
      if (exists("plot_te_source")) {
        titled_print(plot_te_source(sources), "TE Source Plot (affected)")
        ggsave(paste0(output_dir, "source_distribution.png"), width = 7, height = 5)
      }

      # Bar plot of top sources
      top_sources <- sources_summary %>% head(20)
      if (nrow(top_sources) > 0) {
        p_sources <- ggplot(top_sources, aes(x = reorder(source, n), y = n)) +
          geom_bar(stat = "identity", fill = "steelblue") +
          coord_flip() +
          theme_minimal() +
          labs(x = "Source Element", y = "Number of Transductions",
               title = "Top LINE1 Source Elements") +
          theme(axis.text.y = element_text(size = 8))

        titled_print(p_sources, "Top LINE1 source elements")
        ggsave(paste0(output_dir, "source_top20_barplot.png"), plot = p_sources, width = 8, height = 6)
      }

      # Source activity by TP53 status
      sources_by_tp53 <- sources_transduction %>%
        group_by(TP53_status) %>%
        summarise(
          n_transductions = n(),
          n_unique_sources = n_distinct(source),
          n_samples = n_distinct(sample),
          .groups = "drop"
        )
      write_output(quote(sources_by_tp53), "Source activity by TP53 status")

      # Plot transduction count by TP53 status (3-way)
      cat("\n--- Transduction Count by TP53 Status (3-way: WT/Somatic/Germline) ---\n")
      transduction_per_sample <- sources_transduction %>%
        group_by(sample, TP53_status) %>%
        summarise(n_transductions = n(), .groups = "drop")

      # Add samples with 0 transductions (filter to tumour samples only)
      # Note: Get TP53_status from clinical merge, not from te_aff_t, to avoid duplicate columns
      all_tumour_samples <- te_aff_t %>%
        filter(grepl("_T$", sample)) %>%
        select(sample) %>%
        distinct()
      transduction_per_sample <- all_tumour_samples %>%
        left_join(transduction_per_sample %>% select(sample, n_transductions), by = "sample") %>%
        mutate(n_transductions = ifelse(is.na(n_transductions), 0, n_transductions))

      # Merge with clinical and metrics
      transduction_per_sample <- merge_dfs(transduction_per_sample, clinical, include_all_x = FALSE)
      transduction_per_sample <- merge_dfs(transduction_per_sample, metrics, include_all_x = TRUE)

      # Create tumor_type_grouped with >=10 samples threshold
      transduction_per_sample <- transduction_per_sample %>%
        mutate(tumor_type_grouped = ifelse(
          tumor_type %in% names(table(tumor_type))[table(tumor_type) >= 10],
          tumor_type, "Other"
        ))

      # Create 3-way TP53 grouping
      transduction_per_sample <- create_tp53_3way(transduction_per_sample)

      # Rename count column to 'total' for function compatibility
      transduction_per_sample <- transduction_per_sample %>% rename(total = n_transductions)

      # Kruskal test using function
      cat("Plotting transduction count by TP53 status (Kruskal test)...\n")
      p_trans_tp53 <- plot_count_kruskal(transduction_per_sample, type=NA, chr=NA, group="TP53_3way",
                                          log_scale=TRUE, x_lab=x_tp53, y_lab="Transduction count",
                                          x_order=c("WT", "Somatic", "Germline"))
      titled_print(p_trans_tp53, "Transduction count by TP53 status (3-way Kruskal)")
      ggsave(paste0(output_dir, "source_transduction_tp53_kruskal.png"), plot = p_trans_tp53, width = 5, height = 5)

      # Linear model
      cat("Plotting transduction count by TP53 status (Linear model)...\n")
      tryCatch({
        p_trans_tp53_lm <- plot_count_lm(transduction_per_sample, type=NA, chr=NA, group="TP53_3way",
                                          y_lab="Transduction count", covariates=covar_full,
                                          residuals=FALSE, log_scale=FALSE, x_lab=x_tp53, min_samples=3)
        titled_print(p_trans_tp53_lm, "Transduction count by TP53 status (3-way LM)")
        ggsave(paste0(output_dir, "source_transduction_tp53_lm.png"), plot=p_trans_tp53_lm, width = 5, height = 5)
      }, error = function(e) {
        cat("Warning: Could not create LM plot for TP53 status:", e$message, "\n")
      })

      # Source activity by ancestry
      cat("\n--- Source Activity by Ancestry ---\n")
      if ("predicted_ancestry_thres" %in% colnames(transduction_per_sample)) {
        sources_by_ancestry <- sources_transduction %>%
          filter(!is.na(predicted_ancestry_thres) & predicted_ancestry_thres != "Unknown" & predicted_ancestry_thres != "Undefined") %>%
          group_by(predicted_ancestry_thres) %>%
          summarise(
            n_transductions = n(),
            n_unique_sources = n_distinct(source),
            n_samples = n_distinct(sample),
            .groups = "drop"
          )
        write_output(quote(sources_by_ancestry), "Source activity by ancestry")

        # Filter to valid ancestry
        transduction_ancestry <- transduction_per_sample %>%
          filter(!is.na(predicted_ancestry_thres) & predicted_ancestry_thres != "Unknown" & predicted_ancestry_thres != "Undefined")

        if (nrow(transduction_ancestry) > 0 && length(unique(transduction_ancestry$predicted_ancestry_thres)) > 1) {
          # Kruskal test using function
          cat("Plotting transduction count by ancestry (Kruskal test)...\n")
          p_trans_ancestry <- plot_count_kruskal(transduction_ancestry, type=NA, chr=NA, group="predicted_ancestry_thres",
                                                  log_scale=TRUE, x_lab="Ancestry", y_lab="Transduction count")
          titled_print(p_trans_ancestry, "Transduction count by Ancestry (Kruskal)")
          ggsave(paste0(output_dir, "source_transduction_ancestry_kruskal.png"), plot = p_trans_ancestry, width = 6, height = 5)

          # Linear model
          cat("Plotting transduction count by ancestry (Linear model)...\n")
          tryCatch({
            p_trans_ancestry_lm <- plot_count_lm(transduction_ancestry, type=NA, chr=NA, group="predicted_ancestry_thres",
                                                  y_lab="Transduction count", covariates=covar_no_ancestry,
                                                  residuals=FALSE, log_scale=FALSE, x_lab="Ancestry", min_samples=3,
                                                  fill_palette=color_palette_7)
            titled_print(p_trans_ancestry_lm, "Transduction count by Ancestry (LM)")
            ggsave(paste0(output_dir, "source_transduction_ancestry_lm.png"), plot=p_trans_ancestry_lm, width = 6, height = 5)
          }, error = function(e) {
            cat("Warning: Could not create LM plot for ancestry:", e$message, "\n")
          })
        }
      }

      # Source activity by tumour type
      cat("\n--- Source Activity by Tumour Type ---\n")
      if ("tumor_type_grouped" %in% colnames(transduction_per_sample)) {
        sources_by_tt <- sources_transduction %>%
          filter(!is.na(tumor_type_grouped)) %>%
          group_by(tumor_type_grouped) %>%
          summarise(
            n_transductions = n(),
            n_unique_sources = n_distinct(source),
            n_samples = n_distinct(sample),
            .groups = "drop"
          )
        write_output(quote(sources_by_tt), "Source activity by tumour type")

        # Filter to valid tumour type
        transduction_tt <- transduction_per_sample %>% filter(!is.na(tumor_type_grouped))

        if (nrow(transduction_tt) > 0 && length(unique(transduction_tt$tumor_type_grouped)) > 1) {
          # Kruskal test using function
          cat("Plotting transduction count by tumour type (Kruskal test)...\n")
          tryCatch({
            p_trans_tt <- plot_count_kruskal(transduction_tt, type=NA, chr=NA, group="tumor_type_grouped",
                                              log_scale=TRUE, x_lab="Tumour type", y_lab="Transduction count")
            titled_print(p_trans_tt, "Transduction count by Tumour Type (Kruskal)")
            ggsave(paste0(output_dir, "source_transduction_tumour_type_kruskal.png"), plot = p_trans_tt, width = 8, height = 5)
          }, error = function(e) {
            cat("Warning: Could not create Kruskal plot for tumour type:", e$message, "\n")
          })

          # Linear model
          cat("Plotting transduction count by tumour type (Linear model)...\n")
          tryCatch({
            p_trans_tt_lm <- plot_count_lm(transduction_tt, type=NA, chr=NA, group="tumor_type_grouped",
                                            y_lab="Transduction count", covariates=covar_no_tumour_type,
                                            residuals=FALSE, log_scale=FALSE, x_lab="Tumour type", min_samples=3)
            titled_print(p_trans_tt_lm, "Transduction count by Tumour Type (LM)")
            ggsave(paste0(output_dir, "source_transduction_tumour_type_lm.png"), plot=p_trans_tt_lm, width = 8, height = 5)
          }, error = function(e) {
            cat("Warning: Could not create LM plot for tumour type:", e$message, "\n")
          })
        }
      }

      # Top sources by group comparisons
      cat("\n--- Top Sources by Group ---\n")

      # Create 3-way TP53 grouping for sources_transduction
      sources_transduction <- create_tp53_3way(sources_transduction)

      # Top sources by TP53 status (3-way)
      top_sources_tp53 <- sources_transduction %>%
        group_by(TP53_3way, source) %>%
        summarise(n = n(), .groups = "drop") %>%
        group_by(TP53_3way) %>%
        slice_max(n, n = 10) %>%
        ungroup() %>%
        mutate(TP53_3way = factor(TP53_3way, levels = c("WT", "Somatic", "Germline")))

      if (nrow(top_sources_tp53) > 0) {
        p_top_tp53 <- ggplot(top_sources_tp53, aes(x = reorder(source, n), y = n, fill = TP53_3way)) +
          geom_bar(stat = "identity") +
          coord_flip() +
          facet_wrap(~TP53_3way, scales = "free_y") +
          theme_minimal() +
          labs(x = "Source Element", y = "Number of Transductions",
               title = "Top LINE1 Sources by TP53 Status (WT/Somatic/Germline)") +
          theme(axis.text.y = element_text(size = 7), legend.position = "none") +
          scale_fill_manual(values = c("WT" = "#377EB8", "Somatic" = "#FF7F00", "Germline" = "#E41A1C"))
        titled_print(p_top_tp53, "Top LINE1 Sources by TP53 Status (3-way)")
        ggsave(paste0(output_dir, "source_top_by_tp53.png"), plot = p_top_tp53, width = 12, height = 6)
      }

      # Top sources by ancestry (if enough data)
      if ("predicted_ancestry_thres" %in% colnames(sources_transduction)) {
        top_sources_ancestry <- sources_transduction %>%
          filter(!is.na(predicted_ancestry_thres) & predicted_ancestry_thres != "Unknown" & predicted_ancestry_thres != "Undefined") %>%
          group_by(predicted_ancestry_thres, source) %>%
          summarise(n = n(), .groups = "drop") %>%
          group_by(predicted_ancestry_thres) %>%
          slice_max(n, n = 5) %>%
          ungroup()

        if (nrow(top_sources_ancestry) > 0) {
          p_top_ancestry <- ggplot(top_sources_ancestry, aes(x = reorder(source, n), y = n, fill = predicted_ancestry_thres)) +
            geom_bar(stat = "identity") +
            coord_flip() +
            facet_wrap(~predicted_ancestry_thres, scales = "free_y") +
            theme_minimal() +
            labs(x = "Source Element", y = "Number of Transductions",
                 title = "Top LINE1 Sources by Ancestry") +
            theme(axis.text.y = element_text(size = 6), legend.position = "none")
          titled_print(p_top_ancestry, "Top sources by ancestry")
          ggsave(paste0(output_dir, "source_top_by_ancestry.png"), plot = p_top_ancestry, width = 12, height = 8)
        }
      }

      # Top sources by tumour type (if enough data)
      if ("tumor_type_grouped" %in% colnames(sources_transduction)) {
        top_sources_tt <- sources_transduction %>%
          filter(!is.na(tumor_type_grouped)) %>%
          group_by(tumor_type_grouped, source) %>%
          summarise(n = n(), .groups = "drop") %>%
          group_by(tumor_type_grouped) %>%
          slice_max(n, n = 5) %>%
          ungroup()

        if (nrow(top_sources_tt) > 0) {
          p_top_tt <- ggplot(top_sources_tt, aes(x = reorder(source, n), y = n, fill = tumor_type_grouped)) +
            geom_bar(stat = "identity") +
            coord_flip() +
            facet_wrap(~tumor_type_grouped, scales = "free_y") +
            theme_minimal() +
            labs(x = "Source Element", y = "Number of Transductions",
                 title = "Top LINE1 Sources by Tumour Type") +
            theme(axis.text.y = element_text(size = 6), legend.position = "none")
          titled_print(p_top_tt, "Top sources by tumour type")
          ggsave(paste0(output_dir, "source_top_by_tumour_type.png"), plot = p_top_tt, width = 12, height = 8)
        }
      }

    }
  } else {
    cat("\nNo LINE1 transductions found\n")
  }

}, error = function(e) {
  cat("Warning: Could not perform source analysis:", e$message, "\n")
})

################################################################################
#### COMBINED SUMMARY ####
################################################################################
cat("\n===== COMBINED SUMMARY =====\n")

summary_df <- data.frame(
  Analysis = c("Full-length LINE1 (>=5900bp)", "Transductions"),
  N_insertions = c(
    nrow(te_aff_expand_line_fulllength_t),
    if (exists("sources_transduction")) nrow(sources_transduction) else NA
  )
)
write_output(quote(summary_df), "Summary of LINE1 subset analyses")

cat("\n Full-length LINE1 and Source analysis completed successfully\n")

# Close module-specific sink
close_module_sink()
