#!/usr/bin/env Rscript

# Germline TE Visualization - Full-length LINE1, Young TE, and Source Analysis
# Full-length LINE1 (>=5900bp), Young TE subfamilies (L1HS, AluY, SVA_E, SVA_F), and Source/Transduction analysis

# Source common setup and load data
source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/00_viz_common_setup.R")
REQUIRED_DATA <- c("count_matrix", "expand", "clinical", "ancestry")
source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/00_viz_load_data_germline.R")

# Create output directory
output_dir <- paste0(plot_dir, "fulllength_young_source/")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(paste0(output_dir, "files/"), showWarnings = FALSE, recursive = TRUE)

# Initialize module-specific text output
init_module_sink(output_dir, "FULLLENGTH_YOUNG_SOURCE")

cat("Running 02_te_viz_germline_13_fulllength_young_source.R...\n")

# Define x-axis label for TP53 status
x_tp53 <- expression("Germline " * italic("TP53") * " status")

# Define covariates
covar_full <- c("predicted_ancestry_thres", "age_at_diagnosis", "med_cov", "total_reads", "med_read_len", "pct_chimeras", "avg_quality", "tumor_type_grouped")
covar_no_ancestry <- c("age_at_diagnosis", "med_cov", "total_reads", "med_read_len", "pct_chimeras", "avg_quality", "tumor_type_grouped")
covar_no_tumour_type <- c("predicted_ancestry_thres", "age_at_diagnosis", "med_cov", "total_reads", "med_read_len", "pct_chimeras", "avg_quality")

# Filter to LINE1 only
te_aff_expand_line <- te_aff_expand %>% filter(ALT == "LINE1")

################################################################################
#### FULL-LENGTH LINE1 ANALYSIS ####
################################################################################
cat("\n===== FULL-LENGTH LINE1 ANALYSIS =====\n")

# Summary of LINE1 lengths
write_output(quote(summary(te_aff_expand_line$SV_length)), "Summary of LINE1 SV_length (affected)")

# Filter to full-length (>=5900bp)
te_aff_expand_line_fulllength <- te_aff_expand_line %>% filter(SV_length >= 5900)
write_output(quote(nrow(te_aff_expand_line_fulllength)), "Number of full-length LINE1 (SV_length >= 5900)")

if (nrow(te_aff_expand_line_fulllength) > 0) {
  # Process combinations and add nohit samples
  te_fl_processed <- process_all_combinations(te_aff_expand_line_fulllength)
  load(paste0(r_dir, "nohits", ".RData"))
  if (nrow(nohits) > 0) {
    te_fl_processed <- as.data.frame(add_nohit_samples(te_fl_processed, nohits))
  }
  te_fl_processed <- merge_dfs(te_fl_processed, clinical, include_all_x = FALSE)

  # Merge with te_aff to get sequencing metrics and tumor_type_grouped
  metrics_cols <- c("sample", "med_cov", "total_reads", "med_read_len", "pct_chimeras", "avg_quality", "tumor_type", "tumor_type_grouped")
  metrics_cols_available <- metrics_cols[metrics_cols %in% colnames(te_aff)]
  if (length(metrics_cols_available) > 1) {
    te_fl_processed <- te_fl_processed %>%
      left_join(te_aff %>% select(all_of(metrics_cols_available)) %>% distinct(), by = "sample")
  }

  # Fix tumor_type column naming after merge (may have .x/.y suffixes)
  if ("tumor_type.x" %in% colnames(te_fl_processed) && !"tumor_type" %in% colnames(te_fl_processed)) {
    te_fl_processed <- te_fl_processed %>% rename(tumor_type = tumor_type.x)
  }
  if ("tumor_type_grouped.x" %in% colnames(te_fl_processed) && !"tumor_type_grouped" %in% colnames(te_fl_processed)) {
    te_fl_processed <- te_fl_processed %>% rename(tumor_type_grouped = tumor_type_grouped.x)
  }

  # If tumor_type_grouped still doesn't exist, create it with higher threshold for fewer groups
  if (!"tumor_type_grouped" %in% colnames(te_fl_processed) && "tumor_type" %in% colnames(te_fl_processed)) {
    te_fl_processed <- te_fl_processed %>%
      mutate(tumor_type_grouped = ifelse(
        tumor_type %in% names(table(tumor_type))[table(tumor_type) >= 10],
        tumor_type, "Other"
      ))
  }

  cat("\n--- Full-length LINE1 by TP53 Status ---\n")

  # Kruskal test
  cat("Plotting full-length LINE1 count by TP53 status (Kruskal test)...\n")
  p_fl_tp53 <- plot_count_kruskal(te_fl_processed, type=NA, chr=NA, group="TP53_status",
                                   log_scale=TRUE, x_lab=x_tp53, y_lab="Full-length LINE1 count")
  titled_print(p_fl_tp53, "Full-length LINE1 by TP53 status")
  ggsave(paste0(output_dir, "full_length_tp53_kruskal.png"), plot=p_fl_tp53, width = 4, height = 5)

  # Linear model
  cat("Plotting full-length LINE1 count by TP53 status (Linear model)...\n")
  tryCatch({
    p_fl_tp53_lm <- plot_count_lm(te_fl_processed, type=NA, chr=NA, group="TP53_status",
                                   y_lab="Full-length LINE1 count", covariates=covar_full,
                                   residuals=FALSE, log_scale=FALSE, x_lab=x_tp53, min_samples=3)
    titled_print(p_fl_tp53_lm, "Full-length LINE1 by TP53 status (LM)")
    ggsave(paste0(output_dir, "full_length_tp53_lm.png"), plot=p_fl_tp53_lm, width = 4, height = 5)
  }, error = function(e) {
    cat("Warning: Could not create LM plot for TP53 status:", e$message, "\n")
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
} else {
  cat("No full-length LINE1 insertions found\n")
}

################################################################################
#### YOUNG TE ANALYSIS (L1HS, AluY, SVA_E, SVA_F) ####
################################################################################
cat("\n===== YOUNG TE ANALYSIS (L1HS, AluY, SVA_E, SVA_F) =====\n")

# Check if subfamily column exists
if ("subfamily" %in% colnames(te_aff_expand)) {

  # Define young/active TE subfamilies
  young_family <- c("L1HS", "AluY", "SVA_E", "SVA_F")

  # Filter to young TE subfamilies from full expand data (not just LINE1)
  te_aff_expand_young <- te_aff_expand %>% filter(subfamily %in% young_family)

  write_output(quote(table(te_aff_expand_young$subfamily)), "Young TE subfamily counts")
  write_output(quote(nrow(te_aff_expand_young)), "Number of young TE insertions")

  if (nrow(te_aff_expand_young) > 0) {
    # Process combinations and add nohit samples
    te_young_processed <- process_all_combinations(te_aff_expand_young)
    if (nrow(nohits) > 0) {
      te_young_processed <- as.data.frame(add_nohit_samples(te_young_processed, nohits))
    }
    te_young_processed <- merge_dfs(te_young_processed, clinical, include_all_x = FALSE)

    # Merge with te_aff to get sequencing metrics and tumor_type_grouped
    if (length(metrics_cols_available) > 1) {
      te_young_processed <- te_young_processed %>%
        left_join(te_aff %>% select(all_of(metrics_cols_available)) %>% distinct(), by = "sample")
    }

    # Fix tumor_type column naming after merge (may have .x/.y suffixes)
    if ("tumor_type.x" %in% colnames(te_young_processed) && !"tumor_type" %in% colnames(te_young_processed)) {
      te_young_processed <- te_young_processed %>% rename(tumor_type = tumor_type.x)
    }
    if ("tumor_type_grouped.x" %in% colnames(te_young_processed) && !"tumor_type_grouped" %in% colnames(te_young_processed)) {
      te_young_processed <- te_young_processed %>% rename(tumor_type_grouped = tumor_type_grouped.x)
    }

    # If tumor_type_grouped still doesn't exist, create it with higher threshold for fewer groups
    if (!"tumor_type_grouped" %in% colnames(te_young_processed) && "tumor_type" %in% colnames(te_young_processed)) {
      te_young_processed <- te_young_processed %>%
        mutate(tumor_type_grouped = ifelse(
          tumor_type %in% names(table(tumor_type))[table(tumor_type) >= 10],
          tumor_type, "Other"
        ))
    }

    cat("\n--- Young TE by TP53 Status ---\n")

    # Kruskal test
    cat("Plotting young TE count by TP53 status (Kruskal test)...\n")
    p_young_tp53 <- plot_count_kruskal(te_young_processed, type=NA, chr=NA, group="TP53_status",
                                        log_scale=TRUE, x_lab=x_tp53, y_lab="Young TE count")
    titled_print(p_young_tp53, "Young TE by TP53 status")
    ggsave(paste0(output_dir, "young_te_tp53_kruskal.png"), plot=p_young_tp53, width = 4, height = 5)

    # Linear model
    cat("Plotting young TE count by TP53 status (Linear model)...\n")
    tryCatch({
      p_young_tp53_lm <- plot_count_lm(te_young_processed, type=NA, chr=NA, group="TP53_status",
                                        y_lab="Young TE count", covariates=covar_full,
                                        residuals=FALSE, log_scale=FALSE, x_lab=x_tp53, min_samples=3)
      titled_print(p_young_tp53_lm, "Young TE by TP53 status (LM)")
      ggsave(paste0(output_dir, "young_te_tp53_lm.png"), plot=p_young_tp53_lm, width = 4, height = 5)
    }, error = function(e) {
      cat("Warning: Could not create LM plot for TP53 status:", e$message, "\n")
    })

    cat("\n--- Young TE by Ancestry ---\n")

    if ("predicted_ancestry_thres" %in% colnames(te_young_processed)) {
      te_young_ancestry <- te_young_processed %>%
        filter(!is.na(predicted_ancestry_thres) & predicted_ancestry_thres != "Unknown" & predicted_ancestry_thres != "Undefined")

      if (nrow(te_young_ancestry) > 0 && length(unique(te_young_ancestry$predicted_ancestry_thres)) > 1) {
        # Kruskal test
        cat("Plotting young TE count by ancestry (Kruskal test)...\n")
        p_young_ancestry <- plot_count_kruskal(te_young_ancestry, type=NA, chr=NA, group="predicted_ancestry_thres",
                                                log_scale=TRUE, x_lab="Ancestry", y_lab="Young TE count")
        titled_print(p_young_ancestry, "Young TE by ancestry")
        ggsave(paste0(output_dir, "young_te_ancestry_kruskal.png"), plot=p_young_ancestry, width = 6, height = 5)

        # Linear model
        cat("Plotting young TE count by ancestry (Linear model)...\n")
        tryCatch({
          p_young_ancestry_lm <- plot_count_lm(te_young_ancestry, type=NA, chr=NA, group="predicted_ancestry_thres",
                                                y_lab="Young TE count", covariates=covar_no_ancestry,
                                                residuals=FALSE, log_scale=FALSE, x_lab="Ancestry", min_samples=3,
                                                fill_palette=color_palette_7)
          titled_print(p_young_ancestry_lm, "Young TE by ancestry (LM)")
          ggsave(paste0(output_dir, "young_te_ancestry_lm.png"), plot=p_young_ancestry_lm, width = 6, height = 5)
        }, error = function(e) {
          cat("Warning: Could not create LM plot for ancestry:", e$message, "\n")
        })
      }
    }

    cat("\n--- Young TE by Tumour Type ---\n")

    if ("tumor_type_grouped" %in% colnames(te_young_processed)) {
      te_young_tt <- te_young_processed %>% filter(!is.na(tumor_type_grouped))

      if (nrow(te_young_tt) > 0 && length(unique(te_young_tt$tumor_type_grouped)) > 1) {
        # Kruskal test
        cat("Plotting young TE count by tumour type (Kruskal test)...\n")
        tryCatch({
          p_young_tt <- plot_count_kruskal(te_young_tt, type=NA, chr=NA, group="tumor_type_grouped",
                                            log_scale=TRUE, x_lab="Tumour type", y_lab="Young TE count")
          titled_print(p_young_tt, "Young TE by tumour type")
          ggsave(paste0(output_dir, "young_te_tumour_type_kruskal.png"), plot=p_young_tt, width = 8, height = 5)
        }, error = function(e) {
          cat("Warning: Could not create Kruskal plot for tumour type:", e$message, "\n")
        })

        # Linear model
        cat("Plotting young TE count by tumour type (Linear model)...\n")
        tryCatch({
          p_young_tt_lm <- plot_count_lm(te_young_tt, type=NA, chr=NA, group="tumor_type_grouped",
                                          y_lab="Young TE count", covariates=covar_no_tumour_type,
                                          residuals=FALSE, log_scale=FALSE, x_lab="Tumour type", min_samples=3)
          titled_print(p_young_tt_lm, "Young TE by tumour type (LM)")
          ggsave(paste0(output_dir, "young_te_tumour_type_lm.png"), plot=p_young_tt_lm, width = 8, height = 5)
        }, error = function(e) {
          cat("Warning: Could not create LM plot for tumour type:", e$message, "\n")
        })
      }
    }

    # Summary statistics for young TE
    cat("\n--- Young TE Summary Statistics ---\n")
    young_summary <- te_young_processed %>%
      group_by(TP53_status) %>%
      summarise(
        n_samples = n(),
        mean_count = mean(total, na.rm = TRUE),
        median_count = median(total, na.rm = TRUE),
        sd_count = sd(total, na.rm = TRUE),
        .groups = "drop"
      )
    write_output(quote(young_summary), "Young TE summary by TP53 status")

    # Breakdown by subfamily
    cat("\n--- Young TE Breakdown by Subfamily ---\n")
    young_by_subfamily <- te_aff_expand_young %>%
      group_by(subfamily) %>%
      summarise(
        n_insertions = n(),
        n_samples = n_distinct(sample),
        n_TP53_mut = sum(TP53_status == "Mutant", na.rm = TRUE),
        n_TP53_wt = sum(TP53_status == "WT", na.rm = TRUE),
        .groups = "drop"
      ) %>%
      arrange(desc(n_insertions))
    write_output(quote(young_by_subfamily), "Young TE counts by subfamily")

    # Save to CSV
    write.csv(young_by_subfamily, paste0(output_dir, "files/young_te_by_subfamily.csv"), row.names = FALSE)

    ############################################################################
    #### YOUNG VS OLD SUBFAMILY ANALYSIS BY ANCESTRY ####
    ############################################################################
    cat("\n--- Young vs Old Subfamily Analysis by Ancestry ---\n")

    tryCatch({
      # Add young/old classification to all TEs
      te_aff_expand_classified <- te_aff_expand %>%
        mutate(age_class = ifelse(subfamily %in% young_family, "Young", "Old"))

      # Merge with clinical data to get ancestry
      te_aff_expand_classified <- te_aff_expand_classified %>%
        left_join(clinical %>% select(sample, predicted_ancestry_thres) %>% distinct(), by = "sample")

      # Filter for valid ancestry
      te_classified_ancestry <- te_aff_expand_classified %>%
        filter(!is.na(predicted_ancestry_thres) &
               predicted_ancestry_thres != "Unknown" &
               predicted_ancestry_thres != "Undefined")

      if (nrow(te_classified_ancestry) > 0) {

        # --- PIE CHARTS: Young vs Old by Ancestry ---
        cat("Creating pie charts for young vs old TEs by ancestry...\n")

        # Overall pie chart (all ancestries combined)
        overall_counts <- te_classified_ancestry %>%
          count(age_class) %>%
          mutate(pct = n / sum(n) * 100,
                 label = paste0(age_class, "\n", n, " (", round(pct, 1), "%)"))

        p_pie_overall <- ggplot(overall_counts, aes(x = "", y = n, fill = age_class)) +
          geom_bar(stat = "identity", width = 1) +
          coord_polar("y", start = 0) +
          scale_fill_manual(values = c("Young" = "#E41A1C", "Old" = "#377EB8")) +
          labs(title = "Young vs Old TEs - All Samples", fill = "Age Class") +
          theme_void() +
          theme(legend.position = "right") +
          geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 3)

        titled_print(p_pie_overall, "Young vs Old TEs - Overall")
        ggsave(paste0(output_dir, "young_old_pie_overall.png"), plot = p_pie_overall, width = 6, height = 5)

        # Pie chart for each ancestry
        ancestry_counts <- te_classified_ancestry %>%
          count(predicted_ancestry_thres, age_class) %>%
          group_by(predicted_ancestry_thres) %>%
          mutate(pct = n / sum(n) * 100,
                 label = paste0(round(pct, 1), "%"))

        p_pie_ancestry <- ggplot(ancestry_counts, aes(x = "", y = n, fill = age_class)) +
          geom_bar(stat = "identity", width = 1) +
          coord_polar("y", start = 0) +
          facet_wrap(~predicted_ancestry_thres, ncol = 3) +
          scale_fill_manual(values = c("Young" = "#E41A1C", "Old" = "#377EB8")) +
          labs(title = "Young vs Old TEs by Ancestry", fill = "Age Class") +
          theme_void() +
          theme(legend.position = "bottom",
                strip.text = element_text(size = 10, face = "bold"))

        titled_print(p_pie_ancestry, "Young vs Old TEs by Ancestry")
        ggsave(paste0(output_dir, "young_old_pie_by_ancestry.png"), plot = p_pie_ancestry, width = 10, height = 8)

        # --- BOX PLOTS: Young vs Old counts by Ancestry ---
        cat("Creating box plots for young vs old TE counts by ancestry...\n")

        # Calculate counts per sample for young and old
        sample_counts <- te_classified_ancestry %>%
          count(sample, predicted_ancestry_thres, age_class) %>%
          tidyr::pivot_wider(names_from = age_class, values_from = n, values_fill = 0)

        # Add samples with no TEs
        all_samples_ancestry <- clinical %>%
          filter(!is.na(predicted_ancestry_thres) &
                 predicted_ancestry_thres != "Unknown" &
                 predicted_ancestry_thres != "Undefined") %>%
          select(sample, predicted_ancestry_thres) %>%
          distinct()

        sample_counts_full <- all_samples_ancestry %>%
          left_join(sample_counts, by = c("sample", "predicted_ancestry_thres")) %>%
          mutate(Young = ifelse(is.na(Young), 0, Young),
                 Old = ifelse(is.na(Old), 0, Old))

        # Pivot longer for plotting
        sample_counts_long <- sample_counts_full %>%
          tidyr::pivot_longer(cols = c(Young, Old), names_to = "age_class", values_to = "count")

        # Box plot
        p_box_ancestry <- ggplot(sample_counts_long, aes(x = predicted_ancestry_thres, y = count, fill = age_class)) +
          geom_boxplot(outlier.shape = 21) +
          scale_fill_manual(values = c("Young" = "#E41A1C", "Old" = "#377EB8")) +
          scale_y_continuous(trans = scales::pseudo_log_trans(base = 10)) +
          labs(x = "Ancestry", y = "TE count per sample", fill = "Age Class",
               title = "Young vs Old TE counts by Ancestry") +
          theme_bw() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))

        titled_print(p_box_ancestry, "Young vs Old TE counts by Ancestry")
        ggsave(paste0(output_dir, "young_old_boxplot_by_ancestry.png"), plot = p_box_ancestry, width = 8, height = 6)

        # --- PIE CHARTS: All Subfamilies by Ancestry ---
        cat("Creating pie charts for all subfamilies by ancestry...\n")

        # Overall subfamily distribution
        subfamily_overall <- te_classified_ancestry %>%
          count(subfamily) %>%
          arrange(desc(n)) %>%
          mutate(pct = n / sum(n) * 100)

        # Top 10 subfamilies for cleaner visualization
        top_subfamilies <- subfamily_overall %>% head(10) %>% pull(subfamily)

        subfamily_data <- te_classified_ancestry %>%
          mutate(subfamily_grouped = ifelse(subfamily %in% top_subfamilies, subfamily, "Other")) %>%
          count(subfamily_grouped) %>%
          mutate(pct = n / sum(n) * 100,
                 label = ifelse(pct >= 3, paste0(round(pct, 1), "%"), ""))

        p_pie_subfamily_overall <- ggplot(subfamily_data, aes(x = "", y = n, fill = subfamily_grouped)) +
          geom_bar(stat = "identity", width = 1) +
          coord_polar("y", start = 0) +
          labs(title = "Subfamily Distribution - All Samples", fill = "Subfamily") +
          theme_void() +
          theme(legend.position = "right")

        titled_print(p_pie_subfamily_overall, "Subfamily Distribution - Overall")
        ggsave(paste0(output_dir, "subfamily_pie_overall.png"), plot = p_pie_subfamily_overall, width = 8, height = 6)

        # Subfamily distribution by ancestry
        subfamily_by_ancestry <- te_classified_ancestry %>%
          mutate(subfamily_grouped = ifelse(subfamily %in% top_subfamilies, subfamily, "Other")) %>%
          count(predicted_ancestry_thres, subfamily_grouped) %>%
          group_by(predicted_ancestry_thres) %>%
          mutate(pct = n / sum(n) * 100)

        p_pie_subfamily_ancestry <- ggplot(subfamily_by_ancestry, aes(x = "", y = n, fill = subfamily_grouped)) +
          geom_bar(stat = "identity", width = 1) +
          coord_polar("y", start = 0) +
          facet_wrap(~predicted_ancestry_thres, ncol = 3) +
          labs(title = "Subfamily Distribution by Ancestry", fill = "Subfamily") +
          theme_void() +
          theme(legend.position = "bottom",
                strip.text = element_text(size = 10, face = "bold"))

        titled_print(p_pie_subfamily_ancestry, "Subfamily Distribution by Ancestry")
        ggsave(paste0(output_dir, "subfamily_pie_by_ancestry.png"), plot = p_pie_subfamily_ancestry, width = 12, height = 10)

        # --- BOX PLOT: Top subfamilies by Ancestry ---
        cat("Creating box plot for top subfamilies by ancestry...\n")

        # Calculate counts per sample for top subfamilies
        subfamily_sample_counts <- te_classified_ancestry %>%
          filter(subfamily %in% top_subfamilies) %>%
          count(sample, predicted_ancestry_thres, subfamily)

        # Add zeros for missing combinations
        subfamily_sample_full <- expand.grid(
          sample = unique(all_samples_ancestry$sample),
          subfamily = top_subfamilies,
          stringsAsFactors = FALSE
        ) %>%
          left_join(all_samples_ancestry, by = "sample") %>%
          left_join(subfamily_sample_counts, by = c("sample", "predicted_ancestry_thres", "subfamily")) %>%
          mutate(n = ifelse(is.na(n), 0, n))

        p_box_subfamily <- ggplot(subfamily_sample_full, aes(x = predicted_ancestry_thres, y = n, fill = subfamily)) +
          geom_boxplot(outlier.shape = 21) +
          scale_y_continuous(trans = scales::pseudo_log_trans(base = 10)) +
          labs(x = "Ancestry", y = "TE count per sample", fill = "Subfamily",
               title = "Top Subfamilies by Ancestry") +
          theme_bw() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1),
                legend.position = "bottom") +
          guides(fill = guide_legend(nrow = 2))

        titled_print(p_box_subfamily, "Top Subfamilies by Ancestry")
        ggsave(paste0(output_dir, "subfamily_boxplot_by_ancestry.png"), plot = p_box_subfamily, width = 12, height = 8)

        # Save summary statistics
        young_old_summary <- sample_counts_full %>%
          group_by(predicted_ancestry_thres) %>%
          summarise(
            n_samples = n(),
            mean_young = mean(Young),
            median_young = median(Young),
            mean_old = mean(Old),
            median_old = median(Old),
            young_pct = sum(Young) / (sum(Young) + sum(Old)) * 100,
            .groups = "drop"
          )
        write_output(quote(young_old_summary), "Young vs Old TE summary by ancestry")
        write.csv(young_old_summary, paste0(output_dir, "files/young_old_by_ancestry.csv"), row.names = FALSE)

      } else {
        cat("No TEs with valid ancestry information found\n")
      }

    }, error = function(e) {
      cat("Warning: Could not perform young vs old analysis:", e$message, "\n")
    })

  } else {
    cat("No young TE insertions found (L1HS, AluY, SVA_E, SVA_F)\n")
  }

} else {
  cat("Warning: 'subfamily' column not found in data - skipping young TE analysis\n")
}

################################################################################
#### SOURCE/TRANSDUCTION ANALYSIS ####
################################################################################
cat("\n===== SOURCE/TRANSDUCTION ANALYSIS =====\n")

tryCatch({
  te_aff_expand_line <- extract_info_fields(te_aff_expand_line)

  cat("Number of unique LINE1 sources (affected):", length(unique(te_aff_expand_line$source)), "\n")
  cat("Number of LINE1 insertions (affected):", nrow(te_aff_expand_line), "\n")

  # Sources with multiple occurrences
  sources_with_multiple <- table(te_aff_expand_line$source)[table(te_aff_expand_line$source) > 1]
  if (length(sources_with_multiple) > 0) {
    cat("\nLINE1 sources with >1 occurrence (affected):\n")
    print(sources_with_multiple)
  } else {
    cat("\nNo LINE1 sources with >1 occurrence\n")
  }

  # Filter to transductions only
  sources_transduction <- te_aff_expand_line %>% filter(source != "not_transduction")

  # Add tumor_type_grouped to sources_transduction (get from te_aff if available)
  if (nrow(sources_transduction) > 0 && "tumor_type_grouped" %in% colnames(te_aff)) {
    sources_transduction <- sources_transduction %>%
      left_join(te_aff %>% select(sample, tumor_type_grouped) %>% distinct(), by = "sample")
    # Handle .x/.y suffixes if tumor_type_grouped already existed
    if ("tumor_type_grouped.y" %in% colnames(sources_transduction)) {
      sources_transduction <- sources_transduction %>%
        mutate(tumor_type_grouped = coalesce(tumor_type_grouped.y, tumor_type_grouped.x)) %>%
        select(-tumor_type_grouped.x, -tumor_type_grouped.y)
    }
  } else if (nrow(sources_transduction) > 0 && "tumor_type" %in% colnames(sources_transduction)) {
    # Fallback: create with higher threshold
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

      # Plot transduction count by TP53 status
      cat("\n--- Transduction Count by TP53 Status ---\n")
      transduction_per_sample <- sources_transduction %>%
        group_by(sample, TP53_status) %>%
        summarise(n_transductions = n(), .groups = "drop")

      # Add samples with 0 transductions
      all_samples_tp53 <- te_aff %>% select(sample, TP53_status) %>% distinct()
      transduction_per_sample <- all_samples_tp53 %>%
        left_join(transduction_per_sample %>% select(sample, n_transductions), by = "sample") %>%
        mutate(n_transductions = ifelse(is.na(n_transductions), 0, n_transductions))

      # Merge with clinical data for covariates
      transduction_per_sample <- merge_dfs(transduction_per_sample, clinical, include_all_x = FALSE)

      # Merge with te_aff to get sequencing metrics (med_cov, total_reads, etc.)
      if (length(metrics_cols_available) > 1) {
        transduction_per_sample <- transduction_per_sample %>%
          left_join(te_aff %>% select(all_of(metrics_cols_available)) %>% distinct(), by = "sample")
      }

      # Fix tumor_type column naming after merge (may have .x/.y suffixes)
      if ("tumor_type.x" %in% colnames(transduction_per_sample) && !"tumor_type" %in% colnames(transduction_per_sample)) {
        transduction_per_sample <- transduction_per_sample %>% rename(tumor_type = tumor_type.x)
      }
      if ("tumor_type_grouped.x" %in% colnames(transduction_per_sample) && !"tumor_type_grouped" %in% colnames(transduction_per_sample)) {
        transduction_per_sample <- transduction_per_sample %>% rename(tumor_type_grouped = tumor_type_grouped.x)
      }

      # If tumor_type_grouped still doesn't exist, create it with higher threshold for fewer groups
      if (!"tumor_type_grouped" %in% colnames(transduction_per_sample) && "tumor_type" %in% colnames(transduction_per_sample)) {
        transduction_per_sample <- transduction_per_sample %>%
          mutate(tumor_type_grouped = ifelse(
            tumor_type %in% names(table(tumor_type))[table(tumor_type) >= 10],
            tumor_type, "Other"
          ))
      }

      # Rename count column to 'total' for function compatibility
      transduction_per_sample <- transduction_per_sample %>% rename(total = n_transductions)

      # Kruskal test using function
      cat("Plotting transduction count by TP53 status (Kruskal test)...\n")
      p_trans_tp53 <- plot_count_kruskal(transduction_per_sample, type=NA, chr=NA, group="TP53_status",
                                          log_scale=TRUE, x_lab=x_tp53, y_lab="Transduction count")
      titled_print(p_trans_tp53, "Transduction count by TP53 status (Kruskal)")
      ggsave(paste0(output_dir, "source_transduction_tp53_kruskal.png"), plot = p_trans_tp53, width = 4, height = 5)

      # Linear model
      cat("Plotting transduction count by TP53 status (Linear model)...\n")
      tryCatch({
        p_trans_tp53_lm <- plot_count_lm(transduction_per_sample, type=NA, chr=NA, group="TP53_status",
                                          y_lab="Transduction count", covariates=covar_full,
                                          residuals=FALSE, log_scale=FALSE, x_lab=x_tp53, min_samples=3)
        titled_print(p_trans_tp53_lm, "Transduction count by TP53 status (LM)")
        ggsave(paste0(output_dir, "source_transduction_tp53_lm.png"), plot=p_trans_tp53_lm, width = 4, height = 5)
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

      # Top sources by TP53 status
      top_sources_tp53 <- sources_transduction %>%
        group_by(TP53_status, source) %>%
        summarise(n = n(), .groups = "drop") %>%
        group_by(TP53_status) %>%
        slice_max(n, n = 10) %>%
        ungroup()

      if (nrow(top_sources_tp53) > 0) {
        p_top_tp53 <- ggplot(top_sources_tp53, aes(x = reorder(source, n), y = n, fill = TP53_status)) +
          geom_bar(stat = "identity") +
          coord_flip() +
          facet_wrap(~TP53_status, scales = "free_y") +
          theme_minimal() +
          labs(x = "Source Element", y = "Number of Transductions",
               title = "Top LINE1 Sources by TP53 Status") +
          theme(axis.text.y = element_text(size = 7), legend.position = "none") +
          scale_fill_manual(values = c("Mutant" = "#E41A1C", "WT" = "#377EB8"))
        titled_print(p_top_tp53, "Top sources by TP53 status")
        ggsave(paste0(output_dir, "source_top_by_tp53.png"), plot = p_top_tp53, width = 10, height = 6)
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
  Analysis = c("Full-length LINE1 (>=5900bp)", "Young TE (L1HS, AluY, SVA_E, SVA_F)", "Transductions"),
  N_insertions = c(
    if (exists("te_aff_expand_line_fulllength")) nrow(te_aff_expand_line_fulllength) else NA,
    if (exists("te_aff_expand_young")) nrow(te_aff_expand_young) else NA,
    if (exists("sources_transduction")) nrow(sources_transduction) else NA
  )
)
write_output(quote(summary_df), "Summary of LINE1 subset analyses")

cat("\n Full-length, Young TE, and Source analysis completed successfully\n")

# Close module-specific sink
close_module_sink()
