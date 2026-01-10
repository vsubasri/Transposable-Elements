#!/usr/bin/env Rscript

# Germline TE Visualization - Subfamily Analysis
# Subfamily distributions, Fisher's tests by TP53 status, ancestry, and tumour type

# Source common setup and load data
source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/00_viz_common_setup.R")
REQUIRED_DATA <- c("expand", "clinical", "ancestry")
source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/00_viz_load_data_germline.R")

# Initialize module-specific text output
dir.create(paste0(plot_dir, "subfamily/"), showWarnings = FALSE, recursive = TRUE)
init_module_sink(paste0(plot_dir, "subfamily/"), "SUBFAMILY")

cat("Running 02_te_viz_germline_12_subfamily.R...\n")

#### SUBFAMILY PIE CHARTS ####
cat("\n=== Subfamily Pie Charts ===\n")

# Generate pie charts for each TE type using te_all_expand
for (te_type in c("L1", "Alu", "SVA")) {
  cat(paste0("Plotting: subfamily_pie_", te_type, ".pdf\n"))
  tryCatch({
    p <- plot_subfamily_pie(te_all_expand, group = te_type, n = 3)
    p <- p + ggtitle(paste0(te_type, " Subfamily Distribution"))
    print(p)
    ggsave(paste0(plot_dir, "subfamily/subfamily_pie_", te_type, ".pdf"), plot = p, width = 8, height = 6)
  }, error = function(e) {
    cat(paste0("  Warning: Could not generate pie chart for ", te_type, ": ", e$message, "\n"))
  })
}

#### FISHER'S TEST BY TP53 STATUS ####
cat("\n=== Fisher's Test: Subfamily by TP53 Status ===\n")

# Use te_lfs_mut_wt_expand which has both Mutant and WT samples
if ("subfamily" %in% colnames(te_lfs_mut_wt_expand) && "TP53_status" %in% colnames(te_lfs_mut_wt_expand)) {

  fisher_results <- perform_fisher_test_summary(te_lfs_mut_wt_expand, "subfamily", "TP53_status")

  # Add FDR correction
  fisher_results$P_Adjusted <- p.adjust(fisher_results$P_Value, method = "BH")

  # Sort by p-value
  fisher_results <- fisher_results[order(fisher_results$P_Value), ]

  # Print results
  cat("\nFisher's Test Results (Top 20 by p-value):\n")
  print(head(fisher_results, 20), row.names = FALSE)

  # Save full results
  write.csv(fisher_results, paste0(plot_dir, "subfamily/files/fisher_subfamily_tp53.csv"), row.names = FALSE)
  cat(paste0("\nSaved: ", plot_dir, "subfamily/files/fisher_subfamily_tp53.csv\n"))

  # Plot significant subfamilies
  sig_subfamilies <- fisher_results[fisher_results$P_Adjusted < 0.05, ]
  if (nrow(sig_subfamilies) > 0) {
    cat(paste0("\nSignificant subfamilies (FDR < 0.05): ", nrow(sig_subfamilies), "\n"))
    print(sig_subfamilies, row.names = FALSE)

    # Create bar plot of significant subfamilies
    sig_long <- sig_subfamilies %>%
      select(Subfamily, Mutant_Percentage, WT_Percentage, P_Adjusted) %>%
      pivot_longer(cols = c(Mutant_Percentage, WT_Percentage),
                   names_to = "TP53_Status", values_to = "Percentage") %>%
      mutate(TP53_Status = gsub("_Percentage", "", TP53_Status))

    p_sig <- ggplot(sig_long, aes(x = reorder(Subfamily, -Percentage), y = Percentage, fill = TP53_Status)) +
      geom_bar(stat = "identity", position = "dodge") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(x = "Subfamily", y = "Percentage (%)",
           title = "Significant Subfamily Differences by TP53 Status",
           subtitle = "FDR < 0.05") +
      scale_fill_manual(values = c("Mutant" = "#E41A1C", "WT" = "#377EB8"))

    print(p_sig)
    ggsave(paste0(plot_dir, "subfamily/subfamily_fisher_tp53_significant.pdf"), plot = p_sig, width = 10, height = 6)
  } else {
    cat("\nNo subfamilies significant at FDR < 0.05\n")
  }

  # Create volcano-style plot
  fisher_results$neg_log_p <- -log10(fisher_results$P_Value)
  fisher_results$diff <- fisher_results$Mutant_Percentage - fisher_results$WT_Percentage

  p_volcano <- ggplot(fisher_results, aes(x = diff, y = neg_log_p)) +
    geom_point(aes(color = P_Adjusted < 0.05), alpha = 0.6) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    geom_text_repel(data = subset(fisher_results, P_Adjusted < 0.1),
                    aes(label = Subfamily), size = 3, max.overlaps = 15) +
    theme_minimal() +
    labs(x = "Difference (Mutant% - WT%)", y = "-log10(p-value)",
         title = "Subfamily Association with TP53 Status",
         color = "FDR < 0.05") +
    scale_color_manual(values = c("FALSE" = "grey50", "TRUE" = "#E41A1C"))

  print(p_volcano)
  ggsave(paste0(plot_dir, "subfamily/subfamily_fisher_tp53_volcano.pdf"), plot = p_volcano, width = 10, height = 8)

} else {
  cat("Warning: Required columns 'subfamily' or 'TP53_status' not found in te_lfs_mut_wt_expand\n")
}

#### SUBFAMILY BY ANCESTRY ####
cat("\n=== Subfamily Distribution by Ancestry ===\n")

if ("subfamily" %in% colnames(te_all_expand) && "predicted_ancestry_thres" %in% colnames(te_all_expand)) {

  # Count subfamilies by ancestry
  subfamily_ancestry <- te_all_expand %>%
    filter(!is.na(predicted_ancestry_thres) & predicted_ancestry_thres != "Unknown") %>%
    group_by(predicted_ancestry_thres, subfamily) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(predicted_ancestry_thres) %>%
    mutate(percentage = count / sum(count) * 100) %>%
    ungroup()

  # Get top subfamilies overall
  top_subfamilies <- subfamily_ancestry %>%
    group_by(subfamily) %>%
    summarise(total = sum(count), .groups = "drop") %>%
    arrange(desc(total)) %>%
    head(15) %>%
    pull(subfamily)

  # Filter to top subfamilies and lump rest as "Other"
  subfamily_ancestry_plot <- subfamily_ancestry %>%
    mutate(subfamily_grouped = ifelse(subfamily %in% top_subfamilies, subfamily, "Other")) %>%
    group_by(predicted_ancestry_thres, subfamily_grouped) %>%
    summarise(count = sum(count), percentage = sum(percentage), .groups = "drop")

  # Stacked bar plot
  p_ancestry_stack <- ggplot(subfamily_ancestry_plot,
                              aes(x = predicted_ancestry_thres, y = percentage, fill = subfamily_grouped)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "right") +
    labs(x = "Ancestry", y = "Percentage (%)",
         title = "Subfamily Distribution by Ancestry",
         fill = "Subfamily") +
    scale_fill_manual(values = c(colorRampPalette(brewer.pal(12, "Set3"))(length(unique(subfamily_ancestry_plot$subfamily_grouped)))))

  cat("Plotting: subfamily_by_ancestry_stacked.pdf\n")
  print(p_ancestry_stack)
  ggsave(paste0(plot_dir, "subfamily/subfamily_by_ancestry_stacked.pdf"), plot = p_ancestry_stack, width = 12, height = 8)

  # Faceted bar plot by TE type
  for (te_type in c("L1", "Alu", "SVA")) {
    pattern <- switch(te_type,
                      "L1" = "^L1",
                      "Alu" = "^Alu|^FRAM",
                      "SVA" = "^SVA")

    subfamily_ancestry_type <- te_all_expand %>%
      filter(!is.na(predicted_ancestry_thres) & predicted_ancestry_thres != "Unknown") %>%
      filter(grepl(pattern, subfamily)) %>%
      group_by(predicted_ancestry_thres, subfamily) %>%
      summarise(count = n(), .groups = "drop") %>%
      group_by(predicted_ancestry_thres) %>%
      mutate(percentage = count / sum(count) * 100) %>%
      ungroup()

    if (nrow(subfamily_ancestry_type) > 0) {
      # Get top subfamilies for this type
      top_type <- subfamily_ancestry_type %>%
        group_by(subfamily) %>%
        summarise(total = sum(count), .groups = "drop") %>%
        arrange(desc(total)) %>%
        head(10) %>%
        pull(subfamily)

      subfamily_ancestry_type_plot <- subfamily_ancestry_type %>%
        mutate(subfamily_grouped = ifelse(subfamily %in% top_type, subfamily, "Other")) %>%
        group_by(predicted_ancestry_thres, subfamily_grouped) %>%
        summarise(count = sum(count), percentage = sum(percentage), .groups = "drop")

      p_type <- ggplot(subfamily_ancestry_type_plot,
                       aes(x = predicted_ancestry_thres, y = percentage, fill = subfamily_grouped)) +
        geom_bar(stat = "identity") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(x = "Ancestry", y = "Percentage (%)",
             title = paste0(te_type, " Subfamily Distribution by Ancestry"),
             fill = "Subfamily")

      cat(paste0("Plotting: subfamily_", te_type, "_by_ancestry.pdf\n"))
      print(p_type)
      ggsave(paste0(plot_dir, "subfamily/subfamily_", te_type, "_by_ancestry.pdf"), plot = p_type, width = 10, height = 6)
    }
  }

  # Chi-square test for subfamily distribution differences across ancestries
  cat("\nChi-square test for subfamily distribution by ancestry:\n")
  contingency_table <- table(te_all_expand$subfamily, te_all_expand$predicted_ancestry_thres)
  # Filter to subfamilies with enough counts
  contingency_table <- contingency_table[rowSums(contingency_table) >= 10, ]
  if (nrow(contingency_table) > 1 && ncol(contingency_table) > 1) {
    chi_test <- chisq.test(contingency_table)
    cat(paste0("Chi-square statistic: ", round(chi_test$statistic, 2), "\n"))
    cat(paste0("Degrees of freedom: ", chi_test$parameter, "\n"))
    cat(paste0("P-value: ", formatC(chi_test$p.value, format = "e", digits = 2), "\n"))
  }

} else {
  cat("Warning: Required columns not found for ancestry analysis\n")
}

#### SUBFAMILY BY TUMOUR TYPE ####
cat("\n=== Subfamily Distribution by Tumour Type ===\n")

if ("subfamily" %in% colnames(te_all_expand) && "tumor_type_grouped" %in% colnames(te_all_expand)) {

  # Count subfamilies by tumour type
  subfamily_tumour <- te_all_expand %>%
    filter(!is.na(tumor_type_grouped)) %>%
    group_by(tumor_type_grouped, subfamily) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(tumor_type_grouped) %>%
    mutate(percentage = count / sum(count) * 100) %>%
    ungroup()

  # Get top subfamilies overall
  top_subfamilies_tt <- subfamily_tumour %>%
    group_by(subfamily) %>%
    summarise(total = sum(count), .groups = "drop") %>%
    arrange(desc(total)) %>%
    head(15) %>%
    pull(subfamily)

  # Filter to top subfamilies and lump rest as "Other"
  subfamily_tumour_plot <- subfamily_tumour %>%
    mutate(subfamily_grouped = ifelse(subfamily %in% top_subfamilies_tt, subfamily, "Other")) %>%
    group_by(tumor_type_grouped, subfamily_grouped) %>%
    summarise(count = sum(count), percentage = sum(percentage), .groups = "drop")

  # Stacked bar plot
  p_tumour_stack <- ggplot(subfamily_tumour_plot,
                           aes(x = tumor_type_grouped, y = percentage, fill = subfamily_grouped)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "right") +
    labs(x = "Tumour Type", y = "Percentage (%)",
         title = "Subfamily Distribution by Tumour Type",
         fill = "Subfamily") +
    scale_fill_manual(values = c(colorRampPalette(brewer.pal(12, "Set3"))(length(unique(subfamily_tumour_plot$subfamily_grouped)))))

  cat("Plotting: subfamily_by_tumour_type_stacked.pdf\n")
  print(p_tumour_stack)
  ggsave(paste0(plot_dir, "subfamily/subfamily_by_tumour_type_stacked.pdf"), plot = p_tumour_stack, width = 12, height = 8)

  # Faceted bar plot by TE type
  for (te_type in c("L1", "Alu", "SVA")) {
    pattern <- switch(te_type,
                      "L1" = "^L1",
                      "Alu" = "^Alu|^FRAM",
                      "SVA" = "^SVA")

    subfamily_tumour_type <- te_all_expand %>%
      filter(!is.na(tumor_type_grouped)) %>%
      filter(grepl(pattern, subfamily)) %>%
      group_by(tumor_type_grouped, subfamily) %>%
      summarise(count = n(), .groups = "drop") %>%
      group_by(tumor_type_grouped) %>%
      mutate(percentage = count / sum(count) * 100) %>%
      ungroup()

    if (nrow(subfamily_tumour_type) > 0) {
      # Get top subfamilies for this type
      top_type_tt <- subfamily_tumour_type %>%
        group_by(subfamily) %>%
        summarise(total = sum(count), .groups = "drop") %>%
        arrange(desc(total)) %>%
        head(10) %>%
        pull(subfamily)

      subfamily_tumour_type_plot <- subfamily_tumour_type %>%
        mutate(subfamily_grouped = ifelse(subfamily %in% top_type_tt, subfamily, "Other")) %>%
        group_by(tumor_type_grouped, subfamily_grouped) %>%
        summarise(count = sum(count), percentage = sum(percentage), .groups = "drop")

      p_type_tt <- ggplot(subfamily_tumour_type_plot,
                          aes(x = tumor_type_grouped, y = percentage, fill = subfamily_grouped)) +
        geom_bar(stat = "identity") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(x = "Tumour Type", y = "Percentage (%)",
             title = paste0(te_type, " Subfamily Distribution by Tumour Type"),
             fill = "Subfamily")

      cat(paste0("Plotting: subfamily_", te_type, "_by_tumour_type.pdf\n"))
      print(p_type_tt)
      ggsave(paste0(plot_dir, "subfamily/subfamily_", te_type, "_by_tumour_type.pdf"), plot = p_type_tt, width = 10, height = 6)
    }
  }

  # Heatmap of subfamily proportions by tumour type
  cat("Plotting: subfamily_tumour_type_heatmap.pdf\n")

  # Create matrix for heatmap
  subfamily_matrix <- subfamily_tumour %>%
    filter(subfamily %in% top_subfamilies_tt) %>%
    select(tumor_type_grouped, subfamily, percentage) %>%
    pivot_wider(names_from = tumor_type_grouped, values_from = percentage, values_fill = 0) %>%
    column_to_rownames("subfamily") %>%
    as.matrix()

  if (nrow(subfamily_matrix) > 1 && ncol(subfamily_matrix) > 1) {
    pdf(paste0(plot_dir, "subfamily/subfamily_tumour_type_heatmap.pdf"), width = 10, height = 8)
    pheatmap(subfamily_matrix,
             main = "Subfamily Distribution by Tumour Type",
             color = colorRampPalette(c("white", "steelblue", "darkblue"))(100),
             cluster_rows = TRUE,
             cluster_cols = TRUE,
             scale = "none",
             fontsize = 10)
    dev.off()
  }

  # Chi-square test for subfamily distribution differences across tumour types
  cat("\nChi-square test for subfamily distribution by tumour type:\n")
  contingency_table_tt <- table(te_all_expand$subfamily, te_all_expand$tumor_type_grouped)
  # Filter to subfamilies with enough counts
  contingency_table_tt <- contingency_table_tt[rowSums(contingency_table_tt) >= 10, ]
  if (nrow(contingency_table_tt) > 1 && ncol(contingency_table_tt) > 1) {
    chi_test_tt <- chisq.test(contingency_table_tt)
    cat(paste0("Chi-square statistic: ", round(chi_test_tt$statistic, 2), "\n"))
    cat(paste0("Degrees of freedom: ", chi_test_tt$parameter, "\n"))
    cat(paste0("P-value: ", formatC(chi_test_tt$p.value, format = "e", digits = 2), "\n"))
  }

} else {
  cat("Warning: Required columns not found for tumour type analysis\n")
}

#### SUBFAMILY SUMMARY TABLE ####
cat("\n=== Subfamily Summary Statistics ===\n")

if ("subfamily" %in% colnames(te_all_expand)) {
  subfamily_summary <- te_all_expand %>%
    group_by(subfamily) %>%
    summarise(
      total_count = n(),
      n_samples = n_distinct(sample),
      mean_per_sample = n() / n_distinct(sample),
      .groups = "drop"
    ) %>%
    arrange(desc(total_count))

  cat("\nTop 20 Subfamilies by Total Count:\n")
  print(head(subfamily_summary, 20), row.names = FALSE)

  # Save full summary
  dir.create(paste0(plot_dir, "subfamily/files/"), showWarnings = FALSE, recursive = TRUE)
  write.csv(subfamily_summary, paste0(plot_dir, "subfamily/files/subfamily_summary.csv"), row.names = FALSE)
  cat(paste0("\nSaved: ", plot_dir, "subfamily/files/subfamily_summary.csv\n"))
}

cat("\n✓ Subfamily analysis completed successfully\n")

# Close module-specific sink
close_module_sink()
