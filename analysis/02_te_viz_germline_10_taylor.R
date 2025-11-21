#!/usr/bin/env Rscript

# Germline TE Visualization - Taylor Cohort
# Taylor cohort-specific analyses

# Source common setup and load data
source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/00_viz_common_setup.R")
source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/00_viz_load_data_germline.R")

cat("Running 02_te_viz_germline_10_taylor.R...\n")

#### TAYLOR COHORT DATA CHECK ####
if (!exists("te_taylor") || nrow(te_taylor) == 0) {
  cat("⚠ WARNING: Taylor cohort data not available or empty (te_taylor not found or has 0 rows).\n")
  cat("  Skipping Taylor cohort analysis.\n\n")
  cat("✓ Script completed successfully\n")
  quit(save = "no", status = 0)
}

# Define TE types for analysis
types <- c(NA, "LINE1", "ALU", "SVA")

# Load RE path
re_germline_path <- "/Users/briannelaverty/Documents/R_Malkin/TE/data/final/germline_annotSV_output.SV_RE_intersect.report"

#### COHORT SUMMARY ####
write_output(quote({
  cat("Taylor Cohort Summary:\n")
  cat("Total samples:", nrow(te_taylor), "\n")
  cat("\nTumor type distribution:\n")
  print(table(te_taylor$tumor_type))
  cat("\nTumor type subclass distribution:\n")
  print(table(te_taylor$tumor_type_subclass))
  cat("\nAge at diagnosis summary:\n")
  print(summary(te_taylor$age_at_diagnosis))
}), "Taylor Cohort Summary Statistics")

#### ANALYSIS 1: TE FREQUENCY BY TUMOR TYPE SUBCLASS (KRUSKAL-WALLIS) ####
cat("\n===== ANALYSIS 1: TE FREQUENCY BY TUMOR TYPE SUBCLASS =====\n")

for (i in seq_along(types)) {
  type_label <- ifelse(is.na(types[i]), "total", types[i])

  write_output(
    quote(plot_count_kruskal(df = te_taylor, chr = NA, type = types[i],
                            group = "tumor_type_subclass",
                            x_lab = "Tumor Type Subclass",
                            y_lab = "TE Count",
                            log_scale = FALSE)),
    paste0("Taylor: TE Count by Tumor Type Subclass (", type_label, ")")
  )

  p <- plot_count_kruskal(df = te_taylor, chr = NA, type = types[i],
                         group = "tumor_type_subclass",
                         x_lab = "Tumor Type Subclass",
                         y_lab = "TE Count",
                         log_scale = FALSE)
  titled_print(p, paste0("Taylor: TE Count by Tumor Type Subclass (", type_label, ")"))
  ggsave(paste0(plot_dir, "taylor/te_count_by_subclass_", type_label, ".png"),
         plot = p, width = 10, height = 6)
}

#### ANALYSIS 5: SPECIFIC TEs BY TUMOR TYPE SUBCLASS ####
cat("\n===== ANALYSIS 5: SPECIFIC TEs BY TUMOR TYPE SUBCLASS =====\n")

tryCatch({
  for (min_samples in c(3, 5)) {
    cat("\n--- Testing with min_samples =", min_samples, "---\n")
    te_specific_taylor_subclass <- test_specific_tes_by_group(
      te_expand = te_taylor_expand,
      te_count = te_taylor,
      group_column = "tumor_type_subclass",
      min_samples_with = min_samples,
      min_samples_without = min_samples,
      output_dir = r_dir_files,
      output_prefix = paste0("specific_tes_taylor_subclass_min", min_samples)
    )

    if (!is.null(te_specific_taylor_subclass) && nrow(te_specific_taylor_subclass$full_results) > 0) {
      cat("Found", nrow(te_specific_taylor_subclass$significant_tes),
          "significant TEs at min_samples =", min_samples, "\n")
      write_output(quote(head(te_specific_taylor_subclass$full_results, 10)),
                  paste0("Top 10 TEs by subclass (min_samples=", min_samples, ")"))
    }
  }
}, error = function(e) {
  cat("Warning: Could not perform specific TE testing by subclass:", e$message, "\n")
})

#### ANALYSIS 6: CANCER GENE ANALYSIS BY TUMOR TYPE ####
cat("\n===== ANALYSIS 6: CANCER GENE ANALYSIS BY TUMOR TYPE =====\n")

# Process cancer gene overlaps in Taylor cohort
tryCatch({
  # Count cancer genes affected per sample
  te_taylor_split_cancergenes <- te_taylor_split %>%
    filter(Gene_name %in% genes) %>%
    group_by(sample) %>%
    summarise(
      cancer_genes_affected = n_distinct(Gene_name),
      .groups = "drop"
    )

  # Merge with main data
  te_taylor_cancergenes <- te_taylor %>%
    left_join(te_taylor_split_cancergenes, by = "sample") %>%
    mutate(cancer_genes_affected = ifelse(is.na(cancer_genes_affected), 0, cancer_genes_affected))

  write_output(quote({
    cat("Cancer gene summary:\n")
    cat("Total unique cancer genes affected:",
        length(unique(te_taylor_split$Gene_name[te_taylor_split$Gene_name %in% genes])), "\n")
    cat("Samples with cancer gene overlaps:",
        sum(te_taylor_cancergenes$cancer_genes_affected > 0), "\n")
    cat("\nTop 10 most frequently affected cancer genes:\n")
    gene_freq <- table(te_taylor_split$Gene_name[te_taylor_split$Gene_name %in% genes])
    print(head(sort(gene_freq, decreasing = TRUE), 10))
  }), "Taylor: Cancer Gene Summary")

  # Save to file
  write.csv(te_taylor_cancergenes,
           paste0(r_dir_files, "taylor_cancer_genes_per_sample.csv"),
           row.names = FALSE)

  # Test by tumor type subclass
  for (i in seq_along(types)) {
    type_label <- ifelse(is.na(types[i]), "total", types[i])

    p <- plot_count_kruskal(df = te_taylor_cancergenes, chr = NA, type = NA,
                           group = "tumor_type_subclass",
                           column = "cancer_genes_affected",
                           x_lab = "Tumor Type Subclass",
                           y_lab = "Cancer Genes Affected",
                           log_scale = FALSE)
    titled_print(p, "Taylor: Cancer Genes Affected by Tumor Subclass")
    ggsave(paste0(plot_dir, "taylor/cancer_genes_by_subclass.png"),
           plot = p, width = 10, height = 6)
  }
}, error = function(e) {
  cat("Warning: Could not perform cancer gene analysis:", e$message, "\n")
})

#### ANALYSIS 7: AGE AT DIAGNOSIS CORRELATION ####
cat("\n===== ANALYSIS 7: AGE AT DIAGNOSIS CORRELATION =====\n")

write_output(quote(plot_count_age(df = te_taylor, type = "total", chr = NA,
                                 y_lab = "Total TE Count")),
            "Taylor: TE Count by Age at Diagnosis")
p_age <- plot_count_age(df = te_taylor, type = "total", chr = NA,
                       y_lab = "Total TE Count")
titled_print(p_age, "Taylor: TE Count by Age at Diagnosis")
ggsave(paste0(plot_dir, "taylor/te_count_by_age.png"),
       plot = p_age, width = 8, height = 6)

# Age correlation by tumor subclass (if enough samples per subclass)
tryCatch({
  # Create scatter plot colored by subclass
  p_age_subclass <- ggplot(te_taylor %>% filter(!is.na(age_at_diagnosis)),
                          aes(x = age_at_diagnosis, y = total,
                              color = tumor_type_subclass)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_smooth(method = "lm", se = TRUE, alpha = 0.2) +
    labs(x = "Age at Diagnosis (years)",
         y = "Total TE Count",
         color = "Tumor Subclass") +
    theme_bw() +
    theme(legend.position = "right")

  titled_print(p_age_subclass, "Taylor: TE Count vs Age by Tumor Subclass")
  ggsave(paste0(plot_dir, "taylor/te_count_by_age_subclass.png"),
         plot = p_age_subclass, width = 10, height = 6)
}, error = function(e) {
  cat("Warning: Could not create age correlation by subclass:", e$message, "\n")
})

#### ANALYSIS 10: FULL-LENGTH LINE1 ANALYSIS ####
cat("\n===== ANALYSIS 10: FULL-LENGTH LINE1 ANALYSIS =====\n")

tryCatch({
  # Identify full-length LINE1 (>6000 bp)
  te_taylor_expand_line_fulllength <- te_taylor_expand %>%
    filter(ALT == "LINE1" & SV_length > 6000)

  cat("Full-length LINE1 insertions (>6kb):", nrow(te_taylor_expand_line_fulllength), "\n")

  if (nrow(te_taylor_expand_line_fulllength) > 0) {
    # Count per sample
    te_taylor_line_fl_counts <- te_taylor_expand_line_fulllength %>%
      group_by(sample) %>%
      summarise(full_length_line1 = n(), .groups = "drop")

    # Merge with main data
    te_taylor_line_fl <- te_taylor %>%
      left_join(te_taylor_line_fl_counts, by = "sample") %>%
      mutate(full_length_line1 = ifelse(is.na(full_length_line1), 0, full_length_line1))

    # Test by tumor subclass
    p_fl <- plot_count_kruskal(df = te_taylor_line_fl, chr = NA, type = NA,
                               group = "tumor_type_subclass",
                               column = "full_length_line1",
                               x_lab = "Tumor Type Subclass",
                               y_lab = "Full-length LINE1 Count",
                               log_scale = FALSE)
    titled_print(p_fl, "Taylor: Full-length LINE1 by Tumor Subclass")
    ggsave(paste0(plot_dir, "taylor/fulllength_line1_by_subclass.png"),
           plot = p_fl, width = 10, height = 6)

    # Save full-length LINE1 data
    write.csv(te_taylor_expand_line_fulllength,
             paste0(r_dir_files, "taylor_fulllength_line1.csv"),
             row.names = FALSE)
  }
}, error = function(e) {
  cat("Warning: Could not perform full-length LINE1 analysis:", e$message, "\n")
})

#### ANALYSIS 11: CHROMOSOMAL DISTRIBUTION ####
cat("\n===== ANALYSIS 11: CHROMOSOMAL DISTRIBUTION =====\n")

tryCatch({
  # Count TEs per chromosome per tumor subclass
  chr_dist <- te_taylor_expand %>%
    group_by(tumor_type_subclass, SV_chrom) %>%
    summarise(te_count = n(), .groups = "drop") %>%
    pivot_wider(names_from = SV_chrom, values_from = te_count, values_fill = 0)

  write_output(quote(chr_dist), "Chromosomal Distribution by Tumor Subclass")

  # Save to file
  write.csv(chr_dist,
           paste0(r_dir_files, "taylor_chr_distribution_by_subclass.csv"),
           row.names = FALSE)

  # Chi-square test for non-random distribution (if enough samples)
  if (nrow(chr_dist) > 1) {
    chr_matrix <- as.matrix(chr_dist[, -1])
    rownames(chr_matrix) <- chr_dist$tumor_type_subclass

    chi_test <- tryCatch({
      chisq.test(chr_matrix)
    }, error = function(e) NULL)

    if (!is.null(chi_test)) {
      cat("\nChi-square test for chromosomal distribution:\n")
      cat("X-squared =", chi_test$statistic, ", p-value =", chi_test$p.value, "\n")
    }
  }
}, error = function(e) {
  cat("Warning: Could not perform chromosomal distribution analysis:", e$message, "\n")
})

#### ANALYSIS 2: RE (REGULATORY ELEMENT) ANALYSIS - PUT LAST ####
cat("\n===== ANALYSIS 2: RE ANALYSIS =====\n")

#### ANALYSIS 2a: RE BY TUMOR TYPE SUBCLASS (WITHIN TAYLOR) ####
cat("\n=== Analysis 2a: RE by Tumor Type Subclass (within Taylor) ===\n")

tryCatch({
  # Load and join RE data for Taylor cohort
  cat("Loading RE data for Taylor cohort...\n")
  te_taylor_re <- load_and_join_re_data(te_taylor_expand, re_germline_path)
  te_taylor_re_split <- split_re_genes(te_taylor_re)

  cat("Taylor RE data loaded:\n")
  cat("  Samples with RE overlaps:", length(unique(te_taylor_re_split$sample.x)), "\n")
  cat("  Genes in REs:", length(unique(te_taylor_re_split$gene_reg)), "\n")

  # Group genes by tumor_type_subclass
  geneClusters_taylor_subclass <- lapply(
    split(te_taylor_re_split$gene_reg, te_taylor_re_split$tumor_type_subclass),
    unique
  )

  cat("\nGenes in REs per tumor subclass:\n")
  for (subclass in names(geneClusters_taylor_subclass)) {
    cat("  ", subclass, ":", length(geneClusters_taylor_subclass[[subclass]]), "genes\n")
  }

  # Run compareCluster pathway analysis
  if (length(geneClusters_taylor_subclass) > 1 &&
      all(sapply(geneClusters_taylor_subclass, length) >= 5)) {

    ora_re_taylor_subclass <- tryCatch({
      compareCluster(
        geneCluster = geneClusters_taylor_subclass,
        fun = "enrichGO",
        OrgDb = org.Hs.eg.db,
        keyType = "SYMBOL",
        ont = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.1
      )
    }, error = function(e) {
      cat("Error in pathway analysis:", e$message, "\n")
      NULL
    })

    if (!is.null(ora_re_taylor_subclass) && nrow(as.data.frame(ora_re_taylor_subclass)) > 0) {
      cat("✓ Found", nrow(as.data.frame(ora_re_taylor_subclass)),
          "enriched pathways across tumor subclasses\n")

      # Save results
      write.csv(as.data.frame(ora_re_taylor_subclass),
               paste0(r_dir_files, "re_pathway_taylor_subclass.csv"),
               row.names = FALSE)

      # Create plots
      tryCatch({
        p_dot <- dotplot(ora_re_taylor_subclass, showCategory = 10)
        titled_print(p_dot, "RE Pathway Dotplot - Taylor by Tumor Subclass")
        ggsave(paste0(plot_dir, "taylor/re_pathway_subclass_dot.png"),
               plot = p_dot, width = 14, height = 9)
      }, error = function(e) cat("Warning: Could not create dotplot:", e$message, "\n"))

      tryCatch({
        p_cnet <- cnetplot(ora_re_taylor_subclass, showCategory = 5,
                          colorEdge = TRUE, node_label = "category")
        titled_print(p_cnet, "RE Pathway Cnetplot - Taylor by Tumor Subclass")
        ggsave(paste0(plot_dir, "taylor/re_pathway_subclass_cnet.png"),
               plot = p_cnet, width = 14, height = 9)
      }, error = function(e) cat("Warning: Could not create cnetplot:", e$message, "\n"))

      tryCatch({
        ora_re_taylor_pairwise <- pairwise_termsim(ora_re_taylor_subclass)
        p_emap <- emapplot(ora_re_taylor_pairwise, showCategory = 20)
        titled_print(p_emap, "RE Pathway Emapplot - Taylor by Tumor Subclass")
        ggsave(paste0(plot_dir, "taylor/re_pathway_subclass_emap.png"),
               plot = p_emap, width = 14, height = 9)
      }, error = function(e) cat("Warning: Could not create emapplot:", e$message, "\n"))
    } else {
      cat("No significant pathways found for tumor subclass comparison\n")
    }
  } else {
    cat("Not enough genes per subclass for pathway analysis (minimum 5 required)\n")
  }
}, error = function(e) {
  cat("Warning: Could not perform RE analysis by tumor subclass:", e$message, "\n")
})

#### ANALYSIS 2b: RE COMPARISON: TAYLOR vs KICS vs HOSTSEQ ####
cat("\n=== Analysis 2b: RE Comparison - Taylor vs KICS vs HostSeq ===\n")

tryCatch({
  # Load RE data for KICS and HostSeq if not already loaded
  if (!exists("te_kics_hostseq_re_split")) {
    te_kics_hostseq_re <- load_and_join_re_data(te_kics_hostseq_expand, re_germline_path)
    te_kics_hostseq_re_split <- split_re_genes(te_kics_hostseq_re)
  }

  # Add cohort labels
  te_taylor_re_split$cohort <- "Taylor"
  te_kics_hostseq_re_split$cohort <- te_kics_hostseq_re_split$cohort  # Already has cohort column

  # Combine all three
  te_all_cohorts_re <- rbind(
    te_taylor_re_split %>% select(gene_reg, cohort),
    te_kics_hostseq_re_split %>% select(gene_reg, cohort)
  )

  # Group genes by cohort
  geneClusters_cohort <- lapply(
    split(te_all_cohorts_re$gene_reg, te_all_cohorts_re$cohort),
    unique
  )

  cat("\nGenes in REs per cohort:\n")
  for (cohort_name in names(geneClusters_cohort)) {
    cat("  ", cohort_name, ":", length(geneClusters_cohort[[cohort_name]]), "genes\n")
  }

  # Run compareCluster for cohort comparison
  if (all(sapply(geneClusters_cohort, length) >= 5)) {
    ora_re_cohort_comparison <- tryCatch({
      compareCluster(
        geneCluster = geneClusters_cohort,
        fun = "enrichGO",
        OrgDb = org.Hs.eg.db,
        keyType = "SYMBOL",
        ont = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.1
      )
    }, error = function(e) {
      cat("Error in cohort comparison:", e$message, "\n")
      NULL
    })

    if (!is.null(ora_re_cohort_comparison) && nrow(as.data.frame(ora_re_cohort_comparison)) > 0) {
      cat("✓ Found", nrow(as.data.frame(ora_re_cohort_comparison)),
          "enriched pathways across cohorts\n")

      # Save results
      write.csv(as.data.frame(ora_re_cohort_comparison),
               paste0(r_dir_files, "re_pathway_cohort_comparison.csv"),
               row.names = FALSE)

      # Identify unique pathways per cohort
      ora_df <- as.data.frame(ora_re_cohort_comparison)
      taylor_unique <- ora_df %>%
        filter(Cluster == "Taylor" & !Description %in% ora_df$Description[ora_df$Cluster != "Taylor"])
      kics_unique <- ora_df %>%
        filter(Cluster == "KICS" & !Description %in% ora_df$Description[ora_df$Cluster != "KICS"])
      hostseq_unique <- ora_df %>%
        filter(Cluster == "HostSeq" & !Description %in% ora_df$Description[ora_df$Cluster != "HostSeq"])

      cat("\nUnique pathways per cohort:\n")
      cat("  Taylor:", nrow(taylor_unique), "\n")
      cat("  KICS:", nrow(kics_unique), "\n")
      cat("  HostSeq:", nrow(hostseq_unique), "\n")

      # Save unique pathways
      write.csv(taylor_unique,
               paste0(r_dir_files, "re_pathway_taylor_unique.csv"),
               row.names = FALSE)

      # Create comparison plots
      tryCatch({
        p_dot_cohort <- dotplot(ora_re_cohort_comparison, showCategory = 10)
        titled_print(p_dot_cohort, "RE Pathway Dotplot - Taylor vs KICS vs HostSeq")
        ggsave(paste0(plot_dir, "taylor/re_pathway_cohort_comparison_dot.png"),
               plot = p_dot_cohort, width = 14, height = 9)
      }, error = function(e) cat("Warning: Could not create cohort dotplot:", e$message, "\n"))

      tryCatch({
        p_cnet_cohort <- cnetplot(ora_re_cohort_comparison, showCategory = 5,
                                  colorEdge = TRUE, node_label = "category")
        titled_print(p_cnet_cohort, "RE Pathway Cnetplot - Taylor vs KICS vs HostSeq")
        ggsave(paste0(plot_dir, "taylor/re_pathway_cohort_comparison_cnet.png"),
               plot = p_cnet_cohort, width = 14, height = 9)
      }, error = function(e) cat("Warning: Could not create cohort cnetplot:", e$message, "\n"))

      tryCatch({
        ora_re_cohort_pairwise <- pairwise_termsim(ora_re_cohort_comparison)
        p_emap_cohort <- emapplot(ora_re_cohort_pairwise, showCategory = 20)
        titled_print(p_emap_cohort, "RE Pathway Emapplot - Taylor vs KICS vs HostSeq")
        ggsave(paste0(plot_dir, "taylor/re_pathway_cohort_comparison_emap.png"),
               plot = p_emap_cohort, width = 14, height = 9)
      }, error = function(e) cat("Warning: Could not create cohort emapplot:", e$message, "\n"))
    } else {
      cat("No significant pathways found for cohort comparison\n")
    }
  } else {
    cat("Not enough genes per cohort for pathway analysis (minimum 5 required)\n")
  }
}, error = function(e) {
  cat("Warning: Could not perform RE cohort comparison:", e$message, "\n")
})

# Close the PDF device at the end
dev.off()

cat("\n===== SCRIPT COMPLETED SUCCESSFULLY =====\n")
cat("All analysis sections have been processed.\n")
cat("Generated plots saved to:", plot_dir, "\n")
cat("PDF compilation saved to: graph_output.pdf\n")
cat("Text output saved to:", stdout_file, "\n")

# Close sink to stop redirecting output
sink()

cat("✓ Script completed successfully\n")
