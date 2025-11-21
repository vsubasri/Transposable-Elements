#!/usr/bin/env Rscript

# Tumour TE Visualization - Regulatory Elements
# RE pathway analysis and RE-RNA differential expression

# Source common setup and load data
source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/00_viz_common_setup.R")
source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/00_viz_load_data_tumour.R")

cat("Running 02_te_viz_tumour_07_re.R...\n")

#### PARAMETER SWEEP CONFIGURATION ####

# Parameter grid for RE-RNA differential expression analysis
cat("\n*** RE-RNA DIFFERENTIAL EXPRESSION ANALYSIS ***\n")
param_grid_re_rna <- expand.grid(
  min_samples_per_group = c(3),
  p_gene = c(0.05, 0.1),
  p_pathway = c(0.05, 0.1),
  q_pathway = c(0.05, 0.1),
  stringsAsFactors = FALSE
)
cat("Testing", nrow(param_grid_re_rna), "parameter combinations\n\n")

#### REGULATORY ELEMENTS ANALYSIS ####
write_output(quote(NULL), "Regulatory Elements (RE) Pathway Analysis - Tumour")

# Parameter grid for RE pathway analyses
# Filter genes by minimum samples, then run pathway analysis
cat("\n*** RE PATHWAY ANALYSIS ***\n")
param_grid_re_pathway <- expand.grid(
  min_samples_per_group = c(3),
  p_pathway = c(0.05, 0.1),
  q_pathway = c(0.05, 0.1),
  stringsAsFactors = FALSE
)
cat("Testing", nrow(param_grid_re_pathway), "parameter combinations\n\n")

# Load and join RE data for each relevant dataframe (do this once, outside the loop)
re_tumour_path <- "/Users/briannelaverty/Documents/R_Malkin/TE/data/final/tumour_annotSV_output.SV_RE_intersect.report"

cat("\n=== Loading and joining RE data for all cohorts ===\n")
te_aff_re_t <- load_and_join_re_data(te_aff_expand_t, re_tumour_path)
te_aff_re_split_t <- split_re_genes(te_aff_re_t)

# Loop through RE pathway parameter combinations
for (i in 1:nrow(param_grid_re_pathway)) {
  params_re <- param_grid_re_pathway[i, ]

  # Create parameter suffix
  param_suffix_re_pathway <- paste0(
    "_min", params_re$min_samples_per_group,
    "_ppathway", params_re$p_pathway,
    "_qpathway", params_re$q_pathway
  )

  cat("\n\n===== TESTING RE PATHWAY PARAMETERS", i, "/", nrow(param_grid_re_pathway), "=====\n")
  cat("min_samples_per_group =", params_re$min_samples_per_group,
      ", p_pathway =", params_re$p_pathway, ", q_pathway =", params_re$q_pathway, "\n")

# Analysis 1: General affected cohort (te_aff_expand_t)
cat("\n=== Regulatory Elements Analysis: Affected Cohort (Tumour) ===\n")

cat("\nRunning pathway analysis...\n")
n_genes_input_t <- length(unique(te_aff_re_split_t$gene_reg))
n_samples_input_t <- length(unique(te_aff_re_split_t$sample.x))

# Write data summary before analysis
summary_text <- paste0("RE Pathway Analysis Summary - Affected Tumour Cohort\n",
                      "====================================================\n\n",
                      "Total samples: ", n_samples_input_t, "\n",
                      "Total genes: ", n_genes_input_t, "\n",
                      "Parameters: min_samples=", params_re$min_samples_per_group,
                      ", p<", params_re$p_pathway, ", q<", params_re$q_pathway, "\n\n",
                      "Status: Running analysis...")
writeLines(summary_text, paste0(plot_dir, "reg_element/re_pathway_aff_tumour_summary", param_suffix_re_pathway, ".txt"))

ora_re_aff_t <- perform_re_pathway_analysis(te_aff_re_split_t, analysis_type = "general", min_samples = params_re$min_samples_per_group, p_pathway = params_re$p_pathway, q_pathway = params_re$q_pathway)

if (!is.null(ora_re_aff_t) && (inherits(ora_re_aff_t, "enrichResult") || inherits(ora_re_aff_t, "compareClusterResult")) && nrow(as.data.frame(ora_re_aff_t)) > 0) {
  cat("Pathway analysis successful! Found", nrow(as.data.frame(ora_re_aff_t)), "significant pathways\n")
  write.csv(as.data.frame(ora_re_aff_t), paste0(r_dir_files, "re_pathway_aff_tumour", param_suffix_re_pathway, ".csv"), row.names=FALSE)

  # Update summary with results
  summary_text <- paste0("RE Pathway Analysis Summary - Affected Tumour Cohort\n",
                        "====================================================\n\n",
                        "Total samples: ", n_samples_input_t, "\n",
                        "Total genes: ", n_genes_input_t, "\n",
                        "Parameters: min_samples=", params_re$min_samples_per_group,
                        ", p<", params_re$p_pathway, ", q<", params_re$q_pathway, "\n\n",
                        "Status: SUCCESS\n",
                        "Significant pathways found: ", nrow(as.data.frame(ora_re_aff_t)))
  writeLines(summary_text, paste0(plot_dir, "reg_element/re_pathway_aff_tumour_summary", param_suffix_re_pathway, ".txt"))

  cat("Creating barplot...\n")
  p_bar_ora_re_aff_t <- barplot(ora_re_aff_t, showCategory=20)
  titled_print(p_bar_ora_re_aff_t, "ORA barplot (RE - Affected Tumour)")
  ggsave(paste0(plot_dir, "reg_element/re_pathway_aff_tumour_bar", param_suffix_re_pathway, ".png"), plot=p_bar_ora_re_aff_t, width=14, height=9)

  cat("Creating dotplot...\n")
  p_dot_ora_re_aff_t <- dotplot(ora_re_aff_t, showCategory=20)
  titled_print(p_dot_ora_re_aff_t, "ORA dotplot (RE - Affected Tumour)")
  ggsave(paste0(plot_dir, "reg_element/re_pathway_aff_tumour_dot", param_suffix_re_pathway, ".png"), plot=p_dot_ora_re_aff_t, width=14, height=9)

  cat("Creating cnetplot...\n")
  p_cnet_ora_re_aff_t <- cnetplot(ora_re_aff_t, showCategory=10, colorEdge=TRUE, node_label="category")
  titled_print(p_cnet_ora_re_aff_t, "RE ORA cnetplot (Affected Tumour)")
  ggsave(paste0(plot_dir, "reg_element/re_pathway_aff_tumour_cnet", param_suffix_re_pathway, ".png"), plot=p_cnet_ora_re_aff_t, width=14, height=9)

  cat("Creating emapplot...\n")
  tryCatch({
    ora_re_aff_pairwise_t <- pairwise_termsim(ora_re_aff_t)
    p_emap_ora_re_aff_t <- emapplot(ora_re_aff_pairwise_t, showCategory=20)
    titled_print(p_emap_ora_re_aff_t, "ORA emapplot (RE - Affected Tumour)")
    ggsave(paste0(plot_dir, "reg_element/re_pathway_aff_tumour_emap", param_suffix_re_pathway, ".png"), plot=p_emap_ora_re_aff_t, width=14, height=9)
  }, error = function(e) {
    cat("Warning: Could not create emapplot:", e$message, "\n")
  })

  cat("\nTop 5 enriched pathways:\n")
  print(head(as.data.frame(ora_re_aff_t)[, c("Description", "pvalue", "p.adjust", "Count")], 5))
} else {
  # Update summary with failure status
  if (is.null(ora_re_aff_t)) {
    status_msg <- paste0("Status: INSUFFICIENT DATA\n",
                        "Reason: Not enough samples or genes to run pathway analysis\n",
                        "Minimum required: ", params_re$min_samples_per_group, " samples per group")
  } else {
    status_msg <- paste0("Status: NO SIGNIFICANT PATHWAYS\n",
                        "Analysis completed but no pathways met significance thresholds")
  }

  summary_text <- paste0("RE Pathway Analysis Summary - Affected Tumour Cohort\n",
                        "====================================================\n\n",
                        "Total samples: ", n_samples_input_t, "\n",
                        "Total genes: ", n_genes_input_t, "\n",
                        "Parameters: min_samples=", params_re$min_samples_per_group,
                        ", p<", params_re$p_pathway, ", q<", params_re$q_pathway, "\n\n",
                        status_msg)
  writeLines(summary_text, paste0(plot_dir, "reg_element/re_pathway_aff_tumour_summary", param_suffix_re_pathway, ".txt"))

  cat(status_msg, "\n")
}

# Analysis 2: Affected cohort colored by TP53 status (te_aff_expand_t - same data as Analysis 1)
cat("\n\n=== Regulatory Elements Analysis: Affected Cohort (by TP53 Status - Tumour) ===\n")

# Check if TP53_status column exists (using te_aff_re_split_t from Analysis 1)
if ("TP53_status" %in% colnames(te_aff_re_split_t)) {
  cat("TP53_status groups:\n")
  print(table(te_aff_re_split_t$TP53_status))

  cat("\nRunning pathway analysis by TP53 status...\n")
  n_genes_input_tp53_t <- length(unique(te_aff_re_split_t$gene_reg))
  n_samples_input_tp53_t <- length(unique(te_aff_re_split_t$sample.x))
  # Count samples per TP53 group
  samples_per_tp53_t <- table(unique(te_aff_re_split_t[, c("sample.x", "TP53_status")])$TP53_status)
  samples_per_tp53_t_str <- paste(names(samples_per_tp53_t), "=", samples_per_tp53_t, "samples", collapse=", ")

  # Count genes per TP53 group BEFORE running analysis
  genes_per_tp53_t <- te_aff_re_split_t %>%
    group_by(TP53_status) %>%
    summarise(n_genes = n_distinct(gene_reg), .groups = "drop")
  genes_per_tp53_t_str <- paste(genes_per_tp53_t$TP53_status, "=", genes_per_tp53_t$n_genes, "genes", collapse=", ")

  # Write data summary before analysis
  summary_text <- paste0("RE Pathway Analysis Summary - Affected Tumour TP53 Status\n",
                        "=========================================================\n\n",
                        "Total samples: ", n_samples_input_tp53_t, "\n",
                        "Samples per group: ", samples_per_tp53_t_str, "\n",
                        "Total genes: ", n_genes_input_tp53_t, "\n",
                        "Genes per group: ", genes_per_tp53_t_str, "\n",
                        "Parameters: min_samples=", params_re$min_samples_per_group,
                        ", p<", params_re$p_pathway, ", q<", params_re$q_pathway, "\n\n",
                        "Status: Running analysis...")
  writeLines(summary_text, paste0(plot_dir, "reg_element/re_pathway_aff_tp53_tumour_summary", param_suffix_re_pathway, ".txt"))

  ora_re_aff_tp53_t <- perform_re_pathway_analysis(te_aff_re_split_t, analysis_type = "tp53", min_samples = params_re$min_samples_per_group, p_pathway = params_re$p_pathway, q_pathway = params_re$q_pathway)

  if (!is.null(ora_re_aff_tp53_t) && (inherits(ora_re_aff_tp53_t, "enrichResult") || inherits(ora_re_aff_tp53_t, "compareClusterResult")) && nrow(as.data.frame(ora_re_aff_tp53_t)) > 0) {
    cat("Pathway analysis successful! Found", nrow(as.data.frame(ora_re_aff_tp53_t)), "significant pathways\n")
    write.csv(as.data.frame(ora_re_aff_tp53_t), paste0(r_dir_files, "re_pathway_aff_tp53_tumour", param_suffix_re_pathway, ".csv"), row.names=FALSE)

    # Update summary with results
    summary_text <- paste0("RE Pathway Analysis Summary - Affected Tumour TP53 Status\n",
                          "=========================================================\n\n",
                          "Total samples: ", n_samples_input_tp53_t, "\n",
                          "Samples per group: ", samples_per_tp53_t_str, "\n",
                          "Total genes: ", n_genes_input_tp53_t, "\n",
                          "Genes per group: ", genes_per_tp53_t_str, "\n",
                          "Parameters: min_samples=", params_re$min_samples_per_group,
                          ", p<", params_re$p_pathway, ", q<", params_re$q_pathway, "\n\n",
                          "Status: SUCCESS\n",
                          "Significant pathways found: ", nrow(as.data.frame(ora_re_aff_tp53_t)))
    writeLines(summary_text, paste0(plot_dir, "reg_element/re_pathway_aff_tp53_tumour_summary", param_suffix_re_pathway, ".txt"))

    cat("Creating dotplot...\n")
    p_dot_ora_re_tp53_t <- dotplot(ora_re_aff_tp53_t, showCategory=20)
    titled_print(p_dot_ora_re_tp53_t, "ORA dotplot (RE - TP53 Status Tumour)")
    ggsave(paste0(plot_dir, "reg_element/re_pathway_tp53_tumour_dot", param_suffix_re_pathway, ".png"), plot=p_dot_ora_re_tp53_t, width=14, height=9)

    cat("Creating cnetplot...\n")
    p_cnet_ora_re_tp53_t <- cnetplot(ora_re_aff_tp53_t, showCategory=10, colorEdge=TRUE, node_label="category")
    titled_print(p_cnet_ora_re_tp53_t, "RE ORA cnetplot (TP53 Status Tumour)")
    ggsave(paste0(plot_dir, "reg_element/re_pathway_tp53_tumour_cnet", param_suffix_re_pathway, ".png"), plot=p_cnet_ora_re_tp53_t, width=14, height=9)

    cat("Creating emapplot...\n")
    # Check if multiple clusters have results
    cluster_counts <- table(as.data.frame(ora_re_aff_tp53_t)$Cluster)
    cat("Cluster distribution:", paste(names(cluster_counts), "=", cluster_counts, collapse=", "), "\n")

    if (length(cluster_counts) > 1) {
      tryCatch({
        ora_re_tp53_pairwise_t <- pairwise_termsim(ora_re_aff_tp53_t)
        p_emap_ora_re_tp53_t <- emapplot(ora_re_tp53_pairwise_t, showCategory=20,
                                          pie.params = list(pie = "count"),
                                          cluster.params = list(cluster = TRUE, legend = TRUE))
        titled_print(p_emap_ora_re_tp53_t, "ORA emapplot (RE - TP53 Tumour)")
        ggsave(paste0(plot_dir, "reg_element/re_pathway_tp53_tumour_emap", param_suffix_re_pathway, ".png"), plot=p_emap_ora_re_tp53_t, width=14, height=9)
      }, error = function(e) {
        cat("Warning: Could not create emapplot:", e$message, "\n")
      })
    } else {
      cat("Skipping emapplot: Only one cluster has results (", names(cluster_counts), ")\n")
    }
  } else {
    # Update summary with failure status
    if (is.null(ora_re_aff_tp53_t)) {
      status_msg <- paste0("Status: INSUFFICIENT DATA\n",
                          "Reason: Not enough samples per group to run pathway analysis\n",
                          "Minimum required: ", params_re$min_samples_per_group, " samples per group")
    } else {
      status_msg <- paste0("Status: NO SIGNIFICANT PATHWAYS\n",
                          "Analysis completed but no pathways met significance thresholds")
    }

    summary_text <- paste0("RE Pathway Analysis Summary - Affected Tumour TP53 Status\n",
                          "=========================================================\n\n",
                          "Total samples: ", n_samples_input_tp53_t, "\n",
                          "Samples per group: ", samples_per_tp53_t_str, "\n",
                          "Total genes: ", n_genes_input_tp53_t, "\n",
                          "Genes per group: ", genes_per_tp53_t_str, "\n",
                          "Parameters: min_samples=", params_re$min_samples_per_group,
                          ", p<", params_re$p_pathway, ", q<", params_re$q_pathway, "\n\n",
                          status_msg)
    writeLines(summary_text, paste0(plot_dir, "reg_element/re_pathway_aff_tp53_tumour_summary", param_suffix_re_pathway, ".txt"))

    cat(status_msg, "\n")
  }
} else {
  cat("WARNING: TP53_status column not found in te_aff_re_split_t, skipping TP53 analysis\n")
}

} # End of RE pathway parameter loop


#### RE-RNA DIFFERENTIAL EXPRESSION TEST ####
write_output(quote(NULL), "Running RE-RNA Differential Expression Analysis")

# Check if RNA data is available (requires running RNA script first)
if (!exists("rna_filtered")) {
  # Try to load saved rna_filtered from RNA script
  rna_filtered_file <- paste0(r_dir_files, "rna_filtered_tumour.rds")
  if (file.exists(rna_filtered_file)) {
    cat("Loading rna_filtered from previously run RNA script...\n")
    rna_filtered <- readRDS(rna_filtered_file)
    cat("✓ Successfully loaded rna_filtered\n\n")
  } else {
    cat("⚠ WARNING: RNA data not processed (rna_filtered not found).\n")
    cat("  The RE-RNA analysis requires RNA processing from 02_te_viz_tumour_07_rna.R\n")
    cat("  Skipping RE-RNA differential expression analysis.\n")
    cat("  To run this section, first run the RNA script or run the ALL script.\n\n")
  }
}

if (exists("rna_filtered")) {
  # Set up summary tracking for parameter sweep
  sweep_summary <- data.frame()

  re_report_path_tumour <- "/Users/briannelaverty/Documents/R_Malkin/TE/data/final/tumour_annotSV_output.SV_RE_intersect.report"

  # Loop through parameter combinations
  for (i in 1:nrow(param_grid_re_rna)) {
    params <- param_grid_re_rna[i, ]

    # Create parameter suffix for filenames
    param_suffix <- paste0(
      "_min", params$min_samples_per_group,
      "_pgene", params$p_gene,
      "_ppathway", params$p_pathway,
      "_qpathway", params$q_pathway
    )

    cat("\n\n===== TESTING RE-RNA PARAMETERS", i, "/", nrow(param_grid_re_rna), "=====\n")
    cat("min_samples_per_group =", params$min_samples_per_group,
        ", p_gene =", params$p_gene,
        ", p_pathway =", params$p_pathway,
        ", q_pathway =", params$q_pathway, "\n")

    # Run differential expression test
    re_rna_results_t <- test_rna_by_gene_re_status(
      re_report_path = re_report_path_tumour,
      rna_data = rna_filtered,
      sample_type = "tumour",
      min_gene_mentions = 1,
      min_samples_per_group = params$min_samples_per_group
    )

    # Save differential expression results
    if (nrow(re_rna_results_t) > 0) {
      output_file <- paste0(r_dir_files, "re_rna_differential_tumour", param_suffix, ".csv")
      write.csv(re_rna_results_t, output_file, row.names = FALSE)

      # Filter genes by p_gene threshold
      pathway_genes_t <- re_rna_results_t %>%
        filter(p_value < params$p_gene) %>%
        pull(gene) %>%
        unique()

      cat("Genes with p-value <", params$p_gene, ":", length(pathway_genes_t), "\n")

      # Pathway analysis if enough genes
      if (length(pathway_genes_t) >= 5) {
        pathway_df_t <- data.frame(Gene_name = pathway_genes_t)

        ora_re_rna_t <- tryCatch({
          perform_ora_custom_cutoffs(pathway_df_t,
                                    p_pathway = params$p_pathway,
                                    q_pathway = params$q_pathway,
                                    nsample_thresh = 0,
                                    filter_exon = FALSE)
        }, error = function(e) {
          cat("Error in pathway analysis:", e$message, "\n")
          NULL
        })

        if (!is.null(ora_re_rna_t) && nrow(as.data.frame(ora_re_rna_t)) > 0) {
          n_pathways <- nrow(as.data.frame(ora_re_rna_t))
          cat("✓ Found", n_pathways, "enriched pathways\n")

          # Save pathway results
          pathway_file <- paste0(r_dir_files, "re_rna_pathway_tumour", param_suffix, ".csv")
          write.csv(as.data.frame(ora_re_rna_t), pathway_file, row.names = FALSE)

          # Save gene list
          gene_list_file <- paste0(r_dir_files, "re_rna_genes_tumour", param_suffix, ".txt")
          writeLines(pathway_genes_t, gene_list_file)

          # Create plots
            tryCatch({
              p_bar_re_rna_t <- barplot(ora_re_rna_t, showCategory = 20)
              titled_print(p_bar_re_rna_t, "RE-RNA Pathway Barplot (Tumour)")
              ggsave(paste0(plot_dir, "reg_element/re_rna_pathway_bar_tumour", param_suffix, ".png"),
                     plot = p_bar_re_rna_t, width = 14, height = 9)
            }, error = function(e) cat("Warning: Could not create barplot:", e$message, "\n"))

            tryCatch({
              p_dot_re_rna_t <- dotplot(ora_re_rna_t, showCategory = 20)
              titled_print(p_dot_re_rna_t, "RE-RNA Pathway Dotplot (Tumour)")
              ggsave(paste0(plot_dir, "reg_element/re_rna_pathway_dot_tumour", param_suffix, ".png"),
                     plot = p_dot_re_rna_t, width = 14, height = 9)
            }, error = function(e) cat("Warning: Could not create dotplot:", e$message, "\n"))

            tryCatch({
              p_cnet_re_rna_t <- cnetplot(ora_re_rna_t, showCategory = 10, colorEdge = TRUE, node_label = "category")
              titled_print(p_cnet_re_rna_t, "RE-RNA Pathway Cnetplot (Tumour)")
              ggsave(paste0(plot_dir, "reg_element/re_rna_pathway_cnet_tumour", param_suffix, ".png"),
                     plot = p_cnet_re_rna_t, width = 14, height = 9)
            }, error = function(e) cat("Warning: Could not create cnetplot:", e$message, "\n"))

            tryCatch({
              ora_re_rna_pairwise_t <- pairwise_termsim(ora_re_rna_t)
              p_emap_re_rna_t <- emapplot(ora_re_rna_pairwise_t, showCategory = 20)
              titled_print(p_emap_re_rna_t, "RE-RNA Pathway Emapplot (Tumour)")
              ggsave(paste0(plot_dir, "reg_element/re_rna_pathway_emap_tumour", param_suffix, ".png"),
                     plot = p_emap_re_rna_t, width = 14, height = 9)
            }, error = function(e) cat("Warning: Could not create emapplot:", e$message, "\n"))

          # Record sweep summary
          sweep_summary <- rbind(sweep_summary, data.frame(
            min_samples_per_group = params$min_samples_per_group,
            p_gene = params$p_gene,
            p_pathway = params$p_pathway,
            q_pathway = params$q_pathway,
            n_de_genes = nrow(re_rna_results_t),
            n_filtered_genes = length(pathway_genes_t),
            n_pathways = n_pathways,
            stringsAsFactors = FALSE
          ))
        } else {
          cat("No significant pathways found\n")
          sweep_summary <- rbind(sweep_summary, data.frame(
            min_samples_per_group = params$min_samples_per_group,
            p_gene = params$p_gene,
            p_pathway = params$p_pathway,
            q_pathway = params$q_pathway,
            n_de_genes = nrow(re_rna_results_t),
            n_filtered_genes = length(pathway_genes_t),
            n_pathways = 0,
            stringsAsFactors = FALSE
          ))
        }
      } else {
        cat("Not enough genes (< 5) for pathway analysis\n")
        sweep_summary <- rbind(sweep_summary, data.frame(
          min_samples_per_group = params$min_samples_per_group,
          p_gene = params$p_gene,
          p_pathway = params$p_pathway,
          q_pathway = params$q_pathway,
          n_de_genes = nrow(re_rna_results_t),
          n_filtered_genes = length(pathway_genes_t),
          n_pathways = NA,
          stringsAsFactors = FALSE
        ))
      }
    } else {
      cat("No differential expression results for this parameter combination\n")
      sweep_summary <- rbind(sweep_summary, data.frame(
        min_samples_per_group = params$min_samples_per_group,
        p_gene = params$p_gene,
        p_pathway = params$p_pathway,
        q_pathway = params$q_pathway,
        n_de_genes = 0,
        n_filtered_genes = 0,
        n_pathways = NA,
        stringsAsFactors = FALSE
      ))
    }
  }

  # Save sweep summary
  if (nrow(sweep_summary) > 0) {
    summary_file <- paste0(r_dir_files, "re_rna_sweep_tumour_summary.csv")
    write.csv(sweep_summary, summary_file, row.names = FALSE)
    cat("\n✓ Parameter sweep complete! Summary saved to:", summary_file, "\n")
    cat("\nTop 5 parameter combinations by number of pathways:\n")
    print(head(sweep_summary[order(-sweep_summary$n_pathways), ], 5))
  }
}



cat("✓ Script completed successfully\n")
