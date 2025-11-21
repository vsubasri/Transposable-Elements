#!/usr/bin/env Rscript

# Germline TE Visualization - Regulatory Elements
# RE pathway analysis and RE-RNA differential expression

# Source common setup and load data
source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/00_viz_common_setup.R")
source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/00_viz_load_data_germline.R")

cat("Running 02_te_viz_germline_06_re.R...\n")

#### PARAMETER SWEEP CONFIGURATION ####

# Parameter grid for RE-RNA differential expression analysis
cat("\n*** RE-RNA DIFFERENTIAL EXPRESSION ANALYSIS ***\n")
param_grid_re_rna <- expand.grid(
  min_samples_per_group = c(3, 5),
  p_gene = c(0.05, 0.1),
  p_pathway = c(0.05, 0.1),
  q_pathway = c(0.05, 0.1),
  stringsAsFactors = FALSE
)
cat("Testing", nrow(param_grid_re_rna), "parameter combinations\n\n")

#### REGULATORY ELEMENTS ANALYSIS ####
write_output(quote(NULL), "Regulatory Elements (RE) Pathway Analysis")

# Parameter grid for RE pathway analyses
# Filter genes by minimum samples, then run pathway analysis
cat("\n*** RE PATHWAY ANALYSIS ***\n")
param_grid_re_pathway <- expand.grid(
  min_samples_per_group = c(5),
  p_pathway = c(0.05),
  q_pathway = c(0.05),
  stringsAsFactors = FALSE
)
cat("Testing", nrow(param_grid_re_pathway), "parameter combinations\n\n")

# Load and join RE data for each relevant dataframe (do this once, outside the loop)
re_germline_path <- "/Users/briannelaverty/Documents/R_Malkin/TE/data/final/germline_annotSV_output.SV_RE_intersect.report"

cat("\n=== Loading and joining RE data for all cohorts ===\n")
te_aff_re <- load_and_join_re_data(te_aff_expand, re_germline_path)
te_aff_re_split <- split_re_genes(te_aff_re)

te_lfs_re <- load_and_join_re_data(te_lfs_expand, re_germline_path)
te_lfs_re_split <- split_re_genes(te_lfs_re)

te_kics_hostseq_re <- load_and_join_re_data(te_kics_hostseq_expand, re_germline_path)
te_kics_hostseq_re_split <- split_re_genes(te_kics_hostseq_re)

kics_sample_type <- prep_kics_sample_type("/Users/briannelaverty/Documents/R_Malkin/clinical/kics_germline_sample_type.csv")
te_kics_expand_re <- load_and_join_re_data(te_kics_expand, re_germline_path)
te_kics_expand_re_split <- split_re_genes(te_kics_expand_re)
te_kics_re_sampletype <- merge(te_kics_expand_re_split, kics_sample_type, by.x = "sample.x", by.y = "sample", all.x = TRUE)
te_kics_re_sampletype <- te_kics_re_sampletype %>%
  filter(sample_type %in% c("Blood", "Fibroblasts", "Tissue (fresh)"))

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

# Analysis 1: General affected cohort (te_aff_expand)
cat("\n=== Regulatory Elements Analysis: Affected Cohort ===\n")

cat("\nRunning pathway analysis...\n")
n_genes_input <- length(unique(te_aff_re_split$gene_reg))
n_samples_input <- length(unique(te_aff_re_split$sample.x))

# Write data summary before analysis
summary_text <- paste0("RE Pathway Analysis Summary - Affected Cohort\n",
                      "================================================\n\n",
                      "Total samples: ", n_samples_input, "\n",
                      "Total genes: ", n_genes_input, "\n",
                      "Parameters: min_samples=", params_re$min_samples_per_group,
                      ", p<", params_re$p_pathway, ", q<", params_re$q_pathway, "\n\n",
                      "Status: Running analysis...")
writeLines(summary_text, paste0(plot_dir, "reg_element/re_pathway_aff_general_summary", param_suffix_re_pathway, ".txt"))

ora_re_aff <- perform_re_pathway_analysis(te_aff_re_split, analysis_type = "general", min_samples = params_re$min_samples_per_group, p_pathway = params_re$p_pathway, q_pathway = params_re$q_pathway)

if (!is.null(ora_re_aff) && (inherits(ora_re_aff, "enrichResult") || inherits(ora_re_aff, "compareClusterResult")) && nrow(as.data.frame(ora_re_aff)) > 0) {
  cat("Pathway analysis successful! Found", nrow(as.data.frame(ora_re_aff)), "significant pathways\n")
  write.csv(as.data.frame(ora_re_aff), paste0(r_dir_files, "re_pathway_aff_general", param_suffix_re_pathway, ".csv"), row.names=FALSE)

  # Update summary with results
  summary_text <- paste0("RE Pathway Analysis Summary - Affected Cohort\n",
                        "================================================\n\n",
                        "Total samples: ", n_samples_input, "\n",
                        "Total genes: ", n_genes_input, "\n",
                        "Parameters: min_samples=", params_re$min_samples_per_group,
                        ", p<", params_re$p_pathway, ", q<", params_re$q_pathway, "\n\n",
                        "Status: SUCCESS\n",
                        "Significant pathways found: ", nrow(as.data.frame(ora_re_aff)))
  writeLines(summary_text, paste0(plot_dir, "reg_element/re_pathway_aff_general_summary", param_suffix_re_pathway, ".txt"))

  cat("Creating barplot...\n")
  p_bar_ora_re_aff <- barplot(ora_re_aff, showCategory=20)
  titled_print(p_bar_ora_re_aff, "ORA barplot (RE - Affected)")
  ggsave(paste0(plot_dir, "reg_element/re_pathway_aff_general_bar", param_suffix_re_pathway, ".png"), plot=p_bar_ora_re_aff, width=14, height=9)

  cat("Creating dotplot...\n")
  p_dot_ora_re_aff <- dotplot(ora_re_aff, showCategory=20)
  titled_print(p_dot_ora_re_aff, "ORA dotplot (RE - Affected)")
  ggsave(paste0(plot_dir, "reg_element/re_pathway_aff_general_dot", param_suffix_re_pathway, ".png"), plot=p_dot_ora_re_aff, width=14, height=9)

  cat("Creating cnetplot...\n")
  p_cnet_ora_re_aff <- cnetplot(ora_re_aff, showCategory=10, colorEdge=TRUE, node_label="category")
  titled_print(p_cnet_ora_re_aff, "RE ORA cnetplot (Affected)")
  ggsave(paste0(plot_dir, "reg_element/re_pathway_aff_general_cnet", param_suffix_re_pathway, ".png"), plot=p_cnet_ora_re_aff, width=14, height=9)

  cat("Creating emapplot...\n")
  tryCatch({
    ora_re_aff_pairwise <- pairwise_termsim(ora_re_aff)
    p_emap_ora_re_aff <- emapplot(ora_re_aff_pairwise, showCategory=20)
    titled_print(p_emap_ora_re_aff, "ORA emapplot (RE - Affected)")
    ggsave(paste0(plot_dir, "reg_element/re_pathway_aff_general_emap", param_suffix_re_pathway, ".png"), plot=p_emap_ora_re_aff, width=14, height=9)
  }, error = function(e) {
    cat("Warning: Could not create emapplot:", e$message, "\n")
  })

  # Simplified pathway plots
  cat("Creating simplified pathways...\n")
  ora_re_aff_simple <- simplify(ora_re_aff, cutoff=0.5, by="p.adjust", select_fun=min)
  write.csv(as.data.frame(ora_re_aff_simple), paste0(r_dir_files, "re_pathway_aff_general_simple", param_suffix_re_pathway, ".csv"), row.names=FALSE)
  if (!is.null(ora_re_aff_simple) && nrow(as.data.frame(ora_re_aff_simple)) > 0) {
    p_bar_ora_re_aff_simple <- barplot(ora_re_aff_simple, showCategory=20)
    titled_print(p_bar_ora_re_aff_simple, "ORA barplot (RE - Affected, simplified)")
    ggsave(paste0(plot_dir, "reg_element/re_pathway_aff_general_bar_simplified", param_suffix_re_pathway, ".png"), plot=p_bar_ora_re_aff_simple, width=14, height=9)
    p_dot_ora_re_aff_simple <- dotplot(ora_re_aff_simple, showCategory=20)
    titled_print(p_dot_ora_re_aff_simple, "ORA dotplot (RE - Affected, simplified)")
    ggsave(paste0(plot_dir, "reg_element/re_pathway_aff_general_dot_simplified", param_suffix_re_pathway, ".png"), plot=p_dot_ora_re_aff_simple, width=14, height=9)
    p_cnet_ora_re_aff_simple <- cnetplot(ora_re_aff_simple, showCategory=10, colorEdge=TRUE, node_label="category")
    titled_print(p_cnet_ora_re_aff_simple, "RE ORA cnetplot (Affected, simplified)")
    ggsave(paste0(plot_dir, "reg_element/re_pathway_aff_general_cnet_simplified", param_suffix_re_pathway, ".png"), plot=p_cnet_ora_re_aff_simple, width=14, height=9)
    tryCatch({
      ora_re_aff_simple_pairwise <- pairwise_termsim(ora_re_aff_simple)
      p_emap_ora_re_aff_simple <- emapplot(ora_re_aff_simple_pairwise, showCategory=20)
      titled_print(p_emap_ora_re_aff_simple, "ORA emapplot (RE - Affected, simplified)")
      ggsave(paste0(plot_dir, "reg_element/re_pathway_aff_general_emap_simplified", param_suffix_re_pathway, ".png"), plot=p_emap_ora_re_aff_simple, width=14, height=9)
    }, error = function(e) {
      cat("Warning: Could not create simplified emapplot:", e$message, "\n")
    })
  }

  cat("\nTop 5 enriched pathways:\n")
  print(head(as.data.frame(ora_re_aff)[, c("Description", "pvalue", "p.adjust", "Count")], 5))
} else {
  # Update summary with failure reason
  if (is.null(ora_re_aff)) {
    status_msg <- paste0("Status: INSUFFICIENT DATA\n",
                        "Reason: Not enough samples or genes to run pathway analysis")
  } else {
    status_msg <- paste0("Status: NO SIGNIFICANT PATHWAYS\n",
                        "Analysis completed but no pathways met significance thresholds")
  }

  summary_text <- paste0("RE Pathway Analysis Summary - Affected Cohort\n",
                        "================================================\n\n",
                        "Total samples: ", n_samples_input, "\n",
                        "Total genes: ", n_genes_input, "\n",
                        "Parameters: min_samples=", params_re$min_samples_per_group,
                        ", p<", params_re$p_pathway, ", q<", params_re$q_pathway, "\n\n",
                        status_msg)
  writeLines(summary_text, paste0(plot_dir, "reg_element/re_pathway_aff_general_summary", param_suffix_re_pathway, ".txt"))

  cat(status_msg, "\n")
}

# Analysis 2: LFS cohort colored by Cancer status (te_lfs_expand)
cat("\n\n=== Regulatory Elements Analysis: LFS Cohort (by Cancer Status) ===\n")
te_lfs_re <- load_and_join_re_data(te_lfs_expand, re_germline_path)
te_lfs_re_split <- split_re_genes(te_lfs_re)

# Check if Cancer column exists for cancer status analysis
if ("Cancer" %in% colnames(te_lfs_re_split)) {
  cat("Cancer status groups:\n")
  print(table(te_lfs_re_split$Cancer))

  cat("\nRunning pathway analysis by Cancer status (Affected vs Unaffected)...\n")
  n_genes_input_lfs <- length(unique(te_lfs_re_split$gene_reg))
  n_samples_input_lfs <- length(unique(te_lfs_re_split$sample.x))
  # Count samples per Cancer group
  samples_per_cancer <- table(unique(te_lfs_re_split[, c("sample.x", "Cancer")])$Cancer)
  samples_per_cancer_str <- paste(names(samples_per_cancer), "=", samples_per_cancer, "samples", collapse=", ")

  # Write data summary before analysis
  genes_per_cancer <- te_lfs_re_split %>%
    group_by(Cancer) %>%
    summarise(n_genes = n_distinct(gene_reg), .groups = "drop")
  genes_per_cancer_str <- paste(genes_per_cancer$Cancer, "=", genes_per_cancer$n_genes, "genes", collapse=", ")

  summary_text <- paste0("RE Pathway Analysis Summary - LFS Cancer Status\n",
                        "================================================\n\n",
                        "Total samples: ", n_samples_input_lfs, "\n",
                        "Samples per group: ", samples_per_cancer_str, "\n",
                        "Total genes: ", n_genes_input_lfs, "\n",
                        "Genes per group: ", genes_per_cancer_str, "\n",
                        "Parameters: min_samples=", params_re$min_samples_per_group,
                        ", p<", params_re$p_pathway, ", q<", params_re$q_pathway, "\n\n",
                        "Status: Running analysis...")
  writeLines(summary_text, paste0(plot_dir, "reg_element/re_pathway_lfs_cancer_summary", param_suffix_re_pathway, ".txt"))

  ora_re_lfs_cancer <- perform_re_pathway_analysis(te_lfs_re_split, analysis_type = "cancer_status", min_samples = params_re$min_samples_per_group, p_pathway = params_re$p_pathway, q_pathway = params_re$q_pathway)

  if (!is.null(ora_re_lfs_cancer) && (inherits(ora_re_lfs_cancer, "enrichResult") || inherits(ora_re_lfs_cancer, "compareClusterResult")) && nrow(as.data.frame(ora_re_lfs_cancer)) > 0) {
    cat("Pathway analysis successful! Found", nrow(as.data.frame(ora_re_lfs_cancer)), "significant pathways\n")
    write.csv(as.data.frame(ora_re_lfs_cancer), paste0(r_dir_files, "re_pathway_lfs_cancer", param_suffix_re_pathway, ".csv"), row.names=FALSE)

    # Update summary with results
    summary_text <- paste0("RE Pathway Analysis Summary - LFS Cancer Status\n",
                          "================================================\n\n",
                          "Total samples: ", n_samples_input_lfs, "\n",
                          "Samples per group: ", samples_per_cancer_str, "\n",
                          "Total genes: ", n_genes_input_lfs, "\n",
                          "Genes per group: ", genes_per_cancer_str, "\n",
                          "Parameters: min_samples=", params_re$min_samples_per_group,
                          ", p<", params_re$p_pathway, ", q<", params_re$q_pathway, "\n\n",
                          "Status: SUCCESS\n",
                          "Significant pathways found: ", nrow(as.data.frame(ora_re_lfs_cancer)))
    writeLines(summary_text, paste0(plot_dir, "reg_element/re_pathway_lfs_cancer_summary", param_suffix_re_pathway, ".txt"))

    cat("Creating dotplot...\n")
    p_dot_ora_re_lfs <- dotplot(ora_re_lfs_cancer, showCategory=20)
    titled_print(p_dot_ora_re_lfs, "ORA dotplot (RE - LFS by Cancer)")
    ggsave(paste0(plot_dir, "reg_element/re_pathway_lfs_cancer_dot", param_suffix_re_pathway, ".png"), plot=p_dot_ora_re_lfs, width=14, height=9)

    cat("Creating cnetplot...\n")
    p_cnet_ora_re_lfs <- cnetplot(ora_re_lfs_cancer, showCategory=10, colorEdge=TRUE, node_label="category")
    titled_print(p_cnet_ora_re_lfs, "RE ORA cnetplot (LFS by Cancer)")
    ggsave(paste0(plot_dir, "reg_element/re_pathway_lfs_cancer_cnet", param_suffix_re_pathway, ".png"), plot=p_cnet_ora_re_lfs, width=14, height=9)

    cat("Creating emapplot...\n")
    tryCatch({
      # Check if multiple clusters have results
      cluster_counts <- table(as.data.frame(ora_re_lfs_cancer)$Cluster)
      if (length(cluster_counts) > 1) {
        ora_re_lfs_pairwise <- pairwise_termsim(ora_re_lfs_cancer)
        p_emap_ora_re_lfs <- emapplot(ora_re_lfs_pairwise, showCategory=20,
                                       pie.params = list(pie = "count"),
                                       cluster.params = list(cluster = TRUE, legend = TRUE))
        titled_print(p_emap_ora_re_lfs, "ORA emapplot (RE - LFS by Cancer)")
        ggsave(paste0(plot_dir, "reg_element/re_pathway_re_lfs_cancer", param_suffix_re_pathway, ".png"), plot=p_emap_ora_re_lfs, width=14, height=9)
      } else {
        cat("Skipping emapplot: Only one cluster has results (", names(cluster_counts), ")\n")
      }
    }, error = function(e) {
      cat("Warning: Could not create emapplot:", e$message, "\n")
    })

    # Simplified pathway plots
    cat("\nCreating simplified pathway plots...\n")
    ora_re_lfs_cancer_simple <- simplify(ora_re_lfs_cancer, cutoff=0.5, by="p.adjust", select_fun=min)
    write.csv(as.data.frame(ora_re_lfs_cancer_simple), paste0(r_dir_files, "re_pathway_lfs_cancer_simple", param_suffix_re_pathway, ".csv"), row.names=FALSE)

    if (!is.null(ora_re_lfs_cancer_simple) && nrow(as.data.frame(ora_re_lfs_cancer_simple)) > 0) {
      # Dot (no bar for compareCluster)
      p_dot_ora_re_lfs_simple <- dotplot(ora_re_lfs_cancer_simple, showCategory=20)
      titled_print(p_dot_ora_re_lfs_simple, "ORA dotplot simplified (RE - LFS by Cancer)")
      ggsave(paste0(plot_dir, "reg_element/re_pathway_lfs_cancer_dot_simplified", param_suffix_re_pathway, ".png"), plot=p_dot_ora_re_lfs_simple, width=14, height=9)

      # Cnet
      p_cnet_ora_re_lfs_simple <- cnetplot(ora_re_lfs_cancer_simple, showCategory=10, colorEdge=TRUE, node_label="category")
      titled_print(p_cnet_ora_re_lfs_simple, "RE ORA cnetplot simplified (LFS by Cancer)")
      ggsave(paste0(plot_dir, "reg_element/re_pathway_lfs_cancer_cnet_simplified", param_suffix_re_pathway, ".png"), plot=p_cnet_ora_re_lfs_simple, width=14, height=9)

      # Emap
      tryCatch({
        cluster_counts_simple <- table(as.data.frame(ora_re_lfs_cancer_simple)$Cluster)
        if (length(cluster_counts_simple) > 1) {
          ora_re_lfs_simple_pairwise <- pairwise_termsim(ora_re_lfs_cancer_simple)
          p_emap_ora_re_lfs_simple <- emapplot(ora_re_lfs_simple_pairwise, showCategory=20,
                                                pie.params = list(pie = "count"),
                                                cluster.params = list(cluster = TRUE, legend = TRUE))
          titled_print(p_emap_ora_re_lfs_simple, "ORA emapplot simplified (RE - LFS by Cancer)")
          ggsave(paste0(plot_dir, "reg_element/re_pathway_lfs_cancer_emap_simplified", param_suffix_re_pathway, ".png"), plot=p_emap_ora_re_lfs_simple, width=14, height=9)
        } else {
          cat("Skipping simplified emapplot: Only one cluster has results\n")
        }
      }, error = function(e) {
        cat("Warning: Could not create simplified emapplot:", e$message, "\n")
      })
    } else {
      cat("No significant pathways in simplified results\n")
    }
  } else {
    # Update summary with failure reason
    if (is.null(ora_re_lfs_cancer)) {
      status_msg <- paste0("Status: INSUFFICIENT DATA\n",
                          "Reason: Not enough samples or genes to run pathway analysis")
    } else {
      status_msg <- paste0("Status: NO SIGNIFICANT PATHWAYS\n",
                          "Analysis completed but no pathways met significance thresholds")
    }

    summary_text <- paste0("RE Pathway Analysis Summary - LFS Cancer Status\n",
                          "================================================\n\n",
                          "Total samples: ", n_samples_input_lfs, "\n",
                          "Samples per group: ", samples_per_cancer_str, "\n",
                          "Total genes: ", n_genes_input_lfs, "\n",
                          "Genes per group: ", genes_per_cancer_str, "\n",
                          "Parameters: min_samples=", params_re$min_samples_per_group,
                          ", p<", params_re$p_pathway, ", q<", params_re$q_pathway, "\n\n",
                          status_msg)
    writeLines(summary_text, paste0(plot_dir, "reg_element/re_pathway_lfs_cancer_summary", param_suffix_re_pathway, ".txt"))

    cat(status_msg, "\n")
  }
} else {
  cat("WARNING: Cancer column not found in te_lfs_re_split, skipping cancer status analysis\n")
}

# Analysis 3: Affected cohort colored by TP53 status (te_aff_expand - same data as Analysis 1)
cat("\n\n=== Regulatory Elements Analysis: Affected Cohort (by TP53 Status) ===\n")

# Check if TP53_status column exists (using te_aff_re_split from Analysis 1)
if ("TP53_status" %in% colnames(te_aff_re_split)) {
  cat("TP53_status groups:\n")
  print(table(te_aff_re_split$TP53_status))

  cat("\nRunning pathway analysis by TP53 status...\n")
  n_genes_input_tp53 <- length(unique(te_aff_re_split$gene_reg))
  n_samples_input_tp53 <- length(unique(te_aff_re_split$sample.x))
  # Count samples per TP53 group
  samples_per_tp53 <- table(unique(te_aff_re_split[, c("sample.x", "TP53_status")])$TP53_status)
  samples_per_tp53_str <- paste(names(samples_per_tp53), "=", samples_per_tp53, "samples", collapse=", ")

  # Write data summary before analysis
  genes_per_tp53 <- te_aff_re_split %>%
    group_by(TP53_status) %>%
    summarise(n_genes = n_distinct(gene_reg), .groups = "drop")
  genes_per_tp53_str <- paste(genes_per_tp53$TP53_status, "=", genes_per_tp53$n_genes, "genes", collapse=", ")

  summary_text <- paste0("RE Pathway Analysis Summary - Affected TP53 Status\n",
                        "===================================================\n\n",
                        "Total samples: ", n_samples_input_tp53, "\n",
                        "Samples per group: ", samples_per_tp53_str, "\n",
                        "Total genes: ", n_genes_input_tp53, "\n",
                        "Genes per group: ", genes_per_tp53_str, "\n",
                        "Parameters: min_samples=", params_re$min_samples_per_group,
                        ", p<", params_re$p_pathway, ", q<", params_re$q_pathway, "\n\n",
                        "Status: Running analysis...")
  writeLines(summary_text, paste0(plot_dir, "reg_element/re_pathway_aff_tp53_summary", param_suffix_re_pathway, ".txt"))

  ora_re_aff_tp53 <- perform_re_pathway_analysis(te_aff_re_split, analysis_type = "tp53", min_samples = params_re$min_samples_per_group, p_pathway = params_re$p_pathway, q_pathway = params_re$q_pathway)

  if (!is.null(ora_re_aff_tp53) && (inherits(ora_re_aff_tp53, "enrichResult") || inherits(ora_re_aff_tp53, "compareClusterResult")) && nrow(as.data.frame(ora_re_aff_tp53)) > 0) {
    cat("Pathway analysis successful! Found", nrow(as.data.frame(ora_re_aff_tp53)), "significant pathways\n")
    write.csv(as.data.frame(ora_re_aff_tp53), paste0(r_dir_files, "re_pathway_aff_tp53", param_suffix_re_pathway, ".csv"), row.names=FALSE)

    # Update summary with results
    summary_text <- paste0("RE Pathway Analysis Summary - Affected TP53 Status\n",
                          "===================================================\n\n",
                          "Total samples: ", n_samples_input_tp53, "\n",
                          "Samples per group: ", samples_per_tp53_str, "\n",
                          "Total genes: ", n_genes_input_tp53, "\n",
                          "Genes per group: ", genes_per_tp53_str, "\n",
                          "Parameters: min_samples=", params_re$min_samples_per_group,
                          ", p<", params_re$p_pathway, ", q<", params_re$q_pathway, "\n\n",
                          "Status: SUCCESS\n",
                          "Significant pathways found: ", nrow(as.data.frame(ora_re_aff_tp53)))
    writeLines(summary_text, paste0(plot_dir, "reg_element/re_pathway_aff_tp53_summary", param_suffix_re_pathway, ".txt"))

    cat("Creating dotplot...\n")
    p_dot_ora_re_tp53 <- dotplot(ora_re_aff_tp53, showCategory=20)
    titled_print(p_dot_ora_re_tp53, "ORA dotplot (RE - TP53 Status)")
    ggsave(paste0(plot_dir, "reg_element/re_pathway_tp53_dot", param_suffix_re_pathway, ".png"), plot=p_dot_ora_re_tp53, width=14, height=9)

    cat("Creating cnetplot...\n")
    p_cnet_ora_re_tp53 <- cnetplot(ora_re_aff_tp53, showCategory=10, colorEdge=TRUE, node_label="category")
    titled_print(p_cnet_ora_re_tp53, "RE ORA cnetplot (TP53 Status)")
    ggsave(paste0(plot_dir, "reg_element/re_pathway_tp53_cnet", param_suffix_re_pathway, ".png"), plot=p_cnet_ora_re_tp53, width=14, height=9)

    cat("Creating emapplot...\n")
    tryCatch({
      # Check if multiple clusters have results
      cluster_counts <- table(as.data.frame(ora_re_aff_tp53)$Cluster)
      if (length(cluster_counts) > 1) {
        ora_re_tp53_pairwise <- pairwise_termsim(ora_re_aff_tp53)
        p_emap_ora_re_tp53 <- emapplot(ora_re_tp53_pairwise, showCategory=20,
                                        pie.params = list(pie = "count"),
                                        cluster.params = list(cluster = TRUE, legend = TRUE))
        titled_print(p_emap_ora_re_tp53, "ORA emapplot (RE - TP53)")
        ggsave(paste0(plot_dir, "reg_element/re_pathway_re_tp53", param_suffix_re_pathway, ".png"), plot=p_emap_ora_re_tp53, width=14, height=9)
      } else {
        cat("Skipping emapplot: Only one cluster has results (", names(cluster_counts), ")\n")
      }
    }, error = function(e) {
      cat("Warning: Could not create emapplot:", e$message, "\n")
    })

    # Simplified pathway plots
    cat("\nCreating simplified pathway plots...\n")
    ora_re_aff_tp53_simple <- simplify(ora_re_aff_tp53, cutoff=0.5, by="p.adjust", select_fun=min)
    write.csv(as.data.frame(ora_re_aff_tp53_simple), paste0(r_dir_files, "re_pathway_aff_tp53_simple", param_suffix_re_pathway, ".csv"), row.names=FALSE)

    if (!is.null(ora_re_aff_tp53_simple) && nrow(as.data.frame(ora_re_aff_tp53_simple)) > 0) {
      # Dot (no bar for compareCluster)
      p_dot_ora_re_tp53_simple <- dotplot(ora_re_aff_tp53_simple, showCategory=20)
      titled_print(p_dot_ora_re_tp53_simple, "ORA dotplot simplified (RE - TP53 Status)")
      ggsave(paste0(plot_dir, "reg_element/re_pathway_tp53_dot_simplified", param_suffix_re_pathway, ".png"), plot=p_dot_ora_re_tp53_simple, width=14, height=9)

      # Cnet
      p_cnet_ora_re_tp53_simple <- cnetplot(ora_re_aff_tp53_simple, showCategory=10, colorEdge=TRUE, node_label="category")
      titled_print(p_cnet_ora_re_tp53_simple, "RE ORA cnetplot simplified (TP53 Status)")
      ggsave(paste0(plot_dir, "reg_element/re_pathway_tp53_cnet_simplified", param_suffix_re_pathway, ".png"), plot=p_cnet_ora_re_tp53_simple, width=14, height=9)

      # Emap
      tryCatch({
        cluster_counts_simple <- table(as.data.frame(ora_re_aff_tp53_simple)$Cluster)
        if (length(cluster_counts_simple) > 1) {
          ora_re_tp53_simple_pairwise <- pairwise_termsim(ora_re_aff_tp53_simple)
          p_emap_ora_re_tp53_simple <- emapplot(ora_re_tp53_simple_pairwise, showCategory=20,
                                                 pie.params = list(pie = "count"),
                                                 cluster.params = list(cluster = TRUE, legend = TRUE))
          titled_print(p_emap_ora_re_tp53_simple, "ORA emapplot simplified (RE - TP53)")
          ggsave(paste0(plot_dir, "reg_element/re_pathway_tp53_emap_simplified", param_suffix_re_pathway, ".png"), plot=p_emap_ora_re_tp53_simple, width=14, height=9)
        } else {
          cat("Skipping simplified emapplot: Only one cluster has results\n")
        }
      }, error = function(e) {
        cat("Warning: Could not create simplified emapplot:", e$message, "\n")
      })
    } else {
      cat("No significant pathways in simplified results\n")
    }
  } else {
    # Update summary with failure reason
    if (is.null(ora_re_aff_tp53)) {
      status_msg <- paste0("Status: INSUFFICIENT DATA\n",
                          "Reason: Not enough samples or genes to run pathway analysis")
    } else {
      status_msg <- paste0("Status: NO SIGNIFICANT PATHWAYS\n",
                          "Analysis completed but no pathways met significance thresholds")
    }

    summary_text <- paste0("RE Pathway Analysis Summary - Affected TP53 Status\n",
                          "===================================================\n\n",
                          "Total samples: ", n_samples_input_tp53, "\n",
                          "Samples per group: ", samples_per_tp53_str, "\n",
                          "Total genes: ", n_genes_input_tp53, "\n",
                          "Genes per group: ", genes_per_tp53_str, "\n",
                          "Parameters: min_samples=", params_re$min_samples_per_group,
                          ", p<", params_re$p_pathway, ", q<", params_re$q_pathway, "\n\n",
                          status_msg)
    writeLines(summary_text, paste0(plot_dir, "reg_element/re_pathway_aff_tp53_summary", param_suffix_re_pathway, ".txt"))

    cat(status_msg, "\n")
  }
} else {
  cat("WARNING: TP53_status column not found in te_aff_re_split, skipping TP53 analysis\n")
}

# Analysis 4: KICS + HostSeq colored by cohort (te_kics_hostseq_expand)
cat("\n\n=== Regulatory Elements Analysis: KICS + HostSeq (by Cohort) ===\n")
te_kics_hostseq_re <- load_and_join_re_data(te_kics_hostseq_expand, re_germline_path)
te_kics_hostseq_re_split <- split_re_genes(te_kics_hostseq_re)

cat("\nRunning pathway analysis by cohort (KICS vs HostSeq)...\n")
n_genes_input_kics <- length(unique(te_kics_hostseq_re_split$gene_reg))
n_samples_input_kics <- length(unique(te_kics_hostseq_re_split$sample.x))
# Count samples per cohort
samples_per_kics_cancer <- table(unique(te_kics_hostseq_re_split[, c("sample.x", "cohort")])$cohort)
samples_per_kics_cancer_str <- paste(names(samples_per_kics_cancer), "=", samples_per_kics_cancer, "samples", collapse=", ")

# Write data summary before analysis
genes_per_cohort <- te_kics_hostseq_re_split %>%
  group_by(cohort) %>%
  summarise(n_genes = n_distinct(gene_reg), .groups = "drop")
genes_per_cohort_str <- paste(genes_per_cohort$cohort, "=", genes_per_cohort$n_genes, "genes", collapse=", ")

summary_text <- paste0("RE Pathway Analysis Summary - KICS + HostSeq Cohorts\n",
                      "=====================================================\n\n",
                      "Total samples: ", n_samples_input_kics, "\n",
                      "Samples per group: ", samples_per_kics_cancer_str, "\n",
                      "Total genes: ", n_genes_input_kics, "\n",
                      "Genes per group: ", genes_per_cohort_str, "\n",
                      "Parameters: min_samples=", params_re$min_samples_per_group,
                      ", p<", params_re$p_pathway, ", q<", params_re$q_pathway, "\n\n",
                      "Status: Running analysis...")
writeLines(summary_text, paste0(plot_dir, "reg_element/re_pathway_kics_hostseq_cohort_summary", param_suffix_re_pathway, ".txt"))

ora_re_kics_cancer <- perform_re_pathway_analysis(te_kics_hostseq_re_split, analysis_type = "cancer", min_samples = params_re$min_samples_per_group, p_pathway = params_re$p_pathway, q_pathway = params_re$q_pathway)

if (!is.null(ora_re_kics_cancer) && (inherits(ora_re_kics_cancer, "enrichResult") || inherits(ora_re_kics_cancer, "compareClusterResult")) && nrow(as.data.frame(ora_re_kics_cancer)) > 0) {
  cat("Pathway analysis successful! Found", nrow(as.data.frame(ora_re_kics_cancer)), "significant pathways\n")
  write.csv(as.data.frame(ora_re_kics_cancer), paste0(r_dir_files, "re_pathway_kics_hostseq_cohort", param_suffix_re_pathway, ".csv"), row.names=FALSE)

  # Update summary with results
  summary_text <- paste0("RE Pathway Analysis Summary - KICS + HostSeq Cohorts\n",
                        "=====================================================\n\n",
                        "Total samples: ", n_samples_input_kics, "\n",
                        "Samples per group: ", samples_per_kics_cancer_str, "\n",
                        "Total genes: ", n_genes_input_kics, "\n",
                        "Genes per group: ", genes_per_cohort_str, "\n",
                        "Parameters: min_samples=", params_re$min_samples_per_group,
                        ", p<", params_re$p_pathway, ", q<", params_re$q_pathway, "\n\n",
                        "Status: SUCCESS\n",
                        "Significant pathways found: ", nrow(as.data.frame(ora_re_kics_cancer)))
  writeLines(summary_text, paste0(plot_dir, "reg_element/re_pathway_kics_hostseq_cohort_summary", param_suffix_re_pathway, ".txt"))

  cat("Creating dotplot...\n")
  p_dot_ora_re_kics <- dotplot(ora_re_kics_cancer, showCategory=20)
  titled_print(p_dot_ora_re_kics, "ORA dotplot (RE - KICS vs HostSeq)")
  ggsave(paste0(plot_dir, "reg_element/re_pathway_kics_hostseq_cohort_dot", param_suffix_re_pathway, ".png"), plot=p_dot_ora_re_kics, width=14, height=9)

  cat("Creating cnetplot...\n")
  p_cnet_ora_re_kics <- cnetplot(ora_re_kics_cancer, showCategory=10, colorEdge=TRUE, node_label="category")
  titled_print(p_cnet_ora_re_kics, "RE ORA cnetplot (KICS vs HostSeq)")
  ggsave(paste0(plot_dir, "reg_element/re_pathway_kics_hostseq_cohort_cnet", param_suffix_re_pathway, ".png"), plot=p_cnet_ora_re_kics, width=14, height=9)

  cat("Creating emapplot...\n")
  tryCatch({
    # Check if multiple clusters have results
    cluster_counts <- table(as.data.frame(ora_re_kics_cancer)$Cluster)
    if (length(cluster_counts) > 1) {
      ora_re_kics_pairwise <- pairwise_termsim(ora_re_kics_cancer)
      p_emap_ora_re_kics <- emapplot(ora_re_kics_pairwise, showCategory=20,
                                      pie.params = list(pie = "count"),
                                      cluster.params = list(cluster = TRUE, legend = TRUE))
      titled_print(p_emap_ora_re_kics, "ORA emapplot (RE - KICS vs HostSeq)")
      ggsave(paste0(plot_dir, "reg_element/re_pathway_kics_hostseq_cohort_emap", param_suffix_re_pathway, ".png"), plot=p_emap_ora_re_kics, width=14, height=9)
    } else {
      cat("Skipping emapplot: Only one cluster has results (", names(cluster_counts), ")\n")
    }
  }, error = function(e) {
    cat("Warning: Could not create emapplot:", e$message, "\n")
  })

  # Simplified pathway plots
  cat("\nCreating simplified pathway plots...\n")
  ora_re_kics_cancer_simple <- simplify(ora_re_kics_cancer, cutoff=0.5, by="p.adjust", select_fun=min)
  write.csv(as.data.frame(ora_re_kics_cancer_simple), paste0(r_dir_files, "re_pathway_kics_hostseq_cohort_simple", param_suffix_re_pathway, ".csv"), row.names=FALSE)

  if (!is.null(ora_re_kics_cancer_simple) && nrow(as.data.frame(ora_re_kics_cancer_simple)) > 0) {
    # Dot (no bar for compareCluster)
    p_dot_ora_re_kics_simple <- dotplot(ora_re_kics_cancer_simple, showCategory=20)
    titled_print(p_dot_ora_re_kics_simple, "ORA dotplot simplified (RE - KICS vs HostSeq)")
    ggsave(paste0(plot_dir, "reg_element/re_pathway_kics_hostseq_cohort_dot_simplified", param_suffix_re_pathway, ".png"), plot=p_dot_ora_re_kics_simple, width=14, height=9)

    # Cnet
    p_cnet_ora_re_kics_simple <- cnetplot(ora_re_kics_cancer_simple, showCategory=10, colorEdge=TRUE, node_label="category")
    titled_print(p_cnet_ora_re_kics_simple, "RE ORA cnetplot simplified (KICS vs HostSeq)")
    ggsave(paste0(plot_dir, "reg_element/re_pathway_kics_hostseq_cohort_cnet_simplified", param_suffix_re_pathway, ".png"), plot=p_cnet_ora_re_kics_simple, width=14, height=9)

    # Emap
    tryCatch({
      cluster_counts_simple <- table(as.data.frame(ora_re_kics_cancer_simple)$Cluster)
      if (length(cluster_counts_simple) > 1) {
        ora_re_kics_simple_pairwise <- pairwise_termsim(ora_re_kics_cancer_simple)
        p_emap_ora_re_kics_simple <- emapplot(ora_re_kics_simple_pairwise, showCategory=20,
                                               pie.params = list(pie = "count"),
                                               cluster.params = list(cluster = TRUE, legend = TRUE))
        titled_print(p_emap_ora_re_kics_simple, "ORA emapplot simplified (RE - KICS vs HostSeq)")
        ggsave(paste0(plot_dir, "reg_element/re_pathway_kics_hostseq_cohort_emap_simplified", param_suffix_re_pathway, ".png"), plot=p_emap_ora_re_kics_simple, width=14, height=9)
      } else {
        cat("Skipping simplified emapplot: Only one cluster has results\n")
      }
    }, error = function(e) {
      cat("Warning: Could not create simplified emapplot:", e$message, "\n")
    })
  } else {
    cat("No significant pathways in simplified results\n")
  }
} else {
  # Update summary with failure status
  if (is.null(ora_re_kics_cancer)) {
    status_msg <- paste0("Status: INSUFFICIENT DATA\n",
                        "Reason: Not enough samples per group to run pathway analysis\n",
                        "Minimum required: ", params_re$min_samples_per_group, " samples per group")
  } else {
    status_msg <- paste0("Status: NO SIGNIFICANT PATHWAYS\n",
                        "Analysis completed but no pathways met significance thresholds")
  }

  summary_text <- paste0("RE Pathway Analysis Summary - KICS + HostSeq Cohorts\n",
                        "=====================================================\n\n",
                        "Total samples: ", n_samples_input_kics, "\n",
                        "Samples per group: ", samples_per_kics_cancer_str, "\n",
                        "Total genes: ", n_genes_input_kics, "\n",
                        "Genes per group: ", genes_per_cohort_str, "\n",
                        "Parameters: min_samples=", params_re$min_samples_per_group,
                        ", p<", params_re$p_pathway, ", q<", params_re$q_pathway, "\n\n",
                        status_msg)
  writeLines(summary_text, paste0(plot_dir, "reg_element/re_pathway_kics_hostseq_cohort_summary", param_suffix_re_pathway, ".txt"))

  cat(status_msg, "\n")
}

# Analysis 5: KICS by Sample Type (Blood, Fibroblasts, Tissue (fresh))
cat("\n\n=== Regulatory Elements Analysis: KICS (by Sample Type) ===\n")

# Load and prepare KICS sample type data
kics_sample_type <- prep_kics_sample_type(
  "/Users/briannelaverty/Documents/R_Malkin/clinical/kics_germline_sample_type.csv"
)

# Load RE data for KICS expand
te_kics_expand_re <- load_and_join_re_data(te_kics_expand, re_germline_path)
te_kics_expand_re_split <- split_re_genes(te_kics_expand_re)

# Merge sample type with RE data (merge by sample column)
# Note: After RE join, sample column is named 'sample.x' due to duplicate columns
te_kics_re_sampletype <- merge(te_kics_expand_re_split, kics_sample_type, by.x = "sample.x", by.y = "sample", all.x = TRUE)

cat("Merged RE data with sample type:\n")
cat("  Total rows:", nrow(te_kics_re_sampletype), "\n")
cat("  Samples with sample type data:", sum(!is.na(te_kics_re_sampletype$sample_type)), "\n")

# Filter for Blood, Fibroblasts, Tissue (fresh)
te_kics_re_sampletype <- te_kics_re_sampletype %>%
  filter(sample_type %in% c("Blood", "Fibroblasts", "Tissue (fresh)"))

cat("Samples in RE analysis by sample type:",
    length(unique(te_kics_re_sampletype$sample)), "\n")
cat("Sample type distribution:\n")
print(table(te_kics_re_sampletype$sample_type))

cat("\nRunning pathway analysis by sample type...\n")
n_genes_input_sampletype <- length(unique(te_kics_re_sampletype$gene_reg))
n_samples_input_sampletype <- length(unique(te_kics_re_sampletype$sample.x))
# Count samples per sample type
samples_per_sampletype <- table(unique(te_kics_re_sampletype[, c("sample.x", "sample_type")])$sample_type)
samples_per_sampletype_str <- paste(names(samples_per_sampletype), "=", samples_per_sampletype, "samples", collapse=", ")

# Count genes per sample type BEFORE running analysis
genes_per_sampletype <- te_kics_re_sampletype %>%
  group_by(sample_type) %>%
  summarise(n_genes = n_distinct(gene_reg), .groups = "drop")
genes_per_sampletype_str <- paste(genes_per_sampletype$sample_type, "=", genes_per_sampletype$n_genes, "genes", collapse=", ")

# Write data summary before analysis
summary_text <- paste0("RE Pathway Analysis Summary - KICS Sample Type\n",
                      "===============================================\n\n",
                      "Total samples: ", n_samples_input_sampletype, "\n",
                      "Samples per group: ", samples_per_sampletype_str, "\n",
                      "Total genes: ", n_genes_input_sampletype, "\n",
                      "Genes per group: ", genes_per_sampletype_str, "\n",
                      "Parameters: min_samples=", params_re$min_samples_per_group,
                      ", p<", params_re$p_pathway, ", q<", params_re$q_pathway, "\n\n",
                      "Status: Running analysis...")
writeLines(summary_text, paste0(plot_dir, "reg_element/re_pathway_kics_sampletype_summary", param_suffix_re_pathway, ".txt"))

ora_re_sampletype <- perform_ora_sample_type_custom_cutoffs(te_kics_re_sampletype, nsample_thresh = params_re$min_samples_per_group, p_pathway = params_re$p_pathway, q_pathway = params_re$q_pathway, sample_type_column = "sample_type", gene_col = "gene_reg")

if (!is.null(ora_re_sampletype) && (inherits(ora_re_sampletype, "enrichResult") || inherits(ora_re_sampletype, "compareClusterResult")) && nrow(as.data.frame(ora_re_sampletype)) > 0) {
  cat("Pathway analysis successful! Found", nrow(as.data.frame(ora_re_sampletype)), "significant pathways\n")
  write.csv(as.data.frame(ora_re_sampletype), paste0(r_dir_files, "re_pathway_kics_sampletype", param_suffix_re_pathway, ".csv"), row.names=FALSE)

  # Update summary with results
  summary_text <- paste0("RE Pathway Analysis Summary - KICS Sample Type\n",
                        "===============================================\n\n",
                        "Total samples: ", n_samples_input_sampletype, "\n",
                        "Samples per group: ", samples_per_sampletype_str, "\n",
                        "Total genes: ", n_genes_input_sampletype, "\n",
                        "Genes per group: ", genes_per_sampletype_str, "\n",
                        "Parameters: min_samples=", params_re$min_samples_per_group,
                        ", p<", params_re$p_pathway, ", q<", params_re$q_pathway, "\n\n",
                        "Status: SUCCESS\n",
                        "Significant pathways found: ", nrow(as.data.frame(ora_re_sampletype)))
  writeLines(summary_text, paste0(plot_dir, "reg_element/re_pathway_kics_sampletype_summary", param_suffix_re_pathway, ".txt"))

  cat("Creating dotplot...\n")
  p_dot_ora_re_sampletype <- dotplot(ora_re_sampletype, showCategory=20)
  titled_print(p_dot_ora_re_sampletype, "ORA dotplot (RE - KICS by Sample Type)")
  ggsave(paste0(plot_dir, "reg_element/re_pathway_kics_sampletype_dot", param_suffix_re_pathway, ".png"), plot=p_dot_ora_re_sampletype, width=14, height=9)

  cat("Creating cnetplot...\n")
  p_cnet_ora_re_sampletype <- cnetplot(ora_re_sampletype, showCategory=10, colorEdge=TRUE, node_label="category")
  titled_print(p_cnet_ora_re_sampletype, "RE ORA cnetplot (KICS by Sample Type)")
  ggsave(paste0(plot_dir, "reg_element/re_pathway_kics_sampletype_cnet", param_suffix_re_pathway, ".png"), plot=p_cnet_ora_re_sampletype, width=14, height=9)

  cat("Creating emapplot...\n")
  tryCatch({
    # Check if multiple clusters have results
    cluster_counts <- table(as.data.frame(ora_re_sampletype)$Cluster)
    if (length(cluster_counts) > 1) {
      ora_re_sampletype_pairwise <- pairwise_termsim(ora_re_sampletype)
      p_emap_ora_re_sampletype <- emapplot(ora_re_sampletype_pairwise, showCategory=20,
                                           pie.params = list(pie = "count"),
                                           cluster.params = list(cluster = TRUE, legend = TRUE))
      titled_print(p_emap_ora_re_sampletype, "ORA emapplot (RE - KICS by Sample Type)")
      ggsave(paste0(plot_dir, "reg_element/re_pathway_kics_sampletype_emap", param_suffix_re_pathway, ".png"), plot=p_emap_ora_re_sampletype, width=14, height=9)
    } else {
      cat("Skipping emapplot: Only one cluster has results (", names(cluster_counts), ")\n")
    }
  }, error = function(e) {
    cat("Warning: Could not create emapplot:", e$message, "\n")
  })

  # Simplified pathway plots
  cat("\nCreating simplified pathway plots...\n")
  ora_re_sampletype_simple <- simplify(ora_re_sampletype, cutoff=0.5, by="p.adjust", select_fun=min)
  write.csv(as.data.frame(ora_re_sampletype_simple), paste0(r_dir_files, "re_pathway_kics_sampletype_simple", param_suffix_re_pathway, ".csv"), row.names=FALSE)

  if (!is.null(ora_re_sampletype_simple) && nrow(as.data.frame(ora_re_sampletype_simple)) > 0) {
    # Dot (no bar for compareCluster)
    p_dot_ora_re_sampletype_simple <- dotplot(ora_re_sampletype_simple, showCategory=20)
    titled_print(p_dot_ora_re_sampletype_simple, "ORA dotplot simplified (RE - KICS by Sample Type)")
    ggsave(paste0(plot_dir, "reg_element/re_pathway_kics_sampletype_dot_simplified", param_suffix_re_pathway, ".png"), plot=p_dot_ora_re_sampletype_simple, width=14, height=9)

    # Cnet
    p_cnet_ora_re_sampletype_simple <- cnetplot(ora_re_sampletype_simple, showCategory=10, colorEdge=TRUE, node_label="category")
    titled_print(p_cnet_ora_re_sampletype_simple, "RE ORA cnetplot simplified (KICS by Sample Type)")
    ggsave(paste0(plot_dir, "reg_element/re_pathway_kics_sampletype_cnet_simplified", param_suffix_re_pathway, ".png"), plot=p_cnet_ora_re_sampletype_simple, width=14, height=9)

    # Emap
    tryCatch({
      cluster_counts_simple <- table(as.data.frame(ora_re_sampletype_simple)$Cluster)
      if (length(cluster_counts_simple) > 1) {
        ora_re_sampletype_simple_pairwise <- pairwise_termsim(ora_re_sampletype_simple)
        p_emap_ora_re_sampletype_simple <- emapplot(ora_re_sampletype_simple_pairwise, showCategory=20,
                                                     pie.params = list(pie = "count"),
                                                     cluster.params = list(cluster = TRUE, legend = TRUE))
        titled_print(p_emap_ora_re_sampletype_simple, "ORA emapplot simplified (RE - KICS by Sample Type)")
        ggsave(paste0(plot_dir, "reg_element/re_pathway_kics_sampletype_emap_simplified", param_suffix_re_pathway, ".png"), plot=p_emap_ora_re_sampletype_simple, width=14, height=9)
      } else {
        cat("Skipping simplified emapplot: Only one cluster has results\n")
      }
    }, error = function(e) {
      cat("Warning: Could not create simplified emapplot:", e$message, "\n")
    })
  } else {
    cat("No significant pathways in simplified results\n")
  }
} else {
  # Update summary with failure status
  if (is.null(ora_re_sampletype)) {
    status_msg <- paste0("Status: INSUFFICIENT DATA\n",
                        "Reason: Not enough samples per group to run pathway analysis\n",
                        "Minimum required: ", params_re$min_samples_per_group, " samples per group")
  } else {
    status_msg <- paste0("Status: NO SIGNIFICANT PATHWAYS\n",
                        "Analysis completed but no pathways met significance thresholds")
  }

  summary_text <- paste0("RE Pathway Analysis Summary - KICS Sample Type\n",
                        "===============================================\n\n",
                        "Total samples: ", n_samples_input_sampletype, "\n",
                        "Samples per group: ", samples_per_sampletype_str, "\n",
                        "Total genes: ", n_genes_input_sampletype, "\n",
                        "Genes per group: ", genes_per_sampletype_str, "\n",
                        "Parameters: min_samples=", params_re$min_samples_per_group,
                        ", p<", params_re$p_pathway, ", q<", params_re$q_pathway, "\n\n",
                        status_msg)
  writeLines(summary_text, paste0(plot_dir, "reg_element/re_pathway_kics_sampletype_summary", param_suffix_re_pathway, ".txt"))

  cat(status_msg, "\n")
}

} # End of RE pathway parameter loop


#### RE-RNA DIFFERENTIAL EXPRESSION TEST WITH PARAMETER SWEEP ####
write_output(quote(NULL), "Running RE-RNA Differential Expression Analysis")

# Check if RNA data is available (requires running RNA script first)
if (!exists("rna_filtered")) {
  # Try to load saved rna_filtered from RNA script
  rna_filtered_file <- paste0(r_dir_files, "rna_filtered_germline.rds")
  if (file.exists(rna_filtered_file)) {
    cat("Loading rna_filtered from previously run RNA script...\n")
    rna_filtered <- readRDS(rna_filtered_file)
    cat("✓ Successfully loaded rna_filtered\n\n")
  } else {
    cat("⚠ WARNING: RNA data not processed (rna_filtered not found).\n")
    cat("  The RE-RNA analysis requires RNA processing from 02_te_viz_germline_06_rna.R\n")
    cat("  Skipping RE-RNA differential expression analysis.\n")
    cat("  To run this section, first run the RNA script or run the ALL script.\n\n")
  }
}

if (exists("rna_filtered")) {
  # Set up summary tracking for parameter sweep
  sweep_summary <- data.frame()

  re_report_path_germline <- "/Users/briannelaverty/Documents/R_Malkin/te/data/final/germline_annotSV_output.SV_RE_intersect.report"

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

    cat("\n=== Testing parameters", i, "/", nrow(param_grid_re_rna), "===\n")
    cat("min_samples_per_group =", params$min_samples_per_group,
        ", p_gene =", params$p_gene,
        ", p_pathway =", params$p_pathway,
        ", q_pathway =", params$q_pathway, "\n")

    # Run differential expression test
    re_rna_results <- test_rna_by_gene_re_status(
      re_report_path = re_report_path_germline,
      rna_data = rna_filtered,
      sample_type = "germline",
      min_gene_mentions = 1,
      min_samples_per_group = params$min_samples_per_group
    )

    # Save differential expression results
    if (nrow(re_rna_results) > 0) {
      output_file <- paste0(r_dir_files, "re_rna_differential_germline", param_suffix, ".csv")
      write.csv(re_rna_results, output_file, row.names = FALSE)
      cat("✓ Results saved to:", basename(output_file), "\n")

      # Filter genes by p_gene threshold
      pathway_genes <- re_rna_results %>%
        filter(p_value < params$p_gene) %>%
        pull(gene) %>%
        unique()

      cat("Genes with p-value <", params$p_gene, ":", length(pathway_genes), "\n")

      # Pathway analysis if enough genes
      if (length(pathway_genes) >= 5) {
        pathway_df <- data.frame(Gene_name = pathway_genes)

        ora_re_rna <- tryCatch({
          perform_ora_custom_cutoffs(pathway_df,
                                    p_pathway = params$p_pathway,
                                    q_pathway = params$q_pathway,
                                    nsample_thresh = 0,
                                    filter_exon = FALSE)
        }, error = function(e) {
          cat("Error in pathway analysis:", e$message, "\n")
          NULL
        })

        if (!is.null(ora_re_rna) && nrow(as.data.frame(ora_re_rna)) > 0) {
          n_pathways <- nrow(as.data.frame(ora_re_rna))
          cat("✓ Found", n_pathways, "enriched pathways\n")

          # Save pathway results
          pathway_file <- paste0(r_dir_files, "re_rna_pathway_germline", param_suffix, ".csv")
          write.csv(as.data.frame(ora_re_rna), pathway_file, row.names = FALSE)

          # Save gene list
          gene_list_file <- paste0(r_dir_files, "re_rna_genes_germline", param_suffix, ".txt")
          writeLines(pathway_genes, gene_list_file)

          # Create plots
          tryCatch({
            p_bar_re_rna <- barplot(ora_re_rna, showCategory = 20)
            titled_print(p_bar_re_rna, "RE-RNA Pathway Barplot (Germline)")
            ggsave(paste0(plot_dir, "reg_element/re_rna_pathway_bar_germline", param_suffix, ".png"),
                   plot = p_bar_re_rna, width = 14, height = 9)
          }, error = function(e) cat("Warning: Could not create barplot:", e$message, "\n"))

          tryCatch({
            p_dot_re_rna <- dotplot(ora_re_rna, showCategory = 20)
            titled_print(p_dot_re_rna, "RE-RNA Pathway Dotplot (Germline)")
            ggsave(paste0(plot_dir, "reg_element/re_rna_pathway_dot_germline", param_suffix, ".png"),
                   plot = p_dot_re_rna, width = 14, height = 9)
          }, error = function(e) cat("Warning: Could not create dotplot:", e$message, "\n"))

          tryCatch({
            p_cnet_re_rna <- cnetplot(ora_re_rna, showCategory = 10, colorEdge = TRUE, node_label = "category")
            titled_print(p_cnet_re_rna, "RE-RNA Pathway Cnetplot (Germline)")
            ggsave(paste0(plot_dir, "reg_element/re_rna_pathway_cnet_germline", param_suffix, ".png"),
                   plot = p_cnet_re_rna, width = 14, height = 9)
          }, error = function(e) cat("Warning: Could not create cnetplot:", e$message, "\n"))

          tryCatch({
            ora_re_rna_pairwise <- pairwise_termsim(ora_re_rna)
            p_emap_re_rna <- emapplot(ora_re_rna_pairwise, showCategory = 20)
            titled_print(p_emap_re_rna, "RE-RNA Pathway Emapplot (Germline)")
            ggsave(paste0(plot_dir, "reg_element/re_rna_pathway_emap_germline", param_suffix, ".png"),
                   plot = p_emap_re_rna, width = 14, height = 9)
          }, error = function(e) cat("Warning: Could not create emapplot:", e$message, "\n"))

          # Record sweep summary
          sweep_summary <- rbind(sweep_summary, data.frame(
            min_samples_per_group = params$min_samples_per_group,
            p_gene = params$p_gene,
            p_pathway = params$p_pathway,
            q_pathway = params$q_pathway,
            n_de_genes = nrow(re_rna_results),
            n_filtered_genes = length(pathway_genes),
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
            n_de_genes = nrow(re_rna_results),
            n_filtered_genes = length(pathway_genes),
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
          n_de_genes = nrow(re_rna_results),
          n_filtered_genes = length(pathway_genes),
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
    summary_file <- paste0(r_dir_files, "re_rna_sweep_germline_summary.csv")
    write.csv(sweep_summary, summary_file, row.names = FALSE)
    cat("\n✓ Parameter sweep complete! Summary saved to:", summary_file, "\n")
    cat("\nTop 5 parameter combinations by number of pathways:\n")
    print(head(sweep_summary[order(-sweep_summary$n_pathways), ], 5))
  }
}



cat("✓ Script completed successfully\n")
