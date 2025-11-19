#!/usr/bin/env Rscript

# Germline TE Visualization - Pathway Analysis
# All pathway analyses: TP53-specific, Cancer-specific, KICS, TP53, PedCancer

# Source common setup and load data
source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/00_viz_common_setup.R")
source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/00_viz_load_data_germline.R")

cat("Running 02_te_viz_germline_05_pathway.R...\n")

# Define TE types for analysis
types <- c(NA, "LINE1", "ALU", "SVA")
x_tp53 <- expression("Germline " * italic("TP53") * " status")

# Load geneList from cancer genes script (needed for gene-pathway correlation)
geneList_kics_file <- paste0(r_dir_files, "te_top_genes_kics.csv")
if (file.exists(geneList_kics_file)) {
  geneList_kics <- read.csv(geneList_kics_file, stringsAsFactors = FALSE)
  cat("Loaded geneList_kics from:", geneList_kics_file, "\n")
} else {
  cat("Warning: geneList_kics file not found. Run 02_te_viz_germline_04_cancer_genes.R first.\n")
  cat("File expected at:", geneList_kics_file, "\n")
}

#### PARAMETER SWEEP CONFIGURATION ####

# Define parameter grid for pathway analysis
param_grid <- expand.grid(
  min_samples = c(5),
  p_gene = c(0.05, 0.1),
  p_pathway = c(0.05, 0.1),
  q_pathway = c(0.05, 0.1),
  stringsAsFactors = FALSE
)
cat("\n*** PATHWAY ANALYSIS PARAMETER SWEEP ***\n")
cat("Testing", nrow(param_grid), "parameter combinations\n\n")

#### PATHWAY ANALYSIS - TEs SPECIFIC TO TP53 STATUS ####
write_output(quote(NULL), "Pathway Analysis - TEs Specific to TP53 Status")

# Load specific TE files (will be filtered by parameter sweep)
specific_te_files_tp53 <- list()
for (ms in unique(param_grid$min_samples)) {
  file_path <- paste0(r_dir_files, "specific_tes_germline_aff_min", ms, "_full_results.csv")
  if (file.exists(file_path)) {
    specific_te_files_tp53[[paste0("tp53_min", ms)]] <- read.csv(file_path, stringsAsFactors = FALSE)
    cat("Loaded TP53-specific TEs (min_samples=", ms, "):", nrow(specific_te_files_tp53[[paste0("tp53_min", ms)]]), "TEs\n")
  }
}

# Run parameter sweep for TP53-specific TEs
if (length(specific_te_files_tp53) > 0) {
  tp53_sweep_results <- run_pathway_parameter_sweep(
    data_list = specific_te_files_tp53,
    param_grid = param_grid,
    output_base_dir = plot_dir,
    output_subdir = "pathway",
    csv_dir = r_dir_files,
    analysis_name = "pathway_tp53",
    create_plots = TRUE,
    plot_formats = c("png")
  )

  cat("\nTP53 pathway analysis complete\n")
  cat("Plots saved to:", file.path(plot_dir, "pathway/"), "\n")
  cat("CSV files saved to:", r_dir_files, "\n")
}

#### PATHWAY ANALYSIS - TEs SPECIFIC TO CANCER STATUS (LFS) ####
write_output(quote(NULL), "Pathway Analysis - TEs Specific to Cancer Status (LFS)")

# Load specific TE files (will be filtered by parameter sweep)
specific_te_files_cancer <- list()
for (ms in unique(param_grid$min_samples)) {
  file_path <- paste0(r_dir_files, "specific_tes_germline_cancer_min", ms, "_full_results.csv")
  if (file.exists(file_path)) {
    specific_te_files_cancer[[paste0("cancer_min", ms)]] <- read.csv(file_path, stringsAsFactors = FALSE)
    cat("Loaded Cancer-specific TEs (min_samples=", ms, "):", nrow(specific_te_files_cancer[[paste0("cancer_min", ms)]]), "TEs\n")
  }
}

# Run parameter sweep for Cancer-specific TEs
if (length(specific_te_files_cancer) > 0) {
  cancer_sweep_results <- run_pathway_parameter_sweep(
    data_list = specific_te_files_cancer,
    param_grid = param_grid,
    output_base_dir = plot_dir,
    output_subdir = "pathway",
    csv_dir = r_dir_files,
    analysis_name = "pathway_cancer",
    create_plots = TRUE,
    plot_formats = c("png")
  )

  cat("\nCancer pathway analysis complete\n")
  cat("Plots saved to:", file.path(plot_dir, "pathway/"), "\n")
  cat("CSV files saved to:", r_dir_files, "\n")
}

# LEGACY CODE COMMENTED OUT - Now using parameter sweep system above
# Old single-run code retained for reference but not executed


#### PATHWAY PEDIATRIC ####
write_output(quote(NULL), "Pathway Analysis - Pediatric Cohort (KICS)")
# Parameter values: no min_samples, no p_gene filtering (all genes), pvalueCutoff=0.05, qvalueCutoff=0.1
param_suffix_kics <- "_min0_pgene1_ppathway0.05_qpathway0.1"

ora_kics <- perform_ora(te_kics_split, nsample_thresh = 0, filter_exon = FALSE)
write.csv(as.data.frame(ora_kics), paste0(r_dir_files, "pathway_kics", param_suffix_kics, ".csv"), row.names=FALSE)

# Only plot if ora_kics has results
if (!is.null(ora_kics) && nrow(as.data.frame(ora_kics)) > 0) {
  # viz
  p_bar_ora_kics <- barplot(ora_kics, showCategory=20)
  titled_print(p_bar_ora_kics, "ORA barplot (KICS)")
  ggsave(paste0(plot_dir, "pathway/pathway_kics_bar", param_suffix_kics, ".png"), plot=p_bar_ora_kics, width=14, height=9)
  p_dot_ora_kics <- dotplot(ora_kics, showCategory=20)
  titled_print(p_dot_ora_kics, "ORA dotplot (KICS)")
  ggsave(paste0(plot_dir, "pathway/pathway_kics_dot", param_suffix_kics, ".png"), plot=p_dot_ora_kics, width=14, height=9)
  p_cnet_ora_kics <- cnetplot(ora_kics, showCategory=10, color.params=list(edge=TRUE), node_label="category")
  titled_print(p_cnet_ora_kics, "ORA cnetplot (KICS)")
  ggsave(paste0(plot_dir, "pathway/pathway_kics_cnet", param_suffix_kics, ".png"), plot=p_cnet_ora_kics, width=14, height=9)
  ora_kics_pairwise <- pairwise_termsim(ora_kics)
  p_emap_ora_kics <- emapplot(ora_kics_pairwise, pie="count", showCategory=20, group_category=TRUE, group_legend=TRUE)
  titled_print(p_emap_ora_kics, "ORA emapplot (KICS)")
  ggsave(paste0(plot_dir, "pathway/pathway_kics_emap", param_suffix_kics, ".png"), plot=p_emap_ora_kics, width=14, height=9)

  ora_kics_simple <- simplify(ora_kics, cutoff=0.6, by="p.adjust", select_fun=min)
  write.csv(as.data.frame(ora_kics_simple), paste0(r_dir_files, "pathway_kics_simple", param_suffix_kics, ".csv"), row.names=FALSE)
  if (!is.null(ora_kics_simple) && nrow(as.data.frame(ora_kics_simple)) > 0) {
    p_bar_ora_kics_simple <- barplot(ora_kics_simple, showCategory=20)
    titled_print(p_bar_ora_kics_simple, "ORA barplot (KICS, simplified)")
    ggsave(paste0(plot_dir, "pathway/pathway_kics_bar_simplified", param_suffix_kics, ".png"), plot=p_bar_ora_kics_simple, width=14, height=9)
    p_dot_ora_kics_simple <- dotplot(ora_kics_simple, showCategory=20)
    titled_print(p_dot_ora_kics_simple, "ORA dotplot (KICS, simplified)")
    ggsave(paste0(plot_dir, "pathway/pathway_kics_dot_simplified", param_suffix_kics, ".png"), plot=p_dot_ora_kics_simple, width=14, height=9)
    p_cnet_ora_kics_simple <- cnetplot(ora_kics_simple, showCategory=10, color.params=list(edge=TRUE), node_label="category")
    titled_print(p_cnet_ora_kics_simple, "ORA cnetplot (KICS, simplified)")
    ggsave(paste0(plot_dir, "pathway/pathway_kics_cnet_simplified", param_suffix_kics, ".png"), plot=p_cnet_ora_kics_simple, width=14, height=9)
    ora_kics_simple_pairwise <- pairwise_termsim(ora_kics_simple)
    p_emap_ora_kics_simple <- emapplot(ora_kics_simple_pairwise, pie="count", showCategory=20, group_category=TRUE, group_legend=TRUE)
    titled_print(p_emap_ora_kics_simple, "ORA emapplot (KICS, simplified)")
    ggsave(paste0(plot_dir, "pathway/pathway_kics_emap_simplified", param_suffix_kics, ".png"), plot=p_emap_ora_kics_simple, width=14, height=9)
  }

  # pathway gene counts
  pathway_genes <- plot_pathway_gene_counts(ora_kics, n_descriptions=nrow(ora_kics), n_genes=50)
  write.csv(pathway_genes, paste0(r_dir_files, "pathway_genes_kics", param_suffix_kics, ".csv"), row.names=FALSE)

  # correlation with top mutated genes (requires geneList_kics from cancer_genes script)
  p_pathway_corr <- analyze_gene_mutations(geneList_kics, ora_kics, n_descriptions=nrow(ora_kics), column = "num_samples_effected", y_lab="Number of samples with insertion in gene")
  titled_print(p_pathway_corr, "Pathway genes vs number of samples (KICS)")
  ggsave(paste0(plot_dir, "pathway/pathway_kics_genes_vs_numsamples", param_suffix_kics, ".png"), plot=p_pathway_corr, width=9, height=5)
} else {
  warning("No significant enrichment found in ora_kics; skipping pathway plots.")
}



#### PATHWAY TP53 ####
write_output(quote(NULL), "Pathway Analysis - TP53 Cohort (LFS)")
# Parameter values: no min_samples, no p_gene filtering (all genes), pvalueCutoff=0.05, qvalueCutoff=0.1
param_suffix_tp53 <- "_min0_pgene1_ppathway0.05_qpathway0.1"

ora_tp53 <- perform_ora_tp53(te_aff_split)
write.csv(as.data.frame(ora_tp53), paste0(r_dir_files, "pathway_tp53", param_suffix_tp53, ".csv"), row.names=FALSE)

if (!is.null(ora_tp53) && nrow(as.data.frame(ora_tp53)) > 0) {
  p_dot_ora_tp53 <- dotplot(ora_tp53, showCategory=20)
  titled_print(p_dot_ora_tp53, "ORA dotplot (TP53)")
  ggsave(paste0(plot_dir, "pathway/pathway_tp53_dot", param_suffix_tp53, ".png"), plot=p_dot_ora_tp53, width=14, height=9)
  p_cnet_ora_tp53 <- cnetplot(ora_tp53, showCategory=10, color.params=list(edge=TRUE), node_label="category")
  titled_print(p_cnet_ora_tp53, "ORA cnetplot (TP53)")
  ggsave(paste0(plot_dir, "pathway/pathway_tp53_cnet", param_suffix_tp53, ".png"), plot=p_cnet_ora_tp53, width=14, height=9)
  ora_tp53_pairwise <- pairwise_termsim(ora_tp53)
  p_emap_ora_tp53 <- emapplot(ora_tp53_pairwise, pie="count", showCategory=20, group_category=TRUE, group_legend=TRUE)
  titled_print(p_emap_ora_tp53, "ORA emapplot (TP53)")
  ggsave(paste0(plot_dir, "pathway/pathway_tp53_emap", param_suffix_tp53, ".png"), plot=p_emap_ora_tp53, width=14, height=9)

  ora_tp53_simple <- simplify(ora_tp53, cutoff=0.5, by="p.adjust", select_fun=min)
  write.csv(as.data.frame(ora_tp53_simple), paste0(r_dir_files, "pathway_tp53_simple", param_suffix_tp53, ".csv"), row.names=FALSE)
  if (!is.null(ora_tp53_simple) && nrow(as.data.frame(ora_tp53_simple)) > 0) {
    p_dot_ora_tp53_simple <- dotplot(ora_tp53_simple, showCategory=20)
    titled_print(p_dot_ora_tp53_simple, "ORA dotplot (TP53, simplified)")
    ggsave(paste0(plot_dir, "pathway/pathway_tp53_dot_simplified", param_suffix_tp53, ".png"), plot=p_dot_ora_tp53_simple, width=14, height=9)
    p_cnet_ora_tp53_simple <- cnetplot(ora_tp53_simple, showCategory=10, color.params=list(edge=TRUE), node_label="category")
    titled_print(p_cnet_ora_tp53_simple, "ORA cnetplot (TP53, simplified)")
    ggsave(paste0(plot_dir, "pathway/pathway_tp53_cnet_simplified", param_suffix_tp53, ".png"), plot=p_cnet_ora_tp53_simple, width=14, height=9)
    ora_tp53_simple_pairwise <- pairwise_termsim(ora_tp53_simple)
    p_emap_ora_tp53_simple <- emapplot(ora_tp53_simple_pairwise, pie="count", showCategory=20, group_category=TRUE, group_legend=TRUE)
    titled_print(p_emap_ora_tp53_simple, "ORA emapplot (TP53, simplified)")
    ggsave(paste0(plot_dir, "pathway/pathway_tp53_emap_simplified", param_suffix_tp53, ".png"), plot=p_emap_ora_tp53_simple, width=14, height=9)
  } else {
    warning("No significant enrichment found in ora_tp53_simple; skipping simplified pathway plots.")
  }
} else {
  warning("No significant enrichment found in ora_tp53; skipping pathway plots.")
}



#### PATHWAY PED CANCER VS HOST SEQ ####
write_output(quote(NULL), "Pathway Analysis - Cancer vs HostSeq Cohort")
# Parameter values: no min_samples, no p_gene filtering (all genes), pvalueCutoff=0.05, qvalueCutoff=0.1
param_suffix_cancer <- "_min0_pgene1_ppathway0.05_qpathway0.1"

ora_cancer <- perform_ora_cancer(te_kics_hostseq)
write.csv(as.data.frame(ora_cancer), paste0(r_dir_files, "pathway_cancer", param_suffix_cancer, ".csv"), row.names=FALSE)

if (!is.null(ora_cancer) && nrow(as.data.frame(ora_cancer)) > 0) {
  p_dot_ora_cancer <- dotplot(ora_cancer, showCategory=20)
  titled_print(p_dot_ora_cancer, "ORA dotplot (Cancer vs HostSeq)")
  ggsave(paste0(plot_dir, "pathway/pathway_cancer_dot", param_suffix_cancer, ".png"), plot=p_dot_ora_cancer, width=14, height=9)
  p_cnet_ora_cancer <- cnetplot(ora_cancer, showCategory=10, color.params=list(edge=TRUE), node_label="category")
  titled_print(p_cnet_ora_cancer, "ORA cnetplot (Cancer vs HostSeq)")
  ggsave(paste0(plot_dir, "pathway/pathway_cancer_cnet", param_suffix_cancer, ".png"), plot=p_cnet_ora_cancer, width=14, height=9)
  ora_cancer_pairwise <- pairwise_termsim(ora_cancer)
  p_emap_ora_cancer <- emapplot(ora_cancer_pairwise, pie="count", showCategory=20, group_category=TRUE, group_legend=TRUE)
  titled_print(p_emap_ora_cancer, "ORA emapplot (Cancer vs HostSeq)")
  ggsave(paste0(plot_dir, "pathway/pathway_cancer_emap", param_suffix_cancer, ".png"), plot=p_emap_ora_cancer, width=14, height=9)

  # Simplified pathway plots
  ora_cancer_simple <- simplify(ora_cancer, cutoff=0.5, by="p.adjust", select_fun=min)
  write.csv(as.data.frame(ora_cancer_simple), paste0(r_dir_files, "pathway_cancer_simple", param_suffix_cancer, ".csv"), row.names=FALSE)
  if (!is.null(ora_cancer_simple) && nrow(as.data.frame(ora_cancer_simple)) > 0) {
    p_dot_ora_cancer_simple <- dotplot(ora_cancer_simple, showCategory=20)
    titled_print(p_dot_ora_cancer_simple, "ORA dotplot (Cancer vs HostSeq, simplified)")
    ggsave(paste0(plot_dir, "pathway/pathway_cancer_dot_simplified", param_suffix_cancer, ".png"), plot=p_dot_ora_cancer_simple, width=14, height=9)
    p_cnet_ora_cancer_simple <- cnetplot(ora_cancer_simple, showCategory=10, color.params=list(edge=TRUE), node_label="category")
    titled_print(p_cnet_ora_cancer_simple, "ORA cnetplot (Cancer vs HostSeq, simplified)")
    ggsave(paste0(plot_dir, "pathway/pathway_cancer_cnet_simplified", param_suffix_cancer, ".png"), plot=p_cnet_ora_cancer_simple, width=14, height=9)
    ora_cancer_simple_pairwise <- pairwise_termsim(ora_cancer_simple)
    p_emap_ora_cancer_simple <- emapplot(ora_cancer_simple_pairwise, pie="count", showCategory=20, group_category=TRUE, group_legend=TRUE)
    titled_print(p_emap_ora_cancer_simple, "ORA emapplot (Cancer vs HostSeq, simplified)")
    ggsave(paste0(plot_dir, "pathway/pathway_cancer_emap_simplified", param_suffix_cancer, ".png"), plot=p_emap_ora_cancer_simple, width=14, height=9)
  }
} else {
  warning("No significant enrichment found in ora_cancer; skipping pathway plots.")
}



cat("✓ Script completed successfully\n")
