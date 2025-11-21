#!/usr/bin/env Rscript

# Tumour TE Visualization - Pathway Analysis
# Pathway analysis for TP53-specific TEs, KICS cohort, TP53 status

# Source common setup and load data
source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/00_viz_common_setup.R")
source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/00_viz_load_data_tumour.R")

cat("Running 02_te_viz_tumour_06_pathway.R...\n")

# Define TE types for analysis
types <- c(NA, "LINE1", "ALU", "SVA")
x_tp53 <- expression("Somatic " * italic("TP53") * " status")

#### PARAMETER SWEEP CONFIGURATION ####

# Define parameter grid for pathway analysis
param_grid <- expand.grid(
  min_samples = c(3),
  p_gene = c(0.05, 0.1),
  p_pathway = c(0.05, 0.1),
  q_pathway = c(0.05, 0.1),
  stringsAsFactors = FALSE
)
cat("\n*** PATHWAY ANALYSIS PARAMETER SWEEP ***\n")
cat("Testing", nrow(param_grid), "parameter combinations\n\n")

#### PATHWAY ANALYSIS - TEs SPECIFIC TO TP53 STATUS ####
write_output(quote(NULL), "Pathway Analysis - TEs Specific to TP53 STATUS")

# Load specific TE files (will be filtered by parameter sweep)
specific_te_files_tp53_t <- list()
for (ms in unique(param_grid$min_samples)) {
  file_path <- paste0(r_dir_files, "specific_tes_tumour_min", ms, "_full_results.csv")
  if (file.exists(file_path)) {
    specific_te_files_tp53_t[[paste0("tp53_min", ms)]] <- read.csv(file_path, stringsAsFactors = FALSE)
    cat("Loaded TP53-specific TEs (min_samples=", ms, "):", nrow(specific_te_files_tp53_t[[paste0("tp53_min", ms)]]), "TEs\n")
  }
}

# Run parameter sweep for TP53-specific TEs
if (length(specific_te_files_tp53_t) > 0) {
  tp53_sweep_results_t <- run_pathway_parameter_sweep(
    data_list = specific_te_files_tp53_t,
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

# LEGACY CODE COMMENTED OUT - Now using parameter sweep system above
# Old single-run code retained for reference but not executed
#
# min_samples <- 5
# p_thresh <- 0.05
# specific_te_file_t <- paste0(r_dir_files, "specific_tes_tumour_min", min_samples, "_full_results.csv")
# if (file.exists(specific_te_file_t)) {
#   # Visualizations
#   p_dot_ora_tp53_specific_t <- dotplot(ora_tp53_specific_t, showCategory=20)
#   titled_print(p_dot_ora_tp53_specific_t, paste0("ORA dotplot (TP53-specific TEs, min", min_samples, ", FDR<", p_thresh, ")"))
#   ggsave(paste0(plot_dir, "pathway/pathway_tp53_specific_min", min_samples, "_fdr", p_thresh, "_dot.png"), plot=p_dot_ora_tp53_specific_t, width=14, height=9)
#
#   p_cnet_ora_tp53_specific_t <- cnetplot(ora_tp53_specific_t, showCategory=10, colorEdge=TRUE, node_label="category")
#   titled_print(p_cnet_ora_tp53_specific_t, paste0("ORA cnetplot (TP53-specific TEs, min", min_samples, ", FDR<", p_thresh, ")"))
#   ggsave(paste0(plot_dir, "pathway/pathway_tp53_specific_min", min_samples, "_fdr", p_thresh, "_cnet.png"), plot=p_cnet_ora_tp53_specific_t, width=14, height=9)
#
#   ora_tp53_specific_pairwise_t <- pairwise_termsim(ora_tp53_specific_t)
#   p_emap_ora_tp53_specific_t <- emapplot(ora_tp53_specific_pairwise_t, pie="count", showCategory=20, group_category=TRUE, group_legend=TRUE)
#   titled_print(p_emap_ora_tp53_specific_t, paste0("ORA emapplot (TP53-specific TEs, min", min_samples, ", FDR<", p_thresh, ")"))
#   ggsave(paste0(plot_dir, "pathway/pathway_tp53_specific_min", min_samples, "_fdr", p_thresh, ".png"), plot=p_emap_ora_tp53_specific_t, width=14, height=9)
# } else {
#   cat("No significant enrichment found for TP53-specific TEs\n")
# }
# } else {
#   cat("No significant TP53-specific TEs found\n")
# }
# } else {
#   cat("Specific TE results file not found:", specific_te_file_t, "\n")
# }



tp53_fitness <- fread("/Users/briannelaverty/Documents/R_Malkin/TE/data/tp53_variants_predictions.tsv", sep="\t", header=TRUE)
cat("Adding fitness values to affected and LFS data...\n")
te_aff_t_tp53 <- prep_p53_fitness(te_aff_t, tp53_fitness)
te_lfs_t_tp53 <- prep_p53_fitness(te_lfs_t, tp53_fitness)
cat("Scatter plot of TE count vs p53 fitness (affected)...\n")
for (i in seq_along(types)) {
  plot_title <- paste0("TE count vs p53 fitness (", ifelse(is.na(types[i]), "all types", types[i]), ")")
  dependent_var <- ifelse(is.na(types[i]), "total", types[i])
  tryCatch({
    p <- scatter_template_lm_one_variable(te_aff_t_tp53, "p53_fitness", dependent_var, c(0.25, 1), c(0, 1000), "p53 fitness", "Repeat count", shape_column = NULL, colour_column = NULL)
    titled_print(p, plot_title)
    ggsave(paste0(plot_dir, "counts_clinical_lfs/te_count_lfs_fitness_", ifelse(is.na(types[i]), "all", types[i]), ".png"), plot = p, width = 9, height = 5)
  }, error = function(e) {
    cat("Skipping fitness plot for type", ifelse(is.na(types[i]), "all", types[i]), "due to error:", e$message, "\n")
  })
}

cat("\n===== P53 FITNESS BY INHERITANCE =====\n")
p_p53_fitness_inheritance <- plot_box_kruskal(te_lfs_t_tp53, "p53_fitness", "inheritance", "Inheritance", "p53 fitness")
titled_print(p_p53_fitness_inheritance, "p53 fitness by inheritance (LFS)")
ggsave(paste0(plot_dir, "other/p53_fitness_by_inheritance_lfs.png"), plot = p_p53_fitness_inheritance, width = 5, height = 5)

cat("\n===== CANCER GENES ANALYSIS =====\n")
cpg <- read.csv("/Users/briannelaverty/Documents/R_Malkin/te/data/raw/kics_cpg.txt", header=FALSE, sep="\t")
genes <- cpg$V1

# remove samples with too many TEs
#to_remove <- te_all_t$sample[te_all_t$total > 100000]
#te_kics_split_t <- te_kics_split_t %>% filter(!sample %in% to_remove)
#te_aff_split_t <- te_aff_split_t %>% filter(!sample %in% to_remove)

# Cancer genes affected
load(paste0(r_dir, "nohits_final_te_count_t_te_aff_selected_t", ".RData"))
te_aff_split_cancergenes_t <- te_aff_split_t %>% filter(Gene_name %in% genes)
te_aff_split_cancergenes_processed_t <- process_all_combinations(te_aff_split_cancergenes_t)
te_aff_split_cancergenes_processed_t <- as.data.frame(add_nohit_samples(te_aff_split_cancergenes_processed_t, nohits))
write_output(quote(summary(te_aff_split_cancergenes_processed_t$total)), "Summary of cancer genes (affected)")
te_aff_split_cancergenes_processed_t <- merge_dfs(te_aff_split_cancergenes_processed_t, clinical, include_all_x = FALSE) # merge with clinical data

cat("Plotting cancer genes affected by TP53 status (Wilcoxon test)...\n")
for (i in seq_along(types)) {
  plot_title <- paste0(ifelse(is.na(types[i]), "All types", types[i]), " cancer genes affected by TP53 status")
  write_output(
    quote(plot_count_wilcox(te_aff_split_cancergenes_processed_t, type=types[i], group="TP53_status", y_lab="Cancer genes affected", chr=NA, x_lab=x_tp53, log_scale=FALSE)),
    paste0("Cancer genes affected by TP53 status (type=", ifelse(is.na(types[i]), "all", types[i]), ")")
  )
  p <- plot_count_wilcox(te_aff_split_cancergenes_processed_t, type=types[i], group="TP53_status", y_lab="Cancer genes affected", chr=NA, x_lab=x_tp53, log_scale=FALSE)
  titled_print(p, plot_title)
  ggsave(paste0(plot_dir, "cancer_genes/te_cancergenes_", ifelse(is.na(types[i]), "all", types[i]), "_tp53.png"), plot = p, width = 5, height = 5)
}
# ===== TOP MUTATED CANCER GENES =====
# Add gene size to KICS split data
tryCatch({
  te_kics_split_t <- add_gene_size_todf(te_kics_split_t, gene_size)
}, error = function(e) {
  cat("Warning: Could not add gene size data:", e$message, "\n")
})
# Add a dummy gene_size column if it doesn't exist
if (!"gene_size" %in% colnames(te_kics_split_t)) {
  te_kics_split_t$gene_size <- 1
}

geneList_kics_t <- calculate_gene_mut_sample_frequency(te_kics_split_t, filter_exon = FALSE)
geneList_kics_filtered_t <- geneList_kics_t[geneList_kics_t$Gene_name %in% genes, ]

# Plot top mutated cancer genes (KICS)
write_output(quote(plot_top_genes(geneList_kics_filtered_t, column="freq_samples_effected", label_column=NULL, x_lab="Proportion of samples with insertion in gene", top_n=20)), "Top mutated cancer genes (KICS, raw frequency)")
p_cancer1 <- plot_top_genes(geneList_kics_filtered_t, column="freq_samples_effected", label_column=NULL, x_lab="Proportion of samples with insertion in gene", top_n=20)
titled_print(p_cancer1, "Top mutated cancer genes (KICS, raw frequency)")
ggsave(paste0(plot_dir, "cancer_genes/te_top_cancer_genes_affected_kics.png"), plot = p_cancer1, width = 9, height = 5)
write_output(quote(plot_top_genes(geneList_kics_filtered_t, column="freq_samples_effected_normalized", label_column="num_samples_effected", x_lab="Proportion of samples with insertion in gene normalized for gene size", top_n=20)), "Top mutated cancer genes (KICS, normalized for gene size)")
p_cancer2 <- plot_top_genes(geneList_kics_filtered_t, column="freq_samples_effected_normalized", label_column="num_samples_effected", x_lab="Proportion of samples with insertion in gene normalized for gene size", top_n=20)
titled_print(p_cancer2, "Top mutated cancer genes (KICS, normalized for gene size)")
ggsave(paste0(plot_dir, "cancer_genes/te_top_cancer_genes_affected_normalized_kics.png"), plot = p_cancer2, width = 9, height = 5)

# Stacked plot for top genes (KICS)
write_output(
  quote(plot_gene_effects(te_kics_split_t, min_sample_tt = 5, remove_other = FALSE, top_n_genes = 15, SV_type = NULL)),
  "Top Affected Genes (KICS)"
)

# Example: Patients with RAF1 affected
#raf1_patients <- te_kics_split_t %>%
#  filter(Gene_name == "RAF1") %>%
#  group_by(sample, tumor_type, TP53_status) %>%
#  summarise(count = n(), .groups = 'drop') %>%
#  as.data.frame()
#write_output(quote(raf1_patients), "RAF1-affected patients (KICS)")
#
#raf1_aff_samples <- te_aff_split_genes_t %>%
#  filter(Gene_name == "RAF1") %>%
#  select(sample, Gene_name, Location, Location2, GnomAD_pLI, ExAC_pLI)
#write_output(quote(raf1_aff_samples), "RAF1-affected samples (affected)")

# Genes present that meet criteria
location <- c("CDS", "UTR", "3'UTR", "5'UTR", "5'UTR-CDS", "CDS-3'UTR")
genes_aff_t <- identify_te_genes(te_aff_split_genes_t, gene_vector=genes, location="exon", location2=NULL)
genes_kics_t <- identify_te_genes(te_kics_split_t, gene_vector=genes, location=NULL, location2=NULL)
genes_lfs_t <- identify_te_genes(te_lfs_split_t, gene_vector=genes, location="exon", location2=NULL)

write_output(quote(table(genes_aff_t$sample, genes_aff_t$Gene_name)), "samples and genes with exon insertion")
#write_output(quote(table(genes_aff_t$sample, genes_aff_t$TP53_status)), "Table: affected samples x TP53_status")
#write_output(quote(table(genes_aff_t$sample, genes_aff_t$tumor_type)), "Table: affected samples x tumor_type")
#write_output(quote(table(genes_aff_t$TP53_status, genes_aff_t$tumor_type)), "Table: TP53_status x tumor_type")
#write_output(quote(table(genes_aff_t$Gene_name, genes_aff_t$tumor_type)), "Table: gene x tumor_type")
#write_output(quote(table(genes_aff_t$Gene_name, genes_aff_t$TP53_status)), "Table: gene x TP53_status")
#write_output(quote(length(table(genes_kics_t$Gene_name))), "Number of unique cancer genes in KICS")
#write_output(quote(table(genes_lfs_t$Gene_name)), "Table: gene x LFS")

# Wilcox test between affected and controls for cancer genes
write_output(quote(count_location_wilcox(te_aff_split_t, gene=genes, filter_element=NA, group="TP53_status", location=NULL)), "Wilcox test: affected vs controls (cancer genes)")

# Wilcox test between LFS affected and unaffected for cancer genes
#write_output(quote(count_location_wilcox(te_lfs_split_t, gene=genes, filter_element=NA, group="Cancer", location=location)), "Wilcox test: LFS affected vs unaffected (cancer genes)")

cat("\n===== TOP MUTATED GENES =====\n")
# KICS
tryCatch({
  te_kics_split_t_temp <- add_gene_size_todf(te_kics_split_t, gene_size)
  # Only update if the function succeeded
  if ("gene_size" %in% colnames(te_kics_split_t_temp)) {
    te_kics_split_t <- te_kics_split_t_temp
  }
}, error = function(e) {
  cat("Warning: Could not add gene size data to KICS:", e$message, "\n")
})
# Ensure gene_size column exists
if (!"gene_size" %in% colnames(te_kics_split_t)) {
  cat("Adding dummy gene_size column to KICS data\n")
  te_kics_split_t$gene_size <- 1
}
tryCatch({
  geneList_kics_t <- calculate_gene_mut_sample_frequency(te_kics_split_t)
  p1 <- plot_top_genes(geneList_kics_t, column="freq_samples_effected", label_column=NULL, x_lab="Proportion of samples with insertion in gene", top_n=20)
  titled_print(p1, "Top mutated genes (KICS, raw frequency)")
  ggsave(paste0(plot_dir, "cancer_genes/te_top_genes_affected_kics.png"), plot = p1, width = 9, height = 5)
  p2 <- plot_top_genes(geneList_kics_t, column="freq_samples_effected_normalized", label_column="num_samples_effected", x_lab="Proportion of samples with insertion in gene normalized by gene size", top_n=20)
  titled_print(p2, "Top mutated genes (KICS, normalized for gene size)")
  ggsave(paste0(plot_dir, "cancer_genes/te_top_genes_affected_normalized_kics.png"), plot = p2, width = 9, height = 5)
}, error = function(e) {
  cat("Warning: Could not calculate gene frequencies for KICS:", e$message, "\n")
})

# number of genes affecting >50%, 75%, 90% of samples (KICS)
write_output(quote(count_genes_by_threshold(geneList_kics_t, column_name = "freq_samples_effected", thresholds = c(0.5, 0.75, 0.9))), "Number of genes affecting >50%, 75%, 90% of samples (KICS)")

# LFS
tryCatch({
  te_lfs_split_t_temp <- add_gene_size_todf(te_lfs_split_t, gene_size)
  # Only update if the function succeeded
  if ("gene_size" %in% colnames(te_lfs_split_t_temp)) {
    te_lfs_split_t <- te_lfs_split_t_temp
  }
}, error = function(e) {
  cat("Warning: Could not add gene size data to LFS:", e$message, "\n")
})
# Ensure gene_size column exists
if (!"gene_size" %in% colnames(te_lfs_split_t)) {
  cat("Adding dummy gene_size column to LFS data\n")
  te_lfs_split_t$gene_size <- 1
}
tryCatch({
  geneList_lfs_t <- calculate_gene_mut_sample_frequency(te_lfs_split_t, nsample_thresh = 0, filter_exon = FALSE)
  p3 <- plot_top_genes(geneList_lfs_t, column="freq_samples_effected", label_column=NULL, x_lab="Proportion of samples with insertion in gene", top_n=20)
  titled_print(p3, "Top mutated genes (LFS, raw frequency)")
  ggsave(paste0(plot_dir, "cancer_genes/te_top_genes_affected_lfs.png"), plot = p3, width = 9, height = 5)
  p4 <- plot_top_genes(geneList_lfs_t, column="freq_samples_effected_normalized", label_column="freq_samples_effected", x_lab="Proportion of samples with insertion in gene normalized by gene size", top_n=20)
  titled_print(p4, "Top mutated genes (LFS, normalized for gene size)")
  ggsave(paste0(plot_dir, "cancer_genes/te_top_genes_affected_normalized_lfs.png"), plot = p4, width = 9, height = 5)
}, error = function(e) {
  cat("Warning: Could not calculate gene frequencies for LFS:", e$message, "\n")
})

# Control
te_control_split_t <- te_aff_split_t %>% filter(TP53_status == "Control")
if (nrow(te_control_split_t) > 0) {
  # Ensure gene_size column exists
  if (!"gene_size" %in% colnames(te_control_split_t)) {
    cat("Adding dummy gene_size column to Control data\n")
    te_control_split_t$gene_size <- 1
  }
  tryCatch({
    geneList_control_t <- calculate_gene_mut_sample_frequency(te_control_split_t, nsample_thresh = 0, filter_exon = FALSE)
    p5 <- plot_top_genes(geneList_control_t, top_n=20)
    titled_print(p5, "Top mutated genes (Control, raw frequency)")
    ggsave(paste0(plot_dir, "cancer_genes/te_top_genes_affected_control.png"), plot = p5, width = 9, height = 5)
  }, error = function(e) {
    cat("Warning: Could not calculate gene frequencies for Control:", e$message, "\n")
  })
} else {
  cat("No control samples found, skipping control analysis\n")
}

# Specific gene heatmap of insertions (example: RAF1)
#p_raf1_heatmap <- plot_te_insertion_complex_heatmap(te_kics_split_t, gene_name = "RAF1", bin_size = 1000)
#titled_print(p_raf1_heatmap, "TE insertion heatmap for RAF1 (KICS)")

cat("\n===== PATHWAY PEDIATRIC (ORA) =====\n")
# Parameter values: no min_samples, no p_gene filtering (all genes), pvalueCutoff=0.05, qvalueCutoff=0.1
param_suffix_kics_t <- "_min0_pgene1_ppathway0.05_qpathway0.1"

# over representation analysis
tryCatch({
  ora_kics_t <- perform_ora(te_kics_split_t, nsample_thresh = 0, filter_exon = FALSE)
  write.csv(as.data.frame(ora_kics_t), paste0(r_dir_files, "pathway_kics", param_suffix_kics_t, ".csv"), row.names=FALSE)

  # viz
  tryCatch({
    p_ora_bar_kics <- barplot(ora_kics_t, showCategory=20)
    titled_print(p_ora_bar_kics, "ORA Barplot (KICS)")
    ggsave(paste0(plot_dir, "pathway/pathway_kics_bar", param_suffix_kics_t, ".png"), plot = p_ora_bar_kics, width = 14, height = 9)
  }, error = function(e) {
    cat("Warning: Could not create ORA barplot:", e$message, "\n")
  })
  
  tryCatch({
    p_ora_dot_kics <- dotplot(ora_kics_t, showCategory=20)
    titled_print(p_ora_dot_kics, "ORA Dotplot (KICS)")
    ggsave(paste0(plot_dir, "pathway/pathway_kics_dot", param_suffix_kics_t, ".png"), plot = p_ora_dot_kics, width = 14, height = 9)
  }, error = function(e) {
    cat("Warning: Could not create ORA dotplot:", e$message, "\n")
  })

  tryCatch({
    p_ora_cnet_kics <- cnetplot(ora_kics_t, showCategory = 10, color.params = list(edge = TRUE), node_label = "category")
    titled_print(p_ora_cnet_kics, "ORA Cnetplot (KICS)")
    ggsave(paste0(plot_dir, "pathway/pathway_kics_cnet", param_suffix_kics_t, ".png"), plot = p_ora_cnet_kics, width = 14, height = 9)
  }, error = function(e) {
    cat("Warning: Could not create ORA cnetplot:", e$message, "\n")
  })

  tryCatch({
    ora_kics_pairwise_t <- pairwise_termsim(ora_kics_t)
    p_ora_emap_kics <- emapplot(ora_kics_pairwise_t, pie="count", showCategory=20, group_category=TRUE, group_legend=TRUE)
    titled_print(p_ora_emap_kics, "ORA Emapplot (KICS)")
    ggsave(paste0(plot_dir, "pathway/pathway_kics_emap", param_suffix_kics_t, ".png"), plot = p_ora_emap_kics, width = 14, height = 9)
  }, error = function(e) {
    cat("Warning: Could not create ORA emapplot:", e$message, "\n")
  })

  # group descriptions by broader function
  cat("Simplifying ORA (KICS)...\n")
  tryCatch({
    ora_kics_simple_t <- simplify(ora_kics_t, cutoff=0.6, by="p.adjust", select_fun=min)
    write.csv(as.data.frame(ora_kics_simple_t), paste0(r_dir_files, "pathway_kics_simple", param_suffix_kics_t, ".csv"), row.names=FALSE)
    p_ora_simple_bar_kics <- barplot(ora_kics_simple_t, showCategory=20)
    titled_print(p_ora_simple_bar_kics, "ORA Barplot (KICS, simplified)")
    ggsave(paste0(plot_dir, "pathway/pathway_kics_bar_simplified", param_suffix_kics_t, ".png"), plot = p_ora_simple_bar_kics, width = 14, height = 9)
    p_ora_simple_dot_kics <- dotplot(ora_kics_simple_t, showCategory=20)
    titled_print(p_ora_simple_dot_kics, "ORA Dotplot (KICS, simplified)")
    ggsave(paste0(plot_dir, "pathway/pathway_kics_dot_simplified", param_suffix_kics_t, ".png"), plot = p_ora_simple_dot_kics, width = 14, height = 9)
    p_ora_simple_cnet_kics <- cnetplot(ora_kics_simple_t, showCategory = 10, color.params = list(edge = TRUE), node_label = "category")
    titled_print(p_ora_simple_cnet_kics, "ORA Cnetplot (KICS, simplified)")
    ggsave(paste0(plot_dir, "pathway/pathway_kics_cnet_simplified", param_suffix_kics_t, ".png"), plot = p_ora_simple_cnet_kics, width = 14, height = 9)
    ora_kics_simple_pairwise_t <- pairwise_termsim(ora_kics_simple_t)
    p_ora_simple_emap_kics <- emapplot(ora_kics_simple_pairwise_t, pie="count", showCategory=20, group_category=TRUE, group_legend=TRUE)
    titled_print(p_ora_simple_emap_kics, "ORA Emapplot (KICS, simplified)")
    ggsave(paste0(plot_dir, "pathway/pathway_kics_emap_simplified", param_suffix_kics_t, ".png"), plot = p_ora_simple_emap_kics, width = 14, height = 9)
  }, error = function(e) {
    cat("Warning: Could not simplify ORA or create simplified plots:", e$message, "\n")
  })

  # common genes in pathway
  cat("Plotting pathway gene counts (KICS)...\n")
  tryCatch({
    pathway_genes_t <- plot_pathway_gene_counts(ora_kics_t, n_descriptions=nrow(ora_kics_t), n_genes=50)
    titled_print(pathway_genes_t, "Pathway gene counts (KICS)")
    ggsave(paste0(plot_dir, "pathway/pathway_kics_gene_counts", param_suffix_kics_t, ".png"), plot = pathway_genes_t, width = 9, height = 5)
  }, error = function(e) {
    cat("Warning: Could not create pathway gene counts plot:", e$message, "\n")
  })
  
}, error = function(e) {
  cat("Warning: Could not perform ORA analysis:", e$message, "\n")
})

# look at correlation b/w top genes mutated and genes in top pathways
tryCatch({
  p_analyze_gene_mut_kics <- analyze_gene_mutations(geneList_kics_t, ora_kics_t, n_descriptions=nrow(ora_kics_t), column = "num_samples_effected", y_lab="Number of samples with insertion in gene")
  titled_print(p_analyze_gene_mut_kics, "Gene mutations vs pathway genes (KICS)")
  ggsave(paste0(plot_dir, "pathway/pathway_kics_genes_vs_numsamples", param_suffix_kics_t, ".png"), plot = p_analyze_gene_mut_kics, width = 9, height = 5)
}, error = function(e) {
  cat("Warning: Could not analyze gene mutations vs pathways:", e$message, "\n")
})

# gene set enrichment analysis
#cat("Performing GSEA (KICS)...\n")
#te_kics_split <- replace_gene_names(te_kics_split)
#te_kics_split <- add_gene_size_todf(te_kics_split, gene_size)
#gsea_kics <- perform_gsea(te_kics_split, column="freq_samples_effected", nsample_thresh = 0, filter_exon = FALSE)
#
## viz
#p_gsea_dot_kics <- dotplot(gsea_kics, showCategory=20)
#titled_print(p_gsea_dot_kics, "GSEA Dotplot (KICS)")
#ggsave(paste0(plot_dir, "te_gsea_dot_kics.png"), plot = p_gsea_dot_kics, width = 14, height = 9)
#p_gsea_go_kics <- goplot(gsea_kics, showCategory = 5)
#titled_print(p_gsea_go_kics, "GSEA Goplot (KICS)")
#ggsave(paste0(plot_dir, "te_gsea_go_kics.png"), plot = p_gsea_go_kics, width = 9, height = 7)
#p_gsea_cnet_kics <- cnetplot(gsea_kics, showCategory = 10, colorEdge = TRUE, node_label = "category")
#titled_print(p_gsea_cnet_kics, "GSEA Cnetplot (KICS)")
#ggsave(paste0(plot_dir, "te_gsea_cnet_kics.png"), plot = p_gsea_cnet_kics, width = 14, height = 9)
#p_gsea_cnet_circ_kics <- cnetplot(gsea_kics, showCategory = 5, node_label = "category", circular = TRUE, colorEdge = TRUE)
#titled_print(p_gsea_cnet_circ_kics, "GSEA Cnetplot Circular (KICS)")
#ggsave(paste0(plot_dir, "te_gsea_cnet_circular_kics.png"), plot = p_gsea_cnet_circ_kics, width = 9, height = 7)

cat("\n===== PATHWAY TP53 =====\n")
# Parameter values: no min_samples, no p_gene filtering (all genes), pvalueCutoff=0.05, qvalueCutoff=0.1
param_suffix_tp53_t <- "_min0_pgene1_ppathway0.05_qpathway0.1"

tryCatch({
  ora_tp53_t <- perform_ora_tp53(te_aff_split_t)
  ora_tp53_genes_t <- perform_ora_tp53(te_aff_split_genes_t)
  write.csv(as.data.frame(ora_tp53_t), paste0(r_dir_files, "pathway_tp53", param_suffix_tp53_t, ".csv"), row.names=FALSE)

  # viz
  tryCatch({
    p_ora_tp53_dot <- dotplot(ora_tp53_genes_t, showCategory = 20)
    titled_print(p_ora_tp53_dot, "ORA Dotplot (TP53)")
    ggsave(paste0(plot_dir, "pathway/pathway_tp53_dot", param_suffix_tp53_t, ".png"), plot = p_ora_tp53_dot, width = 14, height = 9)
  }, error = function(e) {
    cat("Warning: Could not create TP53 ORA dotplot:", e$message, "\n")
  })

  tryCatch({
    p_ora_tp53_cnet <- cnetplot(ora_tp53_t, showCategory = 10, color.params = list(edge = TRUE), node_label = "category")
    titled_print(p_ora_tp53_cnet, "ORA Cnetplot (TP53)")
    ggsave(paste0(plot_dir, "pathway/pathway_tp53_cnet", param_suffix_tp53_t, ".png"), plot = p_ora_tp53_cnet, width = 14, height = 9)
  }, error = function(e) {
    cat("Warning: Could not create TP53 ORA cnetplot:", e$message, "\n")
  })

  tryCatch({
    ora_tp53_pairwise_t <- pairwise_termsim(ora_tp53_genes_t)
    p_ora_tp53_emap <- emapplot(ora_tp53_pairwise_t, pie="count", showCategory = 20, group_category=TRUE, group_legend=TRUE) + scale_fill_manual(values=colours)
    titled_print(p_ora_tp53_emap, "ORA EMapplot (TP53)")
    ggsave(paste0(plot_dir, "pathway/pathway_tp53_emap", param_suffix_tp53_t, ".png"), plot = p_ora_tp53_emap, width = 14, height = 9)
  }, error = function(e) {
    cat("Warning: Could not create TP53 ORA emapplot:", e$message, "\n")
  })
  
  cat("Simplifying ORA (TP53)...\n")
  tryCatch({
    ora_tp53_simple_t <- simplify(ora_tp53_t, cutoff=0.5, by="p.adjust", select_fun=min)
    write.csv(as.data.frame(ora_tp53_simple_t), paste0(r_dir_files, "pathway_tp53_simple", param_suffix_tp53_t, ".csv"), row.names=FALSE)
    p_ora_tp53_simple_dot <- dotplot(ora_tp53_simple_t, showCategory=20)
    titled_print(p_ora_tp53_simple_dot, "ORA Dotplot (TP53, simplified)")
    ggsave(paste0(plot_dir, "pathway/pathway_tp53_dot_simplified", param_suffix_tp53_t, ".png"), plot = p_ora_tp53_simple_dot, width = 14, height = 9)
    p_ora_tp53_simple_cnet <- cnetplot(ora_tp53_simple_t, showCategory = 10, color.params = list(edge = TRUE), node_label = "category")
    titled_print(p_ora_tp53_simple_cnet, "ORA Cnetplot (TP53, simplified)")
    ggsave(paste0(plot_dir, "pathway/pathway_tp53_cnet_simplified", param_suffix_tp53_t, ".png"), plot = p_ora_tp53_simple_cnet, width = 14, height = 9)
    ora_tp53_simple_pairwise_t <- pairwise_termsim(ora_tp53_simple_t)
    p_ora_tp53_simple_emap <- emapplot(ora_tp53_simple_pairwise_t, pie="count", showCategory=20, group_category=TRUE, group_legend=TRUE)
    titled_print(p_ora_tp53_simple_emap, "ORA Emapplot (TP53, simplified)")
    ggsave(paste0(plot_dir, "pathway/pathway_tp53_emap_simplified", param_suffix_tp53_t, ".png"), plot = p_ora_tp53_simple_emap, width = 14, height = 9)
  }, error = function(e) {
    cat("Warning: Could not simplify TP53 ORA or create simplified plots:", e$message, "\n")
  })
  
  tryCatch({
    cat("Analyzing unique GO terms in TP53 ORA...\n")
    ora_tp53_result_t <- ora_tp53_t@compareClusterResult
    
    # GO terms that were only enriched in KiCS or LFS not both
    ora_tp53_result_unique_t <- ora_tp53_result_t %>%
      group_by(ID) %>%
      filter(n() == 1) %>%
      ungroup()
    
    # split enrichment by group
    ora_tp53_result_control_t <- ora_tp53_result_unique_t %>% filter(Cluster=="Control")
    ora_tp53_result_lfs_t <- ora_tp53_result_unique_t %>% filter(Cluster=="LFS")
    
    if (nrow(ora_tp53_result_control_t) > 0) {
      cat("Plotting pathway gene counts (TP53 Control)...\n")
      pathway_genes_control_t <- plot_pathway_gene_counts(ora_tp53_result_control_t, n_descriptions=nrow(ora_tp53_result_control_t), n_genes=50)
      titled_print(pathway_genes_control_t, "Pathway gene counts (TP53 Control)")
      ggsave(paste0(plot_dir, "pathway/pathway_gene_counts_tp53_control.png"), plot = pathway_genes_control_t, width = 9, height = 5)
    } else {
      cat("No unique control pathways found for TP53 analysis\n")
    }

    if (nrow(ora_tp53_result_lfs_t) > 0) {
      cat("Plotting pathway gene counts (TP53 LFS)...\n")
      pathway_genes_lfs_t <- plot_pathway_gene_counts(ora_tp53_result_lfs_t, n_descriptions=nrow(ora_tp53_result_lfs_t), n_genes=50)
      titled_print(pathway_genes_lfs_t, "Pathway gene counts (TP53 LFS)")
      ggsave(paste0(plot_dir, "pathway/pathway_gene_counts_tp53_lfs.png"), plot = pathway_genes_lfs_t, width = 9, height = 5)
    } else {
      cat("No unique LFS pathways found for TP53 analysis\n")
    }
  }, error = function(e) {
    cat("Warning: Could not analyze unique TP53 pathways:", e$message, "\n")
  })
  
}, error = function(e) {
  cat("Warning: Could not perform TP53 ORA analysis:", e$message, "\n")
})



cat("✓ Script completed successfully\n")
