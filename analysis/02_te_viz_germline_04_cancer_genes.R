#!/usr/bin/env Rscript

# Germline TE Visualization - Cancer Genes
# Cancer genes analysis and top mutated genes

# Source common setup and load data
source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/00_viz_common_setup.R")
source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/00_viz_load_data_germline.R")

cat("Running 02_te_viz_germline_04_cancer_genes.R...\n")

#### CANCER GENES ####
cat("Processing cancer genes...\n")
genes <- cpg$V1

# overall cancer genes effected
te_cancer_hostseq_split <- te_all_split %>%
  filter(cohort == "HostSeq") %>%
  bind_rows(te_aff_split)
te_cancer_hostseq_split_cancergenes  <- te_cancer_hostseq_split %>% filter(Gene_name %in% genes)
te_cancer_hostseq_split_cancergenes_processed <- process_all_combinations(te_cancer_hostseq_split_cancergenes)
te_cancer_hostseq_split_cancergenes_processed <- as.data.frame(add_nohit_samples(te_cancer_hostseq_split_cancergenes_processed, nohits))
write_output(quote(summary(te_cancer_hostseq_split_cancergenes_processed$total)), "Summary: Number of cancer genes affected per sample")

# merge with clinical and metrics
tmp <- merge_dfs(te_cancer_hostseq_split_cancergenes_processed, clinical, include_all_x = FALSE)
te_cancer_hostseq_split_cancergenes_processed <- merge_dfs(tmp, metrics, include_all_x = FALSE)

# sig test
#p_cancergenes_kruskal <- write_output(
#  quote(plot_count_kruskal(te_cancer_hostseq_split_cancergenes_processed, type=NA, group="cancer_cohort", x_lab=x_tp53, y_lab="Number of cancer genes effected", chr=NA, log_scale=TRUE)),
#  "Kruskal test: Number of cancer genes affected by cohort"
#)
#titled_print(p_cancergenes_kruskal, "Kruskal: Number of cancer genes affected by cohort")
#ggsave(paste0(plot_dir, "te_cancergenes_kruskal.png"), plot=p_cancergenes_kruskal, width=9, height=5)

#p_cancergenes_lm <- write_output(
#  quote(plot_count_lm(te_cancer_hostseq_split_cancergenes_processed, type=NA, covariates = "tumor_type", group="cancer_cohort", y_lab="Number of cancer genes effected", residuals=FALSE, log_scale=FALSE, x_lab=x_cancer_tp53, min_samples=5, chr=NA)),
#  "Linear model: Number of cancer genes affected by cohort"
#)
#titled_print(p_cancergenes_lm, "Linear model: Number of cancer genes affected by cohort")
#ggsave(paste0(plot_dir, "te_cancergenes_lm.png"), plot=p_cancergenes_lm, width=9, height=5)


# top genes (KICS)
te_kics_split <- add_gene_size_todf(te_kics_split, gene_size)
geneList_kics <- calculate_gene_mut_sample_frequency(te_kics_split)
geneList_kics_filtered <- geneList_kics[geneList_kics$Gene_name %in% genes, ]
write.csv(geneList_kics_filtered, paste0(r_dir_files, "te_cancer_genes_kics.csv"), row.names=FALSE)

p_top_cancer_genes <- plot_top_genes(geneList_kics_filtered, column="freq_samples_effected", label_column=NULL, x_lab="Proportion of samples with insertion in gene", top_n=20)
titled_print(p_top_cancer_genes, "Top cancer genes affected (KICS)")
ggsave(paste0(plot_dir, "cancer_genes/te_top_cancer_genes_affected_kics.png"), plot=p_top_cancer_genes, width=9, height=5)

p_top_cancer_genes_norm <- plot_top_genes(geneList_kics_filtered, column="freq_samples_effected_normalized", label_column="num_samples_effected", x_lab="Proportion of samples with insertion in gene normalized for gene size", top_n=20)
titled_print(p_top_cancer_genes_norm, "Top cancer genes affected (normalized, KICS)")
ggsave(paste0(plot_dir, "cancer_genes/te_top_cancer_genes_affected_normalized_kics.png"), plot=p_top_cancer_genes_norm, width=9, height=5)

# top genes (LFS)
te_lfs_split <- add_gene_size_todf(te_lfs_split, gene_size)
geneList_lfs <- calculate_gene_mut_sample_frequency(te_lfs_split)
geneList_lfs_filtered <- geneList_lfs[geneList_lfs$Gene_name %in% genes, ]
write.csv(geneList_lfs_filtered, paste0(r_dir_files, "te_cancer_genes_lfs.csv"), row.names=FALSE)

p_top_cancer_genes_lfs <- plot_top_genes(geneList_lfs_filtered, column="freq_samples_effected", label_column=NULL, x_lab="Proportion of samples with insertion in gene", top_n=20)
titled_print(p_top_cancer_genes_lfs, "Top cancer genes affected (LFS)")
ggsave(paste0(plot_dir, "cancer_genes/te_top_cancer_genes_affected_lfs.png"), plot=p_top_cancer_genes_lfs, width=9, height=5)

p_top_cancer_genes_norm_lfs <- plot_top_genes(geneList_lfs_filtered, column="freq_samples_effected_normalized", label_column="num_samples_effected", x_lab="Proportion of samples with insertion in gene normalized for gene size", top_n=20)
titled_print(p_top_cancer_genes_norm_lfs, "Top cancer genes affected (normalized, LFS)")
ggsave(paste0(plot_dir, "cancer_genes/te_top_cancer_genes_affected_normalized_lfs.png"), plot=p_top_cancer_genes_norm_lfs, width=9, height=5)

# genes present that meet criteria
location <- c("CDS", "UTR", "3'UTR", "5'UTR", "5'UTR-CDS","CDS-3'UTR")
genes_aff <- identify_te_genes(te_aff_split_genes, gene_vector=genes, location="exon", location2=NULL)
genes_kics <- identify_te_genes(te_kics_split, gene_vector=genes, location=NULL, location2=NULL)
genes_lfs <- identify_te_genes(te_lfs_split, gene_vector=genes, location=NULL, location2=NULL)
write.csv(genes_aff, paste0(r_dir_files, "te_cancer_genes_aff_exon.csv"), row.names=FALSE)

# plot top cancer genes in germline
#titled_print(plot_top_cancer_genes_germline(genes_aff, top_n=15), "Top cancer genes in germline (affected)")
#ggsave(paste0(plot_dir, "cancer_genes/te_top_cancer_genes_germline_aff.png"), width=9, height=5)
#
#titled_print(plot_top_cancer_genes_germline(genes_kics, top_n=15), "Top cancer genes in germline (KICS)")
#ggsave(paste0(plot_dir, "cancer_genes/te_top_cancer_genes_germline_kics.png"), width=9, height=5)
#
#titled_print(plot_top_cancer_genes_germline(genes_lfs, top_n=15), "Top cancer genes in germline (LFS)")
#ggsave(paste0(plot_dir, "cancer_genes/te_top_cancer_genes_germline_lfs.png"), width=9, height=5)
#
# wilcox test removed - use specific TE tests by gene instead

# wilcox test b/w lfs aff and unaff
#write_output(quote(count_location_wilcox(te_lfs_split_rare, gene=genes, filter_element=NA, group="Cancer", location=location)), "Wilcoxon: Cancer genes in LFS affected vs unaffected")



#### TOP MUTATED GENES ####
# KICS
cat("Calculating top mutated genes in KICS\n")
geneList_kics <- calculate_gene_mut_sample_frequency(te_kics_split)
write.csv(geneList_kics, paste0(r_dir_files, "te_top_genes_kics.csv"), row.names=FALSE)

# Check if required columns exist and are not empty before plotting
if ("freq_samples_effected" %in% colnames(geneList_kics) && 
    "freq_mutations" %in% colnames(geneList_kics) &&
    nrow(geneList_kics[!is.na(geneList_kics$freq_samples_effected) & !is.na(geneList_kics$freq_mutations), ]) > 0) {
  p_top_genes_kics <- plot_top_genes(geneList_kics, column="freq_samples_effected", label_column="freq_mutations", x_lab="Proportion of samples with insertion in gene normalized by gene size", top_n=20)
  titled_print(p_top_genes_kics, "Top mutated genes (KICS)")
  ggsave(paste0(plot_dir, "cancer_genes/te_top_genes_affected_normalized_kics.png"), plot=p_top_genes_kics, width=9, height=5)
} else {
  warning("Required columns for plotting top genes in KICS are missing or empty.")
}

# number of genes affecting > 50%, 75%, 90% of samples
thresh_counts_kics <- count_genes_by_threshold(geneList_kics, column_name = "freq_samples_effected", thresholds = c(0.5, 0.75, 0.9))
write_output(quote(print(thresh_counts_kics)), "Number of genes affecting >50%, 75%, 90% of KICS samples")

# LFS
cat("Calculating top mutated genes in LFS\n")
# Note: te_lfs_split already has gene_size added on line 1369
geneList_lfs <- calculate_gene_mut_sample_frequency(te_lfs_split, nsample_thresh = 0, filter_exon = FALSE)
write.csv(geneList_lfs, paste0(r_dir_files, "te_top_genes_lfs.csv"), row.names=FALSE)

p_top_genes_lfs <- plot_top_genes(geneList_lfs, column="freq_samples_effected_normalized", label_column="freq_samples_effected", x_lab="Proportion of samples with insertion in gene normalized by gene size", top_n=20)
titled_print(p_top_genes_lfs, "Top mutated genes (LFS)")
ggsave(paste0(plot_dir, "cancer_genes/te_top_genes_affected_normalized_lfs.png"), plot=p_top_genes_lfs, width=9, height=5)




cat("✓ Script completed successfully\n")
