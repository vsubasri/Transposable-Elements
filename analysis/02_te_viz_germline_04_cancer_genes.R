#!/usr/bin/env Rscript

# Germline TE Visualization - Cancer Genes
# Cancer genes analysis and top mutated genes

# Source common setup and load data
source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/00_viz_common_setup.R")
REQUIRED_DATA <- c("split", "clinical", "genes")
source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/00_viz_load_data_germline.R")

# Initialize module-specific text output
init_module_sink(paste0(plot_dir, "cancer_genes/"), "CANCER_GENES")

cat("Running 02_te_viz_germline_04_cancer_genes.R...\n")

#### CANCER GENES ####
cat("Processing cancer genes...\n")
genes <- cpg$V1

# top genes (KICS)
te_kics_split <- add_gene_size_todf(te_kics_split, gene_size)
geneList_kics <- calculate_gene_mut_sample_frequency(te_kics_split)
geneList_kics_filtered <- geneList_kics[geneList_kics$Gene_name %in% genes, ]
write.csv(geneList_kics_filtered, paste0(plot_dir, "cancer_genes/te_cancer_genes_kics.csv"), row.names=FALSE)

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
write.csv(geneList_lfs_filtered, paste0(plot_dir, "cancer_genes/te_cancer_genes_lfs.csv"), row.names=FALSE)

p_top_cancer_genes_lfs <- plot_top_genes(geneList_lfs_filtered, column="freq_samples_effected", label_column=NULL, x_lab="Proportion of samples with insertion in gene", top_n=20)
titled_print(p_top_cancer_genes_lfs, "Top cancer genes affected (LFS)")
ggsave(paste0(plot_dir, "cancer_genes/te_top_cancer_genes_affected_lfs.png"), plot=p_top_cancer_genes_lfs, width=9, height=5)

p_top_cancer_genes_norm_lfs <- plot_top_genes(geneList_lfs_filtered, column="freq_samples_effected_normalized", label_column="num_samples_effected", x_lab="Proportion of samples with insertion in gene normalized for gene size", top_n=20)
titled_print(p_top_cancer_genes_norm_lfs, "Top cancer genes affected (normalized, LFS)")
ggsave(paste0(plot_dir, "cancer_genes/te_top_cancer_genes_affected_normalized_lfs.png"), plot=p_top_cancer_genes_norm_lfs, width=9, height=5)

# top genes (HostSeq)
te_hostseq_split <- add_gene_size_todf(te_hostseq_split, gene_size)
geneList_hostseq <- calculate_gene_mut_sample_frequency(te_hostseq_split)
geneList_hostseq_filtered <- geneList_hostseq[geneList_hostseq$Gene_name %in% genes, ]
write.csv(geneList_hostseq_filtered, paste0(plot_dir, "cancer_genes/te_cancer_genes_hostseq.csv"), row.names=FALSE)

p_top_cancer_genes_hostseq <- plot_top_genes(geneList_hostseq_filtered, column="freq_samples_effected", label_column=NULL, x_lab="Proportion of samples with insertion in gene", top_n=20)
titled_print(p_top_cancer_genes_hostseq, "Top cancer genes affected (HostSeq)")
ggsave(paste0(plot_dir, "cancer_genes/te_top_cancer_genes_affected_hostseq.png"), plot=p_top_cancer_genes_hostseq, width=9, height=5)

p_top_cancer_genes_norm_hostseq <- plot_top_genes(geneList_hostseq_filtered, column="freq_samples_effected_normalized", label_column="num_samples_effected", x_lab="Proportion of samples with insertion in gene normalized for gene size", top_n=20)
titled_print(p_top_cancer_genes_norm_hostseq, "Top cancer genes affected (normalized, HostSeq)")
ggsave(paste0(plot_dir, "cancer_genes/te_top_cancer_genes_affected_normalized_hostseq.png"), plot=p_top_cancer_genes_norm_hostseq, width=9, height=5)

# genes present that meet criteria
location <- c("CDS", "UTR", "3'UTR", "5'UTR", "5'UTR-CDS","CDS-3'UTR")
genes_aff <- identify_te_genes(te_aff_split_genes, gene_vector=genes, location="exon", location2=NULL)
genes_kics <- identify_te_genes(te_kics_split, gene_vector=genes, location=NULL, location2=NULL)
genes_lfs <- identify_te_genes(te_lfs_split, gene_vector=genes, location=NULL, location2=NULL)
write.csv(genes_aff, paste0(plot_dir, "cancer_genes/te_cancer_genes_aff_exon.csv"), row.names=FALSE)

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
write.csv(geneList_kics, paste0(plot_dir, "cancer_genes/te_top_genes_kics.csv"), row.names=FALSE)

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
write.csv(geneList_lfs, paste0(plot_dir, "cancer_genes/te_top_genes_lfs.csv"), row.names=FALSE)

p_top_genes_lfs <- plot_top_genes(geneList_lfs, column="freq_samples_effected_normalized", label_column="freq_samples_effected", x_lab="Proportion of samples with insertion in gene normalized by gene size", top_n=20)
titled_print(p_top_genes_lfs, "Top mutated genes (LFS)")
ggsave(paste0(plot_dir, "cancer_genes/te_top_genes_affected_normalized_lfs.png"), plot=p_top_genes_lfs, width=9, height=5)

# HostSeq
cat("Calculating top mutated genes in HostSeq\n")
# Note: te_hostseq_split already has gene_size added earlier
geneList_hostseq <- calculate_gene_mut_sample_frequency(te_hostseq_split, nsample_thresh = 0, filter_exon = FALSE)
write.csv(geneList_hostseq, paste0(plot_dir, "cancer_genes/te_top_genes_hostseq.csv"), row.names=FALSE)

p_top_genes_hostseq <- plot_top_genes(geneList_hostseq, column="freq_samples_effected_normalized", label_column="freq_samples_effected", x_lab="Proportion of samples with insertion in gene normalized by gene size", top_n=20)
titled_print(p_top_genes_hostseq, "Top mutated genes (HostSeq)")
ggsave(paste0(plot_dir, "cancer_genes/te_top_genes_affected_normalized_hostseq.png"), plot=p_top_genes_hostseq, width=9, height=5)

# number of genes affecting > 50%, 75%, 90% of samples
thresh_counts_hostseq <- count_genes_by_threshold(geneList_hostseq, column_name = "freq_samples_effected", thresholds = c(0.5, 0.75, 0.9))
write_output(quote(print(thresh_counts_hostseq)), "Number of genes affecting >50%, 75%, 90% of HostSeq samples")

#### GENE SUMMARY TABLES ####
# Summary of all genes affected by TEs with location information

cat("\n===== GENE SUMMARY TABLES =====\n")

# Create genes_affected directory if it doesn't exist
dir.create(paste0(plot_dir, "genes_affected"), showWarnings = FALSE, recursive = TRUE)

# Generate gene summary table for all affected samples
cat("Generating gene summary table...\n")
gene_summary <- te_aff_split %>%
  group_by(Gene_name) %>%
  summarize(
    total_samples_affected = n_distinct(sample),
    total_te_insertions = n(),
    has_exonic = any(grepl("exon", Location, ignore.case = TRUE), na.rm = TRUE),
    has_intronic = any(grepl("intron", Location, ignore.case = TRUE), na.rm = TRUE),
    location_summary = paste(unique(na.omit(Location)), collapse = "; "),
    gene_features = paste(unique(na.omit(Location2)), collapse = "; "),
    .groups = "drop"
  ) %>%
  arrange(desc(total_samples_affected), desc(total_te_insertions)) %>%
  mutate(
    te_location_type = case_when(
      has_exonic & has_intronic ~ "Both exonic & intronic",
      has_exonic ~ "Exonic only",
      has_intronic ~ "Intronic only",
      TRUE ~ "Other"
    )
  ) %>%
  select(Gene_name, total_samples_affected, total_te_insertions, te_location_type, location_summary, gene_features)

# Save to CSV
gene_summary_file <- paste0(plot_dir, "genes_affected/gene_summary_te_affected_germline.csv")
write.csv(gene_summary, gene_summary_file, row.names = FALSE)
cat(sprintf("✓ Gene summary saved to: %s\n", basename(gene_summary_file)))
cat(sprintf("  Total genes analyzed: %d\n", nrow(gene_summary)))
cat(sprintf("  TEs by location type:\n"))
cat(sprintf("    - Exonic only: %d genes (%.1f%%)\n",
            sum(gene_summary$te_location_type == "Exonic only"),
            100 * sum(gene_summary$te_location_type == "Exonic only") / nrow(gene_summary)))
cat(sprintf("    - Intronic only: %d genes (%.1f%%)\n",
            sum(gene_summary$te_location_type == "Intronic only"),
            100 * sum(gene_summary$te_location_type == "Intronic only") / nrow(gene_summary)))
cat(sprintf("    - Both: %d genes (%.1f%%)\n",
            sum(gene_summary$te_location_type == "Both exonic & intronic"),
            100 * sum(gene_summary$te_location_type == "Both exonic & intronic") / nrow(gene_summary)))
cat(sprintf("  Top gene (most samples): %s (%d samples)\n\n",
            gene_summary$Gene_name[1],
            gene_summary$total_samples_affected[1]))

# Generate gene summary for KICS samples
cat("Generating gene summary table for KICS samples...\n")
gene_summary_kics <- te_kics_split %>%
  group_by(Gene_name) %>%
  summarize(
    total_samples_affected = n_distinct(sample),
    total_te_insertions = n(),
    has_exonic = any(grepl("exon", Location, ignore.case = TRUE), na.rm = TRUE),
    has_intronic = any(grepl("intron", Location, ignore.case = TRUE), na.rm = TRUE),
    location_summary = paste(unique(na.omit(Location)), collapse = "; "),
    gene_features = paste(unique(na.omit(Location2)), collapse = "; "),
    .groups = "drop"
  ) %>%
  arrange(desc(total_samples_affected), desc(total_te_insertions)) %>%
  mutate(
    te_location_type = case_when(
      has_exonic & has_intronic ~ "Both exonic & intronic",
      has_exonic ~ "Exonic only",
      has_intronic ~ "Intronic only",
      TRUE ~ "Other"
    )
  ) %>%
  select(Gene_name, total_samples_affected, total_te_insertions, te_location_type, location_summary, gene_features)

gene_summary_kics_file <- paste0(plot_dir, "genes_affected/gene_summary_te_kics_germline.csv")
write.csv(gene_summary_kics, gene_summary_kics_file, row.names = FALSE)
cat(sprintf("✓ Gene summary for KICS saved to: %s\n", basename(gene_summary_kics_file)))

# Generate gene summary for LFS samples
cat("Generating gene summary table for LFS samples...\n")
gene_summary_lfs <- te_lfs_split %>%
  group_by(Gene_name) %>%
  summarize(
    total_samples_affected = n_distinct(sample),
    total_te_insertions = n(),
    has_exonic = any(grepl("exon", Location, ignore.case = TRUE), na.rm = TRUE),
    has_intronic = any(grepl("intron", Location, ignore.case = TRUE), na.rm = TRUE),
    location_summary = paste(unique(na.omit(Location)), collapse = "; "),
    gene_features = paste(unique(na.omit(Location2)), collapse = "; "),
    .groups = "drop"
  ) %>%
  arrange(desc(total_samples_affected), desc(total_te_insertions)) %>%
  mutate(
    te_location_type = case_when(
      has_exonic & has_intronic ~ "Both exonic & intronic",
      has_exonic ~ "Exonic only",
      has_intronic ~ "Intronic only",
      TRUE ~ "Other"
    )
  ) %>%
  select(Gene_name, total_samples_affected, total_te_insertions, te_location_type, location_summary, gene_features)

gene_summary_lfs_file <- paste0(plot_dir, "genes_affected/gene_summary_te_lfs_germline.csv")
write.csv(gene_summary_lfs, gene_summary_lfs_file, row.names = FALSE)
cat(sprintf("✓ Gene summary for LFS saved to: %s\n", basename(gene_summary_lfs_file)))

# Generate cancer genes only gene summary for all samples
cat("Generating cancer genes only gene summary for all samples...\n")
gene_summary_cancer <- gene_summary %>%
  filter(Gene_name %in% genes)

gene_summary_cancer_file <- paste0(plot_dir, "genes_affected/gene_summary_te_all_cancer_genes_germline.csv")
write.csv(gene_summary_cancer, gene_summary_cancer_file, row.names = FALSE)
cat(sprintf("✓ Cancer genes only gene summary saved to: %s\n", basename(gene_summary_cancer_file)))
cat(sprintf("  Total cancer genes with TEs: %d\n\n", nrow(gene_summary_cancer)))

# Generate cancer genes only gene summary for KICS samples
cat("Generating cancer genes only gene summary for KICS samples...\n")
gene_summary_kics_cancer <- gene_summary_kics %>%
  filter(Gene_name %in% genes)

gene_summary_kics_cancer_file <- paste0(plot_dir, "genes_affected/gene_summary_te_kics_cancer_genes_germline.csv")
write.csv(gene_summary_kics_cancer, gene_summary_kics_cancer_file, row.names = FALSE)
cat(sprintf("✓ Cancer genes only gene summary for KICS saved to: %s\n", basename(gene_summary_kics_cancer_file)))
cat(sprintf("  Total cancer genes with TEs (KICS): %d\n\n", nrow(gene_summary_kics_cancer)))

# Generate cancer genes only gene summary for LFS samples
cat("Generating cancer genes only gene summary for LFS samples...\n")
gene_summary_lfs_cancer <- gene_summary_lfs %>%
  filter(Gene_name %in% genes)

gene_summary_lfs_cancer_file <- paste0(plot_dir, "genes_affected/gene_summary_te_lfs_cancer_genes_germline.csv")
write.csv(gene_summary_lfs_cancer, gene_summary_lfs_cancer_file, row.names = FALSE)
cat(sprintf("✓ Cancer genes only gene summary for LFS saved to: %s\n", basename(gene_summary_lfs_cancer_file)))
cat(sprintf("  Total cancer genes with TEs (LFS): %d\n\n", nrow(gene_summary_lfs_cancer)))

#### DETAILED SAMPLE-LEVEL CANCER GENES OUTPUT ####
cat("\n===== DETAILED SAMPLE-LEVEL CANCER GENES OUTPUT =====\n")

# Define clinical columns to include
clinical_cols <- c("sample", "TP53_status", "Cancer", "tumor_type", "sex",
                   "age_at_diagnosis", "cohort", "mapped_label", "predicted_ancestry_thres")

# Function to create sample-level cancer genes output
create_sample_cancer_genes_output <- function(te_split_df, cancer_genes, output_file, cohort_name) {
  # Filter to cancer genes only
  te_cancer <- te_split_df %>%
    filter(Gene_name %in% cancer_genes)

  if (nrow(te_cancer) == 0) {
    cat(sprintf("  No cancer gene TEs found for %s\n", cohort_name))
    return(NULL)
  }

  # Get available clinical columns
  avail_clinical <- intersect(clinical_cols, colnames(te_cancer))

  te_output <- te_cancer %>%
    select(Gene_name, all_of(avail_clinical),
           SV_chrom, SV_start, SV_end, ALT,
           any_of(c("Location", "Location2", "SV_length"))) %>%
    arrange(Gene_name, sample)

  write.csv(te_output, output_file, row.names = FALSE)
  cat(sprintf("✓ Sample-level cancer genes for %s saved to: %s\n", cohort_name, basename(output_file)))
  cat(sprintf("  Samples: %d, Genes: %d, TEs: %d\n",
              n_distinct(te_output$sample), n_distinct(te_output$Gene_name), nrow(te_output)))

  return(te_output)
}

# All samples
cat("Generating sample-level cancer genes for all samples...\n")
sample_cancer_all <- create_sample_cancer_genes_output(
  te_all_split, genes,
  paste0(plot_dir, "genes_affected/gene_sample_te_all_cancer_genes_germline.csv"),
  "All"
)

# KICS samples
cat("Generating sample-level cancer genes for KICS samples...\n")
sample_cancer_kics <- create_sample_cancer_genes_output(
  te_kics_split, genes,
  paste0(plot_dir, "genes_affected/gene_sample_te_kics_cancer_genes_germline.csv"),
  "KICS"
)

# LFS samples
cat("Generating sample-level cancer genes for LFS samples...\n")
sample_cancer_lfs <- create_sample_cancer_genes_output(
  te_lfs_split, genes,
  paste0(plot_dir, "genes_affected/gene_sample_te_lfs_cancer_genes_germline.csv"),
  "LFS"
)

cat("✓ Script completed successfully\n")

# Close module-specific sink
close_module_sink()
