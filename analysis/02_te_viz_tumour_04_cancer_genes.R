#!/usr/bin/env Rscript

# Tumour TE Visualization - Cancer Genes
# Cancer genes analysis and gene summary tables

# Source common setup and load data
source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/00_viz_common_setup.R")
REQUIRED_DATA <- c("split", "clinical", "genes")
source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/00_viz_load_data_tumour.R")

cat("Running 02_te_viz_tumour_04_cancer_genes.R...\n")

#### GENE SUMMARY TABLES ####
# Summary of all genes affected by TEs with location information

cat("\n===== GENE SUMMARY TABLES =====\n")

# Create genes_affected directory if it doesn't exist
dir.create(paste0(plot_dir, "genes_affected"), showWarnings = FALSE, recursive = TRUE)

# Generate gene summary table for all affected samples
cat("Generating gene summary table...\n")
gene_summary <- te_all_split_t %>%
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
gene_summary_file <- paste0(plot_dir, "genes_affected/gene_summary_te_affected_tumour.csv")
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
gene_summary_kics <- te_kics_split_t %>%
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

gene_summary_kics_file <- paste0(plot_dir, "genes_affected/gene_summary_te_kics_tumour.csv")
write.csv(gene_summary_kics, gene_summary_kics_file, row.names = FALSE)
cat(sprintf("✓ Gene summary for KICS saved to: %s\n", basename(gene_summary_kics_file)))

# Generate gene summary for LFS samples
cat("Generating gene summary table for LFS samples...\n")
gene_summary_lfs <- te_lfs_split_t %>%
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

gene_summary_lfs_file <- paste0(plot_dir, "genes_affected/gene_summary_te_lfs_tumour.csv")
write.csv(gene_summary_lfs, gene_summary_lfs_file, row.names = FALSE)
cat(sprintf("✓ Gene summary for LFS saved to: %s\n", basename(gene_summary_lfs_file)))

# Generate cancer genes only gene summary for all samples
cat("Generating cancer genes only gene summary for all samples...\n")
genes <- cpg$V1
gene_summary_cancer <- gene_summary %>%
  filter(Gene_name %in% genes)

gene_summary_cancer_file <- paste0(plot_dir, "genes_affected/gene_summary_te_all_cancer_genes_tumour.csv")
write.csv(gene_summary_cancer, gene_summary_cancer_file, row.names = FALSE)
cat(sprintf("✓ Cancer genes only gene summary saved to: %s\n", basename(gene_summary_cancer_file)))
cat(sprintf("  Total cancer genes with TEs: %d\n\n", nrow(gene_summary_cancer)))

# Generate cancer genes only gene summary for KICS samples
cat("Generating cancer genes only gene summary for KICS samples...\n")
gene_summary_kics_cancer <- gene_summary_kics %>%
  filter(Gene_name %in% genes)

gene_summary_kics_cancer_file <- paste0(plot_dir, "genes_affected/gene_summary_te_kics_cancer_genes_tumour.csv")
write.csv(gene_summary_kics_cancer, gene_summary_kics_cancer_file, row.names = FALSE)
cat(sprintf("✓ Cancer genes only gene summary for KICS saved to: %s\n", basename(gene_summary_kics_cancer_file)))
cat(sprintf("  Total cancer genes with TEs (KICS): %d\n\n", nrow(gene_summary_kics_cancer)))

# Generate cancer genes only gene summary for LFS samples
cat("Generating cancer genes only gene summary for LFS samples...\n")
gene_summary_lfs_cancer <- gene_summary_lfs %>%
  filter(Gene_name %in% genes)

gene_summary_lfs_cancer_file <- paste0(plot_dir, "genes_affected/gene_summary_te_lfs_cancer_genes_tumour.csv")
write.csv(gene_summary_lfs_cancer, gene_summary_lfs_cancer_file, row.names = FALSE)
cat(sprintf("✓ Cancer genes only gene summary for LFS saved to: %s\n", basename(gene_summary_lfs_cancer_file)))
cat(sprintf("  Total cancer genes with TEs (LFS): %d\n\n", nrow(gene_summary_lfs_cancer)))

#### DETAILED SAMPLE-LEVEL CANCER GENES OUTPUT ####
cat("\n===== DETAILED SAMPLE-LEVEL CANCER GENES OUTPUT =====\n")

# Define clinical columns to include
clinical_cols <- c("sample", "TP53_status", "tumor_type", "sex",
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
  te_all_split_t, genes,
  paste0(plot_dir, "genes_affected/gene_sample_te_all_cancer_genes_tumour.csv"),
  "All"
)

# KICS samples
cat("Generating sample-level cancer genes for KICS samples...\n")
sample_cancer_kics <- create_sample_cancer_genes_output(
  te_kics_split_t, genes,
  paste0(plot_dir, "genes_affected/gene_sample_te_kics_cancer_genes_tumour.csv"),
  "KICS"
)

# LFS samples
cat("Generating sample-level cancer genes for LFS samples...\n")
sample_cancer_lfs <- create_sample_cancer_genes_output(
  te_lfs_split_t, genes,
  paste0(plot_dir, "genes_affected/gene_sample_te_lfs_cancer_genes_tumour.csv"),
  "LFS"
)

cat("✓ Script completed successfully\n")
