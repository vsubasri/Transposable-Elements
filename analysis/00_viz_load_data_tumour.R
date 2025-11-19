#!/usr/bin/env Rscript

# Load Tumour Data for TE Visualization
# This script loads all necessary data for tumour TE visualization
# Source this after 00_viz_common_setup.R

# Check if data is already loaded (skip reloading when running via ALL script)
if (exists("te_all_t") && exists("te_aff_expand_t") && exists("te_all_split_t")) {
  cat("Data already loaded, skipping reload...\n")
  return(invisible(NULL))
}

cat("Loading tumour data...\n")

#### SETUP PATHS ####
plot_dir <- "/Users/briannelaverty/Library/Mobile Documents/com~apple~CloudDocs/Malkin_lab/figures/te/graphs/tumour/"
r_dir_files <- "/Users/briannelaverty/Library/Mobile Documents/com~apple~CloudDocs/Malkin_lab/figures/te/graphs/tumour/files/"

# Set data directory based on TEST_MODE
if (TEST_MODE) {
  r_dir <- "/Users/briannelaverty/Documents/R_Malkin/te/data/R_obj/test_data/tumour/"
  cat("*** TEST MODE ACTIVE: Using test data from", r_dir, "***\n\n")
} else {
  r_dir <- "/Users/briannelaverty/Documents/R_Malkin/te/data/R_obj/tumour/"
}

# TE type filtering options (must match processing script settings)
INCLUDE_ALU <- FALSE  # Set to TRUE to include ALU elements
INCLUDE_SVA <- FALSE  # Set to TRUE to include SVA elements
PROCESS_COMMON_TES <- FALSE  # Set to TRUE to visualize common TEs

# Create directories if they don't exist
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
if (!dir.exists(r_dir_files)) dir.create(r_dir_files, recursive = TRUE)

# Create subdirectories for organized plot storage
subdirs <- c("dataset", "general", "counts_cohort", "counts_tumour_type",
             "counts_clinical_lfs", "counts_clinical_kics", "ancestry",
             "pca_umap", "cancer_genes", "pathway", "reg_element", "rna",
             "survival_burden", "other")
for (subdir in subdirs) {
  subdir_path <- paste0(plot_dir, subdir, "/")
  if (!dir.exists(subdir_path)) dir.create(subdir_path, recursive = TRUE)
}

# Clear output files at start of each run
pdf_output_file <- paste0(plot_dir, "graph_output.pdf")
if (file.exists(pdf_output_file)) file.remove(pdf_output_file)

# Redirect all stdout to text file
stdout_file <- paste0(r_dir_files, "graph_output_text.txt")
if (file.exists(stdout_file)) file.remove(stdout_file)
sink(stdout_file, split = TRUE)

# Open a PDF device at the start to capture all plots in a single PDF
pdf(file = paste0(plot_dir, "graph_output.pdf"), width = 12, height = 8)

# Covariates for linear models
covar_med <- c("med_cov", "total_reads", "med_read_len", "pct_chimeras", "avg_quality", "tumor_type")

#### LOAD ADDITIONAL DATA ####
cat("Loading additional raw data files...\n")
loh_time<- read.delim("/Users/briannelaverty/Documents/R_Malkin/te/data/final/loh_time", sep="\t", header=TRUE)
somatic_conversion <- read.delim("/Users/briannelaverty/Documents/R_Malkin/te/data/final/somatic_variant_conversion.txt", sep="\t", header=TRUE)
kics_tp53_somatic_variants <- fread("/Users/briannelaverty/Documents/R_Malkin/TE/data/kics_tp53_somatic_variants.csv", sep=",", header=TRUE)
tumour_sv <- read.delim("/Users/briannelaverty/Documents/R_Malkin/te/data/sv/vis_sv_data_combined.tsv", sep="\t", header=TRUE)
kics_DOD <- readxl::read_excel("/Users/briannelaverty/Documents/R_Malkin/clinical/kics_DOD.xlsx")

#### LOAD RNA DATA ####
cat("Loading RNA data files...\n")
kics_rna <- read.delim("/Users/briannelaverty/Documents/R_Malkin/te/data/rna/kics_FPKM.tsv", sep="\t", header=TRUE, check.names = FALSE)
lfs_rna <- read.delim("/Users/briannelaverty/Documents/R_Malkin/te/data/rna/lfs_FPKM.tsv", sep="\t", header=TRUE, check.names = FALSE)
lfs_wgs2rna <- read.delim("/Users/briannelaverty/Documents/R_Malkin/te/data/rna/lfs_wgs2rna.csv", sep=",", header=TRUE)
stjude_rna <- read.delim("/Users/briannelaverty/Documents/R_Malkin/te/data/rna/stjude_FPKM.tsv", sep="\t", header=TRUE, check.names = FALSE)
stjude_rna_names <- read.delim("/Users/briannelaverty/Documents/R_Malkin/te/data/rna/stjude_rna_names.txt", sep="\t", header=FALSE)
matched_dna_rna <- read.csv("/Users/briannelaverty/Documents/R_Malkin/te/data/rna/matched_dna_rna.csv", stringsAsFactors = FALSE)

#### LOAD SAVED PROCESSED DATA ####
cat("Loading saved processed tumour data objects...\n")

# Load prepared data objects (use production dir for metadata in test mode)
metadata_dir <- if (TEST_MODE) "/Users/briannelaverty/Documents/R_Malkin/te/data/R_obj/tumour/" else r_dir

load(paste0(metadata_dir, "gene_size.RData"))
load(paste0(metadata_dir, "chr_lengths.RData"))
load(paste0(metadata_dir, "hostseq_cancer.RData"))
load(paste0(metadata_dir, "clinical.RData"))
load(paste0(metadata_dir, "metrics.RData"))

# Load prepared data objects
load(paste0(r_dir, "final_te_count_rare.RData"))

if (TEST_MODE) {
  te_aff_t <- final_te_count_t_test$te_aff_selected
  te_aff_all_t <- final_te_count_t_test$te_aff_all
  te_lfs_t <- final_te_count_t_test$te_lfs_selected
  te_lfs_all_t <- final_te_count_t_test$te_lfs_all
  te_kics_t <- final_te_count_t_test$te_kics_selected
  te_kics_all_t <- final_te_count_t_test$te_kics_all
  te_taylor_t <- final_te_count_t_test$te_taylor_selected
  te_taylor_all_t <- final_te_count_t_test$te_taylor_all
  te_all_t <- final_te_count_t_test$te_all_selected
  te_all_all_t <- final_te_count_t_test$te_all_all
} else {
  te_aff_t <- final_te_count_t$te_aff_selected
  te_aff_all_t <- final_te_count_t$te_aff_all
  te_lfs_t <- final_te_count_t$te_lfs_selected
  te_lfs_all_t <- final_te_count_t$te_lfs_all
  te_kics_t <- final_te_count_t$te_kics_selected
  te_kics_all_t <- final_te_count_t$te_kics_all
  te_taylor_t <- final_te_count_t$te_taylor_selected
  te_taylor_all_t <- final_te_count_t$te_taylor_all
  te_all_t <- final_te_count_t$te_all_selected
  te_all_all_t <- final_te_count_t$te_all_all
}

# Rare TE expanded format
load(paste0(r_dir, "final_te_count_expand_rare.RData"))

if (TEST_MODE) {
  te_aff_expand_t <- final_te_count_expand_t_test$te_aff_selected
  te_aff_expand_all_t <- final_te_count_expand_t_test$te_aff_all
  te_lfs_expand_t <- final_te_count_expand_t_test$te_lfs_selected
  te_lfs_expand_all_t <- final_te_count_expand_t_test$te_lfs_all
  te_kics_expand_t <- final_te_count_expand_t_test$te_kics_selected
  te_kics_expand_all_t <- final_te_count_expand_t_test$te_kics_all
  te_taylor_expand_t <- final_te_count_expand_t_test$te_taylor_selected
  te_taylor_expand_all_t <- final_te_count_expand_t_test$te_taylor_all
  te_all_expand_t <- final_te_count_expand_t_test$te_all_selected
  te_all_expand_all_t <- final_te_count_expand_t_test$te_all_all
} else {
  te_aff_expand_t <- final_te_count_expand_t$te_aff_selected
  te_aff_expand_all_t <- final_te_count_expand_t$te_aff_all
  te_lfs_expand_t <- final_te_count_expand_t$te_lfs_selected
  te_lfs_expand_all_t <- final_te_count_expand_t$te_lfs_all
  te_kics_expand_t <- final_te_count_expand_t$te_kics_selected
  te_kics_expand_all_t <- final_te_count_expand_t$te_kics_all
  te_taylor_expand_t <- final_te_count_expand_t$te_taylor_selected
  te_taylor_expand_all_t <- final_te_count_expand_t$te_taylor_all
  te_all_expand_t <- final_te_count_expand_t$te_all_selected
  te_all_expand_all_t <- final_te_count_expand_t$te_all_all
}

# Rare TE split by gene format (annotSV split mode)
load(paste0(r_dir, "final_te_count_expand_split_split_rare.RData"))

if (TEST_MODE) {
  te_aff_split_t <- final_te_count_expand_split_split_rare_t_test$te_aff_selected
  te_lfs_split_t <- final_te_count_expand_split_split_rare_t_test$te_lfs_selected
  te_kics_split_t <- final_te_count_expand_split_split_rare_t_test$te_kics_selected
  te_taylor_split_t <- final_te_count_expand_split_split_rare_t_test$te_taylor_selected
  te_all_split_t <- final_te_count_expand_split_split_rare_t_test$te_all_selected
} else {
  te_aff_split_t <- final_te_count_expand_split_split_rare_t$te_aff_selected
  te_lfs_split_t <- final_te_count_expand_split_split_rare_t$te_lfs_selected
  te_kics_split_t <- final_te_count_expand_split_split_rare_t$te_kics_selected
  te_taylor_split_t <- final_te_count_expand_split_split_rare_t$te_taylor_selected
  te_all_split_t <- final_te_count_expand_split_split_rare_t$te_all_selected
}

# Create aliases for gene analysis (used by some scripts)
te_aff_split_genes_t <- te_aff_split_t
te_lfs_split_genes_t <- te_lfs_split_t
te_kics_split_genes_t <- te_kics_split_t
te_taylor_split_genes_t <- te_taylor_split_t
te_all_split_genes_t <- te_all_split_t

cat("✓ Tumour data loaded successfully\n")
cat("\n")
cat("========================================\n")
cat("  TUMOUR DATASET SUMMARY\n")
cat("========================================\n")
cat("Cohorts: LFS, KICS, Taylor\n\n")
cat("Available datasets by sample selection:\n\n")
cat("te_all_*_t       : All tumour samples (LFS + KICS)\n")
cat("te_aff_*_t       : Affected samples\n")
cat("te_lfs_*_t       : LFS-related samples\n")
cat("te_kics_*_t      : KICS cohort only\n")
cat("\nData formats:\n")
cat("- *_expand_t     : One row per TE insertion (for TE-level analysis)\n")
cat("- *_split_t      : One row per TE-gene overlap (for gene/pathway analysis)\n")
cat("- (no suffix)    : Count matrix format (one row per sample)\n")
cat("========================================\n\n")
