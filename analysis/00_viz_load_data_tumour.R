#!/usr/bin/env Rscript

# Load Tumour Data for TE Visualization
# This script loads all necessary data for tumour TE visualization
# Source this after 00_viz_common_setup.R
#
# USAGE: Set REQUIRED_DATA before sourcing to load only what you need:
#   REQUIRED_DATA <- c("count_matrix", "expand", "clinical")
#   source("00_viz_load_data_tumour.R")
#
# Available data categories:
#   "count_matrix" - te_*_t count matrix variants
#   "expand"       - te_*_expand_t variants
#   "split"        - te_*_split_t variants (for pathway/gene analysis)
#   "rna"          - kics_rna, lfs_rna, stjude_rna, matched_dna_rna, lfs_wgs2rna
#   "ancestry"     - ancestry data (loaded via clinical)
#   "clinical"     - clinical, metrics
#   "genes"        - cpg, gene_size
#   "loh"          - loh_time, somatic_conversion, kics_tp53_somatic_variants
#   "survival"     - kics_DOD
#   "sv"           - tumour_sv structural variants
#   "common"       - common TE data (if COMMON_MODE=TRUE)
#   "all"          - load everything (default if REQUIRED_DATA not set)

# TE type filtering options (must match processing script settings)
INCLUDE_ALU <- FALSE  # Set to TRUE to include ALU elements
INCLUDE_SVA <- FALSE  # Set to TRUE to include SVA elements

# Analysis mode - determines which TE data to load and output directory
COMMON_MODE <- FALSE  # Set TRUE to use common TEs (no frequency filtering), outputs to /common/ subdir
FULLLENGTH_YOUNG_MODE <- FALSE  # Set TRUE to use full-length young TEs (L1HS ≥5900bp, AluY, SVA_E/F), outputs to /fulllength_young/ subdir
# Default (both FALSE) = rare TEs, outputs to /rare/ subdir



# Helper function to check if data category is needed
needs_data <- function(category) {
  if (!exists("REQUIRED_DATA", envir = .GlobalEnv)) return(TRUE)
  req <- get("REQUIRED_DATA", envir = .GlobalEnv)
  if ("all" %in% req) return(TRUE)
  category %in% req
}

# Check if data is already loaded (skip reloading when running via ALL script)
if (exists("te_all_t") && exists("te_aff_expand_t") && exists("te_all_split_t")) {
  cat("Data already loaded, skipping reload...\n")
  return(invisible(NULL))
}

# Set FULLLENGTH_YOUNG_MODE based on REQUIRED_DATA if not already set
if (!exists("FULLLENGTH_YOUNG_MODE") || !FULLLENGTH_YOUNG_MODE) {
  if (exists("REQUIRED_DATA", envir = .GlobalEnv)) {
    req <- get("REQUIRED_DATA", envir = .GlobalEnv)
    if ("fulllength_young" %in% req) {
      FULLLENGTH_YOUNG_MODE <- TRUE
    }
  }
}

cat("Loading tumour data...\n")

#### SETUP PATHS ####
# Base directory for plots
plot_dir_base <- "/Users/briannelaverty/Library/Mobile Documents/com~apple~CloudDocs/Malkin_lab/figures/te/graphs/tumour/"

# Use subdirectory based on analysis mode (fulllength_young > common > rare)
if (exists("FULLLENGTH_YOUNG_MODE") && FULLLENGTH_YOUNG_MODE) {
  plot_dir <- paste0(plot_dir_base, "fulllength_young/")
  cat("*** FULLLENGTH YOUNG MODE ACTIVE: Using full-length young TEs, output to /fulllength_young/ ***\n")
} else if (exists("COMMON_MODE") && COMMON_MODE) {
  plot_dir <- paste0(plot_dir_base, "common/")
  cat("*** COMMON MODE ACTIVE: Using unfiltered TEs, output to /common/ ***\n")
} else {
  plot_dir <- paste0(plot_dir_base, "rare/")
}
r_dir_files <- paste0(plot_dir, "files/")

# Set data directory based on TEST_MODE
if (TEST_MODE) {
  r_dir <- "/Users/briannelaverty/Documents/R_Malkin/te/data/R_obj/test_data/tumour/"
  cat("*** TEST MODE ACTIVE: Using test data from", r_dir, "***\n\n")
} else {
  r_dir <- "/Users/briannelaverty/Documents/R_Malkin/te/data/R_obj/tumour/"
}


# Create directories if they don't exist
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
if (!dir.exists(r_dir_files)) dir.create(r_dir_files, recursive = TRUE)

# Create subdirectories for organized plot storage
subdirs <- c("dataset", "general", "counts_cohort", "counts_tumour_type",
             "counts_clinical_lfs", "counts_clinical_kics", "ancestry",
             "pca_umap", "cancer_genes", "pathway", "reg_element", "rna",
             "survival_burden", "loh", "re_rna", "other")
for (subdir in subdirs) {
  subdir_path <- paste0(plot_dir, subdir, "/")
  if (!dir.exists(subdir_path)) dir.create(subdir_path, recursive = TRUE)
}

# Clear output files at start of each run
pdf_output_file <- paste0(plot_dir, "graph_output.pdf")
if (file.exists(pdf_output_file)) file.remove(pdf_output_file)

# Open a PDF device at the start to capture all plots in a single PDF
pdf(file = paste0(plot_dir, "graph_output.pdf"), width = 12, height = 8)

# Covariate vectors for linear models
covar_med <- c("predicted_ancestry_thres", "age_at_diagnosis", "med_cov", "total_reads", "med_read_len", "pct_chimeras", "avg_quality", "tumor_type_grouped")
covar_no_ancestry <- c("age_at_diagnosis", "med_cov", "total_reads", "med_read_len", "pct_chimeras", "avg_quality", "tumor_type_grouped")
covar_no_age <- c("predicted_ancestry_thres", "med_cov", "total_reads", "med_read_len", "pct_chimeras", "avg_quality", "tumor_type_grouped")
covar_no_tumour_type <- c("predicted_ancestry_thres", "age_at_diagnosis", "med_cov", "total_reads", "med_read_len", "pct_chimeras", "avg_quality")

# Helper function to add tumor_type_grouped column based on each df's distribution
# Groups tumor types with <min_count samples as "Other"
add_tumor_type_grouped <- function(df, min_count = 3) {
  if ("tumor_type" %in% colnames(df)) {
    tumor_type_counts <- table(df$tumor_type)
    major_types <- names(tumor_type_counts[tumor_type_counts >= min_count])
    df$tumor_type_grouped <- factor(ifelse(
      df$tumor_type %in% major_types,
      as.character(df$tumor_type),
      "Other"
    ))
  }
  return(df)
}

#### LOAD ADDITIONAL DATA ####
if (needs_data("loh")) {
  cat("Loading LOH data files...\n")
  loh_time <- read.delim("/Users/briannelaverty/Documents/R_Malkin/te/data/final/loh_time", sep="\t", header=TRUE)
  somatic_conversion <- read.delim("/Users/briannelaverty/Documents/R_Malkin/te/data/final/somatic_variant_conversion.txt", sep="\t", header=TRUE)
  kics_tp53_somatic_variants <- fread("/Users/briannelaverty/Documents/R_Malkin/TE/data/kics_tp53_somatic_variants.csv", sep=",", header=TRUE)
}

if (needs_data("sv")) {
  cat("Loading structural variant data...\n")
  tumour_sv <- read.delim("/Users/briannelaverty/Documents/R_Malkin/te/data/sv/vis_sv_data_combined.tsv", sep="\t", header=TRUE)
}

if (needs_data("survival")) {
  cat("Loading survival data...\n")
  kics_DOD <- readxl::read_excel("/Users/briannelaverty/Documents/R_Malkin/clinical/kics_DOD.xlsx")
}

#### LOAD RNA DATA ####
if (needs_data("rna")) {
  cat("Loading RNA data files...\n")
  kics_rna <- read.delim("/Users/briannelaverty/Documents/R_Malkin/te/data/rna/kics_FPKM.tsv", sep="\t", header=TRUE, check.names = FALSE)
  lfs_rna <- read.delim("/Users/briannelaverty/Documents/R_Malkin/te/data/rna/lfs_FPKM.tsv", sep="\t", header=TRUE, check.names = FALSE)
  lfs_wgs2rna <- read.delim("/Users/briannelaverty/Documents/R_Malkin/te/data/rna/lfs_wgs2rna.csv", sep=",", header=TRUE)
  stjude_rna <- read.delim("/Users/briannelaverty/Documents/R_Malkin/te/data/rna/stjude_FPKM.tsv", sep="\t", header=TRUE, check.names = FALSE)
  stjude_rna_names <- read.delim("/Users/briannelaverty/Documents/R_Malkin/te/data/rna/stjude_rna_names.txt", sep="\t", header=FALSE)
  matched_dna_rna <- read.csv("/Users/briannelaverty/Documents/R_Malkin/te/data/rna/matched_dna_rna.csv", stringsAsFactors = FALSE)
}

#### LOAD SAVED PROCESSED DATA ####
cat("Loading saved processed tumour data objects...\n")

# Load prepared data objects (use production dir for metadata in test mode)
metadata_dir <- if (TEST_MODE) "/Users/briannelaverty/Documents/R_Malkin/te/data/R_obj/tumour/" else r_dir

if (needs_data("genes")) {
  load(paste0(metadata_dir, "gene_size.RData"))
  cpg <- read.csv("/Users/briannelaverty/Documents/R_Malkin/te/data/raw/kics_cpg.txt", header=FALSE, sep="\t") # cancer predisposition genes
}

# Always load core metadata
load(paste0(metadata_dir, "chr_lengths.RData"))
load(paste0(metadata_dir, "hostseq_cancer.RData"))

if (needs_data("clinical")) {
  load(paste0(metadata_dir, "clinical.RData"))
  load(paste0(metadata_dir, "metrics.RData"))
}

# TE DATA - COUNT MATRIX FORMAT
if (needs_data("count_matrix")) {
  if (FULLLENGTH_YOUNG_MODE) {
    cat("Loading FULLLENGTH YOUNG TE count matrix data...\n")
    load(paste0(r_dir, "final_te_count_fulllength_young.RData"))
    data_source <- final_te_count_fulllength_young_t
  } else if (COMMON_MODE) {
    cat("Loading COMMON TE count matrix data...\n")
    load(paste0(r_dir, "final_te_count_common.RData"))
    data_source <- final_te_count_common_t
  } else {
    cat("Loading rare TE count matrix data...\n")
    load(paste0(r_dir, "final_te_count_rare.RData"))
    data_source <- if (TEST_MODE) final_te_count_t_test else final_te_count_t
  }

  te_aff_t <- data_source$te_aff_selected
  te_aff_all_t <- data_source$te_aff_all
  te_lfs_t <- data_source$te_lfs_selected
  te_lfs_all_t <- data_source$te_lfs_all
  te_kics_t <- data_source$te_kics_selected
  te_kics_all_t <- data_source$te_kics_all
  te_taylor_t <- data_source$te_taylor_selected
  te_taylor_all_t <- data_source$te_taylor_all
  te_all_t <- data_source$te_all_selected
  te_all_all_t <- data_source$te_all_all

  # Add tumor_type_grouped to each dataframe based on its own distribution
  te_aff_t <- add_tumor_type_grouped(te_aff_t)
  te_aff_all_t <- add_tumor_type_grouped(te_aff_all_t)
  te_lfs_t <- add_tumor_type_grouped(te_lfs_t)
  te_lfs_all_t <- add_tumor_type_grouped(te_lfs_all_t)
  te_kics_t <- add_tumor_type_grouped(te_kics_t)
  te_kics_all_t <- add_tumor_type_grouped(te_kics_all_t)
  te_taylor_t <- add_tumor_type_grouped(te_taylor_t)
  te_taylor_all_t <- add_tumor_type_grouped(te_taylor_all_t)
  te_all_t <- add_tumor_type_grouped(te_all_t)
  te_all_all_t <- add_tumor_type_grouped(te_all_all_t)
}

# TE DATA - EXPANDED FORMAT
if (needs_data("expand")) {
  if (FULLLENGTH_YOUNG_MODE) {
    cat("Loading FULLLENGTH YOUNG TE expanded data...\n")
    load(paste0(r_dir, "final_te_count_expand_fulllength_young.RData"))
    data_source <- final_te_count_expand_fulllength_young_t
  } else if (COMMON_MODE) {
    cat("Loading COMMON TE expanded data...\n")
    load(paste0(r_dir, "final_te_count_expand_common.RData"))
    data_source <- final_te_count_expand_common_t
  } else {
    cat("Loading rare TE expanded data...\n")
    load(paste0(r_dir, "final_te_count_expand_rare.RData"))
    data_source <- if (TEST_MODE) final_te_count_expand_t_test else final_te_count_expand_t
  }

  te_aff_expand_t <- data_source$te_aff_selected
  te_aff_expand_all_t <- data_source$te_aff_all
  te_lfs_expand_t <- data_source$te_lfs_selected
  te_lfs_expand_all_t <- data_source$te_lfs_all
  te_kics_expand_t <- data_source$te_kics_selected
  te_kics_expand_all_t <- data_source$te_kics_all
  te_taylor_expand_t <- data_source$te_taylor_selected
  te_taylor_expand_all_t <- data_source$te_taylor_all
  te_all_expand_t <- data_source$te_all_selected
  te_all_expand_all_t <- data_source$te_all_all

  # Add tumor_type_grouped to each dataframe based on its own distribution
  te_aff_expand_t <- add_tumor_type_grouped(te_aff_expand_t)
  te_aff_expand_all_t <- add_tumor_type_grouped(te_aff_expand_all_t)
  te_lfs_expand_t <- add_tumor_type_grouped(te_lfs_expand_t)
  te_lfs_expand_all_t <- add_tumor_type_grouped(te_lfs_expand_all_t)
  te_kics_expand_t <- add_tumor_type_grouped(te_kics_expand_t)
  te_kics_expand_all_t <- add_tumor_type_grouped(te_kics_expand_all_t)
  te_taylor_expand_t <- add_tumor_type_grouped(te_taylor_expand_t)
  te_taylor_expand_all_t <- add_tumor_type_grouped(te_taylor_expand_all_t)
  te_all_expand_t <- add_tumor_type_grouped(te_all_expand_t)
  te_all_expand_all_t <- add_tumor_type_grouped(te_all_expand_all_t)
}

# TE DATA - SPLIT BY GENE FORMAT
if (needs_data("split")) {
  if (FULLLENGTH_YOUNG_MODE) {
    cat("Loading FULLLENGTH YOUNG TE split (gene) data...\n")
    load(paste0(r_dir, "final_te_count_split_fulllength_young.RData"))
    data_source <- final_te_count_expand_split_fulllength_young_t
  } else if (COMMON_MODE) {
    cat("Loading COMMON TE split (gene) data...\n")
    load(paste0(r_dir, "final_te_count_split_common.RData"))
    data_source <- final_te_count_expand_split_common_t
  } else {
    cat("Loading rare TE split (gene) data...\n")
    load(paste0(r_dir, "final_te_count_expand_split_split_rare.RData"))
    data_source <- if (TEST_MODE) final_te_count_expand_split_split_rare_t_test else final_te_count_expand_split_split_rare_t
  }

  te_aff_split_t <- data_source$te_aff_selected
  te_lfs_split_t <- data_source$te_lfs_selected
  te_kics_split_t <- data_source$te_kics_selected
  te_taylor_split_t <- data_source$te_taylor_selected
  te_all_split_t <- data_source$te_all_selected

  # Add tumor_type_grouped to each dataframe based on its own distribution
  te_aff_split_t <- add_tumor_type_grouped(te_aff_split_t)
  te_lfs_split_t <- add_tumor_type_grouped(te_lfs_split_t)
  te_kics_split_t <- add_tumor_type_grouped(te_kics_split_t)
  te_taylor_split_t <- add_tumor_type_grouped(te_taylor_split_t)
  te_all_split_t <- add_tumor_type_grouped(te_all_split_t)

  # Create aliases for gene analysis (used by some scripts)
  te_aff_split_genes_t <- te_aff_split_t
  te_lfs_split_genes_t <- te_lfs_split_t
  te_kics_split_genes_t <- te_kics_split_t
  te_taylor_split_genes_t <- te_taylor_split_t
  te_all_split_genes_t <- te_all_split_t
}

# Print summary of what was loaded
cat("✓ Tumour data loaded successfully\n")
if (exists("REQUIRED_DATA", envir = .GlobalEnv)) {
  req <- get("REQUIRED_DATA", envir = .GlobalEnv)
  cat("  Loaded categories:", paste(req, collapse = ", "), "\n")
}
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
