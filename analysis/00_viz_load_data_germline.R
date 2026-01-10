#!/usr/bin/env Rscript

# Load Germline Data for TE Visualization
# This script loads all necessary data for germline TE visualization
# Source this after 00_viz_common_setup.R
#
# USAGE: Set REQUIRED_DATA before sourcing to load only what you need:
#   REQUIRED_DATA <- c("count_matrix", "expand", "clinical")
#   source("00_viz_load_data_germline.R")
#
# Available data categories:
#   "count_matrix" - te_all, te_aff, te_lfs, te_kics, te_taylor, te_hostseq, te_lfs_mut_wt
#   "expand"       - te_*_expand variants
#   "split"        - te_*_split variants (for pathway/gene analysis)
#   "rna"          - kics_rna, lfs_rna, stjude_rna, matched_dna_rna, lfs_wgs2rna
#   "ancestry"     - ancestry data
#   "location"     - location_100kb_g (genomic windows)
#   "clinical"     - clinical, metrics
#   "genes"        - cpg, gene_size, hg37_genes, num_calls
#   "common"       - common TE data (if PROCESS_COMMON_TES=TRUE)
#   "all"          - load everything (default if REQUIRED_DATA not set)

# Processing options
PROCESS_COMMON_TES <- FALSE # Set to TRUE to visualize common TEs
COMMON_MODE <- FALSE  # Set TRUE to use common TEs (no frequency filtering), outputs to /common/ subdir
FULLLENGTH_YOUNG_MODE <- FALSE  # Set TRUE to use full-length young TEs (L1HS ≥5900bp, AluY, SVA_E/F), outputs to /fulllength_young/ subdir
TEST_HOSTSEQ_SPLITS <- FALSE  # Set to TRUE to test different HostSeq filter/analysis split percentages

# Helper function to check if data category is needed
needs_data <- function(category) {
  if (!exists("REQUIRED_DATA", envir = .GlobalEnv)) return(TRUE)
  req <- get("REQUIRED_DATA", envir = .GlobalEnv)
  if ("all" %in% req) return(TRUE)
  category %in% req
}

# Check if data is already loaded (skip reloading when running via ALL script)
# Note: Check for te_aff which is needed by count_matrix category
if (exists("te_all") && exists("te_aff") && exists("te_aff_expand") && exists("te_all_split")) {
  cat("Data already loaded, skipping reload...\n")
  return(invisible(NULL))
}

cat("Loading germline data...\n")

#### SETUP PATHS ####
# Base directory for plots
plot_dir_base <- "/Users/briannelaverty/Library/Mobile Documents/com~apple~CloudDocs/Malkin_lab/figures/te/graphs/germline/"

# Use subdirectory based on analysis mode (fulllength_young > common > base)
if (exists("FULLLENGTH_YOUNG_MODE") && FULLLENGTH_YOUNG_MODE) {
  plot_dir <- paste0(plot_dir_base, "fulllength_young/")
  cat("*** FULLLENGTH YOUNG MODE ACTIVE: Using full-length young TEs, output to /fulllength_young/ ***\n")
} else if (exists("COMMON_MODE") && COMMON_MODE) {
  plot_dir <- paste0(plot_dir_base, "common/")
  cat("*** COMMON MODE ACTIVE: Using unfiltered TEs, output to /common/ ***\n")
} else {
  plot_dir <- plot_dir_base
}
r_dir_files <- paste0(plot_dir, "files/")

# Set data directory based on TEST_MODE
if (TEST_MODE) {
  r_dir <- "/Users/briannelaverty/Documents/R_Malkin/te/data/R_obj/test_data/germline/"
  cat("*** TEST MODE ACTIVE: Using test data from", r_dir, "***\n\n")
} else {
  r_dir <- "/Users/briannelaverty/Documents/R_Malkin/te/data/R_obj/germline/"
}


# Create directories if they don't exist
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
if (!dir.exists(r_dir_files)) dir.create(r_dir_files, recursive = TRUE)

# Create subdirectories for organized plot storage
subdirs <- c("dataset", "general", "counts_cohort", "counts_tumour_type",
             "counts_clinical_lfs", "counts_clinical_kics", "ancestry",
             "pca_umap", "cancer_genes", "pathway", "reg_element", "rna", "taylor", "other")
for (subdir in subdirs) {
  subdir_path <- paste0(plot_dir, subdir, "/")
  if (!dir.exists(subdir_path)) dir.create(subdir_path, recursive = TRUE)
}

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
if (needs_data("genes")) {
  cat("Loading gene data files...\n")
  num_calls <- read.delim("/Users/briannelaverty/Documents/R_Malkin/te/data/final/count_calls.txt", sep="\t", header=TRUE)
  hg37_genes <- read.delim("/Users/briannelaverty/Documents/R_Malkin/te/data/final/hg37_genes.tsv", sep="\t", header=TRUE)
  cpg <- read.csv("/Users/briannelaverty/Documents/R_Malkin/te/data/raw/kics_cpg.txt", header=FALSE, sep="\t") # cancer predisposition
}

if (needs_data("location")) {
  cat("Loading location data...\n")
  # Location windows with counts of LINE, ALU, SVA with 0 columns removed with clinical info
  # Made by saving te_aff_split then running /hpf/largeprojects/davidm/blaverty/te/ml/scripts/00_ml_pipeline.sh
  location_100kb_g <- read.csv("/Users/briannelaverty/Documents/R_Malkin/te/data/final/100kb_aff_unaff_g.csv", stringsAsFactors = FALSE)
}

#### LOAD RNA DATA ####
if (needs_data("rna")) {
  cat("Loading RNA data files...\n")
  kics_rna <- read.delim("/Users/briannelaverty/Documents/R_Malkin/te/data/rna/kics_FPKM.tsv", sep="\t", header=TRUE, check.names = FALSE)
  lfs_rna <- read.delim("/Users/briannelaverty/Documents/R_Malkin/te/data/rna/lfs_FPKM.tsv", sep="\t", header=TRUE, check.names = FALSE)
  lfs_wgs2rna <- read.delim("/Users/briannelaverty/Documents/R_Malkin/te/data/rna/lfs_wgs2rna.csv", sep=",", header=TRUE)
  stjude_rna <- read.delim("/Users/briannelaverty/Documents/R_Malkin/te/data/rna/stjude_FPKM.tsv", sep="\t", header=TRUE, check.names = FALSE)
  matched_dna_rna <- read.csv("/Users/briannelaverty/Documents/R_Malkin/te/data/rna/matched_dna_rna.csv", stringsAsFactors = FALSE)
}

#### LOAD SAVED PROCESSED DATA ####
cat("Loading saved processed data objects...\n")

# Load prepared data objects (use production dir for metadata in test mode)
metadata_dir <- if (TEST_MODE) "/Users/briannelaverty/Documents/R_Malkin/te/data/R_obj/germline/" else r_dir

if (needs_data("ancestry")) {
  load(paste0(metadata_dir, "ancestry.RData"))
}

if (needs_data("genes")) {
  load(paste0(metadata_dir, "gene_size.RData"))
}

# Always load these core metadata files (needed by most scripts)
load(paste0(metadata_dir, "chr_lengths.RData"))
load(paste0(metadata_dir, "nohits.RData"))
load(paste0(metadata_dir, "hostseq_cancer.RData"))

if (needs_data("clinical")) {
  load(paste0(metadata_dir, "clinical.RData"))
  load(paste0(metadata_dir, "metrics.RData"))
}

if (needs_data("count_matrix") || needs_data("expand") || needs_data("split")) {
  load(paste0(r_dir, "te_raw_prepped.RData"))
}

#### CREATE DATA FRAMES FOR VISUALIZATION ####
cat("Creating data frames for visualization...\n")

# TE DATA - COUNT MATRIX FORMAT
if (needs_data("count_matrix")) {
  # FORMAT: One row per sample, columns for TE type counts (LINE1, ALU, SVA, etc.)
  # USE CASE: Statistical comparisons of total TE counts between groups

  if (FULLLENGTH_YOUNG_MODE) {
    cat("Loading FULLLENGTH YOUNG TE count matrix data...\n")
    load(paste0(r_dir, "final_te_count_fulllength_young.RData"))
    data_source <- final_te_count_fulllength_young
  } else if (COMMON_MODE) {
    cat("Loading COMMON TE count matrix data...\n")
    load(paste0(r_dir, "final_te_count_common.RData"))
    data_source <- final_te_count_common
  } else {
    cat("Loading rare TE count matrix data...\n")
    load(paste0(r_dir, "final_te_count_rare.RData"))
    data_source <- if (TEST_MODE) final_te_count_test else final_te_count
  }

  te_all <- data_source$te_all
  te_aff <- data_source$te_aff
  te_lfs <- data_source$te_lfs
  te_lfs_mut_wt <- data_source$te_lfs_mut_wt
  te_kics <- data_source$te_kics
  te_taylor <- data_source$te_taylor
  te_hostseq <- data_source$te_hostseq

  # Add tumor_type_grouped to each dataframe based on its own distribution
  te_all <- add_tumor_type_grouped(te_all)
  te_aff <- add_tumor_type_grouped(te_aff)
  te_lfs <- add_tumor_type_grouped(te_lfs)
  te_lfs_mut_wt <- add_tumor_type_grouped(te_lfs_mut_wt)
  te_kics <- add_tumor_type_grouped(te_kics)
  te_taylor <- add_tumor_type_grouped(te_taylor)
  te_hostseq <- add_tumor_type_grouped(te_hostseq)
}

# TE DATA - EXPANDED FORMAT
if (needs_data("expand")) {
  # FORMAT: One row per complete TE insertion
  # USE CASE: General TE analysis, counting TEs per sample, TE characteristics
  # GENE INFO: Contains gene overlap information but TE is not split by genes

  if (FULLLENGTH_YOUNG_MODE) {
    cat("Loading FULLLENGTH YOUNG TE expanded data...\n")
    load(paste0(r_dir, "final_te_count_expand_fulllength_young.RData"))
    data_source <- final_te_count_expand_fulllength_young
  } else if (COMMON_MODE) {
    cat("Loading COMMON TE expanded data...\n")
    load(paste0(r_dir, "final_te_count_expand_common.RData"))
    data_source <- final_te_count_expand_common
  } else {
    cat("Loading rare TE expanded data...\n")
    load(paste0(r_dir, "final_te_count_expand_rare.RData"))
    data_source <- if (TEST_MODE) final_te_count_expand_test else final_te_count_expand
  }

  te_all_expand <- data_source$te_all
  te_aff_expand <- data_source$te_aff
  te_lfs_expand <- data_source$te_lfs
  te_lfs_mut_wt_expand <- data_source$te_lfs_mut_wt
  te_kics_expand <- data_source$te_kics
  te_taylor_expand <- data_source$te_taylor
  te_hostseq_expand <- data_source$te_hostseq

  # Add tumor_type_grouped to each dataframe based on its own distribution
  te_all_expand <- add_tumor_type_grouped(te_all_expand)
  te_aff_expand <- add_tumor_type_grouped(te_aff_expand)
  te_lfs_expand <- add_tumor_type_grouped(te_lfs_expand)
  te_lfs_mut_wt_expand <- add_tumor_type_grouped(te_lfs_mut_wt_expand)
  te_kics_expand <- add_tumor_type_grouped(te_kics_expand)
  te_taylor_expand <- add_tumor_type_grouped(te_taylor_expand)
  te_hostseq_expand <- add_tumor_type_grouped(te_hostseq_expand)

  # Create affected + unaffected dataset (excluding Taylor and HostSeq)
  te_aff_unaff_expand <- te_all_expand %>% filter(cohort != "Taylor" & cohort != "HostSeq")

  # Create KICS + HostSeq expand dataset (for RE analysis)
  # Manually combine since te_all_expand doesn't include KICS
  te_kics_hostseq_expand <- rbind(te_kics_expand, te_hostseq_expand)
}

# TE DATA - SPLIT BY GENE FORMAT
if (needs_data("split")) {
  # FORMAT: One row per TE-gene overlap (a single TE can have multiple rows if it overlaps multiple genes)
  # USE CASE: Gene-level and pathway enrichment analysis
  # GENE INFO: Each row represents one TE insertion affecting one specific gene
  # DIFFERENCE FROM EXPAND: If a TE overlaps 3 genes, it appears as 3 separate rows (vs 1 row in expand format)
  # SOURCE: Uses AnnotSV split mode input which pre-identifies gene overlaps

  if (FULLLENGTH_YOUNG_MODE) {
    cat("Loading FULLLENGTH YOUNG TE split (gene) data...\n")
    load(paste0(r_dir, "final_te_count_split_fulllength_young.RData"))
    data_source <- final_te_count_split_fulllength_young
  } else if (COMMON_MODE) {
    cat("Loading COMMON TE split (gene) data...\n")
    load(paste0(r_dir, "final_te_count_split_common.RData"))
    data_source <- final_te_count_split_common
  } else {
    cat("Loading rare TE split (gene) data...\n")
    load(paste0(r_dir, "final_te_count_expand_split.RData"))
    data_source <- if (TEST_MODE) final_te_count_expand_split_test else final_te_count_expand_split
  }

  te_all_split <- data_source$te_all
  te_aff_split <- data_source$te_aff
  te_aff_unaff_split <- data_source$te_aff_unaff
  te_lfs_split <- data_source$te_lfs
  te_lfs_mut_wt_split <- data_source$te_lfs_mut_wt
  te_kics_split <- data_source$te_kics
  te_taylor_split <- data_source$te_taylor
  te_hostseq_split <- data_source$te_hostseq
  te_kics_hostseq <- data_source$te_kics_hostseq

  # Add tumor_type_grouped to each dataframe based on its own distribution
  te_all_split <- add_tumor_type_grouped(te_all_split)
  te_aff_split <- add_tumor_type_grouped(te_aff_split)
  te_aff_unaff_split <- add_tumor_type_grouped(te_aff_unaff_split)
  te_lfs_split <- add_tumor_type_grouped(te_lfs_split)
  te_lfs_mut_wt_split <- add_tumor_type_grouped(te_lfs_mut_wt_split)
  te_kics_split <- add_tumor_type_grouped(te_kics_split)
  te_taylor_split <- add_tumor_type_grouped(te_taylor_split)
  te_hostseq_split <- add_tumor_type_grouped(te_hostseq_split)
  te_kics_hostseq <- add_tumor_type_grouped(te_kics_hostseq)

  # Create aliases for gene analysis (backward compatibility)
  te_all_split_genes <- te_all_split
  te_aff_split_genes <- te_aff_split
  te_lfs_split_genes <- te_lfs_split
  te_kics_split_genes <- te_kics_split
  te_taylor_split_genes <- te_taylor_split
  te_hostseq_split_genes <- te_hostseq_split
}

# COMMON TE DATA
if (PROCESS_COMMON_TES && needs_data("common")) {
  cat("Loading common TE data...\n")
  load(paste0(r_dir, "final_te_count_common.RData"))
  te_aff_common <- final_te_count_common$te_aff
  te_lfs_common <- final_te_count_common$te_lfs
  te_kics_common <- final_te_count_common$te_kics
  te_all_common <- final_te_count_common$te_all

  # Create HostSeq and Taylor subsets from te_all_common
  te_hostseq_common <- te_all_common %>% filter(cohort == "HostSeq")
  te_taylor_common <- te_all_common %>% filter(cohort == "Taylor")

  # Common TE expanded format
  load(paste0(r_dir, "final_te_count_expand_common.RData"))
  te_aff_expand_common <- final_te_count_expand_common$te_aff
  te_lfs_expand_common <- final_te_count_expand_common$te_lfs
  te_kics_expand_common <- final_te_count_expand_common$te_kics
  te_all_expand_common <- final_te_count_expand_common$te_all

  # Common TE split by gene format
  load(paste0(r_dir, "final_te_count_split_common.RData"))
  te_aff_split_common <- final_te_count_split_common$te_aff
  te_lfs_split_common <- final_te_count_split_common$te_lfs
  te_kics_split_common <- final_te_count_split_common$te_kics
  te_all_split_common <- final_te_count_split_common$te_all
} else if (PROCESS_COMMON_TES) {
  cat("Skipping common TE data (not in REQUIRED_DATA)\n")
} else {
  cat("Skipping common TE data (PROCESS_COMMON_TES = FALSE)\n")
}

# Print summary of what was loaded
cat("✓ Germline data loaded successfully\n")
if (exists("REQUIRED_DATA", envir = .GlobalEnv)) {
  req <- get("REQUIRED_DATA", envir = .GlobalEnv)
  cat("  Loaded categories:", paste(req, collapse = ", "), "\n")
}
cat("\n")
cat("========================================\n")
cat("  GERMLINE DATASET SUMMARY\n")
cat("========================================\n")
cat("Cohorts: LFS, KICS, Taylor, HostSeq (analysis group)\n\n")
cat("Available datasets by sample selection:\n\n")
cat("te_all_*         : All germline samples (LFS + KICS + Taylor + HostSeq, affected + unaffected)\n")
cat("te_aff_*         : Affected only (LFS + KICS, excluding Taylor and HostSeq)\n")
cat("te_lfs_*         : All TP53 Mutant samples (from any cohort)\n")
cat("te_lfs_mut_wt_*  : All TP53 Mutant samples + LFS_wt cohort\n")
cat("                   Includes affected + unaffected, with Cancer status indicator\n")
cat("te_kics_*        : KICS cohort only (excluding HostSeq, affected)\n")
cat("te_taylor_*      : Taylor cohort only\n")
cat("te_hostseq_*     : HostSeq analysis group only (cohort == 'HostSeq')\n")
cat("te_aff_unaff_*   : Affected + unaffected (LFS + KICS only, excluding Taylor and HostSeq)\n")
cat("\nNote: HostSeq samples only appear in te_all_* and te_hostseq_* datasets\n")
cat("\nData formats:\n")
cat("- *_expand       : One row per TE insertion (for TE-level analysis)\n")
cat("- *_split        : One row per TE-gene overlap (for gene/pathway analysis)\n")
cat("- (no suffix)    : Count matrix format (one row per sample)\n")
cat("========================================\n\n")
