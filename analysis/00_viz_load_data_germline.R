#!/usr/bin/env Rscript

# Load Germline Data for TE Visualization
# This script loads all necessary data for germline TE visualization
# Source this after 00_viz_common_setup.R

# Check if data is already loaded (skip reloading when running via ALL script)
if (exists("te_all") && exists("te_aff_expand") && exists("te_all_split")) {
  cat("Data already loaded, skipping reload...\n")
  return(invisible(NULL))
}

cat("Loading germline data...\n")

#### SETUP PATHS ####
plot_dir <- "/Users/briannelaverty/Library/Mobile Documents/com~apple~CloudDocs/Malkin_lab/figures/te/graphs/germline/"
r_dir_files <- "/Users/briannelaverty/Library/Mobile Documents/com~apple~CloudDocs/Malkin_lab/figures/te/graphs/germline/files/"

# Set data directory based on TEST_MODE
if (TEST_MODE) {
  r_dir <- "/Users/briannelaverty/Documents/R_Malkin/te/data/R_obj/test_data/germline/"
  cat("*** TEST MODE ACTIVE: Using test data from", r_dir, "***\n\n")
} else {
  r_dir <- "/Users/briannelaverty/Documents/R_Malkin/te/data/R_obj/germline/"
}

# Processing options
PROCESS_COMMON_TES <- FALSE  # Set to TRUE to visualize common TEs
TEST_HOSTSEQ_SPLITS <- FALSE  # Set to TRUE to test different HostSeq filter/analysis split percentages

# Create directories if they don't exist
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
if (!dir.exists(r_dir_files)) dir.create(r_dir_files, recursive = TRUE)

# Create subdirectories for organized plot storage
subdirs <- c("dataset", "general", "counts_cohort", "counts_tumour_type",
             "counts_clinical_lfs", "counts_clinical_kics", "ancestry",
             "pca_umap", "cancer_genes", "pathway", "reg_element", "rna", "other")
for (subdir in subdirs) {
  subdir_path <- paste0(plot_dir, subdir, "/")
  if (!dir.exists(subdir_path)) dir.create(subdir_path, recursive = TRUE)
}

# Redirect all stdout to text file
stdout_file <- paste0(r_dir_files, "graph_output_text.txt")
if (file.exists(stdout_file)) file.remove(stdout_file)
sink(stdout_file, split = TRUE)

# Open a PDF device at the start to capture all plots in a single PDF
pdf(file = paste0(plot_dir, "graph_output.pdf"), width = 12, height = 8)

#### LOAD ADDITIONAL DATA ####
cat("Loading additional raw data files...\n")
num_calls <- read.delim("/Users/briannelaverty/Documents/R_Malkin/te/data/final/count_calls.txt", sep="\t", header=TRUE)
hg37_genes <- read.delim("/Users/briannelaverty/Documents/R_Malkin/te/data/final/hg37_genes.tsv", sep="\t", header=TRUE)
cpg<- read.csv("/Users/briannelaverty/Documents/R_Malkin/te/data/raw/kics_cpg.txt", header=FALSE, sep="\t") # cancer predisposition
# Location windows with counts of LINE, ALU, SVA with 0 columns removed with clinical info
# Made by saving te_aff_split then running /hpf/largeprojects/davidm/blaverty/te/ml/scripts/00_ml_pipeline.sh
location_100kb_g <- read.csv("/Users/briannelaverty/Documents/R_Malkin/te/data/final/100kb_aff_unaff_g.csv", stringsAsFactors = FALSE)

#### LOAD RNA DATA ####
cat("Loading RNA data files...\n")
kics_rna <- read.delim("/Users/briannelaverty/Documents/R_Malkin/te/data/rna/kics_FPKM.tsv", sep="\t", header=TRUE, check.names = FALSE)
lfs_rna <- read.delim("/Users/briannelaverty/Documents/R_Malkin/te/data/rna/lfs_FPKM.tsv", sep="\t", header=TRUE, check.names = FALSE)
lfs_wgs2rna <- read.delim("/Users/briannelaverty/Documents/R_Malkin/te/data/rna/lfs_wgs2rna.csv", sep=",", header=TRUE)
stjude_rna <- read.delim("/Users/briannelaverty/Documents/R_Malkin/te/data/rna/stjude_FPKM.tsv", sep="\t", header=TRUE, check.names = FALSE)
matched_dna_rna <- read.csv("/Users/briannelaverty/Documents/R_Malkin/te/data/rna/matched_dna_rna.csv", stringsAsFactors = FALSE)

#### LOAD SAVED PROCESSED DATA ####
cat("Loading saved processed data objects...\n")

# Load prepared data objects (use production dir for metadata in test mode)
metadata_dir <- if (TEST_MODE) "/Users/briannelaverty/Documents/R_Malkin/te/data/R_obj/germline/" else r_dir

load(paste0(metadata_dir, "ancestry.RData"))
load(paste0(metadata_dir, "gene_size.RData"))
load(paste0(metadata_dir, "chr_lengths.RData"))
load(paste0(metadata_dir, "nohits.RData"))
load(paste0(metadata_dir, "hostseq_cancer.RData"))
load(paste0(metadata_dir, "clinical.RData"))
load(paste0(metadata_dir, "metrics.RData"))
load(paste0(r_dir, "te_raw_prepped.RData"))

#### CREATE DATA FRAMES FOR VISUALIZATION ####
cat("Creating data frames for visualization...\n")

# RARE TE DATA
cat("Loading rare TE data...\n")

# Rare TE count matrix format
# FORMAT: One row per sample, columns for TE type counts (LINE1, ALU, SVA, etc.)
# USE CASE: Statistical comparisons of total TE counts between groups
load(paste0(r_dir, "final_te_count_rare.RData"))

# Handle different variable names in test vs production data
if (TEST_MODE) {
  te_all <- final_te_count_test$te_all
  te_aff <- final_te_count_test$te_aff
  te_lfs <- final_te_count_test$te_lfs
  te_kics <- final_te_count_test$te_kics
  te_taylor <- final_te_count_test$te_taylor
  te_hostseq <- final_te_count_test$te_hostseq
} else {
  te_all <- final_te_count$te_all
  te_aff <- final_te_count$te_aff
  te_lfs <- final_te_count$te_lfs
  te_kics <- final_te_count$te_kics
  te_taylor <- final_te_count$te_taylor
  te_hostseq <- final_te_count$te_hostseq
}

# Rare TE expanded format
# FORMAT: One row per complete TE insertion
# USE CASE: General TE analysis, counting TEs per sample, TE characteristics
# GENE INFO: Contains gene overlap information but TE is not split by genes
load(paste0(r_dir, "final_te_count_expand_rare.RData"))

if (TEST_MODE) {
  te_all_expand <- final_te_count_expand_test$te_all
  te_aff_expand <- final_te_count_expand_test$te_aff
  te_lfs_expand <- final_te_count_expand_test$te_lfs
  te_kics_expand <- final_te_count_expand_test$te_kics
  te_taylor_expand <- final_te_count_expand_test$te_taylor
  te_hostseq_expand <- final_te_count_expand_test$te_hostseq
} else {
  te_all_expand <- final_te_count_expand$te_all
  te_aff_expand <- final_te_count_expand$te_aff
  te_lfs_expand <- final_te_count_expand$te_lfs
  te_kics_expand <- final_te_count_expand$te_kics
  te_taylor_expand <- final_te_count_expand$te_taylor
  te_hostseq_expand <- final_te_count_expand$te_hostseq
}

# Create affected + unaffected dataset (excluding Taylor and HostSeq)
te_aff_unaff_expand <- te_all_expand %>% filter(cohort != "Taylor" & cohort != "HostSeq")

# Create KICS + HostSeq expand dataset (for RE analysis)
# Manually combine since te_all_expand doesn't include KICS
te_kics_hostseq_expand <- rbind(te_kics_expand, te_hostseq_expand)

# Rare TE split by gene format (annotSV split mode)
# FORMAT: One row per TE-gene overlap (a single TE can have multiple rows if it overlaps multiple genes)
# USE CASE: Gene-level and pathway enrichment analysis
# GENE INFO: Each row represents one TE insertion affecting one specific gene
# DIFFERENCE FROM EXPAND: If a TE overlaps 3 genes, it appears as 3 separate rows (vs 1 row in expand format)
# SOURCE: Uses AnnotSV split mode input which pre-identifies gene overlaps
load(paste0(r_dir, "final_te_count_expand_split.RData"))

if (TEST_MODE) {
  te_all_split <- final_te_count_expand_split_test$te_all
  te_aff_split <- final_te_count_expand_split_test$te_aff
  te_lfs_split <- final_te_count_expand_split_test$te_lfs
  te_kics_split <- final_te_count_expand_split_test$te_kics
  te_taylor_split <- final_te_count_expand_split_test$te_taylor
  te_hostseq_split <- final_te_count_expand_split_test$te_hostseq
  te_kics_hostseq <- final_te_count_expand_split_test$te_kics_hostseq
} else {
  te_all_split <- final_te_count_expand_split$te_all
  te_aff_split <- final_te_count_expand_split$te_aff
  te_lfs_split <- final_te_count_expand_split$te_lfs
  te_kics_split <- final_te_count_expand_split$te_kics
  te_taylor_split <- final_te_count_expand_split$te_taylor
  te_hostseq_split <- final_te_count_expand_split$te_hostseq
  te_kics_hostseq <- final_te_count_expand_split$te_kics_hostseq
}

# COMMON TE DATA
if (PROCESS_COMMON_TES) {
  cat("Loading common TE data...\n")
  load(paste0(r_dir, "final_te_count_common.RData"))
  te_aff_common <- final_te_count_common$te_aff
  te_lfs_common <- final_te_count_common$te_lfs
  te_kics_common <- final_te_count_common$te_kics
  te_all_common <- final_te_count_common$te_all

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
} else {
  cat("Skipping common TE data (PROCESS_COMMON_TES = FALSE)\n")
}

# SPLIT TE DATA FOR GENE ANALYSIS
# Now loaded above with te_*_split variables - keeping separate names for backward compatibility
if (TEST_MODE) {
  te_all_split_genes <- final_te_count_expand_split_test$te_all
  te_aff_split_genes <- final_te_count_expand_split_test$te_aff
  te_lfs_split_genes <- final_te_count_expand_split_test$te_lfs
  te_kics_split_genes <- final_te_count_expand_split_test$te_kics
  te_taylor_split_genes <- final_te_count_expand_split_test$te_taylor
  te_hostseq_split_genes <- final_te_count_expand_split_test$te_hostseq
} else {
  te_all_split_genes <- final_te_count_expand_split$te_all
  te_aff_split_genes <- final_te_count_expand_split$te_aff
  te_lfs_split_genes <- final_te_count_expand_split$te_lfs
  te_kics_split_genes <- final_te_count_expand_split$te_kics
  te_taylor_split_genes <- final_te_count_expand_split$te_taylor
  te_hostseq_split_genes <- final_te_count_expand_split$te_hostseq
}

cat("✓ Germline data loaded successfully\n")
cat("\n")
cat("========================================\n")
cat("  GERMLINE DATASET SUMMARY\n")
cat("========================================\n")
cat("Cohorts: LFS, KICS, Taylor, HostSeq (analysis group)\n\n")
cat("Available datasets by sample selection:\n\n")
cat("te_all_*         : All germline samples (LFS + KICS + Taylor + HostSeq, affected + unaffected)\n")
cat("te_aff_*         : Affected only (LFS + KICS, excluding Taylor and HostSeq)\n")
cat("te_lfs_*         : LFS-related samples (TP53 mutant OR cohort in LFS/Nick/SJ)\n")
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
