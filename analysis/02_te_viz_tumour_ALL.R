#!/usr/bin/env Rscript

# Tumour TE Visualization - Run All Analyses
# This script runs all tumour visualization analyses in sequence

cat("======================================\n")
cat("  TUMOUR TE VISUALIZATION - ALL\n")
cat("======================================\n\n")

#### MODULE SELECTION ####
# Set to FALSE to skip a module
RUN_DESCRIPTIVE <- FALSE
RUN_CLINICAL <- FALSE 
RUN_SPECIFIC_TES <- FALSE 
RUN_LOH <- FALSE 
RUN_SURVIVAL <- FALSE 
RUN_PATHWAY <- FALSE 
RUN_RNA <- TRUE
RUN_RE <- TRUE
RUN_PCA_UMAP <- TRUE
RUN_METHYLATION <- TRUE
RUN_SV_OVERLAP <- TRUE
RUN_TAYLOR <- TRUE

# Source common setup
source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/00_viz_common_setup.R")

# Load tumour data
source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/00_viz_load_data_tumour.R")

cat("\n=== Running Analysis Scripts ===\n\n")

# 01 - Descriptive statistics
if (RUN_DESCRIPTIVE) {
  cat("1/12: Descriptive statistics...\n")
  source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/02_te_viz_tumour_01_descriptive.R")
} else {
  cat("1/12: Skipping descriptive statistics\n")
}

# 02 - Clinical associations
if (RUN_CLINICAL) {
  cat("2/12: Clinical associations...\n")
  source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/02_te_viz_tumour_02_clinical.R")
} else {
  cat("2/12: Skipping clinical associations\n")
}

# 03 - Specific TEs
if (RUN_SPECIFIC_TES) {
  cat("3/12: Specific TEs analysis...\n")
  source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/02_te_viz_tumour_03_specific_tes.R")
} else {
  cat("3/12: Skipping specific TEs analysis\n")
}

# 04 - LOH timing
if (RUN_LOH) {
  cat("4/12: LOH timing analysis...\n")
  source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/02_te_viz_tumour_04_loh.R")
} else {
  cat("4/12: Skipping LOH timing analysis\n")
}

# 05 - Survival analysis
if (RUN_SURVIVAL) {
  cat("5/12: Survival analysis...\n")
  source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/02_te_viz_tumour_05_survival.R")
} else {
  cat("5/12: Skipping survival analysis\n")
}

# 06 - Pathway analysis
if (RUN_PATHWAY) {
  cat("6/12: Pathway analysis...\n")
  source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/02_te_viz_tumour_06_pathway.R")
} else {
  cat("6/12: Skipping pathway analysis\n")
}

# 07 - RNA expression
if (RUN_RNA) {
  cat("7/12: RNA expression analysis...\n")
  source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/02_te_viz_tumour_07_rna.R")
} else {
  cat("7/12: Skipping RNA expression analysis\n")
}

# 08 - Regulatory elements
if (RUN_RE) {
  cat("8/12: Regulatory elements analysis...\n")
  source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/02_te_viz_tumour_08_re.R")
} else {
  cat("8/12: Skipping regulatory elements analysis\n")
}

# 09 - PCA/UMAP
if (RUN_PCA_UMAP) {
  cat("9/12: PCA/UMAP analysis...\n")
  source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/02_te_viz_tumour_09_pca_umap.R")
} else {
  cat("9/12: Skipping PCA/UMAP analysis\n")
}

# 10 - Methylation
if (RUN_METHYLATION) {
  cat("10/12: Methylation analysis...\n")
  source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/02_te_viz_tumour_10_methylation.R")
} else {
  cat("10/12: Skipping methylation analysis\n")
}

# 11 - SV overlap
if (RUN_SV_OVERLAP) {
  cat("11/12: SV overlap analysis...\n")
  source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/02_te_viz_tumour_11_sv_overlap.R")
} else {
  cat("11/12: Skipping SV overlap analysis\n")
}

# 12 - Taylor cohort
if (RUN_TAYLOR) {
  cat("12/12: Taylor cohort analysis...\n")
  source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/02_te_viz_tumour_12_taylor.R")
} else {
  cat("12/12: Skipping Taylor cohort analysis\n")
}

# Close PDF device
dev.off()

# Close sink
sink()

cat("\n======================================\n")
cat("  ALL ANALYSES COMPLETED\n")
cat("======================================\n")
cat("Output saved to:", plot_dir, "\n")
cat("PDF: graph_output.pdf\n")
cat("Text log: files/graph_output_text.txt\n")
