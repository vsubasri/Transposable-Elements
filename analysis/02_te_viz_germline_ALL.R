#!/usr/bin/env Rscript

# Germline TE Visualization - Run All Analyses
# This script runs all germline visualization analyses in sequence

cat("======================================\n")
cat("  GERMLINE TE VISUALIZATION - ALL\n")
cat("======================================\n\n")

#### MODULE SELECTION ####
# Set to FALSE to skip a module
RUN_DESCRIPTIVE <- FALSE
RUN_CLINICAL <- FALSE
RUN_SPECIFIC_TES <- FALSE
RUN_CANCER_GENES <- FALSE
RUN_PATHWAY <- FALSE
RUN_RNA <- TRUE
RUN_RE <- TRUE
RUN_PCA_UMAP <- FALSE
RUN_METHYLATION <- FALSE
RUN_TAYLOR <- TRUE 

# Source common setup
source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/00_viz_common_setup.R")

# Load germline data
source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/00_viz_load_data_germline.R")

cat("\n=== Running Analysis Scripts ===\n\n")

# 01 - Descriptive statistics
if (RUN_DESCRIPTIVE) {
  cat("1/10: Descriptive statistics...\n")
  source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/02_te_viz_germline_01_descriptive.R")
} else {
  cat("1/10: Skipping descriptive statistics\n")
}

# 02 - Clinical associations
if (RUN_CLINICAL) {
  cat("2/10: Clinical associations...\n")
  source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/02_te_viz_germline_02_clinical.R")
} else {
  cat("2/10: Skipping clinical associations\n")
}

# 03 - Specific TEs
if (RUN_SPECIFIC_TES) {
  cat("3/10: Specific TEs analysis...\n")
  source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/02_te_viz_germline_03_specific_tes.R")
} else {
  cat("3/10: Skipping specific TEs analysis\n")
}

# 04 - Cancer genes
if (RUN_CANCER_GENES) {
  cat("4/10: Cancer genes analysis...\n")
  source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/02_te_viz_germline_04_cancer_genes.R")
} else {
  cat("4/10: Skipping cancer genes analysis\n")
}

# 05 - Pathway analysis
if (RUN_PATHWAY) {
  cat("5/10: Pathway analysis...\n")
  source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/02_te_viz_germline_05_pathway.R")
} else {
  cat("5/10: Skipping pathway analysis\n")
}

# 06 - RNA expression
if (RUN_RNA) {
  cat("6/10: RNA expression analysis...\n")
  source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/02_te_viz_germline_06_rna.R")
} else {
  cat("6/10: Skipping RNA expression analysis\n")
}

# 07 - Regulatory elements
if (RUN_RE) {
  cat("7/10: Regulatory elements analysis...\n")
  source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/02_te_viz_germline_07_re.R")
} else {
  cat("7/10: Skipping regulatory elements analysis\n")
}

# 08 - PCA/UMAP
if (RUN_PCA_UMAP) {
  cat("8/10: PCA/UMAP analysis...\n")
  source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/02_te_viz_germline_08_pca_umap.R")
} else {
  cat("8/10: Skipping PCA/UMAP analysis\n")
}

# 09 - Methylation
if (RUN_METHYLATION) {
  cat("9/10: Methylation analysis...\n")
  source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/02_te_viz_germline_09_methylation.R")
} else {
  cat("9/10: Skipping methylation analysis\n")
}

# 10 - Taylor cohort
if (RUN_TAYLOR) {
  cat("10/10: Taylor cohort analysis...\n")
  source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/02_te_viz_germline_10_taylor.R")
} else {
  cat("10/10: Skipping Taylor cohort analysis\n")
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
