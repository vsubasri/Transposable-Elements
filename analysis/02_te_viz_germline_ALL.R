#!/usr/bin/env Rscript

# Germline TE Visualization - Run All Analyses
# This script runs all germline visualization analyses in sequence

cat("======================================\n")
cat("  GERMLINE TE VISUALIZATION - ALL\n")
cat("======================================\n\n")

#### TEST MODE ####
# Enable test mode BEFORE sourcing common setup
TEST_MODE <- FALSE 

#### MODULE SELECTION ####
# Set to FALSE to skip a module
RUN_DESCRIPTIVE <- TRUE
RUN_CLINICAL <- TRUE
RUN_SPECIFIC_TES <- TRUE
RUN_CANCER_GENES <- TRUE
RUN_PATHWAY <- TRUE
RUN_RNA <- TRUE
RUN_RE <- TRUE
RUN_RE_RNA <- TRUE
RUN_PCA_UMAP <- TRUE
RUN_METHYLATION <- TRUE
RUN_TAYLOR <- TRUE
RUN_SUBFAMILY <- TRUE
RUN_FULLLENGTH_YOUNG_SOURCE <- TRUE

# Source common setup
source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/00_viz_common_setup.R")

# Load ALL germline data (master script needs everything)
REQUIRED_DATA <- c("all")
source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/00_viz_load_data_germline.R")

cat("\n=== Running Analysis Scripts ===\n\n")

# 01 - Descriptive statistics
if (RUN_DESCRIPTIVE) {
  cat("1/13: Descriptive statistics...\n")
  source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/02_te_viz_germline_01_descriptive.R")
} else {
  cat("1/13: Skipping descriptive statistics\n")
}

# 02 - Clinical associations
if (RUN_CLINICAL) {
  cat("2/13: Clinical associations...\n")
  source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/02_te_viz_germline_02_clinical.R")
} else {
  cat("2/13: Skipping clinical associations\n")
}

# 03 - Specific TEs
if (RUN_SPECIFIC_TES) {
  cat("3/13: Specific TEs analysis...\n")
  source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/02_te_viz_germline_03_specific_tes.R")
} else {
  cat("3/13: Skipping specific TEs analysis\n")
}

# 04 - Cancer genes
if (RUN_CANCER_GENES) {
  cat("4/13: Cancer genes analysis...\n")
  source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/02_te_viz_germline_04_cancer_genes.R")
} else {
  cat("4/13: Skipping cancer genes analysis\n")
}

# 05 - Pathway analysis
if (RUN_PATHWAY) {
  cat("5/13: Pathway analysis...\n")
  source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/02_te_viz_germline_05_pathway.R")
} else {
  cat("5/13: Skipping pathway analysis\n")
}

# 06 - RNA expression
if (RUN_RNA) {
  cat("6/13: RNA expression analysis...\n")
  source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/02_te_viz_germline_06_rna.R")
} else {
  cat("6/13: Skipping RNA expression analysis\n")
}

# 07 - Regulatory elements (ORA)
if (RUN_RE) {
  cat("7/13: Regulatory elements analysis...\n")
  source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/02_te_viz_germline_07_re.R")
} else {
  cat("7/13: Skipping regulatory elements analysis\n")
}

# 08 - RE-RNA differential expression (GSEA) - NEW
if (RUN_RE_RNA) {
  cat("8/13: RE-RNA differential expression analysis...\n")
  source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/02_te_viz_germline_08_re_rna.R")
} else {
  cat("8/13: Skipping RE-RNA differential expression analysis\n")
}

# 09 - PCA/UMAP (file: 08_pca_umap.R)
if (RUN_PCA_UMAP) {
  cat("9/13: PCA/UMAP analysis...\n")
  source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/02_te_viz_germline_08_pca_umap.R")
} else {
  cat("9/13: Skipping PCA/UMAP analysis\n")
}

# 10 - Methylation (file: 09_methylation.R)
if (RUN_METHYLATION) {
  cat("10/13: Methylation analysis...\n")
  source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/02_te_viz_germline_09_methylation.R")
} else {
  cat("10/13: Skipping methylation analysis\n")
}

# 11 - Taylor cohort (file: 10_taylor.R)
if (RUN_TAYLOR) {
  cat("11/13: Taylor cohort analysis...\n")
  source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/02_te_viz_germline_10_taylor.R")
} else {
  cat("11/13: Skipping Taylor cohort analysis\n")
}

# 12 - Subfamily analysis
if (RUN_SUBFAMILY) {
  cat("12/13: Subfamily analysis...\n")
  source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/02_te_viz_germline_12_subfamily.R")
} else {
  cat("12/13: Skipping subfamily analysis\n")
}

# 13 - Full-length, Young L1, Source (file: 13_fulllength_young_source.R)
if (RUN_FULLLENGTH_YOUNG_SOURCE) {
  cat("13/13: Full-length, Young L1, Source analysis...\n")
  source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/02_te_viz_germline_13_fulllength_young_source.R")
} else {
  cat("13/13: Skipping full-length, young L1, source analysis\n")
}

# Close PDF device
dev.off()

cat("\n======================================\n")
cat("  ALL ANALYSES COMPLETED\n")
cat("======================================\n")
cat("Output saved to:", plot_dir, "\n")
cat("PDF: graph_output.pdf\n")
cat("Text logs: See each module's output directory for *_GRAPH_OUTPUT.txt files\n")
