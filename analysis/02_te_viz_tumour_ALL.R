#!/usr/bin/env Rscript

# Tumour TE Visualization - Run All Analyses
# This script runs all tumour visualization analyses in sequence

cat("======================================\n")
cat("  TUMOUR TE VISUALIZATION - ALL\n")
cat("======================================\n\n")

#### TEST MODE ####
# Enable test mode BEFORE sourcing common setup
TEST_MODE <- FALSE

#### MODULE SELECTION ####
# Set to FALSE to skip a module
RUN_DESCRIPTIVE <- TRUE
RUN_CLINICAL <- TRUE
RUN_SPECIFIC_TES <- TRUE
RUN_LOH <- TRUE
RUN_SURVIVAL <- TRUE
RUN_PATHWAY <- TRUE
RUN_RNA <- TRUE
RUN_RE <- TRUE
RUN_RE_RNA <- TRUE
RUN_PCA_UMAP <- TRUE
RUN_METHYLATION <- TRUE
RUN_SV_OVERLAP <- TRUE
RUN_TAYLOR <- TRUE
RUN_CHROMOTHRIPSIS <- TRUE
RUN_FULLLENGTH_YOUNG_SOURCE <- TRUE

# Source common setup
source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/00_viz_common_setup.R")

# Load ALL tumour data (master script needs everything)
REQUIRED_DATA <- c("all")
source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/00_viz_load_data_tumour.R")

cat("\n=== Running Analysis Scripts ===\n\n")

# 01 - Descriptive statistics
if (RUN_DESCRIPTIVE) {
  cat("1/14: Descriptive statistics...\n")
  source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/02_te_viz_tumour_01_descriptive.R")
} else {
  cat("1/14: Skipping descriptive statistics\n")
}

# 02 - Clinical associations
if (RUN_CLINICAL) {
  cat("2/14: Clinical associations...\n")
  source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/02_te_viz_tumour_02_clinical.R")
} else {
  cat("2/14: Skipping clinical associations\n")
}

# 03 - Specific TEs
if (RUN_SPECIFIC_TES) {
  cat("3/14: Specific TEs analysis...\n")
  source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/02_te_viz_tumour_03_specific_tes.R")
} else {
  cat("3/14: Skipping specific TEs analysis\n")
}

# 04 - LOH timing
if (RUN_LOH) {
  cat("4/14: LOH timing analysis...\n")
  source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/02_te_viz_tumour_04_loh.R")
} else {
  cat("4/14: Skipping LOH timing analysis\n")
}

# 05 - Survival analysis
if (RUN_SURVIVAL) {
  cat("5/14: Survival analysis...\n")
  source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/02_te_viz_tumour_05_survival.R")
} else {
  cat("5/14: Skipping survival analysis\n")
}

# 06 - Pathway analysis
if (RUN_PATHWAY) {
  cat("6/14: Pathway analysis...\n")
  source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/02_te_viz_tumour_06_pathway.R")
} else {
  cat("6/14: Skipping pathway analysis\n")
}

# 07 - RNA expression
if (RUN_RNA) {
  cat("7/14: RNA expression analysis...\n")
  source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/02_te_viz_tumour_07_rna.R")
} else {
  cat("7/14: Skipping RNA expression analysis\n")
}

# 08 - Regulatory elements (ORA)
if (RUN_RE) {
  cat("8/14: Regulatory elements analysis...\n")
  source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/02_te_viz_tumour_08_re.R")
} else {
  cat("8/14: Skipping regulatory elements analysis\n")
}

# 09 - RE-RNA differential expression (GSEA) - NEW
if (RUN_RE_RNA) {
  cat("9/14: RE-RNA differential expression analysis...\n")
  source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/02_te_viz_tumour_09_re_rna.R")
} else {
  cat("9/14: Skipping RE-RNA differential expression analysis\n")
}

# 10 - PCA/UMAP (file: 09_pca_umap.R)
if (RUN_PCA_UMAP) {
  cat("10/14: PCA/UMAP analysis...\n")
  source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/02_te_viz_tumour_09_pca_umap.R")
} else {
  cat("10/14: Skipping PCA/UMAP analysis\n")
}

# 11 - Methylation (file: 10_methylation.R)
if (RUN_METHYLATION) {
  cat("11/14: Methylation analysis...\n")
  source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/02_te_viz_tumour_10_methylation.R")
} else {
  cat("11/14: Skipping methylation analysis\n")
}

# 12 - SV overlap (file: 11_sv_overlap.R)
if (RUN_SV_OVERLAP) {
  cat("12/14: SV overlap analysis...\n")
  source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/02_te_viz_tumour_11_sv_overlap.R")
} else {
  cat("12/14: Skipping SV overlap analysis\n")
}

# 13 - Taylor cohort (file: 12_taylor.R)
if (RUN_TAYLOR) {
  cat("13/14: Taylor cohort analysis...\n")
  source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/02_te_viz_tumour_12_taylor.R")
} else {
  cat("13/14: Skipping Taylor cohort analysis\n")
}

# 14 - Chromothripsis (file: 13_chromothripsis.R)
if (RUN_CHROMOTHRIPSIS) {
  cat("14/15: Chromothripsis analysis...\n")
  source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/02_te_viz_tumour_13_chromothripsis.R")
} else {
  cat("14/15: Skipping chromothripsis analysis\n")
}

# 15 - Full-length, Young L1, Source (file: 14_fulllength_young_source.R)
if (RUN_FULLLENGTH_YOUNG_SOURCE) {
  cat("15/15: Full-length, Young L1, Source analysis...\n")
  source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/02_te_viz_tumour_14_fulllength_young_source.R")
} else {
  cat("15/15: Skipping full-length, young L1, source analysis\n")
}

# Close PDF device
dev.off()

cat("\n======================================\n")
cat("  ALL ANALYSES COMPLETED\n")
cat("======================================\n")
cat("Output saved to:", plot_dir, "\n")
cat("PDF: graph_output.pdf\n")
cat("Text logs: See each module's output directory for *_GRAPH_OUTPUT.txt files\n")
