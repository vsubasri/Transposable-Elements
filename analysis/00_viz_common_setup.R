#!/usr/bin/env Rscript

# Common Setup for TE Visualization Scripts
# This file contains shared libraries, paths, and utility functions
# Source this at the beginning of each visualization script

# Check if already loaded (skip when running via ALL script)
if (exists("write_output") && exists("titled_print")) {
  cat("Common setup already loaded, skipping reload...\n")
  return(invisible(NULL))
}

#### CONFIGURATION ####
# Set to TRUE to use small test datasets for faster pipeline testing
# Set to FALSE to use full production datasets
TEST_MODE <- FALSE

#### LIBRARIES ####
suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
  library(ggsignif)
  library(tidyr)
  library(umap)
  library(pheatmap)
  library(purrr)
  library(RColorBrewer)
  library(caret)
  library("org.Hs.eg.db")
  library(clusterProfiler)
  library(AnnotationDbi)
  library(stringr)
  library(scales)
  library(glmnet)
  library(reshape2)
  library(ggrepel)
  library(GO.db)
  library(tibble)
  library(VennDiagram)
  library(eulerr)
  library(enrichplot)
  library(ComplexHeatmap)
  library(circlize)
  library(dplyr)
})

# Source functions
source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/functions_te.R")

#### UTILITY FUNCTIONS ####

# Helper function to print section headers and suppress unwanted output
write_output <- function(expr, section_title) {
  cat("\n")
  cat(paste0("========================================\n"))
  cat(paste0("  ", section_title, "\n"))
  cat(paste0("========================================\n"))
  result <- eval(expr, envir = parent.frame())
  # Only print dataframes if they're small summaries (not full data)
  if (is.data.frame(result) && nrow(result) > 0 && nrow(result) <= 50) {
    print(result)
  }
  cat("\n")
  invisible(result)
}

# Helper function to add a title to a plot before printing to the PDF
titled_print <- function(plot, title) {
  if (inherits(plot, "ggplot")) {
    print(plot + ggtitle(title))
  } else {
    print(plot)
  }
}

cat("✓ Common setup loaded successfully\n")
