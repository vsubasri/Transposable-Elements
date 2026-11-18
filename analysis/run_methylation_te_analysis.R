#!/usr/bin/env Rscript

# ==============================================================================
# Methylation-TE Enrichment Analysis
# ==============================================================================
#
# Purpose: Identify differentially methylated probes (DMPs) between TP53 mutant
#          and WT samples, then test for enrichment of TEs at those sites.
#
# Usage: source("run_methylation_te_analysis.R")
#        results <- run_analysis()
#
# ==============================================================================

# Source functions
source("/Users/vallijahsubasri/Transposable-Elements/analysis/methylation_te.R")

# ==============================================================================
# CONFIGURATION
# ==============================================================================

# Input/Output paths
config <- list(
  methylation_file = "/Users/vallijahsubasri/Transposable-Elements/Noob_beta_mut_wt.rds",
  output_dir = "/Users/vallijahsubasri/Transposable-Elements/analysis/methylation_te_analysis/",

  # Genome
  genome = "hg19",  # or "hg38"

  # DMR thresholds
  delta_beta_threshold = 0.1,  # 10% difference cutoff for bumphunter
  fdr_threshold = 0.05,
  bumphunter_permutations = 100  # increase for final analysis (e.g., 1000)
)

# Create output directory
dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# MAIN ANALYSIS
# ==============================================================================

run_analysis <- function() {

  cat("=======================================================\n")
  cat("    Methylation-TE Enrichment Analysis\n")
  cat("=======================================================\n\n")

  # -------------------------
  # Step 1: Load methylation data
  # -------------------------
  cat("STEP 1: Loading methylation data\n")
  cat("-------------------------------\n")

  meth_data <- readRDS(config$methylation_file)
  cat("Dimensions:", nrow(meth_data), "probes x", ncol(meth_data), "samples\n")
  cat("\nSample names:\n")
  print(colnames(meth_data))

  # -------------------------
  # Step 2: Define groups
  # -------------------------
  cat("\n\nSTEP 2: Defining groups\n")
  cat("------------------------\n")

  sample_names <- colnames(meth_data)

  # Auto-detect mutant/WT from column names
  if (any(grepl("mut|wt", sample_names, ignore.case = TRUE))) {
    group_info <- ifelse(grepl("mut", sample_names, ignore.case = TRUE), "mutant", "WT")
    names(group_info) <- sample_names
  } else {
    cat("ERROR: Could not detect groups from sample names.\n")
    cat("Please manually define group_info.\n")
    return(NULL)
  }

  cat("Group assignments:\n")
  print(table(group_info))

  # -------------------------
  # Step 3: Load RepeatMasker (TE annotations)
  # -------------------------
  cat("\n\nSTEP 3: Loading TE annotations\n")
  cat("-------------------------------\n")

  rmsk_file <- file.path(config$output_dir, paste0("rmsk_", config$genome, ".txt.gz"))
  rmsk_file <- download_repeatmasker(genome = config$genome, output_file = rmsk_file)
  rmsk <- load_repeatmasker(rmsk_file)

  # -------------------------
  # Step 4: Get probe coordinates
  # -------------------------
  cat("\n\nSTEP 4: Getting probe coordinates\n")
  cat("----------------------------------\n")

  array_type <- ifelse(nrow(meth_data) > 500000, "EPIC", "450K")
  cat("Detected array type:", array_type, "\n")

  probe_coords <- get_probe_coordinates(array_type = array_type, genome = config$genome)

  if (is.null(probe_coords)) {
    return(NULL)
  }

  probe_coords <- probe_coords %>% filter(probe_id %in% rownames(meth_data))
  cat("Probes with coordinates:", nrow(probe_coords), "\n")

  # -------------------------
  # Step 5: Find probe-TE overlaps
  # -------------------------
  cat("\n\nSTEP 5: Finding probe-TE overlaps\n")
  cat("----------------------------------\n")

  overlap_df <- find_probe_te_overlaps(probe_coords, rmsk)

  write.csv(overlap_df,
            file.path(config$output_dir, "probe_te_overlaps.csv"),
            row.names = FALSE)

  # -------------------------
  # Step 6: Differential methylation (bumphunter DMRs)
  # -------------------------
  cat("\n\nSTEP 6: Differential methylation analysis (bumphunter)\n")
  cat("-------------------------------------------------------\n")

  dmr_results <- run_differential_methylation(
    meth_data,
    group_info,
    probe_coords,
    cutoff = config$delta_beta_threshold,
    B = config$bumphunter_permutations,
    fdr_thresh = config$fdr_threshold
  )

  write.csv(dmr_results$dmrs,
            file.path(config$output_dir, "dmr_results.csv"),
            row.names = FALSE)
  write.csv(dmr_results$probes,
            file.path(config$output_dir, "probe_dmr_mapping.csv"),
            row.names = FALSE)

  # -------------------------
  # Step 7: TE enrichment
  # -------------------------
  cat("\n\nSTEP 7: TE enrichment analysis\n")
  cat("-------------------------------\n")

  enrichment_results <- run_te_enrichment(dmr_results, overlap_df, rmsk)

  write.csv(enrichment_results$by_class,
            file.path(config$output_dir, "te_class_enrichment.csv"),
            row.names = FALSE)
  if (!is.null(enrichment_results$dmr_te_overlap)) {
    write.csv(enrichment_results$dmr_te_overlap,
              file.path(config$output_dir, "dmr_te_overlaps.csv"),
              row.names = FALSE)
  }

  # -------------------------
  # Step 8: Generate plots
  # -------------------------
  cat("\n\nSTEP 8: Generating plots\n")
  cat("-------------------------\n")

  plot_enrichment_heatmap(
    enrichment_results$by_class,
    file.path(config$output_dir, "te_enrichment_heatmap.pdf")
  )

  plot_enrichment_barplot(
    enrichment_results$by_class,
    file.path(config$output_dir, "te_enrichment_barplot.pdf")
  )

  plot_volcano(
    enrichment_results$probe_annotated,
    config$delta_beta_threshold,
    config$fdr_threshold,
    file.path(config$output_dir, "dmr_probe_summary.pdf")
  )

  plot_dmr_distribution(
    enrichment_results$dmr_table,
    file.path(config$output_dir, "dmr_size_distribution.pdf")
  )

  cat("Plots saved to:", config$output_dir, "\n")

  # -------------------------
  # Save results
  # -------------------------
  saveRDS(enrichment_results, file.path(config$output_dir, "enrichment_results.rds"))

  cat("\n=======================================================\n")
  cat("    Analysis Complete!\n")
  cat("=======================================================\n")
  cat("\nOutput files:\n")
  cat("  - probe_te_overlaps.csv\n")
  cat("  - dmr_results.csv (bumphunter DMRs)\n")
  cat("  - probe_dmr_mapping.csv\n")
  cat("  - dmr_te_overlaps.csv\n")
  cat("  - te_class_enrichment.csv\n")
  cat("  - te_enrichment_heatmap.pdf\n")
  cat("  - te_enrichment_barplot.pdf\n")
  cat("  - dmr_volcano.pdf\n")
  cat("  - enrichment_results.rds\n")
  cat("\nResults directory:", config$output_dir, "\n")

  return(enrichment_results)
}

# ==============================================================================
# RUN
# ==============================================================================

# Uncomment to run:
# results <- run_analysis()
