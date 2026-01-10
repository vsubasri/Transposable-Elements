#!/usr/bin/env Rscript

# Tumour TE Visualization - RNA Expression (Multi-Database GSEA)
# LM-based expression testing with t-statistics for GSEA ranking
# Databases: GO_BP, Reactome, MSigDB_Hallmark, MSigDB_Oncogenic (no KEGG)

# Source common setup and load data
source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/00_viz_common_setup.R")
REQUIRED_DATA <- c("count_matrix", "split", "rna", "clinical")
source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/00_viz_load_data_tumour.R")

# Initialize module-specific text output
init_module_sink(paste0(plot_dir, "rna/"), "RNA")

cat("Running 02_te_viz_tumour_07_rna.R...\n")
cat("LM-based expression testing with t-statistics for GSEA ranking\n")
cat("Databases: GO_BP, Reactome, MSigDB_Hallmark, MSigDB_Oncogenic\n\n")

# Initialize run summary
run_summary <- data.frame(
  analysis_name = character(),
  cohort = character(),
  n_te_samples = integer(),
  n_rna_samples = integer(),
  n_matched_samples = integer(),
  minsample = integer(),
  genes_tested = integer(),
  group_by = character(),
  genes_with_t_stats = integer(),
  gsea_ran = logical(),
  databases_with_results = character(),
  error = character(),
  timestamp = character(),
  stringsAsFactors = FALSE
)

#### CONFIGURATION ####

# Databases to run (no KEGG per plan)
DATABASES_TO_RUN <- c("GO_BP", "Reactome", "MSigDB_Hallmark", "MSigDB_Oncogenic")

# GSEA parameter settings (per plan)
MIN_SAMPLES_VALUES <- c(3, 5)
Q_PATHWAY_VALUES <- c(0.05, 0.1, 0.25)  # For pathway filtering

# Create base output directory
rna_dir <- paste0(plot_dir, "rna/")
dir.create(rna_dir, showWarnings = FALSE, recursive = TRUE)

cat("Base output directory:", rna_dir, "\n\n")

#### ADD 3-LEVEL TP53 CLASSIFICATION ####
cat("Adding 3-level TP53 classification (Germline/Somatic/WT)...\n")
te_aff_split_t <- add_tp53_3level(te_aff_split_t)
cat("TP53_3level distribution:\n")
print(table(te_aff_split_t$TP53_3level, useNA = "always"))
cat("\n")

if (exists("te_kics_split_t")) {
  te_kics_split_t <- add_tp53_3level(te_kics_split_t)
}
if (exists("te_lfs_split_t")) {
  te_lfs_split_t <- add_tp53_3level(te_lfs_split_t)
}
if (exists("te_all_split_t")) {
  te_all_split_t <- add_tp53_3level(te_all_split_t)
}

#### PROCESS RNA DATA ####
write_output(quote(NULL), "Processing Tumour RNA Data")

# Validate, rename, combine and filter RNA datasets
rna_processing_result <- process_tumor_rna_data(kics_rna, lfs_rna, stjude_rna, lfs_wgs2rna)
rna_ready_for_filtering <- rna_processing_result$rna_data

cat("Tumor RNA data prepared:", nrow(rna_ready_for_filtering), "genes\n")

#### MATCH TE AND RNA SAMPLES ####
write_output(quote(NULL), "Matching TE and RNA Samples")

# Save matched samples info
rna_files_dir <- paste0(rna_dir, "files/")
dir.create(rna_files_dir, showWarnings = FALSE, recursive = TRUE)

# Match TE sample names with RNA sample names
rna_matching_result <- match_te_rna_samples(
  rna_ready_for_filtering,
  te_all_t,
  r_dir_files = rna_files_dir,
  dna_rna_mapping = matched_dna_rna,
  original_kics_rna = kics_rna,
  te_all_all_t = te_all_all_t
)
rna_filtered <- rna_matching_result$rna_filtered
sample_overlap <- rna_matching_result$sample_overlap

cat("Final RNA dataset:", nrow(rna_filtered), "genes x", ncol(rna_filtered) - 1, "samples\n")
cat("Matched samples:", length(sample_overlap), "\n\n")

# Create matched_samples dataframe for LM function
# For tumor data, TE samples and RNA samples should have the same names
# sample_overlap contains the matched sample names
rna_sample_names <- colnames(rna_filtered)[-1]  # Exclude gene_name column

# Get all unique TE sample names from split data
all_te_sample_names <- unique(c(
  if(exists("te_all_split_t")) te_all_split_t$sample else character(0),
  if(exists("te_aff_split_t")) te_aff_split_t$sample else character(0),
  if(exists("te_kics_split_t")) te_kics_split_t$sample else character(0),
  if(exists("te_lfs_split_t")) te_lfs_split_t$sample else character(0)
))

cat("Creating TE to RNA sample mapping...\n")
cat("Total unique tumor TE samples:", length(all_te_sample_names), "\n")

# For tumor data, match by exact name or base ID
te_base_ids <- sub("_.*", "", all_te_sample_names)
rna_base_ids <- sub("_.*", "", rna_sample_names)

matched_samples_for_lm <- data.frame(
  dna_sample = all_te_sample_names,
  stringsAsFactors = FALSE
)

# Try exact match first, then base ID match
matched_samples_for_lm$rna_sample <- sapply(seq_along(all_te_sample_names), function(i) {
  te_sample <- all_te_sample_names[i]
  te_base <- te_base_ids[i]

  # Try exact match first
  if (te_sample %in% rna_sample_names) {
    return(te_sample)
  }
  # Try base ID match
  idx <- match(te_base, rna_base_ids)
  if (!is.na(idx)) {
    return(rna_sample_names[idx])
  }
  return(NA_character_)
})

matched_samples_for_lm <- matched_samples_for_lm[!is.na(matched_samples_for_lm$rna_sample), ]
cat("Matched", nrow(matched_samples_for_lm), "tumor TE samples to RNA samples\n")
cat("Example matches (first 5):\n")
print(head(matched_samples_for_lm, 5))

#### RNA EXPRESSION ANALYSIS ####
write_output(quote(NULL), "RNA Expression Analysis")

# Define gene of interest
gene_of_interest <- "MROH7"

# Check if we have enough samples and gene exists
if (ncol(rna_filtered) > 1 && gene_of_interest %in% rna_filtered$gene_name) {
  # Plot 1: Gene Expression by TE Status
  tryCatch({
    result <- plot_gene_expression_te(rna_filtered, te_all_split_t, gene_of_interest, "TP53_status", "TE Status", paste0(gene_of_interest, " Expression (FPKM)"), FALSE, TRUE, FALSE)
    titled_print(result$plot, paste0(gene_of_interest, " Expression by TE Status"))
    ggsave(paste0(rna_dir, tolower(gene_of_interest), "_expression_te_status.png"), width = 9, height = 5)
  }, error = function(e) cat("Warning: Could not create TE status plot:", e$message, "\n"))

  # Plot 2: Gene Expression by TP53 Status
  tryCatch({
    result <- plot_gene_expression_te(rna_filtered, te_all_split_t, gene_of_interest, "TP53_status", "TP53 Status", paste0(gene_of_interest, " Expression (FPKM)"), FALSE, FALSE, FALSE)
    titled_print(result$plot, paste0(gene_of_interest, " Expression by TP53 Status"))
    ggsave(paste0(rna_dir, tolower(gene_of_interest), "_expression_tp53_status.png"), width = 9, height = 5)
  }, error = function(e) cat("Warning: Could not create TP53 status plot:", e$message, "\n"))

  # Plot 3: Gene Expression by Tumor Type
  tryCatch({
    result <- plot_gene_expression_te(rna_filtered, te_all_split_t, gene_of_interest, "tumor_type", "Tumor Type", paste0(gene_of_interest, " Expression (FPKM)"), FALSE, FALSE, FALSE)
    titled_print(result$plot, paste0(gene_of_interest, " Expression by Tumor Type"))
    ggsave(paste0(rna_dir, tolower(gene_of_interest), "_expression_tumor_type.png"), width = 9, height = 5)
  }, error = function(e) cat("Warning: Could not create tumor type plot:", e$message, "\n"))
} else {
  cat("Skipping expression plots - insufficient data or gene not found\n")
}

#### SAME CANCER TYPE TE EXPRESSION ANALYSIS ####
cat("\n===== SAME CANCER TYPE TE EXPRESSION ANALYSIS =====\n")

# Genes to analyze with same-cancer-type filtering
genes_for_same_cancer_analysis <- list(
  list(gene1 = "LRP1B", gene2 = NULL),
  list(gene1 = "ERBB4", gene2 = "ALK"),
  list(gene1 = "ALK", gene2 = NULL)
)

for (gene_config in genes_for_same_cancer_analysis) {
  gene1 <- gene_config$gene1
  gene2 <- gene_config$gene2

  cat("\nAnalyzing", gene1, "with same-cancer-type filtering...\n")
  if (!is.null(gene2)) {
    cat("Also plotting", gene2, "expression by TE in", gene1, "\n")
  }

  tryCatch({
    result <- plot_gene_expression_same_cancer(rna_filtered, te_all_split_t, gene1, gene2 = gene2)

    if (!is.null(result$plot)) {
      titled_print(result$plot, paste0(gene1, " Expression by TE Status (same cancer type)"))
      ggsave(paste0(rna_dir, tolower(gene1), "_expression_te_status_same_cancer.png"),
             width = 9, height = 5)
    }

    if (!is.null(gene2) && !is.null(result$plot2)) {
      titled_print(result$plot2, paste0(gene2, " Expression by TE in ", gene1, " (same cancer type)"))
      ggsave(paste0(rna_dir, tolower(gene2), "_by_te_in_", tolower(gene1), "_same_cancer.png"),
             width = 9, height = 5)
    }
  }, error = function(e) cat("  Warning: Could not create plot:", e$message, "\n"))
}

#### HELPER FUNCTION: Run LM-based GSEA analysis ####
run_lm_gsea_analysis <- function(rna_data,
                                  te_split_data,
                                  analysis_name,
                                  base_dir,
                                  cohort_name = "Unknown",
                                  actual_matched_samples = 0,
                                  databases = DATABASES_TO_RUN,
                                  gene_col = "Gene_name",
                                  group_by_gene = TRUE,
                                  matched_samples = NULL) {

  cat("\n========================================\n")
  cat("  LM-GSEA Analysis:", analysis_name, "\n")
  cat("========================================\n")

  # Get sample counts
  n_te_samples <- length(unique(te_split_data$sample))
  n_rna_samples <- ncol(rna_data) - 1  # Exclude gene_name column

  # Create base files directory
  base_files_dir <- paste0(base_dir, "files/")
  dir.create(base_files_dir, showWarnings = FALSE, recursive = TRUE)

  for (minsample in MIN_SAMPLES_VALUES) {
    cat("\n--- minsample =", minsample, "---\n")

    # Track for run summary
    genes_tested <- 0
    genes_with_t_stats <- 0
    gsea_ran <- FALSE
    dbs_with_results <- c()

    # Run LM-based expression testing
    tryCatch({
      lm_results <- test_te_expression_effects_lm(
        rna_df = rna_data,
        te_split_df = te_split_data,
        min_samples = minsample,
        group_by_gene = group_by_gene,
        matched_samples = matched_samples
      )

      # Always save LM results, even if NULL/empty
      lm_file <- paste0(base_files_dir, "lm_results_", analysis_name, "_minsample", minsample, ".csv")
      if (!is.null(lm_results) && nrow(lm_results) > 0) {
        write.csv(lm_results, lm_file, row.names = FALSE)
        cat("  Saved LM results:", basename(lm_file), "\n")
        genes_tested <- nrow(lm_results)
      } else {
        # Save empty results with header
        empty_df <- data.frame(gene = character(), p_value = numeric(), t_statistic = numeric())
        write.csv(empty_df, lm_file, row.names = FALSE)
        cat("  Saved empty LM results:", basename(lm_file), "\n")
      }

      if (is.null(lm_results) || nrow(lm_results) == 0) {
        cat("No results from LM analysis\n")

        # Add to run summary
        run_summary <<- rbind(run_summary, data.frame(
          analysis_name = analysis_name,
          cohort = cohort_name,
          n_te_samples = n_te_samples,
          n_rna_samples = n_rna_samples,
          n_matched_samples = actual_matched_samples,
          minsample = minsample,
          group_by = if(group_by_gene) "gene" else "coordinate",
          genes_tested = genes_tested,
          genes_with_t_stats = genes_with_t_stats,
          gsea_ran = gsea_ran,
          databases_with_results = "",
          error = "",
          timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
          stringsAsFactors = FALSE
        ))
        next
      }

      cat("LM analysis complete:", nrow(lm_results), "genes tested\n")

      # Get t-statistics for GSEA ranking (extract from LM results)
      t_stats <- setNames(lm_results$t_statistic, lm_results$gene)

      # Remove duplicates (keep max absolute t-stat per gene)
      if (any(duplicated(names(t_stats)))) {
        t_stats <- tapply(t_stats, names(t_stats), function(x) x[which.max(abs(x))])
      }

      # Sort descending for GSEA
      t_stats <- sort(t_stats, decreasing = TRUE)
      genes_with_t_stats <- length(t_stats)

      if (length(t_stats) < 10) {
        cat("Too few genes with valid t-statistics:", length(t_stats), "\n")

        # Add to run summary
        run_summary <<- rbind(run_summary, data.frame(
          analysis_name = analysis_name,
          cohort = cohort_name,
          n_te_samples = n_te_samples,
          n_rna_samples = n_rna_samples,
          n_matched_samples = actual_matched_samples,
          minsample = minsample,
          group_by = if(group_by_gene) "gene" else "coordinate",
          genes_tested = genes_tested,
          genes_with_t_stats = genes_with_t_stats,
          gsea_ran = gsea_ran,
          databases_with_results = "",
          error = "",
          timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
          stringsAsFactors = FALSE
        ))
        next
      }

      cat("Genes with t-statistics:", length(t_stats), "\n")
      cat("t-stat range:", round(min(t_stats), 2), "to", round(max(t_stats), 2), "\n")

      gsea_ran <- TRUE

      # Run GSEA for each database
      for (db in databases) {
        cat("\n  Database:", db, "\n")

        # Create output directories
        analysis_dir <- paste0(base_dir, db, "/", analysis_name, "/")
        files_dir <- paste0(analysis_dir, "files/")
        dir.create(files_dir, showWarnings = FALSE, recursive = TRUE)

        # Save ranked gene list
        ranked_genes_df <- data.frame(
          gene = names(t_stats),
          t_statistic = as.numeric(t_stats),
          stringsAsFactors = FALSE
        )
        write.csv(ranked_genes_df,
                  paste0(files_dir, "rna_gsea_ranked_genes_", db, "_", analysis_name,
                         "_minsample", minsample, ".csv"),
                  row.names = FALSE)

        # Run GSEA
        tryCatch({
          gsea_result <- perform_gsea_single_database(
            gene_list = t_stats,
            database = db,
            pvalueCutoff = 1
          )

          if (is.null(gsea_result) || nrow(as.data.frame(gsea_result)) == 0) {
            cat("    No enriched terms found\n")
            next
          }

          dbs_with_results <- c(dbs_with_results, db)

          # Save full GSEA result object
          gsea_rds_file <- paste0(files_dir, "rna_gsea_result_", db, "_", analysis_name,
                                   "_minsample", minsample, ".rds")
          saveRDS(gsea_result, gsea_rds_file)
          cat("    Saved GSEA result object:", basename(gsea_rds_file), "\n")

          # Save full GSEA CSV (all results, not filtered)
          full_gsea_df <- as.data.frame(gsea_result)
          full_csv_file <- paste0(files_dir, "rna_gsea_full_", db, "_", analysis_name,
                                   "_minsample", minsample, ".csv")
          write.csv(full_gsea_df, full_csv_file, row.names = FALSE)
          cat("    Saved full GSEA CSV:", basename(full_csv_file), "(", nrow(full_gsea_df), "terms)\n")

          # Save unfiltered plots (using all results)
          unfiltered_prefix <- paste0("rna_gsea_", db, "_", analysis_name,
                                      "_minsample", minsample, "_all")
          tryCatch({
            # Dot plot with all results
            p_dot <- dotplot(gsea_result, showCategory = min(20, nrow(full_gsea_df)))
            ggsave(paste0(analysis_dir, unfiltered_prefix, "_dot.png"), p_dot, width = 12, height = 8)

            # Ridge plot
            tryCatch({
              p_ridge <- ridgeplot(gsea_result, showCategory = min(15, nrow(full_gsea_df)))
              ggsave(paste0(analysis_dir, unfiltered_prefix, "_ridge.png"), p_ridge, width = 10, height = 8)
            }, error = function(e) NULL)

            # GSEA plot for top 3 pathways
            if (nrow(full_gsea_df) >= 3) {
              tryCatch({
                p_gsea <- gseaplot2(gsea_result, geneSetID = 1:3)
                ggsave(paste0(analysis_dir, unfiltered_prefix, "_gsea.png"), p_gsea, width = 10, height = 8)
              }, error = function(e) NULL)
            }
            cat("    Saved unfiltered plots:", unfiltered_prefix, "\n")
          }, error = function(e) {
            cat("    Warning: Could not create unfiltered plots:", e$message, "\n")
          })

          # Post-hoc filter by qpathway
          for (qpathway in Q_PATHWAY_VALUES) {
            gsea_df <- as.data.frame(gsea_result)
            filtered_df <- gsea_df %>%
              filter(qvalue < qpathway)

            if (nrow(filtered_df) == 0) next

            # Create output prefix
            prefix <- paste0("rna_gsea_", db, "_", analysis_name,
                             "_minsample", minsample, "_qpathway", qpathway)

            # Save filtered CSV
            csv_file <- paste0(files_dir, prefix, ".csv")
            write.csv(filtered_df, csv_file, row.names = FALSE)

            # Create plots
            tryCatch({
              # Dot plot
              p_dot <- dotplot(gsea_result, showCategory = min(20, nrow(filtered_df)))
              ggsave(paste0(analysis_dir, prefix, "_dot.png"), p_dot, width = 12, height = 8)

              # Ridge plot for GSEA
              tryCatch({
                p_ridge <- ridgeplot(gsea_result, showCategory = min(15, nrow(filtered_df)))
                ggsave(paste0(analysis_dir, prefix, "_ridge.png"), p_ridge, width = 10, height = 8)
              }, error = function(e) NULL)

              # GSEA plot for top pathways
              if (nrow(filtered_df) > 0) {
                tryCatch({
                  p_gsea <- gseaplot2(gsea_result, geneSetID = 1:min(3, nrow(filtered_df)))
                  ggsave(paste0(analysis_dir, prefix, "_gsea.png"), p_gsea, width = 10, height = 8)
                }, error = function(e) NULL)
              }
            }, error = function(e) {
              cat("    Warning: Could not create plots:", e$message, "\n")
            })

            cat("    Saved:", prefix, "(", nrow(filtered_df), "terms)\n")
          }
        }, error = function(e) {
          cat("    Error running GSEA for", db, ":", e$message, "\n")
        })
      }

      # Add to run summary
      run_summary <<- rbind(run_summary, data.frame(
        analysis_name = analysis_name,
        cohort = cohort_name,
        n_te_samples = n_te_samples,
        n_rna_samples = n_rna_samples,
        n_matched_samples = actual_matched_samples,
        minsample = minsample,
        group_by = if(group_by_gene) "gene" else "coordinate",
        genes_tested = genes_tested,
        genes_with_t_stats = genes_with_t_stats,
        gsea_ran = gsea_ran,
        databases_with_results = paste(dbs_with_results, collapse = ","),
        error = "",
        timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
        stringsAsFactors = FALSE
      ))

    }, error = function(e) {
      cat("Error in LM/GSEA analysis:", e$message, "\n")

      # Still record the summary even on error
      run_summary <<- rbind(run_summary, data.frame(
        analysis_name = analysis_name,
        cohort = cohort_name,
        n_te_samples = n_te_samples,
        n_rna_samples = n_rna_samples,
        n_matched_samples = actual_matched_samples,
        minsample = minsample,
        group_by = if(group_by_gene) "gene" else "coordinate",
        genes_tested = genes_tested,
        genes_with_t_stats = genes_with_t_stats,
        gsea_ran = gsea_ran,
        databases_with_results = paste(dbs_with_results, collapse = ","),
        error = e$message,
        timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
        stringsAsFactors = FALSE
      ))
    })
  }
}

#### GSEA ANALYSIS 1: BY GENE (All cohort) ####
write_output(quote(NULL), "GSEA Analysis - By Gene (All Tumour Cohort)")

run_lm_gsea_analysis(
  rna_data = rna_filtered,
  te_split_data = te_all_split_t,
  analysis_name = "all_by_gene",
  base_dir = rna_dir,
  cohort_name = "All",
  actual_matched_samples = length(sample_overlap),
  group_by_gene = TRUE,
  matched_samples = matched_samples_for_lm
)

#### GSEA ANALYSIS 2: BY TE COORDINATE (All cohort) ####
write_output(quote(NULL), "GSEA Analysis - By TE Coordinate (All Tumour Cohort)")

run_lm_gsea_analysis(
  rna_data = rna_filtered,
  te_split_data = te_all_split_t,
  analysis_name = "all_by_coord",
  base_dir = rna_dir,
  cohort_name = "All",
  actual_matched_samples = length(sample_overlap),
  group_by_gene = FALSE,
  matched_samples = matched_samples_for_lm
)

#### GSEA ANALYSIS 3: AFFECTED COHORT ####
write_output(quote(NULL), "GSEA Analysis - Affected Cohort")

if (exists("te_aff_split_t") && nrow(te_aff_split_t) > 0) {
  run_lm_gsea_analysis(
    rna_data = rna_filtered,
    te_split_data = te_aff_split_t,
    analysis_name = "affected_by_gene",
    base_dir = rna_dir,
    cohort_name = "Affected",
    actual_matched_samples = length(sample_overlap),
    group_by_gene = TRUE,
    matched_samples = matched_samples_for_lm
  )
}

#### GSEA ANALYSIS 4: KICS COHORT ####
write_output(quote(NULL), "GSEA Analysis - KICS Cohort")

if (exists("te_kics_split_t") && nrow(te_kics_split_t) > 0) {
  run_lm_gsea_analysis(
    rna_data = rna_filtered,
    te_split_data = te_kics_split_t,
    analysis_name = "kics_by_gene",
    base_dir = rna_dir,
    cohort_name = "KICS",
    actual_matched_samples = length(sample_overlap),
    group_by_gene = TRUE,
    matched_samples = matched_samples_for_lm
  )
}

#### GSEA ANALYSIS 5: LFS COHORT ####
write_output(quote(NULL), "GSEA Analysis - LFS Cohort")

if (exists("te_lfs_split_t") && nrow(te_lfs_split_t) > 0) {
  run_lm_gsea_analysis(
    rna_data = rna_filtered,
    te_split_data = te_lfs_split_t,
    analysis_name = "lfs_by_gene",
    base_dir = rna_dir,
    cohort_name = "LFS",
    actual_matched_samples = length(sample_overlap),
    group_by_gene = TRUE,
    matched_samples = matched_samples_for_lm
  )
}

#### SAVE DATA FOR RE-RNA SCRIPT ####
cat("\n===== SAVING DATA =====\n")

# Save rna_filtered for use by RE-RNA script
cat("Saving rna_filtered for RE-RNA script...\n")
saveRDS(rna_filtered, paste0(rna_files_dir, "rna_filtered_tumour.rds"))
cat("Saved rna_filtered to:", paste0(rna_files_dir, "rna_filtered_tumour.rds\n"))

# Save matched_samples for use by RE-RNA script
cat("Saving matched_samples for RE-RNA script...\n")
saveRDS(matched_samples_for_lm, paste0(rna_files_dir, "matched_samples_tumour.rds"))
cat("Saved matched_samples to:", paste0(rna_files_dir, "matched_samples_tumour.rds\n"))

# Save run summary
run_summary_file <- paste0(rna_dir, "RUN_SUMMARY.csv")
write.csv(run_summary, run_summary_file, row.names = FALSE)
cat("Saved run summary to:", run_summary_file, "\n")

#### SUMMARY ####
write_output(quote(NULL), "RNA Analysis Summary")

cat("Analysis completed:\n")
cat("  1. RNA data processing and sample matching\n")
cat("  2. Gene expression plots for", gene_of_interest, "\n")
cat("  3. Same-cancer-type TE expression analysis\n")
cat("  4. LM-based TE expression analysis (replaces Wilcoxon)\n")
cat("  5. Multi-database GSEA pathway analysis (t-statistics ranking)\n\n")

cat("Analyses run:\n")
cat("  - All tumour cohort by gene\n")
cat("  - All tumour cohort by TE coordinate\n")
cat("  - Affected cohort by gene\n")
cat("  - KICS cohort by gene\n")
cat("  - LFS cohort by gene\n\n")

cat("Databases used:", paste(DATABASES_TO_RUN, collapse = ", "), "\n")
cat("Min samples:", paste(MIN_SAMPLES_VALUES, collapse = ", "), "\n")
cat("Q-pathway thresholds:", paste(Q_PATHWAY_VALUES, collapse = ", "), "\n")

cat("\nOutput structure: {database}/{analysis}/files/ and plots\n")
cat("Output location:", rna_dir, "\n")

cat("\n Script completed successfully\n")

# Close module-specific sink
close_module_sink()
