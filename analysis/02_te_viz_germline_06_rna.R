#!/usr/bin/env Rscript

# Germline TE Visualization - RNA Expression (Multi-Database GSEA)
# LM-based expression testing with t-statistics for GSEA ranking
# Databases: GO_BP, Reactome, MSigDB_Hallmark, MSigDB_Oncogenic (no KEGG)

# Source common setup and load data
source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/00_viz_common_setup.R")
REQUIRED_DATA <- c("count_matrix", "split", "rna", "clinical")
source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/00_viz_load_data_germline.R")

# Initialize module-specific text output
init_module_sink(paste0(plot_dir, "rna/"), "RNA")

cat("Running 02_te_viz_germline_06_rna.R...\n")
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
te_aff_split <- add_tp53_3level(te_aff_split)
cat("TP53_3level distribution:\n")
print(table(te_aff_split$TP53_3level, useNA = "always"))
cat("\n")

if (exists("te_kics_split")) {
  te_kics_split <- add_tp53_3level(te_kics_split)
}
if (exists("te_lfs_split")) {
  te_lfs_split <- add_tp53_3level(te_lfs_split)
}

#### PROCESS RNA DATA ####
write_output(quote(NULL), "Processing Germline RNA Data")
rna_processing_result <- process_germline_rna_data(kics_rna, lfs_rna, stjude_rna, matched_dna_rna, lfs_wgs2rna)
rna_ready_for_filtering <- rna_processing_result$rna_data

write_output(quote({
  cat("Total RNA samples available:", ncol(rna_ready_for_filtering) - 1, "\n")
  cat("Total genes with expression:", nrow(rna_ready_for_filtering), "\n")
}), "RNA Data Summary")

#### MATCH RNA AND TE SAMPLES ####
write_output(quote(NULL), "Matching Germline TE with Tumor RNA Samples")

# Load tumor-matched RNA sample list (from tumour RNA script output)
tumor_rna_file_icloud <- "/Users/briannelaverty/Library/Mobile Documents/com~apple~CloudDocs/Malkin_lab/figures/te/graphs/tumour/rare/rna/files/successfully_matched_samples.csv"
tumor_rna_file_local <- "/tmp/successfully_matched_samples.csv"

# Copy from iCloud to local temp (workaround for R permission issue with iCloud)
system(paste("cp", shQuote(tumor_rna_file_icloud), shQuote(tumor_rna_file_local)))

if (!file.exists(tumor_rna_file_local)) {
  stop("Failed to copy file from iCloud. Check if source exists: ", tumor_rna_file_icloud, "\nRun 02_te_viz_tumour_07_rna.R first.")
}

tumor_rna_samples <- read.csv(tumor_rna_file_local, stringsAsFactors = FALSE)$sample_name
cat("Tumor RNA samples (from tumor analysis):", length(tumor_rna_samples), "\n")

# Direct match - RNA columns already renamed to match tumour sample names
matched_rna_samples <- intersect(tumor_rna_samples, colnames(rna_ready_for_filtering))
cat("Matched RNA samples (available in RNA data):", length(matched_rna_samples), "\n")

if (length(matched_rna_samples) == 0) {
  stop("No matching samples found between tumor RNA file and available RNA data")
}

# Filter RNA data
rna_filtered <- rna_ready_for_filtering %>%
  select(gene_name, all_of(matched_rna_samples))

cat("Final RNA dataset:", nrow(rna_filtered), "genes x", ncol(rna_filtered) - 1, "samples\n")

# Create matched_samples mapping for germline TE -> tumor RNA
# This is critical: germline samples are e.g. "0045_288555_N" while tumor RNA is "0045_288555_T"
# Match using patient base ID (first part before underscore)
cat("\nCreating germline TE to tumor RNA sample mapping...\n")

# Get RNA sample base IDs
rna_sample_base_ids <- sub("_.*", "", matched_rna_samples)

# Get unique TE sample names from ALL cohorts' split data
all_te_sample_names <- unique(c(
  te_aff_split$sample,
  if(exists("te_kics_split")) te_kics_split$sample else character(0),
  if(exists("te_lfs_split")) te_lfs_split$sample else character(0)
))

cat("Total unique germline TE samples across cohorts:", length(all_te_sample_names), "\n")

# Create matched_samples dataframe using vectorized approach
te_base_ids <- sub("_.*", "", all_te_sample_names)

# For each TE sample, find matching RNA sample
matched_samples_for_lm <- data.frame(
  dna_sample = all_te_sample_names,
  stringsAsFactors = FALSE
)

# Add RNA matches
matched_samples_for_lm$rna_sample <- sapply(te_base_ids, function(te_base) {
  idx <- match(te_base, rna_sample_base_ids)
  if (!is.na(idx)) matched_rna_samples[idx] else NA_character_
})

# Keep only matched samples
matched_samples_for_lm <- matched_samples_for_lm[!is.na(matched_samples_for_lm$rna_sample), ]

cat("Matched", nrow(matched_samples_for_lm), "germline TE samples to tumor RNA samples\n")
cat("Example matches (first 5):\n")
print(head(matched_samples_for_lm, 5))

# Save matched samples info
rna_files_dir <- paste0(rna_dir, "files/")
dir.create(rna_files_dir, showWarnings = FALSE, recursive = TRUE)

write.csv(data.frame(rna_sample_id = matched_rna_samples),
          paste0(rna_files_dir, "germline_matched_rna_samples.csv"),
          row.names = FALSE)

#### RNA EXPRESSION ANALYSIS ####
write_output(quote(NULL), "RNA Expression Analysis")

# Define gene of interest for germline analysis
gene_of_interest <- "TP53"

# Check if we have enough samples and gene exists in the RNA data before plotting
if (ncol(rna_filtered) > 1 && gene_of_interest %in% rna_filtered$gene_name) {
  # Plot 1: Gene Expression by TE Status
  tryCatch({
    result <- plot_gene_expression_te(rna_filtered, te_aff_split, gene_of_interest, "TP53_status", "TE Status", paste0(gene_of_interest, " Expression (FPKM)"), FALSE, TRUE, FALSE)
    titled_print(result$plot, paste0(gene_of_interest, " Expression by TE Status"))
    ggsave(paste0(rna_dir, tolower(gene_of_interest), "_expression_te_status.png"), width = 9, height = 5)
  }, error = function(e) cat("Warning: Could not create TE status plot:", e$message, "\n"))

  # Plot 2: Gene Expression by TP53 Status with stats
  tryCatch({
    result <- plot_gene_expression_te(rna_filtered, te_aff_split, gene_of_interest, "TP53_status", "TP53 Status", paste0(gene_of_interest, " Expression (FPKM)"), FALSE, FALSE, FALSE)
    titled_print(result$plot, paste0(gene_of_interest, " Expression by TP53 Status (with stats)"))
    ggsave(paste0(rna_dir, tolower(gene_of_interest), "_expression_tp53_status_stats.png"), width = 9, height = 5)
  }, error = function(e) cat("Warning: Could not create TP53 status plot:", e$message, "\n"))

  # Plot 3: Gene Expression by Tumor Type
  tryCatch({
    result <- plot_gene_expression_te(rna_filtered, te_aff_split, gene_of_interest, "tumor_type", "Tumor Type", paste0(gene_of_interest, " Expression (FPKM)"), FALSE, FALSE, FALSE)
    titled_print(result$plot, paste0(gene_of_interest, " Expression by Tumor Type"))
    ggsave(paste0(rna_dir, tolower(gene_of_interest), "_expression_tumor_type.png"), width = 9, height = 5)
  }, error = function(e) cat("Warning: Could not create tumor type plot:", e$message, "\n"))
} else {
  if (ncol(rna_filtered) <= 1) {
    cat("No matched samples with RNA data; skipping", gene_of_interest, "expression plots\n")
  } else {
    cat(gene_of_interest, "not found in RNA data; skipping expression plots\n")
  }
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

#### GSEA ANALYSIS 1: BY GENE (Affected cohort) ####
write_output(quote(NULL), "GSEA Analysis - By Gene (Affected Cohort)")

run_lm_gsea_analysis(
  rna_data = rna_filtered,
  te_split_data = te_aff_split,
  analysis_name = "affected_by_gene",
  base_dir = rna_dir,
  cohort_name = "Affected",
  actual_matched_samples = length(matched_rna_samples),
  group_by_gene = TRUE,
  matched_samples = matched_samples_for_lm
)

#### GSEA ANALYSIS 2: BY TE COORDINATE (Affected cohort) ####
write_output(quote(NULL), "GSEA Analysis - By TE Coordinate (Affected Cohort)")

run_lm_gsea_analysis(
  rna_data = rna_filtered,
  te_split_data = te_aff_split,
  analysis_name = "affected_by_coord",
  base_dir = rna_dir,
  cohort_name = "Affected",
  actual_matched_samples = length(matched_rna_samples),
  group_by_gene = FALSE,
  matched_samples = matched_samples_for_lm
)

#### GSEA ANALYSIS 3: KICS COHORT ####
write_output(quote(NULL), "GSEA Analysis - KICS Cohort")

if (exists("te_kics_split") && nrow(te_kics_split) > 0) {
  run_lm_gsea_analysis(
    rna_data = rna_filtered,
    te_split_data = te_kics_split,
    analysis_name = "kics_by_gene",
    base_dir = rna_dir,
    cohort_name = "KICS",
    actual_matched_samples = length(matched_rna_samples),
    group_by_gene = TRUE,
    matched_samples = matched_samples_for_lm
  )
}

#### GSEA ANALYSIS 4: LFS COHORT ####
write_output(quote(NULL), "GSEA Analysis - LFS Cohort")

if (exists("te_lfs_split") && nrow(te_lfs_split) > 0) {
  run_lm_gsea_analysis(
    rna_data = rna_filtered,
    te_split_data = te_lfs_split,
    analysis_name = "lfs_by_gene",
    base_dir = rna_dir,
    cohort_name = "LFS",
    actual_matched_samples = length(matched_rna_samples),
    group_by_gene = TRUE,
    matched_samples = matched_samples_for_lm
  )
}

#### GSEA ANALYSIS 5: KICS VS HOSTSEQ ####
write_output(quote(NULL), "GSEA Analysis - KICS vs HostSeq")

if (exists("te_kics_hostseq") && nrow(te_kics_hostseq) > 0) {
  # Need split version - create if not exists
  if (!exists("te_kics_hostseq_split")) {
    te_kics_hostseq_split <- te_kics_hostseq
  }
  te_kics_hostseq_split <- add_tp53_3level(te_kics_hostseq_split)

  cohort_groups <- unique(te_kics_hostseq_split$cohort)
  cohort_groups <- cohort_groups[!is.na(cohort_groups)]

  if (length(cohort_groups) >= 2) {
    for (grp in cohort_groups) {
      grp_data <- te_kics_hostseq_split %>% filter(cohort == grp)
      if (nrow(grp_data) > 0) {
        run_lm_gsea_analysis(
          rna_data = rna_filtered,
          te_split_data = grp_data,
          analysis_name = paste0("kics_hostseq_", tolower(grp)),
          base_dir = rna_dir,
          cohort_name = grp,
          actual_matched_samples = length(matched_rna_samples),
          group_by_gene = TRUE,
          matched_samples = matched_samples_for_lm
        )
      }
    }
  }
}

#### GSEA ANALYSIS 6: TP53 STATUS ####
write_output(quote(NULL), "GSEA Analysis - TP53 Status")

tp53_groups <- unique(te_aff_split$TP53_status)
tp53_groups <- tp53_groups[!is.na(tp53_groups)]

if (length(tp53_groups) >= 2) {
  for (grp in tp53_groups) {
    grp_data <- te_aff_split %>% filter(TP53_status == grp)
    if (nrow(grp_data) > 0) {
      run_lm_gsea_analysis(
        rna_data = rna_filtered,
        te_split_data = grp_data,
        analysis_name = paste0("tp53_status_", tolower(grp)),
        base_dir = rna_dir,
        cohort_name = paste0("TP53_", grp),
        actual_matched_samples = length(matched_rna_samples),
        group_by_gene = TRUE,
        matched_samples = matched_samples_for_lm
      )
    }
  }
}

#### GSEA ANALYSIS 7: TP53 3-LEVEL ####
write_output(quote(NULL), "GSEA Analysis - TP53 3-Level")

tp53_3level_groups <- unique(te_aff_split$TP53_3level)
tp53_3level_groups <- tp53_3level_groups[!is.na(tp53_3level_groups)]

if (length(tp53_3level_groups) >= 2) {
  for (grp in tp53_3level_groups) {
    grp_data <- te_aff_split %>% filter(TP53_3level == grp)
    if (nrow(grp_data) > 0) {
      run_lm_gsea_analysis(
        rna_data = rna_filtered,
        te_split_data = grp_data,
        analysis_name = paste0("tp53_3level_", tolower(grp)),
        base_dir = rna_dir,
        cohort_name = paste0("TP53_3level_", grp),
        actual_matched_samples = length(matched_rna_samples),
        group_by_gene = TRUE,
        matched_samples = matched_samples_for_lm
      )
    }
  }
}

#### GSEA ANALYSIS 8: LFS BY TP53 STATUS ####
write_output(quote(NULL), "GSEA Analysis - LFS by TP53 Status")

if (exists("te_lfs_split") && "TP53_status" %in% colnames(te_lfs_split)) {
  lfs_tp53_groups <- unique(te_lfs_split$TP53_status)
  lfs_tp53_groups <- lfs_tp53_groups[!is.na(lfs_tp53_groups)]

  if (length(lfs_tp53_groups) >= 2) {
    for (grp in lfs_tp53_groups) {
      grp_data <- te_lfs_split %>% filter(TP53_status == grp)
      if (nrow(grp_data) > 0) {
        run_lm_gsea_analysis(
          rna_data = rna_filtered,
          te_split_data = grp_data,
          analysis_name = paste0("lfs_tp53_status_", tolower(grp)),
          base_dir = rna_dir,
          cohort_name = paste0("LFS_TP53_", grp),
          actual_matched_samples = length(matched_rna_samples),
          group_by_gene = TRUE,
          matched_samples = matched_samples_for_lm
        )
      }
    }
  }
}

#### GSEA ANALYSIS 9: KICS BY TUMOR TYPE ####
write_output(quote(NULL), "GSEA Analysis - KICS by Tumor Type")

if (exists("te_kics_split") && "tumor_type" %in% colnames(te_kics_split)) {
  tumor_type_counts <- table(te_kics_split$tumor_type)
  valid_tumor_types <- names(tumor_type_counts[tumor_type_counts >= 3])

  if (length(valid_tumor_types) >= 2) {
    for (grp in valid_tumor_types) {
      grp_data <- te_kics_split %>% filter(tumor_type == grp)
      if (nrow(grp_data) > 0) {
        run_lm_gsea_analysis(
          rna_data = rna_filtered,
          te_split_data = grp_data,
          analysis_name = paste0("kics_tumor_type_", gsub("[^a-zA-Z0-9]", "_", tolower(grp))),
          base_dir = rna_dir,
          cohort_name = grp,
          actual_matched_samples = length(matched_rna_samples),
          group_by_gene = TRUE,
          matched_samples = matched_samples_for_lm
        )
      }
    }
  }
}

#### GSEA ANALYSIS 10: KICS BY ANCESTRY ####
write_output(quote(NULL), "GSEA Analysis - KICS by Ancestry")

if (exists("te_kics_split") && "predicted_ancestry_thres" %in% colnames(te_kics_split)) {
  ancestry_counts <- table(te_kics_split$predicted_ancestry_thres)
  valid_ancestries <- names(ancestry_counts[ancestry_counts >= 10])

  if (length(valid_ancestries) >= 2) {
    for (grp in valid_ancestries) {
      grp_data <- te_kics_split %>% filter(predicted_ancestry_thres == grp)
      if (nrow(grp_data) > 0) {
        run_lm_gsea_analysis(
          rna_data = rna_filtered,
          te_split_data = grp_data,
          analysis_name = paste0("kics_ancestry_", tolower(grp)),
          base_dir = rna_dir,
          cohort_name = paste0("Ancestry_", grp),
          actual_matched_samples = length(matched_rna_samples),
          group_by_gene = TRUE,
          matched_samples = matched_samples_for_lm
        )
      }
    }
  }
}

#### GSEA ANALYSIS 11: LFS BY COHORT ####
write_output(quote(NULL), "GSEA Analysis - LFS by Cohort")

if (exists("te_lfs_split") && "cohort" %in% colnames(te_lfs_split)) {
  cohorts <- unique(te_lfs_split$cohort)
  cohorts <- cohorts[!is.na(cohorts)]

  if (length(cohorts) >= 2) {
    for (grp in cohorts) {
      grp_data <- te_lfs_split %>% filter(cohort == grp)
      if (nrow(grp_data) > 0) {
        run_lm_gsea_analysis(
          rna_data = rna_filtered,
          te_split_data = grp_data,
          analysis_name = paste0("lfs_cohort_", tolower(grp)),
          base_dir = rna_dir,
          cohort_name = grp,
          actual_matched_samples = length(matched_rna_samples),
          group_by_gene = TRUE,
          matched_samples = matched_samples_for_lm
        )
      }
    }
  }
}

#### GSEA ANALYSIS 12: KICS BY SAMPLE TYPE ####
write_output(quote(NULL), "GSEA Analysis - KICS by Sample Type")

sample_type_file <- "/Users/briannelaverty/Documents/R_Malkin/clinical/kics_germline_sample_type.csv"
if (file.exists(sample_type_file) && exists("te_kics_split")) {
  kics_sample_type <- prep_kics_sample_type(sample_type_file)

  te_kics_sampletype <- merge(te_kics_split, kics_sample_type,
                               by = "sample", all.x = TRUE)

  valid_sample_types <- c("Blood", "Fibroblasts", "Tissue (fresh)")
  te_kics_sampletype <- te_kics_sampletype %>%
    filter(sample_type %in% valid_sample_types)

  sample_types <- unique(te_kics_sampletype$sample_type)
  sample_types <- sample_types[!is.na(sample_types)]

  if (length(sample_types) >= 2) {
    for (grp in sample_types) {
      grp_data <- te_kics_sampletype %>% filter(sample_type == grp)
      if (nrow(grp_data) > 0) {
        run_lm_gsea_analysis(
          rna_data = rna_filtered,
          te_split_data = grp_data,
          analysis_name = paste0("kics_sample_type_", gsub("[^a-zA-Z0-9]", "_", tolower(grp))),
          base_dir = rna_dir,
          cohort_name = grp,
          actual_matched_samples = length(matched_rna_samples),
          group_by_gene = TRUE,
          matched_samples = matched_samples_for_lm
        )
      }
    }
  }
}

#### SAVE DATA FOR RE-RNA SCRIPT ####
cat("\n===== SAVING DATA =====\n")

# Save rna_filtered for use by RE-RNA script
cat("Saving rna_filtered for RE-RNA script...\n")
saveRDS(rna_filtered, paste0(rna_files_dir, "rna_filtered_germline.rds"))
cat("Saved rna_filtered to:", paste0(rna_files_dir, "rna_filtered_germline.rds\n"))

# Save matched_samples for use by RE-RNA script
cat("Saving matched_samples for RE-RNA script...\n")
saveRDS(matched_samples_for_lm, paste0(rna_files_dir, "matched_samples_germline.rds"))
cat("Saved matched_samples to:", paste0(rna_files_dir, "matched_samples_germline.rds\n"))

# Save run summary
run_summary_file <- paste0(rna_dir, "RUN_SUMMARY.csv")
write.csv(run_summary, run_summary_file, row.names = FALSE)
cat("Saved run summary to:", run_summary_file, "\n")

#### SUMMARY ####
write_output(quote(NULL), "RNA Analysis Summary")

cat("Analysis completed:\n")
cat("  1. RNA data processing and sample matching\n")
cat("  2. Gene expression plots for", gene_of_interest, "\n")
cat("  3. LM-based TE expression analysis (replaces Wilcoxon)\n")
cat("  4. Multi-database GSEA pathway analysis (t-statistics ranking)\n\n")

cat("Analyses run:\n")
cat("  1. Affected cohort by gene\n")
cat("  2. Affected cohort by TE coordinate\n")
cat("  3. KICS cohort by gene\n")
cat("  4. LFS cohort by gene\n")
cat("  5. KICS vs HostSeq\n")
cat("  6. TP53 status\n")
cat("  7. TP53 3-level\n")
cat("  8. LFS by TP53 status\n")
cat("  9. KICS by tumor type\n")
cat("  10. KICS by ancestry\n")
cat("  11. LFS by cohort\n")
cat("  12. KICS by sample type\n\n")

cat("Databases used:", paste(DATABASES_TO_RUN, collapse = ", "), "\n")
cat("Min samples:", paste(MIN_SAMPLES_VALUES, collapse = ", "), "\n")
cat("Q-pathway thresholds:", paste(Q_PATHWAY_VALUES, collapse = ", "), "\n")

cat("\nOutput structure: {database}/{analysis}/files/ and plots\n")
cat("Output location:", rna_dir, "\n")

cat("\n Script completed successfully\n")

# Close module-specific sink
close_module_sink()
