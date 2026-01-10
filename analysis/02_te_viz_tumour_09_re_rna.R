#!/usr/bin/env Rscript

# Tumour TE Visualization - RE-RNA Differential Expression (Multi-Database GSEA)
# LM-based expression testing with t-statistics for GSEA ranking
# Databases: GO_BP, Reactome, MSigDB_Hallmark, MSigDB_Oncogenic (no KEGG)

# Source common setup and load data
source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/00_viz_common_setup.R")
REQUIRED_DATA <- c("count_matrix", "expand", "split", "rna", "clinical")
source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/00_viz_load_data_tumour.R")

# Create base output directory (top-level re_rna/)
re_rna_dir <- paste0(plot_dir, "re_rna/")
dir.create(re_rna_dir, showWarnings = FALSE, recursive = TRUE)

# Initialize module-specific text output
init_module_sink(re_rna_dir, "RE_RNA")

cat("Running 02_te_viz_tumour_09_re_rna.R...\n")
cat("LM-based expression testing with t-statistics for GSEA ranking\n")
cat("Databases: GO_BP, Reactome, MSigDB_Hallmark, MSigDB_Oncogenic\n\n")

# Initialize run summary
run_summary <- data.frame(
  analysis_name = character(),
  cohort = character(),
  n_re_samples = integer(),
  n_rna_samples = integer(),
  n_matched_samples = integer(),
  minsample = integer(),
  genes_tested = integer(),
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

cat("Base output directory:", re_rna_dir, "\n\n")

#### LOAD RNA DATA ####
write_output(quote(NULL), "Loading RNA Data")

# Try to load rna_filtered from RNA script
rna_filtered_file <- paste0(plot_dir, "rna/files/rna_filtered_tumour.rds")
if (file.exists(rna_filtered_file)) {
  cat("Loading rna_filtered from RNA script...\n")
  rna_filtered <- readRDS(rna_filtered_file)
  cat("Successfully loaded rna_filtered:", nrow(rna_filtered), "genes,", ncol(rna_filtered) - 1, "samples\n")

  # Load matched_samples from RNA script (for proper sample matching)
  matched_samples_file <- paste0(plot_dir, "rna/files/matched_samples_tumour.rds")
  if (file.exists(matched_samples_file)) {
    cat("Loading matched_samples from RNA script...\n")
    matched_samples_for_lm <- readRDS(matched_samples_file)
    cat("Successfully loaded matched_samples:", nrow(matched_samples_for_lm), "sample pairs\n\n")
  } else {
    cat("Warning: matched_samples file not found. Sample matching may fail.\n")
    cat("Run 02_te_viz_tumour_07_rna.R first to generate it.\n\n")
    matched_samples_for_lm <- NULL
  }
} else {
  # Fallback: process RNA data here
  cat("RNA file not found, processing RNA data...\n")
  rna_processing_result <- process_tumor_rna_data(kics_rna, lfs_rna, stjude_rna, lfs_wgs2rna)
  rna_ready_for_filtering <- rna_processing_result$rna_data

  # Match samples
  te_samples <- unique(te_all_t$sample)
  tumour_te_base_ids <- unique(sub("_.*", "", te_samples))
  rna_samples <- setdiff(colnames(rna_ready_for_filtering), "gene_name")
  rna_base_ids <- sub("_.*", "", rna_samples)
  matched_rna_samples <- rna_samples[rna_base_ids %in% tumour_te_base_ids]

  rna_filtered <- rna_ready_for_filtering %>%
    select(gene_name, all_of(intersect(matched_rna_samples, colnames(rna_ready_for_filtering))))

  cat("Processed rna_filtered:", nrow(rna_filtered), "genes,", ncol(rna_filtered) - 1, "samples\n\n")

  # No matched_samples in fallback case
  matched_samples_for_lm <- NULL
}

#### LOAD RE DATA ####
write_output(quote(NULL), "Loading RE Data")

# Try to load RE data from RE script
re_files_dir <- paste0(plot_dir, "reg_element/files/")
te_kics_re_file <- paste0(re_files_dir, "te_kics_re_split_t.rds")
te_aff_re_file <- paste0(re_files_dir, "te_aff_re_split_t.rds")

if (file.exists(te_kics_re_file)) {
  cat("Loading RE data from RE script...\n")
  te_kics_re_split <- readRDS(te_kics_re_file)
  cat("KICS RE data:", nrow(te_kics_re_split), "rows\n")
} else {
  # Fallback: process RE data here
  re_report_path_tumour <- "/Users/briannelaverty/Documents/R_Malkin/te/data/final/tumour_annotSV_output.SV_RE_intersect.report"
  cat("Loading RE data from report...\n")
  re_kics_preprocessed <- preprocess_re_for_rna(
    re_report_path = re_report_path_tumour,
    valid_samples = unique(te_kics_t$sample)
  )
  te_kics_re_split <- re_kics_preprocessed$re_split
  cat("KICS RE data:", nrow(te_kics_re_split), "rows\n")
}

if (file.exists(te_aff_re_file)) {
  te_aff_re_split <- readRDS(te_aff_re_file)
  cat("Affected RE data:", nrow(te_aff_re_split), "rows\n")
} else {
  re_report_path_tumour <- "/Users/briannelaverty/Documents/R_Malkin/te/data/final/tumour_annotSV_output.SV_RE_intersect.report"
  re_aff_preprocessed <- preprocess_re_for_rna(
    re_report_path = re_report_path_tumour,
    valid_samples = unique(te_aff_t$sample)
  )
  te_aff_re_split <- re_aff_preprocessed$re_split
  cat("Affected RE data:", nrow(te_aff_re_split), "rows\n")
}

cat("\n")

# Add TP53 3-level classification
te_all_t <- add_tp53_3level(te_all_t)
te_aff_t <- add_tp53_3level(te_aff_t)
te_kics_re_split <- add_tp53_3level(te_kics_re_split)
te_aff_re_split <- add_tp53_3level(te_aff_re_split)

#### HELPER FUNCTION: Run LM-based GSEA analysis for RE-RNA ####
run_lm_gsea_analysis_re_rna <- function(rna_data,
                                         re_split_data,
                                         analysis_name,
                                         base_dir,
                                         cohort_name = "Unknown",
                                         actual_matched_samples = 0,
                                         databases = DATABASES_TO_RUN,
                                         gene_col = "gene_reg",
                                         matched_samples = NULL) {

  cat("\n========================================\n")
  cat("  LM-GSEA Analysis (RE-RNA):", analysis_name, "\n")
  cat("========================================\n")

  # Get sample counts
  sample_col <- if ("sample" %in% colnames(re_split_data)) "sample" else
                if ("sample.x" %in% colnames(re_split_data)) "sample.x" else "sample_id"
  n_re_samples <- length(unique(re_split_data[[sample_col]]))
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

    # Run LM-based expression testing for RE-affected genes
    tryCatch({
      lm_results <- test_te_expression_effects_lm(
        rna_df = rna_data,
        te_split_df = re_split_data,
        min_samples = minsample,
        group_by_gene = TRUE,
        gene_col = "gene_reg",  # RE data uses gene_reg instead of Gene_name
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
          n_re_samples = n_re_samples,
          n_rna_samples = n_rna_samples,
          n_matched_samples = actual_matched_samples,
          minsample = minsample,
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
          n_re_samples = n_re_samples,
          n_rna_samples = n_rna_samples,
          n_matched_samples = actual_matched_samples,
          minsample = minsample,
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
                  paste0(files_dir, "re_rna_gsea_ranked_genes_", db, "_", analysis_name,
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
          gsea_rds_file <- paste0(files_dir, "re_rna_gsea_result_", db, "_", analysis_name,
                                   "_minsample", minsample, ".rds")
          saveRDS(gsea_result, gsea_rds_file)
          cat("    Saved GSEA result object:", basename(gsea_rds_file), "\n")

          # Save full GSEA CSV (all results, not filtered)
          full_gsea_df <- as.data.frame(gsea_result)
          full_csv_file <- paste0(files_dir, "re_rna_gsea_full_", db, "_", analysis_name,
                                   "_minsample", minsample, ".csv")
          write.csv(full_gsea_df, full_csv_file, row.names = FALSE)
          cat("    Saved full GSEA CSV:", basename(full_csv_file), "(", nrow(full_gsea_df), "terms)\n")

          # Save unfiltered plots (using all results)
          unfiltered_prefix <- paste0("re_rna_gsea_", db, "_", analysis_name,
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
            prefix <- paste0("re_rna_gsea_", db, "_", analysis_name,
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
        n_re_samples = n_re_samples,
        n_rna_samples = n_rna_samples,
        n_matched_samples = actual_matched_samples,
        minsample = minsample,
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
        n_re_samples = n_re_samples,
        n_rna_samples = n_rna_samples,
        n_matched_samples = actual_matched_samples,
        minsample = minsample,
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

#### GSEA ANALYSIS 1: KICS COHORT ####
write_output(quote(NULL), "GSEA Analysis - KICS Cohort (RE-RNA)")

run_lm_gsea_analysis_re_rna(
  rna_data = rna_filtered,
  re_split_data = te_kics_re_split,
  analysis_name = "kics",
  base_dir = re_rna_dir,
  cohort_name = "KICS",
  actual_matched_samples = ncol(rna_filtered) - 1,
  matched_samples = matched_samples_for_lm
)

#### GSEA ANALYSIS 2: AFFECTED COHORT ####
write_output(quote(NULL), "GSEA Analysis - Affected Cohort (RE-RNA)")

run_lm_gsea_analysis_re_rna(
  rna_data = rna_filtered,
  re_split_data = te_aff_re_split,
  analysis_name = "affected",
  base_dir = re_rna_dir,
  cohort_name = "Affected",
  actual_matched_samples = ncol(rna_filtered) - 1,
  matched_samples = matched_samples_for_lm
)

#### TP53-STRATIFIED RE-RNA ANALYSIS ####
write_output(quote(NULL), "TP53-Stratified RE-RNA Analysis")

cat("\n===== TP53-STRATIFIED RE-RNA ANALYSIS =====\n")

cat("TP53_3level distribution:\n")
print(table(te_all_t$TP53_3level, useNA = "always"))
cat("\n")

# Get TP53 groups from te_all_t (tumour data)
tp53_groups <- list(
  "TP53germline" = te_all_t %>% filter(TP53_3level == "Germline") %>% pull(sample) %>% unique(),
  "TP53somatic" = te_all_t %>% filter(TP53_3level == "Somatic") %>% pull(sample) %>% unique(),
  "TP53wt" = te_all_t %>% filter(TP53_3level == "WT") %>% pull(sample) %>% unique()
)

cat("TP53 group sample counts:\n")
for (grp in names(tp53_groups)) {
  cat("  ", grp, ":", length(tp53_groups[[grp]]), "samples\n")
}
cat("\n")

for (tp53_group in names(tp53_groups)) {
  samples_in_group <- tp53_groups[[tp53_group]]

  cat("\n--- Analyzing TP53 group:", tp53_group, "---\n")
  cat("Samples in group:", length(samples_in_group), "\n")

  if (length(samples_in_group) < 5) {
    cat("Skipping", tp53_group, "- too few samples\n")
    next
  }

  # Filter RE data to these samples
  # Need to handle different sample column names
  if ("sample" %in% colnames(te_aff_re_split)) {
    re_tp53_split <- te_aff_re_split %>%
      filter(sample %in% samples_in_group)
  } else if ("sample.x" %in% colnames(te_aff_re_split)) {
    re_tp53_split <- te_aff_re_split %>%
      filter(sample.x %in% samples_in_group)
    re_tp53_split$sample <- re_tp53_split$sample.x
  } else if ("sample_id" %in% colnames(te_aff_re_split)) {
    re_tp53_split <- te_aff_re_split %>%
      filter(sample_id %in% samples_in_group)
    re_tp53_split$sample <- re_tp53_split$sample_id
  } else {
    cat("Could not find sample column in RE data. Skipping.\n")
    next
  }

  if (nrow(re_tp53_split) == 0) {
    cat("No RE data for", tp53_group, ". Skipping.\n")
    next
  }

  cat("RE entries in group:", nrow(re_tp53_split), "\n")

  run_lm_gsea_analysis_re_rna(
    rna_data = rna_filtered,
    re_split_data = re_tp53_split,
    analysis_name = tp53_group,
    base_dir = re_rna_dir,
    cohort_name = tp53_group,
    actual_matched_samples = ncol(rna_filtered) - 1,
    matched_samples = matched_samples_for_lm
  )
}

#### TP53 COMBINED ANALYSIS ####
write_output(quote(NULL), "TP53 Combined Analysis")

cat("\n===== TP53 COMBINED ANALYSIS =====\n")

# Load LM results from saved files for combined boxplot
tp53_de_results <- list()
for (tp53_group in names(tp53_groups)) {
  for (min_samples in MIN_SAMPLES_VALUES) {
    # The run_lm_gsea_analysis_re_rna saves to {base_dir}/files/
    lm_file <- paste0(re_rna_dir, "files/lm_results_", tp53_group, "_minsample", min_samples, ".csv")
    if (file.exists(lm_file)) {
      lm_results <- read.csv(lm_file)
      if (nrow(lm_results) > 0) {
        tp53_de_results[[paste0(tp53_group, "_min", min_samples)]] <- lm_results
        cat("Loaded:", basename(lm_file), "-", nrow(lm_results), "genes\n")
      }
    }
  }
}

# Create combined TP53 x TE boxplots
if (length(tp53_de_results) > 0) {
  cat("\nCreating combined TP53 x TE boxplots...\n")

  for (min_samples in MIN_SAMPLES_VALUES) {
    for (q_thresh in Q_PATHWAY_VALUES) {
      combined_boxplot_path <- paste0(re_rna_dir, "re_rna_de_boxplot_TP53combined_min", min_samples, "_q", q_thresh, ".png")
      tryCatch({
        plot_de_genes_tp53_combined_boxplot(
          de_results_list = tp53_de_results,
          rna_data = rna_filtered,
          tp53_sample_groups = tp53_groups,
          output_path = combined_boxplot_path,
          p_cutoff = q_thresh,
          min_samples_key = paste0("min", min_samples),
          use_raw_p = TRUE,
          max_genes = 10
        )
        cat("  Saved:", basename(combined_boxplot_path), "\n")
      }, error = function(e) {
        cat("  Warning: Could not create combined boxplot for min", min_samples, "q", q_thresh, ":", e$message, "\n")
      })
    }
  }
} else {
  cat("No TP53 DE results available for combined boxplot\n")
}

#### TUMOR TYPE-SPECIFIC RE-RNA ANALYSIS ####
write_output(quote(NULL), "Tumor Type-Specific RE-RNA Analysis")

cat("\n===== TUMOR TYPE-SPECIFIC RE-RNA GSEA =====\n")

# Get tumor types from te_all_t
if ("tumor_type" %in% colnames(te_all_t)) {
  tumor_types <- unique(te_all_t$tumor_type)
  tumor_types <- tumor_types[!is.na(tumor_types) & tumor_types != "" & tumor_types != "U"]
  cat("Tumor types available:", paste(tumor_types, collapse = ", "), "\n\n")

  for (tumor_type in tumor_types) {
    cat("\n--- Analyzing tumor type:", tumor_type, "---\n")

    # Get samples for this tumor type
    tumor_type_samples <- unique(te_all_t$sample[te_all_t$tumor_type == tumor_type])
    cat("Samples:", length(tumor_type_samples), "\n")

    if (length(tumor_type_samples) < 5) {
      cat("Too few samples for tumor type:", tumor_type, ". Skipping.\n")
      next
    }

    # Filter RE data (KICS cohort) for tumor type samples
    if ("sample" %in% colnames(te_kics_re_split)) {
      re_split_tumor_type <- te_kics_re_split %>%
        filter(sample %in% tumor_type_samples)
    } else if ("sample.x" %in% colnames(te_kics_re_split)) {
      re_split_tumor_type <- te_kics_re_split %>%
        filter(sample.x %in% tumor_type_samples)
      re_split_tumor_type$sample <- re_split_tumor_type$sample.x
    } else if ("sample_id" %in% colnames(te_kics_re_split)) {
      re_split_tumor_type <- te_kics_re_split %>%
        filter(sample_id %in% tumor_type_samples)
      re_split_tumor_type$sample <- re_split_tumor_type$sample_id
    } else {
      cat("Could not find sample column in RE data. Skipping.\n")
      next
    }

    if (nrow(re_split_tumor_type) == 0) {
      cat("No RE data for tumor type:", tumor_type, ". Skipping.\n")
      next
    }

    cat("RE entries in group:", nrow(re_split_tumor_type), "\n")

    # Clean tumor type name for file naming
    tumor_type_clean <- gsub(" ", "_", tumor_type)
    tumor_type_clean <- gsub("[^A-Za-z0-9_]", "", tumor_type_clean)

    run_lm_gsea_analysis_re_rna(
      rna_data = rna_filtered,
      re_split_data = re_split_tumor_type,
      analysis_name = paste0("tumor_", tumor_type_clean),
      base_dir = re_rna_dir,
      cohort_name = tumor_type,
      actual_matched_samples = ncol(rna_filtered) - 1,
      matched_samples = matched_samples_for_lm
    )
  }
} else {
  cat("Warning: tumor_type column not found in TE data. Skipping tumor type-specific analysis.\n")
}

#### ANCESTRY-STRATIFIED RE-RNA ANALYSIS ####
write_output(quote(NULL), "Ancestry-Stratified RE-RNA Analysis")

cat("\n===== ANCESTRY-STRATIFIED RE-RNA ANALYSIS =====\n")

# Get ancestry groups from te_kics_t (if available)
if ("ancestry" %in% colnames(te_kics_t)) {
  ancestry_groups <- unique(te_kics_t$ancestry)
  ancestry_groups <- ancestry_groups[!is.na(ancestry_groups) & ancestry_groups != ""]
  cat("Ancestry groups available:", paste(ancestry_groups, collapse = ", "), "\n\n")

  for (ancestry in ancestry_groups) {
    cat("\n--- Analyzing ancestry:", ancestry, "---\n")

    # Get samples for this ancestry
    ancestry_samples <- unique(te_kics_t$sample[te_kics_t$ancestry == ancestry])
    cat("Samples:", length(ancestry_samples), "\n")

    if (length(ancestry_samples) < 5) {
      cat("Too few samples for ancestry:", ancestry, ". Skipping.\n")
      next
    }

    # Filter RE data for ancestry samples
    if ("sample" %in% colnames(te_kics_re_split)) {
      re_split_ancestry <- te_kics_re_split %>%
        filter(sample %in% ancestry_samples)
    } else if ("sample.x" %in% colnames(te_kics_re_split)) {
      re_split_ancestry <- te_kics_re_split %>%
        filter(sample.x %in% ancestry_samples)
      re_split_ancestry$sample <- re_split_ancestry$sample.x
    } else if ("sample_id" %in% colnames(te_kics_re_split)) {
      re_split_ancestry <- te_kics_re_split %>%
        filter(sample_id %in% ancestry_samples)
      re_split_ancestry$sample <- re_split_ancestry$sample_id
    } else {
      cat("Could not find sample column in RE data. Skipping.\n")
      next
    }

    if (nrow(re_split_ancestry) == 0) {
      cat("No RE data for ancestry:", ancestry, ". Skipping.\n")
      next
    }

    cat("RE entries in group:", nrow(re_split_ancestry), "\n")

    # Clean ancestry name for file naming
    ancestry_clean <- gsub(" ", "_", ancestry)
    ancestry_clean <- gsub("[^A-Za-z0-9_]", "", ancestry_clean)

    run_lm_gsea_analysis_re_rna(
      rna_data = rna_filtered,
      re_split_data = re_split_ancestry,
      analysis_name = paste0("ancestry_", ancestry_clean),
      base_dir = re_rna_dir,
      cohort_name = ancestry,
      actual_matched_samples = ncol(rna_filtered) - 1,
      matched_samples = matched_samples_for_lm
    )
  }
} else {
  cat("Warning: ancestry column not found in TE data. Skipping ancestry-specific analysis.\n")
}

#### LFS COHORT RE-RNA ANALYSIS ####
write_output(quote(NULL), "LFS Cohort RE-RNA Analysis")

cat("\n===== LFS COHORT RE-RNA ANALYSIS =====\n")

# Check if we have LFS samples in the tumour data
if (exists("te_lfs_t") && nrow(te_lfs_t) > 0) {
  lfs_samples <- unique(te_lfs_t$sample)
  cat("LFS samples:", length(lfs_samples), "\n")

  if (length(lfs_samples) >= 5) {
    # Filter RE data for LFS samples
    if ("sample" %in% colnames(te_kics_re_split)) {
      re_split_lfs <- te_kics_re_split %>%
        filter(sample %in% lfs_samples)
    } else if ("sample.x" %in% colnames(te_kics_re_split)) {
      re_split_lfs <- te_kics_re_split %>%
        filter(sample.x %in% lfs_samples)
      re_split_lfs$sample <- re_split_lfs$sample.x
    } else if ("sample_id" %in% colnames(te_kics_re_split)) {
      re_split_lfs <- te_kics_re_split %>%
        filter(sample_id %in% lfs_samples)
      re_split_lfs$sample <- re_split_lfs$sample_id
    } else {
      cat("Could not find sample column in RE data.\n")
      re_split_lfs <- data.frame()
    }

    if (nrow(re_split_lfs) > 0) {
      cat("RE entries in LFS:", nrow(re_split_lfs), "\n")

      run_lm_gsea_analysis_re_rna(
        rna_data = rna_filtered,
        re_split_data = re_split_lfs,
        analysis_name = "lfs",
        base_dir = re_rna_dir,
        cohort_name = "LFS",
        actual_matched_samples = ncol(rna_filtered) - 1,
        matched_samples = matched_samples_for_lm
      )
    } else {
      cat("No RE data for LFS samples.\n")
    }
  } else {
    cat("Too few LFS samples for analysis.\n")
  }
} else {
  cat("Warning: No LFS tumour data available. Skipping LFS-specific analysis.\n")
}

# Save run summary
run_summary_file <- paste0(re_rna_dir, "RUN_SUMMARY.csv")
write.csv(run_summary, run_summary_file, row.names = FALSE)
cat("Saved run summary to:", run_summary_file, "\n")

#### SUMMARY ####
write_output(quote(NULL), "RE-RNA Analysis Summary")

cat("Analysis completed:\n")
cat("  1. KICS cohort RE-RNA analysis (LM-based)\n")
cat("  2. Affected cohort RE-RNA analysis (LM-based)\n")
cat("  3. TP53-stratified RE-RNA analysis (Germline/Somatic/WT)\n")
cat("  4. Tumor type-specific RE-RNA analysis\n")
cat("  5. Ancestry-stratified RE-RNA analysis\n")
cat("  6. LFS cohort RE-RNA analysis\n\n")

cat("Databases used:", paste(DATABASES_TO_RUN, collapse = ", "), "\n")
cat("Min samples:", paste(MIN_SAMPLES_VALUES, collapse = ", "), "\n")
cat("Q-pathway thresholds:", paste(Q_PATHWAY_VALUES, collapse = ", "), "\n")

cat("\nOutput structure: {database}/{analysis}/files/ and plots\n")
cat("Output location:", re_rna_dir, "\n")

cat("\n Script completed successfully\n")

# Close module-specific sink
close_module_sink()
