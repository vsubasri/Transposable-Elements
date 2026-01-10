#!/usr/bin/env Rscript

# Germline TE Visualization - Pathway Analysis (Multi-Database)
# GLM-based differential incidence testing with covariate control
# Three analysis types: ora/, glm_ora/, glm_gsea/
# Databases: GO_BP, Reactome, MSigDB_Hallmark, MSigDB_Oncogenic (no KEGG)

# Source common setup and load data
source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/00_viz_common_setup.R")
REQUIRED_DATA <- c("split", "clinical")
source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/00_viz_load_data_germline.R")

# Load emmeans for post-hoc contrasts
library(emmeans)

# Initialize module-specific text output
init_module_sink(paste0(plot_dir, "pathway/"), "PATHWAY")

# Initialize run summary with analysis_type column for ora/glm_ora/glm_gsea
run_summary <- data.frame(
  analysis_name = character(),
  analysis_type = character(),        # "ora", "glm_ora", or "glm_gsea"
  group_col = character(),
  any_pathway_significant = logical(),
  pathways_significant = character(),
  n_samples_per_group = character(),
  genes_tested = integer(),
  genes_significant = integer(),
  enrichment_ran = logical(),
  minsample = integer(),
  qgene = numeric(),
  mingene = integer(),
  qpathway = numeric(),
  groups = character(),
  glm_errors = integer(),
  posthoc_errors = integer(),
  posthoc_skipped = character(),
  stringsAsFactors = FALSE
)

cat("Running 02_te_viz_germline_05_pathway.R...\n")
cat("Multi-database pathway analysis with GLM-based gene filtering\n")
cat("Databases: GO_BP, Reactome, Hallmark, Oncogenic\n\n")

#### CONFIGURATION ####

# Databases to run (no KEGG per plan)
DATABASES_TO_RUN <- c("GO_BP", "Reactome", "Hallmark", "Oncogenic")

# Parameter settings (per plan)
MIN_SAMPLES_VALUES <- c(3, 5)
Q_GENE_VALUES <- c(0.05, 0.1, 0.25)  # For GLM gene filtering
MIN_GENE_VALUES <- c(3, 5)           # Min genes hitting pathway
Q_PATHWAY_VALUES <- c(0.05, 0.1, 0.25)  # For pathway filtering

# Load background genes
background_genes_file <- "/Users/briannelaverty/Documents/R_Malkin/TE/data/final/background_genes.txt"
if (file.exists(background_genes_file)) {
  BACKGROUND_GENES <- readLines(background_genes_file)
  cat("Loaded", length(BACKGROUND_GENES), "background genes for ORA universe\n")
} else {
  BACKGROUND_GENES <- NULL
  cat("Warning: Background genes file not found. ORA will use default universe.\n")
}

# Create base output directories for three analysis types
pathway_dir <- paste0(plot_dir, "pathway/")
ora_dir <- paste0(pathway_dir, "ora/")
glm_ora_dir <- paste0(pathway_dir, "glm_ora/")
glm_gsea_dir <- paste0(pathway_dir, "glm_gsea/")

dir.create(pathway_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(ora_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(glm_ora_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(glm_gsea_dir, showWarnings = FALSE, recursive = TRUE)

cat("Base output directory:", pathway_dir, "\n")
cat("  - ora/: Simple ORA (minsample filter only)\n")
cat("  - glm_ora/: GLM-based ORA\n")
cat("  - glm_gsea/: GLM-based GSEA\n\n")

#### ADD 3-LEVEL TP53 CLASSIFICATION ####
cat("Adding 3-level TP53 classification (Germline/Somatic/WT)...\n")
te_aff_split <- add_tp53_3level(te_aff_split)
cat("TP53_3level distribution:\n")
print(table(te_aff_split$TP53_3level, useNA = "always"))
cat("\n")

if (exists("te_kics_split")) {
  te_kics_split <- add_tp53_3level(te_kics_split)
}
if (exists("te_kics_hostseq") && nrow(te_kics_hostseq) > 0) {
  te_kics_hostseq <- add_tp53_3level(te_kics_hostseq)
}
if (exists("te_lfs_split")) {
  te_lfs_split <- add_tp53_3level(te_lfs_split)
}

#### HELPER FUNCTIONS ####

# Helper: Create ORA visualizations for a filtered result
create_ora_visualizations <- function(ora_result, filtered_df, te_data, db, prefix, db_dir,
                                       gene_col = "Gene_name", sample_col = "sample") {
  n_terms <- nrow(filtered_df)
  if (n_terms == 0) return(NULL)

  # Filter the enrichResult object to only include pathways in filtered_df
  filtered_result <- ora_result
  if ("ID" %in% colnames(filtered_df)) {
    filtered_ids <- unique(filtered_df$ID)
    if (inherits(ora_result, "enrichResult")) {
      filtered_result@result <- filtered_result@result %>%
        filter(ID %in% filtered_ids)
    } else if (inherits(ora_result, "compareClusterResult")) {
      filtered_result@compareClusterResult <- filtered_result@compareClusterResult %>%
        filter(ID %in% filtered_ids)
    }
  }

  tryCatch({
    # 1. Dotplot
    p_dot <- enrichplot::dotplot(filtered_result, showCategory = min(20, n_terms))
    ggsave(paste0(db_dir, prefix, "_dot.png"), p_dot, width = 10, height = 8)

    # 2. Cnetplot (gene-concept network)
    if (inherits(filtered_result, "enrichResult")) {
      tryCatch({
        p_cnet <- enrichplot::cnetplot(filtered_result, showCategory = min(10, n_terms),
                                        categorySize = "pvalue")
        ggsave(paste0(db_dir, prefix, "_cnet.png"), p_cnet, width = 12, height = 10)
      }, error = function(e) cat("    Cnetplot error:", e$message, "\n"))
    }

    # 3. Emapplot (enrichment map) - needs ≥5 terms
    if (n_terms >= 5 && inherits(filtered_result, "enrichResult")) {
      tryCatch({
        ora_pairwise <- enrichplot::pairwise_termsim(filtered_result)
        p_emap <- enrichplot::emapplot(ora_pairwise, showCategory = min(30, n_terms))
        ggsave(paste0(db_dir, prefix, "_emap.png"), p_emap, width = 12, height = 10)
      }, error = function(e) cat("    Emapplot error:", e$message, "\n"))
    }

    # 4. Heatplot
    if (inherits(filtered_result, "enrichResult")) {
      tryCatch({
        p_heat <- enrichplot::heatplot(filtered_result, showCategory = min(20, n_terms))
        ggsave(paste0(db_dir, prefix, "_heatplot.png"), p_heat, width = 14, height = 8)
      }, error = function(e) cat("    Heatplot error:", e$message, "\n"))
    }

    # 5-7. Sample-level visualizations
    if (!is.null(te_data) && n_terms > 0) {
      tryCatch({
        p_balloon <- plot_pathway_balloon(filtered_df, te_data, max_pathways = 20,
                                           gene_col = gene_col)
        if (!is.null(p_balloon)) {
          ggsave(paste0(db_dir, prefix, "_balloon.png"), p_balloon, width = 14, height = 8)
        }
      }, error = function(e) cat("    Balloon plot error:", e$message, "\n"))

      tryCatch({
        ht <- plot_pathway_sample_heatmap(filtered_df, te_data, max_pathways = 30,
                                           gene_col = gene_col)
        if (!is.null(ht)) {
          png(paste0(db_dir, prefix, "_heatmap.png"), width = 12, height = 10, units = "in", res = 150)
          ComplexHeatmap::draw(ht)
          dev.off()
        }
      }, error = function(e) cat("    Heatmap error:", e$message, "\n"))

      tryCatch({
        p_gene <- plot_pathway_gene_heatmap(filtered_df, te_data, max_pathways = 15,
                                             max_genes = 50, gene_col = gene_col)
        if (!is.null(p_gene)) {
          ggsave(paste0(db_dir, prefix, "_gene_pathway.png"), p_gene, width = 12, height = 10)
        }
      }, error = function(e) cat("    Gene-pathway error:", e$message, "\n"))
    }
  }, error = function(e) cat("  Visualization error:", e$message, "\n"))
}

# Helper: Run compareCluster ORA for a specific database
run_compareCluster_ora_by_db <- function(gene_clusters, db, background_genes = NULL, pvalueCutoff = 1) {
  tryCatch({
    if (db == "GO_BP") {
      clusterProfiler::compareCluster(
        geneCluster = gene_clusters,
        fun = "enrichGO",
        OrgDb = org.Hs.eg.db::org.Hs.eg.db,
        keyType = "SYMBOL",
        ont = "BP",
        universe = background_genes,
        pvalueCutoff = pvalueCutoff,
        qvalueCutoff = 1
      )
    } else if (db == "Reactome") {
      gene_clusters_entrez <- lapply(gene_clusters, function(genes) {
        ids <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, genes, "ENTREZID", "SYMBOL")
        ids[!is.na(ids)]
      })
      gene_clusters_entrez <- gene_clusters_entrez[sapply(gene_clusters_entrez, length) > 0]

      if (length(gene_clusters_entrez) >= 2) {
        clusterProfiler::compareCluster(
          geneCluster = gene_clusters_entrez,
          fun = ReactomePA::enrichPathway,
          pvalueCutoff = pvalueCutoff,
          qvalueCutoff = 1
        )
      } else NULL
    } else if (db %in% c("Hallmark", "Oncogenic")) {
      msig_category <- if (db == "Hallmark") "H" else "C6"
      msig_db <- msigdbr::msigdbr(species = "Homo sapiens", category = msig_category)
      msig_t2g <- msig_db %>% dplyr::select(gs_name, gene_symbol)

      clusterProfiler::compareCluster(
        geneCluster = gene_clusters,
        fun = "enricher",
        TERM2GENE = msig_t2g,
        universe = background_genes,
        pvalueCutoff = pvalueCutoff,
        qvalueCutoff = 1
      )
    } else NULL
  }, error = function(e) {
    cat("    compareCluster error for", db, ":", e$message, "\n")
    NULL
  })
}

#### SIMPLE ORA ANALYSIS (no GLM) ####
run_simple_ora_analysis <- function(te_data,
                                     group_col,
                                     groups,
                                     analysis_name,
                                     base_dir,
                                     databases = DATABASES_TO_RUN,
                                     background_genes = BACKGROUND_GENES,
                                     gene_col = "Gene_name",
                                     sample_col = "sample") {

  cat("\n========================================\n")
  cat("  Simple ORA Analysis:", analysis_name, "\n")
  cat("========================================\n")
  cat("Group column:", group_col, "\n")
  cat("Groups:", paste(groups, collapse = ", "), "\n")

  # Calculate sample counts
  actual_sample_col <- if (sample_col %in% colnames(te_data)) sample_col else "sample"
  sample_counts <- te_data %>%
    filter(!is.na(.data[[group_col]])) %>%
    group_by(.data[[group_col]]) %>%
    summarise(n = n_distinct(.data[[actual_sample_col]]), .groups = "drop")
  sample_counts_vec <- setNames(sample_counts$n, as.character(sample_counts[[group_col]]))
  n_samples_str <- paste(names(sample_counts_vec), sample_counts_vec, sep = ":", collapse = ", ")

  # Create directories
  analysis_base_dir <- paste0(base_dir, analysis_name, "/")
  for (db in databases) {
    db_dir <- paste0(analysis_base_dir, db, "/")
    db_files_dir <- paste0(db_dir, "files/")
    dir.create(db_files_dir, showWarnings = FALSE, recursive = TRUE)
  }

  for (minsample in MIN_SAMPLES_VALUES) {
    cat("\n--- minsample =", minsample, "---\n")

    # Get genes with >= minsample samples having TE
    gene_sample_counts <- te_data %>%
      group_by(.data[[gene_col]]) %>%
      summarise(n_samples = n_distinct(.data[[actual_sample_col]]), .groups = "drop") %>%
      filter(n_samples >= minsample)

    genes_passing <- unique(gene_sample_counts[[gene_col]])
    genes_passing <- genes_passing[!is.na(genes_passing) & genes_passing != ""]
    cat("Genes with >=", minsample, "samples:", length(genes_passing), "\n")

    if (length(genes_passing) < 5) {
      cat("Too few genes. Skipping.\n")
      for (mingene in MIN_GENE_VALUES) {
        for (qpathway in Q_PATHWAY_VALUES) {
          run_summary <<- rbind(run_summary, data.frame(
            analysis_name = analysis_name, analysis_type = "ora",
            group_col = group_col, any_pathway_significant = FALSE,
            pathways_significant = paste(sapply(databases, function(db) paste0(db, ":0")), collapse = ", "),
            n_samples_per_group = n_samples_str, genes_tested = length(genes_passing),
            genes_significant = NA_integer_, enrichment_ran = FALSE,
            minsample = minsample, qgene = NA_real_, mingene = mingene, qpathway = qpathway,
            groups = paste(groups, collapse = ","), glm_errors = NA_integer_,
            posthoc_errors = NA_integer_, posthoc_skipped = NA_character_,
            stringsAsFactors = FALSE
          ))
        }
      }
      next
    }

    # Assign genes to groups by highest incidence
    te_filtered <- te_data %>% filter(.data[[gene_col]] %in% genes_passing)
    gene_group_incidence <- te_filtered %>%
      group_by(.data[[gene_col]], .data[[group_col]]) %>%
      summarise(n_samples = n_distinct(.data[[actual_sample_col]]), .groups = "drop") %>%
      filter(!is.na(.data[[group_col]]))

    gene_assignments <- gene_group_incidence %>%
      group_by(.data[[gene_col]]) %>%
      slice_max(n_samples, n = 1, with_ties = FALSE) %>%
      ungroup()

    gene_clusters <- split(gene_assignments[[gene_col]], gene_assignments[[group_col]])
    gene_clusters <- gene_clusters[sapply(gene_clusters, length) > 0]

    if (length(gene_clusters) < 2) {
      cat("Not enough gene clusters. Skipping.\n")
      next
    }

    cat("Gene clusters:", paste(names(gene_clusters), sapply(gene_clusters, length), sep = ":", collapse = ", "), "\n")

    # Run ORA for each database
    for (db in databases) {
      cat("\n  Database:", db, "\n")
      db_dir <- paste0(analysis_base_dir, db, "/")
      db_files_dir <- paste0(db_dir, "files/")

      ora_result <- run_compareCluster_ora_by_db(gene_clusters, db, background_genes, pvalueCutoff = 1)

      if (is.null(ora_result) || nrow(as.data.frame(ora_result)) == 0) {
        cat("    No enriched terms found\n")
        next
      }

      ora_df <- as.data.frame(ora_result)
      cat("    Saved full result:", nrow(ora_df), "terms\n")

      # Save full result
      saveRDS(ora_result, paste0(db_files_dir, "ora_", db, "_", analysis_name,
                                 "_minsample", minsample, "_full.rds"))
      write.csv(ora_df, paste0(db_files_dir, "ora_", db, "_", analysis_name,
                               "_minsample", minsample, "_full.csv"), row.names = FALSE)

      # Post-hoc filter by mingene and qpathway
      for (mingene in MIN_GENE_VALUES) {
        for (q_pathway in Q_PATHWAY_VALUES) {
          filtered_df <- ora_df %>% filter(Count >= mingene, qvalue < q_pathway)
          pathways_sig <- nrow(filtered_df)

          if (pathways_sig > 0) {
            prefix <- paste0("ora_", db, "_", analysis_name, "_minsample", minsample,
                             "_mingene", mingene, "_qpathway", q_pathway)
            write.csv(filtered_df, paste0(db_files_dir, prefix, ".csv"), row.names = FALSE)
            create_ora_visualizations(ora_result, filtered_df, te_data, db, prefix, db_dir,
                                       gene_col, actual_sample_col)
            cat("    Saved:", prefix, "(", pathways_sig, "terms)\n")
          }
        }
      }
    }

    # Add run summary entries
    for (mingene in MIN_GENE_VALUES) {
      for (qpathway in Q_PATHWAY_VALUES) {
        pathway_counts <- sapply(databases, function(db) {
          db_files_dir <- paste0(analysis_base_dir, db, "/files/")
          prefix <- paste0("ora_", db, "_", analysis_name, "_minsample", minsample,
                           "_mingene", mingene, "_qpathway", qpathway)
          csv_file <- paste0(db_files_dir, prefix, ".csv")
          if (file.exists(csv_file)) nrow(read.csv(csv_file)) else 0
        })
        pathways_str <- paste(paste0(databases, ":", pathway_counts), collapse = ", ")

        run_summary <<- rbind(run_summary, data.frame(
          analysis_name = analysis_name, analysis_type = "ora",
          group_col = group_col, any_pathway_significant = any(pathway_counts > 0),
          pathways_significant = pathways_str, n_samples_per_group = n_samples_str,
          genes_tested = length(genes_passing), genes_significant = NA_integer_,
          enrichment_ran = TRUE, minsample = minsample, qgene = NA_real_,
          mingene = mingene, qpathway = qpathway, groups = paste(groups, collapse = ","),
          glm_errors = NA_integer_, posthoc_errors = NA_integer_, posthoc_skipped = NA_character_,
          stringsAsFactors = FALSE
        ))
      }
    }
  }
}

#### GLM-ORA ANALYSIS ####
run_glm_ora_analysis <- function(te_data,
                                  group_col,
                                  groups,
                                  analysis_name,
                                  base_dir,
                                  databases = DATABASES_TO_RUN,
                                  background_genes = BACKGROUND_GENES,
                                  gene_col = "Gene_name",
                                  sample_col = "sample",
                                  covariates = covar_med) {

  cat("\n========================================\n")
  cat("  GLM-ORA Analysis:", analysis_name, "\n")
  cat("========================================\n")
  cat("Group column:", group_col, "\n")
  cat("Groups:", paste(groups, collapse = ", "), "\n")

  n_groups <- length(groups)

  # Calculate sample counts
  actual_sample_col <- if (sample_col %in% colnames(te_data)) sample_col else "sample"
  sample_counts <- te_data %>%
    filter(!is.na(.data[[group_col]])) %>%
    group_by(.data[[group_col]]) %>%
    summarise(n = n_distinct(.data[[actual_sample_col]]), .groups = "drop")
  sample_counts_vec <- setNames(sample_counts$n, as.character(sample_counts[[group_col]]))
  n_samples_str <- paste(names(sample_counts_vec), sample_counts_vec, sep = ":", collapse = ", ")

  # Create directories
  analysis_base_dir <- paste0(base_dir, analysis_name, "/")
  files_dir <- paste0(analysis_base_dir, "files/")
  dir.create(files_dir, showWarnings = FALSE, recursive = TRUE)
  for (db in databases) {
    db_dir <- paste0(analysis_base_dir, db, "/")
    db_files_dir <- paste0(db_dir, "files/")
    dir.create(db_files_dir, showWarnings = FALSE, recursive = TRUE)
  }

  for (minsample in MIN_SAMPLES_VALUES) {
    cat("\n--- minsample =", minsample, "---\n")

    # Run GLM ONCE with p_threshold=1.0
    glm_result <- filter_genes_by_differential_incidence_glm(
      te_data = te_data,
      group_col = group_col,
      groups = groups,
      min_samples = minsample,
      p_threshold = 1.0,
      gene_col = gene_col,
      covariates = covariates
    )

    full_results <- glm_result$full_results
    pairwise_results <- glm_result$pairwise
    genes_tested <- glm_result$genes_attempted
    glm_errors <- glm_result$glm_errors
    posthoc_errors <- glm_result$posthoc_errors
    posthoc_skipped <- glm_result$posthoc_skipped

    # Add comparison info for 2-group
    if (!is.null(full_results) && nrow(full_results) > 0 && n_groups == 2) {
      full_results <- full_results %>%
        mutate(
          comparison = paste0(groups[1], "_vs_", groups[2]),
          enriched_in = case_when(
            z_value > 0 ~ groups[1],
            z_value < 0 ~ groups[2],
            TRUE ~ NA_character_
          )
        )
    }

    # Save full GLM results
    glm_results_file <- paste0(files_dir, "glm_", analysis_name,
                               "_minsample", minsample, "_all.csv")
    if (!is.null(full_results) && nrow(full_results) > 0) {
      write.csv(full_results, glm_results_file, row.names = FALSE)
      cat("  Saved all GLM results:", basename(glm_results_file), "(", nrow(full_results), "genes)\n")
    }

    # Loop over qgene thresholds
    for (qgene in Q_GENE_VALUES) {
      cat("\n--- minsample =", minsample, ", qgene =", qgene, "---\n")

      # Filter by qgene
      if (!is.null(full_results) && nrow(full_results) > 0) {
        p_col <- if ("lrt_padj" %in% colnames(full_results)) "lrt_padj" else "p_adj"
        sig_results <- full_results %>% filter(.data[[p_col]] < qgene)
        sig_genes <- sig_results$gene
      } else {
        sig_results <- data.frame()
        sig_genes <- character(0)
      }

      genes_significant <- length(sig_genes)
      cat("  Significant genes at q <", qgene, ":", genes_significant, "\n")

      # Save filtered GLM results
      glm_results_file <- paste0(files_dir, "glm_", analysis_name,
                                 "_minsample", minsample, "_qgene", qgene, ".csv")
      if (nrow(sig_results) > 0) {
        write.csv(sig_results, glm_results_file, row.names = FALSE)
      }

      pathway_counts_by_db <- list()

      if (length(sig_genes) < 5) {
        cat("  Too few significant genes. Skipping ORA.\n")
        for (db in databases) {
          for (mingene in MIN_GENE_VALUES) {
            for (qpathway in Q_PATHWAY_VALUES) {
              key <- paste(mingene, qpathway, sep = "_")
              if (!db %in% names(pathway_counts_by_db)) pathway_counts_by_db[[db]] <- list()
              pathway_counts_by_db[[db]][[key]] <- 0
            }
          }
        }
      } else if (n_groups >= 3) {
        # 3+ groups: use pairwise results for gene assignment
        if (!is.null(pairwise_results) && nrow(pairwise_results) > 0) {
          cat("  Using post-hoc pairwise contrasts for gene assignment\n")
          gene_assignment_df <- assign_genes_by_posthoc(pairwise_results, p_threshold = qgene)

          # Save post-hoc results
          posthoc_file <- paste0(files_dir, "glm_", analysis_name,
                                 "_posthoc_minsample", minsample, "_qgene", qgene, ".csv")
          write.csv(pairwise_results, posthoc_file, row.names = FALSE)

          if (nrow(gene_assignment_df) > 0) {
            gene_clusters <- gene_assignment_to_clusters(gene_assignment_df)

            # Save gene assignment
            gene_assign_file <- paste0(files_dir, "glm_", analysis_name,
                                       "_gene_assignment_minsample", minsample, "_qgene", qgene, ".csv")
            gene_assign_expanded <- gene_assignment_df %>%
              mutate(assigned_groups_str = sapply(assigned_groups, paste, collapse = ";"))
            write.csv(gene_assign_expanded[, c("gene", "assigned_groups_str")], gene_assign_file, row.names = FALSE)

            # Run compareCluster for each database
            for (db in databases) {
              cat("\n  Database:", db, "(compareCluster)\n")
              db_dir <- paste0(analysis_base_dir, db, "/")
              db_files_dir <- paste0(db_dir, "files/")
              pathway_counts_by_db[[db]] <- list()

              compare_result <- run_compareCluster_ora_by_db(gene_clusters, db, background_genes, pvalueCutoff = 1)

              if (!is.null(compare_result) && nrow(as.data.frame(compare_result)) > 0) {
                # Save full result
                saveRDS(compare_result, paste0(files_dir, "glm_ora_pathway_", db, "_", analysis_name,
                                               "_minsample", minsample, "_qgene", qgene, ".rds"))
                write.csv(as.data.frame(compare_result), paste0(files_dir, "glm_ora_pathway_", db, "_", analysis_name,
                                                                 "_minsample", minsample, "_qgene", qgene, ".csv"), row.names = FALSE)

                # Post-hoc filter
                for (mingene in MIN_GENE_VALUES) {
                  for (qpathway in Q_PATHWAY_VALUES) {
                    compare_df <- as.data.frame(compare_result)
                    filtered_df <- compare_df %>% filter(Count >= mingene, qvalue < qpathway)
                    pathways_sig <- nrow(filtered_df)
                    key <- paste(mingene, qpathway, sep = "_")
                    pathway_counts_by_db[[db]][[key]] <- pathways_sig

                    if (pathways_sig > 0) {
                      prefix <- paste0("glm_ora_pathway_", db, "_", analysis_name,
                                       "_minsample", minsample, "_qgene", qgene,
                                       "_mingene", mingene, "_qpathway", qpathway)
                      write.csv(filtered_df, paste0(db_files_dir, prefix, ".csv"), row.names = FALSE)
                      create_ora_visualizations(compare_result, filtered_df, te_data, db, prefix, db_dir,
                                                gene_col, actual_sample_col)
                      cat("    Saved:", prefix, "(", pathways_sig, "terms)\n")
                    }
                  }
                }
              } else {
                for (mingene in MIN_GENE_VALUES) {
                  for (qpathway in Q_PATHWAY_VALUES) {
                    key <- paste(mingene, qpathway, sep = "_")
                    pathway_counts_by_db[[db]][[key]] <- 0
                  }
                }
              }
            }
          }
        }
      } else {
        # 2-group: assign by enriched_in
        gene_clusters <- split(sig_results$gene, sig_results$enriched_in)
        gene_clusters <- gene_clusters[!is.na(names(gene_clusters))]
        gene_clusters <- gene_clusters[sapply(gene_clusters, length) > 0]

        if (length(gene_clusters) >= 1) {
          for (db in databases) {
            cat("\n  Database:", db, "\n")
            db_dir <- paste0(analysis_base_dir, db, "/")
            db_files_dir <- paste0(db_dir, "files/")
            pathway_counts_by_db[[db]] <- list()

            compare_result <- run_compareCluster_ora_by_db(gene_clusters, db, background_genes, pvalueCutoff = 1)

            if (!is.null(compare_result) && nrow(as.data.frame(compare_result)) > 0) {
              saveRDS(compare_result, paste0(files_dir, "glm_ora_pathway_", db, "_", analysis_name,
                                             "_minsample", minsample, "_qgene", qgene, ".rds"))
              write.csv(as.data.frame(compare_result), paste0(files_dir, "glm_ora_pathway_", db, "_", analysis_name,
                                                               "_minsample", minsample, "_qgene", qgene, ".csv"), row.names = FALSE)

              for (mingene in MIN_GENE_VALUES) {
                for (qpathway in Q_PATHWAY_VALUES) {
                  compare_df <- as.data.frame(compare_result)
                  filtered_df <- compare_df %>% filter(Count >= mingene, qvalue < qpathway)
                  pathways_sig <- nrow(filtered_df)
                  key <- paste(mingene, qpathway, sep = "_")
                  pathway_counts_by_db[[db]][[key]] <- pathways_sig

                  if (pathways_sig > 0) {
                    prefix <- paste0("glm_ora_pathway_", db, "_", analysis_name,
                                     "_minsample", minsample, "_qgene", qgene,
                                     "_mingene", mingene, "_qpathway", qpathway)
                    write.csv(filtered_df, paste0(db_files_dir, prefix, ".csv"), row.names = FALSE)
                    create_ora_visualizations(compare_result, filtered_df, te_data, db, prefix, db_dir,
                                              gene_col, actual_sample_col)
                    cat("    Saved:", prefix, "(", pathways_sig, "terms)\n")
                  }
                }
              }
            } else {
              for (mingene in MIN_GENE_VALUES) {
                for (qpathway in Q_PATHWAY_VALUES) {
                  key <- paste(mingene, qpathway, sep = "_")
                  pathway_counts_by_db[[db]][[key]] <- 0
                }
              }
            }
          }
        }
      }

      # Add run summary entries
      for (mingene in MIN_GENE_VALUES) {
        for (qpathway in Q_PATHWAY_VALUES) {
          key <- paste(mingene, qpathway, sep = "_")
          pathway_counts <- sapply(databases, function(db) {
            count <- pathway_counts_by_db[[db]][[key]]
            if (is.null(count)) count <- 0
            paste0(db, ":", count)
          })
          pathways_str <- paste(pathway_counts, collapse = ", ")

          run_summary <<- rbind(run_summary, data.frame(
            analysis_name = analysis_name, analysis_type = "glm_ora",
            group_col = group_col,
            any_pathway_significant = any(grepl(":[1-9]", pathway_counts)),
            pathways_significant = pathways_str, n_samples_per_group = n_samples_str,
            genes_tested = genes_tested, genes_significant = genes_significant,
            enrichment_ran = length(sig_genes) >= 5, minsample = minsample,
            qgene = qgene, mingene = mingene, qpathway = qpathway,
            groups = paste(groups, collapse = ","), glm_errors = glm_errors,
            posthoc_errors = posthoc_errors, posthoc_skipped = posthoc_skipped,
            stringsAsFactors = FALSE
          ))
        }
      }
    }
  }
}

#### GLM-GSEA ANALYSIS ####
run_glm_gsea_analysis <- function(te_data,
                                   group_col,
                                   groups,
                                   analysis_name,
                                   base_dir,
                                   databases = DATABASES_TO_RUN,
                                   gene_col = "Gene_name",
                                   sample_col = "sample",
                                   covariates = covar_med) {

  cat("\n========================================\n")
  cat("  GLM-GSEA Analysis:", analysis_name, "\n")
  cat("========================================\n")
  cat("Group column:", group_col, "\n")
  cat("Groups:", paste(groups, collapse = ", "), "\n")

  n_groups <- length(groups)

  # Calculate sample counts
  actual_sample_col <- if (sample_col %in% colnames(te_data)) sample_col else "sample"
  sample_counts <- te_data %>%
    filter(!is.na(.data[[group_col]])) %>%
    group_by(.data[[group_col]]) %>%
    summarise(n = n_distinct(.data[[actual_sample_col]]), .groups = "drop")
  sample_counts_vec <- setNames(sample_counts$n, as.character(sample_counts[[group_col]]))
  n_samples_str <- paste(names(sample_counts_vec), sample_counts_vec, sep = ":", collapse = ", ")

  # Create directories
  analysis_base_dir <- paste0(base_dir, analysis_name, "/")
  files_dir <- paste0(analysis_base_dir, "files/")
  dir.create(files_dir, showWarnings = FALSE, recursive = TRUE)
  for (db in databases) {
    db_dir <- paste0(analysis_base_dir, db, "/")
    db_files_dir <- paste0(db_dir, "files/")
    dir.create(db_files_dir, showWarnings = FALSE, recursive = TRUE)
  }

  for (minsample in MIN_SAMPLES_VALUES) {
    cat("\n--- minsample =", minsample, "---\n")

    # Run GLM ONCE with p_threshold=1.0
    glm_result <- filter_genes_by_differential_incidence_glm(
      te_data = te_data,
      group_col = group_col,
      groups = groups,
      min_samples = minsample,
      p_threshold = 1.0,
      gene_col = gene_col,
      covariates = covariates
    )

    full_results <- glm_result$full_results
    pairwise_results <- glm_result$pairwise
    genes_tested <- glm_result$genes_attempted
    glm_errors <- glm_result$glm_errors
    posthoc_errors <- glm_result$posthoc_errors
    posthoc_skipped <- glm_result$posthoc_skipped

    # Add comparison info for 2-group
    if (!is.null(full_results) && nrow(full_results) > 0 && n_groups == 2) {
      full_results <- full_results %>%
        mutate(
          comparison = paste0(groups[1], "_vs_", groups[2]),
          enriched_in = case_when(
            z_value > 0 ~ groups[1],
            z_value < 0 ~ groups[2],
            TRUE ~ NA_character_
          )
        )
    }

    # Save full GLM results
    glm_results_file <- paste0(files_dir, "glm_", analysis_name,
                               "_minsample", minsample, "_all.csv")
    if (!is.null(full_results) && nrow(full_results) > 0) {
      write.csv(full_results, glm_results_file, row.names = FALSE)
      cat("  Saved all GLM results:", basename(glm_results_file), "(", nrow(full_results), "genes)\n")
    }

    pathway_counts_by_db <- list()
    for (db in databases) {
      pathway_counts_by_db[[db]] <- list()
      for (qpathway in Q_PATHWAY_VALUES) {
        pathway_counts_by_db[[db]][[as.character(qpathway)]] <- 0
      }
    }

    gsea_ran <- FALSE

    if (is.null(full_results) || nrow(full_results) == 0) {
      cat("  No GLM results available. Skipping GSEA.\n")
      for (qpathway in Q_PATHWAY_VALUES) {
        pathways_str <- paste(sapply(databases, function(db) paste0(db, ":0")), collapse = ", ")
        run_summary <<- rbind(run_summary, data.frame(
          analysis_name = analysis_name, analysis_type = "glm_gsea",
          group_col = group_col, any_pathway_significant = FALSE,
          pathways_significant = pathways_str, n_samples_per_group = n_samples_str,
          genes_tested = genes_tested, genes_significant = NA_integer_,
          enrichment_ran = FALSE, minsample = minsample, qgene = NA_real_,
          mingene = NA_integer_, qpathway = qpathway,
          groups = paste(groups, collapse = ","), glm_errors = glm_errors,
          posthoc_errors = posthoc_errors, posthoc_skipped = posthoc_skipped,
          stringsAsFactors = FALSE
        ))
      }
      next
    }

    if (n_groups == 2) {
      # 2-group: use z_value directly
      if (!"z_value" %in% colnames(full_results)) {
        cat("  No z_value column found. Skipping GSEA.\n")
        next
      }

      gene_list <- setNames(full_results$z_value, full_results$gene)
      gene_list <- gene_list[!is.na(gene_list)]
      gene_list <- sort(gene_list, decreasing = TRUE)

      cat("  Gene list for GSEA:", length(gene_list), "genes\n")

      for (db in databases) {
        cat("\n  Database:", db, "(GSEA)\n")
        db_dir <- paste0(analysis_base_dir, db, "/")
        db_files_dir <- paste0(db_dir, "files/")

        tryCatch({
          gsea_result <- NULL

          if (db == "GO_BP") {
            gsea_result <- clusterProfiler::gseGO(
              geneList = gene_list,
              OrgDb = org.Hs.eg.db::org.Hs.eg.db,
              keyType = "SYMBOL",
              ont = "BP",
              minGSSize = 10,
              maxGSSize = 500,
              pvalueCutoff = 1,
              eps = 0
            )
          } else if (db == "Reactome") {
            gene_ids <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
                                               keys = names(gene_list),
                                               column = "ENTREZID",
                                               keytype = "SYMBOL",
                                               multiVals = "first")
            gene_list_entrez <- gene_list[!is.na(gene_ids)]
            names(gene_list_entrez) <- gene_ids[!is.na(gene_ids)]

            if (length(gene_list_entrez) >= 10) {
              gsea_result <- ReactomePA::gsePathway(
                geneList = gene_list_entrez,
                organism = "human",
                minGSSize = 10,
                maxGSSize = 500,
                pvalueCutoff = 1,
                eps = 0
              )
            }
          } else if (db %in% c("Hallmark", "Oncogenic")) {
            msig_category <- if (db == "Hallmark") "H" else "C6"
            msig_db <- msigdbr::msigdbr(species = "Homo sapiens", category = msig_category)
            msig_t2g <- msig_db %>% dplyr::select(gs_name, gene_symbol)

            gsea_result <- clusterProfiler::GSEA(
              geneList = gene_list,
              TERM2GENE = msig_t2g,
              minGSSize = 10,
              maxGSSize = 500,
              pvalueCutoff = 1,
              eps = 0
            )
          }

          if (!is.null(gsea_result) && nrow(as.data.frame(gsea_result)) > 0) {
            gsea_ran <- TRUE
            saveRDS(gsea_result, paste0(db_files_dir, "glm_gsea_pathway_", db, "_", analysis_name,
                                        "_minsample", minsample, "_full.rds"))
            write.csv(as.data.frame(gsea_result), paste0(db_files_dir, "glm_gsea_pathway_", db, "_", analysis_name,
                                                          "_minsample", minsample, "_full.csv"), row.names = FALSE)

            for (qpathway in Q_PATHWAY_VALUES) {
              gsea_df <- as.data.frame(gsea_result)
              filtered_df <- gsea_df %>% filter(qvalue < qpathway)
              pathways_sig <- nrow(filtered_df)
              pathway_counts_by_db[[db]][[as.character(qpathway)]] <- pathways_sig

              if (pathways_sig > 0) {
                prefix <- paste0("glm_gsea_pathway_", db, "_", analysis_name,
                                 "_minsample", minsample, "_qpathway", qpathway)
                write.csv(filtered_df, paste0(db_files_dir, prefix, ".csv"), row.names = FALSE)

                # GSEA visualizations
                tryCatch({
                  p_dot <- enrichplot::dotplot(gsea_result, showCategory = min(20, pathways_sig))
                  ggsave(paste0(db_dir, prefix, "_dot.png"), p_dot, width = 10, height = 8)
                }, error = function(e) NULL)

                tryCatch({
                  p_ridge <- enrichplot::ridgeplot(gsea_result, showCategory = min(15, pathways_sig))
                  ggsave(paste0(db_dir, prefix, "_ridge.png"), p_ridge, width = 10, height = 8)
                }, error = function(e) NULL)

                cat("    Saved:", prefix, "(", pathways_sig, "terms)\n")
              }
            }
          }
        }, error = function(e) {
          cat("    Error running GSEA for", db, ":", e$message, "\n")
        })
      }
    } else {
      # 3+ groups: run GSEA for each pairwise contrast
      if (!is.null(pairwise_results) && nrow(pairwise_results) > 0) {
        # Save post-hoc pairwise results
        posthoc_file <- paste0(files_dir, "glm_", analysis_name,
                               "_posthoc_minsample", minsample, ".csv")
        write.csv(pairwise_results, posthoc_file, row.names = FALSE)
        cat("  Saved post-hoc pairwise results:", basename(posthoc_file), "\n")

        contrasts <- unique(pairwise_results$contrast)
        cat("  Running GSEA for", length(contrasts), "pairwise contrasts\n")

        for (contrast_name in contrasts) {
          cat("\n  Contrast:", contrast_name, "\n")
          contrast_results <- pairwise_results %>% filter(contrast == contrast_name)

          if (!"z.ratio" %in% colnames(contrast_results)) {
            cat("    No z.ratio column found. Skipping.\n")
            next
          }

          gene_list <- setNames(contrast_results$z.ratio, contrast_results$gene)
          gene_list <- gene_list[!is.na(gene_list)]
          gene_list <- sort(gene_list, decreasing = TRUE)

          if (length(gene_list) < 10) {
            cat("    Too few genes (", length(gene_list), "). Skipping.\n")
            next
          }

          contrast_clean <- gsub(" - ", "_vs_", contrast_name)
          contrast_clean <- gsub("[^a-zA-Z0-9_]", "", contrast_clean)

          for (db in databases) {
            cat("    Database:", db, "\n")
            db_dir <- paste0(analysis_base_dir, db, "/")
            db_files_dir <- paste0(db_dir, "files/")

            tryCatch({
              gsea_result <- NULL

              if (db == "GO_BP") {
                gsea_result <- clusterProfiler::gseGO(
                  geneList = gene_list,
                  OrgDb = org.Hs.eg.db::org.Hs.eg.db,
                  keyType = "SYMBOL",
                  ont = "BP",
                  minGSSize = 10,
                  maxGSSize = 500,
                  pvalueCutoff = 1,
                  eps = 0
                )
              } else if (db == "Reactome") {
                gene_ids <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
                                                   keys = names(gene_list),
                                                   column = "ENTREZID",
                                                   keytype = "SYMBOL",
                                                   multiVals = "first")
                gene_list_entrez <- gene_list[!is.na(gene_ids)]
                names(gene_list_entrez) <- gene_ids[!is.na(gene_ids)]

                if (length(gene_list_entrez) >= 10) {
                  gsea_result <- ReactomePA::gsePathway(
                    geneList = gene_list_entrez,
                    organism = "human",
                    minGSSize = 10,
                    maxGSSize = 500,
                    pvalueCutoff = 1,
                    eps = 0
                  )
                }
              } else if (db %in% c("Hallmark", "Oncogenic")) {
                msig_category <- if (db == "Hallmark") "H" else "C6"
                msig_db <- msigdbr::msigdbr(species = "Homo sapiens", category = msig_category)
                msig_t2g <- msig_db %>% dplyr::select(gs_name, gene_symbol)

                gsea_result <- clusterProfiler::GSEA(
                  geneList = gene_list,
                  TERM2GENE = msig_t2g,
                  minGSSize = 10,
                  maxGSSize = 500,
                  pvalueCutoff = 1,
                  eps = 0
                )
              }

              if (!is.null(gsea_result) && nrow(as.data.frame(gsea_result)) > 0) {
                gsea_ran <- TRUE
                saveRDS(gsea_result, paste0(db_files_dir, "glm_gsea_pathway_", db, "_", analysis_name,
                                            "_", contrast_clean, "_minsample", minsample, "_full.rds"))
                write.csv(as.data.frame(gsea_result), paste0(db_files_dir, "glm_gsea_pathway_", db, "_", analysis_name,
                                                              "_", contrast_clean, "_minsample", minsample, "_full.csv"), row.names = FALSE)

                for (qpathway in Q_PATHWAY_VALUES) {
                  gsea_df <- as.data.frame(gsea_result)
                  filtered_df <- gsea_df %>% filter(qvalue < qpathway)
                  pathways_sig <- nrow(filtered_df)
                  current <- pathway_counts_by_db[[db]][[as.character(qpathway)]]
                  pathway_counts_by_db[[db]][[as.character(qpathway)]] <- current + pathways_sig

                  if (pathways_sig > 0) {
                    prefix <- paste0("glm_gsea_pathway_", db, "_", analysis_name,
                                     "_", contrast_clean, "_minsample", minsample, "_qpathway", qpathway)
                    write.csv(filtered_df, paste0(db_files_dir, prefix, ".csv"), row.names = FALSE)

                    tryCatch({
                      p_dot <- enrichplot::dotplot(gsea_result, showCategory = min(20, pathways_sig))
                      ggsave(paste0(db_dir, prefix, "_dot.png"), p_dot, width = 10, height = 8)
                    }, error = function(e) NULL)

                    tryCatch({
                      p_ridge <- enrichplot::ridgeplot(gsea_result, showCategory = min(15, pathways_sig))
                      ggsave(paste0(db_dir, prefix, "_ridge.png"), p_ridge, width = 10, height = 8)
                    }, error = function(e) NULL)

                    cat("    Saved:", prefix, "(", pathways_sig, "terms)\n")
                  }
                }
              }
            }, error = function(e) {
              cat("    Error running GSEA for", db, ":", e$message, "\n")
            })
          }
        }
      }
    }

    # Add run summary entries
    for (qpathway in Q_PATHWAY_VALUES) {
      pathway_counts <- sapply(databases, function(db) {
        count <- pathway_counts_by_db[[db]][[as.character(qpathway)]]
        if (is.null(count)) count <- 0
        paste0(db, ":", count)
      })
      pathways_str <- paste(pathway_counts, collapse = ", ")

      run_summary <<- rbind(run_summary, data.frame(
        analysis_name = analysis_name, analysis_type = "glm_gsea",
        group_col = group_col,
        any_pathway_significant = any(grepl(":[1-9]", pathway_counts)),
        pathways_significant = pathways_str, n_samples_per_group = n_samples_str,
        genes_tested = genes_tested, genes_significant = NA_integer_,
        enrichment_ran = gsea_ran, minsample = minsample, qgene = NA_real_,
        mingene = NA_integer_, qpathway = qpathway,
        groups = paste(groups, collapse = ","), glm_errors = glm_errors,
        posthoc_errors = posthoc_errors, posthoc_skipped = posthoc_skipped,
        stringsAsFactors = FALSE
      ))
    }
  }
}

#### PATHWAY ANALYSIS 1: KICS VS HOSTSEQ ####
write_output(quote(NULL), "Pathway Analysis - KICS vs HostSeq")

if (exists("te_kics_hostseq") && nrow(te_kics_hostseq) > 0) {
  cat("Dataset: te_kics_hostseq\n")
  cat("Samples:", length(unique(te_kics_hostseq$sample)), "\n")
  cohort_groups <- unique(te_kics_hostseq$cohort)
  cohort_groups <- cohort_groups[!is.na(cohort_groups)]
  cat("Cohort groups:", paste(cohort_groups, collapse = ", "), "\n")

  if (length(cohort_groups) >= 2) {
    run_simple_ora_analysis(
      te_data = te_kics_hostseq,
      group_col = "cohort",
      groups = cohort_groups,
      analysis_name = "kics_hostseq",
      base_dir = ora_dir
    )
    run_glm_ora_analysis(
      te_data = te_kics_hostseq,
      group_col = "cohort",
      groups = cohort_groups,
      analysis_name = "kics_hostseq",
      base_dir = glm_ora_dir,
      covariates = covar_no_age
    )
    run_glm_gsea_analysis(
      te_data = te_kics_hostseq,
      group_col = "cohort",
      groups = cohort_groups,
      analysis_name = "kics_hostseq",
      base_dir = glm_gsea_dir,
      covariates = covar_no_age
    )
  }
}

#### PATHWAY ANALYSIS 2: TP53 STATUS ####
write_output(quote(NULL), "Pathway Analysis - TP53 Status")

cat("Dataset: te_aff_split\n")
tp53_groups <- unique(te_aff_split$TP53_status)
tp53_groups <- tp53_groups[!is.na(tp53_groups)]
cat("TP53 groups:", paste(tp53_groups, collapse = ", "), "\n")

if (length(tp53_groups) >= 2) {
  run_simple_ora_analysis(
    te_data = te_aff_split,
    group_col = "TP53_status",
    groups = tp53_groups,
    analysis_name = "tp53_status",
    base_dir = ora_dir
  )
  run_glm_ora_analysis(
    te_data = te_aff_split,
    group_col = "TP53_status",
    groups = tp53_groups,
    analysis_name = "tp53_status",
    base_dir = glm_ora_dir
  )
  run_glm_gsea_analysis(
    te_data = te_aff_split,
    group_col = "TP53_status",
    groups = tp53_groups,
    analysis_name = "tp53_status",
    base_dir = glm_gsea_dir
  )
}

#### PATHWAY ANALYSIS 3: TP53 3-LEVEL ####
write_output(quote(NULL), "Pathway Analysis - TP53 3-Level")

te_aff_3level <- te_aff_split %>% filter(!is.na(TP53_3level))
if (nrow(te_aff_3level) > 0) {
  tp53_3level_groups <- unique(te_aff_3level$TP53_3level)
  tp53_3level_groups <- tp53_3level_groups[!is.na(tp53_3level_groups)]
  cat("TP53_3level groups:", paste(tp53_3level_groups, collapse = ", "), "\n")

  if (length(tp53_3level_groups) >= 2) {
    run_simple_ora_analysis(
      te_data = te_aff_3level,
      group_col = "TP53_3level",
      groups = tp53_3level_groups,
      analysis_name = "tp53_3level",
      base_dir = ora_dir
    )
    run_glm_ora_analysis(
      te_data = te_aff_3level,
      group_col = "TP53_3level",
      groups = tp53_3level_groups,
      analysis_name = "tp53_3level",
      base_dir = glm_ora_dir
    )
    run_glm_gsea_analysis(
      te_data = te_aff_3level,
      group_col = "TP53_3level",
      groups = tp53_3level_groups,
      analysis_name = "tp53_3level",
      base_dir = glm_gsea_dir
    )
  }
}

#### PATHWAY ANALYSIS 4: CANCER STATUS (LFS) ####
write_output(quote(NULL), "Pathway Analysis - Cancer Status LFS")

cat("Dataset: te_lfs_split\n")
cancer_groups <- unique(te_lfs_split$Cancer)
cancer_groups <- cancer_groups[!is.na(cancer_groups)]
cat("Cancer groups:", paste(cancer_groups, collapse = ", "), "\n")

if (length(cancer_groups) >= 2) {
  run_simple_ora_analysis(
    te_data = te_lfs_split,
    group_col = "Cancer",
    groups = cancer_groups,
    analysis_name = "cancer_lfs",
    base_dir = ora_dir
  )
  run_glm_ora_analysis(
    te_data = te_lfs_split,
    group_col = "Cancer",
    groups = cancer_groups,
    analysis_name = "cancer_lfs",
    base_dir = glm_ora_dir,
    covariates = covar_no_age
  )
  run_glm_gsea_analysis(
    te_data = te_lfs_split,
    group_col = "Cancer",
    groups = cancer_groups,
    analysis_name = "cancer_lfs",
    base_dir = glm_gsea_dir,
    covariates = covar_no_age
  )
}

#### PATHWAY ANALYSIS 5: LFS BY TP53 STATUS ####
write_output(quote(NULL), "Pathway Analysis - LFS by TP53 Status")

if (exists("te_lfs_split") && "TP53_status" %in% colnames(te_lfs_split)) {
  lfs_tp53_groups <- unique(te_lfs_split$TP53_status)
  lfs_tp53_groups <- lfs_tp53_groups[!is.na(lfs_tp53_groups)]

  if (length(lfs_tp53_groups) >= 2) {
    run_simple_ora_analysis(
      te_data = te_lfs_split,
      group_col = "TP53_status",
      groups = lfs_tp53_groups,
      analysis_name = "lfs_tp53_status",
      base_dir = ora_dir
    )
    run_glm_ora_analysis(
      te_data = te_lfs_split,
      group_col = "TP53_status",
      groups = lfs_tp53_groups,
      analysis_name = "lfs_tp53_status",
      base_dir = glm_ora_dir
    )
    run_glm_gsea_analysis(
      te_data = te_lfs_split,
      group_col = "TP53_status",
      groups = lfs_tp53_groups,
      analysis_name = "lfs_tp53_status",
      base_dir = glm_gsea_dir
    )
  }
}

#### PATHWAY ANALYSIS 6: TAYLOR BY SUBTYPE ####
write_output(quote(NULL), "Pathway Analysis - Taylor by Subtype")

if (exists("te_taylor_split") && nrow(te_taylor_split) > 0) {
  if ("tumor_type_subclass" %in% colnames(te_taylor_split)) {
    # Filter to subtypes with >= 3 samples
    subtype_counts <- table(te_taylor_split$tumor_type_subclass)
    valid_subtypes <- names(subtype_counts[subtype_counts >= 3])
    cat("Subtypes with >= 3 samples:", paste(valid_subtypes, collapse = ", "), "\n")

    if (length(valid_subtypes) >= 2) {
      te_taylor_subtype <- te_taylor_split %>%
        filter(tumor_type_subclass %in% valid_subtypes)

      run_simple_ora_analysis(
        te_data = te_taylor_subtype,
        group_col = "tumor_type_subclass",
        groups = valid_subtypes,
        analysis_name = "taylor_subtype",
        base_dir = ora_dir
      )
      run_glm_ora_analysis(
        te_data = te_taylor_subtype,
        group_col = "tumor_type_subclass",
        groups = valid_subtypes,
        analysis_name = "taylor_subtype",
        base_dir = glm_ora_dir
      )
      run_glm_gsea_analysis(
        te_data = te_taylor_subtype,
        group_col = "tumor_type_subclass",
        groups = valid_subtypes,
        analysis_name = "taylor_subtype",
        base_dir = glm_gsea_dir
      )
    } else {
      cat("Not enough subtypes with >= 3 samples. Skipping.\n")
    }
  }
}

#### PATHWAY ANALYSIS 7: TUMOR TYPE (KICS) ####
write_output(quote(NULL), "Pathway Analysis - KICS by Tumor Type")

if ("tumor_type" %in% colnames(te_kics_split)) {
  # Filter to tumor types with >= 3 samples
  tumor_type_counts <- table(te_kics_split$tumor_type)
  valid_tumor_types <- names(tumor_type_counts[tumor_type_counts >= 3])
  cat("Tumor types with >= 3 samples:", paste(valid_tumor_types, collapse = ", "), "\n")

  if (length(valid_tumor_types) >= 2) {
    te_kics_tumor <- te_kics_split %>%
      filter(tumor_type %in% valid_tumor_types)

    run_simple_ora_analysis(
      te_data = te_kics_tumor,
      group_col = "tumor_type",
      groups = valid_tumor_types,
      analysis_name = "kics_tumor_type",
      base_dir = ora_dir
    )
    run_glm_ora_analysis(
      te_data = te_kics_tumor,
      group_col = "tumor_type",
      groups = valid_tumor_types,
      analysis_name = "kics_tumor_type",
      base_dir = glm_ora_dir,
      covariates = covar_no_tumour_type
    )
    run_glm_gsea_analysis(
      te_data = te_kics_tumor,
      group_col = "tumor_type",
      groups = valid_tumor_types,
      analysis_name = "kics_tumor_type",
      base_dir = glm_gsea_dir,
      covariates = covar_no_tumour_type
    )
  } else {
    cat("Not enough tumor types with >= 3 samples. Skipping.\n")
  }
}

#### PATHWAY ANALYSIS 8: KICS BY ANCESTRY ####
write_output(quote(NULL), "Pathway Analysis - KICS by Ancestry")

if ("predicted_ancestry_thres" %in% colnames(te_kics_split)) {
  ancestries <- unique(te_kics_split$predicted_ancestry_thres)
  ancestries <- ancestries[!is.na(ancestries)]
  cat("Ancestries:", paste(ancestries, collapse = ", "), "\n")

  ancestry_counts <- table(te_kics_split$predicted_ancestry_thres)
  valid_ancestries <- names(ancestry_counts[ancestry_counts >= 10])

  if (length(valid_ancestries) >= 2) {
    te_kics_ancestry <- te_kics_split %>%
      filter(predicted_ancestry_thres %in% valid_ancestries)

    run_simple_ora_analysis(
      te_data = te_kics_ancestry,
      group_col = "predicted_ancestry_thres",
      groups = valid_ancestries,
      analysis_name = "kics_ancestry",
      base_dir = ora_dir
    )
    run_glm_ora_analysis(
      te_data = te_kics_ancestry,
      group_col = "predicted_ancestry_thres",
      groups = valid_ancestries,
      analysis_name = "kics_ancestry",
      base_dir = glm_ora_dir,
      covariates = covar_no_ancestry
    )
    run_glm_gsea_analysis(
      te_data = te_kics_ancestry,
      group_col = "predicted_ancestry_thres",
      groups = valid_ancestries,
      analysis_name = "kics_ancestry",
      base_dir = glm_gsea_dir,
      covariates = covar_no_ancestry
    )
  }
}

#### PATHWAY ANALYSIS 9: LFS BY COHORT ####
write_output(quote(NULL), "Pathway Analysis - LFS by Cohort")

if ("cohort" %in% colnames(te_lfs_split)) {
  cohorts <- unique(te_lfs_split$cohort)
  cohorts <- cohorts[!is.na(cohorts)]
  cat("Cohorts:", paste(cohorts, collapse = ", "), "\n")

  if (length(cohorts) >= 2) {
    run_simple_ora_analysis(
      te_data = te_lfs_split,
      group_col = "cohort",
      groups = cohorts,
      analysis_name = "lfs_cohort",
      base_dir = ora_dir
    )
    run_glm_ora_analysis(
      te_data = te_lfs_split,
      group_col = "cohort",
      groups = cohorts,
      analysis_name = "lfs_cohort",
      base_dir = glm_ora_dir
    )
    run_glm_gsea_analysis(
      te_data = te_lfs_split,
      group_col = "cohort",
      groups = cohorts,
      analysis_name = "lfs_cohort",
      base_dir = glm_gsea_dir
    )
  }
}

#### PATHWAY ANALYSIS 10: KICS BY SAMPLE TYPE ####
write_output(quote(NULL), "Pathway Analysis - KICS by Sample Type")

sample_type_file <- "/Users/briannelaverty/Documents/R_Malkin/clinical/kics_germline_sample_type.csv"
if (file.exists(sample_type_file)) {
  kics_sample_type <- prep_kics_sample_type(sample_type_file)
  cat("Loaded sample type data:", nrow(kics_sample_type), "samples\n")

  te_kics_sampletype <- merge(te_kics_split, kics_sample_type,
                               by = "sample", all.x = TRUE)
  cat("Merged sample type with data:", nrow(te_kics_sampletype), "rows\n")

  valid_sample_types <- c("Blood", "Fibroblasts", "Tissue (fresh)")
  te_kics_sampletype <- te_kics_sampletype %>%
    filter(sample_type %in% valid_sample_types)
  cat("After filtering to valid types:", nrow(te_kics_sampletype), "rows\n")

  sample_types <- unique(te_kics_sampletype$sample_type)
  sample_types <- sample_types[!is.na(sample_types)]
  cat("Sample types:", paste(sample_types, collapse = ", "), "\n")

  if (length(sample_types) >= 2) {
    run_simple_ora_analysis(
      te_data = te_kics_sampletype,
      group_col = "sample_type",
      groups = sample_types,
      analysis_name = "kics_sample_type",
      base_dir = ora_dir
    )
    run_glm_ora_analysis(
      te_data = te_kics_sampletype,
      group_col = "sample_type",
      groups = sample_types,
      analysis_name = "kics_sample_type",
      base_dir = glm_ora_dir
    )
    run_glm_gsea_analysis(
      te_data = te_kics_sampletype,
      group_col = "sample_type",
      groups = sample_types,
      analysis_name = "kics_sample_type",
      base_dir = glm_gsea_dir
    )
  }
}

#### SUMMARY ####
write_output(quote(NULL), "Pathway Analysis Summary")

cat("Analysis completed for the following datasets:\n")
cat("  1. kics_hostseq - KICS vs HostSeq\n")
cat("  2. tp53_status - TP53 status (affected samples)\n")
cat("  3. tp53_3level - TP53 3-level (Germline/Somatic/WT)\n")
cat("  4. cancer_lfs - Cancer status LFS\n")
cat("  5. lfs_tp53_status - LFS by TP53 Status\n")
cat("  6. taylor_subtype - Taylor by subtype\n")
cat("  7. kics_tumor_type - KICS by tumor type\n")
cat("  8. kics_ancestry - KICS by ancestry\n")
cat("  9. lfs_cohort - LFS by cohort\n")
cat("  10. kics_sample_type - KICS by sample type\n\n")

cat("Analysis types run for each:\n")
cat("  - ora/: Simple ORA (minsample filter only, no GLM)\n")
cat("  - glm_ora/: GLM-based ORA with covariate control\n")
cat("  - glm_gsea/: GLM-based GSEA using z-statistic ranking\n\n")

cat("Databases used:", paste(DATABASES_TO_RUN, collapse = ", "), "\n")
cat("Min samples:", paste(MIN_SAMPLES_VALUES, collapse = ", "), "\n")
cat("Q-gene thresholds (glm_ora):", paste(Q_GENE_VALUES, collapse = ", "), "\n")
cat("Min gene counts (ora, glm_ora):", paste(MIN_GENE_VALUES, collapse = ", "), "\n")
cat("Q-pathway thresholds:", paste(Q_PATHWAY_VALUES, collapse = ", "), "\n")

cat("\nOutput structure:\n")
cat("  ora/{analysis}/{database}/ - Simple ORA results and plots\n")
cat("  glm_ora/{analysis}/{database}/ - GLM-ORA results and plots\n")
cat("  glm_gsea/{analysis}/{database}/ - GLM-GSEA results and plots\n")
cat("Output location:", pathway_dir, "\n")

# Sort run summary by analysis_name first to group all same analyses together
run_summary <- run_summary %>%
  arrange(analysis_name, analysis_type, minsample, qgene, mingene, qpathway)

run_summary_file <- paste0(pathway_dir, "RUN_SUMMARY.csv")
write.csv(run_summary, run_summary_file, row.names = FALSE)
cat("\nSaved run summary:", run_summary_file, "\n")

cat("\n✓ Script completed successfully\n")

# Close module-specific sink
close_module_sink()
