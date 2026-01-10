#!/usr/bin/env Rscript

# Tumour TE Visualization - Pathway Analysis (Multi-Database)
# Three analysis types: Simple ORA, GLM-ORA, and GLM-GSEA
# Databases: GO_BP, Reactome, Hallmark, Oncogenic (no KEGG)

# Source common setup and load data
source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/00_viz_common_setup.R")
REQUIRED_DATA <- c("count_matrix", "split", "clinical", "genes")
source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/00_viz_load_data_tumour.R")

# Load emmeans for post-hoc contrasts
library(emmeans)

# Initialize module-specific text output
init_module_sink(paste0(plot_dir, "pathway/"), "PATHWAY")

# Initialize run summary with analysis_type column for ora/glm_ora/glm_gsea
# Sorted by analysis_name first to group all same analyses together
run_summary <- data.frame(
  analysis_name = character(),
  analysis_type = character(),        # "ora", "glm_ora", or "glm_gsea"
  group_col = character(),
  any_pathway_significant = logical(), # TRUE if any pathway count > 0
  pathways_significant = character(),  # "GO_BP:5, Reactome:3, Hallmark:0, Oncogenic:0"
  n_samples_per_group = character(),
  genes_tested = integer(),           # ora: genes with minsample; glm_*: genes tested by GLM
  genes_significant = integer(),      # ora: genes passing minsample; glm_*: genes passing qgene
  enrichment_ran = logical(),         # Did ORA/GSEA run?
  minsample = integer(),
  qgene = numeric(),                  # ora/glm_gsea: NA; glm_ora: qgene threshold
  mingene = integer(),                # ora/glm_ora: mingene; glm_gsea: NA
  qpathway = numeric(),
  groups = character(),
  glm_errors = integer(),             # ora: NA; glm_*: GLM errors
  posthoc_errors = integer(),
  posthoc_skipped = character(),
  stringsAsFactors = FALSE
)

cat("Running 02_te_viz_tumour_06_pathway.R...\n")
cat("Multi-database pathway analysis: Simple ORA, GLM-ORA, and GLM-GSEA\n")
cat("Databases: GO_BP, Reactome, Hallmark, Oncogenic\n\n")

#### CONFIGURATION ####

# Databases to run (no KEGG per plan)
DATABASES_TO_RUN <- c("GO_BP", "Reactome", "Hallmark", "Oncogenic")

# Parameter settings (per plan)
MIN_SAMPLES_VALUES <- c(3, 5)
Q_GENE_VALUES <- c(0.05, 0.1, 0.25)  # For GLM gene filtering
MIN_GENE_VALUES <- c(3, 5)           # Min genes hitting pathway
Q_PATHWAY_VALUES <- c(0.05, 0.1, 0.25)  # For pathway filtering

# Define TE types for legacy analyses
types <- c(NA, "LINE1", "ALU", "SVA")
x_tp53 <- expression("Somatic " * italic("TP53") * " status")

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
pathway_base_dir <- paste0(plot_dir, "pathway/")
ora_dir <- paste0(pathway_base_dir, "ora/")
glm_ora_dir <- paste0(pathway_base_dir, "glm_ora/")
glm_gsea_dir <- paste0(pathway_base_dir, "glm_gsea/")

dir.create(pathway_base_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(ora_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(glm_ora_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(glm_gsea_dir, showWarnings = FALSE, recursive = TRUE)

cat("Base output directory:", pathway_base_dir, "\n")
cat("  - ora/: Simple ORA (minsample filter only)\n")
cat("  - glm_ora/: GLM-based ORA\n")
cat("  - glm_gsea/: GLM-based GSEA\n\n")

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

#### HELPER FUNCTIONS ####

# Helper: Create ORA visualizations for a filtered result
# ora_result: the original enrichResult/compareClusterResult object
# filtered_df: the post-hoc filtered dataframe (by mingene, qpathway)
create_ora_visualizations <- function(ora_result, filtered_df, te_data, db, prefix, db_dir,
                                       gene_col = "Gene_name", sample_col = "sample") {
  n_terms <- nrow(filtered_df)
  if (n_terms == 0) return(NULL)

  # Filter the enrichResult object to only include pathways in filtered_df
  # This ensures plots show only the post-hoc filtered pathways
  filtered_result <- ora_result
  if ("ID" %in% colnames(filtered_df)) {
    filtered_ids <- unique(filtered_df$ID)
    if (inherits(ora_result, "enrichResult")) {
      # enrichResult uses @result slot
      filtered_result@result <- filtered_result@result %>%
        filter(ID %in% filtered_ids)
    } else if (inherits(ora_result, "compareClusterResult")) {
      # compareClusterResult uses @compareClusterResult slot
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

    # 5-7. Sample-level visualizations (balloon, heatmap, gene_pathway)
    # These use filtered_df directly (already filtered)
    if (!is.null(te_data) && n_terms > 0) {
      tryCatch({
        # Balloon plot
        p_balloon <- plot_pathway_balloon(filtered_df, te_data, max_pathways = 20,
                                           gene_col = gene_col)
        if (!is.null(p_balloon)) {
          ggsave(paste0(db_dir, prefix, "_balloon.png"), p_balloon, width = 14, height = 8)
        }
      }, error = function(e) cat("    Balloon plot error:", e$message, "\n"))

      tryCatch({
        # Sample-pathway heatmap
        ht <- plot_pathway_sample_heatmap(filtered_df, te_data, max_pathways = 30,
                                           gene_col = gene_col)
        if (!is.null(ht)) {
          png(paste0(db_dir, prefix, "_heatmap.png"), width = 12, height = 10, units = "in", res = 150)
          ComplexHeatmap::draw(ht)
          dev.off()
        }
      }, error = function(e) cat("    Heatmap error:", e$message, "\n"))

      tryCatch({
        # Gene-pathway heatmap
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
      # Convert to Entrez IDs for Reactome
      gene_clusters_entrez <- lapply(gene_clusters, function(genes) {
        ids <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, genes, "ENTREZID", "SYMBOL")
        ids[!is.na(ids)]
      })
      gene_clusters_entrez <- gene_clusters_entrez[sapply(gene_clusters_entrez, length) > 0]

      if (length(gene_clusters_entrez) >= 2) {
        clusterProfiler::compareCluster(
          geneCluster = gene_clusters_entrez,
          fun = "enrichPathway",
          organism = "human",
          pvalueCutoff = pvalueCutoff,
          qvalueCutoff = 1
        )
      } else {
        NULL
      }
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
    } else {
      NULL
    }
  }, error = function(e) {
    cat("    Error running compareCluster for", db, ":", e$message, "\n")
    NULL
  })
}

# Helper: Assign genes to groups based on which group has highest incidence
assign_genes_to_groups_by_incidence <- function(te_data, genes, group_col, groups,
                                                  gene_col = "Gene_name", sample_col = "sample") {
  gene_clusters <- list()

  for (g in groups) {
    group_data <- te_data %>% filter(.data[[group_col]] == g)
    group_samples <- unique(group_data[[sample_col]])
    n_group <- length(group_samples)

    if (n_group == 0) next

    # For each gene, count samples in this group
    gene_incidence <- group_data %>%
      filter(.data[[gene_col]] %in% genes) %>%
      group_by(.data[[gene_col]]) %>%
      summarise(n_samples = n_distinct(.data[[sample_col]]), .groups = "drop") %>%
      mutate(incidence = n_samples / n_group)

    gene_clusters[[g]] <- gene_incidence
  }

  # Assign each gene to group with highest incidence
  all_genes <- unique(unlist(lapply(gene_clusters, function(x) x[[gene_col]])))
  gene_assignments <- list()

  for (gene in all_genes) {
    max_incidence <- 0
    assigned_group <- groups[1]  # default

    for (g in groups) {
      if (!is.null(gene_clusters[[g]])) {
        gene_row <- gene_clusters[[g]] %>% filter(.data[[gene_col]] == gene)
        if (nrow(gene_row) > 0 && gene_row$incidence[1] > max_incidence) {
          max_incidence <- gene_row$incidence[1]
          assigned_group <- g
        }
      }
    }

    if (!assigned_group %in% names(gene_assignments)) {
      gene_assignments[[assigned_group]] <- c()
    }
    gene_assignments[[assigned_group]] <- c(gene_assignments[[assigned_group]], gene)
  }

  return(gene_assignments)
}

#### SIMPLE ORA FUNCTION (no GLM, just minsample filter) ####
run_simple_ora_analysis <- function(te_data,
                                     group_col,
                                     groups,
                                     analysis_name,
                                     base_dir,  # ora/
                                     databases = DATABASES_TO_RUN,
                                     background_genes = BACKGROUND_GENES,
                                     gene_col = "Gene_name",
                                     sample_col = "sample") {

  cat("\n========================================\n")
  cat("  Simple ORA Analysis:", analysis_name, "\n")
  cat("========================================\n")
  cat("Group column:", group_col, "\n")
  cat("Groups:", paste(groups, collapse = ", "), "\n")

  # Calculate sample counts per group for run summary
  actual_sample_col <- if (sample_col %in% colnames(te_data)) sample_col else "sample"
  sample_counts <- te_data %>%
    filter(!is.na(.data[[group_col]])) %>%
    group_by(.data[[group_col]]) %>%
    summarise(n = n_distinct(.data[[actual_sample_col]]), .groups = "drop")
  sample_counts_vec <- setNames(sample_counts$n, as.character(sample_counts[[group_col]]))
  n_samples_str <- paste(names(sample_counts_vec), sample_counts_vec, sep = ":", collapse = ", ")

  analysis_base_dir <- paste0(base_dir, analysis_name, "/")

  for (minsample in MIN_SAMPLES_VALUES) {
    cat("\n--- minsample =", minsample, "---\n")

    # Get genes with ≥ minsample samples having TE (overall, not per group)
    gene_sample_counts <- te_data %>%
      filter(!is.na(.data[[gene_col]]) & .data[[gene_col]] != "") %>%
      group_by(.data[[gene_col]]) %>%
      summarise(n_samples = n_distinct(.data[[actual_sample_col]]), .groups = "drop") %>%
      filter(n_samples >= minsample)

    genes_passing <- unique(gene_sample_counts[[gene_col]])
    cat("Genes with >=", minsample, "samples:", length(genes_passing), "\n")

    if (length(genes_passing) < 5) {
      cat("Too few genes. Skipping.\n")
      # Add run summary entries for this minsample
      for (mingene in MIN_GENE_VALUES) {
        for (q_pathway in Q_PATHWAY_VALUES) {
          pathways_str <- paste(sapply(databases, function(db) paste0(db, ":0")), collapse = ", ")
          run_summary <<- rbind(run_summary, data.frame(
            analysis_name = analysis_name,
            analysis_type = "ora",
            group_col = group_col,
            any_pathway_significant = FALSE,
            pathways_significant = pathways_str,
            n_samples_per_group = n_samples_str,
            genes_tested = length(genes_passing),
            genes_significant = length(genes_passing),
            enrichment_ran = FALSE,
            minsample = minsample,
            qgene = NA_real_,
            mingene = mingene,
            qpathway = q_pathway,
            groups = paste(groups, collapse = ","),
            glm_errors = NA_integer_,
            posthoc_errors = NA_integer_,
            posthoc_skipped = NA_character_,
            stringsAsFactors = FALSE
          ))
        }
      }
      next
    }

    # Assign genes to groups based on which group has highest incidence
    gene_clusters <- assign_genes_to_groups_by_incidence(
      te_data, genes_passing, group_col, groups, gene_col, actual_sample_col
    )

    # Run compareCluster ORA for each database
    for (db in databases) {
      db_dir <- paste0(analysis_base_dir, db, "/")
      db_files_dir <- paste0(db_dir, "files/")
      dir.create(db_files_dir, showWarnings = FALSE, recursive = TRUE)

      cat("\n  Database:", db, "\n")
      ora_result <- run_compareCluster_ora_by_db(gene_clusters, db, background_genes, pvalueCutoff = 1)

      if (is.null(ora_result) || nrow(as.data.frame(ora_result)) == 0) {
        cat("    No enriched terms found\n")
        next
      }

      # Save full ORA result to {db}/files/
      saveRDS(ora_result, paste0(db_files_dir, "ora_", db, "_", analysis_name,
                                  "_minsample", minsample, "_full.rds"))
      write.csv(as.data.frame(ora_result), paste0(db_files_dir, "ora_", db, "_", analysis_name,
                                  "_minsample", minsample, "_full.csv"), row.names = FALSE)
      cat("    Saved full result:", nrow(as.data.frame(ora_result)), "terms\n")

      # Post-hoc filter by mingene and qpathway
      for (mingene in MIN_GENE_VALUES) {
        for (q_pathway in Q_PATHWAY_VALUES) {
          ora_df <- as.data.frame(ora_result)
          filtered_df <- ora_df %>% filter(Count >= mingene, qvalue < q_pathway)

          pathways_sig <- nrow(filtered_df)

          if (pathways_sig > 0) {
            prefix <- paste0("ora_", db, "_", analysis_name, "_minsample", minsample,
                             "_mingene", mingene, "_qpathway", q_pathway)

            # Save filtered CSV to files/
            write.csv(filtered_df, paste0(db_files_dir, prefix, ".csv"), row.names = FALSE)

            # Create visualizations with filtered data
            create_ora_visualizations(ora_result, filtered_df, te_data, db, prefix, db_dir,
                                       gene_col, actual_sample_col)

            cat("    Saved:", prefix, "(", pathways_sig, "terms)\n")
          }
        }
      }
    }

    # Add run summary entries for this minsample (one row per mingene/qpathway combo)
    pathway_counts_by_db <- list()
    for (db in databases) {
      db_files_dir <- paste0(analysis_base_dir, db, "/files/")
      for (mingene in MIN_GENE_VALUES) {
        for (q_pathway in Q_PATHWAY_VALUES) {
          key <- paste(mingene, q_pathway, sep = "_")
          csv_file <- paste0(db_files_dir, "ora_", db, "_", analysis_name,
                             "_minsample", minsample, "_mingene", mingene, "_qpathway", q_pathway, ".csv")
          if (file.exists(csv_file)) {
            df <- read.csv(csv_file)
            if (!db %in% names(pathway_counts_by_db)) pathway_counts_by_db[[db]] <- list()
            pathway_counts_by_db[[db]][[key]] <- nrow(df)
          } else {
            if (!db %in% names(pathway_counts_by_db)) pathway_counts_by_db[[db]] <- list()
            pathway_counts_by_db[[db]][[key]] <- 0
          }
        }
      }
    }

    for (mingene in MIN_GENE_VALUES) {
      for (q_pathway in Q_PATHWAY_VALUES) {
        key <- paste(mingene, q_pathway, sep = "_")
        pathway_counts <- sapply(databases, function(db) {
          count <- pathway_counts_by_db[[db]][[key]]
          if (is.null(count)) count <- 0
          paste0(db, ":", count)
        })
        pathways_str <- paste(pathway_counts, collapse = ", ")

        # Check if any pathway has count > 0
        any_sig <- any(sapply(databases, function(db) {
          count <- pathway_counts_by_db[[db]][[key]]
          !is.null(count) && count > 0
        }))

        run_summary <<- rbind(run_summary, data.frame(
          analysis_name = analysis_name,
          analysis_type = "ora",
          group_col = group_col,
          any_pathway_significant = any_sig,
          pathways_significant = pathways_str,
          n_samples_per_group = n_samples_str,
          genes_tested = length(genes_passing),
          genes_significant = length(genes_passing),
          enrichment_ran = TRUE,
          minsample = minsample,
          qgene = NA_real_,
          mingene = mingene,
          qpathway = q_pathway,
          groups = paste(groups, collapse = ","),
          glm_errors = NA_integer_,
          posthoc_errors = NA_integer_,
          posthoc_skipped = NA_character_,
          stringsAsFactors = FALSE
        ))
      }
    }
  }
}

#### GLM-BASED ORA FUNCTION ####
# Run GLM once per minsample, filter by qgene afterward
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

  # Calculate sample counts per group for run summary
  actual_sample_col <- if (sample_col %in% colnames(te_data)) sample_col else "sample"
  sample_counts <- te_data %>%
    filter(!is.na(.data[[group_col]])) %>%
    group_by(.data[[group_col]]) %>%
    summarise(n = n_distinct(.data[[actual_sample_col]]), .groups = "drop")

  # Convert to character to avoid factor level issues and build n_samples_str
  sample_counts_vec <- setNames(sample_counts$n, as.character(sample_counts[[group_col]]))
  n_samples_str <- paste(names(sample_counts_vec), sample_counts_vec, sep = ":", collapse = ", ")

  # Create analysis-specific directories: analysis_name/files/ and analysis_name/db/
  analysis_base_dir <- paste0(base_dir, analysis_name, "/")
  files_dir <- paste0(analysis_base_dir, "files/")
  dir.create(files_dir, showWarnings = FALSE, recursive = TRUE)

  for (db in databases) {
    db_plot_dir <- paste0(analysis_base_dir, db, "/")
    dir.create(db_plot_dir, showWarnings = FALSE, recursive = TRUE)
  }

  # Run GLM ONCE per minsample with p_threshold=1.0 to get ALL genes
  for (minsample in MIN_SAMPLES_VALUES) {
    cat("\n--- minsample =", minsample, "---\n")

    # Run GLM ONCE per minsample to get ALL results (no p-value filtering)
    glm_result <- filter_genes_by_differential_incidence_glm(
      te_data = te_data,
      group_col = group_col,
      groups = groups,
      min_samples = minsample,
      p_threshold = 1.0,  # Get ALL genes, filter by qgene later
      gene_col = gene_col,
      covariates = covariates
    )

    # Extract full results
    full_results <- glm_result$full_results
    pairwise_results <- glm_result$pairwise

    # Extract metrics from GLM result
    genes_tested <- glm_result$genes_attempted
    glm_errors <- glm_result$glm_errors
    posthoc_errors <- glm_result$posthoc_errors
    posthoc_skipped <- glm_result$posthoc_skipped

    # Add comparison and enriched_in columns for 2-group comparisons
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

    # Save GLM results ONCE per minsample (all genes, sorted by p-value)
    glm_results_file <- paste0(files_dir, "glm_", analysis_name,
                               "_minsample", minsample, "_all.csv")
    if (!is.null(full_results) && nrow(full_results) > 0) {
      write.csv(full_results, glm_results_file, row.names = FALSE)
      cat("  Saved GLM results:", basename(glm_results_file), "(", nrow(full_results), "genes)\n")
    } else {
      # Save empty file with headers
      if (n_groups >= 3) {
        empty_df <- data.frame(
          gene = character(0), te_type = character(0), te_coordinates = character(0),
          te_location_type = character(0), te_location = character(0), gene_features = character(0),
          samples_with_te = integer(0), samples_with_te_names = character(0),
          samples_without_te = integer(0), effect_direction = character(0),
          lrt_pval = numeric(0), lrt_padj = numeric(0)
        )
      } else {
        empty_df <- data.frame(
          gene = character(0), te_type = character(0), te_coordinates = character(0),
          te_location_type = character(0), te_location = character(0), gene_features = character(0),
          samples_with_te = integer(0), samples_with_te_names = character(0),
          samples_without_te = integer(0), effect_direction = character(0),
          estimate = numeric(0), std_error = numeric(0), z_value = numeric(0),
          p_value = numeric(0), p_adj = numeric(0)
        )
      }
      write.csv(empty_df, glm_results_file, row.names = FALSE)
      cat("  Saved GLM results:", basename(glm_results_file), "(0 genes - empty)\n")
    }

    # Now loop through qgene thresholds and filter the SAME results
    for (qgene in Q_GENE_VALUES) {
      cat("\n  --- qgene =", qgene, "---\n")

      # Filter by qgene from the full results (using p_adj column)
      if (!is.null(full_results) && nrow(full_results) > 0) {
        p_col <- if ("lrt_padj" %in% colnames(full_results)) "lrt_padj" else "p_adj"
        sig_genes_df <- full_results %>% filter(.data[[p_col]] < qgene)
        sig_genes <- sig_genes_df$gene
      } else {
        sig_genes <- character(0)
      }

      genes_significant <- length(sig_genes)
      cat("  Significant genes at q <", qgene, ":", genes_significant, "\n")

      # Storage for pathway counts by database (for Problem 6)
      pathway_counts_by_db <- list()

      if (length(sig_genes) < 5) {
        cat("  Too few significant genes (", length(sig_genes), "). Skipping ORA.\n")

        # Initialize empty counts for all databases
        for (db in databases) {
          pathway_counts_by_db[[db]] <- list()
          for (mingene in MIN_GENE_VALUES) {
            for (qpathway in Q_PATHWAY_VALUES) {
              key <- paste(mingene, qpathway, sep = "_")
              pathway_counts_by_db[[db]][[key]] <- 0
            }
          }
        }

        # Problem 6: Add combined run_summary entries
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
              analysis_name = analysis_name,
              analysis_type = "glm_ora",
              group_col = group_col,
              any_pathway_significant = FALSE,
              pathways_significant = pathways_str,
              n_samples_per_group = n_samples_str,
              genes_tested = genes_tested,
              genes_significant = genes_significant,
              enrichment_ran = FALSE,
              minsample = minsample,
              qgene = qgene,
              mingene = mingene,
              qpathway = qpathway,
              groups = paste(groups, collapse = ","),
              glm_errors = glm_errors,
              posthoc_errors = posthoc_errors,
              posthoc_skipped = posthoc_skipped,
              stringsAsFactors = FALSE
            ))
          }
        }
        next
      }

      cat("  Running ORA on", length(sig_genes), "significant genes\n")

      # For 3+ groups, run compareCluster with post-hoc based gene assignment
      if (n_groups >= 3) {
        tryCatch({
          cat("\n    Running group-stratified compareCluster analysis (3+ groups)...\n")

          # Save pairwise post-hoc results
          if (!is.null(pairwise_results) && nrow(pairwise_results) > 0) {
            pairwise_file <- paste0(files_dir, "pairwise_posthoc_", analysis_name,
                                    "_minsample", minsample, "_qgene", qgene, ".csv")
            write.csv(pairwise_results, pairwise_file, row.names = FALSE)
            cat("    Saved pairwise post-hoc results:", basename(pairwise_file), "\n")

            # Assign genes to groups based on significant post-hoc contrasts
            gene_assignment <- assign_genes_by_posthoc(pairwise_results, p_threshold = 0.05)

            if (!is.null(gene_assignment) && nrow(gene_assignment) > 0) {
              # Save gene assignment
              gene_assign_file <- paste0(files_dir, "gene_group_assignment_", analysis_name,
                                         "_minsample", minsample, "_qgene", qgene, ".csv")
              gene_assign_expanded <- gene_assignment %>%
                tidyr::unnest(assigned_groups) %>%
                rename(assigned_group = assigned_groups)
              write.csv(gene_assign_expanded, gene_assign_file, row.names = FALSE)
              cat("    Saved gene-group assignment:", basename(gene_assign_file), "\n")

              # Create gene clusters for compareCluster
              gene_clusters <- gene_assignment_to_clusters(gene_assignment)

              # Only proceed if we have genes in at least 2 groups
              if (length(gene_clusters) >= 2 && all(sapply(gene_clusters, length) >= 1)) {
                # Run compareCluster for each database and collect counts
                for (db in databases) {
                  cat("\n    compareCluster for", db, "...\n")

                  db_plot_dir <- paste0(analysis_base_dir, db, "/")
                  pathway_counts_by_db[[db]] <- list()

                  tryCatch({
                    if (db == "GO_BP") {
                      compare_result <- clusterProfiler::compareCluster(
                        geneCluster = gene_clusters,
                        fun = "enrichGO",
                        OrgDb = org.Hs.eg.db::org.Hs.eg.db,
                        keyType = "SYMBOL",
                        ont = "BP",
                        universe = background_genes,
                        pvalueCutoff = 0.25,
                        qvalueCutoff = 0.25
                      )
                    } else if (db == "Reactome") {
                      gene_clusters_entrez <- lapply(gene_clusters, function(genes) {
                        ids <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
                                                     keys = genes,
                                                     column = "ENTREZID",
                                                     keytype = "SYMBOL",
                                                     multiVals = "first")
                        ids[!is.na(ids)]
                      })
                      compare_result <- clusterProfiler::compareCluster(
                        geneCluster = gene_clusters_entrez,
                        fun = ReactomePA::enrichPathway,
                        pvalueCutoff = 0.25,
                        qvalueCutoff = 0.25
                      )
                    } else {
                      compare_result <- NULL
                    }

                    if (!is.null(compare_result) && nrow(as.data.frame(compare_result)) > 0) {
                      compare_rds <- paste0(files_dir, "compareCluster_", db, "_", analysis_name,
                                            "_minsample", minsample, "_qgene", qgene, ".rds")
                      saveRDS(compare_result, compare_rds)
                      compare_csv <- paste0(files_dir, "compareCluster_", db, "_", analysis_name,
                                            "_minsample", minsample, "_qgene", qgene, ".csv")
                      compare_df <- as.data.frame(compare_result)
                      write.csv(compare_df, compare_csv, row.names = FALSE)
                      p_compare <- enrichplot::dotplot(compare_result, showCategory = 10)
                      ggsave(paste0(db_plot_dir, "pathway_", db, "_", analysis_name,
                                    "_minsample", minsample, "_qgene", qgene, "_compare_groups.png"),
                             p_compare, width = 14, height = 10)
                      cat("      Saved compareCluster results and plot for", db, "\n")

                      # Collect counts for each mingene/qpathway combo
                      for (mingene in MIN_GENE_VALUES) {
                        for (qpathway in Q_PATHWAY_VALUES) {
                          filtered_df <- compare_df %>% filter(Count >= mingene, qvalue < qpathway)
                          key <- paste(mingene, qpathway, sep = "_")
                          pathway_counts_by_db[[db]][[key]] <- nrow(filtered_df)
                        }
                      }
                    } else {
                      cat("      No enriched terms found for", db, "\n")
                      for (mingene in MIN_GENE_VALUES) {
                        for (qpathway in Q_PATHWAY_VALUES) {
                          key <- paste(mingene, qpathway, sep = "_")
                          pathway_counts_by_db[[db]][[key]] <- 0
                        }
                      }
                    }
                  }, error = function(e) {
                    cat("      Warning: compareCluster failed for", db, ":", e$message, "\n")
                    for (mingene in MIN_GENE_VALUES) {
                      for (qpathway in Q_PATHWAY_VALUES) {
                        key <- paste(mingene, qpathway, sep = "_")
                        pathway_counts_by_db[[db]][[key]] <- 0
                      }
                    }
                  })
                }
              } else {
                cat("    Not enough gene clusters for compareCluster\n")
                for (db in databases) {
                  pathway_counts_by_db[[db]] <- list()
                  for (mingene in MIN_GENE_VALUES) {
                    for (qpathway in Q_PATHWAY_VALUES) {
                      key <- paste(mingene, qpathway, sep = "_")
                      pathway_counts_by_db[[db]][[key]] <- 0
                    }
                  }
                }
              }
            } else {
              cat("    No genes with significant post-hoc contrasts\n")
              for (db in databases) {
                pathway_counts_by_db[[db]] <- list()
                for (mingene in MIN_GENE_VALUES) {
                  for (qpathway in Q_PATHWAY_VALUES) {
                    key <- paste(mingene, qpathway, sep = "_")
                    pathway_counts_by_db[[db]][[key]] <- 0
                  }
                }
              }
            }
          } else {
            cat("    No pairwise post-hoc results available\n")
            for (db in databases) {
              pathway_counts_by_db[[db]] <- list()
              for (mingene in MIN_GENE_VALUES) {
                for (qpathway in Q_PATHWAY_VALUES) {
                  key <- paste(mingene, qpathway, sep = "_")
                  pathway_counts_by_db[[db]][[key]] <- 0
                }
              }
            }
          }
        }, error = function(e) {
          cat("    Warning: Group-stratified analysis failed:", e$message, "\n")
          for (db in databases) {
            pathway_counts_by_db[[db]] <- list()
            for (mingene in MIN_GENE_VALUES) {
              for (qpathway in Q_PATHWAY_VALUES) {
                key <- paste(mingene, qpathway, sep = "_")
                pathway_counts_by_db[[db]][[key]] <- 0
              }
            }
          }
        })

        # Add combined run_summary entries for 3+ group analysis
        for (mingene in MIN_GENE_VALUES) {
          for (qpathway in Q_PATHWAY_VALUES) {
            key <- paste(mingene, qpathway, sep = "_")
            pathway_counts <- sapply(databases, function(db) {
              count <- pathway_counts_by_db[[db]][[key]]
              if (is.null(count)) count <- 0
              paste0(db, ":", count)
            })
            pathways_str <- paste(pathway_counts, collapse = ", ")

            # Check if any pathway has count > 0
            any_sig <- any(sapply(databases, function(db) {
              count <- pathway_counts_by_db[[db]][[key]]
              !is.null(count) && count > 0
            }))

            run_summary <<- rbind(run_summary, data.frame(
              analysis_name = analysis_name,
              analysis_type = "glm_ora",
              group_col = group_col,
              any_pathway_significant = any_sig,
              pathways_significant = pathways_str,
              n_samples_per_group = n_samples_str,
              genes_tested = genes_tested,
              genes_significant = genes_significant,
              enrichment_ran = TRUE,
              minsample = minsample,
              qgene = qgene,
              mingene = mingene,
              qpathway = qpathway,
              groups = paste(groups, collapse = ","),
              glm_errors = glm_errors,
              posthoc_errors = posthoc_errors,
              posthoc_skipped = posthoc_skipped,
              stringsAsFactors = FALSE
            ))
          }
        }
        next
      }

      # Run ORA for each database (2-group analyses only)
      for (db in databases) {
        cat("\n    Database:", db, "\n")

        db_plot_dir <- paste0(analysis_base_dir, db, "/")
        pathway_counts_by_db[[db]] <- list()

        tryCatch({
          ora_result <- perform_ora_single_database(
            genes = sig_genes,
            database = db,
            pvalueCutoff = 1,
            qvalueCutoff = 1,
            universe = background_genes
          )

          if (is.null(ora_result) || nrow(as.data.frame(ora_result)) == 0) {
            cat("      No enriched terms found\n")
            for (mingene in MIN_GENE_VALUES) {
              for (qpathway in Q_PATHWAY_VALUES) {
                key <- paste(mingene, qpathway, sep = "_")
                pathway_counts_by_db[[db]][[key]] <- 0
              }
            }
            next
          }

          # Save full ORA result object
          ora_rds_file <- paste0(files_dir, "pathway_ora_result_", db, "_", analysis_name,
                                  "_minsample", minsample, "_qgene", qgene, ".rds")
          saveRDS(ora_result, ora_rds_file)
          cat("      Saved ORA result object:", basename(ora_rds_file), "\n")

          # Post-hoc filter by mingene and qpathway
          ora_df <- as.data.frame(ora_result)

          for (mingene in MIN_GENE_VALUES) {
            for (qpathway in Q_PATHWAY_VALUES) {
              filtered_df <- ora_df %>% filter(Count >= mingene, qvalue < qpathway)
              key <- paste(mingene, qpathway, sep = "_")
              pathway_counts_by_db[[db]][[key]] <- nrow(filtered_df)

              if (nrow(filtered_df) == 0) next

              prefix <- paste0("pathway_", db, "_", analysis_name,
                               "_minsample", minsample, "_qgene", qgene,
                               "_mingene", mingene, "_qpathway", qpathway)

              csv_file <- paste0(db_plot_dir, prefix, ".csv")
              write.csv(filtered_df, csv_file, row.names = FALSE)

              tryCatch({
                p_dot <- dotplot(ora_result, showCategory = min(20, nrow(filtered_df)))
                ggsave(paste0(db_plot_dir, prefix, "_dot.png"), p_dot, width = 12, height = 8)
                p_bar <- barplot(ora_result, showCategory = min(20, nrow(filtered_df)))
                ggsave(paste0(db_plot_dir, prefix, "_bar.png"), p_bar, width = 10, height = 8)
                if (nrow(filtered_df) > 0 && nrow(filtered_df) <= 10) {
                  tryCatch({
                    p_cnet <- cnetplot(ora_result, showCategory = min(5, nrow(filtered_df)))
                    ggsave(paste0(db_plot_dir, prefix, "_cnet.png"), p_cnet, width = 12, height = 10)
                  }, error = function(e) NULL)
                }
              }, error = function(e) {
                cat("      Warning: Could not create standard plots:", e$message, "\n")
              })

              tryCatch({
                p_balloon <- plot_pathway_balloon(ora_df = filtered_df, te_split_data = te_data,
                                                  max_pathways = 15, gene_col = gene_col)
                if (!is.null(p_balloon)) {
                  ggsave(paste0(db_plot_dir, prefix, "_balloon.png"), p_balloon, width = 14, height = 8)
                }
                ht <- plot_pathway_sample_heatmap(ora_df = filtered_df, te_split_data = te_data,
                                                   max_pathways = 20, gene_col = gene_col)
                if (!is.null(ht)) {
                  png(paste0(db_plot_dir, prefix, "_heatmap.png"), width = 14, height = 8, units = "in", res = 300)
                  draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
                  dev.off()
                }
                p_gene <- plot_pathway_gene_heatmap(ora_df = filtered_df, te_split_data = te_data,
                                                     max_pathways = 10, max_genes = 30, gene_col = gene_col)
                if (!is.null(p_gene)) {
                  ggsave(paste0(db_plot_dir, prefix, "_gene_pathway.png"), p_gene, width = 12, height = 10)
                }
              }, error = function(e) {
                cat("      Warning: Could not create sample-level plots:", e$message, "\n")
              })

              cat("      Saved:", prefix, "(", nrow(filtered_df), "terms)\n")
            }
          }
        }, error = function(e) {
          cat("      Error running ORA for", db, ":", e$message, "\n")
          for (mingene in MIN_GENE_VALUES) {
            for (qpathway in Q_PATHWAY_VALUES) {
              key <- paste(mingene, qpathway, sep = "_")
              pathway_counts_by_db[[db]][[key]] <- 0
            }
          }
        })
      }

      # Add combined run_summary entries for 2-group analysis
      for (mingene in MIN_GENE_VALUES) {
        for (qpathway in Q_PATHWAY_VALUES) {
          key <- paste(mingene, qpathway, sep = "_")
          pathway_counts <- sapply(databases, function(db) {
            count <- pathway_counts_by_db[[db]][[key]]
            if (is.null(count)) count <- 0
            paste0(db, ":", count)
          })
          pathways_str <- paste(pathway_counts, collapse = ", ")

          # Check if any pathway has count > 0
          any_sig <- any(sapply(databases, function(db) {
            count <- pathway_counts_by_db[[db]][[key]]
            !is.null(count) && count > 0
          }))

          run_summary <<- rbind(run_summary, data.frame(
            analysis_name = analysis_name,
            analysis_type = "glm_ora",
            group_col = group_col,
            any_pathway_significant = any_sig,
            pathways_significant = pathways_str,
            n_samples_per_group = n_samples_str,
            genes_tested = genes_tested,
            genes_significant = genes_significant,
            enrichment_ran = TRUE,
            minsample = minsample,
            qgene = qgene,
            mingene = mingene,
            qpathway = qpathway,
            groups = paste(groups, collapse = ","),
            glm_errors = glm_errors,
            posthoc_errors = posthoc_errors,
            posthoc_skipped = posthoc_skipped,
            stringsAsFactors = FALSE
          ))
        }
      }
    }
  }
}

#### GLM-BASED GSEA FUNCTION ####
# Uses z-statistic from GLM as ranking metric for Gene Set Enrichment Analysis
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

  # Calculate sample counts per group for run summary
  actual_sample_col <- if (sample_col %in% colnames(te_data)) sample_col else "sample"
  sample_counts <- te_data %>%
    filter(!is.na(.data[[group_col]])) %>%
    group_by(.data[[group_col]]) %>%
    summarise(n = n_distinct(.data[[actual_sample_col]]), .groups = "drop")
  sample_counts_vec <- setNames(sample_counts$n, as.character(sample_counts[[group_col]]))
  n_samples_str <- paste(names(sample_counts_vec), sample_counts_vec, sep = ":", collapse = ", ")

  # Create analysis-specific directories
  analysis_base_dir <- paste0(base_dir, analysis_name, "/")
  files_dir <- paste0(analysis_base_dir, "files/")
  dir.create(files_dir, showWarnings = FALSE, recursive = TRUE)

  # Create database subdirectories for plots
  for (db in databases) {
    db_plot_dir <- paste0(analysis_base_dir, db, "/")
    dir.create(db_plot_dir, showWarnings = FALSE, recursive = TRUE)
  }

  # Loop over minsample, run GLM ONCE with p_threshold=1.0
  for (minsample in MIN_SAMPLES_VALUES) {
    cat("\n--- minsample =", minsample, "---\n")

    # Run GLM ONCE with p_threshold=1.0 to get ALL results
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

    # Add comparison and enriched_in columns for 2-group comparisons
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

    # Extract metrics from GLM result
    genes_tested <- glm_result$genes_attempted
    glm_errors <- glm_result$glm_errors
    posthoc_errors <- glm_result$posthoc_errors
    posthoc_skipped <- glm_result$posthoc_skipped

    # Initialize pathway_counts_by_db for combined format
    pathway_counts_by_db <- list()
    for (db in databases) {
      pathway_counts_by_db[[db]] <- list()
      for (qpathway in Q_PATHWAY_VALUES) {
        pathway_counts_by_db[[db]][[as.character(qpathway)]] <- 0
      }
    }

    gsea_ran <- FALSE

    # Create gene list for GSEA using z-statistic
    if (is.null(full_results) || nrow(full_results) == 0) {
      cat("  No GLM results available. Skipping GSEA.\n")

      # Add run_summary entries (no qgene, no mingene for GSEA)
      for (qpathway in Q_PATHWAY_VALUES) {
        pathway_counts <- sapply(databases, function(db) paste0(db, ":0"))
        pathways_str <- paste(pathway_counts, collapse = ", ")

        run_summary <<- rbind(run_summary, data.frame(
          analysis_name = analysis_name,
          analysis_type = "glm_gsea",
          group_col = group_col,
          any_pathway_significant = FALSE,
          pathways_significant = pathways_str,
          n_samples_per_group = n_samples_str,
          genes_tested = genes_tested,
          genes_significant = NA_integer_,
          enrichment_ran = FALSE,
          minsample = minsample,
          qgene = NA_real_,
          mingene = NA_integer_,
          qpathway = qpathway,
          groups = paste(groups, collapse = ","),
          glm_errors = glm_errors,
          posthoc_errors = posthoc_errors,
          posthoc_skipped = posthoc_skipped,
          stringsAsFactors = FALSE
        ))
      }
      next
    }

    # For 2-group comparison: use z_value column
    # For 3+ group comparison: need to use pairwise results
    if (n_groups == 2) {
      # Use z_value from main results
      if (!"z_value" %in% colnames(full_results)) {
        cat("  No z_value column found. Skipping GSEA.\n")
        next
      }

      # Create named vector: gene names as names, z-values as values
      gene_list <- setNames(full_results$z_value, full_results$gene)
      gene_list <- gene_list[!is.na(gene_list)]
      gene_list <- sort(gene_list, decreasing = TRUE)

      cat("  Gene list for GSEA:", length(gene_list), "genes\n")

      # Run GSEA for each database
      for (db in databases) {
        cat("\n  Database:", db, "(GSEA)\n")
        db_plot_dir <- paste0(analysis_base_dir, db, "/")
        db_files_dir <- paste0(db_plot_dir, "files/")
        dir.create(db_files_dir, showWarnings = FALSE, recursive = TRUE)

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
            # Convert gene list to Entrez IDs
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

            # Save full GSEA result
            saveRDS(gsea_result, paste0(db_files_dir, "glm_gsea_pathway_", db, "_", analysis_name,
                                        "_minsample", minsample, "_full.rds"))
            write.csv(as.data.frame(gsea_result), paste0(db_files_dir, "glm_gsea_pathway_", db, "_", analysis_name,
                                        "_minsample", minsample, "_full.csv"), row.names = FALSE)

            # Filter by qpathway
            for (qpathway in Q_PATHWAY_VALUES) {
              gsea_df <- as.data.frame(gsea_result)
              filtered_df <- gsea_df %>% filter(qvalue < qpathway)
              pathways_sig <- nrow(filtered_df)

              pathway_counts_by_db[[db]][[as.character(qpathway)]] <- pathways_sig

              if (pathways_sig > 0) {
                prefix <- paste0("glm_gsea_pathway_", db, "_", analysis_name,
                                 "_minsample", minsample, "_qpathway", qpathway)

                # Save filtered CSV
                write.csv(filtered_df, paste0(db_files_dir, prefix, ".csv"), row.names = FALSE)

                # Create GSEA visualizations
                # For binary analysis: groups[1] is coded as 1, groups[2] as 0
                # Positive z-value = higher in groups[1], so "activated" = up in groups[1]
                group1_label <- groups[1]
                group2_label <- groups[2]

                # Dot plot with up/down split
                tryCatch({
                  p_dot <- enrichplot::dotplot(gsea_result,
                                                showCategory = min(20, pathways_sig),
                                                split = ".sign") +
                    ggplot2::facet_grid(~.sign, labeller = ggplot2::labeller(.sign = c(
                      "activated" = paste0("Up in ", group1_label),
                      "suppressed" = paste0("Up in ", group2_label)
                    ))) +
                    ggplot2::theme(strip.text = ggplot2::element_text(size = 10))
                  ggsave(paste0(db_plot_dir, prefix, "_dot.png"), p_dot, width = 14, height = 8)
                }, error = function(e) cat("    Dot plot error:", e$message, "\n"))

                # Ridge plot
                if (pathways_sig >= 3) {
                  tryCatch({
                    p_ridge <- enrichplot::ridgeplot(gsea_result, showCategory = min(20, pathways_sig))
                    ggsave(paste0(db_plot_dir, prefix, "_ridge.png"), p_ridge, width = 10, height = 10)
                  }, error = function(e) cat("    Ridge plot error:", e$message, "\n"))
                }

                # GSEA running score plot for top pathways
                top_pathways <- head(filtered_df$ID, 4)
                for (i in seq_along(top_pathways)) {
                  pathway_id <- top_pathways[i]
                  tryCatch({
                    p_gsea <- enrichplot::gseaplot2(gsea_result, geneSetID = pathway_id)
                    ggsave(paste0(db_plot_dir, prefix, "_gsea_", i, ".png"), p_gsea, width = 10, height = 6)
                  }, error = function(e) cat("    gseaplot2 error for", pathway_id, ":", e$message, "\n"))
                }

                cat("    Saved:", prefix, "(", pathways_sig, "terms)\n")
              }
            }
          } else {
            cat("    No enriched terms found\n")
          }
        }, error = function(e) {
          cat("    Error running GSEA for", db, ":", e$message, "\n")
        })
      }
    } else {
      # 3+ groups: run GSEA for each pairwise contrast
      if (!is.null(pairwise_results) && nrow(pairwise_results) > 0) {
        # Save post-hoc pairwise results (same as GLM ORA does)
        posthoc_file <- paste0(files_dir, "glm_", analysis_name,
                               "_posthoc_minsample", minsample, ".csv")
        write.csv(pairwise_results, posthoc_file, row.names = FALSE)
        cat("  Saved post-hoc pairwise results:", basename(posthoc_file), "\n")

        # Get unique contrasts
        contrasts <- unique(pairwise_results$contrast)
        cat("  Running GSEA for", length(contrasts), "pairwise contrasts\n")

        for (contrast_name in contrasts) {
          cat("\n  Contrast:", contrast_name, "\n")

          contrast_results <- pairwise_results %>% filter(contrast == contrast_name)

          # Use z.ratio as the ranking metric (from emmeans)
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

          # Clean contrast name for file naming
          contrast_clean <- gsub(" - ", "_vs_", contrast_name)
          contrast_clean <- gsub("[^a-zA-Z0-9_]", "", contrast_clean)

          for (db in databases) {
            cat("    Database:", db, "\n")
            db_plot_dir <- paste0(analysis_base_dir, db, "/")
            db_files_dir <- paste0(db_plot_dir, "files/")
            dir.create(db_files_dir, showWarnings = FALSE, recursive = TRUE)

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

                # Save full GSEA result for this contrast
                saveRDS(gsea_result, paste0(db_files_dir, "glm_gsea_pathway_", db, "_", analysis_name,
                                            "_", contrast_clean, "_minsample", minsample, "_full.rds"))
                write.csv(as.data.frame(gsea_result), paste0(db_files_dir, "glm_gsea_pathway_", db, "_", analysis_name,
                                            "_", contrast_clean, "_minsample", minsample, "_full.csv"), row.names = FALSE)

                # Filter by qpathway and count
                for (qpathway in Q_PATHWAY_VALUES) {
                  gsea_df <- as.data.frame(gsea_result)
                  filtered_df <- gsea_df %>% filter(qvalue < qpathway)
                  pathways_sig <- nrow(filtered_df)

                  # Add to counts (sum across contrasts)
                  current <- pathway_counts_by_db[[db]][[as.character(qpathway)]]
                  pathway_counts_by_db[[db]][[as.character(qpathway)]] <- current + pathways_sig

                  if (pathways_sig > 0) {
                    prefix <- paste0("glm_gsea_pathway_", db, "_", analysis_name,
                                     "_", contrast_clean, "_minsample", minsample, "_qpathway", qpathway)

                    write.csv(filtered_df, paste0(db_files_dir, prefix, ".csv"), row.names = FALSE)

                    # Create GSEA visualizations
                    # Parse contrast name for dynamic facet labels
                    contrast_parts <- strsplit(contrast_name, " - ")[[1]]
                    group1_label <- trimws(contrast_parts[1])
                    group2_label <- if(length(contrast_parts) > 1) trimws(contrast_parts[2]) else "Reference"

                    # Dot plot with up/down split
                    tryCatch({
                      p_dot <- enrichplot::dotplot(gsea_result,
                                                    showCategory = min(20, pathways_sig),
                                                    split = ".sign") +
                        ggplot2::facet_grid(~.sign, labeller = ggplot2::labeller(.sign = c(
                          "activated" = paste0("Up in ", group1_label),
                          "suppressed" = paste0("Up in ", group2_label)
                        ))) +
                        ggplot2::theme(strip.text = ggplot2::element_text(size = 10))
                      ggsave(paste0(db_plot_dir, prefix, "_dot.png"), p_dot, width = 14, height = 8)
                    }, error = function(e) cat("      Dot plot error:", e$message, "\n"))

                    # Ridge plot
                    if (pathways_sig >= 3) {
                      tryCatch({
                        p_ridge <- enrichplot::ridgeplot(gsea_result, showCategory = min(20, pathways_sig))
                        ggsave(paste0(db_plot_dir, prefix, "_ridge.png"), p_ridge, width = 10, height = 10)
                      }, error = function(e) cat("      Ridge plot error:", e$message, "\n"))
                    }

                    # GSEA running score plot for top pathways
                    top_pathways <- head(filtered_df$ID, 4)
                    for (i in seq_along(top_pathways)) {
                      pathway_id <- top_pathways[i]
                      tryCatch({
                        p_gsea <- enrichplot::gseaplot2(gsea_result, geneSetID = pathway_id)
                        ggsave(paste0(db_plot_dir, prefix, "_gsea_", i, ".png"), p_gsea, width = 10, height = 6)
                      }, error = function(e) cat("      gseaplot2 error for", pathway_id, ":", e$message, "\n"))
                    }
                  }
                }
              }
            }, error = function(e) {
              cat("      Error running GSEA:", e$message, "\n")
            })
          }
        }
      } else {
        cat("  No pairwise results available for GSEA\n")
      }
    }

    # Add run_summary entries (one per qpathway)
    for (qpathway in Q_PATHWAY_VALUES) {
      pathway_counts <- sapply(databases, function(db) {
        count <- pathway_counts_by_db[[db]][[as.character(qpathway)]]
        if (is.null(count)) count <- 0
        paste0(db, ":", count)
      })
      pathways_str <- paste(pathway_counts, collapse = ", ")

      # Check if any pathway has count > 0
      any_sig <- any(sapply(databases, function(db) {
        count <- pathway_counts_by_db[[db]][[as.character(qpathway)]]
        !is.null(count) && count > 0
      }))

      run_summary <<- rbind(run_summary, data.frame(
        analysis_name = analysis_name,
        analysis_type = "glm_gsea",
        group_col = group_col,
        any_pathway_significant = any_sig,
        pathways_significant = pathways_str,
        n_samples_per_group = n_samples_str,
        genes_tested = genes_tested,
        genes_significant = NA_integer_,
        enrichment_ran = gsea_ran,
        minsample = minsample,
        qgene = NA_real_,
        mingene = NA_integer_,
        qpathway = qpathway,
        groups = paste(groups, collapse = ","),
        glm_errors = glm_errors,
        posthoc_errors = posthoc_errors,
        posthoc_skipped = posthoc_skipped,
        stringsAsFactors = FALSE
      ))
    }
  }
}

#### PATHWAY ANALYSIS 1: KICS COHORT ####
# Note: This is a cohort-level analysis without group comparison
# Only simple ORA is applicable (no GLM-ORA or GLM-GSEA)
write_output(quote(NULL), "Pathway Analysis - KICS Cohort (Simple ORA)")

cat("Dataset: te_kics_split_t\n")
cat("Samples:", length(unique(te_kics_split_t$sample)), "\n")
cat("Genes:", length(unique(te_kics_split_t$Gene_name)), "\n\n")

# For cohort-level analysis, we just run simple ORA on genes in >= minsample samples
for (minsample in MIN_SAMPLES_VALUES) {
  cat("\n--- KICS cohort analysis with minsample =", minsample, "---\n")

  # Get genes present in >= minsample samples
  gene_counts <- te_kics_split_t %>%
    group_by(Gene_name) %>%
    summarise(n_samples = n_distinct(sample), .groups = "drop") %>%
    filter(n_samples >= minsample)

  kics_genes <- unique(gene_counts$Gene_name)
  cat("Genes in >=", minsample, "samples:", length(kics_genes), "\n")

  if (length(kics_genes) >= 10) {
    # Create analysis-level directories for kics_cohort in ora/
    kics_analysis_dir <- paste0(ora_dir, "kics_cohort/")

    for (db in DATABASES_TO_RUN) {
      cat("\n  Database:", db, "\n")

      # Create database-specific directories
      db_dir <- paste0(kics_analysis_dir, db, "/")
      db_files_dir <- paste0(db_dir, "files/")
      dir.create(db_files_dir, showWarnings = FALSE, recursive = TRUE)

      tryCatch({
        ora_result <- perform_ora_single_database(
          genes = kics_genes,
          database = db,
          pvalueCutoff = 1,
          qvalueCutoff = 1,
          universe = BACKGROUND_GENES
        )

        if (is.null(ora_result) || nrow(as.data.frame(ora_result)) == 0) {
          cat("    No enriched terms found\n")
          next
        }

        # Save full ORA result object
        saveRDS(ora_result, paste0(db_files_dir, "ora_", db, "_kics_cohort",
                                "_minsample", minsample, "_full.rds"))
        write.csv(as.data.frame(ora_result), paste0(db_files_dir, "ora_", db, "_kics_cohort",
                                "_minsample", minsample, "_full.csv"), row.names = FALSE)
        cat("    Saved full result:", nrow(as.data.frame(ora_result)), "terms\n")

        # Post-hoc filter by mingene and qpathway
        for (mingene in MIN_GENE_VALUES) {
          for (qpathway in Q_PATHWAY_VALUES) {
            ora_df <- as.data.frame(ora_result)
            filtered_df <- ora_df %>%
              filter(Count >= mingene, qvalue < qpathway)

            pathways_sig <- nrow(filtered_df)
            if (pathways_sig == 0) next

            prefix <- paste0("ora_", db, "_kics_cohort",
                             "_minsample", minsample,
                             "_mingene", mingene, "_qpathway", qpathway)

            # Save filtered CSV
            write.csv(filtered_df, paste0(db_files_dir, prefix, ".csv"), row.names = FALSE)

            # Create plots
            tryCatch({
              p_dot <- dotplot(ora_result, showCategory = min(20, pathways_sig))
              ggsave(paste0(db_dir, prefix, "_dot.png"), p_dot, width = 12, height = 8)

              p_bar <- barplot(ora_result, showCategory = min(20, pathways_sig))
              ggsave(paste0(db_dir, prefix, "_bar.png"), p_bar, width = 10, height = 8)
            }, error = function(e) NULL)

            cat("    Saved:", prefix, "(", pathways_sig, "terms)\n")
          }
        }
      }, error = function(e) {
        cat("    Error running ORA for", db, ":", e$message, "\n")
      })
    }
  }
}

#### PATHWAY ANALYSIS 2: TP53 STATUS ####
write_output(quote(NULL), "Pathway Analysis - TP53 Status")

cat("Dataset: te_aff_split_t\n")
cat("Samples:", length(unique(te_aff_split_t$sample)), "\n")
tp53_groups <- unique(te_aff_split_t$TP53_status)
tp53_groups <- tp53_groups[!is.na(tp53_groups)]
cat("TP53 groups:", paste(tp53_groups, collapse = ", "), "\n\n")

if (length(tp53_groups) >= 2) {
  # Run all 3 analysis types
  run_simple_ora_analysis(
    te_data = te_aff_split_t,
    group_col = "TP53_status",
    groups = tp53_groups,
    analysis_name = "tp53_status",
    base_dir = ora_dir
  )
  run_glm_ora_analysis(
    te_data = te_aff_split_t,
    group_col = "TP53_status",
    groups = tp53_groups,
    analysis_name = "tp53_status",
    base_dir = glm_ora_dir
  )
  run_glm_gsea_analysis(
    te_data = te_aff_split_t,
    group_col = "TP53_status",
    groups = tp53_groups,
    analysis_name = "tp53_status",
    base_dir = glm_gsea_dir
  )
} else {
  cat("Not enough TP53 groups. Skipping.\n")
}

#### PATHWAY ANALYSIS 3: TP53 3-LEVEL ####
write_output(quote(NULL), "Pathway Analysis - TP53 3-Level")

cat("Dataset: te_aff_split_t (TP53_3level)\n")

# Filter to samples with valid TP53_3level
te_aff_3level <- te_aff_split_t %>% filter(!is.na(TP53_3level))
if (nrow(te_aff_3level) > 0) {
  tp53_3level_groups <- unique(te_aff_3level$TP53_3level)
  tp53_3level_groups <- tp53_3level_groups[!is.na(tp53_3level_groups)]
  cat("TP53_3level groups:", paste(tp53_3level_groups, collapse = ", "), "\n")

  if (length(tp53_3level_groups) >= 2) {
    # Run all 3 analysis types
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
  } else {
    cat("Not enough TP53_3level groups. Skipping.\n")
  }
}

#### PATHWAY ANALYSIS 4: TUMOR TYPE (KICS) ####
write_output(quote(NULL), "Pathway Analysis - KICS by Tumor Type")

if ("tumor_type" %in% colnames(te_kics_split_t)) {
  # Filter to tumor types with >= 3 samples
  tumor_type_counts <- table(te_kics_split_t$tumor_type)
  valid_tumor_types <- names(tumor_type_counts[tumor_type_counts >= 3])
  cat("Tumor types with >= 3 samples:", paste(valid_tumor_types, collapse = ", "), "\n")

  if (length(valid_tumor_types) >= 2) {
    te_kics_tumor <- te_kics_split_t %>%
      filter(tumor_type %in% valid_tumor_types)

    # Run all 3 analysis types
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
} else {
  cat("Warning: tumor_type column not found in te_kics_split_t. Skipping tumor type analysis.\n")
}

#### PATHWAY ANALYSIS 5: KICS BY ANCESTRY ####
write_output(quote(NULL), "Pathway Analysis - KICS by Ancestry")

if ("predicted_ancestry_thres" %in% colnames(te_kics_split_t)) {
  ancestries <- unique(te_kics_split_t$predicted_ancestry_thres)
  ancestries <- ancestries[!is.na(ancestries)]
  cat("Ancestries:", paste(ancestries, collapse = ", "), "\n")

  # Only run if we have enough groups with samples
  ancestry_counts <- table(te_kics_split_t$predicted_ancestry_thres)
  valid_ancestries <- names(ancestry_counts[ancestry_counts >= 10])

  if (length(valid_ancestries) >= 2) {
    te_kics_ancestry <- te_kics_split_t %>%
      filter(predicted_ancestry_thres %in% valid_ancestries)

    # Run all 3 analysis types
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
  } else {
    cat("Not enough ancestry groups with >=10 samples. Skipping.\n")
  }
}

#### PATHWAY ANALYSIS 6: LFS BY COHORT ####
write_output(quote(NULL), "Pathway Analysis - LFS by Cohort")

if (exists("te_lfs_split_t") && "cohort" %in% colnames(te_lfs_split_t)) {
  cohorts <- unique(te_lfs_split_t$cohort)
  cohorts <- cohorts[!is.na(cohorts)]
  cat("Cohorts:", paste(cohorts, collapse = ", "), "\n")

  if (length(cohorts) >= 2) {
    # Run all 3 analysis types
    run_simple_ora_analysis(
      te_data = te_lfs_split_t,
      group_col = "cohort",
      groups = cohorts,
      analysis_name = "lfs_cohort",
      base_dir = ora_dir
    )
    run_glm_ora_analysis(
      te_data = te_lfs_split_t,
      group_col = "cohort",
      groups = cohorts,
      analysis_name = "lfs_cohort",
      base_dir = glm_ora_dir
    )
    run_glm_gsea_analysis(
      te_data = te_lfs_split_t,
      group_col = "cohort",
      groups = cohorts,
      analysis_name = "lfs_cohort",
      base_dir = glm_gsea_dir
    )
  }
}

#### PATHWAY ANALYSIS 7: TAYLOR BY SUBTYPE ####
write_output(quote(NULL), "Pathway Analysis - Taylor by Subtype")

if (exists("te_taylor_split_t") && nrow(te_taylor_split_t) > 0) {
  if ("tumor_type_subclass" %in% colnames(te_taylor_split_t)) {
    # Filter to subtypes with >= 3 samples
    subtype_counts <- table(te_taylor_split_t$tumor_type_subclass)
    valid_subtypes <- names(subtype_counts[subtype_counts >= 3])
    cat("Taylor subtypes with >= 3 samples:", paste(valid_subtypes, collapse = ", "), "\n")

    if (length(valid_subtypes) >= 2) {
      te_taylor_subtype <- te_taylor_split_t %>%
        filter(tumor_type_subclass %in% valid_subtypes)

      # Simple ORA
      run_simple_ora_analysis(
        te_data = te_taylor_subtype,
        group_col = "tumor_type_subclass",
        groups = valid_subtypes,
        analysis_name = "taylor_subtype",
        base_dir = ora_dir
      )

      # GLM-ORA (no covariates for Taylor - external cohort)
      run_glm_ora_analysis(
        te_data = te_taylor_subtype,
        group_col = "tumor_type_subclass",
        groups = valid_subtypes,
        analysis_name = "taylor_subtype",
        base_dir = glm_ora_dir,
        covariates = NULL
      )

      # GLM-GSEA
      run_glm_gsea_analysis(
        te_data = te_taylor_subtype,
        group_col = "tumor_type_subclass",
        groups = valid_subtypes,
        analysis_name = "taylor_subtype",
        base_dir = glm_gsea_dir,
        covariates = NULL
      )
    } else {
      cat("Not enough subtypes with >= 3 samples. Skipping.\n")
    }
  } else {
    cat("tumor_type_subclass column not found in te_taylor_split_t\n")
  }
} else {
  cat("te_taylor_split_t not available for Taylor subtype analysis\n")
}

#### P53 FITNESS ANALYSIS ####
write_output(quote(NULL), "P53 Fitness Analysis")

tp53_fitness <- fread("/Users/briannelaverty/Documents/R_Malkin/TE/data/tp53_variants_predictions.tsv", sep="\t", header=TRUE)
cat("Adding fitness values to affected and LFS data...\n")
te_aff_t_tp53 <- prep_p53_fitness(te_aff_t, tp53_fitness)
te_lfs_t_tp53 <- prep_p53_fitness(te_lfs_t, tp53_fitness)

cat("Scatter plot of TE count vs p53 fitness (affected)...\n")
for (i in seq_along(types)) {
  plot_title <- paste0("TE count vs p53 fitness (", ifelse(is.na(types[i]), "all types", types[i]), ")")
  dependent_var <- ifelse(is.na(types[i]), "total", types[i])
  tryCatch({
    p <- scatter_template_lm_one_variable(te_aff_t_tp53, "p53_fitness", dependent_var, c(0.25, 1), c(0, 1000), "p53 fitness", "Repeat count", shape_column = NULL, colour_column = NULL)
    titled_print(p, plot_title)
    ggsave(paste0(plot_dir, "counts_clinical_lfs/te_count_lfs_fitness_", ifelse(is.na(types[i]), "all", types[i]), ".png"), plot = p, width = 9, height = 5)
  }, error = function(e) {
    cat("Skipping fitness plot for type", ifelse(is.na(types[i]), "all", types[i]), "due to error:", e$message, "\n")
  })
}

cat("\n===== P53 FITNESS BY INHERITANCE =====\n")
p_p53_fitness_inheritance <- plot_box_kruskal(te_lfs_t_tp53, "p53_fitness", "inheritance", "Inheritance", "p53 fitness")
titled_print(p_p53_fitness_inheritance, "p53 fitness by inheritance (LFS)")
ggsave(paste0(plot_dir, "counts_clinical_lfs/p53_fitness_by_inheritance_lfs.png"), plot = p_p53_fitness_inheritance, width = 5, height = 5)

#### CANCER GENES ANALYSIS ####
write_output(quote(NULL), "Cancer Genes Analysis")

cpg <- read.csv("/Users/briannelaverty/Documents/R_Malkin/te/data/raw/kics_cpg.txt", header=FALSE, sep="\t")
genes <- cpg$V1

# Cancer genes affected (requires nohits file from descriptive module)
nohits_file <- paste0(r_dir, "nohits_final_te_count_t_te_aff_selected_t", ".RData")
if (file.exists(nohits_file)) {
  load(nohits_file)
  te_aff_split_cancergenes_t <- te_aff_split_t %>% filter(Gene_name %in% genes)
  te_aff_split_cancergenes_processed_t <- process_all_combinations(te_aff_split_cancergenes_t)
  te_aff_split_cancergenes_processed_t <- as.data.frame(add_nohit_samples(te_aff_split_cancergenes_processed_t, nohits))
  write_output(quote(summary(te_aff_split_cancergenes_processed_t$total)), "Summary of cancer genes (affected)")
  te_aff_split_cancergenes_processed_t <- merge_dfs(te_aff_split_cancergenes_processed_t, clinical, include_all_x = FALSE)

  cat("Plotting cancer genes affected by TP53 status (Wilcoxon test)...\n")
  for (i in seq_along(types)) {
    plot_title <- paste0(ifelse(is.na(types[i]), "All types", types[i]), " cancer genes affected by TP53 status")
    write_output(
      quote(plot_count_wilcox(te_aff_split_cancergenes_processed_t, type=types[i], group="TP53_status", y_lab="Cancer genes affected", chr=NA, x_lab=x_tp53, log_scale=FALSE)),
      paste0("Cancer genes affected by TP53 status (type=", ifelse(is.na(types[i]), "all", types[i]), ")")
    )
    p <- plot_count_wilcox(te_aff_split_cancergenes_processed_t, type=types[i], group="TP53_status", y_lab="Cancer genes affected", chr=NA, x_lab=x_tp53, log_scale=FALSE)
    titled_print(p, plot_title)
    ggsave(paste0(plot_dir, "cancer_genes/te_cancergenes_", ifelse(is.na(types[i]), "all", types[i]), "_tp53.png"), plot = p, width = 5, height = 5)
  }
} else {
  cat("Warning: nohits file not found. Skipping cancer genes affected analysis.\n")
  cat("  File:", nohits_file, "\n")
  cat("  Run 02_te_viz_tumour_01_descriptive.R first to generate this file.\n")
}

#### TOP MUTATED CANCER GENES ####
write_output(quote(NULL), "Top Mutated Cancer Genes")

# Add gene size to KICS split data
tryCatch({
  te_kics_split_t <- add_gene_size_todf(te_kics_split_t, gene_size)
}, error = function(e) {
  cat("Warning: Could not add gene size data:", e$message, "\n")
})

# Add a dummy gene_size column if it doesn't exist
if (!"gene_size" %in% colnames(te_kics_split_t)) {
  te_kics_split_t$gene_size <- 1
}

geneList_kics_t <- calculate_gene_mut_sample_frequency(te_kics_split_t, filter_exon = FALSE)
geneList_kics_filtered_t <- geneList_kics_t[geneList_kics_t$Gene_name %in% genes, ]

# Plot top mutated cancer genes (KICS)
if (nrow(geneList_kics_filtered_t) > 0) {
  write_output(quote(plot_top_genes(geneList_kics_filtered_t, column="freq_samples_effected", label_column=NULL, x_lab="Proportion of samples with insertion in gene", top_n=20)), "Top mutated cancer genes (KICS, raw frequency)")
  p_cancer1 <- plot_top_genes(geneList_kics_filtered_t, column="freq_samples_effected", label_column=NULL, x_lab="Proportion of samples with insertion in gene", top_n=20)
  titled_print(p_cancer1, "Top mutated cancer genes (KICS, raw frequency)")
  ggsave(paste0(plot_dir, "cancer_genes/te_top_cancer_genes_affected_kics.png"), plot = p_cancer1, width = 9, height = 5)
  write_output(quote(plot_top_genes(geneList_kics_filtered_t, column="freq_samples_effected_normalized", label_column="num_samples_effected", x_lab="Proportion of samples with insertion in gene normalized for gene size", top_n=20)), "Top mutated cancer genes (KICS, normalized for gene size)")
  p_cancer2 <- plot_top_genes(geneList_kics_filtered_t, column="freq_samples_effected_normalized", label_column="num_samples_effected", x_lab="Proportion of samples with insertion in gene normalized for gene size", top_n=20)
  titled_print(p_cancer2, "Top mutated cancer genes (KICS, normalized for gene size)")
  ggsave(paste0(plot_dir, "cancer_genes/te_top_cancer_genes_affected_normalized_kics.png"), plot = p_cancer2, width = 9, height = 5)
} else {
  cat("Warning: No cancer genes found in KICS data (TEST_MODE may have insufficient data)\n")
}

# Stacked plot for top genes (KICS)
write_output(
  quote(plot_gene_effects(te_kics_split_t, min_sample_tt = 5, remove_other = FALSE, top_n_genes = 15, SV_type = NULL)),
  "Top Affected Genes (KICS)"
)

# Genes present that meet criteria
location <- c("CDS", "UTR", "3'UTR", "5'UTR", "5'UTR-CDS", "CDS-3'UTR")
genes_aff_t <- identify_te_genes(te_aff_split_genes_t, gene_vector=genes, location="exon", location2=NULL)
genes_kics_t <- identify_te_genes(te_kics_split_t, gene_vector=genes, location=NULL, location2=NULL)
genes_lfs_t <- identify_te_genes(te_lfs_split_t, gene_vector=genes, location="exon", location2=NULL)

write_output(quote(table(genes_aff_t$sample, genes_aff_t$Gene_name)), "samples and genes with exon insertion")

# Wilcox test between affected and controls for cancer genes
write_output(quote(count_location_wilcox(te_aff_split_t, gene=genes, filter_element=NA, group="TP53_status", location=NULL)), "Wilcox test: affected vs controls (cancer genes)")

#### TOP MUTATED GENES ####
write_output(quote(NULL), "Top Mutated Genes (All Genes)")

# KICS
tryCatch({
  te_kics_split_t_temp <- add_gene_size_todf(te_kics_split_t, gene_size)
  if ("gene_size" %in% colnames(te_kics_split_t_temp)) {
    te_kics_split_t <- te_kics_split_t_temp
  }
}, error = function(e) {
  cat("Warning: Could not add gene size data to KICS:", e$message, "\n")
})

if (!"gene_size" %in% colnames(te_kics_split_t)) {
  cat("Adding dummy gene_size column to KICS data\n")
  te_kics_split_t$gene_size <- 1
}

tryCatch({
  geneList_kics_t <- calculate_gene_mut_sample_frequency(te_kics_split_t)
  p1 <- plot_top_genes(geneList_kics_t, column="freq_samples_effected", label_column=NULL, x_lab="Proportion of samples with insertion in gene", top_n=20)
  titled_print(p1, "Top mutated genes (KICS, raw frequency)")
  ggsave(paste0(plot_dir, "cancer_genes/te_top_genes_affected_kics.png"), plot = p1, width = 9, height = 5)
  p2 <- plot_top_genes(geneList_kics_t, column="freq_samples_effected_normalized", label_column="num_samples_effected", x_lab="Proportion of samples with insertion in gene normalized by gene size", top_n=20)
  titled_print(p2, "Top mutated genes (KICS, normalized for gene size)")
  ggsave(paste0(plot_dir, "cancer_genes/te_top_genes_affected_normalized_kics.png"), plot = p2, width = 9, height = 5)
}, error = function(e) {
  cat("Warning: Could not calculate gene frequencies for KICS:", e$message, "\n")
})

# Number of genes affecting >50%, 75%, 90% of samples (KICS)
write_output(quote(count_genes_by_threshold(geneList_kics_t, column_name = "freq_samples_effected", thresholds = c(0.5, 0.75, 0.9))), "Number of genes affecting >50%, 75%, 90% of samples (KICS)")

# LFS
tryCatch({
  te_lfs_split_t_temp <- add_gene_size_todf(te_lfs_split_t, gene_size)
  if ("gene_size" %in% colnames(te_lfs_split_t_temp)) {
    te_lfs_split_t <- te_lfs_split_t_temp
  }
}, error = function(e) {
  cat("Warning: Could not add gene size data to LFS:", e$message, "\n")
})

if (!"gene_size" %in% colnames(te_lfs_split_t)) {
  cat("Adding dummy gene_size column to LFS data\n")
  te_lfs_split_t$gene_size <- 1
}

tryCatch({
  geneList_lfs_t <- calculate_gene_mut_sample_frequency(te_lfs_split_t, nsample_thresh = 0, filter_exon = FALSE)
  p3 <- plot_top_genes(geneList_lfs_t, column="freq_samples_effected", label_column=NULL, x_lab="Proportion of samples with insertion in gene", top_n=20)
  titled_print(p3, "Top mutated genes (LFS, raw frequency)")
  ggsave(paste0(plot_dir, "cancer_genes/te_top_genes_affected_lfs.png"), plot = p3, width = 9, height = 5)
  p4 <- plot_top_genes(geneList_lfs_t, column="freq_samples_effected_normalized", label_column="freq_samples_effected", x_lab="Proportion of samples with insertion in gene normalized by gene size", top_n=20)
  titled_print(p4, "Top mutated genes (LFS, normalized for gene size)")
  ggsave(paste0(plot_dir, "cancer_genes/te_top_genes_affected_normalized_lfs.png"), plot = p4, width = 9, height = 5)
}, error = function(e) {
  cat("Warning: Could not calculate gene frequencies for LFS:", e$message, "\n")
})

# Control
te_control_split_t <- te_aff_split_t %>% filter(TP53_status == "Control")
if (nrow(te_control_split_t) > 0) {
  if (!"gene_size" %in% colnames(te_control_split_t)) {
    cat("Adding dummy gene_size column to Control data\n")
    te_control_split_t$gene_size <- 1
  }
  tryCatch({
    geneList_control_t <- calculate_gene_mut_sample_frequency(te_control_split_t, nsample_thresh = 0, filter_exon = FALSE)
    p5 <- plot_top_genes(geneList_control_t, top_n=20)
    titled_print(p5, "Top mutated genes (Control, raw frequency)")
    ggsave(paste0(plot_dir, "cancer_genes/te_top_genes_affected_control.png"), plot = p5, width = 9, height = 5)
  }, error = function(e) {
    cat("Warning: Could not calculate gene frequencies for Control:", e$message, "\n")
  })
} else {
  cat("No control samples found, skipping control analysis\n")
}

#### SUMMARY ####
write_output(quote(NULL), "Pathway Analysis Summary")

cat("Analysis completed for the following datasets:\n")
cat("  1. KICS cohort (general ORA)\n")
cat("  2. TP53 status (GLM 2-group: Mutant vs WT)\n")
cat("  3. TP53 3-level (GLM LRT + emmeans: Germline/Somatic/WT) [NEW]\n")
cat("  4. KICS by tumor type (GLM LRT + emmeans)\n")
cat("  5. KICS by ancestry (GLM LRT + emmeans) [NEW]\n")
cat("  6. LFS by cohort (GLM LRT + emmeans) [NEW]\n\n")

cat("Databases used:", paste(DATABASES_TO_RUN, collapse = ", "), "\n")
cat("Min samples:", paste(MIN_SAMPLES_VALUES, collapse = ", "), "\n")
cat("Q-gene thresholds:", paste(Q_GENE_VALUES, collapse = ", "), "\n")
cat("Min gene counts:", paste(MIN_GENE_VALUES, collapse = ", "), "\n")
cat("Q-pathway thresholds:", paste(Q_PATHWAY_VALUES, collapse = ", "), "\n")

cat("\nOutput structure:\n")
cat("  - ora/: Simple ORA (minsample filter only)\n")
cat("  - glm_ora/: GLM-based ORA\n")
cat("  - glm_gsea/: GLM-based GSEA\n")
cat("Output location:", pathway_base_dir, "\n")

# Save run summary - sort by analysis_name, analysis_type, minsample, qgene, mingene, qpathway
if (nrow(run_summary) > 0) {
  run_summary <- run_summary %>%
    arrange(analysis_name, analysis_type, minsample, qgene, mingene, qpathway)
}
run_summary_file <- paste0(pathway_base_dir, "RUN_SUMMARY.csv")
write.csv(run_summary, run_summary_file, row.names = FALSE)
cat("\nSaved run summary:", run_summary_file, "\n")
cat("  - Total rows:", nrow(run_summary), "\n")
cat("  - Analysis types:", paste(unique(run_summary$analysis_type), collapse = ", "), "\n")

cat("\n===== ADDITIONAL ANALYSES =====\n")
cat("P53 fitness analysis completed\n")
cat("Cancer genes analysis completed\n")
cat("Top mutated genes analysis completed\n")

cat("\nScript completed successfully\n")

# Close module-specific sink
close_module_sink()
