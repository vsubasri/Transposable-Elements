#!/usr/bin/env Rscript

# Germline TE Visualization - Regulatory Elements (Multi-Database ORA)
# GLM-based differential incidence testing with covariate control
# Databases: GO_BP, Reactome, MSigDB_Hallmark, MSigDB_Oncogenic (no KEGG)
# Note: RE-RNA analysis moved to 02_te_viz_germline_08_re_rna.R

# Source common setup and load data
source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/00_viz_common_setup.R")
REQUIRED_DATA <- c("count_matrix", "expand", "split", "clinical")
source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/00_viz_load_data_germline.R")

# Load emmeans for post-hoc contrasts
library(emmeans)

# Initialize module-specific text output
init_module_sink(paste0(plot_dir, "reg_element/"), "RE")

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

cat("Running 02_te_viz_germline_07_re.R...\n")
cat("Multi-database ORA analysis with GLM-based gene filtering\n")
cat("Databases: GO_BP, Reactome, Hallmark, Oncogenic\n\n")

#### CONFIGURATION ####

# Databases to run (no KEGG per plan)
DATABASES_TO_RUN <- c("GO_BP", "Reactome", "Hallmark", "Oncogenic")

# Parameter settings (per plan)
MIN_SAMPLES_VALUES <- c(3, 5)
Q_GENE_VALUES <- c(0.05, 0.1, 0.25)  # For GLM gene filtering
MIN_GENE_VALUES <- c(3, 5)           # Min genes hitting pathway
Q_PATHWAY_VALUES <- c(0.05, 0.1, 0.25)  # For pathway filtering

# Covariate vectors (covar_med, covar_no_ancestry, covar_no_age, covar_no_tumour_type)
# are defined in 00_viz_load_data_germline.R

# Load background regulatory genes
background_genes_file <- "/Users/briannelaverty/Documents/R_Malkin/TE/data/final/background_regulatory_genes.txt"
if (file.exists(background_genes_file)) {
  BACKGROUND_GENES <- readLines(background_genes_file)
  cat("Loaded", length(BACKGROUND_GENES), "background regulatory genes for ORA universe\n")
} else {
  # Fall back to regular background genes
  background_genes_file <- "/Users/briannelaverty/Documents/R_Malkin/TE/data/final/background_genes.txt"
  if (file.exists(background_genes_file)) {
    BACKGROUND_GENES <- readLines(background_genes_file)
    cat("Loaded", length(BACKGROUND_GENES), "background genes for ORA universe\n")
  } else {
    BACKGROUND_GENES <- NULL
    cat("Warning: Background genes file not found. ORA will use default universe.\n")
  }
}

# Create base output directories for three analysis types
re_dir <- paste0(plot_dir, "reg_element/")
ora_dir <- paste0(re_dir, "ora/")
glm_ora_dir <- paste0(re_dir, "glm_ora/")
glm_gsea_dir <- paste0(re_dir, "glm_gsea/")

dir.create(re_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(ora_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(glm_ora_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(glm_gsea_dir, showWarnings = FALSE, recursive = TRUE)

cat("Base output directory:", re_dir, "\n")
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
if (exists("te_lfs_split")) {
  te_lfs_split <- add_tp53_3level(te_lfs_split)
}

#### LOAD AND PROCESS RE DATA ####
write_output(quote(NULL), "Loading Regulatory Elements Data")

re_germline_path <- "/Users/briannelaverty/Documents/R_Malkin/TE/data/final/germline_annotSV_output.SV_RE_intersect.report"

cat("Loading and joining RE data for all cohorts...\n")

# Helper function to fix sample column after join (sample.x becomes sample)
fix_sample_column <- function(df) {
  if ("sample.x" %in% colnames(df) && !"sample" %in% colnames(df)) {
    df <- df %>% rename(sample = sample.x)
  }
  # Also remove sample.y if it exists (from RE report)
  if ("sample.y" %in% colnames(df)) {
    df <- df %>% select(-sample.y)
  }
  return(df)
}

# Affected cohort
te_aff_re <- load_and_join_re_data(te_aff_expand, re_germline_path)
te_aff_re_split <- split_re_genes(te_aff_re) %>% fix_sample_column()
te_aff_re_split <- add_tp53_3level(te_aff_re_split)
cat("Affected RE data:", nrow(te_aff_re_split), "rows,", length(unique(te_aff_re_split$gene_reg)), "genes\n")

# LFS cohort
te_lfs_re <- load_and_join_re_data(te_lfs_expand, re_germline_path)
te_lfs_re_split <- split_re_genes(te_lfs_re) %>% fix_sample_column()
te_lfs_re_split <- add_tp53_3level(te_lfs_re_split)
cat("LFS RE data:", nrow(te_lfs_re_split), "rows,", length(unique(te_lfs_re_split$gene_reg)), "genes\n")

# KICS + HostSeq cohort
te_kics_hostseq_re <- load_and_join_re_data(te_kics_hostseq_expand, re_germline_path)
te_kics_hostseq_re_split <- split_re_genes(te_kics_hostseq_re) %>% fix_sample_column()
te_kics_hostseq_re_split <- add_tp53_3level(te_kics_hostseq_re_split)
cat("KICS+HostSeq RE data:", nrow(te_kics_hostseq_re_split), "rows,", length(unique(te_kics_hostseq_re_split$gene_reg)), "genes\n")

# KICS cohort
te_kics_re <- load_and_join_re_data(te_kics_expand, re_germline_path)
te_kics_re_split <- split_re_genes(te_kics_re) %>% fix_sample_column()
te_kics_re_split <- add_tp53_3level(te_kics_re_split)
cat("KICS RE data:", nrow(te_kics_re_split), "rows,", length(unique(te_kics_re_split$gene_reg)), "genes\n")

# Taylor cohort
if (exists("te_taylor_expand") && nrow(te_taylor_expand) > 0) {
  te_taylor_re <- load_and_join_re_data(te_taylor_expand, re_germline_path)
  te_taylor_re_split <- split_re_genes(te_taylor_re) %>% fix_sample_column()
  cat("Taylor RE data:", nrow(te_taylor_re_split), "rows,", length(unique(te_taylor_re_split$gene_reg)), "genes\n")
} else {
  cat("Taylor RE data: not available\n")
}
cat("\n")

#### HELPER FUNCTIONS ####

# Helper: Create ORA visualizations for a filtered result
# ora_result: the original enrichResult/compareClusterResult object
# filtered_df: the post-hoc filtered dataframe (by mingene, qpathway)
create_ora_visualizations <- function(ora_result, filtered_df, te_data, db, prefix, db_dir,
                                       gene_col = "gene_reg", sample_col = "sample") {
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
                                                  gene_col = "gene_reg", sample_col = "sample") {
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
run_simple_ora_analysis_re <- function(te_data,
                                        group_col,
                                        groups,
                                        analysis_name,
                                        base_dir,  # ora/
                                        databases = DATABASES_TO_RUN,
                                        background_genes = BACKGROUND_GENES,
                                        gene_col = "gene_reg",
                                        sample_col = "sample") {

  cat("\n========================================\n")
  cat("  Simple ORA Analysis (RE):", analysis_name, "\n")
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
# Problem 7: Run GLM ONCE per minsample with p_threshold=1.0, then filter by qgene afterward
# Problem 6: Use combined database format for pathways_significant
run_glm_ora_analysis_re <- function(te_data,
                                     group_col,
                                     groups,
                                     analysis_name,
                                     base_dir,
                                     databases = DATABASES_TO_RUN,
                                     background_genes = BACKGROUND_GENES,
                                     gene_col = "gene_reg",
                                     sample_col = "sample",
                                     covariates = covar_med) {

  cat("\n========================================\n")
  cat("  GLM-ORA Analysis (RE):", analysis_name, "\n")
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

  # Create database subdirectories for plots
  for (db in databases) {
    db_plot_dir <- paste0(analysis_base_dir, db, "/")
    dir.create(db_plot_dir, showWarnings = FALSE, recursive = TRUE)
  }

  # Problem 7: Loop over minsample, run GLM ONCE with p_threshold=1.0
  for (minsample in MIN_SAMPLES_VALUES) {
    cat("\n--- minsample =", minsample, "---\n")

    # Need to use sample.x column for RE data
    te_data_renamed <- te_data
    if (sample_col != "sample_id" && sample_col %in% colnames(te_data_renamed)) {
      te_data_renamed$sample_id <- te_data_renamed[[sample_col]]
    }

    # Run GLM ONCE with p_threshold=1.0 to get ALL results
    glm_result <- filter_genes_by_differential_incidence_glm(
      te_data = te_data_renamed,
      group_col = group_col,
      groups = groups,
      min_samples = minsample,
      p_threshold = 1.0,  # Get ALL genes, filter by qgene later
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

    # Save full GLM results (before qgene filtering)
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

    # Now loop over qgene thresholds (Problem 7 - filter AFTER GLM)
    for (qgene in Q_GENE_VALUES) {
      cat("\n--- minsample =", minsample, ", qgene =", qgene, "---\n")

      # Problem 7: Filter the already-computed full_results by qgene threshold
      # Determine which column to use for filtering based on n_groups
      if (n_groups >= 3) {
        p_col <- "lrt_padj"
      } else {
        p_col <- "p_adj"
      }

      # Filter to significant genes
      if (!is.null(full_results) && nrow(full_results) > 0 && p_col %in% colnames(full_results)) {
        sig_results <- full_results %>% filter(.data[[p_col]] < qgene)
        sig_genes <- unique(sig_results$gene)
      } else {
        sig_results <- data.frame()
        sig_genes <- character(0)
      }

      # Save filtered GLM results for this qgene
      glm_results_file <- paste0(files_dir, "glm_", analysis_name,
                                 "_minsample", minsample, "_qgene", qgene, ".csv")
      if (nrow(sig_results) > 0) {
        write.csv(sig_results, glm_results_file, row.names = FALSE)
        cat("  Saved GLM results:", basename(glm_results_file), "(", nrow(sig_results), "genes )\n")
      } else {
        # Save empty file with headers for 3+ group (LRT) or 2-group case
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

      genes_significant <- length(sig_genes)
      ora_ran <- FALSE

      # Problem 6: Initialize pathway_counts_by_db for combined format
      pathway_counts_by_db <- list()
      for (db in databases) {
        pathway_counts_by_db[[db]] <- list()
        for (mingene in MIN_GENE_VALUES) {
          for (qpathway in Q_PATHWAY_VALUES) {
            key <- paste(mingene, qpathway, sep = "_")
            pathway_counts_by_db[[db]][[key]] <- 0
          }
        }
      }

      if (length(sig_genes) < 5) {
        cat("Too few significant genes (", length(sig_genes), "). Skipping ORA.\n")

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

      cat("Running ORA on", length(sig_genes), "significant genes\n")
      ora_ran <- TRUE

      # Always use compareCluster for all analyses
      cat("\n  === compareCluster ANALYSIS ===\n")

        # Step 1: Assign genes to groups using post-hoc pairwise contrasts
        if (!is.null(pairwise_results) && nrow(pairwise_results) > 0) {
          cat("  Using post-hoc pairwise contrasts for gene assignment\n")
          gene_assignment_df <- assign_genes_by_posthoc(pairwise_results, p_threshold = qgene)

          # Save post-hoc results
          posthoc_file <- paste0(files_dir, "glm_", analysis_name,
                                 "_posthoc_minsample", minsample, "_qgene", qgene, ".csv")
          write.csv(pairwise_results, posthoc_file, row.names = FALSE)
          cat("  Saved post-hoc pairwise results:", basename(posthoc_file), "\n")

          if (nrow(gene_assignment_df) > 0) {
            # Convert to clusters format for compareCluster
            gene_clusters <- gene_assignment_to_clusters(gene_assignment_df)

            # Save gene assignment
            gene_assign_file <- paste0(files_dir, "glm_", analysis_name,
                                       "_gene_assignment_minsample", minsample, "_qgene", qgene, ".csv")
            # Expand list column for saving
            gene_assign_expanded <- gene_assignment_df %>%
              mutate(assigned_groups_str = sapply(assigned_groups, paste, collapse = ";"))
            write.csv(gene_assign_expanded[, c("gene", "assigned_groups_str")], gene_assign_file, row.names = FALSE)
            cat("  Saved gene-group assignment:", basename(gene_assign_file), "(", nrow(gene_assignment_df), "genes)\n")

            # Step 2: Run compareCluster for each database
            for (db in databases) {
              cat("\n  Database:", db, "(compareCluster)\n")
              db_plot_dir <- paste0(analysis_base_dir, db, "/")

              tryCatch({
                # Run compareCluster based on database type
                if (db == "GO_BP") {
                  compare_result <- clusterProfiler::compareCluster(
                    geneCluster = gene_clusters,
                    fun = "enrichGO",
                    OrgDb = org.Hs.eg.db::org.Hs.eg.db,
                    keyType = "SYMBOL",
                    ont = "BP",
                    universe = background_genes,
                    pvalueCutoff = 1,
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
                    compare_result <- clusterProfiler::compareCluster(
                      geneCluster = gene_clusters_entrez,
                      fun = "enrichPathway",
                      organism = "human",
                      pvalueCutoff = 1,
                      qvalueCutoff = 1
                    )
                  } else {
                    compare_result <- NULL
                  }
                } else if (db %in% c("Hallmark", "Oncogenic")) {
                  # MSigDB databases
                  msig_category <- if (db == "Hallmark") "H" else "C6"
                  msig_db <- msigdbr::msigdbr(species = "Homo sapiens", category = msig_category)
                  msig_t2g <- msig_db %>% dplyr::select(gs_name, gene_symbol)

                  compare_result <- clusterProfiler::compareCluster(
                    geneCluster = gene_clusters,
                    fun = "enricher",
                    TERM2GENE = msig_t2g,
                    universe = background_genes,
                    pvalueCutoff = 1,
                    qvalueCutoff = 1
                  )
                } else {
                  compare_result <- NULL
                }

                if (!is.null(compare_result) && nrow(as.data.frame(compare_result)) > 0) {
                  # Save compareCluster result to files_dir (no database subdirectory)
                  compare_rds <- paste0(files_dir, "glm_ora_pathway_", db, "_", analysis_name,
                                        "_minsample", minsample, "_qgene", qgene, ".rds")
                  saveRDS(compare_result, compare_rds)

                  # Save as CSV
                  compare_csv <- paste0(files_dir, "glm_ora_pathway_", db, "_", analysis_name,
                                        "_minsample", minsample, "_qgene", qgene, ".csv")
                  write.csv(as.data.frame(compare_result), compare_csv, row.names = FALSE)

                  # Post-hoc filter by mingene and qpathway - collect counts in pathway_counts_by_db
                  for (mingene in MIN_GENE_VALUES) {
                    for (qpathway in Q_PATHWAY_VALUES) {
                      compare_df <- as.data.frame(compare_result)
                      filtered_df <- compare_df %>%
                        filter(Count >= mingene, qvalue < qpathway)

                      pathways_sig <- nrow(filtered_df)

                      # Problem 6: Collect count in pathway_counts_by_db instead of adding per-db row
                      key <- paste(mingene, qpathway, sep = "_")
                      pathway_counts_by_db[[db]][[key]] <- pathways_sig

                      if (pathways_sig == 0) next

                      # Create output prefix with glm_ora_pathway_ naming
                      prefix <- paste0("glm_ora_pathway_", db, "_", analysis_name,
                                       "_minsample", minsample, "_qgene", qgene,
                                       "_mingene", mingene, "_qpathway", qpathway)

                      # Save filtered CSV to db_plot_dir/files/
                      db_files_dir <- paste0(db_plot_dir, "files/")
                      dir.create(db_files_dir, showWarnings = FALSE, recursive = TRUE)
                      csv_file <- paste0(db_files_dir, prefix, ".csv")
                      write.csv(filtered_df, csv_file, row.names = FALSE)

                      # Create all ORA visualizations (7 types) with filtered data
                      create_ora_visualizations(compare_result, filtered_df, te_data, db, prefix, db_plot_dir,
                                                gene_col, actual_sample_col)

                      cat("    Saved:", prefix, "(", pathways_sig, "terms)\n")
                    }
                  }
                } else {
                  # No results - set counts to 0 in pathway_counts_by_db
                  for (mingene in MIN_GENE_VALUES) {
                    for (qpathway in Q_PATHWAY_VALUES) {
                      key <- paste(mingene, qpathway, sep = "_")
                      pathway_counts_by_db[[db]][[key]] <- 0
                    }
                  }
                  cat("    No enriched terms found\n")
                }
              }, error = function(e) {
                cat("    Error running compareCluster for", db, ":", e$message, "\n")
                # Set counts to 0 in pathway_counts_by_db on error
                for (mingene in MIN_GENE_VALUES) {
                  for (qpathway in Q_PATHWAY_VALUES) {
                    key <- paste(mingene, qpathway, sep = "_")
                    pathway_counts_by_db[[db]][[key]] <- 0
                  }
                }
              })
            }
          } else {
            cat("  No genes could be assigned to groups from post-hoc results\n")
            # pathway_counts_by_db already initialized to 0
          }
        } else {
          cat("  No post-hoc pairwise results available\n")
          # pathway_counts_by_db already initialized to 0
        }

        # Problem 6: Add combined run_summary entries for 3+ groups
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
              enrichment_ran = ora_ran,
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
run_glm_gsea_analysis_re <- function(te_data,
                                      group_col,
                                      groups,
                                      analysis_name,
                                      base_dir,
                                      databases = DATABASES_TO_RUN,
                                      gene_col = "gene_reg",
                                      sample_col = "sample",
                                      covariates = covar_med) {

  cat("\n========================================\n")
  cat("  GLM-GSEA Analysis (RE):", analysis_name, "\n")
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

    # Need to use sample.x column for RE data
    te_data_renamed <- te_data
    if (sample_col != "sample_id" && sample_col %in% colnames(te_data_renamed)) {
      te_data_renamed$sample_id <- te_data_renamed[[sample_col]]
    }

    # Run GLM ONCE with p_threshold=1.0 to get ALL results
    glm_result <- filter_genes_by_differential_incidence_glm(
      te_data = te_data_renamed,
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

                # Create GSEA visualizations - each in its own tryCatch
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

                    # Create GSEA visualizations - each in its own tryCatch
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

#### RE PATHWAY ANALYSIS 1: KICS VS HOSTSEQ ####
write_output(quote(NULL), "RE Pathway Analysis - KICS vs HostSeq")

cat("Dataset: te_kics_hostseq_re_split\n")
cat("Samples:", length(unique(te_kics_hostseq_re_split$sample)), "\n")
cohort_groups <- unique(te_kics_hostseq_re_split$cohort)
cohort_groups <- cohort_groups[!is.na(cohort_groups)]
cat("Cohort groups:", paste(cohort_groups, collapse = ", "), "\n\n")

if (length(cohort_groups) >= 2) {
  # Use covar_no_age because HostSeq samples lack age_at_diagnosis
  # Run all three analysis types: ora, glm_ora, glm_gsea
  run_simple_ora_analysis_re(
    te_data = te_kics_hostseq_re_split,
    group_col = "cohort",
    groups = cohort_groups,
    analysis_name = "kics_hostseq",
    base_dir = ora_dir
  )
  run_glm_ora_analysis_re(
    te_data = te_kics_hostseq_re_split,
    group_col = "cohort",
    groups = cohort_groups,
    analysis_name = "kics_hostseq",
    base_dir = glm_ora_dir,
    covariates = covar_no_age
  )
  run_glm_gsea_analysis_re(
    te_data = te_kics_hostseq_re_split,
    group_col = "cohort",
    groups = cohort_groups,
    analysis_name = "kics_hostseq",
    base_dir = glm_gsea_dir,
    covariates = covar_no_age
  )
} else {
  cat("Not enough cohort groups. Skipping.\n")
}

#### RE PATHWAY ANALYSIS 2: TP53 STATUS ####
write_output(quote(NULL), "RE Pathway Analysis - TP53 Status")

cat("Dataset: te_aff_re_split\n")
cat("Samples:", length(unique(te_aff_re_split$sample)), "\n")
tp53_groups <- unique(te_aff_re_split$TP53_status)
tp53_groups <- tp53_groups[!is.na(tp53_groups)]
cat("TP53 groups:", paste(tp53_groups, collapse = ", "), "\n\n")

if (length(tp53_groups) >= 2) {
  run_simple_ora_analysis_re(
    te_data = te_aff_re_split,
    group_col = "TP53_status",
    groups = tp53_groups,
    analysis_name = "tp53_status",
    base_dir = ora_dir
  )
  run_glm_ora_analysis_re(
    te_data = te_aff_re_split,
    group_col = "TP53_status",
    groups = tp53_groups,
    analysis_name = "tp53_status",
    base_dir = glm_ora_dir
  )
  run_glm_gsea_analysis_re(
    te_data = te_aff_re_split,
    group_col = "TP53_status",
    groups = tp53_groups,
    analysis_name = "tp53_status",
    base_dir = glm_gsea_dir
  )
} else {
  cat("Not enough TP53 groups. Skipping.\n")
}

#### RE PATHWAY ANALYSIS 3: TP53 3-LEVEL ####
write_output(quote(NULL), "RE Pathway Analysis - TP53 3-Level")

cat("Dataset: te_aff_re_split (TP53_3level)\n")

# Filter to samples with valid TP53_3level
te_aff_re_3level <- te_aff_re_split %>% filter(!is.na(TP53_3level))
if (nrow(te_aff_re_3level) > 0) {
  tp53_3level_groups <- unique(te_aff_re_3level$TP53_3level)
  tp53_3level_groups <- tp53_3level_groups[!is.na(tp53_3level_groups)]
  cat("TP53_3level groups:", paste(tp53_3level_groups, collapse = ", "), "\n")

  if (length(tp53_3level_groups) >= 2) {
    run_simple_ora_analysis_re(
      te_data = te_aff_re_3level,
      group_col = "TP53_3level",
      groups = tp53_3level_groups,
      analysis_name = "tp53_3level",
      base_dir = ora_dir
    )
    run_glm_ora_analysis_re(
      te_data = te_aff_re_3level,
      group_col = "TP53_3level",
      groups = tp53_3level_groups,
      analysis_name = "tp53_3level",
      base_dir = glm_ora_dir
    )
    run_glm_gsea_analysis_re(
      te_data = te_aff_re_3level,
      group_col = "TP53_3level",
      groups = tp53_3level_groups,
      analysis_name = "tp53_3level",
      base_dir = glm_gsea_dir
    )
  } else {
    cat("Not enough TP53_3level groups. Skipping.\n")
  }
}

#### RE PATHWAY ANALYSIS 4: CANCER STATUS (LFS) ####
write_output(quote(NULL), "RE Pathway Analysis - Cancer Status LFS")

if ("Cancer" %in% colnames(te_lfs_re_split)) {
  cat("Dataset: te_lfs_re_split\n")
  cat("Samples:", length(unique(te_lfs_re_split$sample)), "\n")
  cancer_groups <- unique(te_lfs_re_split$Cancer)
  cancer_groups <- cancer_groups[!is.na(cancer_groups)]
  cat("Cancer groups:", paste(cancer_groups, collapse = ", "), "\n\n")

  if (length(cancer_groups) >= 2) {
    # Use covar_no_age because unaffected LFS samples lack age_at_diagnosis
    run_simple_ora_analysis_re(
      te_data = te_lfs_re_split,
      group_col = "Cancer",
      groups = cancer_groups,
      analysis_name = "cancer_lfs",
      base_dir = ora_dir
    )
    run_glm_ora_analysis_re(
      te_data = te_lfs_re_split,
      group_col = "Cancer",
      groups = cancer_groups,
      analysis_name = "cancer_lfs",
      base_dir = glm_ora_dir,
      covariates = covar_no_age
    )
    run_glm_gsea_analysis_re(
      te_data = te_lfs_re_split,
      group_col = "Cancer",
      groups = cancer_groups,
      analysis_name = "cancer_lfs",
      base_dir = glm_gsea_dir,
      covariates = covar_no_age
    )
  } else {
    cat("Not enough Cancer groups. Skipping.\n")
  }
}

#### RE PATHWAY ANALYSIS 5: LFS BY TP53 STATUS ####
write_output(quote(NULL), "RE Pathway Analysis - LFS by TP53 Status")

if (exists("te_lfs_re_split") && "TP53_status" %in% colnames(te_lfs_re_split)) {
  lfs_tp53_groups <- unique(te_lfs_re_split$TP53_status)
  lfs_tp53_groups <- lfs_tp53_groups[!is.na(lfs_tp53_groups)]

  if (length(lfs_tp53_groups) >= 2) {
    cat("Dataset: te_lfs_re_split by TP53_status\n")
    cat("Samples:", length(unique(te_lfs_re_split$sample)), "\n")
    cat("TP53 groups:", paste(lfs_tp53_groups, collapse = ", "), "\n")

    run_simple_ora_analysis_re(
      te_data = te_lfs_re_split,
      group_col = "TP53_status",
      groups = lfs_tp53_groups,
      analysis_name = "lfs_tp53_status",
      base_dir = ora_dir
    )
    run_glm_ora_analysis_re(
      te_data = te_lfs_re_split,
      group_col = "TP53_status",
      groups = lfs_tp53_groups,
      analysis_name = "lfs_tp53_status",
      base_dir = glm_ora_dir
    )
    run_glm_gsea_analysis_re(
      te_data = te_lfs_re_split,
      group_col = "TP53_status",
      groups = lfs_tp53_groups,
      analysis_name = "lfs_tp53_status",
      base_dir = glm_gsea_dir
    )
  } else {
    cat("Not enough TP53 groups in LFS data. Skipping.\n")
  }
}

#### RE PATHWAY ANALYSIS 6: KICS TUMOR TYPE ####
write_output(quote(NULL), "RE Pathway Analysis - KICS by Tumor Type")

if ("tumor_type" %in% colnames(te_kics_re_split)) {
  # Filter to tumor types with >= 3 samples
  tumor_type_counts <- table(te_kics_re_split$tumor_type)
  valid_tumor_types <- names(tumor_type_counts[tumor_type_counts >= 3])
  cat("Tumor types with >= 3 samples:", paste(valid_tumor_types, collapse = ", "), "\n")

  if (length(valid_tumor_types) >= 2) {
    te_kics_re_tumor <- te_kics_re_split %>%
      filter(tumor_type %in% valid_tumor_types)

    run_simple_ora_analysis_re(
      te_data = te_kics_re_tumor,
      group_col = "tumor_type",
      groups = valid_tumor_types,
      analysis_name = "kics_tumor_type",
      base_dir = ora_dir
    )
    run_glm_ora_analysis_re(
      te_data = te_kics_re_tumor,
      group_col = "tumor_type",
      groups = valid_tumor_types,
      analysis_name = "kics_tumor_type",
      base_dir = glm_ora_dir,
      covariates = covar_no_tumour_type
    )
    run_glm_gsea_analysis_re(
      te_data = te_kics_re_tumor,
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

#### RE PATHWAY ANALYSIS 7: KICS BY ANCESTRY ####
write_output(quote(NULL), "RE Pathway Analysis - KICS by Ancestry")

if ("predicted_ancestry_thres" %in% colnames(te_kics_re_split)) {
  ancestries <- unique(te_kics_re_split$predicted_ancestry_thres)
  ancestries <- ancestries[!is.na(ancestries)]
  cat("Ancestries:", paste(ancestries, collapse = ", "), "\n")

  # Only run if we have enough groups with samples
  ancestry_counts <- table(te_kics_re_split$predicted_ancestry_thres)
  valid_ancestries <- names(ancestry_counts[ancestry_counts >= 10])

  if (length(valid_ancestries) >= 2) {
    te_kics_re_ancestry <- te_kics_re_split %>%
      filter(predicted_ancestry_thres %in% valid_ancestries)

    run_simple_ora_analysis_re(
      te_data = te_kics_re_ancestry,
      group_col = "predicted_ancestry_thres",
      groups = valid_ancestries,
      analysis_name = "kics_ancestry",
      base_dir = ora_dir
    )
    run_glm_ora_analysis_re(
      te_data = te_kics_re_ancestry,
      group_col = "predicted_ancestry_thres",
      groups = valid_ancestries,
      analysis_name = "kics_ancestry",
      base_dir = glm_ora_dir,
      covariates = covar_no_ancestry
    )
    run_glm_gsea_analysis_re(
      te_data = te_kics_re_ancestry,
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

#### RE PATHWAY ANALYSIS 8: LFS BY COHORT ####
write_output(quote(NULL), "RE Pathway Analysis - LFS by Cohort")

if ("cohort" %in% colnames(te_lfs_re_split)) {
  cohorts <- unique(te_lfs_re_split$cohort)
  cohorts <- cohorts[!is.na(cohorts)]
  cat("Cohorts:", paste(cohorts, collapse = ", "), "\n")

  if (length(cohorts) >= 2) {
    run_simple_ora_analysis_re(
      te_data = te_lfs_re_split,
      group_col = "cohort",
      groups = cohorts,
      analysis_name = "lfs_cohort",
      base_dir = ora_dir
    )
    run_glm_ora_analysis_re(
      te_data = te_lfs_re_split,
      group_col = "cohort",
      groups = cohorts,
      analysis_name = "lfs_cohort",
      base_dir = glm_ora_dir
    )
    run_glm_gsea_analysis_re(
      te_data = te_lfs_re_split,
      group_col = "cohort",
      groups = cohorts,
      analysis_name = "lfs_cohort",
      base_dir = glm_gsea_dir
    )
  }
}

#### RE PATHWAY ANALYSIS 9: KICS BY SAMPLE TYPE ####
write_output(quote(NULL), "RE Pathway Analysis - KICS by Sample Type")

# Load and merge sample type data for KICS
sample_type_file <- "/Users/briannelaverty/Documents/R_Malkin/clinical/kics_germline_sample_type.csv"
if (file.exists(sample_type_file)) {
  kics_sample_type <- prep_kics_sample_type(sample_type_file)
  cat("Loaded sample type data:", nrow(kics_sample_type), "samples\n")

  # Merge with KICS RE data
  te_kics_re_sampletype <- merge(te_kics_re_split, kics_sample_type,
                                  by = "sample", all.x = TRUE)
  cat("Merged sample type with RE data:", nrow(te_kics_re_sampletype), "rows\n")
  cat("Samples with sample_type:", sum(!is.na(te_kics_re_sampletype$sample_type)), "\n")

  # Filter to valid sample types
  valid_sample_types <- c("Blood", "Fibroblasts", "Tissue (fresh)")
  te_kics_re_sampletype <- te_kics_re_sampletype %>%
    filter(sample_type %in% valid_sample_types)
  cat("After filtering to valid types:", nrow(te_kics_re_sampletype), "rows\n")

  sample_types <- unique(te_kics_re_sampletype$sample_type)
  sample_types <- sample_types[!is.na(sample_types)]
  cat("Sample types:", paste(sample_types, collapse = ", "), "\n")

  if (length(sample_types) >= 2) {
    run_simple_ora_analysis_re(
      te_data = te_kics_re_sampletype,
      group_col = "sample_type",
      groups = sample_types,
      analysis_name = "kics_sample_type",
      base_dir = ora_dir
    )
    run_glm_ora_analysis_re(
      te_data = te_kics_re_sampletype,
      group_col = "sample_type",
      groups = sample_types,
      analysis_name = "kics_sample_type",
      base_dir = glm_ora_dir
    )
    run_glm_gsea_analysis_re(
      te_data = te_kics_re_sampletype,
      group_col = "sample_type",
      groups = sample_types,
      analysis_name = "kics_sample_type",
      base_dir = glm_gsea_dir
    )
  } else {
    cat("Not enough sample types. Skipping.\n")
  }
} else {
  cat("Sample type file not found:", sample_type_file, "\n")
}

#### RE PATHWAY ANALYSIS 10: TAYLOR BY SUBTYPE ####
write_output(quote(NULL), "RE Pathway Analysis - Taylor by Subtype")

if (exists("te_taylor_re_split") && nrow(te_taylor_re_split) > 0) {
  if ("tumor_type_subclass" %in% colnames(te_taylor_re_split)) {
    cat("Dataset: te_taylor_re_split\n")

    # Filter to subtypes with >= 3 samples
    subtype_counts <- table(te_taylor_re_split$tumor_type_subclass)
    valid_subtypes <- names(subtype_counts[subtype_counts >= 3])
    cat("Subtypes with >= 3 samples:", paste(valid_subtypes, collapse = ", "), "\n")

    if (length(valid_subtypes) >= 2) {
      te_taylor_re_subtype <- te_taylor_re_split %>%
        filter(tumor_type_subclass %in% valid_subtypes)

      run_simple_ora_analysis_re(
        te_data = te_taylor_re_subtype,
        group_col = "tumor_type_subclass",
        groups = valid_subtypes,
        analysis_name = "taylor_subtype",
        base_dir = ora_dir
      )
      run_glm_ora_analysis_re(
        te_data = te_taylor_re_subtype,
        group_col = "tumor_type_subclass",
        groups = valid_subtypes,
        analysis_name = "taylor_subtype",
        base_dir = glm_ora_dir
      )
      run_glm_gsea_analysis_re(
        te_data = te_taylor_re_subtype,
        group_col = "tumor_type_subclass",
        groups = valid_subtypes,
        analysis_name = "taylor_subtype",
        base_dir = glm_gsea_dir
      )
    } else {
      cat("Not enough subtypes with >= 3 samples. Skipping.\n")
    }
  } else {
    cat("tumor_type_subclass column not found. Skipping.\n")
  }
} else {
  cat("Taylor RE data not available. Skipping.\n")
}

#### SAVE RE DATA FOR RE-RNA SCRIPT ####
cat("\n===== SAVING RE DATA =====\n")

# Save RE split data for use by RE-RNA script
re_files_dir <- paste0(re_dir, "files/")
dir.create(re_files_dir, showWarnings = FALSE, recursive = TRUE)

saveRDS(te_aff_re_split, paste0(re_files_dir, "te_aff_re_split.rds"))
saveRDS(te_lfs_re_split, paste0(re_files_dir, "te_lfs_re_split.rds"))
saveRDS(te_kics_hostseq_re_split, paste0(re_files_dir, "te_kics_hostseq_re_split.rds"))
saveRDS(te_kics_re_split, paste0(re_files_dir, "te_kics_re_split.rds"))
if (exists("te_taylor_re_split") && nrow(te_taylor_re_split) > 0) {
  saveRDS(te_taylor_re_split, paste0(re_files_dir, "te_taylor_re_split.rds"))
}
cat("Saved RE data for RE-RNA script\n")

#### SUMMARY ####
write_output(quote(NULL), "RE Pathway Analysis Summary")

cat("Analysis completed for the following datasets:\n")
cat("  1. kics_hostseq - KICS vs HostSeq\n")
cat("  2. tp53_status - TP53 status (affected samples)\n")
cat("  3. tp53_3level - TP53 3-level (Germline/Somatic/WT)\n")
cat("  4. cancer_lfs - Cancer status LFS\n")
cat("  5. lfs_tp53_status - LFS by TP53 Status\n")
cat("  6. kics_tumor_type - KICS by tumor type\n")
cat("  7. kics_ancestry - KICS by ancestry\n")
cat("  8. lfs_cohort - LFS by cohort\n")
cat("  9. kics_sample_type - KICS by sample type\n")
cat("  10. taylor_subtype - Taylor by subtype\n\n")

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
cat("Output location:", re_dir, "\n")

# Sort run summary by analysis_name first (to group all same analyses together),
# then by analysis_type, then by minsample, qgene, mingene, qpathway
run_summary <- run_summary %>%
  arrange(analysis_name, analysis_type, minsample, qgene, mingene, qpathway)

# Save run summary
run_summary_file <- paste0(re_dir, "RUN_SUMMARY.csv")
write.csv(run_summary, run_summary_file, row.names = FALSE)
cat("\nSaved run summary:", run_summary_file, "\n")
cat("Run summary sorted by: analysis_name, analysis_type, minsample, qgene, mingene, qpathway\n")

cat("\nNote: RE-RNA differential expression analysis is in 02_te_viz_germline_08_re_rna.R\n")

cat("\n Script completed successfully\n")

# Close module-specific sink
close_module_sink()
