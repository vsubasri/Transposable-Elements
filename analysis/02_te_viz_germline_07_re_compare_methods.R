#!/usr/bin/env Rscript

# Germline RE Analysis - GLM vs Fisher Method Comparison
# Purpose: Compare differential gene incidence using:
#   1. GLM with covariates (existing method)
#   2. Simple Fisher's exact test (no covariates)
# Analysis: KICS vs HostSeq only
# Output: Comparison of genes and pathways found by each method

# Source common setup and load data
source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/00_viz_common_setup.R")
REQUIRED_DATA <- c("count_matrix", "expand", "split", "clinical")
source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/00_viz_load_data_germline.R")

library(emmeans)

cat("=============================================================\n")
cat("  RE ANALYSIS: GLM vs FISHER METHOD COMPARISON\n")
cat("  KICS vs HostSeq Analysis Only\n")
cat("=============================================================\n\n")

#### CONFIGURATION ####

DATABASES_TO_RUN <- c("GO_BP", "Reactome", "Hallmark", "Oncogenic")

# Parameter sweep for comprehensive comparison
MIN_SAMPLES_VALUES <- c(1, 3, 5)
Q_GENE_VALUES <- c(0.1, 0.25)
MIN_GENE <- 3
Q_PATHWAY <- 0.25

# Load background genes
background_genes_file <- "/Users/briannelaverty/Documents/R_Malkin/TE/data/final/background_regulatory_genes.txt"
if (file.exists(background_genes_file)) {
  BACKGROUND_GENES <- readLines(background_genes_file)
  cat("Loaded", length(BACKGROUND_GENES), "background regulatory genes\n")
} else {
  BACKGROUND_GENES <- NULL
  cat("Warning: Background genes file not found\n")
}

# Create output directory
compare_dir <- paste0(plot_dir, "reg_element/method_comparison/")
dir.create(compare_dir, showWarnings = FALSE, recursive = TRUE)
files_dir <- paste0(compare_dir, "files/")
dir.create(files_dir, showWarnings = FALSE, recursive = TRUE)

cat("Output directory:", compare_dir, "\n\n")

#### LOAD RE DATA ####
cat("Loading Regulatory Elements Data...\n")

re_germline_path <- "/Users/briannelaverty/Documents/R_Malkin/TE/data/final/germline_annotSV_output.SV_RE_intersect.report"

fix_sample_column <- function(df) {
  if ("sample.x" %in% colnames(df) && !"sample" %in% colnames(df)) {
    df <- df %>% rename(sample = sample.x)
  }
  if ("sample.y" %in% colnames(df)) {
    df <- df %>% select(-sample.y)
  }
  return(df)
}

# Load KICS + HostSeq RE data
te_kics_hostseq_re <- load_and_join_re_data(te_kics_hostseq_expand, re_germline_path)
te_kics_hostseq_re_split <- split_re_genes(te_kics_hostseq_re) %>% fix_sample_column()
cat("KICS+HostSeq RE data:", nrow(te_kics_hostseq_re_split), "rows\n")
cat("Unique genes:", length(unique(te_kics_hostseq_re_split$gene_reg)), "\n")

# Get sample column name
sample_col <- if ("sample" %in% colnames(te_kics_hostseq_re_split)) "sample" else "sample.x"
cat("Sample column:", sample_col, "\n")

# Cohort info
cohort_groups <- unique(te_kics_hostseq_re_split$cohort)
cohort_groups <- cohort_groups[!is.na(cohort_groups)]
cat("Cohort groups:", paste(cohort_groups, collapse = ", "), "\n")

# Sample counts
sample_counts <- te_kics_hostseq_re_split %>%
  group_by(cohort) %>%
  summarise(n = n_distinct(.data[[sample_col]]), .groups = "drop")
cat("\nSample counts:\n")
print(sample_counts)
cat("\n")

#### WILCOXON TEST FUNCTION ####
# Tests whether TE count per gene differs between groups using Wilcoxon rank-sum test
filter_genes_by_wilcox <- function(te_data,
                                   group_col,
                                   groups,
                                   min_samples = 5,
                                   gene_col = "gene_reg",
                                   sample_col = "sample") {

  cat("\n===== FILTERING GENES BY WILCOXON RANK-SUM TEST =====\n")
  cat("Groups:", paste(groups, collapse = ", "), "\n")
  cat("Min samples:", min_samples, "\n")

  if (!sample_col %in% colnames(te_data)) {
    sample_col <- if ("sample" %in% colnames(te_data)) "sample" else "sample.x"
  }

  # Get all samples by group
  samples_by_group <- lapply(groups, function(g) {
    unique(te_data[[sample_col]][te_data[[group_col]] == g & !is.na(te_data[[group_col]])])
  })
  names(samples_by_group) <- groups

  for (g in groups) {
    cat("Group", g, ":", length(samples_by_group[[g]]), "samples\n")
  }

  all_samples <- unlist(samples_by_group)

  # Get genes meeting min_samples threshold
  gene_counts <- te_data %>%
    filter(.data[[sample_col]] %in% all_samples) %>%
    filter(!is.na(.data[[gene_col]]) & .data[[gene_col]] != "") %>%
    group_by(.data[[gene_col]]) %>%
    summarise(n_samples = n_distinct(.data[[sample_col]]), .groups = "drop") %>%
    filter(n_samples >= min_samples)

  cat("Genes passing min_samples filter:", nrow(gene_counts), "\n")

  if (nrow(gene_counts) == 0) {
    return(list(genes = character(0), full_results = NULL, genes_tested = 0))
  }

  # Create sample-level count matrix for each gene
  # Count TEs per sample per gene
  te_counts_per_sample <- te_data %>%
    filter(.data[[sample_col]] %in% all_samples) %>%
    filter(.data[[gene_col]] %in% gene_counts[[gene_col]]) %>%
    group_by(.data[[gene_col]], .data[[sample_col]]) %>%
    summarise(te_count = n(), .groups = "drop")

  # Create complete matrix with zeros for samples without TEs
  all_combinations <- expand.grid(
    gene = unique(gene_counts[[gene_col]]),
    sample = all_samples,
    stringsAsFactors = FALSE
  )
  colnames(all_combinations) <- c(gene_col, sample_col)

  te_counts_complete <- all_combinations %>%
    left_join(te_counts_per_sample, by = c(gene_col, sample_col)) %>%
    mutate(te_count = ifelse(is.na(te_count), 0, te_count))

  # Add group info
  sample_groups <- data.frame(
    sample = all_samples,
    group = sapply(all_samples, function(s) {
      for (g in groups) {
        if (s %in% samples_by_group[[g]]) return(g)
      }
      return(NA)
    }),
    stringsAsFactors = FALSE
  )
  colnames(sample_groups) <- c(sample_col, "group")

  te_counts_complete <- te_counts_complete %>%
    left_join(sample_groups, by = sample_col)

  results_list <- list()
  genes_to_test <- unique(gene_counts[[gene_col]])
  cat("Testing", length(genes_to_test), "genes with Wilcoxon rank-sum test...\n")

  for (gene in genes_to_test) {
    gene_data <- te_counts_complete %>%
      filter(.data[[gene_col]] == gene)

    group1_counts <- gene_data$te_count[gene_data$group == groups[1]]
    group2_counts <- gene_data$te_count[gene_data$group == groups[2]]

    # Wilcoxon rank-sum test
    wilcox_result <- tryCatch({
      wilcox.test(group1_counts, group2_counts, exact = FALSE)
    }, error = function(e) {
      list(p.value = NA, statistic = NA)
    })

    # Calculate summary stats
    mean_group1 <- mean(group1_counts)
    mean_group2 <- mean(group2_counts)
    median_group1 <- median(group1_counts)
    median_group2 <- median(group2_counts)

    # Effect direction based on median
    effect_direction <- if (median_group1 > median_group2) groups[1] else groups[2]

    # Samples with TE (count > 0)
    samples_with <- gene_data[[sample_col]][gene_data$te_count > 0]

    results_list[[length(results_list) + 1]] <- data.frame(
      gene = gene,
      samples_with_te = length(samples_with),
      samples_with_te_names = paste(samples_with, collapse = ";"),
      mean_group1 = mean_group1,
      mean_group2 = mean_group2,
      median_group1 = median_group1,
      median_group2 = median_group2,
      effect_direction = effect_direction,
      W_statistic = as.numeric(wilcox_result$statistic),
      p_value = wilcox_result$p.value,
      stringsAsFactors = FALSE
    )
  }

  if (length(results_list) == 0) {
    return(list(genes = character(0), full_results = NULL, genes_tested = length(genes_to_test)))
  }

  results_df <- bind_rows(results_list)
  results_df <- results_df %>%
    mutate(p_adj = p.adjust(p_value, method = "BH"))

  # Rename columns to include group names
  colnames(results_df)[colnames(results_df) == "mean_group1"] <- paste0("mean_", groups[1])
  colnames(results_df)[colnames(results_df) == "mean_group2"] <- paste0("mean_", groups[2])
  colnames(results_df)[colnames(results_df) == "median_group1"] <- paste0("median_", groups[1])
  colnames(results_df)[colnames(results_df) == "median_group2"] <- paste0("median_", groups[2])

  cat("Genes tested:", nrow(results_df), "\n")
  cat("=============================================\n\n")

  return(list(full_results = results_df, genes_tested = length(genes_to_test)))
}

#### RUN ORA FUNCTION ####
run_ora_for_genes <- function(sig_genes, method_name, output_prefix, databases = DATABASES_TO_RUN) {
  if (length(sig_genes) < 3) {
    cat(method_name, ": Too few genes (", length(sig_genes), ") for ORA\n")
    return(list())
  }

  cat(method_name, ": Running ORA on", length(sig_genes), "genes\n")
  pathway_results <- list()

  for (db in databases) {
    cat("  Database:", db, "... ")

    tryCatch({
      if (db == "GO_BP") {
        ora_result <- clusterProfiler::enrichGO(
          gene = sig_genes,
          OrgDb = org.Hs.eg.db::org.Hs.eg.db,
          keyType = "SYMBOL",
          ont = "BP",
          universe = BACKGROUND_GENES,
          pvalueCutoff = 1,
          qvalueCutoff = 1
        )
      } else if (db == "Reactome") {
        entrez_ids <- AnnotationDbi::mapIds(
          org.Hs.eg.db::org.Hs.eg.db,
          sig_genes, "ENTREZID", "SYMBOL"
        )
        entrez_ids <- entrez_ids[!is.na(entrez_ids)]

        if (length(entrez_ids) >= 3) {
          ora_result <- ReactomePA::enrichPathway(
            gene = entrez_ids,
            organism = "human",
            pvalueCutoff = 1,
            qvalueCutoff = 1
          )
        } else {
          ora_result <- NULL
        }
      } else if (db %in% c("Hallmark", "Oncogenic")) {
        msig_category <- if (db == "Hallmark") "H" else "C6"
        msig_db <- msigdbr::msigdbr(species = "Homo sapiens", category = msig_category)
        msig_t2g <- msig_db %>% dplyr::select(gs_name, gene_symbol)

        ora_result <- clusterProfiler::enricher(
          gene = sig_genes,
          TERM2GENE = msig_t2g,
          universe = BACKGROUND_GENES,
          pvalueCutoff = 1,
          qvalueCutoff = 1
        )
      } else {
        ora_result <- NULL
      }

      if (!is.null(ora_result) && nrow(as.data.frame(ora_result)) > 0) {
        ora_df <- as.data.frame(ora_result)
        ora_sig <- ora_df %>% filter(Count >= MIN_GENE, qvalue < Q_PATHWAY)

        cat(nrow(ora_sig), "significant pathways\n")

        if (nrow(ora_sig) > 0) {
          pathway_results[[db]] <- ora_sig
          write.csv(ora_sig, paste0(files_dir, output_prefix, "_", db, ".csv"), row.names = FALSE)
        }
      } else {
        cat("0 pathways\n")
      }

    }, error = function(e) {
      cat("error:", e$message, "\n")
    })
  }

  return(pathway_results)
}

#### MASTER COMPARISON LOOP ####

# Prepare data for GLM
te_data_glm <- te_kics_hostseq_re_split
if ("sample.x" %in% colnames(te_data_glm) && !"sample_id" %in% colnames(te_data_glm)) {
  te_data_glm$sample_id <- te_data_glm[[sample_col]]
}

# Master results table
master_results <- data.frame(
  min_sample = integer(),
  q_gene = numeric(),
  method = character(),
  genes_tested = integer(),
  genes_significant = integer(),
  pathways_go_bp = integer(),
  pathways_reactome = integer(),
  pathways_hallmark = integer(),
  pathways_oncogenic = integer(),
  pathways_total = integer(),
  stringsAsFactors = FALSE
)

for (min_sample in MIN_SAMPLES_VALUES) {
  for (q_gene in Q_GENE_VALUES) {

    cat("\n")
    cat("##########################################################\n")
    cat("#  PARAMETERS: min_sample =", min_sample, ", q_gene =", q_gene, "\n")
    cat("##########################################################\n")

    prefix <- paste0("minsample", min_sample, "_qgene", q_gene)

    # ===== GLM METHOD =====
    cat("\n--- METHOD 1: GLM WITH COVARIATES ---\n")

    glm_result <- filter_genes_by_differential_incidence_glm(
      te_data = te_data_glm,
      group_col = "cohort",
      groups = cohort_groups,
      min_samples = min_sample,
      p_threshold = 1.0,
      gene_col = "gene_reg",
      covariates = covar_no_age
    )

    glm_full <- glm_result$full_results
    glm_genes_tested <- glm_result$genes_attempted

    if (!is.null(glm_full) && nrow(glm_full) > 0) {
      glm_sig <- glm_full %>% filter(p_adj < q_gene)
      glm_sig_genes <- unique(glm_sig$gene)
      cat("GLM significant genes (q <", q_gene, "):", length(glm_sig_genes), "\n")

      # Save results
      write.csv(glm_full, paste0(files_dir, "glm_", prefix, "_all.csv"), row.names = FALSE)
      if (nrow(glm_sig) > 0) {
        write.csv(glm_sig, paste0(files_dir, "glm_", prefix, "_significant.csv"), row.names = FALSE)
      }
    } else {
      glm_sig <- data.frame()
      glm_sig_genes <- character(0)
    }

    # Run ORA for GLM
    glm_pathways <- run_ora_for_genes(glm_sig_genes, "GLM", paste0("ora_glm_", prefix))

    # Add to master results
    master_results <- rbind(master_results, data.frame(
      min_sample = min_sample,
      q_gene = q_gene,
      method = "GLM",
      genes_tested = glm_genes_tested,
      genes_significant = length(glm_sig_genes),
      pathways_go_bp = if ("GO_BP" %in% names(glm_pathways)) nrow(glm_pathways$GO_BP) else 0,
      pathways_reactome = if ("Reactome" %in% names(glm_pathways)) nrow(glm_pathways$Reactome) else 0,
      pathways_hallmark = if ("Hallmark" %in% names(glm_pathways)) nrow(glm_pathways$Hallmark) else 0,
      pathways_oncogenic = if ("Oncogenic" %in% names(glm_pathways)) nrow(glm_pathways$Oncogenic) else 0,
      pathways_total = if (length(glm_pathways) > 0) sum(sapply(glm_pathways, nrow)) else 0,
      stringsAsFactors = FALSE
    ))

    # ===== WILCOXON METHOD =====
    cat("\n--- METHOD 2: WILCOXON RANK-SUM TEST ---\n")

    wilcox_result <- filter_genes_by_wilcox(
      te_data = te_kics_hostseq_re_split,
      group_col = "cohort",
      groups = cohort_groups,
      min_samples = min_sample,
      gene_col = "gene_reg",
      sample_col = sample_col
    )

    wilcox_full <- wilcox_result$full_results
    wilcox_genes_tested <- wilcox_result$genes_tested

    if (!is.null(wilcox_full) && nrow(wilcox_full) > 0) {
      wilcox_sig <- wilcox_full %>% filter(p_adj < q_gene)
      wilcox_sig_genes <- unique(wilcox_sig$gene)
      cat("Wilcoxon significant genes (q <", q_gene, "):", length(wilcox_sig_genes), "\n")

      # Save results
      write.csv(wilcox_full, paste0(files_dir, "wilcox_", prefix, "_all.csv"), row.names = FALSE)
      if (nrow(wilcox_sig) > 0) {
        write.csv(wilcox_sig, paste0(files_dir, "wilcox_", prefix, "_significant.csv"), row.names = FALSE)
      }
    } else {
      wilcox_sig <- data.frame()
      wilcox_sig_genes <- character(0)
    }

    # Run ORA for Wilcoxon
    wilcox_pathways <- run_ora_for_genes(wilcox_sig_genes, "Wilcoxon", paste0("ora_wilcox_", prefix))

    # Add to master results
    master_results <- rbind(master_results, data.frame(
      min_sample = min_sample,
      q_gene = q_gene,
      method = "Wilcoxon",
      genes_tested = wilcox_genes_tested,
      genes_significant = length(wilcox_sig_genes),
      pathways_go_bp = if ("GO_BP" %in% names(wilcox_pathways)) nrow(wilcox_pathways$GO_BP) else 0,
      pathways_reactome = if ("Reactome" %in% names(wilcox_pathways)) nrow(wilcox_pathways$Reactome) else 0,
      pathways_hallmark = if ("Hallmark" %in% names(wilcox_pathways)) nrow(wilcox_pathways$Hallmark) else 0,
      pathways_oncogenic = if ("Oncogenic" %in% names(wilcox_pathways)) nrow(wilcox_pathways$Oncogenic) else 0,
      pathways_total = if (length(wilcox_pathways) > 0) sum(sapply(wilcox_pathways, nrow)) else 0,
      stringsAsFactors = FALSE
    ))

    # ===== COMPARE =====
    cat("\n--- COMPARISON ---\n")
    genes_both <- intersect(glm_sig_genes, wilcox_sig_genes)
    genes_glm_only <- setdiff(glm_sig_genes, wilcox_sig_genes)
    genes_wilcox_only <- setdiff(wilcox_sig_genes, glm_sig_genes)

    cat("Overlap:", length(genes_both), "genes\n")
    cat("GLM only:", length(genes_glm_only), "genes\n")
    cat("Wilcoxon only:", length(genes_wilcox_only), "genes\n")

    # Save gene comparison
    union_genes <- union(glm_sig_genes, wilcox_sig_genes)
    if (length(union_genes) > 0) {
      gene_comparison <- data.frame(
        gene = union_genes,
        in_glm = union_genes %in% glm_sig_genes,
        in_wilcox = union_genes %in% wilcox_sig_genes,
        in_both = union_genes %in% genes_both,
        stringsAsFactors = FALSE
      )
      write.csv(gene_comparison, paste0(files_dir, "gene_comparison_", prefix, ".csv"), row.names = FALSE)
    }
  }
}

#### FINAL SUMMARY ####
cat("\n")
cat("=============================================================\n")
cat("                    MASTER SUMMARY                           \n")
cat("=============================================================\n\n")

print(master_results)

# Calculate method differences
cat("\n--- METHOD COMPARISON ACROSS PARAMETERS ---\n")

for (min_sample in MIN_SAMPLES_VALUES) {
  for (q_gene in Q_GENE_VALUES) {
    glm_row <- master_results[master_results$min_sample == min_sample &
                               master_results$q_gene == q_gene &
                               master_results$method == "GLM", ]
    wilcox_row <- master_results[master_results$min_sample == min_sample &
                                  master_results$q_gene == q_gene &
                                  master_results$method == "Wilcoxon", ]

    cat("\nmin_sample =", min_sample, ", q_gene =", q_gene, ":\n")
    cat("  Genes:    GLM =", glm_row$genes_significant, ", Wilcoxon =", wilcox_row$genes_significant, "\n")
    cat("  Pathways: GLM =", glm_row$pathways_total, ", Wilcoxon =", wilcox_row$pathways_total, "\n")
  }
}

# Save master results
write.csv(master_results, paste0(files_dir, "master_comparison_summary.csv"), row.names = FALSE)

cat("\n\nOutput location:", compare_dir, "\n")
cat("Master summary:", paste0(files_dir, "master_comparison_summary.csv"), "\n")

cat("\n Script completed successfully\n")
