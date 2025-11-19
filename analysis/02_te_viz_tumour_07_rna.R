#!/usr/bin/env Rscript

# Tumour TE Visualization - RNA Expression
# RNA processing and pathway analysis for genes with expression effects

# Source common setup and load data
source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/00_viz_common_setup.R")
source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/00_viz_load_data_tumour.R")

cat("Running 02_te_viz_tumour_08_rna.R...\n")

#### PARAMETER SWEEP CONFIGURATION ####

# Define parameter grid for RNA pathway analysis
param_grid_rna <- expand.grid(
  min_samples = c(3),
  analysis_type = c("coord", "gene"),
  p_gene = c(0.05, 0.1),
  p_pathway = c(0.05, 0.1),
  q_pathway = c(0.05, 0.1),
  stringsAsFactors = FALSE
)
cat("\n*** RNA PATHWAY ANALYSIS PARAMETER SWEEP ***\n")
cat("Testing", nrow(param_grid_rna), "parameter combinations\n\n")

#### PROCESS RNA DATA ####
# Validate, rename, combine and filter RNA datasets from KICS, LFS, and St. Jude
# KICS now uses matched_dna_rna.csv for direct DNA-RNA mapping
rna_processing_result <- process_tumor_rna_data(kics_rna, lfs_rna, stjude_rna, lfs_wgs2rna)
rna_ready_for_filtering <- rna_processing_result$rna_data

cat("Tumor RNA data prepared for filtering. Will complete after TE data is loaded.\n")

# Match TE sample names directly with RNA sample names (now both use DNA naming)
# Uses loaded matched_dna_rna for mapping, exports unmatched samples to CSV
rna_matching_result <- match_te_rna_samples(rna_ready_for_filtering, te_all_t, r_dir_files = r_dir_files, dna_rna_mapping = matched_dna_rna, original_kics_rna = kics_rna, te_all_all_t = te_all_all_t)
rna_filtered <- rna_matching_result$rna_filtered
sample_overlap <- rna_matching_result$sample_overlap

# RNA gene expression analyses (COMMENTED OUT - replaced by updated section below)
# # Compare MROH7 expression between LFS vs control groups
# write_output(quote(plot_gene_expression_te(rna_filtered, te_all_split_t, "MROH7", "TP53_status", "TP53 Status", "MROH7 Expression (FPKM)", FALSE, FALSE)), "MROH7 Expression by TP53 Status")
# titled_print(plot_gene_expression_te(rna_filtered, te_all_split_t, "MROH7", "TP53_status", "TP53 Status", "MROH7 Expression (FPKM)", FALSE, FALSE), "MROH7 Expression by TP53 Status")
# ggsave(paste0(plot_dir, "mroh7_expression_tp53_status.png"), width = 9, height = 5)
# 
# # Compare MROH7 expression based on TE presence/absence in MROH7 gene
# write_output(quote(plot_gene_expression_te(rna_filtered, te_all_split_t, "MROH7", "TP53_status", "TE Status", "MROH7 Expression (FPKM)", FALSE, TRUE)), "MROH7 Expression by TE Status")
# titled_print(plot_gene_expression_te(rna_filtered, te_all_split_t, "MROH7", "TP53_status", "TE Status", "MROH7 Expression (FPKM)", FALSE, TRUE), "MROH7 Expression by TE Status")
# ggsave(paste0(plot_dir, "mroh7_expression_te_status.png"), width = 9, height = 5)
# 
# # Compare MROH7 expression across different tumor types
# write_output(quote(plot_gene_expression_te(rna_filtered, te_all_split_t, "MROH7", "tumor_type", "Tumor Type", "MROH7 Expression (FPKM)", FALSE, FALSE)), "MROH7 Expression by Tumor Type")
# titled_print(plot_gene_expression_te(rna_filtered, te_all_split_t, "MROH7", "tumor_type", "Tumor Type", "MROH7 Expression (FPKM)", FALSE, FALSE), "MROH7 Expression by Tumor Type")
# ggsave(paste0(plot_dir, "mroh7_expression_tumor_type.png"), width = 9, height = 5)
# 
# # Compare MROH7 expression by TE status within tumor type groups  
# write_output(quote(plot_gene_expression_te(rna_filtered, te_all_split_t, "MROH7", "tumor_type", "TE Status", "MROH7 Expression (FPKM)", FALSE, TRUE)), "MROH7 Expression by TE Status (grouped by tumor type)")
# titled_print(plot_gene_expression_te(rna_filtered, te_all_split_t, "MROH7", "tumor_type", "TE Status", "MROH7 Expression (FPKM)", FALSE, TRUE), "MROH7 Expression by TE Status (grouped by tumor type)")
# ggsave(paste0(plot_dir, "mroh7_expression_te_status_by_tumor_type.png"), width = 9, height = 5)



# Generate RNA expression plots with complete integration

# Plot 1: MROH7 Expression by TE Status (different grouping variable)
write_output(quote(plot_gene_expression_te(rna_filtered, te_all_split_t, "MROH7", "TP53_status", "TE Status", "MROH7 Expression (FPKM)", FALSE, TRUE, FALSE)), "MROH7 Expression by TE Status")
titled_print(plot_gene_expression_te(rna_filtered, te_all_split_t, "MROH7", "TP53_status", "TE Status", "MROH7 Expression (FPKM)", FALSE, TRUE, FALSE)$plot, "MROH7 Expression by TE Status")
ggsave(paste0(plot_dir, "rna/mroh7_expression_te_status.png"), width = 9, height = 5)

# Plot 2: MROH7 Expression by TP53 Status with stats (stats=TRUE, no facet)
write_output(quote(plot_gene_expression_te(rna_filtered, te_all_split_t, "MROH7", "TP53_status", "TP53 Status", "MROH7 Expression (FPKM)", FALSE, FALSE, FALSE)), "MROH7 Expression by TP53 Status (with stats)")
titled_print(plot_gene_expression_te(rna_filtered, te_all_split_t, "MROH7", "TP53_status", "TP53 Status", "MROH7 Expression (FPKM)", FALSE, FALSE, FALSE)$plot, "MROH7 Expression by TP53 Status (with stats)")
ggsave(paste0(plot_dir, "rna/mroh7_expression_tp53_status_stats.png"), width = 9, height = 5)

# Plot 3: MROH7 Expression by Tumor Type (no stats, no facet)
write_output(quote(plot_gene_expression_te(rna_filtered, te_all_split_t, "MROH7", "tumor_type", "Tumor Type", "MROH7 Expression (FPKM)", FALSE, FALSE, FALSE)), "MROH7 Expression by Tumor Type")
titled_print(plot_gene_expression_te(rna_filtered, te_all_split_t, "MROH7", "tumor_type", "Tumor Type", "MROH7 Expression (FPKM)", FALSE, FALSE, FALSE)$plot, "MROH7 Expression by Tumor Type")
ggsave(paste0(plot_dir, "rna/mroh7_expression_tumor_type.png"), width = 9, height = 5)

# Plot 4: MROH7 Expression by Tumor Type with faceting by TP53 Status (no stats, facet=TRUE)
write_output(quote(plot_gene_expression_te(rna_filtered, te_all_split_t, "MROH7", "tumor_type", "Tumor Type", "MROH7 Expression (FPKM)", FALSE, FALSE, FALSE)), "MROH7 Expression by Tumor Type (faceted by TP53)")
titled_print(plot_gene_expression_te(rna_filtered, te_all_split_t, "MROH7", "tumor_type", "Tumor Type", "MROH7 Expression (FPKM)", FALSE, FALSE, FALSE)$plot, "MROH7 Expression by Tumor Type (faceted by TP53)")
ggsave(paste0(plot_dir, "rna/mroh7_expression_tumor_type_faceted.png"), width = 12, height = 5)

cat("\n===== SYSTEMATIC TE EXPRESSION ANALYSIS =====\n")

# Generate gene summary table before systematic analysis
cat("Generating gene summary table...\n")
gene_summary <- te_all_split_t %>%
  group_by(Gene_name) %>%
  summarize(
    total_samples_affected = n_distinct(sample),
    total_te_insertions = n(),
    has_exonic = any(grepl("exon", Location, ignore.case = TRUE), na.rm = TRUE),
    has_intronic = any(grepl("intron", Location, ignore.case = TRUE), na.rm = TRUE),
    location_summary = paste(unique(na.omit(Location)), collapse = "; "),
    gene_features = paste(unique(na.omit(Location2)), collapse = "; "),
    .groups = "drop"
  ) %>%
  arrange(desc(total_samples_affected), desc(total_te_insertions)) %>%
  mutate(
    te_location_type = case_when(
      has_exonic & has_intronic ~ "Both exonic & intronic",
      has_exonic ~ "Exonic only",
      has_intronic ~ "Intronic only",
      TRUE ~ "Other"
    )
  ) %>%
  select(Gene_name, total_samples_affected, total_te_insertions, te_location_type, location_summary, gene_features)

# Save to CSV
gene_summary_file <- paste0(r_dir_files, "gene_summary_te_affected.csv")
write.csv(gene_summary, gene_summary_file, row.names = FALSE)
cat(sprintf("✓ Gene summary saved to: %s\n", basename(gene_summary_file)))
cat(sprintf("  Total genes analyzed: %d\n", nrow(gene_summary)))
cat(sprintf("  TEs by location type:\n"))
cat(sprintf("    - Exonic only: %d genes (%.1f%%)\n",
            sum(gene_summary$te_location_type == "Exonic only"),
            100 * sum(gene_summary$te_location_type == "Exonic only") / nrow(gene_summary)))
cat(sprintf("    - Intronic only: %d genes (%.1f%%)\n",
            sum(gene_summary$te_location_type == "Intronic only"),
            100 * sum(gene_summary$te_location_type == "Intronic only") / nrow(gene_summary)))
cat(sprintf("    - Both: %d genes (%.1f%%)\n",
            sum(gene_summary$te_location_type == "Both exonic & intronic"),
            100 * sum(gene_summary$te_location_type == "Both exonic & intronic") / nrow(gene_summary)))
cat(sprintf("  Top gene (most samples): %s (%d samples)\n\n",
            gene_summary$Gene_name[1],
            gene_summary$total_samples_affected[1]))

# Generate gene summary for KICS samples
cat("Generating gene summary table for KICS samples...\n")
gene_summary_kics <- te_kics_split_t %>%
  group_by(Gene_name) %>%
  summarize(
    total_samples_affected = n_distinct(sample),
    total_te_insertions = n(),
    has_exonic = any(grepl("exon", Location, ignore.case = TRUE), na.rm = TRUE),
    has_intronic = any(grepl("intron", Location, ignore.case = TRUE), na.rm = TRUE),
    location_summary = paste(unique(na.omit(Location)), collapse = "; "),
    gene_features = paste(unique(na.omit(Location2)), collapse = "; "),
    .groups = "drop"
  ) %>%
  arrange(desc(total_samples_affected), desc(total_te_insertions)) %>%
  mutate(
    te_location_type = case_when(
      has_exonic & has_intronic ~ "Both exonic & intronic",
      has_exonic ~ "Exonic only",
      has_intronic ~ "Intronic only",
      TRUE ~ "Other"
    )
  ) %>%
  select(Gene_name, total_samples_affected, total_te_insertions, te_location_type, location_summary, gene_features)

gene_summary_kics_file <- paste0(r_dir_files, "gene_summary_te_kics.csv")
write.csv(gene_summary_kics, gene_summary_kics_file, row.names = FALSE)
cat(sprintf("✓ Gene summary for KICS saved to: %s\n", basename(gene_summary_kics_file)))

# Generate gene summary for LFS samples
cat("Generating gene summary table for LFS samples...\n")
gene_summary_lfs <- te_lfs_split_t %>%
  group_by(Gene_name) %>%
  summarize(
    total_samples_affected = n_distinct(sample),
    total_te_insertions = n(),
    has_exonic = any(grepl("exon", Location, ignore.case = TRUE), na.rm = TRUE),
    has_intronic = any(grepl("intron", Location, ignore.case = TRUE), na.rm = TRUE),
    location_summary = paste(unique(na.omit(Location)), collapse = "; "),
    gene_features = paste(unique(na.omit(Location2)), collapse = "; "),
    .groups = "drop"
  ) %>%
  arrange(desc(total_samples_affected), desc(total_te_insertions)) %>%
  mutate(
    te_location_type = case_when(
      has_exonic & has_intronic ~ "Both exonic & intronic",
      has_exonic ~ "Exonic only",
      has_intronic ~ "Intronic only",
      TRUE ~ "Other"
    )
  ) %>%
  select(Gene_name, total_samples_affected, total_te_insertions, te_location_type, location_summary, gene_features)

gene_summary_lfs_file <- paste0(r_dir_files, "gene_summary_te_lfs.csv")
write.csv(gene_summary_lfs, gene_summary_lfs_file, row.names = FALSE)
cat(sprintf("✓ Gene summary for LFS saved to: %s\n", basename(gene_summary_lfs_file)))

# Generate cancer genes only gene summary for all samples
cat("Generating cancer genes only gene summary for all samples...\n")
# Define cancer predisposition genes list
cpg <- read.csv("/Users/briannelaverty/Documents/R_Malkin/te/data/raw/kics_cpg.txt", header=FALSE, sep="\t")
genes <- cpg$V1
gene_summary_cancer <- gene_summary %>%
  filter(Gene_name %in% genes)

gene_summary_cancer_file <- paste0(r_dir_files, "gene_summary_te_all_cancer_genes.csv")
write.csv(gene_summary_cancer, gene_summary_cancer_file, row.names = FALSE)
cat(sprintf("✓ Cancer genes only gene summary saved to: %s\n", basename(gene_summary_cancer_file)))
cat(sprintf("  Total cancer genes with TEs: %d\n\n", nrow(gene_summary_cancer)))

# Generate cancer genes only gene summary for KICS samples
cat("Generating cancer genes only gene summary for KICS samples...\n")
gene_summary_kics_cancer <- gene_summary_kics %>%
  filter(Gene_name %in% genes)

gene_summary_kics_cancer_file <- paste0(r_dir_files, "gene_summary_te_kics_cancer_genes.csv")
write.csv(gene_summary_kics_cancer, gene_summary_kics_cancer_file, row.names = FALSE)
cat(sprintf("✓ Cancer genes only gene summary for KICS saved to: %s\n", basename(gene_summary_kics_cancer_file)))
cat(sprintf("  Total cancer genes with TEs (KICS): %d\n\n", nrow(gene_summary_kics_cancer)))

# Generate cancer genes only gene summary for LFS samples
cat("Generating cancer genes only gene summary for LFS samples...\n")
gene_summary_lfs_cancer <- gene_summary_lfs %>%
  filter(Gene_name %in% genes)

gene_summary_lfs_cancer_file <- paste0(r_dir_files, "gene_summary_te_lfs_cancer_genes.csv")
write.csv(gene_summary_lfs_cancer, gene_summary_lfs_cancer_file, row.names = FALSE)
cat(sprintf("✓ Cancer genes only gene summary for LFS saved to: %s\n", basename(gene_summary_lfs_cancer_file)))
cat(sprintf("  Total cancer genes with TEs (LFS): %d\n\n", nrow(gene_summary_lfs_cancer)))

# Systematic TE expression analysis - test all TEs for significant expression effects (by coordinates)
systematic_results <- run_systematic_te_expression_analysis(
  rna_filtered = rna_filtered,
  te_split_df = te_all_split_t,
  plot_dir = plot_dir,
  min_samples = 5,
  fdr_cutoff = 0.05,
  create_plots = TRUE,
  max_plots = 3,
  group_by_gene = FALSE,
  files_dir = r_dir_files
)

write_output(quote(systematic_results), "Systematic TE Expression Analysis Results (by coordinates)")

# Additional systematic TE testing with more lenient parameters (by coordinates)
systematic_results_all <- run_systematic_te_expression_analysis(
  rna_filtered = rna_filtered,
  te_split_df = te_all_split_t,
  plot_dir = plot_dir,
  min_samples = 3,  # Lower threshold for tumor to get more results
  fdr_cutoff = 0.1,  # More lenient for discovery
  create_plots = TRUE,
  max_plots = 3,
  group_by_gene = FALSE,
  files_dir = r_dir_files
)

write_output(quote(systematic_results_all), "Comprehensive Systematic TE Expression Analysis Results (by coordinates)")

# Systematic TE expression analysis - grouped by gene (aggregate all TEs per gene)
systematic_results_by_gene <- run_systematic_te_expression_analysis(
  rna_filtered = rna_filtered,
  te_split_df = te_all_split_t,
  plot_dir = plot_dir,
  min_samples = 5,
  fdr_cutoff = 0.05,
  create_plots = TRUE,
  max_plots = 3,
  group_by_gene = TRUE,
  files_dir = r_dir_files
)

write_output(quote(systematic_results_by_gene), "Systematic TE Expression Analysis Results (grouped by gene)")

# Additional gene-grouped testing with more lenient parameters
systematic_results_by_gene_all <- run_systematic_te_expression_analysis(
  rna_filtered = rna_filtered,
  te_split_df = te_all_split_t,
  plot_dir = plot_dir,
  min_samples = 3,  # Lower threshold for tumor to get more results
  fdr_cutoff = 0.1,  # More lenient for discovery
  create_plots = TRUE,
  max_plots = 3,
  group_by_gene = TRUE,
  files_dir = r_dir_files
)

write_output(quote(systematic_results_by_gene_all), "Comprehensive Systematic TE Expression Analysis Results (grouped by gene)")

# Check if RNA expression analysis found any testable TEs
if (nrow(systematic_results) == 0 && nrow(systematic_results_all) == 0 &&
    nrow(systematic_results_by_gene) == 0 && nrow(systematic_results_by_gene_all) == 0) {

  cat("\n=== Creating explanatory note for RNA expression analysis ===\n")

  # Create explanatory note file
  note_content <- paste0(
    "================================================================================\n",
    "TUMOR RNA EXPRESSION ANALYSIS - NO RESULTS\n",
    "================================================================================\n\n",
    "The systematic TE expression analysis did not produce results for the tumor\n",
    "data because there are insufficient recurrent TEs for statistical testing.\n\n",
    "SUMMARY:\n",
    "--------\n",
    "- By coordinates (min 5 samples): ", nrow(systematic_results), " tests\n",
    "- By coordinates (min 3 samples): ", nrow(systematic_results_all), " tests\n",
    "- By gene (min 5 samples): ", nrow(systematic_results_by_gene), " tests\n",
    "- By gene (min 3 samples): ", nrow(systematic_results_by_gene_all), " tests\n\n",
    "REASON:\n",
    "-------\n",
    "Somatic/tumor TEs are mostly unique to individual tumors (patient-specific),\n",
    "unlike germline TEs which are inherited and recurrent across many individuals.\n\n",
    "The analysis requires TEs to be present in multiple samples to test for\n",
    "differential expression effects, but tumor TEs are too rare/unique for this\n",
    "type of analysis.\n\n",
    "EXPECTED OUTPUT FILES THAT DO NOT EXIST:\n",
    "-----------------------------------------\n",
    "- rna_results_5samples_by_coordinates.csv\n",
    "- rna_results_3samples_by_coordinates.csv\n",
    "- rna_results_5samples_by_gene.csv\n",
    "- rna_results_3samples_by_gene.csv\n",
    "- significant_te_expression_results_*.csv\n",
    "- te_genes_with_expression_effects_*.csv\n",
    "- te_pathway_ora_expression_*.csv (pathway analysis)\n\n",
    "These files ARE available for germline analysis where recurrent TEs are common.\n\n",
    "ALTERNATIVE ANALYSIS:\n",
    "--------------------\n",
    "For tumor-specific RNA effects, use:\n",
    "- test_re_rna_tumour.R (tests regulatory element genes)\n",
    "- Results: re_gene_rna_differential_expression_tumour.csv\n",
    "- Pathway: re_rna_pathway_enrichment_tumour_p0.1.csv\n\n",
    "Date: ", Sys.Date(), "\n",
    "Script: 02_te_viz_tumour.R (lines 1924-1982)\n",
    "================================================================================\n"
  )
  writeLines(note_content, paste0(r_dir_files, "RNA_EXPRESSION_ANALYSIS_NOTE.txt"))
  cat("✓ Created explanatory note: RNA_EXPRESSION_ANALYSIS_NOTE.txt\n\n")
}


#### PATHWAY ANALYSIS - TEs WITH EXPRESSION EFFECTS WITH PARAMETER SWEEP ####
write_output(quote(NULL), "Pathway Analysis - TEs with Expression Effects")

# Set up summary tracking for parameter sweep
sweep_summary <- data.frame()

# Loop through parameter combinations
for (i in 1:nrow(param_grid_rna)) {
  params <- param_grid_rna[i, ]

  # Create parameter suffix for filenames
  param_suffix <- paste0(
    "_min", params$min_samples,
    "_", params$analysis_type,
    "_pgene", params$p_gene,
    "_ppathway", params$p_pathway,
    "_qpathway", params$q_pathway
  )

  cat("\n\n===== TESTING RNA PATHWAY PARAMETERS", i, "/", nrow(param_grid_rna), "=====\n")
  cat("min_samples =", params$min_samples,
      ", analysis_type =", params$analysis_type,
      ", p_gene =", params$p_gene,
      ", p_pathway =", params$p_pathway,
      ", q_pathway =", params$q_pathway, "\n")

  # Select the appropriate systematic results based on parameters
  if (params$analysis_type == "coord") {
    if (params$min_samples == 5) {
      systematic_results_to_use <- systematic_results
    }
  } else if (params$analysis_type == "gene") {
    if (params$min_samples == 5) {
      systematic_results_to_use <- systematic_results_by_gene
    }
  }

  # Get genes with expression p < threshold from systematic expression analysis
  if (exists("systematic_results_to_use") && nrow(systematic_results_to_use) > 0 && "p_value" %in% colnames(systematic_results_to_use)) {
    # Get genes with p < threshold directly from systematic results
    significant_genes_t <- systematic_results_to_use %>%
      filter(p_value < params$p_gene) %>%
      pull(gene) %>%
      unique()

    cat("Genes with expression p-value <", params$p_gene, ":", length(significant_genes_t), "\n")

    # Perform pathway analysis if enough genes
    if (length(significant_genes_t) >= 5) {
      genes_df_t <- data.frame(Gene_name = significant_genes_t)

      ora_expression_t <- tryCatch({
        perform_ora_custom_cutoffs(genes_df_t,
                                  p_pathway = params$p_pathway,
                                  q_pathway = params$q_pathway,
                                  nsample_thresh = 0,
                                  filter_exon = FALSE)
      }, error = function(e) {
        cat("Error in pathway analysis:", e$message, "\n")
        NULL
      })

      if (!is.null(ora_expression_t) && nrow(as.data.frame(ora_expression_t)) > 0) {
        n_pathways <- nrow(as.data.frame(ora_expression_t))
        cat("✓ Found", n_pathways, "enriched pathways\n")

        # Save pathway results
        pathway_file <- paste0(r_dir_files, "rna_pathway_tumour", param_suffix, ".csv")
        write.csv(as.data.frame(ora_expression_t), pathway_file, row.names=FALSE)

        # Create plots
        tryCatch({
            p_dot_ora_expression_t <- dotplot(ora_expression_t, showCategory=20)
            titled_print(p_dot_ora_expression_t, "ORA dotplot (Expression TEs)")
            ggsave(paste0(plot_dir, "pathway/rna_pathway_tumour_dot_", param_suffix, ".png"), plot=p_dot_ora_expression_t, width=14, height=9)
          }, error = function(e) cat("Warning: Could not create dotplot:", e$message, "\n"))

          tryCatch({
            p_cnet_ora_expression_t <- cnetplot(ora_expression_t, showCategory=10, colorEdge=TRUE, node_label="category")
            titled_print(p_cnet_ora_expression_t, "ORA cnetplot (Expression TEs)")
            ggsave(paste0(plot_dir, "pathway/rna_pathway_tumour_cnet_", param_suffix, ".png"), plot=p_cnet_ora_expression_t, width=14, height=9)
          }, error = function(e) cat("Warning: Could not create cnetplot:", e$message, "\n"))

          tryCatch({
            ora_expression_pairwise_t <- pairwise_termsim(ora_expression_t)
            p_emap_ora_expression_t <- emapplot(ora_expression_pairwise_t, pie="count", showCategory=20, group_category=TRUE, group_legend=TRUE)
            titled_print(p_emap_ora_expression_t, "ORA emapplot (Expression TEs)")
            ggsave(paste0(plot_dir, "pathway/rna_pathway_tumour_emap_", param_suffix, ".png"), plot=p_emap_ora_expression_t, width=14, height=9)
          }, error = function(e) cat("Warning: Could not create emapplot:", e$message, "\n"))

          tryCatch({
            pathway_genes_expression_t <- plot_pathway_gene_counts(ora_expression_t, n_descriptions=nrow(ora_expression_t), n_genes=50)
            write.csv(pathway_genes_expression_t, paste0(r_dir_files, "pathway_genes_rna_tumour_", param_suffix, ".csv"), row.names=FALSE)
          }, error = function(e) cat("Warning: Could not create pathway gene counts:", e$message, "\n"))

        # Record sweep summary
        sweep_summary <- rbind(sweep_summary, data.frame(
          min_samples = params$min_samples,
          analysis_type = params$analysis_type,
          p_gene = params$p_gene,
          p_pathway = params$p_pathway,
          q_pathway = params$q_pathway,
          n_genes = length(significant_genes_t),
          n_pathways = n_pathways,
          stringsAsFactors = FALSE
        ))
      } else {
        cat("No significant pathways found\n")
        sweep_summary <- rbind(sweep_summary, data.frame(
          min_samples = params$min_samples,
          analysis_type = params$analysis_type,
          p_gene = params$p_gene,
          p_pathway = params$p_pathway,
          q_pathway = params$q_pathway,
          n_genes = length(significant_genes_t),
          n_pathways = 0,
          stringsAsFactors = FALSE
        ))
      }
    } else {
      cat("Not enough genes (< 5) for pathway analysis\n")
      sweep_summary <- rbind(sweep_summary, data.frame(
        min_samples = params$min_samples,
        analysis_type = params$analysis_type,
        p_gene = params$p_gene,
        p_pathway = params$p_pathway,
        q_pathway = params$q_pathway,
        n_genes = length(significant_genes_t),
        n_pathways = NA,
        stringsAsFactors = FALSE
      ))
    }
  } else {
    cat("No systematic expression results available for this combination\n")
    sweep_summary <- rbind(sweep_summary, data.frame(
      min_samples = params$min_samples,
      analysis_type = params$analysis_type,
      p_gene = params$p_gene,
      p_pathway = params$p_pathway,
      q_pathway = params$q_pathway,
      n_genes = 0,
      n_pathways = NA,
      stringsAsFactors = FALSE
    ))
  }
}

# Save sweep summary
if (nrow(sweep_summary) > 0) {
  summary_file <- paste0(r_dir_files, "rna_pathway_sweep_tumour_summary.csv")
  write.csv(sweep_summary, summary_file, row.names = FALSE)
  cat("\n✓ Parameter sweep complete! Summary saved to:", summary_file, "\n")
  cat("\nTop 5 parameter combinations by number of pathways:\n")
  print(head(sweep_summary[order(-sweep_summary$n_pathways), ], 5))
}

##### ANCESTRY ANALYSIS ####
#cat("\n========================================\n")
#cat("ANCESTRY ANALYSIS\n")
#cat("========================================\n\n")
#
## Load ancestry data
#cat("Loading ancestry data...\n")
#load(paste0(r_dir, "ancestry.RData"))
#
## Load location windows data
#cat("Loading location windows data...\n")
#location_100kb_t <- read.csv("/hpf/largeprojects/davidm/blaverty/te/ml/output/tumour/location_csv/100kb_complete_filtered_t.csv", stringsAsFactors = FALSE)
#
## Merge location windows with ancestry using same approach as processing script
#location_ancestry <- write_output(
#  quote(merge_location_ancestry(location_100kb_t, ancestry)),
#  "Merging location windows with ancestry"
#)
#
## Identify clinical columns to exclude from PCA/UMAP
#clinical_cols <- c("sample", "predicted_ancestry_thres", "mapped_label", "age", "sex", "cohort", "tumor_type",
#                   "TP53_germline", "affected", "cancer", "age_diagnosis", "TP53_status")
#exclude_cols <- intersect(clinical_cols, colnames(location_ancestry))
#
## Perform PCA on location windows
#pca_results <- write_output(
#  quote(perform_location_pca(location_ancestry, exclude_cols = exclude_cols)),
#  "Performing PCA on location windows"
#)
#
## Perform UMAP on location windows
#umap_results <- write_output(
#  quote(perform_location_umap(location_ancestry, exclude_cols = exclude_cols)),
#  "Performing UMAP on location windows"
#)
#
## Plot PCA colored by predicted_ancestry_thres
#p_pca_ancestry <- write_output(
#  quote(plot_pca_ancestry(pca_results, color_by = "predicted_ancestry_thres",
#                         output_dir = plot_dir, plot_prefix = "tumour")),
#  "Plotting PCA colored by predicted ancestry"
#)
#titled_print(p_pca_ancestry, "PCA of TE Location Windows by Predicted Ancestry (Tumour)")
#
## Plot PCA colored by mapped_label
#p_pca_mapped <- write_output(
#  quote(plot_pca_ancestry(pca_results, color_by = "mapped_label",
#                         output_dir = plot_dir, plot_prefix = "tumour")),
#  "Plotting PCA colored by mapped label"
#)
#titled_print(p_pca_mapped, "PCA of TE Location Windows by Mapped Label (Tumour)")
#
## Plot UMAP colored by predicted_ancestry_thres
#p_umap_ancestry <- write_output(
#  quote(plot_umap_ancestry(umap_results, color_by = "predicted_ancestry_thres",
#                          output_dir = plot_dir, plot_prefix = "tumour")),
#  "Plotting UMAP colored by predicted ancestry"
#)
#titled_print(p_umap_ancestry, "UMAP of TE Location Windows by Predicted Ancestry (Tumour)")
#
## Plot UMAP colored by mapped_label
#p_umap_mapped <- write_output(
#  quote(plot_umap_ancestry(umap_results, color_by = "mapped_label",
#                          output_dir = plot_dir, plot_prefix = "tumour")),
#  "Plotting UMAP colored by mapped label"
#)
#titled_print(p_umap_mapped, "UMAP of TE Location Windows by Mapped Label (Tumour)")
#
# Plot count LM grouped by predicted_ancestry_thres
# Check if we have enough samples with ancestry data
#if ("predicted_ancestry_thres" %in% colnames(te_aff_t) &&
#    sum(!is.na(te_aff_t$predicted_ancestry_thres)) >= 10) {
#
#  ancestry_palette <- scales::hue_pal()(length(unique(te_aff_t$predicted_ancestry_thres[!is.na(te_aff_t$predicted_ancestry_thres)])))
#
#  p_lm_ancestry_pred <- tryCatch({
#    write_output(
#      quote(plot_count_lm(te_aff_t, min_samples=5, residuals=FALSE, group="predicted_ancestry_thres",
#                          log_scale=TRUE, covariates = covar_med, x_lab="Predicted Ancestry",
#                          y_lab="Total TE count", type=NA, chr=NA, fill_palette=ancestry_palette)),
#      "Linear model by predicted ancestry thres (te_aff_t, all types, log scale)"
#    )
#  }, error = function(e) {
#    cat("Skipping ancestry LM plot - insufficient data or error:", e$message, "\n")
#    NULL
#  })
#
#  if (!is.null(p_lm_ancestry_pred)) {
#    titled_print(p_lm_ancestry_pred, "Linear model by predicted ancestry thres (te_aff_t, all types, log scale)")
#    ggsave(paste0(plot_dir, "te_count_lm_predicted_ancestry_thres_aff_all_t.png"), plot=p_lm_ancestry_pred, width = 9, height = 5)
#  }
#} else {
#  cat("Skipping predicted_ancestry_thres LM plot - insufficient samples with ancestry data\n")
#  cat("Samples with ancestry:", sum(!is.na(te_aff_t$predicted_ancestry_thres)), "\n")
#}
#
## Plot count LM grouped by mapped_label
#if ("mapped_label" %in% colnames(te_aff_t) &&
#    sum(!is.na(te_aff_t$mapped_label)) >= 10) {
#
#  mapped_label_palette <- scales::hue_pal()(length(unique(te_aff_t$mapped_label[!is.na(te_aff_t$mapped_label)])))
#
#  p_lm_ancestry_mapped <- tryCatch({
#    write_output(
#      quote(plot_count_lm(te_aff_t, min_samples=5, residuals=FALSE, group="mapped_label",
#                          log_scale=TRUE, covariates = covar_med, x_lab="Mapped Label",
#                          y_lab="Total TE count", type=NA, chr=NA, fill_palette=mapped_label_palette)),
#      "Linear model by mapped label (te_aff_t, all types, log scale)"
#    )
#  }, error = function(e) {
#    cat("Skipping mapped_label LM plot - insufficient data or error:", e$message, "\n")
#    NULL
#  })
#
#  if (!is.null(p_lm_ancestry_mapped)) {
#    titled_print(p_lm_ancestry_mapped, "Linear model by mapped label (te_aff_t, all types, log scale)")
#    ggsave(paste0(plot_dir, "te_count_lm_mapped_label_aff_all_t.png"), plot=p_lm_ancestry_mapped, width = 9, height = 5)
#  }
#} else {
#  cat("Skipping mapped_label LM plot - insufficient samples with ancestry data\n")
#  cat("Samples with mapped_label:", sum(!is.na(te_aff_t$mapped_label)), "\n")
#}
#
## Plot pie chart for ancestry
#cat("Plotting: clinical_ancestry_tumour.pdf\n")
#p_ancestry_pie_dataset <- plot_ancestry_pie(te_all_t, output_dir = NULL, plot_prefix = "tumour")
#print(p_ancestry_pie_dataset)
#ggsave(paste0(plot_dir, "dataset/clinical_ancestry_tumour.pdf"), width = 6, height = 5)
#
##### PCA/UMAP BY OTHER VARIABLES ####
#cat("\n========================================\n")
#cat("PCA/UMAP COLORED BY OTHER VARIABLES\n")
#cat("========================================\n\n")
#
## Variables to color by (includes both ancestry variables and clinical variables)
#color_variables <- c("mapped_label", "tumor_type", "cohort", "TP53_status", "sex")
#
## Plot PCA by each variable
#for (var in color_variables) {
#  if (var %in% colnames(location_ancestry)) {
#    cat("Plotting PCA colored by", var, "...\n")
#    p_pca <- write_output(
#      quote(plot_pca_by_variable(pca_results, color_by = var,
#                                  output_dir = plot_dir, plot_prefix = "tumour")),
#      paste0("PCA colored by ", var)
#    )
#    if (!is.null(p_pca)) {
#      titled_print(p_pca, paste0("PCA of TE Location Windows by ", var, " (Tumour)"))
#    }
#  } else {
#    cat("Skipping PCA for", var, "- column not found in data\n")
#  }
#}
#
## Plot UMAP by each variable
#for (var in color_variables) {
#  if (var %in% colnames(location_ancestry)) {
#    cat("Plotting UMAP colored by", var, "...\n")
#    p_umap <- write_output(
#      quote(plot_umap_by_variable(umap_results, color_by = var,
#                                   output_dir = plot_dir, plot_prefix = "tumour")),
#      paste0("UMAP colored by ", var)
#    )
#    if (!is.null(p_umap)) {
#      titled_print(p_umap, paste0("UMAP of TE Location Windows by ", var, " (Tumour)"))
#    }
#  } else {
#    cat("Skipping UMAP for", var, "- column not found in data\n")
#  }
#}

#### TAYLOR COHORT ANALYSIS ####
write_output(quote({
  cat("Taylor Cohort Summary:\n")
  cat("Total samples (selected):", nrow(te_taylor_t), "\n")
  cat("Total samples (all):", nrow(te_taylor_all_t), "\n")
  cat("\nTumor type distribution:\n")
  print(table(te_taylor_t$tumor_type))
  cat("\nTumor type subclass distribution:\n")
  print(table(te_taylor_t$tumor_type_subclass))
}), "Taylor Cohort Summary Statistics")

write_output(quote(plot_count_kruskal(df = te_taylor_t, chr = NA, type = "total", group = "tumor_type_subclass", x_lab = "Tumor Type Subclass", y_lab = "Total TE Count", log_scale = FALSE)),
             "Taylor: TE Count by Tumor Type Subclass (Kruskal-Wallis)")
p_kruskal <- plot_count_kruskal(df = te_taylor_t, chr = NA, type = "total", group = "tumor_type_subclass", x_lab = "Tumor Type Subclass", y_lab = "Total TE Count", log_scale = FALSE)
titled_print(p_kruskal, "Taylor: TE Count by Tumor Type Subclass (Kruskal-Wallis)")
ggsave(paste0(plot_dir, "other/taylor_te_count_by_subclass_kruskal.png"), plot = p_kruskal, width = 10, height = 6)

write_output(quote(plot_count_age(df = te_taylor_t, type = "total", chr = NA, y_lab = "Total TE Count")),
             "Taylor: TE Count by Age at Diagnosis")
p_age <- plot_count_age(df = te_taylor_t, type = "total", chr = NA, y_lab = "Total TE Count")
titled_print(p_age, "Taylor: TE Count by Age at Diagnosis")
ggsave(paste0(plot_dir, "other/taylor_te_count_by_age.png"), plot = p_age, width = 8, height = 6)

# Save rna_filtered for use by RE script
cat("\nSaving rna_filtered for RE script...\n")
saveRDS(rna_filtered, paste0(r_dir_files, "rna_filtered_tumour.rds"))
cat("✓ Saved rna_filtered to:", paste0(r_dir_files, "rna_filtered_tumour.rds\n"))

# Close the PDF device at the end
dev.off()

cat("\n===== SCRIPT COMPLETED SUCCESSFULLY =====\n")
cat("All analysis sections have been processed.\n")
cat("Generated plots saved to:", plot_dir, "\n")
cat("PDF compilation saved to: graph_output.pdf\n")
cat("Text output saved to:", stdout_file, "\n")

# Close sink to stop redirecting output
sink()

cat("✓ Script completed successfully\n")
