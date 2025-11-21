#!/usr/bin/env Rscript

# Germline TE Visualization - RNA Expression
# RNA processing, expression analysis, and RNA pathway analysis

# Source common setup and load data
source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/00_viz_common_setup.R")
source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/00_viz_load_data_germline.R")

cat("Running 02_te_viz_germline_07_rna.R...\n")

#### PARAMETER SWEEP CONFIGURATION ####

# Define parameter grid for RNA pathway analysis
param_grid_rna <- expand.grid(
  min_samples = c(5),
  analysis_type = c("coord", "gene"),
  p_gene = c(0.05, 0.1),
  p_pathway = c(0.05, 0.1),
  q_pathway = c(0.05, 0.1),
  stringsAsFactors = FALSE
)
cat("\n*** RNA PATHWAY ANALYSIS PARAMETER SWEEP ***\n")
cat("Testing", nrow(param_grid_rna), "parameter combinations\n\n")

#### PROCESS RNA DATA ####
# Process RNA data using tumor approach (same renaming with matched_dna_rna.csv)
write_output(quote(NULL), "Processing Germline RNA Data")
rna_processing_result <- process_germline_rna_data(kics_rna, lfs_rna, stjude_rna, matched_dna_rna, lfs_wgs2rna)
rna_ready_for_filtering <- rna_processing_result$rna_data

write_output(quote({
  cat("Total RNA samples available:", ncol(rna_ready_for_filtering) - 1, "\n")
  cat("Total genes with expression:", nrow(rna_ready_for_filtering), "\n")
}), "RNA Data Summary")

#### MATCH RNA AND TE SAMPLES ####
write_output(quote(NULL), "Matching Germline TE with Tumor RNA Samples")

# Load tumor-matched RNA sample list
tumor_rna_file <- "/Users/briannelaverty/Library/Mobile Documents/com~apple~CloudDocs/Malkin_lab/figures/te/graphs/tumour/files/successfully_matched_samples.csv"
tumor_rna_samples <- read.csv(tumor_rna_file, stringsAsFactors = FALSE)$sample_name
cat("Tumor RNA samples (from tumor analysis):", length(tumor_rna_samples), "\n")

# Get germline TE samples
te_samples <- unique(te_all$sample)
cat("Germline TE samples:", length(te_samples), "\n")

# Extract base patient IDs
tumor_rna_base_ids <- unique(sub("_.*", "", tumor_rna_samples))
germline_te_base_ids <- unique(sub("_.*", "", te_samples))

cat("Unique tumor RNA patients:", length(tumor_rna_base_ids), "\n")
cat("Unique germline TE patients:", length(germline_te_base_ids), "\n")

# Find which tumor RNA patients are missing from germline TE
missing_patients <- setdiff(tumor_rna_base_ids, germline_te_base_ids)
cat("Tumor patients WITHOUT germline TE:", length(missing_patients), "\n")
if (length(missing_patients) > 0) {
  cat("  First 10 missing:", paste(head(missing_patients, 10), collapse=", "), "\n")
}

# Match tumor RNA samples to germline TE using base patient ID
rna_base_ids <- sub("_.*", "", tumor_rna_samples)
matched_rna_samples <- tumor_rna_samples[rna_base_ids %in% germline_te_base_ids]

cat("Matched RNA samples:", length(matched_rna_samples), "(from", length(unique(sub("_.*", "", matched_rna_samples))), "patients)\n")

# Filter RNA data
rna_filtered <- rna_ready_for_filtering %>%
  select(gene_name, all_of(matched_rna_samples))

cat("Final RNA dataset:", nrow(rna_filtered), "genes ×", ncol(rna_filtered) - 1, "samples\n")

# Save matched samples and missing info
write.csv(data.frame(rna_sample_id = matched_rna_samples),
          paste0(r_dir_files, "germline_matched_rna_samples.csv"),
          row.names = FALSE)

if (length(missing_patients) > 0) {
  write.csv(data.frame(patient_id = missing_patients),
            paste0(r_dir_files, "tumor_patients_missing_germline_te.csv"),
            row.names = FALSE)
  cat("Saved missing patients to tumor_patients_missing_germline_te.csv\n")
}


#### RNA EXPRESSION ANALYSIS ####
# Generate RNA expression plots with complete integration

# Define gene of interest for germline analysis
# Use TP53 since it's a key gene in germline analysis and should be present in RNA data
gene_of_interest <- "TP53"

# Check if we have enough samples and gene exists in the RNA data before plotting
if (ncol(rna_filtered) > 1 && gene_of_interest %in% rna_filtered$gene_name) {
  # Plot 1: Gene Expression by TE Status
  write_output(quote(plot_gene_expression_te(rna_filtered, te_aff_split, gene_of_interest, "TP53_status", "TE Status", paste0(gene_of_interest, " Expression (FPKM)"), FALSE, TRUE, FALSE)), paste0(gene_of_interest, " Expression by TE Status"))
  titled_print(plot_gene_expression_te(rna_filtered, te_aff_split, gene_of_interest, "TP53_status", "TE Status", paste0(gene_of_interest, " Expression (FPKM)"), FALSE, TRUE, FALSE)$plot, paste0(gene_of_interest, " Expression by TE Status"))
  ggsave(paste0(plot_dir, "rna/", tolower(gene_of_interest), "_expression_te_status.png"), width = 9, height = 5)

  # Plot 2: Gene Expression by TP53 Status with stats
  write_output(quote(plot_gene_expression_te(rna_filtered, te_aff_split, gene_of_interest, "TP53_status", "TP53 Status", paste0(gene_of_interest, " Expression (FPKM)"), FALSE, FALSE, FALSE)), paste0(gene_of_interest, " Expression by TP53 Status (with stats)"))
  titled_print(plot_gene_expression_te(rna_filtered, te_aff_split, gene_of_interest, "TP53_status", "TP53 Status", paste0(gene_of_interest, " Expression (FPKM)"), FALSE, FALSE, FALSE)$plot, paste0(gene_of_interest, " Expression by TP53 Status (with stats)"))
  ggsave(paste0(plot_dir, "rna/", tolower(gene_of_interest), "_expression_tp53_status_stats.png"), width = 9, height = 5)

  # Plot 3: Gene Expression by Tumor Type
  write_output(quote(plot_gene_expression_te(rna_filtered, te_aff_split, gene_of_interest, "tumor_type", "Tumor Type", paste0(gene_of_interest, " Expression (FPKM)"), FALSE, FALSE, FALSE)), paste0(gene_of_interest, " Expression by Tumor Type"))
  titled_print(plot_gene_expression_te(rna_filtered, te_aff_split, gene_of_interest, "tumor_type", "Tumor Type", paste0(gene_of_interest, " Expression (FPKM)"), FALSE, FALSE, FALSE)$plot, paste0(gene_of_interest, " Expression by Tumor Type"))
  ggsave(paste0(plot_dir, "rna/", tolower(gene_of_interest), "_expression_tumor_type.png"), width = 9, height = 5)

  # Plot 4: Gene Expression by Tumor Type with faceting by TP53 Status
  write_output(quote(plot_gene_expression_te(rna_filtered, te_aff_split, gene_of_interest, "tumor_type", "Tumor Type", paste0(gene_of_interest, " Expression (FPKM)"), FALSE, FALSE, FALSE)), paste0(gene_of_interest, " Expression by Tumor Type (faceted by TP53)"))
  titled_print(plot_gene_expression_te(rna_filtered, te_aff_split, gene_of_interest, "tumor_type", "Tumor Type", paste0(gene_of_interest, " Expression (FPKM)"), FALSE, FALSE, FALSE)$plot, paste0(gene_of_interest, " Expression by Tumor Type (faceted by TP53)"))
  ggsave(paste0(plot_dir, "rna/", tolower(gene_of_interest), "_expression_tumor_type_faceted.png"), width = 12, height = 5)
} else {
  if (ncol(rna_filtered) <= 1) {
    cat("No matched samples with RNA data; skipping", gene_of_interest, "expression plots\n")
  } else {
    cat(gene_of_interest, "not found in RNA data; skipping expression plots\n")
  }
}

cat("\n===== SYSTEMATIC TE EXPRESSION ANALYSIS =====\n")

# Generate gene summary table before systematic analysis
cat("Generating gene summary table...\n")
gene_summary <- te_aff_split %>%
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
gene_summary_file <- paste0(r_dir_files, "gene_summary_te_affected_germline.csv")
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
gene_summary_kics <- te_kics_split %>%
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

gene_summary_kics_file <- paste0(r_dir_files, "gene_summary_te_kics_germline.csv")
write.csv(gene_summary_kics, gene_summary_kics_file, row.names = FALSE)
cat(sprintf("✓ Gene summary for KICS saved to: %s\n", basename(gene_summary_kics_file)))

# Generate gene summary for LFS samples
cat("Generating gene summary table for LFS samples...\n")
gene_summary_lfs <- te_lfs_split %>%
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

gene_summary_lfs_file <- paste0(r_dir_files, "gene_summary_te_lfs_germline.csv")
write.csv(gene_summary_lfs, gene_summary_lfs_file, row.names = FALSE)
cat(sprintf("✓ Gene summary for LFS saved to: %s\n", basename(gene_summary_lfs_file)))

# Generate cancer genes only gene summary for all samples
cat("Generating cancer genes only gene summary for all samples...\n")
# Define cancer predisposition genes list
genes <- cpg$V1
gene_summary_cancer <- gene_summary %>%
  filter(Gene_name %in% genes)

gene_summary_cancer_file <- paste0(r_dir_files, "gene_summary_te_all_cancer_genes_germline.csv")
write.csv(gene_summary_cancer, gene_summary_cancer_file, row.names = FALSE)
cat(sprintf("✓ Cancer genes only gene summary saved to: %s\n", basename(gene_summary_cancer_file)))
cat(sprintf("  Total cancer genes with TEs: %d\n\n", nrow(gene_summary_cancer)))

# Generate cancer genes only gene summary for KICS samples
cat("Generating cancer genes only gene summary for KICS samples...\n")
gene_summary_kics_cancer <- gene_summary_kics %>%
  filter(Gene_name %in% genes)

gene_summary_kics_cancer_file <- paste0(r_dir_files, "gene_summary_te_kics_cancer_genes_germline.csv")
write.csv(gene_summary_kics_cancer, gene_summary_kics_cancer_file, row.names = FALSE)
cat(sprintf("✓ Cancer genes only gene summary for KICS saved to: %s\n", basename(gene_summary_kics_cancer_file)))
cat(sprintf("  Total cancer genes with TEs (KICS): %d\n\n", nrow(gene_summary_kics_cancer)))

# Generate cancer genes only gene summary for LFS samples
cat("Generating cancer genes only gene summary for LFS samples...\n")
gene_summary_lfs_cancer <- gene_summary_lfs %>%
  filter(Gene_name %in% genes)

gene_summary_lfs_cancer_file <- paste0(r_dir_files, "gene_summary_te_lfs_cancer_genes_germline.csv")
write.csv(gene_summary_lfs_cancer, gene_summary_lfs_cancer_file, row.names = FALSE)
cat(sprintf("✓ Cancer genes only gene summary for LFS saved to: %s\n", basename(gene_summary_lfs_cancer_file)))
cat(sprintf("  Total cancer genes with TEs (LFS): %d\n\n", nrow(gene_summary_lfs_cancer)))

# Systematic TE expression analysis - test all TEs for significant expression effects (by coordinates)
systematic_results <- run_systematic_te_expression_analysis(
  rna_filtered = rna_filtered,
  te_split_df = te_aff_split,
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
  te_split_df = te_aff_split,
  plot_dir = plot_dir,
  min_samples = 3,  # Lower threshold for discovery
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
  te_split_df = te_aff_split,
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
  te_split_df = te_aff_split,
  plot_dir = plot_dir,
  min_samples = 3,  # Lower threshold for discovery
  fdr_cutoff = 0.1,  # More lenient for discovery
  create_plots = TRUE,
  max_plots = 3,
  group_by_gene = TRUE,
  files_dir = r_dir_files
)

write_output(quote(systematic_results_by_gene_all), "Comprehensive Systematic TE Expression Analysis Results (grouped by gene)")

# Systematic TE expression analysis with min_samples = 10 (by coordinates)
systematic_results_10 <- run_systematic_te_expression_analysis(
  rna_filtered = rna_filtered,
  te_split_df = te_aff_split,
  plot_dir = plot_dir,
  min_samples = 10,
  fdr_cutoff = 0.05,
  create_plots = TRUE,
  max_plots = 3,
  group_by_gene = FALSE,
  files_dir = r_dir_files
)

write_output(quote(systematic_results_10), "Systematic TE Expression Analysis Results (by coordinates, min_samples=10)")

# Systematic TE expression analysis with min_samples = 15 (by coordinates)
systematic_results_15 <- run_systematic_te_expression_analysis(
  rna_filtered = rna_filtered,
  te_split_df = te_aff_split,
  plot_dir = plot_dir,
  min_samples = 15,
  fdr_cutoff = 0.05,
  create_plots = TRUE,
  max_plots = 3,
  group_by_gene = FALSE,
  files_dir = r_dir_files
)

write_output(quote(systematic_results_15), "Systematic TE Expression Analysis Results (by coordinates, min_samples=15)")

# Systematic TE expression analysis with min_samples = 20 (by coordinates)
systematic_results_20 <- run_systematic_te_expression_analysis(
  rna_filtered = rna_filtered,
  te_split_df = te_aff_split,
  plot_dir = plot_dir,
  min_samples = 20,
  fdr_cutoff = 0.05,
  create_plots = TRUE,
  max_plots = 3,
  group_by_gene = FALSE,
  files_dir = r_dir_files
)

write_output(quote(systematic_results_20), "Systematic TE Expression Analysis Results (by coordinates, min_samples=20)")

# Systematic TE expression analysis with min_samples = 10 (by gene)
systematic_results_by_gene_10 <- run_systematic_te_expression_analysis(
  rna_filtered = rna_filtered,
  te_split_df = te_aff_split,
  plot_dir = plot_dir,
  min_samples = 10,
  fdr_cutoff = 0.05,
  create_plots = TRUE,
  max_plots = 3,
  group_by_gene = TRUE,
  files_dir = r_dir_files
)

write_output(quote(systematic_results_by_gene_10), "Systematic TE Expression Analysis Results (grouped by gene, min_samples=10)")

# Systematic TE expression analysis with min_samples = 15 (by gene)
systematic_results_by_gene_15 <- run_systematic_te_expression_analysis(
  rna_filtered = rna_filtered,
  te_split_df = te_aff_split,
  plot_dir = plot_dir,
  min_samples = 15,
  fdr_cutoff = 0.05,
  create_plots = TRUE,
  max_plots = 3,
  group_by_gene = TRUE,
  files_dir = r_dir_files
)

write_output(quote(systematic_results_by_gene_15), "Systematic TE Expression Analysis Results (grouped by gene, min_samples=15)")

# Systematic TE expression analysis with min_samples = 20 (by gene)
systematic_results_by_gene_20 <- run_systematic_te_expression_analysis(
  rna_filtered = rna_filtered,
  te_split_df = te_aff_split,
  plot_dir = plot_dir,
  min_samples = 20,
  fdr_cutoff = 0.05,
  create_plots = TRUE,
  max_plots = 3,
  group_by_gene = TRUE,
  files_dir = r_dir_files
)

write_output(quote(systematic_results_by_gene_20), "Systematic TE Expression Analysis Results (grouped by gene, min_samples=20)")


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
    } else if (params$min_samples == 10) {
      systematic_results_to_use <- systematic_results_10
    } else if (params$min_samples == 15) {
      systematic_results_to_use <- systematic_results_15
    } else if (params$min_samples == 20) {
      systematic_results_to_use <- systematic_results_20
    }
  } else if (params$analysis_type == "gene") {
    if (params$min_samples == 5) {
      systematic_results_to_use <- systematic_results_by_gene
    } else if (params$min_samples == 10) {
      systematic_results_to_use <- systematic_results_by_gene_10
    } else if (params$min_samples == 15) {
      systematic_results_to_use <- systematic_results_by_gene_15
    } else if (params$min_samples == 20) {
      systematic_results_to_use <- systematic_results_by_gene_20
    }
  }

  # Get genes with expression p < threshold from systematic expression analysis
  if (exists("systematic_results_to_use") && nrow(systematic_results_to_use) > 0 && "p_value" %in% colnames(systematic_results_to_use)) {
    # Get genes with p < threshold directly from systematic results
    significant_genes <- systematic_results_to_use %>%
      filter(p_value < params$p_gene) %>%
      pull(gene) %>%
      unique()

    cat("Genes with expression p-value <", params$p_gene, ":", length(significant_genes), "\n")

    # Save genes with expression effects
    genes_with_expression_full <- systematic_results_to_use %>%
      filter(p_value < params$p_gene) %>%
      select(gene, p_value, p_adj, mean_with_te, mean_without_te, samples_with_te, samples_without_te) %>%
      arrange(p_value)

    genes_file <- paste0(r_dir_files, "te_genes_with_expression_germline", param_suffix, ".csv")
    write.csv(genes_with_expression_full, genes_file, row.names=FALSE)

    # Perform pathway analysis if enough genes
    if (length(significant_genes) >= 5) {
      genes_df <- data.frame(Gene_name = significant_genes)

      ora_expression <- tryCatch({
        perform_ora_custom_cutoffs(genes_df,
                                  p_pathway = params$p_pathway,
                                  q_pathway = params$q_pathway,
                                  nsample_thresh = 0,
                                  filter_exon = FALSE)
      }, error = function(e) {
        cat("Error in pathway analysis:", e$message, "\n")
        NULL
      })

      if (!is.null(ora_expression) && nrow(as.data.frame(ora_expression)) > 0) {
        n_pathways <- nrow(as.data.frame(ora_expression))
        cat("✓ Found", n_pathways, "enriched pathways\n")

        # Save pathway results
        pathway_file <- paste0(r_dir_files, "rna_pathway_germline", param_suffix, ".csv")
        write.csv(as.data.frame(ora_expression), pathway_file, row.names=FALSE)

        # Create plots
        # Bar (single group)
        tryCatch({
          p_bar_ora_expression <- barplot(ora_expression, showCategory=20)
          titled_print(p_bar_ora_expression, "ORA barplot (Expression TEs)")
          ggsave(paste0(plot_dir, "pathway/rna_pathway_germline_bar", param_suffix, ".png"), plot=p_bar_ora_expression, width=14, height=9)
        }, error = function(e) cat("Warning: Could not create barplot:", e$message, "\n"))

        # Dot
        tryCatch({
          p_dot_ora_expression <- dotplot(ora_expression, showCategory=20)
          titled_print(p_dot_ora_expression, "ORA dotplot (Expression TEs)")
          ggsave(paste0(plot_dir, "pathway/rna_pathway_germline_dot", param_suffix, ".png"), plot=p_dot_ora_expression, width=14, height=9)
        }, error = function(e) cat("Warning: Could not create dotplot:", e$message, "\n"))

        # Cnet
        tryCatch({
          p_cnet_ora_expression <- cnetplot(ora_expression, showCategory=10, colorEdge=TRUE, node_label="category")
          titled_print(p_cnet_ora_expression, "ORA cnetplot (Expression TEs)")
            ggsave(paste0(plot_dir, "pathway/rna_pathway_germline_cnet_", param_suffix, ".png"), plot=p_cnet_ora_expression, width=14, height=9)
          }, error = function(e) cat("Warning: Could not create cnetplot:", e$message, "\n"))

          # Emap
          tryCatch({
            ora_expression_pairwise <- pairwise_termsim(ora_expression)
            p_emap_ora_expression <- emapplot(ora_expression_pairwise, pie="count", showCategory=20, group_category=TRUE, group_legend=TRUE)
            titled_print(p_emap_ora_expression, "ORA emapplot (Expression TEs)")
            ggsave(paste0(plot_dir, "pathway/rna_pathway_germline_emap_", param_suffix, ".png"), plot=p_emap_ora_expression, width=14, height=9)
          }, error = function(e) cat("Warning: Could not create emapplot:", e$message, "\n"))

          # Simplified pathway plots
          cat("\nCreating simplified pathway plots...\n")
          ora_expression_simple <- simplify(ora_expression, cutoff=0.5, by="p.adjust", select_fun=min)
          write.csv(as.data.frame(ora_expression_simple), paste0(r_dir_files, "rna_pathway_germline_simple.csv"), row.names=FALSE)

          if (!is.null(ora_expression_simple) && nrow(as.data.frame(ora_expression_simple)) > 0) {
            # Bar (single group)
            tryCatch({
              p_bar_ora_expression_simple <- barplot(ora_expression_simple, showCategory=20)
              titled_print(p_bar_ora_expression_simple, "ORA barplot simplified (Expression TEs)")
              ggsave(paste0(plot_dir, "pathway/rna_pathway_germline_bar_simplified.png"), plot=p_bar_ora_expression_simple, width=14, height=9)
            }, error = function(e) cat("Warning: Could not create simplified barplot:", e$message, "\n"))

            # Dot
            tryCatch({
              p_dot_ora_expression_simple <- dotplot(ora_expression_simple, showCategory=20)
              titled_print(p_dot_ora_expression_simple, "ORA dotplot simplified (Expression TEs)")
              ggsave(paste0(plot_dir, "pathway/rna_pathway_germline_dot_simplified.png"), plot=p_dot_ora_expression_simple, width=14, height=9)
            }, error = function(e) cat("Warning: Could not create simplified dotplot:", e$message, "\n"))

            # Cnet
            tryCatch({
              p_cnet_ora_expression_simple <- cnetplot(ora_expression_simple, showCategory=10, colorEdge=TRUE, node_label="category")
              titled_print(p_cnet_ora_expression_simple, "ORA cnetplot simplified (Expression TEs)")
              ggsave(paste0(plot_dir, "pathway/rna_pathway_germline_cnet_simplified.png"), plot=p_cnet_ora_expression_simple, width=14, height=9)
            }, error = function(e) cat("Warning: Could not create simplified cnetplot:", e$message, "\n"))

            # Emap
            tryCatch({
              ora_expression_simple_pairwise <- pairwise_termsim(ora_expression_simple)
              p_emap_ora_expression_simple <- emapplot(ora_expression_simple_pairwise, pie="count", showCategory=20, group_category=TRUE, group_legend=TRUE)
              titled_print(p_emap_ora_expression_simple, "ORA emapplot simplified (Expression TEs)")
              ggsave(paste0(plot_dir, "pathway/rna_pathway_germline_emap_simplified.png"), plot=p_emap_ora_expression_simple, width=14, height=9)
            }, error = function(e) cat("Warning: Could not create simplified emapplot:", e$message, "\n"))
          } else {
            cat("No significant pathways in simplified results\n")
          }

          tryCatch({
            pathway_genes_expression <- plot_pathway_gene_counts(ora_expression, n_descriptions=nrow(ora_expression), n_genes=50)
            write.csv(pathway_genes_expression, paste0(r_dir_files, "te_pathway_genes_expression_germline.csv"), row.names=FALSE)
          }, error = function(e) cat("Warning: Could not create pathway gene counts:", e$message, "\n"))

        # Record sweep summary
        sweep_summary <- rbind(sweep_summary, data.frame(
          min_samples = params$min_samples,
          analysis_type = params$analysis_type,
          p_gene = params$p_gene,
          p_pathway = params$p_pathway,
          q_pathway = params$q_pathway,
          n_genes = length(significant_genes),
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
          n_genes = length(significant_genes),
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
        n_genes = length(significant_genes),
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
  summary_file <- paste0(r_dir_files, "rna_pathway_sweep_germline_summary.csv")
  write.csv(sweep_summary, summary_file, row.names = FALSE)
  cat("\n✓ Parameter sweep complete! Summary saved to:", summary_file, "\n")
  cat("\nTop 5 parameter combinations by number of pathways:\n")
  print(head(sweep_summary[order(-sweep_summary$n_pathways), ], 5))
}

#### RNA PATHWAY ANALYSIS BY TP53 STATUS ####
cat("\n\n=== RNA Pathway Analysis: Affected Cohort (by TP53 Status) ===\n")

# Check if TP53_status column exists
if ("TP53_status" %in% colnames(te_aff_split)) {
  cat("TP53_status groups:\n")
  print(table(te_aff_split$TP53_status))

  # Get genes with TEs, grouped by TP53 status
  cat("\nGrouping genes by TP53 status for pathway analysis...\n")
  geneClusters_tp53 <- lapply(split(te_aff_split$Gene_name, te_aff_split$TP53_status), unique)

  cat("Genes per TP53 status:\n")
  n_genes_per_group <- sapply(geneClusters_tp53, length)
  for (status in names(geneClusters_tp53)) {
    cat("  ", status, ":", length(geneClusters_tp53[[status]]), "genes\n")
  }
  n_samples_rna_tp53 <- length(unique(te_aff_split$sample))
  n_genes_rna_tp53 <- length(unique(te_aff_split$Gene_name))

  # Count samples per TP53 group BEFORE running analysis
  samples_per_tp53_rna <- table(unique(te_aff_split[, c("sample", "TP53_status")])$TP53_status)
  samples_per_tp53_rna_str <- paste(names(samples_per_tp53_rna), "=", samples_per_tp53_rna, "samples", collapse=", ")
  genes_per_group_str <- paste(names(n_genes_per_group), "=", n_genes_per_group, "genes", collapse=", ")

  # Write data summary before analysis
  summary_text <- paste0("RNA Pathway Analysis Summary - TP53 Status\n",
                        "===========================================\n\n",
                        "Total samples: ", n_samples_rna_tp53, "\n",
                        "Samples per group: ", samples_per_tp53_rna_str, "\n",
                        "Total genes: ", n_genes_rna_tp53, "\n",
                        "Genes per group: ", genes_per_group_str, "\n",
                        "Parameters: p<0.05, q<0.1\n\n",
                        "Status: Running analysis...")
  writeLines(summary_text, paste0(plot_dir, "rna/rna_pathway_tp53_summary.txt"))

  cat("\nRunning pathway analysis by TP53 status...\n")
  ora_rna_tp53 <- tryCatch({
    compareCluster(
      geneCluster = geneClusters_tp53,
      fun = "enrichGO",
      OrgDb = org.Hs.eg.db,
      keyType = "SYMBOL",
      ont = "BP",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.1
    )
  }, error = function(e) {
    cat("Error in pathway analysis:", e$message, "\n")
    NULL
  })

  if (!is.null(ora_rna_tp53) && nrow(as.data.frame(ora_rna_tp53)) > 0) {
    cat("Pathway analysis successful! Found", nrow(as.data.frame(ora_rna_tp53)), "significant pathways\n")
    write.csv(as.data.frame(ora_rna_tp53), paste0(r_dir_files, "rna_pathway_tp53.csv"), row.names=FALSE)

    # Update summary with results
    summary_text <- paste0("RNA Pathway Analysis Summary - TP53 Status\n",
                          "===========================================\n\n",
                          "Total samples: ", n_samples_rna_tp53, "\n",
                          "Samples per group: ", samples_per_tp53_rna_str, "\n",
                          "Total genes: ", n_genes_rna_tp53, "\n",
                          "Genes per group: ", genes_per_group_str, "\n",
                          "Parameters: p<0.05, q<0.1\n\n",
                          "Status: SUCCESS\n",
                          "Significant pathways found: ", nrow(as.data.frame(ora_rna_tp53)))
    writeLines(summary_text, paste0(plot_dir, "rna/rna_pathway_tp53_summary.txt"))

    # Dot (no bar for compareCluster)
    cat("Creating dotplot...\n")
    p_dot_ora_rna_tp53 <- dotplot(ora_rna_tp53, showCategory=20)
    titled_print(p_dot_ora_rna_tp53, "ORA dotplot (RNA - TP53 Status)")
    ggsave(paste0(plot_dir, "pathway/rna_pathway_tp53_dot.png"), plot=p_dot_ora_rna_tp53, width=14, height=9)

    # Cnet
    cat("Creating cnetplot...\n")
    p_cnet_ora_rna_tp53 <- cnetplot(ora_rna_tp53, showCategory=10, colorEdge=TRUE, node_label="category")
    titled_print(p_cnet_ora_rna_tp53, "RNA ORA cnetplot (TP53 Status)")
    ggsave(paste0(plot_dir, "pathway/rna_pathway_tp53_cnet.png"), plot=p_cnet_ora_rna_tp53, width=14, height=9)

    # Emap
    cat("Creating emapplot...\n")
    tryCatch({
      # Check if multiple clusters have results
      cluster_counts <- table(as.data.frame(ora_rna_tp53)$Cluster)
      if (length(cluster_counts) > 1) {
        ora_rna_tp53_pairwise <- pairwise_termsim(ora_rna_tp53)
        p_emap_ora_rna_tp53 <- emapplot(ora_rna_tp53_pairwise, showCategory=20,
                                        pie.params = list(pie = "count"),
                                        cluster.params = list(cluster = TRUE, legend = TRUE))
        titled_print(p_emap_ora_rna_tp53, "ORA emapplot (RNA - TP53)")
        ggsave(paste0(plot_dir, "pathway/rna_pathway_tp53_emap.png"), plot=p_emap_ora_rna_tp53, width=14, height=9)
      } else {
        cat("Skipping emapplot: Only one cluster has results (", names(cluster_counts), ")\n")
      }
    }, error = function(e) {
      cat("Warning: Could not create emapplot:", e$message, "\n")
    })

    # Simplified pathway plots
    cat("\nCreating simplified pathway plots...\n")
    ora_rna_tp53_simple <- simplify(ora_rna_tp53, cutoff=0.5, by="p.adjust", select_fun=min)
    write.csv(as.data.frame(ora_rna_tp53_simple), paste0(r_dir_files, "rna_pathway_tp53_simple.csv"), row.names=FALSE)

    if (!is.null(ora_rna_tp53_simple) && nrow(as.data.frame(ora_rna_tp53_simple)) > 0) {
      # Dot (no bar for compareCluster)
      p_dot_ora_rna_tp53_simple <- dotplot(ora_rna_tp53_simple, showCategory=20)
      titled_print(p_dot_ora_rna_tp53_simple, "ORA dotplot simplified (RNA - TP53 Status)")
      ggsave(paste0(plot_dir, "pathway/rna_pathway_tp53_dot_simplified.png"), plot=p_dot_ora_rna_tp53_simple, width=14, height=9)

      # Cnet
      p_cnet_ora_rna_tp53_simple <- cnetplot(ora_rna_tp53_simple, showCategory=10, colorEdge=TRUE, node_label="category")
      titled_print(p_cnet_ora_rna_tp53_simple, "RNA ORA cnetplot simplified (TP53 Status)")
      ggsave(paste0(plot_dir, "pathway/rna_pathway_tp53_cnet_simplified.png"), plot=p_cnet_ora_rna_tp53_simple, width=14, height=9)

      # Emap
      tryCatch({
        cluster_counts_simple <- table(as.data.frame(ora_rna_tp53_simple)$Cluster)
        if (length(cluster_counts_simple) > 1) {
          ora_rna_tp53_simple_pairwise <- pairwise_termsim(ora_rna_tp53_simple)
          p_emap_ora_rna_tp53_simple <- emapplot(ora_rna_tp53_simple_pairwise, showCategory=20,
                                                  pie.params = list(pie = "count"),
                                                  cluster.params = list(cluster = TRUE, legend = TRUE))
          titled_print(p_emap_ora_rna_tp53_simple, "ORA emapplot simplified (RNA - TP53)")
          ggsave(paste0(plot_dir, "pathway/rna_pathway_tp53_emap_simplified.png"), plot=p_emap_ora_rna_tp53_simple, width=14, height=9)
        } else {
          cat("Skipping simplified emapplot: Only one cluster has results\n")
        }
      }, error = function(e) {
        cat("Warning: Could not create simplified emapplot:", e$message, "\n")
      })
    } else {
      cat("No significant pathways in simplified results\n")
    }
  } else {
    # Update summary with failure status
    if (is.null(ora_rna_tp53)) {
      status_msg <- paste0("Status: INSUFFICIENT DATA\n",
                          "Reason: Not enough data to run pathway analysis")
    } else {
      status_msg <- paste0("Status: NO SIGNIFICANT PATHWAYS\n",
                          "Analysis completed but no pathways met significance thresholds")
    }

    summary_text <- paste0("RNA Pathway Analysis Summary - TP53 Status\n",
                          "===========================================\n\n",
                          "Total samples: ", n_samples_rna_tp53, "\n",
                          "Samples per group: ", samples_per_tp53_rna_str, "\n",
                          "Total genes: ", n_genes_rna_tp53, "\n",
                          "Genes per group: ", genes_per_group_str, "\n",
                          "Parameters: p<0.05, q<0.1\n\n",
                          status_msg)
    writeLines(summary_text, paste0(plot_dir, "rna/rna_pathway_tp53_summary.txt"))

    cat(status_msg, "\n")
  }
} else {
  cat("WARNING: TP53_status column not found in te_aff_split, skipping TP53 analysis\n")
}

# Save rna_filtered for use by RE script
cat("\nSaving rna_filtered for RE script...\n")
saveRDS(rna_filtered, paste0(r_dir_files, "rna_filtered_germline.rds"))
cat("✓ Saved rna_filtered to:", paste0(r_dir_files, "rna_filtered_germline.rds\n"))

# Close PDF device after all plots
dev.off()

cat("\n===== SCRIPT COMPLETED SUCCESSFULLY =====\n")
cat("All analysis sections have been processed.\n")
cat("Generated plots saved to:", plot_dir, "\n")
cat("PDF compilation saved to: graph_output.pdf\n")
cat("Text output saved to:", stdout_file, "\n")

# Close sink to stop redirecting output
sink()

cat("✓ Script completed successfully\n")
