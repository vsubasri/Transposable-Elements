#!/usr/bin/env Rscript

# Germline TE Visualization - Window Analysis
# Analyze TE insertions within genomic windows from ML feature selection

# Source common setup and load data
source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/00_viz_common_setup.R")
REQUIRED_DATA <- c("expand", "split", "clinical")
source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/00_viz_load_data_germline.R")

cat("Running 02_te_viz_germline_11_window_analysis.R...\n")

#### ANALYSIS PARAMETERS - MODIFY THESE ####
ML_FEATURE_FILE <- "/Users/briannelaverty/Documents/R_Malkin/te/output/ml_analysis/100kb_TP53_status_te_aff/top20_genomic_windows_log_100kb_TP53_status_te_aff_expand_g.csv"
DATASET <- "te_aff"           # Dataset: te_aff, te_lfs, te_kics, te_all, etc.
GROUP <- "TP53_status"        # Grouping column: TP53_status, Cancer_status, etc.
SPLIT <- "expand"             # Format: expand, split, or count
OUTPUT_PREFIX <- "window_analysis"  # Prefix for output files

cat("\n===== CONFIGURATION =====\n")
cat("ML Feature File:", ML_FEATURE_FILE, "\n")
cat("Dataset:", DATASET, "\n")
cat("Group:", GROUP, "\n")
cat("Split:", SPLIT, "\n")
cat("Output Prefix:", OUTPUT_PREFIX, "\n\n")

#### LOAD WINDOW DATA ####
cat("Loading window data...\n")
if (!file.exists(ML_FEATURE_FILE)) {
  stop("ML feature file not found: ", ML_FEATURE_FILE)
}
windows <- read.csv(ML_FEATURE_FILE, stringsAsFactors = FALSE)

# Ensure chr column is character to match TE data format
windows$chr <- as.character(windows$chr)

cat("Loaded", nrow(windows), "windows\n")
cat("Window columns:", paste(colnames(windows), collapse = ", "), "\n")
cat("Sample window chr values:", paste(head(windows$chr, 3), collapse = ", "), "\n\n")

# Check required columns
required_cols <- c("id", "chr", "start", "end", "type")
missing_cols <- setdiff(required_cols, colnames(windows))
if (length(missing_cols) > 0) {
  stop("Missing required columns in window file: ", paste(missing_cols, collapse = ", "))
}

#### SELECT DATASET ####
cat("Selecting dataset...\n")
dataset_name <- paste0(DATASET, "_", SPLIT)
if (!exists(dataset_name)) {
  stop("Dataset not found: ", dataset_name, ". Available datasets: ",
       paste(grep("^te_", ls(), value = TRUE), collapse = ", "))
}
te_data <- get(dataset_name)
cat("Using dataset:", dataset_name, "with", nrow(te_data), "rows\n")
cat("Sample TE SV_chrom values:", paste(head(unique(te_data$SV_chrom), 5), collapse = ", "), "\n")
cat("Sample TE ALT values:", paste(unique(te_data$ALT), collapse = ", "), "\n")

# Check for required columns in TE data
required_te_cols <- c("sample", "SV_chrom", "SV_start", "SV_end", "ALT", "ID", "Gene_name")
missing_te_cols <- setdiff(required_te_cols, colnames(te_data))
if (length(missing_te_cols) > 0) {
  stop("Missing required columns in TE data: ", paste(missing_te_cols, collapse = ", "))
}

# Check for grouping column
if (!GROUP %in% colnames(te_data)) {
  stop("Grouping column '", GROUP, "' not found in dataset. Available columns: ",
       paste(colnames(te_data), collapse = ", "))
}

# Identify groups
groups <- sort(unique(te_data[[GROUP]]))
groups <- groups[!is.na(groups)]
if (length(groups) != 2) {
  warning("Expected 2 groups for ", GROUP, ", found ", length(groups), ": ",
          paste(groups, collapse = ", "))
}
group1 <- groups[1]
group2 <- groups[2]
cat("Groups:", group1, "vs", group2, "\n\n")

#### LOAD RE DATA ####
cat("Loading regulatory element data...\n")
re_report_path <- "/Users/briannelaverty/Documents/R_Malkin/te/data/final/germline_annotSV_output.SV_RE_intersect.report"
if (file.exists(re_report_path)) {
  tryCatch({
    # Load RE report manually to use LEFT JOIN (keeps all TEs)
    re_data <- fread(re_report_path, header = FALSE)
    colnames(re_data) <- c("chr", "start", "end", "ins", "sample", "ID", "ref", "ALT_re", "x",
                           "filter", "info", "y", "genotype", "chr_reg", "start_reg", "end_reg",
                           "type_reg", "gene_reg")

    cat("Loaded", nrow(re_data), "RE records\n")

    # LEFT JOIN to keep all TEs (unlike load_and_join_re_data which uses inner_join)
    te_data_with_re <- te_data %>%
      left_join(re_data %>% select(ID, gene_reg), by = "ID", multiple = "all")

    cat("After LEFT join:", nrow(te_data_with_re), "rows (kept all TEs)\n")
  }, error = function(e) {
    cat("Warning: Could not load RE data:", e$message, "\n")
    cat("Proceeding without RE gene information\n")
    te_data_with_re <- te_data
    te_data_with_re$gene_reg <- NA
  })
} else {
  cat("RE report file not found, proceeding without RE gene information\n")
  te_data_with_re <- te_data
  te_data_with_re$gene_reg <- NA
}
cat("\n")

#### PROCESS WINDOWS ####
cat("Processing windows...\n")

# Initialize output data frames
te_details_list <- list()
window_summary_list <- list()
gene_summary_list <- list()

for (i in 1:nrow(windows)) {
  window <- windows[i, ]
  window_id <- window$id
  window_chr <- window$chr
  window_start <- window$start
  window_end <- window$end
  window_type <- window$type
  window_importance <- ifelse("importance" %in% colnames(window), window$importance, NA)

  if (i %% 5 == 0) {
    cat("  Processing window", i, "of", nrow(windows), "\n")
  }

  # Debug first window
  if (i == 1) {
    cat("\nDEBUG - First window:\n")
    cat("  window_chr:", window_chr, "(class:", class(window_chr), ")\n")
    cat("  window_start:", window_start, "\n")
    cat("  window_end:", window_end, "\n")
    cat("  window_type:", window_type, "\n")
    cat("  Filtering for: SV_chrom ==", window_chr, "& ALT ==", window_type, "\n")

    # Count TEs on this chromosome
    chr_tes <- te_data_with_re %>% filter(SV_chrom == window_chr)
    cat("  TEs on chr", window_chr, ":", nrow(chr_tes), "\n")

    # Count TEs of this type on this chromosome
    type_tes <- chr_tes %>% filter(ALT == window_type)
    cat("  TEs of type", window_type, "on chr", window_chr, ":", nrow(type_tes), "\n")
  }

  # Filter TEs by THREE criteria (matching working script logic):
  # 1. Chromosome match
  # 2. Type match
  # 3. 50% overlap (bedtools -F 0.50 logic)
  matching_tes <- te_data_with_re %>%
    filter(
      SV_chrom == window_chr,
      ALT == window_type
    ) %>%
    mutate(
      overlap_start = pmax(SV_start, window_start),
      overlap_end = pmin(SV_end, window_end),
      overlap_length = pmax(0, overlap_end - overlap_start),
      te_length = SV_end - SV_start,
      overlap_frac = overlap_length / te_length
    ) %>%
    filter(overlap_frac >= 0.50)

  if (i == 1) {
    cat("  Matching TEs in window (≥50% overlap):", nrow(matching_tes), "\n\n")
  }

  if (nrow(matching_tes) == 0) {
    # No TEs in this window - add empty row to summary
    window_summary_list[[i]] <- data.frame(
      window_id = window_id,
      chr = window_chr,
      start = window_start,
      end = window_end,
      type = window_type,
      importance = window_importance,
      n_samples_group1 = 0,
      n_samples_group2 = 0,
      n_unique_tes = 0,
      n_unique_genes = 0,
      n_unique_re_genes = 0,
      samples_group1 = "",
      samples_group2 = "",
      genes_affected = "",
      re_genes = "",
      stringsAsFactors = FALSE
    )
    next
  }

  # Add group information
  matching_tes$group_value <- matching_tes[[GROUP]]

  # Create unique TE identifier (coordinate + type)
  matching_tes$te_coord <- paste0("chr", matching_tes$SV_chrom, ":",
                                   matching_tes$SV_start, "-", matching_tes$SV_end)

  # Process by unique TE coordinate
  unique_tes <- matching_tes %>%
    group_by(te_coord, SV_chrom, SV_start, SV_end, ALT) %>%
    summarise(
      n_samples_total = n_distinct(sample),
      n_group1 = sum(group_value == group1, na.rm = TRUE),
      n_group2 = sum(group_value == group2, na.rm = TRUE),
      samples_group1 = paste(unique(sample[group_value == group1]), collapse = ";"),
      samples_group2 = paste(unique(sample[group_value == group2]), collapse = ";"),
      genes_raw = paste(unique(Gene_name), collapse = ";"),
      re_genes_raw = paste(unique(na.omit(gene_reg)), collapse = ";"),
      .groups = "drop"
    )

  # Split and clean gene lists
  unique_tes$genes_from_this_te <- sapply(unique_tes$genes_raw, function(x) {
    if (is.na(x) || x == "") return("")
    genes <- unlist(strsplit(x, "[;/]"))
    genes <- unique(trimws(genes))
    genes <- genes[genes != "" & !is.na(genes)]
    paste(genes, collapse = ";")
  })

  unique_tes$re_genes_from_this_te <- sapply(unique_tes$re_genes_raw, function(x) {
    if (is.na(x) || x == "") return("")
    genes <- unlist(strsplit(x, "[;/]"))
    genes <- unique(trimws(genes))
    genes <- genes[genes != "" & !is.na(genes)]
    paste(genes, collapse = ";")
  })

  # Build TE details output for this window
  te_details_list[[i]] <- data.frame(
    window_id = window_id,
    window_chr = window_chr,
    window_start = window_start,
    window_end = window_end,
    window_type = window_type,
    importance = window_importance,
    te_coord = unique_tes$te_coord,
    te_chr = unique_tes$SV_chrom,
    te_start = unique_tes$SV_start,
    te_end = unique_tes$SV_end,
    te_type = unique_tes$ALT,
    n_samples_total = unique_tes$n_samples_total,
    n_group1 = unique_tes$n_group1,
    n_group2 = unique_tes$n_group2,
    samples_group1 = unique_tes$samples_group1,
    samples_group2 = unique_tes$samples_group2,
    genes_from_this_te = unique_tes$genes_from_this_te,
    re_genes_from_this_te = unique_tes$re_genes_from_this_te,
    stringsAsFactors = FALSE
  )

  # Build window summary
  all_samples_group1 <- unique(matching_tes$sample[matching_tes$group_value == group1])
  all_samples_group2 <- unique(matching_tes$sample[matching_tes$group_value == group2])

  all_genes <- unique(unlist(strsplit(unique_tes$genes_from_this_te, ";")))
  all_genes <- all_genes[all_genes != "" & !is.na(all_genes)]

  all_re_genes <- unique(unlist(strsplit(unique_tes$re_genes_from_this_te, ";")))
  all_re_genes <- all_re_genes[all_re_genes != "" & !is.na(all_re_genes)]

  window_summary_list[[i]] <- data.frame(
    window_id = window_id,
    chr = window_chr,
    start = window_start,
    end = window_end,
    type = window_type,
    importance = window_importance,
    n_samples_group1 = length(all_samples_group1),
    n_samples_group2 = length(all_samples_group2),
    n_unique_tes = nrow(unique_tes),
    n_unique_genes = length(all_genes),
    n_unique_re_genes = length(all_re_genes),
    samples_group1 = paste(all_samples_group1, collapse = ";"),
    samples_group2 = paste(all_samples_group2, collapse = ";"),
    genes_affected = paste(all_genes, collapse = ";"),
    re_genes = paste(all_re_genes, collapse = ";"),
    stringsAsFactors = FALSE
  )

  # Build gene summary (both Gene_name and RE genes)
  gene_summary_rows <- list()

  # Process Gene_name genes
  if (length(all_genes) > 0) {
    for (gene in all_genes) {
      # Find which TEs affect this gene
      tes_with_gene <- unique_tes[grepl(gene, unique_tes$genes_from_this_te, fixed = TRUE), ]

      samples_with_gene_g1 <- unique(unlist(strsplit(tes_with_gene$samples_group1, ";")))
      samples_with_gene_g1 <- samples_with_gene_g1[samples_with_gene_g1 != ""]

      samples_with_gene_g2 <- unique(unlist(strsplit(tes_with_gene$samples_group2, ";")))
      samples_with_gene_g2 <- samples_with_gene_g2[samples_with_gene_g2 != ""]

      gene_summary_rows[[length(gene_summary_rows) + 1]] <- data.frame(
        window_id = window_id,
        window_type = window_type,
        gene = gene,
        gene_source = "Gene_name",
        n_tes_affecting_gene = nrow(tes_with_gene),
        n_samples_total = length(samples_with_gene_g1) + length(samples_with_gene_g2),
        n_group1 = length(samples_with_gene_g1),
        n_group2 = length(samples_with_gene_g2),
        te_coords_affecting_gene = paste(tes_with_gene$te_coord, collapse = ";"),
        samples_group1 = paste(samples_with_gene_g1, collapse = ";"),
        samples_group2 = paste(samples_with_gene_g2, collapse = ";"),
        stringsAsFactors = FALSE
      )
    }
  }

  # Process RE genes
  if (length(all_re_genes) > 0) {
    for (gene in all_re_genes) {
      # Find which TEs affect this gene
      tes_with_gene <- unique_tes[grepl(gene, unique_tes$re_genes_from_this_te, fixed = TRUE), ]

      samples_with_gene_g1 <- unique(unlist(strsplit(tes_with_gene$samples_group1, ";")))
      samples_with_gene_g1 <- samples_with_gene_g1[samples_with_gene_g1 != ""]

      samples_with_gene_g2 <- unique(unlist(strsplit(tes_with_gene$samples_group2, ";")))
      samples_with_gene_g2 <- samples_with_gene_g2[samples_with_gene_g2 != ""]

      gene_summary_rows[[length(gene_summary_rows) + 1]] <- data.frame(
        window_id = window_id,
        window_type = window_type,
        gene = gene,
        gene_source = "RE",
        n_tes_affecting_gene = nrow(tes_with_gene),
        n_samples_total = length(samples_with_gene_g1) + length(samples_with_gene_g2),
        n_group1 = length(samples_with_gene_g1),
        n_group2 = length(samples_with_gene_g2),
        te_coords_affecting_gene = paste(tes_with_gene$te_coord, collapse = ";"),
        samples_group1 = paste(samples_with_gene_g1, collapse = ";"),
        samples_group2 = paste(samples_with_gene_g2, collapse = ";"),
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(gene_summary_rows) > 0) {
    gene_summary_list[[i]] <- do.call(rbind, gene_summary_rows)
  }
}

cat("Finished processing windows\n")

# Summary statistics
windows_with_tes <- sum(sapply(te_details_list, function(x) !is.null(x) && nrow(x) > 0))
cat("Windows with TEs:", windows_with_tes, "out of", nrow(windows), "\n\n")

#### COMBINE AND SAVE OUTPUT ####
cat("Combining results...\n")

# Combine TE details
te_details_df <- do.call(rbind, te_details_list[!sapply(te_details_list, is.null)])
cat("TE details:", nrow(te_details_df), "rows\n")

# Combine window summary
window_summary_df <- do.call(rbind, window_summary_list)
cat("Window summary:", nrow(window_summary_df), "rows\n")

# Combine gene summary
gene_summary_df <- do.call(rbind, gene_summary_list[!sapply(gene_summary_list, is.null)])
cat("Gene summary:", nrow(gene_summary_df), "rows\n\n")

# Save output files
output_dir <- r_dir_files
cat("Saving output files to:", output_dir, "\n")

te_details_file <- paste0(output_dir, OUTPUT_PREFIX, "_te_details.csv")
window_summary_file <- paste0(output_dir, OUTPUT_PREFIX, "_summary.csv")
gene_summary_file <- paste0(output_dir, OUTPUT_PREFIX, "_gene_summary.csv")

write.csv(te_details_df, te_details_file, row.names = FALSE)
cat("✓ Saved:", te_details_file, "\n")

write.csv(window_summary_df, window_summary_file, row.names = FALSE)
cat("✓ Saved:", window_summary_file, "\n")

write.csv(gene_summary_df, gene_summary_file, row.names = FALSE)
cat("✓ Saved:", gene_summary_file, "\n")

cat("\n✓ Script completed successfully\n")
