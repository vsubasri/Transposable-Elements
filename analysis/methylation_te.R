# ==============================================================================
# Functions for Methylation-TE Enrichment Analysis
# ==============================================================================

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(GenomicRanges)

# ------------------------------------------------------------------------------
# Download RepeatMasker annotations from UCSC
# ------------------------------------------------------------------------------
download_repeatmasker <- function(genome = "hg19", output_file = NULL) {
  url <- paste0("http://hgdownload.soe.ucsc.edu/goldenPath/", genome, "/database/rmsk.txt.gz")

  if (is.null(output_file)) {
    output_file <- paste0("rmsk_", genome, ".txt.gz")
  }

  if (!file.exists(output_file)) {
    cat("Downloading RepeatMasker annotations from UCSC...\n")
    cat("URL:", url, "\n")
    download.file(url, output_file, mode = "wb")
    cat("Download complete:", output_file, "\n")
  } else {
    cat("RepeatMasker file already exists:", output_file, "\n")
  }

  return(output_file)
}

# ------------------------------------------------------------------------------
# Load and parse RepeatMasker data
# ------------------------------------------------------------------------------
load_repeatmasker <- function(rmsk_file, te_classes = c("LINE", "SINE", "LTR", "DNA", "SVA")) {
  cat("Loading RepeatMasker annotations...\n")

  rmsk <- fread(
    rmsk_file,
    header = FALSE,
    col.names = c("bin", "swScore", "milliDiv", "milliDel", "milliIns",
                  "chr", "start", "end", "genoLeft", "strand",
                  "repName", "repClass", "repFamily",
                  "repStart", "repEnd", "repLeft", "id")
  )

  cat("Total repeat elements loaded:", nrow(rmsk), "\n")

  if (!is.null(te_classes)) {
    rmsk <- rmsk %>%
      filter(repClass %in% te_classes | grepl(paste(te_classes, collapse = "|"), repFamily))
    cat("TEs after filtering:", nrow(rmsk), "\n")
  }

  cat("\nTE class distribution:\n")
  print(table(rmsk$repClass))

  return(rmsk)
}

# ------------------------------------------------------------------------------
# Get probe coordinates from Bioconductor annotation
# ------------------------------------------------------------------------------
get_probe_coordinates <- function(array_type = "450K", genome = "hg19") {
  cat("Getting probe coordinates for", array_type, "array...\n")

  if (array_type == "450K") {
    if (requireNamespace("IlluminaHumanMethylation450kanno.ilmn12.hg19", quietly = TRUE)) {
      library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
      anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
      probe_coords <- data.frame(
        probe_id = rownames(anno),
        chr = anno$chr,
        pos = anno$pos,
        strand = anno$strand,
        stringsAsFactors = FALSE
      )
      return(probe_coords)
    }
  } else if (array_type == "EPIC") {
    if (requireNamespace("IlluminaHumanMethylationEPICanno.ilm10b4.hg19", quietly = TRUE)) {
      library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
      anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
      probe_coords <- data.frame(
        probe_id = rownames(anno),
        chr = anno$chr,
        pos = anno$pos,
        strand = anno$strand,
        stringsAsFactors = FALSE
      )
      return(probe_coords)
    }
  }

  cat("Warning: Could not load annotation package.\n")
  cat("Install with: BiocManager::install('IlluminaHumanMethylation450kanno.ilmn12.hg19')\n")
  return(NULL)
}

# ------------------------------------------------------------------------------
# Find overlaps between probes and TEs
# ------------------------------------------------------------------------------
find_probe_te_overlaps <- function(probe_coords, rmsk, flank = 0) {
  cat("Finding overlaps between probes and TEs...\n")

  probe_gr <- GRanges(
    seqnames = probe_coords$chr,
    ranges = IRanges(start = probe_coords$pos, width = 1),
    probe_id = probe_coords$probe_id
  )

  te_gr <- GRanges(
    seqnames = rmsk$chr,
    ranges = IRanges(start = pmax(1, rmsk$start - flank), end = rmsk$end + flank),
    repName = rmsk$repName,
    repClass = rmsk$repClass,
    repFamily = rmsk$repFamily
  )

  overlaps <- findOverlaps(probe_gr, te_gr)

  cat("Found", length(overlaps), "probe-TE overlaps\n")
  cat("Unique probes in TEs:", length(unique(queryHits(overlaps))), "\n")

  overlap_df <- data.frame(
    probe_id = probe_gr$probe_id[queryHits(overlaps)],
    te_name = te_gr$repName[subjectHits(overlaps)],
    te_class = te_gr$repClass[subjectHits(overlaps)],
    te_family = te_gr$repFamily[subjectHits(overlaps)]
  )

  return(overlap_df)
}

# ------------------------------------------------------------------------------
# Differential methylation analysis using bumphunter (DMR detection)
# ------------------------------------------------------------------------------
run_differential_methylation <- function(meth_data, group_info, probe_coords,
                                          cutoff = 0.1, B = 100,
                                          fdr_thresh = 0.05) {
  cat("\n=== Differential Methylation Analysis (bumphunter) ===\n")

  library(bumphunter)

  groups <- unique(group_info)
  cat("Groups:", paste(groups, collapse = " vs "), "\n")

  # Align samples
  samples <- names(group_info)[names(group_info) %in% colnames(meth_data)]
  meth_matrix <- as.matrix(meth_data[, samples])
  group_info <- group_info[samples]

  cat("Group 1 (", groups[1], "):", sum(group_info == groups[1]), "samples\n")
  cat("Group 2 (", groups[2], "):", sum(group_info == groups[2]), "samples\n")

  # Align probes with coordinates
  common_probes <- intersect(rownames(meth_matrix), probe_coords$probe_id)
  meth_matrix <- meth_matrix[common_probes, ]
  probe_coords <- probe_coords %>% filter(probe_id %in% common_probes)
  probe_coords <- probe_coords[match(rownames(meth_matrix), probe_coords$probe_id), ]

  cat("Probes with coordinates:", nrow(meth_matrix), "\n")

  # Design matrix (group2 vs group1)
  design <- model.matrix(~ factor(group_info, levels = groups))

  # Run bumphunter
  cat("Running bumphunter (B =", B, "permutations)...\n")
  cat("This may take a few minutes...\n")

  dmrs <- bumphunter(
    object = meth_matrix,
    design = design,
    chr = probe_coords$chr,
    pos = probe_coords$pos,
    coef = 2,  # coefficient for group comparison
    cutoff = cutoff,
    B = B,
    type = "Beta"
  )

  dmr_table <- dmrs$table

  if (is.null(dmr_table) || nrow(dmr_table) == 0) {
    cat("No DMRs found.\n")
    dmr_table <- data.frame(
      chr = character(), start = integer(), end = integer(),
      value = numeric(), area = numeric(), cluster = integer(),
      indexStart = integer(), indexEnd = integer(), L = integer(),
      clusterL = integer(), p.value = numeric(), fwer = numeric(),
      p.valueArea = numeric(), fwerArea = numeric()
    )
  } else {
    cat("\nDMRs found:", nrow(dmr_table), "\n")
  }

  # Classify DMRs by direction
  dmr_table <- dmr_table %>%
    mutate(
      status = case_when(
        fwer < fdr_thresh & value > 0 ~ "hypermethylated",
        fwer < fdr_thresh & value < 0 ~ "hypomethylated",
        TRUE ~ "not_significant"
      )
    )

  cat("\nDMR summary:\n")
  print(table(dmr_table$status))

  # Also return probe-level results for compatibility
  probe_results <- data.frame(
    probe_id = rownames(meth_matrix),
    chr = probe_coords$chr,
    pos = probe_coords$pos,
    stringsAsFactors = FALSE
  )

  # Map probes to DMRs
  probe_results$dmr_index <- NA
  probe_results$status <- "not_significant"

  if (nrow(dmr_table) > 0) {
    for (i in 1:nrow(dmr_table)) {
      idx <- dmr_table$indexStart[i]:dmr_table$indexEnd[i]
      probe_results$dmr_index[idx] <- i
      if (dmr_table$status[i] != "not_significant") {
        probe_results$status[idx] <- dmr_table$status[i]
      }
    }
  }

  cat("\nProbes in significant DMRs:\n")
  print(table(probe_results$status))

  return(list(
    dmrs = dmr_table,
    probes = probe_results,
    bumphunter_obj = dmrs
  ))
}

# ------------------------------------------------------------------------------
# TE enrichment analysis using Fisher's exact test
# ------------------------------------------------------------------------------
run_te_enrichment <- function(dmr_results, overlap_df, rmsk) {
  cat("\n=== TE Enrichment Analysis ===\n")

  # Extract probe-level results from bumphunter output
  probe_results <- dmr_results$probes
  dmr_table <- dmr_results$dmrs

  te_probes <- unique(overlap_df$probe_id)

  # Annotate probes with TE info
  probe_results <- probe_results %>%
    mutate(in_te = probe_id %in% te_probes)

  # Rename for compatibility with rest of function
  dmp_results <- probe_results

  te_anno <- overlap_df %>%
    group_by(probe_id) %>%
    summarise(
      te_classes = paste(unique(te_class), collapse = ";"),
      te_families = paste(unique(te_family), collapse = ";"),
      .groups = "drop"
    )

  dmp_results <- dmp_results %>%
    left_join(te_anno, by = "probe_id")

  # Overall enrichment
  cat("\n--- Overall TE Enrichment ---\n")

  mat_overall <- matrix(
    c(
      sum(dmp_results$status != "not_significant" & dmp_results$in_te),
      sum(dmp_results$status != "not_significant" & !dmp_results$in_te),
      sum(dmp_results$status == "not_significant" & dmp_results$in_te),
      sum(dmp_results$status == "not_significant" & !dmp_results$in_te)
    ),
    nrow = 2, byrow = TRUE,
    dimnames = list(c("DMP", "Not_DMP"), c("In_TE", "Not_in_TE"))
  )

  cat("Contingency matrix:\n")
  print(mat_overall)

  fisher_overall <- fisher.test(mat_overall)
  cat("\nFisher's exact test:\n")
  cat("  Odds ratio:", round(fisher_overall$estimate, 3), "\n")
  cat("  P-value:", format(fisher_overall$p.value, scientific = TRUE), "\n")

  # By direction (hypo/hyper)
  cat("\n--- Enrichment by Direction ---\n")

  enrichment_by_direction <- list()

  for (direction in c("hypomethylated", "hypermethylated")) {
    dmps <- dmp_results %>% filter(status == direction)
    bg <- dmp_results %>% filter(status == "not_significant")

    mat <- matrix(
      c(sum(dmps$in_te), sum(!dmps$in_te),
        sum(bg$in_te), sum(!bg$in_te)),
      nrow = 2, byrow = TRUE
    )

    fisher_test <- fisher.test(mat)

    cat("\n", toupper(direction), ":\n")
    cat("  N probes:", nrow(dmps), "| In TEs:", sum(dmps$in_te),
        "(", round(100 * mean(dmps$in_te), 1), "%)\n")
    cat("  OR:", round(fisher_test$estimate, 3),
        "| P:", format(fisher_test$p.value, scientific = TRUE), "\n")

    enrichment_by_direction[[direction]] <- list(
      n_probes = nrow(dmps),
      n_in_te = sum(dmps$in_te),
      odds_ratio = fisher_test$estimate,
      pval = fisher_test$p.value
    )
  }

  # By TE class
  cat("\n--- Enrichment by TE Class ---\n")

  te_classes <- c("LINE", "SINE", "LTR", "DNA", "SVA")
  class_results <- list()

  for (te_class in te_classes) {
    class_probes <- overlap_df %>%
      filter(te_class == !!te_class) %>%
      pull(probe_id) %>%
      unique()

    for (direction in c("hypomethylated", "hypermethylated")) {
      dmps <- dmp_results %>% filter(status == direction)
      bg <- dmp_results %>% filter(status == "not_significant")

      mat <- matrix(
        c(sum(dmps$probe_id %in% class_probes),
          nrow(dmps) - sum(dmps$probe_id %in% class_probes),
          sum(bg$probe_id %in% class_probes),
          nrow(bg) - sum(bg$probe_id %in% class_probes)),
        nrow = 2, byrow = TRUE
      )

      fisher_test <- fisher.test(mat)

      class_results[[paste(te_class, direction, sep = "_")]] <- data.frame(
        te_class = te_class,
        direction = direction,
        n_dmp_in_te = sum(dmps$probe_id %in% class_probes),
        n_dmp_total = nrow(dmps),
        pct_in_te = 100 * sum(dmps$probe_id %in% class_probes) / nrow(dmps),
        odds_ratio = fisher_test$estimate,
        pval = fisher_test$p.value
      )
    }
  }

  class_enrichment_df <- bind_rows(class_results)
  class_enrichment_df$fdr <- p.adjust(class_enrichment_df$pval, method = "BH")

  cat("\nTE class enrichment (probe-level):\n")
  print(class_enrichment_df %>% arrange(pval))

  # DMR-level TE overlap analysis
  cat("\n--- DMR-level TE Overlap Analysis ---\n")

  dmr_te_overlap <- NULL
  if (nrow(dmr_table) > 0 && "chr" %in% colnames(dmr_table)) {
    dmr_gr <- GRanges(
      seqnames = dmr_table$chr,
      ranges = IRanges(start = dmr_table$start, end = dmr_table$end),
      status = dmr_table$status,
      value = dmr_table$value
    )

    te_gr <- GRanges(
      seqnames = rmsk$chr,
      ranges = IRanges(start = rmsk$start, end = rmsk$end),
      repClass = rmsk$repClass,
      repFamily = rmsk$repFamily
    )

    dmr_te_hits <- findOverlaps(dmr_gr, te_gr)

    if (length(dmr_te_hits) > 0) {
      dmr_te_overlap <- data.frame(
        dmr_index = queryHits(dmr_te_hits),
        dmr_chr = as.character(seqnames(dmr_gr)[queryHits(dmr_te_hits)]),
        dmr_start = start(dmr_gr)[queryHits(dmr_te_hits)],
        dmr_end = end(dmr_gr)[queryHits(dmr_te_hits)],
        dmr_status = dmr_gr$status[queryHits(dmr_te_hits)],
        dmr_value = dmr_gr$value[queryHits(dmr_te_hits)],
        te_class = te_gr$repClass[subjectHits(dmr_te_hits)],
        te_family = te_gr$repFamily[subjectHits(dmr_te_hits)]
      )

      cat("DMRs overlapping TEs:", length(unique(dmr_te_overlap$dmr_index)), "/", nrow(dmr_table), "\n")

      # Summarize by TE class
      dmr_te_summary <- dmr_te_overlap %>%
        filter(dmr_status != "not_significant") %>%
        group_by(dmr_status, te_class) %>%
        summarise(n_dmrs = n_distinct(dmr_index), .groups = "drop")

      if (nrow(dmr_te_summary) > 0) {
        cat("\nDMR-TE overlaps by class:\n")
        print(dmr_te_summary)
      }
    } else {
      cat("No DMR-TE overlaps found.\n")
    }
  }

  return(list(
    overall = fisher_overall,
    by_direction = enrichment_by_direction,
    by_class = class_enrichment_df,
    probe_annotated = dmp_results,
    dmr_table = dmr_table,
    dmr_te_overlap = dmr_te_overlap
  ))
}

# ------------------------------------------------------------------------------
# Plotting functions
# ------------------------------------------------------------------------------
plot_enrichment_heatmap <- function(class_enrichment_df, output_file) {
  p <- class_enrichment_df %>%
    mutate(
      log_or = log2(odds_ratio),
      sig = ifelse(fdr < 0.05, "*", "")
    ) %>%
    ggplot(aes(x = direction, y = te_class, fill = log_or)) +
    geom_tile(color = "white") +
    geom_text(aes(label = paste0(round(odds_ratio, 2), sig)), size = 4) +
    scale_fill_gradient2(
      low = "#3498DB", mid = "white", high = "#E74C3C",
      midpoint = 0, name = "log2(OR)"
    ) +
    labs(
      title = "TE Enrichment at Differentially Methylated Sites",
      subtitle = "* = FDR < 0.05",
      x = "Methylation Status (Mutant vs WT)",
      y = "TE Class"
    ) +
    theme_minimal() +
    theme(axis.text = element_text(size = 12))

  ggsave(output_file, p, width = 7, height = 5)
  return(p)
}

plot_enrichment_barplot <- function(class_enrichment_df, output_file) {
  p <- class_enrichment_df %>%
    mutate(sig = ifelse(fdr < 0.05, "*", "")) %>%
    ggplot(aes(x = te_class, y = odds_ratio, fill = direction)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
    geom_text(
      aes(label = sig, y = odds_ratio + 0.1),
      position = position_dodge(width = 0.8), size = 5
    ) +
    scale_fill_manual(values = c("hypomethylated" = "#3498DB", "hypermethylated" = "#E74C3C")) +
    labs(title = "TE Enrichment by Class", x = "TE Class", y = "Odds Ratio", fill = "Status") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(output_file, p, width = 8, height = 5)
  return(p)
}

plot_volcano <- function(probe_results, delta_thresh, pval_thresh, output_file) {
  # For bumphunter results, create a summary visualization
  # Show probe status distribution with TE annotation

  if (!"in_te" %in% colnames(probe_results)) {
    probe_results$in_te <- FALSE
  }

  # Summary bar plot of probes in DMRs
  summary_df <- probe_results %>%
    group_by(status, in_te) %>%
    summarise(n = n(), .groups = "drop") %>%
    mutate(te_label = ifelse(in_te, "In TE", "Not in TE"))

  p <- ggplot(summary_df, aes(x = status, y = n, fill = te_label)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = c("In TE" = "#E74C3C", "Not in TE" = "#3498DB")) +
    labs(
      title = "DMR Probe Distribution",
      subtitle = "Probes classified by bumphunter DMRs",
      x = "Methylation Status",
      y = "Number of Probes",
      fill = "TE Overlap"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(output_file, p, width = 8, height = 6)
  return(p)
}

# ------------------------------------------------------------------------------
# Plot DMR genomic distribution
# ------------------------------------------------------------------------------
plot_dmr_distribution <- function(dmr_table, output_file) {
  if (is.null(dmr_table) || nrow(dmr_table) == 0) {
    cat("No DMRs to plot.\n")
    return(NULL)
  }

  dmr_table <- dmr_table %>%
    filter(status != "not_significant") %>%
    mutate(width = end - start)

  p <- ggplot(dmr_table, aes(x = width, fill = status)) +
    geom_histogram(bins = 30, alpha = 0.7, position = "identity") +
    scale_fill_manual(values = c("hypomethylated" = "#3498DB", "hypermethylated" = "#E74C3C")) +
    scale_x_log10() +
    labs(
      title = "DMR Size Distribution",
      x = "DMR Width (bp, log scale)",
      y = "Count",
      fill = "Status"
    ) +
    theme_minimal()

  ggsave(output_file, p, width = 8, height = 5)
  return(p)
}
