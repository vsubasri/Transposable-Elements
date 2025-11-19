#### TRACKING HELPER FUNCTIONS ####

# Split HostSeq samples into filtering and analysis groups
# Args:
#   te_data: TE dataframe with sample column
#   filter_pct: Percentage of HostSeq samples to use for filtering (default 66%)
#   seed: Random seed for reproducibility (default 123)
# Returns:
#   List with:
#     - filter_samples: Vector of sample IDs for filtering group
#     - analysis_samples: Vector of sample IDs for analysis group
#     - te_filter: TE data with only filtering group HostSeq samples (for determining common TEs)
#     - te_analysis: TE data with filtering group HostSeq samples removed (for analysis)
split_hostseq_samples <- function(te_data, filter_pct = 66, seed = 123) {

  cat("\n===== SPLITTING HOSTSEQ SAMPLES (STRATIFIED BY ANCESTRY) =====\n")

  # Set seed for reproducibility
  set.seed(seed)

  # Get all HostSeq sample IDs with their ancestry
  hostseq_data <- te_data %>%
    filter(grepl("^HS_", sample)) %>%
    distinct(sample, .keep_all = TRUE) %>%
    select(sample, predicted_ancestry_thres)

  n_hostseq <- nrow(hostseq_data)

  if (n_hostseq == 0) {
    cat("Warning: No HostSeq samples found. Adding hostseq_group column with NA.\n")
    te_data$hostseq_group <- NA_character_
    return(list(
      filter_samples = character(0),
      analysis_samples = character(0),
      te_data = te_data
    ))
  }

  # Stratified split by ancestry
  filter_samples <- c()
  analysis_samples <- c()

  # Split within each ancestry group
  for (ancestry_group in unique(hostseq_data$predicted_ancestry_thres)) {
    samples_in_group <- hostseq_data %>%
      filter(predicted_ancestry_thres == ancestry_group) %>%
      pull(sample)

    n_in_group <- length(samples_in_group)
    n_filter_group <- round(n_in_group * filter_pct / 100)

    filter_in_group <- sample(samples_in_group, n_filter_group)
    analysis_in_group <- setdiff(samples_in_group, filter_in_group)

    filter_samples <- c(filter_samples, filter_in_group)
    analysis_samples <- c(analysis_samples, analysis_in_group)

    cat("Ancestry:", ancestry_group, "- Total:", n_in_group,
        "| Filter:", length(filter_in_group), "| Analysis:", length(analysis_in_group), "\n")
  }

  cat("\nTotal HostSeq samples:", n_hostseq, "\n")
  cat("Filter group (", filter_pct, "%):", length(filter_samples), "samples\n", sep = "")
  cat("Analysis group (", 100 - filter_pct, "%):", length(analysis_samples), "samples\n", sep = "")

  # Add hostseq_group column to label samples
  te_data$hostseq_group <- case_when(
    te_data$sample %in% filter_samples ~ "filter",
    te_data$sample %in% analysis_samples ~ "analysis",
    TRUE ~ NA_character_
  )

  cat("Labeled HostSeq samples in dataset (stratified by ancestry)\n")
  cat("=====================================\n\n")

  return(list(
    filter_samples = filter_samples,
    analysis_samples = analysis_samples,
    te_data = te_data
  ))
}

# Replace HostSeq filter group samples with analysis group samples
# Used after filtering to swap in the independent analysis group for final analysis
# Args:
#   te_data: Processed TE data that includes filter group HostSeq samples
#   filter_samples: Vector of filter group sample IDs
#   te_analysis: TE data with only analysis group HostSeq + non-HostSeq samples
# Returns:
#   TE data with filter group removed and analysis group added
replace_hostseq_with_analysis_group <- function(te_data, filter_samples, te_analysis, verbose = FALSE) {

  if (verbose) {
    cat("\n===== REPLACING HOSTSEQ FILTER GROUP WITH ANALYSIS GROUP =====\n")
    # Count before
    hostseq_before <- sum(grepl("^HS_", te_data$sample))
    cat("HostSeq rows before replacement:", hostseq_before, "\n")
  }

  # Remove filter group HostSeq samples
  te_no_filter <- te_data %>%
    filter(!sample %in% filter_samples)

  # Get analysis group HostSeq samples from te_analysis
  te_analysis_hostseq <- te_analysis %>%
    filter(grepl("^HS_", sample))

  # Ensure column type compatibility before binding
  # Get common columns
  common_cols <- intersect(names(te_no_filter), names(te_analysis_hostseq))

  # For each common column, ensure matching types
  for (col in common_cols) {
    class_nofilter <- class(te_no_filter[[col]])[1]
    class_analysis <- class(te_analysis_hostseq[[col]])[1]

    if (class_nofilter != class_analysis) {
      # Convert both to character for safety, unless both are numeric-compatible
      if ((class_nofilter %in% c("numeric", "integer", "double")) &&
          (class_analysis %in% c("numeric", "integer", "double"))) {
        te_no_filter[[col]] <- as.numeric(te_no_filter[[col]])
        te_analysis_hostseq[[col]] <- as.numeric(te_analysis_hostseq[[col]])
      } else if (class_nofilter == "logical" || class_analysis == "logical") {
        te_no_filter[[col]] <- as.character(te_no_filter[[col]])
        te_analysis_hostseq[[col]] <- as.character(te_analysis_hostseq[[col]])
      } else {
        te_no_filter[[col]] <- as.character(te_no_filter[[col]])
        te_analysis_hostseq[[col]] <- as.character(te_analysis_hostseq[[col]])
      }
    }
  }

  # Combine - use bind_rows to handle different columns
  te_with_analysis <- bind_rows(te_no_filter, te_analysis_hostseq)

  if (verbose) {
    # Count after
    hostseq_after <- sum(grepl("^HS_", te_with_analysis$sample))
    cat("HostSeq rows after replacement:", hostseq_after, "\n")
    cat("Analysis group HostSeq samples:", length(unique(te_analysis_hostseq$sample)), "\n")
    cat("==============================================================\n\n")
  }

  return(te_with_analysis)
}

# Analyze sensitivity of common TE filtering to HostSeq sample size
# Tests different numbers of HostSeq samples for filtering and plots how many TEs are removed
# Args:
#   te_data: TE data with all HostSeq samples (before splitting)
#   rare_gnomad: gnomAD threshold percentage (e.g., 3)
#   rare_hostseq: HostSeq threshold percentage (e.g., 3)
#   sample_sizes: Vector of sample sizes to test (if NULL, uses seq from 10% to 100% by 10%)
#   output_dir: Directory to save plot
#   output_prefix: Prefix for output file name
#   seed: Random seed for reproducibility
# Returns:
#   ggplot object
analyze_hostseq_filter_sensitivity <- function(te_data, rare_gnomad = 3, rare_hostseq = 3,
                                               sample_sizes = NULL, output_dir = NULL,
                                               output_prefix = "te", seed = 123) {

  cat("\n===== HOSTSEQ FILTER SENSITIVITY ANALYSIS =====\n")

  # Get all HostSeq samples
  all_hostseq <- unique(te_data$sample[grepl("^HS_", te_data$sample)])
  n_total_hostseq <- length(all_hostseq)

  cat("Total HostSeq samples available:", n_total_hostseq, "\n")

  if (n_total_hostseq == 0) {
    cat("Warning: No HostSeq samples found. Cannot perform sensitivity analysis.\n")
    return(NULL)
  }

  # Define sample sizes to test if not provided
  if (is.null(sample_sizes)) {
    # Test 10%, 20%, ..., 100%
    percentages <- seq(10, 100, by = 10)
    sample_sizes <- round(n_total_hostseq * percentages / 100)
    # Make sure we don't exceed total
    sample_sizes <- unique(pmin(sample_sizes, n_total_hostseq))
  }

  cat("Testing filter group sizes:", paste(sample_sizes, collapse = ", "), "\n")

  # Count initial TEs
  n_tes_initial <- nrow(te_data)
  n_unique_tes_initial <- te_data %>%
    distinct(SV_chrom, SV_start, SV_end, SV_length, ALT) %>%
    nrow()

  cat("Initial TEs (total insertions):", n_tes_initial, "\n")
  cat("Initial TEs (unique):", n_unique_tes_initial, "\n")

  # Store results
  results <- data.frame(
    n_samples = integer(),
    n_tes_filtered = integer(),
    n_unique_tes_filtered = integer(),
    n_tes_remaining = integer(),
    n_unique_tes_remaining = integer()
  )

  # Set seed
  set.seed(seed)

  # Test each sample size
  for (n_filter in sample_sizes) {
    cat("\nTesting with", n_filter, "HostSeq samples...\n")

    # Randomly sample HostSeq samples for this test
    filter_samples_test <- sample(all_hostseq, n_filter)

    # Create filter dataset (filter group + non-HostSeq)
    te_filter_test <- te_data %>%
      filter(!grepl("^HS_", sample) | sample %in% filter_samples_test)

    # Apply filtering (use the appropriate filter function based on data type)
    if ("data.table" %in% class(te_data)) {
      te_filtered_test <- filter_common_hostseq_germline_te(
        te = te_filter_test,
        rare_gnomad_threshold = rare_gnomad,
        rare_hostseq_threshold = rare_hostseq
      )
    } else {
      te_filtered_test <- filter_common_hostseq_germline_te(
        te = te_filter_test,
        rare_gnomad_threshold = rare_gnomad,
        rare_hostseq_threshold = rare_hostseq
      )
    }

    # Count TEs after filtering
    n_tes_after <- nrow(te_filtered_test)
    n_unique_tes_after <- te_filtered_test %>%
      distinct(SV_chrom, SV_start, SV_end, SV_length, ALT) %>%
      nrow()

    # Calculate filtered counts
    n_tes_filtered <- n_tes_initial - n_tes_after
    n_unique_tes_filtered <- n_unique_tes_initial - n_unique_tes_after

    cat("  TEs filtered (total):", n_tes_filtered, "\n")
    cat("  TEs filtered (unique):", n_unique_tes_filtered, "\n")

    # Store results
    results <- rbind(results, data.frame(
      n_samples = n_filter,
      n_tes_filtered = n_tes_filtered,
      n_unique_tes_filtered = n_unique_tes_filtered,
      n_tes_remaining = n_tes_after,
      n_unique_tes_remaining = n_unique_tes_after
    ))
  }

  cat("\n--- Sensitivity Analysis Complete ---\n")
  cat("Results summary:\n")
  print(results)

  # Create plot
  p <- ggplot(results, aes(x = n_samples, y = n_unique_tes_filtered)) +
    geom_line(color = "blue", linewidth = 1) +
    geom_point(color = "blue", size = 3) +
    labs(
      title = "Sensitivity of Common TE Filtering to HostSeq Sample Size",
      x = "Number of HostSeq Samples Used for Filtering",
      y = "Number of Unique TEs Filtered Out as Common",
      subtitle = paste0("Total HostSeq samples available: ", n_total_hostseq,
                        " | gnomAD threshold: ", rare_gnomad, "% | HostSeq threshold: ", rare_hostseq, "%")
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5)
    ) +
    scale_x_continuous(breaks = sample_sizes)

  # Save plot
  if (!is.null(output_dir)) {
    ggsave(paste0(output_dir, output_prefix, "_hostseq_filter_sensitivity.png"),
           plot = p, width = 10, height = 6)
    cat("✓ Plot saved to:", paste0(output_prefix, "_hostseq_filter_sensitivity.png"), "\n")
  }

  cat("==============================================\n\n")

  return(p)
}

# Find TEs overlapping with methylation probe regions
# Args:
#   te_expand: TE expanded dataframe (must have SV_chrom, SV_start, SV_length, sample, TP53_status)
#   probe_file: Path to methylation probe CSV file (must have chr, pos, end columns)
#   min_overlap_pct: Minimum percentage of probe region that must overlap with TE (default 70%)
#   output_file: Path to save output CSV (optional)
# Returns:
#   Dataframe with probe regions and their TE overlaps
find_te_probe_overlaps <- function(te_expand, probe_file, min_overlap_pct = 0, output_file = NULL) {

  cat("\n===== FINDING TE-PROBE OVERLAPS =====\n")

  # Load probe data
  if (!file.exists(probe_file)) {
    cat("Error: Probe file not found:", probe_file, "\n")
    return(NULL)
  }

  meth_probes <- read.csv(probe_file, stringsAsFactors = FALSE)
  cat("Loaded", nrow(meth_probes), "methylation probes\n")

  # Check required columns
  required_cols <- c("chr", "pos", "end")
  missing_cols <- setdiff(required_cols, colnames(meth_probes))
  if (length(missing_cols) > 0) {
    cat("Error: Missing required columns in probe file:", paste(missing_cols, collapse = ", "), "\n")
    return(NULL)
  }

  # Calculate TE end from start + length
  te_expand <- te_expand %>%
    mutate(TE_end_calc = SV_start + SV_length)

  # Check required TE columns
  te_required_cols <- c("SV_chrom", "SV_start", "sample", "TP53_status")
  missing_te_cols <- setdiff(te_required_cols, colnames(te_expand))
  if (length(missing_te_cols) > 0) {
    cat("Error: Missing required columns in TE data:", paste(missing_te_cols, collapse = ", "), "\n")
    return(NULL)
  }

  cat("Minimum overlap percentage:", min_overlap_pct, "%\n")

  # Process each probe region
  results <- list()

  for (i in 1:nrow(meth_probes)) {
    probe <- meth_probes[i, ]
    probe_chr <- as.character(probe$chr)
    probe_start <- probe$pos
    probe_end <- probe$end
    probe_length <- probe_end - probe_start + 1

    # Create region identifier
    region_id <- paste0(probe_chr, ":", probe_start, "-", probe_end)

    # Find overlapping TEs
    overlapping_tes <- te_expand %>%
      filter(
        SV_chrom == probe_chr,
        # TE overlaps probe if TE end > probe start AND TE start < probe end
        TE_end_calc > probe_start,
        SV_start < probe_end
      ) %>%
      mutate(
        # Calculate overlap length
        overlap_start = pmax(SV_start, probe_start),
        overlap_end = pmin(TE_end_calc, probe_end),
        overlap_length = overlap_end - overlap_start + 1,
        overlap_pct = (overlap_length / probe_length) * 100
      ) %>%
      filter(overlap_pct >= min_overlap_pct)

    # Calculate statistics
    n_overlaps <- nrow(overlapping_tes)

    if (n_overlaps > 0) {
      samples <- paste(unique(overlapping_tes$sample), collapse = ";")

      # Calculate mean overlap percentage for this probe
      mean_overlap_pct <- mean(overlapping_tes$overlap_pct, na.rm = TRUE)

      # Count TP53 status
      tp53_counts <- overlapping_tes %>%
        distinct(sample, TP53_status) %>%
        count(TP53_status)

      n_tp53_mut <- tp53_counts %>% filter(TP53_status == "Mutant") %>% pull(n) %>% sum()
      n_tp53_wt <- tp53_counts %>% filter(TP53_status == "WT") %>% pull(n) %>% sum()

      if (length(n_tp53_mut) == 0) n_tp53_mut <- 0
      if (length(n_tp53_wt) == 0) n_tp53_wt <- 0
    } else {
      samples <- ""
      mean_overlap_pct <- 0
      n_tp53_mut <- 0
      n_tp53_wt <- 0
    }

    results[[i]] <- data.frame(
      region = region_id,
      chr = probe_chr,
      start = probe_start,
      end = probe_end,
      n_overlaps = n_overlaps,
      mean_overlap_pct = round(mean_overlap_pct, 2),
      samples = samples,
      n_tp53_mut = n_tp53_mut,
      n_tp53_wt = n_tp53_wt,
      stringsAsFactors = FALSE
    )
  }

  # Combine results
  results_df <- bind_rows(results)

  # Summary
  cat("\nResults summary:\n")
  cat("  Total probe regions:", nrow(results_df), "\n")
  cat("  Regions with TE overlaps:", sum(results_df$n_overlaps > 0), "\n")
  cat("  Total TE overlaps:", sum(results_df$n_overlaps), "\n")

  # Filter to only regions with overlaps for output
  results_with_overlaps <- results_df %>%
    filter(n_overlaps > 0)

  # Save to CSV if output file specified - only include rows with overlaps
  if (!is.null(output_file)) {
    write.csv(results_with_overlaps, output_file, row.names = FALSE)
    cat("✓ Results saved to:", output_file, "(", nrow(results_with_overlaps), "probes with overlaps)\n")
  }

  cat("=====================================\n\n")

  return(results_df)
}

# Initialize tracking for a filtering step
init_step_tracking <- function(df, sample_col = "sample") {
  list(
    te_count = nrow(df),
    sample_count = length(unique(df[[sample_col]]))
  )
}

# Calculate loss between two steps
calc_step_loss <- function(prev_count, curr_count) {
  loss <- prev_count - curr_count
  loss_pct <- if (prev_count > 0) (loss / prev_count) * 100 else 0
  list(loss = loss, loss_pct = loss_pct)
}

# Extract count from captured stdout by pattern
extract_count_from_stdout <- function(stdout_lines, pattern) {
  for (line in stdout_lines) {
    if (grepl(pattern, line)) {
      numbers <- as.numeric(unlist(regmatches(line, gregexpr("[0-9]+", line))))
      if (length(numbers) > 0) {
        return(numbers)
      }
    }
  }
  return(NA)
}

# Extract single value from stdout by pattern and regex group
extract_value_from_stdout <- function(stdout_lines, pattern, regex_group = "\\1") {
  for (line in stdout_lines) {
    if (grepl(pattern, line)) {
      value <- as.numeric(gsub(pattern, regex_group, line))
      return(value)
    }
  }
  return(NA)
}

# Count unique TE loci with error handling
count_unique_te_loci <- function(df, primary_cols = c("SV_chrom", "SV_start", "SV_end", "SV_length", "ALT"),
                                  fallback_cols = c("SV_chrom", "SV_start", "SV_end", "ALT"),
                                  fallback_value = NA) {
  tryCatch({
    if (all(primary_cols %in% colnames(df))) {
      result <- nrow(unique(df[, primary_cols]))
      cat("Counted unique TE loci using primary columns:", result, "\n")
      return(result)
    } else if (all(fallback_cols %in% colnames(df))) {
      result <- nrow(unique(df[, fallback_cols]))
      cat("Counted unique TE loci using fallback columns:", result, "\n")
      return(result)
    } else {
      # Silently return fallback value - columns not available
      return(fallback_value)
    }
  }, error = function(e) {
    cat("Error calculating unique TE loci:", e$message, "\n")
    return(fallback_value)
  })
}

# Extract step4 sample metrics from final_te_count
extract_step4_metrics <- function(final_te_count, step3_sample_count, process_stdout) {
  metrics <- list(
    total_samples = NA,
    samples_with_tes = NA,
    samples_no_tes = NA,
    total_te_insertions = NA
  )

  if ("combination" %in% colnames(final_te_count)) {
    total_combo <- final_te_count[final_te_count$combination == "total", ]
    metrics$samples_with_tes <- sum(total_combo$count > 0)
    metrics$samples_no_tes <- sum(total_combo$count == 0)
    metrics$total_samples <- nrow(total_combo)
    metrics$total_te_insertions <- sum(total_combo$count)
  } else {
    # Extract from process output
    no_hits_line <- grep("Adding .* samples with no TE insertions", process_stdout, value = TRUE)
    if (length(no_hits_line) > 0) {
      no_hits_added <- as.numeric(gsub(".*Adding (\\d+) samples.*", "\\1", no_hits_line[1]))
      metrics$total_samples <- step3_sample_count + no_hits_added
      metrics$samples_with_tes <- step3_sample_count
      metrics$samples_no_tes <- no_hits_added
    } else if ("sample" %in% colnames(final_te_count)) {
      metrics$total_samples <- length(unique(final_te_count$sample))
    } else {
      metrics$total_samples <- nrow(final_te_count)
    }
  }

  return(metrics)
}

# Write captured output to file
write_stdout_to_file <- function(stdout_lines, file_path, header, append = TRUE) {
  tryCatch({
    # Ensure the output lines are character vectors
    if (is.null(stdout_lines) || length(stdout_lines) == 0) {
      warning("No output to write for: ", header)
      return(invisible(NULL))
    }

    # Create directory if it doesn't exist
    output_dir <- dirname(file_path)
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }

    # Write output with explicit connection handling
    con <- file(file_path, open = if(append) "a" else "w")
    on.exit(close(con), add = TRUE)

    writeLines(c(paste0("\n--- ", header, " ---"), stdout_lines), con)
    flush(con)

    invisible(TRUE)
  }, error = function(e) {
    cat("ERROR in write_stdout_to_file:", e$message, "\n")
    cat("Attempting fallback write method...\n")
    tryCatch({
      # Fallback: use cat() directly
      cat(c(paste0("\n--- ", header, " ---\n"), stdout_lines, "\n"),
          file = file_path, append = append, sep = "\n")
      invisible(TRUE)
    }, error = function(e2) {
      cat("CRITICAL ERROR: Could not write to file:", e2$message, "\n")
      invisible(FALSE)
    })
  })
}

#### TE FILTERING FUNCTIONS ####

# TE type filtering function
filter_te_types <- function(te_data, include_alu = TRUE, include_sva = TRUE) {
  initial_count <- nrow(te_data)
  te_types_to_remove <- c()

  if (!include_alu) te_types_to_remove <- c(te_types_to_remove, "ALU")
  if (!include_sva) te_types_to_remove <- c(te_types_to_remove, "SVA")

  if (length(te_types_to_remove) > 0) {
    pattern <- paste0("-(", paste(te_types_to_remove, collapse = "|"), ")$")
    te_data <- te_data[!grepl(pattern, ID)]
    removed_count <- initial_count - nrow(te_data)
    cat(sprintf("Removed %d insertions (%s) from analysis\n", removed_count, paste(te_types_to_remove, collapse = ", ")))
  } else {
    cat("Including all TE types (LINE1, ALU, SVA)\n")
  }
  cat(sprintf("Remaining insertions: %d\n", nrow(te_data)))

  return(te_data)
}

# Cohort summary statistics table
create_cohort_summary <- function(te_data, hostseq_cancer = NULL) {
  cohorts <- sort(unique(te_data$cohort))
  cohorts <- cohorts[!is.na(cohorts)]

  cohort_table <- data.frame()

  for (coh in cohorts) {
    cohort_data <- te_data[te_data$cohort == coh, ]
    total <- nrow(cohort_data)
    tp53_wt <- sum(cohort_data$TP53_status == "WT", na.rm = TRUE)
    tp53_mut <- sum(cohort_data$TP53_status == "Mutant", na.rm = TRUE)

    # Determine cancer/control counts
    if (!is.null(hostseq_cancer)) {
      # For germline data with hostseq_cancer mapping
      cancer_yes <- sum(cohort_data$sample %in% hostseq_cancer$sample[hostseq_cancer$Cancer == "Yes"], na.rm = TRUE)
    } else {
      # For tumour data - all samples are cancer
      cancer_yes <- total
    }

    cohort_table <- rbind(cohort_table, data.frame(
      Cohort = coh,
      Total = total,
      TP53_WT = tp53_wt,
      TP53_Mutant = tp53_mut,
      Cancer = cancer_yes
    ))
  }

  # Add TOTAL row
  total_row <- data.frame(
    Cohort = "TOTAL",
    Total = sum(cohort_table$Total),
    TP53_WT = sum(cohort_table$TP53_WT),
    TP53_Mutant = sum(cohort_table$TP53_Mutant),
    Cancer = sum(cohort_table$Cancer)
  )
  cohort_table <- rbind(cohort_table, total_row)

  return(cohort_table)
}

# graph functions
#colours <- c("#DDD8C4", "#5FBFF9")
#colour_palette_3<- c("#DDD8C4", "#5FBFF9", "#235789")
colours <- c("#B1D586", "#0080A3")
colours_3 <- c("#B1D586", "#0080A3", "#AB1368")
colours_sex <- c("#AB1368", "#0080A3")
mutation_colours <- c("#99E9FF", "#0090B8", "#004052")
mutation_colours <- c("#DDD8C4", "#0090B8", "#004052")
colour_palette_3<- c("#AB1368", "#B1D586", "#0080A3")
colour_palette_4_seq<- c("#005066", "#0080A3", "#00B0E0", "#99E9FF")
colour_palette_4<- c("#B1D586", "#0080A3", "#AB1368", '#ffcc66')
colour_palette_5 <- c("#DDD8C4", "#5FBFF9", '#99ff99', '#ffcc66', '#ff6666')
color_palette_6 <- c("#DDD8C4", "#5FBFF9", "#235789", '#99ff99', '#ffcc66', '#ff6666')


# set theme 
theme_set(theme(
  # Axis titles larger than axis labels
  axis.title = element_text(size = 14, color = "black"),
  axis.text = element_text(size = 12, color = "black"),
  
  # Remove grid lines and background
  panel.grid = element_blank(),
  panel.background = element_blank(),
  
  # Black x and y axis
  axis.line = element_line(color = "black"),
  
  # No border around the plot
  plot.background = element_blank(),
  panel.border = element_blank(),
  
  # No legend border
  legend.key = element_blank(),
  
  # Set legend title and text size
  legend.title = element_text(size = 12),
  legend.text = element_text(size = 12)
))



# prep dataframes
prep <- function(df){
  #df$TP53_status <- ifelse(df$TP53_status=="No", "Control", "LFS") # change elements of TP53_status
  df$sample <- ifelse(df$cohort == "KiCS", sprintf("%04d", as.numeric(df$sample)), df$sample) # add leading 0 to kics sample
  df<- df%>% filter(!is.na(sample)) # remove where sample is empty. these are entries where kids id is unknown
  df$tumor_type<- ifelse(df$tumor_type=="ARMS", "RMS", df$tumor_type) # combine RMS samples 
  df$tumor_type<- ifelse(df$tumor_type=="ERMS", "RMS", df$tumor_type) # combine RMS samples 
  return(df)
}

prep_sv <- function(df){
  data <- prep(df)
  sv_aff <- data  
  sv_aff_df <- data %>% filter(tumor_type!="U") # only samples affected with cancer
  sv_lfs_df <- data %>%  # only LFS
    filter(TP53_status=="1") %>%
    mutate(Cancer = ifelse(tumor_type != "U", "Affected", "Unaffected")) # afftected or unaffected
  return(list(sv_aff = sv_aff_df, sv_lfs = sv_lfs_df)) 
}

# prep clinical file
prep_clinical <- function(clinical){
  #clinical <- clinical %>% rename("sample"="germline_sample")
  clinical$TP53_status <- ifelse(clinical$TP53_status==0, "WT", "Mutant") # change elements of TP53_status
  clinical$tumor_type<- ifelse(clinical$tumor_type == "ERMS", "RMS", ifelse(clinical$tumor_type == "ARMS", "RMS", clinical$tumor_type))
  clinical$tumor_type<- ifelse(clinical$tumor_type == "HGG", "G", ifelse(clinical$tumor_type == "LGG", "G", clinical$tumor_type))
  clinical$TP53_status <- as.factor(clinical$TP53_status)
  clinical$cluster<- as.factor(clinical$cluster)
  clinical<- clinical %>%
    mutate(base_sample = str_remove(sample, "(_[-A-Za-z0-9]+)?_T$"))
  clinical <- clinical[!duplicated(clinical), ] # remove duplicates
  
  # cancer cohort column
  clinical <- clinical %>%
    mutate(cancer_cohort = case_when(
      TP53_status == "Mutant" & tumor_type != "U" ~ "Cancer mutant TP53",
      TP53_status == "WT" & tumor_type == "U" ~ "No cancer WT TP53",
      TP53_status == "Mutant" & tumor_type == "U" ~ "No cancer mutant TP53",
      TP53_status == "WT" & tumor_type != "U" ~ "Cancer WT TP53",
      TRUE ~ NA_character_  # Default case when none of the above conditions are met
    )) %>%
    mutate(cancer_cohort = factor(cancer_cohort, levels = c("No cancer WT TP53", "Cancer WT TP53", "Cancer mutant TP53")))
  
  
  return(clinical)
} 

prep_hostseq <- function(df) {
  df %>%
    filter(hostseq_cancer == 1) %>%
    mutate(sample = sub("_N$", "", sample))
}

prep_metrics <- function(df) {

  ## --- 1. Rename columns --------------------------------------------------
  names(df) <- c("sample", "mean_cov", "sd_cov", "med_cov", "total_reads",
                 "mean_read_len", "sd_read_len", "med_read_len",
                 "pct_chimeras", "avg_quality")

  ## --- 2. Normalise sample IDs -------------------------------------------
  # For HostSeq samples (start with 20-, 21-, 22-, or HS_), just add _N if not there
  # For other samples ending with _N, normalize to base_sample_N format
  df$sample <- ifelse(grepl("_N$", df$sample),
                      ifelse(grepl("^(20-|21-|22-|HS_)", df$sample),
                             df$sample,  # HostSeq: keep as is
                             paste0(sub("_.*", "", df$sample), "_N")),  # Non-HostSeq: extract before first underscore
                      paste0(df$sample, "_N"))  # No _N at end: add it

  ## --- 3. Deduplicate -----------------------------------------------------
  setDT(df)                      # convert to data.table *in-place*
  df   <- unique(df, by = "sample")   # keep first row for each sample

  return(df[])
}

prep_metrics_tumour <- function(df) {

  ## --- 1. Rename columns --------------------------------------------------
  names(df) <- c("sample", "mean_cov", "sd_cov", "med_cov", "total_reads",
                 "mean_read_len", "sd_read_len", "med_read_len",
                 "pct_chimeras", "avg_quality")

  ## --- 2. Normalise sample IDs -------------------------------------------
  # For tumor samples ending with _T, keep them as is
  # For HostSeq samples (start with 20-, 21-, 22-, or HS_), just add _T if not there
  # For other samples ending with _N, normalize to base_sample_N format
  df$sample <- ifelse(grepl("_T$", df$sample),
                      df$sample,  # Tumor: keep as is
                      ifelse(grepl("_N$", df$sample),
                             ifelse(grepl("^(20-|21-|22-|HS_)", df$sample),
                                    df$sample,  # HostSeq: keep as is
                                    paste0(sub("_.*", "", df$sample), "_N")),  # Non-HostSeq: extract before first underscore
                             paste0(df$sample, "_T")))  # No suffix: add _T

  ## --- 3. Deduplicate -----------------------------------------------------
  setDT(df)                      # convert to data.table *in-place*
  df   <- unique(df, by = "sample")   # keep first row for each sample

  return(df[])
}
 
prep_ancestry <-function(df) {
  # Ensure the column exists
  if (!"indivID" %in% colnames(df)) {
    stop("The dataframe does not have an 'indivID' column.")
  }

  # Create a new sample column based on indivID with updated rules
  df$sample <- df$indivID

  # Remove 'sorted_fixed_control_' from the beginning
  df$sample <- gsub("^sorted_fixed_control_", "", df$sample)

  # Remove 'sorted_fixed_' from the beginning
  df$sample <- gsub("^sorted_fixed_", "", df$sample)

  # Remove 'KICS_' from the beginning
  df$sample <- gsub("^KICS_", "", df$sample)

  # Remove 'LFS_' from the beginning
  df$sample <- gsub("^LFS_", "", df$sample)

  # Remove 'H_LC-' from the beginning
  df$sample <- gsub("^H_LC-", "", df$sample)

  # If starts with H (flowcell ID) or starts with 20/21/22 and ends with A-02-00, add HS_ at beginning
  df$sample <- ifelse(grepl("^H", df$sample) | grepl("^(20|21|22).*A-02-00", df$sample),
                      paste0("HS_", df$sample),
                      df$sample)

  # Remove '-N_realigned_recalibrated' from the end
  df$sample <- gsub("-N_realigned_recalibrated$", "", df$sample)

  # Remove '-N' from the end
  df$sample <- gsub("-N$", "", df$sample)

  # Remove '_merged' from the end
  df$sample <- gsub("_merged$", "", df$sample)

  # Remove '-G' from the end
  df$sample <- gsub("-G$", "", df$sample)

  # Remove '_G1' from the end
  df$sample <- gsub("_G1$", "", df$sample)

  # Add '_N' at the end only if it's not already there
  df$sample <- ifelse(grepl("_N$", df$sample), df$sample, paste0(df$sample, "_N"))

  # Keep sample, predicted_ancestry_thres, and mapped_label (if available)
  cols_to_keep <- c("sample", "predicted_ancestry_thres")
  if ("mapped_label" %in% colnames(df)) {
    cols_to_keep <- c(cols_to_keep, "mapped_label")
  }
  df <- df[, cols_to_keep]

  return(df)
}

# Load and prepare KICS germline sample type data
prep_kics_sample_type <- function(sample_type_path, add_suffix = "_N") {
  # Load KICS germline sample type data from CSV
  kics_sample_type <- read.csv(sample_type_path, stringsAsFactors = FALSE)

  # Pad sample column with leading zeros (10 -> 0010)
  kics_sample_type$sample <- sprintf("%04d", as.numeric(kics_sample_type$sample))

  # Add suffix if specified (e.g., _N for germline, _T for tumor)
  if (!is.null(add_suffix) && add_suffix != "") {
    kics_sample_type$sample <- paste0(kics_sample_type$sample, add_suffix)
  }

  return(kics_sample_type)
}

# Analyze TE burden by sample type (KICS)
analyze_sample_type <- function(te_data, sample_type_data, types, plot_dir,
                                filter_types = c("Blood", "Fibroblasts", "Tissue (fresh)"),
                                min_samples = 3) {
  # Merge sample type with TE data
  te_sample_type <- merge_dfs(te_data, sample_type_data, include_all_x = TRUE,
                              print_info = TRUE, dataset_name = "sample_type")

  # Filter for specific sample types
  te_sample_type_filtered <- te_sample_type %>%
    filter(sample_type %in% filter_types)

  cat("Samples with sample type data:", nrow(te_sample_type_filtered), "\n")
  cat("Sample type distribution:\n")
  print(table(te_sample_type_filtered$sample_type))

  # Generate plots for all TE types
  plots <- vector("list", length(types))
  for (i in seq_along(types)) {
    type_label <- ifelse(is.na(types[i]), "all", types[i])

    plots[[i]] <- plot_count_kruskal_nogroup(
      te_sample_type_filtered,
      column = "sample_type",
      min = min_samples,
      chr = NA,
      type = types[i],
      x_lab = "Sample Type",
      y_lab = "Repeat count",
      log_scale = TRUE
    )

    # Save plot
    ggsave(
      paste0(plot_dir, "counts_cohort/te_count_kics_sampletype_", type_label, ".pdf"),
      plot = plots[[i]],
      width = 9,
      height = 5
    )
  }

  return(list(
    data = te_sample_type_filtered,
    plots = plots
  ))
}

replace_nohit_samples <- function(nohits, clinical) {
  # Loop through each sample in nohits$V1
  for (i in seq_along(nohits$V1)) {
    # Check if the sample in nohits$V1 is in the form "number_T"
    if (grepl("^[0-9]+_T$", nohits$V1[i])) {
      # Extract the number part
      sample_number <- sub("_T$", "", nohits$V1[i])
      print(sample_number)
      # Adjusted regular expression to match any characters between sample_number and "_T"
      matched_sample <- clinical$sample[grepl(paste0("^", sample_number, ".*_T$"), clinical$sample)]
    
      # If a match is found, replace the nohits sample and print the replacement
      if (length(matched_sample) == 1) {
        cat("Replacing:", nohits$V1[i], "with", matched_sample, "\n")
        nohits$V1[i] <- matched_sample
       }
      }
    }
    
    # Return the modified list
    return(nohits)
  }
  
expand_clinical <- function(df, clinical, nohits) {
  # Remove "sorted_fixed_" prefix if present
  if (any(grepl("^sorted_fixed_", df$Samples_ID))) {
    df$Samples_ID <- sub("^sorted_fixed_", "", df$Samples_ID)
  }
  
  # Prepare nohits samples
  nohits_df <- nohits %>%
    rename(Samples_ID = V1) %>%
    mutate(base_sample = gsub("(_[A-Za-z0-9]+)?_T$", "", Samples_ID))
  
  # Separate df into samples starting with "0" and others
  df_with_zero <- df %>% filter(grepl("^0", Samples_ID))
  df_no_zero <- df %>% filter(!grepl("^0", Samples_ID)) %>%
    mutate(base_sample = gsub("(_[A-Za-z0-9]+)?_T$", "", Samples_ID))
  
  # Combine unique samples from df and nohits
  unique_te_samples <- bind_rows(
    df_no_zero %>% distinct(base_sample, Samples_ID),
    nohits_df
  ) %>%
    rename(sample = Samples_ID)
  
  # Format clinical and unique_te_samples consistently
  clinical <- clinical %>%
    mutate(base_sample = gsub("(_[A-Za-z0-9]+)?_T$", "", sample) %>% tolower() %>% trimws())
  unique_te_samples <- unique_te_samples %>%
    mutate(base_sample = tolower(base_sample) %>% trimws())
  
  # Expand clinical
  expanded_clinical_no_zero <- clinical %>%
    rename(clinical_sample = sample) %>%
    full_join(unique_te_samples, by = "base_sample")  # Use full_join to ensure no samples are lost
  
  # Combine clinical rows for samples starting with "0"
  expanded_clinical <- bind_rows(
    clinical %>% filter(sample %in% df_with_zero$Samples_ID),
    expanded_clinical_no_zero
  )
  
  # Check final size
  cat("Final expanded clinical rows:", nrow(expanded_clinical), "\n")
  return(expanded_clinical)
}

expand_clinical_with_input <- function(input_df, clinical_df) {
  # Step 1: Clean the Samples_ID in input_df
  input_df <- input_df %>%
    mutate(Samples_ID = str_remove(Samples_ID, "^sorted_fixed_"),
           base_sample = str_remove(Samples_ID, "(_[A-Za-z0-9]+)?_T$"))  # Extract base sample
  
  # Step 2: Extract base samples from clinical_df
  clinical_df <- clinical_df %>%
    mutate(base_sample = str_remove(sample, "(_[A-Za-z0-9]+)?_T$"))
  
  # Step 3: Identify all detailed samples from input_df
  detailed_samples <- input_df %>%
    distinct(Samples_ID, base_sample)
  print(head(detailed_samples))
  # Step 4: Join clinical information with detailed samples
  expanded_clinical <- detailed_samples %>%
    left_join(clinical_df, by = "base_sample") %>%  # Match on base_sample
    mutate(sample = Samples_ID) %>%  # Replace `sample` with detailed `Samples_ID`
    select(-Samples_ID, -base_sample)  # Clean up intermediate columns
  
  return(expanded_clinical)
}

filter_by_another_df <- function(df, exclude_df) {
  # Filter out rows where the 'sample' column in df matches any in exclude_df$sample
  filtered_df <- df[!df$base_sample %in% exclude_df$sample, ]
  return(filtered_df)
}

filter_by_metrics <- function(df, metrics_df, column_names, thresholds) {
  # Ensure column_names and thresholds are the same length
  if (length(column_names) != length(thresholds)) {
    stop("The length of column_names and thresholds must be the same.")
  }
  
  for (i in seq_along(column_names)) {
    column_name <- column_names[i]
    threshold <- thresholds[i]
    
    if (column_name == "sd_cov/mean_cov") {
      # Compute the ratio dynamically
      metrics_df$ratio <- metrics_df$sd_cov / metrics_df$mean_cov
      # Identify samples to exclude based on the ratio and threshold
      exclude_samples <- metrics_df$sample[metrics_df$ratio < threshold]
      # Filter out rows from df where the sample is in the exclusion list
      df <- df[!df$sample %in% exclude_samples, ]
    } else if (column_name == "pct_chimeras") {
      # Special logic for pct_chimeras to remove rows where value is > threshold
      exclude_samples <- metrics_df$sample[metrics_df[[column_name]] > threshold]
      df <- df[!df$sample %in% exclude_samples, ]
    } else {
      # Standard filtering for other columns
      exclude_samples <- metrics_df$sample[metrics_df[[column_name]] < threshold]
      df <- df[!df$sample %in% exclude_samples, ]
    }
  }
  
  return(df)
}

filter_by_multiple_criteria <- function(df, df_nonproband, df_noconsent, df_hostseqcancer, metrics_df, column_names, thresholds, type, export_filtered = TRUE, output_dir = NULL) {
  sample_id_col <- "sample"  # Column containing unique sample identifiers

  if (type == "T") {
    # set up base sample
    df$base_sample <- df$ID
    df$base_sample <- sub("_T.*", "", df$base_sample)  # Remove "_T" and everything after
    df$base_sample <- sub("_.*", "", df$base_sample)  # Remove "_" and everything after
    df$base_sample <- sub("-.*", "", df$base_sample) # Remove "-T" and everything after
  } else {
    df$base_sample <- sub("_N.*", "", df$sample)
  }

  # Initialize tracking dataframe for all samples
  all_samples <- unique(df[[sample_id_col]])
  filtered_samples_tracker <- data.frame(
    sample = all_samples,
    filter_reason = "Passed",
    stringsAsFactors = FALSE
  )

  # Step 1: Filter by df_nonproband
  initial_samples <- length(unique(df[[sample_id_col]]))
  nonproband_samples <- unique(df[[sample_id_col]][df[[sample_id_col]] %in% df_nonproband$V1])
  df_filtered <- filter_by_another_df(df, df_nonproband)
  remaining_samples <- length(unique(df_filtered[[sample_id_col]]))
  filtered_out_nonproband <- initial_samples - remaining_samples
  cat("Filtered out", filtered_out_nonproband, "samples due to non proband normal\n")
  filtered_samples_tracker$filter_reason[filtered_samples_tracker$sample %in% nonproband_samples] <- "Non-proband"

  # Step 2: Filter by df_noconsent
  count_after_nonproband <- remaining_samples
  noconsent_samples <- unique(df_filtered[[sample_id_col]][df_filtered[[sample_id_col]] %in% df_noconsent$V1])
  df_filtered <- filter_by_another_df(df_filtered, df_noconsent)
  remaining_samples <- length(unique(df_filtered[[sample_id_col]]))
  filtered_out_noconsent <- count_after_nonproband - remaining_samples
  cat("Filtered out", filtered_out_noconsent, "samples due to no consent\n")
  filtered_samples_tracker$filter_reason[filtered_samples_tracker$sample %in% noconsent_samples] <- "No consent"

  # Step 3: Filter by hostseq cancer
  count_after_noconsent <- remaining_samples
  hostseq_samples <- unique(df_filtered[[sample_id_col]][df_filtered[[sample_id_col]] %in% df_hostseqcancer$V1])
  df_filtered <- filter_by_another_df(df_filtered, df_hostseqcancer)
  remaining_samples <- length(unique(df_filtered[[sample_id_col]]))
  filtered_out_hostseqcancer <- count_after_noconsent - remaining_samples
  cat("Filtered out", filtered_out_hostseqcancer, "samples due to host seq having cancer\n")
  filtered_samples_tracker$filter_reason[filtered_samples_tracker$sample %in% hostseq_samples] <- "HostSeq cancer"

  # Step 3.5: Filter samples with no SNV info
  no_snv_samples <- c("MDT-AP-0224", "MDT-AP-1110", "MDT-AP-1112", "MDT-AP-2006", "MDT-AP-2897")
  count_after_hostseqcancer <- remaining_samples
  # Use pattern matching to handle both tumor (_T) and germline (_N) suffixes
  no_snv_pattern <- paste0("^(", paste(no_snv_samples, collapse="|"), ")(_T|_N)?$")
  nosnv_filtered <- unique(df_filtered[[sample_id_col]][grepl(no_snv_pattern, df_filtered[[sample_id_col]])])
  df_filtered <- df_filtered[!grepl(no_snv_pattern, df_filtered[[sample_id_col]]), ]
  remaining_samples <- length(unique(df_filtered[[sample_id_col]]))
  filtered_out_nosnv <- count_after_hostseqcancer - remaining_samples
  cat("Filtered out", filtered_out_nosnv, "samples due to no SNV info\n")
  filtered_samples_tracker$filter_reason[filtered_samples_tracker$sample %in% nosnv_filtered] <- "No SNV info"

  # Step 4: Filter by metrics_df using multiple columns and thresholds
  count_after_nosnv <- remaining_samples
  samples_before_metrics <- unique(df_filtered[[sample_id_col]])
  df_filtered <- filter_by_metrics(df_filtered, metrics_df, column_names, thresholds)
  samples_after_metrics <- unique(df_filtered[[sample_id_col]])
  metrics_filtered <- setdiff(samples_before_metrics, samples_after_metrics)
  remaining_samples <- length(unique(df_filtered[[sample_id_col]]))
  filtered_out_metrics <- count_after_nosnv - remaining_samples
  cat("Filtered out", filtered_out_metrics, "samples due to metrics with columns\n")
  filtered_samples_tracker$filter_reason[filtered_samples_tracker$sample %in% metrics_filtered] <- "Failed QC metrics"

  # Export filtered samples with metrics
  if (export_filtered && !is.null(output_dir)) {
    filtered_only <- filtered_samples_tracker[filtered_samples_tracker$filter_reason != "Passed", ]

    if (nrow(filtered_only) > 0) {
      # Merge with metrics data
      filtered_with_metrics <- merge(filtered_only, metrics_df, by = "sample", all.x = TRUE)

      # For samples that failed QC metrics, specify which metrics failed
      for (i in seq_along(column_names)) {
        col <- column_names[i]
        thresh <- thresholds[i]

        if (col %in% colnames(filtered_with_metrics)) {
          if (col == "pct_chimeras") {
            filtered_with_metrics[[paste0(col, "_failed")]] <- filtered_with_metrics[[col]] > thresh
          } else {
            filtered_with_metrics[[paste0(col, "_failed")]] <- filtered_with_metrics[[col]] < thresh
          }
        }
      }

      # Update filter_reason for metrics-filtered samples to include which metrics failed
      for (idx in which(filtered_with_metrics$filter_reason == "Failed QC metrics")) {
        failed_metrics <- c()
        for (i in seq_along(column_names)) {
          col <- column_names[i]
          fail_col <- paste0(col, "_failed")
          if (fail_col %in% colnames(filtered_with_metrics) &&
              !is.na(filtered_with_metrics[[fail_col]][idx]) &&
              filtered_with_metrics[[fail_col]][idx]) {
            failed_metrics <- c(failed_metrics, col)
          }
        }
        if (length(failed_metrics) > 0) {
          filtered_with_metrics$filter_reason[idx] <- paste0("Failed QC: ", paste(failed_metrics, collapse = ", "))
        }
      }

      # Write to CSV
      output_file <- paste0(output_dir, "filtered_samples_with_metrics_", type, ".csv")
      write.csv(filtered_with_metrics, output_file, row.names = FALSE)
      cat("Exported filtered samples with metrics to:", output_file, "\n")
    }
  }

  return(df_filtered)
}

rename_alt <- function(df) {
  df <- df %>% rename("ALT" = "SV_type")
  
  # Only rename Samples_ID to sample if sample column doesn't already exist
  if (!"sample" %in% colnames(df) && "Samples_ID" %in% colnames(df)) {
    df <- df %>% rename("sample" = "Samples_ID")
  }
  
  # Extract the final field of each ID after splitting by "-"
  last_field <- sapply(strsplit(as.character(df$ID), "-"), function(x) tail(x, 1))
  
  # Replace the ALT column with the last field for each row
  df$ALT <- last_field
  
  return(df)
}

replace_df_samples <- function(df, type) {
  if (type == "T") {
    # replace sample with patient_sample from ID column 
    df$sample <- ifelse(grepl("_T", df$ID), sub("_T.*", "_T", df$ID), NA)
    
    # a few samples didn't have ID column working so find sample from ID list
    missing_samples <- is.na(df$sample)
    
    if (any(missing_samples)) {
      # Extract IDLIST field and clean it in one vectorized step
      idlist_values <- sub(".*IDLIST=([^;]+).*", "\\1", df$INFO[missing_samples])
      # Extract up to and including _T from IDLIST
      df$sample[missing_samples] <- sub("_T.*", "_T", idlist_values)
    }
    
    # for patients with no sample number get sample number from clinical because those only have one sample
    # Create a dictionary for mapping clinical$base_sample to clinical$sample
    clinical_dict <- setNames(clinical$sample, clinical$base_sample)
    
    # Identify all unique four-digit sample names in df
    unique_four_digit_samples <- unique(sub("_T$", "", df$sample[grepl("^[0-9]{4}_T$", df$sample)]))
    
    # Map the unique sample names to their clinical samples
    mapped_values <- clinical_dict[unique_four_digit_samples]
    
    # Replace all matching samples in df
    # For each unique four-digit sample, replace all occurrences in df
    for (i in seq_along(unique_four_digit_samples)) {
      df$sample[df$sample == paste0(unique_four_digit_samples[i], "_T")] <- mapped_values[i]
    }
    
  } else if (type == "N") {
    # replace sample with patient_sample from ID column
    df <- df %>%
      mutate(
        sample = sub("_N.*", "_N", ID),  # Update sample to just _N if ID contains _N
        sample = sub("^LFS_", "", sample)  # Remove LFS_ from sample
      )

    # Hardcoded sample name fix
    df$sample <- ifelse(df$sample == "830_N", "830_1_N", df$sample)
  }
  
  # Return the modified data frame
  return(df)
}

prep_te<- function(df, df_nonproband, df_noconsent, df_hostseqcancer, df_metrics, column_names, thresholds, type, export_filtered = TRUE, output_dir = NULL) {
  # Step 1: Rename alternative columns
  df <- rename_alt(df)

  # Step 2: Replace sample names to match the clinical data
  df <- replace_df_samples(df, type)

  # Step 3: Filter by nonproband, noconsent, and metrics criteria
  df <- filter_by_multiple_criteria(df, df_nonproband, df_noconsent, df_hostseqcancer, df_metrics, column_names, thresholds, type, export_filtered, output_dir)

  return(df)
}

add_gene_size <- function(df) {
  # Summarize the data to one row per gene
  summarized_df <-df %>%
    filter(grepl("^chr[0-9XYM]+$", chrom)) %>% 
    group_by(name2) %>%
    summarize(
      #chrom = unique(chrom), # dont keep chrom bc some genes are in X and Y
      gene_start = min(txStart),
      gene_end = max(txEnd),
      gene_size = gene_end - gene_start,  # Use calculated gene_start and gene_end
      .groups = "drop" # Ungroup after summarizing
    )
  return(as.data.frame(summarized_df))
}

count_TE_occurrences<- function(te_expand, te_all, te_all_all, nohits_file, type="N") {
  # te is each te and only has samples with a te
  # te_all is sample summary and includes samples without TEs

  # total unique samples
  num_total_samples <- length(unique(te_all_all$sample)) 
  
  # total unique patients
  num_total_patients <- length(unique(te_all$base_sample))

  # Step 1: Count unique samples with at least one te
  num_samples_gt1 <- length(unique(te_expand$sample)) 
  
  # Step 2: Group by TE characteristics and count occurrences per sample
  te_counts <- te_expand %>%
    group_by(SV_chrom, SV_start, SV_end, SV_length, ALT) %>%
    summarize(unique_samples_in_group = n_distinct(sample), .groups = "drop" ) %>%
    ungroup()
  total_te <- nrow(te_counts)
  
  # single sample
  num_single_sample_TEs <- nrow(te_counts %>% filter(unique_samples_in_group == 1))
  perc_single_sample <- num_single_sample_TEs/total_te*100
  
  # Load and count nohits
  load(paste0(r_dir, nohits_file, ".RData"))
  num_no_hits <- nrow(nohits)
  
  # at least 1 te
  gt1 <- num_samples_gt1/ num_total_samples*100
  
  # TP53 counts
  # Count patients (base_sample level)
  tp53_counts_patients <- te_all %>%
    select(base_sample, TP53_status) %>%
    distinct() %>%
    count(TP53_status)
  
  # Count samples (sample level) - for samples actually in the analysis
  tp53_counts_samples <- te_all %>%
    select(sample, TP53_status) %>%
    distinct() %>%
    count(TP53_status)

  # number of patients with multiple samples
  if (type=="T") {
    te_sample_counts <- te_expand %>%
      group_by(base_sample) %>%
      summarise(unique_samples_for_patient = n_distinct(sample), .groups = "drop") %>%
      filter(unique_samples_for_patient > 1) 
    
    # Report the number of such groups and the average count
    num_groups_gt1 <- nrow(te_sample_counts)
    average_count <- mean(te_sample_counts$unique_samples_for_patient)
  } 
  
  # Print the results
  cat("Total number of patients:", num_total_patients, "\n")
  cat("Total number of samples:", num_total_samples, "\n")
  cat("Number of TP53 patients: \n")
  print(tp53_counts_patients)
  cat("Number of TP53 samples: \n")
  print(tp53_counts_samples)
  cat("Number of samples with 0 TEs:", num_no_hits, "\n")
  cat("Percent of samples with at least 1 TE:", gt1, "%\n")
  cat("Number of TEs in only one sample:", num_single_sample_TEs, "(", perc_single_sample, "%)\n")
  cat("Spread of number of TEs shared by samples (ie on average the same TE is in 23 samples:", "\n")
  print(summary(te_counts$unique_samples_in_group))
  cat("Spread of ALUs by sample:", "\n")
  print(summary(te_all$ALU))
  cat("Spread of LINE1s by sample:", "\n")
  print(summary(te_all$LINE1))
  print(summary(te_all$SVA))
}

plot_te_counts_summary <- function(df, y_lab, log_scale=FALSE, breaks=NULL) {
  # Ensure required columns exist
  required_columns <- c("LINE1", "ALU", "SVA")
  if (!all(required_columns %in% colnames(df))) {
    stop("The dataframe must contain the columns: LINE1, ALU, and SVA")
  }
  
  # Pivot data to long format
  data_long <- df %>%
    pivot_longer(cols = all_of(required_columns), names_to = "Element", values_to = "Value")
  
  # Create the box plot with points overlaid
  p<- ggplot(data_long, aes(x = Element, y = Value, fill = Element)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +  # Box plot with reduced opacity
    geom_jitter(position = position_jitterdodge(jitter.width = 0.2), color = "black", size = 1.5) +
    scale_fill_manual(values = colour_palette_3) +  # Apply custom fill colors
    labs(x = "Repeat type", y = y_lab) +
    theme(legend.position = "none") 
  
  if (log_scale) {
    p <- p + scale_y_continuous(
      trans = scales::log1p_trans(),
      breaks = if (!is.null(breaks)) breaks else waiver()
    )
  }
  
  return(p)
}

plot_te_counts_unique<- function(te) {
  # Step 1: Compute unique_samples_in_group
  te_counts <- te %>%
    group_by(SV_chrom, SV_start, SV_end, SV_length, ALT) %>%
    summarize(unique_samples_in_group = n_distinct(sample)) %>%
    ungroup()
  
  # Step 2: Calculate the proportion of samples for each group
  total_samples <- length(unique(te$sample))  # Total number of unique samples
  te_counts <- te_counts %>%
    mutate(proportion_of_samples = unique_samples_in_group / total_samples)  # Calculate proportion
  
  # Step 3: Create the plot
  plot <- ggplot(te_counts, aes(x = ALT, y = proportion_of_samples, fill = ALT)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +  # Boxplot with reduced opacity
    geom_jitter(position = position_jitterdodge(jitter.width = 0.2), color = "black", size = 1.5) +  # Add jittered points
    labs(
      x = "Repeat type",
      y = "Proportion of samples with element"
    ) +
    scale_y_continuous(breaks = scales::pretty_breaks()) +
    theme(legend.position = "none",  # Remove legend
          axis.text.x = element_text(angle = 45, hjust = 1))  # Adjust x-axis labels
  
  return(plot)
}

stacked_bar_plot_num_samples_4 <- function(te, thresholds = c(1, 3, 10, 50)) {
  # Validate thresholds length
  if (length(thresholds) != 4) {
    stop("Thresholds must have exactly 4 elements.")
  }
  print(nrow(te))
  # Calculate the number of unique samples for each TE
  te <- te %>%
    group_by(SV_chrom, SV_start, SV_end, SV_length, ALT) %>%
    summarize(unique_samples_in_group = n_distinct(sample)) %>%
    ungroup()
  print(nrow(te))
  
  # Categorize TEs based on the number of samples they appear in
  te_counts <- te %>%
    mutate(sample_count_category = case_when(
      unique_samples_in_group == thresholds[1] ~ paste0(thresholds[1], ' sample'),
      unique_samples_in_group > thresholds[1] & unique_samples_in_group <= thresholds[2] ~ paste0(thresholds[1] + 1, '-', thresholds[2], ' samples'),
      unique_samples_in_group > thresholds[2] & unique_samples_in_group <= thresholds[3] ~ paste0(thresholds[2] + 1, '-', thresholds[3], ' samples'),
      unique_samples_in_group > thresholds[3] & unique_samples_in_group <= thresholds[4] ~ paste0(thresholds[3] + 1, '-', thresholds[4], ' samples'),
      unique_samples_in_group > thresholds[4] ~ paste0('>', thresholds[4], ' samples')
    ))
  
  # Order the sample count categories
  category_order <- c(paste0(thresholds[1], ' sample'), 
                      paste0(thresholds[1] + 1, '-', thresholds[2], ' samples'), 
                      paste0(thresholds[2] + 1, '-', thresholds[3], ' samples'), 
                      paste0(thresholds[3] + 1, '-', thresholds[4], ' samples'), 
                      paste0('>', thresholds[4], ' samples'))
  
  te_counts <- te_counts %>%
    mutate(sample_count_category = factor(sample_count_category, levels = category_order))
  
  print(head(te_counts))
  
  # Summarize data for plotting
  te_summary <- te_counts %>%
    group_by(ALT, sample_count_category) %>%
    summarize(count = n()) %>%
    ungroup()
  print(head(te_summary)) 
  
  # Order the sample count categories
  category_order <- c(paste0(thresholds[1], ' sample'), 
                      paste0(thresholds[1] + 1, '-', thresholds[2], ' samples'), 
                      paste0(thresholds[2] + 1, '-', thresholds[3], ' samples'), 
                      paste0(thresholds[3] + 1, '-', thresholds[4], ' samples'), 
                      paste0('>', thresholds[4], ' samples'))
  
  te_summary <- te_summary %>%
    mutate(sample_count_category = factor(sample_count_category, levels = rev(category_order)))
  
  # Plot
  p <- ggplot(te_summary, aes(x = ALT, y = count, fill = sample_count_category)) +
    geom_bar(stat = "identity") +
    labs(x = "TE type",  
         y = "Number of TEs") +
    scale_fill_manual(values = colour_palette_5,
                      name = "Samples with TE") 
  
  # Plot with y-axis normalized to proportions
  te_summary <- te_summary %>%
    group_by(ALT) %>%
    mutate(proportion = count / sum(count)) %>%
    ungroup()
  
  p2<- ggplot(te_summary, aes(x = ALT, y = proportion, fill = sample_count_category)) +
    geom_bar(stat = "identity") +
    labs(x = "Repeat type",
         y = "Proportion of repeat type") +
    scale_fill_manual(values = colour_palette_5,
                      name = "Samples with the same repeat")
  
 print(p) 
 print(p2)
}

stacked_bar_plot_num_samples <- function(te, thresholds = c(1, 10, 50)) {
  # Validate thresholds length
  if (length(thresholds) != 3) {
    stop("Thresholds must have exactly 3 elements.")
  }
  # Calculate the number of unique samples for each TE
  te <- te %>%
    group_by(SV_chrom, SV_start, SV_end, SV_length, ALT) %>%
    summarize(unique_samples_in_group = n_distinct(sample)) %>%
    ungroup()
  
  # Categorize TEs based on the number of samples they appear in
  te_counts <- te %>%
    mutate(sample_count_category = case_when(
      unique_samples_in_group == thresholds[1] ~ paste0(thresholds[1], ' sample'),
      unique_samples_in_group > thresholds[1] & unique_samples_in_group <= thresholds[2] ~ paste0(thresholds[1] + 1, '-', thresholds[2], ' samples'),
      unique_samples_in_group > thresholds[2] & unique_samples_in_group <= thresholds[3] ~ paste0(thresholds[2] + 1, '-', thresholds[3], ' samples'),
      unique_samples_in_group > thresholds[3] ~ paste0('>', thresholds[3], ' samples')
    ))
  
  # Order the sample count categories
  category_order <- c(paste0(thresholds[1], ' sample'), 
                      paste0(thresholds[1] + 1, '-', thresholds[2], ' samples'), 
                      paste0(thresholds[2] + 1, '-', thresholds[3], ' samples'), 
                      paste0('>', thresholds[3], ' samples'))
  
  te_counts <- te_counts %>%
    mutate(sample_count_category = factor(sample_count_category, levels = category_order))
  
  # Summarize data for plotting
  te_summary <- te_counts %>%
    group_by(ALT, sample_count_category) %>%
    summarize(count = n()) %>%
    ungroup()
  
  # Order the sample count categories
  category_order <- c(paste0(thresholds[1], ' sample'), 
                      paste0(thresholds[1] + 1, '-', thresholds[2], ' samples'), 
                      paste0(thresholds[2] + 1, '-', thresholds[3], ' samples'), 
                      paste0('>', thresholds[3], ' samples'))
  
  te_summary <- te_summary %>%
    mutate(sample_count_category = factor(sample_count_category, levels = rev(category_order)))
  
  # Plot
  p <- ggplot(te_summary, aes(x = ALT, y = count, fill = sample_count_category)) +
    geom_bar(stat = "identity") +
    labs(x = "Repeat type",  
         y = "Number of repeats") +
    scale_fill_manual(values = colour_palette_4_seq,
                      name = "Samples with the same repeat") 
  
  # Plot with y-axis normalized to proportions
  te_summary <- te_summary %>%
    group_by(ALT) %>%
    mutate(proportion = count / sum(count)) %>%
    ungroup()
  
  p2<- ggplot(te_summary, aes(x = ALT, y = proportion, fill = sample_count_category)) +
    geom_bar(stat = "identity") +
    labs(x = "Repeat type",
         y = "Proportion of repeat type") +
    scale_fill_manual(values = colour_palette_4_seq,
                      name = "Samples with the same repeat")
  
  print(p) 
  print(p2)
}

compare_event_proportions <- function(df, group_column, total_column, x_lab = "Group") {
  # Ensure the columns exist
  if (!(group_column %in% colnames(df)) || !(total_column %in% colnames(df))) {
    stop("Specified columns do not exist in the dataframe.")
  }
  
  # Step 1: Create a binary 'event' column
  df <- df %>%
    mutate(event = ifelse(.data[[total_column]] > 0, "Yes", "No"))
  
  # Step 2: Summarize the data for proportions
  proportions <- df %>%
    group_by(!!sym(group_column), event) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(!!sym(group_column)) %>%
    mutate(proportion = count / sum(count))
  
  # Step 3: Plot proportions
  p <- ggplot(proportions, aes(x = .data[[group_column]], y = proportion, fill = event)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = c("No" = "#0080A3", "Yes" = "#B1D586"), 
                      labels = c("No" = "No repeats", "Yes" = "At least one repeat")) +
    labs(
      x = x_lab,
      y = "Proportion",
      fill= NULL
    ) 
  
  # Step 4: Prepare contingency table for Fisher's Exact Test
  contingency_table <- df %>%
    group_by(!!sym(group_column), event) %>%
    summarise(count = n(), .groups = "drop") %>%
    pivot_wider(names_from = event, values_from = count, values_fill = 0) %>%
    arrange(desc(!!sym(group_column))) # Flip group order
  
  # Convert to matrix
  contingency_matrix <- as.matrix(contingency_table[, c("Yes", "No")])
  rownames(contingency_matrix) <- contingency_table[[group_column]]
  
  # Perform Fisher's Exact Test
  fisher_result <- fisher.test(contingency_matrix)
  
  # Print outputs
  cat("\nContingency Table:\n")
  print(contingency_matrix)
  
  cat("\nFisher's Exact Test Results:\n")
  print(fisher_result)
  
  # Print the plot
  return(p)
}

plot_proportions <- function(data) {
  # Iterate over each row and generate a plot
  plots <- lapply(1:nrow(data), function(i) {
    # Extract row data
    row <- data[i, ]
    
    # Prepare data for ggplot
    proportions <- data.frame(
      tp53_status = c("WT", "Mutant"),
      proportion = c(row$perc_TP53_WT, row$perc_TP53_Mutant),
      event = c("With repeat", "With repeat") # Both values are for "With repeat"
    )
    
    # Add absent proportions
    proportions <- rbind(
      proportions,
      data.frame(
        tp53_status = c("WT", "Mutant"),
        proportion = 100 - proportions$proportion, # Remaining proportion is "Absent"
        event = c("Absent", "Absent")
      )
    )
    
    # Generate the plot
    ggplot(proportions, aes(x = tp53_status, y = proportion, fill = event)) +
      geom_bar(stat = "identity", position = "stack") +
      scale_fill_manual(values = c("Absent" = "#B1D586", "With repeat" = "#0080A3"),
                        labels = c("Absent" = "Absent", "With repeat" = "With repeat")) +
      labs(
        x = "TP53 Status",
        y = "Proportion",
        fill = NULL,
        title = paste("SV Chromosome:", row$SV_chrom, 
                      "| Start:", row$SV_start, "| End:", row$SV_end)
      ) 
  })
  
  # Return the list of plots
  return(plots)
}

plot_proportions<- function(data) {
  # Iterate over each row and generate a plot
  plots <- lapply(1:nrow(data), function(i) {
    # Extract row data
    row <- data[i, ]
    
    # Prepare data for ggplot
    proportions <- data.frame(
      tp53_status = c("WT", "Mutant"),
      proportion = c(row$perc_TP53_WT, row$perc_TP53_Mutant)
    )
    
    x_tp53 <- expression("Germline " * italic("TP53") * " status")
    
    # Generate the plot
    ggplot(proportions, aes(x = tp53_status, y = proportion, fill = tp53_status)) +
      geom_bar(stat = "identity", position = "dodge", fill = "#0080A3") + # Blue color for bars
      labs(
        x = x_tp53,
        y = "Proportion of group with repeat", 
        title = paste("SV Chromosome:", row$SV_chrom, 
                      "\nStart:", row$SV_start, "| End:", row$SV_end)
      ) 
  })
  
  # Return the list of plots
  return(plots)
}

plot_venn <- function(df) {
  # Sum the columns to get total counts
  num_ins <- sum(df$num_ins)
  num_xtea <- sum(df$num_xtea)
  num_melt <- sum(df$num_melt)
  
  num_ins_xtea <- sum(df$num_ins_xtea)
  num_ins_melt <- sum(df$num_ins_melt)
  num_xtea_melt <- sum(df$num_xtea_melt)
  num_ins_xtea_melt <- sum(df$num_ins_xtea_melt)
  
  # Print the totals and overlaps
  cat("INS:", num_ins, "\n")
  cat("XTEA:", num_xtea, "\n")
  cat("MELT:", num_melt, "\n")
  cat("INS & XTEA:", num_ins_xtea, "\n")
  cat("INS & MELT:", num_ins_melt, "\n")
  cat("XTEA & MELT:", num_xtea_melt, "\n")
  cat("INS & XTEA & MELT:", num_ins_xtea_melt, "\n")
  
  # Prepare data for Euler diagram
  fit <- euler(c(
    INS = num_ins,
    XTEA = num_xtea,
    MELT = num_melt,
    "INS&XTEA" = num_ins_xtea,
    "INS&MELT" = num_ins_melt,
    "XTEA&MELT" = num_xtea_melt,
    "INS&XTEA&MELT" = num_ins_xtea_melt
  ))
  
  # Plot the Euler diagram
  p <- plot(
    fit,
    fills = list(fill = colour_palette_3, alpha = 0.5),
    quantities = list(outside=TRUE, size=8))
  
  print(p)
}

filter_and_count <- function(df, filter_element = NA, chr = NA) {
  # Identify sample-level metadata columns (columns that have same value for each sample)
  # These are columns we want to preserve in the count matrix
  sample_cols <- c("sample", "base_sample", "predicted_ancestry_thres", "mapped_label")
  metadata_cols <- intersect(sample_cols, colnames(df))

  # Get unique sample metadata
  complete_data <- df %>%
    dplyr::select(any_of(metadata_cols)) %>%
    distinct()

  # Filter based on 'chr' if applicable
  if (!is.na(chr)) {
    df <- df %>%
      filter(SV_chrom == chr)
  }

  # Perform filtering TE type and count occurence by sample
  if (!is.na(filter_element)) {
    filtered_data <- df %>%
      filter(ALT == filter_element) %>%
      group_by(sample) %>%
      summarise(count = n(), .groups = 'drop')
  } else {
    filtered_data <- df %>%
      group_by(sample) %>%
      summarise(count = n(), .groups = 'drop')
  }

  # Left join to include all samples, with zeros for missing counts
  result <- complete_data %>%
    left_join(filtered_data, by = "sample") %>%
    replace_na(list(count = 0))

  return(result)
}

process_all_combinations <- function(df) {
  # Define the filter elements and chromosome values
  filter_elements <- c(NA, "LINE1", "ALU", "SVA")
  chromosomes <- c(NA, 1:22, "X", "Y", "MT")

  # Initialize an empty list to store results
  results_list <- list()

  # Iterate over each filter_element and chromosome combination
  for (filter_element in filter_elements) {
    for (chr in chromosomes) {
      # Run the filter_and_count function
      result <- filter_and_count(df, filter_element, chr)
      # Create an identifier column for each combination
      combination_label <- if (is.na(chr) && is.na(filter_element)) {
        "total"
      } else if (is.na(chr)) {
        filter_element
      } else if (is.na(filter_element)) {
        paste0("chr", chr)
      } else {
        paste0("chr", chr, "_", filter_element)
      }
      # Safety check before adding combination column
      if (nrow(result) == 0) {
        cat("WARNING: Empty result for combination:", combination_label, "\n")
        # Create properly structured empty result
        result <- data.frame(sample = character(0), count = integer(0), combination = character(0), stringsAsFactors = FALSE)
      } else {
        # save result under specific label
        result$combination <- combination_label
      }

      # Append the result to the list
      results_list <- append(results_list, list(result))
    }
  }

  # Combine all dataframes into one
  combined_result <- bind_rows(results_list)

  # Identify metadata columns (non-count columns to preserve)
  metadata_cols <- setdiff(colnames(combined_result), c("count", "combination"))

  # Pivot the combined result to have combinations as columns, preserving metadata
  wide_result <- combined_result %>%
    pivot_wider(
      id_cols = all_of(metadata_cols),
      names_from = combination,
      values_from = count,
      values_fill = list(count = 0)
    )

  return(wide_result)
}

split_by_gene <- function(df){
  df_processed <- df %>%
    separate_rows(Gene_name, sep = ";") %>%
    filter(Gene_name != "")  # Remove rows where Gene_name is empty
  
  return(df_processed) 
}

# merge with clinical
merge_dfs <- function(df_te, df_info, include_all_x = FALSE, print_info = TRUE, dataset_name = "unknown") {
  # If merging with ancestry_tumour and tumor data has base_sample, remove _N from ancestry
  is_ancestry_tumour <- dataset_name == "ancestry_tumour"

  if (is_ancestry_tumour && "base_sample" %in% colnames(df_te)) {
    # Remove _N from ancestry samples to match tumor base_sample
    df_info <- df_info %>%
      mutate(base_sample = gsub("_N$", "", sample)) %>%
      select(-sample)

    # Merge on base_sample
    all <- merge(df_te, df_info, by = "base_sample", all.x = include_all_x)
    unmatched_merged <- unique(df_te[!df_te$base_sample %in% df_info$base_sample, ]$sample)
  } else {
    # Standard merge by sample
    all <- merge(df_te, df_info, by = "sample", all.x = include_all_x)
    unmatched_merged <- unique(df_te[!df_te$sample %in% df_info$sample, ]$sample)
  }

  if (print_info && length(unmatched_merged) > 0) {
    cat("Samples missing", dataset_name, "info:\n")
    print(unmatched_merged)
  }

  return(all)
}

filter_age <- function(df) {
  # Calculate the number of unique samples before filtering
  original_num_samples <- length(unique(df$sample))
  
  # Apply conditional filtering and keep rows where age_to_use is NA
  df_filt <- df %>%
    dplyr::mutate(age_to_use = ifelse(tumor_type == "U", age_at_enrollment, age_at_diagnosis)) %>%
    dplyr::filter(age_to_use < 30 * 365 | is.na(age_to_use))  # Keep rows where age_to_use is NA
  
  # Print samples with NA for age_to_use
  samples_with_na <- df %>%
    dplyr::mutate(age_to_use = ifelse(tumor_type == "U", age_at_enrollment, age_at_diagnosis)) %>%
    filter(is.na(age_to_use))
  #cat("Samples with NA for age_to_use:\n")
  #print(samples_with_na$sample)
  
  # Calculate the number of unique samples after filtering
  filtered_num_samples <- length(unique(df_filt$sample))
  
  # Calculate the number of unique samples filtered out
  num_filtered_out <- original_num_samples - filtered_num_samples
  cat("Number of samples filtered out for age:", num_filtered_out, "\n")
  
  return(df_filt)
}

six_threshold <- function(df){
  six_years_days <- 6*365 # six years in days
  df$six_years <- ifelse(df$age_at_diagnosis < six_years_days, "Before 6", "After 6")
  return(df)
}

add_nohit_samples <- function(original_data, nohits) {
  # reformat
  new_samples <- as.character(nohits$V1)

  # Print how many samples with no insertions are being added
  cat("Adding", length(new_samples), "samples with no TE insertions\n")

  # Create a new data frame for the new samples with zero counts for TE columns
  new_data <- tibble(sample = new_samples)

  # Add columns with appropriate types for each column in the original data (excluding 'sample')
  for (col in names(original_data)) {
    if (col == "sample") next  # Skip sample column, already added

    # Get the column type from original data
    col_class <- class(original_data[[col]])[1]

    if (col_class %in% c("character", "factor")) {
      # For character/factor columns, use NA
      new_data[[col]] <- NA_character_
    } else if (col_class %in% c("numeric", "double", "integer")) {
      # For numeric columns, use 0
      new_data[[col]] <- 0
    } else {
      # Default to NA for other types
      new_data[[col]] <- NA
    }
  }

  # Combine the original data with the new data
  updated_data <- bind_rows(original_data, new_data)

  return(updated_data)
}

remove_duplicates <- function(df) {
  # Count the number of rows before removal
  initial_rows <- nrow(df)
  
  # Remove duplicate rows
  df_cleaned <- df %>% distinct()
  
  # Count the number of rows after removal
  final_rows <- nrow(df_cleaned)
  
  # Calculate and print the number of rows removed
  rows_removed <- initial_rows - final_rows
  cat("Number of duplicate rows removed:", rows_removed, "\n")
  
  # Return the cleaned data frame
  return(df_cleaned)
}

# Function to process and sort TE data
process_te_data_tumour <- function(te_raw, te_germline=NULL, clinical, complete_samples, metrics_list, ancestry,
                            apply_filter_common = TRUE, rare_gnomad = 3, rare_hostseq = 3,
                            split_by_gene = FALSE,
                            apply_process_combinations = TRUE,
                            select_samples_split = FALSE,
                            nohits_prefix = "default",
                            nohits_output_dir = NULL
                            ) {
  
  # Optionally filter common transposable elements
  if (apply_filter_common) {
    te_filtered <- filter_common_hostseq_tumour_te(te_raw, te_germline, rare_gnomad, rare_hostseq)
  } else {
    te_filtered <- te_raw # dont filter common
  }
  
  # Optionally split by gene so expand version has one row per gene
  if (split_by_gene) {
    te_split_by_gene <- split_by_gene(te_filtered)
  } else{
    te_split_by_gene <- te_filtered
  } 
  
  # Optionally process all combinations
  if (apply_process_combinations) {
    te_processed <- process_all_combinations(te_split_by_gene)
    # Calculate nohits as samples in complete_samples but not in te_processed
    complete_samples_unique <- unique(complete_samples$V1)
    te_processed_samples_unique <- unique(te_processed$sample)
    nohits_samples <- setdiff(complete_samples_unique, te_processed_samples_unique)
    nohits <- data.table(V1 = nohits_samples)
    te_processed <- add_nohit_samples(te_processed, nohits)
  } else {
    te_processed <- te_split_by_gene
  }
  
  # merge ancestry
  te_with_ancestry <- merge_dfs(te_processed, ancestry, include_all_x =TRUE, print_info=TRUE, dataset_name="ancestry_tumour")

  # merge metrics
  te_clinical <- merge_dfs(te_with_ancestry, metrics_list, include_all_x =TRUE, print_info=TRUE, dataset_name="metrics")

  # Merge with clinical data
  te_clinical <- merge_dfs(te_clinical, clinical, include_all_x = FALSE, print_info=TRUE, dataset_name="clinical")
  
  # Remove samples where age at diagnosis is greater than 30
  te_filt <- filter_age(te_clinical)  # Assuming filter_age can accept a max_age parameter
  
  # remove duplicate entries
  te_format <- remove_duplicates(te_filt)
  
  # Set factor order
  te_format$TP53_status<- factor(te_format$TP53_status, levels = c("WT", "Mutant"))
  
  # Sort the data into categories - exclude Taylor from aff and lfs datasets
  te_all_df <- te_format  # All samples including Taylor
  te_aff_df <- te_format %>% filter(tumor_type != "U" & cohort != "Taylor")  # Cancer samples excluding Taylor
  te_lfs_df <- te_format %>%
    filter(TP53_status == "Mutant" & cohort != "Taylor") %>%
    mutate(Cancer = ifelse(tumor_type != "U", "Affected", "Unaffected"))  # LFS samples excluding Taylor
  te_kics_df <- te_format %>% filter(cohort == "KiCS")  # KiCS cohort samples
  te_taylor_df <- te_format %>% filter(cohort == "Taylor")  # Taylor-only dataset

  # Optionally select 1 sample for each patient
  if (select_samples_split) {
    te_all_samples <- select_samples(te_all_df, select_samples_split)
    te_aff_samples <- select_samples(te_aff_df, select_samples_split)
    te_lfs_samples <- select_samples(te_lfs_df, select_samples_split)
    te_kics_samples <- select_samples(te_kics_df, select_samples_split)
    te_taylor_samples <- select_samples(te_taylor_df, select_samples_split)

    result <- list(
      te_all_selected = te_all_samples$selected_samples,
      te_all_all = te_all_samples$all_samples,
      te_aff_selected = te_aff_samples$selected_samples,
      te_aff_all = te_aff_samples$all_samples,
      te_lfs_selected = te_lfs_samples$selected_samples,
      te_lfs_all = te_lfs_samples$all_samples,
      te_kics_selected = te_kics_samples$selected_samples,
      te_kics_all = te_kics_samples$all_samples,
      te_taylor_selected = te_taylor_samples$selected_samples,
      te_taylor_all = te_taylor_samples$all_samples
    )

  } else {
    # Combine the processed data into a list for easy access
    result <- list(
      te_all = te_all_df,
      te_aff = te_aff_df,
      te_lfs = te_lfs_df,
      te_kics = te_kics_df,
      te_taylor = te_taylor_df
    )
  }
  
  # Calculate and save nohits lists for each result dataframe
  if (apply_process_combinations) {
    for (name in names(result)) {
      df <- result[[name]]
      if ("total" %in% colnames(df)) {
        nohits_samples <- df[df$total == 0, ]$sample
        nohits <- data.table(V1 = nohits_samples)
        save(nohits, file = paste0(r_dir, "nohits_", nohits_prefix, "_", name, "_t.RData"))
        # Save only te_aff_selected to nohits_output_dir with simplified name
        if (!is.null(nohits_output_dir) && name == "te_aff_selected") {
          write.table(nohits, file = paste0(nohits_output_dir, "nohits_te_aff_t.csv"), row.names = FALSE, col.names = FALSE, sep = ",", quote = FALSE)
        }
      }
    }
  }
  
  return(result)
}

process_te_data_germline <- function(te_raw, clinical, metrics_list, ancestry_list, 
                            apply_filter_common = TRUE, rare_gnomad = 3, rare_hostseq = 3, 
                            split_by_gene = FALSE,
                            apply_process_combinations = TRUE 
                            ) {
  # Optionally filter common transposable elements
  if (apply_filter_common) {
    te_filtered <- filter_common_hostseq_germline_te(te_raw, rare_gnomad, rare_hostseq)
  } else {
    te_filtered <- te_raw
  }

  # Ensure cohort_hs column exists (needed for filtering later)
  if (!"cohort_hs" %in% colnames(te_filtered)) {
    te_filtered$cohort_hs <- ifelse(grepl("^HS_", te_filtered$sample), "HostSeq", "Other")
  }
  
  # Optionally split by gene
  if (split_by_gene) {
    te_split_by_gene <- split_by_gene(te_filtered)
  } else{
    te_split_by_gene <- te_filtered
  } 
  
  # Optionally process all combinations
  # dont need to worry about no hits because it germline
  if (apply_process_combinations) {
    te_processed <- process_all_combinations(te_split_by_gene)
  } else {
    te_processed <- te_split_by_gene
  }

  # merge ancestry
  te_clinical <- merge_dfs(te_processed, ancestry_list, include_all_x =TRUE, print_info=TRUE, dataset_name="ancestry")

  # merge metrics
  te_clinical <- merge_dfs(te_clinical, metrics_list, include_all_x =TRUE, print_info=TRUE, dataset_name="metrics")

  # Merge with clinical data
  te_clinical <- merge_dfs(te_clinical, clinical, include_all_x = FALSE, print_info=TRUE, dataset_name="clinical")
  
  # Remove samples where age at diagnosis is greater than 30
  te_filt <- filter_age(te_clinical)  # Assuming filter_age can accept a max_age parameter
  
  # remove duplicate entries
  te_format <- remove_duplicates(te_filt)

  # Set factor order
  te_format$TP53_status<- factor(te_format$TP53_status, levels = c("WT", "Mutant"))

  # Sort the data into categories - exclude Taylor and HostSeq from aff, lfs, kics datasets
  te_all_df <- te_format  # All samples including Taylor and HostSeq
  te_aff_df <- te_format %>% filter(tumor_type != "U" & cohort != "Taylor" & cohort != "HostSeq")  # Cancer samples excluding Taylor and HostSeq
  te_aff_unaff_df <- te_format %>% filter(cohort != "Taylor" & cohort != "HostSeq")  # Affected + Unaffected (LFS + KICS only, excluding Taylor and HostSeq)
  te_lfs_df <- te_format %>%
    filter(TP53_status == "Mutant" | cohort %in% c("LFS", "Nick", "SJ")) %>%
    mutate(Cancer = ifelse(tumor_type != "U", "Affected", "Unaffected"))  # LFS-related samples: TP53 mutant OR from LFS/Nick/SJ cohorts (includes Taylor and HostSeq if they match criteria)
  te_kics_df <- te_format %>% filter(cohort == "KiCS")  # KiCS cohort samples only
  te_taylor_df <- te_format %>% filter(cohort == "Taylor")  # Taylor-only dataset
  te_hostseq_df <- te_format %>% filter(cohort == "HostSeq")  # HostSeq analysis group only
  te_kics_hostseq_df <- te_format %>% filter(cohort %in% c("KiCS", "HostSeq"))  # KiCS + HostSeq combined

  # Combine the processed data into a list for easy access
  result <- list(
    te_all = te_all_df,
    te_aff = te_aff_df,
    te_aff_unaff = te_aff_unaff_df,
    te_lfs = te_lfs_df,
    te_kics = te_kics_df,
    te_taylor = te_taylor_df,
    te_hostseq = te_hostseq_df,
    te_kics_hostseq = te_kics_hostseq_df
  )
  
  return(result)
}

filter_common_TEs <- function(te, rare_threshold_percentage) {
  # Common if gnomAD AF is > threshold or TE > threshold in this dataset
  rare_threshold <- rare_threshold_percentage/100 # convert to decimal
  cat("Rare treshold:", rare_threshold, "\n")
  
  # Count unique samples
  total_samples <- length(unique(te$sample))
  cat("Total # samples:", total_samples, "\n")
  
  # Group by TE characteristics and calculate thresholds for filtering
  te_with_flags <- te %>%
    group_by(SV_chrom, SV_start, SV_end, SV_length, ALT) %>%
    mutate(unique_samples_in_group = n_distinct(sample)) %>%
    ungroup() %>%
    # Mark as common if frequency in samples exceeds threshold
    mutate(
      is_common = unique_samples_in_group / total_samples >= rare_threshold,
      exceeds_AF_threshold = GRPMAX_AF >= rare_threshold_percentage
    )  # Mark for GRPMAX_AF filtering
  
  # Calculate the counts for each filter condition
  num_filtered_by_AF <- sum(te_with_flags$exceeds_AF_threshold)
  num_filtered_by_common <- sum(te_with_flags$is_common)
  
  # Calculate the overlap between the two filters
  num_overlap_filters <- sum(te_with_flags$is_common & te_with_flags$exceeds_AF_threshold)
  
  # Filter out TEs that are either common or exceed the AF threshold
  uncommon_TEs <- te_with_flags %>%
    filter(is_common == FALSE, exceeds_AF_threshold == FALSE) %>%
    # Remove helper columns
    dplyr::select(-unique_samples_in_group, -is_common, -exceeds_AF_threshold)
  
  # Calculate the percentage filtered by each condition and by the overlap
  total_TE <- nrow(te)
  percent_filtered_by_AF <- (num_filtered_by_AF / total_TE) * 100
  percent_filtered_by_common <- (num_filtered_by_common / total_TE) * 100
  percent_overlap_filters <- (num_overlap_filters / total_TE) * 100
  
  # Print messages
  cat("Filtered by GRPMAX_AF:", num_filtered_by_AF, "(", round(percent_filtered_by_AF, 2), "%)\n")
  cat("Filtered by is_common:", num_filtered_by_common, "(", round(percent_filtered_by_common, 2), "%)\n")
  cat("Filtered by both (overlap):", num_overlap_filters, "(", round(percent_overlap_filters, 2), "%)\n")
  
  return(uncommon_TEs)
}


filter_common_hostseq_tumour_te <- function(te_df, te_germline, rare_gnomad_threshold, rare_hostseq_threshold) {
  # Convert to data.table
  setDT(te_df)

  # Convert rare thresholds to decimals
  rare_threshold_gnomad <- rare_gnomad_threshold / 100
  cat("Rare gnomad threshold:", rare_threshold_gnomad, "\n")

  rare_threshold_hostseq <- rare_hostseq_threshold / 100
  cat("Rare hostseq threshold:", rare_threshold_hostseq, "\n")

  # Count unique TEs in te_df
  unique_te_df_count <- uniqueN(te_df, by = c("SV_chrom", "SV_start", "SV_end", "SV_length", "ALT"))
  cat("Total # of unique TEs in te_df:", unique_te_df_count, "\n")

  # Check if te_germline is NULL or empty
  if (is.null(te_germline) || nrow(te_germline) == 0) {
    cat("No germline data provided. Filtering only by GRPMAX_AF threshold.\n")
    te_df[, SV_chrom := as.character(SV_chrom)]
    filtered_te_df <- te_df[GRPMAX_AF < rare_threshold_gnomad]
    unique_filtered_te_df_count <- uniqueN(filtered_te_df, by = c("SV_chrom", "SV_start", "SV_end", "SV_length", "ALT"))
    cat("Total # of unique TEs after GRPMAX_AF filtering:", unique_filtered_te_df_count, "\n")
    return(filtered_te_df)
  }

  setDT(te_germline)

  # Ensure SV_chrom has the same data type in both datasets
  te_df[, SV_chrom := as.character(SV_chrom)]
  te_germline[, SV_chrom := as.character(SV_chrom)]

  # Add cohort column to te_germline
  te_germline[, cohort_hs := ifelse(grepl("^HS_", sample), "HostSeq", "Other")]

  # Count unique samples in te_germline for the "HostSeq" cohort
  total_samples_hostseq <- uniqueN(te_germline[cohort_hs == "HostSeq", sample])
  cat("Total # of samples in HostSeq cohort:", total_samples_hostseq, "\n")

  # If no HostSeq samples, only filter by gnomAD
  if (total_samples_hostseq == 0) {
    cat("No HostSeq samples found. Filtering only by GRPMAX_AF threshold.\n")
    filtered_te_df <- te_df[GRPMAX_AF < rare_threshold_gnomad]
    unique_filtered_te_df_count <- uniqueN(filtered_te_df, by = c("SV_chrom", "SV_start", "SV_end", "SV_length", "ALT"))
    cat("Total # of unique TEs after GRPMAX_AF filtering:", unique_filtered_te_df_count, "\n")
    return(filtered_te_df)
  }

  # Ensure ALT is character not list
  te_germline[, ALT := as.character(ALT)]

  # Group te_germline to calculate filtering flags
  te_with_flags_germline <- te_germline[
    cohort_hs == "HostSeq",
    .(
      unique_samples_in_group_hostseq = uniqueN(sample) # Count unique samples in HostSeq
    ),
    by = list(SV_chrom, SV_start, SV_end, SV_length, ALT)
  ]
  
  # Add the is_common_hostseq flag
  te_with_flags_germline[, is_common_hostseq := unique_samples_in_group_hostseq / total_samples_hostseq >= rare_threshold_hostseq]
  
  # Filter to only include common TEs
  te_with_flags_germline_common <- te_with_flags_germline[is_common_hostseq == TRUE]
  
  # Calculate end positions for te_germline and te_df
  te_with_flags_germline_common[, SV_end_calc := SV_start + SV_length]  # Calculate end position of TEs in te_germline
  te_df[, SV_end_calc := SV_start + SV_length]        # Calculate end position of TEs in te_df
  
  # Add a 100 bp buffer to the start and end positions for both te_df and te_with_flags_germline_common
  te_df[, `:=`(SV_start_buff = SV_start - 100, SV_end_buff = SV_end_calc + 100)]
  te_with_flags_germline_common[, `:=`(common_SV_start_buff = SV_start - 100, common_SV_end_buff = SV_end_calc + 100)]
  
  # Filter te_df by removing TEs that overlap with common TEs within 100 bp
  overlapping_indices <- te_df[
    te_with_flags_germline_common, 
    on = .(SV_chrom, 
           ALT,
           SV_start_buff <= common_SV_end_buff, 
           SV_end_buff >= common_SV_start_buff),
    which = TRUE
  ]
  
  # Remove NA indices and duplicates
  overlapping_indices <- unique(na.omit(overlapping_indices))

  filtered_by_common_te_df <- te_df[-overlapping_indices]
  
  # Filter te_df by GRPMAX_AF threshold
  filtered_te_df <- filtered_by_common_te_df[GRPMAX_AF < rare_threshold_gnomad]
  
  # Calculate unique TEs filtered by common overlap and GRPMAX threshold
  unique_filtered_by_common_te_df_count <- uniqueN(filtered_by_common_te_df, by = c("SV_chrom", "SV_start", "SV_end", "SV_length", "ALT"))
  unique_filtered_te_df_count <- uniqueN(filtered_te_df, by = c("SV_chrom", "SV_start", "SV_end", "SV_length", "ALT"))
  
  # Calculate stats
  num_filtered_by_common <- unique_te_df_count - unique_filtered_by_common_te_df_count
  num_filtered_by_AF <- unique_filtered_by_common_te_df_count - unique_filtered_te_df_count
  percent_filtered_by_common <- (num_filtered_by_common / unique_te_df_count) * 100
  percent_filtered_by_AF <- (num_filtered_by_AF / unique_te_df_count) * 100
  
  # Calculate total insertions removed
  total_te_insertions_input <- nrow(te_df)
  total_te_insertions_output <- nrow(filtered_te_df)
  total_insertions_removed <- total_te_insertions_input - total_te_insertions_output
  percent_insertions_removed <- (total_insertions_removed / total_te_insertions_input) * 100
  
  # Calculate insertions removed by each step
  insertions_after_common_filter <- nrow(filtered_by_common_te_df)
  insertions_removed_by_common <- total_te_insertions_input - insertions_after_common_filter
  insertions_removed_by_AF <- insertions_after_common_filter - total_te_insertions_output
  
  # Print detailed filtering stats
  cat("FILTERING SUMMARY (Unique TE Loci):\n")
  cat("Filtered by common TEs:", num_filtered_by_common, "(", round(percent_filtered_by_common, 2), "%)\n")
  cat("Filtered by GRPMAX_AF:", num_filtered_by_AF, "(", round(percent_filtered_by_AF, 2), "%)\n")
  cat("FILTERING SUMMARY (Total TE Insertions):\n")
  cat("Total insertions removed:", total_insertions_removed, "of", total_te_insertions_input, "(", round(percent_insertions_removed, 2), "%)\n")
  cat("Insertions removed by common overlap:", insertions_removed_by_common, "\n")
  cat("Insertions removed by GRPMAX_AF:", insertions_removed_by_AF, "\n")
  
  return(filtered_te_df)
}

filter_common_hostseq_germline_te <- function(te, rare_gnomad_threshold, rare_hostseq_threshold) {
  # make cohort column
  te <- te %>%
    mutate(cohort_hs = case_when(
      grepl("^HS_", sample) ~ "HostSeq",
      TRUE ~ "Other"
    ))
  
  # Common if gnomAD AF is > threshold or TE > threshold in hostseq 
  rare_threshold_gnomad <- rare_gnomad_threshold/100 # convert to decimal
  cat("Rare gnomad threshold:", rare_threshold_gnomad, "\n")
  
  rare_threshold_hostseq<- rare_hostseq_threshold/100 # convert to decimal
  cat("Rare hostseq threshold:", rare_threshold_hostseq, "\n")
  
  # Count unique samples
  total_samples <- length(unique(te$sample))
  cat("Total # samples:", total_samples, "\n")
  
  # total samples for hostseq - ONLY use filter group for calculations
  # Check if hostseq_group column exists (new method)
  if ("hostseq_group" %in% colnames(te)) {
    total_samples_hostseq <- length(unique(te$sample[te$hostseq_group == "filter"]))
    cat("Total # samples for HostSeq FILTER group:", total_samples_hostseq, "\n")
  } else {
    # Fallback to old method (all HostSeq)
    total_samples_hostseq <- length(unique(te$sample[te$cohort_hs == "HostSeq"]))
    cat("Total # samples for cohort 'HostSeq' (all):", total_samples_hostseq, "\n")
    cat("Warning: hostseq_group column not found, using all HostSeq samples\n")
  }

  # total tes
  total_te_df <- nrow(te)
  cat("Total # TEs called in df:", total_te_df, "\n")

  # Group by TE characteristics and calculate thresholds for filtering
  # ONLY count samples from HostSeq filter group
  te_with_flags <- te %>%
    group_by(SV_chrom, SV_start, SV_end, SV_length, ALT) %>%
    summarise(
      # Count unique samples from HostSeq FILTER cohort only
      unique_samples_in_group_hostseq = if ("hostseq_group" %in% colnames(cur_data())) {
        n_distinct(sample[hostseq_group == "filter"])
      } else {
        n_distinct(sample[cohort_hs == "HostSeq"])
      },

      # Calculate the filtering conditions
      exceeds_AF_threshold = any(GRPMAX_AF >= rare_threshold_gnomad, na.rm = TRUE),  # At least one sample exceeds the AF threshold in the group
      is_common_hostseq = unique_samples_in_group_hostseq / total_samples_hostseq >= rare_threshold_hostseq,  # Common condition based on HostSeq FILTER samples only
      .groups = "keep"
    ) 
  
  # total unique TEs
  total_unique_te <- nrow(te_with_flags)
  cat("Total # of unique TEs:", total_unique_te, "\n")
  
  # Calculate the counts for each filter condition (unique loci)
  num_filtered_by_AF <- sum(te_with_flags$exceeds_AF_threshold)
  num_filtered_by_common <- sum(te_with_flags$is_common_hostseq)
  num_overlap_filters <- sum(te_with_flags$is_common_hostseq & te_with_flags$exceeds_AF_threshold) # overlap between two filters
  
  # Filter out rows that are common or exceed the AF threshold
  uncommon_TEs <- te %>%
    left_join(te_with_flags %>% select(SV_chrom, SV_start, SV_end, SV_length, ALT, is_common_hostseq, exceeds_AF_threshold), 
              by = c("SV_chrom", "SV_start", "SV_end", "SV_length", "ALT")) %>%  # Join to keep the flags in the original df
    filter(is_common_hostseq == FALSE, exceeds_AF_threshold == FALSE) %>%  # Keep only uncommon and below threshold
    dplyr::select(-is_common_hostseq, -exceeds_AF_threshold)  # Remove helper columns
  
  # Calculate insertions removed (before vs after filtering)
  total_insertions_removed <- total_te_df - nrow(uncommon_TEs)
  
  # Calculate insertions removed by each filter condition
  insertions_filtered_by_AF <- te %>%
    left_join(te_with_flags %>% select(SV_chrom, SV_start, SV_end, SV_length, ALT, exceeds_AF_threshold), 
              by = c("SV_chrom", "SV_start", "SV_end", "SV_length", "ALT")) %>%
    filter(exceeds_AF_threshold == TRUE) %>% nrow()
  
  insertions_filtered_by_common <- te %>%
    left_join(te_with_flags %>% select(SV_chrom, SV_start, SV_end, SV_length, ALT, is_common_hostseq), 
              by = c("SV_chrom", "SV_start", "SV_end", "SV_length", "ALT")) %>%
    filter(is_common_hostseq == TRUE) %>% nrow()
  
  # Calculate the percentage filtered by each condition (unique loci)
  percent_filtered_by_AF <- (num_filtered_by_AF / total_unique_te) * 100
  percent_filtered_by_common <- (num_filtered_by_common / total_unique_te) * 100
  percent_overlap_filters <- (num_overlap_filters / total_unique_te) * 100
  
  # Calculate percentage of total insertions removed
  percent_insertions_removed <- (total_insertions_removed / total_te_df) * 100
  
  # Print messages with both unique loci and total insertions
  cat("FILTERING SUMMARY (Unique TE Loci):\n")
  cat("Filtered by common TEs:", num_filtered_by_common, "(", round(percent_filtered_by_common, 2), "%)\n")
  cat("Filtered by GRPMAX_AF:", num_filtered_by_AF, "(", round(percent_filtered_by_AF, 2), "%)\n")
  cat("Filtered by both (overlap):", num_overlap_filters, "(", round(percent_overlap_filters, 2), "%)\n")
  cat("FILTERING SUMMARY (Total TE Insertions):\n")
  cat("Total insertions removed:", total_insertions_removed, "of", total_te_df, "(", round(percent_insertions_removed, 2), "%)\n")
  cat("Insertions at loci with high GRPMAX_AF:", insertions_filtered_by_AF, "\n")
  cat("Insertions at common HostSeq loci:", insertions_filtered_by_common, "\n")

  # Remove HostSeq filter group samples (they were only used for frequency calculations)
  if ("hostseq_group" %in% colnames(uncommon_TEs)) {
    n_before <- nrow(uncommon_TEs)
    uncommon_TEs <- uncommon_TEs %>% filter(hostseq_group != "filter" | is.na(hostseq_group))
    n_removed <- n_before - nrow(uncommon_TEs)
    cat("Removed", n_removed, "HostSeq filter group TEs (kept analysis group)\n")
  }

  # Check missing samples in uncommon_TEs
  original_samples <- unique(te$sample)  # List of original samples
  remaining_samples <- unique(uncommon_TEs$sample)  # Samples left after filtering
  missing_samples <- setdiff(original_samples, remaining_samples)  # Find missing samples

  if (length(missing_samples) > 0) {
    cat("Samples lost after filtering common TEs:", missing_samples)
  } else {
    cat("No samples were lost after filtering common TE.\n")
  }
  
  return(uncommon_TEs)
}

select_samples <- function(df, select_samples_split = FALSE) {
  # Check if "base_sample.y" exists in column names and rename it
  if ("base_sample.y" %in% colnames(df)) {
    colnames(df)[colnames(df) == "base_sample.y"] <- "base_sample"
  }
  
  # Define ranking for lesion_type and disease_state
  lesion_rank <- c("primary" = 1, "metastasis" = 2, "relapse" = 2, "unknown" = 3)
  disease_rank <- c("initial" = 1, "progressive" = 2, "relapsed" = 3, "relapse" = 3, "unknown" = 4)
  
  # Handle NA and empty string as the lowest priority
  lesion_rank <- c(lesion_rank, "NA" = max(lesion_rank) + 1)
  disease_rank <- c(disease_rank, "NA" = max(disease_rank) + 1)
  
  # Assign lesion_rank
  df$lesion_rank <- ifelse(
    is.na(df$lesion_type) | df$lesion_type == "NA" | df$lesion_type == "" | !(df$lesion_type %in% names(lesion_rank)),
    lesion_rank["NA"],
    lesion_rank[df$lesion_type]
  )
  # Assign disease_rank
  df$disease_rank <- ifelse(
    is.na(df$disease_state) | df$disease_state == "NA" | df$disease_state == "" | !(df$disease_state %in% names(disease_rank)),
    disease_rank["NA"],
    disease_rank[df$disease_state]
  )
  
  # Identify invalid lesion_type and disease_state values and their rows
  # Only include values that are not NA, not empty, and not in the ranking
  lesion_mask <- !is.na(df$lesion_type) & df$lesion_type != "" & !(df$lesion_type %in% names(lesion_rank))
  disease_mask <- !is.na(df$disease_state) & df$disease_state != "" & !(df$disease_state %in% names(disease_rank))
  
  invalid_lesion_rows <- df[lesion_mask, ]
  invalid_disease_rows <- df[disease_mask, ]
  
  # Print warning for invalid lesion_type values only if there are actual problematic values
  if (nrow(invalid_lesion_rows) > 0) {
    unique_invalid_lesion <- unique(invalid_lesion_rows$lesion_type)
    if (length(unique_invalid_lesion) > 0) {
      warning("The following lesion_type values do not match any keys in lesion_rank: ", 
              paste(unique_invalid_lesion, collapse = ", "))
    }
  }
  
  # Print warning for invalid disease_state values only if there are actual problematic values
  if (nrow(invalid_disease_rows) > 0) {
    unique_invalid_disease <- unique(invalid_disease_rows$disease_state)
    if (length(unique_invalid_disease) > 0) {
      warning("The following disease_state values do not match any keys in disease_rank: ", 
              paste(unique_invalid_disease, collapse = ", "))
    }
  }
  
  # Select appropriate logic based on `select_samples_split`
  if (!select_samples_split) {
    # Version 1: Process the entire dataframe, grouped by `base_sample`
    df_selected <- df %>%
      group_by(base_sample) %>%
      arrange(age_at_enrollment, lesion_rank, disease_rank) %>%
      slice(1) %>%
      ungroup() %>%
      as.data.frame()  # Convert to base R data frame
    
    # Initialize warnings
    warnings <- list()
    
    unique_base_samples <- unique(df$base_sample)
    for (base in unique_base_samples) {
      patient_samples <- df[df$base_sample == base, ]
      selected_sample <- df_selected[df_selected$base_sample == base, ]
      
      # Check rankings
      max_lesion_rank <- min(patient_samples$lesion_rank)
      max_disease_rank <- min(patient_samples$disease_rank)
      
      # Debugging: Print problematic rows
     # cat("Checking base_sample:", base, "\n")
     # print(selected_sample)
     # print(patient_samples)
     # cat("Max lesion rank:", max(patient_samples$lesion_rank, na.rm = TRUE), "\n")
     # cat("Max disease rank:", max(patient_samples$disease_rank, na.rm = TRUE), "\n")
      
      if (selected_sample$lesion_rank > max_lesion_rank || selected_sample$disease_rank > max_disease_rank) {
        primary_initial_sample <- patient_samples %>%
          filter(lesion_type == "primary", disease_state == "initial") %>%
          slice(1)
        
        if (nrow(primary_initial_sample) > 0) {
          df_selected <- df_selected %>%
            filter(!(base_sample == base)) %>%
            bind_rows(primary_initial_sample) %>%
            as.data.frame()
          
          warning_message <- paste(
            "Warning: The selected sample for patient", base,
            "does not have the highest ranking. A primary lesion and initial disease state sample was selected instead.\n"
          )
          sample_details <- paste0(
            capture.output(print(patient_samples[, c("sample", "age_at_enrollment", "lesion_type", "disease_state")])),
            collapse = "\n"
          )
          warnings[[length(warnings) + 1]] <- paste0(warning_message, sample_details)
        } else {
          warning_message <- paste(
            "Warning: The selected sample for patient", base,
            "does not have the highest ranking, and no primary lesion and initial disease state sample was found.\n"
          )
          sample_details <- paste0(
            capture.output(print(patient_samples[, c("sample", "age_at_enrollment", "lesion_type", "disease_state")])),
            collapse = "\n"
          )
          warnings[[length(warnings) + 1]] <- paste0(warning_message, sample_details)
        }
      }
    }
  } else {
    # Version 2: Summarize at the sample level for selection purposes
    df_sampled <- df %>%
      group_by(sample, base_sample, age_at_enrollment, lesion_type, disease_state) %>%
      summarise(
        lesion_rank = dplyr::first(lesion_rank),
        disease_rank = dplyr::first(disease_rank),
        .groups = "drop"
      )
    
    df_selected_summary <- df_sampled %>%
      group_by(base_sample) %>%
      arrange(age_at_enrollment, lesion_rank, disease_rank) %>%
      slice(1) %>%
      ungroup()
    
    df_selected <- df %>%
      filter(sample %in% df_selected_summary$sample)
    
    # Initialize warnings
    warnings <- list()
    
    unique_base_samples <- unique(df_sampled$base_sample)
    for (base in unique_base_samples) {
      patient_samples <- df_sampled[df_sampled$base_sample == base, ]
      selected_sample <- df_selected_summary[df_selected_summary$base_sample == base, ]
      
      # Check rankings
      max_lesion_rank <- min(patient_samples$lesion_rank)
      max_disease_rank <- min(patient_samples$disease_rank)
      
      # Debugging: Print problematic rows
     # cat("Checking base_sample:", base, "\n")
     # print(selected_sample)
     # print(patient_samples)
     # cat("Max lesion rank:", max(patient_samples$lesion_rank, na.rm = TRUE), "\n")
     # cat("Max disease rank:", max(patient_samples$disease_rank, na.rm = TRUE), "\n")
      
      if (selected_sample$lesion_rank > max_lesion_rank || selected_sample$disease_rank > max_disease_rank) {
        primary_initial_sample <- patient_samples %>%
          filter(lesion_type == "primary", disease_state == "initial") %>%
          slice(1)
        
        if (nrow(primary_initial_sample) > 0) {
          df_selected_summary <- df_selected_summary %>%
            filter(!(base_sample == base)) %>%
            bind_rows(primary_initial_sample)
          
          df_selected <- df %>%
            filter(sample %in% df_selected_summary$sample)
          
          warning_message <- paste(
            "Warning: The selected sample for patient", base,
            "does not have the highest ranking. A primary lesion and initial disease state sample was selected instead.\n"
          )
          warnings[[length(warnings) + 1]] <- warning_message
        } else {
          warning_message <- paste(
            "Warning: The selected sample for patient", base,
            "does not have the highest ranking, and no primary lesion and initial disease state sample was found.\n"
          )
          warnings[[length(warnings) + 1]] <- warning_message
        }
      }
    }
  }
  
  # Print warnings
  if (length(warnings) > 0) {
    cat(paste(warnings, collapse = "\n\n"))
  }
  
  # Return the full dataframe and the selected samples
  return(list(all_samples = as.data.frame(df), selected_samples = as.data.frame(df_selected)))
}

nsample <- function(df, group) {
  df %>%
    group_by(!!sym(group)) %>%
    summarise(unique_samples = n_distinct(sample))
}

compare_proportions <- function(df, group_column, total_column) {
  # Step 1: Create a binary column for the event
  df <- df %>%
    mutate(event = ifelse(.data[[total_column]] > 0, "Yes", "No"))
  
  # Step 2: Create a contingency table for the event and group column
  contingency_table <- table(df$event, df[[group_column]])
  
  # Step 3: Print the contingency table
  cat("Contingency Table:\n")
  print(contingency_table)
  
  # Step 4: Perform Chi-Square Test
  chi_test_result <- chisq.test(contingency_table)
  
  # Step 5: Return results
  cat("\nChi-Square Test Results:\n")
  print(chi_test_result)
  
  # Return the test result for further use
  return(chi_test_result)
}

compare_proportions_fisher <- function(df, group_column, total_column) {
  # Step 1: Create a binary column for the event
  df <- df %>%
    mutate(event = ifelse(.data[[total_column]] > 0, "Yes", "No"))
  
  # Step 2: Create a contingency table for the event and group column
  contingency_table <- table(df$event, df[[group_column]])
  
  # Step 3: Print the contingency table
  cat("Contingency Table:\n")
  print(contingency_table)
  
  # Step 4: Perform Fisher's Exact Test
  fisher_test_result <- fisher.test(contingency_table)
  
  # Step 5: Print and return results
  cat("\nFisher's Exact Test Results:\n")
  print(fisher_test_result)
  
  # Return the test result for further use
  return(fisher_test_result)
}

# Helper function to construct the combination label
construct_combination_label <- function(chr, type) {
  if (is.na(chr) && is.na(type)) {
    return("total")
  } else if (is.na(chr)) {
    return(type)
  } else if (is.na(type)) {
    return(paste0("chr", chr))
  } else {
    return(paste0("chr", chr, "_", type))
  }
}

# plot box plot
plot_box <- function(df, count_y, group, xlab, y_lab, colours, log_scale = FALSE) {
  p <- ggplot(df, aes(x = !!sym(group), y = !!sym(count_y), fill = !!sym(group))) +
    geom_boxplot(outlier.shape = NA) +
    scale_fill_manual(values = colours) +
    geom_jitter(position = position_jitterdodge(jitter.width = 0.2), color = "black", size = 1.5) +
    labs(x = xlab, y = y_lab) + 
    guides(fill = "none")
  
  # Apply logarithmic scale if log_scale is TRUE
  if (log_scale) {
    p <- p + scale_y_continuous(
      trans = scales::log1p_trans()
    )
  }
  
  
  return(p)
}

# plot total te calls
plot_te_sum <- function(df) {
  # Step 1: Summing each TE column across all samples
  sum_data <- df %>%
    summarise(
      ALU = sum(ALU, na.rm = TRUE),
      LINE = sum(LINE1, na.rm = TRUE),
      SVA = sum(SVA, na.rm = TRUE)
    )
  
  # Step 2: Reshape the data for plotting
  plot_data <- sum_data %>%
    pivot_longer(cols = everything(), names_to = "TE_type", values_to = "Total")
  
  # Step 3: Create the bar plot
  plot <- ggplot(plot_data, aes(x = TE_type, y = Total, fill = TE_type)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = colour_palette_3) +
    labs(
      x = "TE Type",
      y = "Total TE count of all samples"
    ) 
  
  return(plot)
}

# plot number of tes total per group
# wilcoxon test
plot_count_wilcox <- function(df, chr, type, group, x_lab, y_lab, log_scale = FALSE) {
  # Use the helper function to construct the combination label
  combination <- construct_combination_label(chr, type)
  
  # Check if there are any NA values in the group column
  if (any(is.na(df[[group]]))) {
    cat("Group contains NA values. Excluding NA group from the plot and Wilcoxon test.\n")
    df <- df %>% filter(!is.na(!!sym(group)))
  }
  
  # Check if grouping factor has exactly 2 levels for Wilcoxon test
  unique_groups <- unique(df[[group]])
  if (length(unique_groups) != 2) {
    cat(paste("Warning: Grouping factor has", length(unique_groups), "levels:", paste(unique_groups, collapse=", "), "\n"))
    cat("Wilcoxon test requires exactly 2 groups. Skipping statistical test.\n")
    p_value <- NA
    p_value_formatted <- "N/A"
  } else {
    # Print p-value from Wilcoxon test
    print("dependent ~ independent")
    formula <- reformulate(group, combination)
    print(formula)
    p_value <- wilcox.test(formula, data = df)$p.value
    print(paste("Wilcoxon test p-value:", p_value))
    # Format the p-value for display
    p_value_formatted <- formatC(p_value, format = "e", digits = 2)
  }
  
  # Calculate and print medians for each group
  medians <- df %>%
    group_by(!!sym(group)) %>%
    summarise(median_count = median(!!sym(combination), na.rm = TRUE))
  print(medians)
  
  # Choose appropriate color palette based on number of groups
  n_groups <- length(unique_groups)
  if (n_groups <= 2) {
    color_palette <- colours
  } else if (n_groups == 3) {
    color_palette <- colours_3
  } else {
    # For more than 3 groups, use a default palette
    color_palette <- rainbow(n_groups)
  }
  
  # Plot with optional log scale
  plot <- plot_box(df, combination, group, x_lab, y_lab, color_palette, log_scale)
  
  # Add statistical annotations only if we have exactly 2 groups
  if (length(unique_groups) == 2) {
    plot <- plot + 
      geom_signif(test = "wilcox.test", 
                  comparisons = list(levels(factor(df[[group]]))),
                  map_signif_level = TRUE,
                  textsize = 5) + 
      coord_cartesian(clip = 'off') # dont cut off annotation
  }
  
  return(plot)
}

generate_plots <- function(plot_function, df, types, y_end = "count", ...) {
  # Use lapply to iterate over each type and print the corresponding plot
  invisible(lapply(types, function(type) {
    # Set the y-axis label based on type and y_end
    y_lab <- if (is.na(type)) {
      paste("TE", y_end)
    } else {
      paste(type, y_end)
    }
    
    # Call the provided plotting function and print the plot
    print(plot_function(df, type = type, y_lab = y_lab, ...))
  }))
}

stacked_bar_plot_num_te <- function(data, min_samples = 3, te_thresholds = c(0, 10, 100), te_type = NULL) {
  # Step 1: Select a specific TE type column if provided
  if (!is.null(te_type) && te_type %in% colnames(data)) {
    data <- data %>%
      mutate(total = !!sym(te_type))  # Set `total` to the selected TE type column
  }
  
  # Step 2: Define TE ranges based on thresholds and assign categories
  data <- data %>%
    mutate(
      TE_category = case_when(
        total == te_thresholds[1] ~ "0",
        total > te_thresholds[1] & total <= te_thresholds[2] ~ "1-10",
        total > te_thresholds[2] & total <= te_thresholds[3] ~ "10-100",
        total > te_thresholds[3] ~ ">100"
      ),
      TE_category = factor(TE_category, levels = c(">100", "10-100", "1-10", "0")) # Set factor levels for order
    )
  
  # Step 3: Filter out tumor types with fewer than the minimum sample count
  data <- data %>%
    group_by(tumor_type) %>%
    filter(n() >= min_samples) %>%
    ungroup()
  
  # Step 4: Calculate the proportion of samples within each TE range for each tumor type
  data_summary <- data %>%
    group_by(tumor_type, TE_category) %>%
    summarise(count = n(), .groups = 'drop') %>%
    group_by(tumor_type) %>%
    mutate(proportion = count / sum(count)) %>%
    ungroup()
  
  # Step 5: Sort tumor types by the proportion of >100 (red), then 10-100 (orange), then 1-10 (yellow)
  tumor_type_order <- data_summary %>%
    group_by(tumor_type) %>%
    summarise(
      prop_red = sum(proportion[TE_category == ">100"]),
      prop_orange = sum(proportion[TE_category == "10-100"]),
      prop_yellow = sum(proportion[TE_category == "1-10"]),
      .groups = 'drop'
    ) %>%
    arrange(desc(prop_red), desc(prop_orange), desc(prop_yellow)) %>%
    pull(tumor_type)
  
  # Step 6: Add sample count to tumor type labels
  tumor_type_labels <- data %>%
    group_by(tumor_type) %>%
    summarise(num_samples = n(), .groups = 'drop') %>%
    mutate(label = paste0(tumor_type, " (", num_samples, ")")) %>%
    dplyr::select(tumor_type, label) %>%
    deframe()  # Convert to named vector for ggplot2
  
  # Set tumor_type as a factor with the sorted levels
  data_summary <- data_summary %>%
    mutate(tumor_type = factor(tumor_type, levels = tumor_type_order))
  
  # Step 7: Plot the data as a stacked bar chart
  ggplot(data_summary, aes(x = tumor_type, y = proportion, fill = TE_category)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(
      x = "Tumor Type",
      y = "Proportion of Samples",
      fill = "TE Range"
    ) +
    scale_x_discrete(labels = tumor_type_labels) +  # Apply labels with sample counts
    scale_fill_manual(
      values = c("0" = "grey", "1-10" = "yellow", "10-100" = "orange", ">100" = "red")
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

stacked_bar_plot_per_chromosome <- function(data, te_thresholds = c(0, 10, 100), type = NA) {
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  
  # Step 1: Determine the type of columns to select
  if (is.na(type)) {
    chr_columns <- grep("^chr[0-9XY]+$", colnames(data), value = TRUE)
  } else if (type == "LINE1") {
    chr_columns <- grep("^chr[0-9XY]+_LINE1$", colnames(data), value = TRUE)
  } else if (type == "ALU") {
    chr_columns <- grep("^chr[0-9XY]+_ALU$", colnames(data), value = TRUE)
  } else if (type == "SVA") {
    chr_columns <- grep("^chr[0-9XY]+_SVA$", colnames(data), value = TRUE)
  } else {
    stop("Invalid type. Use NA, 'LINE1', 'ALU', or 'SVA'.")
  }
  
  if (length(chr_columns) == 0) {
    stop("No matching chromosome-specific columns found in the dataframe.")
  }
  
  # Step 2: Pivot data to long format for easier processing
  data_long <- data %>%
    select(sample, all_of(chr_columns)) %>%
    pivot_longer(cols = -sample, names_to = "Chromosome", values_to = "Total") %>%
    mutate(Chromosome = gsub("_.*$", "", Chromosome))  # Remove suffix (e.g., _LINE1)
  
  # Step 3: Categorize chromosome counts based on thresholds
  categories <- c(
    "0",
    paste0(te_thresholds[1] + 1, "-", te_thresholds[2]),
    paste0(te_thresholds[2] + 1, "-", te_thresholds[3]),
    paste0(">", te_thresholds[length(te_thresholds)])
  )
  
  data_long <- data_long %>%
    mutate(
      TE_category = case_when(
        Total == te_thresholds[1] ~ "0",
        Total > te_thresholds[1] & Total <= te_thresholds[2] ~ paste0(te_thresholds[1] + 1, "-", te_thresholds[2]),
        Total > te_thresholds[2] & Total <= te_thresholds[3] ~ paste0(te_thresholds[2] + 1, "-", te_thresholds[3]),
        Total > te_thresholds[3] ~ paste0(">", te_thresholds[3])
      ),
      TE_category = factor(TE_category, levels = rev(categories))  # Reverse levels to put "0" closest to the X-axis
    )
  
  # Step 4: Summarize data for the plot
  data_summary <- data_long %>%
    group_by(Chromosome, TE_category) %>%
    summarize(Count = n(), .groups = "drop") %>%
    group_by(Chromosome) %>%
    mutate(Proportion = Count / sum(Count)) %>%
    ungroup()
  
  # order chromosomes
  data_summary <- data_summary %>%
    mutate(
      Chromosome = factor(Chromosome, 
                          levels = c(paste0("chr", 1:22), "chrX", "chrY"))  # Correct numeric order
    )
  
  # Step 5: Dynamically generate the color palette
  dynamic_palette <- setNames(
    c("grey", "yellow", "orange", "red"),  # Fixed order to ensure "0" is always gray
    categories
  )
  
  # Step 6: Create the stacked bar plot
  plot <- ggplot(data_summary, aes(x = Chromosome, y = Proportion, fill = TE_category)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_y_continuous(expand = c(0, 0)) + 
    labs(
      x = "Chromosome",
      y = "Proportion of Samples",
      fill = "TE Count"
    ) +
    scale_fill_manual(values = dynamic_palette) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.5),
          axis.title.x = element_text(margin = margin(t = 22)))
  
  return(plot)
}

stats_and_plot <- function(data, group = "TP53_status", count_col = "total") {
  # Perform Wilcoxon test
  p_value <- wilcox.test(reformulate(group, count_col), data = data)$p.value
  print(paste("Wilcoxon test p-value:", p_value))
  
  # Calculate and print medians for each group
  medians <- data %>%
    group_by(!!sym(group)) %>%
    summarise(median_count = median(!!sym(count_col), na.rm = TRUE))
  print(medians)
  
  # Use your custom plot_box function for plotting
  plot <- plot_box(data, count_y = count_col, group = group, xlab = "TP53 Status", y_lab = "Count", 
                   colours = colours) +
    geom_signif(comparisons = list(c("Control", "LFS")), test = "wilcox.test", 
                map_signif_level = TRUE, textsize = 5) +
    coord_cartesian(clip = 'off') # To prevent clipping of annotations
  
  # Return both the plot and statistics as a list
  list(plot = plot, p_value = p_value, medians = medians)
}

plot_count_lm_nomerge_tt <- function(df, chr, type, group, covariates, x_lab, y_lab, log_scale = FALSE) {
  # Use the helper function to construct the combination label
  combination <- construct_combination_label(chr, type)
  
  # Check if there are any NA values in the group or covariates columns
  all_vars <- c(group, covariates)
  if (anyNA(df[all_vars])) {
    cat("Data contains NA values in group or covariates. Excluding rows with NA values.\n")
    df <- df %>% filter(complete.cases(df[all_vars]))
  }
  
  # Construct the formula for the linear model
  formula <- as.formula(
    paste(combination, "~", group, "+", paste(covariates, collapse = "+"))
  )
  
  # Fit the linear model
  lm_model <- lm(formula, data = df)
  lm_summary <- summary(lm_model)
  
  # Print the model summary
  print(lm_summary)
  
  # Extract the p-value for the group effect
  dummy_var <- grep(paste0("^", group), rownames(coef(lm_summary)), value = TRUE)
  group_p_value <- coef(lm_summary)[dummy_var, "Pr(>|t|)"]
  print(paste("Linear model p-value for", group, ":", group_p_value))
  
  # Format the p-value for display
  p_value_formatted <- formatC(group_p_value, format = "e", digits = 2)
  
  # Calculate and print medians for each group
  medians <- df %>%
    group_by(!!sym(group)) %>%
    summarise(median_count = median(!!sym(combination), na.rm = TRUE))
  print(medians)
  
  # Plot with optional log scale
  plot <- plot_box(df, combination, group, x_lab, y_lab, colours, log_scale) + 
    annotate("text", x = 1.5, y = max(df[[combination]], na.rm = TRUE), 
             label = paste("p =", p_value_formatted), size = 5, hjust = 0) +
    coord_cartesian(clip = 'off') # Avoid cutting off annotations
  
  return(plot)
}

plot_count_lasso <- function(df, chr, type, group, covariates, x_lab, y_lab, log_scale = FALSE, alpha = 1, lambda = NULL) {
  # Use the helper function to construct the combination label
  combination <- construct_combination_label(chr, type)
  
  # Check if there are any NA values in the group or covariates columns
  all_vars <- c(group, covariates)
  if (anyNA(df[all_vars])) {
    cat("Data contains NA values in group or covariates. Excluding rows with NA values.\n")
    df <- df %>% filter(complete.cases(df[all_vars]))
  }
  
  # Ensure group is a factor
  df[[group]] <- factor(df[[group]])
  
  # Create the model matrix
  X <- model.matrix(as.formula(paste("~", group, "+", paste(covariates, collapse = "+"))), data = df)[, -1]
  y <- df[[combination]]
  
  # Fit Lasso regression
  lasso_model <- glmnet(X, y, alpha = alpha, lambda = lambda)
  
  # Cross-validation to select lambda if not provided
  if (is.null(lambda)) {
    cv_lasso <- cv.glmnet(X, y, alpha = alpha)
    lambda <- cv_lasso$lambda.min
    print(paste("Selected Lambda:", lambda))
    lasso_model <- glmnet(X, y, alpha = alpha, lambda = lambda)
  }
  
  # Extract coefficients
  coef_lasso <- coef(lasso_model, s = lambda)
  print("Lasso Coefficients:")
  print(coef_lasso)
  
  # Plot with optional log scale
  plot <- plot_box(df, combination, group, x_lab, y_lab, colours, log_scale) +
    annotate("text", x = 1.5, y = max(df[[combination]], na.rm = TRUE),
             label = paste("Lasso Regularization Applied"), size = 5, hjust = 0) +
    coord_cartesian(clip = 'off') # Avoid cutting off annotations
  
  return(plot)
}

plot_count_lm <- function(df, chr, type, group, covariates, x_lab, y_lab, residuals = FALSE, log_scale = FALSE, breaks=NULL, min_samples = 5, fill_palette = NULL) {
  # Use the helper function to construct the combination label
  combination <- construct_combination_label(chr, type)
  print(table(df[[group]]))
  # Check if there are any NA values in the group or covariates columns
  all_vars <- c(group, covariates)
  if (anyNA(df[all_vars])) {
    cat("Data contains NA values in group or covariates. Excluding rows with NA values.\n")
    df <- df %>% filter(complete.cases(df[all_vars]))
  }
  
  # Aggregate rare tumor types
  df <- df %>%
    group_by(tumor_type) %>%
    mutate(tumor_type = ifelse(n() < min_samples, "Other", tumor_type)) %>%
    ungroup()
  
  # Ensure group is a factor
  df[[group]] <- factor(df[[group]])
  
  # Construct the formula for the linear model
  formula <- as.formula(
    paste(combination, "~", if (residuals) paste(covariates, collapse = "+") else paste(group, "+", paste(covariates, collapse = "+")))
  )
  print(table(df[[group]]))
  # Fit the linear model
  lm_model <- lm(formula, data = df)
  lm_summary <- summary(lm_model)
  
  # Print the model summary
  print(lm_summary)
  print(formula)
  
  if (residuals) {
    # Approach 1: Extract residuals and perform Wilcoxon test
    df$residuals <- residuals(lm_model)
    
    # Perform Wilcoxon test on residuals between groups
    wilcox_result <- wilcox.test(reformulate(group, response = "residuals"), data = df, exact = FALSE)
    wilcox_p_value <- wilcox_result$p.value
    cat("p-value for residuals between groups:", wilcox_p_value, "\n")
    
    # Format the p-value for display
    p_value_formatted <- formatC(wilcox_p_value, format = "e", digits = 2)
    
    # Calculate and print medians for residuals
    medians <- df %>%
      group_by(!!sym(group)) %>%
      summarise(median_residuals = median(residuals, na.rm = TRUE))
    print(medians)
    
    # Plot residuals
    p <- ggplot(df, aes(x = !!sym(group), y = residuals)) +
      geom_boxplot(outlier.shape = NA, color = "black", fill = "#5FBFF9") +
      geom_jitter(color = "black", size = 1.5, width = 0.2) +
      labs(x = x_lab, y = "Residuals") +
      annotate("text", x = 1.5, y = max(df$residuals, na.rm = TRUE),
               label = paste("p =", p_value_formatted), size = 5, hjust = 0) 
    
    if (log_scale) {
      p <- p + scale_y_continuous(
        trans = scales::log1p_trans(),
        breaks = if (!is.null(breaks)) breaks else waiver()
      )
    }
    
    return(p)
    #return(list(lm_model = lm_model, wilcox_result = wilcox_result, plot = plot))
    
  } else {
    # Approach 2: Analyze group effect in the linear model
    # Print the full linear model summary to see pairwise comparisons vs reference group
    cat("\n=== LINEAR MODEL SUMMARY ===\n")
    cat("Coefficients (each group compared to reference group):\n")
    print(coef(lm_summary))
    cat("\n")

    # Print ANOVA table to see overall F-test for each variable
    cat("=== ANOVA TABLE (F-tests for overall effect of each variable) ===\n")
    anova_result <- anova(lm_model)
    print(anova_result)
    cat("\n")

    # Use ANOVA F-test to get overall group effect p-value for plot
    # This tests: "Do ANY of the groups differ?" (works for both 2-level and multi-level groups)
    group_p_value <- anova_result[group, "Pr(>F)"]
    cat("Using overall ANOVA F-test p-value for", group, ":", group_p_value, "\n")
    cat("(This tests whether ANY groups differ, adjusting for all covariates)\n\n")

    # Format the p-value for display
    p_value_formatted <- formatC(group_p_value, format = "e", digits = 2)
    
    # Calculate and print medians for the group
    medians <- df %>%
      group_by(.data[[group]]) %>%
      summarise(median_count = median(!!sym(combination), na.rm = TRUE))
    print(medians)
    
    # Plot the group effect
    # Use provided palette or default to colour_palette_4
    if (is.null(fill_palette)) {
      fill_palette <- colour_palette_4
    }

    p <- ggplot(df, aes(x = .data[[group]], y = !!sym(combination), fill = .data[[group]])) +
      geom_boxplot(outlier.shape = NA, color = "black") +
      scale_fill_manual(values = fill_palette, guide = "none") +
      geom_jitter(color = "black", size = 1.5, width = 0.2) +
      labs(x = x_lab, y = y_lab) +
      annotate("text", x = 1.7, y = max(df[[combination]], na.rm = TRUE),
               label = paste("p =", p_value_formatted), hjust = 0) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
    
    if (log_scale) {
      p <- p + scale_y_continuous(
        trans = scales::log1p_trans(),
        breaks = if (!is.null(breaks)) breaks else waiver()
      )
    }
    
    return(p)
    #return(list(lm_model = lm_model, plot = plot, group_p_value = group_p_value))
  }
}

bootstrap_test <- function(df, column, value_col, test_function, n_bootstraps = 1000, step_size = 100, ...) {
  # `step_size` determines intervals at which to calculate p-values
  # Remove rows with NA values in the specified columns (both grouping and value columns)
  df_clean <- df[!is.na(df[[column]]) & !is.na(df[[value_col]]), ]
  
  # Get the two unique groups in the column (after removing NA)
  groups <- unique(df_clean[[column]])
  
  # Check if there are exactly two unique groups
  if (length(groups) != 2) {
    stop("The column must contain exactly two unique groups, excluding NAs.")
  }
  
  # Subset the groups after cleaning the data
  group1 <- df_clean[[value_col]][df_clean[[column]] == groups[1]]
  group2 <- df_clean[[value_col]][df_clean[[column]] == groups[2]]
  
  # Calculate the observed test statistic
  if (identical(test_function, wilcox.test)) {
    observed_stat <- test_function(group1, group2, exact = FALSE, ...)$statistic
  } else {
    observed_stat <- test_function(group1) - test_function(group2)
  }
  
  # Combine both groups
  combined <- c(group1, group2)
  
  # Initialize vectors for bootstrap test statistics and p-value tracking
  bootstrap_stats <- numeric(n_bootstraps)
  p_value_convergence <- numeric(floor(n_bootstraps / step_size))
  bootstrap_sizes <- seq(step_size, n_bootstraps, by = step_size)
  
  for (i in 1:n_bootstraps) {
    # Resample from the combined population
    resampled_group1 <- sample(combined, length(group1), replace = TRUE)
    resampled_group2 <- sample(combined, length(group2), replace = TRUE)
    
    # Calculate the test statistic for this bootstrap sample
    if (identical(test_function, wilcox.test)) {
      bootstrap_stats[i] <- test_function(resampled_group1, resampled_group2, exact = FALSE, ...)$statistic
    } else {
      bootstrap_stats[i] <- test_function(resampled_group1) - test_function(resampled_group2)
    }
    
    # Track p-value at intervals of step_size
    if (i %% step_size == 0) {
      current_p_value <- mean(abs(bootstrap_stats[1:i]) >= abs(observed_stat))
      p_value_convergence[i / step_size] <- current_p_value
    }
  }
  
  # Calculate the final p-value (two-sided)
  p_value <- mean(abs(bootstrap_stats) >= abs(observed_stat))
  
  # Print the final p-value
  cat("Final P-value:", p_value, "\n")
  
  # Plot the bootstrap distribution with ggplot2
  plot_data <- data.frame(bootstrap_stats = bootstrap_stats)
  p_bootstrap <- ggplot(plot_data, aes(x = bootstrap_stats)) +
    geom_histogram(bins = 30, fill = "#5FBFF9", color = "black", alpha = 0.7) +
    geom_vline(xintercept = observed_stat, color = "red", linewidth = 1, linetype = "solid") +
    labs(
         x = "Test Statistic",
         y = "Frequency")
  
  print(p_bootstrap)
  
  # Plot the p-value convergence
  convergence_data <- data.frame(bootstrap_sizes = bootstrap_sizes, p_values = p_value_convergence)
  p_convergence <- ggplot(convergence_data, aes(x = bootstrap_sizes, y = p_values)) +
    geom_line(color = "blue", linewidth = 1) +
    geom_point(color = "blue") +
    labs(
         x = "Number of Bootstraps",
         y = "P-value") 
  
  print(p_convergence)
  
  # Return the results
  list(
    observed_stat = observed_stat,
    bootstrap_stats = bootstrap_stats,
    p_value = p_value,
    p_value_convergence = p_value_convergence
  )
}

# plot number of tes total per group
# wilcoxon test
plot_count_location_wilcox<- function(sv_df, filter_var="SV.type", filter_element, group, chr=NA, location=NA, gene=NA, x_lab, y_lab){
  if (!is.na(gene)){
    sv_df <- sv_df %>% filter(Gene_name==gene)
  }
  if (!is.na(chr)){
    sv_df <- sv_df %>% filter(SV.chrom==chr)
  }
  if (is.character(location) && !any(is.na(location))){
    if (length(location) == 1){
      sv_df <- sv_df %>% filter(Location2 == location)
    } else if (length(location) == 2){
      sv_df <- sv_df %>% filter(Location2 == location[1] | Location2 == location[2])
    }
  }
  
  if (!is.na(filter_element)){ # if need to filter by SV type
    data <- sv_df %>%
      filter(!!sym(filter_var) == {{ filter_element}} ) %>% # filter for deletions
      group_by(!!sym(group), sample) %>% # group
      summarise(count = n(), .groups = 'drop') # count
  }
  else { # if doing overall TE count
    data <- sv_df %>%
      group_by(!!sym(group), sample) %>% # group
      summarise(count = n(), .groups = 'drop') # count by group
  }
  
  # Calculate medians
  medians <- data %>%
    group_by(!!sym(group)) %>%
    summarise(median_count = median(count))
  
  # print p value
  print(wilcox.test(reformulate(group, "count"), data = data)$p.value)
  print(medians)
  
  # Plot
  plot <- plot_box(data, count, !!sym(group), x_lab, y_lab) + 
    geom_signif(test = "wilcox.test", # sig test
                comparisons = list(levels(factor(sv_df[[group]]))),
                map_signif_level = TRUE,
                textsize = 5) + # stars 
    coord_cartesian(clip = 'off') # dont cut off annotation
  return(plot)
}

plot_count_mult_location_wilcox <- function(sv_df, filter_var="SV.type", filter_element, group, chr=NA, location=NA, gene=NA, x_lab, y_lab){
  if (!is.na(gene)){
    sv_df <- sv_df %>% filter(Gene_name == gene)
  }
  if (!is.na(chr)){
    sv_df <- sv_df %>% filter(SV.chrom == chr)
  }
  if (is.character(location) && !any(is.na(location))){
    if (length(location) == 1){
      sv_df <- sv_df %>% filter(Location2 == location)
    } else if (length(location) > 1){
      sv_df <- sv_df %>% filter(Location2 %in% location)
    }
  }
  
  if (!is.na(filter_element)){ # if need to filter by SV type
    data <- sv_df %>%
      filter(!!sym(filter_var) == {{ filter_element }} ) %>%
      group_by(!!sym(group), sample, Location2) %>%
      summarise(count = n(), .groups = 'drop')
  } else { # if doing overall TE count
    data <- sv_df %>%
      group_by(!!sym(group), sample, Location2) %>%
      summarise(count = n(), .groups = 'drop')
  }
  
  # Calculate medians
  medians <- data %>%
    group_by(Location2, !!sym(group)) %>%
    summarise(median_count = median(count), .groups = 'drop')
  
  # Print p-value and medians
  p_value <- wilcox.test(reformulate(group, "count"), data = data)$p.value
  print(paste("P-value:", p_value))
  print(medians)
  
  # Calculate p-values for each Location2
  p_values <- sapply(split(data, data$Location2), function(x) {
    wilcox.test(count ~ get(group), data = x)$p.value
  })
  
  # Create labels
  labels <- symnum(p_values, corr = FALSE, cutpoints = c(0, .001, .01, .05, 1), symbols = c("***", "**", "*", "n.s."))
  
  # Calculate y values for significance annotations
  y_values <- sapply(split(data, data$Location2), function(x) {
    max(sapply(split(x, x[[group]]), function(xx) {
      boxplot(xx$count, plot = FALSE)$stats[5, ]
    })) + 2
  })
  
  # Create position vectors
  data$interaction_var <- interaction(data$Location2, data[[group]], drop = TRUE)
  levels_interaction_var <- levels(data$interaction_var)
  positions <- data.frame(
    Location2 = rep(unique(data$Location2), each = 2),
    xmin = rep(seq(1, length(unique(data$Location2))) - 0.2, each = 2),
    xmax = rep(seq(1, length(unique(data$Location2))) + 0.2, each = 2),
    y = rep(y_values, each = 2)
  )
  
  # Plot
  plot <- ggplot(data, aes(x = interaction(Location2, !!sym(group)), y = count, fill = !!sym(group))) + 
    geom_boxplot(position = position_dodge(width = 0.75)) +
    labs(x = x_lab, y = y_lab) +
    theme_minimal() +
    geom_signif(
      y_position = positions$y,
      xmin = positions$xmin,
      xmax = positions$xmax,
      annotations = rep(labels, each = 2),
      textsize = 5
    ) +
    coord_cartesian(clip = 'off') +
    scale_x_discrete(labels = function(x) gsub("\\..*", "", x)) +
    theme(axis.text.x = element_text(angle = 45, vjust=5))
  
  return(plot)
}

plot_count_kruskal <- function(df, chr=NA, type=NA, group, x_lab, y_lab, x_order = NULL, log_scale = FALSE, breaks=NULL) {
  # Use the helper function to construct the combination label
  combination <- construct_combination_label(chr, type)
  
  # Check if there are any NA values in the group column
  if (any(is.na(df[[group]]))) {
    cat("Group contains NA values. Excluding NA group from the plot and Wilcoxon test.\n")
    df <- df %>% filter(!is.na(!!sym(group)))
  }
  
  # Check if there are at least 2 groups for Kruskal-Wallis test
  unique_groups <- unique(df[[group]])
  unique_groups <- unique_groups[!is.na(unique_groups)]
  
  if (length(unique_groups) < 2) {
    cat("Kruskal-Wallis test requires at least 2 groups. Column:", group, "has", length(unique_groups), "unique value(s):", paste(unique_groups, collapse = ", "), ". Skipping statistical test.\n")
    p_value <- NA
    p_value_formatted <- "NA"
  } else {
    # Perform Kruskal-Wallis test
    test_result <- kruskal.test(reformulate(group, combination), data = df)
    p_value <- test_result$p.value
    p_value_formatted <- formatC(p_value, format = "e", digits = 2)
    print(p_value)
  }

  # Calculate and print medians for each group
  medians <- df %>%
    group_by(!!sym(group)) %>%
    summarise(median_count = median(!!sym(combination), na.rm = TRUE))
  print(medians)
  
  # Plot
  plot <- plot_box_kruskal(df, combination, group = group, xlab = x_lab, y_lab = y_lab, x_order = x_order, log_scale) 
  plot <- plot +  
    annotate("text", x = Inf, y = Inf, label = paste("p =", p_value_formatted), vjust = 2, hjust = 1, size = 4)
  
  # Apply logarithmic scale if log_scale is TRUE
  if (log_scale) {
    plot <- plot + scale_y_continuous(
      trans = scales::log1p_trans(),
      breaks = if (!is.null(breaks)) breaks else waiver()
    )
  }
  
  return(plot)
}

plot_box_kruskal<- function(df, count_y, group, xlab, y_lab, x_order = NULL, log_scale = FALSE){
  # Filter out non-finite values
  df_clean <- df %>% filter(is.finite(!!sym(count_y)))
  
  plot <- ggplot(df_clean, aes(x = !!sym(group), y = !!sym(count_y), fill = !!sym(group))) +
    geom_boxplot(outlier.shape = NA) + 
    scale_fill_manual(values =  color_palette_6) + #mutation_colours
    geom_jitter(position=position_jitterdodge(jitter.width=0.2), color = "black", size = 1.5) +
    labs(x = xlab, y = y_lab) + 
    guides(fill="none") + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  
  # Apply custom order for x-axis if provided
  if (!is.null(x_order)) {
    plot <- plot + scale_x_discrete(limits = x_order)
  }
  
  # Apply logarithmic scale if log_scale is TRUE
  if (log_scale) {
    plot <- plot + scale_y_continuous(
      trans = scales::log1p_trans()
    )
  }
  return(plot)
}

# Faceted version of plot_count_kruskal by ancestry
# Creates separate panels for each ancestry group with individual Kruskal-Wallis tests
plot_count_kruskal_facet_ancestry <- function(df, chr=NA, type=NA, group, x_lab, y_lab,
                                              x_order = NULL, log_scale = FALSE, breaks=NULL,
                                              ancestry_col = "predicted_ancestry_thres") {
  # Use the helper function to construct the combination label
  combination <- construct_combination_label(chr, type)

  # Check if ancestry column exists
  if (!ancestry_col %in% colnames(df)) {
    cat("Warning: Ancestry column '", ancestry_col, "' not found in dataframe. Cannot create faceted plot.\n", sep="")
    return(NULL)
  }

  # Remove rows with NA in group or ancestry columns
  df_clean <- df %>%
    filter(!is.na(!!sym(group)), !is.na(!!sym(ancestry_col)))

  if (nrow(df_clean) == 0) {
    cat("Warning: No data remaining after removing NA values. Cannot create faceted plot.\n")
    return(NULL)
  }

  # Calculate Kruskal-Wallis test p-values for each ancestry group
  ancestry_groups <- unique(df_clean[[ancestry_col]])
  p_values_df <- data.frame()

  for (anc in ancestry_groups) {
    df_anc <- df_clean %>% filter(!!sym(ancestry_col) == anc)
    unique_groups <- unique(df_anc[[group]])
    unique_groups <- unique_groups[!is.na(unique_groups)]

    if (length(unique_groups) >= 2) {
      test_result <- kruskal.test(reformulate(group, combination), data = df_anc)
      p_value <- test_result$p.value
      p_value_formatted <- formatC(p_value, format = "e", digits = 2)
    } else {
      p_value_formatted <- "NA"
    }

    p_values_df <- rbind(p_values_df, data.frame(
      ancestry = anc,
      p_value = p_value_formatted,
      stringsAsFactors = FALSE
    ))
  }

  # Print p-values for each ancestry group
  cat("\nKruskal-Wallis p-values by ancestry:\n")
  print(p_values_df)

  # Rename ancestry column in p_values_df to match the actual column name
  # This ensures geom_text properly matches p-values to facets
  colnames(p_values_df)[1] <- ancestry_col

  # Calculate and print medians for each group x ancestry combination
  medians <- df_clean %>%
    group_by(!!sym(ancestry_col), !!sym(group)) %>%
    summarise(median_count = median(!!sym(combination), na.rm = TRUE), .groups = "drop")
  cat("\nMedians by ancestry and group:\n")
  print(medians)

  # Create base plot with faceting
  plot <- plot_box_kruskal_facet_ancestry(df_clean, combination, group = group,
                                          ancestry_col = ancestry_col,
                                          xlab = x_lab, y_lab = y_lab,
                                          x_order = x_order, log_scale = log_scale)

  # Add p-values as text annotations for each facet (only for non-NA values)
  # Position at top right of each panel
  p_values_df_valid <- p_values_df[p_values_df$p_value != "NA", ]
  if (nrow(p_values_df_valid) > 0) {
    plot <- plot +
      geom_text(data = p_values_df_valid,
                aes(x = Inf, y = Inf, label = paste("p =", p_value)),
                inherit.aes = FALSE,
                vjust = 2, hjust = 1.1, size = 3)
  }

  # Apply logarithmic scale if log_scale is TRUE
  if (log_scale) {
    plot <- plot + scale_y_continuous(
      trans = scales::log1p_trans(),
      breaks = if (!is.null(breaks)) breaks else waiver()
    )
  }

  return(plot)
}

# Helper function for faceted boxplot by ancestry
plot_box_kruskal_facet_ancestry <- function(df, count_y, group, ancestry_col, xlab, y_lab,
                                           x_order = NULL, log_scale = FALSE) {
  # Filter out non-finite values
  df_clean <- df %>% filter(is.finite(!!sym(count_y)))

  # Rename ancestry column for plotting (assign within aes)
  plot <- ggplot(df_clean, aes(x = !!sym(group), y = !!sym(count_y), fill = !!sym(group))) +
    geom_boxplot(outlier.shape = NA) +
    scale_fill_manual(values = color_palette_6) +
    geom_jitter(position=position_jitterdodge(jitter.width=0.2), color = "black", size = 1.5) +
    labs(x = xlab, y = y_lab) +
    guides(fill="none") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    facet_wrap(as.formula(paste("~", ancestry_col)), scales = "free_x")

  # Apply custom order for x-axis if provided
  if (!is.null(x_order)) {
    plot <- plot + scale_x_discrete(limits = x_order)
  }

  # Apply logarithmic scale if log_scale is TRUE
  if (log_scale) {
    plot <- plot + scale_y_continuous(
      trans = scales::log1p_trans()
    )
  }

  return(plot)
}

# Faceted version of plot_count_lm by ancestry
# Creates separate panels for each ancestry group with individual linear model tests
plot_count_lm_facet_ancestry <- function(df, chr, type, group, covariates, x_lab, y_lab,
                                         residuals = FALSE, log_scale = FALSE, breaks = NULL,
                                         min_samples = 5, fill_palette = NULL,
                                         ancestry_col = "predicted_ancestry_thres") {
  # Use the helper function to construct the combination label
  combination <- construct_combination_label(chr, type)

  # Check if ancestry column exists
  if (!ancestry_col %in% colnames(df)) {
    cat("Warning: Ancestry column '", ancestry_col, "' not found in dataframe. Cannot create faceted plot.\n", sep="")
    return(NULL)
  }

  # Check if there are any NA values in the required columns
  all_vars <- c(group, covariates, ancestry_col)
  if (anyNA(df[all_vars])) {
    cat("Data contains NA values in group, covariates, or ancestry. Excluding rows with NA values.\n")
    df <- df %>% filter(complete.cases(df[all_vars]))
  }

  if (nrow(df) == 0) {
    cat("Warning: No data remaining after removing NA values. Cannot create faceted plot.\n")
    return(NULL)
  }

  # Aggregate rare tumor types
  df <- df %>%
    group_by(tumor_type) %>%
    mutate(tumor_type = ifelse(n() < min_samples, "Other", tumor_type)) %>%
    ungroup()

  # Ensure group is a factor
  df[[group]] <- factor(df[[group]])

  # Calculate linear model p-values for each ancestry group
  ancestry_groups <- unique(df[[ancestry_col]])
  p_values_df <- data.frame()

  cat("\n=== LINEAR MODEL RESULTS BY ANCESTRY ===\n")
  for (anc in ancestry_groups) {
    df_anc <- df %>% filter(!!sym(ancestry_col) == anc)

    cat("\n--- Ancestry:", anc, "---\n")
    cat("Sample sizes:\n")
    print(table(df_anc[[group]]))

    # Initialize p_value_formatted
    p_value_formatted <- "NA"

    # Construct the formula for the linear model
    formula <- as.formula(
      paste(combination, "~", if (residuals) paste(covariates, collapse = "+") else paste(group, "+", paste(covariates, collapse = "+")))
    )

    # Fit the linear model
    tryCatch({
      lm_model <- lm(formula, data = df_anc)
      lm_summary <- summary(lm_model)
      anova_result <- anova(lm_model)

      # Get overall ANOVA F-test p-value for the group
      if (!residuals && group %in% rownames(anova_result)) {
        group_p_value <- anova_result[group, "Pr(>F)"]
        p_value_formatted <<- formatC(group_p_value, format = "e", digits = 2)
        cat("ANOVA p-value for", group, ":", p_value_formatted, "\n")
      } else {
        p_value_formatted <<- "NA"
        cat("Could not calculate p-value\n")
      }

      # Calculate medians
      medians <- df_anc %>%
        group_by(.data[[group]]) %>%
        summarise(median_count = median(!!sym(combination), na.rm = TRUE))
      cat("Medians:\n")
      print(medians)

    }, error = function(e) {
      cat("Error fitting model:", e$message, "\n")
    })

    p_values_df <- rbind(p_values_df, data.frame(
      ancestry = anc,
      p_value = p_value_formatted,
      stringsAsFactors = FALSE
    ))
  }

  # Rename ancestry column in p_values_df to match the actual column name
  colnames(p_values_df)[1] <- ancestry_col

  # Use provided palette or default to color_palette_6 (same as Kruskal plots)
  if (is.null(fill_palette)) {
    fill_palette <- color_palette_6
  }

  # Create faceted plot
  p <- ggplot(df, aes(x = .data[[group]], y = !!sym(combination), fill = .data[[group]])) +
    geom_boxplot(outlier.shape = NA, color = "black") +
    scale_fill_manual(values = fill_palette, guide = "none") +
    geom_jitter(color = "black", size = 1.5, width = 0.2) +
    labs(x = x_lab, y = y_lab) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    facet_wrap(as.formula(paste("~", ancestry_col)), scales = "free_x")

  # Add p-values as text annotations for each facet (only for non-NA values)
  p_values_df_valid <- p_values_df[p_values_df$p_value != "NA", ]
  if (nrow(p_values_df_valid) > 0) {
    p <- p +
      geom_text(data = p_values_df_valid,
                aes(x = Inf, y = Inf, label = paste("p =", p_value)),
                inherit.aes = FALSE,
                vjust = 2, hjust = 1.1, size = 3)
  }

  # Apply logarithmic scale if log_scale is TRUE
  if (log_scale) {
    p <- p + scale_y_continuous(
      trans = scales::log1p_trans(),
      breaks = if (!is.null(breaks)) breaks else waiver()
    )
  }

  return(p)
}

calc_location_wilcox <- function(sv_df, filter_var="SV.type", filter_element, group, location, gene){
  # Check if the combination is present in the dataframe
  if (is.character(location) && !any(is.na(location))){
    if (!any(sv_df$Gene_name == gene & sv_df$Location2 %in% location)) {
      #print(paste("Combination of gene:", gene, "and location:", location, "not present in dataframe"))
      return(NA)  # Return NA if the combination is not found
    }
  }else{
    if (!any(sv_df$Gene_name == gene)) {
      #print(paste("Combination of gene:", gene, "not present in dataframe"))
      return(NA)  # Return NA if the combination is not found
    }
  }
  # filter for gene name
  if (!is.na(gene)){
    sv_df <- sv_df %>% filter(Gene_name==gene)
  }
  # filter for 3' UTR
  if (is.character(location) && !any(is.na(location))){
    if (length(location) == 1){
      sv_df <- sv_df %>% filter(Location2 == location)
    } else if (length(location) == 2){
      sv_df <- sv_df %>% filter(Location2 == location[1] | Location2 == location[2])
    }
  }
  # Filter for TE type
  if (!is.na(filter_element)){ 
    sv_df <- sv_df %>% filter(!!sym(filter_var) == {{ filter_element}} )
  }
  
  # group by group and sample and count TE
  data <- sv_df %>%
    group_by(!!sym(group), sample) %>%
    summarise(count = n(), .groups = 'drop')
  
  # Check if there are at least two samples
  if (nrow(data) == 1) {
    #print(paste("Only one sample for gene:", gene, "and location:", location))
    return(NA)  # Return NA if not enough groups
  }

  # Check if there are at least two groups with actual values (not NA)
  non_na_groups <- unique(data[[group]][!is.na(data[[group]])])
  if (length(non_na_groups) < 2) {
    print(paste("Not enough groups for gene:", gene, "and location:", location, "num samples:", nrow(data)))
    print(data[[group]])
    return(NA)  # Return NA if not enough groups
  }

  # Only return p-value from Wilcox test
  return(wilcox.test(reformulate(group, "count"), data = data, exact=FALSE)$p.value)
}

count_location_wilcox <- function(te_df, gene_vector, filter_var="SV.type", filter_element, group, location) {
  # Initialize a dataframe to store the results
  results <- data.frame(gene = character(), p_value = numeric(),
                       mutant_affected = integer(), mutant_unaffected = integer(),
                       wt_affected = integer(), wt_unaffected = integer(),
                       stringsAsFactors = FALSE)

  # Get unique samples per group
  all_samples_per_group <- te_df %>%
    distinct(sample, !!sym(group)) %>%
    group_by(!!sym(group)) %>%
    summarise(all_samples = list(sample), .groups = "drop")

  # Loop through each gene
  for (i in 1:length(gene_vector)) {
    p_value <- calc_location_wilcox(te_df, filter_element=filter_element, group=group, location=location, gene=gene_vector[i])

    if (is.na(p_value)) {
      next  # Skip this iteration if the combination is not present
    }

    # Get samples affected by this gene
    gene_df <- te_df %>% filter(Gene_name == gene_vector[i])

    # Apply location filter if specified
    if (is.character(location) && !any(is.na(location))) {
      if (length(location) == 1) {
        gene_df <- gene_df %>% filter(Location2 == location)
      } else if (length(location) == 2) {
        gene_df <- gene_df %>% filter(Location2 == location[1] | Location2 == location[2])
      }
    }

    # Apply filter_element if specified
    if (!is.na(filter_element)) {
      gene_df <- gene_df %>% filter(!!sym(filter_var) == filter_element)
    }

    # Get affected samples per group
    affected_samples_per_group <- gene_df %>%
      distinct(sample, !!sym(group)) %>%
      group_by(!!sym(group)) %>%
      summarise(affected_samples = list(sample), .groups = "drop")

    # Calculate counts for each group
    group_levels <- unique(te_df[[group]])
    mutant_affected <- 0
    mutant_unaffected <- 0
    wt_affected <- 0
    wt_unaffected <- 0

    for (grp in group_levels) {
      all_samples <- all_samples_per_group %>% filter(!!sym(group) == grp) %>% pull(all_samples) %>% unlist()
      affected_samples <- affected_samples_per_group %>% filter(!!sym(group) == grp) %>% pull(affected_samples)

      if (length(affected_samples) > 0) {
        affected_samples <- unlist(affected_samples)
      } else {
        affected_samples <- character(0)
      }

      n_affected <- length(affected_samples)
      n_unaffected <- length(all_samples) - n_affected

      if (grepl("Mutant", grp, ignore.case = TRUE)) {
        mutant_affected <- n_affected
        mutant_unaffected <- n_unaffected
      } else {
        wt_affected <- n_affected
        wt_unaffected <- n_unaffected
      }
    }

    results <- rbind(results, data.frame(
      gene = gene_vector[i],
      mutant_affected = mutant_affected,
      mutant_unaffected = mutant_unaffected,
      wt_affected = wt_affected,
      wt_unaffected = wt_unaffected,
      p_value = p_value
    ))
  }

  # Correct p-values for multiple testing
  results$fdr <- p.adjust(results$p_value, method = "fdr")

  # Return the results dataframe
  return(results)
}

downsample_df <- function(df) {
  # Splitting the dataframe into LFS and Control groups
  lfs_group <- df %>% filter(TP53_status == "LFS")
  control_group <- df %>% filter(TP53_status == "Control")
  
  # Downsampling the Control group
  # Adjust the size parameter as needed. Here it matches the size of the LFS group
  downsampled_control <- control_group %>% sample_n(size = nrow(lfs_group), replace = FALSE)
  
  # Combining back the LFS and the downsampled Control groups
  combined_df <- bind_rows(lfs_group, downsampled_control)
  
  return(combined_df)
}

identify_te_genes <- function(gene_df, gene_vector = NULL, location2 = NULL, location = NULL) {
  
  # Filter by gene_vector if specified
  if (!is.null(gene_vector)) {
    gene_df <- gene_df %>%
      filter(Gene_name %in% gene_vector)
  }
  
  # Filter by location2 if specified
  if (!is.null(location2)) {
    gene_df <- gene_df %>%
      filter(Location2 %in% location2)
  }
  
  # Filter by location (string match) if specified
  if (!is.null(location)) {
    gene_df <- gene_df %>%
      filter(grepl(location, Location, ignore.case = TRUE))
  }
  
  # Return the filtered data frame
  return(data.frame(gene_df))
}

gene_fisher<- function(sv_df, gene_names){
  results <- list()
  
  # Perform Fisher's test for each gene
  for (gene_name in gene_names) {
    # Filter for the specific gene
    gene_data <- sv_df %>% filter(Gene_name == gene_name)
    
    # Skip gene if not present in dataframe
    if (nrow(gene_data) == 0) {
      next
    }
    
    # Create a summary for each sample: is the gene mutated (at least once)?
    sample_summary <- gene_data %>%
      group_by(sample) %>%
      summarise(Gene_Mutated = n() >= 1, .groups = 'drop') %>%
      right_join(sv_df %>% distinct(sample, TP53_status), by = 'sample')
    
    # Replace NA with FALSE for Gene_Mutated (assumes NA means not mutated)
    sample_summary$Gene_Mutated[is.na(sample_summary$Gene_Mutated)] <- FALSE
    
    # Check if there are any mutated instances of the gene
    if (sum(sample_summary$Gene_Mutated) == 0) {
      next  # Skip to the next gene if there are no mutated instances
    }
    
    # Create the contingency table
    contingency_table <- table(sample_summary$TP53_status, sample_summary$Gene_Mutated)
    #print(contingency_table)
    # Check if the table has the right dimensions for Fisher's test
    if (all(dim(contingency_table) >= 2)) {
      # Perform Fisher's Exact Test
      fisher_test_result <- fisher.test(contingency_table)
      # Store results in the list
      results[[gene_name]] <- list(contingency_table = contingency_table, 
                                   p_value = fisher_test_result$p.value)
    }else{
      results[[gene_name]] <- list(contingency_table = contingency_table, 
                                   p_value = NA)
    }
  }
  
  # Create a dataframe from the results
  results_df <- do.call(rbind, lapply(names(results), function(gene_name) {
    data.frame(Gene = gene_name, 
               Contingency_Table = I(list(results[[gene_name]]$contingency_table)), 
               P_Value = results[[gene_name]]$p_value)
  }))
  
  # Adjust p-values for multiple testing using Benjamini-Hochberg method
  results_df$FDR <- p.adjust(results_df$P_Value, method = "BH")
  
  return(results_df)
}

generate_plots_perchrom <- function(plot_function, df, types, group) {
  # Use lapply to iterate over each type and print the corresponding plot
  invisible(lapply(types, function(type) {
    # Dynamically set y_lab based on the type
    y_lab <- if (is.na(type)) "Total TE count normalized by chromosome length" else paste(type, "count normalized by chromosome length")
    
    # Call the provided plotting function and print the plot
    print(plot_function(df, type = type, group = group, y_lab = y_lab))
  }))
}

# plot box plot per chromosome
plot_box_perchr <- function(df, group, y_lab){
  ggplot(df, aes(chr, normalized_count, fill = !!sym(group))) +
    geom_boxplot(position = position_dodge(width = 0.8), outlier.color = "grey") +
    scale_fill_manual(values =  colours) +
    #geom_jitter(position=position_jitterdodge(jitter.width=0.2), color = "black", size = 1.5) +
    labs(x = "Chromosome", y = y_lab) 
}

# plot number of TE per chromsomes
plot_count_perchr <- function(df, type = NA, group, y_lab, chr_length = chr_lengths) {
  # Determine the columns to select based on the presence of type
  if (!is.na(type)) {
    columns_to_select <- grep(paste0("_", type, "$"), colnames(df), value = TRUE)
  } else {
    columns_to_select <- grep("^chr[0-9XY]+$", colnames(df), value = TRUE)
  }
  
  if (length(columns_to_select) == 0) {
    stop("Specified combination columns not found in the data frame.")
  }
  
  # Select the relevant columns and pivot to long format
  columns_to_select <- c("sample", group, columns_to_select)
  data <- df %>%
    select(all_of(columns_to_select)) %>%
    pivot_longer(cols = starts_with("chr"), names_to = "chr", values_to = "count") %>%
    mutate(chr = str_remove(chr, "^chr")) %>%
    mutate(chr = str_remove(chr, "_.*$"))

  # Normalize the counts
  data <- as.data.frame(data)
  
  # Convert chr_length to data.frame if it's a named vector
  if (is.vector(chr_length)) {
    chr_length <- data.frame(
      chr = names(chr_length),
      length = as.numeric(chr_length),
      stringsAsFactors = FALSE
    )
  } else {
    chr_length <- as.data.frame(chr_length)
    if (!"length" %in% colnames(chr_length)) {
      # Try to find the length column by other names
      length_col <- grep("length|Length|LENGTH", colnames(chr_length), value = TRUE)[1]
      if (!is.na(length_col)) {
        chr_length$length <- chr_length[[length_col]]
      } else {
        stop("Cannot find length column in chr_length")
      }
    }
    chr_length$length <- as.numeric(chr_length$length)
  }
  
  data$chr <- as.character(data$chr)
  chr_length$chr <- as.character(chr_length$chr)
  
  data <- data %>%
    left_join(chr_length, by = "chr") %>%
    mutate(normalized_count = count / length)
  
  # Set factor levels for chromosome ordering
  data$chr <- factor(data$chr, levels = chr_length$chr)
  
  # Collecting p-values for each comparison
  p_values <- data.frame(chromosome = character(), group1 = character(), group2 = character(), p_value = numeric())
  for (chrom in unique(data$chr)) {
    sub_data <- subset(data, chr == chrom)
    group_levels <- levels(factor(sub_data[[group]]))
    if (length(group_levels) >= 2) {
      combn(group_levels, 2, function(x) {
        group1_data <- sub_data %>% 
          filter(!!sym(group) == x[1]) %>% 
          pull(normalized_count)
        group2_data <- sub_data %>% 
          filter(!!sym(group) == x[2]) %>% 
          pull(normalized_count)
        test_result <- wilcox.test(group1_data, group2_data)
        p_values <<- rbind(p_values, data.frame(chromosome = chrom, group1 = x[1], group2 = x[2], p_value = test_result$p.value))
      }, simplify = FALSE)
    }
  }
  
  # Apply BH correction
  p_values$p_value_corrected <- p.adjust(p_values$p_value, method = "BH")
  
  # Plot
  plot_data <-plot_box_perchr(data, group, y_lab)
  
  # Find y position for each annotation
  max_values <- aggregate(normalized_count ~ chr, data, max)
  names(max_values) <- c("chromosome", "y_max")
  
  # Add annotations
  for (row in 1:nrow(p_values)) {
    ann <- p_values[row, ]
    y_pos <- max_values[max_values$chromosome == ann$chromosome, "y_max"] * 1.05
    p_label <- as.character(ifelse(ann$p_value_corrected < 0.001, "***", ifelse(ann$p_value_corrected < 0.01, "**", ifelse(ann$p_value_corrected < 0.05, "*", ""))))
    
    # Calculate the x-position for the annotation based on the chromosome
    x_pos <- which(levels(data$chr) == ann$chromosome)
    
    # Only add annotation if chromosome is found in the levels and has significant p-value
    if (length(x_pos) > 0 && length(y_pos) > 0 && !is.na(y_pos) && !is.na(p_label) && p_label != "") {
      plot_data <- plot_data + annotate("text", x = x_pos, y = y_pos, label = p_label, size = 16/.pt, vjust = -0.5)
    }
  }
  
  return(plot_data)
}

# plot box plot per chromosome
plot_box_perchr_notest <- function(df, chrom, count_y, y_lab){
  ggplot(df, aes(!!sym(chrom), !!sym(count_y))) +
    geom_boxplot(fill = "#5FBFF9", outlier.color = "grey", position = position_dodge(width = 0.8)) +
    #geom_jitter(alpha=0.8, color = "grey50", size = 1.5) +
    labs(x = "Chromosome", y = y_lab) 
}

# plot number of TE per chromsomes
plot_count_perchr_notest <- function(df, chr_length = chr_lengths, type = NA, y_lab, log_scale = FALSE) {
  # Determine the columns to select based on the presence of filter_element
  if (!is.na(type)) {
    columns_to_select <- grep(paste0("_", type, "$"), colnames(df), value = TRUE)
  } else {
    columns_to_select <- grep("^chr[0-9XY]+$", colnames(df), value = TRUE)
  }

  if (length(columns_to_select) == 0) {
    stop("Specified combination columns not found in the data frame.")
  }

  # Select the relevant columns and pivot to long format
  data <- df %>%
    select(sample, all_of(columns_to_select)) %>%
    pivot_longer(-sample, names_to = "chr", values_to = "count") %>%
    mutate(
      chr = str_remove(chr, "^chr"),
      chr = str_remove(chr, "_.*$")
    )

  # Normalize the counts
  # Convert chr_length to data.frame if it's a named vector
  if (is.vector(chr_length)) {
    chr_length <- data.frame(
      chr = names(chr_length),
      length = as.numeric(chr_length),
      stringsAsFactors = FALSE
    )
  } else {
    chr_length <- as.data.frame(chr_length)
    if (!"length" %in% colnames(chr_length)) {
      # Try to find the length column by other names
      length_col <- grep("length|Length|LENGTH", colnames(chr_length), value = TRUE)[1]
      if (!is.na(length_col)) {
        chr_length$length <- chr_length[[length_col]]
      } else {
        stop("Cannot find length column in chr_length")
      }
    }
    chr_length$length <- as.numeric(chr_length$length)
  }
  
  data$chr <- as.character(data$chr)
  chr_length$chr <- as.character(chr_length$chr)
  data$count <- as.numeric(data$count)
  
  data <- data %>%
    left_join(chr_length, by = "chr") %>%
    mutate(normalized_count = count / length)

  # Set factor levels for chromosome ordering
  data$chr <- factor(data$chr, levels = chr_length$chr)

  # Handle log scale
  if (log_scale) {
    epsilon <- 1e-6
    data <- data %>%
      mutate(normalized_count = normalized_count + epsilon)
  }

  # Plot using the custom function
  plot_data <- plot_box_perchr_notest(data, "chr", "normalized_count", y_lab)

  # Add log scale if needed
  if (log_scale) {
    plot_data <- plot_data +
      scale_y_log10(labels = scales::label_log(base = 10)) +
      labs(y = paste("Log10(", y_lab, " + ", epsilon, ")"))
  }

  return(plot_data)
}


# plot box plot per tumor type
plot_box_tt <- function(df, count_y, group, y_lab, legend_lab=NULL, log_scale=FALSE){
  plot_data <- ggplot(df, aes(tumor_type, !!sym(count_y), fill = !!sym(group))) +
    geom_boxplot(outlier.shape = NA) +
    scale_fill_manual(values = colours) +
    geom_jitter(position=position_jitterdodge(jitter.width=0.2), color = "black", size = 1.5) +
    labs(x = "Tumor type", y = y_lab)
  
  # Add log scale to the plot if requested
  # Apply logarithmic scale if log_scale is TRUE
  if (log_scale) {
    plot_data <- plot_data + scale_y_continuous(
      trans = scales::log1p_trans(),
    )
  }
  
  # Set custom legend title only if specified
  if (!is.null(legend_lab)) {
    plot_data <- plot_data + labs(fill = legend_lab)
  }
  
  return(plot_data)

}

# plot number of tes total per group per tumor type
plot_count_tt <- function(df, chr, type, group, min, y_lab, legend_title=NULL, log_scale=FALSE){
  # Use the helper function to construct the combination label
  combination <- construct_combination_label(chr, type)
  
  # Check if there are any NA values in the group column
  if (any(is.na(df[[group]]))) {
    cat("Group contains NA values. Excluding NA group from the plot and Wilcoxon test.\n")
    df <- df %>% filter(!is.na(!!sym(group)))
  }
  
  # Identify the two unique groups in the 'group' variable
  unique_groups <- df %>%
    distinct(!!sym(group)) %>%
    pull()
  
  # Extract samples for each group and filter by minimum sample size per tumor type
  group1_samples <- df %>%
    filter(!!sym(group) == unique_groups[1]) %>%
    group_by(tumor_type) %>%
    filter(n() >= min) %>%
    ungroup()
  
  group2_samples <- df %>%
    filter(!!sym(group) == unique_groups[2]) %>%
    group_by(tumor_type) %>%
    filter(n() >= min) %>%
    ungroup()
  print(unique(group1_samples$tumor_type))
  print(unique(group2_samples$tumor_type))
  
  # Find tumor types present in both control_samples and lfs_samples
  valid_tumor_types <- intersect(group1_samples$tumor_type, group2_samples$tumor_type)
  
  # Filter original data for these valid tumor types
  df_filtered <- df %>%
    filter(tumor_type %in% valid_tumor_types) 
  
  # Collecting p-values for each tumor type
  p_values <- data.frame(tumor_type = character(), p_value = numeric())
  
  for (tumor in unique(df_filtered$tumor_type)) {
    sub_data <- filter(df_filtered, tumor_type == tumor)
    
    # Check if TP53_status has exactly two levels
    if (length(unique(sub_data[[group]])) == 2) {
      print("dependent ~ independent")
      formula <- reformulate(group, combination)
      print(formula)
      p_value <- wilcox.test(formula, data = sub_data)$p.value
    } else {
      p_value = NA  # Assign NA if not exactly two levels
    }
    
    p_values <- rbind(p_values, data.frame(tumor_type = tumor, p_value = p_value))
  }
  
  # Apply BH correction
  p_values$p_value_corrected <- p.adjust(p_values$p_value, method = "BH")
  
  # Print the p-values
  print(p_values)
  
  # Plot
  plot<- plot_box_tt(df_filtered, combination, group, y_lab, legend_lab=legend_title, log_scale=log_scale)
  
  # Find y position for each annotation
  max_values <- df_filtered %>%
    group_by(tumor_type) %>%
    summarise(y_max = max(!!sym(combination), na.rm = TRUE)) %>%
    ungroup()
  
  # Add annotations
  for (row in seq_len(nrow(p_values))) {
    ann <- p_values[row, ]
    y_pos <- max_values %>% filter(tumor_type == ann$tumor_type) %>% pull(y_max) * 1.05
    label <- dplyr::case_when(
      ann$p_value_corrected < 0.001 ~ "***",
      ann$p_value_corrected < 0.01 ~ "**",
      ann$p_value_corrected < 0.05 ~ "*",
      TRUE ~ "n.s."
    )
    x_pos <- which(levels(factor(df_filtered$tumor_type)) == ann$tumor_type)    
    
    plot<- plot+ 
      annotate("text", x = x_pos, y = y_pos, label = label, vjust = -0.5, size = 5) +
      coord_cartesian(clip = 'off')
  }
  
  return(plot)
}

# plot number of tes total per group per tumor type
plot_count_kruskal_nogroup <- function(df,  column, min, chr, type, x_lab, y_lab, log_scale = FALSE, breaks=NULL) {
  # Use the helper function to construct the combination label
  combination <- construct_combination_label(chr, type)
  
  df_filtered <- df %>%
    dplyr::select(sample, !!sym(column), !!sym(combination)) %>%
    filter(!is.na(!!sym(column)), !!sym(column) != "NANANA") %>%  # Exclude rows with "NANANA"
    group_by(!!sym(column)) %>%
    filter(n() >= min) %>%  # Keep groups with at least n samples
    ungroup()
  
  # Calculate medians for each group and reorder
  medians <- df_filtered %>%
    group_by(!!sym(column)) %>%
    summarize(median_value = median(!!sym(combination), na.rm = TRUE)) %>%
    arrange(desc(median_value))
  
  df_filtered <- df_filtered %>%
    mutate(!!sym(column) := factor(!!sym(column), levels = medians[[column]]))  # Reorder column
  
  # Check if there are at least 2 groups for Kruskal-Wallis test
  unique_groups <- unique(df_filtered[[column]])
  unique_groups <- unique_groups[!is.na(unique_groups)]
  
  if (length(unique_groups) < 2) {
    cat("Kruskal-Wallis test requires at least 2 groups. Column:", column, "has", length(unique_groups), "unique value(s):", paste(unique_groups, collapse = ", "), ". Skipping statistical test.\n")
    p_value <- NA
  } else {
    # Perform Kruskal-Wallis test
    formula <- reformulate(column, combination)
    print("dependent ~ independent")
    print(formula)
    p_value <- kruskal.test(formula, data = df_filtered)$p.value
  }
  
  # Format the p-value for display
  if (is.na(p_value)) {
    p_value_formatted <- "NA"
  } else {
    p_value_formatted <- formatC(p_value, format = "e", digits = 2)
  }
  
  # Plot
  p <- ggplot(df_filtered, aes_string(x = column, y = combination, fill = column)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(color = "black", size = 1.5) +
    labs(x = x_lab, y = y_lab) +
    guides(fill = "none")
  
  # Annotate the plot with the p-value
  p <- p +
    annotate("text", x = Inf, y = Inf, label = paste("p =", p_value_formatted), vjust = 2, hjust = 1, size = 4)
  
  # Apply logarithmic scale if log_scale is TRUE
  if (log_scale) {
    p <- p + scale_y_continuous(
      trans = scales::log1p_trans(),
      breaks = if (!is.null(breaks)) breaks else waiver()
    )
  }
  
  return(p)
}

# plot number of tes total per group
plot_size_box <- function(sv_df, filter_var, filter_element, group, x_lab, y_lab){
  if (!is.na(filter_element)){ # if need to filter by SV type
    data <- sv_df %>%
      filter(!!sym(filter_var) == {{ filter_element }} ) %>% # filter for deletions
      group_by(!!sym(group), sample) # group
  }
  else { # if doing overall TE count
    data <- sv_df %>%
      group_by(!!sym(group), sample) # group
  }
  
  # print p value
  print(wilcox.test(reformulate(group, "SV.length"), data = data)$p.value)
  
  # Plot
  plot <- plot_box(data, SV.length, !!sym(group), x_lab, y_lab) + 
    geom_signif(test = "wilcox.test", # sig test
                comparisons = list(levels(factor(sv_df[[group]]))),
                map_signif_level = TRUE) + # stars 
    coord_cartesian(clip = 'off') # dont cut off annotation
  return(plot)
}

# plot number of tes total per group
plot_size_box_avg <- function(sv_df, filter_var, filter_element, group, x_lab, y_lab){
  if (!is.na(filter_element)){ # if need to filter by SV type
    data <- sv_df %>%
      filter(!!sym(filter_var) == {{ filter_element }} ) %>% # filter for deletions
      group_by(!!sym(group), sample) %>% # group
      summarise(avg_length = mean(SV.length), .groups = 'drop') # count
  }
  else { # if doing overall TE count
    data <- sv_df %>%
      group_by(!!sym(group), sample) %>% # group
      summarise(avg_length = mean(SV.length), .groups = 'drop') # count
  }
  
  # print p value
  print(wilcox.test(reformulate(group, "avg_length"), data = data)$p.value)
  
  # Plot
  plot <- plot_box(data, avg_length, !!sym(group), x_lab, y_lab) + 
    geom_signif(test = "wilcox.test", # sig test
                comparisons = list(levels(factor(sv_df[[group]]))),
                map_signif_level = TRUE) + # stars 
    coord_cartesian(clip = 'off') # dont cut off annotation
  return(plot)
}

# plot number of sv total for LFS affected and unaffected
plot_count_lfs_affected <- function(sv_df, filter_var, x_lab,  y_lab){
  if (!is.na(filter_var)){ # if need to filter by SV type
    data <- sv_df %>%
      filter(SV.type == {{ filter_var}} ) %>% # filter for deletions
      group_by(Cancer, sample) %>% # group
      summarise(count = n(), .groups = 'drop') # count
  }
  else { # if doing overall TE count
    data <- sv_df %>%
      group_by(Cancer, sample) %>% # group
      summarise(count = n(), .groups = 'drop') # count by group
  }
  
  # print p value
  print(wilcox.test(reformulate("Cancer", "count"), data = data)$p.value)
  
  # Plot
  plot <- plot_box(data, count, Cancer, x_lab, y_lab) + 
    geom_signif(test = "wilcox.test", # sig test
                comparisons = list(levels(factor(data$Cancer))),
                map_signif_level = TRUE) + # stars 
    coord_cartesian(clip = 'off') # dont cut off annotation
  return(plot)
}

# plot number of sv by age at diagnosis
plot_count_age <- function(df, type, chr, y_lab){
  # Ensure age_at_diagnosis is numeric and filter out NA
  df <- df %>%
    mutate(age_at_diagnosis = as.numeric(as.character(age_at_diagnosis))) %>%
    filter(!is.na(age_at_diagnosis)) %>%
    mutate(age_at_diagnosis = age_at_diagnosis/365)
  
  # Use the helper function to construct the combination label
  combination <- construct_combination_label(chr, type)
  
  # Fit linear model
  formula <- reformulate("age_at_diagnosis", combination)
  print("dependent ~ independent")
  print(formula)
  model <- lm(formula, data = df)
  summary_model <- summary(model)
  print(summary_model)
  
  # Extract R-squared value, p-value, and coefficient
  r_squared <- sprintf("%.2e", summary_model$r.squared)
  p_value <- sprintf("%.2e", summary_model$coefficients[2, "Pr(>|t|)"])
  coefficient <- sprintf("%.2e", summary_model$coefficients[2, "Estimate"])
  
  # Plotting
  p <- ggplot(df, aes(x = age_at_diagnosis, y = !!sym(combination))) +
    geom_point() +
    geom_smooth(method="lm", se=FALSE, color="#5FBFF9") + 
    scale_x_continuous(breaks = seq(0, 30, by = 10), limits = c(0, 30)) + 
   # scale_y_continuous(limits = c(0, 4000)) + 
    labs(x = "Age at Diagnosis",
         y = y_lab)  
  
  # Add annotation with R-squared, p-value, and coefficient
  p <- p +
    annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.1,
             label = sprintf("R²: %s\np-value: %s\nCoefficient: %s", r_squared, p_value, coefficient),
             size = 3)
  return(p)
}

# plot number of sv by age at diagnosis
plot_count_age_enrollment <- function(sv_df, filter_var, filter_element, y_lab){
  # Ensure age_at_enrollment is numeric
  sv_df <- sv_df %>%
    mutate(age_at_enrollment = as.numeric(as.character(age_at_enrollment)))
  
  if (!is.na(filter_element)){ # if need to filter by SV type
    data <- sv_df %>%
      filter(!!sym(filter_var) == {{ filter_element}} ) %>% # filter for deletions
      group_by(age_at_enrollment, sample) %>% # group
      summarise(count = n(), .groups = 'drop') # count
  }
  else { # if doing overall TE count
    data <- sv_df %>%
      group_by(age_at_enrollment, sample) %>% # group
      summarise(count = n(), .groups = 'drop') # count by group
  }
  
  # Fit linear model
  model <- lm(count ~ age_at_enrollment, data = data)
  summary_model <- summary(model)
  print(summary_model)
  
  # Extract R-squared value, p-value, and coefficient
  r_squared <- sprintf("%.2e", summary_model$r.squared)
  p_value <- sprintf("%.2e", summary_model$coefficients[2, "Pr(>|t|)"])
  coefficient <- sprintf("%.2e", summary_model$coefficients[2, "Estimate"])
  
  # Plotting
  p <- ggplot(data, aes(x = age_at_enrollment, y = count)) +
    geom_point() +
    geom_smooth(method="lm", se=FALSE, color="#5FBFF9") + 
    scale_x_continuous(breaks = seq(0, 30, by = 10)) + 
    labs(x = "Age at Diagnosis",
         y = y_lab) 
  
  # Add annotation with R-squared, p-value, and coefficient
  p <- p +
    annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.1,
             label = sprintf("R²: %s\np-value: %s\nCoefficient: %s", r_squared, p_value, coefficient),
             size = 3)
  return(p)
}

# plot number of sv by age at diagnosis
plot_length_age <- function(df, type, y_lab){
  df <- df %>%
    mutate(age_at_diagnosis = as.numeric(as.character(age_at_diagnosis))) %>%
    filter(!is.na(age_at_diagnosis)) %>%
    mutate(age_at_diagnosis = age_at_diagnosis/365)
  
  if (!is.na(type)){ # if need to filter by SV type
    df_filtered <- df %>%
      filter(ALT == type) %>% 
      group_by(sample, age_at_diagnosis) %>% # group
      summarise(avg_length = mean(SV_length), .groups = 'drop') # count
  }
  else { # if doing overall TE count
    df_filtered <- df %>%
      group_by(sample, SV_length, age_at_diagnosis) %>% # group
      summarise(avg_length = mean(SV_length), .groups = 'drop') # count
  }
  
  # Fit linear model
  model <- lm(avg_length ~ age_at_diagnosis, data = df_filtered)
  summary_model <- summary(model)
  print(summary_model)
  
  # Extract R-squared value, p-value, and coefficient
  r_squared <- sprintf("%.2e", summary_model$r.squared)
  p_value <- sprintf("%.2e", summary_model$coefficients[2, "Pr(>|t|)"])
  coefficient <- sprintf("%.2e", summary_model$coefficients[2, "Estimate"])
  
  
  # Plotting
  p <- ggplot(df_filtered, aes(x = age_at_diagnosis, y = avg_length)) +
    geom_point() +
    geom_smooth(method="lm", se=FALSE, color="#5FBFF9") + 
    scale_x_continuous(breaks = seq(0, 30, by = 10), limits = c(0, 30)) + 
    labs(x = "Age at Diagnosis",
         y = y_lab)
  
  # Add annotation with R-squared, p-value, and coefficient
  p <- p +
    annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.1,
             label = sprintf("R²: %s\np-value: %s\nCoefficient: %s", r_squared, p_value, coefficient),
             size = 3)
  return(p)
  
}

# size of sv
plot_size_density <- function(df, n, filter_var, filter_element, group, x_lab) {
  if (!is.na(filter_element)){ # if need to filter by SV type
    data <- df %>%
      filter(!!sym(filter_var) == {{ filter_element }} ) %>% # filter for deletions
      filter(SV.length<n) # filter size 
  }
  else { # if doing overall TE count
    data <- df %>%
      filter(SV.length<n) # filter size 
  }
  
  ggplot(data, aes(x = SV.length, fill = !!sym(group))) +
    geom_density(alpha = 0.5) +
    labs(x = x_lab,
         y = "Density") +
    scale_fill_manual(values = colours)
}

# size of sv no group for fill
plot_size_density_nogroup <- function(df, n, filter_var, filter_element, x_lab) {
  if (!is.na(filter_element)){ # if need to filter by SV type
    data <- df %>%
      filter(!!sym(filter_var) == {{ filter_element }} ) %>% # filter for deletions
      filter(SV.length<n) # filter size 
  }
  else { # if doing overall TE count
    data <- df %>%
      filter(SV.length<n) # filter size 
  }
  
  ggplot(data, aes(x = SV.length)) +
    geom_density(alpha = 0.5) +
    labs(x = x_lab,
         y = "Density") +
    scale_fill_manual(values = "#5FBFF9")
}

# size of te line graph
# DOESNT CONTROL FOR NUMBER OF SAMPLES IN GROUP SO USE DENSITY
plot_size_line <- function(df, n, filter_var, filter_element, xlab){
  if (!is.na(filter_element)){ # if need to filter by SV type
    df <- df %>%
      filter(!!sym(filter_var) == {{ filter_element }} ) # filter for deletions
  }
  
  df_summary <- df %>%
    filter(SV.length<n) %>%
    group_by(SV.length) %>%
    summarise(count = n())
  
  # Plot the data
  ggplot(df_summary, aes(x = SV.length, y = count)) +
    geom_line() +
    labs(x = xlab,
         y = "Count")
}

# size of te line graph
plot_size_line_allsv <- function(df, filter_length = NULL, log_x = FALSE) {
  if (!is.null(filter_length)) {
    df <- df %>%
      filter(SV_length < filter_length) %>%
      group_by(SV_length) %>%
      summarise(count = n(), .groups = 'drop')
  } else {
    df <- df %>% 
      group_by(SV_length) %>%
      summarise(count = n(), .groups = 'drop')
  }
  
  # Create the base plot
  p <- ggplot(df, aes(x = SV_length, y = count)) +
    geom_line(color="#5FBFF9") +
    labs(x = ifelse(log_x, "TE Length (log scale)", "TE Length"),
         y = "Count")
  
  # Apply logarithmic scale to the x-axis if log_x is TRUE
  if (log_x) {
    p <- p + scale_x_log10(
      labels = trans_format("log10", math_format(10^.x)),
      breaks = function(x) {
        exp(seq(ceiling(log10(min(x))), floor(log10(max(x))), by = 1) * log(10))
      }
    ) +
      annotation_logticks(
        sides = "b", 
        outside=TRUE) + 
      coord_cartesian(clip = "off")
  }
  
  return(p)
}

plot_count_line_allchr<- function(df, type, y_lab){
  if (!is.na(type)){ # if need to filter by SV type
    df <- df %>%
      filter(ALT == type) # filter for deletions
  }
  
  # Bin the SV_start values into 1 million bp bins
  df_summary <- df %>%
    mutate(SV_start_bin = cut_width(SV_start, 1e7, boundary=0)) %>%
    group_by(SV_chrom, SV_start_bin, ALT) %>%
    summarise(count = n(), .groups = 'drop')
  
  # Create bin labels that represent the start of each bin with chromosome numbers
  df_summary <- df_summary %>%
    mutate(SV.chrom_num = as.numeric(as.character(SV_chrom)), # Convert chromosome to numeric
           SV_start_bin_num = as.numeric(gsub(",.*", "", gsub("\\(|\\[", "", SV_start_bin)))) %>%
    arrange(SV.chrom_num, SV_start_bin_num) %>%
    mutate(bin_label = paste(SV_chrom, SV_start_bin_num, sep = "-"))
  
  # Add a row number column
  df_summary <- df_summary %>%
    mutate(row_number = row_number())
  
  # Calculate middle positions for each SV.chrom
  middle_positions <- df_summary %>%
    group_by(SV_chrom) %>%
    summarise(middle = mean(row_number))
  
  # Plot line graph 
  line_plot <- ggplot(df_summary, aes(x = row_number, y = count, group = ALT, color = ALT)) +
    geom_line() +
    scale_x_continuous(breaks = middle_positions$middle, labels = middle_positions$SV_chrom, expand = c(0, 0)) +
    scale_fill_manual(values = colour_palette_3, name = NULL) +
    labs(x = "Chromosome",
         y = y_lab)  + 
    theme(axis.text.x = element_text(hjust = 0.5)) # Rotate labels for readability
  
  # Plot bar graph 
  bar_plot <- ggplot(df_summary, aes(x = row_number, y = count, group = ALT, fill= ALT)) +
    geom_bar(stat = "identity", position = position_dodge(width = 1.0)) +
    scale_fill_manual(values = colour_palette_3, name = NULL) + # comment out if multiple colours
    scale_x_continuous(breaks = middle_positions$middle, labels = middle_positions$SV_chrom, expand = c(0, 0)) +
    labs(x = "Chromosome",
         y = y_lab)  + 
    theme(axis.text.x = element_text(hjust = 0.5)) # Rotate labels for readability
  
  if (!is.na(type)){
    line_plot <- line_plot + guides(color = "none")
    bar_plot <- bar_plot + guides(fill = "none") # Remove legend if filter_element is specified
  }
  
  return(list(line_plot = line_plot, bar_plot = bar_plot))
}

plot_count_line_allchr_group <- function(df, group, y_lab){
  # Calculate group sizes
  group_sizes <- df %>%
    group_by(!!sym(group)) %>%
    summarise(group_size = n_distinct(sample), .groups = 'drop')
  
  # If the group is 'tumor_type', filter for types with at least 5 samples
  if(group == "tumor_type") {
    group_sizes <- group_sizes %>%
      filter(group_size >= 5)
  }
  
  # Proceed only with filtered groups
  df <- df %>%
    inner_join(group_sizes, by = group)
  
  # Bin the SV_start values into 1 million bp bins
  df_summary <- df %>%
    mutate(SV_start_bin = cut_width(SV_start, 1e7, boundary=0)) %>%
    group_by(!!sym(group), SV.chrom, SV_start_bin) %>%
    summarise(count = n(), .groups = 'drop') 
  
  # Join group sizes to df_summary
  df_summary <- df_summary %>%
    left_join(group_sizes, by = group)
  
  # Ensure count and group_size are numeric
  df_summary <- df_summary %>%
    mutate(count = as.numeric(count),
           group_size = as.numeric(group_size),
           normalized_count = count / group_size)
  
  # Create bin labels that represent the start of each bin with chromosome numbers
  df_summary <- df_summary %>%
    mutate(SV.chrom_num = as.numeric(as.character(SV.chrom)), # Convert chromosome to numeric
           SV_start_bin_num = as.numeric(gsub(",.*", "", gsub("\\(|\\[", "", SV_start_bin)))) %>%
    arrange(SV.chrom_num, SV_start_bin_num) %>%
    mutate(bin_label = paste(SV.chrom, SV_start_bin_num, sep = "-"))
  
  # Add a row number column
  df_summary <- df_summary %>%
    mutate(row_number = row_number())
  
  # Calculate middle positions for each SV.chrom
  middle_positions <- df_summary %>%
    group_by(SV.chrom) %>%
    summarise(middle = mean(row_number))
  
  # Plot the data
  line_plot <- ggplot(df_summary, aes(x = row_number, y = normalized_count, group = !!sym(group), color = !!sym(group))) +
    geom_line() +
    scale_x_continuous(breaks = middle_positions$middle, labels = middle_positions$SV.chrom, expand = c(0, 0)) +
    labs(x = "Chromosome",
         y = "Normalized TE count per 10 Mb")  + 
    theme(axis.text.x = element_text(hjust = 0.5)) # Rotate labels for readability
  
  # Plot bar graph 
  bar_plot <- ggplot(df_summary, aes(x = row_number, y = normalized_count, group = !!sym(group), fill= !!sym(group))) +
    geom_bar(stat = "identity", position = position_dodge(width = 1.0)) +
    scale_x_continuous(breaks = middle_positions$middle, labels = middle_positions$SV.chrom, expand = c(0, 0)) +
    labs(x = "Chromosome",
         y = "TE count per 10 Mb")  + 
    theme(axis.text.x = element_text(hjust = 0.5)) # Rotate labels for readability
  
  return(list(line_plot = line_plot, bar_plot = bar_plot))
}

# number across genome of te line graph
plot_count_line_chr<- function(df, type, chr, chr_length){
  if (!is.na(type)){ # if need to filter by SV type
    df <- df %>%
      filter(ALT == type) # filter for deletions
  }
  
  # Bin the SV_start values into 1 million bp bins
  df_summary <- df %>%
    filter(SV_chrom == chr) %>%
    mutate(SV_start_bin = cut_width(SV_start, 1e6, boundary=0)) %>%
    group_by(SV_start_bin, ALT) %>%
    summarise(count = n(), .groups = 'drop')
  
  # Define specific x-axis labels
  x_labels <- seq(1, chr_length, length.out = 5) / 1e6 # Convert to Mb
  
  # Plot the data
  p <- ggplot(df_summary, aes(x = SV_start_bin, y = count, group = ALT, color = ALT)) +
    geom_line() +
    scale_x_discrete(labels = x_labels, breaks = x_labels) +
    labs(x = paste0("Chromosome ", chr, " (Mb)"),
         y = "TE count per Mb")  + 
    theme(axis.text.x = element_text(angle = 45, hjust = 0.5)) # Rotate labels for readability
  
  if (!is.na(type)){
    p <- p + guides(color="none")
  }
  
  return(p)
}



# PCA
# format
count_overlaps <- function(df){
  df$window_chr<- gsub("chr", "", df$window_chr)
  
  # Create a unique identifier for each window
  bed_data <- df %>%
    mutate(window = paste(window_chr, window_start, window_end, sep = "_"))
  
  # Count the number of overlaps for each sample and window
  overlap_counts <- bed_data %>%
    group_by(sample, window) %>%
    summarise(overlap_count = n(), .groups = 'drop')
  
  # Reshape the data to have one row per sample and one column per window
  final_df <- pivot_wider(overlap_counts, names_from = window, values_from = overlap_count, values_fill = 0)
  
  return(final_df)
}

# sort columns numerically
sort_numeric_columns <- function(cols) {
  # Filter out columns that start with 'X_' or 'Y_'
  numeric_cols <- cols[!grepl("^(X|Y)_", cols)]
  
  # Split the column names and convert to numeric for sorting
  numeric_parts <- lapply(strsplit(numeric_cols, "_"), function(x) c(as.numeric(x[1]), as.numeric(x[2])))
  numeric_parts <- do.call(rbind, numeric_parts)
  
  # Return the sorted numeric columns based on the numeric parts
  numeric_cols[order(numeric_parts[, 1], numeric_parts[, 2])]
}

# sort x y columns
sort_X_Y_columns <- function(cols, prefix) {
  # Filter and sort X or Y columns
  xy_cols <- cols[grepl(paste0("^", prefix, "_"), cols)]
  xy_parts <- lapply(strsplit(xy_cols, "_"), function(x) as.numeric(x[2]))
  xy_parts <- do.call(rbind, xy_parts)
  xy_cols[order(xy_parts)]
}

sort_columns <- function(df){
  # Extract and sort the column names (excluding 'sample')
  col_names <- names(df)[-1]
  numeric_cols_sorted <- sort_numeric_columns(col_names)
  x_cols_sorted <- sort_X_Y_columns(col_names, "X")
  y_cols_sorted <- sort_X_Y_columns(col_names, "Y")
  
  # Combine sorted columns and add 'sample' at the beginning
  final_col_order <- c("sample", numeric_cols_sorted, x_cols_sorted, y_cols_sorted)
  
  # Rearrange the columns
  final_df <- df[final_col_order]
  
  return(final_df)
}	

# prep kics  file
prep_kics_location <- function(df, id){
  # change column names in te location df
  colnames(df) <- c("window_chr", "window_start", "window_end", "overlap_chr", "overlap_start", "overlap_end", "sample", "full_overlap_name")
  
  # format kics id df
  id$kics_id <- sprintf("%04d", id$kics_id) # add leading zeros to kics id
  
  # print unmatched
  unmatched <- unique(df[!df$sample %in% id$sample, ]$sample) # samples not in kicds id 
  print("no kics id for these samples")
  print(unmatched)
  
  # merge
  df <- merge(df, id, by = "sample") # add kics sample
  
  # format columns
  df <- df %>%
    select(-sample) %>%                  
    rename(sample = kics_id) %>%         
    select(sample, everything())    
  
  # remove cmmrd samples
  df <- remove_lynch(df) # remove lynch samples
  return(df)
}

# remove lynch samples
remove_lynch <- function(df){
  lynch_ids <- c("0063", "0083", "0120", "0141", "0156", "0171", "0232", "0219")
  df <- df[!df$sample %in% lynch_ids, ]
  return(df)
}

# prep lfs file
prep_lfs_location <- function(df){
  # change column names in te location df
  colnames(df) <- c("window_chr", "window_start", "window_end", "overlap_chr", "overlap_start", "overlap_end", "sample", "full_overlap_name")
  
  # format sample names so find match in clinical
  df$sample <- gsub("_.*", "", df$sample)
  #df$sample <- gsub(".realigned-recalibrated", "", df$sample)
  return(df)
}

# prep clinical file
prep_clinical2 <- function(clinical){
  
  clinical$ID <- ifelse(clinical$cohort == "KiCS", sprintf("%04d", as.numeric(clinical$ID)), clinical$ID)
  clinical <- clinical %>%
    select(ID, cohort, tumor_type, tumor_class, sex, age_at_diagnosis, TP53_status, total_tumor_samples, tumor_site, disease_state, treatment_status) %>%
    rename(sample=ID)
  return(clinical)
}

merge_all <- function(kics_sv, lfs_sv, clinical_df){
  merged <- rbind(kics_sv, lfs_sv) # merge kics and lfs
  all <- merge(merged, clinical_df, by.x="kics_id", by.y="ID", all.x=TRUE) # merge with clinical
  unmatched_merged <- unique(merged[!merged$kics_id %in% clinical_df$ID, ]$kics_id) # samples not in clinical 
  print("no clinical info for these samples")
  print(unmatched_merged)
  return(all)
}

# merge te location and clinical 
merge_df_clinical <- function(te_df, clinical_df){
  unmatched_merged <- unique(te_df[!te_df$sample %in% clinical_df$sample, ]$sample) # samples not in clinical 
  print("no clinical info for these samples")
  print(unmatched_merged)
  
  # merge with clinical
  merged_clinical <- merge(te_df, clinical_df, by="sample") # merge 
  
  return(merged_clinical)
}

scale_keep_order <- function(df){
  # Copy the dataframe to keep the original order of columns
  df_scaled <- df 
  
  # Scale only the numeric columns, excluding 'sample'
  numeric_cols <- which(sapply(df, is.numeric))
  df_scaled[numeric_cols] <- scale(df[numeric_cols])
  
  return(df_scaled) 
}

# Function to identify columns to scale
get_cols_to_scale <- function(df) {
  grep("^[^_]+_\\d+_\\d+$", names(df), value = TRUE)
}

# Function to scale columns
scale_columns <- function(df, cols) {
  scale_factors <- lapply(df[, cols, drop = FALSE], function(x) {
    list(mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE))
  })
  
  df[, cols] <- Map(function(x, y) (x - y$mean) / y$sd, df[, cols, drop = FALSE], scale_factors)
  list(df = df, scale_factors = scale_factors)
}

# Function to apply scaling factors to another dataframe
apply_scaling_to_other_df <- function(df, cols, scale_factors) {
  df[, cols] <- Map(function(x, y) (x - y$mean) / y$sd, df[, cols, drop = FALSE], scale_factors)
  df
}

filter_df <- function(df, filter){
  if (filter=="aff"){ # only affected with cancer
    df_filt<- df %>%
      filter(tumor_type!="U") # only samples affected with cancer
  } else if (filter=="lfs"){ # only lfs
    df_filt<- df %>%  
      filter(TP53_status=="Yes") %>%
      mutate(Cancer = ifelse(tumor_type != "U", "Affected", "Unaffected")) # afftected or unaffected
  } else{ # all samples
    df_filt <- df
  }
  df_filt <- df_filt %>% filter(as.numeric(age_at_diagnosis)/365<30)
  return(df_filt)
}

filter_count_clin_lfs <- function(df){
  return(df %>%  
           filter(TP53_status=="Yes") %>%
           mutate(Cancer = ifelse(tumor_type != "U", "Affected", "Unaffected")))
}

get_samples_below_age <- function(df, n) {
  # Ensure that age_at_diagnosis is a numeric column
  df$age_at_diagnosis <- as.numeric(df$age_at_diagnosis)
  
  # Filter the dataframe where age_at_diagnosis is less than n
  samples_above_age <- df[df$age_at_diagnosis > n, "sample"]
  
  return(samples_above_age)
}

remove_samples <- function(df, samples_to_remove) {
  # Ensure that the 'sample' column is of the same type as samples_to_remove
  df$sample <- as.character(df$sample)
  samples_to_remove <- as.character(samples_to_remove)
  
  # Filter out the rows where the 'sample' column matches any value in samples_to_remove
  filtered_df <- df[!(df$sample %in% samples_to_remove), ]
  
  return(filtered_df)
}

pca_location <- function(df_scaled, colour, filter, samples_to_exclude){
  # pca dimensions
  pca_result <- prcomp(df_scaled[, -which(names(df_scaled) == "sample")], center = FALSE, scale. = FALSE)
  
  # make df for graphing 
  pca_df <- data.frame(PC1 = pca_result$x[, 1],
                       PC2 = pca_result$x[, 2],
                       sample = df_scaled$sample)
  
  # merge with clinical 
  pca_df <- merge(pca_df, clinical, by="sample")
  
  # filter for affected, lfs or NA 
  pca_df <- filter_df(pca_df, filter)
  
  # remove PCA outlier
  if (length(samples_to_exclude) > 0) {
    pca_df <- pca_df %>% 
      filter(!sample %in% samples_to_exclude)
  }
  
  # Plotting the first two principal components
  ggplot(pca_df, aes(x = PC1, y = PC2, color = {{ colour }})) +
    geom_point(alpha=0.8, size=5) +  # Scatter plot
    geom_text(aes(label=sample), vjust=2, hjust=-2, size=3, color="blue") +  # Labels
    labs(x = "Principal Component 1",
         y = "Principal Component 2",
         title = "PCA of Scaled Data") 
}

umap_location <- function(df_scaled, colour, filter, n_neighbors = 15, min_dist = 0.1, n_components = 2, samples_to_exclude, add_labels=FALSE){
  # Perform UMAP
  umap_result <- umap::umap(df_scaled[, -which(names(df_scaled) == "sample")],
                            n_neighbors = n_neighbors, 
                            min_dist = min_dist, 
                            n_components = n_components)
  
  # Create a dataframe for plotting
  umap_df <- data.frame(UMAP1 = umap_result$layout[,1],
                        UMAP2 = umap_result$layout[,2],
                        sample = df_scaled$sample)
  
  # merge with clinical
  umap_df <- merge(umap_df, clinical, by="sample")
  
  # filter for affected, lfs or NA 
  umap_df <- filter_df(umap_df, filter)
  
  # Remove rows with samples to exclude if provided
  if (length(samples_to_exclude) > 0) {
    umap_df <- umap_df %>%
      filter(!sample %in% samples_to_exclude)
  }
  
  # Start the ggplot
  p <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = {{ colour }})) +
    geom_point()
  
  # Add sample labels if add_labels is TRUE
  if (add_labels) {
    p <- p + geom_text(aes(label = sample), vjust = 2, hjust = 2, size = 3, color = "blue")
  }
  
  # Finalize the plot
  p <- p + labs(x = "UMAP Dimension 1",
                y = "UMAP Dimension 2")
  
  return(p)
}

create_heatmap <- function(data, chr, scales, cluster_rows = TRUE, cluster_cols = TRUE, cluster_dist_cols, cluster_dist_rows, cluster_method, annotation_cols = NULL) {
  
  # Filter chromosome
  if (!is.na(chr)){
    pattern <- paste0("^", chr, "_")
    selected_columns <- grep(pattern, names(data), value = TRUE)
    all_selected_columns <- union(selected_columns, annotation_cols)
    data <- data[, all_selected_columns]
  } 
  
  # separate data and annotations into separate dataframes 
  annotations <- data[, annotation_cols, drop = FALSE]
  data <- data[, !(colnames(data) %in% annotation_cols)]
  
  # make age continuous
  if ("age_at_diagnosis" %in% names(annotations)){
    annotations$age_at_diagnosis <- as.numeric(annotations$age_at_diagnosis)
  }
  if ("age_diff" %in% names(annotations)){
    annotations$age_diff <- as.numeric(annotations$age_diff)
  }
  
  # set rownames
  rownames(data) <- data$sample
  rownames(annotations) <- rownames(data)
  
  # Remove columns that do not fit the format 'text_digit_digit' because other clinical data that wasnt used for annotations
  pattern <- "^[^_]+_\\d+_\\d+$" 
  data <- data.frame(data[, grepl(pattern, colnames(data))])
  
  # set colour palette
  # Define the color palette
  color_palette <- colorRampPalette(c("#5FBFF9", "white", "#931621"))(100)
  
  # Calculate the max absolute value to ensure symmetry around 0
  max_abs_value <- max(abs(c(min(data), max(data))))
  
  # Create symmetric breaks around 0
  breaks <- seq(-max_abs_value, max_abs_value, length.out = length(color_palette) + 1)
  
  # Create the heatmap
  pheatmap_result <- pheatmap::pheatmap(data, 
                                        scale = scales, 
                                        cluster_rows = cluster_rows, 
                                        cluster_cols = cluster_cols, 
                                        clustering_distance_rows = cluster_dist_rows,
                                        clustering_distance_cols = cluster_dist_cols,  
                                        clustering_method = cluster_method,  
                                        annotation_row = annotations,
                                        show_rownames = FALSE,
                                        show_colnames = FALSE,
                                        color = color_palette, 
                                        breaks = breaks)
  
  return(pheatmap_result)
}

filter_te_count <- function(df, n_te, filter){
  # filter to at least 5 TEs in window
  numeric_columns <-  df%>%
    select(where(is.numeric)) %>% 
    select_if(~sum(.) > n_te)
  
  # Adding specific non-numeric columns
  other_columns <- c("sample", "TP53_status", "tumor_type", "tumor_class", "age_at_diagnosis", "sex", "cluster")
  
  # Combine numeric columns with the non-numeric columns
  df_filt<- df%>% 
    select(all_of(other_columns), all_of(names(numeric_columns)))
  
  # filter for affected or lfs	
  df_filt<- filter_df(df_filt, filter) # filter can be "aff" or "lfs" of "none"
  
  return(df_filt)
}

# wilcox test for location bins
wilcox_location <- function(df, n_te, chr, group){
  # filter for chromosome
  if (!is.na(chr)){
    df <- df %>% select(matches(paste0("^", chr, "_")))  # select specific chromosome 
  }
  
  wilcox_results <- df %>%
    select(-cluster) %>%
    select(where(~ is.numeric(.) && sum(.) >= n_te)) %>% # filter for numeric columns with at least n TE in window
    names() %>%
    
    # Perform Wilcoxon test and calculate medians for each column
    map_df(~{
      column_name <- .x
      test_result <- wilcox.test(df[[column_name]] ~ df[[group]])
      
      # Get unique levels of the group variable
      group_levels <- levels(factor(df[[group]]))
      
      # Calculate medians for each group
      group1_median <- median(df[[column_name]][df[[group]] == group_levels[1]], na.rm = TRUE)
      group2_median <- median(df[[column_name]][df[[group]] == group_levels[2]], na.rm = TRUE)
      
      # make final tibble
      tibble(
        window = column_name,
        group1_median = group1_median,
        group2_median = group2_median,
        p_val = test_result$p.value
      )
    }) %>%
    # Adjust the p-values using BH method
    mutate(fdr = p.adjust(p_val, method = "BH"))
  
  return(wilcox_results)
}

# identify significant windows
significant_windows <- function(df){
  results <- df %>%
    filter(fdr < 0.05) %>%
    pull(window)
  return(results)
}

# summary of significant columns
significant_windows_summary <- function(df, sig_window){
  for (column in sig_window){
    print(column)
    print(summary(df %>%  select(all_of(column)) %>% unlist))
  }
}

# plot significant windows
significant_windows_plot <- function(df, sig_window, group, xlab){
  for (column in sig_window){
    p <- ggplot(df, aes(x = {{ group }}, y = .data[[column]], fill = {{ group }})) +
      geom_boxplot(outlier.shape = NA) +
      scale_fill_manual(values =  colours) +
      geom_jitter(width = 0.4, height=0, color = "black", size = 1.5) +
      labs(title = paste("TE counts for window", column),
           x = xlab,
           y = "TE count in window")
    print(p)
  }
}

# top windows from results
get_top_results <- function(sig_window, n_windows){
  top_results <- sig_window %>%
    arrange(fdr, p_val) %>%
    slice_head(n = n_windows)
  return(top_results)
}

extract_columns <- function(df, top_results, annotations){
  # filter cleaned data for top results columns
  # Identify matching column names
  matching_col_names <- top_results$window[which(top_results$window %in% names(df))]
  
  #  column names to keep
  col_names_to_keep <- c("sample", annotations, matching_col_names)
  
  # Subset cleaned_data to get only the desired columns
  extracted_data <- select(df, all_of(col_names_to_keep))
  
  return(extracted_data)
}

# heatmap with top windows
significant_windows_prep_heatmap <- function(df, sig_window, n_windows, annotations){
  # get top n tes
  top_results <- get_top_results(sig_window, n_windows)
  
  # get top results columns from dataframe and annotations
  final_data <- extract_columns(df, top_results, annotations)
  
  return(final_data)
}

check_constant_columns <- function(data) {
  constant_columns <- sapply(data, function(column) {
    length(unique(column)) <= 1
  })
  constant_column_names <- names(constant_columns[constant_columns])
  
  if (length(constant_column_names) > 0) {
    cat("Constant or near-constant columns found:\n")
    print(constant_column_names)
  } else {
    cat("No constant or near-constant columns found.\n")
  }
}

remove_constant_columns <- function(data) {
  # Identify constant or near-constant columns
  constant_columns <- sapply(data, function(column) {
    length(unique(column)) <= 1
  })
  
  # Columns to keep (non-constant columns)
  columns_to_keep <- names(constant_columns[!constant_columns])
  
  # Subset the data to keep only non-constant columns
  cleaned_data <- data[, columns_to_keep]
  
  return(cleaned_data)
}

split_data <- function(df, stratify_by, train_size = 0.8, set_seed = 123) {
  # Set seed for reproducibility
  set.seed(set_seed)
  
  # Initialize vectors for train and test indices
  train_indices <- c()
  test_indices <- c()
  
  # Get the unique levels of the stratify_by column
  levels_stratify_by <- levels(factor(df[[stratify_by]]))
  
  # Loop through each level and perform stratified sampling
  for (level in levels_stratify_by) {
    # Get indices of rows for the current level
    indices <- which(df[[stratify_by]] == level)
    
    # Determine the number of samples for the training set
    train_size_level <- round(length(indices) * train_size)
    
    # Sample indices for the training set
    train_indices_level <- sample(indices, train_size_level)
    
    # The rest of the indices form the test set
    test_indices_level <- setdiff(indices, train_indices_level)
    
    # Append indices
    train_indices <- c(train_indices, train_indices_level)
    test_indices <- c(test_indices, test_indices_level)
  }
  
  # Subset the original dataframe to create the training and test sets
  train_set <- df[train_indices, , drop = FALSE]
  test_set <- df[test_indices, , drop = FALSE]
  
  # Return the training and test sets
  list(train = train_set, test = test_set)
}

# filter out columns where n% of rows are not 0
filter_columns <- function(df, n_percent) {
  threshold <- n_percent / 100 * nrow(df)
  cols_to_keep <- sapply(df[, !names(df) %in% "sample"], function(col) sum(col != 0) >= threshold)
  df[, c("sample", names(df)[!names(df) %in% "sample"][cols_to_keep])]
}

# count number of columns where n% of rows are not 0
count_columns <- function(df, n_percent) {
  threshold <- n_percent / 100 * nrow(df)
  count <- sum(sapply(df, function(col) sum(col != 0) >= threshold))
  cat("Number of columns with at least", n_percent, "% of rows not 0:", count, "\n")
}

# Function to plot histogram for age_at_diagnosis
plot_histogram_age <- function(data, bin) {
  # convert age to years
  data$age_at_diagnosis <- data$age_at_diagnosis / 365.25
  
  # Plot
  ggplot(data, aes(x = age_at_diagnosis)) +
    geom_histogram(binwidth = bin, fill = "#AB1368", color = "black") +
    labs(x = "Age at Diagnosis (years)", y = "Frequency")
}

# Function to plot a pie chart for a given column with 'Other' category and labels
plot_pie_chart <- function(data, column_name, ylab) {
  # Calculate frequencies
  data_to_plot <- as.data.frame(table(data[[column_name]]))
  names(data_to_plot) <- c("Category", "Freq")

  # Calculate total
  total <- sum(data_to_plot$Freq)

  # Determine threshold for grouping into 'Other'
  threshold <- 0.02 * total

  # Group small categories into 'Other'
  data_to_plot$Category <- ifelse(data_to_plot$Freq < threshold, 'Other', as.character(data_to_plot$Category))

  # Aggregate frequencies by category
  data_to_plot <- aggregate(Freq ~ Category, data_to_plot, sum)

  if (column_name == "TP53_status") {
    data_to_plot$Category <- factor(data_to_plot$Category, levels = c("WT", "Mutant"))
  }

  if (column_name == "sex") {
    category_mapping <- c("F" = "Female", "M" = "Male")
    # Map abbreviations to full names
    data_to_plot$Category <- category_mapping[as.character(data_to_plot$Category)]
    # Factorize with desired levels
    data_to_plot$Category <- factor(data_to_plot$Category, levels = c("Female", "Male"))
  }

  if (column_name == "tumor_type") {
    data_to_plot <- data_to_plot %>% dplyr::filter(Category != "U") # remove no cancer
  }

  # Calculate percentage for labels
  data_to_plot$Percentage <- round((data_to_plot$Freq / total) * 100, 1)
  data_to_plot$Label <- paste0(data_to_plot$Category, ": ", data_to_plot$Percentage, "%")

  # Dynamically generate a color palette for the number of categories
  n_cats <- nrow(data_to_plot)
  # Use multiple qualitative palettes for more variety
  palette_list <- c("Set3", "Paired", "Dark2", "Pastel1", "Pastel2", "Set1", "Set2", "Accent")
  all_colours <- unlist(lapply(palette_list, function(pal) RColorBrewer::brewer.pal(RColorBrewer::brewer.pal.info[pal,]$maxcolors, pal)))
  # Remove duplicates and ensure enough colors
  all_colours <- unique(all_colours)
  if (n_cats > length(all_colours)) {
    pie_colours <- colorRampPalette(all_colours)(n_cats)
  } else {
    pie_colours <- all_colours[1:n_cats]
  }

  # Plot
  p <- ggplot(data_to_plot, aes(x = "", y = Freq, fill = Category)) +
    geom_bar(width = 1, stat = "identity") +
    coord_polar("y", start = 0) +
    theme_void() +
    theme(legend.position = "right")

  # Always use the dynamically generated palette
  color_mapping <- setNames(pie_colours, data_to_plot$Category)
  p <- p + scale_fill_manual(values = color_mapping, name = ylab, labels = data_to_plot$Label)

  # Adjust legend title and text size
  p <- p + theme(legend.title = element_text(size = 18),
                 legend.text = element_text(size = 16))

  return(p)
}

make_ml_matrix <- function(df) {
  # Select columns based on the pattern and include 'sample' and 'TP53_status'
  cols_to_keep <- c('sample', 'TP53_status', grep("^[^_]+_\\d+_\\d+$", names(df), value = TRUE))
  df <- df[, cols_to_keep]
  
  # Rename 'TP53_status' to 'label'
  names(df)[names(df) == 'TP53_status'] <- 'label'
  
  # Recode 'label': 'LFS' to 1 and 'Control' to 0
  df$label <- ifelse(df$label == 'LFS', 1, ifelse(df$label == 'Control', 0, df$label))
  
  # Return the processed dataframe
  return(df)
}

identify_sv_in_te <- function(df_sv, df_te) {
  setDT(df_sv)
  setDT(df_te)
  
  # Rename columns in df_te for clarity
  setnames(df_te, old = c("SV_start", "SV_end", "SV.chrom", "SV.type"), new = c("TE.start", "TE.end", "TE.chrom", "TE.type"))
  
  results_list <- list()
  samples <- unique(df_sv$sample)
  
  for (sample_id in samples) {
    sv_sample <- df_sv[sample == sample_id]
    te_sample <- df_te[sample == sample_id]
    
    # Find SVs with start positions within TEs
    start_matches <- sv_sample[te_sample, on = .(SV.chrom = TE.chrom), nomatch = 0, allow.cartesian = TRUE]
    start_matches <- start_matches[TE.start <= SV.start & TE.end >= SV.start, 
                                   .(sample, SV.chrom, SV.start, SV.end, SV.type, TE.start, TE.end, TE.type)]
    
    # Find SVs with end positions within TEs
    end_matches <- sv_sample[te_sample, on = .(SV.chrom = TE.chrom), nomatch = 0, allow.cartesian = TRUE]
    end_matches <- end_matches[TE.start <= SV.end & TE.end >= SV.end, 
                               .(sample, SV.chrom, SV.start, SV.end, SV.type, TE.start, TE.end, TE.type)]
    
    # Merge start and end matches on common columns
    combined_matches <- merge(start_matches, end_matches, 
                              by = c("sample", "SV.chrom", "SV.start", "SV.end", "SV.type", "TE.type"),
                              suffixes = c(".one", ".two"), allow.cartesian = TRUE)
    
    # Filter for different TE instances with the same TE type
    final_matches <- combined_matches[TE.start.one!= TE.start.two & TE.end.one != TE.end.two, ]
    
    if (nrow(final_matches) > 0) {
      results_list[[sample_id]] <- final_matches
    }
  }
  
  # Combine all results into a single data.table
  results <- rbindlist(results_list, use.names = TRUE, fill = TRUE)
  
  return(results)
}

identify_sv_in_te_onebreak <- function(df_sv, df_te) {
  setDT(df_sv)
  setDT(df_te)
  
  # Rename columns in df_te for clarity
  setnames(df_te, old = c("SV_start", "SV_end", "SV.chrom", "SV.type"), new = c("TE.start", "TE.end", "TE.chrom", "TE.type"))
  
  results_list <- list()
  samples <- unique(df_sv$sample)
  
  for (sample_id in samples) {
    sv_sample <- df_sv[sample == sample_id]
    te_sample <- df_te[sample == sample_id]
    
    # Find SVs with start positions within TEs
    start_matches <- sv_sample[te_sample, on = .(SV.chrom = TE.chrom), nomatch = 0, allow.cartesian = TRUE]
    start_matches <- start_matches[TE.start <= SV.start & TE.end >= SV.start, 
                                   .(sample, SV.chrom, SV.start, SV.end, SV.type, TE.start, TE.end, TE.type)]
    
    # Find SVs with end positions within TEs
    end_matches <- sv_sample[te_sample, on = .(SV.chrom = TE.chrom), nomatch = 0, allow.cartesian = TRUE]
    end_matches <- end_matches[TE.start <= SV.end & TE.end >= SV.end, 
                               .(sample, SV.chrom, SV.start, SV.end, SV.type, TE.start, TE.end, TE.type)]
    
    # Combine start and end matches for the current sample, removing duplicates
    combined_matches <- unique(rbindlist(list(start_matches, end_matches), use.names = TRUE, fill = TRUE))
    
    if (nrow(combined_matches) > 0) {
      results_list[[sample_id]] <- combined_matches
    }
  }
  
  # Combine all results from all samples into a single data.table
  results <- rbindlist(results_list, use.names = TRUE, fill = TRUE)
  
  return(results)
}

identify_sv_in_te_linx <- function(df_sv, df_te) {
  setDT(df_sv)
  setDT(df_te)

  # Rename columns in df_te for clarity
  setnames(df_te, old = c("SV_start", "SV_chrom", "SV_length", "ALT"), new = c("TE.start", "TE.chrom", "TE.length", "TE.type"))
  setnames(df_sv, old = c("SampleId", "PosStart", "PosEnd", "ChrStart", "Type"), new = c("sample", "SV.start", "SV.end", "SV.chrom", "SV.type"))

  # end of TE is start plus length
  df_te$TE.end = df_te$TE.start + df_te$TE.length

  results_list <- list()
  samples_sv <- unique(df_sv$sample)

  # Get TE samples that have at least one TE (count > 0)
  # df_te has been filtered so each row is a TE insertion, meaning all samples here have count > 0
  samples_te_with_count <- unique(df_te$sample)
  samples_te_base <- sub("_T$", "", samples_te_with_count)  # Remove _T suffix

  # Find overlap - also check with underscore to hyphen conversion
  samples_with_both <- intersect(samples_sv, samples_te_base)

  # Also check if SV samples with _ converted to - match TE samples
  samples_sv_alt <- gsub("_", "-", samples_sv)
  samples_with_both_alt <- intersect(samples_sv_alt, samples_te_base)

  # Combine matches
  all_matched_te_samples <- unique(c(samples_te_base[samples_te_base %in% samples_sv],
                                      samples_te_base[samples_te_base %in% samples_sv_alt]))
  samples_te_no_sv <- setdiff(samples_te_base, all_matched_te_samples)

  # Hardcoded sample mappings (SV sample -> TE sample)
  hardcoded_mappings <- c(
    "U02H2D_5524A" = "U02H2D_A_T",
    "5009_4856_1" = "5009_1_T",
    "5471_5787_2" = "5471_2_T",
    "621_3311A_1" = "621_1_T",
    "PD13489_PD13489a" = "PD13489_T"
  )

  for (sample_id in samples_sv) {
    sv_sample <- df_sv[sample == sample_id]

    # Strategy 0: Check hardcoded mappings first
    if (sample_id %in% names(hardcoded_mappings)) {
      te_sample_id <- hardcoded_mappings[[sample_id]]
      te_sample <- df_te[sample == te_sample_id]
      if (nrow(te_sample) > 0) {
        cat("  Note: Matched SV sample", sample_id, "to TE sample", te_sample_id, "(hardcoded mapping)\n")
      }
    } else {
      # Strategy 1: Try exact match with _T suffix
      te_sample_id <- paste0(sample_id, "_T")
      te_sample <- df_te[sample == te_sample_id]
    }

    # Strategy 2: Try converting SV all underscores to hyphens (e.g., 0218_18_751 -> 0218-18-751_T)
    if (nrow(te_sample) == 0) {
      te_sample_id_alt <- paste0(gsub("_", "-", sample_id), "_T")
      te_sample <- df_te[sample == te_sample_id_alt]
      if (nrow(te_sample) > 0) {
        cat("  Note: Matched SV sample", sample_id, "to TE sample", te_sample_id_alt, "(SV _ to -)\n")
      }
    }

    # Strategy 3: Try keeping first underscore, convert rest to hyphens (e.g., 0296_19_6855 -> 0296_19-6855_T)
    if (nrow(te_sample) == 0) {
      # Split on first underscore, then convert remaining underscores to hyphens
      parts <- strsplit(sample_id, "_", fixed = TRUE)[[1]]
      if (length(parts) > 2) {
        te_sample_id_alt3 <- paste0(parts[1], "_", paste(parts[-1], collapse = "-"), "_T")
        te_sample <- df_te[sample == te_sample_id_alt3]
        if (nrow(te_sample) > 0) {
          cat("  Note: Matched SV sample", sample_id, "to TE sample", te_sample_id_alt3, "(SV keep first _, rest to -)\n")
        }
      }
    }

    # Strategy 4: Try matching with flexible -/_ replacement
    if (nrow(te_sample) == 0) {
      # Create pattern by replacing each _ or - with a regex that matches either
      sample_pattern <- gsub("[-_]", "[-_]", sample_id)
      matching_samples <- grep(paste0("^", sample_pattern, "_T$"), df_te$sample, value = TRUE)
      if (length(matching_samples) > 0) {
        te_sample <- df_te[sample == matching_samples[1]]
        if (nrow(te_sample) > 0) {
          cat("  Note: Matched SV sample", sample_id, "to TE sample", matching_samples[1], "(flexible -/_)\n")
        }
      }
    }

    # Skip if still no TE data for this sample
    if (nrow(te_sample) == 0) {
      next
    }

    # Find SVs with start positions within TEs
    start_matches <- sv_sample[te_sample, on = .(SV.chrom = TE.chrom), nomatch = 0, allow.cartesian = TRUE]
    start_matches <- start_matches[TE.start <= SV.start & TE.end >= SV.start, 
                                   .(sample, SV.chrom, SV.start, SV.end, SV.type, TE.start, TE.end, TE.type)]
    
    # Find SVs with end positions within TEs
    end_matches <- sv_sample[te_sample, on = .(SV.chrom = TE.chrom), nomatch = 0, allow.cartesian = TRUE]
    end_matches <- end_matches[TE.start <= SV.end & TE.end >= SV.end, 
                               .(sample, SV.chrom, SV.start, SV.end, SV.type, TE.start, TE.end, TE.type)]
    
    # Merge start and end matches on common columns
    combined_matches <- merge(start_matches, end_matches, 
                              by = c("sample", "SV.chrom", "SV.start", "SV.end", "SV.type", "TE.type"),
                              suffixes = c(".one", ".two"), allow.cartesian = TRUE)
    
    # Filter for different TE instances with the same TE type
    final_matches <- combined_matches[TE.start.one!= TE.start.two & TE.end.one != TE.end.two, ]
    
    if (nrow(final_matches) > 0) {
      results_list[[sample_id]] <- final_matches
    }
  }
  
  # Combine all results into a single data.table
  results <- rbindlist(results_list, use.names = TRUE, fill = TRUE)

  cat("  Number of samples with SV data:", length(samples_sv), "\n")
  cat("  Number of samples with TE count > 0:", length(samples_te_base), "\n")
  cat("  Number of samples with both SV and TE data:", length(samples_with_both), "\n")
  cat("  Number of samples with TE count > 0 but no SV data:", length(samples_te_no_sv), "\n")
  cat("  SVs found with both breakpoints in TEs:", nrow(results), "\n")

  # Return both results and sample lists
  attr(results, "samples_te_no_sv") <- samples_te_no_sv
  return(results)
}

identify_sv_in_te_linx_onebreak <- function(df_sv, df_te) {
  setDT(df_sv)
  setDT(df_te)

  # Rename columns in df_te for clarity
  setnames(df_te, old = c("SV_start", "SV_chrom", "SV_length", "ALT"), new = c("TE.start", "TE.chrom", "TE.length", "TE.type"))
  setnames(df_sv, old = c("SampleId", "PosStart", "PosEnd", "ChrStart", "Type"), new = c("sample", "SV.start", "SV.end", "SV.chrom", "SV.type"))

  # end of TE is start plus length
  df_te$TE.end = df_te$TE.start + df_te$TE.length

  results_list <- list()
  samples_sv <- unique(df_sv$sample)

  # Get TE samples that have at least one TE (count > 0)
  # df_te has been filtered so each row is a TE insertion, meaning all samples here have count > 0
  samples_te_with_count <- unique(df_te$sample)
  samples_te_base <- sub("_T$", "", samples_te_with_count)  # Remove _T suffix

  # Find overlap - also check with underscore to hyphen conversion
  samples_with_both <- intersect(samples_sv, samples_te_base)

  # Also check if SV samples with _ converted to - match TE samples
  samples_sv_alt <- gsub("_", "-", samples_sv)
  samples_with_both_alt <- intersect(samples_sv_alt, samples_te_base)

  # Combine matches
  all_matched_te_samples <- unique(c(samples_te_base[samples_te_base %in% samples_sv],
                                      samples_te_base[samples_te_base %in% samples_sv_alt]))
  samples_te_no_sv <- setdiff(samples_te_base, all_matched_te_samples)

  # Hardcoded sample mappings (SV sample -> TE sample)
  hardcoded_mappings <- c(
    "U02H2D_5524A" = "U02H2D_A_T",
    "5009_4856_1" = "5009_1_T",
    "5471_5787_2" = "5471_2_T",
    "621_3311A_1" = "621_1_T",
    "PD13489_PD13489a" = "PD13489_T"
  )

  for (sample_id in samples_sv) {
    sv_sample <- df_sv[sample == sample_id]

    # Strategy 0: Check hardcoded mappings first
    if (sample_id %in% names(hardcoded_mappings)) {
      te_sample_id <- hardcoded_mappings[[sample_id]]
      te_sample <- df_te[sample == te_sample_id]
      if (nrow(te_sample) > 0) {
        cat("  Note: Matched SV sample", sample_id, "to TE sample", te_sample_id, "(hardcoded mapping)\n")
      }
    } else {
      # Strategy 1: Try exact match with _T suffix
      te_sample_id <- paste0(sample_id, "_T")
      te_sample <- df_te[sample == te_sample_id]
    }

    # Strategy 2: Try converting SV all underscores to hyphens (e.g., 0218_18_751 -> 0218-18-751_T)
    if (nrow(te_sample) == 0) {
      te_sample_id_alt <- paste0(gsub("_", "-", sample_id), "_T")
      te_sample <- df_te[sample == te_sample_id_alt]
      if (nrow(te_sample) > 0) {
        cat("  Note: Matched SV sample", sample_id, "to TE sample", te_sample_id_alt, "(SV _ to -)\n")
      }
    }

    # Strategy 3: Try keeping first underscore, convert rest to hyphens (e.g., 0296_19_6855 -> 0296_19-6855_T)
    if (nrow(te_sample) == 0) {
      # Split on first underscore, then convert remaining underscores to hyphens
      parts <- strsplit(sample_id, "_", fixed = TRUE)[[1]]
      if (length(parts) > 2) {
        te_sample_id_alt3 <- paste0(parts[1], "_", paste(parts[-1], collapse = "-"), "_T")
        te_sample <- df_te[sample == te_sample_id_alt3]
        if (nrow(te_sample) > 0) {
          cat("  Note: Matched SV sample", sample_id, "to TE sample", te_sample_id_alt3, "(SV keep first _, rest to -)\n")
        }
      }
    }

    # Strategy 4: Try matching with flexible -/_ replacement
    if (nrow(te_sample) == 0) {
      # Create pattern by replacing each _ or - with a regex that matches either
      sample_pattern <- gsub("[-_]", "[-_]", sample_id)
      matching_samples <- grep(paste0("^", sample_pattern, "_T$"), df_te$sample, value = TRUE)
      if (length(matching_samples) > 0) {
        te_sample <- df_te[sample == matching_samples[1]]
        if (nrow(te_sample) > 0) {
          cat("  Note: Matched SV sample", sample_id, "to TE sample", matching_samples[1], "(flexible -/_)\n")
        }
      }
    }

    # Skip if still no TE data for this sample
    if (nrow(te_sample) == 0) {
      next
    }

    # Find SVs with start positions within TEs
    start_matches <- sv_sample[te_sample, on = .(SV.chrom = TE.chrom), nomatch = 0, allow.cartesian = TRUE]
    start_matches <- start_matches[TE.start <= SV.start & TE.end >= SV.start,
                                   .(sample, SV.chrom, SV.start, SV.end, SV.type, TE.start, TE.end, TE.type)]

    # Find SVs with end positions within TEs
    end_matches <- sv_sample[te_sample, on = .(SV.chrom = TE.chrom), nomatch = 0, allow.cartesian = TRUE]
    end_matches <- end_matches[TE.start <= SV.end & TE.end >= SV.end,
                               .(sample, SV.chrom, SV.start, SV.end, SV.type, TE.start, TE.end, TE.type)]

    # Combine start and end matches for the current sample, removing duplicates
    combined_matches <- unique(rbindlist(list(start_matches, end_matches), use.names = TRUE, fill = TRUE))

    if (nrow(combined_matches) > 0) {
      results_list[[sample_id]] <- combined_matches
    }
  }

  # Combine all results from all samples into a single data.table
  results <- rbindlist(results_list, use.names = TRUE, fill = TRUE)

  cat("  Number of samples with SV data:", length(samples_sv), "\n")
  cat("  Number of samples with TE count > 0:", length(samples_te_base), "\n")
  cat("  Number of samples with both SV and TE data:", length(samples_with_both), "\n")
  cat("  Number of samples with TE count > 0 but no SV data:", length(samples_te_no_sv), "\n")
  cat("  SVs found with at least one breakpoint in TE:", nrow(results), "\n")
  return(results)
}
  
extract_location <- function(df, mutation_column) {
  df %>%
    mutate(location = ifelse(str_detect(.data[[mutation_column]], "\\s"),
                             NA_integer_,  # Assign NA of type integer if there's a space
                             as.integer(str_extract(.data[[mutation_column]], "\\d+"))))  # Extract digits otherwise
}

# plot TE frequency by location
plot_te_locations <- function(df, chr, type, hotspots, log_scale=FALSE){
  # Use the helper function to construct the combination label
  combination <- construct_combination_label(chr, type)
  
  # Basic plot setup
  p <- ggplot(df, aes(x = protein.codon.num, y = !!sym(combination))) +
    geom_point(alpha=0.7, size=4, color = "#5FBFF9") +  # Apply custom colors
    scale_x_continuous(limits = c(0, 350)) +
    labs(x = "p53 protein residues", y = "TE count") 
  
  # Add arrows at hotspots, pointing downwards
  if (!is.null(hotspots) && length(hotspots) > 0) {
    p <- p + geom_segment(data = data.frame(x = hotspots), 
                          aes(x = x, xend = x, y = -0.5, yend = -1),  # Adjust y positions for arrows
                          arrow = arrow(type = "closed", length = unit(0.1, "inches"), ends = "first"),  # "ends = first" makes the arrow point downwards
                          colour = "black", inherit.aes = FALSE)
  }
  
  # Apply logarithmic scale if log_scale is TRUE
  if (log_scale) {
    p <- p + scale_y_continuous(
      trans = scales::log1p_trans()
    )
  }
  
  return(p)
}

# make shelly te count df
make_count_df <- function(te_df, filter_var="SV.type", filter_element, chr=NA, name){
  if (!is.na(chr)){
    te_df <- te_df %>% filter(SV.chrom == chr)
  }
  
  # Get all unique samples
  all_samples <- unique(te_df$sample)
  
  if (!is.na(filter_element)){ # if need to filter by TE type
    data <- te_df %>%
      filter(!!sym(filter_var) == {{ filter_element}} ) %>% # filter for deletions
      group_by(sample) %>% # group
      summarise(count = n(), .groups = 'drop') # count
  }
  else { # if doing overall TE count
    data <- te_df %>%
      group_by(sample) %>% # group
      summarise(count = n(), .groups = 'drop') # count by group
  }
  
  # Include all samples, setting count to 0 for samples not in the data
  data <- data %>%
    right_join(tibble(sample = all_samples), by = "sample") %>%
    replace_na(list(count = 0))
  
  colnames(data)[2] <- name
  return(data)
}

# group by gene with start and end within thrshold
group_sv <- function(data, threshold = 50) {
  data <- data %>% arrange(SV_start, SV_end)
  group <- numeric(nrow(data))
  group_id <- 1
  
  for (i in 1:nrow(data)) {
    if (group[i] == 0) {
      group[i] <- group_id
      for (j in (i+1):nrow(data)) {
        if (group[j] == 0 &&
            !is.na(data$SV_start[i]) && !is.na(data$SV_start[j]) &&
            abs(data$SV_start[i] - data$SV_start[j]) <= threshold &&
            !is.na(data$SV_end[i]) && !is.na(data$SV_end[j]) &&
            abs(data$SV_end[i] - data$SV_end[j]) <= threshold) {
          group[j] <- group_id
        }
      }
      group_id <- group_id + 1
    }
  }
  
  data$group <- group
  return(data)
}

# tes grouped by gene and start and end within 50bp. output summary df and df for each group
summarize_sv_groups <- function(data, group_by_cols, threshold = 50) {
  grouped_data <- data %>%
    group_by(across(all_of(group_by_cols))) %>%
    group_modify(~group_sv(.x, threshold)) %>%
    ungroup()
  
  summary_data <- grouped_data %>%
    group_by(across(all_of(group_by_cols)), group) %>%
    summarize(
      SV_start = min(SV_start),
      SV_end = max(SV_end),
      count = n(),
      .groups = 'drop'
    )
  
  result_list <- lapply(1:nrow(summary_data), function(i) {
    row <- summary_data[i, ]
    details <- grouped_data %>%
      filter(
        !!sym(group_by_cols[1]) == row[[group_by_cols[1]]] &
          !!sym(group_by_cols[2]) == row[[group_by_cols[2]]] &
          !!sym(group_by_cols[3]) == row[[group_by_cols[3]]] &
          group == row$group
      )
    list(summary = row, details = details)
  })
  
  return(list(summary = summary_data, details = result_list))
}

group_variants <- function(df, threshold = 50) {
  # Initial sorting
  df <- df %>%
    arrange(SV.chrom, SV.type, SV_start, SV_end)
  
  # Create initial groups based on proximity
  df <- df %>%
    group_by(SV.chrom, SV.type) %>%
    mutate(
      group_number = cumsum(
        (abs(SV_start - lag(SV_start, default = SV_start[1])) > threshold) |
          (abs(SV_end - lag(SV_end, default = SV_end[1])) > threshold)
      ) + 1,
      group = paste(SV.chrom, SV.type, group_number, sep = "_")
    ) %>%
    ungroup()
  
  # Function to check if all elements in a group are within the threshold
  check_group <- function(group_df, threshold) {
    starts_within_threshold <- all(as.matrix(dist(group_df$SV_start)) <= threshold)
    ends_within_threshold <- all(as.matrix(dist(group_df$SV_end)) <= threshold)
    return(starts_within_threshold & ends_within_threshold)
  }
  
  # Apply the check and create a consistency dataframe
  consistency_df <- df %>%
    group_by(SV.chrom, SV.type, group) %>%
    summarize(consistent = check_group(cur_data(), threshold), .groups = 'drop')
  
  # Merge the consistency information back to the original dataframe
  df <- df %>%
    left_join(consistency_df, by = c("SV.chrom", "SV.type", "group"))
  
  # Select the specified columns
  df <- df %>%
    select(sample, chr = SV.chrom, start = SV_start, end = SV_end, SV.type = SV.type, group, consistent)
  
  return(df)
}

plot_top_genes <- function(df, column, label_column=NULL, top_n = 10, x_lab) {
  if ("Gene_name" %in% colnames(df)) {
    colnames(df)[colnames(df) == "Gene_name"] <- "geneID"
  }
  
  # Check if required columns exist
  if (!"geneID" %in% colnames(df)) {
    stop("The DataFrame must contain a 'geneID' column.")
  }
  
  # Remove rows with NA in critical columns
  df <- df[complete.cases(df[, c(column, label_column, "geneID")]), ]
  
  # Sort the DataFrame by the specified column in descending order
  df <- df[order(-df[[column]]), ]
  
  # Take the top n genes
  top_genes_df <- head(df, n = top_n)
  
  # Check if top_genes_df is non-empty
  if (nrow(top_genes_df) == 0) {
    stop("No data available for plotting after filtering.")
  }
  
  plot <- ggplot(top_genes_df, aes(x = .data[[column]], y = reorder(geneID, .data[[column]]))) +
    geom_bar(stat = "identity", fill = "#0080A3") +  # Bar plot
    labs(
      x = x_lab,
      y = "Gene"
    )
  
  # Add labels only if label_column is not NULL
  if (!is.null(label_column)) {
    plot <- plot + 
      geom_text(aes(label = round(.data[[label_column]], 2)), hjust = -0.1, size = 4)
  }
  return(plot)
}

plot_top_cancer_genes_germline <- function(df, top_n = 10) {
  # Create a frequency table for the "Gene_name" column
  gene_freq <- as.data.frame(table(df$Gene_name))
  print(head(gene_freq, n = 10))  # Print the first 10 rows of the frequency table for debugging
  # Rename the columns for clarity
  colnames(gene_freq) <- c("Gene", "Frequency")
  
  # Sort the frequency table in descending order
  gene_freq <- gene_freq[order(-gene_freq$Frequency), ]
  
  # Take the top n genes
  top_genes_df <- head(gene_freq, n = top_n)
  
  # Create the bar plot
  ggplot(top_genes_df, aes(x = Frequency, y = reorder(Gene, Frequency))) +
    geom_bar(stat = "identity", fill = "#0080A3") +  # Bar plot
    labs(
      x = "Frequency",
      y = "Gene",
    ) 
}

plot_te_insertion_heatmap_line <- function(df, gene_name) {
  # Filter data for the specified gene
  gene_data <- df %>%
    filter(Gene_name == gene_name) %>%
    mutate(SV_end = SV_start + SV_length - 1) # Dynamically calculate SV_end
  
  if (nrow(gene_data) == 0) {
    stop(paste("No insertions found for gene:", gene_name))
  }
  
  # Dynamically determine x-axis limits
  x_min <- min(gene_data$SV_start) - 100
  x_max <- max(gene_data$SV_end) + 100
  
  # Create the heatmap plot showing the start-to-end range of insertions
  heatmap_plot <- ggplot(gene_data, aes(y = sample)) +
    geom_tile(aes(x = (SV_start + SV_end) / 2, width = SV_end - SV_start + 1, height = 0.9), fill = "blue", alpha = 0.7) +
    scale_x_continuous(limits = c(x_min, x_max), expand = c(0, 0)) +
    labs(
      title = paste("TE Insertion Heatmap for Gene:", gene_name),
      x = "Genomic Coordinate",
      y = "Sample"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank()
    )
  
  # Print the plot
  print(heatmap_plot)
  
  return(heatmap_plot)
}

plot_te_insertion_binned_heatmap <- function(df, gene_name, bin_size = 1000) {  
  # Filter data for the specified gene
  gene_data <- df %>%
    filter(Gene_name == gene_name) %>%
    mutate(SV_end = SV_start + SV_length - 1)  # Dynamically calculate SV_end
  
  if (nrow(gene_data) == 0) {
    stop(paste("No insertions found for gene:", gene_name))
  }
  
  # Dynamically determine x-axis limits
  x_min <- min(gene_data$SV_start) - 30000
  x_max <- max(gene_data$SV_end) + 30000
  
  # Bin genomic coordinates
  binned_data <- gene_data %>%
    mutate(
      Bin_start = floor(SV_start / bin_size) * bin_size,
      Bin_end = floor(SV_end / bin_size) * bin_size + bin_size - 1
    ) %>%
    dplyr::select(sample, SV_start, SV_end, Bin_start, Bin_end, tumor_type) %>%  # Ensure tumor_type is included
    group_by(sample, Bin_start) %>%
    reframe(
      SV_start = min(SV_start),  # Minimum SV_start within the bin
      SV_end = max(SV_end),      # Maximum SV_end within the bin
      Count = n(),
      tumor_type = tumor_type  # Use the first tumor type in the group
    ) %>%
    arrange(Bin_start)  # Ensure bins are ordered numerically
  
  print(binned_data, n = 50)
  
  # Create the heatmap with tumor type annotation
  heatmap_plot <- ggplot() +
    # Heatmap for insertion counts
    geom_tile(data = binned_data, aes(x = Bin_start, y = sample, fill = Count), color = "white") +
    scale_fill_gradient(low = "lightblue", high = "darkblue", name = "Insertions") +
    scale_x_continuous(
      limits = c(x_min, x_max), 
      expand = c(0, 0), 
      breaks = seq(x_min, x_max, by = 10000),
      labels = function(x) x  # Ensure numeric labels are clean
    ) +
    labs(
      x = "Binned genomic coordinate",
      y = "Sample"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),  # Enable and angle x-axis labels
      axis.ticks.x = element_line(),  # Ensure x-axis ticks are displayed
      panel.grid = element_blank(),
      axis.text.y = element_blank(),  # Remove y-axis text
      axis.ticks.y = element_blank()  # Remove y-axis ticks
    )
  
  print(heatmap_plot)
}

plot_te_insertion_complex_heatmap <- function(df, gene_name, bin_size = 1000, annotations = NULL, lwd=3) {
  # Filter data for the specified gene
  gene_data <- df %>%
    filter(Gene_name == gene_name) %>%
    mutate(SV_end = SV_start + SV_length - 1)  # Dynamically calculate SV_end
  
  if (nrow(gene_data) == 0) {
    stop(paste("No insertions found for gene:", gene_name))
  }
  
  # Dynamically determine x-axis limits
  x_min <- min(gene_data$SV_start) - 20000
  x_max <- max(gene_data$SV_end) + 20000
  
  # Bin genomic coordinates
  binned_data <- gene_data %>%
    mutate(
      Bin_start = floor(SV_start / bin_size) * bin_size,
      Bin_end = floor(SV_end / bin_size) * bin_size + bin_size - 1,
      Bin_label = paste(Bin_start, Bin_end, sep = "-")
    ) %>%
    group_by(sample, Bin_start) %>%
    reframe(
      Count = n()
    ) %>%
    ungroup()
  
  # Handle annotations
  row_annotation <- NULL
  if (!is.null(annotations)) {
    if (is.character(annotations)) {
      # Extract specified columns for annotations
      annotation_data <- gene_data %>%
        dplyr::select(sample, all_of(annotations)) %>%
        distinct() %>%
        column_to_rownames("sample")
      
      # Define annotation colors
      annotation_colors <- list()
      for (col in colnames(annotation_data)) {
        if (is.numeric(annotation_data[[col]])) {
          # Continuous variable
          annotation_colors[[col]] <- circlize::colorRamp2(
            c(min(annotation_data[[col]], na.rm = TRUE), max(annotation_data[[col]], na.rm = TRUE)),
            c("white", "#AB1368")
          )
        } else {
          # Categorical variable
          unique_vals <- unique(annotation_data[[col]])
          annotation_colors[[col]] <- setNames(
            scales::hue_pal()(length(unique_vals)),  # Generate ggplot-like default colors
            unique_vals
          )
        }
      }
      
      # Create row annotations
      row_annotation <- do.call(rowAnnotation, c(annotation_data, list(col = annotation_colors)))
      
    } else {
      stop("The 'annotations' argument must be a character vector of column names or NULL.")
    }
  }
  
  # Create a matrix for the heatmap
  heatmap_data <- binned_data %>%
    pivot_wider(names_from = Bin_start, values_from = Count, values_fill = list(Count = 0)) %>%
    column_to_rownames("sample") %>%
    as.data.frame()
  
  # Ensure the column names (bins) are ordered numerically
  col_order <- as.numeric(colnames(heatmap_data))  # Convert column names to numeric
  heatmap_data <- heatmap_data[, order(col_order)]  # Reorder columns by numeric bin order
  
  # Convert to matrix
  heatmap_matrix <- as.matrix(heatmap_data)
  
  # Determine if the insertion legend should be shown
  max_insertion <- max(heatmap_matrix, na.rm = TRUE)
  show_legend <- ifelse(max_insertion > 1, TRUE, FALSE)
  
  # Create the heatmap
  heatmap <- Heatmap(
    heatmap_matrix,
    name = "Insertions",
    col = colorRamp2(c(0, max_insertion), c("white",  "#0080A3")),
    show_row_names = FALSE,
    show_column_names = TRUE,
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    rect_gp = gpar(col = "white", lwd = lwd),
    left_annotation = row_annotation,  # Add the annotations if provided
    heatmap_legend_param = list(
      title = if (show_legend) "Number of insertions" else NULL
    ),
    show_heatmap_legend = show_legend  # Disable the heatmap legend entirely if needed
  )
  
  # Draw the heatmap
  return(heatmap)
  
  #return(binned_data)
}

plot_te_insertion_complex_heatmap_allbins <- function(df, gene_name, bin_size = 1000, annotations = NULL) {
  # Filter data for the specified gene
  gene_data <- df %>%
    filter(Gene_name == gene_name) %>%
    mutate(SV_end = SV_start + SV_length - 1)  # Dynamically calculate SV_end
  
  if (nrow(gene_data) == 0) {
    stop(paste("No insertions found for gene:", gene_name))
  }
  
  # Bin genomic coordinates
  binned_data <- gene_data %>%
    mutate(
      Bin_start = floor(SV_start / bin_size) * bin_size,
      Bin_end = floor(SV_end / bin_size) * bin_size + bin_size - 1,
      Bin_label = paste(Bin_start, Bin_end, sep = "-")
    ) %>%
    group_by(sample, Bin_start) %>%
    reframe(
      Count = n()
    ) %>%
    ungroup()
  
  # Generate a complete sequence of bins between the first and last filled bins
  min_bin <- min(binned_data$Bin_start)
  max_bin <- max(binned_data$Bin_start)
  all_bins <- seq(min_bin, max_bin, by = bin_size)
  
  # Expand the dataset to include all bins for each sample
  expanded_data <- expand.grid(
    sample = unique(binned_data$sample),
    Bin_start = all_bins
  ) %>%
    left_join(binned_data, by = c("sample", "Bin_start")) %>%
    mutate(Count = replace_na(Count, 0))  # Fill missing counts with 0
  
  # Handle annotations
  row_annotation <- NULL
  if (!is.null(annotations)) {
    if (is.character(annotations)) {
      # Extract specified columns for annotations
      annotation_data <- gene_data %>%
        dplyr::select(sample, all_of(annotations)) %>%
        distinct() %>%
        column_to_rownames("sample")
      
      # Define annotation colors
      annotation_colors <- list()
      for (col in colnames(annotation_data)) {
        if (is.numeric(annotation_data[[col]])) {
          # Continuous variable
          annotation_colors[[col]] <- circlize::colorRamp2(
            c(min(annotation_data[[col]], na.rm = TRUE), max(annotation_data[[col]], na.rm = TRUE)),
            c("lightyellow", "red")
          )
        } else {
          # Categorical variable
          unique_vals <- unique(annotation_data[[col]])
          annotation_colors[[col]] <- setNames(
            RColorBrewer::brewer.pal(n = length(unique_vals), "Set1"),
            unique_vals
          )
        }
      }
      
      # Create row annotations
      row_annotation <- do.call(rowAnnotation, c(annotation_data, list(col = annotation_colors)))
    } else {
      stop("The 'annotations' argument must be a character vector of column names or NULL.")
    }
  }
  
  # Create a matrix for the heatmap
  heatmap_data <- expanded_data %>%
    pivot_wider(names_from = Bin_start, values_from = Count, values_fill = list(Count = 0)) %>%
    column_to_rownames("sample") %>%
    as.data.frame()
  
  # Ensure the column names (bins) are ordered numerically
  col_order <- as.numeric(colnames(heatmap_data))  # Convert column names to numeric
  heatmap_data <- heatmap_data[, order(col_order)]  # Reorder columns by numeric bin order
  
  # Convert to matrix
  heatmap_matrix <- as.matrix(heatmap_data)
  
  # Determine if the insertion legend should be shown
  max_insertion <- max(heatmap_matrix, na.rm = TRUE)
  show_legend <- ifelse(max_insertion > 1, TRUE, FALSE)
  
  # Create the heatmap
  heatmap <- Heatmap(
    heatmap_matrix,
    name = "Insertions",
    col = colorRamp2(c(0, max_insertion), c("white", "darkblue")),
    show_row_names = TRUE,
    show_column_names = FALSE,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    left_annotation = row_annotation,  # Add the annotations if provided
    heatmap_legend_param = list(
      title = if (show_legend) "Insertions" else NULL
    ),
    show_heatmap_legend = show_legend  # Disable the heatmap legend entirely if needed
  )
  
  # Draw the heatmap
  draw(heatmap)
  
  return(expanded_data)
}

# Function to calculate the number of genes exceeding thresholds
count_genes_by_threshold <- function(df, column_name = "freq_samples_effected", thresholds = c(0.5, 0.75, 0.9)) {
  # Ensure the column exists in the DataFrame
  if (!column_name %in% colnames(df)) {
    stop(paste("Column", column_name, "not found in the DataFrame"))
  }
  
  # Calculate counts for each threshold
  results <- sapply(thresholds, function(threshold) {
    sum(df[[column_name]] > threshold)
  })
  
  # Create a DataFrame of results
  results_df <- data.frame(Threshold = thresholds, Count = results)
  
  return(results_df)
}

plot_pathway_gene_counts <- function(df, n_descriptions=nrow(df), n_genes=10) {
  # Subset to the top n_descriptions pathways
  top_df <- df[1:n_descriptions, ]
  
  # Extract and split the geneIDs into individual genes
  gene_list <- unlist(strsplit(top_df$geneID, split = "/"))
  
  # Count the frequency of each gene
  gene_counts <- table(gene_list)
  
  # Sort the genes by frequency in decreasing order
  gene_counts_sorted <- sort(gene_counts, decreasing = TRUE)
  
  # Convert to a data frame for plotting
  gene_counts_df <- as.data.frame(gene_counts_sorted)
  colnames(gene_counts_df) <- c("Gene", "Count")
  
  # If n_genes is specified, select the top genes
  if (!is.null(n_genes)) {
    n_genes <- min(n_genes, nrow(gene_counts_df))
    gene_counts_df <- gene_counts_df[1:n_genes, ]
  }
  
  # Plot using ggplot2 with specific fill color
  p <- ggplot(gene_counts_df, aes(x = Count, y = reorder(Gene, Count))) +
    geom_bar(stat = "identity", fill = "#5FBFF9") +
    labs(
      title = paste("Top", ifelse(is.null(n_genes), nrow(gene_counts_df), n_genes),
                    "genes in ", n_descriptions, "significant pathways"),
      x = "Frequency",
      y = "Genes"
    ) +
    theme(
      axis.text.y = element_text(size = 8)
    )
  
  print(p)
  
  return(gene_counts_df)
} 

analyze_gene_mutations <- function(sample_genes, ora, n_descriptions=nrow(ora), column, y_lab) {
  if ("Gene_name" %in% colnames(sample_genes)) {
    colnames(sample_genes)[colnames(sample_genes) == "Gene_name"] <- "geneID"
  }
  
  # Subset to the top n_descriptions pathways
  top_df <- ora[1:n_descriptions, ]
  
  # Extract and split the geneIDs into individual genes
  gene_list <- unlist(strsplit(top_df$geneID, split = "/"))
  
  # Count the frequency of each gene
  gene_counts <- table(gene_list)
  
  # Sort the genes by frequency in decreasing order
  gene_counts_sorted <- sort(gene_counts, decreasing = TRUE)
  
  # Convert to a data frame for plotting
  gene_counts_df <- as.data.frame(gene_counts_sorted)
  colnames(gene_counts_df) <- c("geneID", "Count")
  
  # Merge the data frames on 'geneID'
  merged_df <- merge(gene_counts_df, sample_genes, by = "geneID", all.x = TRUE)
  
  # Replace NA values with 0 for genes not present in mutation data
  merged_df[is.na(merged_df)] <- 0
  
  # View the merged data
  print("Merged Data:")
  print(head(merged_df))
  
  # Ensure numeric data
  merged_df$Count <- as.numeric(merged_df$Count)
  merged_df[[column]] <- as.numeric(merged_df[[column]])
  
  # Scatter Plot with gene labels
  p1 <- ggplot(merged_df, aes_string(x = "Count", y = column, label = "geneID")) +
    geom_point(color = "#0080A3", size = 3) +
    geom_text_repel(vjust = -1, size = 3, box.padding = 0.3, max.overlaps = 20 ) +  # Non-overlapping labels
    labs(
      x = "Frequency of gene in significant pathways",
      y = y_lab 
    )
  
  # Bar Plot
  # Prepare plot data
  plot_data <- merged_df[, c("geneID", "Count", column)]
  
  # Melt the data for plotting
  plot_data_melt <- melt(plot_data, id.vars = "geneID")
  
  # Set the factor levels for 'variable' to ensure correct legend and ordering
  plot_data_melt$variable <- factor(plot_data_melt$variable, levels = c("Count", column))
  
  # Create a named vector for 'values' where the names match the levels in 'variable'
  fill_values <- setNames(c("#5FBFF9", "#F8766D"), c("Count", column))
  
  # Plotting
  p2 <- ggplot(plot_data_melt, aes(x = geneID, y = value, fill = variable)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(
      x = "Genes",
      y = "Counts"
    ) +
    scale_fill_manual(
      values = fill_values,  # Use the named vector here
      name = "Metric",
      labels = c("Pathway Count", column)
    ) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1)
    )
  
  # Return plots
  #list(scatter_plot = p1, bar_plot = p2, merged_data = merged_df)
  return(p1)
}

filter_cilia_genes <- function(gene_vector, data_frame) {
  # Split the vector into individual gene names
  gene_list <- unlist(strsplit(gene_vector, "/"))
  
  # Get unique gene names
  unique_genes <- unique(gene_list)
  
  # Filter the data frame for matching gene names in the Gene_name column
  filtered_data <- data_frame[data_frame$Gene_name %in% unique_genes, ]
  
  # Sort the filtered data by the num_samples_effected column in descending order
  sorted_data <- filtered_data[order(-filtered_data$num_samples_effected), ]
  
  return(sorted_data)
}

# gene list input for compareCluster 
# frequency of samples with variant 
calculate_gene_frequency <- function(df, nsample_thresh, filter_exon = FALSE) {
  # Extract number of unique samples
  num_samples <- length(unique(df$sample))  # Total number of unique samples
  
  # Optionally filter for rows where Location contains 'exon'
  if (filter_exon) {
    df <- df[grepl("exon", df$Location), ]
  }
  
  # Count the number of unique samples affected for each gene
  gene_sample_count <- df %>%
    group_by(Gene_name) %>%
    summarise(
      num_samples_mutated = n_distinct(sample),
      .groups = "drop"
    )
  
  # Calculate frequency of samples affected
  gene_sample_count$freq_samples_mutated <- gene_sample_count$num_samples_mutated / num_samples
  
  # Filter for genes with at least `nsample_thresh` samples affected
  gene_sample_count <- subset(gene_sample_count, num_samples_mutated >= nsample_thresh)
  
  # Sort by frequency of samples affected in descending order
  gene_sample_count <- gene_sample_count[order(-gene_sample_count$freq_samples_mutated), ]
  
  # Convert to a named vector with frequency as values and gene names as names
  gene_vec <- setNames(
    gene_sample_count$freq_samples_mutated, 
    gene_sample_count$Gene_name
  )
  
  return(gene_vec)
}

calculate_gene_mut_sample_frequency <- function(df, nsample_thresh = 0, filter_exon = FALSE) {
  # Extract the total number of unique samples
  total_samples <- length(unique(df$sample))
  
  # Optionally filter for rows where Location contains 'exon'
  if (filter_exon) {
    df <- df[grepl("exon", df$Location), ]
  }
  
  # Count the number of unique samples affected for each gene
  gene_sample_count <- df %>%
    group_by(Gene_name, gene_size) %>%
    summarise(
      num_samples_effected = n_distinct(sample),
      .groups = "drop"
    )
  
  # Calculate mutation count (total occurrences of the gene)
  gene_total_count <- as.data.frame(table(df$Gene_name), stringsAsFactors = FALSE)
  colnames(gene_total_count) <- c("geneID", "mutation_count")
  
  # Merge sample count and total mutation count
  gene_df <- merge(gene_sample_count, gene_total_count, by.x = "Gene_name", by.y = "geneID")
  
  # Calculate frequencies
  gene_df$freq_mutations <- gene_df$mutation_count / gene_df$num_samples_effected 
  gene_df$freq_samples_effected <- gene_df$num_samples_effected / total_samples
  
  # Normalize by gene size
  gene_df$freq_mutations_normalized <- gene_df$freq_mutations / gene_df$gene_size
  gene_df$freq_samples_effected_normalized <- gene_df$freq_samples_effected / gene_df$gene_size
  
  # Filter for genes with at least `nsample_thresh` samples affected
  gene_df <- subset(gene_df, num_samples_effected >= nsample_thresh)
  
  return(gene_df)
}

# Function to create a named vector
create_named_vector <- function(df, column_name) {
  # Ensure the specified column exists in the DataFrame
  if (!(column_name %in% colnames(df))) {
    stop("Specified column does not exist in the DataFrame.")
  }
  
  # Create the named vector
  named_vector <- setNames(df[[column_name]], df$Gene_name)
  
  # Sort the named vector in decreasing order
  named_vector <- sort(named_vector, decreasing = TRUE)
  
  return(named_vector)
}

replace_gene_names <- function(df) {
  # Hardcoded mapping as a dataframe
  mapping <- data.frame(
    Current = c("C1orf100", "C12orf40", "CASTOR3", "C18orf25", "C6orf201", 
                "RPSAP58", "C4orf47", "C5orf49", "DDX58", "CLECL1", 
                "MRPS36", "ZNRD1ASP", "KIAA1522", "CSNKA2IP", "ARNTL2", 
                "FAM102B", "C11orf53", "ZBED9", "BMT2", "FAM102A", 
                "ZC3H12A-DT", "CARD17"),
    Replacement = c("SPMIP3", "REDIC1", "CASTOR3P", "ARK2N", "TEX56P", 
                    "RPSA2", "CFAP96", "CFAP90", "RIGI", "CLECL1P", 
                    "KGD4", "POLR1HASP", "NHSL3", "CSNK2A2IP", "BMAL2", 
                    "EEIG2", "POU2AF2", "SCAND3", "SAMTOR", "EEIG1", 
                    "LITATS1", "CARD17P")
  )
  
  # Replace Gene_name values based on the mapping
  df$Gene_name <- ifelse(
    df$Gene_name %in% mapping$Current,  # Check if Gene_name matches a value in the 'Current' column
    mapping$Replacement[match(df$Gene_name, mapping$Current)],  # Replace with corresponding value
    df$Gene_name  # Keep the original value if no match
  )
  
  return(df)
}

add_gene_size_todf <- function(df_te, df_size) {
  # Ensure gene_size has unique entries per gene (take first occurrence if duplicates)
  df_size_unique <- df_size %>%
    distinct(name2, .keep_all = TRUE)

  # Merge te_aff_split with gene_size where Gene_name matches name2
  merged_df <- df_te %>%
    left_join(df_size_unique, by = c("Gene_name" = "name2"))

  # Filter rows where gene_size is NA
  no_size_rows <- merged_df %>%
    filter(is.na(gene_size))

  # Count and list gene names with no size
  count_no_size <- nrow(no_size_rows)
  gene_names_no_size <- unique(no_size_rows$Gene_name)

  # Print information
  cat("Number of genes with no size information:", count_no_size, "\n")
  cat("Gene names with no size information:\n", paste(gene_names_no_size, collapse = ", "), "\n")

  # Return the merged dataframe
  return(merged_df)
}

perform_ora <- function(df, nsample_thresh = 0, filter_exon = FALSE) {
  # Step 2: Get unique genes for ORA
  geneList <- unique(df$Gene_name)
  
  # Step 3: Perform ORA using enrichGO
  ora<- enrichGO(
      gene          = geneList,
      OrgDb         = org.Hs.eg.db,
      keyType       = "SYMBOL",
      ont           = "BP",               
      pAdjustMethod = "BH",               
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.1
  )
  
  # Return ORA results
  return(ora)
}

# Run pathway analysis with multiple parameter combinations and organize results
run_pathway_parameter_sweep <- function(data_list, param_grid, output_base_dir,
                                       output_subdir = NULL,
                                       csv_dir = NULL,
                                       analysis_name = "pathway_sweep",
                                       create_plots = TRUE,
                                       plot_formats = c("png", "pdf")) {
  # data_list: named list of input data (e.g., list(tp53_specific = df1, cancer_specific = df2))
  # param_grid: data frame with columns for each parameter to vary
  # output_base_dir: base directory for plots
  # output_subdir: subdirectory within output_base_dir for plots (e.g., "pathway")
  # csv_dir: directory for CSV files (if NULL, uses output_base_dir)

  cat("\n========================================\n")
  cat("PATHWAY PARAMETER SWEEP\n")
  cat("========================================\n")
  cat("Analysis:", analysis_name, "\n")
  cat("Parameter combinations:", nrow(param_grid), "\n")
  cat("Data sources:", length(data_list), "\n")
  cat("Plot directory:", output_base_dir, "\n")
  if (!is.null(csv_dir)) cat("CSV directory:", csv_dir, "\n")
  cat("\n")

  # Set up plot directory
  if (!is.null(output_subdir)) {
    plot_dir <- file.path(output_base_dir, output_subdir)
  } else {
    plot_dir <- output_base_dir
  }
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

  # Set up CSV directory
  if (is.null(csv_dir)) {
    csv_dir <- plot_dir
  }
  dir.create(csv_dir, recursive = TRUE, showWarnings = FALSE)

  # Initialize results tracking
  all_results <- list()
  summary_table <- data.frame()

  # Loop through each parameter combination
  for (i in 1:nrow(param_grid)) {
    params <- param_grid[i, ]

    # Create parameter suffix for file naming with new naming convention
    param_names <- names(params)
    param_values <- as.character(params)
    param_suffix <- paste(
      ifelse(param_names == "min_samples", paste0("min", param_values),
      ifelse(param_names == "p_gene", paste0("pgene", param_values),
      ifelse(param_names == "p_pathway", paste0("ppathway", param_values),
      ifelse(param_names == "q_pathway", paste0("qpathway", param_values),
      paste0(param_names, param_values))))),
      collapse = "_"
    )

    cat("\n--- Parameter set", i, "of", nrow(param_grid), "---\n")
    print(params)

    # Loop through each data source
    for (data_name in names(data_list)) {
      # Extract dataset name (e.g., "tp53" from "tp53_min5")
      dataset_name <- sub("_min.*", "", data_name)

      # Create run_id in format: pathway_dataset_paramvalues
      run_id <- paste0(analysis_name, "_", dataset_name, "_", param_suffix)
      cat("\nProcessing:", run_id, "\n")

      tryCatch({
        # Get data and apply filters based on parameters
        data <- data_list[[data_name]]

        # Apply min_samples filter if parameter exists
        if ("min_samples" %in% names(params)) {
          # Assume data has 'sample' column to count
          if ("sample" %in% colnames(data)) {
            gene_sample_counts <- data %>%
              group_by(Gene_name) %>%
              summarise(n_samples = n_distinct(sample), .groups = "drop") %>%
              filter(n_samples >= params$min_samples)

            data <- data %>% filter(Gene_name %in% gene_sample_counts$Gene_name)
            cat("  Genes passing min_samples filter:", length(unique(data$Gene_name)), "\n")
          }
        }

        # Apply p_gene filter if parameter exists and p_value/fdr column present
        if ("p_gene" %in% names(params)) {
          if ("fdr" %in% colnames(data)) {
            data <- data %>% filter(fdr < params$p_gene)
            cat("  Genes passing p_gene (fdr) filter:", length(unique(data$Gene_name)), "\n")
          } else if ("p_value" %in% colnames(data)) {
            data <- data %>% filter(p_value < params$p_gene)
            cat("  Genes passing p_gene (p_value) filter:", length(unique(data$Gene_name)), "\n")
          }
        }

        if (nrow(data) == 0) {
          cat("  WARNING: No data remaining after filters. Skipping.\n")
          next
        }

        # Perform ORA with custom cutoffs if provided
        p_pathway <- ifelse("p_pathway" %in% names(params), params$p_pathway, 0.05)
        q_pathway <- ifelse("q_pathway" %in% names(params), params$q_pathway, 0.1)

        ora_result <- perform_ora_custom_cutoffs(data,
                                                 p_pathway = p_pathway,
                                                 q_pathway = q_pathway)

        # Save results
        if (!is.null(ora_result) && nrow(as.data.frame(ora_result)) > 0) {
          # Save full ORA table to CSV directory
          write.csv(as.data.frame(ora_result),
                   file.path(csv_dir, paste0(run_id, ".csv")),
                   row.names = FALSE)

          # Save gene list to CSV directory
          write.table(unique(data$Gene_name),
                     file.path(csv_dir, paste0(run_id, "_genes.txt")),
                     row.names = FALSE, col.names = FALSE, quote = FALSE)

          # Record summary metrics
          summary_row <- data.frame(
            run_id = run_id,
            data_source = data_name,
            params,
            n_genes_input = length(unique(data$Gene_name)),
            n_pathways_enriched = nrow(as.data.frame(ora_result)),
            top_pathway = as.data.frame(ora_result)$Description[1],
            top_pathway_pval = as.data.frame(ora_result)$pvalue[1],
            stringsAsFactors = FALSE
          )
          summary_table <- rbind(summary_table, summary_row)

          # Create plots if requested
          if (create_plots && nrow(as.data.frame(ora_result)) > 0) {
            for (fmt in plot_formats) {
              # Dotplot - format: pathway_dataset_graphtype_params
              p_dot <- dotplot(ora_result, showCategory = 20) +
                ggtitle(paste0(run_id, "\n", nrow(as.data.frame(ora_result)), " enriched pathways"))
              ggsave(file.path(plot_dir, paste0(analysis_name, "_", dataset_name, "_dot_", param_suffix, ".", fmt)),
                     plot = p_dot, width = 10, height = 8)

              # Enrichment map (if enough pathways)
              if (nrow(as.data.frame(ora_result)) >= 5) {
                ora_pairwise <- pairwise_termsim(ora_result)
                p_emap <- emapplot(ora_pairwise, showCategory = 20, pie = "count")
                ggsave(file.path(plot_dir, paste0(analysis_name, "_", dataset_name, "_emap_", param_suffix, ".", fmt)),
                       plot = p_emap, width = 12, height = 10)
              }
            }
          }

          all_results[[run_id]] <- ora_result
          cat("  SUCCESS:", nrow(as.data.frame(ora_result)), "pathways enriched\n")
        } else {
          cat("  No significant enrichment found\n")
          summary_row <- data.frame(
            run_id = run_id,
            data_source = data_name,
            params,
            n_genes_input = length(unique(data$Gene_name)),
            n_pathways_enriched = 0,
            top_pathway = NA,
            top_pathway_pval = NA,
            stringsAsFactors = FALSE
          )
          summary_table <- rbind(summary_table, summary_row)
        }

      }, error = function(e) {
        cat("  ERROR:", e$message, "\n")
      })
    }
  }

  # Save master summary table to CSV directory
  summary_file <- file.path(csv_dir, paste0(analysis_name, "_summary.csv"))
  write.csv(summary_table, summary_file, row.names = FALSE)

  cat("\n========================================\n")
  cat("PARAMETER SWEEP COMPLETE\n")
  cat("========================================\n")
  cat("Total runs:", nrow(summary_table), "\n")
  cat("Successful runs:", sum(summary_table$n_pathways_enriched > 0, na.rm = TRUE), "\n")
  cat("Plots saved to:", plot_dir, "\n")
  cat("CSV files saved to:", csv_dir, "\n")
  cat("Summary table:", summary_file, "\n\n")

  return(list(
    summary = summary_table,
    results = all_results,
    plot_dir = plot_dir,
    csv_dir = csv_dir
  ))
}

# Perform ORA with custom p-value and q-value cutoffs
perform_ora_custom_cutoffs <- function(df, p_pathway = 0.05, q_pathway = 0.1,
                                      nsample_thresh = 0, filter_exon = FALSE, gene_col = "Gene_name") {
  # Get unique genes for ORA
  geneList <- unique(df[[gene_col]])

  if (length(geneList) < 3) {
    cat("Too few genes (n =", length(geneList), ") for pathway analysis\n")
    return(NULL)
  }

  # Perform ORA using enrichGO
  ora <- enrichGO(
    gene          = geneList,
    OrgDb         = org.Hs.eg.db,
    keyType       = "SYMBOL",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = p_pathway,
    qvalueCutoff  = q_pathway
  )

  return(ora)
}

perform_ora_group <- function(df, group_column, n = 5) {
  # Step 1: Count unique samples per group
  group_sample_counts <- df %>%
    distinct(sample, !!sym(group_column)) %>%
    count(!!sym(group_column), name = "sample_count")
  
  # Step 2: Identify groups with at least `n` samples
  valid_groups <- group_sample_counts %>%
    filter(sample_count >= n) %>%
    pull(!!sym(group_column))
  
  print(valid_groups)
  
  # Step 3: Filter the original dataframe to include only valid groups
  df_filtered <- df %>%
    filter(!!sym(group_column) %in% valid_groups)
  
  if (nrow(df_filtered) == 0) {
    stop("No groups have at least the required number of samples.")
  }
  
  print(table(df_filtered$tumor_type))
  
  # Step 4: Create gene clusters based on the filtered group
  geneClusters_group <- df_filtered %>%
    group_by(!!sym(group_column)) %>%
    summarise(Genes = list(unique(Gene_name))) %>%
    deframe()  # Convert to named list
  
  # Perform ORA with compareCluster
  ora_result <- compareCluster(
    geneCluster = geneClusters_group,
    fun = "enrichGO",
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.1
  )
  
  # Return the ORA results
  return(ora_result)
}

perform_gsea <- function(df, column, nsample_thresh = 0, filter_exon = FALSE) {
  # Calculate gene frequency
  frequencies <- calculate_gene_mut_sample_frequency(df,
                                                     nsample_thresh = nsample_thresh, 
                                                     filter_exon = filter_exon)
  print(head(frequencies))
  gene_vector <- create_named_vector(frequencies, column) 
  print(length(gene_vector))
  
  # Perform GSEA
  gsea_results <- gseGO(
      geneList      = gene_vector,
      OrgDb         = org.Hs.eg.db,
      keyType       = "SYMBOL",
      ont           = "BP",              # Biological Process
      minGSSize     = 100,
      maxGSSize     = 500,
      pvalueCutoff  = 1,
      pAdjustMethod = "BH",
      scoreType     = "pos",
      verbose       = TRUE
  )
  
  # Return GSEA results
  return(gsea_results)
}

perform_gsea_tp53 <- function(df, nsample_thresh = 0, filter_exon = FALSE) {
  # Split the dataframe by TP53 status
  te_aff_tp3_control <- df %>% filter(TP53_status == "Control")
  te_aff_tp3_lfs <- df %>% filter(TP53_status == "LFS")
  
  # Calculate gene frequency for Control and LFS
  geneList_tumour_tp53_gsea_control <- calculate_gene_frequency(te_aff_tp3_control, 
                                                                nsample_thresh = nsample_thresh, 
                                                                filter_exon = filter_exon)
  geneList_tumour_tp53_gsea_lfs <- calculate_gene_frequency(te_aff_tp3_lfs, 
                                                            nsample_thresh = nsample_thresh, 
                                                            filter_exon = filter_exon)
  print(length(geneList_tumour_tp53_gsea_control))
  print(length(geneList_tumour_tp53_gsea_lfs))
  
  # Create the input list for compareCluster
  geneClusters_tumour_tp53_gsea <- list(
    Control = geneList_tumour_tp53_gsea_control, 
    LFS = geneList_tumour_tp53_gsea_lfs
  )
  
  # Perform enrichment analysis with compareCluster
  gse_tp53 <- tryCatch(
    compareCluster(
      geneCluster = geneClusters_tumour_tp53_gsea,
      fun = "gseGO",          
      OrgDb = org.Hs.eg.db,
      keyType = "SYMBOL",             
      ont = "BP",                 
      minGSSize = 100,
      maxGSSize = 5000,
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      scoreType = "pos",
      verbose = TRUE
    ),
    error = function(e) NULL  # Return NULL in case of error
  )

  # Return the result
  return(gse_tp53)
}

# Perform ORA with custom cutoffs for TP53 status
perform_ora_tp53_custom_cutoffs <- function(df, p_pathway = 0.05, q_pathway = 0.1, nsample_thresh = 0, gene_col = "Gene_name") {
  # Filter groups by sample threshold if specified
  if (nsample_thresh > 0) {
    # First, identify TP53 status groups that meet the sample threshold
    samples_per_tp53 <- df %>%
      group_by(TP53_status) %>%
      summarise(n_samples = n_distinct(sample.x), .groups = "drop")

    tp53_with_enough_samples <- samples_per_tp53 %>%
      filter(n_samples >= nsample_thresh) %>%
      pull(TP53_status)

    cat("TP53 status groups with >=", nsample_thresh, "samples:", paste(tp53_with_enough_samples, collapse=", "), "\n")

    # If fewer than 2 groups meet threshold, return NULL
    if (length(tp53_with_enough_samples) < 2) {
      cat("Insufficient TP53 status groups (need at least 2 groups with >=", nsample_thresh, "samples)\n")
      return(NULL)
    }

    # Filter to only groups that meet the threshold
    df <- df %>% filter(TP53_status %in% tp53_with_enough_samples)

    cat("After filtering groups:", length(unique(df[[gene_col]])), "genes across", length(tp53_with_enough_samples), "TP53 groups\n")
  }

  # Over-representation analysis
  geneClusters_tumour_tp53_ora <- lapply(split(df[[gene_col]], df$TP53_status), unique)

  # Perform ORA with compareCluster
  ora_tp53 <- compareCluster(
    geneCluster = geneClusters_tumour_tp53_ora,
    fun = "enrichGO",
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    pvalueCutoff = p_pathway,
    qvalueCutoff = q_pathway
  )

  # Return the ORA results
  return(ora_tp53)
}

perform_ora_tp53 <- function(df) {
  # Over-representation analysis
  geneClusters_tumour_tp53_ora <- lapply(split(df$Gene_name, df$TP53_status), unique)

  # Perform ORA with compareCluster
  ora_tp53 <- compareCluster(
    geneCluster = geneClusters_tumour_tp53_ora,
    fun = "enrichGO",
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.1
  )

  # Return the ORA results
  return(ora_tp53)
}

# Perform ORA with custom cutoffs for cancer cohort
perform_ora_cancer_custom_cutoffs <- function(df, p_pathway = 0.05, q_pathway = 0.1, nsample_thresh = 0, gene_col = "Gene_name") {
  # Filter genes by sample threshold if specified
  if (nsample_thresh > 0) {
    # First, identify cohorts that meet the sample threshold
    samples_per_cohort <- df %>%
      group_by(cohort) %>%
      summarise(n_samples = n_distinct(sample.x), .groups = "drop")

    cohorts_with_enough_samples <- samples_per_cohort %>%
      filter(n_samples >= nsample_thresh) %>%
      pull(cohort)

    cat("Cohorts with >=", nsample_thresh, "samples:", paste(cohorts_with_enough_samples, collapse=", "), "\n")

    # If fewer than 2 cohorts meet threshold, return NULL
    if (length(cohorts_with_enough_samples) < 2) {
      cat("Insufficient cohorts (need at least 2 cohorts with >=", nsample_thresh, "samples)\n")
      return(NULL)
    }

    # Filter to only cohorts that meet the threshold
    df <- df %>% filter(cohort %in% cohorts_with_enough_samples)

    cat("After filtering groups:", length(unique(df[[gene_col]])), "genes across", length(cohorts_with_enough_samples), "cohorts\n")
  }

  # Over-representation analysis
  geneClusters_ora <- lapply(split(df[[gene_col]], df$cohort), unique)

  # Perform ORA with compareCluster
  ora_cancer <- compareCluster(
    geneCluster = geneClusters_ora,
    fun = "enrichGO",
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    pvalueCutoff = p_pathway,
    qvalueCutoff = q_pathway
  )

  # Return the ORA results
  return(ora_cancer)
}

perform_ora_cancer <- function(df) {
  # Over-representation analysis
  geneClusters_ora <- lapply(split(df$Gene_name, df$cohort), unique)

  # Perform ORA with compareCluster
  ora_cancer<- compareCluster(
    geneCluster = geneClusters_ora,
    fun = "enrichGO",
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.1
  )

  # Return the ORA results
  return(ora_cancer)
}

# Perform ORA by sample type with custom cutoffs (KICS)
perform_ora_sample_type_custom_cutoffs <- function(df, p_pathway = 0.05, q_pathway = 0.1, nsample_thresh = 0, sample_type_column = "sample_type", gene_col = "Gene_name") {
  # Filter genes by sample threshold if specified
  if (nsample_thresh > 0) {
    # First, identify sample types that meet the sample threshold
    samples_per_type <- df %>%
      group_by(!!sym(sample_type_column)) %>%
      summarise(n_samples = n_distinct(sample.x), .groups = "drop")

    types_with_enough_samples <- samples_per_type %>%
      filter(n_samples >= nsample_thresh) %>%
      pull(!!sym(sample_type_column))

    cat("Sample types with >=", nsample_thresh, "samples:", paste(types_with_enough_samples, collapse=", "), "\n")

    # If fewer than 2 sample types meet threshold, return NULL
    if (length(types_with_enough_samples) < 2) {
      cat("Insufficient sample types (need at least 2 types with >=", nsample_thresh, "samples)\n")
      return(NULL)
    }

    # Filter to only sample types that meet the threshold
    df <- df %>% filter(!!sym(sample_type_column) %in% types_with_enough_samples)

    cat("After filtering groups:", length(unique(df[[gene_col]])), "genes across", length(types_with_enough_samples), "sample types\n")
  }

  # Over-representation analysis by sample type
  geneClusters_sample_type <- lapply(split(df[[gene_col]], df[[sample_type_column]]), unique)

  # Perform ORA with compareCluster
  ora_sample_type <- compareCluster(
    geneCluster = geneClusters_sample_type,
    fun = "enrichGO",
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    pvalueCutoff = p_pathway,
    qvalueCutoff = q_pathway
  )

  # Return the ORA results
  return(ora_sample_type)
}

# Perform ORA by sample type (KICS)
perform_ora_sample_type <- function(df, sample_type_column = "sample_type") {
  # Over-representation analysis by sample type
  geneClusters_sample_type <- lapply(split(df$Gene_name, df[[sample_type_column]]), unique)

  # Perform ORA with compareCluster
  ora_sample_type <- compareCluster(
    geneCluster = geneClusters_sample_type,
    fun = "enrichGO",
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.1
  )

  # Return the ORA results
  return(ora_sample_type)
}

perform_gsea_tp53_tumour_types <- function(df, tumour_types) {
  # Initialize a list to store the comparison results for each tumor type
  results <- list()
  
  # Loop through each tumor type
  for (tumour in tumour_types) {
    # Filter for Control and LFS within the specific tumor type
    te_aff_tp3_control_tt <- df %>% 
      filter(TP53_status == "Control" & tumor_type == tumour)
    te_aff_tp3_lfs_tt <- df %>% 
      filter(TP53_status == "LFS" & tumor_type == tumour)
    
    # Calculate gene frequency for Control and LFS
    geneList_tumour_tp53_gsea_control_tt <- calculate_gene_frequency(te_aff_tp3_control_tt, nsample_thresh = 0, filter_exon = FALSE)
    geneList_tumour_tp53_gsea_lfs_tt <- calculate_gene_frequency(te_aff_tp3_lfs_tt, nsample_thresh = 0, filter_exon = FALSE)
    
    # print length
    print(length(geneList_tumour_tp53_gsea_control_tt))
    print(length(geneList_tumour_tp53_gsea_lfs_tt))
    
    # Create input for compareCluster
    geneClusters_tumour_tp53_gsea_tt <- list(
      Control = geneList_tumour_tp53_gsea_control_tt,
      LFS = geneList_tumour_tp53_gsea_lfs_tt
    )
    
    # Compare enriched GO terms for Biological Process (BP)
    gse_tp53_tt <- tryCatch(
      compareCluster(
        geneCluster = geneClusters_tumour_tp53_gsea_tt,
        fun = "gseGO",          
        OrgDb = org.Hs.eg.db,
        keyType = "SYMBOL",             
        ont = "BP",                 
        minGSSize = 100,
        maxGSSize = 5000,
        pvalueCutoff = 0.05,
        pAdjustMethod = "BH",
        scoreType = "pos",
        verbose = TRUE
      ),
      error = function(e) NULL  # Handle errors gracefully
    )
    
    if (is.null(gse_tp53_tt)) {
      # Save a message indicating no enrichment was found
      results[[tumour]] <- "no enrichment found"
    } else {
      # Visualize the dotplot
      dotplot(gse_tp53_tt, showCategory = 30)
      
      # Save the results in the list
      results[[tumour]] <- gse_tp53_tt
    }
  }
  
  # Return the results
  return(results)
}

analyze_grouped_variants_output <- function(df) {
  # Filter to get only the inconsistent groups
  inconsistent_groups <- df %>%
    filter(!consistent) %>%
    select(group) %>%
    distinct()
  
  # Count the number of inconsistent groups
  inconsistent_count <- nrow(inconsistent_groups)
  
  # Get the details of inconsistent groups
  inconsistent_group_details <- df %>%
    filter(group %in% inconsistent_groups$group)
  
  # Calculate the size of each group
  group_sizes <- df %>%
    group_by(group) %>%
    summarize(group_size = max(end) - min(start)) %>%
    ungroup()
  
  # Merge the group sizes with the inconsistent group details
  inconsistent_group_details <- inconsistent_group_details %>%
    left_join(group_sizes, by = "group")
  
  # Print the number of inconsistent groups
  print(paste("Number of inconsistent groups:", inconsistent_count))
  
  # Print the details of inconsistent groups
  print("Details of inconsistent groups:")
  print(inconsistent_group_details)
  
  # Return the dataframe with inconsistent group details
  return(list(
    inconsistent_count = inconsistent_count,
    inconsistent_groups = inconsistent_group_details
  ))
}

prep_intersect <- function(df, clinical){
  df_lynch <- remove_lynch(df) # remove lynch
  df_clin <- merge_df_clinical(df_lynch, clinical) # add clincial
  df_age <- filter_age(df_clin) # filter by age
  return(df_age)
}

plot_te_location_distribution_absolute <- function(data) {
  
  # Summarize the data: count the number of TEs for each combination of ALT and Location2
  summarized_data <- data %>%
    group_by(ALT, Location2) %>%
    summarise(TE_count = n()) %>%
    ungroup()
  
  # Create the absolute count plot
  g <- ggplot(summarized_data, aes(x = ALT, y = TE_count, fill = Location2)) +
    geom_bar(stat = "identity") +
    labs(x = "TE type", y = "TE count") +
    scale_fill_manual(values = color_palette_6) +
    theme(legend.title = element_blank())  # Remove the legend title
  
  return(g)
}

plot_te_location_distribution_proportion <- function(data) {
  
  # Summarize the data: count the number of TEs for each combination of ALT and Location2
  summarized_data <- data %>%
    group_by(ALT, Location2) %>%
    summarise(TE_count = n()) %>%
    ungroup()
  
  # Calculate the proportion within each ALT group
  proportion_data <- summarized_data %>%
    group_by(ALT) %>%
    mutate(Proportion = TE_count / sum(TE_count)) %>%
    ungroup()
  
  # Create the proportion plot
  g <- ggplot(proportion_data, aes(x = ALT, y = Proportion, fill = Location2)) +
    geom_bar(stat = "identity") +
    labs(x = "TE type", y = "Proportion of TEs") +
    scale_fill_manual(values = color_palette_6) +
    theme(legend.title = element_blank())  # Remove the legend title
  
  return(g)
}

prep_p53_fitness <- function(te_aff_input, fitness_data_input) {
  # Rename columns in fitness_data_input
  names(fitness_data_input) <- c("mutation", "p53_fitness")
  
  # Remove "p." prefix from the mutation column in fitness_data_input
  fitness_data_input$mutation <- gsub("^p\\.", "", fitness_data_input$mutation)
  
  # Merge te_aff_input with fitness_data_input by "mutation"
  merged_data <- merge(te_aff_input, fitness_data_input, by = "mutation", all.x = TRUE)
  
  # Count the number of samples where TP53_status is "LFS" and p53_fitness is NA
  lfs_na_count <- sum(merged_data$TP53_status == "LFS" & is.na(merged_data$p53_fitness))
  
  # Print the count with the message
  message("Number of LFS samples with NA p53_fitness: ", lfs_na_count)
  
  # Return the updated data frame
  return(merged_data)
}

scatter_template <- function(data, independent, dependent, x_limit, y_limit, x_lab, y_lab, shape_column = NULL, colour_column = NULL) {
  
  # Build the base ggplot object with optional shape and color mappings
  g <- ggplot(data, aes(x = !!sym(independent), y = !!sym(dependent)))
  
  # Add shape if the shape_column is provided
  if (!is.null(shape_column)) {
    g <- g + aes(shape = !!sym(shape_column))
  }
  
  # Add color if the color_column is provided
  if (!is.null(colour_column)) {
    g <- g + aes(color = !!sym(colour_column))
  }
  
  # Add the point plot and the smoothing line
  g <- g +
    geom_point() +
    geom_smooth(method = "lm", col = "#5FBFF9") +
    coord_cartesian(xlim = x_limit, ylim = y_limit) +
    labs(x = x_lab, y = y_lab)
  
  return(g)
}

linear_model_one_variable <- function(data, independent, dependent){
  # fit model
  formula <- as.formula(paste(dependent, "~", independent)) # make formula because column names given to function in quotes
  fit <- lm(formula, data)
  # Get the model summary 
  model_summary <- summary(fit)
  print(model_summary) # linear model summary
  # Extract R squared 
  r_squared <- round(model_summary$r.squared, 3)
  # Extract p value
  p_value <- round(model_summary$coefficients[2, 4], 3) # extract specific p value of interest
  # Get correlation
  c <- cor(data[[independent]], data[[dependent]]) # access columns
  # Return the results as a list
  return(list(correlation = c, r_squared = r_squared, p_value = p_value))
}

scatter_template_lm_one_variable <- function(data, independent, dependent, x_limit, y_limit, x_lab, y_lab, shape_column = NULL, colour_column = NULL) {
  data <- data %>%
    filter(is.finite(!!sym(independent)) & is.finite(!!sym(dependent)))
  
  # Compute metrics using the linear model
  metrics <- linear_model_one_variable(data, independent, dependent)
  correlation <- metrics$correlation
  r_squared <- metrics$r_squared
  p_value <- metrics$p_value
  
  # Create the scatter plot with optional shape and color
  g <- scatter_template(data, independent, dependent, x_limit, y_limit, x_lab, y_lab, shape_column = shape_column, colour_column = colour_column)
  
  # Add annotations for correlation, R-squared, and p-value
  g <- g + annotate("text", x = 0.85, y = y_limit[2]- 100, # coordinates of label
                    label = paste("Correlation =", format(correlation, digits = 3), 
                                  "\nR-squared =", format(r_squared, digits = 3), 
                                  "\np-value =", format(p_value, digits = 3)),
                    size = 4, hjust = 0, vjust = 0)
  
  return(g)
}

generate_scatterplots <- function(plot_function, df, independent, types, x_limit, y_limit, x_lab, y_end = "count", ...) {
  # Use lapply to iterate over each independent variable and create the corresponding scatter plot
  invisible(lapply(types, function(type) {
    # Set the y-axis label based on the independent variable and y_end
    y_lab <- if (is.na(type)) {
      paste("TE", y_end)
    } else {
      paste(type, y_end)
    }
    
    # Set the dependent variable; use "total" if type is NA
    dependent_var <- if (is.na(type)) {
      "total"
    } else {
      type
    }
    
    # Call function for each independent variable
    print(plot_function(df, independent = independent, dependent = dependent_var,
                        x_limit = x_limit, y_limit = y_limit, x_lab = x_lab, y_lab= y_lab, ...))
  }))
}

generate_plots <- function(plot_function, df, types, y_end = "count", ...) {
  # Use lapply to iterate over each type and print the corresponding plot
  invisible(lapply(types, function(type) {
    # Set the y-axis label based on type and y_end
    y_lab <- if (is.na(type)) {
      paste("TE", y_end)
    } else {
      paste(type, y_end)
    }
    
    
    # Call the provided plotting function and print the plot
    print(plot_function(df, type = type, y_lab = y_lab, ...))
  }))
}

prep_somatic_fitness <- function(somatic_variants_input, fitness_data_input) {
  # Rename columns in fitness_data_input
  names(fitness_data_input) <- c("aa", "p53_fitness")
  
  # Merge somatic_variants_input with fitness_data_input by "aa"
  merged_variants <- merge(somatic_variants_input, fitness_data_input, by = "aa", all.x = TRUE)
  
  # Count the number of variants without fitness value
  num_variants_without_fitness <- sum(is.na(merged_variants$p53_fitness))
  
  # Print the count with the message
  message("Number of variants without fitness value: ", num_variants_without_fitness)
  
  # Return the updated data frame
  return(merged_variants)
}

prep_variants_tumour <- function(variants_data, conversion, sample_data) {
  # Ensure variants_data$sample is formatted as 4-digit strings
  variants_data <- variants_data %>%
    mutate(variants_sample = sprintf("%04d", as.numeric(KiCS_ID)))
  
  # Filter for KiCS_ID in sample_data$sample
  variants_data <- subset(variants_data, variants_sample %in% sample_data$base_sample)
  
  # Ensure conversion$sample is also 4 digits
  conversion <- conversion %>%
    mutate(sample = sprintf("%04d", as.numeric(sample)))
  
  # Merge the two data frames by sample
  merged_data <- merge(variants_data, conversion, by = "ccp_sample_id", all.x = TRUE)
  
  # Remove original sample columns
  merged_data <- merged_data %>%
    select(-sample) %>%
    mutate(sample = new_sample)  %>% # Create new column `sample` 
    filter(!is.na(sample)) # remove smaples where ccp id didnt line up bc a few dont have in clinical
  
  # Return the processed data
  return(merged_data)
}

prep_variants_germline <- function(variants_data, sample_data) {
  # Pad KiCS_ID with leading zeros to ensure they are four digits
  variants_data$KiCS_ID <- sprintf("%04d", as.numeric(variants_data$KiCS_ID))
  
  # Filter for KiCS_ID in sample_data$sample
  filtered_data <- subset(variants_data, sample %in% sample_data$sample)
  
  # Rename geneSymbol to gene if it exists
  if ("geneSymbol" %in% colnames(filtered_data)) {
    filtered_data <- filtered_data %>% rename("gene" = "geneSymbol")
  }
  
  return(filtered_data)
}

variants_gene_table <- function(variants_data, sample_data) {
  # Further filter for Pathogenic or Likely Pathogenic interpretations
  filtered_data <- subset(filtered_data, interpretation %in% c("Pathogenic", "Likely Pathogenic"))
  
  # Create a table of the geneSymbol column
  gene_symbol_table <- table(filtered_data$gene)
  
  # Return the table
  return(gene_symbol_table)
}

add_variant_columns <- function(data, variant_data, variant_type) {
  # Set the prefix based on the variant type (g for germline, s for somatic)
  prefix <- ifelse(variant_type == "germline", "g_", "s_")
  
  # If germline, filter by P/LP interpretation
  if (variant_type == "germline") {
    unique_genes <- unique(variant_data$gene[variant_data$interpretation %in% c("Pathogenic", "Likely Pathogenic")])
  } else {
    # For somatic, include all unique genes regardless of interpretation
    unique_genes <- unique(variant_data$gene)
  }
  
  # Loop through each gene and create columns for Pathogenic/Likely Pathogenic variants (for germline)
  for (gene in unique_genes) {
    if (variant_type == "germline") {
      # Create a column for Pathogenic/Likely Pathogenic variants for germline
      data[[paste0(prefix, gene)]] <- ifelse(
        data$sample %in% variant_data$sample[variant_data$gene == gene & variant_data$interpretation %in% c("Pathogenic", "Likely Pathogenic")], 
        1, 
        0
      )
    } else {
      # Create a column for all variants for somatic, no interpretation filtering
      data[[paste0(prefix, gene)]] <- ifelse(
        data$sample %in% variant_data$sample[variant_data$gene == gene], 
        1, 
        0
      )
    }
  }
  
  # Return the updated data
  return(data)
}

add_variant_columns <- function(data, variant_data, variant_type) {
  # Set the prefix based on the variant type (g for germline, s for somatic)
  prefix <- ifelse(variant_type == "germline", "g_", "s_")
  
  # If germline, filter by P/LP interpretation
  if (variant_type == "germline") {
    unique_genes <- unique(variant_data$gene[variant_data$interpretation %in% c("Pathogenic", "Likely Pathogenic")])
  } else {
    # For somatic, include all unique genes regardless of interpretation
    unique_genes <- unique(variant_data$gene)
  }
  
  # Loop through each gene and create columns for Pathogenic/Likely Pathogenic variants (for germline)
  for (gene in unique_genes) {
    if (variant_type == "germline") {
      # Create a column for Pathogenic/Likely Pathogenic variants for germline
      data[[paste0(prefix, gene)]] <- ifelse(
        data$sample %in% variant_data$sample[variant_data$gene == gene & variant_data$interpretation %in% c("Pathogenic", "Likely Pathogenic")], 
        1, 
        0
      )
    } else {
      # Create a column for all variants for somatic, no interpretation filtering
      data[[paste0(prefix, gene)]] <- ifelse(
        data$sample %in% variant_data$sample[variant_data$gene == gene], 
        1, 
        0
      )
    }
  }
  
  # Return the updated data
  return(data)
}

add_variant_columns <- function(data, variant_data, variant_type) {
  # Set the prefix based on the variant type (g for germline, s for somatic)
  prefix <- ifelse(variant_type == "germline", "g_", "s_")
  
  # If germline, filter by P/LP interpretation
  if (variant_type == "germline") {
    unique_genes <- unique(variant_data$gene[variant_data$interpretation %in% c("Pathogenic", "Likely Pathogenic")])
  } else {
    # For somatic, include all unique genes regardless of interpretation
    unique_genes <- unique(variant_data$gene)
  }
  
  # Loop through each gene and create columns for Pathogenic/Likely Pathogenic variants (for germline)
  for (gene in unique_genes) {
    if (variant_type == "germline") {
      # Create a column for Pathogenic/Likely Pathogenic variants for germline
      data[[paste0(prefix, gene)]] <- ifelse(
        data$base_sample %in% variant_data$sample[variant_data$gene == gene & variant_data$interpretation %in% c("Pathogenic", "Likely Pathogenic")], 
        1, 
        0
      )
    } else {
      # Create a column for all variants for somatic, no interpretation filtering
      data[[paste0(prefix, gene)]] <- ifelse(
        data$sample %in% variant_data$sample[variant_data$gene == gene], 
        1, 
        0
      )
    }
  }
  
  # Return the updated data
  return(data)
}

prep_variant_df <- function(te_kics_df, te_aff_df) {
  
  # Load variants
  kics_germline_variants <- fread("/Users/briannelaverty/Documents/R_Malkin/TE/data/kics_germline_variants_panel.csv", sep = ",", header = TRUE)
  kics_tp53_somatic_variants <- fread("/Users/briannelaverty/Documents/R_Malkin/TE/data/kics_tp53_somatic_variants.csv", sep = ",", header = TRUE)
  
  # Prep variants: Add leading 0 and filter out samples without TE info
  kics_germline_variants <- prep_variants_germline(kics_germline_variants, te_kics_df)
  kics_tp53_somatic_variants <- prep_variants_tumour(kics_tp53_somatic_variants, te_kics_df)
  
  # Number of samples with P/LP variants by gene
  variants_gene_table(kics_germline_variants, te_kics_df)
  
  # Add variant columns
  te_variant_df <- add_variant_columns(data = te_aff_df, variant_data = kics_germline_variants, variant_type = "germline")
  te_variant_df <- add_variant_columns(data = te_variant_df, variant_data = kics_tp53_somatic_variants, variant_type = "somatic")
  
  # Add filler column for germline TP53 (currently absent)
  te_variant_df$g_TP53 <- NA
  
  # Plot preparation
  te_variant_df <- te_variant_df %>%
    mutate(TP53 = case_when(
      TP53_status == "LFS" ~ "Germline",
      g_TP53 == 1 ~ "Germline",
      s_TP53 == 1 ~ "Somatic",
      TRUE ~ "None" # Catch-all case
    ))
  
  return(te_variant_df)
}

plot_tp53_variants <- function(df, chr = NA, type = NA, group, x_lab, y_lab, x_order = NULL, log_scale = FALSE, breaks = NULL) {
  # Define custom colors for None, Somatic, and Germline
  custom_colors <- c("None" = "#DDD8C4", "Somatic" = "#0090B8", "Germline" = "#004052")
  
  # Use helper function to construct the combination label
  combination <- construct_combination_label(chr, type)
  
  # Check if there are any NA values in the group column
  if (any(is.na(df[[group]]))) {
    cat("Group contains NA values. Excluding NA group from the plot and Kruskal-Wallis test.\n")
    df <- df %>% filter(!is.na(!!sym(group)))
  }
  
  # Check if there are at least 2 groups for Kruskal-Wallis test
  unique_groups <- unique(df[[group]])
  unique_groups <- unique_groups[!is.na(unique_groups)]
  
  if (length(unique_groups) < 2) {
    cat("Kruskal-Wallis test requires at least 2 groups. Column:", group, "has", length(unique_groups), "unique value(s):", paste(unique_groups, collapse = ", "), ". Skipping statistical test.\n")
    p_value <- NA
    p_value_formatted <- "NA"
  } else {
    # Perform Kruskal-Wallis test
    test_result <- kruskal.test(reformulate(group, combination), data = df)
    p_value <- test_result$p.value
    p_value_formatted <- formatC(p_value, format = "e", digits = 2)
  }
  
  # Calculate and print medians for each group
  medians <- df %>%
    group_by(!!sym(group)) %>%
    summarise(median_count = median(!!sym(combination), na.rm = TRUE))
  print(medians)
  
  # Plot
  df_clean <- df %>% filter(is.finite(!!sym(combination)))  # Remove non-finite values
  
  plot <- ggplot(df_clean, aes(x = !!sym(group), y = !!sym(combination), fill = !!sym(group))) +
    geom_boxplot(outlier.shape = NA) + 
    geom_jitter(position = position_jitterdodge(jitter.width = 0.2), color = "black", size = 1.5) +
    scale_fill_manual(values = custom_colors) +  # Apply custom colors
    labs(x = x_lab, y = y_lab) + 
    guides(fill = "none") +
    annotate("text", x = Inf, y = Inf, label = paste("p =", p_value_formatted), vjust = 2, hjust = 1, size = 4)
  
  # Apply custom order for x-axis if provided
  if (!is.null(x_order)) {
    plot <- plot + scale_x_discrete(limits = x_order)
  }
  
  # Apply logarithmic scale if log_scale is TRUE
  if (log_scale) {
    plot <- plot + scale_y_continuous(
      trans = scales::log1p_trans(),
      breaks = if (!is.null(breaks)) breaks else waiver()
    )
  }
  
  return(plot)
}

extract_info_fields <- function(df) {
  df <- df %>%
    mutate(
      source = str_extract(INFO, "TD_SRC=[^;]+") %>% str_replace("TD_SRC=", ""),
      subtype = str_extract(INFO, "SUBTYPE=[^;]+") %>% str_replace("SUBTYPE=", ""),
      length = str_extract(INFO, "AVG_LEN=[^;]+") %>% str_replace("AVG_LEN=", ""),
      seq = str_extract(INFO, "SVINSSEQ=[^;]+") %>% str_replace("SVINSSEQ=", "")
    )
  return(df)
}

plot_multisample_scatter_kics <- function(data, te_count_col = "total") {
  # filter for kics 
  data <- data %>%
    filter(cohort=="KiCS")
  
  # Step 1: Filter to keep only patients with multiple samples
  df_multiple_samples <- data %>%
    group_by(base_sample) %>%
    filter(n() > 1) %>%  # Keep only patients with more than one sample
    ungroup()
  
  # Step 2: Plot the TE counts
  p <- ggplot(df_multiple_samples, aes(x = factor(base_sample), y = !!sym(te_count_col), color=disease_state, shape=lesion_type)) +
    geom_jitter(size = 5, alpha = 0.7, position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.4)) +  # Scatter plot with black outline
    scale_y_continuous(trans = scales::log1p_trans(), breaks = c(1, 10, 100, 1000, 2000)) +  # Log1p transformation with specified breaks
    expand_limits(y = 0.1) +  # Ensure y-axis includes small values and avoids cutoff at zero
    labs(
      x = "Patient",
      y = "TE count"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  
  return(p)
}

plot_multisample_scatter_nick <- function(data, te_count_col = "total") {
  # Filter for LFS cohort
  data <- data %>%
    filter(cohort == "LFS")
  
  # Step 1: Filter to keep only patients with multiple samples
  df_multiple_samples <- data %>%
    group_by(base_sample) %>%
    filter(n() > 1) %>%
    ungroup()
  
  # Define a color palette that can handle more tumor types
  tumor_type_colors <- c("lightgreen", "lightcoral", "skyblue", "gold", "purple", "orange", "pink")
  
  # Step 2: Plot the TE counts with tumor type as color
  p <- ggplot(df_multiple_samples, aes(x = factor(base_sample), y = !!sym(te_count_col), fill = tumor_type)) +
    geom_jitter(shape = 21, size = 5, alpha = 0.7, 
                position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.4)) +
    scale_color_manual(values = tumor_type_colors, name = "Tumor Type") +
    expand_limits(y = 0.1) +  # Ensure y-axis includes small values and avoids cutoff at zero
    labs(
      x = "Patient",
      y = "TE Count",
      color = "Tumor Type"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(p)
}

plot_multisample_scatter_clonality <- function(data, clonal_df, clonal_column, te_count_col = "total", legend_lab) {
  # Filter for LFS cohort
  data <- data %>%
    filter(cohort == "LFS")
  
  # Add a filter to exclude samples containing "merge" in the "sample" column**
  data <- data %>%
    filter(!grepl("merge", sample, ignore.case = TRUE))
  
  # Step 1: Filter to keep only patients with multiple samples
  df_multiple_samples <- data %>%
    group_by(base_sample) %>%
    filter(n() > 1) %>%
    ungroup()
  
  # Join with clonality data on "sample"
  df_multiple_samples <- df_multiple_samples %>%
    left_join(clonal_df, by = "base_sample") %>%
    mutate(
      ssm_prop_clonal = as.numeric(ssm_prop_clonal),  
      cnv_prop_clonal = as.numeric(cnv_prop_clonal)
    )
  
  # Step 2: Plot the TE counts with clonality as color (using ssm_prop_clonal as an example)
  p <- ggplot(df_multiple_samples, aes(x = factor(base_sample), y = !!sym(te_count_col), fill= !!sym(clonal_column))) +
    geom_jitter(shape = 21, size = 5, alpha = 0.7, 
                position = position_jitter(width = 0.4, height = 0)) +
    scale_fill_gradient(low = "white", high = "#007FA3", name=legend_lab) +
    expand_limits(y = 0.1) +  # Ensure y-axis includes small values and avoids cutoff at zero
    labs(
      x = "Patient",
      y = "Repeat count"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(p)
}

plot_multisample_line <- function(data) {
  # Filter for patients with multiple samples
  df_multiple_samples <- data %>%
    group_by(base_sample) %>%
    filter(n() > 1) %>%
    ungroup()
  
  # Plot age at enrollment vs. total with lines connecting samples for each patient
  p <- ggplot(df_multiple_samples, aes(x = age_at_enrollment, y = total, group = base_sample)) +
    geom_line(aes(color = disease_state), size = 1, alpha = 0.6) +  # Line connecting points for each patient
    geom_point(aes(color = disease_state, shape = lesion_type), size = 3) +  # Points for each sample
    scale_y_continuous(trans = scales::log1p_trans(), breaks = c(10, 100, 1000, 2000)) +  # Log1p transformation with specified breaks
    scale_color_brewer(palette = "Set1") +  # Use a color palette for disease states
    labs(
      x = "Age at Enrollment",
      y = "Total TEs",
      color = "Disease State",
      shape = "Lesion Type"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(p)
}

plot_multisample_line_sample <- function(data) {
  # Filter for patients with multiple samples and add a sample order column
  df_multiple_samples <- data %>%
    group_by(base_sample) %>%
    filter(n() > 1) %>%
    mutate(sample_order = row_number()) %>%  # Assign sample order within each patient
    ungroup()
  
  # Plot sample order vs. total with lines connecting samples for each patient
  p <- ggplot(df_multiple_samples, aes(x = sample_order, y = total, group = base_sample)) +
    geom_line(aes(color = disease_state), size = 1, alpha = 0.6) +  # Line connecting points for each patient
    geom_point(aes(color = disease_state, shape = lesion_type), size = 3) +  # Points for each sample
    scale_x_continuous(breaks = unique(df_multiple_samples$sample_order),
                       labels = paste("Sample", unique(df_multiple_samples$sample_order))) +  # Custom x-axis labels
    scale_y_continuous(trans = scales::log1p_trans(), breaks = c(10, 100, 1000, 2000)) +  # Log1p transformation with specified breaks
    scale_color_brewer(palette = "Set1") +  # Use a color palette for disease states
    labs(
      x = "Sample Order",
      y = "TE count",
      color = "Disease State",
      shape = "Lesion Type"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(p)
}

plot_multisample_sametime <- function(data) {
  # Filter for kics
  data <- data %>%
    filter(cohort=="KiCS")
  
  # Filter for patients with two samples taken at the same time
  paired_samples <- data %>%
    group_by(base_sample, age_at_enrollment) %>%
    filter(n() == 2) %>%  # Keep only cases where there are exactly two samples at the same time
    ungroup()
  print(paired_samples[,c("base_sample", "total", "disease_state", "lesion_type")], n=Inf)
  
  # Plot total TEs with color by disease state and shape by lesion type
  ggplot(paired_samples, aes(x = factor(base_sample), y = total, color = disease_state, shape = lesion_type)) +
    geom_jitter(size = 3, alpha = 0.7, position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.4)) +  # Jitter dodge for both axes
    scale_y_continuous(trans = scales::log1p_trans(), breaks = c(10, 50, 100, 200)) +  # Log1p transformation with specified breaks
    labs(
      x = "Patient",
      y = "TE count",
      color = "Disease State",
      shape = "Lesion Type"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

plot_individual_patient_samples <- function(data) {
  # Define consistent mappings for disease_state and lesion_type
  disease_state_levels <- c("initial", "progressive", "relapsed")
  disease_state_colors <- c("initial" = "blue", "progressive" = "green", "relapsed" = "red")
  
  lesion_type_levels <- c("primary", "metastasis")
  lesion_type_shapes <- c("primary" = 16, "metastasis" = 17)
  
  # Filter for KiCS cohort and convert to factors with specified levels
  data <- data %>%
    filter(cohort == "KiCS") %>%
    mutate(
      disease_state = factor(disease_state, levels = disease_state_levels),
      lesion_type = factor(lesion_type, levels = lesion_type_levels)
    )
  
  # Filter for patients with multiple samples
  multiple_samples <- data %>%
    group_by(base_sample) %>%
    filter(n() > 1) %>%  # Keep only patients with more than one sample
    ungroup()
  
  # Print the filtered data for verification
  print(multiple_samples[, c("base_sample", "total",  "sample_topography", "treatment", "age_at_enrollment", "disease_state", "lesion_type", "tumor_type")], n = Inf)
  
  # Get a list of unique patients
  unique_patients <- unique(multiple_samples$base_sample)
  
  # Create a list of individual plots for each patient
  patient_plots <- map(unique_patients, function(patient) {
    # Filter data for the current patient
    patient_data <- multiple_samples %>% filter(base_sample == patient)
    
    # Handle missing or NA values in total
    patient_data <- patient_data %>% filter(!is.na(total))
    
    # Calculate dynamic nudges based on y-axis range
    y_range <- range(patient_data$total, na.rm = TRUE)
    y_nudge <- diff(y_range) * 0.002  # Adjust label upward by 2% of the y-axis range
    
    # Create the plot for the current patient
    p <- ggplot(patient_data, aes(x = age_at_enrollment, y = total, color = disease_state, shape = lesion_type)) +
      geom_point(size = 5, alpha = 0.7) +
      geom_text(
        aes(label = sample_topography),
        hjust = 0,       # Align label to the left of the x position
        vjust = 0,       # Align label to the bottom of the y position
        nudge_x = 0.5,   # Move label slightly to the right
        nudge_y = y_nudge,   # Move label slightly upward
        size = 4,
        check_overlap = TRUE,
        color = "black" 
      ) +
      scale_color_manual(values = disease_state_colors, drop = FALSE) +
      scale_shape_manual(values = lesion_type_shapes, drop = FALSE) +
      scale_y_continuous(trans = scales::log1p_trans()) +
      labs(
        x = "Age at Enrollment",
        y = "TE count",
        color = "Disease State",
        shape = "Lesion Type",
        title = paste0("Patient: ", patient, ", Tumor type: ", paste(unique(patient_data$tumor_type), collapse = ", "))
      ) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
    
    print(p)
  })
  
  # Set names of the list elements to the patient IDs for easy identification
  names(patient_plots) <- unique_patients
}

plot_individual_patient_samples_facet <- function(data) {
  # Define consistent mappings for disease_state and lesion_type
  disease_state_levels <- c("initial", "progressive", "relapsed")
  disease_state_colors <- c("initial" = "blue", "progressive" = "green", "relapsed" = "red")
  
  lesion_type_levels <- c("primary", "metastasis")
  lesion_type_shapes <- c("primary" = 16, "metastasis" = 17)
  
  # Filter for KiCS cohort and convert to factors with specified levels
  data <- data %>%
    filter(cohort == "KiCS") %>%
    mutate(
      disease_state = factor(disease_state, levels = disease_state_levels),
      lesion_type = factor(lesion_type, levels = lesion_type_levels)
    )
  
  # Filter for patients with multiple samples
  multiple_samples <- data %>%
    group_by(base_sample) %>%
    filter(n() > 1) %>%  # Keep only patients with more than one sample
    ungroup()
  
  # Print the filtered data for verification
  print(multiple_samples[, c("base_sample", "total", "sample_topography", "treatment", "age_at_enrollment", "disease_state", "lesion_type", "tumor_type")], n = Inf)
  
  # Create a faceted plot for all patients
  p <- ggplot(multiple_samples, aes(x = age_at_enrollment, y = total, color = disease_state, shape = lesion_type)) +
    geom_point(size = 5, alpha = 0.7) +
#    geom_text(
#      aes(label = sample_topography),
#      hjust = 0,       # Align label to the left of the x position
#      vjust = 0,       # Align label to the bottom of the y position
#      nudge_x = 0.5,   # Move label slightly to the right
#      nudge_y = 0.3,   # Dynamically calculated nudge
#      size = 4,
#      check_overlap = TRUE,
#      color = "black"  # Explicitly set label color
#    ) +
    scale_color_manual(values = disease_state_colors, drop = FALSE) +
    scale_shape_manual(values = lesion_type_shapes, drop = FALSE) +
    scale_y_continuous(trans = scales::log1p_trans()) +
    labs(
      x = "Age at Enrollment",
      y = "TE count",
      color = "Disease State",
      shape = "Lesion Type"
    ) +
    theme(
      axis.text.x = element_blank(),   # Hide x-axis text
      axis.ticks.x = element_blank(),  # Hide x-axis ticks
      axis.title.x = element_blank(),  # Hide x-axis title
      strip.text = element_text(size = 12, face = "bold")
    ) +
    facet_wrap(~ base_sample, scales = "free", ncol = 3)  # Allow both free x and y axes
  
  return(p)
}

plot_specific_sample <- function(data, sample_id, x_nudge, y_nudge, log_scale=FALSE) {
  # Define consistent mappings for disease_state and lesion_type
  disease_state_levels <- c("initial", "progressive", "relapsed")
  disease_state_colors <- c("initial" = "#5FBFF9", "progressive" = '#99ff99', "relapsed" = '#ffcc66')
  
  lesion_type_levels <- c("primary", "metastasis")
  lesion_type_shapes <- c("primary" = 16, "metastasis" = 17)
  
  # Filter for KiCS cohort and convert to factors with specified levels
  data <- data %>%
    filter(cohort == "KiCS") %>%
    mutate(
      disease_state = factor(disease_state, levels = disease_state_levels),
      lesion_type = factor(lesion_type, levels = lesion_type_levels)
    )
  
  # Filter for the specific sample
  specific_sample <- data %>%
    filter(base_sample == sample_id)
  
  # Handle missing or NA values in total
  specific_sample <- specific_sample %>% filter(!is.na(total))
  
  # If no data for the specific sample, return a message
  if (nrow(specific_sample) == 0) {
    stop(paste("No data found for sample ID:", sample_id))
  }
  
  # Create the plot for the specific sample
  p <- ggplot(specific_sample, aes(x = age_at_enrollment, y = total, color = disease_state, shape = lesion_type)) +
    geom_point(size = 5, alpha = 1) +
    geom_text(
      aes(label = sample_topography),
      hjust = 0,       # Align label to the left of the x position
      vjust = 0,       # Align label to the bottom of the y position
      nudge_x = x_nudge,   # Move label slightly to the right
      nudge_y = y_nudge,   # Move label slightly upward
      size = 3,
      check_overlap = FALSE,
      color = "black" 
    ) +
    scale_color_manual(values = disease_state_colors, drop = FALSE) +
    scale_shape_manual(values = lesion_type_shapes, drop = FALSE) +
    labs(
      x = "Age at sample collection (days)",
      y = "TE count",
      color = "Disease state",
      shape = "Lesion type",
      title = paste0("Patient: ", sample_id, ", Tumor type: ", paste(unique(specific_sample$tumor_type), collapse = ", "))
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Apply logarithmic scale if log_scale is TRUE
  if (log_scale) {
    p <- p + scale_y_continuous(
      trans = scales::log1p_trans()
    )
  }
  
  return(p)
}

plot_te_change <- function(data, log_scale = FALSE) {
  # Step 1: Filter for the "KiCS" cohort
  data <- data %>%
    filter(cohort == "KiCS")
  
  # Step 2: Keep only patients with multiple samples
  df_multiple_samples <- data %>%
    group_by(base_sample) %>%
    filter(n() > 1) %>%
    ungroup()
  
  if (nrow(df_multiple_samples) == 0) {
    stop("No patients with more than one sample found.")
  }
  
  # Step 3: Calculate absolute changes for each patient
  change_data <- df_multiple_samples %>%
    arrange(base_sample, age_at_enrollment) %>%  # Order samples by age_at_enrollment
    group_by(base_sample) %>%
    mutate(
      sample_order = row_number(),
      absolute_change = total - lag(total),  # Absolute change
      transition = case_when(
        sample_order == 2 ~ "1-2",
        sample_order == 3 ~ "2-3",
        sample_order == 4 ~ "3-4"
      )
    ) %>%
    filter(!is.na(absolute_change)) %>%
    ungroup()
  
  # Step 4: Plot absolute change with color for transitions
  plot <- ggplot(change_data, aes(
    x = base_sample, y = absolute_change, 
    color = tumor_type, shape = transition  # Color by tumor type, shape by transition
  )) +
    geom_hline(yintercept = 0, color = "black", linetype = "solid") +  # Horizontal line at 0
    geom_point(size = 5, position = position_dodge(width = 0.5)) +  # Dodge to prevent overlap
    scale_shape_manual(
      values = c("1-2" = 16, "2-3" = 17, "3-4" = 15),  # Different shapes for transitions
      name = "Sample Transition"
    ) +
    labs(x = "Patient", y = "Absolute change") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  
  # Step 5: Apply log scale if required
  if (log_scale) {
    plot <- plot + scale_y_continuous(
      trans = scales::pseudo_log_trans(sigma = 10),  # Handles negative and positive values smoothly
      breaks = c(-300, -100, -10, 0, 10, 100, 300)
    )
  }
  
  # Print the plot
  print(plot)
  
  # Return the summarized data for inspection
  return(change_data)
}

generate_te_insertion_heatmap_genes <- function(te_df, clinical_df, genes) {
  # Split the concatenated gene string into a vector
  genes_vector <- unlist(strsplit(genes, "/"))
  
  # Filter the TE insertion data by the specified genes
  te_df_filtered <- te_df %>% 
    filter(Gene_name %in% genes_vector)
  
  # Count TE insertions per gene for each sample
  insertion_counts <- te_df_filtered %>%
    group_by(sample, Gene_name) %>%
    summarise(insertion_count = n(), .groups = 'drop') %>%
    pivot_wider(names_from = Gene_name, values_from = insertion_count, values_fill = 0)
  
  # Merge with clinical data
  merged_df <- clinical_df %>%
    left_join(insertion_counts, by = "sample")
  
  # Extract numeric data for heatmap
  numeric_data <- as.matrix(merged_df[, 144:287])
  numeric_data[is.na(numeric_data) | is.nan(numeric_data) | is.infinite(numeric_data)] <- 0 
  
  # Extract clinical annotations
  clinical_annotations <- merged_df %>%
    dplyr::select(tumor_class, sex, age_at_diagnosis, total)
  
  # Create row annotations
  row_annot <- rowAnnotation(
    `Tumor class`= clinical_annotations$tumor_class,
    Sex = clinical_annotations$sex,
    `Age of onset` = clinical_annotations$age_at_diagnosis,
    `Total insertions` = clinical_annotations$total,
    annotation_name_gp = gpar(fontsize = 10),
    col = list(
      Sex = c("M" = "#0080A3", "F" = "#AB1368"),
      `Age of onset` = colorRamp2(c(min(clinical_annotations$age_at_diagnosis, na.rm = TRUE), 
                                  max(clinical_annotations$age_at_diagnosis, na.rm = TRUE)), 
                                c("white", "#301934")),
      `Total insertions` = colorRamp2(c(min(clinical_annotations$total, na.rm = TRUE), 
                                      max(clinical_annotations$total, na.rm = TRUE)), 
                                    c("white", "#004052"))
    )
  )
  
  # Generate the heatmap
  heatmap <- Heatmap(
    numeric_data,
    name = "Insertions",
    col = colorRamp2(
      c(min(numeric_data, na.rm = TRUE), max(numeric_data, na.rm = TRUE)), 
      c("white", "#004052")
    ),
    show_row_names = FALSE,
    show_column_names = TRUE,
    row_title = NULL,  # Ensures no row title
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    clustering_distance_rows = "euclidean",
    clustering_distance_columns = "euclidean",
    clustering_method_rows = "ward.D2",
    clustering_method_columns = "ward.D2",
    heatmap_legend_param = list(title = "Insertions"),
    left_annotation = row_annot  
  )
  
  draw(heatmap, annotation_legend_side = "right")
}

generate_te_insertion_heatmap_genes_cut <- function(te_df, clinical_df, genes, top_percent) {
  # Split the concatenated gene string into a vector
  genes_vector <- unlist(strsplit(genes, "/"))
  
  # Filter the TE insertion data by the specified genes
  te_df_filtered <- te_df %>% 
    filter(Gene_name %in% genes_vector)
  
  # Count TE insertions per gene for each sample
  insertion_counts <- te_df_filtered %>%
    group_by(sample, Gene_name) %>%
    summarise(insertion_count = n(), .groups = 'drop') %>%
    pivot_wider(names_from = Gene_name, values_from = insertion_count, values_fill = 0)
  
  # Merge with clinical data
  merged_df <- clinical_df %>%
    left_join(insertion_counts, by = "sample")
  
  # Filter to keep only the top n% of samples based on total insertions
  cutoff <- quantile(merged_df$total, probs = 1 - top_percent / 100, na.rm = TRUE)
  filtered_df <- merged_df %>% filter(total >= cutoff)
  
  # Extract numeric data for heatmap
  numeric_data <- as.matrix(filtered_df[, 144:287])
  numeric_data[is.na(numeric_data) | is.nan(numeric_data) | is.infinite(numeric_data)] <- 0 
  
  # Extract clinical annotations
  clinical_annotations <- filtered_df %>%
    dplyr::select(tumor_type, sex, age_at_diagnosis, TP53_status, total)
  
  # Create row annotations
  row_annot <- rowAnnotation(
    `Tumor type` = clinical_annotations$tumor_type,
    `Age of onset` = clinical_annotations$age_at_diagnosis,
    `Germline TP53 status` = clinical_annotations$TP53_status,
    Sex = clinical_annotations$sex,
    `Total insertions` = clinical_annotations$total,
    annotation_name_gp = gpar(fontsize = 10),
    col = list(
      Sex = c("M" = "#0080A3", "F" = "#AB1368"),
      `Germline TP53 status`= c("WT" = "lightgrey", "Mutant" = "darkgreen"),
      `Age of onset` = colorRamp2(c(min(clinical_annotations$age_at_diagnosis, na.rm = TRUE), 
                                  max(clinical_annotations$age_at_diagnosis, na.rm = TRUE)), 
                                c("white", "#301934")),
      `Total insertions` = colorRamp2(
        log1p(c(0, 1, max(clinical_annotations$total, na.rm = TRUE))),  # Log-transform the color mapping range
        c("white", "lightgrey", "#004052")
      )
    ),
    annotation_legend_param = list(
      `Total insertions` = list(
        at = c(0, 1, max(clinical_annotations$total, na.rm = TRUE)),  # Original values for legend
        labels = c("0", "1", as.character(max(clinical_annotations$total, na.rm = TRUE))),
        title = "Total insertions"
      )
    )
  )
  
  # Generate the heatmap
  heatmap <- Heatmap(
    numeric_data,
    name = "Insertions",
    col = colorRamp2(
      c(0, 1, max(numeric_data, na.rm = TRUE)), 
      c("white", "lightgrey", "#004052")  
    ),
    show_row_names = FALSE,
    show_column_names = FALSE,
    row_title = NULL,  # Ensures no row title
    column_title = "Cilia genes",
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    clustering_distance_rows = "euclidean",
    clustering_distance_columns = "euclidean",
    clustering_method_rows = "ward.D2",
    clustering_method_columns = "ward.D2",
    heatmap_legend_param = list(title = "Insertions"),
    left_annotation = row_annot
  )
  
  draw(heatmap, annotation_legend_side = "right")
}


generate_te_insertion_heatmap_genes_cut <- function(te_df, clinical_df, genes, top_percent = NULL, at_least_one_insertion = FALSE) {
  ht_opt$message = FALSE
  
  # Split the concatenated gene string into a vector
  genes_vector <- unlist(strsplit(genes, "/"))
  
  # Filter the TE insertion data by the specified genes
  te_df_filtered <- te_df %>% 
    filter(Gene_name %in% genes_vector)
  
  # Count TE insertions per gene for each sample
  insertion_counts <- te_df_filtered %>%
    group_by(sample, Gene_name) %>%
    summarise(insertion_count = n(), .groups = 'drop') %>%
    pivot_wider(names_from = Gene_name, values_from = insertion_count, values_fill = 0)
  
  # Merge with clinical data
  merged_df <- clinical_df %>%
    left_join(insertion_counts, by = "sample")
  
  # Apply filtering based on top percentage or at least one insertion
  if (!is.null(top_percent)) {
    # Filter to keep only the top n% of samples based on total insertions
    cutoff <- quantile(merged_df$total, probs = 1 - top_percent / 100, na.rm = TRUE)
    filtered_df <- merged_df %>% filter(total >= cutoff)
  } else if (at_least_one_insertion) {
    # Filter to keep only samples with at least one gene insertion
    filtered_df <- merged_df %>% filter(rowSums(as.matrix(merged_df[, 144:287])) > 0)
  } else {
    filtered_df <- merged_df
  }
  
  # Extract numeric data for heatmap
  numeric_data <- as.matrix(filtered_df[, 144:287])
  numeric_data[is.na(numeric_data) | is.nan(numeric_data) | is.infinite(numeric_data)] <- 0 
  
  # Extract clinical annotations
  clinical_annotations <- filtered_df %>%
    dplyr::select(tumor_type, sex, age_at_diagnosis, TP53_status, total)
  
  # Create row annotations
  row_annot <- rowAnnotation(
    `Tumor type` = clinical_annotations$tumor_type,
    `Age of onset` = clinical_annotations$age_at_diagnosis,
    `Germline TP53 status` = clinical_annotations$TP53_status,
    Sex = clinical_annotations$sex,
    `Total insertions` = clinical_annotations$total,
    annotation_name_gp = gpar(fontsize = 10),
    col = list(
      Sex = c("M" = "#0080A3", "F" = "#AB1368"),
      `Germline TP53 status`= c("WT" = "lightgrey", "Mutant" = "darkgreen"),
      `Age of onset` = colorRamp2(c(min(clinical_annotations$age_at_diagnosis, na.rm = TRUE), 
                                    max(clinical_annotations$age_at_diagnosis, na.rm = TRUE)), 
                                  c("white", "#301934")),
      `Total insertions` = colorRamp2(
        c(
          0,
          log1p(1),
          log1p(quantile(clinical_annotations$total, 0.95, na.rm = TRUE)),
          log1p(quantile(clinical_annotations$total, 0.95, na.rm = TRUE)),
          max(clinical_annotations$total, na.rm = TRUE)
        ),
        c("white", "lightgrey", "#004052", "red", "red")  # Ensure no transition from blue to red
      )
    ),
    annotation_legend_param = list(
      `Total insertions` = list(
        at = c(0, 1, quantile(clinical_annotations$total, 0.95, na.rm = TRUE), max(clinical_annotations$total, na.rm = TRUE)),
        labels = c("0", "1", "95%", "Max"),
        col = c("white", "lightgrey", "#004052", "red", "red"),
        title = "Total insertions"
      )
    )
  )
  
  # Generate the heatmap
  heatmap <- Heatmap(
    numeric_data,
    name = "Insertions",
    col = colorRamp2(
      c(0, 1, max(numeric_data, na.rm = TRUE)), 
      c("white", "lightgrey", "#004052")  
    ),
    show_row_names = FALSE,
    show_column_names = FALSE,
    row_title = NULL,  # Ensures no row title
    column_title = "Cilia genes",
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    clustering_distance_rows = "euclidean",
    clustering_distance_columns = "euclidean",
    clustering_method_rows = "ward.D2",
    clustering_method_columns = "ward.D2",
    heatmap_legend_param = list(title = "Insertions"),
    left_annotation = row_annot
  )
  
  draw(heatmap, annotation_legend_side = "right")
}

generate_te_insertion_heatmap_genes_cut_quantiles <- function(te_df, clinical_df, genes, top_percent = NULL, at_least_one_insertion = FALSE) {
  ht_opt$message = FALSE
  
  # Split the concatenated gene string into a vector
  genes_vector <- unlist(strsplit(genes, "/"))
  
  # Filter the TE insertion data by the specified genes
  te_df_filtered <- te_df %>%
    filter(Gene_name %in% genes_vector)
  
  # Count TE insertions per gene for each sample
  insertion_counts <- te_df_filtered %>%
    group_by(sample, Gene_name) %>%
    summarise(insertion_count = n(), .groups = 'drop') %>%
    pivot_wider(names_from = Gene_name, values_from = insertion_count, values_fill = 0)
  
  # Merge with clinical data
  merged_df <- clinical_df %>%
    left_join(insertion_counts, by = "sample")
  
  # Apply filtering based on top percentage or at least one insertion
  if (!is.null(top_percent)) {
    # Filter to keep only the top n% of samples based on total insertions
    cutoff <- quantile(merged_df$total, probs = 1 - top_percent / 100, na.rm = TRUE)
    filtered_df <- merged_df %>% filter(total >= cutoff)
  } else if (at_least_one_insertion) {
    # Filter to keep only samples with at least one gene insertion
    numeric_columns <- setdiff(names(merged_df), names(clinical_df))
    filtered_df <- merged_df %>%
      filter(rowSums(as.matrix(merged_df[, numeric_columns])) > 0)
  } else {
    filtered_df <- merged_df
  }
  
  # Extract numeric data for heatmap
  numeric_columns <- setdiff(names(filtered_df), names(clinical_df))
  numeric_data <- as.matrix(filtered_df[, numeric_columns])
  numeric_data[is.na(numeric_data) | is.nan(numeric_data) | is.infinite(numeric_data)] <- 0 
  
  # Extract clinical annotations
  clinical_annotations <- filtered_df %>%
    dplyr::select(tumor_type, sex, age_at_diagnosis, TP53_status, total)
  
  # Define quantile-based breaks and explicitly include the full range of values
  quantiles <- quantile(clinical_annotations$total, probs = c(0, 0.25, 0.5, 0.75, 0.95, 1), na.rm = TRUE)
  breaks <- c(0, quantiles[2], quantiles[3], quantiles[4], quantiles[5], max(clinical_annotations$total, na.rm = TRUE))
  
  # Define labels for each interval (one less than breaks)
  labels <- c("0", "0-25%", "25-50%", "50-75%", "75-95%", "95-100%")
  
  # Ensure no missing values in clinical_annotations$total
  clinical_annotations$total[is.na(clinical_annotations$total)] <- 0
  
  # Map total values to bins using cut()
  total_insertions_col <- cut(
    clinical_annotations$total,
    breaks = breaks,
    labels = labels[-1],  # Exclude the first label since it's for non-zero bins
    include.lowest = TRUE
  )
  
  # Convert to character to allow direct assignment
  total_insertions_col <- as.character(total_insertions_col)
  
  # Assign "0" to rows where total is exactly 0
  total_insertions_col[clinical_annotations$total == 0] <- "0"
  
  # Convert back to factor and ensure levels include "0"
  total_insertions_col <- factor(total_insertions_col, levels = labels)
  
  # Define discrete colors for each bin
  discrete_colors <- setNames(
    c("white", "lightgrey", "#A9A9A9", "#707070", "#004052", "red"),
    labels
  )
  
  # Debugging print statements
  print("Debugging Breaks:")
  print(breaks)
  
  print("Debugging Labels:")
  print(labels)
  
  print("Debugging Discrete Colors:")
  print(discrete_colors)
  
  print("Debugging Total Insertions Column Levels:")
  print(levels(total_insertions_col))
  
  # Create row annotations
  row_annot <- rowAnnotation(
    `Tumor type` = clinical_annotations$tumor_type,
    `Age of onset` = clinical_annotations$age_at_diagnosis,
    `Germline TP53 status` = clinical_annotations$TP53_status,
    Sex = clinical_annotations$sex,
    `Total insertions` = total_insertions_col,
    annotation_name_gp = gpar(fontsize = 10),
    col = list(
      Sex = c("M" = "#0080A3", "F" = "#AB1368"),
      `Germline TP53 status`= c("WT" = "lightgrey", "Mutant" = "darkgreen"),
      `Age of onset` = colorRamp2(
        c(min(clinical_annotations$age_at_diagnosis, na.rm = TRUE), 
          max(clinical_annotations$age_at_diagnosis, na.rm = TRUE)), 
        c("white", "#301934")
      ),
      `Total insertions` = discrete_colors
    ),
    annotation_legend_param = list(
      `Total insertions` = list(
        at = labels,
        labels = labels,
        col = discrete_colors,
        title = "Total insertions"
      )
    )
  )
  
  # Generate the heatmap
  heatmap <- Heatmap(
    numeric_data,
    name = "Insertions",
    col = colorRamp2(
      c(0, 1, max(numeric_data, na.rm = TRUE)), 
      c("white", "lightgrey", "#004052")   
    ),
    show_row_names = FALSE,
    show_column_names = FALSE,
    row_title = NULL,  # Ensures no row title
    column_title = "Cilia genes",
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    clustering_distance_rows = "euclidean",
    clustering_distance_columns = "euclidean",
    clustering_method_rows = "ward.D2",
    clustering_method_columns = "ward.D2",
    heatmap_legend_param = list(title = "Insertions"),
    left_annotation = row_annot
  )
  
  draw(heatmap, annotation_legend_side = "right")
}

generate_te_insertion_heatmap_genes_cut_quantiles_normalized<- function(te_df, clinical_df, genes, top_percent = NULL, at_least_one_insertion = FALSE) {
  ht_opt$message = FALSE
  
  # Split the concatenated gene string into a vector
  genes_vector <- unlist(strsplit(genes, "/"))
  
  # Filter the TE insertion data by the specified genes
  te_df_filtered <- te_df %>%
    filter(Gene_name %in% genes_vector)
  
  # Count TE insertions per gene for each sample
  insertion_counts <- te_df_filtered %>%
    group_by(sample, Gene_name) %>%
    summarise(insertion_count = n(), .groups = 'drop') %>%
    pivot_wider(names_from = Gene_name, values_from = insertion_count, values_fill = 0)
  
  # Merge with clinical data
  merged_df <- clinical_df %>%
    left_join(insertion_counts, by = "sample")
  
  # Apply filtering based on top percentage or at least one insertion
  if (!is.null(top_percent)) {
    # Filter to keep only the top n% of samples based on total insertions
    cutoff <- quantile(merged_df$total, probs = 1 - top_percent / 100, na.rm = TRUE)
    filtered_df <- merged_df %>% filter(total >= cutoff)
  } else if (at_least_one_insertion) {
    # Filter to keep only samples with at least one gene insertion
    numeric_columns <- setdiff(names(merged_df), names(clinical_df))
    filtered_df <- merged_df %>%
      filter(rowSums(as.matrix(merged_df[, numeric_columns])) > 0)
  } else {
    filtered_df <- merged_df
  }
  
  # Extract numeric data for heatmap
  numeric_columns <- setdiff(names(filtered_df), names(clinical_df))
  numeric_data <- as.matrix(filtered_df[, numeric_columns])
  numeric_data[is.na(numeric_data) | is.nan(numeric_data) | is.infinite(numeric_data)] <- 0 
  
  # Extract clinical annotations
  clinical_annotations <- filtered_df %>%
    dplyr::select(tumor_type, sex, age_at_diagnosis, TP53_status, total)
  
  # Compute rowSum and normalize by total
  row_sums <- rowSums(numeric_data)
  clinical_annotations$Normalized_Insertion_Percentage <- row_sums / filtered_df$total
  
  # Quantile-based binning for Normalized Insertion %
  quantile_breaks <- quantile(clinical_annotations$Normalized_Insertion_Percentage, probs = seq(0, 1, length.out = 6), na.rm = TRUE)
  quantile_colors <- c("lightgrey", "lightblue", "lightgreen", "yellow", "orange")
  
  # Compute midpoints for the intervals
  interval_midpoints <- (head(quantile_breaks, -1) + tail(quantile_breaks, -1)) / 2
  
  # Define quantile labels
  quantile_labels <- paste0("Q", 1:5)  
  
  # Map values to quantile bins
  clinical_annotations$Quantile_Normalized <- cut(
    clinical_annotations$Normalized_Insertion_Percentage,
    breaks = quantile_breaks,
    labels = quantile_labels,
    include.lowest = TRUE
  )

  # Print debugging information
  print("Quantile Breaks:")
  print(quantile_breaks)
  print("Quantile Labels:")
  print(quantile_labels)
  print("Quantile Assignments:")
  print(table(clinical_annotations$Quantile_Normalized))
  print("set names:")
  print(setNames(quantile_colors, quantile_labels))
  
  # Define discrete colors for total insertions
  discrete_colors <- setNames(
    c("white", "lightgrey", "#A9A9A9", "#707070", "#004052", "red"),
    c("0", "0-25%", "25-50%", "50-75%", "75-95%", "95-100%")
  )
  
  # Create row annotations
  row_annot <- rowAnnotation(
    `Tumor type` = clinical_annotations$tumor_type,
    `Age of onset` = clinical_annotations$age_at_diagnosis,
    `Germline TP53 status` = clinical_annotations$TP53_status,
    Sex = clinical_annotations$sex,
    `Normalized Insertion %` = clinical_annotations$Quantile_Normalized,  # Use binned values
    col = list(
      `Normalized Insertion %` = setNames(quantile_colors, quantile_labels),  # Map colors to labels
      Sex = c("M" = "#0080A3", "F" = "#AB1368"),
      `Germline TP53 status` = c("WT" = "lightgrey", "Mutant" = "darkgreen"),
      `Age of onset` = colorRamp2(
        c(
          min(clinical_annotations$age_at_diagnosis, na.rm = TRUE), 
          max(clinical_annotations$age_at_diagnosis, na.rm = TRUE)
        ), 
        c("white", "#301934")
      )
    )
  )
  
  # Generate the heatmap
  heatmap <- Heatmap(
    numeric_data,
    name = "Insertions",
    col = colorRamp2(
      c(0, 1),                # Scaled range
      c("lightgrey", "orange")  # Light grey to orange
    ),
    show_row_names = FALSE,
    show_column_names = FALSE,
    row_title = NULL,         # Ensures no row title
    column_title = "Cilia genes",
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    clustering_distance_rows = "euclidean",
    clustering_distance_columns = "euclidean",
    clustering_method_rows = "ward.D2",
    clustering_method_columns = "ward.D2",
    heatmap_legend_param = list(title = "Insertions"),
    left_annotation = row_annot
  )
  
  draw(heatmap, annotation_legend_side = "right")
}

plot_te_gene_insertions <- function(te_df, clinical_df, genes, xlim, ylim) {
  # Split the concatenated gene string into a vector
  genes_vector <- unlist(strsplit(genes, "/"))
  
  # Filter the TE insertion data by the specified genes
  te_df_filtered <- te_df %>%
    filter(Gene_name %in% genes_vector)
  
  # Count the number of mutated genes per sample
  mutated_genes_count <- te_df_filtered %>%
    group_by(sample) %>%
    summarise(mutated_genes_count = n_distinct(Gene_name), .groups = 'drop')
  
  # Merge mutated genes count with clinical data
  merged_df <- clinical_df %>%
    left_join(mutated_genes_count, by = "sample")
  
  # Replace NA values with 0 for samples with no mutations
  merged_df <- merged_df %>%
    mutate(
      mutated_genes_count = replace_na(mutated_genes_count, 0)
    )
  
  # Fit a linear model
  fit <- lm(total ~ mutated_genes_count, data = merged_df)
  fit_summary <- summary(fit)
  r_squared <- fit_summary$r.squared
  p_value <- coef(fit_summary)[2, 4]
  print(paste("R-squared:", round(r_squared, 3)))
  print(paste("P-value:", signif(p_value, 3)))
  
  # Plot total TE insertions vs. number of mutated genes
  plot <- ggplot(merged_df, aes(x = mutated_genes_count, y = total)) +
    geom_point() +
    (if (!is.null(ylim)) ylim(ylim) else NULL) +
    (if (!is.null(xlim)) xlim(xlim) else NULL) +
    geom_smooth(method = "lm", se = FALSE, color = "blue") +
    labs(
      x = "Number of mutated genes",
      y = "Total repeats"
    ) 
  
  return(plot)
}

# Function to calculate and print percentages
calculate_shared_diff <- function(df1, df2, col_name) {
  # Extract columns by name
  col1 <- df1[[col_name]]
  col2 <- df2[[col_name]]
  
  # Find shared and different elements
  shared_elements <- intersect(col1, col2)
  different_elements <- setdiff(union(col1, col2), shared_elements)
  
  # Calculate percentages
  total_elements <- length(union(col1, col2))
  percent_shared <- (length(shared_elements) / total_elements) * 100
  percent_different <- (length(different_elements) / total_elements) * 100
  
  # Print results
  cat("Percentage of shared elements:", percent_shared, "%\n")
  cat("Percentage of different elements:", percent_different, "%\n")
  
  # Optionally return results as a list
  return(list(percent_shared = percent_shared, percent_different = percent_different))
}

# Function to perform Fisher's test for each subfamily by sample
perform_fisher_test_summary <- function(df, subfamily_col, tp53_status_col) {
  unique_subfamilies <- unique(df[[subfamily_col]])
  
  # Initialize an empty dataframe
  results_df <- data.frame(
    Subfamily = character(),
    Mutant_Percentage = numeric(),
    WT_Percentage = numeric(),
    P_Value = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (subfamily in unique_subfamilies) {
    # Create a contingency table for the current subfamily and TP53_status
    contingency_table <- table(
      Subfamily = df[[subfamily_col]] == subfamily,
      TP53_Status = df[[tp53_status_col]]
    )
    
    # Ensure the table has both rows and columns
    if (nrow(contingency_table) < 2 || ncol(contingency_table) < 2) {
      # Skip this subfamily if the table is invalid
      next
    }
    
    # Perform Fisher's exact test
    fisher_test <- fisher.test(contingency_table)
    
    # Calculate percentages
    mutant_cases <- ifelse("Mutant" %in% colnames(contingency_table), contingency_table[2, "Mutant"], 0)
    wt_cases <- ifelse("WT" %in% colnames(contingency_table), contingency_table[2, "WT"], 0)
    total_cases <- mutant_cases + wt_cases
    
    if (total_cases == 0) {
      mutant_percentage <- 0
      wt_percentage <- 0
    } else {
      mutant_percentage <- (mutant_cases / total_cases) * 100
      wt_percentage <- (wt_cases / total_cases) * 100
    }
    
    # Add a row to the results dataframe
    results_df <- rbind(results_df, data.frame(
      Subfamily = subfamily,
      Mutant_Percentage = mutant_percentage,
      WT_Percentage = wt_percentage,
      P_Value = fisher_test$p.value
    ))
  }
  
  return(results_df)
}

plot_subfamily_pie <- function(df, group, n) {
  # Filter the dataframe based on the group
  if (group == "Alu") {
    relevant_df <- df[grepl("^Alu|^FRAM", df$subfamily), ]
  } else if (group == "L1") {
    relevant_df <- df[grepl("^L1", df$subfamily), ]
  } else if (group == "SVA") {
    relevant_df <- df[grepl("^SVA", df$subfamily), ]
  } else {
    stop("Invalid group specified. Choose from 'Alu', 'L1', or 'SVA'.")
  }
  
  # Count occurrences of each subfamily
  subfamily_counts <- as.data.frame(table(relevant_df$subfamily))
  colnames(subfamily_counts) <- c("Subfamily", "Count")
  
  # Calculate percentages
  total_count <- sum(subfamily_counts$Count)
  subfamily_counts$Percentage <- (subfamily_counts$Count / total_count) * 100
  
  # Lump subfamilies with < n% into "Other"
  subfamily_counts$Subfamily <- ifelse(
    subfamily_counts$Percentage < n,
    "Other",
    as.character(subfamily_counts$Subfamily)
  )
  
  # Recalculate percentages after lumping
  lumped_counts <- aggregate(Count ~ Subfamily, data = subfamily_counts, sum)
  lumped_counts$Percentage <- (lumped_counts$Count / total_count) * 100
  
  # Compute positions for labels
  lumped_counts$Cumulative <- cumsum(lumped_counts$Count)
  lumped_counts$Midpoint <- lumped_counts$Cumulative - lumped_counts$Count / 2
  lumped_counts$Angle <- (lumped_counts$Midpoint / total_count) * 360
  lumped_counts$Label <- paste0(lumped_counts$Subfamily, " (", round(lumped_counts$Percentage, 1), "%)")
  
  # Create the pie chart with labels outside
  ggplot(lumped_counts, aes(x = 2, y = Count, fill = Subfamily)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar(theta = "y") +
    theme_void() +
    geom_text(
      aes(
        label = Label,
        x = 3, # Position outside the pie
        angle = ifelse(Angle > 90 & Angle < 270, Angle + 180, Angle) # Flip angle for better readability
      ),
      size = 3.5,
      hjust = 0
    ) +
    xlim(0.5, 3.5) # Add space for labels
}

# Function to perform Fisher test and plot proportions
fisher_test_and_plot <- function(group1_yes, group1_no, group2_yes, group2_no) {
  # Create a contingency table
  contingency_table <- matrix(c(group1_yes, group1_no, group2_yes, group2_no), 
                              nrow = 2, 
                              byrow = TRUE,
                              dimnames = list(Group = c("Group 1", "Group 2"),
                                              Event = c("Yes", "No")))
  
  # Print the contingency table
  print("Contingency Table:")
  print(contingency_table)
  
  # Perform Fisher's exact test
  fisher_result <- fisher.test(contingency_table)
  print("Fisher Test Result:")
  print(fisher_result)
  
  # Create a data frame for plotting
  total1 <- group1_yes + group1_no
  total2 <- group2_yes + group2_no
  
  proportions <- data.frame(
    Group = rep(c("Adult", "Pediatric"), each = 2),
    Event = rep(c("Yes", "No"), times = 2),
    Count = c(group1_yes, group1_no, group2_yes, group2_no)
  )
  
  proportions$Proportion <- proportions$Count / 
    c(rep(total1, 2), rep(total2, 2))
  
  # Plot proportions
  p <- ggplot(proportions, aes(x = Group, y = Proportion, fill = Event)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = c("No" = "#B1D586", "Yes" = "#0080A3"), 
                      labels = c("No" = "No repeats", "Yes" = "At least one repeat")) +
    labs(
      x = "Age group",
      y = "Proportion",
      fill = NULL
    ) 
  
  print(p)
}

filter_and_print_g_columns <- function(data) {
  # Identify columns starting with "g_"
  g_columns <- grep("^g_", colnames(data), value = TRUE)
  
  if (length(g_columns) == 0) {
    stop("No columns starting with 'g_' found in the dataset.")
  }
  
  # Check if g_columns contain numeric data
  if (!all(sapply(data[, g_columns], is.numeric))) {
    stop("Some 'g_' columns are not numeric. Ensure all 'g_' columns are numeric.")
  }
  
  # Check for non-zero values in any of the "g_" columns
  row_has_nonzero_g <- rowSums(data[, g_columns], na.rm = TRUE) > 0
  
  # Debugging: Print summary of g_ columns and row sums
  print("Column sums for g_ columns:")
  print(colSums(data[, g_columns], na.rm = TRUE))
  
  print("Row sums for g_ columns:")
  print(rowSums(data[, g_columns], na.rm = TRUE))
  
  # Filter rows where any g_ column is non-zero
  filtered_data <- data[row_has_nonzero_g, ]
  
  # Select specific columns: sample, total, tumour type, and the "g_" columns
  selected_columns <- c("sample", "total", "tumor_type", g_columns)
  selected_columns <- selected_columns[selected_columns %in% colnames(filtered_data)]
  
  # Print the resulting filtered and selected data
  if (nrow(filtered_data) == 0) {
    warning("No rows with non-zero 'g_' column values found.")
    return(NULL)
  }
  
  return(filtered_data[, selected_columns])
}

match_samples <- function(loh_time, te_aff_t) {
  
  # Create the base_sample column by removing everything after the first "_"
  loh_time <- loh_time %>%
    mutate(base_sample = sub("_.*$", "", sample))
  
  # Step 1: Exact match between loh_time$sample and te_aff_t$sample
  matched_samples <- loh_time %>%
    filter(sample %in% te_aff_t$sample)
  
  # Print exact matches
  print("Exact matches:")
  print(matched_samples %>% select(sample))
  
  # Step 2: For unmatched samples, find base_sample in te_aff_t$base_sample
  base_sample_matches <- loh_time %>%
    filter(!(sample %in% te_aff_t$sample)) %>%  # Filter unmatched samples
    inner_join(te_aff_t %>% mutate(base_sample = sub("_.*$", "", sample)), 
               by = "base_sample")  # Match on base_sample
  
  # Print base_sample matches
  print("Base sample matches:")
  print(base_sample_matches %>% select(sample.x, sample.y, base_sample))  # Adjust to the appropriate column names
  
  # Step 3: Print samples with no matches
  no_matches <- loh_time %>%
    filter(!(sample %in% te_aff_t$sample) & !(base_sample %in% te_aff_t$base_sample))
  
  print("Samples with no matches:")
  print(no_matches %>% select(sample, base_sample))  # Print the columns with no matches
}

fisher_test_unique_te <- function(data, min) {
  
  # Calculate total counts for unique TP53_WT and TP53_Mutant samples
  total_TP53_WT <- data %>%
    filter(TP53_status == "WT") %>%
    distinct(sample) %>%
    nrow()  # Count the number of unique TP53_WT samples
  
  total_TP53_Mutant <- data %>%
    filter(TP53_status == "Mutant") %>%
    distinct(sample) %>%
    nrow()  # Count the number of unique TP53_Mutant samples
  
  # Perform grouping and Fisher's test
  results <- data %>%
    group_by(SV_chrom, SV_start, SV_end, SV_length, ALT) %>%
    summarize(
      count_TP53_WT = sum(TP53_status == "WT"),        # Number of TP53_WT with the variant
      absent_TP53_WT = total_TP53_WT - count_TP53_WT,  # Calculate "absent" TP53_WT
      count_TP53_Mutant = sum(TP53_status == "Mutant"), # Number of TP53_Mutant with the variant
      absent_TP53_Mutant = total_TP53_Mutant - count_TP53_Mutant,  # Calculate "absent" TP53_Mutant
      perc_TP53_WT = count_TP53_WT/total_TP53_WT*100, # percentage of WT with TE 
      perc_TP53_Mutant = count_TP53_Mutant/total_TP53_Mutant*100, # percentage of mutant with TE
      fisher_p_value = {
        # Use the pre-calculated values in the contingency table
        contingency_table <- matrix(
          c(
            count_TP53_WT, 
            count_TP53_Mutant, 
            absent_TP53_WT, 
            absent_TP53_Mutant
          ),
          nrow = 2,
          byrow = TRUE
        )
        
        # Perform Fisher's Exact Test if the table is valid
        if (all(contingency_table >= min)) {
          fisher.test(contingency_table)$p.value
        } else {
          # Check if any 3/4 values in the table are > 3 and one is 0
          values <- as.vector(contingency_table)
          if (sum(values > 3) == 3 && sum(values == 0) == 1) {
            # Identify the index with 0 and the ones > 3
            zero_index <- which(values == 0)
            warning_indices <- which(values > 3)
            
            # Print a specific warning for this condition
            print(paste("Warning: 3/4 of the values in the table are >3, and the index with 0 is", zero_index))
            print("Details of all columns:")
            print(contingency_table)
          }
          
          NA  # Return NA if the test cannot be performed
        }
      },
      .groups = "drop"  # Ungroup after summarizing
    ) %>%
    filter(!is.na(fisher_p_value))  # Remove rows with NA p-values
  
  # Apply Benjamini-Hochberg correction to the raw p-values for multiple comparisons
  results$fisher_p_value_BH <- p.adjust(results$fisher_p_value, method = "BH")
  
  # Sort results in ascending order by adjusted p-value
  results <- results %>%
    arrange(fisher_p_value_BH)
  
  # Return the results with all calculated columns
  return(results)
}

fisher_test_by_tumor_type <- function(data, min_samples_tt = 3, min_samples_te=5) {
  
  # Step 1: Collapse rare tumor types into "Other"
  tumor_type_counts <- data %>%
    distinct(sample, tumor_type) %>%
    count(tumor_type, name = "n")
  
  common_tumors <- tumor_type_counts %>%
    filter(n >= min_samples_tt) %>%
    pull(tumor_type)
  
  data <- data %>%
    filter(!is.na(tumor_type)) %>%
    mutate(tumor_type_grouped = ifelse(tumor_type %in% common_tumors, tumor_type, "Other"))

  # Step 2: Group and run Fisher's test for each TE
  # Filter TEs that appear in fewer than min_samples_te unique samples
  data <- data %>%
    group_by(SV_chrom, SV_start, SV_end, SV_length, ALT) %>%
    filter(n_distinct(sample) >= min_samples_te) %>%
    ungroup()

  # Check if data is empty after filtering
  if (nrow(data) == 0) {
    cat("Warning: No TEs remain after filtering. Returning empty results.\n")
    return(data.frame(SV_chrom = character(0), SV_start = numeric(0), SV_end = numeric(0),
                      SV_length = numeric(0), ALT = character(0), fisher_p_value = numeric(0),
                      fisher_p_value_BH = numeric(0)))
  }

  # Create sample-tumor map for the entire dataset
  all_sample_tumor_map <- data %>%
    distinct(sample, tumor_type_grouped)

  results <- data %>%
    group_by(SV_chrom, SV_start, SV_end, SV_length, ALT) %>%
    do({
      this_data <- .

      # Determine which samples in this group have the current TE
      has_te_samples <- this_data %>% distinct(sample) %>% pull(sample)

      # Merge and create binary presence
      contingency_data <- all_sample_tumor_map %>%
        mutate(has_variant = ifelse(sample %in% has_te_samples, "Present", "Absent")) %>%
        count(tumor_type_grouped, has_variant) %>%
        tidyr::pivot_wider(names_from = has_variant, values_from = n, values_fill = 0)

      # Ensure both columns exist for Fisher test BEFORE converting to rownames
      if (!"Present" %in% colnames(contingency_data)) contingency_data$Present <- 0
      if (!"Absent" %in% colnames(contingency_data)) contingency_data$Absent <- 0

      # Only convert to rownames if we have valid data
      if (nrow(contingency_data) == 0 || any(is.na(contingency_data$tumor_type_grouped))) {
        return(data.frame(fisher_p_value = NA))
      }

      contingency_data <- contingency_data %>% column_to_rownames("tumor_type_grouped")

      # get p value
      fisher_p <- tryCatch({
        fisher.test(contingency_data)$p.value
      }, error = function(e) NA)

      # Get gene names and other annotations from the first row
      gene_names <- paste(unique(this_data$Gene_name), collapse = ";")
      location <- paste(unique(this_data$Location), collapse = ";")
      n_samples <- n_distinct(this_data$sample)

      data.frame(
        fisher_p_value = fisher_p,
        gene_names = gene_names,
        location = location,
        n_samples = n_samples
      )
    }) %>%
    ungroup() %>%
    filter(!is.na(fisher_p_value))

  # Step 3: Adjust p-values (only if there are results)
  if (nrow(results) > 0) {
    results$fisher_p_value_BH <- p.adjust(results$fisher_p_value, method = "BH")
  } else {
    cat("Warning: No valid Fisher test results. Returning empty dataframe.\n")
    return(data.frame(SV_chrom = character(0), SV_start = numeric(0), SV_end = numeric(0),
                      SV_length = numeric(0), ALT = character(0), fisher_p_value = numeric(0),
                      fisher_p_value_BH = numeric(0)))
  }
  
  # Step 4: Return sorted results
  results <- results %>% arrange(fisher_p_value_BH)
  
  return(results)
}


find_sig_te_samples <- function(sig_df, te_df) {
  # Define the columns you want to retain in the results
  selected_columns <- c("sample", "AnnotSV_ID", "SV_chrom", "SV_start", "SV_end", 
                        "SV_length", "ALT", "Gene_name", "Location", "Location2", 
                        "Overlapped_CDS_length", "Overlapped_CDS_percent", "Frameshift", 
                        "Dist_nearest_SS", "Nearest_SS_type", 
                        "Closest_left", "Closest_right", "Gene_count", "Exon_count", 
                        "RE_gene", "P_gain_phen", "P_gain_hpo", "P_gain_source", 
                        "P_gain_coord", "P_loss_phen", "P_loss_hpo", "P_loss_source", 
                        "P_loss_coord", "P_ins_phen", "P_ins_hpo", "P_ins_source", 
                        "P_ins_coord", "po_P_gain_phen", "po_P_gain_hpo", "po_P_gain_source", 
                        "po_P_gain_coord", "po_P_gain_percent", "po_P_loss_phen", 
                        "po_P_loss_hpo", "po_P_loss_source", "po_P_loss_coord", 
                        "po_P_loss_percent", "P_snvindel_nb", "P_snvindel_phen", "AN", 
                        "AC", "AF", "N_BI_GENOS", "N_HOMREF", "N_HET", "N_HOMALT", 
                        "FREQ_HOMREF", "FREQ_HET", "FREQ_HOMALT", "GRPMAX_AF", 
                        "GC_content_left", "GC_content_right", "Repeat_coord_left", 
                        "Repeat_type_left", "Repeat_coord_right", "Repeat_type_right", 
                        "ACMG", "HI", "TS", "DDD_HI_percent", "ExAC_delZ", "ExAC_dupZ", 
                        "ExAC_cnvZ", "ExAC_synZ", "ExAC_misZ", "OMIM_ID", "OMIM_phenotype", 
                        "OMIM_inheritance", "OMIM_morbid", "OMIM_morbid_candidate", 
                        "GnomAD_pLI", "ExAC_pLI", "AnnotSV_ranking_score", 
                        "AnnotSV_ranking_criteria", "subfamily", "tumor_type", "tumor_class", "age_at_diagnosis", 
                        "sex", "TP53_status", "Variant_location", "Variant_function", 
                        "inheritance", "protein.codon.change", "protein.codon.num", 
                        "Variant_Classification", "Variant_Type", "sample_centre", 
                        "sequence_centre", "cohort", "age_at_enrollment", "vital_status", 
                        "treatment", "lesion_type", "disease_state", "germline_source", 
                        "sample_topography", "hostseq_cancer", "first_aa", "second_aa", 
                        "first_aa_letter", "second_aa_letter", "mutation", "Domain", 
                        "cluster", "base_sample.y", "age_to_use")
  
  # Create an empty list to store dataframes for each row of sig_df
  result_list <- list()
  
  # Loop through each row of sig_df
  for (i in 1:nrow(sig_df)) {
    # Get the row from sig_df
    row <- sig_df[i, ]
    
    # Filter te_df for matching rows
    filtered_data <- te_df %>%
      filter(SV_chrom == row$SV_chrom & 
               SV_start == row$SV_start & 
               SV_end == row$SV_end & 
               SV_length == row$SV_length & 
               ALT == row$ALT) %>%
      select(all_of(selected_columns))  # Only keep selected columns
    
    # Store the filtered data in the result list
    result_list[[i]] <- filtered_data
  }
  
  # Return the list of dataframes
  return(result_list)
}

print_summary_sig_te_samples<- function(df) {
  # Loop through each unique TP53_status (WT and Mutant)
  for (tp53_status in unique(df$TP53_status)) {
    
    # Filter the dataframe for the current TP53_status
    filtered_df <- df %>% filter(TP53_status == tp53_status)
    
    # Print the sample and AnnotSV_ID for the filtered dataframe
    print(filtered_df %>% select(sample, AnnotSV_ID))
    
    # Print a table of tumor_type for the current TP53_status group
    cat("\nTable of tumor_type for TP53_status:", tp53_status, "\n")
    print(table(filtered_df$tumor_type))
    
    # Print a table of cohort for the current TP53_status group
    cat("\nTable of cohort for TP53_status:", tp53_status, "\n")
    print(table(filtered_df$cohort))
    
    # Print a table of sex for the current TP53_status group
    cat("\nTable of sex for TP53_status:", tp53_status, "\n")
    print(table(filtered_df$sex))
    
    # Print a summary of age_at_diagnosis for the current TP53_status group
    cat("\nSummary of age_at_diagnosis for TP53_status:", tp53_status, "\n")
    print(summary(filtered_df$age_at_diagnosis))
    
    cat("\n########\n")
  }
}

get_sig_te_contingency_tables <- function(data, sig_results, min_samples_tt) {
  # Group rare tumor types into "Other"
  tumor_type_counts <- data %>%
    distinct(sample, tumor_type) %>%
    count(tumor_type, name = "n")

  common_tumors <- tumor_type_counts %>%
    filter(n >= min_samples_tt) %>%
    pull(tumor_type)

  data <- data %>%
    filter(!is.na(tumor_type)) %>%
    mutate(tumor_type_grouped = ifelse(tumor_type %in% common_tumors, tumor_type, "Other"))

  # Build contingency table for each row (TE) in sig_results
  sig_results_with_tables <- sig_results %>%
    rowwise() %>%
    mutate(contingency_table = list({
      # Extract current row as a list
      te <- cur_data()

      # Filter matching rows in full dataset
      this_data <- data %>%
        filter(SV_chrom == te$SV_chrom,
               SV_start == te$SV_start,
               SV_end == te$SV_end,
               SV_length == te$SV_length,
               ALT == te$ALT)

      has_te_samples <- this_data %>%
        distinct(sample) %>%
        pull(sample)

      sample_tumor_map <- data %>%
        distinct(sample, tumor_type_grouped)

      contingency_data <- sample_tumor_map %>%
        filter(!is.na(tumor_type_grouped)) %>%
        mutate(has_variant = ifelse(sample %in% has_te_samples, "Present", "Absent")) %>%
        count(tumor_type_grouped, has_variant) %>%
        tidyr::pivot_wider(names_from = has_variant, values_from = n, values_fill = 0)

      # Ensure both columns exist BEFORE converting to rownames
      if (!"Present" %in% colnames(contingency_data)) contingency_data$Present <- 0
      if (!"Absent" %in% colnames(contingency_data)) contingency_data$Absent <- 0

      # Only convert to rownames if valid
      if (nrow(contingency_data) > 0 && !any(is.na(contingency_data$tumor_type_grouped))) {
        contingency_data <- contingency_data %>% column_to_rownames("tumor_type_grouped")
      } else {
        return(matrix(NA, nrow = 0, ncol = 0))
      }

      as.matrix(contingency_data)
    })) %>%
    ungroup()

  return(sig_results_with_tables)
}

# Convert contingency tables to long format suitable for CSV export
# Returns a data frame with one row per TE-tumor_type combination
get_sig_te_contingency_tables_long <- function(data, sig_results, min_samples_tt) {
  # Group rare tumor types into "Other"
  tumor_type_counts <- data %>%
    distinct(sample, tumor_type) %>%
    count(tumor_type, name = "n")

  common_tumors <- tumor_type_counts %>%
    filter(n >= min_samples_tt) %>%
    pull(tumor_type)

  data <- data %>%
    filter(!is.na(tumor_type)) %>%
    mutate(tumor_type_grouped = ifelse(tumor_type %in% common_tumors, tumor_type, "Other"))

  # Get all unique tumor types for consistency
  all_tumor_types <- unique(data$tumor_type_grouped)

  # For each significant TE, get counts by tumor type
  result_list <- list()

  for (i in 1:nrow(sig_results)) {
    te <- sig_results[i, ]

    # Filter matching rows in full dataset
    this_data <- data %>%
      filter(SV_chrom == te$SV_chrom,
             SV_start == te$SV_start,
             SV_end == te$SV_end,
             SV_length == te$SV_length,
             ALT == te$ALT)

    has_te_samples <- this_data %>%
      distinct(sample) %>%
      pull(sample)

    sample_tumor_map <- data %>%
      distinct(sample, tumor_type_grouped)

    # Create contingency data in long format
    contingency_long <- sample_tumor_map %>%
      filter(!is.na(tumor_type_grouped)) %>%
      mutate(has_variant = ifelse(sample %in% has_te_samples, "Present", "Absent")) %>%
      count(tumor_type_grouped, has_variant) %>%
      tidyr::pivot_wider(names_from = has_variant, values_from = n, values_fill = 0)

    # Ensure both columns exist
    if (!"Present" %in% colnames(contingency_long)) contingency_long$Present <- 0
    if (!"Absent" %in% colnames(contingency_long)) contingency_long$Absent <- 0

    # Add TE information and stats
    contingency_long <- contingency_long %>%
      mutate(
        SV_chrom = te$SV_chrom,
        SV_start = te$SV_start,
        SV_end = te$SV_end,
        SV_length = te$SV_length,
        ALT = te$ALT,
        fisher_p_value = te$fisher_p_value,
        fisher_p_value_BH = te$fisher_p_value_BH,
        total_samples = Present + Absent,
        percent_present = round(100 * Present / (Present + Absent), 2)
      ) %>%
      select(SV_chrom, SV_start, SV_end, SV_length, ALT,
             tumor_type = tumor_type_grouped,
             present = Present, absent = Absent, total_samples, percent_present,
             fisher_p_value, fisher_p_value_BH)

    result_list[[i]] <- contingency_long
  }

  # Combine all results
  result_df <- dplyr::bind_rows(result_list)

  return(result_df)
}


get_samples_with_sig_te <- function(data, sig_results) {
  sig_results_with_samples <- sig_results %>%
    rowwise() %>%
    mutate(samples_with_te = list({
      chr    <- SV_chrom
      start  <- SV_start
      end    <- SV_end
      length <- SV_length
      alt    <- ALT
      
      data %>%
        filter(SV_chrom == chr,
               SV_start == start,
               SV_end == end,
               SV_length == length,
               ALT == alt) %>%
        distinct(sample) %>%
        pull(sample)
    })) %>%
    ungroup()
  
  return(sig_results_with_samples)
}

find_sig_te_genes <- function(sig_df, te_df) {
  # Define the columns you want to retain in the results
  selected_columns <- c("sample", "AnnotSV_ID", "SV_chrom", "SV_start", "SV_end", 
                        "SV_length", "ALT", "Gene_name", "Location", "Location2", 
                        "Overlapped_CDS_length", "Overlapped_CDS_percent", "Frameshift", 
                        "Dist_nearest_SS", "Nearest_SS_type", 
                        "Closest_left", "Closest_right", "Gene_count", "Exon_count", 
                        "RE_gene", "P_gain_phen", "P_gain_hpo", "P_gain_source", 
                        "P_gain_coord", "P_loss_phen", "P_loss_hpo", "P_loss_source", 
                        "P_loss_coord", "P_ins_phen", "P_ins_hpo", "P_ins_source", 
                        "P_ins_coord", "po_P_gain_phen", "po_P_gain_hpo", "po_P_gain_source", 
                        "po_P_gain_coord", "po_P_gain_percent", "po_P_loss_phen", 
                        "po_P_loss_hpo", "po_P_loss_source", "po_P_loss_coord", 
                        "po_P_loss_percent", "P_snvindel_nb", "P_snvindel_phen", "AN", 
                        "AC", "AF", "N_BI_GENOS", "N_HOMREF", "N_HET", "N_HOMALT", 
                        "FREQ_HOMREF", "FREQ_HET", "FREQ_HOMALT", "GRPMAX_AF", 
                        "GC_content_left", "GC_content_right", "Repeat_coord_left", 
                        "Repeat_type_left", "Repeat_coord_right", "Repeat_type_right", 
                        "ACMG", "HI", "TS", "DDD_HI_percent", "ExAC_delZ", "ExAC_dupZ", 
                        "ExAC_cnvZ", "ExAC_synZ", "ExAC_misZ", "OMIM_ID", "OMIM_phenotype", 
                        "OMIM_inheritance", "OMIM_morbid", "OMIM_morbid_candidate", 
                        "GnomAD_pLI", "ExAC_pLI", "AnnotSV_ranking_score", 
                        "AnnotSV_ranking_criteria", "subfamily")
  
  # Create an empty list to store the results
  result_list <- list()
  
  # Loop through each row of sig_df
  for (i in 1:nrow(sig_df)) {
    # Get the row from sig_df
    row <- sig_df[i, ]
    
    # Filter te_df for matching rows and select only the selected columns
    filtered_data <- te_df %>%
      filter(SV_chrom == row$SV_chrom & 
               SV_start == row$SV_start & 
               SV_end == row$SV_end & 
               SV_length == row$SV_length & 
               ALT == row$ALT) %>%
      select(all_of(selected_columns)) %>%
      slice(1)  # Only get the first match (if there are multiple matches)
    
    # Append the filtered data to the result list
    if (nrow(filtered_data) > 0) {
      result_list[[i]] <- filtered_data
    }
  }
  
  # Combine all results into a single dataframe
  final_result <- bind_rows(result_list)
  
  # Return the combined dataframe
  return(final_result)
}

combine_gene_expression <- function(df1, df2) {
  # Merge the two data frames by 'gene_name'
  combined_df <- merge(df1, df2, by = "gene_name", all = TRUE)
  
  # Replace NA values with 0
  combined_df[is.na(combined_df)] <- 0
  
  return(combined_df)
  
}

rename_lfs_rna_columns <- function(rna_df, names_df) {
  # Create a named vector for renaming
  rename_map <- setNames(names_df$sample, names_df$rna_sample)
  
  # Rename columns in rna_df using the mapping
  colnames(rna_df) <- ifelse(colnames(rna_df) %in% names(rename_map), 
                             rename_map[colnames(rna_df)], 
                             colnames(rna_df))
  
  return(rna_df)
}

rename_stjude_rna <- function(df) {
  colnames(df)[-1] <- paste0(sub("_.*", "", colnames(df)[-1]), "_T")
  return(df)
}

rename_kics_rna<- function(rna_df, naming_df) {
  # Extract current column names from the RNA dataframe
  rna_names <- colnames(rna_df)
  
  # Separate gene_name column and other RNA columns
  gene_name_col <- rna_names[1]  # Assuming gene_name is the first column
  other_rna_cols <- rna_names[-1]  # This gets all columns EXCEPT the first one
  
  # Find which RNA names from the columns have mappings in the naming dataframe
  # Hardcoded to use "rna_name" column in naming_df
  matched_indices <- match(other_rna_cols, naming_df[["rna_name"]])
  
  # Identify which columns to keep (those with a match in naming_df)
  columns_to_keep <- !is.na(matched_indices)
  
  # Create a vector of column indices to keep (including gene_name)
  columns_to_keep_indices <- c(1, which(columns_to_keep) + 1)
  
  # Filter the dataframe to keep only columns with matches and the gene_name column
  filtered_df <- rna_df[, columns_to_keep_indices, drop = FALSE]
  
  # Get the new column names from the mapping
  # Hardcoded to use "sample" column in naming_df
  new_colnames <- c(gene_name_col, naming_df[["sample"]][matched_indices[columns_to_keep]])
  
  # Rename the columns
  colnames(filtered_df) <- new_colnames
  
  # Return the modified dataframe
  return(filtered_df)
}

filter_columns_by_sample <- function(main_df, reference_df) {
  # Get column names from main_df
  col_names <- colnames(main_df)
  
  # Check which column names are in the reference_df$sample
  cols_to_keep <- col_names %in% reference_df$sample
  
  # Return filtered dataframe
  return(main_df[, cols_to_keep, drop = FALSE])
}

split_re_gene_rows <- function(df, gene_col = "RE_gene", sep = ";") {
  df %>%
    separate_rows(!!sym(gene_col), sep = sep) %>%
    mutate(!!gene_col := str_trim(!!sym(gene_col)))
}

plot_gene_effects <- function(df, min_sample_tt = 5, remove_other = FALSE, top_n_genes = 10, SV_type = NULL, RE = FALSE, cancer_genes=FALSE) {
  # Filter by SV_type if provided
  if (!is.null(SV_type)) {
    df <- df %>% filter(ALT == SV_type)
  }

  # Add a new column `tumor_type_other` to group tumor types into "Other"
  df <- df %>%
    group_by(tumor_type) %>%
    mutate(tumor_type_other = ifelse(n_distinct(sample) < min_sample_tt, "Other", tumor_type)) %>%
    ungroup()
  
  # If RE is TRUE, preprocess RE_gene and use it instead of Gene_name
  if (RE) {
    df <- df %>%
      mutate(RE_gene = sub(" \\(.*$", "", RE_gene)) # Extract everything before " ("
  }
  
  # If cancer_genes is provided, filter by RE_gene_symbol
  if (!isFALSE(cancer_genes)) {
    df <- df %>% filter(RE_gene %in% cancer_genes)
  }
  
  # Determine the column to use for grouping genes
  gene_column <- if (RE) "RE_gene" else "Gene_name"
  
  # Count unique samples per gene and tumor type, filter top genes, and optionally remove "Other"
  gene_counts <- df %>%
    distinct(sample, .data[[gene_column]], tumor_type_other) %>%
    group_by(.data[[gene_column]], tumor_type_other) %>%
    summarise(sample_count = n_distinct(sample), .groups = "drop") %>%
    group_by(.data[[gene_column]]) %>%
    mutate(total_samples = sum(sample_count)) %>%
    ungroup() %>%
    arrange(desc(total_samples)) %>%
    filter(.data[[gene_column]] %in% head(unique(.data[[gene_column]]), top_n_genes)) %>%
    filter(!(remove_other & tumor_type_other == "Other"))
  
  # Rename the gene column for consistent plotting
  colnames(gene_counts)[1] <- "Gene"
  
  print(head(gene_counts))
  # Plot the stacked bar plot
  p <- ggplot(gene_counts, aes(x = reorder(Gene, -total_samples), y = sample_count, fill = tumor_type_other)) +
    geom_bar(stat = "identity") +
    labs(x = "Gene", y = "Number of samples affected", fill = "Tumor type") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(p)
}

# Filter false positive transposable element insertions using manual IGV review data
# - Removes all two-caller FPs listed in l1_merge_fp
# - For samples with Filtered_true="Y": keeps only TPs from unique_l1_tp
# - For samples with Filtered_true=blank: removes FPs from unique_l1_fp
# - Excludes samples not present in unique_l1_master (unreviewed)
filter_false_positives <- function(te_data, l1_merge_fp, unique_l1_tp, unique_l1_fp, unique_l1_master, check_missing_ids = TRUE) {
  
  # Use r_dir from global environment if available, otherwise current directory
  output_dir <- if(exists("r_dir")) r_dir else "."
  
  # l1_merge_fp should already be loaded as data.frame from fread
  cat("Filtering false positives using manual review data...\n")
  
  # Clear validation warnings file for fresh start
  if (exists("r_dir_files")) {
    validation_warnings_file <- file.path(r_dir_files, "validation_warnings.txt")
    if (file.exists(validation_warnings_file)) {
      file.remove(validation_warnings_file)
    }
  }
  
  # Create sample column if it doesn't exist (extract from ID)
  if (!"sample" %in% colnames(te_data)) {
    if ("ID" %in% colnames(te_data)) {
      # Extract sample from ID format: "0453_20-10584-A-02-00_T-1-56405386-1537-LINE1"
      # Pattern: everything up to and including "_T", then followed by a dash and numbers/letters
      te_data$sample <- sub("^(.+_T)-.*", "\\1", te_data$ID)
    } else {
      stop("Cannot create sample column - no ID column found")
    }
  }

  # Check for missing IDs BEFORE any filtering (only if requested)
  if (check_missing_ids) {
    # Only check TP IDs since FPs will be intentionally removed
    tp_ids <- unique_l1_tp$ID[!is.na(unique_l1_tp$ID)]
    missing_tp_ids <- tp_ids[!(tp_ids %in% te_data$ID)]

    # Filter out short IDs (length < 50) - these are likely truncated or malformed
    if (length(missing_tp_ids) > 0) {
      missing_tp_ids_long <- missing_tp_ids[nchar(missing_tp_ids) >= 50]
    } else {
      missing_tp_ids_long <- character(0)
    }

    # Write missing TP IDs to file if any found
    if (length(missing_tp_ids_long) > 0) {
      missing_file <- file.path(r_dir_files, "missing_ids_in_te_data.txt")

      # Add diagnostic info to the file
      diagnostic_info <- c(
        paste0("# Missing True Positive IDs (", length(missing_tp_ids_long), " total)"),
        paste0("# These IDs are in unique_l1_tp.csv but not found in the main te_data file"),
        paste0("# Total TPs in review file: ", length(tp_ids)),
        paste0("# Total IDs in te_data: ", nrow(te_data)),
        paste0("# Missing IDs (length >= 50):"),
        "",
        missing_tp_ids_long
      )
      writeLines(diagnostic_info, missing_file)
      cat("WARNING:", length(missing_tp_ids_long), "true positive IDs from review files not found in te_data. Written to:", missing_file, "\n")
    } else {
      cat("SUCCESS: All reviewed true positive IDs found in te_data\n")
    }
  }

  # Check for samples that haven't been reviewed yet
  te_samples <- unique(trimws(as.character(te_data$sample)))
  te_samples <- te_samples[!is.na(te_samples) & te_samples != ""]
  reviewed_samples <- unique(trimws(as.character(unique_l1_master$sample)))
  unreviewed_samples <- te_samples[!(te_samples %in% reviewed_samples)]
  
  cat("Sample review status: Found", length(te_samples), "TE samples,", length(reviewed_samples), "reviewed,", length(unreviewed_samples), "unreviewed\n")
  
  # Print and exclude unreviewed samples
  if (length(unreviewed_samples) > 0) {
    # Write unreviewed samples to file
    unreviewed_file <- file.path(r_dir_files, "unreviewed_samples_excluded.txt")
    writeLines(unreviewed_samples, unreviewed_file)
    cat("WARNING:", length(unreviewed_samples), "samples not reviewed and removed from analysis. Written to:", unreviewed_file, "\n")

    # Remove unreviewed samples from te_data
    te_data <- te_data[!(trimws(as.character(sample)) %in% unreviewed_samples)]
  } else {
    cat("SUCCESS: All samples have been reviewed - no unreviewed samples excluded\n")
  }
  
  # Calculate masks AFTER removing unreviewed samples
  # One-caller: contains "xtea" or "totalrecall" in ID
  # Two-caller: does not contain these terms
  one_caller_mask <- grepl("xtea|totalrecall", te_data$ID, ignore.case = TRUE)
  two_caller_mask <- !one_caller_mask

  cat("Caller types:", sum(one_caller_mask), "one-caller,", sum(two_caller_mask), "two-caller insertions\n")
  cat("Samples after unreviewed removal:", length(unique(te_data$sample)), "samples\n")
  
  remove_ids <- character(0)
  
  # 1. Filter ALL FP calls by two callers listed in l1_merge_fp
  two_caller_fps <- l1_merge_fp$ID[!is.na(l1_merge_fp$ID)]
  remove_ids <- c(remove_ids, as.character(two_caller_fps))
  cat("Two-caller FP removal:", sum(te_data$ID %in% two_caller_fps), "of", length(two_caller_fps), "FPs found in data\n")
  
  # 2. Filter one-caller calls based on sample-specific rules from unique_l1_master
  te_samples_clean <- trimws(as.character(te_data$sample))
  master_samples_clean <- trimws(as.character(unique_l1_master$sample))
  reviewed_samples_in_data <- unique(te_samples_clean[te_samples_clean %in% master_samples_clean])
  
  cat("One-caller sample-specific filtering for", length(reviewed_samples_in_data), "reviewed samples\n")
  
  for (sample_id in reviewed_samples_in_data) {
    # Get filtering rule for this sample from unique_l1_master
    sample_master_row <- unique_l1_master[trimws(as.character(unique_l1_master$sample)) == sample_id, ]
    
    # Handle case where sample might appear multiple times - take first occurrence
    if (nrow(sample_master_row) > 1) {
      sample_master_row <- sample_master_row[1, ]
    }
    
    # Get all one-caller TE IDs for this sample (don't touch two-caller calls)
    sample_one_caller_ids <- te_data$ID[trimws(as.character(te_data$sample)) == sample_id & one_caller_mask]
    
    if (nrow(sample_master_row) > 0) {
      filtered_true_status <- sample_master_row$Filtered_true
      
      # Additional validation logic based on number false/true columns
      if ("number_false" %in% colnames(sample_master_row) && "number_true" %in% colnames(sample_master_row) && "sum_unique" %in% colnames(sample_master_row)) {
        num_false <- as.numeric(sample_master_row$number_false)
        num_true <- as.numeric(sample_master_row$number_true) 
        sum_unique <- as.numeric(sample_master_row$sum_unique)
        
        # Skip filtering if no false positives
        if (!is.na(num_false) && num_false == 0) {
          #cat("Sample", sample_id, ": No false positives, skipping filtering\n")
          next
        }
        
        # Filter out all if all are false positives  
        if (!is.na(num_false) && !is.na(sum_unique) && num_false == sum_unique) {
          #cat("Sample", sample_id, ": All TEs are false positives, removing all\n")
          sample_ids_to_remove <- sample_one_caller_ids
          remove_ids <- c(remove_ids, as.character(sample_ids_to_remove))
          next
        }
        
        # Validate counts match review files
        if (!is.na(filtered_true_status) && trimws(as.character(filtered_true_status)) == "Y") {
          sample_tps_count <- sum(trimws(as.character(unique_l1_tp$sample)) == sample_id & !is.na(unique_l1_tp$ID))
          if (!is.na(num_true) && sample_tps_count != num_true) {
            warning_msg <- paste("Sample", sample_id, ": Mismatch in true positives. Master file:", num_true, "vs unique_l1_tp:", sample_tps_count)
            cat("WARNING:", warning_msg, "\n")
            warning_file <- file.path(r_dir_files, "validation_warnings.txt")
            write(warning_msg, warning_file, append = TRUE)
          }
        } else if (is.na(filtered_true_status) || trimws(as.character(filtered_true_status)) == "") {
          sample_fps_count <- sum(trimws(as.character(unique_l1_fp$sample)) == sample_id & !is.na(unique_l1_fp$ID))
          if (!is.na(num_false) && sample_fps_count != num_false) {
            warning_msg <- paste("Sample", sample_id, ": Mismatch in false positives. Master file:", num_false, "vs unique_l1_fp:", sample_fps_count)
            cat("WARNING:", warning_msg, "\n")
            warning_file <- file.path(r_dir_files, "validation_warnings.txt")
            write(warning_msg, warning_file, append = TRUE)
          }
        }
      }
      
      if (!is.na(filtered_true_status) && trimws(as.character(filtered_true_status)) == "Y") {
        # For Filtered_true = Y: Only keep TPs listed in unique_l1_tp for this sample
        sample_tps <- as.character(unique_l1_tp$ID[trimws(as.character(unique_l1_tp$sample)) == sample_id & !is.na(unique_l1_tp$ID)])
        
        # Remove all other one-caller calls from this sample (except TPs)
        sample_ids_to_remove <- sample_one_caller_ids[!(sample_one_caller_ids %in% sample_tps)]
        remove_ids <- c(remove_ids, as.character(sample_ids_to_remove))
        
      } else if (is.na(filtered_true_status) || trimws(as.character(filtered_true_status)) == "") {
        # For Filtered_true = blank: Remove FPs listed in unique_l1_fp for this sample
        sample_fps <- as.character(unique_l1_fp$ID[trimws(as.character(unique_l1_fp$sample)) == sample_id & !is.na(unique_l1_fp$ID)])
        
        # Remove only the specific one-caller FPs
        sample_fps_to_remove <- sample_one_caller_ids[sample_one_caller_ids %in% sample_fps]
        remove_ids <- c(remove_ids, as.character(sample_fps_to_remove))
        
      }
    }
  }
  
  # Debug: Check if remove_ids actually match anything in te_data
  if (length(remove_ids) > 0) {
    matching_ids <- sum(te_data$ID %in% remove_ids)
    cat("Final removal summary:\n")
    cat("- Unique IDs marked for removal:", length(unique(remove_ids)), "\n")
    cat("- IDs in te_data that match remove_ids:", matching_ids, "\n")
    
    if (matching_ids == 0) {
      cat("WARNING: No matching IDs found! Checking ID formats...\n")
      cat("- First few remove_ids:", paste(head(remove_ids, 3), collapse = ", "), "\n")
      cat("- First few te_data IDs:", paste(head(te_data$ID, 3), collapse = ", "), "\n")
      cat("- te_data ID class:", class(te_data$ID), "\n")
      cat("- remove_ids class:", class(remove_ids), "\n")
    }
  }
  
  # Filter te_data: remove all IDs in remove_ids
  te_filtered <- te_data[!(ID %in% remove_ids)]
  
  # Print detailed filtering summary by caller type
  if ("caller_cat" %in% colnames(te_data)) {
    # Original counts by caller type
    original_one_caller <- nrow(te_data[caller_cat == "one_caller"])
    original_two_caller <- nrow(te_data[caller_cat == "two_caller"])
    
    # Filtered counts by caller type
    filtered_one_caller <- nrow(te_filtered[caller_cat == "one_caller"])
    filtered_two_caller <- nrow(te_filtered[caller_cat == "two_caller"])
    
    cat("Detailed filtering summary by caller type:\n")
    cat("- One-caller: ", original_one_caller, " → ", filtered_one_caller, " (removed ", original_one_caller - filtered_one_caller, ", ", round(100*filtered_one_caller/original_one_caller, 1), "% retention)\n", sep = "")
    cat("- Two-caller: ", original_two_caller, " → ", filtered_two_caller, " (removed ", original_two_caller - filtered_two_caller, ", ", round(100*filtered_two_caller/original_two_caller, 1), "% retention)\n", sep = "")
  }
  
  # Print overall filtering summary
  cat("Filtering complete:", nrow(te_data), "→", nrow(te_filtered), "insertions (removed", nrow(te_data) - nrow(te_filtered), ")\n")
  
  # Safety check for empty dataset
  if (nrow(te_filtered) == 0) {
    cat("WARNING: All data was filtered out! This may cause downstream processing to fail.\n")
  }
  
  return(te_filtered)
}

test_all_te_expression_effects <- function(rna_df, te_split_df, min_samples = 5, fdr_cutoff = 0.05, group_by_gene = FALSE) {
  cat("Systematic TE expression analysis\n")
  cat("Minimum samples per TE:", min_samples, "\n")
  cat("FDR cutoff:", fdr_cutoff, "\n")

  # Check input types
  if (!is.data.frame(te_split_df)) {
    cat("Error: te_split_df is not a data frame. Class:", class(te_split_df), "\n")
    return(data.frame())
  }
  if (!is.data.frame(rna_df)) {
    cat("Error: rna_df is not a data frame. Class:", class(rna_df), "\n")
    return(data.frame())
  }

  # Clean sample names
  te_split_df$sample_clean <- sub("_.*", "", te_split_df$sample)

  # Identify unique TEs by genomic coordinates and type (optionally grouped by gene)
  if (group_by_gene) {
    cat("Identifying unique TE insertions grouped by gene...\n")
    te_sample_counts <- te_split_df %>%
      group_by(Gene_name, ALT) %>%
      summarise(
        num_samples = n_distinct(sample_clean),
        te_coordinates = paste(unique(paste(SV_chrom, SV_start, SV_end, sep=":")), collapse=";"),
        samples_with_te = paste(unique(sample_clean), collapse=";"),
        te_location = paste(unique(na.omit(Location)), collapse="; "),
        gene_features = paste(unique(na.omit(Location2)), collapse="; "),
        .groups = "drop"
      ) %>%
      filter(num_samples >= min_samples) %>%
      arrange(desc(num_samples))
  } else {
    cat("Identifying unique TE insertions...\n")
    te_sample_counts <- te_split_df %>%
      group_by(SV_chrom, SV_start, SV_end, SV_length, ALT) %>%
      summarise(
        num_samples = n_distinct(sample_clean),
        genes_affected = paste(unique(Gene_name), collapse=";"),
        samples_with_te = paste(unique(sample_clean), collapse=";"),
        te_location = paste(unique(na.omit(Location)), collapse="; "),
        gene_features = paste(unique(na.omit(Location2)), collapse="; "),
        .groups = "drop"
      ) %>%
      filter(num_samples >= min_samples) %>%
      arrange(desc(num_samples))
  }
  
  cat("TEs with ≥", min_samples, "samples:", nrow(te_sample_counts), "\n")
  
  if (nrow(te_sample_counts) == 0) {
    cat("No TEs meet minimum sample threshold. Returning empty results.\n")
    return(data.frame())
  }

  # Get all unique samples from RNA data (not just TE samples)
  # This includes samples with NO TEs as valid controls
  rna_samples <- colnames(rna_df)[-1]  # Exclude gene_name column
  rna_base_ids <- unique(sub("_.*", "", rna_samples))
  cat("Total samples available in RNA data:", length(rna_base_ids), "\n")
  cat("Samples with gene-associated TEs:", length(unique(te_split_df$sample_clean)), "\n")
  
  # Prepare results dataframe
  results <- data.frame()
  
  cat("Testing TE effects on gene expression...\n")
  
  for (i in 1:nrow(te_sample_counts)) {
    te_info <- te_sample_counts[i, ]

    samples_with_te <- unlist(strsplit(te_info$samples_with_te, ";"))
    # Use ALL RNA samples as potential controls, not just those with gene-associated TEs
    samples_without_te <- setdiff(rna_base_ids, samples_with_te)
    
    if (group_by_gene) {
      # When grouping by gene, each row is a gene-TE type combination
      genes_affected <- te_info$Gene_name
      # Suppressed verbose TE iteration output
    } else {
      # Original mode: split genes affected by this TE
      genes_affected <- unlist(strsplit(te_info$genes_affected, ";"))
      # Suppressed verbose TE iteration output
    }
    
    # Test each gene affected by this TE
    for (gene in genes_affected) {
      gene <- trimws(gene)  # Remove whitespace
      
      # Check if gene exists in RNA data
      first_col_name <- names(rna_df)[1]
      if (!gene %in% rna_df[[first_col_name]]) {
        # Suppressed verbose gene lookup messages
        next
      }
      
      # Get expression data for this gene
      gene_expr_row <- rna_df[rna_df[[first_col_name]] == gene, ]
      
      # Match RNA samples with TE samples
      rna_samples <- colnames(rna_df)[-1]  # Exclude gene_name column
      
      # Extract base IDs from RNA sample names for matching with TE data
      rna_base_ids <- sub("_.*", "", rna_samples)
      
      # Find RNA samples that match TE samples (by base ID)
      matched_with_te_indices <- which(rna_base_ids %in% samples_with_te)
      matched_without_te_indices <- which(rna_base_ids %in% samples_without_te)
      
      matched_with_te <- rna_samples[matched_with_te_indices]
      matched_without_te <- rna_samples[matched_without_te_indices]
      
      
      if (length(matched_with_te) < min_samples || length(matched_without_te) < min_samples) {
        # Suppressed insufficient samples message
        next
      }
      
      # Extract expression values
      expr_with_te <- as.numeric(gene_expr_row[, matched_with_te])
      expr_without_te <- as.numeric(gene_expr_row[, matched_without_te])
      
      # Remove any NA values
      expr_with_te <- expr_with_te[!is.na(expr_with_te)]
      expr_without_te <- expr_without_te[!is.na(expr_without_te)]
      
      if (length(expr_with_te) < min_samples || length(expr_without_te) < min_samples) {
        # Suppressed insufficient non-NA samples message
        next
      }
      
      # Perform Wilcoxon test
      tryCatch({
        wilcox_result <- wilcox.test(expr_with_te, expr_without_te)

        # Calculate medians and means
        median_with_te <- median(expr_with_te)
        median_without_te <- median(expr_without_te)
        mean_with_te <- mean(expr_with_te)
        mean_without_te <- mean(expr_without_te)
        fold_change <- median_with_te / (median_without_te + 0.001)  # Add small value to avoid division by zero

        # Determine effect direction
        effect_direction <- ifelse(median_with_te > median_without_te, "Up", "Down")

        # Determine TE location type
        te_location_str <- te_info$te_location
        te_location_type <- case_when(
          grepl("exon", te_location_str, ignore.case = TRUE) & grepl("intron", te_location_str, ignore.case = TRUE) ~ "Both exonic & intronic",
          grepl("exon", te_location_str, ignore.case = TRUE) ~ "Exonic",
          grepl("intron", te_location_str, ignore.case = TRUE) ~ "Intronic",
          TRUE ~ "Other"
        )

        # Store results
        if (group_by_gene) {
          result_row <- data.frame(
            gene = gene,
            te_type = te_info$ALT,
            te_coordinates = te_info$te_coordinates,
            te_location_type = te_location_type,
            te_location = te_info$te_location,
            gene_features = te_info$gene_features,
            samples_with_te = length(matched_with_te),
            samples_without_te = length(matched_without_te),
            median_with_te = median_with_te,
            median_without_te = median_without_te,
            mean_with_te = mean_with_te,
            mean_without_te = mean_without_te,
            fold_change = fold_change,
            effect_direction = effect_direction,
            p_value = wilcox_result$p.value,
            stringsAsFactors = FALSE
          )
        } else {
          result_row <- data.frame(
            te_chrom = te_info$SV_chrom,
            te_start = te_info$SV_start,
            te_end = te_info$SV_end,
            te_length = te_info$SV_length,
            te_type = te_info$ALT,
            gene = gene,
            te_location_type = te_location_type,
            te_location = te_info$te_location,
            gene_features = te_info$gene_features,
            samples_with_te = length(matched_with_te),
            samples_without_te = length(matched_without_te),
            median_with_te = median_with_te,
            median_without_te = median_without_te,
            mean_with_te = mean_with_te,
            mean_without_te = mean_without_te,
            fold_change = fold_change,
            effect_direction = effect_direction,
            p_value = wilcox_result$p.value,
            stringsAsFactors = FALSE
          )
        }

        results <- rbind(results, result_row)

        # Suppressed verbose gene test output

      }, error = function(e) {
        # Suppressed verbose error output
      })
    }
  }
  
  if (nrow(results) == 0) {
    cat("No successful tests completed.\n")
    return(data.frame())
  }
  
  # Apply FDR correction
  cat("Applying FDR correction...\n")
  results$p_adj <- p.adjust(results$p_value, method = "fdr")
  
  # Sort by adjusted p-value
  results <- results[order(results$p_adj), ]
  
  # Identify significant results
  significant_results <- results[results$p_adj < fdr_cutoff, ]
  
  cat("Results summary:\n")
  cat("Total tests performed:", nrow(results), "\n")
  cat("Significant results (FDR <", fdr_cutoff, "):", nrow(significant_results), "\n")
  
  if (nrow(significant_results) > 0) {
    cat("\nSignificant TE-gene expression associations:\n")
    for (i in 1:min(10, nrow(significant_results))) {
      row <- significant_results[i, ]
      cat(sprintf("  %s (chr%d:%d-%d) → %s: FC=%.2f, p_adj=%.2e\n", 
                  row$te_type, row$te_chrom, row$te_start, row$te_end,
                  row$gene, row$fold_change, row$p_adj))
    }
  }
  
  return(results)
}

# Helper function to extract full sample name from TE ID column
extract_full_sample_from_te_id <- function(te_id) {
  # For TE IDs like "0197_20-15931-A-02-00_T-totalrecall-1-245334-844-LINE1"
  # Extract the sample part: "0197_20-15931-A-02-00_T"
  # Split by "-" and take parts before the method name
  
  if (is.na(te_id) || te_id == "") return(te_id)
  
  # Look for pattern: sample_ID followed by method (totalrecall, xtea, etc.)
  # Pattern: everything before "-totalrecall-" or "-xtea-" or similar method indicators
  sample_part <- gsub("-(totalrecall|xtea|transurveyor)-.*$", "", te_id)
  
  return(sample_part)
}

# Helper function to extract base sample ID from sample names (keep for compatibility)
extract_te_sample_id <- function(te_sample) {
  # For samples like "0198", "3425", "3872", return as is
  # For samples with suffixes like "0002_321321_T", extract the base ID "0002"
  base_id <- gsub("_.*$", "", te_sample)  # Remove everything after first underscore
  return(base_id)
}

# Wrapper function for systematic TE expression analysis with plotting
run_systematic_te_expression_analysis <- function(rna_filtered, te_split_df, plot_dir, min_samples = 5, fdr_cutoff = 0.05, create_plots = TRUE, max_plots = 3, group_by_gene = FALSE, files_dir = NULL) {
  cat("=== Systematic TE Expression Analysis ===\n")

  # Check inputs
  if (!is.data.frame(te_split_df)) {
    cat("ERROR: te_split_df is not a data frame!\n")
    cat("  Class:", class(te_split_df), "\n")
    cat("  Type:", typeof(te_split_df), "\n")
    return(data.frame())
  }
  if (!is.data.frame(rna_filtered)) {
    cat("ERROR: rna_filtered is not a data frame!\n")
    cat("  Class:", class(rna_filtered), "\n")
    return(data.frame())
  }

  cat("Parameters:\n")
  cat("  - Minimum samples per TE:", min_samples, "\n")
  cat("  - FDR cutoff:", fdr_cutoff, "\n")
  cat("  - Group by gene:", group_by_gene, "\n")
  cat("  - RNA samples available:", ncol(rna_filtered) - 1, "\n")
  cat("  - TE samples available:", length(unique(te_split_df$sample)), "\n")
  cat("  - Genes in RNA data:", nrow(rna_filtered), "\n")
  cat("\n")

  systematic_results <- test_all_te_expression_effects(
    rna_df = rna_filtered,
    te_split_df = te_split_df,
    min_samples = min_samples,
    fdr_cutoff = fdr_cutoff,
    group_by_gene = group_by_gene
  )

  if (nrow(systematic_results) > 0) {
    # Determine output directory (use files_dir if provided, otherwise plot_dir)
    output_dir <- ifelse(!is.null(files_dir), files_dir, plot_dir)

    # Save complete results
    mode_suffix <- ifelse(group_by_gene, "_by_gene", "_by_coordinates")
    results_file <- paste0(output_dir, "rna_results_", min_samples, "samples", mode_suffix, ".csv")
    write.csv(systematic_results, results_file, row.names=FALSE)
    cat("Complete results (", nrow(systematic_results), "tests) saved to", basename(results_file), "\n")

    # Display significant results
    significant <- systematic_results[systematic_results$p_adj < fdr_cutoff, ]
    if (nrow(significant) > 0) {
      cat("Found", nrow(significant), "significant TE-gene expression associations at FDR <", fdr_cutoff, "!\n")

      # Save significant results to separate file
      sig_file <- paste0(output_dir, "significant_te_expression_results_", min_samples, "samples_fdr", fdr_cutoff, mode_suffix, ".csv")
      write.csv(significant, sig_file, row.names=FALSE)
      cat("Significant results saved to", basename(sig_file), "\n")
      
      # Print top 10 significant results
      cat("\nTop", min(10, nrow(significant)), "significant TE-gene associations:\n")
      for (i in 1:min(10, nrow(significant))) {
        row <- significant[i, ]
        if (group_by_gene) {
          cat(sprintf("  %d. %s + %s: FC=%.2f, p_adj=%.2e (n_with=%d, n_without=%d)\n", 
                      i, row$gene, row$te_type, row$fold_change, row$p_adj, 
                      row$samples_with_te, row$samples_without_te))
        } else {
          cat(sprintf("  %d. %s (chr%d:%d-%d) → %s: FC=%.2f, p_adj=%.2e (n_with=%d, n_without=%d)\n", 
                      i, row$te_type, row$te_chrom, row$te_start, row$te_end,
                      row$gene, row$fold_change, row$p_adj, row$samples_with_te, row$samples_without_te))
        }
      }
      
      # Create plots for top significant results if requested
      if (create_plots) {
        cat("\nCreating plots for top", min(max_plots, nrow(significant)), "significant associations...\n")
        for (i in 1:min(max_plots, nrow(significant))) {
          result_row <- significant[i, ]
          tryCatch({
            cat("  Creating plot", i, "for", result_row$gene, "TE association...\n")
            plot_result <- plot_gene_expression_te(
              rna_df = rna_filtered,
              te_split_df = te_split_df,
              gene_of_interest = result_row$gene,
              group_column = "TP53_status",
              x_lab = "TE Status",
              y_lab = paste(result_row$gene, "Expression (FPKM)"),
              log_scale = FALSE,
              te_based_grouping = TRUE
            )
            plot_title <- paste("Significant:", result_row$gene, "by TE presence")
            titled_print(plot_result$plot, plot_title)
            
            plot_filename <- paste0(plot_dir, "significant_", result_row$gene, "_te_effect.png")
            ggsave(plot_filename, plot=plot_result$plot, width=8, height=6)
            cat("    Plot saved to", basename(plot_filename), "\n")
          }, error = function(e) {
            cat("    Error creating plot for", result_row$gene, ":", e$message, "\n")
          })
        }
      }
    } else {
      cat("No significant TE-gene expression associations found at FDR <", fdr_cutoff, "\n")
      cat("Summary of all results:\n")
      if (nrow(systematic_results) > 0) {
        cat("  - Total tests performed:", nrow(systematic_results), "\n")
        cat("  - Minimum p-value:", min(systematic_results$p_value), "\n")
        cat("  - Minimum adjusted p-value:", min(systematic_results$p_adj), "\n")
        cat("  - Tests with p_adj < 0.1:", sum(systematic_results$p_adj < 0.1), "\n")
        cat("  - Tests with p_adj < 0.2:", sum(systematic_results$p_adj < 0.2), "\n")
      }
    }
  } else {
    cat("No TE-gene pairs could be tested.\n")
    cat("Possible reasons:\n")
    cat("  - No TEs present in >=", min_samples, "samples\n")
    cat("  - No overlap between RNA and TE samples\n")
    cat("  - No genes affected by TEs found in RNA data\n")
  }
  
  cat("\n")
  return(systematic_results)
}

# Function to clean RNA samples for tumor data (keep RNA format)
rename_rna_samples_kics_tumor <- function(kics_rna_data, matched_dna_rna_file = "/Users/briannelaverty/Documents/R_Malkin/te/data/rna/matched_dna_rna.csv") {
  cat("Processing KICS tumor RNA samples using matched_dna_rna.csv...\n")
  
  # Load the DNA-RNA mapping file
  dna_rna_mapping <- read.csv(matched_dna_rna_file, stringsAsFactors = FALSE)
  cat("Loaded DNA-RNA mapping with", nrow(dna_rna_mapping), "entries\n")
  
  # Get current RNA sample names (excluding gene_name column)
  rna_sample_names <- colnames(kics_rna_data)[-1]
  
  # Remove "X" prefix from RNA sample names if present (R adds this to numeric column names)
  rna_sample_names_clean <- gsub("^X", "", rna_sample_names)
  
  # Create mapping dictionary: RNA name -> DNA name
  mapping_dict <- setNames(dna_rna_mapping$sample_name_with_t, dna_rna_mapping$kics_rna_name)
  
  # Apply mapping to rename RNA samples to DNA sample names
  new_names <- sapply(rna_sample_names_clean, function(rna_name) {
    if (rna_name %in% names(mapping_dict)) {
      return(mapping_dict[rna_name])
    } else {
      return(rna_name)  # Keep original if not found
    }
  })
  
  # Update column names with DNA names
  colnames(kics_rna_data)[-1] <- new_names
  
  # Count successful mappings
  mapped_count <- sum(new_names != rna_sample_names_clean)
  
  cat("KICS tumor RNA samples processed. Sample count:", length(new_names), "\n")
  cat("Successfully mapped", mapped_count, "RNA samples to DNA names\n")
  cat("Example mapped samples (first 5):", paste(head(new_names, 5), collapse=", "), "\n")
  
  return(kics_rna_data)
}

# Function to rename RNA samples to match TE sample IDs for germline data
rename_rna_samples_kics_germline <- function(kics_rna_data, kics_dna2rna_mapping) {
  cat("Renaming KICS germline RNA samples to match TE sample IDs...\n")
  
  # Create a reverse mapping from RNA name to FULL WGS sample name
  reverse_mapping <- setNames(kics_dna2rna_mapping$sample, kics_dna2rna_mapping$rna_name)
  
  # Get current RNA sample names (excluding gene_name column)
  rna_sample_names <- colnames(kics_rna_data)[-1]
  
  # Remove "X" prefix from RNA sample names if present
  rna_sample_names_clean <- gsub("^X", "", rna_sample_names)
  
  # Apply the reverse mapping to get full sample names
  new_names <- sapply(rna_sample_names_clean, function(rna_name) {
    if (rna_name %in% names(reverse_mapping)) {
      return(reverse_mapping[rna_name])
    } else {
      # If not found in mapping, keep the original cleaned name
      return(rna_name)
    }
  })
  
  # Update column names with full sample names
  colnames(kics_rna_data)[-1] <- new_names
  
  cat("KICS germline RNA samples renamed. Sample count:", length(new_names), "\n")
  cat("Unique samples after renaming:", length(unique(new_names)), "\n")
  cat("Example mappings (first 5):\n")
  for (i in 1:min(5, length(rna_sample_names_clean))) {
    cat("  ", rna_sample_names_clean[i], "->", new_names[i], "\n")
  }
  
  return(kics_rna_data)
}

# Function to rename LFS RNA samples to match TE sample IDs  
rename_rna_samples_lfs <- function(lfs_rna_data, lfs_wgs2rna_mapping) {
  cat("Renaming LFS RNA samples to match TE sample IDs...\n")
  
  # Create mapping from RNA sample to WGS sample
  mapping_dict <- setNames(lfs_wgs2rna_mapping$sample, lfs_wgs2rna_mapping$rna_sample)
  
  # Get current RNA sample names (excluding gene_name column)
  rna_sample_names <- colnames(lfs_rna_data)[-1]
  
  # Apply mapping
  new_names <- sapply(rna_sample_names, function(rna_name) {
    if (rna_name %in% names(mapping_dict)) {
      return(mapping_dict[rna_name])
    } else {
      return(rna_name)  # Keep original if not found
    }
  })
  
  # Update column names
  colnames(lfs_rna_data)[-1] <- new_names
  
  cat("LFS RNA samples renamed. Sample count:", length(new_names), "\n")
  return(lfs_rna_data)
}

# Function to rename St. Jude RNA samples (keep full sample names)
rename_rna_samples_stjude <- function(stjude_rna_data) {
  cat("Processing St. Jude RNA samples (removing duplicates)...\n")
  
  # Get St. Jude sample names
  rna_sample_names <- colnames(stjude_rna_data)[-1]
  
  # Remove "X" prefix if R added it to numeric column names
  clean_names <- gsub("^X", "", rna_sample_names)
  
  # Remove duplicate parts (e.g., SJACT071_SJACT071 -> SJACT071_T)
  new_names <- sapply(clean_names, function(name) {
    parts <- strsplit(name, "_")[[1]]
    if (length(parts) >= 2 && parts[1] == parts[2]) {
      # If first two parts are identical, use first part + _T
      return(paste0(parts[1], "_T"))
    } else {
      return(name)  # Keep original if no duplicate pattern
    }
  })
  
  colnames(stjude_rna_data)[-1] <- new_names
  
  cat("St. Jude RNA samples processed. Sample count:", length(new_names), "\n")
  cat("Example samples (first 5):", paste(head(new_names, 5), collapse=", "), "\n")
  return(stjude_rna_data)
}

# Function to plot gene expression vs TE status
plot_gene_expression_te <- function(rna_df, te_split_df, gene_of_interest, group_column,
                                   x_lab, y_lab, log_scale = FALSE,
                                   te_based_grouping = FALSE, combined_grouping = FALSE) {

  # Validate inputs (duplicate function - add same checks)
  if (!is.data.frame(te_split_df)) {
    cat("ERROR in plot_gene_expression_te (v2): te_split_df is not a data frame!\n")
    cat("  Class:", class(te_split_df), "\n")
    cat("  Type:", typeof(te_split_df), "\n")
    return(list(plot = NULL, merged_data = data.frame()))
  }
  if (!is.data.frame(rna_df)) {
    cat("ERROR in plot_gene_expression_te (v2): rna_df is not a data frame!\n")
    cat("  Class:", class(rna_df), "\n")
    return(list(plot = NULL, merged_data = data.frame()))
  }

  # Filter for the gene of interest
  if (!gene_of_interest %in% rna_df$gene_name) {
    stop(paste("Gene", gene_of_interest, "not found in RNA data"))
  }
  
  gene_expr <- rna_df[rna_df$gene_name == gene_of_interest, ]
  
  # Reshape to long format
  expr_long <- gene_expr %>%
    select(-gene_name) %>%
    pivot_longer(everything(), names_to = "Sample", values_to = "Expression") %>%
    filter(!is.na(Expression), Expression > 0)  # Remove missing/zero values
  
  cat("Gene expression data for", gene_of_interest, "- samples with data:", nrow(expr_long), "\n")
  
  if (te_based_grouping) {
    # Use TE data for grouping
    te_gene_data <- te_split_df %>%
      filter(Gene_name == gene_of_interest) %>%
      select(sample, any_of(group_column)) %>%
      rename(Sample = sample)  # Rename to match RNA data format
    
    if (nrow(te_gene_data) == 0) {
      cat("No TE data found for gene", gene_of_interest, ". Creating all samples as 'No TE'.\n")
      # Create a dataframe with all RNA samples marked as "No TE"
      te_gene_data <- expr_long %>%
        select(Sample) %>%
        distinct() %>%
        mutate(TE_Status = "No TE")
    } else {
      te_gene_data <- te_gene_data %>%
        mutate(TE_Status = "Has TE")
      
      # Add samples without TEs
      samples_without_te <- expr_long %>%
        filter(!Sample %in% te_gene_data$Sample) %>%
        select(Sample) %>%
        distinct() %>%
        mutate(TE_Status = "No TE")
      
      # Add clinical data if available
      if (group_column %in% colnames(te_gene_data)) {
        samples_without_te[[group_column]] <- NA
      }
      
      te_gene_data <- bind_rows(te_gene_data, samples_without_te)
    }
    
    # Merge with expression data
    plot_data <- expr_long %>%
      left_join(te_gene_data, by = "Sample")
    
    # Set up grouping variable
    if (combined_grouping && group_column %in% colnames(plot_data)) {
      plot_data <- plot_data %>%
        mutate(Group = paste(TE_Status, get(group_column), sep = " & "))
      x_var <- "Group"
    } else {
      x_var <- "TE_Status"
    }
    
  } else {
    # Use clinical data for grouping
    clinical_data <- te_split_df %>%
      select(sample, any_of(group_column)) %>%
      distinct() %>%
      rename(Sample = sample)  # Rename to match RNA data format
    
    plot_data <- expr_long %>%
      left_join(clinical_data, by = "Sample")
    
    x_var <- group_column
  }
  
  # Apply log transformation if requested
  if (log_scale) {
    plot_data <- plot_data %>%
      mutate(Expression = log2(Expression + 1))
    y_lab <- paste("log2(", y_lab, " + 1)")
  }
  
  # Create the plot
  p <- ggplot(plot_data, aes_string(x = x_var, y = "Expression")) +
    geom_boxplot() +
    geom_jitter(width = 0.2, alpha = 0.6) +
    labs(x = x_lab, y = y_lab, title = paste("Expression of", gene_of_interest)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Add statistical test if there are exactly 2 groups
  if (length(unique(plot_data[[x_var]])) == 2) {
    p <- p + geom_signif(comparisons = list(unique(plot_data[[x_var]])), 
                        map_signif_level = TRUE, test = "wilcox.test")
  }
  
  # Print summary statistics
  summary_stats <- plot_data %>%
    group_by(!!sym(x_var)) %>%
    summarise(
      n = n(),
      median = median(Expression, na.rm = TRUE),
      mean = mean(Expression, na.rm = TRUE),
      .groups = 'drop'
    )
  
  cat("Summary statistics for", gene_of_interest, "expression:\n")
  print(summary_stats)
  
  return(list(plot = p, data = plot_data, summary = summary_stats))
}

validate_rna_mapping <- function(kics_rna, lfs_rna, kics_dna2rna, lfs_wgs2rna) {
  cat("Checking KICS mapping file for RNA sample availability...\n")
  
  # Check KICS mapping entries that don't have corresponding RNA data
  kics_rna_samples <- colnames(kics_rna)[-1]  # Exclude gene_name column
  mapping_col <- NULL
  if ("rna_name" %in% colnames(kics_dna2rna)) {
    mapping_col <- "rna_name"
  } else if ("rna" %in% colnames(kics_dna2rna)) {
    mapping_col <- "rna"
  }
  
  if (!is.null(mapping_col)) {
    # Get valid RNA IDs from mapping file (remove NA and empty values)
    valid_mapping_rows <- kics_dna2rna[!is.na(kics_dna2rna[[mapping_col]]) & kics_dna2rna[[mapping_col]] != "", ]
    
    if (nrow(valid_mapping_rows) > 0) {
      # Check which mapping entries don't have corresponding RNA samples
      unmatched_mapping_rows <- c()
      for (i in 1:nrow(valid_mapping_rows)) {
        rna_id <- valid_mapping_rows[[mapping_col]][i]
        # Look for RNA sample ending with this RNA ID
        matching_rna <- kics_rna_samples[grepl(paste0("_", rna_id, "$"), kics_rna_samples)]
        if (length(matching_rna) == 0) {
          unmatched_mapping_rows <- c(unmatched_mapping_rows, i)
        }
      }
      
      if (length(unmatched_mapping_rows) > 0) {
        cat("KICS mapping rows without corresponding RNA samples:\n")
        for (i in head(unmatched_mapping_rows, 10)) {
          row <- valid_mapping_rows[i, ]
          cat(sprintf("  Sample %s (DNA: %s, RNA: %s) - no RNA sample ending with _%s\n", 
                      row$sample, row$dna, row[[mapping_col]], row[[mapping_col]]))
        }
        if (length(unmatched_mapping_rows) > 10) {
          cat("  ... and", length(unmatched_mapping_rows) - 10, "more mapping rows\n")
        }
      } else {
        cat("All KICS mapping entries have corresponding RNA samples.\n")
      }
      
      cat("KICS mapping summary:", nrow(valid_mapping_rows), "total entries,", 
          length(unmatched_mapping_rows), "without RNA matches\n")
    } else {
      cat("No valid KICS mapping entries found.\n")
    }
  } else {
    cat("Warning: Neither 'rna_name' nor 'rna' column found in kics_dna2rna mapping file.\n")
    cat("Available columns:", paste(colnames(kics_dna2rna), collapse=", "), "\n")
  }
  
  # LFS and St. Jude have separate processing - just note they're handled differently
  cat("LFS and St. Jude samples processed separately with their own mapping logic.\n")
  
  invisible(NULL)
}

# ====== TUMOR RNA PROCESSING FUNCTIONS ======

#' Create unified RNA dataframe structure - only include samples with actual RNA data
#' @param rna_df Input RNA dataframe
#' @param all_samples Vector of all sample names (UNUSED - kept for compatibility)
#' @return RNA dataframe with only actual RNA samples (no zero-filled columns)
create_unified_rna <- function(rna_df, all_samples = NULL) {
  # Only return columns that actually exist (no zero-filled columns)
  # Keep gene_name, all existing sample columns, and dataset
  existing_sample_cols <- setdiff(colnames(rna_df), c("gene_name", "dataset"))
  rna_df <- rna_df[, c("gene_name", existing_sample_cols, "dataset")]
  return(rna_df)
}

#' Process and combine tumor RNA datasets
#' @param kics_rna KICS RNA FPKM data
#' @param lfs_rna LFS RNA FPKM data  
#' @param stjude_rna St. Jude RNA FPKM data
#' @param lfs_wgs2rna LFS WGS to RNA mapping file
#' @return List with combined RNA data and sample information
process_tumor_rna_data <- function(kics_rna, lfs_rna, stjude_rna, lfs_wgs2rna) {
  cat("Processing tumor RNA data with corrected sample matching...\n")
  
  # Validate RNA mapping before processing (KICS now uses matched_dna_rna.csv)
  cat("KICS RNA processing now uses matched_dna_rna.csv for direct DNA-RNA mapping.\n")
  
  # Rename RNA samples to match TE sample IDs (using matched_dna_rna.csv)
  kics_rna_renamed <- rename_rna_samples_kics_tumor(kics_rna)
  lfs_rna_renamed <- rename_rna_samples_lfs(lfs_rna, lfs_wgs2rna)
  stjude_rna_renamed <- rename_rna_samples_stjude(stjude_rna)
  
  # Add dataset identifiers
  kics_rna_renamed$dataset <- "KICS"
  lfs_rna_renamed$dataset <- "LFS"
  stjude_rna_renamed$dataset <- "StJude"
  
  # Find common genes across all datasets
  common_genes <- intersect(intersect(kics_rna_renamed$gene_name, lfs_rna_renamed$gene_name), stjude_rna_renamed$gene_name)
  cat("Common genes across all RNA datasets:", length(common_genes), "\n")
  
  # Filter datasets to common genes and combine
  kics_common <- kics_rna_renamed[kics_rna_renamed$gene_name %in% common_genes, ]
  lfs_common <- lfs_rna_renamed[lfs_rna_renamed$gene_name %in% common_genes, ]
  stjude_common <- stjude_rna_renamed[stjude_rna_renamed$gene_name %in% common_genes, ]
  
  # No need for unified structure - just combine datasets as-is
  # Remove dataset column for binding and keep only actual samples
  kics_clean <- kics_common %>% select(-dataset)
  lfs_clean <- lfs_common %>% select(-dataset)
  stjude_clean <- stjude_common %>% select(-dataset)
  
  # Combine datasets - dplyr will handle mismatched columns automatically
  rna_combined <- bind_rows(kics_clean, lfs_clean, stjude_clean, .id = "source")
  
  # Aggregate by gene (sum expression values across datasets)
  # Get all sample columns (excluding gene_name and source)
  sample_cols <- setdiff(colnames(rna_combined), c("gene_name", "source"))
  
  rna_aggregated <- rna_combined %>%
    select(-source) %>%
    group_by(gene_name) %>%
    summarise(across(all_of(sample_cols), \(x) sum(x, na.rm = TRUE)), .groups = 'drop')
  
  # Filter for genes with reasonable expression (remove all-zero genes)
  genes_with_expression <- rna_aggregated %>%
    filter(if_any(-gene_name, ~ .x > 0))
  
  cat("Genes with expression data:", nrow(genes_with_expression), "\n")
  
  # Get final sample list from actual data
  final_samples <- setdiff(colnames(genes_with_expression), "gene_name")
  
  return(list(
    rna_data = genes_with_expression,
    all_samples = final_samples,
    common_genes = common_genes
  ))
}

#' Process germline RNA data - combine KICS, LFS, and St. Jude datasets
#' @param kics_rna KICS RNA expression data
#' @param lfs_rna LFS RNA expression data
#' @param stjude_rna St. Jude RNA expression data
#' @param matched_dna_rna DNA-RNA mapping file
#' @param lfs_wgs2rna LFS WGS-RNA mapping file
#' @return List with combined RNA data and metadata
process_germline_rna_data <- function(kics_rna, lfs_rna, stjude_rna, matched_dna_rna, lfs_wgs2rna) {
  # Rename RNA samples to match TE sample IDs (same as tumor)
  kics_rna_renamed <- rename_rna_samples_kics_tumor(kics_rna)
  lfs_rna_renamed <- rename_rna_samples_lfs(lfs_rna, lfs_wgs2rna)
  stjude_rna_renamed <- rename_rna_samples_stjude(stjude_rna)

  # Find common genes across all datasets
  common_genes <- intersect(intersect(kics_rna_renamed$gene_name, lfs_rna_renamed$gene_name), stjude_rna_renamed$gene_name)

  # Filter datasets to common genes
  kics_common <- kics_rna_renamed[kics_rna_renamed$gene_name %in% common_genes, ]
  lfs_common <- lfs_rna_renamed[lfs_rna_renamed$gene_name %in% common_genes, ]
  stjude_common <- stjude_rna_renamed[stjude_rna_renamed$gene_name %in% common_genes, ]

  # Combine datasets
  rna_combined <- bind_rows(kics_common, lfs_common, stjude_common, .id = "source")

  # Aggregate by gene (sum expression values across datasets)
  sample_cols <- setdiff(colnames(rna_combined), c("gene_name", "source"))

  rna_aggregated <- rna_combined %>%
    select(-source) %>%
    group_by(gene_name) %>%
    summarise(across(all_of(sample_cols), \(x) sum(x, na.rm = TRUE)), .groups = 'drop')

  # Filter for genes with expression (remove all-zero genes)
  genes_with_expression <- rna_aggregated %>%
    filter(if_any(-gene_name, ~ .x > 0))

  # Get final sample list
  final_samples <- setdiff(colnames(genes_with_expression), "gene_name")

  return(list(
    rna_data = genes_with_expression,
    all_samples = final_samples,
    common_genes = common_genes
  ))
}

#' Match germline TE samples with RNA samples and create filtered dataset
#' @param rna_ready_for_filtering Prepared RNA data ready for TE matching
#' @param te_aff_split Germline TE split data with sample information
#' @param r_dir_files Output directory for CSV files (optional)
#' @return List with filtered RNA data and matching statistics
match_te_rna_samples_germline <- function(rna_ready_for_filtering, te_aff_split, r_dir_files = NULL) {
  cat("\n=== MATCHING RNA AND TE SAMPLES (GERMLINE) ===\n")

  # Get unique TE samples
  te_samples <- unique(te_aff_split$sample)
  cat("Total unique TE samples:", length(te_samples), "\n")

  # Get RNA sample IDs (already renamed to match TE IDs)
  rna_samples_available <- colnames(rna_ready_for_filtering)[-1]  # Exclude gene_name
  cat("Total RNA samples available:", length(rna_samples_available), "\n")

  # Step 1: Direct full ID matching
  matched_rna_samples <- intersect(te_samples, rna_samples_available)
  cat("Direct matches:", length(matched_rna_samples), "\n")

  # Step 2: Base ID matching for unmatched TE patients
  remaining_rna_samples <- setdiff(rna_samples_available, matched_rna_samples)

  te_base_ids <- unique(sub("_.*", "", te_samples))
  matched_base_ids <- unique(sub("_.*", "", matched_rna_samples))
  unmatched_te_base_ids <- setdiff(te_base_ids, matched_base_ids)

  cat("TE patients without direct RNA matches:", length(unmatched_te_base_ids), "\n")

  if (length(unmatched_te_base_ids) > 0 && length(remaining_rna_samples) > 0) {
    rna_base_ids <- sub("_.*", "", remaining_rna_samples)
    matching_base_ids <- intersect(rna_base_ids, unmatched_te_base_ids)
    base_matched_rna <- remaining_rna_samples[rna_base_ids %in% matching_base_ids]

    matched_rna_samples <- c(matched_rna_samples, base_matched_rna)
    cat("Additional base ID matches:", length(base_matched_rna), "\n")
  }

  cat("Total matched RNA samples:", length(matched_rna_samples), "\n")
  cat("  From", length(unique(sub("_.*", "", matched_rna_samples))), "unique patients\n")

  # Filter RNA data to matched samples only
  rna_filtered <- rna_ready_for_filtering %>%
    select(gene_name, all_of(matched_rna_samples))

  # Report unmatched samples
  unmatched_te <- setdiff(te_samples, matched_rna_samples)
  unmatched_rna <- setdiff(rna_samples_available, matched_rna_samples)

  if (length(unmatched_te) > 0) {
    cat("TE samples without RNA:", length(unmatched_te), "\n")
    cat("  First few:", paste(head(unmatched_te, 5), collapse=", "), "\n")
  }
  if (length(unmatched_rna) > 0) {
    cat("RNA samples without TE:", length(unmatched_rna), "\n")
    cat("  First few:", paste(head(unmatched_rna, 5), collapse=", "), "\n")
  }

  cat("\nFinal RNA dataset:", nrow(rna_filtered), "genes x", ncol(rna_filtered)-1, "samples\n")

  # Export matched/unmatched samples to CSV if directory provided
  if (!is.null(r_dir_files)) {
    # Export matched samples
    if (length(matched_rna_samples) > 0) {
      matched_df <- data.frame(
        sample_id = matched_rna_samples,
        status = "Matched (TE + RNA)",
        stringsAsFactors = FALSE
      )
      write.csv(matched_df, paste0(r_dir_files, "matched_samples_germline.csv"), row.names = FALSE)
    }

    # Export unmatched TE samples
    if (length(unmatched_te) > 0) {
      unmatched_te_df <- data.frame(
        sample_id = unmatched_te,
        status = "TE only (no RNA)",
        stringsAsFactors = FALSE
      )
      write.csv(unmatched_te_df, paste0(r_dir_files, "unmatched_te_samples_germline.csv"), row.names = FALSE)
    }

    # Export unmatched RNA samples
    if (length(unmatched_rna) > 0) {
      unmatched_rna_df <- data.frame(
        sample_id = unmatched_rna,
        status = "RNA only (no TE)",
        stringsAsFactors = FALSE
      )
      write.csv(unmatched_rna_df, paste0(r_dir_files, "unmatched_rna_samples_germline.csv"), row.names = FALSE)
    }

    cat("\nMatching results exported to:", r_dir_files, "\n")
  }

  return(list(
    rna_filtered = rna_filtered
  ))
}

#' Match TE samples with RNA samples and create filtered dataset (using matched_dna_rna.csv)
#' @param rna_ready_for_filtering Prepared RNA data ready for TE matching
#' @param te_all_t TE data with sample information  
#' @param matched_dna_rna_file Path to matched_dna_rna.csv file (used if dna_rna_mapping not provided)
#' @param r_dir_files Output directory for CSV files
#' @param dna_rna_mapping Pre-loaded DNA-RNA mapping dataframe (optional)
#' @return List with filtered RNA data and matching statistics
match_te_rna_samples <- function(rna_ready_for_filtering, te_all_t, matched_dna_rna_file = "/Users/briannelaverty/Documents/R_Malkin/te/data/rna/matched_dna_rna.csv", r_dir_files, dna_rna_mapping = NULL, original_kics_rna = NULL, te_all_all_t = NULL) {
  cat("Completing tumor RNA data filtering with TE sample overlap...\n")
  
  # Get primary TE sample names (selected samples - use full ID matching)
  te_primary_samples <- unique(te_all_t$sample)
  cat("Primary TE sample names from te_all_t (first 10):", paste(head(te_primary_samples, 10), collapse=", "), "\n")
  
  # Get extended TE sample names for fallback matching
  te_extended_samples <- if (!is.null(te_all_all_t)) unique(te_all_all_t$sample) else c()
  if (!is.null(te_all_all_t)) {
    cat("Extended TE sample names from te_all_all_t (first 10):", paste(head(te_extended_samples, 10), collapse=", "), "\n")
  }
  
  # Load matched DNA-RNA mapping (use pre-loaded if available)
  if (is.null(dna_rna_mapping)) {
    dna_rna_mapping <- read.csv(matched_dna_rna_file, stringsAsFactors = FALSE)
    cat("Loaded DNA-RNA mapping from file:", matched_dna_rna_file, "\n")
  } else {
    cat("Using pre-loaded DNA-RNA mapping data\n")
  }
  
  # Get available RNA samples (now they should have DNA names)
  rna_samples_available <- colnames(rna_ready_for_filtering)[-1]  # Exclude gene_name
  
  # Step 1: Direct full ID matching with te_all_t (selected samples)
  matched_rna_samples <- intersect(te_primary_samples, rna_samples_available)
  cat("Direct matches with te_all_t:", length(matched_rna_samples), "\n")
  
  # Step 2: Base ID matching for patients in te_all_t that didn't get direct matches
  if (!is.null(te_all_all_t)) {
    remaining_rna_samples <- setdiff(rna_samples_available, matched_rna_samples)
    
    # Only consider base IDs from patients in te_all_t (selected samples)
    te_primary_base_ids <- unique(sub("_.*", "", te_primary_samples))
    matched_base_ids <- unique(sub("_.*", "", matched_rna_samples))
    unmatched_te_base_ids <- setdiff(te_primary_base_ids, matched_base_ids)
    
    cat("TE patients without direct RNA matches:", length(unmatched_te_base_ids), "\n")
    
    if (length(unmatched_te_base_ids) > 0) {
      # Extract base IDs for remaining RNA samples
      rna_base_ids <- sub("_.*", "", remaining_rna_samples)
      
      # Only match base IDs from unmatched TE patients
      matching_base_ids <- intersect(rna_base_ids, unmatched_te_base_ids)
      base_matched_rna <- remaining_rna_samples[rna_base_ids %in% matching_base_ids]
      
      # Add base ID matches to final list
      matched_rna_samples <- c(matched_rna_samples, base_matched_rna)
      cat("Additional base ID matches for unmatched TE patients:", length(base_matched_rna), "\n")
    }
  }
  
  # Special case: Force match for specific LFS samples with same base sample
  special_lfs_samples <- c("5009_4_T", "2921_5_T")
  for (lfs_sample in special_lfs_samples) {
    if (lfs_sample %in% rna_samples_available && !lfs_sample %in% matched_rna_samples) {
      # Extract base sample (e.g., "5009" from "5009_4_T")
      base_sample <- gsub("_.*$", "", lfs_sample)
      # Find any TE sample with same base sample from either dataset
      all_te_samples <- c(te_primary_samples, te_extended_samples)
      matching_te_samples <- all_te_samples[grepl(paste0("^", base_sample, "_"), all_te_samples)]
      if (length(matching_te_samples) > 0) {
        # Force this RNA sample to be considered "matched"
        matched_rna_samples <- c(matched_rna_samples, lfs_sample)
        cat(sprintf("Special case: Matched %s to TE samples with base %s\n", lfs_sample, base_sample))
      }
    }
  }
  
  cat("Total TE samples processed:", length(c(te_primary_samples, te_extended_samples)), "\n")
  cat("Available RNA samples:", length(rna_samples_available), "\n") 
  cat("Matched RNA samples:", length(matched_rna_samples), "\n")
  
  # Show matching details for primary TE samples
  for (te_sample in head(te_primary_samples, 10)) {
    if (te_sample %in% matched_rna_samples) {
      cat("✓ TE sample", te_sample, "-> direct match with RNA sample\n")
    } else {
      cat("✗ TE sample", te_sample, "-> no RNA sample match\n")
    }
  }
  if (length(te_primary_samples) > 10) {
    cat("... and", length(te_primary_samples) - 10, "more TE samples\n")
  }
  
  # Generate unmatched TE samples report  
  unmatched_te_samples <- setdiff(te_primary_samples, matched_rna_samples)
  
  # Get original RNA sample names for filtering stale mappings
  original_rna_samples <- NULL
  if (!is.null(original_kics_rna)) {
    original_rna_samples <- colnames(original_kics_rna)[-1]  # Exclude gene_name
  }
  
  unmatched_te_results <- export_unmatched_te_samples_new(unmatched_te_samples, dna_rna_mapping, r_dir_files, original_rna_samples)
  
  # Generate unmatched RNA samples report
  te_dataset_for_rna_filtering <- if (!is.null(te_all_all_t)) te_all_all_t else te_all_t
  unmatched_rna_samples <- setdiff(rna_samples_available, matched_rna_samples)
  unmatched_rna_results <- export_unmatched_rna_samples_new(unmatched_rna_samples, dna_rna_mapping, r_dir_files, original_kics_rna, te_dataset_for_rna_filtering, matched_rna_samples)
  
  # Export successful matches to CSV
  if (length(matched_rna_samples) > 0) {
    matched_samples_df <- data.frame(
      sample_name = matched_rna_samples,
      status = "Successfully matched TE and RNA data",
      stringsAsFactors = FALSE
    )
    write.csv(matched_samples_df, paste0(r_dir_files, "successfully_matched_samples.csv"), row.names = FALSE)

    # Print summary
    cat("\n=== FINAL MATCHING SUMMARY ===\n")
    cat(sprintf("✓ %d samples successfully matched (TE + RNA)\n", length(matched_rna_samples)))
    cat(sprintf("  - From %d unique patients\n", length(unique(gsub("_.*", "", matched_rna_samples)))))
    if (unmatched_te_results$unmatched_count > 0) {
      cat(sprintf("⚠ %d TE samples without RNA (see unmatched_te_samples_new.csv)\n", unmatched_te_results$unmatched_count))
    }
    if (unmatched_rna_results$unmatched_count > 0) {
      cat(sprintf("⚠ %d RNA samples without TE (see unmatched_rna_samples_new.csv)\n", unmatched_rna_results$unmatched_count))
    }
    cat(sprintf("\nDetails saved to: successfully_matched_samples.csv\n"))
  }
  
  # Apply same sample selection logic as TE data (one sample per patient)
  if (length(matched_rna_samples) > 0) {
    # Load clinical data to get patient info for sample selection
    if (!is.null(te_all_t)) {
      # Create a sample-to-patient mapping from TE data
      te_clinical_info <- te_all_t %>%
        select(sample, base_sample, age_at_enrollment, lesion_type, disease_state) %>%
        distinct()
      
      # Filter to only matched RNA samples that have clinical info
      rna_samples_with_clinical <- matched_rna_samples[matched_rna_samples %in% te_clinical_info$sample]
      
      if (length(rna_samples_with_clinical) > 0) {
        # Apply sample selection logic (same as TE data)
        selected_rna_samples_df <- te_clinical_info %>%
          filter(sample %in% rna_samples_with_clinical) %>%
          group_by(base_sample) %>%
          arrange(age_at_enrollment, 
                  ifelse(lesion_type == "primary", 1, 
                         ifelse(lesion_type %in% c("metastasis", "relapse"), 2, 3)),
                  ifelse(disease_state == "initial", 1,
                         ifelse(disease_state == "progressive", 2,
                                ifelse(disease_state %in% c("relapsed", "relapse"), 3, 4)))) %>%
          slice(1) %>%
          ungroup()
        
        selected_rna_samples <- selected_rna_samples_df$sample
        cat("Applied sample selection: ", length(rna_samples_with_clinical), "->", length(selected_rna_samples), "samples (one per patient)\n")
      } else {
        selected_rna_samples <- matched_rna_samples
        cat("No clinical info available for sample selection, using all matched samples\n")
      }
    } else {
      selected_rna_samples <- matched_rna_samples
      cat("No TE data available for sample selection, using all matched samples\n")
    }
    
    rna_filtered <- rna_ready_for_filtering %>%
      select(gene_name, all_of(selected_rna_samples))
    
    cat("Final selected RNA samples (first 10):", paste(head(selected_rna_samples, 10), collapse=", "), "\n")
  } else {
    # Fallback to no filtering if no matches found
    rna_filtered <- rna_ready_for_filtering
    cat("Warning: No RNA-TE sample matches found. Using all RNA data.\n")
  }
  
  # Check final results
  rna_samples_final <- colnames(rna_filtered)[-1]  # Exclude gene_name
  cat("Final tumor RNA dataset dimensions:", nrow(rna_filtered), "genes x", length(rna_samples_final), "samples\n")
  
  return(list(
    rna_filtered = rna_filtered,
    sample_overlap = matched_rna_samples,
    matched_samples = length(matched_rna_samples),
    total_te_samples = length(te_primary_samples),
    unmatched_te_info = unmatched_te_results,
    unmatched_rna_info = unmatched_rna_results
  ))
}

#' Export unmatched WGS to RNA mapping rows to CSV
#' @param te_sample_names Vector of TE sample names
#' @param rna_ready_for_filtering RNA data ready for filtering
#' @param kics_dna2rna KICS WGS to RNA mapping file
#' @param r_dir_files Output directory for CSV files
#' @return List with unmatched sample information
export_unmatched_wgs2rna_rows <- function(te_sample_names, rna_ready_for_filtering, kics_dna2rna, r_dir_files) {
  # Find unmatched TE samples
  unmatched_te_samples <- c()
  for (te_sample in te_sample_names) {
    mapping_row <- kics_dna2rna[kics_dna2rna$sample == te_sample, ]
    if (nrow(mapping_row) > 0 && !is.na(mapping_row$rna[1])) {
      rna_id <- mapping_row$rna[1]
      rna_samples_available <- colnames(rna_ready_for_filtering)[-1]
      matching_rna <- rna_samples_available[grepl(paste0("_", rna_id, "$"), rna_samples_available)]
      if (length(matching_rna) == 0) {
        unmatched_te_samples <- c(unmatched_te_samples, te_sample)
      }
    }
  }
  
  if (length(unmatched_te_samples) > 0) {
    cat("\nkics_wgs2rna_tumour rows that don't have matching RNA samples:\n")
    cat("Format: sample | dna | rna | (reason)\n")
    
    # Create dataframe for unmatched rows
    unmatched_rows <- data.frame()
    
    for (te_sample in head(unmatched_te_samples, 20)) {
      mapping_row <- kics_dna2rna[kics_dna2rna$sample == te_sample, ]
      if (nrow(mapping_row) > 0) {
        cat(sprintf("  %s | %s | %s | (RNA sample ending with %s not found)\n", 
                    mapping_row$sample[1], mapping_row$dna[1], mapping_row$rna[1], mapping_row$rna[1]))
      }
    }
    
    if (length(unmatched_te_samples) > 20) {
      cat("  ... and", length(unmatched_te_samples) - 20, "more unmapped rows\n")
    }
    
    # Create complete unmatched rows dataframe for CSV export
    for (te_sample in unmatched_te_samples) {
      mapping_row <- kics_dna2rna[kics_dna2rna$sample == te_sample, ]
      if (nrow(mapping_row) > 0) {
        unmatched_rows <- rbind(unmatched_rows, data.frame(
          sample = mapping_row$sample[1],
          dna = mapping_row$dna[1], 
          rna = mapping_row$rna[1],
          reason = paste("RNA sample ending with", mapping_row$rna[1], "not found"),
          stringsAsFactors = FALSE
        ))
      }
    }
    
    # Export unmatched rows to CSV
    write.csv(unmatched_rows, paste0(r_dir_files, "unmatched_wgs2rna_rows.csv"), row.names = FALSE)
    cat("Unmatched rows exported to:", paste0(r_dir_files, "unmatched_wgs2rna_rows.csv"), "\n")
    
    return(list(
      unmatched_count = length(unmatched_te_samples),
      unmatched_rows = unmatched_rows
    ))
  } else {
    cat("\nAll TE samples have matching RNA samples.\n")
    return(list(unmatched_count = 0, unmatched_rows = data.frame()))
  }
}

#' Export unmatched TE samples using new matched_dna_rna.csv approach
#' @param unmatched_te_samples Vector of TE sample names that don't have RNA matches
#' @param dna_rna_mapping DNA-RNA mapping dataframe
#' @param r_dir_files Output directory for CSV files
#' @return List with unmatched sample information
export_unmatched_te_samples_new <- function(unmatched_te_samples, dna_rna_mapping, r_dir_files, rna_samples_available = NULL) {
  if (length(unmatched_te_samples) > 0) {
    # Filter to only include samples that have mappings but missing RNA data
    samples_with_mapping <- unmatched_te_samples[unmatched_te_samples %in% dna_rna_mapping$sample_name_with_t]
    samples_without_mapping <- length(unmatched_te_samples) - length(samples_with_mapping)

    cat("\n=== TE-RNA MATCHING SUMMARY ===\n")

    # Create dataframe for unmatched TE samples (only those with mappings)
    unmatched_rows <- data.frame()

    if (length(samples_with_mapping) > 0) {
      cat(sprintf("⚠ %d TE samples expected RNA but none found:\n", length(samples_with_mapping)))
      for (te_sample in head(samples_with_mapping, 5)) {
        # Get mapping row (we know it exists)
        mapping_row <- dna_rna_mapping[dna_rna_mapping$sample_name_with_t == te_sample, ]

        cat(sprintf("  • %s (expected RNA: %s)\n",
                    te_sample, mapping_row$kics_rna_name[1]))
        reason <- "RNA sample processed but not available in expression data"

        unmatched_rows <- rbind(unmatched_rows, data.frame(
          te_sample = te_sample,
          dna_name = te_sample,
          rna_name = mapping_row$kics_rna_name[1],
          reason = reason,
          stringsAsFactors = FALSE
        ))
      }

      if (length(samples_with_mapping) > 5) {
        cat(sprintf("  ... and %d more (see CSV file)\n", length(samples_with_mapping) - 5))
      }
    } else {
      cat("✓ All TE samples with RNA mappings have RNA data\n")
    }

    # Report samples without mapping (but don't include in unmatched list)
    if (samples_without_mapping > 0) {
      cat(sprintf("✓ %d TE samples are DNA-only (no RNA expected)\n", samples_without_mapping))
    }
    
    # Create complete unmatched dataframe for remaining samples with mappings
    if (length(samples_with_mapping) > 5) {
      for (te_sample in tail(samples_with_mapping, -5)) {
        mapping_row <- dna_rna_mapping[dna_rna_mapping$sample_name_with_t == te_sample, ]
        reason <- "RNA sample processed but not available in expression data"

        unmatched_rows <- rbind(unmatched_rows, data.frame(
          te_sample = te_sample,
          dna_name = te_sample,
          rna_name = mapping_row$kics_rna_name[1],
          reason = reason,
          stringsAsFactors = FALSE
        ))
      }
    }

    # Export unmatched rows to CSV
    if (nrow(unmatched_rows) > 0) {
      write.csv(unmatched_rows, paste0(r_dir_files, "unmatched_te_samples_new.csv"), row.names = FALSE)
      cat(sprintf("\nDetails saved to: unmatched_te_samples_new.csv\n"))
    }
    
    return(list(
      unmatched_count = length(unmatched_te_samples),
      unmatched_rows = unmatched_rows
    ))
  } else {
    cat("\nAll TE samples have matching RNA samples.\n")
    return(list(unmatched_count = 0, unmatched_rows = data.frame()))
  }
}

#' Export RNA samples that don't have corresponding TE/DNA data
#' @param unmatched_rna_samples Vector of RNA sample names that don't have TE matches
#' @param dna_rna_mapping DNA-RNA mapping dataframe
#' @param r_dir_files Output directory for CSV files
#' @return List with unmatched RNA sample information
export_unmatched_rna_samples_new <- function(unmatched_rna_samples, dna_rna_mapping, r_dir_files, original_kics_rna = NULL, te_all_t = NULL, matched_rna_samples = NULL) {
  if (length(unmatched_rna_samples) > 0) {
    # Get base samples that already have matched RNA
    matched_base_samples <- c()
    if (!is.null(matched_rna_samples)) {
      matched_base_samples <- unique(gsub("_.*$", "", matched_rna_samples))
    }

    # Filter to only include RNA samples where NO RNA from that base sample matched to te_all_t
    filtered_rna_samples <- unmatched_rna_samples
    if (!is.null(te_all_t)) {
      # Get TE sample names
      te_sample_names <- unique(te_all_t$sample)

      # Extract base sample IDs from TE samples
      te_base_samples <- unique(gsub("_.*$", "", te_sample_names))

      # Filter unmatched RNA samples:
      # 1. Base sample must exist in TE data
      # 2. Base sample must NOT have any RNA already matched
      filtered_rna_samples <- c()
      for (rna_sample in unmatched_rna_samples) {
        base_sample <- gsub("_.*$", "", rna_sample)  # Extract base sample ID
        if (base_sample %in% te_base_samples && !base_sample %in% matched_base_samples) {
          filtered_rna_samples <- c(filtered_rna_samples, rna_sample)
        }
      }

    }

    if (length(filtered_rna_samples) > 0) {
      cat(sprintf("\n⚠ %d RNA samples (from %d patients) have no matching TE data:\n",
                  length(filtered_rna_samples), length(unique(gsub("_.*", "", filtered_rna_samples)))))

      # Create dataframe for unmatched RNA samples
      unmatched_rna_rows <- data.frame()

      for (rna_sample in head(filtered_rna_samples, 5)) {
      # Check if RNA sample is a DNA-format name that should have a reverse mapping
      mapping_row <- dna_rna_mapping[dna_rna_mapping$sample_name_with_t == rna_sample, ]

      if (nrow(mapping_row) > 0) {
        # This RNA sample has a mapping but no TE data
        cat(sprintf("  • %s (original: %s)\n",
                    rna_sample, mapping_row$kics_rna_name[1]))
        reason <- "Mapped TE sample not found"
        original_rna_name <- mapping_row$kics_rna_name[1]
      } else {
        # Check if this might be an LFS or St. Jude sample (different naming scheme)
        cat(sprintf("  • %s (LFS/StJude)\n", rna_sample))
        reason <- "RNA sample not in mapping file"
        original_rna_name <- rna_sample
      }

      unmatched_rna_rows <- rbind(unmatched_rna_rows, data.frame(
        rna_sample = rna_sample,
        original_rna_name = original_rna_name,
        reason = reason,
        stringsAsFactors = FALSE
      ))
    }

      if (length(filtered_rna_samples) > 5) {
        cat(sprintf("  ... and %d more (see CSV file)\n", length(filtered_rna_samples) - 5))
      }

      # Create complete unmatched dataframe for remaining RNA samples
      if (length(filtered_rna_samples) > 5) {
        for (rna_sample in tail(filtered_rna_samples, -5)) {
        mapping_row <- dna_rna_mapping[dna_rna_mapping$sample_name_with_t == rna_sample, ]

        if (nrow(mapping_row) > 0) {
          reason <- "Mapped TE sample not found"
          original_rna_name <- mapping_row$kics_rna_name[1]
        } else {
          reason <- "RNA sample not in mapping file"
          original_rna_name <- rna_sample
        }

        unmatched_rna_rows <- rbind(unmatched_rna_rows, data.frame(
          rna_sample = rna_sample,
          original_rna_name = original_rna_name,
          reason = reason,
          stringsAsFactors = FALSE
        ))
      }
    }

      # Export unmatched RNA rows to CSV
      write.csv(unmatched_rna_rows, paste0(r_dir_files, "unmatched_rna_samples_new.csv"), row.names = FALSE)
      cat(sprintf("\nDetails saved to: unmatched_rna_samples_new.csv\n"))

      return(list(
        unmatched_count = length(filtered_rna_samples),
        unmatched_rows = unmatched_rna_rows
      ))
    } else {
      cat("\n✓ All patients with TE data have matched RNA samples\n")
      return(list(unmatched_count = 0, unmatched_rows = data.frame()))
    }
  } else {
    cat("\n✓ All RNA samples have corresponding TE data\n")
    return(list(unmatched_count = 0, unmatched_rows = data.frame()))
  }
}

# Prepare survival data for Kaplan-Meier analysis
# Merges TE burden data with survival outcomes and creates survival time/event variables
prepare_survival_data <- function(te_data, dod_data) {
  # Prepare TE data with base sample (extract everything before first underscore)
  te_survival <- te_data %>%
    mutate(base_sample = sub("_.*", "", sample))

  # Clean DOD data - pad KiCS ID with leading zeros to 4 digits (e.g., 2 -> "0002")
  kics_DOD_clean <- dod_data %>%
    rename(kics_id = `KiCS ID`) %>%
    mutate(base_sample = sprintf("%04d", kics_id))

  cat("DEBUG: First 10 base_sample from TE data:\n")
  print(head(unique(te_survival$base_sample), 10))
  cat("DEBUG: First 10 base_sample from DOD data:\n")
  print(head(kics_DOD_clean$base_sample, 10))

  # Merge TE burden with survival data
  survival_data <- te_survival %>%
    select(base_sample, sample, total, LINE1, ALU, SVA, TP53_status, tumor_type, age_at_diagnosis) %>%
    left_join(kics_DOD_clean, by = "base_sample") %>%
    filter(!is.na(base_sample))

  cat("Survival data merged:\n")
  cat("  Total samples with TE data:", nrow(te_survival), "\n")
  cat("  Samples matched with survival data:", sum(!is.na(survival_data$base_sample)), "\n")
  cat("  Samples with Vital Status available:", sum(!is.na(survival_data$`Vital Status`)), "\n")

  # Define high vs low TE burden (median split)
  # If median is 0, use samples with total > 0 as "High"
  median_te <- median(survival_data$total, na.rm = TRUE)

  if (median_te == 0) {
    cat("  Note: Median TE count is 0, using presence/absence split instead\n")
    survival_data <- survival_data %>%
      mutate(te_burden = ifelse(total > 0, "High", "Low"))
  } else {
    survival_data <- survival_data %>%
      mutate(te_burden = ifelse(total >= median_te, "High", "Low"))
  }

  cat("TE burden groups:\n")
  cat("  Median TE count:", median_te, "\n")
  cat("  Number of samples with High TE burden:", sum(survival_data$te_burden == "High"), "\n")
  cat("  Number of samples with Low TE burden:", sum(survival_data$te_burden == "Low"), "\n")

  # Create time and event columns from kics_DOD data
  # Handle "Not applicable" and "UNK" values in Age at Death

  # Check if required columns exist
  if (!"Vital Status" %in% colnames(survival_data)) {
    cat("ERROR: 'Vital Status' column not found in merged data\n")
    cat("Available columns:", paste(colnames(survival_data), collapse=", "), "\n")
    return(data.frame())
  }

  if (!"Age at Death (days)" %in% colnames(survival_data)) {
    cat("ERROR: 'Age at Death (days)' column not found in merged data\n")
    cat("Available columns:", paste(colnames(survival_data), collapse=", "), "\n")
    return(data.frame())
  }

  survival_data <- survival_data %>%
    mutate(
      # Convert Age at Death to numeric (will turn "Not applicable" and "UNK" to NA)
      age_at_death_days = suppressWarnings(as.numeric(`Age at Death (days)`))
    )

  # Create event indicator
  survival_data <- survival_data %>%
    mutate(
      event = case_when(
        `Vital Status` == "Dead" ~ 1,
        `Vital Status` == "Alive" ~ 0,
        TRUE ~ NA_real_
      )
    )

  cat("After creating event indicator:", nrow(survival_data), "rows\n")
  cat("  Dead:", sum(survival_data$event == 1, na.rm=TRUE), "\n")
  cat("  Alive:", sum(survival_data$event == 0, na.rm=TRUE), "\n")
  cat("  NA events:", sum(is.na(survival_data$event)), "\n")

  # Create time variable (time from diagnosis to death/censoring)
  # For dead patients: age at death - age at diagnosis
  # For alive patients: We need to estimate time from diagnosis to today
  # Assuming diagnosis happened at age_at_diagnosis, we need current date to calculate follow-up
  # Use today's date as censoring date for alive patients

  current_year <- as.numeric(format(Sys.Date(), "%Y"))
  current_day_of_year <- as.numeric(format(Sys.Date(), "%j"))

  survival_data <- survival_data %>%
    mutate(
      age_at_diagnosis_days = age_at_diagnosis * 365.25,  # convert years to days
      # For alive patients, we need to estimate years since diagnosis
      # This requires knowing diagnosis year/date, which we may not have
      # As a proxy, assume diagnosis happened recently and use a conservative estimate
      time = case_when(
        event == 1 & !is.na(age_at_death_days) & !is.na(age_at_diagnosis_days) ~
          age_at_death_days - age_at_diagnosis_days,  # survival time from diagnosis to death
        event == 0 & !is.na(age_at_diagnosis_days) ~
          # For alive: use 5 years as default follow-up time (conservative estimate)
          # This should be replaced with actual date of last contact if available
          5 * 365.25,
        TRUE ~ NA_real_
      )
    )

  cat("WARNING: For alive patients, using 5-year default follow-up time.\n")
  cat("         Ideally, date of last contact should be available in the data.\n")

  cat("After creating time variable:", nrow(survival_data), "rows\n")
  cat("  Rows with valid time:", sum(!is.na(survival_data$time)), "\n")
  cat("  Rows with time > 0:", sum(survival_data$time > 0, na.rm=TRUE), "\n")
  cat("  Rows with NA time:", sum(is.na(survival_data$time)), "\n")
  cat("  Rows with time <= 0:", sum(survival_data$time <= 0, na.rm=TRUE), "\n")

  survival_data <- survival_data %>%
    filter(!is.na(time), !is.na(event), time > 0)  # Remove rows with missing or invalid survival data

  cat("Samples with complete survival data:", nrow(survival_data), "\n")
  if (nrow(survival_data) > 0) {
    cat("  Event counts - Dead:", sum(survival_data$event == 1),
        ", Alive:", sum(survival_data$event == 0), "\n")

    # Check if we have both groups
    n_high <- sum(survival_data$te_burden == "High")
    n_low <- sum(survival_data$te_burden == "Low")
    if (n_high == 0 || n_low == 0) {
      cat("  WARNING: Only one TE burden group present. Cannot perform survival comparison.\n")
    }
  }

  return(survival_data)
}

# Plot Kaplan-Meier survival curves comparing high vs low TE burden
# Returns list with plot object, fit object, and test results
plot_survival_curves <- function(survival_data, output_dir = NULL) {
  if (nrow(survival_data) == 0 || sum(!is.na(survival_data$time)) == 0) {
    cat("Warning: No samples with complete survival data available for plotting.\n")
    return(NULL)
  }

  # Check if we have both groups
  n_high <- sum(survival_data$te_burden == "High")
  n_low <- sum(survival_data$te_burden == "Low")
  if (n_high == 0 || n_low == 0) {
    cat("Warning: Only one TE burden group present (High:", n_high, ", Low:", n_low, "). Cannot perform survival comparison.\n")
    return(NULL)
  }

  # Calculate Kaplan-Meier estimates manually for ggplot
  km_data <- survival_data %>%
    arrange(te_burden, time) %>%
    group_by(te_burden) %>%
    mutate(
      n_risk = n():1,  # number at risk (reverse order)
      n_event = event,
      surv_prob = cumprod(1 - n_event / n_risk)
    ) %>%
    ungroup()

  # Create step data for plotting
  km_steps <- km_data %>%
    group_by(te_burden) %>%
    arrange(time) %>%
    mutate(
      time_end = lead(time, default = max(time) * 1.1),
      surv_prob_carry = surv_prob
    ) %>%
    select(te_burden, time, time_end, surv_prob = surv_prob_carry) %>%
    ungroup()

  # Create ggplot
  p_survival <- ggplot(km_steps, aes(x = time, y = surv_prob, color = te_burden)) +
    geom_step(linewidth = 1.2) +
    scale_color_manual(
      values = c("High" = "#E7B800", "Low" = "#2E9FDF"),
      labels = c("High TE burden", "Low TE burden")
    ) +
    labs(
      title = "Kaplan-Meier Survival Curves by TE Burden",
      x = "Time (days)",
      y = "Survival Probability",
      color = "TE Burden"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5, face = "bold")
    ) +
    ylim(0, 1)

  # Add sample size info
  n_high_text <- paste0("High TE burden: n=", n_high)
  n_low_text <- paste0("Low TE burden: n=", n_low)

  p_survival <- p_survival +
    annotate("text", x = max(km_steps$time) * 0.7, y = 0.15,
             label = n_high_text, color = "#E7B800", size = 4, hjust = 0) +
    annotate("text", x = max(km_steps$time) * 0.7, y = 0.08,
             label = n_low_text, color = "#2E9FDF", size = 4, hjust = 0)

  # Perform log-rank test if survival package is available
  logrank_test <- NULL
  fit <- NULL
  if (requireNamespace("survival", quietly = TRUE)) {
    tryCatch({
      library(survival)
      surv_obj <- Surv(time = survival_data$time, event = survival_data$event)
      fit <- survfit(surv_obj ~ te_burden, data = survival_data)
      logrank_test <- survdiff(surv_obj ~ te_burden, data = survival_data)

      # Extract p-value
      pval <- 1 - pchisq(logrank_test$chisq, df = 1)
      pval_text <- sprintf("Log-rank p = %.4f", pval)

      # Add p-value to plot
      p_survival <- p_survival +
        annotate("text", x = max(km_steps$time) * 0.7, y = 0.95,
                 label = pval_text, size = 4.5, fontface = "bold")

      cat("\nLog-rank test results:\n")
      print(logrank_test)
      cat("\nSurvival summary:\n")
      print(summary(fit))
    }, error = function(e) {
      cat("Warning: Could not perform log-rank test:", e$message, "\n")
    })
  }

  # Save plot if output directory provided
  if (!is.null(output_dir)) {
    ggsave(paste0(output_dir, "survival_te_burden.png"),
           plot = p_survival, width = 10, height = 8)
    cat("✓ Survival plot saved to:", paste0(output_dir, "survival_te_burden.png"), "\n")
  }

  return(list(
    plot = p_survival,
    fit = fit,
    logrank_test = logrank_test
  ))
}

# Test specific TE insertions for differential representation between groups
# Uses Fisher's exact test with Benjamini-Hochberg multiple testing correction
# Only tests TEs with sufficient sample representation (min_samples_with and min_samples_without)
#
# Args:
#   te_expand: Expanded TE dataframe (one row per TE insertion)
#   te_count: Count matrix (one row per sample) - used for total sample counts
#   group_column: Column name for grouping (e.g., "TP53_status")
#   min_samples_with: Minimum samples that must have the TE (default: 5)
#   min_samples_without: Minimum samples that must not have the TE (default: 5)
#   output_dir: Directory to save CSV files (optional)
#   output_prefix: Prefix for output file names (default: "specific_tes")
#
# Returns:
#   List with full results and significant results dataframes
test_specific_tes_by_group <- function(te_expand, te_count, group_column,
                                       min_samples_with = 5, min_samples_without = 5,
                                       output_dir = NULL, output_prefix = "specific_tes") {

  cat("Testing individual TE insertions for differential representation by", group_column, "...\n")

  # Create a unique TE identifier for each insertion
  te_expand_annotated <- te_expand %>%
    mutate(te_id = paste(SV_chrom, SV_start, ALT, sep = "_"))

  # Get unique group values
  groups <- unique(te_count[[group_column]])
  groups <- groups[!is.na(groups)]

  if (length(groups) != 2) {
    cat("Error: Expected exactly 2 groups, found", length(groups), "\n")
    return(NULL)
  }

  group1 <- groups[1]
  group2 <- groups[2]

  cat("Comparing groups:", group1, "vs", group2, "\n")

  # Get TE counts by group for each unique TE
  te_by_group <- te_expand_annotated %>%
    filter(!is.na(.data[[group_column]])) %>%
    group_by(te_id, SV_chrom, SV_start, SV_end, ALT, .data[[group_column]]) %>%
    summarise(
      n_samples = n_distinct(sample),
      samples = paste(unique(sample), collapse = ";"),
      .groups = "drop"
    ) %>%
    pivot_wider(
      names_from = all_of(group_column),
      values_from = c(n_samples, samples),
      values_fill = list(n_samples = 0, samples = "")
    )

  # Calculate total samples per group
  n_group1 <- te_count %>% filter(.data[[group_column]] == group1) %>% nrow()
  n_group2 <- te_count %>% filter(.data[[group_column]] == group2) %>% nrow()

  cat("Total samples:", group1, "=", n_group1, ",", group2, "=", n_group2, "\n")

  # Get column names dynamically
  n_samples_col1 <- paste0("n_samples_", group1)
  n_samples_col2 <- paste0("n_samples_", group2)

  # Perform Fisher's exact test for each TE
  te_test_results <- te_by_group %>%
    mutate(
      total_samples = .data[[n_samples_col1]] + .data[[n_samples_col2]],
      freq_group1 = .data[[n_samples_col1]] / n_group1,
      freq_group2 = .data[[n_samples_col2]] / n_group2,
      fold_change = ifelse(freq_group1 == 0, Inf, freq_group2 / freq_group1),
      samples_with_te = total_samples,
      samples_without_te = (n_group1 + n_group2) - total_samples
    ) %>%
    # Filter for TEs with sufficient representation in EACH group
    filter(.data[[n_samples_col1]] >= min_samples_with & .data[[n_samples_col2]] >= min_samples_with) %>%
    rowwise() %>%
    mutate(
      # Fisher's exact test for count data
      p_value = tryCatch({
        fisher_matrix <- matrix(c(.data[[n_samples_col2]], n_group2 - .data[[n_samples_col2]],
                                   .data[[n_samples_col1]], n_group1 - .data[[n_samples_col1]]),
                                nrow = 2, byrow = TRUE)
        fisher.test(fisher_matrix)$p.value
      }, error = function(e) NA_real_)
    ) %>%
    ungroup() %>%
    filter(!is.na(p_value)) %>%
    mutate(
      p_adj_BH = p.adjust(p_value, method = "BH"),
      significant_BH = p_adj_BH < 0.05,
      enriched_in = case_when(
        !significant_BH ~ "Not significant",
        freq_group2 > freq_group1 ~ group2,
        freq_group1 > freq_group2 ~ group1,
        TRUE ~ "Equal"
      )
    ) %>%
    arrange(p_adj_BH)

  # Print summary
  cat("\nTotal TEs tested:", nrow(te_test_results), "\n")
  cat("  (Filtered to TEs with >=", min_samples_with, "samples with TE in EACH group)\n")
  cat("Significant TEs (BH < 0.05):", sum(te_test_results$significant_BH), "\n")
  cat("  Enriched in", group2, ":", sum(te_test_results$enriched_in == group2), "\n")
  cat("  Enriched in", group1, ":", sum(te_test_results$enriched_in == group1), "\n")

  # Save to CSV if output directory provided
  if (!is.null(output_dir)) {
    # Full results
    te_by_group_file <- paste0(output_dir, output_prefix, "_by_", group_column, ".csv")
    write.csv(te_test_results, te_by_group_file, row.names = FALSE)
    cat("✓ Full results saved to:", basename(te_by_group_file), "\n")

    # Count significant results
    te_test_results_sig <- te_test_results %>% filter(significant_BH)
    if (nrow(te_test_results_sig) == 0) {
      cat("  No significant results found\n")
    } else {
      cat("  Found", nrow(te_test_results_sig), "significant results\n")
    }
  }

  return(list(
    full_results = te_test_results,
    significant_results = te_test_results %>% filter(significant_BH)
  ))
}

# Test specific TE insertions grouped by gene for differential representation between groups
# Aggregates all TEs within each gene before testing
# Uses Fisher's exact test with Benjamini-Hochberg multiple testing correction
#
# Args:
#   te_expand: Expanded TE dataframe (one row per TE insertion) - must have Gene_name column
#   te_count: Count matrix (one row per sample) - used for total sample counts
#   group_column: Column name for grouping (e.g., "TP53_status")
#   min_samples_with: Minimum samples that must have TEs in the gene (default: 5)
#   min_samples_without: Minimum samples that must not have TEs in the gene (default: 5)
#   gene_filter: Optional vector of genes to test (e.g., cancer genes only)
#   output_dir: Directory to save CSV files (optional)
#   output_prefix: Prefix for output file names (default: "specific_tes_by_gene")
#
# Returns:
#   List with full results and significant results dataframes
test_specific_tes_by_gene <- function(te_expand, te_count, group_column,
                                      min_samples_with = 5, min_samples_without = 5,
                                      gene_filter = NULL,
                                      output_dir = NULL, output_prefix = "specific_tes_by_gene") {

  cat("Testing gene-level TE burden for differential representation by", group_column, "...\n")

  # Filter to specific genes if requested
  if (!is.null(gene_filter)) {
    te_expand <- te_expand %>% filter(Gene_name %in% gene_filter)
    cat("Filtered to", length(gene_filter), "genes\n")
  }

  # Get unique group values
  groups <- unique(te_count[[group_column]])
  groups <- groups[!is.na(groups)]

  if (length(groups) != 2) {
    cat("Error: Expected exactly 2 groups, found", length(groups), "\n")
    return(NULL)
  }

  group1 <- groups[1]
  group2 <- groups[2]

  cat("Comparing groups:", group1, "vs", group2, "\n")

  # Get samples with TEs in each gene by group
  gene_by_group <- te_expand %>%
    filter(!is.na(.data[[group_column]])) %>%
    group_by(Gene_name, .data[[group_column]]) %>%
    summarise(
      n_samples = n_distinct(sample),
      samples = paste(unique(sample), collapse = ";"),
      n_insertions = n(),
      .groups = "drop"
    ) %>%
    pivot_wider(
      names_from = all_of(group_column),
      values_from = c(n_samples, samples, n_insertions),
      values_fill = list(n_samples = 0, samples = "", n_insertions = 0)
    )

  # Calculate total samples per group
  n_group1 <- te_count %>% filter(.data[[group_column]] == group1) %>% nrow()
  n_group2 <- te_count %>% filter(.data[[group_column]] == group2) %>% nrow()

  cat("Total samples:", group1, "=", n_group1, ",", group2, "=", n_group2, "\n")

  # Get column names dynamically
  n_samples_col1 <- paste0("n_samples_", group1)
  n_samples_col2 <- paste0("n_samples_", group2)
  n_insertions_col1 <- paste0("n_insertions_", group1)
  n_insertions_col2 <- paste0("n_insertions_", group2)

  # Perform Fisher's exact test for each gene
  gene_test_results <- gene_by_group %>%
    mutate(
      total_samples = .data[[n_samples_col1]] + .data[[n_samples_col2]],
      total_insertions = .data[[n_insertions_col1]] + .data[[n_insertions_col2]],
      freq_group1 = .data[[n_samples_col1]] / n_group1,
      freq_group2 = .data[[n_samples_col2]] / n_group2,
      fold_change = ifelse(freq_group1 == 0, Inf, freq_group2 / freq_group1)
    ) %>%
    # Filter for genes with sufficient representation
    filter(.data[[n_samples_col1]] >= min_samples_with & .data[[n_samples_col2]] >= min_samples_with) %>%
    rowwise() %>%
    mutate(
      # Fisher's exact test
      p_value = tryCatch({
        fisher_matrix <- matrix(c(.data[[n_samples_col2]], n_group2 - .data[[n_samples_col2]],
                                   .data[[n_samples_col1]], n_group1 - .data[[n_samples_col1]]),
                                nrow = 2, byrow = TRUE)
        fisher.test(fisher_matrix)$p.value
      }, error = function(e) NA_real_)
    ) %>%
    ungroup() %>%
    filter(!is.na(p_value)) %>%
    mutate(
      p_adj = p.adjust(p_value, method = "BH"),
      significant = p_adj < 0.05,
      enriched_in = case_when(
        !significant ~ "Not significant",
        freq_group2 > freq_group1 ~ group2,
        freq_group1 > freq_group2 ~ group1,
        TRUE ~ "Equal"
      )
    ) %>%
    arrange(p_value)

  cat("✓ Tested", nrow(gene_test_results), "genes\n")
  cat("  Found", sum(gene_test_results$significant), "significant genes (FDR < 0.05)\n")

  # Save results if output directory provided
  if (!is.null(output_dir)) {
    gene_file <- paste0(output_dir, output_prefix, "_full_results.csv")
    write.csv(gene_test_results, gene_file, row.names = FALSE)
    cat("✓ Full results saved to:", basename(gene_file), "\n")
  }

  return(list(
    full_results = gene_test_results,
    significant_results = gene_test_results %>% filter(significant)
  ))
}

# Plot TE frequency distributions for gnomAD and HostSeq before common filtering
# Creates histograms showing how common/rare TEs are in reference populations
#
# Args:
#   te_data: TE dataframe after prep_te (must have GRPMAX_AF column)
#   output_dir: Directory to save plots
#   output_prefix: Prefix for output file names (e.g., "tumour" or "germline")
#
# Returns:
#   List with gnomad_plot and hostseq_plot
plot_te_frequency_distributions <- function(te_data, output_dir = NULL, output_prefix = "te") {

  cat("\n===== TE FREQUENCY DISTRIBUTIONS =====\n")

  # Check if GRPMAX_AF column exists
  if (!"GRPMAX_AF" %in% colnames(te_data)) {
    cat("Warning: GRPMAX_AF column not found in data. Cannot plot gnomAD frequencies.\n")
    return(NULL)
  }

  # Filter to unique TEs (by genomic position)
  te_unique <- te_data %>%
    distinct(SV_chrom, SV_start, ALT, .keep_all = TRUE)

  cat("Total unique TEs:", nrow(te_unique), "\n")

  # ===== PLOT 1: gnomAD Frequency Distribution =====
  cat("\n--- gnomAD Frequency Distribution ---\n")

  # Count TEs by gnomAD AF bins
  gnomad_summary <- te_unique %>%
    mutate(
      gnomad_AF = as.numeric(GRPMAX_AF),
      gnomad_category = case_when(
        is.na(gnomad_AF) | gnomad_AF == 0 ~ "Not in gnomAD (AF=0)",
        gnomad_AF < 0.001 ~ "Very rare (AF<0.1%)",
        gnomad_AF < 0.01 ~ "Rare (0.1-1%)",
        gnomad_AF < 0.05 ~ "Uncommon (1-5%)",
        gnomad_AF >= 0.05 ~ "Common (AF≥5%)",
        TRUE ~ "Unknown"
      )
    )

  # Print summary
  gnomad_counts <- table(gnomad_summary$gnomad_category)
  cat("gnomAD frequency categories:\n")
  print(gnomad_counts)

  # Create histogram
  gnomad_plot <- ggplot(gnomad_summary %>% filter(!is.na(gnomad_AF), gnomad_AF > 0),
                         aes(x = gnomad_AF)) +
    geom_histogram(bins = 50, fill = "blue", color = "black", alpha = 0.7) +
    scale_x_log10(labels = scales::percent) +
    labs(
      title = "TE Frequency Distribution in gnomAD",
      x = "gnomAD Allele Frequency (GRPMAX_AF)",
      y = "Number of TEs"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.grid.minor = element_blank()
    ) +
    geom_vline(xintercept = 0.03, linetype = "dashed", color = "red", size = 1) +
    annotate("text", x = 0.03, y = Inf, label = "3% threshold",
             vjust = -0.5, hjust = -0.1, color = "red", size = 3.5)

  # Save plot
  if (!is.null(output_dir)) {
    ggsave(paste0(output_dir, output_prefix, "_gnomad_frequency.png"),
           plot = gnomad_plot, width = 10, height = 6)
    cat("✓ gnomAD plot saved to:", paste0(output_prefix, "_gnomad_frequency.png"), "\n")
  }

  cat("\n")

  return(gnomad_plot)
}

# Plot TEs in dataset and their frequency in HostSeq
# Shows each unique TE and what percentage of HostSeq samples have it
#
# Args:
#   te_expand: TE expanded dataframe (one row per TE insertion)
#   output_dir: Directory to save plot
#   output_prefix: Prefix for output file name (e.g., "tumour" or "germline")
#
# Returns:
#   ggplot object
plot_te_hostseq_frequency <- function(te_expand, output_dir = NULL, output_prefix = "te") {

  cat("\n===== TEs IN DATASET AND THEIR HOSTSEQ FREQUENCY =====\n")

  # Identify HostSeq samples
  hostseq_samples <- unique(te_expand$sample[grepl("^HS_", te_expand$sample)])
  n_hostseq <- length(hostseq_samples)
  cat("Number of HostSeq samples:", n_hostseq, "\n")

  if (n_hostseq == 0) {
    cat("Warning: No HostSeq samples found. Cannot calculate HostSeq frequencies.\n")
    return(NULL)
  }

  # Get all unique TEs in the dataset
  te_unique <- te_expand %>%
    distinct(SV_chrom, SV_start, ALT, .keep_all = TRUE) %>%
    mutate(te_id = paste(SV_chrom, SV_start, ALT, sep = "_"))

  cat("Total unique TEs in dataset:", nrow(te_unique), "\n")

  # Calculate frequency in HostSeq for each TE
  te_hostseq_freq <- te_expand %>%
    filter(grepl("^HS_", sample)) %>%
    group_by(SV_chrom, SV_start, ALT) %>%
    summarise(
      n_hostseq_samples = n_distinct(sample),
      freq_in_hostseq = n_hostseq_samples / n_hostseq,
      .groups = "drop"
    ) %>%
    mutate(te_id = paste(SV_chrom, SV_start, ALT, sep = "_"))

  # Merge with all TEs (TEs not in HostSeq will have freq = 0)
  te_freq_complete <- te_unique %>%
    select(te_id, SV_chrom, SV_start, ALT) %>%
    left_join(te_hostseq_freq %>% select(te_id, n_hostseq_samples, freq_in_hostseq),
              by = "te_id") %>%
    mutate(
      n_hostseq_samples = ifelse(is.na(n_hostseq_samples), 0, n_hostseq_samples),
      freq_in_hostseq = ifelse(is.na(freq_in_hostseq), 0, freq_in_hostseq)
    ) %>%
    arrange(desc(freq_in_hostseq))

  cat("TEs found in HostSeq:", sum(te_freq_complete$n_hostseq_samples > 0), "\n")
  cat("TEs not in HostSeq:", sum(te_freq_complete$n_hostseq_samples == 0), "\n")

  # Summary statistics
  cat("\nHostSeq frequency summary:\n")
  cat("  Min:", min(te_freq_complete$freq_in_hostseq), "\n")
  cat("  Median:", median(te_freq_complete$freq_in_hostseq), "\n")
  cat("  Mean:", mean(te_freq_complete$freq_in_hostseq), "\n")
  cat("  Max:", max(te_freq_complete$freq_in_hostseq), "\n")

  # Create plot
  p <- ggplot(te_freq_complete, aes(x = freq_in_hostseq)) +
    geom_histogram(bins = 50, fill = "blue", color = "black", alpha = 0.7) +
    scale_x_continuous(labels = scales::percent,
                       breaks = seq(0, max(te_freq_complete$freq_in_hostseq), by = 0.1)) +
    labs(
      title = paste0("TE Frequency in HostSeq Cohort\n(", nrow(te_freq_complete),
                     " unique TEs, ", n_hostseq, " HostSeq samples)"),
      x = "Frequency in HostSeq",
      y = "Number of TEs"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.grid.minor = element_blank()
    ) +
    geom_vline(xintercept = 0.03, linetype = "dashed", color = "red", size = 1) +
    annotate("text", x = 0.03, y = Inf, label = "3% threshold",
             vjust = -0.5, hjust = -0.1, color = "red", size = 3.5)

  # Save plot
  if (!is.null(output_dir)) {
    ggsave(paste0(output_dir, output_prefix, "_te_hostseq_frequency.png"),
           plot = p, width = 10, height = 6)
    cat("✓ Plot saved to:", paste0(output_prefix, "_te_hostseq_frequency.png"), "\n")
  }

  return(p)
}

#### ANCESTRY PCA FUNCTIONS ####

# Merge location windows with ancestry data
# Args:
#   location_df: Data frame with location window counts (from 100kb_complete_filtered_g.csv or _t.csv)
#   ancestry_df: Ancestry data frame with predicted_ancestry_thres column
# Returns:
#   Merged data frame
merge_location_ancestry <- function(location_df, ancestry_df) {
  cat("\n===== MERGING LOCATION WINDOWS WITH ANCESTRY =====\n")

  # Check for sample column
  if (!"sample" %in% colnames(location_df)) {
    stop("location_df must have a 'sample' column")
  }
  if (!"sample" %in% colnames(ancestry_df)) {
    stop("ancestry_df must have a 'sample' column")
  }

  # Determine if this is tumor data (has _T suffix) or germline data (has _N suffix)
  sample_example <- location_df$sample[1]
  is_tumor <- grepl("_T$", sample_example)

  cat("Data type detected:", ifelse(is_tumor, "TUMOR", "GERMLINE"), "\n")
  cat("Example sample:", sample_example, "\n")

  # Prepare ancestry data for merging
  if (is_tumor) {
    # For tumor data: Remove _N from ancestry to create base_sample
    # Tumor samples like "0074_20-10579-A-02-00_T" have base_sample "0074"
    # Ancestry samples like "0074_N" become base_sample "0074"
    cat("Using base_sample matching for tumor data\n")

    # Create base_sample in ancestry (remove _N)
    ancestry_for_merge <- ancestry_df %>%
      mutate(base_sample = gsub("_N$", "", sample)) %>%
      select(-sample)

    # Create base_sample in location data if not present
    if (!"base_sample" %in% colnames(location_df)) {
      location_df <- location_df %>%
        mutate(base_sample = gsub("_.*_T$", "", sample))  # Extract base ID before tissue/flowcell info
    }

    # Merge on base_sample
    merged <- location_df %>%
      left_join(ancestry_for_merge, by = "base_sample")

  } else {
    # For germline data: Direct sample matching (both end in _N)
    cat("Using direct sample matching for germline data\n")

    # Merge - include both predicted_ancestry_thres and mapped_label
    ancestry_cols <- c("sample", "predicted_ancestry_thres")
    if ("mapped_label" %in% colnames(ancestry_df)) {
      ancestry_cols <- c(ancestry_cols, "mapped_label")
    }

    merged <- location_df %>%
      left_join(ancestry_df %>% select(all_of(ancestry_cols)), by = "sample")
  }

  cat("Location windows samples:", length(unique(location_df$sample)), "\n")
  cat("Ancestry samples:", length(unique(ancestry_df$sample)), "\n")
  cat("Merged samples with ancestry:", sum(!is.na(merged$predicted_ancestry_thres)), "\n")
  cat("Merged samples missing ancestry:", sum(is.na(merged$predicted_ancestry_thres)), "\n")

  # Show examples of unmatched samples
  if (sum(is.na(merged$predicted_ancestry_thres)) > 0) {
    unmatched <- unique(merged$sample[is.na(merged$predicted_ancestry_thres)])
    cat("\nFirst 10 unmatched samples:\n")
    print(head(unmatched, 10))
  }

  cat("====================================================\n\n")

  return(merged)
}


# Perform PCA on location window counts
# Args:
#   location_df: Data frame with location window counts
#   exclude_cols: Column names to exclude from PCA (clinical, ancestry, sample ID, etc.)
# Returns:
#   List with PCA results and transformed data
perform_location_pca <- function(location_df, exclude_cols = c("sample", "predicted_ancestry_thres", "mapped_label")) {
  cat("\n===== PERFORMING PCA ON LOCATION WINDOWS =====\n")

  # Identify location window columns (exclude clinical/ancestry/sample columns)
  all_cols <- colnames(location_df)
  location_cols <- setdiff(all_cols, exclude_cols)

  # Select only these columns
  location_data <- location_df[, location_cols, drop = FALSE]

  # Keep only numeric columns
  numeric_cols <- sapply(location_data, is.numeric)
  location_data <- location_data[, numeric_cols, drop = FALSE]

  # Remove columns that are all NA or constant
  location_data <- location_data[, sapply(location_data, function(x) {
    length(unique(x[!is.na(x)])) > 1
  }), drop = FALSE]

  cat("Total columns:", length(all_cols), "\n")
  cat("Numeric columns after filtering:", ncol(location_data), "\n")
  cat("Samples:", nrow(location_data), "\n")

  # Handle missing values - impute with column mean
  for (col in colnames(location_data)) {
    if (any(is.na(location_data[[col]]))) {
      location_data[[col]][is.na(location_data[[col]])] <- mean(location_data[[col]], na.rm = TRUE)
    }
  }

  # Perform PCA
  pca_result <- prcomp(location_data, center = TRUE, scale. = TRUE)

  # Create data frame with PCA results
  # Start with PC coordinates
  pca_df <- data.frame(
    sample = location_df$sample,
    PC1 = pca_result$x[, 1],
    PC2 = pca_result$x[, 2]
  )

  # Add back all excluded columns from original data (for coloring plots)
  # This includes both numeric and non-numeric clinical/ancestry columns
  excluded_cols_to_add <- setdiff(exclude_cols, "sample")  # Don't duplicate sample
  for (col in excluded_cols_to_add) {
    if (col %in% colnames(location_df)) {
      pca_df[[col]] <- location_df[[col]]
    }
  }

  # Calculate variance explained
  var_explained <- summary(pca_result)$importance[2, 1:2] * 100

  cat("PC1 variance explained:", round(var_explained[1], 2), "%\n")
  cat("PC2 variance explained:", round(var_explained[2], 2), "%\n")
  cat("===============================================\n\n")

  return(list(
    pca_result = pca_result,
    pca_df = pca_df,
    var_explained = var_explained
  ))
}

# Plot PCA colored by predicted ancestry or mapped label
# Args:
#   pca_result: List returned from perform_location_pca()
#   color_by: Column name to color by (e.g., "predicted_ancestry_thres" or "mapped_label")
#   output_dir: Directory to save plot
#   plot_prefix: Prefix for plot filename
# Returns:
#   ggplot object
plot_pca_ancestry <- function(pca_result, color_by = "predicted_ancestry_thres",
                              output_dir = NULL, plot_prefix = "germline") {
  pca_df <- pca_result$pca_df
  var_explained <- pca_result$var_explained

  # Check if color_by column exists
  if (!color_by %in% colnames(pca_df)) {
    stop(paste0("Column '", color_by, "' not found in PCA data frame"))
  }

  # Remove samples with missing values in color_by column
  pca_df_filtered <- pca_df %>% filter(!is.na(.data[[color_by]]))

  # Set title based on color_by
  plot_title <- if (color_by == "predicted_ancestry_thres") {
    "PCA of TE Location Windows Colored by Predicted Ancestry"
  } else if (color_by == "mapped_label") {
    "PCA of TE Location Windows Colored by Mapped Label"
  } else {
    paste0("PCA of TE Location Windows Colored by ", color_by)
  }

  legend_title <- if (color_by == "predicted_ancestry_thres") {
    "Predicted Ancestry"
  } else if (color_by == "mapped_label") {
    "Mapped Label"
  } else {
    color_by
  }

  p <- ggplot(pca_df_filtered, aes(x = PC1, y = PC2, color = .data[[color_by]])) +
    geom_point(size = 3, alpha = 0.7) +
    labs(
      title = plot_title,
      x = paste0("PC1 (", round(var_explained[1], 2), "% variance)"),
      y = paste0("PC2 (", round(var_explained[2], 2), "% variance)"),
      color = legend_title
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      legend.position = "right"
    )

  # Save plot
  if (!is.null(output_dir)) {
    filename <- paste0(plot_prefix, "_pca_", color_by, ".png")
    ggsave(paste0(output_dir, "pca_umap/", filename), plot = p, width = 10, height = 8)
    cat("✓ PCA plot saved to:", filename, "\n")
  }

  return(p)
}

# Plot pie chart of predicted ancestry
# Args:
#   ancestry_df: Data frame with predicted_ancestry_thres column (after merging with samples)
#   output_dir: Directory to save plot
#   plot_prefix: Prefix for plot filename
# Returns:
#   ggplot object
plot_ancestry_pie <- function(ancestry_df, output_dir = NULL, plot_prefix = "germline") {
  # Count ancestry categories
  ancestry_counts <- ancestry_df %>%
    filter(!is.na(predicted_ancestry_thres)) %>%
    group_by(predicted_ancestry_thres) %>%
    summarise(count = n(), .groups = "drop") %>%
    mutate(
      percentage = round(count / sum(count) * 100, 1),
      label = paste0(predicted_ancestry_thres, "\n", count, " (", percentage, "%)")
    )

  cat("\nAncestry distribution:\n")
  print(ancestry_counts)

  p <- ggplot(ancestry_counts, aes(x = "", y = count, fill = predicted_ancestry_thres)) +
    geom_bar(stat = "identity", width = 1, color = "white") +
    coord_polar("y", start = 0) +
    geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 4) +
    labs(
      title = "Distribution of Predicted Ancestry",
      fill = "Predicted Ancestry"
    ) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      legend.position = "right"
    )

  # Save plot
  if (!is.null(output_dir)) {
    ggsave(paste0(output_dir, "dataset/", plot_prefix, "_ancestry_pie.png"),
           plot = p, width = 10, height = 8)
    cat("✓ Ancestry pie chart saved to:", paste0(plot_prefix, "_ancestry_pie.png"), "\n")
  }

  return(p)
}

# Plot pie charts for HostSeq ancestry distribution by filter/analysis group
# Creates 3 pie charts showing ancestry distribution for:
#   1. All HostSeq samples
#   2. HostSeq filter group (used for determining common TEs)
#   3. HostSeq analysis group (used for final analysis)
#
# Args:
#   te_data: TE dataframe with hostseq_group and predicted_ancestry_thres columns
#   output_dir: Directory to save plots (optional)
#   plot_prefix: Prefix for plot filenames (default: "germline")
# Returns:
#   List with three ggplot objects (all, filter, analysis)
plot_hostseq_ancestry_pies <- function(te_data, output_dir = NULL, plot_prefix = "germline", ancestry_col = "predicted_ancestry_thres") {

  cat("\n===== PLOTTING HOSTSEQ ANCESTRY PIE CHARTS =====\n")
  cat("Using ancestry column:", ancestry_col, "\n")

  # Filter for HostSeq samples only
  hostseq_data <- te_data %>%
    filter(grepl("^HS_", sample)) %>%
    distinct(sample, .keep_all = TRUE)

  if (nrow(hostseq_data) == 0) {
    cat("Warning: No HostSeq samples found in data\n")
    return(NULL)
  }

  cat("Total HostSeq samples:", nrow(hostseq_data), "\n")

  # Helper function to create a single pie chart
  create_ancestry_pie <- function(data, title_suffix) {
    # Use the specified ancestry column dynamically
    ancestry_counts <- data %>%
      filter(!is.na(!!sym(ancestry_col))) %>%
      group_by(!!sym(ancestry_col)) %>%
      summarise(count = n(), .groups = "drop") %>%
      mutate(
        percentage = round(count / sum(count) * 100, 1),
        label = paste0(!!sym(ancestry_col), "\n", count, " (", percentage, "%)")
      )

    cat("\n", title_suffix, "- Ancestry distribution:\n", sep = "")
    print(ancestry_counts)

    p <- ggplot(ancestry_counts, aes(x = "", y = count, fill = !!sym(ancestry_col))) +
      geom_bar(stat = "identity", width = 1, color = "white") +
      coord_polar("y", start = 0) +
      geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 3.5) +
      labs(
        title = paste0("HostSeq Ancestry Distribution\n(", title_suffix, ")"),
        fill = ancestry_col
      ) +
      theme_void() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
        legend.position = "right"
      )

    return(p)
  }

  # Create pie chart for all HostSeq samples
  p_all <- create_ancestry_pie(hostseq_data, "All HostSeq")

  # Create pie chart for filter group
  filter_data <- hostseq_data %>% filter(hostseq_group == "filter")
  p_filter <- if (nrow(filter_data) > 0) {
    create_ancestry_pie(filter_data, "Filter Group")
  } else {
    cat("Warning: No HostSeq filter group samples found\n")
    NULL
  }

  # Create pie chart for analysis group
  analysis_data <- hostseq_data %>% filter(hostseq_group == "analysis")
  p_analysis <- if (nrow(analysis_data) > 0) {
    create_ancestry_pie(analysis_data, "Analysis Group")
  } else {
    cat("Warning: No HostSeq analysis group samples found\n")
    NULL
  }

  # Save plots
  if (!is.null(output_dir)) {
    dir.create(paste0(output_dir, "dataset/"), showWarnings = FALSE, recursive = TRUE)

    # Add suffix to filename if using mapped_label instead of default
    file_suffix <- if (ancestry_col == "mapped_label") "_mapped" else ""

    if (!is.null(p_all)) {
      ggsave(paste0(output_dir, "dataset/", plot_prefix, "_hostseq_ancestry_all", file_suffix, ".pdf"),
             plot = p_all, width = 7, height = 5)
      cat("✓ All HostSeq ancestry pie saved\n")
    }

    if (!is.null(p_filter)) {
      ggsave(paste0(output_dir, "dataset/", plot_prefix, "_hostseq_ancestry_filter", file_suffix, ".pdf"),
             plot = p_filter, width = 7, height = 5)
      cat("✓ Filter group ancestry pie saved\n")
    }

    if (!is.null(p_analysis)) {
      ggsave(paste0(output_dir, "dataset/", plot_prefix, "_hostseq_ancestry_analysis", file_suffix, ".pdf"),
             plot = p_analysis, width = 7, height = 5)
      cat("✓ Analysis group ancestry pie saved\n")
    }
  }

  cat("=====================================\n\n")

  return(list(
    all = p_all,
    filter = p_filter,
    analysis = p_analysis
  ))
}

# Perform UMAP on location window counts
# Args:
#   location_df: Data frame with location window counts
#   exclude_cols: Column names to exclude from UMAP (clinical, ancestry, sample ID, etc.)
# Returns:
#   List with UMAP results and transformed data
perform_location_umap <- function(location_df, exclude_cols = c("sample", "predicted_ancestry_thres", "mapped_label")) {
  cat("\n===== PERFORMING UMAP ON LOCATION WINDOWS =====\n")

  # Identify location window columns (exclude clinical/ancestry/sample columns)
  all_cols <- colnames(location_df)
  location_cols <- setdiff(all_cols, exclude_cols)

  # Select only these columns
  location_data <- location_df[, location_cols, drop = FALSE]

  # Keep only numeric columns
  numeric_cols <- sapply(location_data, is.numeric)
  location_data <- location_data[, numeric_cols, drop = FALSE]

  # Remove columns that are all NA or constant
  location_data <- location_data[, sapply(location_data, function(x) {
    length(unique(x[!is.na(x)])) > 1
  }), drop = FALSE]

  cat("Total columns:", length(all_cols), "\n")
  cat("Numeric columns after filtering:", ncol(location_data), "\n")
  cat("Samples:", nrow(location_data), "\n")

  # Handle missing values - impute with column mean
  for (col in colnames(location_data)) {
    if (any(is.na(location_data[[col]]))) {
      location_data[[col]][is.na(location_data[[col]])] <- mean(location_data[[col]], na.rm = TRUE)
    }
  }

  # Perform UMAP
  umap_result <- umap(location_data, n_neighbors = 25, min_dist = 0.3)

  # Create data frame with UMAP results
  # Start with UMAP coordinates
  umap_df <- data.frame(
    sample = location_df$sample,
    UMAP1 = umap_result$layout[, 1],
    UMAP2 = umap_result$layout[, 2]
  )

  # Add back all excluded columns from original data (for coloring plots)
  # This includes both numeric and non-numeric clinical/ancestry columns
  excluded_cols_to_add <- setdiff(exclude_cols, "sample")  # Don't duplicate sample
  for (col in excluded_cols_to_add) {
    if (col %in% colnames(location_df)) {
      umap_df[[col]] <- location_df[[col]]
    }
  }

  cat("UMAP completed successfully\n")
  cat("===============================================\n\n")

  return(list(
    umap_result = umap_result,
    umap_df = umap_df
  ))
}

# Plot UMAP colored by predicted ancestry or mapped label
# Args:
#   umap_result: List returned from perform_location_umap()
#   color_by: Column name to color by (e.g., "predicted_ancestry_thres" or "mapped_label")
#   output_dir: Directory to save plot
#   plot_prefix: Prefix for plot filename
# Returns:
#   ggplot object
plot_umap_ancestry <- function(umap_result, color_by = "predicted_ancestry_thres",
                               output_dir = NULL, plot_prefix = "germline") {
  umap_df <- umap_result$umap_df

  # Check if color_by column exists
  if (!color_by %in% colnames(umap_df)) {
    stop(paste0("Column '", color_by, "' not found in UMAP data frame"))
  }

  # Remove samples with missing values in color_by column
  umap_df_filtered <- umap_df %>% filter(!is.na(.data[[color_by]]))

  # Set title based on color_by
  plot_title <- if (color_by == "predicted_ancestry_thres") {
    "UMAP of TE Location Windows Colored by Predicted Ancestry"
  } else if (color_by == "mapped_label") {
    "UMAP of TE Location Windows Colored by Mapped Label"
  } else {
    paste0("UMAP of TE Location Windows Colored by ", color_by)
  }

  legend_title <- if (color_by == "predicted_ancestry_thres") {
    "Predicted Ancestry"
  } else if (color_by == "mapped_label") {
    "Mapped Label"
  } else {
    color_by
  }

  p <- ggplot(umap_df_filtered, aes(x = UMAP1, y = UMAP2, color = .data[[color_by]])) +
    geom_point(size = 3, alpha = 0.7) +
    labs(
      title = plot_title,
      x = "UMAP1",
      y = "UMAP2",
      color = legend_title
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      legend.position = "right"
    )

  # Save plot
  if (!is.null(output_dir)) {
    filename <- paste0(plot_prefix, "_umap_", color_by, ".png")
    ggsave(paste0(output_dir, "pca_umap/", filename), plot = p, width = 10, height = 8)
    cat("✓ UMAP plot saved to:", filename, "\n")
  }

  return(p)
}


# Generic wrapper to plot PCA colored by any variable
# Args:
#   pca_result: List returned from perform_location_pca()
#   color_by: Column name to color by (e.g., "tumor_type", "total", "cohort")
#   plot_title: Custom title for the plot (optional)
#   output_dir: Directory to save plot
#   plot_prefix: Prefix for plot filename
# Returns:
#   ggplot object
plot_pca_by_variable <- function(pca_result, color_by, plot_title = NULL,
                                 output_dir = NULL, plot_prefix = "germline") {
  # Check if color_by column exists in PCA data
  if (!color_by %in% colnames(pca_result$pca_df)) {
    cat("Warning: Column '", color_by, "' not found in PCA data frame - skipping plot\n")
    return(NULL)
  }

  # Use the existing plot_pca_ancestry function which already handles any color_by variable
  if (is.null(plot_title)) {
    # Generate title based on color_by if not provided
    plot_title <- paste0("PCA of TE Location Windows Colored by ", gsub("_", " ", color_by))
  }

  # Call the existing function but override the title
  p <- plot_pca_ancestry(pca_result, color_by = color_by,
                        output_dir = output_dir, plot_prefix = plot_prefix)

  # Update title if custom one provided
  p <- p + labs(title = plot_title)

  return(p)
}

# Generic wrapper to plot UMAP colored by any variable
# Args:
#   umap_result: List returned from perform_location_umap()
#   color_by: Column name to color by (e.g., "tumor_type", "total", "cohort")
#   plot_title: Custom title for the plot (optional)
#   output_dir: Directory to save plot
#   plot_prefix: Prefix for plot filename
# Returns:
#   ggplot object
plot_umap_by_variable <- function(umap_result, color_by, plot_title = NULL,
                                  output_dir = NULL, plot_prefix = "germline") {
  # Check if color_by column exists in UMAP data
  if (!color_by %in% colnames(umap_result$umap_df)) {
    cat("Warning: Column '", color_by, "' not found in UMAP data frame - skipping plot\n")
    return(NULL)
  }

  # Use the existing plot_umap_ancestry function which already handles any color_by variable
  if (is.null(plot_title)) {
    # Generate title based on color_by if not provided
    plot_title <- paste0("UMAP of TE Location Windows Colored by ", gsub("_", " ", color_by))
  }

  # Call the existing function but override the title
  p <- plot_umap_ancestry(umap_result, color_by = color_by,
                         output_dir = output_dir, plot_prefix = plot_prefix)

  # Update title if custom one provided
  p <- p + labs(title = plot_title)

  return(p)
}

# ============================================================================
# Regulatory Elements (RE) Analysis Functions
# ============================================================================

# Load RE report data and join to TE dataframe by ID
#
# This function loads regulatory elements data from an external report file
# and joins it to the input dataframe based on the ID column.
#
# @param df A dataframe with an ID column (typically TE data)
# @param re_report_path Path to the RE report file
# @return A dataframe with RE annotations joined
load_and_join_re_data <- function(df, re_report_path) {
  # Load RE report
  re_data <- fread(re_report_path, header = FALSE)

  # Add column names
  colnames(re_data) <- c("chr", "start", "end", "ins", "sample", "ID", "ref", "ALT", "x",
                         "filter", "info", "y", "genotype", "chr_reg", "start_reg", "end_reg",
                         "type_reg", "gene_reg")

  # Join to input dataframe by ID
  # Use left_join to keep all rows from df, or inner_join to keep only matches
  # Setting multiple = "all" to handle cases where one TE maps to multiple RE regions
  df_re <- df %>%
    inner_join(re_data, by = "ID", multiple = "all")

  cat("Loaded RE data with", nrow(re_data), "rows\n")
  cat("Joined to split dataframe:", nrow(df_re), "rows out of", nrow(df),
      "(", round(100 * nrow(df_re) / nrow(df), 1), "%)\n")

  return(df_re)
}

# Split gene_reg column into one gene per row
#
# This function takes a dataframe with a gene_reg column containing semicolon-
# separated gene names and expands it so that each gene gets its own row.
# This is necessary for pathway analysis which requires one gene per row.
#
# @param df A dataframe with a gene_reg column
# @return A dataframe with one row per gene
split_re_genes <- function(df) {
  # Count rows before splitting
  n_before <- nrow(df)

  # Split genes using separate_rows (same approach as existing Gene_name splitting)
  df_split <- df %>%
    separate_rows(gene_reg, sep = ";") %>%
    filter(!is.na(gene_reg) & gene_reg != "")

  # Count rows after splitting
  n_after <- nrow(df_split)

  cat("Split gene_reg column:", n_before, "rows ->", n_after, "rows\n")
  cat("Expansion factor:", round(n_after / n_before, 2), "x\n")
  cat("Unique genes:", length(unique(df_split$gene_reg)), "\n")

  return(df_split)
}

# Perform ORA comparing Affected vs Unaffected cancer status
#
# This function performs over-representation analysis to compare gene sets
# between samples with Affected vs Unaffected cancer status.
#
# @param df A dataframe with Gene_name and Cancer columns
# @return compareCluster ORA results object
# Perform ORA with custom cutoffs for cancer status (Affected vs Unaffected)
perform_ora_cancer_status_custom_cutoffs <- function(df, p_pathway = 0.05, q_pathway = 0.1, nsample_thresh = 0, gene_col = "Gene_name") {
  # Filter genes by sample threshold if specified
  if (nsample_thresh > 0) {
    # First, identify groups that meet the sample threshold
    samples_per_group <- df %>%
      group_by(Cancer) %>%
      summarise(n_samples = n_distinct(sample.x), .groups = "drop")

    groups_with_enough_samples <- samples_per_group %>%
      filter(n_samples >= nsample_thresh) %>%
      pull(Cancer)

    cat("Groups with >=", nsample_thresh, "samples:", paste(groups_with_enough_samples, collapse=", "), "\n")

    # If fewer than 2 groups meet threshold, return NULL
    if (length(groups_with_enough_samples) < 2) {
      cat("Insufficient groups (need at least 2 groups with >=", nsample_thresh, "samples)\n")
      return(NULL)
    }

    # Filter to only groups that meet the threshold
    df <- df %>% filter(Cancer %in% groups_with_enough_samples)

    cat("After filtering groups:", length(unique(df[[gene_col]])), "genes across", length(groups_with_enough_samples), "groups\n")
  }

  # Over-representation analysis comparing groups
  geneClusters_ora <- lapply(split(df[[gene_col]], df$Cancer), unique)

  # Perform ORA with compareCluster
  ora_cancer_status <- compareCluster(
    geneCluster = geneClusters_ora,
    fun = "enrichGO",
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    pvalueCutoff = p_pathway,
    qvalueCutoff = q_pathway
  )

  # Return the ORA results
  return(ora_cancer_status)
}

# Perform pathway analysis on regulatory elements genes
#
# This function performs over-representation analysis (ORA) on genes from
# regulatory elements, using the existing pathway analysis functions.
# It handles renaming gene_reg to Gene_name for compatibility.
#
# @param df_re_split A split dataframe with gene_reg column (one gene per row)
# @param analysis_type Type of analysis: "general", "tp53", "cancer", or "cancer_status"
# @return ORA results object
perform_re_pathway_analysis <- function(df_re_split, analysis_type = "general", min_samples = 0, p_pathway = 0.05, q_pathway = 0.1) {
  # Convert to data.frame if it's a tibble
  df_for_ora <- as.data.frame(df_re_split)

  # Determine which gene column to use
  if ("gene_reg" %in% colnames(df_for_ora)) {
    cat("Using gene_reg column from RE data\n")
    gene_column <- "gene_reg"
  } else if ("Gene_name" %in% colnames(df_for_ora)) {
    cat("Using Gene_name column\n")
    gene_column <- "Gene_name"
  } else {
    stop("Neither gene_reg nor Gene_name column found in dataframe. Available columns: ",
         paste(colnames(df_for_ora), collapse=", "))
  }

  cat("Performing", analysis_type, "pathway analysis on RE genes\n")
  cat("Unique genes for analysis:", length(unique(df_for_ora[[gene_column]])), "\n")

  # Perform appropriate ORA based on analysis type, passing gene column name
  ora_result <- switch(analysis_type,
    "general" = perform_ora_custom_cutoffs(df_for_ora, p_pathway = p_pathway, q_pathway = q_pathway, nsample_thresh = min_samples, filter_exon = FALSE, gene_col = gene_column),
    "tp53" = perform_ora_tp53_custom_cutoffs(df_for_ora, p_pathway = p_pathway, q_pathway = q_pathway, nsample_thresh = min_samples, gene_col = gene_column),
    "cancer" = perform_ora_cancer_custom_cutoffs(df_for_ora, p_pathway = p_pathway, q_pathway = q_pathway, nsample_thresh = min_samples, gene_col = gene_column),
    "cancer_status" = perform_ora_cancer_status_custom_cutoffs(df_for_ora, p_pathway = p_pathway, q_pathway = q_pathway, nsample_thresh = min_samples, gene_col = gene_column),
    "sample_type" = perform_ora_sample_type_custom_cutoffs(df_for_ora, p_pathway = p_pathway, q_pathway = q_pathway, nsample_thresh = min_samples, sample_type_column = "sample_type", gene_col = gene_column),
    stop("Unknown analysis_type. Use 'general', 'tp53', 'cancer', 'cancer_status', or 'sample_type'")
  )

  return(ora_result)
}

perform_ora_cancer_status <- function(df) {
  # Over-representation analysis comparing Affected vs Unaffected
  geneClusters_ora <- lapply(split(df$Gene_name, df$Cancer), unique)

  # Perform ORA with compareCluster
  ora_cancer_status <- compareCluster(
    geneCluster = geneClusters_ora,
    fun = "enrichGO",
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.1
  )

  # Return the ORA results
  return(ora_cancer_status)
}

# Test RNA expression differences by TE status at gene from RE report
#
# This function takes an RE report, identifies genes with sufficient TE occurrences,
# and tests whether RNA expression differs between samples with/without TEs at each gene.
#
# @param re_report_path Path to the RE report file
# @param rna_data RNA expression data (genes × samples with gene_name column)
# @param sample_type "germline" or "tumour" - affects sample ID parsing
# @param min_gene_mentions Minimum number of times a gene must appear in RE report (default: 3)
# @param min_samples_per_group Minimum samples required in each group for testing (default: 3)
# @return Dataframe with differential expression results
test_rna_by_gene_re_status <- function(re_report_path, rna_data, sample_type = "germline",
                                       min_gene_mentions = 3, min_samples_per_group = 3) {

  cat("\n=== Testing RNA Expression by TE Status from RE Report ===\n")
  cat("Sample type:", sample_type, "\n")
  cat("Minimum gene mentions:", min_gene_mentions, "\n")
  cat("Minimum samples per group:", min_samples_per_group, "\n\n")

  # Load RE report
  cat("Loading RE report from:", re_report_path, "\n")
  re_data <- fread(re_report_path, header = FALSE)

  # Add column names based on load_and_join_re_data function
  colnames(re_data) <- c("chr", "start", "end", "ins", "sample", "ID", "ref", "ALT", "x",
                         "filter", "info", "y", "genotype", "chr_reg", "start_reg", "end_reg",
                         "type_reg", "gene_reg")

  cat("Loaded", nrow(re_data), "rows from RE report\n")

  # Extract TE type from ALT column (e.g., <INS:ME:ALU> -> ALU)
  re_data$te_type <- gsub(".*:([^>]+)>", "\\1", re_data$ALT)

  # Create TE coordinates
  re_data$te_coordinates <- paste0(re_data$chr, ":", re_data$start, "-", re_data$end)

  # Create RE location
  re_data$te_location <- paste0(re_data$chr_reg, ":", re_data$start_reg, "-", re_data$end_reg)

  # Split gene_reg column (semicolon-separated) into one row per gene
  cat("\nSplitting gene_reg column...\n")
  re_split <- re_data %>%
    separate_rows(gene_reg, sep = ";") %>%
    filter(!is.na(gene_reg) & gene_reg != "")

  cat("After splitting:", nrow(re_split), "rows\n")

  # Count gene occurrences
  gene_counts <- re_split %>%
    count(gene_reg, name = "n_occurrences") %>%
    filter(n_occurrences >= min_gene_mentions) %>%
    arrange(desc(n_occurrences))

  cat("\nGenes with >=", min_gene_mentions, "mentions:", nrow(gene_counts), "\n")
  cat("Top 10 genes by occurrence:\n")
  print(head(gene_counts, 10))

  if (nrow(gene_counts) == 0) {
    cat("\nNo genes meet the minimum occurrence threshold.\n")
    return(data.frame())
  }

  # Extract sample IDs from ID column
  # Germline: HS_21-8063-A-02-00_N-1-825009-291-ALU -> HS_21-8063-A-02-00_N
  # Tumor: 0400_20-1851-A-02-00_T-xtea-1-827186-247-ALU -> 0400_20-1851-A-02-00_T
  cat("\nExtracting sample IDs from ID column...\n")

  # Pattern: sample_id-[tool]-chr-pos-len-type (tool is optional)
  # Remove: -[tool]-chr-pos-len-type from the end
  # Use perl=TRUE for extended regex with optional group
  re_split$sample_id <- gsub("-([a-z]+-)?[0-9XY]+-[0-9]+-[0-9]+-[A-Z0-9]+$", "", re_split$ID, perl = TRUE)

  cat("Extracted", length(unique(re_split$sample_id)), "unique samples from RE report\n")
  cat("Example sample IDs (first 5):\n")
  print(head(unique(re_split$sample_id), 5))

  # Get RNA sample names (excluding gene_name column)
  rna_samples <- setdiff(colnames(rna_data), "gene_name")
  cat("\nRNA data:", nrow(rna_data), "genes ×", length(rna_samples), "samples\n")

  # For germline RE report matched to tumor RNA: extract base patient IDs
  if (sample_type == "germline") {
    cat("\nMatching germline TE samples to tumor RNA samples by base patient ID...\n")

    # Extract base patient IDs (remove _N or _T suffix)
    re_base_ids <- unique(sub("_.*", "", re_split$sample_id))
    rna_base_ids <- sub("_.*", "", rna_samples)

    cat("Unique RE patient IDs:", length(re_base_ids), "\n")
    cat("Unique RNA patient IDs:", length(unique(rna_base_ids)), "\n")

    # Find overlap
    common_base_ids <- intersect(re_base_ids, unique(rna_base_ids))
    cat("Overlapping patient IDs:", length(common_base_ids), "\n")

    # Filter RNA samples to those with germline data
    matched_rna_samples <- rna_samples[rna_base_ids %in% re_base_ids]
    cat("Matched RNA samples:", length(matched_rna_samples), "\n")

    # Create mapping from base ID to RNA sample name
    base_to_rna <- setNames(matched_rna_samples, sub("_.*", "", matched_rna_samples))

  } else {
    # Tumor: direct matching only (no base ID matching for somatic samples)
    cat("\nMatching tumor TE samples to tumor RNA samples...\n")

    # Direct matching only
    matched_rna_samples <- intersect(unique(re_split$sample_id), rna_samples)
    cat("Matched samples:", length(matched_rna_samples), "\n")
  }

  # Test each gene
  results_list <- list()

  cat("\nTesting genes...\n")
  pb <- txtProgressBar(min = 0, max = nrow(gene_counts), style = 3)

  for (i in 1:nrow(gene_counts)) {
    setTxtProgressBar(pb, i)

    gene <- gene_counts$gene_reg[i]

    # Check if gene exists in RNA data
    if (!gene %in% rna_data$gene_name) {
      next
    }

    # Get samples with TE affecting this gene
    samples_with_te <- re_split %>%
      filter(gene_reg == gene) %>%
      pull(sample_id) %>%
      unique()

    # Map to RNA samples
    if (sample_type == "germline") {
      # Convert TE sample IDs to base IDs, then to RNA sample names
      base_ids_with_te <- unique(sub("_.*", "", samples_with_te))
      rna_samples_with_te <- base_to_rna[base_ids_with_te]
      rna_samples_with_te <- rna_samples_with_te[!is.na(rna_samples_with_te)]
    } else {
      # Tumor: direct matching only (no base ID matching for somatic samples)
      rna_samples_with_te <- intersect(samples_with_te, matched_rna_samples)
    }

    # Samples without TE = all matched RNA samples except those with TE
    rna_samples_without_te <- setdiff(matched_rna_samples, rna_samples_with_te)

    # Check minimum sample size
    if (length(rna_samples_with_te) < min_samples_per_group ||
        length(rna_samples_without_te) < min_samples_per_group) {
      next
    }

    # Extract RNA expression for this gene
    gene_expr <- rna_data %>%
      filter(gene_name == gene) %>%
      select(-gene_name) %>%
      unlist()

    expr_with_te <- gene_expr[rna_samples_with_te]
    expr_without_te <- gene_expr[rna_samples_without_te]

    # Remove zeros and NAs
    expr_with_te <- expr_with_te[!is.na(expr_with_te) & expr_with_te > 0]
    expr_without_te <- expr_without_te[!is.na(expr_without_te) & expr_without_te > 0]

    # Recheck sample sizes after filtering
    if (length(expr_with_te) < min_samples_per_group ||
        length(expr_without_te) < min_samples_per_group) {
      next
    }

    # Calculate statistics
    median_with <- median(expr_with_te, na.rm = TRUE)
    median_without <- median(expr_without_te, na.rm = TRUE)
    mean_with <- mean(expr_with_te, na.rm = TRUE)
    mean_without <- mean(expr_without_te, na.rm = TRUE)

    # Calculate fold change
    fold_change <- mean_with / mean_without
    effect_direction <- ifelse(mean_with > mean_without, "up", "down")

    # Perform Wilcoxon test
    test_result <- wilcox.test(expr_with_te, expr_without_te, alternative = "two.sided")

    # Get TE-specific information for this gene
    gene_te_info <- re_split %>%
      filter(gene_reg == gene) %>%
      summarise(
        te_types = paste(unique(te_type), collapse = ";"),
        te_coords = paste(unique(te_coordinates), collapse = ";"),
        te_location_types = paste(unique(type_reg), collapse = ";"),
        te_locations = paste(unique(te_location), collapse = ";"),
        n_unique_tes = n_distinct(ID)
      )

    # Store result
    results_list[[length(results_list) + 1]] <- data.frame(
      gene = gene,
      te_type = gene_te_info$te_types,
      te_coordinates = gene_te_info$te_coords,
      te_location_type = gene_te_info$te_location_types,
      te_location = gene_te_info$te_locations,
      gene_features = "",  # Placeholder for future annotations
      n_unique_tes = gene_te_info$n_unique_tes,
      samples_with_te = length(expr_with_te),
      samples_without_te = length(expr_without_te),
      median_with_te = median_with,
      median_without_te = median_without,
      mean_with_te = mean_with,
      mean_without_te = mean_without,
      fold_change = fold_change,
      effect_direction = effect_direction,
      p_value = test_result$p.value,
      stringsAsFactors = FALSE
    )
  }

  close(pb)

  # Combine results
  if (length(results_list) == 0) {
    cat("\nNo genes could be tested.\n")
    return(data.frame())
  }

  results_df <- bind_rows(results_list)

  # Calculate adjusted p-values
  results_df$p_adj <- p.adjust(results_df$p_value, method = "BH")

  # Sort by p-value
  results_df <- results_df %>%
    arrange(p_value)

  cat("\n\nResults summary:\n")
  cat("Total genes tested:", nrow(results_df), "\n")
  cat("Significant at p < 0.05:", sum(results_df$p_value < 0.05), "\n")
  cat("Significant at p_adj < 0.05:", sum(results_df$p_adj < 0.05), "\n")
  cat("Significant at p_adj < 0.1:", sum(results_df$p_adj < 0.1), "\n")

  if (nrow(results_df) > 0) {
    cat("\nTop 10 results by p-value:\n")
    print(results_df[1:min(10, nrow(results_df)), c("gene", "samples_with_te", "samples_without_te",
                                                      "fold_change", "effect_direction", "p_value", "p_adj")])
  }

  return(results_df)
}
