#!/usr/bin/env Rscript

# Tumour TE Visualization - SV Overlap
# Structural variant - TE overlap analysis

# Source common setup and load data
source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/00_viz_common_setup.R")
source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/00_viz_load_data_tumour.R")

cat("Running 02_te_viz_tumour_11_sv_overlap.R...\n")

cat("\n===== STRUCTURAL VARIANT - TE OVERLAP ANALYSIS =====\n")
tryCatch({
  # Print SV dataset info
  cat("SV dataset summary:\n")
  cat("  Total SVs:", nrow(tumour_sv), "\n")
  cat("  Unique samples in SV data:", length(unique(tumour_sv$SampleId)), "\n")
  write_output(quote(length(unique(tumour_sv$SampleId))), "Number of unique samples in SV dataset")

  # Identify SVs with both breakpoints in TEs (connecting two TEs)
  cat("\nIdentifying SVs with both breakpoints in TEs...\n")
  sv_te_both <- identify_sv_in_te_linx(tumour_sv, te_aff_expand_all_t)

  write_output(quote(nrow(sv_te_both)), "Number of SVs with both breakpoints in TEs")

  # Identify samples with TE count > 0 but no SV data
  # Get samples with any TE count > 0 from count matrix
  te_count_cols <- grep("^(LINE|SINE|SVA|ERV|DNA|Retroposon|RNA|Other)", names(te_aff_t), value = TRUE)
  te_aff_t$total_te_count <- rowSums(te_aff_t[, te_count_cols, drop = FALSE], na.rm = TRUE)
  samples_with_te_count <- te_aff_t %>%
    filter(total_te_count > 0) %>%
    pull(sample) %>%
    sub("_T$", "", .)  # Remove _T suffix to match SV sample IDs

  # Get samples in SV data
  samples_sv <- unique(tumour_sv$SampleId)

  # Find samples with TE count > 0 but no SV data
  # Apply conversion strategies to check if TE samples match SV samples

  # Hardcoded sample mappings (TE sample -> SV sample)
  hardcoded_te_to_sv <- c(
    "U02H2D_A" = "U02H2D_5524A",
    "5009_1" = "5009_4856_1",
    "5471_2" = "5471_5787_2",
    "621_1" = "621_3311A_1",
    "PD13489" = "PD13489_PD13489a"
  )

  samples_te_matched <- sapply(samples_with_te_count, function(te_sample) {
    # Strategy 0: Check hardcoded mappings first
    if (te_sample %in% names(hardcoded_te_to_sv)) {
      sv_sample <- hardcoded_te_to_sv[[te_sample]]
      if (sv_sample %in% samples_sv) return(TRUE)
    }

    # Strategy 1: Exact match
    if (te_sample %in% samples_sv) return(TRUE)

    # Strategy 2: Convert all _ to -
    te_all_hyphen <- gsub("_", "-", te_sample)
    if (te_all_hyphen %in% samples_sv) return(TRUE)

    # Strategy 3: Keep first _, convert rest to -
    parts <- strsplit(te_sample, "_", fixed = TRUE)[[1]]
    if (length(parts) > 2) {
      te_first_underscore <- paste0(parts[1], "_", paste(parts[-1], collapse = "-"))
      if (te_first_underscore %in% samples_sv) return(TRUE)
    }

    # Strategy 4: Flexible regex - check if any SV sample matches the pattern
    sample_pattern <- gsub("[-_]", "[-_]", te_sample)
    if (any(grepl(paste0("^", sample_pattern, "$"), samples_sv))) return(TRUE)

    return(FALSE)
  })

  samples_te_no_sv <- samples_with_te_count[!samples_te_matched]

  if (length(samples_te_no_sv) > 0) {
    # Get TE counts for these samples
    samples_te_no_sv_with_T <- paste0(samples_te_no_sv, "_T")
    samples_te_no_sv_df <- te_aff_t %>%
      filter(sample %in% samples_te_no_sv_with_T) %>%
      select(sample, total_te_count, all_of(te_count_cols)) %>%
      arrange(desc(total_te_count))

    samples_te_no_sv_file <- paste0(r_dir_files, "samples_te_no_sv.csv")
    write.csv(samples_te_no_sv_df, samples_te_no_sv_file, row.names = FALSE)
    cat("✓ Samples with TE count > 0 but no SV data saved to:", basename(samples_te_no_sv_file), "\n")
    cat("\nSamples with TE count > 0 but no SV data (n=", length(samples_te_no_sv), "):\n", sep="")
    write_output(quote(samples_te_no_sv_df), "Samples with TE count > 0 but no SV data")
  } else {
    cat("  All samples with TE count > 0 have SV data available\n")
  }

  if (nrow(sv_te_both) > 0) {
    write_output(quote(head(sv_te_both, 20)), "SVs with both breakpoints in TEs (first 20)")
    # Create detailed summary with all SV and both TE info
    # Note: TE.type is the same for both TEs (from merge), need to get separate types
    sv_te_both_summary <- sv_te_both %>%
      mutate(
        TE1.start = TE.start.one,
        TE1.end = TE.end.one,
        TE1.type = TE.type,  # This is the shared TE type from the merge
        TE2.start = TE.start.two,
        TE2.end = TE.end.two,
        TE2.type = TE.type   # Same type for both since merge was on TE.type
      ) %>%
      select(sample, SV.chrom, SV.start, SV.end, SV.type,
             TE1.start, TE1.end, TE1.type,
             TE2.start, TE2.end, TE2.type) %>%
      arrange(sample, SV.chrom, SV.start)

    # Add clinical info - match by base sample
    sv_te_both_summary <- sv_te_both_summary %>%
      mutate(base_sample = sub("_.*", "", sample)) %>%
      left_join(clinical %>%
                  mutate(base_sample = sub("_.*", "", sample)) %>%
                  select(base_sample, TP53_status, tumor_type, sex, age_at_diagnosis) %>%
                  distinct(base_sample, .keep_all = TRUE),
                by = "base_sample") %>%
      select(-base_sample)  # Remove the temporary column

    # Add n_samples column: count how many samples have the same TE1 coordinates
    sv_te_both_summary <- sv_te_both_summary %>%
      group_by(SV.chrom, TE1.start, TE1.end, TE1.type) %>%
      mutate(n_samples = n_distinct(sample)) %>%
      ungroup()

    sv_te_both_file <- paste0(r_dir_files, "sv_te_both_breakpoints.csv")
    write.csv(sv_te_both_summary, sv_te_both_file, row.names = FALSE)
    cat("✓ SVs with both breakpoints in TEs saved to:", basename(sv_te_both_file), "\n")
  } else {
    cat("  No SVs found with both breakpoints in TEs\n")
  }

  # Identify SVs with at least one breakpoint in TE
  cat("\nIdentifying SVs with at least one breakpoint in TE...\n")
  sv_te_onebreak <- identify_sv_in_te_linx_onebreak(tumour_sv, te_aff_expand_all_t)

  write_output(quote(nrow(sv_te_onebreak)), "Number of SVs with at least one breakpoint in TE")

  if (nrow(sv_te_onebreak) > 0) {
    write_output(quote(head(sv_te_onebreak, 20)), "SVs with at least one breakpoint in TE (first 20)")
    # Create summary for one breakpoint in TE
    sv_te_onebreak_summary <- sv_te_onebreak %>%
      select(sample, SV.chrom, SV.start, SV.end, SV.type,
             TE1.start = TE.start, TE1.end = TE.end, TE1.type = TE.type) %>%
      arrange(sample, SV.chrom, SV.start)

    # Add clinical info - match by base sample
    sv_te_onebreak_summary <- sv_te_onebreak_summary %>%
      mutate(base_sample = sub("_.*", "", sample)) %>%
      left_join(clinical %>%
                  mutate(base_sample = sub("_.*", "", sample)) %>%
                  select(base_sample, TP53_status, tumor_type, sex, age_at_diagnosis) %>%
                  distinct(base_sample, .keep_all = TRUE),
                by = "base_sample") %>%
      select(-base_sample)  # Remove the temporary column

    # Add n_samples column: count how many samples have the same TE coordinates
    sv_te_onebreak_summary <- sv_te_onebreak_summary %>%
      group_by(SV.chrom, TE1.start, TE1.end, TE1.type) %>%
      mutate(n_samples = n_distinct(sample)) %>%
      ungroup()

    write_output(quote(head(sv_te_onebreak_summary, 20)), "SV-TE overlaps with one breakpoint (first 20)")

    sv_te_onebreak_file <- paste0(r_dir_files, "sv_te_onebreak.csv")
    write.csv(sv_te_onebreak_summary, sv_te_onebreak_file, row.names = FALSE)
    cat("✓ SVs with one breakpoint in TE saved to:", basename(sv_te_onebreak_file), "\n")
  } else {
    cat("  No SVs found with one breakpoint in TE\n")
  }
}, error = function(e) {
  cat("Warning: Could not perform SV-TE overlap analysis:", e$message, "\n")
})



cat("✓ Script completed successfully\n")
