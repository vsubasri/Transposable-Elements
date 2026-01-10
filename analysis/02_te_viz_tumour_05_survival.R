#!/usr/bin/env Rscript

# Tumour TE Visualization - Survival Analysis
# Survival analysis and TE burden

# Source common setup and load data
source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/00_viz_common_setup.R")
REQUIRED_DATA <- c("count_matrix", "survival", "clinical")
source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/viz/00_viz_load_data_tumour.R")

# Initialize module-specific text output
init_module_sink(paste0(plot_dir, "survival_burden/"), "SURVIVAL")

cat("Running 02_te_viz_tumour_05_survival.R...\n")

cat("\n===== SURVIVAL ANALYSIS =====\n")
tryCatch({
  # Prepare survival data using custom function
  cat("Preparing survival data...\n")

  # Debug: check te_kics_t sample names
  cat("DEBUG: First 10 te_kics_t samples:\n")
  print(head(te_kics_t$sample, 10))
  cat("DEBUG: First 10 kics_DOD IDs:\n")
  print(head(kics_DOD$`KiCS ID`, 10))

  survival_data <- prepare_survival_data(te_kics_t, kics_DOD)

  write_output(quote(head(survival_data)), "Head of merged survival data")
  write_output(quote(table(survival_data$te_burden)), "TE burden group counts")
  write_output(quote(table(survival_data$event)), "Event counts (0=censored, 1=dead)")
  write_output(quote(summary(survival_data$time)), "Summary of survival times (days)")

  # Plot survival curves using custom function
  if (nrow(survival_data) > 0) {
    survival_results <- plot_survival_curves(survival_data, output_dir = plot_dir)

    if (!is.null(survival_results)) {
      write_output(quote(survival_results$logrank_test), "Log-rank test: High vs Low TE burden")
      write_output(quote(summary(survival_results$fit)), "Survival summary by TE burden")

      # Print survival plot to PDF
      cat("\nAdding survival plot to PDF...\n")
      titled_print(survival_results$plot, "Survival by TE Burden")

      # Save survival data with explanation
      survival_data_export <- survival_data %>%
        select(sample, base_sample, total, LINE1, ALU, SVA, te_burden,
               TP53_status, tumor_type, age_at_diagnosis,
               event, time, age_at_diagnosis_days) %>%
        mutate(
          time_calculation = case_when(
            event == 1 ~ "Age at death - Age at diagnosis",
            event == 0 ~ "Fixed 5-year follow-up (conservative estimate)",
            TRUE ~ NA_character_
          ),
          time_years = round(time / 365.25, 2)
        )

      survival_data_file <- paste0(plot_dir, "survival/survival_burden_data.csv")
      write.csv(survival_data_export, survival_data_file, row.names = FALSE)
      cat("✓ Survival data saved to:", basename(survival_data_file), "\n")
    }
  } else {
    cat("Warning: No samples with complete survival data available for analysis.\n")
  }
}, error = function(e) {
  cat("Warning: Could not perform survival analysis:", e$message, "\n")
})



cat("✓ Script completed successfully\n")

# Close module-specific sink
close_module_sink()
