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
color_palette_5 <- c("#DDD8C4", "#5FBFF9", '#99ff99', '#ffcc66', '#ff6666')
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
  # Rename the columns of the dataframe
  names(df) <- c("sample", "mean_cov", "sd_cov", "med_cov", "total_reads", 
                 "mean_read_len", "sd_read_len", "med_read_len", "pct_chimeras", "avg_quality")
  
  # Modify the 'sample' column for entries ending in _N
  df$sample <- ifelse(grepl("_N$", df$sample),
                      paste0(sub("_.*", "", df$sample), "_N"),
                      df$sample)
  
  return(df)
}
 
prep_ancestry_kics <-function(df) {
  # Ensure the column exists
  if (!"indivID" %in% colnames(df)) {
    stop("The dataframe does not have an 'indivID' column.")
  }
  
  # Create a new sample column based on indivID
  df$sample <- gsub("^KICS_|-N$", "", df$indivID) # Remove 'KICS_' and '-N'
  df$sample <- paste0(df$sample, "_N")            # Add '_N' at the end
  
  # Keep only sample and predicted_ancestry_thres columns
  df <- df[, c("sample", "predicted_ancestry_thres")]
  
  return(df)
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

filter_by_multiple_criteria <- function(df, df_nonproband, df_noconsent, df_hostseqcancer, metrics_df, column_names, thresholds, type) {
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
  
  # Step 1: Filter by df_nonproband
  initial_samples <- length(unique(df[[sample_id_col]]))
  df_filtered <- filter_by_another_df(df, df_nonproband)
  remaining_samples <- length(unique(df_filtered[[sample_id_col]]))
  filtered_out_nonproband <- initial_samples - remaining_samples
  cat("Filtered out", filtered_out_nonproband, "samples due to non proband normal\n")
  
  # Step 2: Filter by df_noconsent
  count_after_nonproband <- remaining_samples
  df_filtered <- filter_by_another_df(df_filtered, df_noconsent)
  remaining_samples <- length(unique(df_filtered[[sample_id_col]]))
  filtered_out_noconsent <- count_after_nonproband - remaining_samples
  cat("Filtered out", filtered_out_noconsent, "samples due to no consent\n")
  
  # Step 3: Filter by hostseq cancer
  count_after_noconsent <- remaining_samples
  df_filtered <- filter_by_another_df(df_filtered, df_hostseqcancer)
  remaining_samples <- length(unique(df_filtered[[sample_id_col]]))
  filtered_out_hostseqcancer <- count_after_noconsent - remaining_samples
  cat("Filtered out", filtered_out_hostseqcancer, "samples due to host seq having cancer\n")
  
  # Step 4: Filter by metrics_df using multiple columns and thresholds
  count_after_hostseqcancer <- remaining_samples
  df_filtered <- filter_by_metrics(df_filtered, metrics_df, column_names, thresholds)
  remaining_samples <- length(unique(df_filtered[[sample_id_col]]))
  filtered_out_metrics <- count_after_hostseqcancer - remaining_samples
  cat("Filtered out", filtered_out_metrics, "samples due to metrics with columns\n")
  
  return(df_filtered)
}

rename_alt <- function(df) {
  df <- df %>% rename("ALT" = "SV_type") %>% rename("sample" = "Samples_ID")
  
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
  }
  
  # Return the modified data frame
  return(df)
}

prep_te<- function(df, df_nonproband, df_noconsent, df_hostseqcancer, df_metrics, column_names, thresholds, type) {
  # Step 1: Rename alternative columns
  df <- rename_alt(df)
  
  # Step 2: Replace sample names to match the clinical data
  df <- replace_df_samples(df, type)
  
  # Step 3: Filter by nonproband, noconsent, and metrics criteria
  df <- filter_by_multiple_criteria(df, df_nonproband, df_noconsent, df_hostseqcancer, df_metrics, column_names, thresholds, type)
  
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

count_TE_occurrences <- function(te, rare_threshold_percentage) {
  # Step 1: Count unique samples
  total_samples <- length(unique(te$sample))
  
  # Step 2: Group by TE characteristics and count occurrences per sample
  te_counts <- te %>%
    group_by(SV_chrom, SV_start, SV_end, SV_length, SV_type) %>%
    summarize(unique_samples_in_group = n_distinct(sample)) %>%
    ungroup()
  
  # Step 3: Calculate thresholds for common and rare TEs
  rare_threshold <- (rare_threshold_percentage / 100) * total_samples
  common_threshold <- (1 - rare_threshold_percentage / 100) * total_samples
  
  # Step 4: Classify TEs
  num_common_TEs <- nrow(te_counts %>% filter(unique_samples_in_group >= common_threshold))
  num_rare_TEs <- nrow(te_counts %>% filter(unique_samples_in_group < rare_threshold))
  num_single_sample_TEs <- nrow(te_counts %>% filter(unique_samples_in_group == 1))
  
  # Print the results
  cat("Number of common TEs (>=", (1 - rare_threshold_percentage / 100) * 100, "% of samples):", num_common_TEs, "\n")
  cat("Number of rare TEs (<", rare_threshold_percentage, "% of samples):", num_rare_TEs, "\n")
  cat("Number of TEs in only one sample:", num_single_sample_TEs, "\n")
  cat("Spread of TEs in samples:", "\n")
  print(summary(te_counts$unique_samples_in_group))
}

count_TE_occurrences<- function(te, te_all, nohits, type="N") {
  # Step 1: Count unique samples
  num_samples_gt1 <- length(unique(te$sample))
  
  # Step 2: Group by TE characteristics and count occurrences per sample
  te_counts <- te %>%
    group_by(SV_chrom, SV_start, SV_end, SV_length, ALT) %>%
    summarize(unique_samples_in_group = n_distinct(sample), .groups = "drop" ) %>%
    ungroup()
  total_te <- nrow(te_counts)
  
  # single sample
  num_single_sample_TEs <- nrow(te_counts %>% filter(unique_samples_in_group == 1))
  perc_single_sample <- num_single_sample_TEs/total_te*100
  
  # nohits
  num_no_hits <- nrow(nohits)
  
  # at least 1 te
  num_total_samples <- num_samples_gt1 + num_no_hits 
  gt1 <- num_samples_gt1/ num_total_samples*100
  
  # number of patients with multiple samples
  if (type=="T") {
    te_sample_counts <- te %>%
      group_by(base_sample) %>%
      summarise(unique_samples_for_patient = n_distinct(sample), .groups = "drop") %>%
      filter(unique_samples_for_patient > 1) 
    
    # Report the number of such groups and the average count
    num_groups_gt1 <- nrow(te_sample_counts)
    average_count <- mean(te_sample_counts$unique_samples_for_patient)
  } 
  
  # Print the results
  cat("Total number of samples before filtering:", num_total_samples, "\n")
  cat("Number of samples with 0 TEs:", num_no_hits, "\n")
  cat("Percent of samples with at least 1 TE:", gt1, "\n")
  cat("Number of TEs in only one sample:", num_single_sample_TEs, "(", perc_single_sample, "%)\n")
  cat("Spread of number of TEs shared by samples (ie on average the same TE is in 23 samples:", "\n")
  print(summary(te_counts$unique_samples_in_group))
  cat("Spread of ALUs by sample:", "\n")
  print(summary(te_all$ALU))
  cat("Spread of LINE1s by sample:", "\n")
  print(summary(te_all$LINE1))
  cat("Spread of SVAs by sample:", "\n")
  print(summary(te_all$SVA))
}

plot_te_counts_summary <- function(df, log_scale=FALSE, breaks=NULL) {
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
    labs(x = "Repeat element", y = "Count") +
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
    scale_fill_manual(values = color_palette_5,
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
    scale_fill_manual(values = color_palette_5,
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
  # Create a complete list of unique samples
  complete_data <- data.frame(sample = unique(df$sample))
  
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
  chromosomes <- c(NA, 1:22, "X", "Y")
  
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
      # save result under specific label
      result$combination <- combination_label
      
      # Append the result to the list
      results_list <- append(results_list, list(result))
    }
  }
  
  # Combine all dataframes into one
  combined_result <- bind_rows(results_list)
  
  # Pivot the combined result to have combinations as columns
  wide_result <- combined_result %>%
    dplyr::select(sample, combination, count) %>%
    pivot_wider(names_from = combination, values_from = count, values_fill = list(count = 0))
  
  return(wide_result)
}

split_by_gene <- function(df){
  df_processed <- df %>%
    separate_rows(Gene_name, sep = ";") %>%
    filter(Gene_name != "")  # Remove rows where Gene_name is empty
  
  return(df_processed) 
}

# merge with clinical
merge_dfs <- function(df_te, df_info, include_all_x = FALSE, print_info = TRUE) {
  # Perform the merge with or without all.x based on the argument
  all <- merge(df_te, df_info, by = "sample", all.x = include_all_x)
  
  # Find samples in df_te that are not in df_info
  unmatched_merged <- unique(df_te[!df_te$sample %in% df_info$sample, ]$sample)
  
  if (print_info) {
    # Print unmatched samples 
    print("No info for these samples:")
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
  
  # Create a new data frame for the new samples with zero counts for TE columns
  new_data <- tibble(sample = new_samples)
  
  # Add columns with zero counts for each TE column in the original data (excluding 'sample')
  for (col in names(original_data)[-1]) {  # Skip the first column ('sample')
    new_data[[col]] <- 0
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
process_te_data_tumour <- function(te_raw, te_germline=NULL, clinical, nohits_list, metrics_list, 
                            apply_filter_common = TRUE, rare_gnomad = 3, rare_hostseq = 3, 
                            split_by_gene = FALSE,
                            apply_process_combinations = TRUE, 
                            select_samples_split = FALSE) {
  
  # Optionally filter common transposable elements
  if (apply_filter_common) {
    result <- filter_common_hostseq_tumour_te(te_raw, te_germline, rare_gnomad, rare_hostseq)
    
    # Check if the result is a list (meaning there are two values returned)
    if (is.list(result)) {
      te_filtered <- result$filtered_te_df  # Extract filtered TE dataframe
      missing_samples <- data.table(V1 = result$missing_samples)  # Extract missing samples, if any
      nohits_list <- rbind(nohits_list, missing_samples)
    } else {
      te_filtered <- result  # If only one value is returned, it's just the filtered TE dataframe
    }
  } else {
    te_filtered <- te_raw # dont filter common
  }
  
  # Optionally split by gene
  if (split_by_gene) {
    te_split_by_gene <- split_by_gene(te_filtered)
  } else{
    te_split_by_gene <- te_filtered
  } 
  
  # Optionally process all combinations
  if (apply_process_combinations) {
    te_processed <- process_all_combinations(te_split_by_gene)
    te_processed <- add_nohit_samples(te_processed, nohits_list)
  } else {
    te_processed <- te_split_by_gene
  }
  
  # merge ancestry
  #print("merging with ancestry")
  #te_clinical <- merge_dfs(te_processed, ancestry_list, include_all_x =TRUE, print_info=FALSE)
  
  # merge metrics
  print("merging with metrics")
  te_clinical <- merge_dfs(te_processed, metrics_list, include_all_x =TRUE, print_info=FALSE)
  
  # Merge with clinical data
  print("merging with clinical")
  te_clinical <- merge_dfs(te_clinical, clinical, include_all_x = FALSE, print_info=TRUE)
  
  # Remove samples where age at diagnosis is greater than 30
  te_filt <- filter_age(te_clinical)  # Assuming filter_age can accept a max_age parameter
  
  # remove duplicate entries
  te_format <- remove_duplicates(te_filt)
  
  # Set factor order
  te_format$TP53_status<- factor(te_format$TP53_status, levels = c("WT", "Mutant"))
  
  # Sort the data into categories
  te_aff_df <- te_format %>% filter(tumor_type != "U")  # Cancer samples
  te_lfs_df <- te_format %>%
    filter(TP53_status == "Mutant") %>%
    mutate(Cancer = ifelse(tumor_type != "U", "Affected", "Unaffected"))  # LFS samples, with cancer status
  te_kics_df <- te_aff_df %>% filter(cohort == "KiCS")  # KiCS cohort samples
  
  # Optionally select 1 sample for each patient 
  if (select_samples_split) {
    te_aff_samples <- select_samples(te_aff_df, select_samples_split)
    te_lfs_samples <- select_samples(te_lfs_df, select_samples_split)
    te_kics_samples <- select_samples(te_kics_df, select_samples_split)
    te_all_sampels <- select_samples(te_format, select_samples_split)
    
    result <- list(
      te_aff_selected = te_aff_samples$selected_samples,
      te_aff_all = te_aff_samples$all_samples,
      te_lfs_selected = te_lfs_samples$selected_samples,
      te_lfs_all = te_lfs_samples$all_samples,
      te_kics_selected = te_kics_samples$selected_samples,
      te_kics_all = te_kics_samples$all_samples,
      te_all_selected = te_all_sampels$selected_samples,
      te_all_all = te_all_sampels$all_samples
    )
    
  } else {
    # Combine the processed data into a list for easy access
    result <- list(
      te_all = te_format,
      te_aff= te_aff_df,
      te_lfs= te_lfs_df,
      te_kics= te_kics_df
    )
  }
  
  return(result)
}

process_te_data_germline <- function(te_raw, clinical, nohits_list, metrics_list, ancestry_list, 
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
  
  # Optionally split by gene
  if (split_by_gene) {
    te_split_by_gene <- split_by_gene(te_filtered)
  } else{
    te_split_by_gene <- te_filtered
  } 
  
  # Optionally process all combinations
  if (apply_process_combinations) {
    te_processed <- process_all_combinations(te_split_by_gene)
    te_processed <- add_nohit_samples(te_processed, nohits_list)
  } else {
    te_processed <- te_split_by_gene
  }
  
  # merge ancestry
  print("merging with ancestry")
  te_clinical <- merge_dfs(te_processed, ancestry_list, include_all_x =TRUE, print_info=FALSE)
  
  # merge metrics
  print("merging with metrics")
  te_clinical <- merge_dfs(te_clinical, metrics_list, include_all_x =TRUE, print_info=FALSE)
  
  # Merge with clinical data
  print("merging with clinical")
  te_clinical <- merge_dfs(te_clinical, clinical, include_all_x = FALSE, print_info=TRUE)
  
  # Remove samples where age at diagnosis is greater than 30
  te_filt <- filter_age(te_clinical)  # Assuming filter_age can accept a max_age parameter
  
  # remove duplicate entries
  te_format <- remove_duplicates(te_filt)
  
  # Set factor order
  te_format$TP53_status<- factor(te_format$TP53_status, levels = c("WT", "Mutant"))
  
  # Sort the data into categories
  te_aff_df <- te_format %>% filter(tumor_type != "U")  # Cancer samples
  te_lfs_df <- te_format %>%
    filter(TP53_status == "Mutant") %>%
    mutate(Cancer = ifelse(tumor_type != "U", "Affected", "Unaffected"))  # LFS samples, with cancer status
  te_kics_df <- te_aff_df %>% filter(cohort == "KiCS")  # KiCS cohort samples
  
  # Combine the processed data into a list for easy access
  result <- list(
    te_all = te_format,
    te_aff= te_aff_df,
    te_lfs= te_lfs_df,
    te_kics= te_kics_df
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
  setDT(te_germline)
  
  # Add cohort column to te_germline
  te_germline[, cohort_hs := ifelse(grepl("^HS_", sample), "HostSeq", "Other")]
  
  # Convert rare thresholds to decimals
  rare_threshold_gnomad <- rare_gnomad_threshold / 100
  cat("Rare gnomad threshold:", rare_threshold_gnomad, "\n")
  
  rare_threshold_hostseq <- rare_hostseq_threshold / 100
  cat("Rare hostseq threshold:", rare_threshold_hostseq, "\n")
  
  # Count unique TEs in te_df
  unique_te_df_count <- uniqueN(te_df, by = c("SV_chrom", "SV_start", "SV_end", "SV_length", "ALT"))
  cat("Total # of unique TEs in te_df:", unique_te_df_count, "\n")
  
  # Count unique samples in te_germline for the "HostSeq" cohort
  total_samples_hostseq <- uniqueN(te_germline[cohort_hs == "HostSeq", sample])
  cat("Total # of samples in HostSeq cohort:", total_samples_hostseq, "\n")
  
  # Group te_germline to calculate filtering flags
  te_with_flags_germline <- te_germline[
    cohort_hs == "HostSeq", 
    .(
      unique_samples_in_group_hostseq = uniqueN(sample) # Count unique samples in HostSeq
    ),
    by = .(SV_chrom, SV_start, SV_end, SV_length, ALT)
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
  
  # Remove NA indices
  overlapping_indices <- na.omit(overlapping_indices)
  
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
  
  # Print filtered stats
  cat("Filtered by common TEs:", num_filtered_by_common, "(", round(percent_filtered_by_common, 2), "%)\n")
  cat("Filtered by GRPMAX_AF:", num_filtered_by_AF, "(", round(percent_filtered_by_AF, 2), "%)\n")
  
  # Check for missing samples after filtering
  original_samples <- unique(te_df$sample)
  remaining_samples <- unique(filtered_te_df$sample)
  missing_samples <- setdiff(original_samples, remaining_samples)
  
  if (length(missing_samples) > 0) {
    cat("Samples lost after filtering common TEs:", missing_samples, "\n")
    
    return(list(filtered_te_df = filtered_te_df, missing_samples = missing_samples))
    
  } else {
    cat("No samples were lost after filtering common TEs.\n")
    
    return(filtered_te_df)
  }
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
  
  # total samples for hostseq 
  total_samples_hostseq <- length(unique(te$sample[te$cohort_hs == "HostSeq"]))
  cat("Total # samples for cohort 'HostSeq':", total_samples_hostseq, "\n")
  
  # total tes
  total_te_df <- nrow(te)
  cat("Total # TEs called in df:", total_te_df, "\n")
    
  # Group by TE characteristics and calculate thresholds for filtering
  te_with_flags <- te %>%
    group_by(SV_chrom, SV_start, SV_end, SV_length, ALT) %>%
    summarise(
      # Count unique samples from HostSeq cohort
      unique_samples_in_group_hostseq = n_distinct(sample[cohort_hs == "HostSeq"]),
      
      # Calculate the filtering conditions
      exceeds_AF_threshold = any(GRPMAX_AF >= rare_threshold_gnomad),  # At least one sample exceeds the AF threshold in the group
      is_common_hostseq = unique_samples_in_group_hostseq / total_samples_hostseq >= rare_threshold_hostseq,  # Common condition based on HostSeq samples
      .groups = "keep"
    ) 
  
  # total unique TEs
  total_unique_te <- nrow(te_with_flags)
  cat("Total # of unique TEs:", total_unique_te, "\n")
  
  # Calculate the counts for each filter condition
  num_filtered_by_AF <- sum(te_with_flags$exceeds_AF_threshold)
  num_filtered_by_common <- sum(te_with_flags$is_common_hostseq)
  num_overlap_filters <- sum(te_with_flags$is_common_hostseq & te_with_flags$exceeds_AF_threshold) # overlap between two filters
  
  # Filter out rows that are common or exceed the AF threshold
  uncommon_TEs <- te %>%
    left_join(te_with_flags %>% select(SV_chrom, SV_start, SV_end, SV_length, ALT, is_common_hostseq, exceeds_AF_threshold), 
              by = c("SV_chrom", "SV_start", "SV_end", "SV_length", "ALT")) %>%  # Join to keep the flags in the original df
    filter(is_common_hostseq == FALSE, exceeds_AF_threshold == FALSE) %>%  # Keep only uncommon and below threshold
    dplyr::select(-is_common_hostseq, -exceeds_AF_threshold)  # Remove helper columns
  
  # Calculate the percentage filtered by each condition and by the overlap
  percent_filtered_by_AF <- (num_filtered_by_AF / total_unique_te) * 100
  percent_filtered_by_common <- (num_filtered_by_common / total_unique_te) * 100
  percent_overlap_filters <- (num_overlap_filters / total_unique_te) * 100
  
  # Print messages
  cat("Filtered by GRPMAX_AF:", num_filtered_by_AF, "(", round(percent_filtered_by_AF, 2), "%)\n")
  cat("Filtered by is_common_hostseq:", num_filtered_by_common, "(", round(percent_filtered_by_common, 2), "%)\n")
  cat("Filtered by both (overlap):", num_overlap_filters, "(", round(percent_overlap_filters, 2), "%)\n")
  
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
  lesion_rank <- c("primary" = 1, "metastasis" = 2)
  disease_rank <- c("initial" = 1, "progressive" = 2, "relapsed" = 3)
  
  # Handle NA as the lowest priority
  lesion_rank <- c(lesion_rank, "NA" = max(lesion_rank) + 1)
  disease_rank <- c(disease_rank, "NA" = max(disease_rank) + 1)
  
  # Assign lesion_rank
  df$lesion_rank <- ifelse(
    is.na(df$lesion_type) | df$lesion_type == "NA" | !(df$lesion_type %in% names(lesion_rank)),
    lesion_rank["NA"],
    lesion_rank[df$lesion_type]
  )
  # Assign disease_rank
  df$disease_rank <- ifelse(
    is.na(df$disease_state) | df$disease_state == "NA" | !(df$disease_state %in% names(disease_rank)),
    disease_rank["NA"],
    disease_rank[df$disease_state]
  )
  
  # Identify invalid lesion_type and disease_state values and their rows
  invalid_lesion_rows <- df[!(is.na(df$lesion_type)) & !(df$lesion_type %in% names(lesion_rank)), ]
  invalid_disease_rows <- df[!(is.na(df$disease_state)) & !(df$disease_state %in% names(disease_rank)), ]
  
  # Print warning for invalid lesion_type values
  if (nrow(invalid_lesion_rows) > 0) {
    warning("The following lesion_type values and corresponding samples do not match any keys in lesion_rank:\n")
    print(invalid_lesion_rows[, c("sample", "lesion_type")])  # Print only sample and lesion_type columns
  }
  
  # Print warning for invalid disease_state values
  if (nrow(invalid_disease_rows) > 0) {
    warning("The following disease_state values and corresponding samples do not match any keys in disease_rank:\n")
    print(invalid_disease_rows[, c("sample", "disease_state")])  # Print only sample and disease_state columns
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
  
  # Print p-value from Wilcoxon test
  print("dependent ~ independent")
  formula <- reformulate(group, combination)
  print(formula)
  p_value <- wilcox.test(formula, data = df)$p.value
  print(paste("Wilcoxon test p-value:", p_value))
  
  # Format the p-value for display
  p_value_formatted <- formatC(p_value, format = "e", digits = 2)
  
  # Calculate and print medians for each group
  medians <- df %>%
    group_by(!!sym(group)) %>%
    summarise(median_count = median(!!sym(combination), na.rm = TRUE))
  print(medians)
  
  # Plot with optional log scale
  plot <- plot_box(df, combination, group, x_lab, y_lab, colours, log_scale) + 
    geom_signif(test = "wilcox.test", 
                comparisons = list(levels(factor(df[[group]]))),
                map_signif_level = TRUE,
                textsize = 5) + 
    coord_cartesian(clip = 'off') # dont cut off annotation
  
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

plot_count_lm <- function(df, chr, type, group, covariates, x_lab, y_lab, residuals = FALSE, log_scale = FALSE, breaks=NULL, min_samples = 5) {
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
    dummy_var <- grep(paste0("^", group), rownames(coef(lm_summary)), value = TRUE)
    group_p_value <- coef(lm_summary)[dummy_var, "Pr(>|t|)"]
    print(paste("p-value for", group, ":", group_p_value))
    
    # Format the p-value for display
    p_value_formatted <- formatC(group_p_value, format = "e", digits = 2)
    
    # Calculate and print medians for the group
    medians <- df %>%
      group_by(.data[[group]]) %>%
      summarise(median_count = median(!!sym(combination), na.rm = TRUE))
    print(medians)
    
    # Plot the group effect
    p <- ggplot(df, aes(x = .data[[group]], y = !!sym(combination), fill = .data[[group]])) +
      geom_boxplot(outlier.shape = NA, color = "black") +
      scale_fill_manual(values = colours_3, guide = "none") +
      geom_jitter(color = "black", size = 1.5, width = 0.2) +
      labs(x = x_lab, y = y_lab) +
      annotate("text", x = 1.7, y = max(df[[combination]], na.rm = TRUE),
               label = paste("p =", p_value_formatted), hjust = 0) 
    
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
    geom_vline(xintercept = observed_stat, color = "red", size = 1, linetype = "solid") +
    labs(title = "Bootstrap Distribution of Test Statistics",
         x = "Test Statistic",
         y = "Frequency")
  
  print(p_bootstrap)
  
  # Plot the p-value convergence
  convergence_data <- data.frame(bootstrap_sizes = bootstrap_sizes, p_values = p_value_convergence)
  p_convergence <- ggplot(convergence_data, aes(x = bootstrap_sizes, y = p_values)) +
    geom_line(color = "blue", size = 1) +
    geom_point(color = "blue") +
    labs(title = "Convergence of P-value with Number of Bootstraps",
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
  
  # Perform Kruskal-Wallis test
  test_result <- kruskal.test(reformulate(group, combination), data = df)
  p_value <- test_result$p.value
  p_value_formatted <- formatC(p_value, format = "e", digits = 2)
  
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
    scale_fill_manual(values =  colour_palette_3) + #mutation_colours
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
  
  # Check if there are at least two groups
  if (length(unique(data[[group]])) < 2) {
    print(paste("Not enough groups for gene:", gene, "and location:", location, "num samples:", nrow(data)))
    print(data[[group]])
    return(NA)  # Return NA if not enough groups
  }
  
  # Only return p-value from Wilcox test
  return(wilcox.test(reformulate(group, "count"), data = data, exact=FALSE)$p.value)
}

count_location_wilcox <- function(te_df, gene_vector, filter_var="SV.type", filter_element, group, location) {
  # Initialize a dataframe to store the results
  results <- data.frame(gene = character(), p_value = numeric(), stringsAsFactors = FALSE)
  
  # Loop through each gene
  for (i in 1:length(gene_vector)) {
    p_value <- calc_location_wilcox(te_df, filter_element=filter_element, group=group, location=location, gene=gene_vector[i])
    
    if (is.na(p_value)) {
      next  # Skip this iteration if the combination is not present
    }
    
    results <- rbind(results, data.frame(gene = gene_vector[i], p_value = p_value))
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
plot_count_perchr <- function(df, type = NA, group, y_lab) {
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
    
    plot_data <- plot_data + annotate("text", x = x_pos, y = y_pos, label = p_label, size = 16/.pt, vjust = -0.5)
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
plot_count_perchr_notest <- function(df, type=NA, y_lab, log_scale=FALSE) {
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
    mutate(chr = str_remove(chr, "^chr")) %>%
    mutate(chr = str_remove(chr, "_.*$"))
  
  # Normalize the counts
  data <- data %>%
    left_join(chr_length, by = "chr") %>%
    mutate(normalized_count = count / length)
  
  # Set factor levels for chromosome ordering
  data$chr <- factor(data$chr, levels = chr_length$chr)
  
  # Handle log scale by adding epsilon to normalized_count
  if (log_scale == TRUE) {
    epsilon <- 1e-6  # Small constant to handle zeros
    data <- data %>%
      mutate(normalized_count = normalized_count + epsilon)  # Add epsilon to avoid log(0)
  }
  
  # Plot using the custom plot function
  plot_data <- plot_box_perchr_notest(data, "chr", "normalized_count", y_lab)
  
  # Add log scale to the plot if requested
  if (log_scale == TRUE) {
    plot_data <- plot_data +
      scale_y_log10(labels = scales::label_log(base = 10)) +
      labs(y = paste("Log10(", y_lab, " + ", epsilon, ")"))  # Update Y-axis label
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
  
  # Perform Kruskal-Wallis test
  formula <- reformulate(combination, column)
  print("dependent ~ independent")
  print(formula)
  p_value <- kruskal.test(formula, data = df_filtered)$p.value
  
  # Format the p-value for display
  p_value_formatted <- formatC(p_value, format = "e", digits = 2)
  
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
    data_to_plot$Category<- factor(data_to_plot$Category, levels = c("WT", "Mutant"))
  }
  
  if (column_name == "sex") {
    category_mapping <- c("F" = "Female", "M" = "Male")
    
    # Map abbreviations to full names
    data_to_plot$Category <- category_mapping[as.character(data_to_plot$Category)]
    
    # Factorize with desired levels
    data_to_plot$Category <- factor(data_to_plot$Category, levels = c("Female", "Male"))
  }
  
  if (column_name == "tumor_type") {
    data_to_plot <- data_to_plot %>% filter(Category!="U") # remove no cancer
  }
    
  # Calculate percentage for labels
  data_to_plot$Percentage <- round((data_to_plot$Freq / total) * 100, 1)
  data_to_plot$Label <- paste0(data_to_plot$Category, ": ", data_to_plot$Percentage, "%")
  
  # Plot
  p <- ggplot(data_to_plot, aes(x = "", y = Freq, fill = Category)) +
    geom_bar(width = 1, stat = "identity") +
    coord_polar("y", start = 0) + 
    theme_void() +
    theme(legend.position = "right")
  
  # Define color scheme based on column_name
  if (column_name == "TP53_status") {
    p <- p +scale_fill_manual(values = colours, name = ylab, labels = data_to_plot$Label)
  } else if (column_name == "tumor_type") {
    default_colors <- scales::hue_pal()(nrow(data_to_plot) - 1) # Generate default colors
    color_mapping <- setNames(c(default_colors, "grey"), 
                              c(data_to_plot$Category[data_to_plot$Category != "Other"], "Other"))
    p <- p + scale_fill_manual(values = color_mapping, name = ylab, labels = data_to_plot$Label)
  } else if (column_name == "sex") {
    p <- p + scale_fill_manual(values = colours_sex, name = ylab, labels = data_to_plot$Label)
  } else if (column_name == "new_cohort") {
    p <- p + scale_fill_manual(values = colour_palette_4, name=ylab)
  }
  
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
  # Merge te_aff_split with gene_size where Gene_name matches name2
  merged_df <- df_te %>%
    left_join(df_size, by = c("Gene_name" = "name2"))
  
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
  
  # Perform Kruskal-Wallis test
  test_result <- kruskal.test(reformulate(group, combination), data = df)
  p_value <- test_result$p.value
  p_value_formatted <- formatC(p_value, format = "e", digits = 2)
  
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
      length = str_extract(INFO, "AVG_LEN=[^;]+") %>% str_replace("AVG_LEN=", "")
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


