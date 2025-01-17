#### LIBRARIES ####
library(ggplot2)
library(data.table)
library(ggsignif)
library(tidyr)
library(umap)
library(pheatmap)
library(purrr)
library(RColorBrewer)
library(caret)
library("org.Hs.eg.db")
library(clusterProfiler)
library(AnnotationDbi)
library(stringr)
library(scales)
library(glmnet)
library(stringr)
library(reshape2)
library(ggrepel)
library(GO.db)
library(tibble)
library(VennDiagram)
library(eulerr)
library(enrichplot)
library(ComplexHeatmap)
library(circlize)
library(data.table)
library(dplyr)

source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/functions_te.R")

#### LOAD DATA #### 
te_raw_t <- fread("/Users/briannelaverty/Documents/R_Malkin/te/data/final/tumour_annotSV_full.tsv", sep="\t", header=TRUE)
te_split_t <- fread("/Users/briannelaverty/Documents/R_Malkin/te/data/final/tumour_annotSV_split.tsv", sep="\t", header=TRUE)
chr_length <- read.csv("/Users/briannelaverty/Documents/R_Malkin/te/data/raw/chromosome_length.csv", header=TRUE)
clinical <- read.delim("/Users/briannelaverty/Documents/R_Malkin/clinical/te_clinical.csv", sep=",", header=TRUE)
nohits_t_input <- fread("/Users/briannelaverty/Documents/R_Malkin/te/data/final/tumour_nohit_samples.txt", sep="\t", header=FALSE)
metrics <- read.delim("/Users/briannelaverty/Documents/R_Malkin/te/data/final/combined_metrics.txt", sep="\t", header=TRUE)
num_calls <- read.delim("/Users/briannelaverty/Documents/R_Malkin/te/data/final/count_calls.txt", sep="\t", header=TRUE)
nonproband <- read.delim("/Users/briannelaverty/Documents/R_Malkin/te/data/final/non_proband_normals", sep="\t", header=TRUE, colClasses = c("character"))  
noconsent <- read.delim("/Users/briannelaverty/Documents/R_Malkin/te/data/final/kics_samples_exclude", sep="\t", header=TRUE, colClasses = c("character"))
hg37_genes <- read.delim("/Users/briannelaverty/Documents/R_Malkin/te/data/final/hg37_genes.tsv", sep="\t", header=TRUE)
loh_time<- read.delim("/Users/briannelaverty/Documents/R_Malkin/te/data/final/loh_time", sep="\t", header=TRUE)
somatic_conversion <- read.delim("/Users/briannelaverty/Documents/R_Malkin/te/data/final/somatic_variant_conversion.txt", sep="\t", header=TRUE)
kics_tp53_somatic_variants <- fread("/Users/briannelaverty/Documents/R_Malkin/TE/data/kics_tp53_somatic_variants.csv", sep=",", header=TRUE)
plot_dir_tumour <- "/Users/briannelaverty/Library/Mobile Documents/com~apple~CloudDocs/Malkin_lab/figures/te/graphs/tumour/"

#### PREP ####
# metrics
metrics <- prep_metrics(metrics)

# clinical
clinical <- prep_clinical(clinical) 

# host seq with cancer
hostseq_cancer <- prep_hostseq(clinical)

# nohits
nohits_t <- replace_nohit_samples(nohits_t_input, clinical) # change sample names in nohits to match

# chr length
chr_length$chr <- factor(chr_length$chr, levels=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"X","Y")) # order factor
chr_lengths <- setNames(chr_length$length, chr_length$chr)

# gene size
gene_size <- add_gene_size(hg37_genes)

# te_raw
# rename alt column, replace sample with full sample id, remove samples based on non proband and no consent, filter by quality metrics
te_raw_prepped_t <- prep_te(te_raw_t, nonproband, noconsent, hostseq_cancer, metrics, c("mean_cov", "avg_quality", "pct_chimeras"), c(20, 30, 2), type="T") 


# remove outlier samples
# do this after count_TE_occurrences
samples_to_remove <- te_all_t$sample[te_all_t$total > 100000]
te_raw_prepped_t <- te_raw_prepped_t %>% 
  filter(!sample %in% samples_to_remove)
te_split_prepped_t <- te_split_prepped_t %>% 
  filter(!sample %in% samples_to_remove)

  
# filter out common TE, create count matrix, merge with clinical, filter for age, sort into 3 tables. only select one sample per patient
final_te_count_t <- process_te_data_tumour(te_raw_prepped_t, te_raw_prepped, clinical, nohits_t, metrics, apply_filter_common = TRUE, rare_gnomad = 3, rare_hostseq = 8, split_by_gene = FALSE, apply_process_combinations = TRUE, select_samples_split = TRUE)
load(paste0(plot_dir_tumour, "final_te_count_rare.RData"))
te_aff_t <- final_te_count_t$te_aff_selected
te_aff_all_t <- final_te_count_t$te_aff_all
te_lfs_t <- final_te_count_t$te_lfs_selected
te_lfs_all_t <- final_te_count_t$te_lfs_all
te_kics_t <- final_te_count_t$te_kics_selected
te_kics_all_t <- final_te_count_t$te_kics_all
te_all_t <- final_te_count_t$te_all_selected # all patients and selected samples 
te_all_all_t <- final_te_count_t$te_all_all # all samples and all patients 
save(final_te_count_t, file = paste0(plot_dir_tumour, "final_te_count_rare.RData"))

# filter by common, dont split by gene, dont create count. one row per full TE (not split by gene). doesn't add no hit samples, select one sample per patient
final_te_count_expand_t <- process_te_data_tumour(te_raw_prepped_t, te_raw_prepped, clinical, nohits_t, metrics, apply_filter_common = TRUE, rare_gnomad = 3, rare_hostseq = 8, split_by_gene = FALSE, apply_process_combinations = FALSE, select_samples_split = TRUE)
load(paste0(plot_dir_tumour, "final_te_count_expand_rare.RData"))
te_aff_expand_t <- final_te_count_expand_t$te_aff_selected
te_aff_expand_all_t <- final_te_count_expand_t$te_aff_all
te_lfs_expand_t <- final_te_count_expand_t$te_lfs_selected
te_lfs_expand_all_t <- final_te_count_expand_t$te_lfs_all
te_kics_expand_t <- final_te_count_expand_t$te_kics_selected
te_kics_expand_all_t <- final_te_count_expand_t$te_kics_all
te_all_expand_t <- final_te_count_expand_t$te_all_selected # all patients and selected samples 
te_all_expand_all_t <- final_te_count_expand_t$te_all_all # all samples and all patients 
save(final_te_count_expand_t, file = paste0(plot_dir_tumour, "final_te_count_expand_rare.RData"))

# filter by common, split by gene, no count. doesnt add no hit samples 
final_te_count_expand_split_t <- process_te_data_tumour(te_raw_prepped_t, te_raw_prepped, clinical, nohits_t, metrics, apply_filter_common = TRUE, rare_gnomad = 3, rare_hostseq = 8, split_by_gene = TRUE, apply_process_combinations = FALSE, select_samples_split = TRUE)
load(paste0(plot_dir_tumour, "final_te_count_split_rare.RData"))
te_aff_split_t <- final_te_count_expand_split_t$te_aff_selected
te_aff_split_all_t <- final_te_count_expand_split_t$te_aff_all
te_lfs_split_t <- final_te_count_expand_split_t$te_lfs_selected
te_lfs_split_all_t <- final_te_count_expand_split_t$te_lfs_all
te_kics_split_t <- final_te_count_expand_split_t$te_kics_selected
te_kics_split_all_t <- final_te_count_expand_split_t$te_kics_all
te_all_split_t <- final_te_count_expand_split_t$te_all_selected
te_all_split_all_t <- final_te_count_expand_split_t$te_all_all
save(final_te_count_expand_split_t, file = paste0(plot_dir_tumour, "final_te_count_split_rare.RData"))


# common
# keep common, dont split by gene, create count, select one sample per patient
final_te_count_common_t <- process_te_data_tumour(te_raw_prepped_t, te_raw_prepped, clinical, nohits_t, metrics, apply_filter_common = FALSE, split_by_gene = FALSE, apply_process_combinations = TRUE, select_samples_split = TRUE)
load(paste0(plot_dir_tumour, "final_te_count_common.RData"))
te_aff_common_t <- final_te_count_common_t$te_aff_selected
te_aff_all_common_t <- final_te_count_common_t$te_aff_all
te_lfs_common_t <- final_te_count_common_t$te_lfs_selected
te_lfs_all_common_t <- final_te_count_common_t$te_lfs_all
te_kics_common_t <- final_te_count_common_t$te_kics_selected
te_kics_all_common_t <- final_te_count_common_t$te_kics_all
te_all_common_t <- final_te_count_common_t$te_all_selected # all patients and selected samples 
te_all_all_common_t <- final_te_count_common_t$te_all_all # all samples and all patients 
save(final_te_count_common_t, file = paste0(plot_dir_tumour, "final_te_count_common.RData"))

# keep common, dont split by gene, dont create count. one row per full TE (not split by gene). doesn't add no hit samples
final_te_count_expand_split_common_t <- process_te_data_tumour(te_split_prepped_t, te_raw_prepped, clinical, nohits_t, metrics, apply_filter_common = FALSE, split_by_gene = FALSE, apply_process_combinations = FALSE, select_samples_split = TRUE)
#load(paste0(plot_dir_tumour, "final_te_count_expand_split_common.RData"))
te_aff_split_common_t <- final_te_count_expand_split_common_t$te_aff_selected
te_aff_split_all_common_t <- final_te_count_expand_split_common_t$te_aff_all
te_lfs_split_common_t <- final_te_count_expand_split_common_t$te_lfs_selected
te_aff_split_all_common_t <- final_te_count_expand_split_common_t$te_lfs_all
te_kics_split_common_t <- final_te_count_expand_split_common_t$te_kics_selected
te_aff_split_all_common_t <- final_te_count_expand_split_common_t$te_kics_all
save(final_te_count_expand_split_common_t, file = paste0(plot_dir_tumour, "final_te_count_expand_split_common.RData"))



## ADD LAST TWO FROM GERMLINE ###


# te_split
# rename alt column, replace sample with full sample id, remove samples based on non proband and no consent, filter by quality metrics
te_split_prepped_t <- prep_te(te_split_t, nonproband, noconsent, hostseq_cancer, metrics, c("mean_cov", "avg_quality", "pct_chimeras"), c(20, 30, 2), type="T") 

# dont use split directly because could have one te over 5 genes for one sample so counts inflated for common calculation
# only use to look at gene characteristics
# te input is split by gene, dont create count, filter by common
final_te_count_expand_split_t <- process_te_data_tumour(te_split_prepped_t, te_raw_prepped, clinical, nohits_t, metrics, apply_filter_common = TRUE, rare_gnomad = 3, rare_hostseq = 8, apply_process_combinations = FALSE, select_samples_split = TRUE)
te_aff_split_genes_t <- final_te_count_expand_split_t$te_aff_selected


#### GENERAL STATS ####
count_TE_occurrences(te_all_expand_t, nohits=nohits_t, te_all=te_all_t) # doesnt include no hit samples for summary 

# zero L1
total <- nrow(te_aff_t)
zero <- nrow(te_aff_t[te_aff_t$LINE1==0,])
one <- nrow(te_aff_t[te_aff_t$LINE1>0,])
one/total


# te calls
plot_te_counts_summary(te_all_t)
ggsave(paste0(plot_dir_tumour, "te_count_type.png"), width = 9, height = 5)

# tes shared by samples
plot_te_counts_unique(te_all_expand_t)
stacked_bar_plot_num_samples(te_all_expand_t, c(1, 3, 5)) # plot of # TEs in range of samples
ggsave(paste0(plot_dir_tumour, "te_type_sample.png"), width = 9, height = 5)

# percent of samples with at least 1 te 
x_tp53 <- expression("Germline " * italic("TP53") * " status")
compare_event_proportions(te_all_t, group_column="TP53_status", total_column="total", x_lab=x_tp53) 
ggsave(paste0(plot_dir_tumour, "te_tumour_odds_norepeats.png"), width = 9, height = 5)

# median number of TEs per sample
summary(te_kics_t$total)
summary(te_lfs_t$total)

# total te calls
plot_te_sum(te_kics_t)

te_aff_t[te_aff_t$total>10000, "sample"]
te_aff_t[te_aff_t$LINE1>100, "sample"]
te_aff_t[te_aff_t$LINE1>100, c("sample", "tumor_type", "sex", "age_at_diagnosis", "lesion_type", "disease_state", "sample_topography")]


#### TYPES OF COUNT ####
# wilcoxon test with all samples
plot_count_wilcox(te_aff_t, type=NA, chr=NA, group="TP53_status", x_lab="Disease state", y_lab="log(Total TE count + 1)", log_scale=TRUE)

# linear model controlling for multiple covariates. get p value for TP53 status
# group some tumour types into other
covar_mean <- c("mean_cov", "total_reads", "mean_read_len", "pct_chimeras", "avg_quality", "tumor_type")
covar_med <- c("med_cov", "total_reads", "med_read_len", "pct_chimeras", "avg_quality", "tumor_type")

# use this one
# model with TP53 status, quality metrics and tumour type. p value is from model  
model <- plot_count_lm(te_aff_t, min_samples=5, residuals=FALSE, log_scale=TRUE, group="TP53_status", covariates = covar_med, x_lab="Disease state", y_lab="Total TE count", type=NA, chr=NA)

# model with quality metrics and tumour type then wilcox with residuals. res are effects not caused by variables in model 
model <- plot_count_lm(te_aff_t, min_samples=5, residuals=TRUE, group="TP53_status", log_scale=FALSE, covariates = covar_med, x_lab="Disease state", y_lab="Total TE count", type=NA, chr=NA)

# bootstrap to test difference of test 
te_aff_t_bootstrap <- bootstrap_test(te_aff_t, "TP53_status", "total", mean, n_bootstraps = 10000, step_size=100) 


#### TE COUNT ####
types <- c(NA, "LINE1", "ALU", "SVA")
x_tp53 <- expression("Germline " * italic("TP53") * " status")

# affected
generate_plots(plot_count_wilcox, te_aff_t, chr=NA, type=types, group="TP53_status", x_lab="Disease state")
plot_count_wilcox(te_aff_t, type=NA, group="TP53_status", y_lab="LINE1 count", chr=NA, x_lab=x_tp53, log_scale=FALSE)
ggsave(paste0(plot_dir_tumour, "te_count_wilcox_aff_line.png"), width = 3, height = 5)

plot_count_lm(te_aff_common_t, type=NA, group="TP53_status", y_lab="Somatic repeat count", breaks=c(10, 100, 1000, 10000), covariates = covar_med, residuals=FALSE, log_scale=TRUE, x_lab=x_tp53, min_samples=5, chr=NA)
ggsave(paste0(plot_dir_tumour, "te_count_lm_aff_common_all.png"), width = 3, height = 5)
bootstrap_test(te_aff_t, "TP53_status", "total", mean, n_bootstraps = 10000, step_size=100) 

plot_count_lm(te_aff_common_t, type="LINE1", group="TP53_status", y_lab="Somatic LINE1 count", breaks=c(10, 100, 500), covariates = covar_med, residuals=FALSE, log_scale=TRUE, x_lab=x_tp53, min_samples=5, chr=NA)
ggsave(paste0(plot_dir_tumour, "te_count_lm_aff_common_line.png"), width = 3, height = 5)
plot_count_lm(te_aff_common_t, type="ALU", group="TP53_status", y_lab="Somatic ALU count", breaks=c(10, 100, 1000, 10000), covariates = covar_med, residuals=FALSE, log_scale=TRUE, x_lab=x_tp53, min_samples=5, chr=NA)
ggsave(paste0(plot_dir_tumour, "te_count_lm_aff_common_alu.png"), width = 3, height = 5)
plot_count_lm(te_aff_common_t, type="SVA", group="TP53_status", y_lab="Somatic SVA count", breaks=c(1, 10, 50), covariates = covar_med, residuals=FALSE, log_scale=TRUE, x_lab=x_tp53, min_samples=5, chr=NA)
ggsave(paste0(plot_dir_tumour, "te_count_lm_aff_common_sva.png"), width = 3, height = 5)

# affected primary initial
te_aff_primary_initial_t <- te_aff_t %>%
  subset(lesion_type=="primary") %>% 
  subset(disease_state=="initial") 

generate_plots(plot_count_wilcox, te_aff_primary_initial_t, chr=NA, type=types, group="TP53_status", x_lab="Disease state")
plot_count_wilcox(te_aff_primary_initial_t, type="SVA", group="TP53_status", y_lab="LINE1 count", chr=NA, x_lab=x_tp53, log_scale=TRUE)
ggsave(paste0(plot_dir_tumour, "te_count_wilcox_aff_line.png"), width = 3, height = 5)

plot_count_lm(te_aff_t, type="LINE1", group="TP53_status", y_lab="LINE1 count", covariates = covar_med, residuals=FALSE, log_scale=FALSE, x_lab=x_tp53, min_samples=5, chr=NA)
ggsave(paste0(plot_dir_tumour, "te_count_lm_aff_sva.png"), width = 3, height = 5)
bootstrap_test(te_aff_t, "TP53_status", "total", mean, n_bootstraps = 10000, step_size=100) 

#lfs
generate_plots(plot_count_wilcox, te_lfs_t, chr=NA, type=types, group="Cancer", x_lab="Cancer")
plot_count_wilcox(te_lfs_t, type=NA, chr=NA, group="Cancer", x_lab="Cancer", y_lab="Total TE count")

# by cluster
generate_plots(plot_count_kruskal, te_lfs_t, chr=NA, type=types, group="cluster", x_lab="Cluster")
plot_count_kruskal(te_lfs_t, type=NA, chr=NA, group="cluster", x_lab="Cluster", y_lab="Repeat count")
ggsave(paste0(plot_dir_tumour, "te_count_lfs_cluster.png"), width = 9, height = 5)

# by inheritance
te_lfs_t <- te_lfs_t[!grepl("DELETE", te_lfs_t$inheritance), ]
generate_plots(plot_count_kruskal, te_lfs_t, chr=NA, type=types, group="inheritance", x_lab="Inheritance")
plot_count_kruskal(te_lfs_t, type=NA, chr=NA, group="inheritance", x_lab="Inheritancee", y_lab="Total TE count")

# by variant location
generate_plots(plot_count_kruskal, te_lfs_t, chr=NA, type=types, group="Variant_location", x_lab="TP53 variant location")
plot_count_kruskal(te_lfs_t, type=NA, chr=NA, group="Variant_location", x_lab="TP53 variant location", y_lab="Total TE count")

# by variant classification 
generate_plots(plot_count_kruskal, te_lfs_t, chr=NA, type=types, group="Variant_Classification", x_lab="TP53 variant classification")
plot_count_kruskal(te_lfs_t, type=NA, chr=NA, group="Variant_Classification", x_lab="TP53 variant classification", y_lab="Total TE count")


#### COMPARE COUNT TO ADULT ####
te_gt1 <- nrow(te_aff_t[te_aff_t$total>0,])
te_0 <- nrow(te_aff_t[te_aff_t$total==0,])

fisher_test_and_plot(group1_yes = 1046, group1_no = 2954, group2_yes = 143, group2_no = 23)
ggsave(paste0(plot_dir_tumour, "te_tumour_atleast1.png"), width = 9, height = 5)

#### COUNT PER CHROMOSOME ####
# count per chromosome for one group
plot_count_perchr_notest(te_kics_t, type=NA, "Total TE count normalized by chromosome length", log_scale=FALSE)

# comparing by tp53 status
generate_plots_perchrom(plot_count_perchr, te_aff_t, type=types, group="TP53_status")
plot_count_perchr(te_aff_t, type=NA, "TP53_status", "Count")

# lfs  
generate_plots_perchrom(plot_count_perchr, te_lfs_t, type=types, group="Cancer")
plot_count_perchr(Cancer, type=NA, "Cancer", "Count")


#### COUNT BY TT ####
# count by tumour type
generate_plots(plot_count_kruskal_nogroup, te_aff_t, column="tumor_type", x_lab="Tumour type", type=types, chr=NA, min=3) 
plot_count_kruskal_nogroup(te_kics_t, column="tumor_type", type=NA, breaks=c(10, 100, 1000, 10000, 100000), min=3, x_lab="Tumour type", y_lab="Somatic repeat count", chr=NA, log_scale = TRUE) 
ggsave(paste0(plot_dir_tumour, "te_count_tt_kics_all.png"), width = 9, height = 5)

# count by tumour type and grouped by group
generate_plots(plot_count_tt, te_aff_t, chr=NA, type=types, group="TP53_status", min=2) 
plot_count_tt(te_aff_t, chr=NA, type=NA, group="TP53_status", min=3, y_lab="Repeat count", legend_title=x_tp53, log_scale=TRUE) 
ggsave(paste0(plot_dir_tumour, "te_count_tt_tp53.png"), width = 9, height = 5)


#### COUNT BY CLINCAL VARIABLES####
# count by age
generate_plots(plot_count_age, te_kics_t, chr=NA, type=types) 
plot_count_age(te_kics_t, chr=NA, type=NA, y_lab="Total TE count")

# treatment
generate_plots(plot_count_wilcox, te_aff_t, chr=NA, type=types, group="treatment", x_lab="Treatment status")
plot_count_wilcox(te_aff_t, type=NA, chr=NA, group="treatment", x_lab="Treatment status", y_lab="log(Total TE count + 1)", log_scale=TRUE)

plot_count_lm(te_aff_t, min_samples=3, residuals=FALSE, log_scale=TRUE, group="treatment", covariates = covar_med, x_lab="Disease state", y_lab="Total TE count", type=NA, chr=NA)
bootstrap_test(te_aff_t, "treatment", "total", mean, n_bootstraps = 10000, step_size=100) 

# lesion type
generate_plots(plot_count_kruskal, te_aff_t, chr=NA, type=types, group="lesion_type", x_lab="Lesion type", log_scale=TRUE) 
plot_count_kruskal(te_aff_t, type=NA, chr=NA, group="lesion_type", x_lab="Treatment status", y_lab="log(Total TE count + 1)", log_scale=TRUE)

plot_count_lm(te_aff_t, min_samples=5, residuals=FALSE, log_scale=TRUE, group="lesion_type", covariates = covar_med, x_lab="Disease state", y_lab="Total TE count", type=NA, chr=NA)
bootstrap_test(te_aff_t, "lesion_type", "total", mean, n_bootstraps = 10000, step_size=100) 

# disease state
generate_plots(plot_count_kruskal, te_aff_t, chr=NA, type=types, group="disease_state", x_lab="Disease state", log_scale=TRUE) 
plot_count_kruskal(te_aff_t, type=NA, chr=NA, group="disease_state", x_lab="Disease state", y_lab="TE count", log_scale=TRUE)

plot_count_lm(te_aff_t, min_samples=5, residuals=FALSE, log_scale=TRUE, group="disease_state", covariates = covar_med, x_lab="Disease state", y_lab="Total TE count", type=NA, chr=NA)
bootstrap_test(te_aff_t, "disease_state", "total", mean, n_bootstraps = 10000, step_size=100) 

# sex
plot_count_lm(te_aff_t, min_samples=5, residuals=FALSE, group="sex", log_scale=TRUE, covariates = covar_med, x_lab="Sex", y_lab="Total TE count", type=NA, chr=NA)
bootstrap_test(te_aff_t, "sex", "total", mean, n_bootstraps = 10000, step_size=100) 


#### FULL LENGTH ####
te_aff_expand_t_line <- te_aff_expand_t %>% filter(ALT=="LINE1")

summary(te_aff_expand_t_line$SV_length)

te_aff_expand_t_line_fulllength <- te_aff_expand_t_line %>% filter(SV_length >= 5900)
nrow(te_aff_expand_t_line_fulllength)

# process combinations
te_aff_expand_t_line_fulllength_processed <- process_all_combinations(te_aff_expand_t_line_fulllength)
te_aff_expand_t_line_fulllength_processed <- as.data.frame(add_nohit_samples(te_aff_expand_t_line_fulllength_processed, nohits_t))

# merge with clinical
te_aff_expand_t_line_fulllength_processed <- merge_dfs(te_aff_expand_t_line_fulllength_processed, clinical, include_all_x = FALSE)

# samples with lots of full length
te_aff_expand_t_line_fulllength_processed_lots <- te_aff_expand_t_line_fulllength_processed %>% filter(total>15)
table(te_aff_expand_t_line_fulllength_processed_lots$TP53_status, te_aff_expand_t_line_fulllength_processed_lots$tumor_type) 

# plot
plot_count_kruskal(te_aff_expand_t_line_fulllength_processed, type=NA, chr=NA, group="TP53_status", log_scale=TRUE, x_lab=x_tp53, y_lab="Full-length LINE1 count")
ggsave(paste0(plot_dir_tumour, "te_count_full_length.png"), width = 3, height = 5)

#### P53 MUTATION and FITNESS ####
# p53 mutation
# count by mutaiton
plot_count_kruskal_nogroup(te_lfs_t, column="mutation", min=3, chr=NA, type=NA, x_lab="p53 mutation", y_lab="Total TE count", log_scale = TRUE) 

# hotspot
hotspot <- c(175, 213, 245, 248, 273, 282, 337)

# plot te frequency with mutation along p53 gene
plot_te_locations(te_lfs_t, chr=NA, type=NA, hotspot=hotspot, log_scale=TRUE)
plot_te_locations(te_lfs_t, "ALU", hotspot)
plot_te_locations(te_lfs_t, "LINE1")


# fitness
tp53_fitness <- fread("/Users/briannelaverty/Documents/R_Malkin/TE/data/tp53_variants_predictions.tsv", sep="\t", header=TRUE)

# add fitness values
te_aff_t_tp53 <- prep_p53_fitness(te_aff_t, tp53_fitness)
te_lfs_t_tp53 <- prep_p53_fitness(te_lfs_t, tp53_fitness)

# scatter plot of count vs fitness
generate_scatterplots(scatter_template_lm_one_variable, te_aff_t_tp53, independent="p53_fitness", types=types, x_limit=c(0,1), y_limit=c(0,75), x_lab="p53 fitness")    
scatter_template_lm_one_variable(te_aff_t_tp53, "p53_fitness", "LINE1", c(0.25,1), c(0,1000), "p53 fitness", "Repeat count", shape_column = NULL, colour_column = NULL)
ggsave(paste0(plot_dir_tumour, "te_count_lfs_fitness.png"), width = 9, height = 5)

# inheritance  
plot_box_kruskal(te_lfs_t_tp53, "p53_fitness", "inheritance", "Inheritance", "p53 fitness")


#### CANCER GENES####
cpg<- read.csv("/Users/briannelaverty/Documents/R_Malkin/te/data/raw/kics_cpg.txt", header=FALSE, sep="\t") # cancer predisposition
genes <- cpg$V1

# remove outlier
samples_to_remove <- te_all_t$sample[te_all_t$total > 100000]
te_kics_split_t <- te_kics_split_t %>% 
  filter(!sample %in% samples_to_remove)
te_aff_split_t <- te_aff_split_t %>% 
  filter(!sample %in% samples_to_remove)


# overall cancer genes effected
# process combinations
te_aff_split_cancergenes_t  <- te_aff_split_t %>% filter(Gene_name %in% genes) # this is counting cancer genes effected not number of TEs that effect cancer genes as one TE may effect multiple cancer genes
te_aff_split_cancergenes_processed_t <- process_all_combinations(te_aff_split_cancergenes_t) 
te_aff_split_cancergenes_processed_t <- as.data.frame(add_nohit_samples(te_aff_split_cancergenes_processed_t, nohits_t))
summary(te_aff_split_cancergenes_processed_t$total)

# merge with clinical
te_aff_split_cancergenes_processed_t <- merge_dfs(te_aff_split_cancergenes_processed_t, clinical, include_all_x = FALSE)

# sig test
plot_count_wilcox(te_aff_split_cancergenes_processed_t, type="SVA", group="TP53_status", y_lab="Cancer genes effected", chr=NA, x_lab=x_tp53, log_scale=FALSE)


# top mutated cancer genes
# add gene size
te_kics_split_t <- add_gene_size_todf(te_kics_split_t, gene_size)
geneList_kics_t <- calculate_gene_mut_sample_frequency(te_kics_split_t, filter_exon = FALSE)
geneList_kics_filtered_t <- geneList_kics_t[geneList_kics_t$Gene_name %in% genes, ]
names(geneList_kics_t)

# plot top genes
plot_top_genes(geneList_kics_filtered_t, column="freq_samples_effected", label_column=NULL, x_lab="Proportion of samples with insertion in gene", top_n=20)
ggsave(paste0(plot_dir_tumour, "te_top_cancer_genes_affected_kics.png"), width = 9, height = 5)
plot_top_genes(geneList_kics_filtered_t, column="freq_samples_effected_normalized", label_column="num_samples_effected", x_lab="Proportion of samples with insertion in gene normalized for gene size", top_n=20)
ggsave(paste0(plot_dir_tumour, "te_top_cancer_genes_affected_normalized_kics.png"), width = 9, height = 5)

# patients with particular gene affected
te_kics_split_t %>%
  filter(Gene_name == "RAF1") %>%
  group_by(sample, tumor_type, TP53_status) %>%
  summarise(count = n(), .groups = 'drop') %>%
  as.data.frame()
te_aff_split_genes_t %>% 
  filter(Gene_name == "RAF1") %>%
  select(sample, Gene_name, Location, Location2, GnomAD_pLI, ExAC_pLI)


# genes present that meet critiera
location=c("CDS", "UTR", "3'UTR", "5'UTR", "5'UTR-CDS","CDS-3'UTR")

genes_aff_t <- identify_te_genes(te_aff_split_genes_t, gene_vector=genes, location="exon", location2=NULL)
genes_kics_t <- identify_te_genes(te_kics_split_t, gene_vector=genes, location=NULL, location2=NULL)
genes_lfs_t <- identify_te_genes(te_lfs_split_t, gene_vector=genes, location="exon", location2=NULL)

table(genes_aff_t$sample, genes_aff_t$Gene_name)
table(genes_aff_t$sample, genes_aff_t$TP53_status)
table(genes_aff_t$sample, genes_aff_t$tumor_type)
table(genes_aff_t$TP53_status, genes_aff_t$tumor_type)
table(genes_aff_t$Gene_name, genes_aff_t$tumor_type)
table(genes_aff_t$Gene_name, genes_aff_t$TP53_status)

length(table(genes_kics_t$Gene_name))
table(genes_lfs_t$Gene_name)


# wilcox test b/w lfs and controls
count_location_wilcox(te_aff_split_t, gene=genes, filter_element=NA, group="TP53_status", location=NULL)

# wilcox test b/w lfs aff and unaff
count_location_wilcox(te_lfs_split_t, gene=genes, filter_element=NA, group="Cancer", location=location)


#### TOP MUTATED GENES ####
# kics
te_kics_split_t <- add_gene_size_todf(te_kics_split_t, gene_size)
geneList_kics_t <- calculate_gene_mut_sample_frequency(te_kics_split_t)
names(geneList_kics_t)

# plot top genes
plot_top_genes(geneList_kics_t, column="freq_samples_effected", label_column=NULL, x_lab="Proportion of samples with insertion in gene", top_n=20)
ggsave(paste0(plot_dir_tumour, "te_top_genes_affected_kics.png"), width = 9, height = 5)
plot_top_genes(geneList_kics_t, column="freq_samples_effected_normalized", label_column="num_samples_effected", x_lab="Proportion of samples with insertion in gene normalized by gene size", top_n=20)
ggsave(paste0(plot_dir_tumour, "te_top_genes_affected_normalized_kics.png"), width = 9, height = 5)

# number of genes affecting > 50% of samples
count_genes_by_threshold(geneList_kics_t, column_name = "freq_samples_effected", thresholds = c(0.5, 0.75, 0.9))


# lfs
te_lfs_split_t <- add_gene_size_todf(te_lfs_split_t, gene_size)
geneList_lfs_t<- calculate_gene_mut_sample_frequency(te_lfs_split_t, nsample_thresh = 0, filter_exon = FALSE)

# plot top genes
plot_top_genes(geneList_lfs_t, column="freq_samples_effected", label_column=NULL, x_lab="Proportion of samples with insertion in gene normalized by gene size", top_n=20)
plot_top_genes(geneList_lfs_t, column="freq_samples_effected_normalized", label_column="freq_samples_effected", x_lab="Proportion of samples with insertion in gene normalized by gene size", top_n=20)
ggsave(paste0(plot_dir_tumour, "te_top_genes_affected_normalized_lfs.png"), width = 9, height = 5)


# control
te_control_split_t <- te_aff_split_t %>% filter(TP53_status == "Control")
geneList_control_t<- calculate_gene_mut_sample_frequency(te_control_split_t, nsample_thresh = 0, filter_exon = FALSE)
plot_top_genes(geneList_control_t, 20)


# specific gene heatmap of insertions
plot_te_insertion_complex_heatmap(te_kics_split_t, gene_name="RAF1", bin_size=1000, annotations=c("tumor_class", "age_at_diagnosis", "TP53_status"), lwd=0)
#ggsave(paste0(plot_dir_tumour, "te_location_kics_gpc5.png"), width = 9, height = 5)


#### PATHWAY PEDIATRIC ####
# over representation analysis
ora_kics_t <- perform_ora(te_kics_split_t, nsample_thresh = 0, filter_exon = FALSE)

# viz 
barplot(ora_kics_t, showCategory=20) 
dotplot(ora_kics_t, showCategory=20) 
cnetplot(ora_kics_t, showCategory = 10, colorEdge = TRUE, node_label = "category") 
ggsave(paste0(plot_dir_tumour, "te_pathway_ora_kics.png"), width = 14, height = 9)
#ora_kics_pairwise_t<- pairwise_termsim(ora_kics_t)
#emapplot(ora_kics_pairwise_t, group_category = TRUE, group_legend = TRUE) +
#  scale_fill_identity(values = "#0080A3") +
#  guides(fill = "none")

# group descriptions by broader function
ora_kics_simple_t <- simplify(ora_kics_t, cutoff=0.6, by="p.adjust", select_fun=min)
barplot(ora_kics_simple_t, showCategory=20) 
cnetplot(ora_kics_simple_t, showCategory = 10, colorEdge = TRUE, node_label = "category")
ggsave(paste0(plot_dir_tumour, "te_pathway_kics_simplified.png"), width = 14, height = 9)


# common genes in pathway
# freq of genes in top pathways
pathway_genes_t <- plot_pathway_gene_counts(ora_kics_t, n_descriptions=nrow(ora_kics_t), n_genes=50)

# look at correlation b/w top genes mutated and genes in top pathways
# column is num_samples_effected or mutation_count
analyze_gene_mutations(geneList_kics_t, ora_kics_t, n_descriptions=nrow(ora_kics_t), column = "num_samples_effected", y_lab="Number of samples with insertion in gene")
ggsave(paste0(plot_dir_tumour, "te_pathwaygenes_vs_numsamples_kics.png"), width = 9, height = 5)


# gene set enrichment analysis
# add giene size
te_kics_split <- replace_gene_names(te_kics_split)
te_kics_split <- add_gene_size_todf(te_kics_split, gene_size)

# ensure max genes still 5000
gsea_kics <- perform_gsea(te_kics_split, column="freq_samples_effected", nsample_thresh = 0, filter_exon = FALSE)

# viz 
dotplot(gsea_kics, showCategory=20) 
goplot(gsea_kics, showCategory = 5)
cnetplot(gsea_kics, showCategory = 10, colorEdge = TRUE, node_label = "category")
cnetplot(gsea_kics, showCategory = 5, node_label = "category", circular = TRUE, colorEdge = TRUE)   


#### CILIA GENES ####
te_kics_split_t <- add_gene_size_todf(te_kics_split_t, gene_size)
geneList_kics_t <- calculate_gene_mut_sample_frequency(te_kics_split_t)

cilia_genes <- ora_kics_t@result[ora_kics_t@result$Description %in% c("cilium assembly", "cilium organization"), "geneID"]
cilia_list <- filter_cilia_genes(cilia_genes, geneList_kics_t)
head(cilia_list)

generate_te_insertion_heatmap_genes_cut(te_kics_split_t, te_aff_variant_t, cilia_genes, top_percent=100, at_least_one_insertion = FALSE)
generate_te_insertion_heatmap_genes_cut_quantiles(te_kics_split_t, te_aff_t, cilia_genes, top_percent=NULL, at_least_one_insertion = FALSE)
generate_te_insertion_heatmap_genes_cut_quantiles_normalized(te_kics_split_t, te_aff_variant_t, cilia_genes, top_percent=NULL, at_least_one_insertion = TRUE)
plot_te_gene_insertions(te_kics_split_t, te_aff_t, cilia_genes, xlim=NULL, ylim=NULL)
plot_te_gene_insertions(te_kics_split_t, te_aff_t, cilia_genes, xlim=c(0,75), ylim=c(0,6000))

neuron_genes<- ora_kics_t@result[ora_kics_t@result$Description %in% c("regulation of small GTPase mediated signal transduction"), "geneID"]
neuron_list <- filter_cilia_genes(neuron_genes, geneList_kics_t)
generate_te_insertion_heatmap_genes_cut_quantiles_normalized(te_kics_split_t, te_aff_variant_t, neuron_genes, top_percent=NULL, at_least_one_insertion = TRUE)
plot_te_gene_insertions(te_kics_split_t, te_aff_t, neuron_genes, xlim=NULL, ylim=NULL)

calculate_shared_diff(cilia_list, neuron_list, "Gene_name")

# samples with specific gene mutated
s <- unique(te_kics_split_t[te_kics_split_t$Gene_name=="SPAG16", "sample"])
te_aff_t[te_aff_t$sample %in% s, c("sample", "total", "TP53_status", "tumor_type")]

#### PATHWAY TP53 ####
ora_tp53_t <- perform_ora_tp53(te_aff_split_t)
ora_tp53_genes_t <- perform_ora_tp53(te_aff_split_genes_t)

# viz 
dotplot(ora_tp53_genes_t, showCategory = 20)
cnetplot(ora_tp53_t, showCategory = 10, colorEdge = TRUE, node_label = "category")
ora_tp53_pairwise_t<- pairwise_termsim(ora_tp53_genes_t)
emapplot(ora_tp53_pairwise_t, pie="count", showCategory = 20, group_category=TRUE, group_legend=TRUE) + scale_fill_manual(values=colours)
ggsave(paste0(plot_dir_tumour, "te_pathway_tp53.png"), width = 14, height = 9)

# simplify
ora_tp53_simple_t <- simplify(ora_tp53_t, cutoff=0.5, by="p.adjust", select_fun=min)
dotplot(ora_tp53_simple_t, showCategory=20) 
cnetplot(ora_tp53_simple_t, showCategory = 10, colorEdge = TRUE, node_label = "category")


# common genes
# enrichment result
ora_tp53_result_t <- ora_tp53_t@compareClusterResult

# GO terms that were only enriched in KiCS or LFS not both
ora_tp53_result_unique_t <- ora_tp53_result_t %>%
  group_by(ID) %>%
  filter(n() == 1) %>%
  ungroup()

# split enrichment by group
ora_tp53_result_control_t <- ora_tp53_result_unique_t %>% filter(Cluster=="Control")
ora_tp53_result_lfs_t <- ora_tp53_result_unique_t %>% filter(Cluster=="LFS")

# freq of genes in top pathways
pathway_genes_control_t <- plot_pathway_gene_counts(ora_tp53_result_control_t, n_descriptions=nrow(ora_tp53_result_control_t), n_genes=50)
pathway_genes_lfs_t<- plot_pathway_gene_counts(ora_tp53_result_lfs_t, n_descriptions=nrow(ora_tp53_result_lfs_t), n_genes=50)

# look at correlation b/w top genes mutated and genes in top pathways
# column is num_samples_effected or mutation_count
analyze_gene_mutations(geneList_control_t, pathway_genes_control_t, column = "freq_samples_effected") # this is different than the kics one because it is only pathways unique to controls not all sig enriched pathways
analyze_gene_mutations(geneList_lfs_t, pathway_genes_lfs_t, column = "freq_samples_effected")


#### UNIQUE TE IN LFS ####
specific_te_fisher <- as.data.frame(fisher_test_unique_te(te_aff_expand_t, min=3))
head(specific_te_fisher)
sig_te <- specific_te_fisher %>% filter(fisher_p_value_BH < 0.1)

# samples with significant TE
sig_te_samples <- find_sig_te_samples(sig_te, te_aff_expand)

# sample characteristics
for (i in 1:length(sig_te_samples)) {
  df <- sig_te_samples[[i]]
  print_summary_sig_te_samples(df)
}
one <- sig_te_samples[[1]]

# genes effected
# two TE dont effect gene
# one te in intron of NUMB, no frameshift, and gnomad and exac pLI close to 0 meaning tolerant of LOF variants 
sig_te_genes<- find_sig_te_genes(sig_te, te_aff_split_genes)
te_aff_split_genes


#### GERMLINE VARIANTS ####
# load variants, add variants to te_aff_t, make overall TP53 column
kics_germline_variants<- fread("/Users/briannelaverty/Documents/R_Malkin/TE/data/kics_germline_variants_panel.csv", sep = ",", header = TRUE)

# Pad KiCS_ID with leading zeros to ensure they are four digits
kics_germline_variants$sample<- sprintf("%04d", as.numeric(kics_germline_variants$KiCS_ID))

# Filter for KiCS_ID in sample_data$sample
kics_germline_variants <- subset(kics_germline_variants, sample %in% te_kics_t$base_sample)
kics_germline_variants <- subset(kics_germline_variants, interpretation %in% c("Pathogenic", "Likely Pathogenic"))
  
# Rename geneSymbol to gene if it exists
kics_germline_variants <- kics_germline_variants %>% rename("gene" = "geneSymbol")
  
# Number of samples with P/LP variants by gene
table(kics_germline_variants$gene)

# Add variant columns
te_variant_df <- add_variant_columns(data = te_kics_t, variant_data = kics_germline_variants, variant_type = "germline")

filter_and_print_g_columns(te_variant_df) 
# 2/4 cancers with >10000 TE have gemrline MPO variant
# Meningioma and NBL

high <- subset(te_variant_df, sample %in% te_aff_t[te_aff_t$LINE1>100, "sample"])
high[,c(1,144:162)]



te_aff_variant_t <- prep_variant_df(te_kics_t, te_aff_t)

generate_plots(plot_count_kruskal, te_aff_variant_t, chr=NA, type=types, group="TP53", x_lab="Cluster")
plot_count_kruskal(te_aff_variant_t, group="TP53", x_order=c("None", "Somatic", "Germline"), x_lab="TP53 variant", y_lab="Total TE count", type=NA, chr=NA, log_scale=TRUE)


#### TP53 VARIANTS COUNT ####
kics_tp53_somatic_variants <- fread("/Users/briannelaverty/Documents/R_Malkin/TE/data/kics_tp53_somatic_variants.csv", sep=",", header=TRUE)
somatic_conversion <- read.delim("/Users/briannelaverty/Documents/R_Malkin/te/data/final/somatic_variant_conversion.txt", sep="\t", header=TRUE)

# samples in kics 
kics_tp53_somatic_variants_processed <- prep_variants_tumour(kics_tp53_somatic_variants, somatic_conversion, te_kics_t)

# add variant columns 
te_aff_variant_t <- add_variant_columns(data=te_aff_t, variant_data=kics_tp53_somatic_variants_processed, variant_type="somatic")

# check got all 
out <- te_aff_variant_t[te_aff_variant_t$s_TP53==1, "sample"]
var <- unique(kics_tp53_somatic_variants_processed$sample)
var[!var%in% out] # samples with variant and base sample match but dont have te calls for that specific sample
te_aff_t[te_aff_t$base_sample=="0010", "sample"] # see specific sample that we have

# change TP53 status column
te_aff_variant_t <- te_aff_variant_t %>%
  mutate(TP53 = case_when(
    TP53_status == "Mutant" ~ "Germline", 
#    g_TP53 == 1 ~ "Germline",
    s_TP53 == 1 ~ "Somatic",
    TRUE ~ "None"  # Catch-all case for "else is None"
  ))

generate_plots(plot_count_kruskal, te_aff_variant_t, chr=NA, type=types, group="TP53", x_lab="Cluster")
plot_tp53_variants(te_aff_variant_t, y_lab="LINE1 count", type="SVA", group="TP53", breaks=c(10, 100, 1000, 10000), x_order=c("None", "Somatic", "Germline"), x_lab="TP53 variant", chr=NA, log_scale=TRUE)
ggsave(paste0(plot_dir_tumour, "te_tumour_tp53_byvariant_line1.png"), width = 9, height = 5)


#### MERGED TUMOURS ####
# load lfs clonality
clonality <- read.csv("/Users/briannelaverty/Documents/R_Malkin/clinical/lfs_clonality.csv", header=TRUE)

# nicks multi region samples. scatter with patient on x and te count on y.
# colour based on tumour type
plot_multisample_scatter_nick(te_all_all_t, te_count_col ="total")
# colour based on clonality
plot_multisample_scatter_clonality(te_all_all_t, clonal_df=clonality, clonal_column="ssm_prop_clonal", legend_lab = "SSM clonality", te_count_col ="total")
ggsave(paste0(plot_dir, "te_tumour_lfs_clonality.png"), width = 9, height = 5)


#### MULTISAMPLE ####

# kics multi sample. patinet on x and count on y. colour and shape differentiates disease adn lesion 
plot_multisample_scatter_kics(te_all_all_t, te_count_col ="total")

# plot for each patient with multiple samples
plot_individual_patient_samples(te_all_all_t)
plot_individual_patient_samples_facet(te_all_all_t)
plot_specific_sample(te_all_all_t, "0074", x_nudge=-100, y_nudge=1, log_scale=FALSE) 

# scatter with patient on x for SAMPLES TAKEN AT SAME TIMEPOINT 
plot_multisample_sametime(te_all_all_t)

# percent change by sample
plot_te_change(te_all_all, log_scale=TRUE)
ggsave(paste0(plot_dir, "te_tumour_kics_change.png"), width = 9, height = 5)


#### TE SOURCE ####
te_aff_expand_line_t <- te_aff_expand_t %>% filter(ALT=="LINE1")
te_aff_expand_line_t <- extract_info_fields(te_aff_expand_line_t)

# number of unqiue sources
nrow(table(te_aff_expand_line_t$source)) 
nrow(te_aff_expand_line_t)
423/2790
22/2790 # not siblings

# specific sources
table(te_aff_expand_line_t$source)[table(te_aff_expand_line_t$source) > 1]

# samples with sources
table(te_aff_expand_line_t[te_aff_expand_line_t$source != "not_transduction", c("sample", "source")])
sources <- te_aff_expand_line_t %>%
  filter(source != "not_transduction") %>% # Exclude unwanted rows
  count(sample, source) %>% # Group by sample and source, count occurrences
  filter(n > 0) %>% # Keep only rows where count is greater than 0
  arrange(desc(n))






# te aff and te lfs
te_aff_expand_t <- extract_info_fields(te_aff_expand_t)
te_lfs_expand_t <- extract_info_fields(te_lfs_expand_t)

nrow(table(te_aff_expand_t$source)) # number of unique sources
nrow(te_aff_expand_t)

# sepcfiic source
table(te_aff_expand_t$source)[table(te_aff_expand_t$source) > 2]
table(te_lfs_expand_t$source)[table(te_lfs_expand_t$source) > 1]

# samples with sources
table(te_aff_expand_t[te_aff_expand_t$source != "not_transduction", c("sample", "source")])
sources <- te_aff_expand_t %>%
  filter(source != "not_transduction") %>% # Exclude unwanted rows
  count(sample, source) %>% # Group by sample and source, count occurrences
  filter(n > 0) %>% # Keep only rows where count is greater than 0
  arrange(desc(n))

#### FULL LENGTH ####
te_aff_expand_line_t <- te_aff_expand_t %>% filter(ALT=="LINE1")

summary(te_aff_expand_line_t$SV_length)

te_aff_expand_line_fulllength_t <- te_aff_expand_line_t %>% filter(SV_length >= 5900)
nrow(te_aff_expand_line_fulllength_t)

# process combinations
te_aff_expand_line_fulllength_processed_t <- process_all_combinations(te_aff_expand_line_fulllength_t)
te_aff_expand_line_fulllength_processed_t <- as.data.frame(add_nohit_samples(te_aff_expand_line_fulllength_processed_t, nohits_t))

# merge with clinical
te_aff_expand_line_fulllength_processed_t <- merge_dfs(te_aff_expand_line_fulllength_processed_t, clinical, include_all_x = FALSE)

# plot
plot_count_wilcox(te_aff_expand_line_fulllength_processed_t, type="LINE1", group="TP53_status", y_lab="LINE1 count", chr=NA, x_lab=x_tp53, log_scale=FALSE)
ggsave(paste0(plot_dir, "te_count_wilcox_aff_line.png"), width = 3, height = 5)


#### TP53 LOH TIMING ####
loh_time <- loh_time %>% 
  select(sample, time) %>%
  mutate(sample = ifelse(grepl("^KiCS", sample), 
                         paste0(sub("^KiCS", "", sample), "_T"), 
                         sample))

timing <- merge(te_aff_t, loh_time, by="sample")

# checking samples 
loh_time$sample
timing$sample
te_aff_t[te_aff_t$base_sample=="0366", "sample"]
match_samples(loh_time, te_aff_t)

# wilcox for time 
wilcox.test(time ~ TP53_status, data = timing)

# save
write.table(timing, file = "/Users/briannelaverty/Documents/R_Malkin/te/data/final/te_loh.csv", quote = FALSE, sep = ',', row.names = FALSE, col.names = TRUE) 


#### DESCRIBE DATASET ####
# extract clinical from samples with have
clin<- te_all_t %>% dplyr::select(sample, tumor_type, sex, age_at_diagnosis, TP53_status)
clin <- clin %>% distinct()

# Plot pie chart for tumor_type
plot_pie_chart(clin, "tumor_type", "Tumour type")
ggsave(paste0(plot_dir_tumour, "clinical_germline_tt.png"), width = 6, height = 5)

# Plot pie chart for TP53_status
plot_pie_chart(clin, "TP53_status", x_tp53)
ggsave(paste0(plot_dir_tumour, "clinical_germline_tp53.png"), width = 6, height = 5)

# Plot pie chart for sex 
plot_pie_chart(clin, "sex", "Sex")
ggsave(paste0(plot_dir_tumour, "clinical_sex.png"), width = 6, height = 5)

# Plot histogram for age_at_diagnosis
plot_histogram_age(clin, bin=2)
ggsave(paste0(plot_dir_tumour, "clinical_germline_age.png"), width = 9, height = 5)



### end ####

