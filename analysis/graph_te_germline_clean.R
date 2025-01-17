plot_dir <- "/Users/briannelaverty/Library/Mobile Documents/com~apple~CloudDocs/Malkin_lab/figures/te/graphs/germline/"

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
library(reshape2)
library(ggrepel)
library(GO.db)
library(tibble)
library(VennDiagram)
library(eulerr)
library(enrichplot)
library(ComplexHeatmap)
library(circlize)
library(dplyr)

source("/Users/briannelaverty/Documents/R_Malkin/te/scripts/functions_te.R")

#### LOAD DATA #### 
te_raw <- fread("/Users/briannelaverty/Documents/R_Malkin/te/data/final/germline_annotSV_full.tsv", sep="\t", header=TRUE)
te_split<- fread("/Users/briannelaverty/Documents/R_Malkin/te/data/final/germline_annotSV_split.tsv", sep="\t", header=TRUE)
chr_length <- read.csv("/Users/briannelaverty/Documents/R_Malkin/te/data/raw/chromosome_length.csv", header=TRUE)
clinical <- read.delim("/Users/briannelaverty/Documents/R_Malkin/clinical/te_clinical.csv", sep=",", header=TRUE)
nohits <- fread("/Users/briannelaverty/Documents/R_Malkin/te/data/final/germline_nohit_samples.txt", sep="\t", header=FALSE)
metrics <- read.delim("/Users/briannelaverty/Documents/R_Malkin/te/data/final/combined_metrics.txt", sep="\t", header=TRUE)
num_calls <- read.delim("/Users/briannelaverty/Documents/R_Malkin/te/data/final/count_calls.txt", sep="\t", header=TRUE)
nonproband <- read.delim("/Users/briannelaverty/Documents/R_Malkin/te/data/final/non_proband_normals", sep="\t", header=TRUE, colClasses = c("character"))  
noconsent <- read.delim("/Users/briannelaverty/Documents/R_Malkin/te/data/final/kics_samples_exclude", sep="\t", header=TRUE, colClasses = c("character"))
hg37_genes <- read.delim("/Users/briannelaverty/Documents/R_Malkin/te/data/final/hg37_genes.tsv", sep="\t", header=TRUE)
ancestry <- read.delim("/Users/briannelaverty/Documents/R_Malkin/te/data/final/kics_predicted_ancestry.txt", sep="\t", header=TRUE)

#### PREP ####
# metrics
metrics <- prep_metrics(metrics)

# clinical
clinical <- prep_clinical(clinical) 

# host seq with cancer
hostseq_cancer <- prep_hostseq(clinical)

# nohits
nohits <- replace_nohit_samples(nohits, clinical) # change sample names in nohits to match

# chr length
chr_length$chr <- factor(chr_length$chr, levels=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"X","Y")) # order factor
chr_lengths <- setNames(chr_length$length, chr_length$chr)

# gene size
gene_size <- add_gene_size(hg37_genes)

# ancestry
ancestry <- prep_ancestry_kics(ancestry)

# te_raw
# rename alt column, replace sample with full sample id, remove samples based on non proband and no consent, filter by quality metrics
te_raw_prepped <- prep_te(te_raw, nonproband, noconsent, hostseq_cancer, metrics, c("mean_cov", "avg_quality", "pct_chimeras"), c(20, 30, 2), type="N") 
save(te_raw_prepped, file = paste0(plot_dir, "te_raw_prepped.RData"))

# filter out common TE, create count matrix,  merge with clinical, filter for age, prep formatting, sort into 3 tables
final_te_count <- process_te_data_germline(te_raw_prepped, clinical, nohits, metrics, ancestry, apply_filter_common = TRUE, rare_gnomad = 3, rare_hostseq = 8, split_by_gene = FALSE, apply_process_combinations = TRUE)
load(paste0(plot_dir, "final_te_count_rare.RData"))
te_aff <- final_te_count$te_aff
te_lfs <- final_te_count$te_lfs
te_kics <- final_te_count$te_kics
te_all<- final_te_count$te_all
save(final_te_count, file = paste0(plot_dir, "final_te_count_rare.RData"))

# filter by common, dont split by gene, dont create count. one row per full TE (not split by gene). doesn't add no hit samples
final_te_count_expand <- process_te_data_germline(te_raw_prepped, clinical, nohits, metrics, ancestry, apply_filter_common = TRUE, rare_gnomad = 3, rare_hostseq = 8, split_by_gene = FALSE, apply_process_combinations = FALSE)
load(paste0(plot_dir, "final_te_count_expand_rare.RData"))
te_aff_expand <- final_te_count_expand$te_aff
te_lfs_expand <- final_te_count_expand$te_lfs
te_kics_expand <- final_te_count_expand$te_kics
te_all_expand <- final_te_count_expand$te_all
save(final_te_count_expand, file = paste0(plot_dir, "final_te_count_expand_rare.RData"))

# filter by common, split by gene, no count. doesnt add no hit samples 
final_te_count_split <- process_te_data_germline(te_raw_prepped, clinical, nohits, metrics, ancestry, apply_filter_common = TRUE, rare_gnomad = 3, rare_hostseq = 8, split_by_gene = TRUE, apply_process_combinations = FALSE)
load(paste0(plot_dir, "final_te_count_split_rare.RData"))
te_aff_split <- final_te_count_split$te_aff
te_lfs_split <- final_te_count_split$te_lfs
te_kics_split <- final_te_count_split$te_kics
te_all_split <- final_te_count_split$te_all
save(final_te_count_split, file = paste0(plot_dir, "final_te_count_split_rare.RData"))


# common
# keep common, dont split by gene, create count
final_te_count_common <- process_te_data_germline(te_raw_prepped, clinical, nohits, metrics, ancestry, apply_filter_common = FALSE, split_by_gene = FALSE, apply_process_combinations = TRUE)
#load(paste0(plot_dir, "final_te_count_common.RData"))
te_aff_common<- final_te_count_common$te_aff
te_lfs_common<- final_te_count_common$te_lfs
te_kics_common<- final_te_count_common$te_kics
te_all_common<- final_te_count_common$te_all
save(final_te_count_common, file = paste0(plot_dir, "final_te_count_common.RData"))

# keep common, dont split by gene, dont create count. one row per full TE (not split by gene). doesn't add no hit samples
final_te_count_expand_common <- process_te_data_germline(te_raw_prepped, clinical, nohits, metrics, ancestry, apply_filter_common = FALSE, split_by_gene = FALSE, apply_process_combinations = FALSE)
load(paste0(plot_dir, "final_te_count_expand_common.RData"))
te_aff_expand_common <- final_te_count_expand_common$te_aff
te_lfs_expand_common <- final_te_count_expand_common$te_lfs
te_kics_expand_common <- final_te_count_expand_common$te_kics
te_all_expand_common <- final_te_count_expand_common$te_all
save(final_te_count_expand_common, file = paste0(plot_dir, "final_te_count_expand_common.RData"))

# keep common, split by gene, no count. doesnt add no hit samples 
final_te_count_split_common <- process_te_data_germline(te_raw_prepped, clinical, nohits, metrics, ancestry, apply_filter_common = FALSE, split_by_gene = TRUE, apply_process_combinations = FALSE)
#load(paste0(plot_dir, "final_te_count_split_common.RData"))
te_aff_split_common<- final_te_count_split_common$te_aff
te_lfs_split_common<- final_te_count_split_common$te_lfs
te_kics_split_common<- final_te_count_split_common$te_kics
te_all_split_common<- final_te_count_split_common$te_all
save(final_te_count_split_common, file = paste0(plot_dir, "final_te_count_split_common.RData"))


# te_split
# rename alt column, replace sample with full sample id, remove samples based on non proband and no consent, filter by quality metrics
te_split_prepped <- prep_te(te_split, nonproband, noconsent, hostseq_cancer, metrics, c("mean_cov", "avg_quality", "pct_chimeras"), c(20, 30, 2), type="N") 

# dont use split directly because could have one te over 5 genes for one sample so counts inflated for common calculation
# only use to look at gene characteristics
# te input is split by gene, dont create count, filter by common
final_te_count_expand_split <- process_te_data_germline(te_split_prepped, clinical, nohits, metrics, ancestry, apply_filter_common = TRUE, rare_gnomad = 3, rare_hostseq = 8, apply_process_combinations = FALSE)
load(paste0(plot_dir, "final_te_count_expand_split_nocommon.RData"))
te_aff_split_genes <- final_te_count_expand_split$te_aff
te_lfs_split_genes <- final_te_count_expand_split$te_lfs
te_kics_split_genes <- final_te_count_expand_split$te_kics
te_all_split_genes <- final_te_count_expand_split$te_all
save(final_te_count_expand_split, file = paste0(plot_dir, "final_te_count_expand_split.RData"))


#### GENERAL STATS ####
count_TE_occurrences(te_all_expand, nohits=nohits, te_all=te_all) # doesnt include no hit samples for summary 

# te calls all common
plot_te_counts_summary(te_all_common, log_scale=TRUE, breaks=c(10,100,1000))
ggsave(paste0(plot_dir, "te_count_type_all_common.png"), width = 9, height = 5)

# te calls all rare
plot_te_counts_summary(te_all, log_scale=TRUE, breaks=c(10,100,500))
ggsave(paste0(plot_dir, "te_count_type_all_rare.png"), width = 9, height = 5)

# te calls aff rare
plot_te_counts_summary(te_aff, log_scale=TRUE, breaks=c(10,100,500))
ggsave(paste0(plot_dir, "te_count_type_aff_rare.png"), width = 9, height = 5)

# tes shared by samples
plot_te_counts_unique(te_aff_expand)
stacked_bar_plot_num_samples(te_aff_expand_common, c(1, 5, 20)) # plot of # TEs in range of samples
ggsave(paste0(plot_dir, "te_type_sample_common.png"), width = 9, height = 5)

# venn diagram of callers
plot_venn(num_calls)
#ggsave(paste0(plot_dir, "te_callers.png"), width = 9, height = 5) # doesnt work bc grob

# percent of samples with at least 1 te 
nrow(te_kics[te_kics$total > 0,])/nrow(te_kics) * 100
nrow(te_lfs[te_lfs$total > 0,])/nrow(te_lfs) * 100

# median number of TEs per sample
summary(te_kics$total)
summary(te_lfs$total)

# total te calls
plot_te_sum(te_kics)


#### TYPES OF COUNT ####
# wilcoxon test with all samples
plot_count_wilcox(te_aff, type=NA, chr=NA, group="TP53_status", x_lab="Disease state", y_lab="log(Total TE count + 1)", log_scale=TRUE)

# linear model controlling for multiple covariates. get p value for TP53 status
# group some tumour types into other
covar_mean <- c("mean_cov", "total_reads", "mean_read_len", "pct_chimeras", "avg_quality", "tumor_type")#, "predicted_ancestry_thres")
covar_med <- c("med_cov", "total_reads", "med_read_len", "pct_chimeras", "avg_quality", "tumor_type")#, "predicted_ancestry_thres")

# use this one
# model with TP53 status, quality metrics and tumour type. p value is from model  
model <- plot_count_lm(te_aff, min_samples=5, residuals=FALSE, log_scale=TRUE, group="TP53_status", covariates = covar_med, x_lab="Disease state", y_lab="Total TE count", type=NA, chr=NA)

# model with quality metrics and tumour type then wilcox with residuals. res are effects not caused by variables in model 
model <- plot_count_lm(te_aff, min_samples=5, residuals=TRUE, group="TP53_status", log_scale=FALSE, covariates = covar_med, x_lab="Disease state", y_lab="Total TE count", type=NA, chr=NA)

# bootstrap to test difference of test 
te_aff_bootstrap <- bootstrap_test(te_aff, "TP53_status", "total", mean, n_bootstraps = 10000, step_size=100) 


#### TE COUNT ####
types <- c(NA, "LINE1", "ALU", "SVA")
x_tp53 <- expression("Cancer status and germline " * italic("TP53") * " genotype")

# all
# df with kics, lfs with cancer, and hostseq
te_cancer_hostseq <- te_all %>%
  filter(cohort == "HostSeq") %>%  # Filter rows where cohort is HostSeq
  bind_rows(te_aff)  # Combine with te_aff dataframe

generate_plots(plot_count_kruskal, te_cancer_hostseq, type=types, log_scale=TRUE, group="cancer_cohort", x_lab="Cluster", chr=NA)
plot_count_kruskal(te_cancer_hostseq, type=NA, y_lab="Rare repeat count", chr=NA, group="cancer_cohort", log_scale=TRUE, x_lab=x_tp53)
ggsave(paste0(plot_dir, "te_count_cancer_cluster_all.png"), width = 3, height = 5)

# need metrics
plot_count_lm(te_cancer_hostseq, type=NA, group="cancer_cohort", y_lab="LINE1 count", covariates = covar_med, residuals=FALSE, log_scale=FALSE, x_lab=x_tp53, min_samples=5, chr=NA)

# affected
generate_plots(plot_count_wilcox, te_aff_common, chr=NA, type=types, group="TP53_status", x_lab="Disease state")
plot_count_wilcox(te_aff, type=NA, group="TP53_status", y_lab="LINE1 count", chr=NA, x_lab=x_tp53, log_scale=FALSE)
ggsave(paste0(plot_dir, "te_count_wilcox_aff_line.png"), width = 3, height = 5)

plot_count_lm(te_aff, type=NA, group="TP53_status", y_lab="LINE1 count", covariates = covar_med, residuals=FALSE, log_scale=FALSE, x_lab=x_tp53, min_samples=5, chr=NA)
ggsave(paste0(plot_dir, "te_count_lm_aff_sva.png"), width = 3, height = 5)
bootstrap_test(te_aff, "TP53_status", "total", mean, n_bootstraps = 10000, step_size=100) 

#lfs
generate_plots(plot_count_wilcox, te_lfs, chr=NA, type=types, group="Cancer", x_lab="Cancer")
plot_count_wilcox(te_lfs, type=NA, chr=NA, group="Cancer", x_lab="Cancer", y_lab="Total TE count")

# sequencing cohort
te_lfs_sequencing <- te_aff %>%
  filter( # filter only row rows with it
    (cohort == "KiCS" & TP53_status == "Mutant") |
      (cohort == "LFS") |
      (cohort == "SJ")
  ) %>%
  mutate( # specify new column
    lfs_center = ifelse(
      cohort == "KiCS" & TP53_status == "Mutant", "KiCS",
      ifelse(cohort == "LFS", "LFS", "SJ")
    )
  ) 
generate_plots(plot_count_kruskal, te_lfs_sequencing, type=types, log_scale=TRUE, group="lfs_center", x_lab="Cluster", chr=NA)
plot_count_kruskal(te_lfs_sequencing, type=NA, y_lab="Rare repeat count", chr=NA, group="lfs_center", log_scale=TRUE, x_lab=x_tp53)
plot_count_lm(te_lfs_sequencing, type=NA, group="lfs_center", y_lab="LINE1 count", covariates = covar_med, residuals=FALSE, log_scale=FALSE, x_lab=x_tp53, min_samples=5, chr=NA)

# seuqence center
  
# sample center


#### COUNT BY TT ####
# count by tumour type
generate_plots(plot_count_kruskal_nogroup, te_aff, type=types, chr=NA, min=5) 
plot_count_kruskal_nogroup(te_kics, column="tumor_type", chr=NA, type=NA, y_lab="Rare repeat count", min=5,x_lab="Tumour type", log_scale = TRUE) 
ggsave(paste0(plot_dir, "te_count_tt_kics_all.png"), width = 9, height = 5)

# count by tumour type and grouped by group
generate_plots(plot_count_tt, te_aff, chr=NA, type=types, group="TP53_status", min=2) 
plot_count_tt(te_aff, chr=NA, type=NA, group="TP53_status", min=3, y_lab="Repeat count", legend_title=x_tp53, log_scale=TRUE) 
ggsave(paste0(plot_dir, "te_count_tt_tp53.png"), width = 9, height = 5)


#### COUNT BY CLINCAL VARIABLES####
# by cluster
generate_plots(plot_count_kruskal, te_lfs, chr=NA, type=types, group="cluster", x_lab="Cluster")
plot_count_kruskal(te_lfs, type=NA, chr=NA, group="cluster", x_lab="Cluster", y_lab="Repeat count", log_scale=TRUE)
ggsave(paste0(plot_dir, "te_count_lfs_cluster.png"), width = 9, height = 5)

# by inheritance
te_lfs <- te_lfs[!grepl("DELETE", te_lfs$inheritance), ]
generate_plots(plot_count_kruskal, te_lfs, chr=NA, type=types, group="inheritance", x_lab="Inheritance")
plot_count_kruskal(te_lfs, type="ALU", chr=NA, group="inheritance", x_lab="Inheritancee", y_lab="Total TE count", log_scale=TRUE)

# by variant location
generate_plots(plot_count_kruskal, te_lfs, chr=NA, type=types, group="Variant_location", x_lab="TP53 variant location")
plot_count_kruskal(te_lfs, type=NA, chr=NA, group="Variant_location", x_lab="TP53 variant location", y_lab="Total TE count")

# by variant classification 
generate_plots(plot_count_kruskal, te_lfs, chr=NA, type=types, group="Variant_Classification", x_lab="TP53 variant classification")
plot_count_kruskal(te_lfs, type=NA, chr=NA, group="Variant_Classification", x_lab="TP53 variant classification", y_lab="Total TE count")

# count by age
generate_plots(plot_count_age, te_kics, chr=NA, type=types) 
plot_count_age(te_kics, chr=NA, type=NA, y_lab="Total TE count")

# treatment
generate_plots(plot_count_wilcox, te_aff, chr=NA, type=types, group="treatment", x_lab="Treatment status")
plot_count_wilcox(te_aff, type=NA, chr=NA, group="treatment", x_lab="Treatment status", y_lab="log(Total TE count + 1)", log_scale=TRUE)

plot_count_lm(te_aff, min_samples=5, residuals=FALSE, log_scale=TRUE, group="treatment", covariates = covar_med, x_lab="Disease state", y_lab="Total TE count", type=NA, chr=NA)
bootstrap_test(te_aff, "treatment", "total", mean, n_bootstraps = 10000, step_size=100) 

# lesion type
generate_plots(plot_count_kruskal, te_aff, chr=NA, type=types, group="lesion_type", x_lab="Lesion type", log_scale=TRUE) 
plot_count_kruskal(te_aff, type=NA, chr=NA, group="lesion_type", x_lab="Treatment status", y_lab="log(Total TE count + 1)", log_scale=TRUE)

plot_count_lm(te_aff, min_samples=5, residuals=FALSE, log_scale=TRUE, group="lesion_type", covariates = covar_med, x_lab="Disease state", y_lab="Total TE count", type=NA, chr=NA)
bootstrap_test(te_aff, "lesion_type", "total", mean, n_bootstraps = 10000, step_size=100) 

# disease state
generate_plots(plot_count_kruskal, te_aff, chr=NA, type=types, group="disease_state", x_lab="Disease state", log_scale=TRUE) 
plot_count_kruskal(te_aff, type=NA, chr=NA, group="disease_state", x_lab="Disease state", y_lab="TE count", log_scale=TRUE)

plot_count_lm(te_aff, min_samples=5, residuals=FALSE, log_scale=TRUE, group="disease_state", covariates = covar_med, x_lab="Disease state", y_lab="Total TE count", type=NA, chr=NA)
bootstrap_test(te_aff, "disease_state", "total", mean, n_bootstraps = 10000, step_size=100) 

# sex
plot_count_lm(te_aff, min_samples=5, residuals=FALSE, group="sex", log_scale=TRUE, covariates = covar_med, x_lab="Sex", y_lab="Total TE count", type=NA, chr=NA)
bootstrap_test(te_aff, "sex", "total", mean, n_bootstraps = 10000, step_size=100) 


#### TE SUBFAMILY ####
te_cancer_hostseq_split <- te_all_split %>%
  filter(cohort == "HostSeq") %>%  # Filter rows where cohort is HostSeq
  bind_rows(te_aff_split)  # Combine with te_aff dataframe

# plot
plot_subfamily_pie(te_aff_split, "Alu", 1)

# fisher test
r <- perform_fisher_test_summary(te_aff_split, "subfamily", "TP53_status")
r %>% filter(P_Value < 0.05)
r <- perform_fisher_test_summary(te_cancer_hostseq_split, "subfamily", "cancer_cohort")

# compare number of young
unique(te_cancer_hostseq_split$subfamily)

young_family <- c("L1HS", "AluY", "SVA_E", "SVA_F")

# Filter the data frame for subfamilies in the young_family list
te_cancer_hostseq_split_young <- te_cancer_hostseq_split %>%
  filter(subfamily %in% young_family)

# process combinations
te_cancer_hostseq_split_young_processed <- process_all_combinations(te_cancer_hostseq_split_young)

# merge with clinical
te_cancer_hostseq_split_young_processed <- merge_dfs(te_cancer_hostseq_split_young_processed, clinical, include_all_x = FALSE)

# plot
generate_plots(plot_count_wilcox, te_aff_split_young_processed, chr=NA, type=types, group="TP53_status", x_lab="Disease state")
plot_count_wilcox(te_aff, type="LINE1", group="TP53_status", y_lab="LINE1 count", chr=NA, x_lab=x_tp53, log_scale=FALSE)
ggsave(paste0(plot_dir, "te_count_wilcox_aff_line.png"), width = 3, height = 5)

# plot
plot_count_kruskal(te_cancer_hostseq_split_young_processed, type=NA, chr=NA, group="cancer_cohort", log_scale=TRUE, x_lab="Cluster", y_lab="Young element count")
ggsave(paste0(plot_dir, "te_count_young.png"), width = 3, height = 5)

#### FULL LENGTH ####
# df with kics, lfs with cancer, and hostseq
te_cancer_hostseq_expand<- te_all_expand %>%
  filter(cohort == "HostSeq") %>%  # Filter rows where cohort is HostSeq
  bind_rows(te_aff_expand)  # Combine with te_aff dataframe

te_cancer_hostseq_expand_line <- te_cancer_hostseq_expand %>% filter(ALT=="LINE1")

summary(te_cancer_hostseq_expand_line$SV_length)

te_cancer_hostseq_expand_line_fulllength <- te_cancer_hostseq_expand_line %>% filter(SV_length >= 5900)
nrow(te_cancer_hostseq_expand_line_fulllength)

# process combinations
te_cancer_hostseq_expand_line_fulllength_processed <- process_all_combinations(te_cancer_hostseq_expand_line_fulllength)
te_cancer_hostseq_expand_line_fulllength_processed <- as.data.frame(add_nohit_samples(te_cancer_hostseq_expand_line_fulllength_processed, nohits))

# merge with clinical
te_cancer_hostseq_expand_line_fulllength_processed <- merge_dfs(te_cancer_hostseq_expand_line_fulllength_processed, clinical, include_all_x = FALSE)

# plot
plot_count_kruskal(te_cancer_hostseq_expand_line_fulllength_processed, type=NA, chr=NA, group="cancer_cohort", log_scale=TRUE, x_lab="Cluster", y_lab="Full-length LINE1 count")
ggsave(paste0(plot_dir, "te_count_full_length.png"), width = 3, height = 5)

#### P53 MUTATION and FITNESS ####
# p53 mutation
# count by mutaiton
plot_count_kruskal_nogroup(te_lfs, column="mutation", min=3, chr=NA, type=NA, x_lab="p53 mutation", y_lab="Total TE count", log_scale = TRUE) 

# hotspot
hotspot <- c(175, 213, 245, 248, 273, 282, 337)

# plot te frequency with mutation along p53 gene
plot_te_locations(te_lfs, chr=NA, type=NA, hotspot=hotspot, log_scale=FALSE)
plot_te_locations(te_lfs, "ALU", hotspot)
plot_te_locations(te_lfs, "LINE1")


# fitness
tp53_fitness <- fread("/Users/briannelaverty/Documents/R_Malkin/TE/data/tp53_variants_predictions.tsv", sep="\t", header=TRUE)

# add fitness values
te_aff_tp53 <- prep_p53_fitness(te_aff, tp53_fitness)
te_lfs_tp53 <- prep_p53_fitness(te_lfs, tp53_fitness)

# scatter plot of count vs fitness
generate_scatterplots(scatter_template_lm_one_variable, te_aff_tp53, independent="p53_fitness", types=types, x_limit=c(0,1), y_limit=c(0,75), x_lab="p53 fitness")    
scatter_template_lm_one_variable(te_aff_tp53, "p53_fitness", "total", c(0.35,1), c(800,1700), "p53 fitness", "Repeat count", shape_column = NULL, colour_column = NULL)
ggsave(paste0(plot_dir, "te_count_lfs_fitness.png"), width = 9, height = 5)

# inheritance  
plot_box_kruskal(te_lfs_tp53, "p53_fitness", "inheritance", "Inheritance", "p53 fitness")


#### CANCER GENES####
cpg<- read.csv("/Users/briannelaverty/Documents/R_Malkin/te/data/raw/kics_cpg.txt", header=FALSE, sep="\t") # cancer predisposition
genes <- cpg$V1


# overall cancer genes effected
# process combinations
te_cancer_hostseq_split <- te_all_split %>%
  filter(cohort == "HostSeq") %>%  # Filter rows where cohort is HostSeq
  bind_rows(te_aff_split)  # Combine with te_aff dataframe

te_cancer_hostseq_split_cancergenes  <- te_cancer_hostseq_split %>% filter(Gene_name %in% genes) # this is counting cancer genes effected not number of TEs that effect cancer genes as one TE may effect multiple cancer genes
te_cancer_hostseq_split_cancergenes_processed <- process_all_combinations(te_cancer_hostseq_split_cancergenes) 
te_cancer_hostseq_split_cancergenes_processed <- as.data.frame(add_nohit_samples(te_cancer_hostseq_split_cancergenes_processed, nohits))
summary(te_cancer_hostseq_split_cancergenes_processed$total)

# merge with clinical and metrics
te_cancer_hostseq_split_cancergenes_processed <- merge_dfs(te_cancer_hostseq_split_cancergenes_processed, clinical, include_all_x = FALSE)
te_cancer_hostseq_split_cancergenes_processed <- merge_dfs(te_cancer_hostseq_split_cancergenes_processed, metrics, include_all_x = FALSE)

# sig test
plot_count_kruskal(te_cancer_hostseq_split_cancergenes_processed, type=NA, group="cancer_cohort", x_lab=x_tp53, y_lab="Number of cancer genes effected", chr=NA, log_scale=TRUE)
plot_count_lm(te_cancer_hostseq_split_cancergenes_processed, type=NA, covariates = "tumor_type", group="cancer_cohort", y_lab="Number of cancer genes effected", residuals=FALSE, log_scale=FALSE, x_lab=x_tp53, min_samples=5, chr=NA)
bootstrap_test(te_cancer_hostseq_split_cancergenes_processed, "TP53_status", "total", mean, n_bootstraps = 10000, step_size=100) 


# top genes
# add gene size
te_kics_split <- add_gene_size_todf(te_kics_split, gene_size)
geneList_kics <- calculate_gene_mut_sample_frequency(te_kics_split)
geneList_kics_filtered <- geneList_kics[geneList_kics$Gene_name %in% genes, ]
names(geneList_kics)

# plot top genes
plot_top_genes(geneList_kics_filtered, column="freq_samples_effected", label_column=NULL, x_lab="Proportion of samples with insertion in gene", top_n=20)
ggsave(paste0(plot_dir, "te_top_cancer_genes_affected_kics.png"), width = 9, height = 5)
plot_top_genes(geneList_kics_filtered, column="freq_samples_effected_normalized", label_column="num_samples_effected", x_lab="Proportion of samples with insertion in gene normalized for gene size", top_n=20)
ggsave(paste0(plot_dir, "te_top_cancer_genes_affected_normalized_kics.png"), width = 9, height = 5)

# patients with particular gene affected
patients_mutyh <- as.data.frame(te_kics_split %>% 
                                  filter(Gene_name=="LRP1B") %>%
                                  group_by("sample"))
table(patients_mutyh$tumor_type)
table(patients_mutyh$TP53_status)


# genes present that meet critiera
# comparing gene instances
location=c("CDS", "UTR", "3'UTR", "5'UTR", "5'UTR-CDS","CDS-3'UTR")

genes_aff <- identify_te_genes(te_aff_split_genes, gene_vector=genes, location="exon", location2=NULL)
genes_kics <- identify_te_genes(te_kics_split, gene_vector=genes, location=NULL, location2=NULL)
genes_lfs <- identify_te_genes(te_lfs_split, gene_vector=genes, location="exon", location2=NULL)

table(genes_aff$Gene_name)
table(genes_kics$Gene_name)
table(genes_lfs$Gene_name)

plot_top_cancer_genes_germline(genes_aff, top_n=15)
plot_top_cancer_genes_germline(genes_kics, top_n=15)
plot_top_cancer_genes_germline(genes_lfs, top_n)


# wilcox test b/w lfs and controls
count_location_wilcox(te_aff_split, gene=genes, filter_element=NA, group="TP53_status", location=NULL)

# wilcox test b/w lfs aff and unaff
count_location_wilcox(te_lfs_split_rare, gene=genes, filter_element=NA, group="Cancer", location=location)


#### TOP MUTATED GENES ####
# kics
te_kics_split <- add_gene_size_todf(te_kics_split, gene_size)
geneList_kics <- calculate_gene_mut_sample_frequency(te_kics_split)
names(geneList_kics)

# plot top genes
plot_top_genes(geneList_kics, column="freq_samples_effected", label_column="freq_mutations", x_lab="Proportion of samples with insertion in gene normalized by gene size", top_n=20)
ggsave(paste0(plot_dir, "te_top_genes_affected_normalized_kics.png"), width = 9, height = 5)

# number of genes affecting > 50% of samples
count_genes_by_threshold(geneList_kics, column_name = "freq_samples_effected", thresholds = c(0.5, 0.75, 0.9))


# lfs
te_lfs_split <- add_gene_size_todf(te_lfs_split, gene_size)
geneList_lfs<- calculate_gene_mut_sample_frequency(te_lfs_split, nsample_thresh = 0, filter_exon = FALSE)

# plot top genes
plot_top_genes(geneList_lfs, column="freq_samples_effected_normalized", label_column="freq_samples_effected", x_lab="Proportion of samples with insertion in gene normalized by gene size", top_n=20)
ggsave(paste0(plot_dir, "te_top_genes_affected_normalized_lfs.png"), width = 9, height = 5)


# control
te_control_split <- te_aff_split %>% filter(TP53_status == "Control")
geneList_tumour_control<- calculate_gene_mut_sample_frequency(te_control_split, nsample_thresh = 0, filter_exon = FALSE)
plot_top_genes(geneList_tumour_control, 20)


# specific gene heatmap of insertions
plot_te_insertion_complex_heatmap(te_kics_split, gene_name="LRP1B", bin_size=1000, annotations=c("tumor_class", "age_at_diagnosis", "TP53_status"), lwd=0)
#ggsave(paste0(plot_dir, "te_location_kics_BCL2.png"), width = 9, height = 5)

te_aff_split_genes %>% filter(Gene_name=="BCL2") # intron, no frameshift for BCL2 and LRP1B


#### PATHWAY PEDIATRIC ####
# over representation analysis
ora_kics <- perform_ora(te_kics_split, nsample_thresh = 0, filter_exon = FALSE)

# viz 
barplot(ora_kics, showCategory=20) 
dotplot(ora_kics, showCategory=20) 
goplot(ora_kics, showCategory=5)
cnetplot(ora_kics, showCategory = 10, colorEdge = TRUE, node_label = "category") 
ggsave(paste0(plot_dir, "te_pathway_ora_kics.png"), width = 14, height = 9)
#ora_kics_pairwise<- pairwise_termsim(ora_kics)
#emapplot(ora_kics_pairwise, group_category = TRUE, group_legend = TRUE) +
#  scale_fill_identity(values = "#0080A3") +
#  guides(fill = "none")

# group descriptions by broader function
ora_kics_simple <- simplify(ora_kics, cutoff=0.6, by="p.adjust", select_fun=min)
barplot(ora_kics_simple, showCategory=20) 
cnetplot(ora_kics_simple, showCategory = 10, colorEdge = TRUE, node_label = "category")
ggsave(paste0(plot_dir, "te_pathway_kics_simplified.png"), width = 14, height = 9)


# common genes in pathway
# freq of genes in top pathways
pathway_genes <- plot_pathway_gene_counts(ora_kics, n_descriptions=nrow(ora_kics), n_genes=50)

# look at correlation b/w top genes mutated and genes in top pathways
# column is num_samples_effected or mutation_count
analyze_gene_mutations(geneList_kics, ora_kics, n_descriptions=nrow(ora_kics), column = "num_samples_effected", y_lab="Number of samples with insertion in gene")
ggsave(paste0(plot_dir, "te_pathwaygenes_vs_numsamples_kics.png"), width = 9, height = 5)


#### PATHWAY TP53 ####
ora_tp53 <- perform_ora_tp53(te_aff_split)

# viz 
dotplot(ora_tp53, showCategory = 20)
cnetplot(ora_tp53, showCategory = 10, colorEdge = TRUE, node_label = "category")
ora_tp53_pairwise<- pairwise_termsim(ora_tp53)
emapplot(ora_tp53_pairwise, pie="count", showCategory = 20, group_category=TRUE, group_legend=TRUE) + scale_fill_manual(values=colours)
ggsave(paste0(plot_dir, "te_pathway_tp53.png"), width = 14, height = 9)

# simplify
ora_tp53_simple <- simplify(ora_tp53, cutoff=0.5, by="p.adjust", select_fun=min)
dotplot(ora_tp53_simple, showCategory=20) 
cnetplot(ora_tp53_simple, showCategory = 10, colorEdge = TRUE, node_label = "category")


# common genes
# enrichment result
ora_tp53_result <- ora_tp53@compareClusterResult

# GO terms that were only enriched in KiCS or LFS not both
ora_tp53_result_unique <- ora_tp53_result %>%
  group_by(ID) %>%
  filter(n() == 1) %>%
  ungroup()

# split enrichment by group
ora_tp53_result_control <- ora_tp53_result_unique %>% filter(Cluster=="Control")
ora_tp53_result_lfs <- ora_tp53_result_unique %>% filter(Cluster=="LFS")

# freq of genes in top pathways
pathway_genes_control <- plot_pathway_gene_counts(ora_tp53_result_control, n_descriptions=nrow(ora_tp53_result_control), n_genes=50)
pathway_genes_lfs <- plot_pathway_gene_counts(ora_tp53_result_lfs, n_descriptions=nrow(ora_tp53_result_lfs), n_genes=50)

# look at correlation b/w top genes mutated and genes in top pathways
# column is num_samples_effected or mutation_count
analyze_gene_mutations(geneList_tumour_control, pathway_genes_control, column = "freq_samples_effected") # this is different than the kics one because it is only pathways unique to controls not all sig enriched pathways
analyze_gene_mutations(geneList_tumour_lfs, pathway_genes_lfs, column = "freq_samples_effected")


#### PATHWAY PED CANCER VS HOST SEQ ####
# make df
te_cancer_hostseq_split <- te_all_split %>%
  filter(cohort == "HostSeq") %>%  # Filter rows where cohort is HostSeq
  bind_rows(te_kics_split)  # Combine with te_aff dataframe

ora_cancer <- perform_ora_cancer(te_cancer_hostseq_split)

# viz 
dotplot(ora_cancer, showCategory = 20)
cnetplot(ora_cancer, showCategory = 10, colorEdge = TRUE, node_label = "category")
ora_cancer_pairwise<- pairwise_termsim(ora_cancer)
emapplot(ora_cancer_pairwise, pie="count", showCategory = 20, group_category=TRUE, group_legend=TRUE) + scale_fill_manual(values=colours)
ggsave(paste0(plot_dir, "te_pathway_cancer.png"), width = 14, height = 9)




#### PATHWAY LFS AFFECTED ####
#### UNIQUE TE IN LFS ####
specific_te_fisher <- as.data.frame(fisher_test_unique_te(te_aff_expand, min=3))
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
table(one$sample, one$tumor_type)
two <- sig_te_samples[[2]]
table(two$sample, two$tumor_type)
three <- sig_te_samples[[3]]
table(three$sample, three$tumor_type)


# genes effected
# two TE dont effect gene
# one te in intron of NUMB, no frameshift, and gnomad and exac pLI close to 0 meaning tolerant of LOF variants 
sig_te_genes<- find_sig_te_genes(sig_te, te_aff_split_genes)

# viz
out <- plot_proportions(sig_te)
out[[1]]
ggsave(paste0(plot_dir, "te_unique_1.png"), width = 3, height = 5)
out[[2]]
ggsave(paste0(plot_dir, "te_unique_2.png"), width = 3, height = 5)
out[[3]]
ggsave(paste0(plot_dir, "te_unique_3.png"), width = 3, height = 5)


#### GERMLINE AND SOMATIC VARIANTS ####
# load variants, add variants to te_aff, make overall TP53 column
te_aff_variant <- prep_variant_df(te_kics, te_aff)

generate_plots(plot_count_kruskal, te_aff_variant, chr=NA, type=types, group="TP53", x_lab="Cluster")
plot_count_kruskal(te_aff_variant, group="TP53", x_order=c("None", "Somatic", "Germline"), x_lab="TP53 variant", y_lab="Total TE count", type=NA, chr=NA, log_scale=TRUE)


#### DESCRIBE DATASET ####
# extract clinical from samples with have
clin<- te_all %>% dplyr::select(sample, tumor_type, sex, age_at_diagnosis, TP53_status, cohort)
clin <- clin %>% distinct()

# Plot pie chart for tumor_type
plot_pie_chart(clin, "tumor_type", "Tumour type")
ggsave(paste0(plot_dir, "clinical_germline_tt.png"), width = 6, height = 5)

# Plot pie chart for TP53_status
tp53_title <- expression("Germline " * italic("TP53") * " status")
plot_pie_chart(clin, "TP53_status", tp53_title)
ggsave(paste0(plot_dir, "clinical_germline_tp53.png"), width = 6, height = 5)

# Plot pie chart for sex 
plot_pie_chart(clin, "sex", "Sex")
ggsave(paste0(plot_dir, "clinical_sex.png"), width = 6, height = 5)

# Plot histogram for age_at_diagnosis
plot_histogram_age(clin, bin=2)
ggsave(paste0(plot_dir, "clinical_germline_age.png"), width = 9, height = 5)

# cohort
clin <- clin %>%
  mutate(new_cohort = case_when(
    cohort == "LFS_mut" ~ "LFS",  
    cohort == "nick" ~ "LFS",    
    TRUE ~ cohort               
  ))
plot_pie_chart(clin, "new_cohort", "Dataset")
ggsave(paste0(plot_dir, "clinical_dataset.png"), width = 9, height = 5)


#### ANCESTRY VIZ ###
# overall viz ####
# te aff
# pca
te_aff_viz <- te_aff[,c(28:51, 53:76, 78:101)]
te_aff_viz_clean <- te_aff_viz[, apply(te_aff_viz, 2, var) != 0] # remove all 0 columns

pca_te_aff <- prcomp(te_aff_viz_clean, scale=TRUE)

# make df for graphing 
pca_df <- data.frame(PC1 = pca_te_aff$x[, 1],
                     PC2 = pca_te_aff$x[, 2],
                     sample = te_aff$sample)

# merge with clinical 
pca_df <- merge(pca_df, te_aff, by="sample")

# plot
ggplot(pca_df, aes(x = PC1, y = PC2, color = predicted_ancestry_thres)) +
  geom_point(alpha=0.8, size=5) +  # Scatter plot
  geom_text(aes(label=sample), vjust=2, hjust=-2, size=3, color="blue") +  # Labels
  labs(x = "Principal Component 1",
       y = "Principal Component 2")


# umap
te_aff_viz_clean_scaled <- scale(te_aff_viz_clean)

umap_result <- umap::umap(te_aff_viz_clean_scaled,
                          n_neighbors = 25, 
                          min_dist = 0.2, 
                          n_components = 2)

# Create a dataframe for plotting
umap_df <- data.frame(UMAP1 = umap_result$layout[,1],
                      UMAP2 = umap_result$layout[,2],
                      sample = te_aff$sample)

# merge with clinical
umap_df <- merge(umap_df, clinical, by="sample")

# plot
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = TP53_status)) +
  geom_point(alpha=0.8, size=5) +  # Scatter plot
  geom_text(aes(label=sample), vjust=2, hjust=-2, size=3, color="blue") +  # Labels
  labs(x = "Component 1",
       y = "Component 2")


# te split
te_viz <- te_aff_split %>%
  group_by()


### end ####


