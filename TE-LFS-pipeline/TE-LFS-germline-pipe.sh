#!/bin/bash

# Transposable Elements Analysis Pipeline GERMLINE SAMPLES

# Usage: ./TE-LFS-germline-pipe.sh /path/to/bamdir

#########################################
# Step 0: Setup and Configuration
#########################################

# filepaths.txt has paths for storing original bam file data, in case of germline pipeline it should have only 1 line for one folder. lets do bamfilepaths_NT.txt and bamfilepaths.txt

# Define the directory containing BAM files
BAM_DIR="${1:-/hpf/largeprojects/davidm/data/te_test/germline_bam/}"
# default BAM_DIR is also in case

# output directories
OUTPUT_DIR=$PWD/results
BAM_fixed=$OUTPUT_DIR/fixed_bams
XTEA_DIR=$OUTPUT_DIR/xtea_runs
MELT_DIR=$OUTPUT_DIR/melt_runs
INS_DIR=$OUTPUT_DIR/ins_runs

mkdir -p $XTEA_DIR $MELT_DIR $INS_DIR $BAM_fixed $XTEA_DIR/logs $MELT_DIR/logs $INS_DIR/logs $BAM_fixed/logs "logs"

# Create bamfilepaths.txt bamfile_desc.txt and sample_id.txt
ls $BAM_DIR/*.bam > data/bamfilepaths.txt
awk -F"/" '{print $NF" "$0}' data/bamfilepaths.txt >data/bamfile_desc.txt
awk -F".bam" '{print $1}' data/bamfile_desc.txt >data/sample_id.txt

#########################################
# Step 1: File Preprocessing
#########################################
cp data/bamfilepaths.txt $BAM_fixed
(
cd $BAM_fixed || exit
while read -r bamfilepath; do
    sbatch submit_preprocessbam.sh "$bamfilepath"
done < bamfilepaths.txt
)
# Wait for preprocessing to complete
while squeue -u $USER | grep -q 'bam_preprocess'; do
    sleep 100
done

#########################################
# Step 2: xTea
#########################################

cp data/bamfile_desc.txt data/sample_id.txt scripts/run_gnrt_pipeline_hg19.sh scripts/interm_set_prep_sbatch.sh scripts/slurm_header.txt $XTEA_DIR
(
cd $XTEA_DIR ||exit
	bash run_gnrt_pipeline_hg19.sh
	#prepating files for slurm job
	bash interm_set_prep_sbatch.sh
	bash submit_scripts.sh
)

#run_gnrt uses bamfile_desc.txt and creates folders corresponding to their sample_id

#########################################
# Step 3: MELT
#########################################

process_melt() {
    local bamfilepath=$1
    local sample_id=$(basename "$bamfilepath" | awk -F"." '{print $1}')
    local sample_output_dir=$MELT_DIR/$sample_id
    mkdir -p $sample_output_dir
    sbatch submit_meltrun.sh "$bamfilepath" "$sample_output_dir"
}
cp scripts/submit_meltrun.sh data/bamfilepaths.txt $MELT_DIR
(
cd $MELT_DIR || exit
while read -r bamfilepath; do
    process_melt $bamfilepath
done < bamfilepaths.txt
)

##########################################
# Step 4: INSurVeyor
##########################################

process_ins() {
    local bamfilepath=$1
    #local sample_id=$(awk -F".bam" '{print $1}' $bamfilepath)
    local sample_id=$(basename "$bamfilepath" | awk -F"." '{print $1}')
    local sample_output_dir=$INS_DIR/$sample_id
    mkdir -p $sample_output_dir
    sbatch submit_meltrun.sh "$bamfilepath" "$sample_output_dir"
}
cp scripts/submit_singu_run_insurveyor.sh data/bamfilepaths.txt $INS_DIR

(
cd $INS_DIR || exit
while read -r bamfilepath; do
	process_ins $bamfilepath
done < bamfilepaths.txt

)
#########################################
# Step 5: new tool
#########################################


#########################################
# Step 6: SURVIVOR
#########################################

#wait for TE-runs to complete
#while squeue -u $USER | grep -q 'xtea_runs | melt_runs| ins_runs'; do
#    sleep 100
#done

# List VCF files
#ls *vcf > sample_files.txt

# Merge VCF files with SURVIVOR
#./SURVIVOR merge sample_files.txt 1000 2 1 1 0 30 sample_merged.vcf
#   SURVIVOR merge sample_files 100 2 1 0 0 30 sample_merged.vcf
#ALU SVA L1
#########################################
# Step 7: Annotations
#########################################

# Load necessary modules
#module load AnnotSV
#module load bcftools
#module load bedtools

# Annotate merged VCF file
#AnnotSV -SVinputFile sample_merged.vcf

# This should create an output folder with annotated files.

#########################################

