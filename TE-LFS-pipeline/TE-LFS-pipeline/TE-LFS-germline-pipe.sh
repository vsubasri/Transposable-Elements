#!/bin/bash

# Transposable Elements Analysis Pipeline GERMLINE SAMPLES

# Usage: ./TE-LFS-germline-pipe.sh /path/to/bamdir (this would be used as $BAM_DIR)

#########################################
# Step 0: Setup and Configuration
#########################################

# Define the directory containing BAM files
BAM_DIR="${1:-/hpf/largeprojects/davidm/data/te_test/germline_bam/}"
# default BAM_DIR is also in case

# output directories
TEPIPE_DIR="/hpf/largeprojects/davidm/shilpa/TE-LFS-pipeline"
DATA=$TEPIPE_DIR/data
OUTPUT_DIR=$PWD/results
BAM_fixed=$PWD/fixed_bams
XTEA_DIR=$OUTPUT_DIR/xtea_runs
MELT_DIR=$OUTPUT_DIR/melt_runs
INS_DIR=$OUTPUT_DIR/ins_runs

mkdir -p $XTEA_DIR $MELT_DIR $INS_DIR $BAM_fixed $XTEA_DIR/logs $MELT_DIR/logs $INS_DIR/logs $BAM_fixed/logs "logs"

# Create bamfilepaths.txt bamfile_desc.txt and sample_id.txt
ls $BAM_DIR/*.bam > $DATA/rawbamfilepaths.txt

#########################################
# Step 1: File Preprocessing
#########################################
cp $DATA/rawbamfilepaths.txt $TEPIPE_DIR/scripts/submit_preprocessbam.sh $BAM_fixed
(
cd $BAM_fixed || exit
while read -r bamfilepath; do
    sbatch submit_preprocessbam.sh "$bamfilepath"
done < rawbamfilepaths.txt
)
# Wait for preprocessing to complete
while squeue -u $USER | grep -q 'bam_preprocess'; do
    echo "**waiting for the bam preprocessing steps ..."
    sleep 1000
done
#########################################

ls $BAM_fixed/sorted_fixed*.bam > $DATA/bamfilepaths.txt
awk -F"/" '{print $NF" "$0}' $DATA/bamfilepaths.txt >$DATA/bamfile_desc.txt
awk -F".bam" '{print $1}' $DATA/bamfile_desc.txt >$DATA/sample_id.txt

#########################################
# Step 2: xTea
#########################################

cp $DATA/bamfile_desc.txt $DATA/sample_id.txt $TEPIPE_DIR/scripts/run_gnrt_pipeline_hg19.sh $TEPIPE_DIR/scripts/interm_set_prep_sbatch.sh $TEPIPE_DIR/scripts/slurm_header_xtea.txt $XTEA_DIR
(
cd $XTEA_DIR || exit
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
cp $TEPIPE_DIR/scripts/submit_meltrun.sh $DATA/bamfilepaths.txt $MELT_DIR
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
cp $TEPIPE_DIR/scripts/submit_singu_run_insurveyor.sh $DATA/bamfilepaths.txt $INS_DIR

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

