#!/bin/bash

# Transposable Elements prediction Pipeline for GERMLINE SAMPLES

# Usage: ./TE-LFS-germline-pipe.sh /path/to/bamdir (this would be used as $BAM_DIR)

#########################################
# Step 0: Setup and Configuration
#########################################

# Check if BAM_DIR argument is provided
if [ -z "$1" ]; then
  echo "Error: Missing argument for BAM_DIR."
  echo "Usage: $0 /path/to/bamdir"
  echo "Please provide the directory containing input BAM files as the first argument."
  exit 1
fi

BAM_DIR=$1 # the directory containing BAM files

TEPIPE_DIR="/hpf/largeprojects/davidm/shilpa/TE-LFS-pipeline" # the directory where pipeline is located

base_dir=$PWD # this is where all the output files will be created

echo "The output files will be created in the current path: $base_dir"

OUTPUT_DIR=$PWD/results_hg38
BAM_fixed=$PWD/fixed_bams_hg38
XTEA_DIR=$OUTPUT_DIR/xtea_runs
MELT_DIR=$OUTPUT_DIR/melt_runs
INS_DIR=$OUTPUT_DIR/ins_runs

mkdir -p $XTEA_DIR $MELT_DIR $INS_DIR $BAM_fixed $XTEA_DIR/logs $MELT_DIR/logs $INS_DIR/logs $BAM_fixed/logs $INS_DIR/logs

# Create bamfilepaths.txt bamfile_desc.txt and sample_id.txt
find "$BAM_DIR" -name "*.bam" > rawbamfilepaths.txt

#########################################
# Step 1: File Preprocessing
#########################################
cp rawbamfilepaths.txt $TEPIPE_DIR/scripts/submit_preprocessbam_hg38.sh $TEPIPE_DIR/scripts/process_metrics.sh $BAM_fixed
cd $BAM_fixed
bam_preprocess_job_ids=()
while read -r bamfilepath; do
	job_id=$(sbatch --parsable submit_preprocessbam_hg38.sh "$bamfilepath")
	bam_preprocess_job_ids+=("$job_id")
done < rawbamfilepaths.txt
cd ..

echo "Processing bam files ..."

job_ids_str=$(IFS=,; echo "${bam_preprocess_job_ids[*]}")
echo " job ids corresponding are: $job_ids_str"
#########################################

sed "s|.*/|$BAM_fixed/sorted_fixed_|" rawbamfilepaths.txt > bamfilepaths.txt #ensure the filenames are correct _N and _T would be used in tumor pipeline
awk -F"/" '{print substr($NF, 1, length($NF)-4) " " $0}' bamfilepaths.txt >bamfile_desc.txt
awk '{print $1}' bamfile_desc.txt >sample_id.txt

#########################################
# TE prediction tools
#########################################
# Step 2: xTea
#########################################
cp $base_dir/bamfile_desc.txt $base_dir/sample_id.txt $TEPIPE_DIR/scripts/run_gnrt_pipeline_hg38.sh $TEPIPE_DIR/scripts/interm_set_prep_sbatch.sh $TEPIPE_DIR/scripts/slurm_header_xtea.txt $XTEA_DIR

cd $XTEA_DIR
   bash run_gnrt_pipeline_hg38.sh
   #prepating files for slurm job
   bash interm_set_prep_sbatch.sh
   sed -i "s/sbatch < /sbatch --dependency=afterok:"\$1" /g" submit_scripts.sh #the TE calls will be submitted if the bam preprocessing runs are successful (afterok ensures this)
   bash submit_scripts.sh $job_ids_str
cd ../..

#run_gnrt uses bamfile_desc.txt and creates folders corresponding to their sample_id

echo "xTea runs submitted"

#########################################
# Step 3: MELT
#########################################

cp $TEPIPE_DIR/scripts/submit_meltrun_hg38.sh $base_dir/bamfilepaths.txt $MELT_DIR

cd $MELT_DIR
while read -r bamfilepath; do
    sbatch --dependency=afterok:$job_ids_str submit_meltrun_hg38.sh "$bamfilepath"
done < bamfilepaths.txt
cd ../..

echo "Melt runs submitted"

##########################################
# Step 4: INSurVeyor
##########################################

cp $TEPIPE_DIR/scripts/submit_singu_run_insurveyor_hg38.sh $base_dir/bamfilepaths.txt $INS_DIR

cd $INS_DIR
while read -r bamfilepath; do
	sbatch --dependency=afterok:$job_ids_str submit_singu_run_insurveyor_hg38.sh "$bamfilepath"
done < bamfilepaths.txt
cd ../..

echo "INSurVeyor runs submitted"
##########################################
squeue -u $USER
##########################################
