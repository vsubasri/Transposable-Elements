#!/bin/bash

# Transposable Elements Analysis Pipeline GERMLINE SAMPLES

# Usage: ./TE-LFS-germline-pipe.sh /path/to/bamdir (this would be used as $BAM_DIR)

#########################################
# Step 0: Setup and Configuration
#########################################

# Define the directory containing BAM files
BAM_DIR=$1
# default BAM_DIR path is mentioned

TEPIPE_DIR="/hpf/largeprojects/davidm/shilpa/TE-LFS-pipeline"

# Output directories
OUTPUT_DIR=$PWD/results
BAM_fixed=$PWD/fixed_bams
XTEA_DIR=$OUTPUT_DIR/xtea_runs
MELT_DIR=$OUTPUT_DIR/melt_runs
INS_DIR=$OUTPUT_DIR/ins_runs
ANALY_DIR=$OUTPUT_DIR/analysis
SURVIVOR_PATH="/hpf/largeprojects/davidm/shilpa/TE-tools/SURVIVOR-master/Debug/SURVIVOR"

mkdir -p $XTEA_DIR $MELT_DIR $INS_DIR $BAM_fixed $XTEA_DIR/logs $MELT_DIR/logs $INS_DIR/logs $BAM_fixed/logs $INS_DIR/logs

# Create bamfilepaths.txt bamfile_desc.txt and sample_id.txt
find "$BAM_DIR" -name "*.bam" > rawbamfilepaths.txt

#########################################
# Step 1: File Preprocessing
#########################################

#########################################
# Step 1: File Preprocessing
#########################################
cp rawbamfilepaths.txt $TEPIPE_DIR/scripts/submit_preprocessbam.sh $TEPIPE_DIR/scripts/process_metrics.sh $BAM_fixed

(
cd $BAM_fixed || exit
bam_preprocess_job_ids=()
while read -r bamfilepath; do
        output=$(sbatch --parsable $PWD/submit_preprocessbam.sh "$bamfilepath")
        echo "sbatch output: $output"  # Debug output
        job_id=$(echo "$output" | awk '{print $1}')
	echo "job_id is $job_id"
        if [ -n "$job_id" ]; then
            bam_preprocess_job_ids+=("$job_id")
        else
            echo "Error: Failed to capture job ID for $bamfilepath"
        fi
done < rawbamfilepaths.txt
)

job_ids_str=$(IFS=,; echo "${bam_preprocess_job_ids[*]}")

if [ -z "$job_ids_str" ]; then
    echo "Error: No job IDs captured, exiting."
fi

echo "Captured Job IDs: ${bam_preprocess_job_ids[@]}"
echo "Job ID string: $job_ids_str"
#########################################

find "$BAM_fixed" -name "sorted_fixed_*_N.bam" > bamfilepaths.txt
awk -F"/" '{print $NF" "$0}' bamfilepaths.txt | sed 's/.bam//' > bamfile_desc.txt
awk '{print $1}' bamfile_desc.txt > sample_id.txt
#########################################
# Step 2: MELT
#########################################

process_melt() {
    local bamfilepath=$1
    local sample_id=$(basename "$bamfilepath" | awk -F"." '{print $1}')
    local sample_output_dir=$MELT_DIR/$sample_id
    mkdir -p $sample_output_dir
    local job_id=$(sbatch --parsable --dependency=afterok:$job_ids_str submit_meltrun.sh "$bamfilepath" "$sample_output_dir")
    if [ $? -eq 0 ]; then
        melt_job_ids+=("$job_id")
    else
        echo "Error: Failed to submit MELT job for $bamfilepath" >&2
    fi
}

cp $TEPIPE_DIR/scripts/submit_meltrun.sh bamfilepaths.txt $MELT_DIR

(
cd $MELT_DIR || exit
while read -r bamfilepath; do
    process_melt "$bamfilepath"
done < bamfilepaths.txt
)

if [ ${#melt_job_ids[@]} -eq 0 ]; then
    echo "Error: No MELT jobs were submitted, exiting." >&2
    exit 1
fi

echo "MELT jobs submitted: ${melt_job_ids[@]}"
##########################################

