#!/bin/bash
#SBATCH --job-name=totalrecall       # name of the job array
#SBATCH --output=logs/%A_%a.tmp.out  # Temporary log placeholder
#SBATCH --error=logs/%A_%a.tmp.err   # Temporary log placeholder
#SBATCH --mem=64G
#SBATCH -t 47:00:00
#SBATCH --qos=shilpa_q
#SBATCH -N 1 -c 16
#SBATCH --array=1-4 #get num from wc -l PAIR_FILE

module load Singularity/3.11.3

# Paths and files
PAIR_FILE="NT-pairs.csv"

IMAGE_FILE="/hpf/largeprojects/davidm/shilpa/TE-tools/totalrecall-hg19.sif"

# Extract sample id, and tumor file, germline info for this array task, racall config uses names as case control
bn=$(awk "NR==${SLURM_ARRAY_TASK_ID}" $PAIR_FILE | awk '{print $1}')
ft=$(awk "NR==${SLURM_ARRAY_TASK_ID}" $PAIR_FILE | awk '{print $2}')
fn=$(awk "NR==${SLURM_ARRAY_TASK_ID}" $PAIR_FILE | awk '{print $3}')

echo "Patient ID: $bn"
echo "Case BAM Path: $ft"
echo "Control BAM Path: $fn"

# Redirect logs to dynamically named files based on sample ID
exec 1>logs/${bn}.out
exec 2>logs/${bn}.err

# Create and move into the directory for this specific patient
mkdir -p "$bn" "logs"
cd "$bn"

# Link case and control BAM files as expected by Snakemake
ln -s "${ft}" case.bam
ln -s "${ft}.bai" case.bam.bai
ln -s "${fn}" control.bam
ln -s "${fn}.bai" control.bam.bai

# Run Snakemake using the singularity image
singularity exec --no-home --bind /hpf:/hpf "$IMAGE_FILE" snakemake -j16 -p -s /opt/totalrecall/Snakefile results.tgz
