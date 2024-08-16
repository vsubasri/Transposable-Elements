#!/bin/bash

#SBATCH -J traficrun
#SBATCH -e %x.e%j
#SBATCH -o %x.o%j
#SBATCH -N 1 -c 16
#SBATCH --mem=64G
#SBATCH -t 24:59:00

set -euxo 

module load Singularity/3.11.3

# Define your data and singularity image paths
SINGULARITY_IMAGE="/hpf/largeprojects/davidm/shilpa/TE-tools/trafic_dockerimage.sif"
CONFIG_FILE=$1
OUTPUT_RESULTS=$2
INPUTBAMFIXED="/hpf/largeprojects/davidm/shilpa/TE_runs/fixed_bams"

singularity exec --bind ${INPUTBAMFIXED}:/trafic_input --bind $OUTPUT_RESULTS:/trafic_output --home $SCRATCH:/root ${SINGULARITY_IMAGE} /usr/bin/snakemake --snakefile /trafic/workflows/SnakefileTraficHg19 --configfile /trafic_input/${CONFIG_FILE} --verbose

