#!/bin/bash

#SBATCH -J ins_run
#SBATCH -o logs/%x_output_%j.log
#SBATCH -e logs/%x_error_%j.log
#SBATCH -N 1 -c 16
#SBATCH --mem=64G
#SBATCH -t 47:59:00

set -euxo 

module load Singularity/3.11.3

SINGULARITY_IMAGE="/hpf/largeprojects/davidm/shilpa/TE-tools/insurveyor.sif"
bamfilepath=$1
outputfoldername=$2
#reference_fasta="/hpf/largeprojects/davidm/resources/GRCh38.primary_assembly.genome.fa"
reference_fasta="/hpf/largeprojects/davidm/resources/te_scripts/GRCh38_full_analysis_set_plus_decoy_hla.fa"
export TMPDIR=$SCRATCH
singularity run -B /hpf:/hpf $SINGULARITY_IMAGE --threads 16 $bamfilepath $outputfoldername $reference_fasta
