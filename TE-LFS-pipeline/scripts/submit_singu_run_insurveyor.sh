#!/bin/bash

#SBATCH -J ins_runs
#SBATCH -e %x.e%j
#SBATCH -o %x.o%j
#SBATCH -N 1 -c 16
#SBATCH --mem=64G
#SBATCH -t 47:59:00

set -euxo 

module load Singularity/3.11.3

SINGULARITY_IMAGE="/hpf/largeprojects/davidm/shilpa/TE-tools/insurveyor.sif"
bamfilepath=$1
outputfoldername=$2
reference_fasta="/hpf/largeprojects/davidm/shilpa/TE-tools/MELTv2.2.2/Demo/hs37d5.fa"

mkdir -p $outputfoldername
export TMPDIR=$SCRATCH

singularity run -B /hpf:/hpf $SINGULARITY_IMAGE --threads 16 $bamfilepath $outputfoldername $reference_fasta
