#!/bin/bash

#SBATCH -J melt_run
#SBATCH -o logs/%x_output_%j.log
#SBATCH -e logs/%x_error_%j.log
#SBATCH -N 1 -c 16
#SBATCH --mem=64G
#SBATCH -t 47:59:00

set -euxo
module load java/1.8.0_91 bowtie2 samtools

melt_full_path="/hpf/largeprojects/davidm/shilpa/TE-tools/MELTv2.2.2"

bamfilepath=$1
wfolder=$2

bam=$(basename "$bamfilepath")

#ls $melt_full_path/me_refs/1KGP_Hg19/*.zip > mei_list.txt #need to be done only once

# Run MELT with Single mode
java -jar $melt_full_path/MELT.jar Single \
    -a \
    -b hs37d5/NC_007605 \
    -h $melt_full_path/Demo/hs37d5.fa \
    -bamfile $bamfilepath \
    -n $melt_full_path/add_bed_files/1KGP_Hg19/hg19.genes.bed \
    -t $melt_full_path/mei_list.txt \
    -w $wfolder
