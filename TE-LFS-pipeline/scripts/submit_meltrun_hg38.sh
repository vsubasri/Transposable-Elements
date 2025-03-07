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
#reference_fasta="/hpf/largeprojects/davidm/resources/hs37d5.fa"
#reference_fasta="/hpf/largeprojects/davidm/resources/GRCh38.primary_assembly.genome.fa"
#reference_fasta="/hpf/largeprojects/davidm/resources/HG38/hg38.fa"
reference_fasta="/hpf/largeprojects/davidm/resources/te_scripts/GRCh38_full_analysis_set_plus_decoy_hla.fa"

bam=$(basename "$bamfilepath")

ls $melt_full_path/me_refs/Hg38/*.zip > $melt_full_path/mei_list_hg38.txt #need to be done only once

# Run MELT with Single mode
java -jar $melt_full_path/MELT.jar Single \
    -ac \
    -h $reference_fasta \
    -bamfile $bamfilepath \
    -n $melt_full_path/add_bed_files/Hg38/Hg38.genes.bed \
    -t $melt_full_path/mei_list_hg38.txt \
    -w $wfolder

