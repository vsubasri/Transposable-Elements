#!/bin/bash

#SBATCH -J bam_preprocess
#SBATCH -o logs/%x_output_%j.log
#SBATCH -e logs/%x_error_%j.log
#SBATCH -N 1 -c 16
#SBATCH --mem=32G
#SBATCH -t 24:00:00

bamfilepath=$1

if [ ! -f "$bamfilepath" ]; then
    echo "BAM file does not exist: $bamfilepath"
    exit 1
fi

module load java/1.8.0_91 bowtie2 samtools

output_dir="$PWD/fixed_bams"
mkdir -p "$output_dir"

bam=$(basename "$bamfilepath")

# Define the sorted BAM file names
sorted_bam="sorted_${bam}"
sorted_fixed="sorted_fixed_${bam}"

# checking if sorted BAM file already exists
if [ ! -f "$sorted_bam" ]; then
    echo "Sorting BAM file: $bamfilepath"
    samtools sort -@ 16 -o "${output_dir}/${sorted_bam}" "$bamfilepath"

    # Index sorted BAM file (index will be saved in the same directory as sorted BAM file)
    samtools index -@ 16 "${output_dir}/${sorted_bam}"

    # Run Picard to fix mate information
    java -jar -Xmx32G /hpf/largeprojects/davidm/shilpa/TE-tools/picard.jar FixMateInformation \
             I="${output_dir}/${sorted_bam}" \
	     O="${output_dir}/${sorted_fixed}" \
	     ADD_MATE_CIGAR=true

else
    echo "Sorted BAM file already exists: $sorted_bam"
fi

