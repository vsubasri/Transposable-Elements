#!/bin/bash

#SBATCH -J bam_preprocess
#SBATCH -o logs/%x_output_%j.log
#SBATCH -e logs/%x_error_%j.log
#SBATCH -N 1 -c 16
#SBATCH --mem=64G
#SBATCH -t 24:00:00

bamfilepath=$1

if [ ! -f "$bamfilepath" ]; then
    echo "BAM file does not exist: $bamfilepath"
    exit 1
fi

module load java/1.8.0_91 bowtie2 samtools
picardpath="/hpf/largeprojects/davidm/shilpa/TE-tools/picard.jar"
#reference_fasta="/hpf/largeprojects/davidm/resources/GRCh38.primary_assembly.genome.fa"
reference_fasta="/hpf/largeprojects/davidm/resources/te_scripts/GRCh38_full_analysis_set_plus_decoy_hla.fa"

output_dir=$PWD
#mkdir -p "$output_dir"

bam=$(basename "$bamfilepath")

# Define the sorted BAM file names
sorted_bam="sorted_${bam}"
sorted_fixed="sorted_fixed_${bam}"

# checking if sorted BAM file already exists
if [ ! -f "$sorted_bam" ] || [ -f "$sorted_fixed" ]; then
    echo "Sorting BAM file: $bamfilepath"
    samtools sort -@ 16 -o "${output_dir}/${sorted_bam}" "$bamfilepath"

    # Index sorted BAM file (index will be saved in the same directory as sorted BAM file)
    samtools index -@ 16 "${output_dir}/${sorted_bam}"

    # Run Picard to fix mate information
    java -jar -Xmx32G $picardpath FixMateInformation \
             I="${output_dir}/${sorted_bam}" \
	     O="${output_dir}/${sorted_fixed}" \
	     ADD_MATE_CIGAR=true

    #Index sorted_fix BAM file (index should be saved in the same directory)
    samtools index -@ 16 "${output_dir}/${sorted_fixed}"

else
    echo "Sorted BAM file already exists: $sorted_bam"
fi

echo "calculating coverage and summary metrics calculations"
# coverage metrics
java -jar -Xmx32G $picardpath CollectWgsMetrics I="$sorted_fixed" O="${output_dir}/${bam}_wgs_metrics.txt" R="$reference_fasta"

# alignment summary metrics
java -jar -Xmx32G $picardpath CollectAlignmentSummaryMetrics I="$sorted_fixed" O="${output_dir}/${bam}_alignment_metrics.txt" R="$reference_fasta"

# Calculate base quality statistics using samtools stats
samtools stats -@ 16 "${sorted_fixed}" > "${output_dir}/${bam}_bamfilestats.txt"

combined_metrics="${output_dir}/${bam}_combined_metrics.txt"
{
echo "Coverage Metrics:"
cat "${output_dir}/${bam}_wgs_metrics.txt"
echo ""
echo "Alignment Summary:"
cat "${output_dir}/${bam}_alignment_metrics.txt"
echo ""
echo "Base quality Stat:"
cat "${output_dir}/${bam}_bamfilestats.txt"
} > $combined_metrics

processed_combined_metrics="${output_dir}/${bam}_processed_combined_metrics.txt"
bash process_metrics.sh $combined_metrics $processed_combined_metrics 

#rm "${output_dir}/${bam}_wgs_metrics.txt" "${output_dir}/${bam}_alignment_metrics.txt" "${output_dir}/${bam}_bamfilestats.txt"

echo "Preprocessing and metrics calculation completed for: $bamfilepath"

