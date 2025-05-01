#!/bin/bash

#SBATCH -J bam_stat
#SBATCH -o logs/%x_output_%j.log
#SBATCH -e logs/%x_error_%j.log
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 47:00:00

module load java/1.8.0_91 bowtie2 samtools

picardpath="/hpf/largeprojects/davidm/shilpa/TE-tools/picard.jar"
reference_fasta="/hpf/largeprojects/davidm/resources/hs37d5.fa"

sorted_fixed=$1
bam=$(basename "$sorted_fixed" | sed 's/^sorted_fixed_//')

mkdir -p bamstat

# Coverage metrics
java -jar -Xmx64G $picardpath CollectWgsMetrics \
    I="$sorted_fixed" \
    O="bamstat/${bam}_wgs_metrics.txt" \
    R="$reference_fasta"

# Alignment summary metrics
java -jar -Xmx64G $picardpath CollectAlignmentSummaryMetrics \
    I="$sorted_fixed" \
    O="bamstat/${bam}_alignment_metrics.txt" \
    R="$reference_fasta"

# Calculate base quality statistics using samtools stats
samtools stats "$sorted_fixed" > "bamstat/${bam}_bamfilestats.txt"

# Combine metrics into one file
combined_metrics="bamstat/${bam}_combined_metrics.txt"
{
    echo "Coverage Metrics:"
    cat "bamstat/${bam}_wgs_metrics.txt"
    echo ""
    echo "Alignment Summary:"
    cat "bamstat/${bam}_alignment_metrics.txt"
    echo ""
    echo "Base quality Stat:"
    cat "bamstat/${bam}_bamfilestats.txt"
} > "$combined_metrics"

# Process combined metrics
processed_combined_metrics="bamstat/${bam}_processed_combined_metrics.txt"
cp process_metrics.sh bamstat/
bash bamstat/process_metrics.sh "$combined_metrics" "$processed_combined_metrics"

# Clean up intermediate files
rm "bamstat/${bam}_wgs_metrics.txt" "bamstat/${bam}_alignment_metrics.txt" "bamstat/${bam}_bamfilestats.txt"

# Completion message
echo "Preprocessing and metrics calculation completed for: $sorted_fixed"
