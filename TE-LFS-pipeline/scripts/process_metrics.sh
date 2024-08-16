#!/bin/bash

# Define the input and output files
input_file=$1
output_file=$2

# Initialize variables
mean_coverage=""
sd_coverage=""
median_coverage=""
total_reads=""
mean_read_length=""
sd_read_length=""
median_read_length=""
pct_chimeras=""
average_quality=""
sample=""

# Process the input file once
while IFS= read -r line; do
  # Extract sample name
  if [[ $line == "# CollectWgsMetrics --INPUT sorted_fixed_"*".bam"* ]]; then
    sample=$(echo $line | sed -n 's/.*--INPUT sorted_fixed_\([^ ]*\).bam.*/\1/p')
  fi

  # Extract Coverage Metrics
  if [[ $line == GENOME_TERRITORY* ]]; then
    read -r line
    mean_coverage=$(echo $line | awk '{print $2}')
    sd_coverage=$(echo $line | awk '{print $3}')
    median_coverage=$(echo $line | awk '{print $4}')
  fi
  
  # Extract Alignment Summary for PAIR
  if [[ $line == PAIR* ]]; then
    total_reads=$(echo $line | awk '{print $2}')
    mean_read_length=$(echo $line | awk '{print $16}')
    sd_read_length=$(echo $line | awk '{print $17}')
    median_read_length=$(echo $line | awk '{print $18}')
    pct_chimeras=$(echo $line | awk '{print $29}')
  fi
  
  # Extract Base Quality Stat
  if [[ $line == *"average quality"* ]]; then
    average_quality=$(echo $line | awk '{print $4}')
  fi

done < "$input_file"

# Check if the output file exists; if not, write the header
if [ ! -f "$output_file" ]; then
  echo -e "Sample\tMean Coverage\tSD Coverage\tMedian Coverage\tTotal Reads\tMean Read Length\tSD Read Length\tMedian Read Length\tPCT Chimeras\tAverage Quality" > "$output_file"
fi

# Append the extracted values to the output file
echo -e "$sample\t$mean_coverage\t$sd_coverage\t$median_coverage\t$total_reads\t$mean_read_length\t$sd_read_length\t$median_read_length\t$pct_chimeras\t$average_quality" >> "$output_file"
