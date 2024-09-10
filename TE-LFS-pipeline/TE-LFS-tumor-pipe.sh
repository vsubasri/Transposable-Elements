#!/bin/bash

# Transposable Elements Analysis Pipeline TUMOR SAMPLES

# Usage: ./TE-LFS-tumor-pipe.sh /path/to/bamdir (this would be used as $BAM_DIR)

#########################################
# Step 0: Setup and Configuration
#########################################

# Define the directory containing BAM files
BAM_DIR="${1:-/hpf/largeprojects/davidm/data/te_test/tumor_bam/}"
# default BAM_DIR

# output directories
TEPIPE_DIR="/hpf/largeprojects/davidm/shilpa/TE-LFS-pipeline"
OUTPUT_DIR=$PWD/results
BAM_fixed=$PWD/fixed_bams
XTEA_DIR=$OUTPUT_DIR/xtea_casectrl_runs
TRAF_DIR=$OUTPUT_DIR/trafic_runs
SURVIVOR_PATH="/hpf/largeprojects/davidm/shilpa/TE-tools/SURVIVOR-master/Debug/SURVIVOR"

mkdir -p $XTEA_DIR $BAM_fixed $XTEA_DIR/logs $BAM_fixed/logs $TRAF_DIR $TRAF_DIR/logs

# Create bamfilepaths.txt bamfile_desc.txt and sample_id.txt
ls $BAM_DIR/*.bam > rawbamfilepaths.txt

#########################################
# Step 1: File Preprocessing
#########################################

cp rawbamfilepaths.txt $TEPIPE_DIR/scripts/submit_preprocessbam.sh $TEPIPE_DIR/scripts/process_metrics.sh $BAM_fixed
(
cd $BAM_fixed || exit
while read -r bamfilepath; do
    sbatch submit_preprocessbam.sh "$bamfilepath"
done < rawbamfilepaths.txt
)
# Wait for preprocessing to complete
#while squeue -u $USER | grep -q 'bam_prep'; do
#    echo "**waiting for the bam preprocessing steps ..."
#    sleep 1000
#done

#########################################

ls $BAM_fixed/sorted_fixed_*_T.bam > bamfilepathsNT.txt
awk -F"/" '{print $NF" "$0}' bamfilepathsNT.txt | sed 's/_T.bam//' >bamfile_desc.txt
awk '{print $1}' bamfile_desc.txt >sample_id.txt

#########################################
# Step 2: xTea
#########################################

cp bamfilepathsNT.txt $TEPIPE_DIR/data/NT-pairs.csv sample_id.txt $TEPIPE_DIR/scripts/run_gnrt_pipeline_case_ctrl_hg19.sh $TEPIPE_DIR/scripts/interm_set_prep_sbatch.sh $TEPIPE_DIR/scripts/slurm_header_xtea.txt $XTEA_DIR
(
cd $XTEA_DIR || exit
	sed -i '1d' NT-pairs.csv # This step should be done only once, run_gnrt script from xtea needs this file without header
	bash run_gnrt_pipeline_case_ctrl_hg19.sh
	#prepating files for slurm job
	bash interm_set_prep_sbatch.sh
	bash submit_scripts.sh
)

#run_gnrt uses bamfile_desc.txt and creates folders corresponding to their sample_id

#########################################
## MELT and Insurveyor don't have an option to tumor samples
#########################################
# Step 3: TraFiC
#########################################
cp bamfilepathsNT.txt $TEPIPE_DIR/data/NT-pairs.csv sample_id.txt $TEPIPE_DIR/scripts/make-config-trafic.py $TEPIPE_DIR/scripts/submit_singu_run_trafic.sh $TRAF_DIR

process_trafic() {
    local sample_id=$1
    local sample_output_dir=$TRAF_DIR/$sample_id
    mkdir -p $sample_output_dir
    sbatch submit_singu_run_trafic.sh "config_trafic_$sample_id.yaml" "$sample_output_dir"
}
module load python/3.10.9

(
cd $TRAF_DIR || exit
python3 make-config-trafic.py ##the config files should be created at fixed_bam because that folder will be used to --bind for TraFic container
while read -r sampleid; do
    process_trafic $sampleid
done < sample_id.txt
)

#########################################
# Step 4: SURVIVOR
#########################################
# Step 5: Annotations
#########################################
#Merging and annotations currently under discussion

