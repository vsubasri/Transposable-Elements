#!/bin/bash

# Transposable Elements Analysis Pipeline GERMLINE SAMPLES

# Usage: ./TE-LFS-germline-pipe.sh /path/to/bamdir (this would be used as $BAM_DIR)

#########################################
# Step 0: Setup and Configuration
#########################################

# Define the directory containing BAM files
BAM_DIR="${1:-/hpf/largeprojects/davidm/data/te_test/germline_bam/}"
# default BAM_DIR is also in case

# output directories
TEPIPE_DIR="/hpf/largeprojects/davidm/shilpa/TE-LFS-pipeline"
OUTPUT_DIR=$PWD/results
BAM_fixed=$PWD/fixed_bams
XTEA_DIR=$OUTPUT_DIR/xtea_runs
MELT_DIR=$OUTPUT_DIR/melt_runs
INS_DIR=$OUTPUT_DIR/ins_runs
SURVIVOR_PATH="/hpf/largeprojects/davidm/shilpa/TE-tools/SURVIVOR-master/Debug/SURVIVOR"

mkdir -p $XTEA_DIR $MELT_DIR $INS_DIR $BAM_fixed $XTEA_DIR/logs $MELT_DIR/logs $INS_DIR/logs $BAM_fixed/logs "logs"

# Create bamfilepaths.txt bamfile_desc.txt and sample_id.txt
ls $BAM_DIR/*.bam > rawbamfilepaths.txt

: <<'COMMENT_BLOCK'
#########################################
# Step 1: File Preprocessing
#########################################
cp rawbamfilepaths.txt $TEPIPE_DIR/scripts/submit_preprocessbam.sh $BAM_fixed
(
cd $BAM_fixed || exit
while read -r bamfilepath; do
    sbatch submit_preprocessbam.sh "$bamfilepath"
done < rawbamfilepaths.txt
)
# Wait for preprocessing to complete
while squeue -u $USER | grep -q 'bam_preprocess'; do
    echo "**waiting for the bam preprocessing steps ..."
    sleep 1000
done
#########################################
COMMENT_BLOCK

ls $BAM_fixed/sorted_fixed*.bam > bamfilepaths.txt
awk -F"/" '{print $NF" "$0}' bamfilepaths.txt >bamfile_desc.txt
awk -F".bam" '{print $1}' bamfile_desc.txt >sample_id.txt

#########################################
# Step 2: xTea
#########################################

cp bamfile_desc.txt sample_id.txt $TEPIPE_DIR/scripts/run_gnrt_pipeline_hg19.sh $TEPIPE_DIR/scripts/interm_set_prep_sbatch.sh $TEPIPE_DIR/scripts/slurm_header_xtea.txt $XTEA_DIR
(
cd $XTEA_DIR || exit
	bash run_gnrt_pipeline_hg19.sh
	#prepating files for slurm job
	bash interm_set_prep_sbatch.sh
	bash submit_scripts.sh
)

#run_gnrt uses bamfile_desc.txt and creates folders corresponding to their sample_id

#########################################
# Step 3: MELT
#########################################

process_melt() {
    local bamfilepath=$1
    local sample_id=$(basename "$bamfilepath" | awk -F"." '{print $1}')
    local sample_output_dir=$MELT_DIR/$sample_id
    mkdir -p $sample_output_dir
    sbatch submit_meltrun.sh "$bamfilepath" "$sample_output_dir"
}
cp $TEPIPE_DIR/scripts/submit_meltrun.sh bamfilepaths.txt $MELT_DIR

(
cd $MELT_DIR || exit
while read -r bamfilepath; do
    process_melt $bamfilepath
done < bamfilepaths.txt
)

echo "Melt runs submitted\n"

##########################################
# Step 4: INSurVeyor
##########################################

process_ins() {
    local bamfilepath=$1
    local sample_id=$(basename "$bamfilepath" | awk -F"." '{print $1}')
    local sample_output_dir=$INS_DIR/$sample_id
    mkdir -p $sample_output_dir
    sbatch submit_singu_run_insurveyor.sh "$bamfilepath" "$sample_output_dir"
}
cp $TEPIPE_DIR/scripts/submit_singu_run_insurveyor.sh bamfilepaths.txt $INS_DIR

(
cd $INS_DIR || exit
while read -r bamfilepath; do
	process_ins $bamfilepath
done < bamfilepaths.txt

)

echo "INSurVeyor runs submitted\n"

#########################################
# Step 5: SURVIVOR
#########################################

#wait for TE-runs to complete
#while squeue -u $USER | grep -q 'xtea_runs | melt_runs| ins_runs'; do
#    sleep 100
#done

# Merge VCF files with SURVIVOR

merge_vcfs() {
	local bamfilepath=$1
        local sample_id=$(basename "$bamfilepath" | awk -F"." '{print $1}')
	local alu_list=${OUTPUT_DIR}/${sample_id}_alu_vcfs.list
	local line1_list=${OUTPUT_DIR}/${sample_id}_l1_vcfs.list
	local sva_list=${OUTPUT_DIR}/${sample_id}_sva_vcfs.list

        ls ${XTEA_DIR}/${sample_id}/Alu/sorted_${sample_id}_ALU.vcf ${MELT_DIR}/${sample_id}/ALU.final_comp.vcf ${INS_DIR}/${sample_id}/out.pass.vcf.gz > $alu_list
	ls ${XTEA_DIR}/${sample_id}/L1/sorted_${sample_id}_LINE1.vcf ${MELT_DIR}/${sample_id}/LINE1.final_comp.vcf ${INS_DIR}/${sample_id}/out.pass.vcf.gz > $line1_list
        ls ${XTEA_DIR}/${sample_id}/SVA/sorted_${sample_id}_SVA.vcf ${MELT_DIR}/${sample_id}/SVA.final_comp.vcf ${INS_DIR}/${sample_id}/out.pass.vcf.gz > $sva_list
        #need to confirm the use of Insurveyor outfile out.pass.vcf.gz that those can be added to all the files above
	echo "Running: $survivor_cmd"
        ${SURVIVOR_PATH} merge ${alu_list} 100 2 1 0 0 30 $OUTPUT_DIR/merged_ALU_${sample_id}.vcf
        ${SURVIVOR_PATH} merge ${line1_list} 100 2 1 0 0 30 $OUTPUT_DIR/merged_LINE_${sample_id}.vcf
        ${SURVIVOR_PATH} merge ${sva_list} 100 2 1 0 0 30 $OUTPUT_DIR/merged_SVA_${sample_id}.vcf

}
echo "Running: SURVIVOR\n"

while read -r bamfilepath; do
        merge_vcfs $bamfilepath
done < bamfilepaths.txt

#########################################
# Step 7: Annotations
#########################################

# Load necessary modules
module load AnnotSV bcftool bedtools bedtools

# Annotate merged VCF file
#AnnotSV -SVinputFile sample_merged.vcf
echo "Running AnnotSV\n"

for $sample_id in $(cat sample_id.txt)
do
AnnotSV -SVinputFile $OUTPUT_DIR/merged_ALU_${sample_id}.vcf -outputFile annoSVoutput${sample_id}.tsv -genomBuild GRCh37 -overlap 70 -REreport -includeCI -hpo -candidateGenesFiltering -candidateGenesFile /hpf/largeprojects/davidm/data/te_test/kics_cpg.txt
done

# place /hpf/largeprojects/davidm/data/te_test/gnomad.v4.1.sv.sites_te_hg37_final.bed into SVincludedInFt not FtIncludedInSV in the AnnotSV directory for the custom annotations
#########################################

