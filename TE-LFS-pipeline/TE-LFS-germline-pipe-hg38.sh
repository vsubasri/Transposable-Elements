#!/bin/bash

# Transposable Elements Analysis Pipeline GERMLINE SAMPLES

# Usage: ./TE-LFS-germline-pipe.sh /path/to/bamdir (this would be used as $BAM_DIR)

#########################################
# Step 0: Setup and Configuration
#########################################

# Define the directory containing BAM files
BAM_DIR="${1:-/hpf/largeprojects/davidm/data/te_test/germline_bam_hg38/}"
# default BAM_DIR is mentioned in case (this this currently does not have hg38 data)

# output directories
TEPIPE_DIR="/hpf/largeprojects/davidm/shilpa/TE-LFS-pipeline"
OUTPUT_DIR=$PWD/results_hg38
BAM_fixed=$PWD/fixed_bams_hg38
XTEA_DIR=$OUTPUT_DIR/xtea_runs
MELT_DIR=$OUTPUT_DIR/melt_runs
INS_DIR=$OUTPUT_DIR/ins_runs
ANALY_DIR=$OUTPUT_DIR/analysis
SURVIVOR_PATH="/hpf/largeprojects/davidm/shilpa/TE-tools/SURVIVOR-master/Debug/SURVIVOR"

mkdir -p $XTEA_DIR $MELT_DIR $INS_DIR $BAM_fixed $XTEA_DIR/logs $MELT_DIR/logs $INS_DIR/logs $BAM_fixed/logs $INS_DIR/logs

# Create bamfilepaths.txt bamfile_desc.txt and sample_id.txt
ls $BAM_DIR/*.bam > rawbamfilepaths.txt

#########################################
# Step 1: File Preprocessing
#########################################
cp rawbamfilepaths.txt $TEPIPE_DIR/scripts/submit_preprocessbam_hg38.sh $TEPIPE_DIR/scripts/process_metrics.sh $BAM_fixed
(
cd $BAM_fixed || exit
while read -r bamfilepath; do
    sbatch submit_preprocessbam_hg38.sh "$bamfilepath"
done < rawbamfilepaths.txt
)
# Wait for preprocessing to complete
#while squeue -u $USER | grep -q 'bam_preprocess'; do
#    echo "**waiting for the bam preprocessing steps ..."
#    sleep 1000
#done
#########################################

ls $BAM_fixed/sorted_fixed_*_N.bam > bamfilepaths.txt
awk -F"/" '{print $NF" "$0}' bamfilepaths.txt | sed 's/.bam//' >bamfile_desc.txt
awk '{print $1}' bamfile_desc.txt >sample_id.txt

#########################################
# Step 2: xTea
#########################################
#
cp bamfile_desc.txt sample_id.txt $TEPIPE_DIR/scripts/run_gnrt_pipeline_hg38.sh $TEPIPE_DIR/scripts/interm_set_prep_sbatch.sh $TEPIPE_DIR/scripts/slurm_header_xtea.txt $XTEA_DIR
(
cd $XTEA_DIR || exit
	bash run_gnrt_pipeline_hg38.sh
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
    sbatch submit_meltrun_hg38.sh "$bamfilepath" "$sample_output_dir"
}
cp $TEPIPE_DIR/scripts/submit_meltrun_hg38.sh bamfilepaths.txt $MELT_DIR

(
cd $MELT_DIR || exit
while read -r bamfilepath; do
    process_melt $bamfilepath
done < bamfilepaths.txt
)

echo "Melt runs submitted"

##########################################
# Step 4: INSurVeyor
##########################################

process_ins() {
    local bamfilepath=$1
    local sample_id=$(basename "$bamfilepath" | awk -F"." '{print $1}')
    local sample_output_dir=$INS_DIR/$sample_id
    mkdir -p $sample_output_dir
    sbatch submit_singu_run_insurveyor_hg38.sh "$bamfilepath" "$sample_output_dir"
}
cp $TEPIPE_DIR/scripts/submit_singu_run_insurveyor_hg38.sh bamfilepaths.txt $INS_DIR

(
cd $INS_DIR || exit
while read -r bamfilepath; do
	process_ins $bamfilepath
done < bamfilepaths.txt

)

echo "INSurVeyor runs submitted"

#########################################
# Step 5: SURVIVOR
# Merge VCF files with SURVIVOR
#########################################

#wait for TE-runs to complete
#while squeue -u $USER | grep -q 'xtea_runs | melt_runs| ins_runs'; do
#    sleep 100
#done

cd $INS_DIR
for sample_id in $(cat sample_id.txt ); do cd $sample_id; gunzip out.pass.vcf.gz; cd ..; done
cd ..

merge_vcfs() {
        local bamfilepath=$1
        local sample_id=$(basename "$bamfilepath" | awk -F"." '{print $1}')
        local alu_list=${OUTPUT_DIR}/${sample_id}_alu_vcfs.list
        local line1_list=${OUTPUT_DIR}/${sample_id}_l1_vcfs.list
        local sva_list=${OUTPUT_DIR}/${sample_id}_sva_vcfs.list

        ls ${XTEA_DIR}/${sample_id}/Alu/${sample_id}_ALU.vcf ${MELT_DIR}/${sample_id}/ALU.final_comp.vcf ${INS_DIR}/${sample_id}/out.pass.vcf > $alu_list
        ls ${XTEA_DIR}/${sample_id}/L1/${sample_id}_LINE1.vcf ${MELT_DIR}/${sample_id}/LINE1.final_comp.vcf ${INS_DIR}/${sample_id}/out.pass.vcf > $line1_list
        ls ${XTEA_DIR}/${sample_id}/SVA/${sample_id}_SVA.vcf ${MELT_DIR}/${sample_id}/SVA.final_comp.vcf ${INS_DIR}/${sample_id}/out.pass.vcf > $sva_list
        #need to confirm the use of Insurveyor outfile out.pass.vcf.gz that those can be added to all the files above
        echo "Running: $survivor_cmd"
        ${SURVIVOR_PATH} merge ${alu_list} 100 2 1 1 0 0 $ANALY_DIR/merged_ALU_${sample_id}.vcf
        ${SURVIVOR_PATH} merge ${line1_list} 100 2 1 1 0 0 $ANALY_DIR/merged_LINE_${sample_id}.vcf
        ${SURVIVOR_PATH} merge ${sva_list} 100 2 1 1 0 0 $ANALY_DIR/merged_SVA_${sample_id}.vcf
        ${SURVIVOR_PATH} genComp merged_ALU_${sample_id}.vcf 0 merged_ALU_${sample_id}.mat.txt
        ${SURVIVOR_PATH} genComp merged_LINE_${sample_id}.vcf 0 merged_LINE_${sample_id}.mat.txt
        ${SURVIVOR_PATH} genComp merged_SVA_${sample_id}.vcf 0 merged_SVA_${sample_id}.mat.txt
}
echo "Running: SURVIVOR"
cp bamfilepaths.txt $ANALY_DIR
(
cd $ANALY_DIR || exit
while read -r bamfilepath; do
        merge_vcfs $bamfilepath
done < bamfilepaths.txt
)

#########################################
# Step 6: Annotations
# Annotate merged VCF file
#########################################

# Load necessary modules
module load bcftools bedtools
source /hpf/largeprojects/davidm/shilpa/TE-tools/annotsv_env/setup_annotsv_env.sh

echo "Running AnnotSV"

cd $ANALY_DIR
ls *.vcf >allvcf.list

for id in $(cat allvcf.list)
do
        echo "$id"
        temp=$(basename "$id" .vcf)
        AnnotSV -SVinputFile $id -outputFile "$temp/SVoutput_$temp.tsv" -genomeBuild GRCh37 -overlap 70 -REreport 1 -includeCI 1 -candidateGenesFiltering 1 -candidateGenesFile /hpf/largeprojects/davidm/data/te_test/kics_cpg.txt
done
cd ..

# debug pending: CUSTOM Annotations added: /hpf/largeprojects/davidm/data/te_test/gnomad.v4.1.sv.sites_te_hg37_final.bed file is placed under /hpf/largeprojects/davidm/shilpa/TE-tools/annotsv_env/share/AnnotSV/Annotations_Human/Users/GRCh37/SVincludedInFt 
#########################################

