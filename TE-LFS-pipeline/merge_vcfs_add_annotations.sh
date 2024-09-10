#!/bin/bash

# Transposable Elements Analysis Pipeline for GERMLINE SAMPLES

TEPIPE_DIR="/hpf/largeprojects/davidm/shilpa/TE-LFS-pipeline" # the directory where pipeline is located

base_dir=$PWD # this is where all the output files will be created

OUTPUT_DIR=$PWD/results
BAM_fixed=$PWD/fixed_bams # if any edits are needed in sample_id.txt, one can go to this folder
XTEA_DIR=$OUTPUT_DIR/xtea_runs
MELT_DIR=$OUTPUT_DIR/melt_runs
INS_DIR=$OUTPUT_DIR/ins_runs
ANALY_DIR=$OUTPUT_DIR/analysis

mkdir -p $ANALY_DIR

echo "Analysis will be carried out at $ANALY_DIR"

module load java/11

JASMINE_PATH="/hpf/largeprojects/davidm/resources/jasmine/Jasmine"

#SURVIVOR_PATH="/hpf/largeprojects/davidm/shilpa/TE-tools/SURVIVOR-master/Debug/SURVIVOR"

#########################################
# Step 5: SURVIVOR
# Merge VCF files with JASMINE
#########################################

cp $base_dir/sample_id.txt $TEPIPE_DIR/scripts/clean*.py $XTEA_DIR
cp $base_dir/sample_id.txt $TEPIPE_DIR/scripts/clean*.py $MELT_DIR
cp $base_dir/sample_id.txt $TEPIPE_DIR/scripts/clean*.py $INS_DIR
cp $base_dir/sample_id.txt $TEPIPE_DIR/scripts/clean*.py $ANALY_DIR

echo "Cleaning VCFs ..."

####### Cleaning vcfs ###################
cd $INS_DIR
   for sample_id in $(cat sample_id.txt ); do cd $sample_id; gunzip out.pass.vcf.gz; cd ..; done
cd ../..

cd $XTEA_DIR
   for sample_id in $(cat sample_id.txt)
   do
   cp clean_xtea.py $sample_id
   cd $sample_id
   python clean_xtea.py Alu/${sample_id}_ALU.vcf $sample_id > Alu/${sample_id}_ALU_edited.vcf
   python clean_xtea.py L1/${sample_id}_LINE1.vcf $sample_id > L1/${sample_id}_LINE1_edited.vcf
   python clean_xtea.py SVA/${sample_id}_SVA.vcf $sample_id > SVA/${sample_id}_SVA_edited.vcf
   cd ..
   done
cd ../..

cd $MELT_DIR
	for sample_id in $(cat sample_id.txt)
        do
	cp clean_melt.py $sample_id
        cd $sample_id
        python clean_melt.py ALU.final_comp.vcf $sample_id > ALU.final_comp_edited.vcf
	python clean_melt.py LINE1.final_comp.vcf $sample_id > LINE1.final_comp_edited.vcf
	python clean_melt.py SVA.final_comp.vcf $sample_id > SVA.final_comp_edited.vcf
	cd ..
done
cd ../..

cd $INS_DIR
	for sample_id in $(cat sample_id.txt)
        do
	cp clean_ins.py $sample_id
        cd $sample_id
        python clean_ins.py out.pass.vcf $sample_id > out.pass_edited.vcf
	cd ..
done
cd ../..

#########################################

merge_vcfs() {
        ##variables
        #local bamfilepath=$1
        #local sample_id=$(basename "$bamfilepath" | awk -F"." '{print $1}')
        local sample_id=$1
        local alu_list=${sample_id}_alu_vcfs.list
        local line1_list=${sample_id}_l1_vcfs.list
        local sva_list=${sample_id}_sva_vcfs.list

        ls ${INS_DIR}/${sample_id}/out.pass_edited.vcf ${XTEA_DIR}/${sample_id}/Alu/${sample_id}_ALU_edited.vcf ${MELT_DIR}/${sample_id}/ALU.final_comp_edited.vcf > $alu_list #INS DIR first in the list needed for Jasmine
        ls ${INS_DIR}/${sample_id}/out.pass_edited.vcf ${XTEA_DIR}/${sample_id}/L1/${sample_id}_LINE1_edited.vcf ${MELT_DIR}/${sample_id}/LINE1.final_comp_edited.vcf  > $line1_list
        ls ${INS_DIR}/${sample_id}/out.pass_edited.vcf ${XTEA_DIR}/${sample_id}/SVA/${sample_id}_SVA_edited.vcf ${MELT_DIR}/${sample_id}/SVA.final_comp_edited.vcf  > $sva_list

        echo "Running: $JASMINE_PATH"
        
        ${JASMINE_PATH}/jasmine -max_dist=100 -min_support=2 --ignore_type --allow_intrasample file_list=${alu_list} out_file=merged_ALU_${sample_id}.vcf
        ${JASMINE_PATH}/jasmine -max_dist=100 -min_support=2 --ignore_type --allow_intrasample file_list=${line1_list} out_file=merged_LINE_${sample_id}.vcf
        ${JASMINE_PATH}/jasmine -max_dist=100 -min_support=2 --ignore_type --allow_intrasample file_list=${sva_list} out_file=merged_SVA_${sample_id}.vcf

	python clean_jasmine_caller.py merged_ALU_${sample_id}.vcf ${sample_id} > merged_ALU_edited_${sample_id}.vcf
	python clean_jasmine_caller.py merged_LINE_${sample_id}.vcf ${sample_id} > merged_LINE_edited_${sample_id}.vcf
	python clean_jasmine_caller.py merged_SVA_${sample_id}.vcf ${sample_id} > merged_SVA_edited_${sample_id}.vcf
        
	cat merged_ALU_edited_${sample_id}.vcf > concatenated_${sample_id}.vcf
	grep -v '^#' merged_LINE_edited_${sample_id}.vcf >> concatenated_${sample_id}.vcf
	grep -v '^#' merged_SVA_edited_${sample_id}.vcf >> concatenated_${sample_id}.vcf
        #cat merged_ALU_edited_${sample_id}.vcf merged_LINE_edited_${sample_id}.vcf merged_SVA_edited_${sample_id}.vcf > concatenated_${sample_id}.vcf
}

(
cd $ANALY_DIR || exit
while read -r sample_id; do
        merge_vcfs $sample_id
done < sample_id.txt
)

echo "JASMINE RUN completed ..."
#########################################
# Step 6: Annotations
# Annotate merged VCF file
#########################################

# Load necessary modules
module load bcftools bedtools
#source /hpf/largeprojects/davidm/shilpa/TE-tools/annotsv_env/setup_annotsv_env.sh
ANNOSV_PATH="/hpf/largeprojects/davidm/shilpa/TE-tools/annotsv_env/AnnotSV/bin/"

cd $ANALY_DIR
#ls merged_*edited*.vcf > all_merged_vcf.list
ls concatenated*.vcf > all_concat_vcf.list

for id in $(cat all_concat_vcf.list)
do
        echo "$id"
	echo "**************"
        temp=$(basename "$id" .vcf)
        temp=${temp#concatenated_}
        $ANNOSV_PATH/AnnotSV -SVinputFile $id -outputFile "$temp/SVoutput_$temp.tsv" -genomeBuild GRCh37 -overlap 70 -REreport 1 -includeCI 1 -candidateGenesFiltering 1 -candidateGenesFile /hpf/largeprojects/davidm/data/te_test/kics_cpg.txt
done
cd ../..

echo "AnnotSV RUN completed..."

##debug pending: CUSTOM Annotations added: /hpf/largeprojects/davidm/data/te_test/gnomad.v4.1.sv.sites_te_hg37_final.bed file is placed under /hpf/largeprojects/davidm/shilpa/TE-tools/annotsv_env/share/AnnotSV/Annotations_Human/Users/GRCh37/SVincludedInFt 
###########################################
