#!/bin/bash

module load samtools/1.9 bwa/0.7.17 anaconda/4.6.14

source /hpf/tools/centos6/anaconda/4.6.14/bin/activate

conda activate TransposableElements

PREFIX=/hpf/largeprojects/davidm/vsubasri/transposable_elements/xTea/test/NA12878/Alu/
############
############
ANNOTATION=/hpf/largeprojects/davidm/vsubasri/transposable_elements/xTea/rep_lib_annotation/Alu/hg19/hg19_Alu.out
ANNOTATION1=/hpf/largeprojects/davidm/vsubasri/transposable_elements/xTea/rep_lib_annotation/Alu/hg19/hg19_Alu.out
REF=/hpf/largeprojects/davidm/resources/Homo_sapiens.GRCh37.75.dna_rm.primary_assembly.fa
GENE=/hpf/largeprojects/davidm/vsubasri/transposable_elements/xTea/gencode.v33lift37.annotation.gff3
BLACK_LIST=null
L1_COPY_WITH_FLANK=/hpf/largeprojects/davidm/vsubasri/transposable_elements/xTea/rep_lib_annotation/Alu/hg19/hg19_AluJabc_copies_with_flank.fa
SF_FLANK=null
L1_CNS=/hpf/largeprojects/davidm/vsubasri/transposable_elements/xTea/rep_lib_annotation/consensus/ALU.fa
XTEA_PATH=/hpf/largeprojects/davidm/vsubasri/transposable_elements/xTea/xtea/
BAM_LIST=${PREFIX}"bam_list.txt"
BAM1=${PREFIX}"10X_phased_possorted_bam.bam"
BARCODE_BAM=${PREFIX}"10X_barcode_indexed.sorted.bam"
TMP=${PREFIX}"tmp/"
TMP_CLIP=${PREFIX}"tmp/clip/"
TMP_CNS=${PREFIX}"tmp/cns/"
TMP_TNSD=${PREFIX}"tmp/transduction/"
############
############
python ${XTEA_PATH}"x_TEA_main.py" -C -i ${BAM_LIST} --lc 3 --rc 3 --cr 1  -r ${L1_COPY_WITH_FLANK}  -a ${ANNOTATION} --cns ${L1_CNS} --ref ${REF} -p ${TMP} -o ${PREFIX}"candidate_list_from_clip.txt"  -n 8 --cp /hpf/largeprojects/davidm/vsubasri/transposable_elements/xTea/test/NA12878/pub_clip/     
python ${XTEA_PATH}"x_TEA_main.py"  -D -i ${PREFIX}"candidate_list_from_clip.txt" --nd 5 --ref ${REF} -a ${ANNOTATION} -b ${BAM_LIST} -p ${TMP} -o ${PREFIX}"candidate_list_from_disc.txt" -n 8    
python ${XTEA_PATH}"x_TEA_main.py" -N --cr 3 --nd 5 -b ${BAM_LIST} -p ${TMP_CNS} --fflank ${SF_FLANK} --flklen 3000 -n 8 -i ${PREFIX}"candidate_list_from_disc.txt" -r ${L1_CNS} --ref ${REF} -a ${ANNOTATION} -o ${PREFIX}"candidate_disc_filtered_cns.txt"    
python ${XTEA_PATH}"x_TEA_main.py" --transduction --cr 3 --nd 5 -b ${BAM_LIST} -p ${TMP_TNSD} --fflank ${SF_FLANK} --flklen 3000 -n 8 -i ${PREFIX}"candidate_disc_filtered_cns.txt" -r ${L1_CNS} --ref ${REF} --input2 ${PREFIX}"candidate_list_from_disc.txt.clip_sites_raw_disc.txt" --rtype 2 -a ${ANNOTATION1}    -o ${PREFIX}"candidate_disc_filtered_cns2.txt"
python ${XTEA_PATH}"x_TEA_main.py" --sibling --cr 3 --nd 5 -b ${BAM_LIST} -p ${TMP_TNSD} --fflank ${SF_FLANK} --flklen 3000 -n 8 -i ${PREFIX}"candidate_disc_filtered_cns2.txt" -r ${L1_CNS} --ref ${REF} --input2 ${PREFIX}"candidate_list_from_disc.txt.clip_sites_raw_disc.txt" --rtype 2 -a ${ANNOTATION1} --blacklist ${BLACK_LIST}    -o ${PREFIX}"candidate_sibling_transduction2.txt"
python ${XTEA_PATH}"x_TEA_main.py" --postF --rtype 2 -p ${TMP_CNS} -n 8 -i ${PREFIX}"candidate_disc_filtered_cns2.txt" -a ${ANNOTATION1}  -o ${PREFIX}"candidate_disc_filtered_cns_post_filtering.txt"
python ${XTEA_PATH}"x_TEA_main.py" --postF --rtype 2 -p ${TMP_CNS} -n 8 -i ${PREFIX}"candidate_disc_filtered_cns2.txt.high_confident" -a ${ANNOTATION1} --blacklist ${BLACK_LIST}  -o ${PREFIX}"candidate_disc_filtered_cns.txt.high_confident.post_filtering.txt"
python ${XTEA_PATH}"x_TEA_main.py" --gene -a ${GENE} -i ${PREFIX}"candidate_disc_filtered_cns.txt.high_confident.post_filtering.txt"  -n 8 -o ${PREFIX}"candidate_disc_filtered_cns.txt.high_confident.post_filtering_with_gene.txt"
python ${XTEA_PATH}"x_TEA_main.py" --gntp_classify -i ${PREFIX}"candidate_disc_filtered_cns.txt.high_confident.post_filtering_with_gene.txt"  -n 1 --model ${XTEA_PATH}"genotyping/trained_model_ssc_py2_random_forest_two_category.pkl"  -o ${PREFIX}"candidate_disc_filtered_cns.txt.high_confident.post_filtering_with_gene_gntp.txt"
python ${XTEA_PATH}"x_TEA_main.py" --gVCF -i ${PREFIX}"candidate_disc_filtered_cns.txt.high_confident.post_filtering_with_gene_gntp.txt"  -o ${PREFIX} -b ${BAM_LIST} --ref ${REF} --rtype 2
