#!/bin/bash

SAMPLE_ID=sample_id.txt
BAMS=NT-pairs.csv
X10_BAM=null
WFOLDER=/hpf/largeprojects/davidm/shilpa/TE_runs/results/xtea_casectrl_runs
OUT_SCRTP=submit_scripts.sh
TIME=0-05:00
REP_LIB=/hpf/largeprojects/davidm/shilpa/TE-tools/xTea-new/xTea/rep_lib_annotation
REF=/hpf/largeprojects/davidm/resources/te_scripts/GRCh38_full_analysis_set_plus_decoy_hla.fa
GENE=/hpf/largeprojects/davidm/resources/gencode.v33.annotation.gff3_hg38
XTEA=/hpf/largeprojects/davidm/shilpa/TE-tools/xTea-new/xTea/
BLK_LIST=/hpf/largeprojects/davidm/shilpa/TE-tools/xTea-new/xTea/rep_lib_annotation/blacklist/hg38/centromere.bed

python3 ${XTEA}bin/xtea \
	 --case_control \
	 --tumor \
         -i ${SAMPLE_ID} \
         -b ${BAMS} \
         -x null \
         -p ${WFOLDER} \
         -o ${OUT_SCRTP} \
         -l ${REP_LIB} \
         -r ${REF} \
         -g ${GENE} \
         --xtea ${XTEA}/xtea \
         -f 5907 \
         -y 15 \
         --blacklist ${BLK_LIST}

