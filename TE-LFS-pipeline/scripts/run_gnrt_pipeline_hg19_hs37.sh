#!/bin/bash

SAMPLE_ID=sample_id.txt
BAMS=bamfile_desc.txt
X10_BAM=null
WFOLDER=$PWD
OUT_SCRTP=submit_scripts.sh
TIME=0-05:00
REP_LIB=/hpf/largeprojects/davidm/shilpa/TE-tools/xTea-new/xTea/rep_lib_annotation
REF=/hpf/largeprojects/davidm/resources/hs37d5.fa
GENE=/hpf/largeprojects/davidm/resources/gencode.v33lift37.annotation.gff3
XTEA=/hpf/largeprojects/davidm/shilpa/TE-tools/xTea-new/xTea/
BLK_LIST=/hpf/largeprojects/davidm/shilpa/TE-tools/xTea-new/xTea/rep_lib_annotation/blacklist/hg19/centromere.bed

python3 ${XTEA}bin/xtea_hg19 \
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
	 -y 7 \
	 --blacklist ${BLK_LIST}

#python3 ${XTEA}xtea/gnrt_pipeline_local.py --case_control --tumor -i ${SAMPLE_ID} -b ${BAMS} -x null -p ${WFOLDER} -o ${OUT_SCRTP} -l ${REP_LIB} -r ${REF} -g ${GENE} --xtea ${XTEA}/xtea -f 5907 -y 7

