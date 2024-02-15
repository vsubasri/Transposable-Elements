#!/bin/bash

outdir=/hpf/largeprojects/davidm/vsubasri/transposable_elements/xTea/KiCS
header=/hpf/largeprojects/davidm/vsubasri/transposable_elements/xTea/header
bamlist=/hpf/largeprojects/davidm/data/kics/archived_bams/test.list

python parse_xTea_input.py -b $bamlist -o $outdir

/hpf/largeprojects/davidm/vsubasri/transposable_elements/xTea/bin/xtea_hg19 \
	-i $outdir/sample_id.txt \
	-b $outdir/illumina_bam_list.txt \
	-x null \
	-p $outdir/ \
	-o $outdir/submit_jobs.sh \
	-l /hpf/largeprojects/davidm/vsubasri/transposable_elements/xTea/rep_lib_annotation/ \
	-r /hpf/largeprojects/davidm/resources/Homo_sapiens.GRCh37.75.dna_rm.primary_assembly.fa \
	-g /hpf/largeprojects/davidm/vsubasri/transposable_elements/xTea/gencode.v33lift37.annotation.gff3 \
	--xtea /hpf/largeprojects/davidm/vsubasri/transposable_elements/xTea/xtea/ \
	-f 5907 -y 15

for file in ${outdir}/*/*/run_xTEA_pipeline.sh
do
	cat $header $file > tmp
	mv tmp $file
done

