#!/bin/bash

module load bedtools/2.27.1 java/1.8.0_91 bcftools/1.9 samtools/1.9

ref_fasta=/hpf/largeprojects/davidm/resources/human_g1k_v37_decoy.fasta
#ref_fasta=/hpf/largeprojects/davidm/resources/Homo_sapiens.GRCh37.75.dna_rm.primary_assembly.fa
gatk3=/hpf/largeprojects/davidm/resources/gatk-3.8/GenomeAnalysisTK.jar
annotate=/hpf/largeprojects/davidm/vsubasri/transposable_elements/vep.sh
bam_input=/hpf/largeprojects/davidm/vsubasri/transposable_elements/xTea/lfs.bam.list
#bam_input=/hpf/largeprojects/davidm/data/kics/archived_bams/test.list
tedir=/hpf/largeprojects/davidm/vsubasri/transposable_elements
breakpoint_overlap=10000
benign_overlap=50
outdir=${tedir}/final_results/lfs_te
#outdir=${tedir}/final_results/kics_te
group=LFS
#group=KiCS

if [[ -d $bam_input ]]; then
        bam_files=( $(find $bam_input -type f -name "*.bam") )
else
        mapfile -t bam_files < $bam_input
fi

if [ ! -d $outdir ] ; then
	mkdir $outdir
fi

for bam_file in ${bam_files[@]}
do
	sample=$(basename $bam_file .realigned-recalibrated.bam)
	melt_dir=${tedir}/MELT/${group}
	xtea_dir=${tedir}/xTea/${group}/${sample}

	if [ ! -d $outdir/$sample ] ; then
		mkdir $outdir/$sample
	fi

#	if [ -f $outdir/$sample/${sample}.SVA.vcf ] ; then

		echo "[ SAMPLE ]: $sample"

		melt_te=("ALU" "LINE1" "SVA")
		xtea_te=("Alu" "LINE1" "SVA")

		for i in "${!xtea_te[@]}"; do
			xtea_vcf=${xtea_dir}/${xtea_te[i]}/${sample}.realigned-recalibrated_${xtea_te[i]^^}.vcf
			if [[ ${group} == "KiCS" ]]; then
				melt_vcf=${melt_dir}/${melt_te[i]}/${sample}.vcf
			else
				melt_vcf=${melt_dir}/${sample}/${melt_te[i]}.final_comp.vcf
			fi

			echo "[ TE ]: ${melt_te[i]}"

			if [[ ${xtea_te[i]} == "LINE1" ]] ; then
				xtea_vcf=${xtea_dir}/L1/${sample}.realigned-recalibrated_${xtea_te[i]^^}.vcf
			fi

			printf '%s\n' $xtea_vcf $melt_vcf > $outdir/$sample/${xtea_te[i]}.input.list

			if [ ! -f $outdir/$sample/${sample}.${xtea_te[i]}_${breakpoint_overlap}.vcf ] ; then
				/hpf/largeprojects/davidm/vsubasri/code/SURVIVOR/Debug/SURVIVOR merge \
					$outdir/$sample/${xtea_te[i]}.input.list ${breakpoint_overlap} 2 1 1 1 0 \
					$outdir/$sample/${sample}.${xtea_te[i]}_${breakpoint_overlap}.vcf 
				/hpf/largeprojects/davidm/vsubasri/code/SURVIVOR/Debug/SURVIVOR genComp $outdir/$sample/${sample}.${xtea_te[i]}_${breakpoint_overlap}.vcf $outdir/$sample/${sample}.${xtea_te[i]}_${breakpoint_overlap}.mat.txt
				/hpf/largeprojects/davidm/vsubasri/code/SURVIVOR/Debug/SURVIVOR stats $outdir/$sample/${sample}.${xtea_te[i]}_${breakpoint_overlap}.vcf -1 -1 -1 $outdir/$sample/${sample}.${xtea_te[i]}_${breakpoint_overlap}.stats.txt
			fi

			if [ ! -f $outdir/$sample/${sample}.${xtea_te[i]}_${breakpoint_overlap}_${benign_overlap}.tsv ] ; then
				$ANNOTSV/bin/AnnotSV -SVinputFile $outdir/$sample/${sample}.${xtea_te[i]}_${breakpoint_overlap}.vcf -SVinputInfo 1 -outputFile $outdir/$sample/${sample}.${xtea_te[i]}_${breakpoint_overlap}_${benign_overlap} -genomeBuild GRCh37 -overlap ${benign_overlap}
			fi

#			input=$outdir/$sample/${sample}.${xtea_te[i]}.vcf;output_vep=$outdir/$sample/${sample}.${xtea_te[i]}.vep.vcf;output_tab=$outdir/$sample/${sample}.${xtea_te[i]}.vep.tab
#			qsub -v input=$input,output_vep=$output_vep,output_tab=$output_tab $annotate 

		done
#	fi
#	cd $outdir/$sample/
#	bgzip -c ${sample}.vcf > ${sample}.vcf.gz
#	tabix -p vcf ${sample}.vcf.gz
#	bcftools +fill-from-fasta ${sample}.vcf.gz -- -c REF -f $ref_fasta > ${sample}_final.vcf

done

