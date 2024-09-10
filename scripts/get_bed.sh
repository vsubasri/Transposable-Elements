#!/bin/bash
te_input=/hpf/largeprojects/davidm/vsubasri/transposable_elements/final_results/kics_te/
if [[ -d $te_input ]]; then
        te_files=( $(find $te_input -type f -name "*.final.txt") )
else
        mapfile -t te_files < te_input
fi

for file in ${te_files[@]}
do
	outdir=$(dirname $file)
	id=$(basename $outdir)
	echo $id
	awk '{ print $2"\t"$3"\t"$4"\t"$7"\t"$1}' $file | sed "s/$/\t"$id"/;1d;s/^/chr/;s/.realigned-recalibrated//g;s/[ \t]\+$//" > $outdir/${id}.bed
done

cat $outdir/*.bed > $outdir/all.bed
