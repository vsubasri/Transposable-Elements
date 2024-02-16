#!/usr/bin/env python
import subprocess
import pandas as pd
import numpy as np
import glob
import argparse
import re
from functools import reduce
from collections import Counter
import os
import csv

def read_largefile(fname, chunk_size=10000):
    reader = pd.read_csv(fname, header=0, delimiter='\t',low_memory=False, iterator=True)
    chunks = []
    loop = True
    while loop:
        try:
            chunk = reader.get_chunk(chunk_size)
            chunks.append(chunk)
        except StopIteration:
            loop = False
            print("Iteration is stopped")
    df_ac = pd.concat(chunks, ignore_index=True)
    return(df_ac)

def combine(tedir, suffix):
	pd.options.mode.chained_assignment = None
	files = glob.glob(os.path.join(tedir, "*/*"+suffix+".final.txt"))	
	te_all = []
	for file in files:
		te = pd.read_csv(file, delimiter='\t')
		
		te.columns.values[14:16] = ["xTea","MELT"]
		te['Samples_ID'] = [samples[0] for samples in te['Samples_ID'].str.split(',')]
		#te = te[te['Annotation_mode'] == "full"]
		### Filter MEIs @@@
		#svs_outfile = os.path.splitext(file)[0]+".fpfiltered.txt"

		#sv.to_csv(outfile,sep='\t', index=False)
		te_all.append(te) ## remove
	te_all = pd.concat(te_all, ignore_index=True)
	te_all.to_csv(os.path.join(tedir,"all_te_"+suffix+".txt"), sep='\t', index=False)
	return(te_all)

def te_intersect_meth(te, tedir, window):
	te = pd.read_csv(os.path.join(tedir,"all_te.txt"),delimiter='\t')
	te['SV_chrom'] = "chr"+ te['SV_chrom'].astype(str)
	bed = te[['SV_chrom','SV_start','SV_end','AnnotSV_ID']]
	bed.to_csv(os.path.join(tedir,"all_te.bed"), sep='\t', header=False,index=False)
	os.system("bedtools window -a " + os.path.join(tedir,"all_te.bed") + " -b /hpf/largeprojects/davidm/vsubasri/transposable_elements/meth_manifest.bed -w "+ window + " > "+ os.path.join(tedir,"te_meth_"+window+".bed"))
	te_meth = pd.read_csv(os.path.join(tedir,"te_meth_"+window+".bed"),delimiter='\t',header=None)
	te_meth = te_meth[[3,7]]; te_meth.columns = ["AnnotSV_ID","probe"]
	te_meth = te_meth.groupby("AnnotSV_ID")["probe"].apply(list).reset_index(name="probes")
	te_meth['probes'] = [','.join(set(map(str, l))) for l in te_meth['probes']] 
#	te_meth = pd.concat([pd.DataFrame(te_meth["AnnotSV_ID"],columns=["AnnotSV_ID"]),te_meth.probes.apply(pd.Series).add_prefix('probe_')])
	te_meth.to_csv(os.path.join(tedir,"te_meth_"+window+".txt"), sep='\t',index=False,quoting = csv.QUOTE_NONE)
	return(te_meth)

#def sample_probes(n):
	

def check_cpg(gene):
        cp_genes_list= open("/hpf/largeprojects/davidm/resources/actionable_cpgs.txt", "r")
        cp_genes = [x for x in cp_genes_list.read().split('\n') if x]
        if gene in cp_genes:
                cpg='CPG'
        else:
                cpg='.'
        return(cpg)

def annotate_cpg(tedir,suffix):
        files = glob.glob(os.path.join(tedir, "*/*"+suffix+".tsv"))
        for file in files:
                outfile= os.path.splitext(file)[0]+".final.txt"
                if os.path.exists(outfile):
                    continue
                else:
                        print(outfile)
                        with open(file,'r') as in_f, open(outfile, 'w') as out_f:
                                fl_indices=in_f.readline().strip().split("\t") + ['CPG']
                                firstline='\t'.join(fl_indices)+"\n"
                                out_f.write(firstline)
                                for line in in_f:
                                        line_indices = line.strip().split("\t")
                                        genes = line_indices[16].split('/')
                                        cpgs = []
                                        for gene in genes:
                                                cpg_tmp = check_cpg(gene)
                                                cpgs.append(cpg_tmp)
                                        if "CPG" in cpgs:
                                                newline = '\t'.join(line_indices + ['CPG']) + '\n'
                                        else:
                                                newline = '\t'.join(line_indices + ['.']) + '\n'
                                        out_f.write(newline)


def main():
        parser = argparse.ArgumentParser(
                description='Filter structural variants',
                formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        )
        parser.add_argument("-s", '--suffix', action="store", dest="suffix", help='file ending')
        parser.add_argument("-t", '--tedir', action="store", dest="tedir", help='path to structural variants')
        parser.add_argument("-i", '--intersect', action = "store", type=str, dest="intersect",help='one of the following: meth')
        args = parser.parse_args()
        print("Annotate with CPG...")
        annotate_cpg(args.tedir,args.suffix)
        print("Combine...")
        te = combine(args.tedir,args.suffix)
        if args.intersect:
                print("Interect Transposable Elements with Methylation...")
                te_intersect_meth(te,args.tedir, args.intersect)

if __name__ == '__main__':
	main()

