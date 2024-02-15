#!/usr/bin/env python

import argparse
import csv
import pandas as pd
import os
import subprocess
import re 


def parse_bams(bams_input):
        if os.path.isfile(bams_input):
                bams=pd.read_csv(bams_input, sep='\n', header=None)[0]
        elif os.path.isdir(bams_input):
                bams=glob.glob(os.path.join(bams_input, "*.bam"))
        else:
                raise IOError("Wrong bam input")
        return(bams)

def create_input_file(bams_input,outdir):
	bams = parse_bams(bams_input)
	samples = [] ; map = []
	for bam in bams:
		sample= regex(bam)
		print("Sample: ",sample)
		samples.append(sample)
		map.append([sample, bam])
	samples = pd.DataFrame(samples)
	samples.to_csv(os.path.join(outdir,"sample_id.txt"), sep='\t', header=False, index=False)
	map = pd.DataFrame(map)
	map.to_csv(os.path.join(outdir,"illumina_bam_list.txt"), sep='\t', header=False, index=False)	

def regex(path):
	if 'merged.bam' in path:
		found = re.search('(.*?).merged.bam', os.path.basename(path))
	elif 'realigned-recalibrated.bam' in path:
		found = re.search('(.*?).realigned-recalibrated.bam', os.path.basename(path))
	newname = found.group(1)  
	return newname

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Create map file for combining vcfs')
	parser.add_argument('-b','--bams_input', required=True, help='Name of file or directory with bam files')
	parser.add_argument('-o','--outdir', required=True, help='Name of output directory ending with a slash ie. "/Users/vsubasri/testdir/')
	args = parser.parse_args()

	create_input_file(args.bams_input,args.outdir)
