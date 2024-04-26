#!/bin/bash

#SBATCH --job-name conda
#SBATCH -t 50:00:00
#SBATCH -N 1 -c 10
#SBATCH --mem=40g


#conda create -n my_snakemake python=3.9.13 snakemake -y
conda env create --quiet --file "/hpf/largeprojects/davidm/blaverty/te/scripts/Transposable-Elements/te_discovery/envs/melt.yaml" --prefix "/desired/prefix/for/environment"



