#!/bin/bash
#SBATCH --job-name=melt
#SBATCH --time=100:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --mem=4G
#SBATCH --output=%x-%j.out
#SBATCH --constraint=CentOS7
#SBATCH --mail-type=ALL

SF=Snakefile
CP="/hpf/largeprojects/davidm/resources/snakemake"
CONFIG=config.yaml

source /hpf/largeprojects/davidm/resources/etc/profile.d/conda.sh
conda activate snakemake

snakemake --use-conda -s ${SF} --cores 4 --conda-prefix ${CP} --configfile ${CONFIG} --printshellcmds

