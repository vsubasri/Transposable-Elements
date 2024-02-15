#!/bin/bash
#SBATCH --job-name=melt
#SBATCH --time=100:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --mem=4G
#SBATCH --output=%x-%j.out
#SBATCH --constraint=CentOS7
#SBATCH --mail-type=ALL

SF=Snakefile
CP="/home/vsubasri/.conda/"
SLURM=~/slurm_profile/
CONFIG=config.yaml

source ~/anaconda3/etc/profile.d/conda.sh

module load miniconda/4.4.10 snakemake/20230426

snakemake --use-conda -s ${SF} --cores 4 --conda-prefix ${CP} --configfile ${CONFIG} --profile ${SLURM} --printshellcmds
