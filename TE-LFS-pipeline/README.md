# TE-LFS-pipeline

Transposable Elements prediction pipeline.

The TE-LFS-pipeline is designed for the prediction and analysis of transposable elements (TEs) in genomic data. The pipeline currently supports germline variant calling with the script `TE-LFS-germline-pipe.sh`. 

## Usage

```bash
cd $PWD
bash /pathtopipeline/TE-LFS-germline-pipe.sh /path/to/bamfilefolder

```

If no path is provided, the default path on HPC will be used as mentioned in the bash script.

## Output directory Structure

When you run the pipeline, the following directories will be created in the $PWD

1. **fixed_bams/**: Contains preprocessed BAM files, including steps such as sorting and adding MC/MQ tags.
2. **results/**: Contains output from each tool, organized into subdirectories for each BAM file.

The `results/` directory will have subdirectories for each tool:
- `xtea_runs/`
- `melt_runs/`
- `ins_runs/`

Logs of individual runs are inside the _runs_ folder.

Each of these tool directories will have subdirectories corresponding to each BAM file processed.

Other pipeline scripts have a very similar setup:

TE-LFS-germline-pipe.sh : with hg19 reference genome
TE-LFS-germline-pipe-hg38.sh : with hg38 reference genome

TE-LFS-tumor-pipe.sh : with hg19 reference genome : This is to process germline and tumor samples. Melt is not available here. xTea's case-control feature has been used. New tool added is: TraFiC. 
TE-LFS-tumor-pipe-hg38.sh : with hg38 reference genome : Same as tumor-pipe.sh, TraFiC not available for hg38 since its an old tool.


