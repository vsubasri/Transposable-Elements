# TE-LFS-pipeline

Transposable Elements prediction pipeline.

The TE-LFS-pipeline is designed for the prediction and analysis of transposable elements (TEs) in genomic data. The pipeline currently supports germline variant calling with the script `TE-LFS-germline-pipe.sh`. 

## Usage

```bash
bash TE-LFS-germline-pipe.sh /path/to/bamfilefolder
```

If no path is provided, the default path on HPC will be used.

## Output directory Structure

When you run the pipeline, the following directories will be created in the CWD:

1. **fixed_bams/**: Contains preprocessed BAM files, including steps such as sorting and adding MC/MQ tags.
2. **results/**: Contains output from each tool, organized into subdirectories for each BAM file.

The `results/` directory will have subdirectories for each tool:
- `xtea_runs/`
- `melt_runs/`
- `ins_runs/`

Each of these tool directories will have subdirectories corresponding to each BAM file processed.

## Future Plans

1. A similar script will be created for tumor variant calling: `TE-LFS-tumour-pipe.sh`.
2. New scripts will be developed for hg38 variant calling: for germline variant calling and for tumor variant calling.

