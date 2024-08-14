#!/usr/bin/env python
import pandas as pd
import os

df = pd.read_csv('NT-pairs.csv', sep='\t')
output_dir = "/hpf/largeprojects/davidm/shilpa/TE_runs/fixed_bams/"

for _,row in df.iterrows():
    sample_id = row['sample_id']
    tumor_bam = row['tumor_samples']
    germline_bam = row['germline_samples']

    #output_file = f'config_trafic_{sample_id}.yaml'
    #base_sample_id = full_sample_id.replace('sorted_fixed_', '').replace('.bam', '')
    output_file = os.path.join(output_dir, f'config_trafic_{sample_id}.yaml')
    yaml_content = f"""user_id : shilpa
project: TE-project
sample_name: {sample_id}
tumour_name: {sample_id}_T
tumour_bam: /trafic_input/{tumor_bam.split('/')[-1]}
normal_name: {sample_id}_N
normal_bam: /trafic_input/{germline_bam.split('/')[-1]}
output_dir: outputs
bwa_threads: 16
"""
    with open(output_file, 'w') as file:
        file.write(yaml_content)
        
    print(f'Config file created for {sample_id}')
