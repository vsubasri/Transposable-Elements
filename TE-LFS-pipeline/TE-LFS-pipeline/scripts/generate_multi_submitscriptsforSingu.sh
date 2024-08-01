#!/bin/bash

if [[ $# -ne 1 ]]; then
  echo "Usage: $0 <script_names.txt>"
  exit 1
fi

input_file=$1

if [[ ! -f $input_file ]]; then
  echo "File $input_file not found!"
  exit 1
fi

counter=1
while IFS= read -r script_name; do
  if [[ -n $script_name ]]; then
    submit_script_content=$(cat << EOF
#!/bin/bash

#SBATCH -J xTea_test2
#SBATCH -e %x.e%j
#SBATCH -o %x.o%j
#SBATCH -N 4 -c 16
#SBATCH --mem=16G
#SBATCH -t 47:59:00

set -euxo 
cd $PWD

module load Singularity/3.11.3

# Define your data and singularity image paths
SINGULARITY_IMAGE="/home/shilpa/xteam_image.sif"

SCRIPT="$script_name"

singularity exec -e -B /hpf:/hpf \$SINGULARITY_IMAGE bash \$SCRIPT
EOF
    )

    # Extract script name without extension
    script_name_no_ext="${script_name##*/}"
    script_name_no_ext="${script_name_no_ext%.*}"

    output_file="submit_singu_${script_name_no_ext}_${counter}.sh"
    ((counter++))

    echo "$submit_script_content" > "$output_file"
    chmod +x "$output_file"
    echo "Created $output_file"
  fi
done < "$input_file"
