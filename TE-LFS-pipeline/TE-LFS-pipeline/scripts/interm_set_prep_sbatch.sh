while read -r line; do
    if [[ $line == sbatch* ]]; then
        script_path=$(echo "$line" | awk '{print $3}')
    	{ cat slurm_header_xtea.txt; cat "$script_path"; } > tmp && mv tmp "$script_path"
    fi
done < submit_scripts.sh
