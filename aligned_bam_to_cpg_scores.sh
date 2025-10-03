#!/bin/bash

# --- Configuration ---
palladium_bam_dir="/scratch/ucgd/lustre-labs/quinlan/data-shared/read-backed-phasing"
base_output_dir="/scratch/ucgd/lustre-labs/quinlan/data-shared/dna-methylation"
bin_dir="/uufs/chpc.utah.edu/common/HIPAA/u6018199/pb-CpG-tools-v3.0.0-x86_64-unknown-linux-gnu/bin"

# --- Add pb-CpG-tools to the PATH ---
export PATH=${bin_dir}:$PATH

# aligned_bam_to_cpg_scores --help

# --- Get Sample Prefixes ---
echo "Getting Palladium prefixes..."
prefixes=$(python src/util/get_palladium_prefixes.py)
# prefixes="200081" 
echo "... Got prefixes"

# --- Main Loop for Each Mode ---
for mode in "model" "count"; do
    
    # Define the output directory for the current mode
    output_dir="${base_output_dir}/CEPH1463.GRCh38.hifi.${mode}.read-backed-phased"
    echo "----------------------------------------------------"
    echo "Starting mode: '${mode}'"
    echo "Output will be written to: ${output_dir}"
    echo "----------------------------------------------------"

    # Create the directory if it doesn't exist
    mkdir -p ${output_dir}

    # Loop through each sample prefix and submit the job
    for prefix in $prefixes; do
        echo "Submitting job for sample '${prefix}' in '${mode}' mode..."
        
        nohup aligned_bam_to_cpg_scores \
            --bam "${palladium_bam_dir}/${prefix}.GRCh38.haplotagged.bam" \
            --output-prefix "${output_dir}/${prefix}.GRCh38.haplotagged" \
            --threads 8 \
            --min-coverage 10 \
            --min-mapq 1 \
            --pileup-mode "${mode}" \
            > /dev/null 2>&1 &
    done
done

echo "All jobs have been submitted."
