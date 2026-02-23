#!/bin/bash

# --- Configuration ---
whatshap_dir="/scratch/ucgd/lustre-labs/quinlan/data-shared/whatshap-phasing"
kid_sample_id="NA12883"
kid_haplotagged_bam="${whatshap_dir}/${kid_sample_id}.haplotagged.bam"
base_output_dir="/scratch/ucgd/lustre-labs/quinlan/data-shared/dna-methylation"
bin_dir="/uufs/chpc.utah.edu/common/HIPAA/u6018199/pb-CpG-tools-v3.0.0-x86_64-unknown-linux-gnu/bin"

# --- Add pb-CpG-tools to the PATH ---
export PATH=${bin_dir}:$PATH

# aligned_bam_to_cpg_scores --help

for mode in "model" "count"; do
    
    # Define the output directory for the current mode
    output_dir="${base_output_dir}/CEPH1463.GRCh38.hifi.${mode}.trio.iht-phased"
    echo "----------------------------------------------------"
    echo "Starting mode: '${mode}'"
    echo "Output will be written to: ${output_dir}"
    echo "----------------------------------------------------"

    # Create the directory if it doesn't exist
    mkdir -p ${output_dir}

    # Submit the job
	echo "Submitting job for sample '${kid_sample_id}' in '${mode}' mode..."
	nohup aligned_bam_to_cpg_scores \
		--bam $kid_haplotagged_bam \
		--output-prefix "${output_dir}/${kid_sample_id}.GRCh38.haplotagged" \
		--threads 8 \
		--min-coverage 10 \
		--min-mapq 1 \
		--pileup-mode "${mode}" \
		> output_dir/aligned_bam_to_cpg_scores.$kid_sample_id.log 2>&1 &
	
done

echo "All aligned_bam_to_cpg_scores jobs have been submitted."
