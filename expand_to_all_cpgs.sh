export PYTHONPATH="src/util" 

pb_cpg_tool_mode='model' # mode of aligned_bam_to_cpg_scores

# INPUT DIRS 
reference="/scratch/ucgd/lustre-labs/quinlan/data-shared/constraint-tools/reference/grch38/hg38.analysisSet.fa" # /scratch/ucgd/lustre-labs/quinlan/u6018199/constraint-tools/download-process-data/download-reference-grch38.sh
meth_founder_phased_dir="/scratch/ucgd/lustre-labs/quinlan/data-shared/dna-methylation/CEPH1463.GRCh38.hifi.${pb_cpg_tool_mode}.founder-phased" # output dir of phase_meth_to_founder_haps.py

# OUTPUT DIRS 
output_dir="/scratch/ucgd/lustre-labs/quinlan/data-shared/dna-methylation/CEPH1463.GRCh38.hifi.${pb_cpg_tool_mode}.founder-phased.all-cpgs" 

mkdir -p ${output_dir}

echo "Writing all CpG sites to ${output_dir} ..." 

# OUTPUT FILE
bed_all_cpgs="${output_dir}/all_cpg_sites.bed" # output of src/write_all_cpgs.py

# TODO: uncomment 
# python src/write_all_cpgs.py \
# 	--reference ${reference} \
# 	--bed_all_cpgs ${bed_all_cpgs} \
# 	> ${output_dir}/write_all_cpgs.log 2>&1

# # TODO: use a ped file to get the prefixes for non-founder individuals in the pedigree
# prefixes=$(python src/util/get_palladium_prefixes.py)
# echo "Done loading prefixes..."

prefixes="200081" # TESTING 

for prefix in $prefixes; do
    uid="${prefix}" # sample ID in joint-called multi-sample vcf

    # INPUT FILES 
	bed_hap_map="${meth_founder_phased_dir}/${uid}.hap-map-blocks.bed" # bed file of hap-map blocks from src/phase_meth_to_founder_haps.py
	bed_meth_founder_phased="${meth_founder_phased_dir}/${uid}.dna-methylation.founder-phased.bed" # bed file of founder-phased methylation levels from src/phase_meth_to_founder_haps.py

	# OUTPUT FILES 
	bed_meth_founder_phased_all_cpgs="${output_dir}/${uid}.dna-methylation.founder-phased.all_cpgs.bed"

    python src/expand_to_all_cpgs.py \
		--bed_all_cpgs ${bed_all_cpgs} \
		--bed_hap_map ${bed_hap_map} \
        --bed_meth_founder_phased ${bed_meth_founder_phased} \
        --bed_meth_founder_phased_all_cpgs ${bed_meth_founder_phased_all_cpgs} \
        > ${output_dir}/${uid}.log 2>&1

    echo "Started ${uid} ..."
done

echo "All expand_to_all_cpgs.py processes started."
echo "Check ${output_dir} for logs and outputs."
