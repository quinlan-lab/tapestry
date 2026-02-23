export PYTHONPATH="src/util"

# INPUT DIRS
reference="/path/to/reference.fa"
meth_trio_phased_dir="/path/to/trio-founder-phased"           # output dir of trio_phase_meth_to_parent_haps.py
meth_count_haplotagged_dir="/path/to/meth/count/haplotagged"  # output dir of aligned_bam_to_cpg_scores (count-based unphased)
meth_model_haplotagged_dir="/path/to/meth/model/haplotagged"  # output dir of aligned_bam_to_cpg_scores (model-based unphased)
vcf_trio_phased="/path/to/whatshap-phasing/trio_phased.vcf.gz"

# OUTPUT DIR
output_dir="/path/to/trio-founder-phased-all-cpgs"

mkdir -p ${output_dir}

echo "Writing all CpG sites in reference genome to '${output_dir}' ..."

bed_all_cpgs_in_reference="${output_dir}/all_cpg_sites_in_reference.bed"

python -m memory_profiler src/write_all_cpgs_in_reference.py \
    --reference ${reference} \
    --bed_all_cpgs_in_reference ${bed_all_cpgs_in_reference} \
    > ${output_dir}/write_all_cpgs_in_reference.log 2>&1

# Child sample ID
uid="kid"

# INPUT FILES
bed_meth_trio_phased="${meth_trio_phased_dir}/${uid}.dna-methylation.founder-phased.bed"
bed_meth_count_unphased="${meth_count_haplotagged_dir}/${uid}.GRCh38.haplotagged.combined.bed.gz"
bed_meth_model_unphased="${meth_model_haplotagged_dir}/${uid}.GRCh38.haplotagged.combined.bed.gz"

# OUTPUT FILES
bed_meth_trio_phased_all_cpgs="${output_dir}/${uid}.dna-methylation.founder-phased.all_cpgs.bed"

echo "Started ${uid} ..."

python src/expand_to_all_cpgs.py \
    --bed_all_cpgs_in_reference ${bed_all_cpgs_in_reference} \
    --bed_meth_founder_phased ${bed_meth_trio_phased} \
    --bed_meth_count_unphased ${bed_meth_count_unphased} \
    --bed_meth_model_unphased ${bed_meth_model_unphased} \
    --bed_meth_founder_phased_all_cpgs ${bed_meth_trio_phased_all_cpgs} \
    --uid ${uid} \
    --vcf_joint_called ${vcf_trio_phased} \
    > ${output_dir}/${uid}.log 2>&1

echo "... Finished ${uid}"
echo "Check ${output_dir} for logs and outputs."
echo "Done running expand_to_all_cpgs_trio.sh"
