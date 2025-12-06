export PYTHONPATH="src/util" 

# INPUT DIRS 
reference="/scratch/ucgd/lustre-labs/quinlan/data-shared/constraint-tools/reference/grch38/hg38.analysisSet.fa" # /scratch/ucgd/lustre-labs/quinlan/u6018199/constraint-tools/download-process-data/download-reference-grch38.sh
meth_founder_phased_dir="/scratch/ucgd/lustre-labs/quinlan/data-shared/dna-methylation/CEPH1463.GRCh38.hifi.founder-phased" # output dir of phase_meth_to_founder_haps.py (containing count-based and model-based founder-phased meth)
meth_count_read_phased_dir="/scratch/ucgd/lustre-labs/quinlan/data-shared/dna-methylation/CEPH1463.GRCh38.hifi.count.read-backed-phased" # output dir of aligned_bam_to_cpg_scores (containing count-based unphased meth)
meth_model_read_phased_dir="/scratch/ucgd/lustre-labs/quinlan/data-shared/dna-methylation/CEPH1463.GRCh38.hifi.model.read-backed-phased" # output dir of aligned_bam_to_cpg_scores (countaining model-based unphased meth)

# OUTPUT DIRS 
output_dir="/scratch/ucgd/lustre-labs/quinlan/data-shared/dna-methylation/CEPH1463.GRCh38.hifi.founder-phased.all-cpgs" 

mkdir -p ${output_dir}

echo "Writing all CpG sites in reference genome to ${output_dir} ..." 

# OUTPUT FILE
bed_all_cpgs_in_reference="${output_dir}/all_cpg_sites_in_reference.bed" # output of src/write_all_cpgs_in_reference.py

python src/write_all_cpgs_in_reference.py \
	--reference ${reference} \
	--bed_all_cpgs_in_reference ${bed_all_cpgs_in_reference} \
	> ${output_dir}/write_all_cpgs_in_reference.log 2>&1

# TODO: use a ped file to get the prefixes for non-founder individuals in the pedigree
prefixes=$(python src/util/get_palladium_prefixes.py)
echo "Done loading prefixes..."

# prefixes="200081" # TESTING

for prefix in $prefixes; do
    uid="${prefix}" # sample ID in joint-called multi-sample vcf

    # INPUT FILES 
	bed_meth_founder_phased="${meth_founder_phased_dir}/${uid}.dna-methylation.founder-phased.bed" # bed file of founder-phased methylation levels from src/phase_meth_to_founder_haps.py
	bed_meth_count_unphased="${meth_count_read_phased_dir}/${uid}.GRCh38.haplotagged.combined.bed.gz" # bed file from aligned_bam_to_cpg_scores (unphased count-based meth)
	bed_meth_model_unphased="${meth_model_read_phased_dir}/${uid}.GRCh38.haplotagged.combined.bed.gz" # bed file from aligned_bam_to_cpg_scores (unphased model-based meth)
	bed_het_site_mismatches="${meth_founder_phased_dir}/${uid}.bit-vector-sites-mismatches.bed" # bed file of heterozygous sites at which bit-vectors are mismatched, from src/phase_meth_to_founder_haps.py
	vcf_joint_called="/scratch/ucgd/lustre-labs/quinlan/data-shared/datasets/Palladium/deepvariant/CEPH-1463.joint.GRCh38.deepvariant.glnexus.phased.vcf.gz" # joint-called multi-sample vcf 

	# OUTPUT FILES 
	bed_meth_founder_phased_all_cpgs="${output_dir}/${uid}.dna-methylation.founder-phased.all_cpgs.bed"

    nohup python src/expand_to_all_cpgs.py \
		--bed_all_cpgs_in_reference ${bed_all_cpgs_in_reference} \
        --bed_meth_founder_phased ${bed_meth_founder_phased} \
        --bed_meth_count_unphased ${bed_meth_count_unphased} \
        --bed_meth_model_unphased ${bed_meth_model_unphased} \
        --bed_meth_founder_phased_all_cpgs ${bed_meth_founder_phased_all_cpgs} \
		--bed_het_site_mismatches ${bed_het_site_mismatches} \
		--uid ${uid} \
        --vcf_joint_called ${vcf_joint_called} \
        > ${output_dir}/${uid}.log 2>&1 & 

    echo "Started ${uid} ..."
done

echo "All expand_to_all_cpgs.py processes started."
echo "Check ${output_dir} for logs and outputs."
