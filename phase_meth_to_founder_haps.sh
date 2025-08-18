pb_cpg_tool_mode='model' # mode of aligned_bam_to_cpg_scores

# INPUT DIRS 
read_phased_dir='/scratch/ucgd/lustre-labs/quinlan/data-shared/read-backed-phasing' # output dir of hiphase
iht_phased_dir='/scratch/ucgd/lustre-labs/quinlan/data-shared/haplotype-maps/CEPH1463.GRCh38' # output dir of gtg-ped-map/gtg-concordance
meth_read_phased_dir="/scratch/ucgd/lustre-labs/quinlan/data-shared/dna-methylation/CEPH1463.GRCh38.hifi.${pb_cpg_tool_mode}.read-backed-phased" # output dir of aligned_bam_to_cpg_scores

# OUTPUT DIRS 
output_dir="/scratch/ucgd/lustre-labs/quinlan/data-shared/dna-methylation/CEPH1463.GRCh38.hifi.${pb_cpg_tool_mode}.founder-phased" # output dir of phase_meth_to_founder_haps.py
# output_dir="/scratch/ucgd/lustre-labs/quinlan/data-shared/dna-methylation/CEPH1463.GRCh38.hifi.${pb_cpg_tool_mode}.founder-phased.test" # testing 

mkdir -p ${output_dir}

# TODO: use a ped file to get the prefixes for non-founder individuals in the pedigree
prefixes=$(python src/util/get_palladium_prefixes.py)
echo "Done loading prefixes..."

# prefixes="200080 200081" # testing 

for prefix in $prefixes; do
    uid="${prefix}" # sample ID in joint-called multi-sample vcf (see below)

    # INPUT FILES 
    vcf_read_phased="${read_phased_dir}/${uid}.GRCh38.deepvariant.glnexus.phased.vcf.gz" # single-sample vcf from hiphase
    tsv_read_phase_blocks="${read_phased_dir}/${uid}.GRCh38.hiphase.blocks.tsv" # single-sample tsv from hiphase
    vcf_iht_phased="${iht_phased_dir}/CEPH1463.GRCh38.pass.sorted.vcf.gz" # joint-called multi-sample vcf from gtg-ped-map/gtg-concordance
    txt_iht_blocks="${iht_phased_dir}/CEPH1463.GRCh38.iht.sorted.txt" # multi-sample iht blocks file from gtg-ped-map/gtg-concordance
    bed_meth_hap1="${meth_read_phased_dir}/${uid}.GRCh38.haplotagged.hap1.bed.gz" # bed file from aligned_bam_to_cpg_scores for hap1
    bed_meth_hap2="${meth_read_phased_dir}/${uid}.GRCh38.haplotagged.hap2.bed.gz" # bed file from aligned_bam_to_cpg_scores for hap2

    nohup python src/phase_meth_to_founder_haps.py \
        --pb_cpg_tool_mode ${pb_cpg_tool_mode} \
        --uid ${uid} \
        --vcf_read_phased ${vcf_read_phased} \
        --tsv_read_phase_blocks ${tsv_read_phase_blocks} \
        --vcf_iht_phased ${vcf_iht_phased} \
        --txt_iht_blocks ${txt_iht_blocks} \
        --bed_meth_hap1 ${bed_meth_hap1} \
        --bed_meth_hap2 ${bed_meth_hap2} \
        --output_dir ${output_dir} \
        > ${output_dir}/${uid}.log 2>&1 &

    echo "Started ${uid} ..."
done

echo "All phase_meth_to_founder_haps.py jobs started."
echo "Check ${output_dir} for logs and outputs."
