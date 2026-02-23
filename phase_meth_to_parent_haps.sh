export PYTHONPATH="src:src/util"

# INPUT DIRS
whatshap_phasing_dir="/scratch/ucgd/lustre-labs/quinlan/data-shared/whatshap-phasing"  # output dir of run-whatshap-phase.sh
dna_meth_dir="/scratch/ucgd/lustre-labs/quinlan/data-shared/dna-methylation"
meth_count_haplotagged_dir="${dna_meth_dir}/CEPH1463.GRCh38.hifi.count.trio.iht-phased"  # output dir of aligned_bam_to_cpg_scores in count mode
meth_model_haplotagged_dir="${dna_meth_dir}/CEPH1463.GRCh38.hifi.model.trio.iht-phased"  # output dir of aligned_bam_to_cpg_scores in model mode

# OUTPUT DIR
output_dir="/scratch/ucgd/lustre-labs/quinlan/data-shared/dna-methylation/CEPH1463.GRCh38.hifi.parent-phased.trio" 

mkdir -p ${output_dir}

# Pedigree file
ped="trio.ped"

# Child sample ID
uid="NA12883"

# INPUT FILES
vcf_trio_phased="${whatshap_phasing_dir}/trio_phased.vcf.gz"
bed_meth_count_hap1="${meth_count_haplotagged_dir}/${uid}.GRCh38.haplotagged.hap1.bed.gz"
bed_meth_count_hap2="${meth_count_haplotagged_dir}/${uid}.GRCh38.haplotagged.hap2.bed.gz"
bed_meth_model_hap1="${meth_model_haplotagged_dir}/${uid}.GRCh38.haplotagged.hap1.bed.gz"
bed_meth_model_hap2="${meth_model_haplotagged_dir}/${uid}.GRCh38.haplotagged.hap2.bed.gz"

python src/phase_meth_to_parent_haps.py \
    --uid ${uid} \
    --vcf_trio_phased ${vcf_trio_phased} \
    --ped ${ped} \
    --bed_meth_count_hap1 ${bed_meth_count_hap1} \
    --bed_meth_count_hap2 ${bed_meth_count_hap2} \
    --bed_meth_model_hap1 ${bed_meth_model_hap1} \
    --bed_meth_model_hap2 ${bed_meth_model_hap2} \
    --output_dir ${output_dir} \
    > ${output_dir}/${uid}.log 2>&1

echo "Done. Check ${output_dir}/${uid}.log for output."
