export PYTHONPATH="src:src/util"

# INPUT DIRS
whatshap_phasing_dir="/path/to/whatshap-phasing"           # output dir of run-whatshap-phase.sh
meth_count_haplotagged_dir="/path/to/meth/count/haplotagged"  # output dir of aligned_bam_to_cpg_scores in count mode
meth_model_haplotagged_dir="/path/to/meth/model/haplotagged"  # output dir of aligned_bam_to_cpg_scores in model mode

# OUTPUT DIR
output_dir="/path/to/trio-founder-phased"

mkdir -p ${output_dir}

# Pedigree file
ped="/path/to/trio.ped"

# Child sample ID
uid="kid"

# INPUT FILES
vcf_trio_phased="${whatshap_phasing_dir}/trio_phased.vcf.gz"
bed_meth_count_hap1="${meth_count_haplotagged_dir}/${uid}.GRCh38.haplotagged.hap1.bed.gz"
bed_meth_count_hap2="${meth_count_haplotagged_dir}/${uid}.GRCh38.haplotagged.hap2.bed.gz"
bed_meth_model_hap1="${meth_model_haplotagged_dir}/${uid}.GRCh38.haplotagged.hap1.bed.gz"
bed_meth_model_hap2="${meth_model_haplotagged_dir}/${uid}.GRCh38.haplotagged.hap2.bed.gz"

python src/trio_phase_meth_to_parent_haps.py \
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
