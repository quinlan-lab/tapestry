#!/bin/bash

# Usage:
#   Production:  ./expand_to_all_cpgs.trio.sh
#   Dev mode:    ./expand_to_all_cpgs.trio.sh --dev-dir trio_dev_data

source src/util/logging.sh

export PYTHONPATH="src/util"

DEV_DIR=""

# --- Argument Parsing ---
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --dev-dir|-d)
            DEV_DIR="${2%/}"
            shift 2
            ;;
        *)
            echo "Error: Unknown parameter: $1"
            echo "Usage: $0 [--dev-dir <DEV_DATA_DIR>]"
            exit 1
            ;;
    esac
done

# --- Configuration ---
kid_id="NA12878"
dad_id="NA12891"
mom_id="NA12892"

# --- Default Configurations (Production) ---

# INPUT DIRS
reference="/scratch/ucgd/lustre-labs/quinlan/data-shared/constraint-tools/reference/grch38/hg38.analysisSet.fa" # /scratch/ucgd/lustre-labs/quinlan/u6018199/constraint-tools/download-process-data/download-reference-grch38.sh
meth_parent_phased_dir="/scratch/ucgd/lustre-labs/quinlan/data-shared/dna-methylation/CEPH1463.GRCh38.hifi.parent-phased" # output dir of phase_meth_to_parent_haps.py (containing count-based and model-based parent-phased meth)
meth_count_pedmec_phased_dir="/scratch/ucgd/lustre-labs/quinlan/data-shared/dna-methylation/CEPH1463.GRCh38.hifi.count.pedmec-phased" # output dir of aligned_bam_to_cpg_scores (containing count-based unphased meth)
meth_model_pedmec_phased_dir="/scratch/ucgd/lustre-labs/quinlan/data-shared/dna-methylation/CEPH1463.GRCh38.hifi.model.pedmec-phased" # output dir of aligned_bam_to_cpg_scores (containing model-based unphased meth)

# OUTPUT DIRS
# This is passed as an arg to src/expand_to_all_cpgs_trio.py so it knows where to write the unphased bigwig files 
output_dir="/scratch/ucgd/lustre-labs/quinlan/data-shared/dna-methylation/CEPH1463.GRCh38.hifi.parent-phased.all-cpgs"

# --- Optional Dev Data Overrides ---
if [ -n "$DEV_DIR" ]; then
    log_info "DEV MODE ENABLED: Reading from and writing to ${DEV_DIR}"

    # INPUT DIRS
    meth_dir="${DEV_DIR}/output/dna-methylation"
    meth_parent_phased_dir="${meth_dir}/CEPH1463.GRCh38.hifi.parent-phased" # output dir of phase_meth_to_parent_haps.py
    meth_count_pedmec_phased_dir="${meth_dir}/CEPH1463.GRCh38.hifi.count.pedmec-phased" # output dir of aligned_bam_to_cpg_scores (containing count-based unphased meth)
    meth_model_pedmec_phased_dir="${meth_dir}/CEPH1463.GRCh38.hifi.model.pedmec-phased" # output dir of aligned_bam_to_cpg_scores (containing model-based unphased meth)

    # Use dev chromosome reference so that write_all_cpgs_in_reference.py only outputs CpGs on that chromosome
    source "$(dirname "$0")/trio_dev_data_config.sh"
    dev_chrom="${DEV_REGION%%:*}"
    reference="${DEV_DIR}/input/${dev_chrom}.fa"

    # OUTPUT DIRS
    output_dir="${meth_dir}/CEPH1463.GRCh38.hifi.parent-phased.all-cpgs"
fi

mkdir -p ${output_dir}

log_info "Writing all CpG sites in reference genome to '${output_dir}' ..."

# OUTPUT FILE
bed_all_cpgs_in_reference="${output_dir}/all_cpg_sites_in_reference.bed" # output of src/write_all_cpgs_in_reference.py

python -m memory_profiler src/write_all_cpgs_in_reference.py \
	--reference ${reference} \
	--bed_all_cpgs_in_reference ${bed_all_cpgs_in_reference} \
	> ${output_dir}/write_all_cpgs_in_reference.log 2>&1

meth_bed() {
    # Usage: meth_bed <mode> <uid> <hap>
    local mode="$1" uid="$2" hap="$3"
    if [ "$mode" = "count" ]; then
        echo "${meth_count_pedmec_phased_dir}/${uid}.GRCh38.haplotagged.${hap}.bed.gz"
    else
        echo "${meth_model_pedmec_phased_dir}/${uid}.GRCh38.haplotagged.${hap}.bed.gz"
    fi
}

# INPUT FILES
bed_meth_parent_phased="${meth_parent_phased_dir}/trio.dna-methylation.bed" # parent-phased methylation levels from phase_meth_to_parent_haps.py
bed_het_site_mismatches_pat="${meth_parent_phased_dir}/${kid_id}.bit-vector-sites-mismatches.paternal.bed"
bed_het_site_mismatches_mat="${meth_parent_phased_dir}/${kid_id}.bit-vector-sites-mismatches.maternal.bed"
vcf_joint_called="/scratch/ucgd/lustre-labs/quinlan/data-shared/datasets/Palladium/deepvariant/CEPH-1463.joint.GRCh38.deepvariant.glnexus.vcf.gz" # joint-called multi-sample vcf 

# Optional Dev Data Overrides for input files
if [ -n "$DEV_DIR" ]; then
    vcf_joint_called="${DEV_DIR}/input/CEPH-1463.joint.GRCh38.deepvariant.glnexus.vcf.gz"
fi

# OUTPUT FILE
bed_meth_parent_phased_all_cpgs="${output_dir}/trio.dna-methylation.all_cpgs.bed"

log_info "Expanding methylation to all CpG sites ..."

PYTHONPATH=src:src/util python src/expand_to_all_cpgs_trio.py \
    --bed_all_cpgs_in_reference "$bed_all_cpgs_in_reference" \
    --bed_meth_parent_phased "$bed_meth_parent_phased" \
    --bed_meth_count_unphased_kid "$(meth_bed count "$kid_id" combined)" \
    --bed_meth_model_unphased_kid "$(meth_bed model "$kid_id" combined)" \
    --bed_meth_count_unphased_dad "$(meth_bed count "$dad_id" combined)" \
    --bed_meth_model_unphased_dad "$(meth_bed model "$dad_id" combined)" \
    --bed_meth_count_unphased_mom "$(meth_bed count "$mom_id" combined)" \
    --bed_meth_model_unphased_mom "$(meth_bed model "$mom_id" combined)" \
    --bed_meth_parent_phased_all_cpgs "$bed_meth_parent_phased_all_cpgs" \
    --bed_het_site_mismatches_pat "$bed_het_site_mismatches_pat" \
    --bed_het_site_mismatches_mat "$bed_het_site_mismatches_mat" \
    --kid_id "$kid_id" \
    --dad_id "$dad_id" \
    --mom_id "$mom_id" \
    --vcf_joint_called "$vcf_joint_called" \
    --output_dir "$output_dir" \
    > ${output_dir}/expand_to_all_cpgs_trio.log 2>&1

log_info "Done expanding methylation to all CpG sites"
