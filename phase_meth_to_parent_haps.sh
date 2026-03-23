#!/bin/bash

# Usage:
#   Production:  ./phase_meth_to_parent_haps.sh
#   Dev mode:    ./phase_meth_to_parent_haps.sh --dev-dir trio_dev_data

source src/util/logging.sh

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
pedmec_phasing_dir="/scratch/ucgd/lustre-labs/quinlan/data-shared/pedmec-phasing" # output dir of run-whatshap.sh
meth_count_dir="/scratch/ucgd/lustre-labs/quinlan/data-shared/dna-methylation/CEPH1463.GRCh38.hifi.count.pedmec-phased" # output dir of aligned_bam_to_cpg_scores.sh in count mode
meth_model_dir="/scratch/ucgd/lustre-labs/quinlan/data-shared/dna-methylation/CEPH1463.GRCh38.hifi.model.pedmec-phased" # output dir of aligned_bam_to_cpg_scores.sh in model mode

# OUTPUT DIRS
output_dir="/scratch/ucgd/lustre-labs/quinlan/data-shared/dna-methylation/CEPH1463.GRCh38.hifi.parent-phased" # output dir of phase_meth_to_parent_haps.py

# --- Optional Dev Data Overrides ---
if [ -n "$DEV_DIR" ]; then
    log_info "DEV MODE ENABLED: Reading from and writing to ${DEV_DIR}"

    # INPUT DIRS
    pedmec_phasing_dir="${DEV_DIR}/output/pedmec-phasing" # output dir of run-whatshap.sh
    meth_dir="${DEV_DIR}/output/dna-methylation"
    meth_count_dir="${meth_dir}/CEPH1463.GRCh38.hifi.count.pedmec-phased" # output dir of aligned_bam_to_cpg_scores.sh in count mode
    meth_model_dir="${meth_dir}/CEPH1463.GRCh38.hifi.model.pedmec-phased" # output dir of aligned_bam_to_cpg_scores.sh in model mode

    # OUTPUT DIRS
    output_dir="${meth_dir}/CEPH1463.GRCh38.hifi.parent-phased" # output dir of phase_meth_to_parent_haps.py
fi

# INPUT FILES
vcf_pedmec_phased="${pedmec_phasing_dir}/CEPH-1463.joint.GRCh38.deepvariant.glnexus.phased.vcf.gz" # pedmec-phased trio VCF from run-whatshap.sh

blocks_tsv() {
    echo "${pedmec_phasing_dir}/CEPH-1463.joint.GRCh38.deepvariant.glnexus.phased.${1}.blocks.tsv" # pedmec-phasing blocks TSV from run-whatshap.sh
}

meth_bed() {
    # Usage: meth_bed <mode> <uid> <hap>
    local mode="$1" uid="$2" hap="$3"
    if [ "$mode" = "count" ]; then
        echo "${meth_count_dir}/${uid}.GRCh38.haplotagged.${hap}.bed.gz" # bed file of count-based methylation from aligned_bam_to_cpg_scores.sh
    else
        echo "${meth_model_dir}/${uid}.GRCh38.haplotagged.${hap}.bed.gz" # bed file of model-based methylation from aligned_bam_to_cpg_scores.sh
    fi
}

log_info "Phasing methylation to parent haplotypes for trio: ${kid_id}, ${dad_id}, ${mom_id}"

mkdir -p "${output_dir}"

PYTHONPATH=src:src/util .venv/bin/python src/phase_meth_to_parent_haps.py \
    --kid_id "$kid_id" \
    --dad_id "$dad_id" \
    --mom_id "$mom_id" \
    --vcf_pedmec_phased "$vcf_pedmec_phased" \
    --blocks_tsv_kid "$(blocks_tsv "$kid_id")" \
    --blocks_tsv_dad "$(blocks_tsv "$dad_id")" \
    --blocks_tsv_mom "$(blocks_tsv "$mom_id")" \
    --bed_meth_count_hap1_kid "$(meth_bed count "$kid_id" hap1)" \
    --bed_meth_count_hap2_kid "$(meth_bed count "$kid_id" hap2)" \
    --bed_meth_model_hap1_kid "$(meth_bed model "$kid_id" hap1)" \
    --bed_meth_model_hap2_kid "$(meth_bed model "$kid_id" hap2)" \
    --bed_meth_count_hap1_dad "$(meth_bed count "$dad_id" hap1)" \
    --bed_meth_count_hap2_dad "$(meth_bed count "$dad_id" hap2)" \
    --bed_meth_model_hap1_dad "$(meth_bed model "$dad_id" hap1)" \
    --bed_meth_model_hap2_dad "$(meth_bed model "$dad_id" hap2)" \
    --bed_meth_count_hap1_mom "$(meth_bed count "$mom_id" hap1)" \
    --bed_meth_count_hap2_mom "$(meth_bed count "$mom_id" hap2)" \
    --bed_meth_model_hap1_mom "$(meth_bed model "$mom_id" hap1)" \
    --bed_meth_model_hap2_mom "$(meth_bed model "$mom_id" hap2)" \
    --bed_meth_count_combined_kid "$(meth_bed count "$kid_id" combined)" \
    --bed_meth_model_combined_kid "$(meth_bed model "$kid_id" combined)" \
    --bed_meth_count_combined_dad "$(meth_bed count "$dad_id" combined)" \
    --bed_meth_model_combined_dad "$(meth_bed model "$dad_id" combined)" \
    --bed_meth_count_combined_mom "$(meth_bed count "$mom_id" combined)" \
    --bed_meth_model_combined_mom "$(meth_bed model "$mom_id" combined)" \
    --output_dir "$output_dir"

log_info "Done phasing methylation to parent haplotypes"
