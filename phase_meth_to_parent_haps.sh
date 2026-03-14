#!/bin/bash

# Usage:
#   ./phase_meth_to_parent_haps.sh --dev-dir trio-dev-data

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

if [ -n "$DEV_DIR" ]; then
    log_info "DEV MODE ENABLED: Reading from and writing to ${DEV_DIR}"

    ped="${DEV_DIR}/input/trio.ped"
    vcf_trio_phased="${DEV_DIR}/output/trio-phasing/CEPH-1463.joint.GRCh38.deepvariant.glnexus.phased.vcf.gz"
    output_dir="${DEV_DIR}/output/dna-methylation/CEPH1463.GRCh38.hifi.parent-phased"

    trio_phasing_dir="${DEV_DIR}/output/trio-phasing"
    blocks_tsv_kid="${trio_phasing_dir}/CEPH-1463.joint.GRCh38.deepvariant.glnexus.phased.${kid_id}.blocks.tsv"
    blocks_tsv_dad="${trio_phasing_dir}/CEPH-1463.joint.GRCh38.deepvariant.glnexus.phased.${dad_id}.blocks.tsv"
    blocks_tsv_mom="${trio_phasing_dir}/CEPH-1463.joint.GRCh38.deepvariant.glnexus.phased.${mom_id}.blocks.tsv"

    meth_count_dir="${DEV_DIR}/output/dna-methylation/CEPH1463.GRCh38.hifi.count.read-backed-phased"
    meth_model_dir="${DEV_DIR}/output/dna-methylation/CEPH1463.GRCh38.hifi.model.read-backed-phased"
else
    ped="XXX" # TODO: fill this in for production
    vcf_trio_phased="/scratch/ucgd/lustre-labs/quinlan/data-shared/whatshap-phasing/CEPH-1463.joint.GRCh38.deepvariant.glnexus.phased.vcf.gz"
    output_dir="/scratch/ucgd/lustre-labs/quinlan/data-shared/dna-methylation/CEPH1463.GRCh38.hifi.parent-phased"

    trio_phasing_dir="/scratch/ucgd/lustre-labs/quinlan/data-shared/whatshap-phasing"
    blocks_tsv_kid="${trio_phasing_dir}/CEPH-1463.joint.GRCh38.deepvariant.glnexus.phased.${kid_id}.blocks.tsv"
    blocks_tsv_dad="${trio_phasing_dir}/CEPH-1463.joint.GRCh38.deepvariant.glnexus.phased.${dad_id}.blocks.tsv"
    blocks_tsv_mom="${trio_phasing_dir}/CEPH-1463.joint.GRCh38.deepvariant.glnexus.phased.${mom_id}.blocks.tsv"

    meth_count_dir="/scratch/ucgd/lustre-labs/quinlan/data-shared/dna-methylation/CEPH1463.GRCh38.hifi.count.read-backed-phased"
    meth_model_dir="/scratch/ucgd/lustre-labs/quinlan/data-shared/dna-methylation/CEPH1463.GRCh38.hifi.model.read-backed-phased"
fi

meth_bed() {
    local uid=$1
    local mode=$2
    local hap=$3
    if [ "$mode" == "count" ]; then
        echo "${meth_count_dir}/${uid}.GRCh38.haplotagged.${hap}.bed.gz"
    else
        echo "${meth_model_dir}/${uid}.GRCh38.haplotagged.${hap}.bed.gz"
    fi
}

log_info "Phasing methylation to parent haplotypes for trio: ${kid_id}, ${dad_id}, ${mom_id}"

mkdir -p ${output_dir}

python src/phase_meth_to_parent_haps.py \
    --kid_id ${kid_id} \
    --dad_id ${dad_id} \
    --mom_id ${mom_id} \
    --ped ${ped} \
    --vcf_trio_phased ${vcf_trio_phased} \
    --blocks_tsv_kid ${blocks_tsv_kid} \
    --blocks_tsv_dad ${blocks_tsv_dad} \
    --blocks_tsv_mom ${blocks_tsv_mom} \
    --bed_meth_count_hap1_kid $(meth_bed ${kid_id} count hap1) \
    --bed_meth_count_hap2_kid $(meth_bed ${kid_id} count hap2) \
    --bed_meth_model_hap1_kid $(meth_bed ${kid_id} model hap1) \
    --bed_meth_model_hap2_kid $(meth_bed ${kid_id} model hap2) \
    --bed_meth_count_hap1_dad $(meth_bed ${dad_id} count hap1) \
    --bed_meth_count_hap2_dad $(meth_bed ${dad_id} count hap2) \
    --bed_meth_model_hap1_dad $(meth_bed ${dad_id} model hap1) \
    --bed_meth_model_hap2_dad $(meth_bed ${dad_id} model hap2) \
    --bed_meth_count_hap1_mom $(meth_bed ${mom_id} count hap1) \
    --bed_meth_count_hap2_mom $(meth_bed ${mom_id} count hap2) \
    --bed_meth_model_hap1_mom $(meth_bed ${mom_id} model hap1) \
    --bed_meth_model_hap2_mom $(meth_bed ${mom_id} model hap2) \
    --output_dir ${output_dir}

log_info "Done phasing methylation to parent haplotypes"
