#!/usr/bin/env bash
set -euo pipefail

KID_ID="NA12878"
DAD_ID="NA12891"
MOM_ID="NA12892"

TRIO_PHASING_DIR="trio-dev-data/output/trio-phasing"
METH_DIR="trio-dev-data/output/dna-methylation"
COUNT_DIR="${METH_DIR}/CEPH1463.GRCh38.hifi.count.pedmec-phased"
MODEL_DIR="${METH_DIR}/CEPH1463.GRCh38.hifi.model.pedmec-phased"
OUTPUT_DIR="${METH_DIR}/CEPH1463.GRCh38.hifi.parent-phased"

VCF_TRIO_PHASED="${TRIO_PHASING_DIR}/CEPH-1463.joint.GRCh38.deepvariant.glnexus.phased.vcf.gz"

blocks_tsv() {
    echo "${TRIO_PHASING_DIR}/CEPH-1463.joint.GRCh38.deepvariant.glnexus.phased.${1}.blocks.tsv"
}

meth_bed() {
    # Usage: meth_bed <mode> <uid> <hap>
    local mode="$1" uid="$2" hap="$3"
    if [ "$mode" = "count" ]; then
        echo "${COUNT_DIR}/${uid}.GRCh38.haplotagged.${hap}.bed.gz"
    else
        echo "${MODEL_DIR}/${uid}.GRCh38.haplotagged.${hap}.bed.gz"
    fi
}

PYTHONPATH=src:src/util .venv/bin/python src/phase_meth_to_parent_haps.py \
    --kid_id "$KID_ID" \
    --dad_id "$DAD_ID" \
    --mom_id "$MOM_ID" \
    --vcf_trio_phased "$VCF_TRIO_PHASED" \
    --blocks_tsv_kid "$(blocks_tsv "$KID_ID")" \
    --blocks_tsv_dad "$(blocks_tsv "$DAD_ID")" \
    --blocks_tsv_mom "$(blocks_tsv "$MOM_ID")" \
    --bed_meth_count_hap1_kid "$(meth_bed count "$KID_ID" hap1)" \
    --bed_meth_count_hap2_kid "$(meth_bed count "$KID_ID" hap2)" \
    --bed_meth_model_hap1_kid "$(meth_bed model "$KID_ID" hap1)" \
    --bed_meth_model_hap2_kid "$(meth_bed model "$KID_ID" hap2)" \
    --bed_meth_count_hap1_dad "$(meth_bed count "$DAD_ID" hap1)" \
    --bed_meth_count_hap2_dad "$(meth_bed count "$DAD_ID" hap2)" \
    --bed_meth_model_hap1_dad "$(meth_bed model "$DAD_ID" hap1)" \
    --bed_meth_model_hap2_dad "$(meth_bed model "$DAD_ID" hap2)" \
    --bed_meth_count_hap1_mom "$(meth_bed count "$MOM_ID" hap1)" \
    --bed_meth_count_hap2_mom "$(meth_bed count "$MOM_ID" hap2)" \
    --bed_meth_model_hap1_mom "$(meth_bed model "$MOM_ID" hap1)" \
    --bed_meth_model_hap2_mom "$(meth_bed model "$MOM_ID" hap2)" \
    --output_dir "$OUTPUT_DIR"
