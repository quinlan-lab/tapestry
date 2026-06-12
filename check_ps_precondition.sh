#!/bin/bash

# Verify the PS-tag precondition for the phase-set assignment fix in
# phasing_trio.py before regenerating hap-map output.
#
# Each SNV is assigned to its phase block by the per-sample PS (phase-set)
# FORMAT tag, which must equal phase_block_id in the WhatsHap blocks TSV. This
# driver runs src/util/check_ps_matches_blocks.py for the kid, dad, and mom
# against the same pedmec-phased VCF + blocks TSVs that phase_meth_to_parent_haps.sh
# consumes. If any sample fails, the fix would silently fall back to span-based
# assignment, so do NOT trust regenerated output until all three pass.
#
# Usage:
#   Production:  ./check_ps_precondition.sh
#   Dev mode:    ./check_ps_precondition.sh --dev-dir trio_dev_data

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
pedmec_phasing_dir="/scratch/ucgd/lustre-labs/quinlan/data-shared/pedmec-phasing" # output dir of run-whatshap.sh

# --- Optional Dev Data Overrides ---
if [ -n "$DEV_DIR" ]; then
    log_info "DEV MODE ENABLED: Reading from ${DEV_DIR}"
    pedmec_phasing_dir="${DEV_DIR}/output/pedmec-phasing" # output dir of run-whatshap.sh
fi

# INPUT FILES
vcf_pedmec_phased="${pedmec_phasing_dir}/CEPH-1463.joint.GRCh38.deepvariant.glnexus.phased.vcf.gz" # pedmec-phased trio VCF from run-whatshap.sh

blocks_tsv() {
    echo "${pedmec_phasing_dir}/CEPH-1463.joint.GRCh38.deepvariant.glnexus.phased.${1}.blocks.tsv" # pedmec-phasing blocks TSV from run-whatshap.sh
}

log_info "Verifying PS precondition for trio: ${kid_id}, ${dad_id}, ${mom_id}"

failed=0
for uid in "$kid_id" "$dad_id" "$mom_id"; do
    log_info "=== ${uid} ==="
    PYTHONPATH=src:src/util .venv/bin/python src/util/check_ps_matches_blocks.py \
        --vcf "$vcf_pedmec_phased" \
        --blocks_tsv "$(blocks_tsv "$uid")" \
        --sample "$uid" || failed=1
done

if [ "$failed" -ne 0 ]; then
    log_info "PS precondition FAILED for at least one sample; do not trust regenerated output."
    exit 1
fi
log_info "PS precondition holds for all samples; safe to regenerate hap-map output."
