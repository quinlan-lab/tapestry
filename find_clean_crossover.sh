#!/bin/bash

# Rank candidate meiotic crossovers in the kid by single-crossover cleanliness,
# from the hap-map blocks and bit-vector mismatch sites produced by
# phase_meth_to_parent_haps.sh. Writes a ranked TSV plus IGV navigation/batch
# outputs, and logs the cleanest single-parent candidate to confirm in IGV
# before recutting the dev set around it.
#
# Usage:
#   Production:  ./find_clean_crossover.sh
#   Dev mode:    ./find_clean_crossover.sh --dev-dir trio_dev_data

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

# --- Default Configurations (Production) ---
# INPUT/OUTPUT DIR: the parent-phased dir written by phase_meth_to_parent_haps.py
output_dir="/scratch/ucgd/lustre-labs/quinlan/data-shared/dna-methylation/CEPH1463.GRCh38.hifi.parent-phased"

# --- Optional Dev Data Overrides ---
if [ -n "$DEV_DIR" ]; then
    log_info "DEV MODE ENABLED: Reading from ${DEV_DIR}"
    output_dir="${DEV_DIR}/output/dna-methylation/CEPH1463.GRCh38.hifi.parent-phased"
fi

# INPUT FILES (outputs of phase_meth_to_parent_haps.py)
blocks_paternal="${output_dir}/${kid_id}.hap-map-blocks.paternal.sorted.bed.gz"
blocks_maternal="${output_dir}/${kid_id}.hap-map-blocks.maternal.sorted.bed.gz"
mismatch_paternal="${output_dir}/${kid_id}.bit-vector-sites-mismatches.paternal.bed"
mismatch_maternal="${output_dir}/${kid_id}.bit-vector-sites-mismatches.maternal.bed"

# OUTPUT FILES
out_tsv="${output_dir}/${kid_id}.crossover-candidates.tsv"
out_igv_bed="${output_dir}/${kid_id}.crossover-candidates.igv.bed"
out_igv_batch="${output_dir}/${kid_id}.crossover-candidates.igv.bat"

log_info "Ranking candidate crossovers for ${kid_id} from ${output_dir}"

PYTHONPATH=src:src/util .venv/bin/python src/util/find_clean_crossover.py \
    --paternal-blocks   "$blocks_paternal" \
    --maternal-blocks   "$blocks_maternal" \
    --paternal-mismatch "$mismatch_paternal" \
    --maternal-mismatch "$mismatch_maternal" \
    --out       "$out_tsv" \
    --igv-bed   "$out_igv_bed" \
    --igv-batch "$out_igv_batch"

log_info "Done; ranked candidates in ${out_tsv}"
