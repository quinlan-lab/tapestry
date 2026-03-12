#!/bin/bash

# Usage:
#   Pedigree mode (default): ./aligned_bam_to_cpg_scores.sh
#   Trio mode:                ./aligned_bam_to_cpg_scores.sh -t <kid_id> <dad_id> <mom_id>
#   Trio mode with dev data:  ./aligned_bam_to_cpg_scores.sh -t <kid_id> <dad_id> <mom_id> --dev-dir dev-data
#   Example trio usage:       ./aligned_bam_to_cpg_scores.sh -t NA12878 NA12891 NA12892 --dev-dir dev-data

source src/util/logging.sh

MODE=""
DEV_DIR=""
TRIO_IDS=()

# --- Argument Parsing ---
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -t)
            MODE="trio"
            shift
            # Collect the three IDs
            if [[ "$#" -lt 3 ]]; then
                echo "Error: -t requires three arguments: <kid_id> <dad_id> <mom_id>"
                exit 1
            fi
            TRIO_IDS=("$1" "$2" "$3")
            shift 3
            ;;
        --dev-dir|-d)
            DEV_DIR="${2%/}"
            shift 2
            ;;
        *)
            echo "Error: Unknown parameter: $1"
            echo "Usage: $0 [-t <kid_id> <dad_id> <mom_id>] [--dev-dir <DEV_DATA_DIR>]"
            exit 1
            ;;
    esac
done

# --- Configuration ---
bin_dir="/uufs/chpc.utah.edu/common/HIPAA/u6018199/pb-CpG-tools-v3.0.0-x86_64-unknown-linux-gnu/bin"
export PATH=${bin_dir}:$PATH

if [ "$MODE" == "trio" ]; then
    # --- Trio Mode ---
    kid_id="${TRIO_IDS[0]}"
    dad_id="${TRIO_IDS[1]}"
    mom_id="${TRIO_IDS[2]}"

    if [ -n "$DEV_DIR" ]; then
        log_info "DEV MODE ENABLED: Reading from and writing to ${DEV_DIR}"
        bam_dir="${DEV_DIR}/output"
        base_output_dir="${DEV_DIR}/output/dna-methylation"

        if [ ! -f "${bam_dir}/${kid_id}.GRCh38.haplotagged.bam" ]; then
            echo "Error: Dev BAM not found: ${bam_dir}/${kid_id}.GRCh38.haplotagged.bam"
            echo "Did you run run-whatshap.sh --dev-dir ${DEV_DIR}?"
            exit 1
        fi
    else
        bam_dir="/scratch/ucgd/lustre-labs/quinlan/data-shared/whatshap-phasing"
        base_output_dir="/scratch/ucgd/lustre-labs/quinlan/data-shared/dna-methylation"
    fi

    for mode in "model" "count"; do
        output_dir="${base_output_dir}/CEPH1463.GRCh38.hifi.${mode}.trio-phased"
        log_info "Starting mode: '${mode}'"
        log_info "Output will be written to: ${output_dir}"
        mkdir -p ${output_dir}

        for id in ${kid_id} ${dad_id} ${mom_id}; do
            log_info "Submitting job for sample '${id}' in '${mode}' mode..."

            aligned_bam_to_cpg_scores \
                --bam "${bam_dir}/${id}.GRCh38.haplotagged.bam" \
                --output-prefix "${output_dir}/${id}.GRCh38.haplotagged" \
                --threads 8 \
                --min-coverage 10 \
                --min-mapq 1 \
                --pileup-mode "${mode}"
        done
    done

else
    # --- Pedigree Mode (original behavior) ---
    palladium_bam_dir="/scratch/ucgd/lustre-labs/quinlan/data-shared/read-backed-phasing"
    base_output_dir="/scratch/ucgd/lustre-labs/quinlan/data-shared/dna-methylation"

    echo "Getting Palladium prefixes..."
    prefixes=$(python src/util/get_palladium_prefixes.py)
    echo "... Got prefixes"

    for mode in "model" "count"; do
        output_dir="${base_output_dir}/CEPH1463.GRCh38.hifi.${mode}.read-backed-phased"
        echo "----------------------------------------------------"
        echo "Starting mode: '${mode}'"
        echo "Output will be written to: ${output_dir}"
        echo "----------------------------------------------------"

        mkdir -p ${output_dir}

        for prefix in $prefixes; do
            echo "Submitting job for sample '${prefix}' in '${mode}' mode..."

            nohup aligned_bam_to_cpg_scores \
                --bam "${palladium_bam_dir}/${prefix}.GRCh38.haplotagged.bam" \
                --output-prefix "${output_dir}/${prefix}.GRCh38.haplotagged" \
                --threads 8 \
                --min-coverage 10 \
                --min-mapq 1 \
                --pileup-mode "${mode}" \
                > /dev/null 2>&1 &
        done
    done
fi

log_info "All jobs have been submitted."
