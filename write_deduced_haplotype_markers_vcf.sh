#!/bin/bash

# Write a VCF of the markers at which a kid has >=1 deduced haplotype, from the
# CEPH1463 haplotype-map markers table, for display as an IGV track. The output
# filename is derived from the sample ID(s) by the Python script and written
# next to the markers file.
#
# Usage:
#   ./write_deduced_haplotype_markers_vcf.sh                  # default kid
#   ./write_deduced_haplotype_markers_vcf.sh NA12879 NA12881  # one or more kids

source src/util/logging.sh

# --- Configuration ---
markers="/scratch/ucgd/lustre-labs/quinlan/data-shared/haplotype-maps/CEPH1463.GRCh38/CEPH1463.GRCh38.markers.sorted.txt"

# Kid sample ID(s): override by passing them as positional arguments.
if [ "$#" -gt 0 ]; then
    samples=("$@")
else
    samples=("NA12879")
fi

# Build the repeated --sample flags the Python script expects.
sample_flags=()
for s in "${samples[@]}"; do
    sample_flags+=(--sample "$s")
done

log_info "Selecting deduced-haplotype markers for: ${samples[*]}"

PYTHONPATH=src:src/util .venv/bin/python src/util/write_deduced_haplotype_markers_vcf.py \
    --markers "$markers" \
    "${sample_flags[@]}"

log_info "Done: write_deduced_haplotype_markers_vcf.sh"
