#!/bin/bash

# Download the full chromosome reference FASTA for the dev region.
# This is separate from trio_dev_data_create.sh because it requires downloading
# the entire chromosome (~250MB for chr1) to preserve correct coordinates,
# and only needs to be run once.

if [ "$#" -ne 1 ]; then
    echo "Error: Target directory required."
    echo "Usage: $0 <DEV_DATA_OUTPUT_DIR>"
    exit 1
fi

DEV_DIR="${1%/}"
mkdir -p "${DEV_DIR}/input"

source "$(dirname "$0")/dev_data_config.sh"

dev_ref="${DEV_DIR}/input/dev_reference.fa"
dev_chrom="${DEV_REGION%%:*}"

echo "Downloading reference for ${dev_chrom} from UCSC..."
curl -sS "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/${dev_chrom}.fa.gz" | gunzip > "${dev_ref}"
samtools faidx "${dev_ref}"

echo "Reference saved to ${dev_ref}"
