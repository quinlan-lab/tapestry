#!/bin/bash
set -euo pipefail

usage() { echo "Usage: $0 <DEV_DATA_OUTPUT_DIR>"; }
if [[ $# == 1 && ( $1 == -h || $1 == --help ) ]]; then usage; exit 0; fi
if [[ $# != 1 || -z $1 || $1 == -* || $1 == / ]]; then usage >&2; exit 1; fi
for tool in curl gunzip samtools; do
    command -v "$tool" >/dev/null || { echo "Missing prerequisite: $tool" >&2; exit 1; }
done
source "$(dirname "$0")/trio_dev_data_config.sh"
DEV_DIR="${1%/}"
mkdir -p "$DEV_DIR/input"
dev_chrom="${DEV_REGION%%:*}"
dev_ref="$DEV_DIR/input/$dev_chrom.fa"
staging=$(mktemp -d "$DEV_DIR/input/.reference.XXXXXX")
trap 'rm -rf -- "$staging"' EXIT
url="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/${dev_chrom}.fa.gz"
echo "Downloading full hg38 $dev_chrom reference (coordinates preserved)..."
curl --fail --location --show-error --silent --retry 3 "$url" | gunzip > "$staging/$dev_chrom.fa"
samtools faidx "$staging/$dev_chrom.fa"
# Expected GRCh38 chr22 length; reject a truncated/wrong reference.
awk -v chrom="$dev_chrom" '$1 == chrom && $2 == 50818468 {ok=1} END {exit !ok}' "$staging/$dev_chrom.fa.fai"
mv "$staging/$dev_chrom.fa" "$dev_ref"
mv "$staging/$dev_chrom.fa.fai" "$dev_ref.fai"
printf '%s\n' "$url" > "$dev_ref.source.txt"
echo "Reference saved to $dev_ref"
