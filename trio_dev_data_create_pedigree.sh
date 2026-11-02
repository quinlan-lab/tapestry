#!/bin/bash
set -euo pipefail

usage() {
    echo "Usage: $0 <DEV_DATA_DIR>"
    echo "Extract matching original pedigree outputs on CHPC into input/pedigree/."
    echo "Overrides: TAPESTRY_READ_PHASED_DIR, TAPESTRY_INHERITANCE_DIR."
}
if [[ $# == 1 && ( $1 == -h || $1 == --help ) ]]; then usage; exit 0; fi
if [[ $# != 1 || -z $1 || $1 == -* || $1 == / ]]; then usage >&2; exit 1; fi
repo_dir=$(cd "$(dirname "$0")" && pwd)
source "$repo_dir/trio_dev_data_config.sh"
source "$repo_dir/src/util/logging.sh"
source "$repo_dir/src/util/trio_ped.sh"
dev_dir="${1%/}"
read_trio_ids "$dev_dir/input/trio.ped"
[[ "$kid_id,$dad_id,$mom_id" == "$DEV_KID_ID,$DEV_DAD_ID,$DEV_MOM_ID" ]] || {
    echo "Fixture PED does not match trio_dev_data_config.sh" >&2; exit 1;
}
[[ "$DEV_REGION" =~ ^([^:]+):([0-9]+)-([0-9]+)$ ]] || {
    echo "Expected a bounded region chrom:start-end in trio_dev_data_config.sh" >&2; exit 1;
}
chrom="${BASH_REMATCH[1]}"; first="${BASH_REMATCH[2]}"; last="${BASH_REMATCH[3]}"
(( first >= 1 && last >= first )) || { echo "Invalid fixture coordinates" >&2; exit 1; }
read_dir="${TAPESTRY_READ_PHASED_DIR:-/scratch/ucgd/lustre-labs/quinlan/data-shared/read-backed-phasing}"
iht_dir="${TAPESTRY_INHERITANCE_DIR:-/scratch/ucgd/lustre-labs/quinlan/data-shared/haplotype-maps/CEPH1463.GRCh38}"
prefix="CEPH1463.GRCh38"
for tool in bcftools awk; do
    command -v "$tool" >/dev/null || { echo "Missing prerequisite: $tool" >&2; exit 1; }
done
sources=("$iht_dir/$prefix.pass.sorted.vcf.gz" "$iht_dir/$prefix.iht.sorted.txt")
for id in "$kid_id" "$dad_id" "$mom_id"; do
    sources+=("$read_dir/$id.GRCh38.deepvariant.glnexus.phased.vcf.gz" "$read_dir/$id.GRCh38.hiphase.blocks.tsv")
done
for file in "${sources[@]}"; do
    [[ -r "$file" ]] || { echo "Missing source: $file" >&2; exit 1; }
done
staging=$(mktemp -d "$dev_dir/input/.pedigree.XXXXXX")
trap 'rm -rf -- "$staging"' EXIT
for id in "$kid_id" "$dad_id" "$mom_id"; do
    echo "Extracting original HiPhase outputs: $id $DEV_REGION"
    vcf="$staging/$id.GRCh38.deepvariant.glnexus.phased.vcf.gz"
    bcftools view -r "$DEV_REGION" -s "$id" "$read_dir/$id.GRCh38.deepvariant.glnexus.phased.vcf.gz" -Oz -o "$vcf"
    bcftools index -t "$vcf"
    # HiPhase block coordinates are one-based inclusive. Preserve original spans
    # and phase-set IDs; only select blocks overlapping the fixture interval.
    awk -F '\t' -v OFS='\t' -v chrom="$chrom" -v first="$first" -v last="$last" -v sample="$id" '
        NR == 1 {
            for (i=1;i<=NF;i++) col[$i]=i
            if (!("chrom" in col) || !("start" in col) || !("end" in col) || !("sample_name" in col)) exit 2
            print; next
        }
        $(col["chrom"]) == chrom && $(col["sample_name"]) == sample && $(col["start"]) <= last && $(col["end"]) >= first {print}
    ' "$read_dir/$id.GRCh38.hiphase.blocks.tsv" > "$staging/$id.GRCh38.hiphase.blocks.tsv"
done
# Subset AFTER full-pedigree inference. Never rerun gtg on the reduced trio.
vcf="$staging/$prefix.pass.sorted.vcf.gz"
bcftools view -r "$DEV_REGION" -s "$kid_id,$dad_id,$mom_id" "$iht_dir/$prefix.pass.sorted.vcf.gz" -Oz -o "$vcf"
bcftools index -t "$vcf"
# IHT coordinates are zero-based half-open. Keep original spans/labels but
# publish only this trio's columns. Omit the unused full-family markers field.
awk -v OFS=' ' -v chrom="$chrom" -v first="$first" -v last="$last" -v kid="$kid_id" -v dad="$dad_id" -v mom="$mom_id" '
    NR == 1 {
        for (i=1;i<=NF;i++) col[$i]=i
        if (!(kid in col) || !(dad in col) || !(mom in col) || !("marker_count" in col) || !("len" in col)) exit 2
        print "#chrom", "start", "end", kid, dad, mom, "marker_count", "len", "markers"; next
    }
    $1 == chrom && $2 < last && $3 > first-1 {
        print $1,$2,$3,$(col[kid]),$(col[dad]),$(col[mom]),$(col["marker_count"]),$(col["len"]),"."
    }
' "$iht_dir/$prefix.iht.sorted.txt" > "$staging/$prefix.iht.sorted.txt"
printf 'region\t%s\nchild\t%s\nfather\t%s\nmother\t%s\nread_phased_source\t%s\ninheritance_source\t%s\n' \
    "$DEV_REGION" "$kid_id" "$dad_id" "$mom_id" "$read_dir" "$iht_dir" > "$staging/fixture.tsv"
if [[ -e "$dev_dir/input/pedigree" ]]; then
    mkdir -p "$dev_dir/archive"
    archive=$(mktemp -d "$dev_dir/archive/pedigree-$(date -u +%Y%m%dT%H%M%SZ).XXXXXX")
    mv "$dev_dir/input/pedigree" "$archive/input-pedigree"
    echo "Previous pedigree subset archived: $archive"
fi
mv "$staging" "$dev_dir/input/pedigree"
echo "Pedigree subset ready: $dev_dir/input/pedigree"
echo "Transfer with: git add '$dev_dir/input/pedigree' && git commit -m 'Add matching chr22 pedigree inputs' && git push"
