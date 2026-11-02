#!/bin/bash
set -euo pipefail

usage() {
    echo "Usage: $0 <DEV_DATA_OUTPUT_DIR>"
    echo "Subset the configured chr22 trio on CHPC; archive previous inputs/outputs."
    echo "Optional source overrides: TAPESTRY_SOURCE_VCF, TAPESTRY_SOURCE_BAM_DIR."
}
if [[ $# == 1 && ( $1 == -h || $1 == --help ) ]]; then usage; exit 0; fi
if [[ $# != 1 || -z $1 || $1 == -* || $1 == / ]]; then usage >&2; exit 1; fi
source "$(dirname "$0")/trio_dev_data_config.sh"
DEV_DIR="${1%/}"
vcf_joint_called="${TAPESTRY_SOURCE_VCF:-/scratch/ucgd/lustre-labs/quinlan/data-shared/datasets/Palladium/deepvariant/CEPH-1463.joint.GRCh38.deepvariant.glnexus.phased.vcf.gz}"
palladium_bam_dir="${TAPESTRY_SOURCE_BAM_DIR:-/scratch/ucgd/lustre-labs/quinlan/data-shared/datasets/Palladium/hifi-bams/GRCh38}"
for tool in bcftools samtools; do
    command -v "$tool" >/dev/null || { echo "Missing prerequisite: $tool" >&2; exit 1; }
done
for file in "$vcf_joint_called" "$palladium_bam_dir/$DEV_KID_ID.GRCh38.haplotagged.bam" "$palladium_bam_dir/$DEV_DAD_ID.GRCh38.haplotagged.bam" "$palladium_bam_dir/$DEV_MOM_ID.GRCh38.haplotagged.bam"; do
    [[ -r "$file" ]] || { echo "Missing source: $file (run on CHPC or set source overrides)" >&2; exit 1; }
done
# Build off to the side: failed extraction cannot publish a partial fixture.
mkdir -p "$DEV_DIR"
staging=$(mktemp -d "${DEV_DIR}/.input.XXXXXX")
trap 'rm -rf -- "$staging"' EXIT
cat > "$staging/trio.ped" <<'PED_COMMENTS'
# Public-sharing provenance (Platinum Pedigree/Palladium release):
# NA12879, NA12877, and NA12878 are explicitly listed as consented for public access:
# https://pmc.ncbi.nlm.nih.gov/articles/PMC11326147/ (Data availability)
# Consortium access guidance (distinguishes open and controlled samples):
# https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Datasets#accessing-controlled-samples
# Public dataset license: CC BY 4.0; attribute the source and identify subsetting/modifications.
# https://registry.opendata.aws/platinum-pedigree/
# Applies to the public release; exact CHPC source-file provenance is not yet verified.
# Subset to these three individuals before sharing; do not publish the full source family VCF.
# Sex follows the consortium sample metadata: NA12879 = 2 (female).
# family child father mother sex phenotype
PED_COMMENTS
printf '1463\t%s\t%s\t%s\t2\t0\n' \
    "$DEV_KID_ID" "$DEV_DAD_ID" "$DEV_MOM_ID" >> "$staging/trio.ped"
dev_vcf="$staging/CEPH-1463.joint.GRCh38.deepvariant.glnexus.vcf.gz"
echo "Extracting $DEV_REGION for $DEV_KID_ID / $DEV_DAD_ID / $DEV_MOM_ID"
bcftools view -r "$DEV_REGION" -s "$DEV_KID_ID,$DEV_DAD_ID,$DEV_MOM_ID" \
    "$vcf_joint_called" -Oz -o "$dev_vcf"
bcftools index "$dev_vcf"
for id in "$DEV_KID_ID" "$DEV_DAD_ID" "$DEV_MOM_ID"; do
    in_bam="$palladium_bam_dir/$id.GRCh38.haplotagged.bam"
    out_bam="$staging/$id.GRCh38.haplotagged.bam"
    echo "Subsetting $in_bam"
    samtools view -b -o "$out_bam" "$in_bam" "$DEV_REGION"
    samtools quickcheck "$out_bam"
    samtools index "$out_bam"
done
printf 'region\t%s\nchild\t%s\nfather\t%s\nmother\t%s\nsource_vcf\t%s\nsource_bam_dir\t%s\n' \
    "$DEV_REGION" "$DEV_KID_ID" "$DEV_DAD_ID" "$DEV_MOM_ID" "$vcf_joint_called" "$palladium_bam_dir" > "$staging/fixture.tsv"
# Preserve the complete old run together, including its reference and PED.
if [[ -e "$DEV_DIR/input" || -e "$DEV_DIR/output" ]]; then
    mkdir -p "$DEV_DIR/archive"
    archive=$(mktemp -d "$DEV_DIR/archive/replaced-$(date -u +%Y%m%dT%H%M%SZ).XXXXXX")
    for name in input output; do
        if [[ -e "$DEV_DIR/$name" ]]; then mv "$DEV_DIR/$name" "$archive/$name"; fi
    done
    echo "Previous fixture archived: $archive"
fi
mv "$staging" "$DEV_DIR/input"
echo "Fixture ready. Run ./trio_dev_data_get_ref.sh '$DEV_DIR', then ./run-whatshap.sh --dev-dir '$DEV_DIR' --crossover-distance-bp 0"
