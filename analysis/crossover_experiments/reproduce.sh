#!/bin/bash
# Reproduce the chr22 slide-deck experiments from the checked-in trio fixture.
set -euo pipefail
cd "$(dirname "$0")/../.."

usage() {
    echo "Usage: $0 UCSC_GENOMIC_SUPER_DUPS.txt.gz ENCODE_HG38_BLACKLIST.bed.gz [OUTPUT_DIR]"
}
if [[ ${1:-} == -h || ${1:-} == --help ]]; then usage; exit 0; fi
if [[ $# -lt 2 || $# -gt 3 ]]; then usage >&2; exit 1; fi
for tool in bcftools bgzip gzip; do
    command -v "$tool" >/dev/null || { echo "Missing prerequisite: $tool" >&2; exit 1; }
done

python=.venv/bin/python
whatshap=.venv/bin/whatshap
[[ "$("$whatshap" --version)" == 2.8 ]] || { echo "Install WhatsHap 2.8 in .venv" >&2; exit 1; }
export PYTHONPATH=src:src/util
source trio_dev_data_config.sh
source src/util/logging.sh
source src/util/trio_ped.sh

input=trio_dev_data/input
out=${3:-trio_dev_data/output/crossover-experiments}
scripts=analysis/crossover_experiments
ped=$input/trio.ped
read_trio_ids "$ped"
chrom=${DEV_REGION%%:*}
reference=$input/$chrom.fa
joint=$input/CEPH-1463.joint.GRCh38.deepvariant.glnexus.vcf.gz
bams=("$input/$kid_id.GRCh38.haplotagged.bam" "$input/$dad_id.GRCh38.haplotagged.bam" "$input/$mom_id.GRCh38.haplotagged.bam")
for file in "$ped" "$reference" "$reference.fai" "$joint" "$1" "$2" "${bams[@]}" "${bams[@]/%/.bai}"; do
    [[ -s "$file" ]] || { echo "Missing input: $file" >&2; exit 1; }
done
if [[ -e "$out" ]]; then
    echo "Output already exists. Choose a new output directory: $out" >&2
    exit 1
fi
mkdir -p "$out"

# Extract the public annotations on the fixture chromosome. UCSC table has a bin column.
gzip -dc "$1" | awk -v OFS='\t' -v chrom="$chrom" '$2 == chrom {print $2,$3,$4}' > "$out/segdup.bed"
gzip -dc "$2" | awk -v OFS='\t' -v chrom="$chrom" '$1 == chrom {print $1,$2,$3}' > "$out/blacklist.bed"
for bed in "$out/segdup.bed" "$out/blacklist.bed"; do
    [[ -s "$bed" ]] || { echo "No $chrom annotations in $bed" >&2; exit 1; }
done
sort -k1,1 -k2,2n -k3,3n -u "$out/segdup.bed" "$out/blacklist.bed" > "$out/public-mask.bed"

quality="$out/quality-input.vcf.gz"
"$python" src/prefilter_whatshap_vcf.py --input "$joint" --output "$quality" \
    --samples "$kid_id" "$dad_id" "$mom_id"

# Mask case 1 uses independent public tracks. Case 2 is a targeted perturbation
# chosen after inspecting the reads, using exactly four one-based VCF positions.
bcftools view -T "^$out/public-mask.bed" --targets-overlap 0 -Oz -o "$out/public-mask-input.vcf.gz" "$quality"
bcftools view -e 'CHROM="chr22" && (POS=20503731 || POS=20503741 || POS=20503743 || POS=20503755)' \
    -Oz -o "$out/four-snp-mask-input.vcf.gz" "$quality"
# Test only Recombinator's depth floor on top of the Tapestry quality filter.
bcftools view -i 'N_PASS(FMT/DP>=18)=3' -Oz -o "$out/depth18-input.vcf.gz" "$quality"

phase() {
    local name=$1 vcf=$2 dest
    dest="$out/$name"
    mkdir -p "$dest"
    printf 'WhatsHap %s\n' "$("$whatshap" --version)" > "$dest/versions.txt"
    bcftools --version >> "$dest/versions.txt"
    bcftools view -H "$vcf" | wc -l > "$dest/input-records.txt"
    "$whatshap" unphase "$vcf" | bgzip -c > "$dest/unphased.vcf.gz"
    local args=(phase --ped "$ped" --reference "$reference" --chromosome "$chrom"
        --sample "$kid_id" --sample "$dad_id" --sample "$mom_id"
        --recombination-list "$dest/recombination.txt" --output "$dest/phased.vcf.gz"
        "$dest/unphased.vcf.gz" "${bams[@]}")
    printf '%q ' "$whatshap" "${args[@]}" > "$dest/phase.command.txt"
    printf '\n' >> "$dest/phase.command.txt"
    log_info "Phasing $name ($DEV_REGION)"
    "$whatshap" "${args[@]}" > "$dest/phase.log" 2>&1
    bcftools index -t "$dest/phased.vcf.gz"
    "$python" "$scripts/summarize.py" "$dest" "$out/segdup.bed" "$out/blacklist.bed"
}

for branch in public-mask four-snp-mask depth18; do
    bcftools index -t "$out/$branch-input.vcf.gz"
done

phase baseline "$joint"
phase quality "$quality"
phase public-mask "$out/public-mask-input.vcf.gz"
phase four-snp-mask "$out/four-snp-mask-input.vcf.gz"
phase depth18 "$out/depth18-input.vcf.gz"

"$python" src/filter_whatshap_crossovers_by_reads.py \
    --crossovers "$out/quality/crossovers.tsv" --phased-vcf "$out/quality/phased.vcf.gz" \
    --ped "$ped" --bam "$kid_id" "${bams[0]}" --bam "$dad_id" "${bams[1]}" --bam "$mom_id" "${bams[2]}" \
    --output-dir "$out/quality/read-evidence"

log_info "Completed chr22 crossover experiments: $out"
