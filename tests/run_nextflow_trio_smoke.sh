#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
nextflow_bin="${NEXTFLOW:-nextflow}"
export PATH="$repo_root/.venv/bin:$PATH"
export NXF_OFFLINE=true

usage() {
    echo "Usage: $0"
    echo "Check local Nextflow trio phasing and crossover screening on the chr22 fixture."
}
if [[ ${1:-} == -h || ${1:-} == --help ]]; then
    usage
    exit 0
fi
if [[ $# -ne 0 ]]; then
    usage >&2
    exit 2
fi
for tool in "$nextflow_bin" bcftools bgzip python3; do
    command -v "$tool" >/dev/null || { echo "Missing prerequisite: $tool" >&2; exit 1; }
done
if [[ "$(python3 -m whatshap --version)" != 2.8 ]]; then
    echo "The fixture check requires WhatsHap 2.8" >&2
    exit 1
fi

fixture="$repo_root/trio_dev_data/input"
if [[ ! -s "$fixture/chr22.fa" || ! -s "$fixture/chr22.fa.fai" ]]; then
    echo "Prepare the chr22 reference with ./trio_dev_data_get_ref.sh trio_dev_data" >&2
    exit 1
fi
mkdir -p "$repo_root/.test-work"
root="$(mktemp -d "$repo_root/.test-work/trio-smoke.XXXXXX")"
cleanup() {
    if [[ "$?" -eq 0 ]]; then
        rm -rf "$root"
    else
        echo "Trio smoke test diagnostics retained in $root" >&2
    fi
}
trap cleanup EXIT
cd "$root"

common=(
    run "$repo_root" -entry trio -ansi-log false
    --ped "$fixture/trio.ped"
    --joint-vcf "$fixture/CEPH-1463.joint.GRCh38.deepvariant.glnexus.vcf.gz"
    --reference-fasta "$fixture/chr22.fa"
    --child-bam "$fixture/NA12879.GRCh38.haplotagged.bam"
    --father-bam "$fixture/NA12877.GRCh38.haplotagged.bam"
    --mother-bam "$fixture/NA12878.GRCh38.haplotagged.bam"
    --chromosomes chr22 --region chr22:20500000-21500000
    --outdir results
)
"$nextflow_bin" "${common[@]}"

results=results/per_chrom_phased
vcf="$results/chr22.phased.vcf.gz"
test -s "$vcf.tbi"
diff -u <(printf '%s\n' NA12879 NA12877 NA12878) <(bcftools query -l "$vcf")
bcftools query -r chr22:20500000-21500000 -f '%CHROM\t%POS[\t%GT]\n' "$vcf" > records.tsv
awk '
    $1 != "chr22" || $2 < 20500000 || $2 > 21500000 { exit 1 }
    $3 ~ /[|]/ { phased = 1 }
    END { if (NR != 2303 || !phased) exit 1 }
' records.tsv
[[ "$(bcftools index -n "$vcf")" == 2303 ]]

# Compare the raw candidate brackets and transmitted states from the quality-filtered experiment.
awk '!/^#/ { print $1, $2, $3, $4, $5, $6, $7, $8 }' "$results/chr22.recombination.txt" > recombination.tsv
diff -u - recombination.tsv <<'EXPECTED'
NA12879 chr22 20502763 20503731 1 0 1 0
NA12879 chr22 20503755 20510653 0 1 0 0
NA12879 chr22 20959352 20961970 1 0 0 0
EXPECTED

# The dual-parent bracket contributes two candidates with zero nearest gap.
cut -f 2-6,10-11 "$results/chr22.crossovers.tsv" > candidates.tsv
diff -u - candidates.tsv <<'EXPECTED'
#child_id	chromosome	parent	position1	position2	nearest_gap_bp	clustered
NA12879	chr22	father	20502763	20503731	0	1
NA12879	chr22	mother	20502763	20503731	0	1
NA12879	chr22	father	20503755	20510653	24	1
NA12879	chr22	father	20959352	20961970	448699	0
EXPECTED
diff -u - "$results/chr22.crossovers.bed" <<'EXPECTED'
chr22	20502762	20503731	NA12879:chr22:20502763-20503731:father	0	.	20502762	20503731	0	0	1
chr22	20502762	20503731	NA12879:chr22:20502763-20503731:mother	0	.	20502762	20503731	0	0	1
chr22	20503754	20510653	NA12879:chr22:20503755-20510653:father	0	.	20503754	20510653	0	24	1
chr22	20959351	20961970	NA12879:chr22:20959352-20961970:father	0	.	20959351	20961970	0	448699	0
EXPECTED

for cutoff in 24 23; do
    "$nextflow_bin" "${common[@]}" -resume --crossover-distance-bp "$cutoff"
    expected=1
    if [[ "$cutoff" -eq 23 ]]; then
        expected=0
    fi
    [[ "$(cut -f 11 "$results/chr22.crossovers.tsv" | paste -sd ',' -)" == "clustered,1,1,$expected,0" ]]
done

echo "PASS: indexed trio phasing, raw crossover parity, and crossover screening on the chr22 fixture"
