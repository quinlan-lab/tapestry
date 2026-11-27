#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
nextflow_bin="${NEXTFLOW:-nextflow}"
container_image="${TAPESTRY_TRIO_CONTAINER:-}"
export PATH="$repo_root/.venv/bin:$PATH"
export NXF_OFFLINE=true

usage() {
    echo "Usage: $0"
    echo "Check local chromosome-parallel trio phasing and model methylation on the chr22 fixture."
    echo "Set TAPESTRY_TRIO_CONTAINER=ghcr.io/quinlan-lab/tapestry:dev to test the default Docker image."
}
if [[ ${1:-} == -h || ${1:-} == --help ]]; then
    usage
    exit 0
fi
if [[ $# -ne 0 ]]; then
    usage >&2
    exit 2
fi
if [[ -n "$container_image" ]]; then
    required_tools=("$nextflow_bin" docker)
else
    required_tools=("$nextflow_bin" bcftools bgzip tabix samtools aligned_bam_to_cpg_scores python3)
fi
for tool in "${required_tools[@]}"; do
    command -v "$tool" >/dev/null || { echo "Missing prerequisite: $tool" >&2; exit 1; }
done
# Use the same runtime for workflow tasks and scientific assertions.
runtime() {
    if [[ -n "$container_image" ]]; then
        docker run --rm --user "$(id -u):$(id -g)" \
            --volume "$repo_root:$repo_root" --workdir "$PWD" \
            --env PYTHONDONTWRITEBYTECODE=1 \
            --env PYTHONPATH="$repo_root/src:$repo_root/src/util" \
            "$container_image" "$@"
    else
        "$@"
    fi
}
if [[ "$(runtime python3 -m whatshap --version)" != 2.8 ]]; then
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

# The 1 Mb fixture fits small local Docker VMs and CI runners.
cat > fixture.config <<'CONFIG'
process {
    withLabel: /trio_.*/ {
        cpus = 2
        memory = '2 GB'
    }
}
CONFIG
common=(
    run "$repo_root" -ansi-log false -c "$root/fixture.config"
    --ped "$fixture/trio.ped"
    --joint-vcf "$fixture/CEPH-1463.joint.GRCh38.deepvariant.glnexus.vcf.gz"
    --reference-fasta "$fixture/chr22.fa"
    --child-bam "$fixture/NA12879.GRCh38.haplotagged.bam"
    --father-bam "$fixture/NA12877.GRCh38.haplotagged.bam"
    --mother-bam "$fixture/NA12878.GRCh38.haplotagged.bam"
    --chromosomes chr22 --region chr22:20500000-21500000
    --outdir results
)
if [[ -n "$container_image" ]]; then
    common+=(-profile docker)
    if [[ "$container_image" != ghcr.io/quinlan-lab/tapestry:dev ]]; then
        common+=(--container "$container_image")
    fi
fi
"$nextflow_bin" "${common[@]}"

results=results/per_chrom_phased
vcf="$results/chr22.phased.vcf.gz"
test -s "$vcf.tbi"
diff -u <(printf '%s\n' NA12879 NA12877 NA12878) <(runtime bcftools query -l "$vcf")
runtime bcftools query -r chr22:20500000-21500000 -f '%CHROM\t%POS[\t%GT]\n' "$vcf" > records.tsv
awk '
    $1 != "chr22" || $2 < 20500000 || $2 > 21500000 { exit 1 }
    $3 ~ /[|]/ { phased = 1 }
    END { if (NR != 2303 || !phased) exit 1 }
' records.tsv
[[ "$(runtime bcftools index -n "$vcf")" == 2303 ]]

# Read selection differs from the macOS/Python 3.13 host fixture. Keep strict
# expectations for each tested runtime instead of accepting either result.
awk '!/^#/ { print $1, $2, $3, $4, $5, $6, $7, $8 }' "$results/chr22.recombination.txt" > recombination.tsv
cut -f 2-6,10-11 "$results/chr22.crossovers.tsv" > candidates.tsv
if [[ -n "$container_image" ]]; then
    diff -u - recombination.tsv <<'EXPECTED'
NA12879 chr22 20502763 20503731 1 1 1 0
NA12879 chr22 20959352 20961970 1 0 0 0
EXPECTED
    diff -u - candidates.tsv <<'EXPECTED'
#child_id	chromosome	parent	position1	position2	nearest_gap_bp	clustered
NA12879	chr22	mother	20502763	20503731	455621	0
NA12879	chr22	father	20959352	20961970	455621	0
EXPECTED
    diff -u - "$results/chr22.crossovers.bed" <<'EXPECTED'
chr22	20502762	20503731	NA12879:chr22:20502763-20503731:mother	0	.	20502762	20503731	0	455621	0
chr22	20959351	20961970	NA12879:chr22:20959352-20961970:father	0	.	20959351	20961970	0	455621	0
EXPECTED
else
    # Compare the raw candidate brackets and transmitted states from the quality-filtered experiment.
    diff -u - recombination.tsv <<'EXPECTED'
NA12879 chr22 20502763 20503731 1 0 1 0
NA12879 chr22 20503755 20510653 0 1 0 0
NA12879 chr22 20959352 20961970 1 0 0 0
EXPECTED

    # The dual-parent bracket contributes two candidates with zero nearest gap.
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
fi

hap_maps=results/per_chrom_hap_maps/chr22
paternal="$hap_maps/NA12879.hap-map-blocks.paternal.bed.gz"
# Query the bases immediately outside the one-based inclusive isolated bracket.
left="$(runtime tabix "$paternal" chr22:20959351-20959351 | cut -f 4 | cut -d, -f 1-2)"
right="$(runtime tabix "$paternal" chr22:20961971-20961971 | cut -f 4 | cut -d, -f 1-2)"
[[ "$left" == A,1.0 && "$right" == B,1.0 ]]
assignments="$(runtime tabix "$paternal" chr22:20959352-20961970)"
[[ -z "$assignments" ]]

for parent in paternal maternal; do
    map="$hap_maps/NA12879.hap-map-blocks.$parent.bed"
    mismatch="$hap_maps/NA12879.bit-vector-sites-mismatches.$parent"
    test -s "$map.header"
    test -s "$map.gz.tbi"
    if [[ -z "$container_image" ]]; then
        # Include the gap between the two clustered brackets and both outside markers.
        assignments="$(runtime tabix "$map.gz" chr22:20502763-20510653)"
        [[ -z "$assignments" ]]
    elif [[ "$parent" == maternal ]]; then
        assignments="$(runtime tabix "$map.gz" chr22:20502763-20503731)"
        [[ -z "$assignments" ]]
    fi
    # This fixture has no mismatches in the resolved blocks.
    test -f "$mismatch.bed"
    test ! -s "$mismatch.bed"
    test -s "$mismatch.bed.header"
    test -s "$mismatch.vcf.gz.tbi"
    mismatches="$(runtime bcftools view -H "$mismatch.vcf.gz")"
    [[ -z "$mismatches" ]]
done

checker_args=()
cutoffs=(24 23)
if [[ -n "$container_image" ]]; then
    checker_args+=(--docker-runtime)
    cutoffs=(455621 455620)
fi
PYTHONPATH="$repo_root/src:$repo_root/src/util" runtime python3 "$repo_root/tests/check_trio_model_fixture.py" results "${checker_args[@]}"

for cutoff in "${cutoffs[@]}"; do
    "$nextflow_bin" "${common[@]}" -resume --crossover-distance-bp "$cutoff"
    expected=1
    if [[ "$cutoff" -eq "${cutoffs[1]}" ]]; then
        expected=0
    fi
    if [[ -n "$container_image" ]]; then
        [[ "$(cut -f 11 "$results/chr22.crossovers.tsv" | paste -sd ',' -)" == "clustered,$expected,$expected" ]]
        for parent in paternal maternal; do
            assignments="$(runtime tabix "$hap_maps/NA12879.hap-map-blocks.$parent.bed.gz" chr22:20510653-20510653)"
            if [[ "$expected" == 1 ]]; then
                [[ -z "$assignments" ]]
            else
                [[ -n "$assignments" ]]
            fi
        done
    else
        [[ "$(cut -f 11 "$results/chr22.crossovers.tsv" | paste -sd ',' -)" == "clustered,1,1,$expected,0" ]]
    fi
done

echo "PASS: indexed trio phasing, screened parental maps, and model methylation on the chr22 fixture"
