#!/usr/bin/env bash
# Reproduce the second deck's 1 Mb trio workflow and its four data figures.
set -euo pipefail
repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$repo_root"

usage() {
    echo "Usage: $0 [OUTPUT_DIR]"
    echo "Run the chr22 fixture into a fresh directory, then generate slide figures."
    echo "Default: trio_dev_data/output/nextflow-trio-model"
}
if [[ ${1:-} == -h || ${1:-} == --help ]]; then usage; exit 0; fi
if [[ $# -gt 1 || ${1:-} == -* ]]; then usage >&2; exit 2; fi
export PATH="$repo_root/.venv/bin:$PATH"
export PYTHONPATH="$repo_root/src:$repo_root/src/util"
export PYTHONDONTWRITEBYTECODE=1
python="$repo_root/.venv/bin/python"
nextflow_bin="${NEXTFLOW:-nextflow}"
for tool in "$python" "$nextflow_bin" bcftools bgzip tabix samtools aligned_bam_to_cpg_scores; do
    command -v "$tool" >/dev/null || { echo "Missing prerequisite: $tool" >&2; exit 1; }
done
[[ "$("$python" -m whatshap --version)" == 2.8 ]] || { echo "Install WhatsHap 2.8 in .venv" >&2; exit 1; }
"$python" -c 'import matplotlib, polars, pysam, pyBigWig'
fixture="$repo_root/trio_dev_data/input"
source src/util/logging.sh
source src/util/trio_ped.sh
read_trio_ids "$fixture/trio.ped"
for file in "$fixture/chr22.fa" "$fixture/chr22.fa.fai" \
    "$fixture/CEPH-1463.joint.GRCh38.deepvariant.glnexus.vcf.gz" \
    "$fixture/$kid_id.GRCh38.haplotagged.bam" \
    "$fixture/$dad_id.GRCh38.haplotagged.bam" "$fixture/$mom_id.GRCh38.haplotagged.bam"; do
    [[ -s "$file" ]] || { echo "Missing input: $file (see analysis/nextflow-trio-model/README.md)" >&2; exit 1; }
done
out=${1:-trio_dev_data/output/nextflow-trio-model}
if [[ -e "$out" ]]; then echo "Output already exists; choose a fresh directory: $out" >&2; exit 1; fi
mkdir -p "$out"
out="$(cd "$out" && pwd)"
mkdir "$out/pipeline"
cd "$out/pipeline"
args=(run "$repo_root" -entry trio -ansi-log false
    --ped "$fixture/trio.ped"
    --joint-vcf "$fixture/CEPH-1463.joint.GRCh38.deepvariant.glnexus.vcf.gz"
    --reference-fasta "$fixture/chr22.fa"
    --child-bam "$fixture/$kid_id.GRCh38.haplotagged.bam"
    --father-bam "$fixture/$dad_id.GRCh38.haplotagged.bam"
    --mother-bam "$fixture/$mom_id.GRCh38.haplotagged.bam"
    --chromosomes chr22 --region chr22:20500000-21500000
    --crossover-distance-bp 10000 --min-coverage 10 --trio-min-mapq 1
    --outdir "$out" -with-trace "$out/trace.tsv")
printf '%q ' "$nextflow_bin" "${args[@]}" > "$out/run.command.txt"
printf '\n' >> "$out/run.command.txt"
{ "$python" -m whatshap --version; bcftools --version; samtools --version;
  aligned_bam_to_cpg_scores --version; "$nextflow_bin" -version; } > "$out/versions.txt"
log_info "Running chr22 trio model experiment: $out"
"$nextflow_bin" "${args[@]}" 2>&1 | tee "$out/run.log"
mkdir "$out/phase-blocks"
for sample in "$kid_id" "$dad_id" "$mom_id"; do
    "$python" -m whatshap stats --sample "$sample" \
        --block-list "$out/phase-blocks/$sample.blocks.tsv" \
        "$out/per_chrom_phased/chr22.phased.vcf.gz" \
        > "$out/phase-blocks/$sample.stats.txt" 2> "$out/phase-blocks/$sample.stats.log"
done
"$python" "$repo_root/tests/check_trio_model_fixture.py" "$out"
"$python" "$repo_root/analysis/nextflow-trio-model/summarize.py" --run-dir "$out"
log_info "Completed chr22 trio model experiment: $out"
