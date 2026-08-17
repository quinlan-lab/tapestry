#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
nextflow_bin="${NEXTFLOW:-nextflow}"
container_image="${TAPESTRY_CONTAINER:-tapestry:dev}"

mkdir -p "${repo_root}/.test-work"
fixture_root="$(mktemp -d -p "${repo_root}/.test-work" e2e.XXXXXX)"
case "${fixture_root}" in
    "${repo_root}/.test-work"/e2e.*) ;;
    *) echo "refusing unsafe test directory: ${fixture_root}" >&2; exit 2 ;;
esac

test_succeeded=false
cleanup() {
    if [[ "${test_succeeded}" == true || -z "${CI:-}" ]]; then
        rm -rf "${fixture_root}"
    else
        echo "CI failure diagnostics retained in ${fixture_root}" >&2
    fi
}
trap cleanup EXIT

mkdir -p "${fixture_root}/fixture"
fixture_name="$(basename "${fixture_root}")"
docker run --rm \
    --user "$(id -u):$(id -g)" \
    --volume "${repo_root}:/work" \
    --workdir /work \
    --env PYTHONDONTWRITEBYTECODE=1 \
    --env PYTHONPATH=/work/src:/work/src/util \
    "${container_image}" python -c \
    'from pathlib import Path; from tests.test_tapestry_validate import make_fixture; import sys; make_fixture(Path(sys.argv[1]), inheritance_informative=True)' \
    "/work/.test-work/${fixture_name}/fixture"

common=(
    run "${repo_root}"
    -profile docker
    --container "${container_image}"
    --outputs-json "${fixture_root}/fixture/miniwdl-outputs.json"
    --ped "${fixture_root}/fixture/family.ped"
    --reference-fasta "${fixture_root}/fixture/data/reference.fa"
    --outdir "${fixture_root}/fixture/results/fixture"
    --project-id fixture
    --samples CHILD
    --regions chr1
    --min-run-markers 1
    -work-dir "${fixture_root}/work"
)

NXF_OFFLINE=true "${nextflow_bin}" -log "${fixture_root}/nextflow.first.log" \
    "${common[@]}" \
    -with-report "${fixture_root}/report.first.html" \
    -with-trace "${fixture_root}/trace.first.tsv" \
    -with-timeline "${fixture_root}/timeline.first.html" \
    -with-dag "${fixture_root}/dag.first.html"

NXF_OFFLINE=true "${nextflow_bin}" -log "${fixture_root}/nextflow.resume.log" \
    "${common[@]}" -resume \
    -with-report "${fixture_root}/report.resume.html" \
    -with-trace "${fixture_root}/trace.resume.tsv" \
    -with-timeline "${fixture_root}/timeline.resume.html" \
    -with-dag "${fixture_root}/dag.resume.html"

results="${fixture_root}/fixture/results/fixture"
test -s "${results}/results-manifest.json"
test -s "${results}/pipeline_info/versions.json"
test -s "${results}/inheritance/fixture.pass.vcf.gz.tbi"
test -s "${results}/samples/CHILD/CHILD.dna-methylation.all-cpgs.bed.gz.tbi"

python3 -c '
import json
from pathlib import Path
import sys
root = Path(sys.argv[1])
document = json.loads((root / "results-manifest.json").read_text())
listed = {"results-manifest.json"}
stack = [document]
while stack:
    value = stack.pop()
    if isinstance(value, dict):
        listed.update(value[key] for key in ("path", "index") if key in value)
        stack.extend(value.values())
    elif isinstance(value, list):
        stack.extend(value)
published = {str(path.relative_to(root)) for path in root.rglob("*") if path.is_file()}
if listed != published:
    raise SystemExit(f"manifest mismatch: missing={sorted(listed-published)} unlisted={sorted(published-listed)}")
' "${results}"

awk -F '\t' '
    NR > 1 {
        process_count++
        if ($5 != "CACHED") exit 1
    }
    END {
        if (process_count != 9) exit 1
    }
' "${fixture_root}/trace.resume.tsv"
test_succeeded=true
echo "PASS: full Docker workflow and -resume cache"
