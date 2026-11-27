# Developer tests

These checks verify changes to pipeline code and container dependencies. Run
them from the repository root with Docker and Nextflow installed. The
`TAPESTRY_*_CONTAINER` variables select the image to test; substitute a local
image such as `tapestry:dev` when testing a local build.

## Trio smoke test

Prepare the checked-in chr22 fixture's reference once:

```bash
docker run --rm --user "$(id -u):$(id -g)" \
  -v "$PWD:/work" -w /work ghcr.io/quinlan-lab/tapestry:dev \
  bash ./trio_dev_data_get_ref.sh trio_dev_data
```

Run the trio workflow and check crossover-distance boundaries, haplotype
labels, unresolved nulls, BED/BAM indexes, and nine BigWigs:

```bash
TAPESTRY_TRIO_CONTAINER=ghcr.io/quinlan-lab/tapestry:dev tests/run_nextflow_trio_smoke.sh
```

Workflow tasks and scientific assertions use the selected container. A failed
run retains diagnostics under `.test-work/trio-smoke.*`.

## Pedigree end-to-end test

Generate a synthetic informative family and run the pedigree workflow twice.
The second run must reuse every process through Nextflow's `-resume` cache:

```bash
TAPESTRY_CONTAINER=ghcr.io/quinlan-lab/tapestry:dev tests/run_nextflow_e2e.sh
```

The test generates its own inputs inside the container. Failed CI runs retain
diagnostics under `.test-work/e2e.*`; ordinary local runs clean up automatically.
