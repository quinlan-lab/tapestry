# Completed work

This file records implemented Tapestry milestones and settled decisions. It is
historical evidence, not the public input specification. The user-facing
command lives in `README.md`.
Remaining work is tracked in `ROADMAP.md`.

## Generic autosomal pedigree/model-only MVP

Completed 2026-08-15.

Tapestry now consumes completed PacBio HiFi human WGS WDL `family` outputs
through a Tapestry-owned manifest. It does not rerun HiPhase or pb-CpG-tools.
The implemented path supports:

- GRCh38 whole autosomes `chr1` through `chr22`.
- Pedigrees with one or more eligible children.
- Model methylation from the WDL combined/hap1/hap2 CpG BEDs.
- WDL v3.3.0 and v3.3.1 tag/commit pairs.
- Local, Docker, Apptainer, and Slurm Nextflow profiles.

Count methylation, sex chromosomes, other references, arbitrary intervals,
other workflow-engine adapters, and generic trio orchestration were explicitly
left for later phases.

## Settled interface decisions

- The entry point accepts miniwdl `outputs.json`, PED, reference FASTA, and
  output directory directly as Nextflow parameters.
- Validation detects the WDL release and internally generates the canonical
  Tapestry manifest and normalized run. Users do not author either contract.
- Miniwdl adaptation and contract validation are one validator operation; there
  is no preparatory command or intermediate user-facing YAML.
- Scientific overrides are optional Nextflow parameters. Executor, queue,
  container, account, and site resource policy live in Nextflow profiles.
- Omitted `--samples` selects every eligible PED member; duplicate, unknown, or
  ineligible selections are rejected.
- WDL v3.3.0 and v3.3.1 map to their audited commits. Other releases are
  rejected.
- Reference-dependent filenames derive from `reference.name`; they are not
  hard-coded in workflow or manifest-generation code.

## Joint-VCF compatibility decision

The WDL `joint_small_variants_vcf` is accepted provisionally for inheritance
mapping and final CpG annotation. Source inspection established that it contains
the genotype and depth fields needed by pinned `gtg`, is bgzip-compressed and
tabix-indexed, and that input HiPhase orientation does not affect `gtg`'s
inheritance inference.

The WDL artifact is nevertheless transformed rather than a lossless GLnexus
callset. A v3.3.1 split/merge reproduction changed 5,198 records to 5,167:
all-no-call records were removed, some allele keys were normalized or merged,
and depth fields for excluded missing calls were lost. This led to the following
implemented normalization contract:

1. Require exact VCF/PED sample-set equality.
2. Validate indexes, unique samples, `GT`, numeric `QUAL`, compatible contigs,
   and at least one of `DP`, `AD`, or `SD`.
3. Remove HiPhase `PS` and `PF` from the inheritance-oriented copy.
4. Preserve an all-site VCF for concordance and annotation.
5. Derive a map-only VCF containing complete-family genotypes.
6. Apply no implicit PASS-only, biallelic-only, or SNV-only filter.
7. Sort, bgzip, index, and publish both normalized views and their QC.

Inheritance QC includes per-sample called/missing genotype counts and streaming
depth summaries. Called-genotype depth uses `DP`, then summed `AD`, then summed
`SD`; QC records which source supplied each observation.

Real CEPH comparisons for v3.3.0 and v3.3.1 remain release work and therefore
stay in `ROADMAP.md`.

## Implemented workflow

The DSL2 workflow now performs nine stages:

1. Validate and normalize the run, manifest, PED, selected samples, reference,
   VCFs, HiPhase blocks, model BEDs, output location, and release policy.
2. Filter model BEDs while preserving upstream metadata and recording the
   effective coverage threshold.
3. Capture container, executable, Python-package, Nextflow, and upstream WDL
   provenance.
4. Produce cleaned all-site and complete-family map VCFs.
5. Run pinned `gtg-ped-map` and `gtg-concordance`, then sort and index outputs.
6. Generate the configured reference-autosome CpG set from the supplied FASTA.
7. Reconcile HiPhase and inheritance phase and assign model methylation to
   founder haplotypes.
8. Expand each selected sample to reference and sample-created CpGs with variant
   and mismatch annotations.
9. Publish a work-directory-free results manifest containing every artifact,
   index, disabled mode, and per-sample completion status.

Samples without an informative phase intersection produce valid null-founder
outputs and are reported as `no_inheritance_phase`, not falsely as `complete`.

## Runtime and reproducibility

- Python dependencies are pinned in `requirements-pipeline.txt`, including
  `pyBigWig`.
- `gtg` is pinned to commit
  `e12aca6b49ee7208952467db4a2a9e2f79b98efb`; its Cargo lockfile is committed.
- The container includes bcftools/tabix and records tool versions in run
  provenance.
- Production BigWigs use FAI-backed `pyBigWig`, eliminating hard-coded hg38
  chromosome lengths.
- The checksum-pinned UCSC `bedGraphToBigWig` v2.10 binary remains in the image
  only as a parity oracle. The deterministic comparison passes for chromosome
  metadata, intervals, missing values, and values to six decimal places. Its
  SHA-256 is
  `1a1527cf364e1e572a81c7284fc9ccd2b3690b5896baa5b57399864f85ad7771`.

## Completed milestones

- Direct miniwdl JSON parsing, resolved provenance records, and validate-only
  entry.
- Deterministic joint-VCF normalization and inheritance slice.
- Portable model-only founder phasing with empty/no-phase handling.
- Reference CpG generation, all-CpG expansion, BigWigs, QC, and complete results
  manifest.
- Docker runtime and local/Docker/Apptainer/Slurm profiles.
- File-relative input resolution, autosome enforcement, release allowlist, and
  output-collision protection.

## Verification evidence

As of 2026-08-17:

- All 23 containerized unit/regression tests pass with no skips.
- The informative synthetic Docker workflow completes all nine processes.
- A full `-resume` rerun reports all nine processes cached.
- Every published synthetic artifact is represented exactly once by
  `results-manifest.json`.
- The BigWig migration parity test passes against the pinned UCSC binary.
- Targeted Pyright reports zero errors and warnings.
- Nextflow profile configuration, shell syntax, and `git diff --check` pass.

## Continuous integration

Implemented 2026-08-16.

One GitHub Actions workflow now runs on pull requests, pushes to `main`, and
manual dispatch. It deliberately has no runtime matrix and no publication or
write permissions:

- `Static checks` creates a Python 3.11 virtual environment, installs pinned
  runtime and CI dependencies, type-checks the maintained generic modules, and
  checks the E2E shell syntax.
- `Workflow` uses Java 17 and Nextflow 24.04.2, builds `tapestry:ci` with cached
  Docker layers, runs all containerized unit/regression tests, and executes the
  informative workflow twice.
- The resume assertion requires exactly nine trace rows and requires every row
  to be `CACHED`.
- CI fixture generation uses the pipeline image rather than host Python. Failed
  CI runs retain only compact synthetic logs, lifecycle reports, QC, and the
  results manifest for short-lived artifact upload; successful and local runs
  clean up their fixture.
- Third-party actions are pinned to full commit SHAs. The workflow neither
  consumes secrets nor uses `pull_request_target`, and it never pushes an image.

The implementation was reproduced locally with the same Nextflow version: all
23 tests passed, the nine-stage synthetic workflow completed, and all nine
stages were cached on resume. The first hosted green run and repository branch
protection remain release-administration tasks in `ROADMAP.md`.

## Compact user interface

Implemented 2026-08-17.

- The primary command is one resumable Nextflow run; validation remains the
  mandatory first stage, while `-entry validate` is an optional preflight.
- Nextflow accepts miniwdl `outputs.json`, PED, reference, and output directory
  directly. Validation generates the internal canonical manifest and normalized
  run.
- Validation prints the resolved family, samples, WDL release, reference,
  regions, inheritance and methylation thresholds, BigWig choice, and output
  directory.
- Completion prints the results-manifest path and per-sample status counts.
- Startup validation checks contracts, headers, indexes, samples, contigs, and
  metadata without scanning genome-wide VCF/BED records. Record-level checks
  remain with the scientific stages that already read those records.
- Founder methylation phasing reads tabix-indexed hap1/hap2 BEDs one autosome at
  a time, appends a stable-schema BED in canonical order, and writes BigWigs in
  bounded chunks. A whole-genome-versus-chromosome regression covers reversed
  region order, empty chromosomes, and haplotype-specific CpGs.

These tests establish the implemented contract on deterministic fixtures. They
do not replace the real WDL/CEPH, cluster-runtime, hosted-CI, and
release-publication gates in `ROADMAP.md`.
