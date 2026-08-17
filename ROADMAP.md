# Tapestry roadmap

The generic autosomal pedigree/model-only workflow is implemented. See
[`completed.md`](completed.md) for settled decisions, implementation history,
and verification evidence. This roadmap contains only remaining work and future
scope.

The public Nextflow parameters are documented in `README.md`. Resolved JSON
records are pipeline outputs for provenance and downstream stages, not a second
user workflow.

## Current product boundary

Tapestry starts after a completed miniwdl run of the PacBio HiFi human WGS WDL
`family` workflow:

```bash
# Validate, then run inheritance, founder phasing, and all-CpG expansion.
nextflow run . -profile docker \
  --outputs-json /path/to/outputs.json \
  --ped /path/to/family.ped \
  --reference-fasta /path/to/GRCh38.fa \
  --outdir /path/to/results \
  -resume

# Optional validation-only preflight.
nextflow run . -entry validate -profile docker \
  --outputs-json /path/to/outputs.json \
  --ped /path/to/family.ped \
  --reference-fasta /path/to/GRCh38.fa \
  --outdir /path/to/results
```

The workflow supports GRCh38 whole autosomes, pedigree mode, model methylation,
and WDL v3.3.0/v3.3.1. It consumes WDL-produced HiPhase and
model-CpG artifacts;
it does not rerun those upstream tools. Count methylation, sex chromosomes,
other references, arbitrary intervals, other engine adapters, and generic trio
orchestration are outside the current boundary.

## Immediate priority: release readiness

Keep this phase narrow. Do not add new scientific modes until these gates are
complete.

### 1. Real WDL/CEPH parity

Run the same CEPH family through PacBio WDL v3.3.0 and v3.3.1 and then through
Tapestry. For each release, retain a small de-identified contract fixture plus a
reproducible comparison report.

Compare the new path with the legacy path using identical `gtg` versions and
thresholds:

- normalized and complete-family marker positions;
- semantic inheritance assignments and block boundaries;
- concordance pass, fail, and no-call counts;
- phased genotype alleles and mismatch sites;
- founder model-methylation values;
- final all-CpG coordinates, null patterns, and allele flags;
- record counts, per-sample genotype/depth QC, and effective thresholds.

Explain differences caused by WDL split/merge behavior, complete-family map
normalization, or model-coverage filtering. Do not require byte-identical VCFs
when the semantic result agrees.

Exit criteria:

- Both v3.3.0 and v3.3.1 have recorded fixtures, commands, versions, deltas, and
  dispositions.
- No unexplained scientific difference remains.
- The support statement in `README.md` names the releases with completed parity
  evidence.

### 2. Runtime smoke tests

Run the same informative synthetic fixture with:

- Docker/local execution;
- Apptainer on the target cluster;
- Slurm with the intended queue/account profile overlay.

Record Nextflow and runtime versions, commands, resource traces, and any
site-local configuration needed outside the repository. Revise the current
all-CpG memory defaults only from observed traces.

### 3. Release and distribution

- Confirm the first hosted GitHub Actions run is green, then require the
  `Static checks` and `Workflow` jobs in branch protection and enable weekly
  Dependabot updates for pinned CI dependencies.
- Publish the OCI image and record its immutable digest.
- Decide whether supported runs start from a tagged GitHub revision, a local
  bundle, or both.
- Tag the repository only after parity, CI, and runtime gates pass.
- Add the released command and digest to `README.md`.

## Next functional phase

After release readiness, the recommended next addition is optional count-based
methylation. It uses WDL haplotagged BAMs without changing the model inputs.

Before implementation, settle the smallest interface extension that answers:

- whether count mode is explicitly enabled or inferred from named inputs;
- where minimum coverage and MAPQ live;
- how count-disabled and count-enabled outputs share one stable table schema;
- whether count generation is a pipeline stage or an externally supplied
  artifact contract.

Model-only runs must remain valid, keep typed-null count columns, and avoid
paying count-mode compute or I/O costs.

## Later decisions

These are intentionally deferred. Resolve them from concrete fixtures rather
than expanding the interface speculatively.

1. **Additional WDL releases.** Add another stable v3.x release only after
   source audit and parity evidence.
2. **Workflow-engine adapters.** Add Cromwell or Terra adapters only when
   needed. Every adapter must emit the existing internal canonical Tapestry
   manifest; engine dialects must not leak into scientific stages.
3. **Trio support.** Decide whether the generic product retains WhatsHap/pedMEC,
   converges trio and pedigree outputs, or versions trio separately.
4. **Reference and region expansion.** Add sex chromosomes, non-GRCh38
   references, arbitrary intervals, or alternate contigs one policy at a time,
   with reference-aware fixtures.
5. **Threshold harmonization.** Decide whether future products expose WDL model
   data at coverage 4, the Tapestry-filtered coverage-10 view, or both. Preserve
   the upstream threshold in metadata.

## Guardrails

- Keep the primary interface to one Nextflow invocation using miniwdl
  `outputs.json`, PED, reference, output directory, and optional scientific
  parameter overrides.
- Record pedigree mode, GRCh38, `gtg`, model mode, and all effective scientific
  parameters in the normalized run.
- Keep execution policy in Nextflow profiles.
- Treat the canonical manifest as the only internal WDL-output contract.
- Never derive sample IDs, roles, reference labels, or output names from
  hard-coded CEPH/GRCh38 constants in generic workflow code.
- Preserve all-site and complete-family VCF branches; do not add undocumented
  PASS, biallelic, or SNV filters.
- Preserve zero-based half-open BED and one-based VCF coordinate conventions.
- Record configuration, versions, reference, QC, and per-sample completion
  status with every run.
- Prefer extending existing validation and process boundaries over adding
  wrappers or alternate orchestration paths.
