# Generic autosomal pedigree/model-only MVP

**Archived implementation plan.** The planned MVP is implemented. See
`completed.md` for the completion record and `ROADMAP.md` for remaining release
work. This file retains detailed design rationale and acceptance criteria.
Its interface examples are historical; `README.md` is the current command
contract.

## Goal

Build a Nextflow DSL2 workflow that consumes a Tapestry-defined manifest for
outputs from the PacBio HiFi human WGS WDL, derives pedigree inheritance, phases
WDL-produced model methylation by founder, expands the result to all CpGs, and
publishes reproducible QC and provenance.

This phase deliberately supports one narrow, complete path:

- GRCh38.
- Pedigree mode.
- Model methylation only.
- PacBio HiFi human WGS WDL v3.3.0 or v3.3.1 output represented by Tapestry's
  canonical manifest.
- Autosomes only (chr1-chr22).

Count methylation, trio-only operation, other references, arbitrary intervals,
and direct Cromwell/miniwdl metadata adapters are deferred. Existing command-line
entry points remain compatible while their implementation is refactored for
model-only operation.

## User interface

The normative schema-v1 contracts live in `schemas/run.schema.json` and
`schemas/upstream-manifest.schema.json`; `examples/generic/` contains portable
examples. This historical plan does not define a second manifest dialect.

The public entry point is `--run-config`, rather than Nextflow's
`-params-file`. This lets Tapestry resolve paths relative to the YAML or JSON
file that contains them.

```bash
# Validate inputs without running scientific processes.
nextflow run . -entry validate -profile docker \
  --run-config family.yaml

# Run locally with Docker.
nextflow run . -profile docker \
  --run-config family.yaml \
  -resume

# Run a released pipeline on a cluster.
nextflow run quinlan-lab/tapestry -r <tag> \
  -profile slurm,apptainer \
  --run-config family.yaml \
  -work-dir /path/to/work \
  -resume
```

Scientific settings belong in the run configuration. Executor, queue, resource,
and container-runtime settings belong in Nextflow profiles.

### Run configuration

YAML is the documented format because it is easiest for people to maintain.
JSON is accepted as an equivalent machine-generated representation. Both are
validated against the same JSON Schema Draft 2020-12 schema.

```yaml
schema_version: 1
mode: pedigree

project:
  id: family_001
  outdir: results/family_001

pedigree:
  ped: inputs/family.ped

reference:
  name: GRCh38
  fasta: references/hg38.analysisSet.fa
  fasta_index: references/hg38.analysisSet.fa.fai
  # fasta_gzi: references/hg38.analysisSet.fa.gz.gzi

upstream:
  manifest: inputs/pacbio-wdl.tapestry.json
  # Required only to evaluate another stable v3.x release.
  allow_unaudited_release: false

samples:
  # Omit to process every eligible PED member.
  include: [child_01, child_02]

inheritance:
  method: gtg
  map:
    min_qual: 20
    min_depth: 10
    min_run_markers: 10
  concordance:
    min_qual: 20
    min_depth: 5

methylation:
  modes: [model]
  min_coverage: 10
  mismatch_window_bp: 50

regions:
  # Omit for chr1-chr22. Schema v1 accepts whole autosomes only.
  include: [chr1, chr2]

outputs:
  bigwig: true
```

Rules for schema version 1:

- Reject unknown keys, custom YAML tags, and unsupported schema versions.
- Accept `mode: pedigree` and `methylation.modes: [model]` only.
- Accept `reference.name: GRCh38` only.
- Accept whole autosomes `chr1` through `chr22` only; reject chrX, chrY,
  mitochondrial contigs, alternate contigs, and interval strings.
- Resolve paths in the run configuration relative to that file.
- Resolve paths in the upstream manifest relative to that manifest.
- Treat an omitted `samples.include` as all eligible PED members. An explicitly
  empty list is invalid.

### Canonical upstream manifest

Tapestry owns a small, versioned manifest contract rather than binding the
pipeline to one workflow engine's metadata format.

```json
{
  "schema_version": 1,
  "provider": "pacbio_hifi_human_wgs_wdl",
  "workflow": {
    "name": "humanwgs_family",
    "release": "v3.3.1",
    "commit": "477ef39ad69e86e90897ea7e313b86bfc12a2a96"
  },
  "family_id": "family_001",
  "joint_small_variants": {
    "vcf": "outputs/family.small_variants.phased.vcf.gz",
    "index": "outputs/family.small_variants.phased.vcf.gz.tbi"
  },
  "samples": [
    {
      "id": "child_01",
      "phased_small_variants": {
        "vcf": "outputs/child_01.small_variants.phased.vcf.gz",
        "index": "outputs/child_01.small_variants.phased.vcf.gz.tbi"
      },
      "phase_blocks": "outputs/child_01.hiphase.blocks.tsv",
      "cpg_model": {
        "combined": {
          "bed": "outputs/child_01.cpg.combined.bed.gz",
          "index": "outputs/child_01.cpg.combined.bed.gz.tbi"
        },
        "hap1": {
          "bed": "outputs/child_01.cpg.hap1.bed.gz",
          "index": "outputs/child_01.cpg.hap1.bed.gz.tbi"
        },
        "hap2": {
          "bed": "outputs/child_01.cpg.hap2.bed.gz",
          "index": "outputs/child_01.cpg.hap2.bed.gz.tbi"
        }
      }
    }
  ]
}
```

The manifest may also carry provenance-only paths such as phase statistics,
haplotags, haplotagged BAMs, and the original workflow-engine metadata. Those
fields are recorded but are not required scientific inputs for this phase.

Support WDL v3.3.0 and v3.3.1 by default and run both through their own parity
fixtures before release. Reject any other release by default. A user may set
`upstream.allow_unaudited_release: true` to deliberately evaluate another stable
v3.x release; validation must emit a prominent warning and record the opt-in in
provenance. The opt-in never permits a prerelease or a different major version.

## PED and sample contract

Normalize PED missing-parent tokens `0`, `.`, and `NA` to one internal missing
value. Validation must enforce:

- Exactly six PED fields per record.
- One family ID per run.
- Unique individual IDs.
- Parent IDs that either exist in the PED or are missing.
- An acyclic pedigree graph.
- Parent/child sex assignments compatible with PED roles where sex is known.
- Exact equality among the PED, joint VCF, and manifest sample sets. Ordering may
  differ.

An eligible output sample is a PED member whose two parents are both present in
the run. With no include-list, process every eligible member. A supplied
include-list must be nonempty and contain eligible members only.

Artifacts needed only by selected output samples may be null for unselected
members. Required joint and selected-sample artifacts may not be null.

## Implementation

### 1. Validation and workflow foundation

Use plugin-free Nextflow DSL2 compatible with Nextflow 24.04.2 or newer. Provide
both a `validate` workflow and the default full workflow.

A small Python validator should:

1. Parse enough of the YAML or JSON to locate the upstream manifest.
2. Validate both documents against committed schemas.
3. Resolve and normalize paths.
4. Validate the PED, sample relationships, and selected samples.
5. Inspect required input files and their indexes.
6. Inspect VCF headers and sample names.
7. Validate expected HiPhase phase-block columns and model BED headers.
8. Verify compatible contig names and lengths, sorted coordinates, and a safe
   output location.
9. Emit normalized artifacts consumed by every later process:
   - resolved run configuration JSON;
   - resolved upstream manifest JSON;
   - normalized `gtg` PED;
   - selected-sample table;
   - validation report;
   - validation-success sentinel.

All scientific processes must depend on the validation sentinel. This prevents
partial scientific work from starting while validation is still underway.

Compute a configuration fingerprint from the normalized run configuration,
manifest, reference identity, and relevant tool versions. Refuse to publish into
an existing run directory with a different fingerprint. Permit `-resume` when
the fingerprint is unchanged.

Do not depend on `nf-schema` in this phase: current releases require a newer
Nextflow baseline than the repository's 24.04.2 target.

### 2. Reproducible runtime

Create one pinned, multi-stage OCI image based on Python 3.11. It contains:

- Pinned Python dependencies.
- `bcftools`, HTSlib, `bgzip`, and `tabix`.
- `gtg` built from exact commit
  `e12aca6b49ee7208952467db4a2a9e2f79b98efb`.
- A committed Cargo lockfile with `hts-sys` pinned to 2.2.0; 2.2.1 is
  incompatible with the selected `rust-htslib` release.

Publish the image as `ghcr.io/quinlan-lab/tapestry` and pin releases by digest in
the production Nextflow configuration. Record image digest and tool versions in
provenance.

Replace hard-coded hg38 BigWig construction with the already pinned `pyBigWig`
dependency using contig lengths from the supplied FASTA index. This changes the
generator from legacy `bedGraphToBigWig`, so the migration requires chromosome,
coordinate, missing-value, and numeric-value equivalence tests against legacy
BigWigs, not merely valid output headers.

The deterministic migration fixture passes against UCSC `bedGraphToBigWig`
v2.10 (SHA-256
`1a1527cf364e1e572a81c7284fc9ccd2b3690b5896baa5b57399864f85ad7771`).
It checks chromosome names and lengths, interval placement, missing positions,
and values to six decimal places. The checksum-pinned binary remains in the
container as a test oracle; production tracks continue to use `pyBigWig`.

Provide Docker, Apptainer, and Slurm profiles. Keep resources configurable in
profiles. A reasonable initial allocation for the CpG expansion stage is 4 CPUs
and 48 GB RAM, to be revised using traces from the CEPH acceptance run.

### 3. Joint VCF normalization and inheritance

The WDL family workflow produces its final joint phased VCF by running GLnexus,
splitting to per-sample VCFs with `--exclude-uncalled`, phasing each sample with
HiPhase, and merging those VCFs. That representation is suitable for
concordance and final annotation, but it changes missing calls and can change
the sites and depth statistics used by `gtg-ped-map`.

Use two explicit branches:

```text
WDL joint VCF
  -> cleaned all-site VCF
       -> gtg-concordance
       -> final variant annotation
  -> complete-family map VCF
       -> gtg-ped-map
```

The cleaned all-site VCF must:

- Keep selected GRCh38 contigs and the exact family samples.
- Preserve records, QUAL, FILTER, INFO, alleles, GT, and depth fields.
- Remove `PS` and `PF` FORMAT fields conditionally because they may describe
  phasing before inheritance concordance.
- Apply no implicit PASS-only, biallelic-only, or SNV-only filters.
- Be deterministically sorted, bgzip-compressed, and tabix-indexed.

The map VCF is the cleaned VCF further restricted to sites with a called GT in
every PED sample. This is a map-only representation. The all-site VCF remains
the source for concordance and final variant annotation.

This rule is based on a regression against the upstream split/merge behavior:
the split/merge path dropped all-no-call records, rewrote some variant keys, and
erased `DP`/`AD` for missing calls. `gtg-ped-map` then produced different depth
statistics and inheritance maps. Restricting map construction to complete-family
sites made marker positions and the semantic inheritance map agree. The final
marker-description ordering was nondeterministic, so regression comparisons must
compare semantic fields rather than raw lines.

Run `gtg-ped-map` on the complete-family VCF. Deterministically sort its IHT,
marker, and recombination outputs. Run `gtg-concordance` against the cleaned
all-site VCF, then sort, compress, and index the pass and fail VCFs.

Family QC must include original, cleaned, complete-map, marker, inheritance,
concordance-pass, and concordance-fail record counts; per-sample called and
missing genotype counts; depth summaries; and map/concordance parameters.
Depth summaries cover called genotypes and use the first available
`DP`, `AD`, or `SD` field in that order, summing vector-valued allele depths.

### 4. Model-only founder phasing

Refactor pedigree phasing functions around enabled methylation modes. Existing
legacy CLI usage with both count and model inputs remains valid, but count inputs
become optional when only model mode is enabled.

For each selected sample:

1. Filter WDL combined, haplotype 1, and haplotype 2 model BEDs to coverage 10 or
   greater.
2. Preserve their metadata header and add the effective Tapestry threshold.
3. Record pre-filter, post-filter, and upstream-declared threshold information in
   QC.
4. Combine the sample HiPhase VCF and phase-block table with the family
   inheritance pass VCF and IHT.
5. Build phase intersections and mismatch outputs.
6. Assign haplotype model methylation values to founders.

Retain full and left joins used by the current implementation so absence of a
value is not confused with a zero value.

If a selected sample has no informative inheritance/phase intersection, publish
valid outputs with founder fields null and set QC status to
`no_inheritance_phase`. An empty family IHT is a fatal family-level error.

Produce model BigWigs only. Keep the established final table schema stable:
legacy count-mode columns remain present as typed null columns, and output
metadata records `enabled_modes: [model]`.

### 5. All-CpG expansion and publication

Generate the reference CpG set once from the configured FASTA and FAI for
chr1-chr22, or the configured whole-autosome subset.

Full-join reference CpGs with the filtered combined model BED so sample-created
CpGs are retained. Then join founder-phased model values, mismatch proximity, and
SNV annotations from the cleaned joint VCF.

Replace shell-command strings in Python with argument arrays or direct library
APIs. Write process outputs atomically before Nextflow publishes them.

Use this output layout:

```text
<outdir>/
  pipeline_info/
  reference/
  inheritance/
  samples/
    <sample-id>/
  results-manifest.json
```

`pipeline_info/` contains resolved inputs, validation results, versions, and the
configuration fingerprint. Nextflow trace/report/timeline/DAG files live in
`.nextflow-reports/` under the launch directory by default because they are
finalized after the scientific DAG; users may redirect them with Nextflow's
standard `-with-*` options.

`inheritance/` contains normalized VCFs, IHT, markers, recombinations,
concordance pass/fail VCFs, indexes, and family QC.

Each sample directory contains haplotype/founder maps, mismatch BED/VCF outputs,
founder model BED and BigWigs, the all-CpG table, and sample QC.

The all-CpG table header records the pipeline version, configuration
fingerprint, reference, coverage threshold, and enabled modes. Upstream model
BED metadata is preserved and augmented with the Tapestry threshold. Native
`gtg`/legacy-compatible tables retain their required headers; their versions and
parameters are linked through `pipeline_info/versions.json` and QC sidecars. No
published metadata exposes ephemeral Nextflow work paths.

`results-manifest.json` lists every published result, its index where relevant,
its completion status, selected samples, and explicitly disabled output modes.

## Delivery milestones

### Milestone 1: Contracts and validate-only workflow

- Add run and manifest schemas.
- Add valid and invalid example configurations.
- Implement path resolution, PED/sample validation, file/header inspection, and
  normalized validation artifacts.
- Add `validate` workflow, Docker development profile, and validation tests.

This milestone establishes the public contract before scientific refactoring.

### Milestone 2: Inheritance slice

- Build cleaned all-site and complete-family map VCFs.
- Package pinned `gtg` and supporting tools.
- Run map and concordance stages.
- Publish deterministic inheritance outputs and family QC.
- Add split/merge compatibility regressions.

### Milestone 3: Per-sample model phasing

- Refactor current pedigree code for model-only operation.
- Filter WDL model BEDs to the effective coverage threshold.
- Phase each eligible selected sample by founder.
- Publish stable-schema model outputs, model BigWigs, mismatches, and sample QC.

### Milestone 4: All-CpG output

- Generate reusable reference CpGs.
- Join model, founder, mismatch, and variant annotations.
- Preserve sample-created CpGs and typed-null count columns.
- Publish indexed/compressed outputs and the results manifest.

### Milestone 5: Operational release

- Add Slurm and Apptainer profiles and tune resources.
- Verify restartability and publication collision behavior.
- Complete README, ROADMAP, AGENTS, schema, and migration documentation.
- Add CI and run the CEPH parity gate.
- Pin the released image digest and tag the first generic MVP release.

## Test and acceptance plan

### Schema and validation tests

Cover YAML/JSON normalization, relative paths, unknown keys, unsupported schema
versions, default acceptance of v3.3.0 and v3.3.1, rejection of other releases,
explicit stable-v3.x opt-in warnings, PED graph errors, mismatched sample sets,
ineligible selections, required null/missing artifacts, missing indexes,
non-autosomal regions, contig length mismatches, and malformed HiPhase TSV/model
BED headers.

### VCF regression tests

Assert that normalization:

- Removes `PS`/`PF` when present and tolerates their absence.
- Preserves all other required fields and records.
- Uses complete-family sites for map construction only.
- Produces sorted, valid, indexed VCFs.
- Applies no undocumented record filters.

Use a synthetic split/merge fixture to prove that normalized raw and WDL-style
remerged inputs produce the same marker positions and semantic IHT result.

### Methylation unit tests

Cover coverage filtering, preservation of headers, typed-null count columns,
full/left join behavior, no-inheritance-phase samples, empty mismatch sets,
multiallelic and missing variant annotations, BigWig chromosome headers, and the
fixed final column contract. On a deterministic fixture, compare `pyBigWig`
tracks with legacy `bedGraphToBigWig` tracks for chromosome order and lengths,
interval coordinates, missing intervals, and numeric values (with a documented
floating-point tolerance).

### End-to-end tests

Create a tiny synthetic GRCh38-like reference, pedigree, joint and per-sample
VCFs, HiPhase block files, and model BEDs. Run both validation and the full Docker
workflow. Assert output contents, indexes, provenance, results manifest, and
successful `-resume` behavior.

Retain the existing recombination-deduplication regression and run Python type
checks and shell syntax checks in CI.

### CEPH parity gate

Run PacBio WDL v3.3.0 and v3.3.1 CEPH family outputs through the new workflow,
using a version-specific fixture and the same scientific thresholds as the
legacy path. Both releases must pass before they are advertised as supported.
Compare:

- Marker positions and semantic IHT assignments.
- Concordance pass/fail classifications.
- Mismatch positions.
- Founder model values.
- Final all-CpG coordinates, null patterns, and allele flags.
- `pyBigWig` output against legacy `bedGraphToBigWig` for chromosome metadata,
  interval coordinates, missing intervals, and methylation values within the
  documented floating-point tolerance.

Block release on unexplained differences. Document expected differences caused
by the explicit complete-family map normalization, model coverage filter, or
stable typed-null output columns.

## Completion criteria

The MVP is complete when:

- One YAML file plus one canonical WDL-output manifest can validate and run a
  multi-child pedigree without source edits.
- Every selected eligible PED member receives model-only founder-phased and
  all-CpG outputs.
- Inheritance mapping is stable across raw versus WDL split/merge joint-VCF
  representations under the complete-family-site policy.
- Every VCF and BED requiring random access is sorted, compressed, and indexed.
- Outputs contain sufficient configuration, version, reference, and QC metadata
  to reproduce and audit the run.
- Docker execution passes the synthetic end-to-end suite and Slurm/Apptainer
  profiles pass a representative smoke test.
- The CEPH parity gate has no unexplained scientific differences.
- BigWigs generated with `pyBigWig` are value-equivalent to the legacy
  `bedGraphToBigWig` tracks under the documented tolerance.
- No production code contains hard-coded sample IDs, input paths, chromosome
  loops, or hg38 downloads.

## Open decisions for later phases

The following decisions are intentionally deferred until the pedigree/model-only
MVP provides evidence for them:

- Whether to add direct Cromwell, miniwdl, or WOMtool metadata adapters. Any
  adapter should translate to the canonical manifest rather than create a second
  internal contract.
- How to represent and validate count-based methylation inputs alongside WDL
  model outputs.
- Whether schema version 2 should support arbitrary BED intervals, non-GRCh38
  references, or both.
- Whether a future workflow should permit subsets of the PED in the joint VCF;
  schema version 1 requires exact sample-set equality.
- Which additional stable WDL v3.x releases should move from explicit opt-in to
  the supported compatibility matrix after automated fixtures are available.
- Whether the final large all-CpG table should additionally be published in a
  columnar format such as Parquet after the TSV/BED-compatible contract is
  stable.
