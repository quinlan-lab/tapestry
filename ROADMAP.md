# Tapestry roadmap: generic WDL-driven workflow

**Status: implementation in progress.** The schema-v1 contract, validate-only
entry point, inheritance slice, model-only founder phasing, all-CpG expansion,
Docker runtime, and synthetic end-to-end/resume path are implemented. The real
CEPH v3.3.0/v3.3.1 parity gate, released image digest, CI, and representative
Slurm/Apptainer runs remain release blockers. Synthetic BigWig generator parity
passes against the checksum-pinned UCSC `bedGraphToBigWig` v2.10 binary.

## Goal

Make Tapestry a generic downstream consumer of the
[PacBio HiFi human WGS WDL](https://github.com/PacificBiosciences/HiFi-human-WGS-WDL).
Tapestry should ingest the WDL's published outputs instead of re-running HiPhase
or other upstream work. One versioned YAML file describes a run; a workflow
engine validates it and runs the remaining Tapestry stages.

## Consuming the PacBio WDL

Target the `family` entry point first. Despite its name, that WDL joint-calls a
family and runs HiPhase separately per sample; it does **not** perform Tapestry's
founder inheritance mapping and does not replace the trio path's pedigree-aware
WhatsHap/pedMEC step.

Outputs (in sample order) and their Tapestry use:

| WDL output | Tapestry use | Legacy work replaced |
| --- | --- | --- |
| `sample_ids` | Key per-sample arrays; reconcile with PED | Hard-coded sample discovery |
| `phased_small_variant_vcf` + index | Read-backed phased variants | `run-hiphase.sh` VCF work |
| `phase_blocks` | HiPhase block TSV for hap-map construction | `run-hiphase.sh` block work |
| `merged_haplotagged_bam` + index | Optional count pileup and visualization | HiPhase haplotagging |
| `cpg_combined_bed`, `cpg_hap1_bed`, `cpg_hap2_bed` + indexes | Upstream model methylation | Tapestry model-mode pb-CpG-tools |
| `joint_small_variants_vcf` + index | Inheritance mapping and final variant annotation; provisionally supported after validation/normalization below | Site-local joint-VCF path |
| `phase_stats`, `phase_haplotags` | Provenance and optional QC | Ad hoc HiPhase logs |

Canonical-manifest rules:

- Do not glob a WDL output directory. The MVP accepts a small, versioned,
  Tapestry-owned manifest that names joint artifacts and per-sample artifacts
  explicitly.
- Validate sample IDs, optional/null outputs, indexes, PED membership, WDL
  release, and artifact metadata before analysis.
- Direct Cromwell, miniwdl, Terra, or localized-output adapters are deferred.
  Future adapters translate engine output JSON into the canonical manifest and
  do not create additional internal contracts.
- Pin and record the exact WDL release and commit per run. Support v3.3.0 and
  v3.3.1 by default. Reject other releases unless a stable v3.x run explicitly
  opts into unaudited-release evaluation; never permit another major version or
  a prerelease through that opt-in.
  The supported tag commits are v3.3.0 at
  `db06f0af2354d847b971b0548eaade9ff145c912` and v3.3.1 at
  `477ef39ad69e86e90897ea7e313b86bfc12a2a96`; validation requires the matching
  tag/commit pair.

### Joint small-variant VCF compatibility

**Provisional decision (2026-08-14):** accept the PacBio `family` WDL v3.3.0 and
v3.3.1 `joint_small_variants_vcf` as MVP inputs for inheritance mapping and
final CpG variant annotation, subject to version-specific normalization and CEPH
parity gates below. The v3.3.1 artifact is structurally compatible, but it is not
a lossless copy of the original GLnexus joint callset; v3.3.0 must be checked
with the same source audit and fixture tests before release.

The source-level contract is suitable for Tapestry:

- The WDL's custom GLnexus configuration retains `GT`, `DP`, `AD`, `GQ`, and
  `PL`; `gtg` requires genotypes plus depth from `DP`, `AD`, or `SD` and a
  numeric `QUAL`.
- The published joint VCF is bgzip-compressed and has the `.tbi` index that
  `gtg-ped-map` and `gtg-concordance` specifically require.
- HiPhase changes genotype phase and adds `PS`/`PF`, but both `gtg` programs
  deliberately discard input phase and sort the allele indices before building
  or applying the inheritance map. Existing HiPhase phase orientation therefore
  does not affect inheritance inference.
- Tapestry's CpG annotation readers only need the requested sample columns,
  coordinates, `REF`/`ALT`, and `GT`; they do not depend on INFO annotations or
  the input phase orientation.

The adapter must nevertheless treat this as a transformed joint callset. In
v3.3.1, the WDL splits the GLnexus VCF into single-sample VCFs with
`--exclude-uncalled`, runs HiPhase per sample, and uses `bcftools merge` to
construct `joint_small_variants_vcf`. A split/merge reproduction on the
official `gtg` fixture changed 5,198 input records to 5,167 output records: 30
all-no-call records were removed, one additional net record was removed through
same-position merging, 32 variant keys were rewritten through allele trimming
or multiallelic merging, and `DP`/`AD` values for excluded `./.` sample calls
became missing. Most such records cannot pass `gtg`'s complete-family genotype
filter, but the changes can still affect depth statistics, marginal marker
selection, and CpG-overlap annotations.

Normalize and validate the joint VCF as follows:

1. Require the VCF sample set to equal the PED sample set exactly. Do not accept
   mere containment: extra samples participate in `gtg-concordance`'s missing
   genotype/depth filter.
2. Require a readable `.vcf.gz.tbi`, unique sample IDs, `GT`, numeric `QUAL`, and
   at least one supported depth field (`DP`, `AD`, or `SD`). Validate contigs
   against the configured reference and the supported chromosome policy.
3. Preserve the original WDL artifact for read-backed phasing and provenance.
   Derive a separate `gtg` input with HiPhase `FORMAT/PS` and `FORMAT/PF`
   removed. `gtg-concordance` replaces `GT` but otherwise preserves input FORMAT
   fields, so leaving `PS` would attach the old read-backed phase-set IDs to the
   new inheritance-oriented genotypes.
4. Do not silently add `FILTER=PASS`, biallelic-only, or SNV-only filtering.
   Current `gtg` behavior filters by numeric `QUAL` and allele length but ignores
   the VCF `FILTER` column; changing that policy could alter legacy results and
   requires a separate, measured decision.
5. Record preflight counts for records, samples, missing genotypes, non-PASS
   records, multiallelic sites, and records lacking usable depth. Publish them
   with the resolved config and tool versions.

The compatibility gate remains open until real CEPH family outputs from both
v3.3.0 and v3.3.1 are compared with the legacy joint VCF using identical `gtg`
versions and thresholds. Compare marker positions, inheritance-block boundaries,
concordance pass/fail/no-call counts, phased genotype alleles, and final
CpG-overlap/allele-specific annotations. Document expected differences rather
than requiring byte-identical VCFs.

Source audit: PacBio WDL
[`family.wdl`](https://github.com/PacificBiosciences/HiFi-human-WGS-WDL/blob/v3.3.1/workflows/family.wdl),
[`joint.wdl`](https://github.com/PacificBiosciences/HiFi-human-WGS-WDL/blob/v3.3.1/workflows/joint/joint.wdl),
[`glnexus.wdl`](https://github.com/PacificBiosciences/HiFi-human-WGS-WDL/blob/v3.3.1/workflows/wdl-common/wdl/tasks/glnexus.wdl), and
[`bcftools.wdl`](https://github.com/PacificBiosciences/HiFi-human-WGS-WDL/blob/v3.3.1/workflows/wdl-common/wdl/tasks/bcftools.wdl);
`gtg` commit
[`e12aca6`](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance/tree/e12aca6b49ee7208952467db4a2a9e2f79b98efb).

Methylation, count vs model:

- The WDL CpG task runs pb-CpG-tools with a model and publishes one
  combined/hap1/hap2 BED set; it does not publish Tapestry's count mode. Consume
  those BEDs as `model` data and preserve their header metadata.
- A later schema may make count-based methylation an explicit optional stage
  using the WDL haplotagged BAM; schema v1 is model-only and never manufactures
  or renames a count value.
- Thresholds differ: the WDL CpG task uses min MAPQ 1, min coverage 4; legacy
  Tapestry uses min coverage 10. Record effective upstream/downstream thresholds
  in provenance and test the effect on joins/QC.

BigWig generation:

- Use the already pinned `pyBigWig` dependency with chromosome lengths from the
  configured FASTA index instead of the legacy hard-coded hg38/
  `bedGraphToBigWig` path.
- Treat this as a generator migration. Compare chromosome metadata, intervals,
  missing intervals, and numeric values against legacy `bedGraphToBigWig`
  output with a documented floating-point tolerance before release. This gate
  passes on the deterministic fixture: the test checks chromosome names and
  lengths, interval coordinates, missing positions, and values to six decimal
  places against UCSC v2.10 (SHA-256
  `1a1527cf364e1e572a81c7284fc9ccd2b3690b5896baa5b57399864f85ad7771`).
  The pinned binary is a containerized parity oracle; production output uses
  `pyBigWig`.

## Pipeline boundary

The first generic implementation starts *after* a completed WDL run. Do not nest
miniwdl/Cromwell inside Nextflow in the MVP; launching the WDL is a separate
execution/provenance problem for later.

```text
Tapestry canonical WDL-output manifest
        v
VALIDATE_INPUTS
        v
inheritance/pedMEC phase
        v
build hap maps and phase methylation
        v
expand/annotate all CpGs
        v
publish BED, BigWig, QC, provenance
```

For pedigree, the WDL replaces read-backed phasing and model methylation
generation; Tapestry still runs `gtg-ped-map`/`gtg-concordance`, founder phasing,
and all-CpG expansion. For trio, explicitly decide whether to keep
WhatsHap/pedMEC — do not substitute per-sample HiPhase output for pedigree-aware
phasing merely because both carry `PS`/`HP` labels.

## Run format and interface

Canonical human-authored format is **YAML** validated by a versioned JSON Schema
(passed to Tapestry via `--run-config`; JSON is accepted as a machine
equivalent). Tapestry uses its own parameter instead of Nextflow's
`-params-file` so relative paths can be resolved against the configuration
file's directory rather than the caller's working directory.
Use ordinary YAML types only — no custom tags, executable expressions, env
lookups, or secrets. Keep two concerns separate:

- the run YAML holds scientific inputs, upstream manifest, reference, pedigree,
  requested outputs, and reproducibility-affecting thresholds;
- `nextflow.config` and profiles hold executor, queue, container, retry, and
  site resource policy.

MVP schema version 1:

```yaml
schema_version: 1
mode: pedigree
project: {id: family_001, outdir: results/family_001}
pedigree: {ped: inputs/family.ped}
reference:
  name: GRCh38
  fasta: references/hg38.analysisSet.fa
  fasta_index: references/hg38.analysisSet.fa.fai
upstream:
  manifest: inputs/pacbio-wdl.tapestry.json
  # Required only to evaluate another stable v3.x release.
  allow_unaudited_release: false
samples:
  include: [child_01, child_02] # optional; omit for all eligible PED members
inheritance:
  method: gtg
  map: {min_qual: 20, min_depth: 10, min_run_markers: 10}
  concordance: {min_qual: 20, min_depth: 5}
methylation:
  modes: [model]
  min_coverage: 10
  mismatch_window_bp: 50
regions:
  include: [chr1, chr2, chr3] # optional; omit for chr1-chr22
outputs:
  bigwig: true
```

Schema-v1 rules:

- Reject unknown keys, custom YAML tags, and unsupported schema versions.
- Accept pedigree mode, model methylation, GRCh38, and whole autosomes `chr1`
  through `chr22` only. Reject chrX, chrY, mitochondrial/alternate contigs, and
  interval strings.
- Resolve run-config paths relative to the run-config file. Resolve canonical
  manifest paths relative to the manifest file.
- Treat an omitted `samples.include` as all eligible PED members; reject an
  explicitly empty list or an ineligible member.
- Support WDL v3.3.0 and v3.3.1 by default. Reject other releases unless
  `upstream.allow_unaudited_release` is true and the manifest names another
  stable v3.x release. Always reject prereleases and other major versions.

The file named by `upstream.manifest` uses this Tapestry-owned contract; it is
not raw Cromwell, miniwdl, or Terra output JSON:

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

The manifest may additionally carry provenance-only phase statistics,
haplotags, haplotagged BAMs, and original workflow-engine metadata. Required
joint artifacts and artifacts for selected samples may not be null. Artifacts
needed only by unselected samples may be null.

Single entry point is a Nextflow DSL2 pipeline:

```bash
nextflow run quinlan-lab/tapestry -r <tag> -profile slurm,apptainer \
  --run-config family.yaml -work-dir /path/to/work -resume
# validate-only:
nextflow run . -entry validate --run-config family.yaml
```

Single-dash args are Nextflow's (`-profile`, `-resume`, `-work-dir`, `-r`,
and `-entry`); `--run-config` is Tapestry's pipeline parameter. Output selection
and scientific thresholds live in YAML; queue, memory, container runtime, and
account live in profiles. The pipeline publishes resolved inputs, the canonical
sample manifest, tool/container versions, QC, and scientific outputs under
`project.outdir`. Nextflow writes its trace/report/timeline/DAG to
`.nextflow-reports/` under the launch directory by default; those lifecycle
reports are not entries in `results-manifest.json` and may be redirected with
Nextflow's standard `-with-*` options.

This section is the source of truth for the public schema and invocation. The
implementation plan in `impl.md`, examples, and eventual JSON Schema must remain
consistent with it.

## Design rules

1. The config (and/or PED) is the source of truth for sample lists and roles.
   Never add hard-coded CEPH/`NA128*` IDs, `GRCh38` labels, `/scratch` paths, or
   `chr{1..22}` lists to workflow code.
2. Resolve relative input paths relative to the config file, not the caller's
   CWD. Keep run products under the configured output root unless overridden.
3. Validate early, before expensive jobs: schema/version, supported WDL
   provider/release, required keys and equal array lengths, unique sample IDs,
   PED consistency, VCF sample presence, input/index existence, reference-contig
   compatibility, autosome-only regions, output collisions, supported modes, and
   the WDL release policy. Report errors with field name and offending value; do
   not partially start on an invalid config.
4. Keep algorithms in importable Python functions on explicit
   arguments/dataframes. Nextflow owns channels, retries, containers, publishing;
   small Python components own schema validation, canonical-manifest validation,
   and transformations.
5. Prefer argv lists with `subprocess.run(..., check=True)`; avoid building shell
   strings from config values. `src/util/shell.py` uses `shell=True` and must not
   receive untrusted config text without refactoring.
6. Preserve existing Python CLIs as low-level tools; the config-driven entry may
   call them or their extracted functions rather than duplicate their logic.
7. Keep current production/dev behavior working, but mark `run-hiphase.sh`,
   model-mode `aligned_bam_to_cpg_scores.sh`, and the WIP `Snakefile` as legacy
   rather than wrapping them into the new DAG.
8. Put a portable example config and schema docs under `examples/`. Do not commit
   credentials, private data locations, or generated genome-wide outputs.
9. Record the normalized config, command/tool versions, and run metadata with
   outputs.

Add tests proving equivalent YAML and JSON produce the same run model and that a
realistic canonical manifest resolves into the expected sample model. Cover
failure cases including missing/null artifacts, sample-set errors, inconsistent
pedigrees, reference mismatch, non-autosomal regions, unsafe/colliding outputs,
unsupported WDL releases, and unaudited-release opt-in behavior.

## Open decisions

These decisions are deferred beyond schema v1 or remain release gates. Resolve
them through fixtures and small end-to-end tests before expanding the public
schema or declaring the MVP complete:

1. **Additional PacBio WDL releases.** v3.3.0 and v3.3.1 are the initial support
   targets. Decide which other stable v3.x releases can move from explicit
   opt-in to the default compatibility matrix after version-specific fixtures
   and parity runs exist.
2. **Joint-VCF parity gate.** Complete the real CEPH comparisons for v3.3.0 and
   v3.3.1. Record fixtures, commands, WDL and `gtg` versions, thresholds, result
   deltas, and final dispositions here before advertising either release as
   supported.
3. **Engine adapters.** Decide which direct Cromwell, miniwdl, or Terra adapters
   are worth adding after the canonical-manifest MVP. Every adapter must emit the
   existing Tapestry manifest contract.
4. **Count-mode follow-on.** Decide how count-based methylation inputs are added
   after the model-only MVP. The MVP retains legacy count columns as typed nulls.
5. **Threshold harmonization beyond the MVP.** Determine whether future modes
   expose WDL model data at minimum coverage 4, the Tapestry-filtered view at
   coverage 10, or both. Never obscure the upstream threshold recorded in the
   BED header.
6. **Trio support.** Decide whether the generic product retains the
   WhatsHap/pedMEC trio path, converges trio and pedigree results on one output
   schema, or releases pedigree support first and treats trio as a separately
   versioned follow-on.
7. **Reference and region expansion.** The MVP is GRCh38 chr1-chr22 only. Decide
   when to add sex chromosomes, non-GRCh38 references, arbitrary intervals, and
   alternate contigs.
8. **Execution and distribution.** Confirm Nextflow DSL2 as the public
   orchestrator after the first vertical slice, choose the supported container
   runtime/profile matrix, and decide whether releases are run from GitHub,
   local bundles, or both.

Record each resolution here with its rationale, date, and the fixture/test that
supports it; then move durable conclusions into `AGENTS.md`, the schema, or
user documentation as appropriate.

## Implementation status

1. **Freeze the upstream contract — partial.** The canonical manifest and
   supported tag/commit pairs are frozen. De-identified real-output fixtures
   from WDL `family`
   v3.3.0 and v3.3.1; document each consumed canonical-manifest field and its
   optionality. Include the original joint VCF, its derived PS/PF-free `gtg`
   input, and their preflight summaries.
2. **Define and validate inputs — implemented.** YAML example, JSON Schema, file-relative path
   resolution, PED/reference checks, release policy, and canonical-manifest
   validation. The `validate` entry point works before any scientific stage is
   ported.
3. **Make CLIs portable — implemented for the generic path.** Remove hard-coded `hg38`, sample naming, and
   output-path assumptions while preserving dataframe logic and CLI
   compatibility.
4. **Build the pedigree/model-only slice — implemented.** Ingest WDL phased VCFs, phase blocks,
   haplotagged BAM metadata, model CpG BEDs; run inheritance mapping, founder
   phasing, all-CpG expansion; publish provenance and QC.
5. **Prove equivalence — release blocker.** The deterministic BigWig-generator
   comparison is complete. Run the remaining scientific comparisons on the CEPH fixture/data; compare joint-VCF marker
   sites, hap-map blocks, concordance pass/fail/no-call counts, phased genotype
   alleles, mismatch sites, CpG variant annotations, model methylation, row
   counts, null patterns, coordinates, QC, and BigWig chromosome/interval/value
   content to legacy. Explain differences from WDL version, VCF
   split/merge behavior, coverage thresholds, or generator semantics.
6. **Add optional count mode.** Generate count BEDs from WDL haplotagged BAMs,
   merge with model BEDs; verify disabling count yields a valid model-only schema.
7. **Port or retire the trio branch deliberately.** Decide whether to keep
   WhatsHap/pedMEC and the trio output schema. Do not let this block the pedigree
   MVP.
8. **Release — pending.** Local/test and CHPC Slurm+Apptainer profiles, containers,
   `-resume` tests, user docs, a tagged release.

**Pedigree MVP is complete** when a user can point one YAML file at a pinned WDL
output manifest, run validate then the main Nextflow command, and get
founder-phased model methylation plus all-CpG/QC outputs without invoking
`run-hiphase.sh` or editing repo scripts. Count mode, WDL launching, and trio
support are follow-ons.
