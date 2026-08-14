# Tapestry roadmap: generic WDL-driven workflow

**Status: proposal, not yet implemented.** This is a design target, not a
description of current behavior. See [AGENTS.md](AGENTS.md) for the repo as it is.

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

Adapter rules:

- Do not glob a WDL output directory. Read the engine's output JSON/manifest,
  read `sample_ids`, zip arrays by index, emit a small canonical per-sample
  manifest. Validate array lengths, sample order, optional/null outputs, indexes,
  and PED membership before analysis.
- Keep the adapter isolated so Cromwell, miniwdl, Terra, or a localized copy work
  without changing scientific code.
- Pin and record the exact WDL release and output schema per run (start with a
  supported release such as `v3.3.1`); keep release-specific field aliases in the
  adapter, not throughout Tapestry.

### Joint small-variant VCF compatibility

**Provisional decision (2026-08-14):** accept the PacBio `family` WDL v3.3.1
`joint_small_variants_vcf` as the MVP input for inheritance mapping and final
CpG variant annotation, subject to the normalization and CEPH parity gate below.
It is structurally compatible, but it is not a lossless copy of the original
GLnexus joint callset.

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

The compatibility gate remains open until one real CEPH family output from the
pinned WDL release is compared with the legacy joint VCF using identical `gtg`
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
- Make count-based methylation an explicit optional stage using the WDL
  haplotagged BAM; never manufacture or rename a count column.
- Thresholds differ: the WDL CpG task uses min MAPQ 1, min coverage 4; legacy
  Tapestry uses min coverage 10. Record effective upstream/downstream thresholds
  in provenance and test the effect on joins/QC.

## Pipeline boundary

The first generic implementation starts *after* a completed WDL run. Do not nest
miniwdl/Cromwell inside Nextflow in the MVP; launching the WDL is a separate
execution/provenance problem for later.

```text
PacBio family-WDL output JSON
        v
INGEST_WDL_OUTPUTS + VALIDATE_INPUTS
        |
        +----------------+---------------- (optional) count CpG pileup
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
(Nextflow loads it via `-params-file`; JSON accepted as a machine equivalent).
Use ordinary YAML types only — no custom tags, executable expressions, env
lookups, or secrets. Keep two concerns separate:

- the run YAML holds scientific inputs, upstream manifest, reference, pedigree,
  requested outputs, and reproducibility-affecting thresholds;
- `nextflow.config` and profiles hold executor, queue, container, retry, and
  site resource policy.

Proposed minimum (not a settled interface):

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
  provider: pacbio_hifi_human_wgs_wdl
  release: v3.3.1
  outputs_manifest: inputs/pacbio-wdl.outputs.json
inheritance: {method: gtg}
methylation:
  model: {source: upstream}
  count: {enabled: true, min_coverage: 10, min_mapq: 1}
regions:
  include: [chr1, chr2, chr3]   # optional; omit to use input-supported contigs
```

Single entry point is a Nextflow DSL2 pipeline:

```bash
nextflow run quinlan-lab/tapestry -r <tag> -profile slurm,apptainer \
  -params-file family.yaml -work-dir /path/to/work -resume
# validate-only:
nextflow run main.nf -entry validate -params-file family.yaml
```

Single-dash args are Nextflow's (`-profile`, `-resume`, `-work-dir`, `-r`,
`-params-file`). Output selection and scientific thresholds live in YAML; queue,
memory, container runtime, and account live in profiles. The pipeline publishes a
resolved YAML/JSON snapshot, the canonical sample manifest, tool/container
versions, a trace/timeline, and final outputs under `project.outdir`.

## Design rules

1. The config (and/or PED) is the source of truth for sample lists and roles.
   Never add hard-coded CEPH/`NA128*` IDs, `GRCh38` labels, `/scratch` paths, or
   `chr{1..22}` lists to workflow code.
2. Resolve relative input paths relative to the config file, not the caller's
   CWD. Keep run products under the configured output root unless overridden.
3. Validate early, before expensive jobs: schema/version, supported WDL
   provider/release, required keys and equal array lengths, unique sample IDs,
   PED consistency, VCF sample presence, input/index existence, reference-contig
   compatibility, output collisions, supported modes. Report errors with field
   name and offending value; do not partially start on an invalid config.
4. Keep algorithms in importable Python functions on explicit
   arguments/dataframes. Nextflow owns channels, retries, containers, publishing;
   small Python adapters own schema validation, WDL-manifest normalization, and
   transformations.
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

Add tests proving equivalent YAML and JSON produce the same run model, that
realistic WDL output JSON normalizes into the expected sample manifest, and for
failure cases: missing/null artifacts, mismatched array lengths, sample-order
errors, inconsistent pedigrees, reference mismatch, unsafe/colliding outputs.

## Open decisions

Resolve these through fixtures and small end-to-end tests before freezing the
public schema or declaring the MVP complete:

1. **Initial PacBio WDL release.** Choose the first supported release, capture a
   representative output manifest, and pin its exact tag/commit and submodule
   versions. `v3.3.1` is currently an example, not yet the compatibility floor.
2. **Joint-VCF parity gate.** The source audit and `gtg` fixture split/merge test
   support `joint_small_variants_vcf` as the provisional v3.3.1 MVP input; the
   remaining decision is whether the real CEPH comparison meets the compatibility
   gate defined above. Record the fixture, commands, `gtg` version, thresholds,
   result deltas, and final disposition here before calling it fully supported.
3. **Manifest dialect.** Decide which engine output JSON is canonical for the
   MVP (for example Cromwell or miniwdl) and whether other engines are handled
   by adapter modes or first normalized into a small Tapestry-owned manifest.
4. **Count-mode default and schema.** Decide whether count mode defaults off for
   the model-only MVP, how model-only outputs represent absent count columns,
   and whether the legacy minimum coverage of 10 remains the default when count
   mode is enabled.
5. **Threshold harmonization.** Determine whether WDL model data at minimum
   coverage 4 is consumed unchanged, filtered to a Tapestry threshold, or
   accompanied by both raw and harmonized QC. Never obscure the upstream
   threshold recorded in the BED header.
6. **Trio support.** Decide whether the generic product retains the
   WhatsHap/pedMEC trio path, converges trio and pedigree results on one output
   schema, or releases pedigree support first and treats trio as a separately
   versioned follow-on.
7. **Reference and region policy.** Define supported reference requirements,
   contig-name compatibility checks, behavior for non-GRCh38 references, and
   whether region-restricted runs are development-only or a supported public
   feature.
8. **Execution and distribution.** Confirm Nextflow DSL2 as the public
   orchestrator after the first vertical slice, choose the supported container
   runtime/profile matrix, and decide whether releases are run from GitHub,
   local bundles, or both.

Record each resolution here with its rationale, date, and the fixture/test that
supports it; then move durable conclusions into `AGENTS.md`, the schema, or
user documentation as appropriate.

## Implementation sequence

1. **Freeze the upstream contract.** De-identified fixture from one pinned WDL
   `family` release; document each consumed field and its optionality. Include
   the original joint VCF, its derived PS/PF-free `gtg` input, and their preflight
   summaries.
2. **Define and validate inputs.** YAML example, JSON Schema, path resolution,
   PED/reference checks, WDL-output adapter writing the canonical manifest. The
   `validate` entry point works before any scientific stage is ported.
3. **Make CLIs portable.** Remove hard-coded `hg38`, sample naming, and
   output-path assumptions while preserving dataframe logic and CLI
   compatibility.
4. **Build the pedigree/model-only slice.** Ingest WDL phased VCFs, phase blocks,
   haplotagged BAM metadata, model CpG BEDs; run inheritance mapping, founder
   phasing, all-CpG expansion; publish provenance and QC.
5. **Prove equivalence.** Run on the CEPH fixture/data; compare joint-VCF marker
   sites, hap-map blocks, concordance pass/fail/no-call counts, phased genotype
   alleles, mismatch sites, CpG variant annotations, model methylation, row
   counts, null patterns, coordinates, and QC to legacy. Explain differences
   from WDL version, VCF split/merge behavior, or coverage thresholds.
6. **Add optional count mode.** Generate count BEDs from WDL haplotagged BAMs,
   merge with model BEDs; verify disabling count yields a valid model-only schema.
7. **Port or retire the trio branch deliberately.** Decide whether to keep
   WhatsHap/pedMEC and the trio output schema. Do not let this block the pedigree
   MVP.
8. **Release.** Local/test and CHPC Slurm+Apptainer profiles, containers,
   `-resume` tests, user docs, a tagged release.

**Pedigree MVP is complete** when a user can point one YAML file at a pinned WDL
output manifest, run validate then the main Nextflow command, and get
founder-phased model methylation plus all-CpG/QC outputs without invoking
`run-hiphase.sh` or editing repo scripts. Count mode, WDL launching, and trio
support are follow-ons.
