# Tapestry repository guide

Tapestry phases PacBio HiFi CpG methylation onto inherited haplotypes. It
combines variant phasing, haplotagged reads, and per-CpG methylation calls to
produce BED/BigWig outputs in which methylation is traced to a founder or parent
haplotype.

The legacy scripts target the CEPH 1463/Palladium data set and GRCh38. Many
top-level shell scripts contain CHPC `/scratch/...` paths, fixed filename
patterns, tool locations, sample IDs, chromosome lists, and resource settings.
Treat those as a worked deployment, not portable API contracts.

Tapestry now includes a generic downstream workflow for the
[PacBio HiFi human WGS WDL](https://github.com/PacificBiosciences/HiFi-human-WGS-WDL),
driven by one versioned YAML run file and a canonical WDL-output manifest. See
[`ROADMAP.md`](ROADMAP.md) for the schema and remaining release gates and
[`impl.md`](impl.md) for the implementation/acceptance plan. The legacy
site-specific shell paths remain available but are not the generic API.

### Generic pedigree/model-only path

`main.nf` validates the contract, normalizes the family VCF, runs pinned
`gtg-ped-map`/`gtg-concordance`, filters WDL model BEDs, reconciles HiPhase and
inheritance phase per selected sample, generates reference-autosome CpGs, and
publishes all-CpG tables, BigWigs, QC, provenance, and a results manifest. The
generic path is GRCh38 `chr1`-`chr22` only and supports WDL v3.3.0/v3.3.1 by
default. Another stable v3.x release requires an explicit unaudited-release
opt-in; prereleases and other major versions are rejected.

## Current workflows

The full-pedigree and trio paths share the methylation-processing concepts but
use different phasing sources.

### Pedigree-wise

1. `run-hiphase.sh` — read-backed phasing and haplotagging per sample (HiPhase).
2. `build-iht-based-haplotype-map-and-phase-variants.sh` — `gtg-ped-map`/
   `gtg-concordance` derive inheritance blocks and a pedigree-wide,
   inheritance-phased VCF.
3. `aligned_bam_to_cpg_scores.sh` — pb-CpG-tools in `count` and `model` modes on
   haplotagged BAMs.
4. `phase_meth_to_founder_haps.sh` → `src/phase_meth_to_founder_haps.py`.
   Reconciles read-backed and inheritance-backed phase, builds haplotype-map
   blocks, records mismatch sites, and re-labels hap1/hap2 methylation with
   parental and founder haplotypes.
5. `expand_to_all_cpgs.sh` → `src/expand_to_all_cpgs.py`. Joins phased and
   unphased methylation to the full reference CpG set, annotates CpGs affected by
   variants, computes QC, writes final BED/BigWig.

### Trio-wise

1. `run-whatshap.sh` — pedigree-aware WhatsHap/pedMEC phasing, then haplotags
   child, father, mother BAMs.
2. `aligned_bam_to_cpg_scores.sh -t ...` — count/model methylation for the three
   samples.
3. `phase_meth_to_parent_haps.sh` → `src/phase_meth_to_parent_haps.py`. Compares
   the child's phased alleles with each parent's, builds paternal/maternal
   haplotype maps, assigns methylation to parent haplotypes.
4. `expand_to_all_cpgs.trio.sh` → `src/expand_to_all_cpgs_trio.py`. Adds all
   reference/sample CpGs, per-member unphased methylation, allele-specific
   annotations, QC, visualization tracks.

`trio_dev_data/` is a small checked-in example used by the trio `--dev-dir`
paths; its data-creation scripts still depend on site-local source data.

## Repository map

- `main.nf`, `nextflow.config`: generic DSL2 workflow and local/Docker/
  Apptainer/Slurm profiles. Scientific inputs come from `--run-config`; executor
  and resource policy belong in profiles.
- `schemas/`, `examples/generic/`: authoritative machine contracts and portable
  human-authored examples. `ROADMAP.md` remains the source of truth for schema
  intent.
- `src/tapestry_validate.py`: strict YAML/JSON, manifest, PED, reference, VCF,
  HiPhase-table, model-BED, release, and output-collision validation.
- `src/normalize_joint_vcf.py`, `src/run_gtg_inheritance.py`: deterministic
  all-site/complete-family VCF branches and pinned `gtg` orchestration.
- `src/filter_model_beds.py`, `src/generate_reference_cpgs.py`,
  `src/expand_model_to_all_cpgs.py`: generic model-only filtering and all-CpG
  publication.
- `Dockerfile`, `requirements-pipeline.txt`, `containers/gtg.Cargo.lock`: pinned
  generic runtime; do not update the `gtg` commit or Cargo lock independently.
- `README.md`: authoritative user-facing description, workflow details, and
  final output-column definitions.
- `docs/pedigree_workflow.mmd`: pedigree workflow/data-flow diagram.
- `src/phasing_pedigree.py`, `src/hap_map_pedigree.py`: reconcile HiPhase and
  inheritance phasing; build founder haplotype maps.
- `src/phasing_trio.py`, `src/hap_map_trio.py`: parse per-sample phase sets and
  build paternal/maternal maps for a trio.
- `src/phase_meth_to_founder_haps.py`, `src/phase_meth_to_parent_haps.py`: main
  phase-to-methylation programs.
- `src/expand_to_all_cpgs.py`, `src/expand_to_all_cpgs_trio.py`: final expansion,
  variant annotation, QC, output generation.
- `src/get_meth_hap1_hap2.py`: combines hap1/hap2 pb-CpG-tools BEDs.
- `src/write_all_cpgs_in_reference.py`: enumerates reference CpGs.
- `src/util/`: shared BED/VCF I/O, hap-map output, sorting/indexing, logging
  (`logging.sh` for shell, `logging_util.py` for Python), phase helpers,
  diagnostics.
- Top-level `*.sh`: current end-to-end orchestration and deployment settings.
- `tests/test_recombination_dedup.py`: regression test for nested phase sets and
  duplicate hap-map blocks.
- `wiki/`: code-linked conceptual walkthrough with reproducible figures.
- `manuscript/`: paper text, analyses, figure sources; not pipeline code.
- `Snakefile`: an incomplete experiment, not the working workflow or a source of
  truth.

## Data contracts and invariants

- BED coordinates are zero-based, half-open; VCF positions are one-based.
  Preserve explicit conversions at file boundaries.
- Chromosome ordering must be version/natural, not lexical; use the existing
  version-sort and sort/compress/index helpers.
- bgzipped BED/VCF outputs generally need tabix/CSI indexes. Do not silently
  change established output columns or filenames without updating downstream
  scripts and `README.md`.
- A variant belongs to the phase set named by that sample's VCF `FORMAT/PS`
  value. Do not infer membership only by overlapping block spans: block spans can
  nest. Each informative SNV must contribute to exactly one `(child phase set,
  parent phase set)` group. The recombination regression test protects this.
- Trio labels: child hap1 = paternal, child hap2 = maternal; father hap1/hap2 =
  A/B; mother hap1/hap2 = C/D. A/B/C/D identify fixed parental haplotypes, not
  transmitted/non-transmitted status.
- Pedigree path: hap1/hap2 are read-backed labels that must be reconciled with
  paternal/maternal inheritance labels before attaching founder IDs.
- Retain unphasable/uncovered CpGs with null measurements where the all-CpG data
  model calls for them. A left/full join is often biologically meaningful, not an
  accident.
- Bit-vector disagreement sites are QC outputs feeding the `within 50 bp of
  mismatch` annotations. Keep paternal and maternal mismatch tracks separate in
  the trio workflow.
- Count-based and model-based methylation are distinct measurements. Never
  conflate them.
- Large inputs are normal. Prefer Polars/lazy or streaming operations, avoid
  unnecessary pandas conversions, and do not materialize genome-wide Cartesian
  joins.

## Development conventions

- Python 3.11; the repo expects a local `.venv` with deps pinned in
  `requirements.txt`. This repo predates uv adoption and has no `pyproject.toml`;
  use the `.venv` form here rather than migrating it.
- Imports rely on both `src` and `src/util` being on `PYTHONPATH`:

  ```bash
  PYTHONPATH=src:src/util .venv/bin/python <script> ...
  ```

- Use Polars for pipeline tables. Keep any bioframe/pandas conversions (interval
  or BigWig APIs) localized.
- Follow existing logging: announce inputs/arguments, dataframe sizes, output
  paths, and completion.
- Reuse `src/util` readers/writers, hap-map writers, and sorting/indexing helpers
  instead of introducing subtly different formats.
- Shell entry points: add `-h/--help`, quote paths, fail on unknown arguments,
  fail fast on missing prerequisites.
- Keep generated data, indexes, logs, and large genomic outputs out of git.
  Small deterministic fixtures belong in `tests/` or `examples/` (note
  `examples/` already holds `apex_trio_relabeling/`).
- Update `README.md` and example configs whenever user-facing inputs, output
  schemas, labels, or workflow commands change. Update wiki/manuscript only when
  the corresponding scientific explanation or figure source changes.

## Verification

Run the checks that match the change. The lightweight regression test needs no
production data:

```bash
PYTHONPATH=src:src/util .venv/bin/python tests/test_recombination_dedup.py
```

The generic suite is designed for the pipeline container:

```bash
docker run --rm -v "$PWD:/work" -w /work \
  -e PYTHONDONTWRITEBYTECODE=1 \
  -e PYTHONPATH=/work/src:/work/src/util \
  tapestry:dev python -m unittest discover -v -s tests -p 'test_*.py'

docker run --rm -v "$PWD:/work" -w /work \
  -e PYTHONDONTWRITEBYTECODE=1 \
  -e PYTHONPATH=/work/src:/work/src/util \
  tapestry:dev python tests/test_recombination_dedup.py
```

For workflow changes, run both `-entry validate` and the synthetic informative
family through the Docker profile, then repeat with `-resume` and confirm every
task is cached. Real WDL v3.3.0 and v3.3.1 CEPH fixtures and the
scientific parity comparison remain mandatory release gates even when synthetic
tests pass. The deterministic `pyBigWig`/`bedGraphToBigWig` migration test runs
against the checksum-pinned UCSC binary included in the container and must not
be skipped in container test runs.

For Python changes, also run `.venv/bin/pyright`. For changed shell entry points,
run `bash -n <script>`. For workflow changes, prefer the smallest
`trio_dev_data` region that exercises the modified stages (external tools and the
reference fixture must be installed first). Do not launch the hard-coded
production scripts or large background jobs merely as a test.

For scientific transformations, test invariants rather than only row snapshots:
unique phase-set assignment, no duplicate hap-map blocks, stable coordinate
conventions, expected null preservation, correct haplotype labels, and
sorted/indexable output.
