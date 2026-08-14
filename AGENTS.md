# Tapestry repository guide

Tapestry phases PacBio HiFi CpG methylation onto inherited haplotypes. It
combines variant phasing, haplotagged reads, and per-CpG methylation calls to
produce BED/BigWig outputs in which methylation is traced to a founder or parent
haplotype.

The repo currently targets the CEPH 1463/Palladium data set and GRCh38. Many
top-level shell scripts contain CHPC `/scratch/...` paths, fixed filename
patterns, tool locations, sample IDs, chromosome lists, and resource settings.
Treat those as the current worked deployment, not portable API contracts.

The planned next step is to make Tapestry a generic downstream consumer of the
[PacBio HiFi human WGS WDL](https://github.com/PacificBiosciences/HiFi-human-WGS-WDL),
driven by one versioned YAML run file. That work is a **proposal, not yet
implemented** — see [ROADMAP.md](ROADMAP.md). This guide describes the repo as it
is today.

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

For Python changes, also run `.venv/bin/pyright`. For changed shell entry points,
run `bash -n <script>`. For workflow changes, prefer the smallest
`trio_dev_data` region that exercises the modified stages (external tools and the
reference fixture must be installed first). Do not launch the hard-coded
production scripts or large background jobs merely as a test.

For scientific transformations, test invariants rather than only row snapshots:
unique phase-set assignment, no duplicate hap-map blocks, stable coordinate
conventions, expected null preservation, correct haplotype labels, and
sorted/indexable output.
