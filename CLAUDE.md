# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Tapestry is a bioinformatics pipeline that phases DNA methylation from PacBio HiFi reads in a human pedigree to the haplotypes of the pedigree's founders. It supports two workflows:

- **Pedigree-wise workflow** (large pedigrees): Uses HiPhase for read-backed phasing + `gtg-ped-map`/`gtg-concordance` for inheritance-based phasing, then combines both to build a haplotype map that phases methylation to founder haplotypes.
- **Trio workflow**: Uses WhatsHap for pedigree-aware phasing of a trio (child + parents), directly mapping hap1=paternal, hap2=maternal.

## Environment Setup

```bash
python3.11 -m venv .venv
source .venv/bin/activate

# cyvcf2 must be installed from source (see README.md for details)
pip install -r requirements.txt
```

Python 3.11 is required. Pyright is configured via `pyrightconfig.json` with `extraPaths: ["src", "src/util"]`.

## External CLI Dependencies

`bedGraphToBigWig`, `bgzip`, `tabix`, `bcftools`, `gtg-ped-map`, `gtg-concordance`, `hiphase`, `aligned_bam_to_cpg_scores`, `whatshap`

## Architecture

### Pipeline Steps (shell scripts at repo root)

1. **Variant phasing**: `run-hiphase.sh` (pedigree) or `run-whatshap-phase.sh` (trio) — phase variants from BAMs
2. **Inheritance-based haplotype map** (pedigree only): `build-iht-based-haplotype-map-and-phase-variants.sh`
3. **Methylation extraction**: `aligned_bam_to_cpg_scores.sh` / `aligned_bam_to_cpg_scores.trio.sh` — generate count- and model-based methylation from haplotagged BAMs
4. **Phase methylation to founder/parent haplotypes**: `phase_meth_to_founder_haps.sh` / `phase_meth_to_parent_haps.sh`
5. **Expand to all CpGs**: `expand_to_all_cpgs.sh` / `expand_to_all_cpgs.trio.sh` — add reference CpGs, unphased methylation, allele-specific CpG flags

Shell scripts accept `--dev-dir` for local development with small test data in `trio-dev-data/`.

### Core Python Modules (`src/`)

- `phase_meth_to_founder_haps.py` — Main entry point for pedigree workflow step 4. Orchestrates all phasing logic and bigwig output.
- `phase_meth_to_parent_haps.py` — Main entry point for trio workflow step 4. Uses WhatsHap phase blocks.
- `expand_to_all_cpgs.py` — Main entry point for step 5. Joins CpG reference sites with phased/unphased methylation, labels allele-specific CpGs via SNV overlap. Memory-intensive; uses jemalloc tuning (`MALLOC_CONF`) and explicit `gc.collect()`.
- `get_all_phasing.py` — Reads read-backed (cyvcf2) and inheritance-based phasing data, joins them via bioframe overlap.
- `get_hap_map.py` — Builds the haplotype map by comparing bit vectors (read-backed vs inheritance-backed allele sequences) to determine which founder haplotype corresponds to hap1/hap2.
- `get_meth_hap1_hap2.py` — Reads pb-CpG-tools BED output for hap1/hap2 methylation (count and model modes).

### Utilities (`src/util/`)

- `shell.py` — Subprocess wrapper (uses `/usr/bin/bash`)
- `read_data.py` / `write_data.py` — BED file I/O with header convention (`##source=...` metadata, `#`-prefixed column headers)
- `version_sort.py` — Chromosome-aware sorting (chr1, chr2, ..., chr22, chrX, chrY)
- `remove_funky_chromosomes.py` — Filters out non-standard contigs (`_random`, `chrEBV`, `chrUn_`)
- `compress-index-vcf`, `sort-compress-index-bed` — Shell utilities for bgzip/tabix

## Key Libraries and Patterns

- **Polars** (`pl`) is the primary DataFrame library (not pandas). Bioframe (`bf.overlap`, `bf.closest`) is used for genomic interval operations but requires pandas conversion; there is a TODO to migrate to polars-bio.
- **cyvcf2** is used to parse VCF files in `get_all_phasing.py`.
- Data flows through Polars DataFrames. Conversions to/from pandas happen only at bioframe boundaries.
- BED output files use a two-line header: `##source='...'` and `#column_names`.

## Development Data

`trio-dev-data/` contains small BAMs and VCF for the NA12878 trio (chr21 subset). Run `create_dev_data.sh` to regenerate. Shell scripts accept `--dev-dir trio-dev-data` to use this data instead of production paths.

## Snakemake

`Snakefile` is a WIP stub for automating the pipeline. Not yet functional.
