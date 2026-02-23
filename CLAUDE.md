# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**tapestry** is a bioinformatics pipeline that phases DNA methylation from HiFi reads in a human pedigree to the haplotypes of the pedigree's founders. It is designed for the CEPH1463 Platinum Pedigree dataset on the UCGD cluster.

## Environment Setup

```bash
/path/to/python3.11 -m venv .venv
source .venv/bin/activate

# cyvcf2 must be built from source with external htslib
cd $HOME
git clone --recursive https://github.com/brentp/cyvcf2
cd cyvcf2
pip install -r requirements.txt
CYVCF2_HTSLIB_MODE=EXTERNAL python setup.py install

cd /path/to/tapestry
pip install -r requirements.txt
```

Type checking uses Pyright (basic mode, Python 3.11, venv at `.venv`).

External CLI tools required in PATH: `bedGraphToBigWig`, `bgzip`, `tabix`, `bcftools`, `gtg-ped-map`, `gtg-concordance`, `hiphase`, `aligned_bam_to_cpg_scores`, `whatshap`, `samtools`.

## Pipeline Workflow

Five sequential shell scripts drive the pipeline:

1. `run-hiphase.sh` — read-backed phasing with HiPhase
2. `build-iht-based-haplotype-map-and-phase-variants.sh` — inheritance-based haplotype map
3. `aligned_bam_to_cpg_scores.sh` — count- and model-based methylation from haplotagged BAMs
4. `phase_meth_to_founder_haps.sh` — calls `src/phase_meth_to_founder_haps.py` to phase methylation to founder haplotypes
5. `expand_to_all_cpgs.sh` — calls `src/expand_to_all_cpgs.py` to include all CpG sites and flag allele-specific sites

There are no automated tests in this repository.

## Architecture

### Core Python Scripts (`src/`)

- **`phase_meth_to_founder_haps.py`** — main pipeline orchestrator. Accepts 8+ CLI args (`--uid`, `--vcf_read_phased`, `--tsv_read_phase_blocks`, `--vcf_iht_phased`, `--txt_iht_blocks`, 4 methylation BED files, `--output_dir`). Outputs founder-phased BED and BigWig files.
- **`expand_to_all_cpgs.py`** — largest script (580 lines). Generalizes output to all CpG sites in reference and sample genome, flags allele-specific CpGs, generates BigWig files, and computes QC stats. Uses Jemalloc for memory optimization.
- **`get_all_phasing.py`** — combines read-backed (HiPhase) and inheritance-based phasing; uses cyvcf2 for VCF parsing; constructs bit vectors from het SNPs.
- **`get_hap_map.py`** — creates hap-map blocks (genomic intervals for founder phasing); computes concordance between phasing methods.
- **`get_meth_hap1_hap2.py`** — separates methylation by haplotype.

### Utility Modules (`src/util/`)

- **`read_data.py`** — `read_dataframe_from_bed()` reads gzipped BED files with `##key=value` metadata headers; `read_tapestry()` reads tapestry output with proper type casting.
- **`write_data.py`** — `write_dataframe_to_bed()` writes Polars DataFrames to BED with source metadata in `##source=...` header.
- **`shell.py`** — subprocess wrapper for bash execution.
- **`get_id_to_paths.py`**, **`get_palladium_prefixes.py`** — map sample IDs to file paths for the Palladium dataset.

### Data Conventions

- Primary data library is **Polars** (not Pandas).
- BED output files have a two-part header: `##key=value` metadata lines followed by a `#col1 col2 ...` column header line.
- The final expanded output (step 4) has 22 columns: chrom, start_cpg, end_cpg, total_read_count, methylation_level_count, methylation_level_model, start_hap_map_block, end_hap_map_block, haplotype_concordance_in_hap_map_block, num_het_SNVs_in_hap_map_block, total_read_count_pat, total_read_count_mat, founder_haplotype_pat, founder_haplotype_mat, methylation_level_pat_count, methylation_level_mat_count, methylation_level_pat_model, methylation_level_mat_model, cpg_is_within_50bp_of_mismatch_site, cpg_overlaps_at_least_one_snv, snv_genotypes, cpg_is_allele_specific.
- Jupyter notebooks (`.ipynb`) in `src/` serve as worked examples and analysis documentation.

## Trio-wise Workflow

A simpler alternative pipeline for a single parent–child trio. WhatsHap's `--ped` option phases hap1=paternal and hap2=maternal directly, so no inheritance-based phasing step is needed.

Five sequential shell scripts:

1. `run-whatshap-phase.sh` — unphase input VCF with `whatshap unphase`, then pedigree-aware phasing with `whatshap phase --ped` (unphasing first prevents mixed phasing: `whatshap phase --ped` copies over variants it cannot phase, so pre-existing phase assignments on those variants would otherwise be retained)
2. `run-whatshap-haplotag.sh` — tag child's reads with `HP:Z:1/2` using `whatshap haplotag`
3. `aligned_bam_to_cpg_scores.sh` — count- and model-based methylation from the haplotagged BAM (reused unchanged)
4. `trio_phase_meth_to_parent_haps.sh` — calls `src/trio_phase_meth_to_parent_haps.py`; maps hap1→pat and hap2→mat using phase blocks from `whatshap stats --block-list`; writes the same BED/BigWig format as `phase_meth_to_founder_haps.py` (`haplotype_concordance_in_hap_map_block` is null)
5. `expand_to_all_cpgs_trio.sh` — calls `src/expand_to_all_cpgs.py` without `--bed_het_site_mismatches` (which is now optional); `cpg_is_within_50bp_of_mismatch_site` is set to False for all rows

Mock test data for verification: `tests/mock_data/` (generated by `python tests/create_mock_data.py`).

### In-Progress Work

- `Snakefile` exists but is incomplete — intended to replace the manual shell-script workflow.
