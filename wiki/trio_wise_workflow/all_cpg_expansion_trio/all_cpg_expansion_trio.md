# All-CpG expansion (Step 4, trio)

This page is part of the [trio-wise workflow](../index.md). It covers
Step 4 of the trio-wise pipeline:
[`expand_to_all_cpgs.trio.sh`](https://github.com/quinlan-lab/tapestry/blob/af81a6bfe3fe2d66df9212d5df1e55fd5129d395/expand_to_all_cpgs.trio.sh#L1)
and its Python driver
[`src/expand_to_all_cpgs_trio.py`](https://github.com/quinlan-lab/tapestry/blob/af81a6bfe3fe2d66df9212d5df1e55fd5129d395/src/expand_to_all_cpgs_trio.py#L1).

By the end of Step 3, parent-phased methylation is reported only at
the subset of CpGs where the kid's haplotagged reads supplied a
phaseable measurement. Step 4 expands that table to every CpG in the
reference, attaches per-sample unphased methylation as a fallback,
adds two QC labels (proximity to bit-vector mismatches; allele-specific
flags from the joint-called VCF), and writes a sorted, tabix-indexed
BED.

The script proceeds in four stages.

## 1. Enumerate every CpG in the reference

[`write_all_cpgs_in_reference.py`](https://github.com/quinlan-lab/tapestry/blob/af81a6bfe3fe2d66df9212d5df1e55fd5129d395/src/write_all_cpgs_in_reference.py#L1)
([`scan_and_write_cpgs`](https://github.com/quinlan-lab/tapestry/blob/af81a6bfe3fe2d66df9212d5df1e55fd5129d395/src/write_all_cpgs_in_reference.py#L8))
scans the GRCh38 FASTA and emits a BED of every `CG` dinucleotide. It
is invoked from
[`expand_to_all_cpgs.trio.sh`](https://github.com/quinlan-lab/tapestry/blob/af81a6bfe3fe2d66df9212d5df1e55fd5129d395/expand_to_all_cpgs.trio.sh#L71)
and produces `all_cpg_sites_in_reference.bed`, the spine that the
remaining stages join against.

## 2. Expand parent-phased methylation to all CpGs and attach unphased meth

[`expand_to_all_cpgs_trio.py`](https://github.com/quinlan-lab/tapestry/blob/af81a6bfe3fe2d66df9212d5df1e55fd5129d395/src/expand_to_all_cpgs_trio.py#L630)
reads three inputs: the all-CpG BED from stage 1, the parent-phased
methylation BED from Step 3, and the per-sample (kid/dad/mom)
count-based and model-based unphased methylation BEDs produced by
`aligned_bam_to_cpg_scores.sh` in trio mode.
[`read_meth_unphased_trio`](https://github.com/quinlan-lab/tapestry/blob/af81a6bfe3fe2d66df9212d5df1e55fd5129d395/src/expand_to_all_cpgs_trio.py#L83)
joins the count and model BEDs per individual; then
[`expand_meth_to_all_cpgs`](https://github.com/quinlan-lab/tapestry/blob/af81a6bfe3fe2d66df9212d5df1e55fd5129d395/src/expand_to_all_cpgs_trio.py#L120)
left-joins the parent-phased and unphased frames onto the all-CpG
spine, so every reference CpG appears as a row, with parent-phased
levels (`A`/`B`/`C`/`D`) where available and unphased per-sample
levels alongside. CpGs with no measurement carry nulls.

In parallel,
[`write_combined_bigwig`](https://github.com/quinlan-lab/tapestry/blob/af81a6bfe3fe2d66df9212d5df1e55fd5129d395/src/expand_to_all_cpgs_trio.py#L361)
emits one bigWig per (sample × `count`/`model`) combination of the
unphased input BEDs, for IGV display.

## 3. Within-50bp-of-mismatch QC flag

For each side of the trio (paternal, maternal), Step 3 wrote a BED of
heterozygous sites at which the bit-vector concordance disagreed
between the kid and that parent (see the
[parent-phased-methylation page](../parent_phased_methylation/parent_phased_methylation.md)).
[`compute_proximity_to_mismatched_heterozygous_sites`](https://github.com/quinlan-lab/tapestry/blob/af81a6bfe3fe2d66df9212d5df1e55fd5129d395/src/expand_to_all_cpgs_trio.py#L202)
labels each CpG with its distance to the nearest mismatch site on each
side and sets a within-50bp flag. The fraction of CpGs flagged is
logged by
[`compute_fraction_of_cpgs_that_are_close_to_mismatches`](https://github.com/quinlan-lab/tapestry/blob/af81a6bfe3fe2d66df9212d5df1e55fd5129d395/src/expand_to_all_cpgs_trio.py#L233).
This is the trio analogue of the same QC step in the pedigree-wise
all-CpG expansion: phase-decision quality near locally-discordant
sites is suspect, and downstream analyses can drop or down-weight
these CpGs.

## 4. Allele-specific labelling against the joint-called VCF

[`get_joint_called_variants`](https://github.com/quinlan-lab/tapestry/blob/af81a6bfe3fe2d66df9212d5df1e55fd5129d395/src/expand_to_all_cpgs_trio.py#L421)
reads the trio's joint-called multi-sample VCF and pulls SNV
genotypes for kid, dad, and mom.
[`label_with_variants`](https://github.com/quinlan-lab/tapestry/blob/af81a6bfe3fe2d66df9212d5df1e55fd5129d395/src/expand_to_all_cpgs_trio.py#L484)
annotates each CpG with whether it overlaps 0/1/2 SNVs, and
[`label_cpgs_as_allele_specific`](https://github.com/quinlan-lab/tapestry/blob/af81a6bfe3fe2d66df9212d5df1e55fd5129d395/src/expand_to_all_cpgs_trio.py#L557)
flags CpGs as allele-specific (per family member) when they overlap a
heterozygous SNV — that is, the methylation signal at those CpGs may
reflect cis-acting allele-specific effects rather than parental
inheritance per se.

The expanded, fully-labelled frame is then written as a sorted,
bgzipped, tabixed BED at
`bed_meth_parent_phased_all_cpgs` — the canonical per-CpG output of
the trio-wise workflow.
