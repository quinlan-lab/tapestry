# pedMEC phasing (Step 1)

This page is part of the [trio-wise workflow](../index.md). It covers
Step 1 of the trio-wise pipeline:
[`run-whatshap.sh`](https://github.com/quinlan-lab/tapestry/blob/97c13602dbbaea12a3aad781310c1b9a3a5359a4/run-whatshap.sh#L1) — the
script that takes a joint-called trio VCF and a set of HiFi BAMs and
produces, for each sample, a pedMEC-phased VCF, a per-sample
phase-block table, and a haplotagged BAM. The downstream
[parent-phased-methylation page](../parent_phased_methylation/parent_phased_methylation.md)
consumes the per-sample phase-block tables and the multi-sample
phased VCF;
[`aligned_bam_to_cpg_scores.sh`](../../../aligned_bam_to_cpg_scores.sh)
in *trio mode* (`-t kid_id dad_id mom_id`) consumes the haplotagged
BAMs.

The script proceeds in four steps.

## 1. Unphase the input VCF

The input joint-called VCF arrives with whatever phasing the upstream
caller emitted. To avoid mixed phasing,
[`whatshap unphase`](https://github.com/quinlan-lab/tapestry/blob/97c13602dbbaea12a3aad781310c1b9a3a5359a4/run-whatshap.sh#L95) strips
all existing phase information. The unphased VCF is then bgzipped and
tabixed for use as input to the per-chromosome pedMEC step.

## 2. Pedigree-aware (pedMEC) phasing

Per-chromosome,
[`whatshap phase --ped`](https://github.com/quinlan-lab/tapestry/blob/97c13602dbbaea12a3aad781310c1b9a3a5359a4/run-whatshap.sh#L128)
runs the pedMEC algorithm jointly on the three samples (kid, dad,
mom) using the trio's PED file and all three HiFi BAMs as evidence.
pedMEC adds an inheritance-consistency constraint to single-sample
read-based MEC, which resolves phasing ambiguities that single-sample
whatshap cannot. By tapestry's convention the kid is ordered as
`pat|mat`, so the kid's `hap1` is paternal and `hap2` is maternal;
dad is `A|B` and mom is `C|D`.

The 25 chromosomes (chr1–22, X, Y, M) are phased in parallel
([`xargs -P PHASE_THREADS`](https://github.com/quinlan-lab/tapestry/blob/97c13602dbbaea12a3aad781310c1b9a3a5359a4/run-whatshap.sh#L142))
and the per-chromosome phased VCFs are concatenated with
[`bcftools concat --naive`](https://github.com/quinlan-lab/tapestry/blob/97c13602dbbaea12a3aad781310c1b9a3a5359a4/run-whatshap.sh#L154)
into one genome-wide phased VCF.

## 3. Per-sample block stats

For each of the three samples,
[`whatshap stats --block-list`](https://github.com/quinlan-lab/tapestry/blob/97c13602dbbaea12a3aad781310c1b9a3a5359a4/run-whatshap.sh#L166)
emits a TSV listing every phase block (chromosome, span, phase-set
ID, number of variants). This is the file the next step of the
pipeline reads via
[`get_phase_blocks`](https://github.com/quinlan-lab/tapestry/blob/97c13602dbbaea12a3aad781310c1b9a3a5359a4/src/phasing_trio.py#L7)
(consumed at
[`phase_meth_to_parent_haps.py:359`](https://github.com/quinlan-lab/tapestry/blob/97c13602dbbaea12a3aad781310c1b9a3a5359a4/src/phase_meth_to_parent_haps.py#L359)).

## 4. Haplotag the BAMs

For each sample, prior `HP`/`PS` tags are stripped from the input
BAM and
[`whatshap haplotag`](https://github.com/quinlan-lab/tapestry/blob/97c13602dbbaea12a3aad781310c1b9a3a5359a4/run-whatshap.sh#L200)
re-tags the reads from the pedMEC-phased VCF, producing a
`*.haplotagged.bam` whose reads carry trio-aware `HP`/`PS` tags.
This haplotagged BAM is the input to Step 2
([`aligned_bam_to_cpg_scores.sh`](https://github.com/quinlan-lab/tapestry/blob/97c13602dbbaea12a3aad781310c1b9a3a5359a4/aligned_bam_to_cpg_scores.sh#L1)
run in *trio mode* via `-t kid_id dad_id mom_id`), which calls
pb-CpG-tools once per HP-tagged BAM and
produces the per-haplotype methylation BEDs that Step 3 relabels
onto the parental-letter alphabet.
