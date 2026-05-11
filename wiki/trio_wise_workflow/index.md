# Trio-wise workflow

This section of the wiki walks through tapestry's trio-wise workflow,
in which pedMEC / whatshap phasing across a parent–parent–child trio
is used to label every CpG-level methylation measurement with one of
the parents' haplotypes — `A`/`B` in dad and `C`/`D` in mom.
`A`/`B`/`C`/`D` are *fixed as dad's hap1/hap2 and mom's hap1/hap2*;
they are not defined as transmitted vs non-transmitted.

For motivation — *why* phase methylation in the first place — see
[the shared motivation page](../motivation/motivation.md).

## Pages

| Page | What it covers |
|---|---|
| [pedMEC phasing](pedmec_phasing/pedmec_phasing.md) | Step 1 — `run-whatshap.sh`: trio-aware pedMEC phasing of the joint-called VCF, per-sample phase-block stats, and haplotagging of each sample's BAM. |
| [Parent-phased methylation](parent_phased_methylation/parent_phased_methylation.md) | Step 3 — `phase_meth_to_parent_haps.py` + `hap_map_trio.py`: the conceptual centre of the trio-wise workflow. Within each hap-map block (intersection of the kid's whatshap phase block and one parent's whatshap phase block), a bit-vector concordance decides which of the kid's two haplotypes descends from which parental homolog (`A`/`B` for dad, `C`/`D` for mom), and per-CpG methylation re-buckets mechanically. |
| [All-CpG expansion (trio)](all_cpg_expansion_trio/all_cpg_expansion_trio.md) | Step 4 — `expand_to_all_cpgs.trio.sh` + `src/expand_to_all_cpgs_trio.py`: expand the parent-phased methylation BED to every CpG in the reference, attach unphased per-sample methylation, flag CpGs within 50 bp of a bit-vector-mismatch heterozygous site, and label CpGs as allele-specific against the joint-called trio VCF. |

The trio-wise side is shorter than the pedigree-wise side because the
bit-vector concordance machinery is established in general form on the
[pedigree-wise founder-phased-methylation
page](../pedigree_wise_workflow/founder_phased_methylation/founder_phased_methylation.md);
the trio-wise version only has to explain the two differences:
(a) the phase source is pedMEC / whatshap rather than
`gtg-concordance`, and (b) the kid's haplotypes are labelled by
parental letters (A/B in dad, C/D in mom) rather than by
founder-of-the-pedigree letters.
