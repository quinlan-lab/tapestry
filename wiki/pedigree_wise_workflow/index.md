# Pedigree-wise workflow

This section of the wiki walks through tapestry's pedigree-wise
workflow, in which read-backed phasing (hiphase) is combined with
inheritance-based phasing across a full pedigree (`gtg-ped-map` +
`gtg-concordance`) to label every CpG-level methylation measurement
with a founder haplotype.

For motivation — *why* phase methylation in the first place — see
[the shared motivation page](../motivation/motivation.md).

## Pages

| Page | What it covers |
|---|---|
| [Overview](overview/overview.md) | CEPH1463 pedigree diagram, end-to-end pipeline DAG, one paragraph per step. |
| [Read-backed phasing](read_backed_phasing/read_backed_phasing.md) | Step 1A — `run-hiphase.sh`: haplotagging reads with HP/PS tags, why two sibling phased VCFs (DeepVariant / pbsv) get produced. |
| [Inheritance mapping](inheritance_mapping/README.md) | Step 1B — `build-iht-based-haplotype-map-and-phase-variants.sh`: vendored from the upstream `Platinum-Pedigree-Inheritance` wiki; covers `gtg-ped-map` (nuclear family + three generations) and `gtg-concordance`. |
| [Founder-phased methylation](founder_phased_methylation/founder_phased_methylation.md) | Step 3 — `phase_meth_to_founder_haps.py`: the conceptual centre of tapestry. Bit-vector matching across the intersection of a hiphase phase block and an iht block, founder-haplotype labelling of haplotagged reads, and mechanical relabelling of per-CpG methylation. |
| [All-CpG expansion](all_cpg_expansion/all_cpg_expansion.md) | Step 4 — `expand_to_all_cpgs.py`: reference CpGs vs sample CpGs vs measured CpGs, allele-specific CpGs, and the within-50bp-of-mismatch QC flag. |
| [Output format](output_format/output_format.md) | BED schema and column dictionary for the pedigree-wise output. |

