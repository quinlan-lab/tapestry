Read-backed phasing of long reads splits methylation into two
tracks, for each each read-backed label (`hap1` / `hap2`; Fig 1B). But read-backed labels are limited: they carry no founder identity, flip independently
between adjacent read-backed phase blocks, and cannot be aligned
across generations. To 
support the parent-of-origin queries
that motivate this work, we must route methylation to founder haplotypes, 
viz., identify, at every position, which of each founder's two
homologs every descendant carries (Fig 2A).
Tapestry recovers this purely from inheritance structure, 
using unphased joint-called variants (not reads). 
At *mom-informative* sites
(mom heterozygous, dad homozygous) the maternally
transmitted allele is deduced from each kid's genotype, and, strung together along the chromosome,
these alleles reveal each kid's
maternal haplotype (Fig 2B). The same procedure on the
paternal side produces, in this particular case, putative crossovers co-located
across independent meioses (Fig 2C iii) — vanishingly unlikely at
the per-meiosis crossover rate of ≪1/Mb. We resolve this 
by swapping the conflicting labels (Fig 2C iv), yielding 
the far more likely event of a single crossover (Fig 2C v).
Intersecting the paternal and maternal haplotype segments gives *IBD segments*
(Fig 2D).
The same machinery generalises to a multi-generation
pedigree, labelling each descendant's homologs by the *root* founders'
homologs. Fig 2 thus attaches *founder labels* (A-D in a trio; A–H in a
four-grandparent pedigree) to each descendant's two homologs within
an IBD segment. To map read-backed labels (`hap1`/`hap2`; Fig 1B) onto these founder labels, we will shortly match their *allele sequences*. The allele sequences corresponding to each read-backed label come from the phased
VCF emitted by the read-backed phasing software, but the allele sequences corresponding to each founder label
remain to be
reconstructed. We now turn to that task.

Though the alleles at *informative* sites were deduced in the process of
assigning founder labels, such a deduction is not possible at a
*non-informative* site as both parents are heterozygous. However, 
in IBD segments, we do know each kid's two founder labels — one paternal, one maternal. Fig 3A–D shows that information is enough to phase an uninformative site by enumerating all ways to assign parental alleles (0/1) to founder labels (A/B and C/D), followed by scoring each assignment against the kids' observed genotypes. Applying this phasing algorithm to every non-informative
site and composing with the informative-site deduction in Fig 2 yields,
for each kid, the full allele sequence corresponding to each founder
label that the kid carries (Fig 3E).

The conceptual centre of the workflow is the *hap-map block*
(Fig 4A): the intersection of a read-backed phase block (which
defines `hap1`/`hap2` contiguously, but without founder identity)
with an inheritance-based linkage block (which defines the founder
labels, but says nothing about which reads belong to which homolog).
Inside a hap-map block, both pieces of information are simultaneously
defined, and a single bit-vector concordance — the per-site agreement
between `hap1`'s allele sequence and each founder's reconstructed
bitstring at heterozygous sites in the block — decides which of
`hap1` and `hap2` descends from which founder. Per-CpG methylation
values are then *relabelled* (not recomputed) by routing each
`hap1`/`hap2` track into its founder-labelled track within each
hap-map block (Fig 4B). The relabelling is what removes the
apparent block-to-block flips in the raw per-haplotype tracks: the
`hap1`/`hap2` ↔ founder mapping can swap between adjacent read-backed
phase blocks inside the same IBD segment, and a consistent
biological pattern (here founder `A` hyper-methylated, founder `C`
hypo-methylated) only emerges after the relabelling. The
concordance value is retained alongside each block as a per-block
quality score.

If tapestry's founder labels are correct, two predictions follow.
First, methylation assigned to a given founder homolog should be
*similar across every descendant* that carries that homolog —
identity-by-descent of the homolog implies identity of its
methylation state up to measurement noise. Second, at a *methylation
QTL* (meQTL), founder-grouped methylation should track the founder's
allele at the QTL site (Fig 5A). We tested both predictions in the
four-generation CEPH1463 pedigree (K1463) at three meQTLs reported
by Rosenski et al. 2025
([doi:10.1038/s41467-025-57433-1](https://doi.org/10.1038/s41467-025-57433-1),
[doi:10.1101/2025.09.15.675351](https://doi.org/10.1101/2025.09.15.675351)).
At all three loci — rs9330298 (Fig 5B), rs12499263 (Fig 5C), and
rs12636296 at the *ALG1L* promoter (Fig 5D) — the methylation
levels that tapestry assigns to each of the eight founder homologs
across the pedigree group tightly within a homolog, split cleanly
between the two alleles at the meQTL, and recover the published
direction of the methylation–allele association (e.g. C-allele
homologs hypo-methylated and A-allele homologs hyper-methylated at
rs9330298). The clean within-homolog clustering across multiple
descendants is the stronger of the two checks: it would not appear
under any mislabelling of homolog identity.

The two scenarios that motivate the workflow reduce to short queries
against the per-CpG haplotype-resolved BED that tapestry emits
(Supp Fig 1), and — because the BED reports every CpG covered by
the HiFi reads — both queries run genome-wide rather than at a
pre-selected list of loci. A homolog-specific LOM at an imprinted
locus surfaces as a contiguous run of CpGs at which a kid's
methylation on one transmitted founder departs sharply from the
level seen on the *same* founder in the transmitting parent
(Fig 4C); restricting the run to kids who share the disease
phenotype, and to a parent of origin consistent with the
locus's imprint, recovers the segregation pattern expected for
an imprinting disorder. A compound genetic-epigenetic
heterozygote surfaces as the inner join of two outlier queries — a
methylation outlier on one parental founder and a deleterious-variant
outlier on the other — restricted to loci where the founder labels
place the two hits in trans (Fig 4D).

The primary output of each workflow is a per-CpG BED file with a
self-describing header (run metadata, full column schema; see
README for the exact columns). For each CpG site it reports total
and per-founder read counts, count-based and model-based methylation
levels on each founder homolog
([pb-CpG-tools](https://github.com/PacificBiosciences/pb-CpG-tools),
[doi placeholder](https://doi.org/PLACEHOLDER)), the enclosing
hap-map block coordinates with its bit-vector concordance and
heterozygous-SNV count, and the founder labels themselves. The
workflow also emits a small set of per-CpG QC flags — proximity to
sites where read-backed and inheritance-based bit vectors disagree
(a local phasing-disagreement marker), overlap with SNVs, and
whether the CpG is *allele-specific* (present in the reads of only
one homolog, e.g. created or destroyed by a heterozygous variant) —
and summary statistics over the run (e.g. the fraction of CpGs
successfully phased to at least one founder, and the fraction within
50 bp of a bit-vector-mismatch site). Alongside the BED, the
workflow writes a small set of files for downstream inspection in
IGV: per-founder methylation bigwig tracks, BED files of hap-map
block boundaries and of bit-vector-mismatch sites, and the
haplotagged BAMs they overlay.
