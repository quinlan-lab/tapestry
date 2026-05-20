Read-backed phasing of long reads splits methylation into two
tracks, for each each read-backed label (`hap1` / `hap2`; Fig 1B). But read-backed labels are limited: they carry no founder identity, flip independently
between adjacent read-backed phase blocks, and cannot be aligned
across generations. To 
support the parent-of-origin queries
that motivate this work, we must route methylation to founder haplotypes, 
viz., identify, at every position, which of each founder's two
homologs every descendant carries (Fig 2A).
Tapestry recovers this haplotype transmission purely from inheritance structure, 
using unphased joint-called variants (not reads) — an instance of *genetic phasing*, 
the assignment of alleles to haplotypes from the Mendelian transmission pattern in a 
family ([Roach et al. 2011, Am J Hum Genet](https://doi.org/10.1016/j.ajhg.2011.07.023)). 
At *mom-informative* sites
(mom heterozygous, dad homozygous) the maternally
transmitted allele is deduced from each kid's genotype, and, strung together along the chromosome,
these alleles reveal each kid's
maternal haplotype (Fig 2B). The same procedure on the
paternal side produces, in this particular case, putative crossovers co-located
across independent meioses (Fig 2C iii). Since that is vanishingly unlikely at
the per-meiosis crossover rate of ≪1/Mb, these are "switch errors". As is standard in pedigree 
haplotyping, we favour the reconstruction implying the fewest recombinants 
([Abecasis et al. 2002, Nat Genet](https://doi.org/10.1038/ng786)), which we enforce  
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

Though the alleles at *informative* sites are deduced in the process of
assigning founder labels, such a deduction is not possible at a
*non-informative* site as both parents are heterozygous. However, 
in IBD segments, we do know each kid's two founder labels — one paternal, one maternal — from Fig 2. Fig 3A–D shows that information is enough to phase an uninformative site by enumerating all ways to assign parental alleles (0/1) to founder labels (A/B and C/D), followed by scoring each assignment against the kids' observed genotypes (prior art: [Eberle et al. 2017, Genome Res](https://doi.org/10.1101/gr.210500.116); [Kronenberg et al. 2025](https://www.nature.com/articles/s41592-025-02750-y)). Applying this phasing algorithm to every non-informative
site and composing with the informative-site deduction in Fig 2 yields,
for each kid, the full allele sequence corresponding to each founder
label that the kid carries (Fig 3E).

The steps so far — genetic phasing (Fig 2) and the reconstruction of 
founder allele sequences (Fig 3) — adapt established genetic-phasing 
methodology. We now turn to the unsolved problem of 
routing methylation onto those founder haplotypes. 

Phasing methylation proceeds in two steps, both operating
on the *hap-map block* (Fig 4A), defined as the intersection of a read-backed
phase block (Fig 1B) with an IBD segment (Fig 2D). Inside hap-map blocks, the
allele sequences are known for both the read-backed
labels (from the phased VCF) and the founder labels (from Fig 3).
Step 1 is to map read-backed labels to founder labels. At
the heterozygous sites of the block, `hap1`'s allele sequence is
compared with that of the founder labels (Fig 4A). The founder
whose sequence agrees with `hap1`'s more often is declared `hap1`'s
founder, and the larger concordance is retained as a per-block
quality score in tapestry's output. In step 2, per-CpG methylation levels stratified
by read-backed labels are *rebucketed* in each hap-map block onto
founder-labelled tracks (Fig 4B). Rebucketing can absorb meaningless block-to-block
flips in the raw read-backed tracks, exposing true biological
pattern across an IBD segment (Fig 4B). The result — *founder-phased
methylation*, tapestry's primary product — is emitted as a per-CpG
BED (full column schema in repository's README) with per-founder read counts and
methylation levels,
hap-map block coordinates and allele-sequence concordance, founder
labels, and per-CpG QC flags. Tapestry also writes founder-phased
methylation bigwig tracks and haplotagged BAMs for IGV inspection.
In the four-generation CEPH K1463 pedigree, Tapestry phases ~95% of CpG sites to founder haplotypes, with <2% of CpG sites lying within 50bp of a discordant heterozygous site (c.f., Fig 4A), consistent with the finding that genetic phasing approaches completeness once three or more children are available ([Roach et al. 2011, Am J Hum Genet](https://doi.org/10.1016/j.ajhg.2011.07.023)). 

With a BED of founder-phased methylation in hand, the two discovery scenarios
described in the introduction become short queries against this
BED (Supp Fig 1) — and, because the BED reports every CpG covered by
the HiFi reads, both queries run genome-wide rather than at a
pre-selected list of loci. A homolog-specific LOM
at an imprinted locus surfaces as a contiguous run of CpGs at which
a kid's methylation on one transmitted founder haplotype falls sharply below
the level seen on the *same* founder haploptype in the transmitting parent
(Fig 4C). Restricting the genome-wide scan to kids who share the disease
phenotype, and to a parent of origin consistent with the locus's
imprint, should recover the segregation pattern expected for an imprinting
disorder. A compound genetic-epigenetic heterozygote surfaces as the
inner join of two outlier queries — a methylation outlier on one
parental founder and a coding loss-of-function variant outlier on
the other — restricted to loci where the founder labels place the
two hits in trans (Fig 4D).

To validate tapestry's founder labels we turned to *methylation
QTLs* (meQTLs). If the founder
labels are correct, we expect methylation at such loci to differ only between haplotypes
carrying different alleles at the QTL site, and to otherwise be independent of haplotype, all regardless of who carries the haplotypes (Fig 5A). 
We tested these expectations in the
four-generation CEPH K1463 pedigree at three published meQTLs (Fig 5B-D). 
At all three loci, the methylation
levels that tapestry assigns to each of the eight founder homologs
across the pedigree (i) group tightly within a homolog,
(ii) split cleanly between two classes of homologs depending on the meQTL allele they harbor,
and (iii) recover the published direction of the methylation–allele
association (e.g. C-allele homologs hypo-methylated and A-allele
homologs hyper-methylated at rs9330298). These results support the correctness of tapestry's founder-phasing of methylation levels. 
