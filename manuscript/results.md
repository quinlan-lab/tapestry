Read-backed phasing of long reads splits methylation into two
tracks, for each each read-backed label (`hap1` / `hap2`; Fig 1B). But read-backed labels are limited: they carry no founder identity, flip independently
between adjacent read-backed phase blocks, and cannot be aligned
across generations. To 
support the parent-of-origin queries
that motivate this work, we must route methylation to founder haplotypes, 
viz., identify, at every position, which of each founder's two
homologs every descendant carries (Fig 2A).
Tapestry recovers this haplotype transmission purely from inheritance structure, 
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
in IBD segments, we do know each kid's two founder labels — one paternal, one maternal — from Fig 2. Fig 3A–D shows that information is enough to phase an uninformative site by enumerating all ways to assign parental alleles (0/1) to founder labels (A/B and C/D), followed by scoring each assignment against the kids' observed genotypes. Applying this phasing algorithm to every non-informative
site and composing with the informative-site deduction in Fig 2 yields,
for each kid, the full allele sequence corresponding to each founder
label that the kid carries (Fig 3E).

Having reconstructed haplotype sequences in IBD segments, 
we can finally route methylation to founder haplotypes in two steps, both operating
on the *hap-map block* (Fig 4A), defined as the intersection of a read-backed
phase block (Fig 1B) with an IBD segment (Fig 2D). Inside hap-map blocks, the
allele sequences are known for both the read-backed
labels (from the phased VCF) and the founder labels (from Fig 3).
**Step 1: map read-backed labels to founder labels (Fig 4A).** At
the heterozygous sites of the block, `hap1`'s allele sequence is
compared with that of the founder labels. The founder
whose sequence agrees with `hap1`'s more often is declared `hap1`'s
founder, and the larger concordance is retained as a per-block
quality score in tapestry's output. **Step 2: rebucket methylation from read-backed labels
to founder labels (Fig 4B).** Per-CpG methylation levels stratified
by read-backed labels are *rebucketed* in each hap-map block onto
founder-labelled tracks. Rebucketing can absorb meaningless block-to-block
flips in the raw read-backed tracks, exposing true biological
pattern across an IBD segment. The result — *founder-phased
methylation*, tapestry's primary product — is emitted as a per-CpG
BED (full column schema in repository's README) with per-founder read counts and
methylation levels,
hap-map block coordinates and allele-sequence concordance, founder
labels, and per-CpG QC flags. Tapestry also writes founder-phased
methylation bigwig tracks and haplotagged BAMs for IGV inspection.

With a BED of founder-phased methylation in hand, the two discovery scenarios
described in the introduction reduce to short queries against this
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
QTLs* (meQTLs) — loci where methylation tracks the local SNP
allele, so methylation is expected to differ between haplotypes
carrying different alleles at the QTL site (Fig 5A). If the founder
labels are correct, two predictions follow at such a locus.
First, methylation assigned to a given founder homolog should be
*similar across every descendant* that carries that homolog. Second, founder-grouped
methylation should split between the two alleles at the QTL site in
the direction reported for the meQTL. We tested both in the
four-generation CEPH1463 pedigree (K1463) at three published meQTLs (Fig 5B-D). 
At all three loci the methylation
levels that tapestry assigns to each of the eight founder homologs
across the pedigree group tightly within a homolog (prediction 1),
split cleanly between the two alleles at the meQTL (prediction 2),
and recover the published direction of the methylation–allele
association (e.g. C-allele homologs hypo-methylated and A-allele
homologs hyper-methylated at rs9330298). 
