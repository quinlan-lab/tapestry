# Why phase DNA methylation?

This page is part of the [tapestry wiki](../index.md). It motivates
the rest of the wiki by answering one question: *why is it worth
phasing per-CpG methylation onto each haplotype, instead of leaving it
as a single number per CpG?* Two cartoon figures take a small synthetic
locus, show what it looks like before phasing, then redraw the same
reads after phasing to make two biological phenomena visible.

The motivation does not depend on which phasing strategy gets you to
parental-haplotype resolution, so this page is shared between
tapestry's [pedigree-wise](../pedigree_wise_workflow/index.md) and
[trio-wise](../trio_wise_workflow/README.md) workflows.

## Definitions

- **Epimutation** — a change in the methylation state of a specific
  physical homolog across a single meiosis. Detectable only when a
  parent and child can be compared on the *same* haplotype, which in
  turn requires phased methylation.
- **Compound genetic-epigenetic heterozygote** — an individual who is
  simultaneously heterozygous at a SNV *and* heterozygous in
  methylation state at a nearby CpG, with the methylation state
  tracking the SNV genotype across the two haplotypes. Invisible to
  either an unphased genotype call or an unphased methylation call,
  but trivial to read off paired phased tracks.

## Figure 1 — before phasing

![Figure 1 — before](fig1_before_unphased.png)

A small synthetic locus with two SNVs (`SNV1`, `SNV2`) bracketing two
CpGs (`CpG1`, `CpG2`). Twelve reads span the whole locus. Each read
is annotated at the four marked positions: a `0` or `1` at each SNV
(REF/ALT bit, matching the convention used throughout the wiki and in
the vendored
[`inheritance_mapping/`](../pedigree_wise_workflow/inheritance_mapping/README.md)
section), and a filled (`●`) or open (`○`) circle at each CpG
indicating methylation status on that read.

In this toy, six of the twelve reads carry `(0, 0)` at the two SNVs
and the other six carry `(1, 1)` — informative phasing anchors that
the unphased pile-up is *not yet using*. The two-bar plot to the
right reports the pooled methylation level at each CpG: both come
out at 50 %, the least-informative possible value. A 50 % per-CpG
number is consistent with many biological truths — both haplotypes
at 50 %, one fully methylated and the other unmethylated, asymmetric
methylation between maternal and paternal homologs — and the
unphased pile-up cannot distinguish them.

## Figure 2 — after phasing

![Figure 2 — after](fig2_after_phased.png)

The same twelve reads, now redrawn under the assumption that phasing
has been done.

- **Panel A** regroups the reads by haplotype (paternal above
  maternal, with a faint divider) and recolours them. Both SNVs
  agree on the assignment, so phasing is robust against a single
  sequencing error at either site.
- **Panel B** replaces each of Figure 1's two pooled bars with two
  per-haplotype bars, drawn on the same 0–100 % y-axis. CpG1 comes
  out at 0 % paternal / 100 % maternal; CpG2 at 100 % paternal /
  0 % maternal. Each Figure 1 bar was hiding a maximally-far-apart
  pair of haplotype-level values, and the two CpGs point in
  *opposite* directions — an unphased aggregate erases both.
- **Panel C** zooms out across ~30 CpGs in a short genomic window.
  Most positions are concordant between the two haplotypes (both
  methylated, both unmethylated). Two positions diverge: CpG1, at
  which the haplotype-resolved methylation co-segregates with the
  SNV1 genotype (a *compound genetic-epigenetic heterozygote*), and
  CpG2, at which the child's paternal homolog is methylated but
  dad's same physical homolog is not (an *epimutation* — a
  methylation change on a specific physical homolog across one
  meiosis).

## What's next

The rest of the wiki is concerned with how tapestry actually performs
the phasing summarised by Figure 2:

- The [pedigree-wise workflow](../pedigree_wise_workflow/index.md)
  combines hiphase's read-backed phasing with `gtg-ped-map` /
  `gtg-concordance`'s inheritance-based phasing across a full
  pedigree, labelling each measurement with a founder haplotype.
- The [trio-wise workflow](../trio_wise_workflow/README.md) uses
  whatshap / pedMEC trio phasing to label each measurement with one
  of the parents' haplotypes — an option for the trio-only setting.

For a worked example of how to query tapestry's output for the two
phenomena above (epimutations, compound genetic-epigenetic
heterozygotes), see
[`trio_wise_workflow/output_format_trio/output_format_trio.md`](../trio_wise_workflow/output_format_trio/output_format_trio.md).
