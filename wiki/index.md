# tapestry — wiki

A reference for tapestry, a pipeline that phases DNA methylation from
PacBio HiFi reads onto the haplotypes of a human pedigree's founders
(or, in the trio-wise special case, onto the haplotypes of each
parent). Tapestry combines read-backed phasing of single-sample
variants with inheritance-based phasing across the pedigree, then uses
the intersection of the two block structures to label every CpG-level
methylation measurement with its founder (or parental) haplotype of
origin.

This wiki is a worked-out, line-numbered exposition of the method and
the pipeline. It is meant to be read alongside the source — every
figure caption and prose reference points at a specific commit-pinned
line in the tapestry repo (or, for the inheritance-mapping section, in
the upstream `Platinum-Pedigree-Inheritance` repo). The wiki is
deliberately *not* the source of figures for the Bioinformatics
manuscript: paper figures are composed separately in Illustrator and
optimise for a different audience. See
[`claude_plan.md`](claude_plan.md) for the design rationale and the
phased build plan.

## Two workflows

Tapestry exposes two phasing strategies that converge on the same kind
of output (haplotype-resolved per-CpG methylation):

- **[Pedigree-wise workflow](pedigree_wise_workflow/index.md)** —
  read-backed phasing (hiphase) plus inheritance-based phasing
  (`gtg-ped-map` + `gtg-concordance`) across a full pedigree. Labels
  every methylation measurement with one of the pedigree's founder
  haplotypes. This is the workflow built out first; it is the focus of
  the current wiki.
- **[Trio-wise workflow](trio_wise_workflow/README.md)** — pedMEC /
  whatshap phasing across a parent–parent–child trio, labelling each
  measurement with one of the parents' haplotypes. Most of this section
  is currently stubbed; the
  [output-format page](trio_wise_workflow/output_format_trio/output_format_trio.md)
  is built out early because it hosts the searchable-output worked
  example used to motivate downstream analyses.

## Why phase methylation?

Both workflows are motivated by the same observation: an unphased
per-CpG methylation level averages over both haplotypes and erases the
two phenomena tapestry exists to surface — *epimutations* (a
methylation change on a specific physical homolog across a single
meiosis) and *compound genetic-epigenetic heterozygotes* (a haplotype
where a SNV genotype and a methylation state co-segregate). The
shared [motivation page](motivation/motivation.md) walks through both
with a small worked toy.

## How to reproduce the figures in this wiki

Every figure in the tapestry wiki (everything *except* the vendored
`pedigree_wise_workflow/inheritance_mapping/` section) is regenerated
deterministically by

```
python wiki/generate_wiki.py
```

The script takes an optional `--page <name>` to regenerate a single
page's figures and markdown. All toy simulations are hard-coded; no
RNG, no real VCFs, no external data. Re-running is byte-reproducible.

The vendored `inheritance_mapping/` section is regenerated from its
own upstream repo (see
[`pedigree_wise_workflow/inheritance_mapping/README.md`](pedigree_wise_workflow/inheritance_mapping/README.md)
for the sync procedure and the pinned upstream SHA).

## Acknowledgements

The `pedigree_wise_workflow/inheritance_mapping/` section is vendored
from the
[Platinum-Pedigree-Inheritance wiki](https://github.com/petermchale/Platinum-Pedigree-Inheritance/tree/main/wiki),
which documents the `gtg-ped-map` and `gtg-concordance` Rust binaries
that tapestry's pedigree-wise workflow invokes. The structure of both
wikis follows Andrej Karpathy's
[LLM-wiki pattern](https://gist.github.com/karpathy/442a6bf555914893e9891c11519de94f):
a thin, hand-curated catalog page links out to self-contained topic
pages, each bundled with its own assets and regenerable from source.
