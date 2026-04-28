"""Tapestry wiki figure generator.

Each page is produced by a `page_*` function that writes its narrative
`*.md` and figure assets into its own folder under `wiki/`. The function
for a given page is the single source of truth for that page; markdown
is written from inside the function so the prose and the figure
filenames stay in lockstep.

Run:
    python wiki/generate_wiki.py                     # all pages
    python wiki/generate_wiki.py --page motivation   # one page

Every wiki page that is *not* the vendored
`pedigree_wise_workflow/inheritance_mapping/` section is regenerated
by this script — including the catalog pages (`wiki/README.md`,
`pedigree_wise_workflow/index.md`, `trio_wise_workflow/index.md`).
The vendored section is regenerated separately; see that folder's
`README.md` for the sync procedure.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from _helpers import head_sha, permalink  # noqa: F401  (permalink reserved for later pages)

WIKI_DIR = Path(__file__).resolve().parent
SHA = head_sha()
sys.path.insert(0, str(WIKI_DIR))

from motivation import single_indiv_phasing, trio_discovery  # noqa: E402
from pedigree_wise_workflow.founder_phased_methylation import (  # noqa: E402
    founder_phased_methylation,
)
from trio_wise_workflow.parent_phased_methylation import (  # noqa: E402
    parent_phased_methylation,
)


# ---------------------------------------------------------------------------
# Page: motivation (shared between pedigree-wise and trio-wise workflows)
# ---------------------------------------------------------------------------

MOTIVATION_MD = """\
# Why phase DNA methylation?

This page is part of the [tapestry wiki](../README.md). It motivates
the rest of the wiki by answering one question: *why is it worth
phasing per-CpG methylation onto each haplotype, instead of leaving it
as a single number per CpG?* Two pairs of cartoon figures build the
answer on a small synthetic locus: the first pair takes a single
individual through a before/after of phasing to show what the read
partition does to the methylation profile; the second pair places that
individual into a trio so the two flagship use cases — *de novo
epimutations* and *compound genetic-epigenetic heterozygotes* — fall
out directly, together with the polars query a user runs against
tapestry's BED output to surface each one.

## Definitions

- **De novo epimutation** — a change in the methylation state of a
  specific physical homolog across a single meiosis. Detectable only
  when a parent and child can be compared on the *same* haplotype,
  which in turn requires phased methylation.
- **Compound genetic-epigenetic heterozygote** — an individual who
  inherits two different parental haplotypes such that one carries a
  genetic variant and the other carries an aberrant methylation state
  in *trans*. Invisible to either an unphased genotype call or an
  unphased methylation call, but detectable directly from per-haplotype
  tracks plus a trio.

## Single individual: before vs. after phasing

![Before phasing](single_indiv_before_phasing.png)

A small synthetic locus with mixed SNVs and CpGs. Ten ragged reads
cover overlapping windows; each read is one monospace row carrying a
`0`/`1` glyph at SNV columns and a small **red** (unmethylated) /
**blue** (methylated) box at CpG columns. Above the pile-up, a single
gray bigwig-style methylation track reports the *pooled* per-CpG
methylation level. With reads unpartitioned, the two haplotypes'
methylation states get averaged into one number per CpG — so per-CpG
heterozygous methylation is hidden inside an intermediate value.

![After phasing](single_indiv_after_phasing.png)

The same ten reads, now partitioned by haplotype: rows are grouped
paternal-above-maternal (one `Pat hap` / `Mat hap` label per group on
the left) and the SNV bits themselves carry the partition. With the
partition in hand, the pooled profile splits into two stacked
per-haplotype bigwig tracks. What Figure 2 shows is what tapestry computes mechanically for
every CpG once the haplotype partition is in hand. The two payoffs of
having those per-haplotype profiles only become visible once the
individual is placed into a pedigree.

## Trio use case 1: de novo epimutation

![Trio de novo](trio_denovo.svg)

The kid, now in a trio. Each individual carries two haplotypes drawn
as horizontal coloured bars; each bar is decorated with methylation
lollipops at CpG positions (filled = methylated, open = unmethylated).
At the right CpG cluster, two de novo events sit on opposite
parental haplotypes. The kid's paternal homolog (`A`, inherited from
dad) is methylated there, but dad's same physical homolog is
unmethylated — a de novo *gain* of methylation. Symmetrically, the
kid's maternal homolog (`C`, inherited from mom) is unmethylated at
the same cluster, but mom's same physical homolog is methylated — a
de novo *loss*. The left cluster is concordant on both sides
(unmethylated on `A` from dad to kid; methylated on `C` from mom to
kid).

![Trio de novo BED + polars](trio_denovo_bed.png)

Once tapestry produces a per-CpG haplotype-resolved BED, surfacing
the de novo gain of methylation on the kid's paternal hap is a
one-pass polars filter: pick the dad haplotype that was transmitted
to the kid (using the `pat_hap` phasing column), then keep CpGs where
`kid_pat - dad_transmitted` is *positive* on average over a small
window. Mirror the filter sign to surface a de novo *loss*; mirror
both the column choice (`mat_hap`/`mom_C`/`mom_D`) and the sign to
surface either de novo direction on the maternal side.

## Trio use case 2: compound genetic-epigenetic heterozygote

![Trio compound het](trio_compound_het.svg)

Same trio, different scenario: the kid inherits **dad's hap A** —
which carries a SNV variant (red star) — *and* **mom's hap C** —
which carries an aberrantly hyper-methylated promoter region (filled
lollipops at the leftmost three CpGs). The two hits sit on opposite
parental haplotypes, so the kid is a compound genetic-epigenetic
heterozygote *in trans* — silenced on the maternal allele by aberrant
methylation, broken on the paternal allele by a coding variant. Each
parent is a silent carrier of one of the two hits.

![Trio compound het BED + polars](trio_compound_het_bed.png)

Discovery from tapestry's BED is again a polars chain over two
inputs: a per-CpG methylation table and a per-SNV genotype table,
both with one column per parental haplotype (`dad_A`, `dad_B`,
`mom_C`, `mom_D`) plus the kid's per-side phased values
(`kid_pat`, `kid_mat`) and the phasing labels (`pat_hap`, `mat_hap`).
Step 1 finds methylation-outlier regions where exactly one parental
hap is the outlier and the kid inherits that hap with ~unchanged
meth. Step 2 finds genotype-outlier SNVs with the analogous
inheritance check. Step 3 joins the two and keeps loci where the
meth-outlier and geno-outlier come from *different* parents — i.e.,
the two hits are in trans.

## What's next

The rest of the wiki is concerned with how tapestry actually performs
the phasing summarised above:

- The [pedigree-wise workflow](../pedigree_wise_workflow/index.md)
  combines hiphase's read-backed phasing with `gtg-ped-map` /
  `gtg-concordance`'s inheritance-based phasing across a full
  pedigree, labelling each measurement with a founder haplotype.
- The [trio-wise workflow](../trio_wise_workflow/index.md) uses
  whatshap / pedMEC trio phasing to label each measurement with one
  of the parents' haplotypes — an option for the trio-only setting.
"""


def page_motivation(outdir: Path) -> None:
    page_dir = outdir / "motivation"
    page_dir.mkdir(parents=True, exist_ok=True)
    single_indiv_phasing.main()
    trio_discovery.main()
    md_path = page_dir / "motivation.md"
    md_path.write_text(MOTIVATION_MD)
    print(f"[wiki] Wrote {md_path}")


# ---------------------------------------------------------------------------
# Page: founder-phased methylation (Step 3)
# ---------------------------------------------------------------------------

FOUNDER_PHASED_METHYLATION_MD = f"""\
# Founder-phased methylation (Step 3)

This page is part of the
[pedigree-wise workflow](../index.md). It covers Step 3 of the
pipeline — `phase_meth_to_founder_haps.py` and `hap_map_pedigree.py` —
which is the conceptual centre of tapestry. Step 3 takes the two
phasings produced upstream and uses them to relabel each per-CpG
methylation measurement with one of the four founder homologs of
the pedigree.

## The two phasings, and where they meet

Two independent phasings reach Step 3:

- **Read-backed phasing** (Step 1A, hiphase). Within each *phase
  block*, hiphase emits two anonymous read partitions, `hap1` and
  `hap2`, and writes the partition into the BAM (`HP`/`PS` tags).
  Inside a block the two partitions are linked across heterozygous
  sites; across blocks the partition labels are independent.
- **Inheritance-based phasing** (Step 1B, `gtg-ped-map` +
  `gtg-concordance`). Within each *linkage block*,
  inheritance-based phasing reaches Step 3 as two parallel data
  structures:
  - The **iht-phased VCF**, which phases each het variant as `p|m`
    where `p` is the allele on the paternal homolog and `m` the
    allele on the maternal homolog. Tapestry parses this in
    [`get_iht_phasing`]({permalink('src/phasing_pedigree.py', 99, SHA)})
    (called from
    [`phase_meth_to_founder_haps.py:29`]({permalink('src/phase_meth_to_founder_haps.py', 29, SHA)})).
  - The **iht-blocks table** ([`df_iht_blocks`]({permalink('src/phase_meth_to_founder_haps.py', 32, SHA)})), which assigns one of the pedigree's
    *founder-haplotype letters* to the paternal and maternal slot in
    each block. The founders of the pedigree are the individuals at
    its root — for the CEPH1463 pedigree, the four great-grandparents
    — and `gtg-ped-map` labels their eight homologs with letters
    `A`, `B`, `C`, `D`, `E`, `F`, `G`, `H`. The kid's paternal and
    maternal homologs in each linkage block are therefore labelled
    with one of *those* root-founder letters, not with letters tied
    to dad and mom (unless the parents themselves happen to be
    founders, as in a strict nuclear-family setting). Within a block
    the assignment is fixed; across blocks it is not. Each linkage
    block is bounded by a crossover event in a transmitting
    ancestor — for example, a crossover on the paternal transmission
    path swaps the founder letter carried in the kid's paternal slot
    from one block to the next.

The two block partitions of the genome do not align. The natural
unit on which to relate them is the **hap-map block**: the
intersection of one read-backed phase block with one linkage block.
Inside such an intersection the read partition is consistent and
the founder labels are fixed, so a single decision relates one to
the other.

## Bit-vector match

Inside one hap-map block, each het SNV contributes one bit to three
parallel vectors:

- `hap1` — the read-partition allele on hiphase's `hap1` side.
- `pat` — the allele on the paternal founder homolog.
- `mat` — the allele on the maternal founder homolog.

![Bit-vector match](bit_vector_match.png)

`hap1` is compared against `pat` and `mat` in turn. Since every
contributing site is heterozygous, `pat` and `mat` carry opposite
alleles at each site; consequently the per-site match is mutually
exclusive, and `concordance(hap1, pat) + concordance(hap1, mat) = 1`
exactly. The decision is a single threshold:

- If `concordance(hap1, pat) > 0.5` → hiphase's `hap1` reads descend
  from the paternal founder homolog (and `hap2` reads from the
  maternal one).
- Otherwise → the assignment is flipped.

This decision is made once per hap-map block; the relevant code is
the bit-vector comparison in [`get_hap_map`]({permalink('src/hap_map_pedigree.py', 13, SHA)}),
specifically the `similarity_hap1_pat > 0.5` branch at
[line 58]({permalink('src/hap_map_pedigree.py', 58, SHA)}); the
function is called from
[`phase_meth_to_founder_haps.py:208`]({permalink('src/phase_meth_to_founder_haps.py', 208, SHA)}).

## Relabelling per-CpG methylation

Per-CpG methylation has already been stratified by hiphase's `hap1`
and `hap2` partition by the time Step 3 runs: Step 2
([`aligned_bam_to_cpg_scores.sh`]({permalink('aligned_bam_to_cpg_scores.sh', 1, SHA)}))
calls pb-CpG-tools once on each HP-tagged BAM, so each per-CpG
methylation level is computed only over reads carrying the same HP
tag. Step 3 does not recompute those per-haplotype levels — it
*relabels* them. Inside
[`phase_meth_to_founder_haps`]({permalink('src/phase_meth_to_founder_haps.py', 44, SHA)})
(called once per pb-CpG-tools mode at
[line 236]({permalink('src/phase_meth_to_founder_haps.py', 236, SHA)})
and
[line 238]({permalink('src/phase_meth_to_founder_haps.py', 238, SHA)})),
two stages of relabelling happen within each hap-map block:

1. **`hap1`/`hap2` → `pat`/`mat`.** The bit-vector decision recorded
   in `df_hap_map` is read off and used to copy the pre-computed
   `methylation_level_hap{1,2}` (and the matching read counts) onto
   `methylation_level_pat` and `methylation_level_mat` columns
   ([lines 71–103]({permalink('src/phase_meth_to_founder_haps.py', 71, SHA)})).
2. **`pat`/`mat` → founder letter.** The founder-haplotype letter
   on each parental slot, carried by `df_hap_map` from the
   iht-blocks table, is split off into `founder_haplotype_pat` and
   `founder_haplotype_mat`
   ([lines 104–105]({permalink('src/phase_meth_to_founder_haps.py', 104, SHA)})).

## Mismatch sites and Step 4 QC

Most hap-map blocks have a small number of sites where `hap1`
disagrees with the chosen parental vector — concordance < 1. Those
sites are recorded as `df_sites_mismatch` alongside the hap-map
block — see the mismatch-collection branches at
[lines 63–66]({permalink('src/hap_map_pedigree.py', 63, SHA)})
(when `hap1` matches `pat`) and
[lines 72–75]({permalink('src/hap_map_pedigree.py', 72, SHA)})
(when it matches `mat`), assembled into a single DataFrame at
[line 109]({permalink('src/hap_map_pedigree.py', 109, SHA)}) — and
propagated to Step 4
([all-CpG expansion](../all_cpg_expansion/all_cpg_expansion.md)),
where
[`compute_proximity_to_mismatched_heterozygous_sites`]({permalink('src/expand_to_all_cpgs.py', 115, SHA)})
annotates every CpG within 50 bp of one of them with the
[`cpg_is_within_50bp_of_mismatch_site`]({permalink('src/expand_to_all_cpgs.py', 378, SHA)})
QC flag in tapestry's BED output. Downstream consumers can choose
whether to keep or filter those CpGs.
"""


def page_founder_phased_methylation(outdir: Path) -> None:
    page_dir = outdir / "pedigree_wise_workflow" / "founder_phased_methylation"
    page_dir.mkdir(parents=True, exist_ok=True)
    founder_phased_methylation.main()
    md_path = page_dir / "founder_phased_methylation.md"
    md_path.write_text(FOUNDER_PHASED_METHYLATION_MD)
    print(f"[wiki] Wrote {md_path}")


# ---------------------------------------------------------------------------
# Catalog pages
# ---------------------------------------------------------------------------

WIKI_README_MD = """\
# tapestry — wiki

[![Why phase methylation?](motivation/single_indiv_before_phasing.png)](motivation/motivation.md)

*Bars in discordant region at 50 %: an unphased per-CpG methylation level cannot tell
two haplotypes apart. See the [motivation page](motivation/motivation.md)
for what the same locus looks like once phasing is done — and for the
two biological phenomena (de novo epimutations, compound genetic-epigenetic
heterozygotes) that phasing exposes.*

Tapestry is a pipeline that phases DNA methylation from
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
- **[Trio-wise workflow](trio_wise_workflow/index.md)** — pedMEC /
  whatshap phasing across a parent–parent–child trio, labelling each
  measurement with one of the parents' haplotypes. Most of this section
  is currently stubbed.

## Why phase methylation?

Both workflows are motivated by the same observation: an unphased
per-CpG methylation level averages over both haplotypes and erases the
two phenomena tapestry exists to surface — *de novo epimutations* (a
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
that tapestry's pedigree-wise workflow invokes.

The structure of both
wikis follows Andrej Karpathy's
[LLM-wiki pattern](https://gist.github.com/karpathy/442a6bf555914893e9891c11519de94f):
a thin, hand-curated catalog page links out to self-contained topic
pages, each bundled with its own assets and regenerable from source.
"""


PEDIGREE_WISE_INDEX_MD = """\
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
| [Read-backed phasing](read_backed_phasing/read_backed_phasing.md) | Step 1A — `run-hiphase.sh`: haplotagging reads with HP/PS tags, why two sibling phased VCFs (DeepVariant / pbsv) get produced. |
| [Inheritance mapping](inheritance_mapping/README.md) | Step 1B — `build-iht-based-haplotype-map-and-phase-variants.sh`: vendored from the upstream `Platinum-Pedigree-Inheritance` wiki; covers `gtg-ped-map` (nuclear family + three generations) and `gtg-concordance`. |
| [Founder-phased methylation](founder_phased_methylation/founder_phased_methylation.md) | Step 3 — `phase_meth_to_founder_haps.py` + `hap_map_pedigree.py`: the conceptual centre of tapestry. Within each hap-map block (intersection of a read-backed phase block and a linkage block), a bit-vector concordance decides which of hiphase's read partitions descends from which founder homolog, and per-CpG methylation re-buckets mechanically. |
| [All-CpG expansion](all_cpg_expansion/all_cpg_expansion.md) | Step 4 — `expand_to_all_cpgs.py`: reference CpGs vs sample CpGs vs measured CpGs, allele-specific CpGs, and the within-50bp-of-mismatch QC flag. |
"""


TRIO_WISE_INDEX_MD = """\
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
| [All-CpG expansion (trio)](all_cpg_expansion_trio/all_cpg_expansion_trio.md) | Step 4 — `expand_to_all_cpgs.trio.sh`: trio analogue of the pedigree-wise all-CpG expansion (reference CpGs vs sample CpGs vs measured CpGs, allele-specific CpGs, within-50bp-of-mismatch QC flag). *(page TODO; for now see [`expand_to_all_cpgs.trio.sh`](../../expand_to_all_cpgs.trio.sh).)* |

The trio-wise side is shorter than the pedigree-wise side because the
bit-vector concordance machinery is established in general form on the
[pedigree-wise founder-phased-methylation
page](../pedigree_wise_workflow/founder_phased_methylation/founder_phased_methylation.md);
the trio-wise version only has to explain the two differences:
(a) the phase source is pedMEC / whatshap rather than
`gtg-concordance`, and (b) the kid's haplotypes are labelled by
parental letters (A/B in dad, C/D in mom) rather than by
founder-of-the-pedigree letters.
"""


PEDMEC_PHASING_MD = f"""\
# pedMEC phasing (Step 1)

This page is part of the [trio-wise workflow](../index.md). It covers
Step 1 of the trio-wise pipeline:
[`run-whatshap.sh`]({permalink('run-whatshap.sh', 1, SHA)}) — the
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
[`whatshap unphase`]({permalink('run-whatshap.sh', 95, SHA)}) strips
all existing phase information. The unphased VCF is then bgzipped and
tabixed for use as input to the per-chromosome pedMEC step.

## 2. Pedigree-aware (pedMEC) phasing

Per-chromosome,
[`whatshap phase --ped`]({permalink('run-whatshap.sh', 128, SHA)})
runs the pedMEC algorithm jointly on the three samples (kid, dad,
mom) using the trio's PED file and all three HiFi BAMs as evidence.
pedMEC adds an inheritance-consistency constraint to single-sample
read-based MEC, which resolves phasing ambiguities that single-sample
whatshap cannot. By tapestry's convention the kid is ordered as
`pat|mat`, so the kid's `hap1` is paternal and `hap2` is maternal;
dad is `A|B` and mom is `C|D`.

The 25 chromosomes (chr1–22, X, Y, M) are phased in parallel
([`xargs -P PHASE_THREADS`]({permalink('run-whatshap.sh', 142, SHA)}))
and the per-chromosome phased VCFs are concatenated with
[`bcftools concat --naive`]({permalink('run-whatshap.sh', 154, SHA)})
into one genome-wide phased VCF.

## 3. Per-sample block stats

For each of the three samples,
[`whatshap stats --block-list`]({permalink('run-whatshap.sh', 166, SHA)})
emits a TSV listing every phase block (chromosome, span, phase-set
ID, number of variants). This is the file the next step of the
pipeline reads via
[`get_phase_blocks`]({permalink('src/phasing_trio.py', 7, SHA)})
(consumed at
[`phase_meth_to_parent_haps.py:359`]({permalink('src/phase_meth_to_parent_haps.py', 359, SHA)})).

## 4. Haplotag the BAMs

For each sample, prior `HP`/`PS` tags are stripped from the input
BAM and
[`whatshap haplotag`]({permalink('run-whatshap.sh', 200, SHA)})
re-tags the reads from the pedMEC-phased VCF, producing a
`*.haplotagged.bam` whose reads carry trio-aware `HP`/`PS` tags.
This haplotagged BAM is the input to Step 2
([`aligned_bam_to_cpg_scores.sh`]({permalink('aligned_bam_to_cpg_scores.sh', 1, SHA)})
run in *trio mode* via `-t kid_id dad_id mom_id`), which calls
pb-CpG-tools once per HP-tagged BAM and
produces the per-haplotype methylation BEDs that Step 3 relabels
onto the parental-letter alphabet.
"""


PARENT_PHASED_METHYLATION_MD = f"""\
# Parent-phased methylation (Step 3)

This page is part of the
[trio-wise workflow](../index.md). It covers Step 3 of the trio-wise
pipeline — `phase_meth_to_parent_haps.py` and `hap_map_trio.py` —
which is the conceptual centre of tapestry's trio-wise side. Step 3
takes the single phasing produced upstream by pedMEC / whatshap and
uses it to relabel each per-CpG methylation measurement with one of
the four parental homologs of the trio (`A`/`B` in dad; `C`/`D` in
mom).

This page mirrors the structure of the pedigree-wise
[founder-phased-methylation page](../../pedigree_wise_workflow/founder_phased_methylation/founder_phased_methylation.md).
Read that page first for the general bit-vector machinery; this page
only spells out the two differences specific to the trio-wise case.

## The phasing, and where blocks meet

Unlike the pedigree-wise workflow — which combines a read-backed
phasing (hiphase) with a separate inheritance-based phasing
(`gtg-ped-map` + `gtg-concordance`) — the trio-wise workflow has a
single phasing source: **pedMEC / whatshap**, run once on the
trio's joint-called VCF
([`run-whatshap.sh:128`]({permalink('run-whatshap.sh', 128, SHA)})).
That single run produces, for each individual, a set of *whatshap
phase blocks*, where the alleles inside one block are linked across
heterozygous sites and across blocks the labels are independent.
Tapestry parses those blocks per individual via
[`get_phase_blocks`]({permalink('src/phasing_trio.py', 7, SHA)})
(called from
[`phase_meth_to_parent_haps.py:359–361`]({permalink('src/phase_meth_to_parent_haps.py', 359, SHA)}))
and parses paired kid-parent allele sequences from the pedMEC-phased
VCF in
[`get_all_phasing`]({permalink('src/phasing_trio.py', 48, SHA)})
(called from
[`phase_meth_to_parent_haps.py:366`]({permalink('src/phase_meth_to_parent_haps.py', 366, SHA)})).

The natural unit on which to relate the kid's phasing to a parent's
phasing is the **hap-map block**: the intersection of one of the
kid's whatshap phase blocks with one of the parent's whatshap phase
blocks. There are *two* such intersections per kid block — one with
dad and one with mom — handled independently and labelled
*paternal* and *maternal* downstream.

## Bit-vector match

Inside one paternal-side hap-map block, each het SNV in dad
contributes one bit to three parallel vectors:

- `kid_pat` — the kid's allele on the paternal homolog (kid's
  whatshap `hap1`, by tapestry's convention that pedMEC orders the
  kid as `pat|mat`).
- `dad_A` — dad's allele on hap1 (= `A` by definition).
- `dad_B` — dad's allele on hap2 (= `B` by definition).

![Bit-vector match (trio)](bit_vector_match_trio.png)

`kid_pat` is compared against `dad_A` and `dad_B` in turn. Since
every contributing site is heterozygous in dad, `dad_A` and `dad_B`
carry opposite alleles at each site; consequently the per-site match
is mutually exclusive, and
`concordance(kid_pat, dad_A) + concordance(kid_pat, dad_B) = 1`
exactly. The decision is a single threshold:

- If `concordance(kid_pat, dad_A) > 0.5` → kid's paternal homolog
  in this block is `A` (and would have been `B` otherwise).

This decision is made once per paternal-side hap-map block; the
relevant code is the bit-vector comparison in
[`_build_hap_map`]({permalink('src/hap_map_trio.py', 10, SHA)}),
specifically the `similarity > 0.5` branch at
[line 71]({permalink('src/hap_map_trio.py', 71, SHA)}); the
`(A, B)` label pair is wired in at
[`get_hap_map`]({permalink('src/hap_map_trio.py', 113, SHA)})
(called from
[`phase_meth_to_parent_haps.py:375`]({permalink('src/phase_meth_to_parent_haps.py', 375, SHA)})).
The maternal side is symmetric: replace `kid_pat`/`dad_A`/`dad_B`
with `kid_mat`/`mom_C`/`mom_D` and intersect the kid's phase block
with mom's instead of dad's.

## Relabelling per-CpG methylation

Per-CpG methylation has already been stratified by each individual's
whatshap `hap1`/`hap2` partition by the time Step 3 runs: Step 2
([`aligned_bam_to_cpg_scores.sh`]({permalink('aligned_bam_to_cpg_scores.sh', 1, SHA)}))
calls pb-CpG-tools once on each HP-tagged BAM, so each per-CpG
methylation level is computed only over reads carrying the same HP
tag. Step 3 does not recompute those per-haplotype levels — it
*relabels* them. Inside
[`phase_meth_to_parent_haplotypes`]({permalink('src/phase_meth_to_parent_haps.py', 45, SHA)})
(called once at
[line 407]({permalink('src/phase_meth_to_parent_haps.py', 407, SHA)})),
the kid's `hap1`/`hap2` and each parent's `hap1`/`hap2` are renamed
to the parent-letter alphabet: kid `hap1`/`hap2` → `kid_pat`/`kid_mat`,
dad `hap1`/`hap2` → `dad_A`/`dad_B`, mom `hap1`/`hap2` →
`mom_C`/`mom_D` ([lines 80–110]({permalink('src/phase_meth_to_parent_haps.py', 80, SHA)})).
The paternal and maternal hap maps are then overlaid onto each CpG
([lines 117–155]({permalink('src/phase_meth_to_parent_haps.py', 117, SHA)})),
and CpGs that fall outside any hap-map block on a given parental
side are nulled out on that side
([lines 157–260]({permalink('src/phase_meth_to_parent_haps.py', 157, SHA)})).
Note that, unlike the pedigree-wise workflow, there is no second
relabelling stage from `pat`/`mat` to founder letters: the trio-wise
output stops at the parental-letter alphabet `A`/`B`/`C`/`D` because
the kid's parents *are* the most distal individuals in this workflow.

## Mismatch sites

Most hap-map blocks have a small number of sites where `kid_pat`
disagrees with the chosen `dad_A`/`dad_B` vector — concordance < 1.
Those sites are recorded as `df_mismatch_pat` (and `df_mismatch_mat`
for the maternal side) alongside the hap-map block — see the
mismatch-collection branch at
[line 86]({permalink('src/hap_map_trio.py', 86, SHA)}) of
`hap_map_trio.py`, populated whichever way the
[`similarity > 0.5`]({permalink('src/hap_map_trio.py', 71, SHA)})
branch went. Both mismatch DataFrames are returned by
[`get_hap_map`]({permalink('src/hap_map_trio.py', 113, SHA)}) and
written out to BED + VCF by
[`phase_meth_to_parent_haps.py:385–388`]({permalink('src/phase_meth_to_parent_haps.py', 385, SHA)})
for the same downstream proximity-flag QC that the pedigree-wise
workflow performs in Step 4.
"""


def page_wiki_readme(outdir: Path) -> None:
    md_path = outdir / "README.md"
    md_path.write_text(WIKI_README_MD)
    print(f"[wiki] Wrote {md_path}")


def page_pedigree_wise_index(outdir: Path) -> None:
    page_dir = outdir / "pedigree_wise_workflow"
    page_dir.mkdir(parents=True, exist_ok=True)
    md_path = page_dir / "index.md"
    md_path.write_text(PEDIGREE_WISE_INDEX_MD)
    print(f"[wiki] Wrote {md_path}")


def page_trio_wise_index(outdir: Path) -> None:
    page_dir = outdir / "trio_wise_workflow"
    page_dir.mkdir(parents=True, exist_ok=True)
    md_path = page_dir / "index.md"
    md_path.write_text(TRIO_WISE_INDEX_MD)
    print(f"[wiki] Wrote {md_path}")


def page_pedmec_phasing(outdir: Path) -> None:
    page_dir = outdir / "trio_wise_workflow" / "pedmec_phasing"
    page_dir.mkdir(parents=True, exist_ok=True)
    md_path = page_dir / "pedmec_phasing.md"
    md_path.write_text(PEDMEC_PHASING_MD)
    print(f"[wiki] Wrote {md_path}")


def page_parent_phased_methylation(outdir: Path) -> None:
    page_dir = outdir / "trio_wise_workflow" / "parent_phased_methylation"
    page_dir.mkdir(parents=True, exist_ok=True)
    parent_phased_methylation.main()
    md_path = page_dir / "parent_phased_methylation.md"
    md_path.write_text(PARENT_PHASED_METHYLATION_MD)
    print(f"[wiki] Wrote {md_path}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

PAGES = {
    "wiki_readme": page_wiki_readme,
    "motivation": page_motivation,
    "pedigree_wise_index": page_pedigree_wise_index,
    "founder_phased_methylation": page_founder_phased_methylation,
    "trio_wise_index": page_trio_wise_index,
    "pedmec_phasing": page_pedmec_phasing,
    "parent_phased_methylation": page_parent_phased_methylation,
}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--page", choices=sorted(PAGES.keys()),
        help="Render only one page (default: render every page).",
    )
    parser.add_argument(
        "--outdir", type=Path, default=WIKI_DIR,
        help="Directory to write the wiki into.",
    )
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    targets = [PAGES[args.page]] if args.page else list(PAGES.values())
    for fn in targets:
        fn(args.outdir)


if __name__ == "__main__":
    main()
