# Plan: building the tapestry wiki

Living planning document. Records design choices, page inventory, and
phased implementation order. Detailed figure specs that have been
implemented are deliberately *not* in this file — the code in the
relevant page folder is the spec. Specs are kept here only for figures
not yet built.

## 1. Goals

1. Stand up a wiki for the **pedigree-wise** tapestry workflow
   (`run-hiphase.sh` →
   `build-iht-based-haplotype-map-and-phase-variants.sh` →
   `aligned_bam_to_cpg_scores.sh` → `phase_meth_to_founder_haps.sh` →
   `expand_to_all_cpgs.sh`).
2. Re-use — not re-write — the three walkthrough pages from the
   Platinum-Pedigree-Inheritance wiki
   (`nuclear_family/`, `three_generations/`, `concordance/`), since
   `build-iht-based-haplotype-map-and-phase-variants.sh` is the tapestry
   step that runs `gtg-ped-map` + `gtg-concordance`.
3. Leave a clean seam to add a **trio-wise** workflow section later
   (pedMEC/whatshap, `phase_meth_to_parent_haps.py`,
   `expand_to_all_cpgs_trio.py`).
4. Every figure is **regenerable from `wiki/generate_wiki.py`**. The
   wiki is a reference users consult while drafting the Bioinformatics
   manuscript in Illustrator; PNGs are not piped into the manuscript
   verbatim.

## 2. Design principles (inherited from the upstream wiki)

The upstream `Platinum-Pedigree-Inheritance/wiki` follows the
"LLM-wiki" pattern (Karpathy
[tweet](https://x.com/karpathy/status/2040572272944324650),
[gist](https://gist.github.com/karpathy/442a6bf555914893e9891c11519de94f)):
a thin catalog `index.md` links to self-contained topic pages, each
bundled with its assets. Concretely:

- **One directory per page.** `<page>.md` + `fig*.png` co-located.
- **Descriptive folder names** (`read_backed_phasing/`, not `step1/`).
- **Commit-SHA-pinned source permalinks** via `_helpers.py`.
- **Reproducible figures.** Hard-coded toy simulations; no RNG seeds,
  no real data. `python wiki/generate_wiki.py [--page X]`.
- **Monospace matplotlib panels** for ASCII-style walkthroughs.
- **Hand-written SVG (Python f-strings + stdlib) for pixel-perfect,
  non-data-driven cartoons.** Where a figure is a fixed schematic
  (pedigree shapes, haplotype bars, lollipops, callouts) rather than
  a render of toy data, prefer emitting SVG directly from Python
  rather than going through matplotlib. Vector output, exact control
  over geometry, no Agg/dpi dance. See `motivation/trio_discovery.py`
  for the pattern (`rect`/`circle`/`line`/`star` primitives plus a
  `build()` that composes a `Scenario` dataclass into one SVG
  document).
- **No cross-page implicit state.** Every page is independently
  readable.

## 3. Vendoring the upstream Rust wiki

Three options were considered; chosen: **vendor** the three upstream
pages into `wiki/pedigree_wise_workflow/inheritance_mapping/`, marked
"synced from upstream at commit `<sha>`", and resync via the procedure
documented in that folder's `README.md`. Self-contained site, single
upstream source of truth, easy resync.

Concretely, on each sync:
- Pull upstream + run its `generate_wiki.py`.
- `rsync` the three subfolders into
  `pedigree_wise_workflow/inheritance_mapping/`.
- Apply one local edit: rewrite `[wiki](../index.md)` →
  `[wiki](../../../index.md)` so "up" lands on the tapestry catalog,
  not the upstream one.
- Upstream permalinks (Rust source) are left alone.

Until the user's open PR against the upstream repo merges, "upstream"
is the user's fork (`petermchale/Platinum-Pedigree-Inheritance`).

## 4. Directory layout

```
wiki/
  index.md
  claude_plan.md            # this file
  generate_wiki.py
  _helpers.py

  motivation/               # shared between pedigree-wise + trio-wise
    motivation.md
    single_indiv_phasing.py
    trio_discovery.py
    single_indiv_before_phasing.png
    single_indiv_after_phasing.png
    trio_denovo.svg
    trio_denovo_bed.png
    trio_compound_het.svg
    trio_compound_het_bed.png

  pedigree_wise_workflow/
    index.md
    overview/                  # CEPH1463 pedigree + DAG
    read_backed_phasing/       # Step 1A — hiphase
    inheritance_mapping/       # Step 1B — VENDORED (see §3)
      README.md
      nuclear_family/
      three_generations/
      concordance/
    founder_phased_methylation/  # Step 3 — heart of the workflow (§8)
    all_cpg_expansion/         # Step 4
    output_format/             # BED schema + column dictionary

  trio_wise_workflow/          # mostly stubbed
    index.md
    output_format/        # built early — column dictionary
```

Step 2 (`aligned_bam_to_cpg_scores.sh`) is a pure pass-through of the
upstream pb-CpG-tools binary; both workflows describe it in one
paragraph inside their `overview.md`.

## 5. Page-by-page manifest

### 5.1 `wiki/index.md`
One-paragraph elevator pitch + two workflow links + reproduction
command. Hero image: `motivation/single_indiv_before_phasing.png`.

### 5.2 `motivation/motivation.md` — IMPLEMENTED
Shared hook page. Six figures: a single-individual before/after pair
showing the read-partition payoff, then a trio scenario per use case
(de novo epimutation; compound genetic-epigenetic heterozygote), each
with the SVG cartoon plus a matplotlib panel showing the polars
discovery query against tapestry's BED. Generated by
`single_indiv_phasing.py` and `trio_discovery.py` in the same folder
(both expose a `main()` invoked from `generate_wiki.py`).

### 5.3 `pedigree_wise_workflow/index.md`
Catalog scoped to this workflow. One row per page below.

### 5.4 `pedigree_wise_workflow/overview/overview.md`
CEPH1463 pedigree (re-drawn from `images/pedigree.jpg`) + workflow DAG
(re-rendered from `docs/pedigree_workflow.mmd` via `mmdc`) + a
paragraph per step linking forward.

### 5.5 `pedigree_wise_workflow/read_backed_phasing/`
Step 1A (hiphase). What hiphase does, two toy figures (phase block
formation, BAM HP/PS tags), pointers to upstream docs.

### 5.6 `pedigree_wise_workflow/inheritance_mapping/` — VENDORED (§3)
Step 1B. Three sub-folders, vendored. A short top-level
`inheritance_mapping/README.md` notes the upstream SHA, the resync
procedure, and the canonical upstream URL.

### 5.7 `pedigree_wise_workflow/founder_phased_methylation/`
Step 3 — `phase_meth_to_founder_haps.py`. Conceptual heart of the
workflow. Hero figure: bit-vector matching cartoon (§8). Sections:
two phasing sources → bit-vector matching across block intersection →
founder-haplotype labeling of haplotagged reads → mechanical
re-bucketing of methylation into pat/mat. Plus an IGV screenshot reused
from `images/tapestry.pedigree.png`. QC: mismatch sites annotated in
the cartoon, become the "within 50 bp of mismatch site" flag in
Step 4's output.

### 5.8 `pedigree_wise_workflow/all_cpg_expansion/`
Step 4 — `expand_to_all_cpgs.py`. Three figures: ref vs. sample vs.
measured CpG universes (Venn); allele-specific CpG toggle (SNV
creates/destroys a CpG); mismatch-site proximity QC.

### 5.9 `pedigree_wise_workflow/output_format/`
Lifts the column dictionary from `README.md` lines 121–156, expanded
with worked examples and back-pointers to the page that explains each
column. No new figures.

### 5.10 `trio_wise_workflow/output_format/` — built early
Trio-output BED column dictionary lifted from `README.md` lines
189–239. The discovery worked example is **already covered** by the
matplotlib panels in `motivation/trio_discovery.py`
(`trio_denovo_bed.png` + `trio_compound_het_bed.png`); this page
links back to those rather than re-emitting them.

## 6. `generate_wiki.py` — script architecture

- Top-level CLI: `--page <name>`; with no arg, regenerate everything
  *except* the vendored `inheritance_mapping/` pages.
- One `page_*()` function per page; each writes both `*.md` and
  `fig*.png` atomically into its folder.
- Per-page rendering code lives **inside the page's folder** as a
  module (e.g. `motivation/single_indiv_phasing.py`,
  `motivation/trio_discovery.py`), exposing a `main()` that
  `generate_wiki.py` calls. Keeps each page a self-contained unit and
  keeps `generate_wiki.py` as a thin orchestrator.
- Stdlib + matplotlib only. Hard-coded deterministic toys.
- Captions use the `_helpers.permalink(path, line, sha)` helper.

## 7. Trio-wise workflow — future extension

The sibling directory `wiki/trio_wise_workflow/` exists as a stub
except for `output_format/`, which is built early. Pending:
`pedmec_phasing/` (whatshap/pedMEC), `parent_haplotype_phasing/`
(`phase_meth_to_parent_haps.py`), `all_cpg_expansion_trio/`
(`expand_to_all_cpgs.trio.sh`).

The pedigree-wise `founder_phased_methylation/` page already explains
bit-vector concordance in the general setting, so
`parent_haplotype_phasing/` will only need to highlight the
differences: (a) phase source is pedMEC/whatshap rather than
`gtg-concordance`; (b) kid haplotypes are labelled by *parental*
letters (A/B in dad, C/D in mom) — fixed as dad's hap1/hap2 and mom's
hap1/hap2, **not** transmitted/non-transmitted.

## 8. Bit-vector matching cartoon (founder_phased_methylation hero figure)

Single end-to-end cartoon at
`pedigree_wise_workflow/founder_phased_methylation/fig1_bit_vector_matching.png`,
illustrating `phase_meth_to_founder_haps.py` + `hap_map_pedigree.py`.
Four stacked panels sharing a horizontal genomic axis:

- **Panel A.** Two phasing sources side by side. Left: hiphase output —
  ~6 haplotagged reads (4 hap1 above, 2 hap2 below) at ~5 het-SNV
  columns; collapsed bit vectors `hap1: 0 1 0 1 1`, `hap2: 1 0 1 0 0`.
  Right: gtg-concordance output — same five sites with phased
  genotypes `p|m`; collapsed `pat: 0 1 0 1 0`, `mat: 1 0 1 0 1`.
- **Panel B.** Block intersection shaded. Below, align `hap1` vs `pat`
  character-by-character with ✓/✗ → 4/5 = 0.8; complementary `hap1`
  vs `mat` → 1/5 = 0.2. The single mismatch site is circled and
  labeled "flagged for QC — `cpg_is_within_50bp_of_mismatch_site`".
- **Panel C.** Panel A's hiphase pileup re-drawn with founder labels
  per read (`dad_hap1` for the four `hap1` reads, `mom_hap2` for the
  two `hap2` reads), inferred from Panel B.
- **Panel D.** A CpG column added; each read carries `●`/`○`. Because
  reads are bucketed into `hap1.bed.gz` / `hap2.bed.gz` by HP tag and
  Panel C has now attached a founder label to each bucket, methylation
  rebuckets mechanically: 4/4 methylated on pat → 100 %; 0/2 on mat →
  0 %. Render on a 0–100 % bar pair.

Toy data hard-coded; no RNG. `motivation/`'s 0–100 % bar style is
reused for Panel D so it visually anchors back to the motivation page.

## 9. Open questions

1. **Upstream for `inheritance_mapping/`.** Currently
   `petermchale/Platinum-Pedigree-Inheritance` (the user's fork);
   switch to `Platinum-Pedigree-Consortium/...` once the user's PR
   merges?
2. **Permalink commit pinning.** Pin to `HEAD` at generate-time, or to
   a named release tag? Current behaviour: `HEAD`.
3. **Re-drawing `images/pedigree.jpg`.** The existing JPG is 3.7 MB
   and looks photographed. For the `pedigree_wise_workflow/overview/` page it can start as a
   direct embed; the manuscript figure will need a clean vector
   re-draw — base on a published pedigree style or draw from scratch?
4. **GitHub Pages vs. raw markdown.** Upstream uses raw markdown
   browsing; default to that unless the user asks otherwise.

---

*Owner: plan authored by Claude; implementation interactive against
this plan, component at a time.*
