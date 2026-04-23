# Plan: building the tapestry wiki

This is a planning document. No wiki pages or scripts have been written yet;
this file records design choices, an inventory of pages to build, and a phased
implementation order so the work can be picked up later without re-discovering
context.

## 1. Goals

1. Stand up a wiki for the **pedigree-wise** tapestry workflow (the four-step
   pipeline already documented in `README.md`: `run-hiphase.sh` →
   `build-iht-based-haplotype-map-and-phase-variants.sh` →
   `aligned_bam_to_cpg_scores.sh` → `phase_meth_to_founder_haps.sh` →
   `expand_to_all_cpgs.sh`).
2. Re-use — not re-write — the three existing walkthrough pages from the
   Platinum-Pedigree-Inheritance wiki
   (`nuclear_family/`, `three_generations/`, `concordance/`), since
   `build-iht-based-haplotype-map-and-phase-variants.sh` is the tapestry
   step that actually runs `gtg-ped-map` and `gtg-concordance`. Those pages
   should appear inside the tapestry wiki as a dedicated section so a reader
   of the tapestry manuscript does not have to leave the site to understand
   inheritance-based phasing.
3. Leave a clean seam to later add a **trio-wise** workflow section
   (pedMEC/whatshap phasing, `phase_meth_to_parent_haps.py`,
   `expand_to_all_cpgs_trio.py`).
4. Keep every figure in the wiki **regenerable from a single `generate_wiki.py`
   script**. The wiki's role vis-à-vis the Bioinformatics manuscript is as a
   *reference the user will consult while drafting the paper's figures* — a
   worked-out, reproducible, line-numbered exposition of the method and the
   pipeline. The paper's figures will be composed separately (in Adobe
   Illustrator); wiki PNGs are not piped into the manuscript verbatim. The
   wiki therefore does not need a manuscript-assembly step or a
   figure-numbering convention that mirrors the paper — it only needs to be
   clear, complete, and easy to re-read.

## 2. Design principles (inherited from the upstream wiki)

The upstream `Platinum-Pedigree-Inheritance/wiki` follows Andrej Karpathy's
"LLM-wiki" pattern (see his
[tweet](https://x.com/karpathy/status/2040572272944324650) and accompanying
[gist](https://gist.github.com/karpathy/442a6bf555914893e9891c11519de94f)):
a thin, hand-curated catalog page (`index.md`) links out to self-contained
topic pages, each bundled with its own assets. The tapestry wiki will use
the same pattern, with these concrete rules:

- **One directory per page.** Each page has its own folder containing
  `<page>.md` + `fig*.png`. No shared image pool; co-locating assets keeps
  pages moveable.
- **Descriptive names, not numeric.** Folder names describe content
  (`read_backed_phasing/`, not `step1/`).
- **Commit-SHA-pinned source permalinks.** Every figure caption and prose
  reference to source code uses
  `https://github.com/quinlan-lab/tapestry/blob/<sha>/...#L<line>` with
  `<sha>` captured from `git rev-parse HEAD` at generate time. For the
  inheritance-mapping section the permalinks point to the
  `Platinum-Pedigree-Consortium/Platinum-Pedigree-Inheritance` repo and its
  pinned commit, exactly as the upstream pages do today.
- **Reproducible figures.** Every figure in the wiki is emitted by
  `wiki/generate_wiki.py`, driven by hard-coded toy simulations (no VCFs,
  no real data, no RNG seeds). `python wiki/generate_wiki.py --page X`
  regenerates one page's assets; running it with no args regenerates
  everything.
- **Monospace matplotlib panels.** Continue the upstream convention of
  text-in-matplotlib so every panel reads as ASCII but exports as
  vector-friendly PNG.
- **No cross-page implicit state.** A reader landing on any single page gets
  enough context (or a one-paragraph "see also" link) to read it
  independently.

## 3. How the rust wiki will be merged in

Three options were considered:

| Option | Pros | Cons |
|---|---|---|
| A. Link out to the upstream wiki from a short "Step 1B" stub. | Zero duplication. | Network-dependent; breaks if upstream moves; tapestry readers have to click away to follow the method. |
| B. Copy the three upstream pages into tapestry wholesale, then maintain both copies. | Self-contained site. | Two sources of truth; guaranteed to drift. |
| C. **Hybrid (chosen):** vendor the upstream pages into `wiki/pedigree_wise_workflow/inheritance_mapping/`, marked as "synced from upstream at commit `<sha>`," and sync by running the upstream `generate_wiki.py` (which is idempotent) and checking the result into tapestry. | Self-contained site; single upstream source of truth; easy to re-sync. | Requires a documented sync procedure. |

Option C is the plan. Concretely:

- Inside the pedigree-wise section (see §4), create
  `wiki/pedigree_wise_workflow/inheritance_mapping/` as a sibling of the
  other pedigree-wise page folders (because it is specifically the
  Step 1B documentation for the pedigree-wise workflow; the trio-wise
  workflow has its own Step 1, which uses whatshap/pedMEC instead).
- Inside it, vendor three sub-folders mirroring upstream exactly:
  `nuclear_family/`, `three_generations/`, `concordance/`. Each sub-folder
  contains the full `*.md` + `fig*.png` bundle already produced by the
  upstream `generate_wiki.py`.
- Put a short `inheritance_mapping/README.md` at the top of that folder
  stating: (a) this section is vendored from
  `petermchale/Platinum-Pedigree-Inheritance/wiki/` at commit `<sha>`;
  (b) re-sync instructions ("`git -C ~/Platinum-Pedigree-Inheritance pull
  && python wiki/generate_wiki.py && rsync -a ~/Platinum-Pedigree-Inheritance/wiki/{nuclear_family,three_generations,concordance} wiki/pedigree_wise_workflow/inheritance_mapping/`",
  or equivalent); (c) link to the upstream URL for the canonical copy.
- Update the relative links inside the vendored pages only where they
  cross the section boundary (the `[wiki](../index.md)` links at the top
  of each upstream page should now resolve to `../../../index.md` — the
  top-level tapestry wiki index — so a reader navigating "up" lands on
  the tapestry catalog, not the upstream one). All other internal links
  (`../three_generations/three_generations.md` from `nuclear_family.md`,
  etc.) continue to work because the three sub-folders retain their
  upstream sibling structure.
- Leave the upstream GitHub-permalink line numbers alone; they should
  continue to point at the upstream Rust commit SHA, since that is where
  `gtg-ped-map` / `gtg-concordance` actually live.

Note: the user has an open PR against the upstream repo that is not yet
merged. Until it lands, the "upstream" for sync purposes is the user's
fork (`petermchale/Platinum-Pedigree-Inheritance`). The sync-command
example in `inheritance_mapping/README.md` should reflect that, and be
updated once the PR is merged upstream.

## 4. Proposed directory layout

The two workflows live in **sibling directories** under `wiki/`. Neither
workflow owns the top-level wiki; `wiki/index.md` is a thin catalog that
routes a reader to one or the other. The pedigree-wise side is the
focus for now; the trio-wise side is a stub (see §8).

```
wiki/
  index.md                         # top-level catalog; routes to the two workflow sections
  claude_plan.md                   # this file
  generate_wiki.py                 # one script, --page <path> dispatches to one component

  motivation/                      # shared motivation page — see §5.2 and §11
    motivation.md                  # why phase methylation at all?
    fig1_before_unphased.png       # toy reads with CpG + SNV; one averaged meth level
    fig2_after_phased.png          # same reads partitioned by haplotype; two profiles
                                   # revealing epimutations and compound het examples

  pedigree_wise_workflow/
    index.md                       # catalog for this workflow, sibling to trio's catalog
    overview/
      overview.md                  # pedigree diagram, pipeline DAG
      fig1_pedigree.png            # CEPH1463 pedigree diagram
      fig2_workflow_dag.png        # re-rendered from docs/pedigree_workflow.mmd
    read_backed_phasing/
      read_backed_phasing.md       # Step 1A — hiphase + haplotagging
      fig1_reads_to_phase.png      # toy: 3 reads, 2 hets, phase-block formation
      fig2_haplotag.png            # HP/PS tags on a read alignment
    inheritance_mapping/           # Step 1B — VENDORED from upstream wiki (see §3)
      README.md                    # sync procedure, upstream SHA
      nuclear_family/
        nuclear_family.md
        fig1.png ... fig6_3.png
      three_generations/
        three_generations.md
        fig1.png fig2.png
      concordance/
        concordance.md
        fig1.png fig2.png fig3.png
    founder_phased_methylation/
      founder_phased_methylation.md  # Step 3 — phase_meth_to_founder_haps.py
      fig1_bit_vector_matching.png   # hero cartoon — see §13 — 4 stacked panels:
                                     # (A) haplotagged reads + two bit vectors from the
                                     #     two phasing sources,
                                     # (B) bit-vector matching in the block intersection,
                                     # (C) reads relabeled with founder haplotype,
                                     # (D) methylation re-bucketed into pat/mat bars
      fig2_igv_screenshot.png        # real IGV screenshot reused from images/tapestry.pedigree.png
    all_cpg_expansion/
      all_cpg_expansion.md         # Step 4 — expand_to_all_cpgs.py
      fig1_cpg_universe.png        # ref-genome CpGs vs sample CpGs vs measured CpGs
      fig2_allele_specific.png     # SNV creates/destroys CpG
      fig3_qc_mismatch.png         # mismatch-site proximity flag
    output_format/
      output_format.md             # BED schema + header + column dictionary (expanded from README)

  trio_wise_workflow/              # mostly stubbed — see §8
    README.md                      # stub: "most pages to be populated later; see README.md §Trio-wise workflow"
    output_format_trio/            # BUILT EARLY — holds the searchable-output mock-up, see §12
      output_format_trio.md        # column dictionary + searchable-output figure with mock data
      fig1_search_epi_compound.png # annotated mock BED table + polars queries
```

Step 2 (`aligned_bam_to_cpg_scores.sh`) is a pure pass-through of the
upstream pb-CpG-tools binary and does not get its own page; both
workflows can describe it in one paragraph inside their
`overview.md`. If later it becomes useful to dedicate a page to it
(e.g. to discuss count vs. model mode tradeoffs), that page would live
inside whichever workflow section is using it.

## 5. Page-by-page manifest

Each page below lists: its audience purpose, the tapestry code it
explains, and the figure(s) it contains. The wiki is a standalone
reference (see §1 goal 4); there is no explicit mapping from wiki
pages to manuscript figures, although §7 lists the pages the user is
likely to consult when drafting each paper figure.

### 5.1 `wiki/index.md` (top-level catalog)

- One-paragraph elevator pitch for tapestry.
- Two sections: "Pedigree-wise workflow" (→
  `pedigree_wise_workflow/index.md`) and "Trio-wise workflow" (→
  `trio_wise_workflow/README.md`, stub).
- "How to reproduce the figures in this wiki" —
  `python wiki/generate_wiki.py`.
- Cite the upstream Platinum-Pedigree-Inheritance wiki as the source of
  the `inheritance_mapping/` section.

### 5.2 `motivation/motivation.md` (shared across both workflows)

The hook page. Answers a reader's first question — *why phase DNA
methylation at all?* — with two cartoon figures built by
`generate_wiki.py` (details in §11). Both workflow sections link to
this page from their respective `overview.md`, and the top-level
`wiki/index.md` embeds `fig1_before_unphased.png` as its hero image.
Living in its own top-level folder rather than inside either workflow
keeps it shared and keeps the pedigree-wise and trio-wise sections
truly symmetric — the motivation is the same regardless of which
phasing strategy gets you to parental-haplotype resolution.

### 5.3 `pedigree_wise_workflow/index.md`

Catalog page scoped to the pedigree-wise workflow. Linked table of
contents, one row per page below, with a one-sentence hook each.

### 5.4 `pedigree_wise_workflow/overview/overview.md`

- CEPH1463 pedigree diagram (`fig1_pedigree.png`, re-drawn from
  `images/pedigree.jpg`).
- End-to-end workflow DAG (`fig2_workflow_dag.png`, re-rendered from
  `docs/pedigree_workflow.mmd` by piping it through `mmdc`).
- Text walk-through of the four steps at the top level, one paragraph
  per step, each linking to its dedicated page.

### 5.5 `pedigree_wise_workflow/read_backed_phasing/read_backed_phasing.md`

Covers Step 1A of the pedigree-wise workflow (`run-hiphase.sh`).

- What hiphase does: reads carrying multiple hets tie the hets' phase
  together into phase blocks; the HP/PS BAM tags propagate that phasing
  to every read overlapping the block.
- Toy figure 1 — 3 reads spanning 2 hets → phase block formation with a
  simulated read-to-allele assignment.
- Toy figure 2 — a BAM row with HP=1, PS=12345 annotations.
- Why two sibling phased VCFs (DeepVariant / pbsv) both get produced.
- Pointers to `run-hiphase.sh` and to the hiphase upstream docs.

### 5.6 `pedigree_wise_workflow/inheritance_mapping/*` (vendored)

Covers Step 1B of the pedigree-wise workflow
(`build-iht-based-haplotype-map-and-phase-variants.sh`). Nothing to
author here — see §3 for the sync strategy. A short
`inheritance_mapping/README.md` introduces the section and links to the
three sub-pages:

- `nuclear_family/nuclear_family.md` — `gtg-ped-map` on a two-founder
  three-child pedigree.
- `three_generations/three_generations.md` — extension to
  outside-marriage + grandchildren, depth-ordered ancestor-first
  processing.
- `concordance/concordance.md` — `gtg-concordance`'s `2^F` founder-phase
  enumeration for phasing every VCF variant inside an IHT block.

### 5.7 `pedigree_wise_workflow/founder_phased_methylation/founder_phased_methylation.md`

Covers Step 3 of the pedigree-wise workflow
(`phase_meth_to_founder_haps.py`). This is the **conceptual heart** of
tapestry and will be the largest new page.

Sections:

1. **Two kinds of phase, each local to its own block type.** Step 1A
   gives us `hap1 | hap2` on each sample's reads — locally phased within a
   hiphase phase block but unnamed (we don't know which of `hap1`/`hap2` is
   paternal vs. maternal). Step 1B gives us `pat | mat` on each sample's
   variants inside an iht block — locally directed by inheritance but only
   at variant sites, and the two halves (pat and mat) trace back to the
   founder letters. Tapestry's job is to cross-reference the two inside
   their *intersection*.
2. **Bit-vector matching across the block intersection — the hero
   cartoon.** Rather than three separate small figures, the core
   mechanism is shown as one end-to-end cartoon
   (`fig1_bit_vector_matching.png`, spec in §13) covering bit-vector
   extraction from both phasing sources, bit-vector comparison inside
   the `(hiphase-block, iht-block)` intersection, founder-haplotype
   labeling of the haplotagged reads, and phasing of methylation as a
   mechanical consequence. Source references:
   `src/hap_map_pedigree.py:13` (`get_hap_map`) and
   `src/phase_meth_to_founder_haps.py:44`
   (`phase_meth_to_founder_haps`). The detailed prose in this page
   walks the reader through the cartoon panel by panel (bit-vector
   extraction → concordance → founder label → methylation phasing).
3. **IGV visualization.** Reuse `images/tapestry.pedigree.png` as
   `fig2_igv_screenshot.png` to show what the output looks like in
   practice after the cartoon's abstract pipeline has been applied to
   real data.
4. **QC: mismatch sites.** Heterozygous SNV sites inside the intersection
   at which the bit-vectors disagree — these are exactly the positions
   where read-backed and inheritance-backed phasing conflict, and become
   the "within 50bp of a mismatch site" flag in Step 4's output. Annotated
   in the hero cartoon (Panel B) rather than given its own figure.

### 5.8 `pedigree_wise_workflow/all_cpg_expansion/all_cpg_expansion.md`

Covers Step 4 of the pedigree-wise workflow (`expand_to_all_cpgs.py`).

Sections:

1. **Why expand.** `aligned_bam_to_cpg_scores` only emits CpGs that had
   ≥ `min-coverage` reads. A full epigenetic atlas has to also
   acknowledge reference CpGs with no data and sample-only CpGs created
   by variants.
2. **Three populations of CpG sites.** Reference-genome CpGs (from
   `write_all_cpgs_in_reference.py`), sample-genome CpGs (inferred from
   the sample's variants), measured CpGs (emitted by pb-CpG-tools).
   Figure — `fig1_cpg_universe.png`, a three-circle Venn.
3. **Allele-specific CpGs.** A CpG dinucleotide that exists on only one
   of the two sample haplotypes (the other haplotype carries a SNV that
   destroys the CpG). Figure — `fig2_allele_specific.png`, a toy site
   showing how a single SNV toggles allele-specificity.
4. **QC: mismatch-site proximity.** Flag any CpG within 50 bp of a
   bit-vector mismatch site. Figure — `fig3_qc_mismatch.png`.

### 5.9 `pedigree_wise_workflow/output_format/output_format.md`

Lifts the column dictionary from `README.md` lines 121-156 and expands
each entry with a worked example and a pointer back to the page that
explains the column. No new figures.

### 5.10 `trio_wise_workflow/output_format_trio/output_format_trio.md` (built early)

Although the rest of the trio-wise section is stubbed (§8), this one
page is built as part of the first-pass wiki because its content is
already fully specified by `README.md` §Trio-wise-workflow and because
it hosts the searchable-output figure (§12), which is a high-value
asset for the manuscript. The page has two parts:

1. A column dictionary lifted from `README.md` lines 189–239, each
   entry expanded with a one-sentence "what you'd use this for"
   annotation.
2. A worked example — the §12 figure — showing a ~6-row mock BED
   table alongside polars-style query snippets that find candidate
   epimutations and compound genetic-epigenetic heterozygotes
   directly in the output.

## 6. `generate_wiki.py` — script architecture

Inherit the structure from the upstream `wiki/generate_wiki.py`:

- Top-level CLI: `--page <name>` to regenerate one page. Page names
  reflect the nested layout, e.g.
  `pedigree/overview`, `pedigree/read_backed_phasing`,
  `pedigree/founder_phased_methylation`,
  `pedigree/all_cpg_expansion`. (Short aliases are fine; the key point
  is that the page name is independent of any filesystem refactor of
  sibling workflow directories.) With no arg, regenerate everything
  *except* the vendored `inheritance_mapping/` pages (those are owned
  by upstream; see §3).
- One top-level function per page: `page_pedigree_overview()`,
  `page_pedigree_read_backed_phasing()`,
  `page_pedigree_founder_phased_methylation()`,
  `page_pedigree_all_cpg_expansion()`. Each function writes both the
  page's `*.md` narrative and its `fig*.png` bundle atomically into its
  own folder under `pedigree_wise_workflow/`.
- Shared helpers (`text_panel`, monospace layout, character-span
  highlighting) ported from upstream `generate_wiki.py`. If that module
  remains available in a sibling checkout of
  `Platinum-Pedigree-Inheritance`, consider vendoring those helpers into
  `wiki/_helpers.py` rather than re-deriving them.
- Stdlib + matplotlib only. No numpy dependency beyond what matplotlib
  pulls in. Simulations are hard-coded and deterministic — the same
  discipline as upstream — so figures remain byte-reproducible across
  reruns.
- Captions use the permalink helper that writes
  `https://github.com/quinlan-lab/tapestry/blob/<sha>/<path>#L<line>`
  (where `<sha>` is captured from `git rev-parse HEAD` at generate-time).
  For the `founder_phased_methylation` page specifically, the permalinks
  target `src/phase_meth_to_founder_haps.py` and
  `src/hap_map_pedigree.py`; for `all_cpg_expansion`, they target
  `src/expand_to_all_cpgs.py` and `src/write_all_cpgs_in_reference.py`.
- When the trio-wise section is built out (§8), add parallel
  `page_trio_*()` functions that write into
  `trio_wise_workflow/`. No existing page needs to move.

## 7. Using the wiki as a reference while drafting manuscript figures

The wiki is a reference, not a pipeline for manuscript assets. When the
user sits down to compose a figure in Illustrator, they will open the
relevant wiki pages in one window and use them to decide what the
figure needs to show and how to lay it out. The table below is a
suggested routing — "if you're drafting figure X, consult pages Y" —
not a machine-readable manifest:

| If drafting… | Consult |
|---|---|
| Abstract / introduction / graphical-abstract figure motivating the method | `motivation/` |
| Study design + pipeline overview figure | `pedigree_wise_workflow/overview/` |
| Method figure (bit-vector concordance, block intersection) | `pedigree_wise_workflow/founder_phased_methylation/` |
| Output-example figure (IGV-style phased methylation) | `pedigree_wise_workflow/founder_phased_methylation/` + `images/tapestry.pedigree.png` |
| CpG-coverage & allele-specific-CpG figure | `pedigree_wise_workflow/all_cpg_expansion/` |
| Supplementary inheritance-mapping methods | `pedigree_wise_workflow/inheritance_mapping/` (vendored) |
| Supplementary "how to query the output for epimutations / compound hets" figure | `trio_wise_workflow/output_format_trio/` (see §12) |
| Supplementary QC figures | `pedigree_wise_workflow/all_cpg_expansion/` + `images/QC.png` |

Consequences of this framing:

- No `manuscript_figures/` assembly step in `generate_wiki.py` and no
  renaming convention tying wiki PNGs to paper figure numbers (both
  would have rotted the moment the paper's layout shifted).
- The wiki does not need to be visually consistent with the paper — it
  optimises for pedagogical clarity (monospace ASCII-style panels,
  line-numbered source permalinks). The paper figures optimise for a
  different audience and will look different.
- If a wiki panel happens to be exactly what a paper figure needs, it
  is fine to open the PNG in Illustrator and trace/copy it. But that
  is an authoring-time decision, not a build-time contract.

## 8. Trio-wise workflow — future extension

The sibling directory `wiki/trio_wise_workflow/` exists from day one
(see §4) and is *mostly* stubbed. The one exception is
`output_format_trio/`, which is built as part of the first-pass wiki
— its content (the trio output BED column dictionary) is already
fully specified by `README.md`, and it is the natural home for the
searchable-output mock-up figure (§12) that the user wants available
for manuscript drafting immediately. The remaining folders below are
stubs until the trio-wise side is properly built out:

```
wiki/trio_wise_workflow/
  index.md
  pedmec_phasing/              # run-whatshap.sh
  parent_haplotype_phasing/    # phase_meth_to_parent_haps.py
  all_cpg_expansion_trio/      # expand_to_all_cpgs.trio.sh
  output_format_trio/          # columns dictionary lifted from README §trio workflow
```

Since the pedigree-wise `founder_phased_methylation` page already
establishes bit-vector concordance in the general setting,
`parent_haplotype_phasing` will be shorter — it just has to explain the
two key differences: (a) the phase source is pedMEC/whatshap rather
than `gtg-concordance`, and (b) the kid's haplotypes are labelled by
*parental* letters (A/B in dad, C/D in mom) rather than by
*founder-of-the-pedigree* letters. Note that A/B/C/D here are fixed as
dad's hap1/hap2 and mom's hap1/hap2 — not defined as transmitted vs.
non-transmitted.

The top-level `wiki/index.md` catalog already has slots for both
workflows. Populating the trio-wise side will not move or rename
anything on the pedigree-wise side.

## 9. Phased implementation plan

Phased so each phase leaves a site that renders cleanly in a GitHub web
view — useful for showing intermediate progress to collaborators without
waiting for the whole wiki to be ready.

**Phase 1 — scaffolding and vendored section (half day).**

- Create `wiki/index.md` with two section stubs (pedigree-wise →
  `pedigree_wise_workflow/index.md`, trio-wise →
  `trio_wise_workflow/README.md`).
- Create `wiki/pedigree_wise_workflow/index.md` as a stub with links to
  empty page folders.
- Create `wiki/trio_wise_workflow/README.md` as a stub per §8.
- Vendor the upstream `inheritance_mapping/` section per §3 (into
  `pedigree_wise_workflow/inheritance_mapping/`).
- Verify every upstream cross-link resolves on GitHub's web view after
  the one-level-deeper nesting.
- Commit.

**Phase 2 — `motivation/` page (half day).**

- Build the two cartoon figures per the spec in §11 via
  `generate_wiki.py`.
- Write the `motivation.md` narrative (before/after framing; brief
  definitions of epimutation and compound genetic-epigenetic
  heterozygote; link forward to each workflow's `overview.md`).
- Wire `wiki/index.md` to embed `fig1_before_unphased.png` as its hero
  image with a caption that links to the full motivation page.
- Commit.

**Phase 3 — `overview/` and `output_format/` (half day).**

- `pedigree_wise_workflow/overview/overview.md`: re-render
  `docs/pedigree_workflow.mmd` to PNG and write the narrative. Pedigree
  figure derives from `images/pedigree.jpg` (can start as a direct
  embed of that image; replace with a cleaner re-drawing later if
  needed). Link back to `motivation/motivation.md` in the first
  paragraph.
- `pedigree_wise_workflow/output_format/output_format.md`: lift the
  column dictionary from `README.md` and expand it.
- Commit.

**Phase 4 — `read_backed_phasing/` page (half day).**

- Narrative + two toy figures from `generate_wiki.py`.
- Minimal math content — hiphase is well-documented upstream, so this
  page mostly anchors a reader and hands off.
- Commit.

**Phase 5 — `founder_phased_methylation/` page (one to two days).**

This is the content-heaviest page. Build it component-at-a-time per
the user's documented preference for incremental review, where the
components map to the four panels of the hero cartoon described in
§13.

- Component 1 (Panel A): two haplotagged read pileups + their bit
  vectors (hiphase hap1/hap2 on the left, iht pat/mat on the right).
  Sets up the two phasing sources side by side on a shared genomic
  axis.
- Component 2 (Panel B): block intersection + bit-vector comparison
  with ✓/✗ marks and a similarity readout; annotate any mismatch
  site as the QC flag for Step 4.
- Component 3 (Panel C): read pileup re-coloured and relabeled with
  founder-haplotype tags (dad_hap1, mom_hap2, etc.) derived from the
  similarity winner.
- Component 4 (Panel D): methylation bars on the shared 0–100 %
  y-axis (same helper as §11.5) showing how the HP-bucketed
  methylation streams are mechanically relabeled into `pat`/`mat`
  once the read→founder map is in hand.
- Component 5: IGV screenshot panel (`fig2_igv_screenshot.png`,
  reuses `images/tapestry.pedigree.png`) closing the page with a
  real-data example.
- Narrative wired around all four cartoon panels + the IGV screenshot,
  with permalinks to `src/phase_meth_to_founder_haps.py` and
  `src/hap_map_pedigree.py`.
- Commit.

**Phase 6 — `all_cpg_expansion/` page (one day).**

- Three CpG-universe figures (Venn, allele-specific toggle, mismatch
  proximity).
- Narrative wired around them, with permalinks to
  `src/expand_to_all_cpgs.py` and `src/write_all_cpgs_in_reference.py`.
- Commit.

**Phase 7 — `trio_wise_workflow/output_format_trio/` page (half day).**

- Lift the trio BED column dictionary from `README.md` lines 189–239.
- Build the §12 figure (annotated mock BED table + polars queries for
  epimutations and compound hets) via `generate_wiki.py`.
- Write the narrative that frames the figure as "once tapestry
  produces its output, here is how you go looking for the phenomena
  motivated in `motivation/motivation.md`."
- Cross-link from `motivation/motivation.md`'s closing paragraph to
  this page, so the motivation → output story arc is walkable.
- Commit.

**Deferred — rest of `trio_wise_workflow/`.** All other trio-wise
pages remain stubs per §8. (The original "Phase 6 — manuscript
assembly script" is dropped; see §1 and §7 — the wiki is a reference,
not a source feed for the paper.)

## 10. Open questions to confirm before starting Phase 1

1. **Is the user's fork (`petermchale/Platinum-Pedigree-Inheritance`) the
   right upstream for the `inheritance_mapping/` section right now**, or
   should the sync point at `Platinum-Pedigree-Consortium/...` once the
   PR merges?
2. **Commit-SHA pinning for tapestry permalinks**: pin to the current
   `HEAD` at generate time, or to a named release tag? Upstream uses a
   specific commit SHA — matching that convention is probably right.
3. **Monospace matplotlib panels vs. diagram-style panels.** Upstream
   uses monospace for the ASCII-style walkthroughs and that is the
   right default. But Bioinformatics readers may expect a more
   conventional schematic style for Fig 1 (pedigree + DAG). Plan: use
   monospace for the inheritance-mapping section (already vendored) and
   for any step-by-step tapestry computations; use conventional
   boxes-and-arrows (via mmd → PNG) for Fig 1 and the block-intersection
   illustration. Confirm.
4. **Do we want the wiki browsable as a GitHub Pages site** (needs
   `_config.yml`, URL-prefix rewriting in links, maybe a Jekyll theme)
   or purely as navigable markdown on github.com (nothing extra
   required)? The upstream wiki uses the latter; default to that unless
   the user asks otherwise.
5. **Re-drawing `images/pedigree.jpg`.** The existing JPG is 3.7 MB and
   looks like a photograph of a printed pedigree. For a manuscript
   figure it will need to be a clean vector re-draw. Question: can we
   base a re-draw on an existing crowdsource pedigree diagram
   (e.g. the one in Eberle 2017), or do we need to draw from scratch?

## 11. Motivation figure spec

Three cartoon figures, all produced by `generate_wiki.py` into
`wiki/motivation/`. All three are fully synthetic — no real reads, no
real sample. The goal is to make a reader, on first contact with the
wiki, understand *why* haplotype-level methylation is worth computing.
The page is shared between the pedigree-wise and trio-wise workflows
(see §5.2) because the motivation does not depend on which phasing
strategy you use to reach parental-haplotype resolution.

The three figures play distinct roles:

1. **Figure 1 (single individual, before phasing).** Read pileup with
   one *pooled* methylation profile on top, illustrating that an
   unphased per-CpG number cannot resolve haplotype-level structure.
2. **Figure 2 (same individual, after phasing).** Same reads now
   coloured by haplotype, with *two* methylation profiles on top
   (one per haplotype). Sells the mechanical claim of tapestry: the
   read partition splits one pooled profile into two per-haplotype
   profiles.
3. **Figure 3 (trio, no reads).** Pedigree-style layout — dad, mom,
   kid — with each individual's two haplotypes drawn as horizontal
   coloured lines carrying methylation lollipops above and SNV allele
   bits below. Four founder colours; the kid's two haplotype lines
   inherit the colour of the parental homolog they descend from. This
   is where the *epimutation* and *compound genetic-epigenetic
   heterozygote* use cases become concrete, because both require a
   trio (epimutation: kid vs same-physical-homolog parent;
   compound-het: visible from a single trio member but only when
   phased into the parents' haplotype labels).

The visual style for the read pileups in Figs 1 and 2 borrows
deliberately from
`/Users/petermchale/phasing_simulations/simulate_reads_single_sample.py`
— monospace glyphs at fixed columns, ragged read starts/ends rendered
as `-` for "no information at this site". The lollipop-style trio in
Fig 3 borrows from the multi-generation methylation-tracking diagram
the user attached when redesigning this page (open / filled circles on
sticks above each haplotype line).

### 11.1 Figure 1 — "Before" (single individual, unphased)

`wiki/motivation/fig1_before_unphased.png`. Two stacked subplots
sharing the genomic x-axis:

- **Top subplot — pooled methylation profile, 0–1 axis.** A single
  series of vertical stems with markers at the top, one stem per CpG.
  Both stems sit at exactly 0.5 in the toy — the least-informative
  possible value.
- **Bottom subplot — read pileup.** ~10 ragged reads, each a row of
  monospace glyphs at fixed x = site-index columns. Sites inside the
  read's window carry the read's allele bit (`0`/`1`) at SNV columns
  and methylation glyph (`●`/`○`) at CpG columns; sites outside the
  window are rendered as a faint `-` so the ragged extents read
  cleanly. Every glyph is drawn in a single uniform grey — the
  haplotype partition is not yet known. Reads are interleaved so the
  unphased picture does not visually pre-suggest the partition.

Locus layout, left-to-right: `SNV1`, `CpG1`, `SNV2`, `CpG2`, `SNV3`.
Three SNVs (instead of two) so the partition in Fig 2 is visibly
robust against a single sequencing error at any one of them. CpG1
sits immediately next to SNV1 to set up the compound-het framing
in Fig 3; CpG2 is mid-locus with no adjacent SNV so it is a pure
methylation phenomenon (the epimutation in Fig 3).

### 11.2 Figure 2 — "After" (same individual, phased)

`wiki/motivation/fig2_after_phased.png`. Same two-subplot layout as
Fig 1, but:

- **Top subplot — two methylation profiles**, one per haplotype, on
  the same 0–1 axis. Paternal stems in dark blue and maternal in dark
  red, slightly x-offset so they don't overlap. In the toy: paternal
  CpG1 = 0 / CpG2 = 1; maternal CpG1 = 1 / CpG2 = 0. Each Fig 1 stem
  splits into two stems on opposite extremes, in *opposite directions*
  across the two CpGs.
- **Bottom subplot — partitioned read pileup.** Same reads as Fig 1,
  now coloured by haplotype (paternal = dark blue, maternal = dark
  red) and grouped vertically (pat above mat) with a faint divider.
  Side labels "paternal reads" / "maternal reads" on the left.

The colours in Fig 2 deliberately match the colours of the kid's two
haplotypes in Fig 3 (paternal = dark blue = dad hap A; maternal = dark
red = mom hap C), so a reader moving from Fig 2 to Fig 3 can trace the
kid's two homologs by colour into the trio context.

### 11.3 Figure 3 — Trio (epimutations and compound-het use cases)

`wiki/motivation/fig3_trio_methylation.png`. A small pedigree drawn on
a single matplotlib axes:

- **Pedigree symbols.** Square for dad (top-left), circle for mom
  (top-right), square for kid (centred below), connected by a couple
  line + descent line in standard pedigree convention.
- **Per-individual haplotype displays.** Each individual has two thick
  horizontal coloured lines (the two haplotypes) stacked vertically
  below their pedigree symbol. Each line carries:
    - **Methylation lollipops above the line** at CpG positions.
      Filled circle = methylated; open circle = unmethylated.
    - **SNV allele bits below the line** at SNV positions, rendered
      in the haplotype's colour.
- **Four-colour founder palette.** Dad hap A (transmitted to kid) =
  dark blue; dad hap B (not transmitted) = light blue; mom hap C
  (transmitted to kid) = dark red; mom hap D (not transmitted) =
  light orange. The kid's pat haplotype line uses the same dark blue
  as dad hap A; the kid's mat haplotype line uses the same dark red as
  mom hap C. The colour-match is the visual cue that says "this is the
  *same physical homolog*, traced from parent to child."
- **Epimutation annotation at CpG2.** Two dashed boxes — one around
  dad hap A's CpG2 lollipop, one around the kid's pat CpG2 lollipop —
  plus a text label in the open space below mom's haplotype display
  with an arrow pointing to the kid's box. Kid pat at CpG2 is
  methylated; dad hap A at CpG2 is unmethylated; same physical
  homolog (matched colours) → de novo gain of methylation in meiosis.
- **Compound-het annotation at CpG1.** A single dashed box covering
  both kid haplotype lines around the SNV1 + CpG1 region (encompasses
  the two SNV1 allele bits and the two CpG1 lollipops on the kid's
  two haplotype lines). Text label in the open space to the left of
  the kid: methylation state at CpG1 co-segregates with SNV1 genotype
  across the two kid haplotypes → compound genetic-epigenetic
  heterozygote.

### 11.4 Toy simulation details

Fully deterministic; no RNG. Site index 0..4 = SNV1, CpG1, SNV2, CpG2,
SNV3.

- **Founder haplotypes** (per-site values; for SNVs `0` = REF, `1` =
  ALT; for CpGs `0` = unmethylated, `1` = methylated):
    - `dad_A = (0, 0, 0, 0, 0)` — transmitted to kid
    - `dad_B = (1, 1, 1, 0, 1)`
    - `mom_C = (1, 1, 1, 0, 1)` — transmitted to kid
    - `mom_D = (0, 0, 0, 0, 0)`
- **Kid haplotypes:**
    - `kid_pat = (0, 0, 0, 1, 0)` — equals `dad_A` except CpG2 is now
      methylated (the epimutation: de novo gain of methylation in the
      meiosis from dad to kid)
    - `kid_mat = (1, 1, 1, 0, 1)` — equals `mom_C`
- **Read simulation (Figs 1 & 2, kid only).** 10 ragged reads, 5 from
  pat and 5 from mat, with windows chosen so both CpGs have at least
  3 reads from each haplotype after partitioning:
  ```
  pat: (0,2), (1,4), (0,3), (0,4), (2,4)
  mat: (0,3), (0,4), (2,4), (1,4), (0,2)
  ```
  Pooled methylation: CpG1 = 0.5, CpG2 = 0.5 (Fig 1's two stems).
  Phased: pat CpG1 = 0, pat CpG2 = 1, mat CpG1 = 1, mat CpG2 = 0
  (Fig 2's four stems).
- **Compound-het visibility.** Kid is heterozygous at every SNV
  (pat = `0`, mat = `1`); the methylation state at the
  immediately-adjacent CpG1 co-segregates with the SNV1 genotype
  (pat = `○`, mat = `●`). This is the dashed-box annotation in
  Fig 3.
- **Epimutation visibility.** Dad's hap A and the kid's pat homolog
  are the same physical homolog (transmitted unchanged in the
  meiosis), so they share a colour in Fig 3. They differ at CpG2: dad
  A is `○`, kid pat is `●`. This is the dashed-box + arrow annotation
  in Fig 3.

### 11.5 Where the figures live in the reader's path

- `wiki/index.md` embeds `fig1_before_unphased.png` as a hero image at
  the top of the page, with a one-line caption and a link to
  `motivation/motivation.md` for the full story.
- `motivation/motivation.md` contains all three figures in order
  (Fig 1 → Fig 2 → Fig 3) surrounded by the narrative above. Its
  closing paragraph links forward to both workflow overviews and to
  the trio-wise output_format searchable example (§12) for "once you
  have these per-haplotype profiles in hand, here is how you query
  for the two phenomena".
- Each workflow's `overview.md` has a back-reference in its opening
  paragraph: *"For motivation, see [why phase DNA
  methylation?](../../motivation/motivation.md)."*

### 11.6 Implementation notes

- All three figures are emitted by `page_motivation()` in
  `wiki/generate_wiki.py`, selectable via
  `python wiki/generate_wiki.py --page motivation`.
- The read pileup in Figs 1 and 2 is rendered as monospace
  `pyplot.text` glyphs at integer site-index x-coordinates with a
  faint coloured backbone line for each read's window (so out-of-
  window cells, drawn as `-`, are visually distinguishable from
  in-window cells without overloading the glyph alphabet).
- The methylation-profile subplots above the read pileups use a
  shared 0–1 y-axis with identical ticks (`0`, `0.5`, `1`).
  `_draw_methylation_profile()` in `generate_wiki.py` is the helper
  that enforces the axis style.
- The trio in Fig 3 is rendered on a single axes with manual
  positioning of the three "person cells" (pedigree symbol +
  haplotype display). Lollipops are `pyplot.scatter` markers above a
  short stem line; the haplotype backbone is a single `plot()` call
  per haplotype with `linewidth=6.5` so the line reads as a clean
  bar even with the lollipop sticks crossing it.
- Because the simulation is hard-coded, no unit tests are required;
  correctness is checked by eye on the rendered PNG.
- The colour palette for Figs 2 and 3 is intentionally consistent
  (paternal = dark blue, maternal = dark red, with kid's haplotype
  lines in Fig 3 reusing those exact colours), so a reader's eye can
  trace a single haplotype across figures.

## 12. Searchable-output figure (trio-wise BED mock-up)

A figure that closes the conceptual loop opened in §11: once the
reader believes that epimutations and compound genetic-epigenetic
heterozygotes are interesting, they will want to know how to *find*
them in a tapestry output BED. This figure shows that the relevant
trio-wise columns in the output make both searches a small filter
expression, not a custom pipeline.

The figure lives at
`wiki/trio_wise_workflow/output_format_trio/fig1_search_epi_compound.png`
and is referenced from `output_format_trio.md` (see §5.10) and from
the closing paragraph of `motivation/motivation.md`.

### 12.1 What the figure shows

Two stacked panels on a single PNG:

**Panel A — mock BED table.** A rendered table with a header row and
6–7 data rows of hard-coded mock data. Column widths are set so the
table fits a single figure width without wrapping. Only the subset of
trio-output columns relevant to the two searches is shown; the rest
are elided with a trailing "…" column so the reader knows more exist
and can find them in the dictionary above.

Columns included (abbreviated headers in the figure; full names are
the real BED column names from `README.md` lines 189–239):

| Abbreviated header | Full BED column name |
|---|---|
| `chrom` / `start_cpg` | `chrom`, `start_cpg` |
| `kid_pat_m` | `methylation_level_kid_pat_count` |
| `kid_mat_m` | `methylation_level_kid_mat_count` |
| `pat_hap` | `paternal_haplotype` (`A` or `B`) |
| `dad_A_m` / `dad_B_m` | `methylation_level_dad_A_count`, `methylation_level_dad_B_count` |
| `mat_hap` | `maternal_haplotype` (`C` or `D`) |
| `mom_C_m` / `mom_D_m` | `methylation_level_mom_C_count`, `methylation_level_mom_D_count` |
| `kid_SNV` | `snv_genotypes_kid` (`hom` or `het`) |
| `AS_kid` | `cpg_is_allele_specific_kid` (`T`/`F`) |

Mock rows, each with a coloured left-gutter label tagging what
phenomenon (if any) the row illustrates:

1. *(no tag)* `chr1 1000` — boring: `kid_pat_m=0.05`, `kid_mat_m=0.03`,
   `pat_hap=A`, `dad_A_m=0.04`, `mat_hap=C`, `mom_C_m=0.03`,
   `kid_SNV=hom`, `AS_kid=F`. Both haplotypes unmethylated and
   concordant with parents. Included to show that the vast majority
   of rows are unremarkable.
2. *(no tag)* `chr1 1500` — boring: all methylation levels ≈ 0.9, no
   SNV overlap, not allele-specific. The other "boring" baseline.
3. **Paternal epimutation** `chr1 2000` — `kid_pat_m=0.95`,
   `pat_hap=A`, `dad_A_m=0.05` (kid's paternal haplotype is fully
   methylated but dad's *same physical homolog* is not — methylation
   appears to have been acquired during transmission). Maternal
   side concordant: `kid_mat_m=0.10`, `mat_hap=C`, `mom_C_m=0.08`.
4. **Maternal epimutation** `chr1 2500` — symmetric: `kid_mat_m=0.05`,
   `mat_hap=D`, `mom_D_m=0.95`. Paternal side concordant.
5. **Compound genetic-epigenetic heterozygote (kid)** `chr1 3000` —
   `kid_pat_m=0.98`, `kid_mat_m=0.03`, `kid_SNV=het`, `AS_kid=F` (the
   CpG exists on *both* of the kid's haplotypes — it is not destroyed
   by the SNV; the two haplotypes just carry opposite methylation
   states, and that state tracks the SNV genotype). Dad/mom methylation
   values consistent with transmission, so this is not also an
   epimutation — just a stable hap-specific methylation state coupled
   to a heterozygous SNV.
6. **Allele-specific CpG (kid)** `chr1 3500` — `kid_pat_m=0.90`,
   `kid_mat_m=` *(null, shown as `NA` or a grey dash)*, `kid_SNV=het`,
   `AS_kid=T`. The SNV destroys the CpG dinucleotide on the kid's
   maternal haplotype, so no methylation level can be reported there.
   Included to make the distinction from row 5 concrete: row 5 has
   the CpG on both haps with different methylation states; row 6 has
   the CpG on only one hap.

Each "tagged" row (3, 4, 5, 6) gets a coloured left-gutter stripe
plus a small coloured text label (e.g. "pat epi", "mat epi", "cGEH",
"AS CpG"). Rows 1 and 2 get no stripe. The colour palette is
consistent with the rest of the wiki's matplotlib style.

**Panel B — search queries that find the tagged rows.** Two short
polars snippets, each with a one-line comment naming what it finds.
The snippets are real, copy-pasteable polars filter expressions
(the user can lift them into their own notebook). Approximate form:

```python
# Candidate epimutations (one haplotype in kid ≠ same physical hap in parent)
epi_pat = df.filter(
    ((pl.col("paternal_haplotype") == "A") &
     ((pl.col("methylation_level_kid_pat_count")
       - pl.col("methylation_level_dad_A_count")).abs() > 0.8))
    | ((pl.col("paternal_haplotype") == "B") &
       ((pl.col("methylation_level_kid_pat_count")
         - pl.col("methylation_level_dad_B_count")).abs() > 0.8))
)
# A symmetric epi_mat expression uses maternal_haplotype with mom_C / mom_D.

# Candidate compound genetic-epigenetic heterozygote in the kid
compound_het_kid = df.filter(
    (pl.col("snv_genotypes_kid") == "het")
    & (pl.col("cpg_is_allele_specific_kid") == False)
    & ((pl.col("methylation_level_kid_pat_count")
        - pl.col("methylation_level_kid_mat_count")).abs() > 0.8)
)
```

The thresholds (`> 0.8`) are placeholders; the figure caption notes
that real analyses will tune these against a realistic mismatch /
noise budget and a multiple-testing correction.

### 12.2 Why this figure belongs to the trio-wise workflow (not the pedigree-wise one)

The pedigree-wise output BED has one sample's point of view — it
records `methylation_level_pat` and `methylation_level_mat` for a
sample and carries founder-haplotype labels (`founder_haplotype_pat`,
`founder_haplotype_mat`) but does *not* include the corresponding
parent's per-haplotype methylation values in the same row. To detect
an epimutation from the pedigree-wise output alone you would have to
join a kid's BED with the corresponding parent's BED on
`founder_haplotype_*`. Doable, but it adds a join step that the
figure would have to document.

The trio-wise output already bakes that join in: each row reports
the kid's `pat`/`mat` methylation, both of dad's `A`/`B` methylation
values, both of mom's `C`/`D` methylation values, and the pointer
columns (`paternal_haplotype` ∈ {A, B}, `maternal_haplotype` ∈ {C,
D}) that say which of dad's or mom's haplotypes the kid inherited in
that block. So a trio-wise row is a direct expression of "here is
everything needed to compare the kid's methylation on a haplotype to
the parent's methylation on the *same physical homolog*." That makes
the figure a perfect natural fit for the trio-wise output schema and
the cleanest possible worked example.

A follow-up figure for the pedigree-wise output page (§5.9) can be
added later if needed: its content would be the same mock BED style
but with the explicit join step visible.

### 12.3 Narrative placement and cross-linking

- `trio_wise_workflow/output_format_trio/output_format_trio.md`
  hosts the figure in its "Worked example: finding epimutations and
  compound hets" section. Above the figure, the page walks through
  one row of the mock table in prose, reading out each column so the
  reader sees what the numbers mean before the figure condenses them.
- `motivation/motivation.md`'s closing paragraph links forward: *"For
  a worked example of how to query tapestry's output for the
  phenomena shown above (epimutations, compound genetic-epigenetic
  heterozygotes), see
  [`trio_wise_workflow/output_format_trio/output_format_trio.md`](../trio_wise_workflow/output_format_trio/output_format_trio.md)."*
- `wiki/index.md` adds a one-line bullet under "What can I do with
  this?" pointing at the same page.

### 12.4 Toy simulation / implementation notes

- The mock table is hard-coded as a list of dicts inside
  `page_output_format_trio()` in `generate_wiki.py` and rendered as
  a `matplotlib.table.Table` (or as monospace text in a text panel,
  matching the rest of the wiki's style). Bake the data in — no RNG,
  no external sample — so the figure is byte-reproducible.
- The left-gutter stripes + tag text are drawn as matplotlib patches
  and `pyplot.text` sitting in a narrow column to the left of the
  table.
- Panel B's polars snippets are rendered as text in a monospace
  matplotlib panel, *not* as live-executed Python — the snippets are
  illustrative; the page's role is to show the query shape, not to
  run it against real data.
- Written as `page_output_format_trio()` inside `generate_wiki.py`,
  selectable via `python wiki/generate_wiki.py --page output_format_trio`.

## 13. Bit-vector matching cartoon (pedigree-wise Step 3 hero figure)

A single end-to-end cartoon illustrating the mechanism implemented by
`phase_meth_to_founder_haps.py` + `src/hap_map_pedigree.py`: how the
sequence of alleles ("bit vector") recovered from hiphase's
single-sample read-backed phasing is matched against the sequence of
alleles recovered from the Rust inheritance-based phasing, inside the
intersection of a hiphase phase block and an iht block, and how that
match lets every haplotagged read (and therefore every per-CpG
methylation measurement) be re-labeled with the sample's founder
haplotype.

This is the pedagogical centre of the pedigree-wise workflow. It
replaces the three previously-separated figures in §5.7
(`fig1_bit_vector_intuition`, `fig2_hap_map_block`,
`fig3_founder_phased_cpg`) with a single coherent artifact so a
reader sees the full chain at once rather than three disconnected
fragments.

File: `wiki/pedigree_wise_workflow/founder_phased_methylation/fig1_bit_vector_matching.png`.

### 13.1 What the cartoon shows

Four stacked panels on a single PNG, sharing a horizontal genomic
coordinate axis so the panels are visually aligned.

**Panel A — two phasing sources, each with its own bit vector at the
same het-SNV sites.** Left half renders the hiphase output: ~6
haplotagged reads grouped by HP tag (hap1 reads on top, hap2 reads
below, with a faint horizontal divider matching §11's pat-above-mat
convention in spirit — here hiphase's `hap1`/`hap2` are *unnamed*
haplotypes). At each of ~5 het-SNV columns, each read carries a `0`
or `1` bit under the REF/ALT convention (matching the wiki's
site-wide symbology, per §11). Below the pileup, the two collapsed
bit vectors appear as
`hap1: 0 1 0 1 1` and `hap2: 1 0 1 0 0`. Above the left half, a
labeled bracket marks the hiphase phase block extent.

Right half renders the gtg-concordance output: the same five
het-SNV sites, each with a phased genotype `p|m` where `p` is the
paternal allele and `m` is the maternal allele. Below, two
collapsed bit vectors `pat: 0 1 0 1 0` and `mat: 1 0 1 0 1`
(deliberately chosen to agree with hap1 at four of the five sites
and disagree at one site — see Panel B). Above, a labeled bracket
marks the iht block extent and lists its founder-haplotype labels
(e.g. `founders: dad_hap1 dad_hap2 mom_hap1 mom_hap2`).

**Panel B — block intersection + bit-vector match.** A single
genomic axis shows the hiphase phase block and the iht block overlaid,
with their intersection shaded. The het-SNV sites inside the
intersection are the ones the comparison runs over.

Below, align `hap1` against `pat` character-by-character with ✓/✗
marks: `0 1 0 1 1` vs `0 1 0 1 0` → ✓✓✓✓✗ → 4/5 = 0.8 similarity.
Also show the complementary comparison `hap1` vs `mat` → 1/5 = 0.2
similarity, to make clear that the two options are mutually
exclusive. The one mismatching site is circled with a small tag
"flagged for QC — becomes `cpg_is_within_50bp_of_mismatch_site` in
Step 4's output."

Panel caption: "Because similarity(`hap1`, `pat`) > 0.5, the
sample's `hap1` reads descend from the paternal founder homolog
tagged `A` (or `dad_hap1`, in the example above) for this iht
block; `hap2` reads descend from the maternal one. If the
similarity had come out < 0.5, the assignment flips. This is the
single binary decision that ties read-backed phasing to
inheritance-based phasing."

**Panel C — founder-haplotype labeling of the haplotagged reads.**
Redraw Panel A's hiphase read pileup, but each read now carries an
additional tag: `dad_hap1` for the four former-`hap1` reads,
`mom_hap2` for the former-`hap2` reads (or whatever founder labels
the iht block's header specified). Use the paternal/maternal colour
convention already established by §11 (blue = paternal, red =
maternal). The reads and their allele bits are otherwise unchanged
from Panel A — the only new information is the founder label per
read, inferred via Panel B.

**Panel D — methylation phasing, mechanical.** Add a single CpG
column to the right of the alignment (or to a nearby position, as
long as it is covered by all the reads) with a `●`/`○` glyph on
each read showing methylation status on that read. Because the
reads are already bucketed into `hap1.bed.gz` and `hap2.bed.gz`
files by `aligned_bam_to_cpg_scores` based on the BAM's HP tag, and
Panel C has now attached a founder label to each HP bucket, the
methylation levels are trivially relabeled: count methylated reads
per founder haplotype and divide by coverage. Render the two
resulting methylation levels as two bars on the **same 0–100 %
y-axis used in §11's bars**, labeled `methylation_level_pat` and
`methylation_level_mat` in that top-to-bottom order (matching §11's
pat-above-mat convention). In the toy, say 4/4 methylated on `pat`
and 0/2 on `mat`, giving bars at 100 % and 0 %.

**Overall caption (under Panel D).** "Read-backed phasing
(hiphase) locally phases each sample's reads into `hap1`/`hap2` but
does not name which is paternal vs. maternal. Inheritance-based
phasing (`gtg-ped-map` + `gtg-concordance`) globally phases variant
sites into `pat`/`mat` using pedigree inheritance. Matching the two
bit vectors inside the intersection of their blocks labels every
read with its founder haplotype — and, because the methylation
stream is already bucketed by HP, this re-labels every per-CpG
methylation level as `pat` or `mat` at no extra cost. This is the
mechanism at the heart of Step 3 of tapestry's pedigree-wise
workflow."

### 13.2 Toy simulation details

Fully deterministic; no RNG. Concretely:

- 6 reads from hiphase, split 4 / 2 between the two HP tags: 4 reads
  tagged `hap1`, 2 reads tagged `hap2`. (An unequal split is chosen
  deliberately so the paternal methylation bar in Panel D reads 100 %
  — i.e., all four `pat` reads methylated — rather than a round 50 %
  that could be confused with an unphased aggregate.)
- 5 het-SNV sites inside the chosen `(hiphase-block, iht-block)`
  intersection, plus 1 CpG site covered by all 6 reads.
- hiphase bit vectors (extracted from the 4+2 read pileup): `hap1 = 0
  1 0 1 1`, `hap2 = 1 0 1 0 0`.
- iht (gtg-concordance) bit vectors: `pat = 0 1 0 1 0`, `mat = 1 0 1
  0 1`.
- similarity(`hap1`, `pat`) = 4/5 = 0.8, similarity(`hap1`, `mat`) =
  1/5 = 0.2. Decision: `hap1` reads → paternal founder, `hap2` reads
  → maternal founder.
- Site 5 is the deliberate mismatch; it is circled in Panel B and
  labeled as the QC-flagged site.
- At the CpG, the 4 paternal-assigned reads are all `●` (methylated);
  the 2 maternal-assigned reads are both `○` (unmethylated). Phased
  methylation levels: `pat = 4/4 = 100 %`, `mat = 0/2 = 0 %`.
- Founder labels for the iht block are hard-coded as
  `dad_hap1 dad_hap2 mom_hap1 mom_hap2`, and the sample's paternal
  founder label for this block is set to `dad_hap1`; Panel C's read
  annotation reads `dad_hap1` on the four paternal-assigned reads and
  `mom_hap2` on the two maternal-assigned reads.

### 13.3 Implementation notes

- Single `page_founder_phased_methylation()` function in
  `generate_wiki.py` that emits the PNG plus the
  `founder_phased_methylation.md` narrative.
- Reuses the monospace text-panel helpers from upstream
  `generate_wiki.py` for all four panels. Read bars are matplotlib
  patches; bit-vector rows are `pyplot.text`; methylation bars use
  the shared `draw_methylation_bars(ax, labels, fractions)` helper
  from §11.5 so the Panel D bars are on the same y-axis as Figs 1
  and 2 of the motivation page (directly comparable across the
  whole wiki).
- The shared genomic axis across the four panels is produced by
  laying out the panels with `gridspec` on a single matplotlib
  figure, with the axis extent passed in from a single
  `(genomic_start, genomic_end)` tuple.
- Written as a component-at-a-time build per §9 Phase 5 so the user
  can review each panel in isolation before the four are stitched
  into the final PNG.

### 13.4 Relationship to the existing `docs/bit_vector_toy_example.mmd`

The existing mermaid file `docs/bit_vector_toy_example.mmd` is a
smaller, earlier sketch of exactly Panel B of this cartoon (the
bit-vector match + concordance calculation, without the read pileups
on either side or the downstream founder-labeling / methylation
phasing). It can be kept as-is in `docs/` for anyone who wants the
condensed version — or retired after the cartoon ships, since the
cartoon supersedes it. Default plan: keep both (the mermaid file is
tiny, costs nothing) until the user explicitly decides to retire it.

---

*Owner: plan authored by Claude; implementation to be performed
interactively against this plan, component at a time, per the user's
documented preference for incremental review.*
