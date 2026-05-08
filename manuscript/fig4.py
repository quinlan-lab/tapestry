"""Render Fig 4 for the tapestry manuscript.

Fig 4 collects the founder-phased-methylation machinery and the trio
discovery tables. These panels were previously B–E of Fig 3 and have
been relettered A–D for Fig 4. The exhaustive-enumeration table at
clean site N1 remains in Fig 3 (panel A there).

Panel A — `manuscript/fig4/fig4_panelA.pdf`
    The bit-vector match cartoon (hap1 vs pat/mat over a hap-map block).
    Wiki source:
      `wiki/pedigree_wise_workflow/founder_phased_methylation/founder_phased_methylation.md`
      section "Bit-vector match", which embeds `bit_vector_match.png`.
    Rendering is delegated to that wiki page's own generator
    (`wiki/.../founder_phased_methylation.py:render_match`), so panel A
    is byte-equivalent to the wiki figure modulo the PDF container.

Panel B — `manuscript/fig4/fig4_panelB.pdf`
    The rebucketing moment: the same per-CpG methylation values for one
    kid (Kid1) shown first under hiphase's arbitrary hap1/hap2 labels,
    then relabelled to founder-aware (pat/mat + founder letter) labels.
    The two sub-tables share their numerical content; only the column
    names change — this is the literal "relabel, don't recompute"
    message of the wiki section.
    Wiki source:
      `wiki/pedigree_wise_workflow/founder_phased_methylation/founder_phased_methylation.md`
      section "Relabelling per-CpG methylation".
    Rendered by `wiki/motivation/trio_discovery.py:render_rebucket_panel`.

Panel C — `manuscript/fig4/fig4_panelC.pdf`
    De novo gain-of-methylation discovery table (use case 1). Moved from
    Fig 1D so the reader sees it after the rebucketing machinery has been
    introduced. Trio context (parental cartoon) remains in Fig 1C.
    Wiki source: this is the manuscript-side table, simulated in
    `wiki/motivation/trio_discovery.py` (DENOVO_TABLE_*); the polars-query
    discovery code lives alongside (DENOVO_CODE) and is rendered into
    Supp Fig 1.
    Rendered by `wiki/motivation/trio_discovery.py:render_trio_denovo_meth_table`.
    The row data was extended (10 CpGs, 2 highlighted) so highlighted
    rows read as a discovery rather than as the entire region.

Panel D — `manuscript/fig4/fig4_panelD.pdf`
    Compound genetic-epigenetic heterozygote discovery tables (use case
    2). Moved from Fig 1F. Trio context remains in Fig 1E.
    Wiki source: simulated in `wiki/motivation/trio_discovery.py`
    (COMPOUND_METH_*, COMPOUND_GENO_*); the polars-query discovery code
    is rendered into Supp Fig 1.
    Rendered by `wiki/motivation/trio_discovery.py:render_trio_compound_het_tables`.
    The methylation row data was extended (7 CpGs, 3 highlighted) so the
    highlighted run reads as a discovery within a baseline window.

Run:
    .venv/bin/python manuscript/fig4.py
"""
from __future__ import annotations

import sys
from pathlib import Path

OUT = Path(__file__).resolve().parent / "fig4"
OUT.mkdir(parents=True, exist_ok=True)


def _ensure_wiki_on_path() -> None:
    repo_root = Path(__file__).resolve().parent.parent
    if str(repo_root) not in sys.path:
        sys.path.insert(0, str(repo_root))


def _render_panelA(out_path: Path) -> None:
    """Panel A — bit-vector match cartoon. Reuses the wiki figure's
    rendering at `wiki/.../founder_phased_methylation.py:render_match`
    (single source of truth) but emits PDF for Illustrator."""
    _ensure_wiki_on_path()
    from wiki.pedigree_wise_workflow.founder_phased_methylation import (  # noqa: E402
        founder_phased_methylation as fpm,
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fpm.render_match(out_path=out_path)


def _render_panelB(out_path: Path) -> None:
    """Panel B — visual rebucketing cartoon. Two read-backed phase
    blocks within one IBD segment; the hap1↔hap2 assignment flips
    between the blocks. Bottom shows the same bars re-routed to
    founder-aware tracks (pat_A / mat_C). Reuses
    `wiki/motivation/trio_discovery.py:render_rebucket_visual`."""
    _ensure_wiki_on_path()
    from wiki.motivation import trio_discovery  # noqa: E402
    out_path.parent.mkdir(parents=True, exist_ok=True)
    trio_discovery.render_rebucket_visual(
        out_path=out_path,
        show_title=False,
        trim_whitespace=True,
    )


def _render_panelC(out_path: Path) -> None:
    """Panel C — de novo discovery table (moved from Fig 1D).
    Reuses `wiki/motivation/trio_discovery.py:render_trio_denovo_meth_table`."""
    _ensure_wiki_on_path()
    from wiki.motivation import trio_discovery  # noqa: E402
    out_path.parent.mkdir(parents=True, exist_ok=True)
    trio_discovery.render_trio_denovo_meth_table(
        out_path=out_path,
        show_title=False,
        col_dx=0.105,
        cell_fontsize=9,
        trim_whitespace=True,
    )


def _render_panelD(out_path: Path) -> None:
    """Panel D — compound het discovery tables (moved from Fig 1F).
    Reuses `wiki/motivation/trio_discovery.py:render_trio_compound_het_tables`."""
    _ensure_wiki_on_path()
    from wiki.motivation import trio_discovery  # noqa: E402
    out_path.parent.mkdir(parents=True, exist_ok=True)
    trio_discovery.render_trio_compound_het_tables(
        out_path=out_path,
        show_title=False,
        trim_whitespace=True,
    )


def main() -> None:
    _render_panelA(OUT / "fig4_panelA.pdf")
    _render_panelB(OUT / "fig4_panelB.pdf")
    _render_panelC(OUT / "fig4_panelC.pdf")
    _render_panelD(OUT / "fig4_panelD.pdf")


if __name__ == "__main__":
    main()
