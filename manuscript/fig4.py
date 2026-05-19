"""Render Fig 4 for the tapestry manuscript.

Fig 4 collects the founder-phased-methylation machinery and the trio
discovery tables. These panels were previously B–E of Fig 3 and have
been relettered A–D for Fig 4. The exhaustive-enumeration table at
clean site N1 remains in Fig 3 (panel A there).

Panels A + B — `manuscript/fig4/fig4_panelAB.pdf`
    Combined into one PDF so row labels and rectangle left edges are
    vertically registered across both panels (no Illustrator-side
    margin / font tweaks needed). Panel A is the bit-vector match
    cartoon (hap1 vs pat/mat over a hap-map block); Panel B is the
    rebucketing cartoon — two read-backed phase blocks inside one IBD
    segment, the hap1↔hap2 assignment flipping between blocks, with
    bars re-routed to founder-aware tracks below.
    Source data for Panel A's bit vectors and ranges is imported from
    the wiki cartoon's module
    (`wiki/.../founder_phased_methylation.py`), so this PDF and the
    wiki PNG share one source of truth.
    Rendered by `wiki/motivation/trio_discovery.py:render_panel_AB_combined`.

Panel C — `manuscript/fig4/fig4_panelC.pdf`
    Homolog-specific LOM discovery table (use case 1): the canonical
    maternal-LOM case at an imprinted locus, where the kid's methylation
    on the maternally transmitted founder (kid_mat) drops sharply below
    mom's methylation on the same founder (mom_C). Trio context
    (parental cartoon) remains in Fig 1C.
    Wiki source: this is the manuscript-side table, simulated in
    `wiki/motivation/trio_discovery.py` (LOM_TABLE_*); the polars-query
    discovery code lives alongside (LOM_CODE) and is rendered into
    Supp Fig 1.
    Rendered by `wiki/motivation/trio_discovery.py:render_trio_lom_meth_table`.

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


def _render_panelC(out_path: Path) -> None:
    """Panel C — homolog-specific LOM discovery table (Fig 1C scenario).
    Reuses `wiki/motivation/trio_discovery.py:render_trio_lom_meth_table`."""
    _ensure_wiki_on_path()
    from wiki.motivation import trio_discovery  # noqa: E402
    out_path.parent.mkdir(parents=True, exist_ok=True)
    trio_discovery.render_trio_lom_meth_table(
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


def _render_panelAB(out_path: Path) -> None:
    """Combined Panels A+B as one PDF with row labels and rectangle
    left edges vertically registered across both panels — so the
    manuscript figure can drop in one PDF instead of stitching A and B
    separately and tweaking font sizes / margins in Illustrator."""
    _ensure_wiki_on_path()
    from wiki.motivation import trio_discovery  # noqa: E402
    out_path.parent.mkdir(parents=True, exist_ok=True)
    trio_discovery.render_panel_AB_combined(
        out_path=out_path,
        trim_whitespace=True,
    )


def main() -> None:
    _render_panelAB(OUT / "fig4_panelAB.pdf")
    _render_panelC(OUT / "fig4_panelC.pdf")
    _render_panelD(OUT / "fig4_panelD.pdf")


if __name__ == "__main__":
    main()
