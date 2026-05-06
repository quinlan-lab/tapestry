"""Render Fig 3 for the tapestry manuscript.

Each panel mirrors a specific wiki section so the manuscript figure and
its caption can be traced back to a single source of truth. When writing
the figure caption, pull context from the linked sections.

Panel A — `manuscript_figures/fig3/fig3_panelA.pdf`
    The exhaustive-enumeration table at clean site N1.
    Wiki source:
      `wiki/pedigree_wise_workflow/inheritance_mapping/concordance/concordance.md`
      section 3 ("Deducing variant phase by exhaustive enumeration at a
      clean site"), which embeds `concordance/fig2.png`.
    The orientation-enumeration logic is replicated from upstream
    `Platinum-Pedigree-Inheritance/wiki/generate_wiki.py` at the SHA
    pinned in
    `wiki/pedigree_wise_workflow/inheritance_mapping/README.md`
    (`7448e5e946adbc7969ad5fd5e0730d7cace23a8d`), so the row contents
    stay in lockstep with the wiki figure rather than being hand-typed.
    Differences from the wiki rendering: the table is laid out with a
    two-row header so A/B/C/D sit under "Founder assignment" and
    K1/K2/K3 sit under each "Kids (...)" group, the "<- winner"
    annotation is dropped, and the output is a vector PDF for Illustrator
    placement.

Panel B — `manuscript_figures/fig3/fig3_panelB.pdf`
    The bit-vector match cartoon (hap1 vs pat/mat over a hap-map block).
    Wiki source:
      `wiki/pedigree_wise_workflow/founder_phased_methylation/founder_phased_methylation.md`
      section "Bit-vector match", which embeds `bit_vector_match.png`.
    Rendering is delegated to that wiki page's own generator
    (`wiki/.../founder_phased_methylation.py:render_match`), so panel B
    is byte-equivalent to the wiki figure modulo the PDF container.

Panel C — `manuscript_figures/fig3/fig3_panelC.pdf`
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

Panel D — `manuscript_figures/fig3/fig3_panelD.pdf`
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

Panel E — `manuscript_figures/fig3/fig3_panelE.pdf`
    Compound genetic-epigenetic heterozygote discovery tables (use case
    2). Moved from Fig 1F. Trio context remains in Fig 1E.
    Wiki source: simulated in `wiki/motivation/trio_discovery.py`
    (COMPOUND_METH_*, COMPOUND_GENO_*); the polars-query discovery code
    is rendered into Supp Fig 1.
    Rendered by `wiki/motivation/trio_discovery.py:render_trio_compound_het_tables`.
    The methylation row data was extended (7 CpGs, 3 highlighted) so the
    highlighted run reads as a discovery within a baseline window.

Run:
    .venv/bin/python manuscript_figures/fig3.py
"""
from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

OUT = Path(__file__).resolve().parent / "fig3"
OUT.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# Replicated from upstream generate_wiki.py at the pinned SHA.
# ---------------------------------------------------------------------------
# Site N1: clean pass at a non-informative site (dad=0/1, mom=0/1).
N1_TRUTH = {"A": 0, "B": 1, "C": 0, "D": 1}

# Block letter labels, kid -> (paternal letter, maternal letter).
BLOCK_LABELS: Dict[str, Tuple[str, str]] = {
    "Kid1": ("A", "C"),
    "Kid2": ("B", "D"),
    "Kid3": ("A", "C"),
}

# All 2^2 = 4 founder-letter orientations.
ORIENTATIONS: List[Tuple[Tuple[str, str], Tuple[str, str]]] = [
    (("A", "B"), ("C", "D")),
    (("B", "A"), ("C", "D")),
    (("A", "B"), ("D", "C")),
    (("B", "A"), ("D", "C")),
]


def _sorted_unphased(pair: Tuple[int, int]) -> Tuple[int, int]:
    return (min(pair), max(pair))


def _site_observed() -> Dict[str, Tuple[int, int]]:
    """Observed unphased genotypes at site N1 (no injected error)."""
    a, b, c, d = N1_TRUTH["A"], N1_TRUTH["B"], N1_TRUTH["C"], N1_TRUTH["D"]
    truth = {
        "Dad": (a, b),
        "Mom": (c, d),
        "Kid1": (a, c),
        "Kid2": (b, d),
        "Kid3": (a, c),
    }
    return {k: _sorted_unphased(v) for k, v in truth.items()}


def _expected_under_orientation(
    dad_letters: Tuple[str, str],
    mom_letters: Tuple[str, str],
    dad_vcf: Tuple[int, int],
    mom_vcf: Tuple[int, int],
) -> Dict[str, Tuple[int, int]]:
    letter_to_allele: Dict[str, int] = {
        dad_letters[0]: dad_vcf[0],
        dad_letters[1]: dad_vcf[1],
        mom_letters[0]: mom_vcf[0],
        mom_letters[1]: mom_vcf[1],
    }
    expected: Dict[str, Tuple[int, int]] = {
        "Dad": (letter_to_allele[dad_letters[0]], letter_to_allele[dad_letters[1]]),
        "Mom": (letter_to_allele[mom_letters[0]], letter_to_allele[mom_letters[1]]),
    }
    for kid, (lp, lm) in BLOCK_LABELS.items():
        expected[kid] = (letter_to_allele[lp], letter_to_allele[lm])
    return expected


def _count_mismatches(
    observed: Dict[str, Tuple[int, int]],
    expected: Dict[str, Tuple[int, int]],
) -> int:
    return sum(
        1 for k in observed
        if _sorted_unphased(observed[k]) != _sorted_unphased(expected[k])
    )


# ---------------------------------------------------------------------------
# Build the rows shown in the table.
# ---------------------------------------------------------------------------
def _build_rows() -> List[List[str]]:
    observed = _site_observed()
    dad_vcf = observed["Dad"]
    mom_vcf = observed["Mom"]

    rows: List[Tuple[Tuple[int, int, int, int], List[str]]] = []
    for dad_letters, mom_letters in ORIENTATIONS:
        expected = _expected_under_orientation(dad_letters, mom_letters, dad_vcf, mom_vcf)
        letter_to_allele = {
            dad_letters[0]: dad_vcf[0], dad_letters[1]: dad_vcf[1],
            mom_letters[0]: mom_vcf[0], mom_letters[1]: mom_vcf[1],
        }
        a, b, c, d = (letter_to_allele[L] for L in "ABCD")
        n_mis = _count_mismatches(observed, expected)
        cells: List[str] = [str(a), str(b), str(c), str(d)]
        for k in ("Kid1", "Kid2", "Kid3"):
            cells.append(f"{expected[k][0]} {expected[k][1]}")
        for k in ("Kid1", "Kid2", "Kid3"):
            s = _sorted_unphased(expected[k])
            cells.append(f"{s[0]} {s[1]}")
        for k in ("Kid1", "Kid2", "Kid3"):
            s = observed[k]
            cells.append(f"{s[0]} {s[1]}")
        cells.append(str(n_mis))
        rows.append(((a, b, c, d), cells))

    rows.sort(key=lambda r: r[0])
    return [cells for _, cells in rows]


# ---------------------------------------------------------------------------
# Rendering — two-row header with grouped sub-columns.
# ---------------------------------------------------------------------------
GROUPS: List[Tuple[str, List[str]]] = [
    ("Hap-allele assignment in parents", ["A", "B", "C", "D"]),
    ("Kids (expected, phased)",   ["K1", "K2", "K3"]),
    ("Kids (expected, unphased)", ["K1", "K2", "K3"]),
    ("Kids (observed, unphased)", ["K1", "K2", "K3"]),
    ("# Mis",                      [""]),
]


def _render(out_path: Path, rows: List[List[str]]) -> None:
    sub_cols = [s for _, subs in GROUPS for s in subs]
    n_cols = len(sub_cols)
    n_rows = len(rows)
    header_rows = 2

    col_w = 0.85
    row_h = 0.36
    left_pad = right_pad = 0.2
    top_pad = bot_pad = 0.15

    fig_w = left_pad + right_pad + n_cols * col_w
    fig_h = top_pad + bot_pad + (header_rows + n_rows) * row_h

    fig = plt.figure(figsize=(fig_w, fig_h))
    ax = fig.add_axes([
        left_pad / fig_w, bot_pad / fig_h,
        (n_cols * col_w) / fig_w, ((header_rows + n_rows) * row_h) / fig_h,
    ])
    ax.set_xlim(0, n_cols)
    ax.set_ylim(0, header_rows + n_rows)
    ax.set_axis_off()

    font_size = 14
    group_font_size = 11
    body_font_size = 15

    def y_for_row(idx_from_top: int) -> float:
        return header_rows + n_rows - idx_from_top - 1

    # Group label row.
    cursor = 0
    for label, subs in GROUPS:
        span = len(subs)
        y0 = y_for_row(0)
        ax.add_patch(Rectangle(
            (cursor, y0), span, 1,
            facecolor="white", edgecolor="none",
        ))
        ax.text(
            cursor + span / 2, y0 + 0.5, label,
            ha="center", va="center", fontsize=group_font_size, family="monospace",
        )
        cursor += span

    # Sub-column row.
    for j, sub in enumerate(sub_cols):
        y0 = y_for_row(1)
        ax.add_patch(Rectangle(
            (j, y0), 1, 1,
            facecolor="#f4f4f4", edgecolor="none",
        ))
        if sub:
            ax.text(
                j + 0.5, y0 + 0.5, sub,
                ha="center", va="center", fontsize=font_size, family="monospace",
            )

    # Body rows.
    for i, row in enumerate(rows):
        for j, cell in enumerate(row):
            y0 = y_for_row(header_rows + i)
            ax.add_patch(Rectangle(
                (j, y0), 1, 1,
                facecolor="white", edgecolor="none",
            ))
            ax.text(
                j + 0.5, y0 + 0.5, cell,
                ha="center", va="center", fontsize=body_font_size, family="monospace",
            )

    # Inner separators only (no outer table border).
    line_kw = dict(color="black", linewidth=0.5)
    # Horizontal: only under the sub-column header row.
    ax.plot([0, n_cols], [n_rows, n_rows], **line_kw)
    # Vertical separators only at group boundaries (full height).
    cursor = 0
    boundaries = []
    for _, subs in GROUPS:
        cursor += len(subs)
        boundaries.append(cursor)
    for j in boundaries[:-1]:  # skip the rightmost (outer edge)
        ax.plot([j, j], [0, header_rows + n_rows], **line_kw)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)
    print(f"wrote {out_path}")


def _ensure_wiki_on_path() -> None:
    import sys
    repo_root = Path(__file__).resolve().parent.parent
    if str(repo_root) not in sys.path:
        sys.path.insert(0, str(repo_root))


def _render_panelB(out_path: Path) -> None:
    """Panel B — bit-vector match cartoon. Reuses the wiki figure's
    rendering at `wiki/.../founder_phased_methylation.py:render_match`
    (single source of truth) but emits PDF for Illustrator."""
    _ensure_wiki_on_path()
    from wiki.pedigree_wise_workflow.founder_phased_methylation import (  # noqa: E402
        founder_phased_methylation as fpm,
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fpm.render_match(out_path=out_path)


def _render_panelC(out_path: Path) -> None:
    """Panel C — rebucketing cartoon. Reuses
    `wiki/motivation/trio_discovery.py:render_rebucket_panel`."""
    _ensure_wiki_on_path()
    from wiki.motivation import trio_discovery  # noqa: E402
    out_path.parent.mkdir(parents=True, exist_ok=True)
    trio_discovery.render_rebucket_panel(
        out_path=out_path,
        show_title=False,
        cell_fontsize=9,
        trim_whitespace=True,
    )


def _render_panelD(out_path: Path) -> None:
    """Panel D — de novo discovery table (moved from Fig 1D).
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


def _render_panelE(out_path: Path) -> None:
    """Panel E — compound het discovery tables (moved from Fig 1F).
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
    rows = _build_rows()
    _render(OUT / "fig3_panelA.pdf", rows)
    _render_panelB(OUT / "fig3_panelB.pdf")
    _render_panelC(OUT / "fig3_panelC.pdf")
    _render_panelD(OUT / "fig3_panelD.pdf")
    _render_panelE(OUT / "fig3_panelE.pdf")


if __name__ == "__main__":
    main()
