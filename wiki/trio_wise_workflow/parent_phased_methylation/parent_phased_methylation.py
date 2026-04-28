"""Parent-phased methylation: bit-vector matching cartoon (trio).

Step 3 of tapestry's trio-wise workflow
(`phase_meth_to_parent_haps.py` + `hap_map_trio.py`). One panel showing
how the kid's paternal-allele bit vector is matched against dad's
hap1/hap2 (= A/B) bit vectors over the hap-map block (intersection of
the kid's whatshap phase block and dad's whatshap phase block). The
maternal side is symmetric: replace `kid_pat`/`dad_A`/`dad_B` with
`kid_mat`/`mom_C`/`mom_D` and intersect the kid's phase block with
mom's.

Run:
    .venv/bin/python wiki/trio_wise_workflow/parent_phased_methylation/parent_phased_methylation.py

Render kwargs mirror the pedigree-wise sibling
(`wiki/pedigree_wise_workflow/founder_phased_methylation/founder_phased_methylation.py`)
so a `manuscript_figures/*.py` driver can reuse this function with
manuscript-specific styling overrides.

Outputs:
    bit_vector_match_trio.png
"""
from __future__ import annotations
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rcParams["font.family"] = "Arial"

OUT = Path(__file__).resolve().parent

# ---------------------------------------------------------------------------
# Toy data — paternal-side match. Kid's phase block spans sites 0..6 (closed),
# dad's phase block spans sites 2..9; their intersection (the hap-map block)
# is sites 2..6. Inside the intersection, kid_pat agrees with dad_A at
# 4 of 5 sites → concordance(kid_pat, dad_A) = 0.80, so kid's paternal
# haplotype is labelled `A`.
# ---------------------------------------------------------------------------

N_SITES = 10
KID_RANGE = (0, 6)
DAD_RANGE = (2, 9)
INT_RANGE = (max(KID_RANGE[0], DAD_RANGE[0]),
             min(KID_RANGE[1], DAD_RANGE[1]))

KID_PAT = [1, 0, 0, 1, 0, 1, 1, ".", ".", "."]
DAD_A   = [".", ".", 0, 1, 0, 1, 0, 1, 0, 1]
DAD_B   = [".", ".", 1, 0, 1, 0, 1, 0, 1, 0]

# ---------------------------------------------------------------------------
# Render
# ---------------------------------------------------------------------------

LABEL_X_DEFAULT = -0.1
CONCORD_X_DEFAULT = -2.4
SITE_DX_DEFAULT = 0.22
BAR_COLOR = "#222222"
TEXT_COLOR = "#222222"
GLYPH_FONTSIZE = 12


def _similarity(a, b, lo, hi):
    n_match = sum(1 for i in range(lo, hi + 1) if a[i] == b[i])
    return n_match, hi - lo + 1


def _draw_bar(ax, y, lo, hi, site_dx, drop=0.18):
    x0, x1 = lo * site_dx, hi * site_dx
    ax.plot([x0, x1], [y, y], color=BAR_COLOR, linewidth=1.2)
    ax.plot([x0, x0], [y - drop, y + drop], color=BAR_COLOR, linewidth=1.2)
    ax.plot([x1, x1], [y - drop, y + drop], color=BAR_COLOR, linewidth=1.2)


def _label(ax, y, text, label_x, fontsize):
    ax.text(label_x, y, text, ha="right", va="center",
            fontsize=fontsize, color=TEXT_COLOR)


def _bit_row(ax, y, items, site_dx, fontsize):
    for i, v in enumerate(items):
        ax.text(i * site_dx, y, str(v), ha="center", va="center",
                fontsize=fontsize, color=TEXT_COLOR)


def render_match(
    out_path: Path | None = None,
    fig_width: float = 6.4,
    fig_height: float = 4.6,
    glyph_fontsize: int = GLYPH_FONTSIZE,
    label_x: float = LABEL_X_DEFAULT,
    concord_x: float = CONCORD_X_DEFAULT,
    site_dx: float = SITE_DX_DEFAULT,
) -> Path:
    n_A, n_total = _similarity(KID_PAT, DAD_A, *INT_RANGE)
    n_B, _ = _similarity(KID_PAT, DAD_B, *INT_RANGE)

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    fig.subplots_adjust(left=0.02, right=0.98, top=0.96, bottom=0.04)
    ax.set_axis_off()

    # Y layout (top → bottom).
    y_kid = 6.0
    y_dad = 5.2
    y_int = 4.4
    y_kid_pat = 3.0
    y_A = 2.2
    y_B = 1.4
    y_concord_A = -0.7
    y_concord_B = -1.4
    y_concord_note = -2.4
    y_concord_decision = -3.3

    ax.set_xlim(concord_x - 0.1, (N_SITES - 1) * site_dx + 0.3)
    ax.set_ylim(y_concord_decision - 0.6, y_kid + 0.6)

    # Block-extent bars.
    _label(ax, y_kid, "kid phase block:", label_x, glyph_fontsize)
    _draw_bar(ax, y_kid, *KID_RANGE, site_dx)
    _label(ax, y_dad, "dad phase block:", label_x, glyph_fontsize)
    _draw_bar(ax, y_dad, *DAD_RANGE, site_dx)
    _label(ax, y_int, "hap-map block (intersection):", label_x, glyph_fontsize)
    _draw_bar(ax, y_int, *INT_RANGE, site_dx)

    # Bit rows.
    _label(ax, y_kid_pat, "kid_pat:", label_x, glyph_fontsize)
    _bit_row(ax, y_kid_pat, KID_PAT, site_dx, glyph_fontsize)
    _label(ax, y_A, "dad_A:", label_x, glyph_fontsize)
    _bit_row(ax, y_A, DAD_A, site_dx, glyph_fontsize)
    _label(ax, y_B, "dad_B:", label_x, glyph_fontsize)
    _bit_row(ax, y_B, DAD_B, site_dx, glyph_fontsize)

    # Haplotype concordance.
    ax.text(concord_x, y_concord_A,
            f"concordance(kid_pat, dad_A)  =  {n_A}/{n_total}  =  {n_A / n_total:.2f}",
            ha="left", va="center", fontsize=glyph_fontsize, color=TEXT_COLOR)
    ax.text(concord_x, y_concord_B,
            f"concordance(kid_pat, dad_B)  =  {n_B}/{n_total}  =  {n_B / n_total:.2f}",
            ha="left", va="center", fontsize=glyph_fontsize, color=TEXT_COLOR)
    ax.text(concord_x, y_concord_note,
            "(sites are heterozygous in dad, so the two concordances sum to 1)",
            ha="left", va="center", fontsize=glyph_fontsize - 1, color=TEXT_COLOR,
            style="italic")
    ax.text(concord_x, y_concord_decision,
            "concordance(kid_pat, dad_A) > 0.5  →  kid_pat = A",
            ha="left", va="center", fontsize=glyph_fontsize, color=TEXT_COLOR,
            weight="bold")

    out_path = out_path or (OUT / "bit_vector_match_trio.png")
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
    print(f"[wiki] Wrote {out_path}")
    return out_path


def main() -> None:
    render_match()


if __name__ == "__main__":
    main()
