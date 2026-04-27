"""Founder-phased methylation: bit-vector matching cartoon.

Step 3 of tapestry's pedigree-wise workflow
(`phase_meth_to_founder_haps.py` + `hap_map_pedigree.py`). One panel
showing how a hap1 bit-vector from hiphase is matched against the
pat/mat bit-vectors from gtg-concordance over the hap-map block
(intersection of the read-backed phase block and the linkage block).

Run:
    .venv/bin/python wiki/pedigree_wise_workflow/founder_phased_methylation/founder_phased_methylation.py

Outputs:
    bit_vector_match.png
"""
from __future__ import annotations
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rcParams["font.family"] = "Arial"

OUT = Path(__file__).resolve().parent

# ---------------------------------------------------------------------------
# Toy data
# ---------------------------------------------------------------------------

N_SITES = 10
RBP_RANGE = (0, 6)
LNK_RANGE = (2, 9)
INT_RANGE = (max(RBP_RANGE[0], LNK_RANGE[0]),
             min(RBP_RANGE[1], LNK_RANGE[1]))

HAP1 = [1, 0, 0, 1, 0, 1, 1, ".", ".", "."]
PAT  = [".", ".", 0, 1, 0, 1, 0, 1, 0, 1]
MAT  = [".", ".", 1, 0, 1, 0, 1, 0, 1, 0]

# ---------------------------------------------------------------------------
# Render
# ---------------------------------------------------------------------------

LABEL_X = -0.1     # right edge of label column, in site-x units
CONCORD_X = -2.0   # left edge of the haplotype-concordance block
SITE_DX = 0.22
BAR_COLOR = "#222222"
TEXT_COLOR = "#222222"
FS = 12


def _similarity(a, b, lo, hi):
    n_match = sum(1 for i in range(lo, hi + 1) if a[i] == b[i])
    return n_match, hi - lo + 1


def _draw_bar(ax, y, lo, hi, drop=0.18):
    x0, x1 = lo * SITE_DX, hi * SITE_DX
    ax.plot([x0, x1], [y, y], color=BAR_COLOR, linewidth=1.2)
    ax.plot([x0, x0], [y - drop, y + drop], color=BAR_COLOR, linewidth=1.2)
    ax.plot([x1, x1], [y - drop, y + drop], color=BAR_COLOR, linewidth=1.2)


def _label(ax, y, text):
    ax.text(LABEL_X, y, text, ha="right", va="center",
            fontsize=FS, color=TEXT_COLOR)


def _bit_row(ax, y, items):
    for i, v in enumerate(items):
        ax.text(i * SITE_DX, y, str(v), ha="center", va="center",
                fontsize=FS, color=TEXT_COLOR)


def render_match(out_path: Path | None = None) -> Path:
    n_pat, n_total = _similarity(HAP1, PAT, *INT_RANGE)
    n_mat, _ = _similarity(HAP1, MAT, *INT_RANGE)

    fig, ax = plt.subplots(figsize=(6.0, 4.6))
    fig.subplots_adjust(left=0.02, right=0.98, top=0.96, bottom=0.04)
    ax.set_axis_off()

    # Y layout (top → bottom).
    y_rbp = 6.0
    y_lnk = 5.2
    y_int = 4.4
    y_hap1 = 3.0
    y_pat = 2.2
    y_mat = 1.4
    y_concord_title = 0.0
    y_concord_pat = -0.7
    y_concord_mat = -1.4
    y_concord_note = -2.4
    y_concord_decision = -3.3

    ax.set_xlim(CONCORD_X - 0.1, (N_SITES - 1) * SITE_DX + 0.3)
    ax.set_ylim(y_concord_decision - 0.6, y_rbp + 0.6)

    # Block-extent bars.
    _label(ax, y_rbp, "read-backed phase block:")
    _draw_bar(ax, y_rbp, *RBP_RANGE)
    _label(ax, y_lnk, "linkage block:")
    _draw_bar(ax, y_lnk, *LNK_RANGE)
    _label(ax, y_int, "hap-map block (intersection):")
    _draw_bar(ax, y_int, *INT_RANGE)

    # Bit rows.
    _label(ax, y_hap1, "hap1:")
    _bit_row(ax, y_hap1, HAP1)
    _label(ax, y_pat, "pat:")
    _bit_row(ax, y_pat, PAT)
    _label(ax, y_mat, "mat:")
    _bit_row(ax, y_mat, MAT)

    # Haplotype concordance.
    ax.text(CONCORD_X, y_concord_pat,
            f"concordance(hap1, pat)  =  {n_pat}/{n_total}  =  {n_pat / n_total:.2f}",
            ha="left", va="center", fontsize=FS, color=TEXT_COLOR)
    ax.text(CONCORD_X, y_concord_mat,
            f"concordance(hap1, mat)  =  {n_mat}/{n_total}  =  {n_mat / n_total:.2f}",
            ha="left", va="center", fontsize=FS, color=TEXT_COLOR)
    ax.text(CONCORD_X, y_concord_note,
            "(sites are heterozygous, so the two concordances sum to 1)",
            ha="left", va="center", fontsize=FS - 1, color=TEXT_COLOR,
            style="italic")
    ax.text(CONCORD_X, y_concord_decision,
            "concordance(hap1, pat) > 0.5  →  hap1 = pat",
            ha="left", va="center", fontsize=FS, color=TEXT_COLOR,
            weight="bold")

    out_path = out_path or (OUT / "bit_vector_match.png")
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
    print(f"[wiki] Wrote {out_path}")
    return out_path


def main() -> None:
    render_match()


if __name__ == "__main__":
    main()
