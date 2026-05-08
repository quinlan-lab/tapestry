"""Founder-phased methylation: bit-vector matching cartoon.

Step 3 of tapestry's pedigree-wise workflow
(`phase_meth_to_founder_haps.py` + `hap_map_pedigree.py`). One panel
showing how a hap1 bit-vector from hiphase is matched against the
pat/mat bit-vectors from gtg-concordance over the hap-map block
(intersection of the read-backed phase block and an IBD segment).

Rendering goes through the shared `wiki/_panel_grid.py` row-tuple model
so this panel and Fig 2 use the same drawing code (palette, monospace
font, two-stripe rectangles).

Run:
    .venv/bin/python wiki/pedigree_wise_workflow/founder_phased_methylation/founder_phased_methylation.py

Outputs:
    bit_vector_match.png
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.patheffects as pe
import matplotlib.pyplot as plt

OUT = Path(__file__).resolve().parent

_REPO_ROOT = Path(__file__).resolve().parents[3]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from wiki._panel_grid import (  # noqa: E402
    SPACER,
    Row,
    draw_panel,
    row_y,
)

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

# Match Fig 3's tick styling: heavy glyph (✔/✘), bold, ~1.75× body
# font, with a thin stroke path-effect for emphasis.
TICK_COLOR = "#1b8a1b"
CROSS_COLOR = "#c0392b"
TICK_FONT = 19  # body font is 11 → 11 × 1.75 ≈ 19

# ---------------------------------------------------------------------------
# Row builders
# ---------------------------------------------------------------------------


def _similarity(a, b, lo, hi):
    n_match = sum(1 for i in range(lo, hi + 1) if a[i] == b[i])
    return n_match, hi - lo + 1


def _block_row(label, lo, hi, top_key, bot_key, in_label, in_label_fs=None) -> Row:
    cells = [""] * N_SITES
    cells[lo] = in_label
    colors = [None] * N_SITES
    bottom_colors = [None] * N_SITES
    for i in range(lo, hi + 1):
        colors[i] = top_key
        bottom_colors[i] = bot_key
    if in_label_fs is None:
        return (label, cells, colors, True, bottom_colors)
    return (label, cells, colors, True, bottom_colors, None, in_label_fs)


def _bit_row(label, items, color_key, color_lo, color_hi) -> Row:
    cells = [str(v) for v in items]
    colors = [None] * N_SITES
    for i in range(color_lo, color_hi + 1):
        if items[i] != ".":
            colors[i] = color_key
    return (label, cells, colors, False)


def _tick_spec(hap1, other, lo, hi):
    """Return a list of (column_index, glyph, color) for a tick row.
    Drawn as overlays after `draw_panel` so we can use a larger font
    and stroke path-effect (Fig 3's body-tick styling)."""
    out = []
    for i in range(lo, hi + 1):
        if hap1[i] == "." or other[i] == ".":
            continue
        if hap1[i] == other[i]:
            out.append((i, "✔", TICK_COLOR))
        else:
            out.append((i, "✘", CROSS_COLOR))
    return out


# Empty placeholder row that reserves the same vertical pitch as a
# regular row, so a tick overlay sits cleanly between bit rows.
_TICK_ROW: Row = ("", [""] * N_SITES, [None] * N_SITES, False)


# ---------------------------------------------------------------------------
# Render
# ---------------------------------------------------------------------------

PITCH = 0.045
CELLS_X0 = 0.32


def render_match(out_path: Path | None = None) -> Path:
    n_pat, n_total = _similarity(HAP1, PAT, *INT_RANGE)
    n_mat, _ = _similarity(HAP1, MAT, *INT_RANGE)

    rows = [
        _block_row("Read-backed phase block:", *RBP_RANGE,
                   "neutral_top", "neutral_bot", "hap1  /  hap2"),
        _block_row("IBD segment:", *LNK_RANGE,
                   "A", "C", "A  |  C"),
        SPACER,
        _bit_row("hap1:", HAP1, "neutral_top", *RBP_RANGE),
        _bit_row("pat_A:", PAT, "A", *LNK_RANGE),
        _TICK_ROW,
        _bit_row("mat_C:", MAT, "C", *LNK_RANGE),
        _TICK_ROW,
        SPACER,
        _block_row("Deduced hap-map block:", *INT_RANGE,
                   "A", "C", "hap1 → A  |  hap2 → C", in_label_fs=9),
    ]

    # Match Fig 2's row pitch (0.30 in/row → cell_h ≈ 0.255 in at the
    # shared 85 % cell-to-row ratio), so the rectangle:text proportion
    # in Fig 4A reads the same as Fig 2.
    fig, ax = plt.subplots(figsize=(10.0, 0.30 * len(rows)))
    fig.subplots_adjust(left=0.02, right=0.98, top=0.96, bottom=0.04)
    draw_panel(ax, rows, pitch=PITCH, cells_x0=CELLS_X0,
               font_size=11, label_x=0.02)

    # Tick / cross overlays — drawn after `draw_panel` so we can use a
    # larger font and a stroke path-effect (matching Fig 3's body ticks).
    n = len(rows)
    pat_bit_idx = rows.index(_bit_row("pat_A:", PAT, "A", *LNK_RANGE))
    mat_bit_idx = rows.index(_bit_row("mat_C:", MAT, "C", *LNK_RANGE))
    tick_pat_y = row_y(pat_bit_idx + 1, n)
    tick_mat_y = row_y(mat_bit_idx + 1, n)
    for ty, spec in (
        (tick_pat_y, _tick_spec(HAP1, PAT, *INT_RANGE)),
        (tick_mat_y, _tick_spec(HAP1, MAT, *INT_RANGE)),
    ):
        for col, glyph, color in spec:
            t = ax.text(CELLS_X0 + col * PITCH, ty, glyph,
                        ha="center", va="center",
                        fontsize=TICK_FONT, family="monospace",
                        color=color, weight="bold")
            t.set_path_effects([
                pe.Stroke(linewidth=0.8, foreground=color),
                pe.Normal(),
            ])

    # Concordance overlays — placed to the right of the corresponding
    # bit rows, using the panel's row-y formula so they stay aligned
    # if the row order changes.
    concord_x = CELLS_X0 + N_SITES * PITCH + 0.01
    for idx, text in (
        (pat_bit_idx, f"concord(pat_A, hap1) = {n_pat}/{n_total}"),
        (mat_bit_idx, f"concord(mat_C, hap1) = {n_mat}/{n_total}"),
    ):
        ax.text(concord_x, row_y(idx, n), text,
                ha="left", va="center", fontsize=11, family="monospace")

    out_path = out_path or (OUT / "bit_vector_match.png")
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"[wiki] Wrote {out_path}")
    return out_path


def main() -> None:
    render_match()


if __name__ == "__main__":
    main()
