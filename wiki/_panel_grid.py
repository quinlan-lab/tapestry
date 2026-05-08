"""Shared panel-grid renderer for the manuscript.

Single source of truth for the row-tuple model used by `manuscript/fig2.py`
and the bit-vector match panel in
`wiki/.../founder_phased_methylation.py`.

A panel is a list of `Row` tuples drawn into one matplotlib axes in
0..1 coordinates. Rows occupy equal vertical pitch.

Row tuple layout (positional, trailing fields optional):
    (label, cells, colors, merge_runs)
    (..., bottom_colors)
    (..., bottom_colors, text_colors)
    (..., bottom_colors, text_colors, label_fontsize)

- label:           row label drawn on the left.
- cells:           per-column string content.
- colors:          per-column palette key or None for fill.
- merge_runs:      collapse adjacent equal-color cells into one
                   rectangle with the run-start string as a single
                   centered label.
- bottom_colors:   per-column palette key for the bottom half of a
                   merged-run rectangle (top stripe = colors[j], bottom
                   stripe = bottom_colors[j]).
- text_colors:     per-column hex color overriding the default black
                   for cell text. Used by tick rows.
- label_fontsize:  font size for the *centered* merged-run label,
                   when the long-form label needs to fit inside a
                   small rectangle.
"""
from __future__ import annotations

from pathlib import Path
from typing import List, Tuple

import matplotlib.pyplot as plt

HAP_PALETTE = {
    "A": "#a6cee3",  # light blue   — dad hap1
    "B": "#fdbf6f",  # light orange — dad hap2
    "C": "#b2df8a",  # light green  — mom hap1
    "D": "#fb9a99",  # light pink   — mom hap2
    "neutral_top": "#c8c8c8",
    "neutral_bot": "#969696",
}

Row = Tuple
SPACER: Row = ("", [], [], False)


def plain(label: str, cells: List[str]) -> Row:
    return (label, cells, [None] * len(cells), False)


def _field(row, idx, default=None):
    return row[idx] if len(row) > idx else default


def draw_two_stripe_block(ax, x0, x1, y_center, height,
                          top_key, bot_key,
                          label=None, label_fontsize=11):
    """Two-stripe rectangle (top half = top_key, bottom half = bot_key)
    used by `draw_panel`'s merge_runs+bottom_colors path AND by
    standalone callers that want the same visual idiom in their own
    axes (e.g. Fig 4B's phase-block / hap-map-block rows over a bar
    chart). Keeping the primitive in one place ensures Fig 4A and
    Fig 4B render their blocks identically."""
    width = x1 - x0
    ax.add_patch(plt.Rectangle((x0, y_center), width, height / 2,
                               facecolor=HAP_PALETTE[top_key],
                               edgecolor="none", zorder=1))
    ax.add_patch(plt.Rectangle((x0, y_center - height / 2), width, height / 2,
                               facecolor=HAP_PALETTE[bot_key],
                               edgecolor="none", zorder=1))
    if label:
        ax.text((x0 + x1) / 2, y_center, label,
                ha="center", va="center", fontsize=label_fontsize,
                family="monospace", color="black", zorder=2)


def draw_block(ax, x0, x1, y_center, height, color_key,
               label=None, label_fontsize=11):
    """Single-color rectangle counterpart to `draw_two_stripe_block`."""
    width = x1 - x0
    ax.add_patch(plt.Rectangle((x0, y_center - height / 2), width, height,
                               facecolor=HAP_PALETTE[color_key],
                               edgecolor="none", zorder=1))
    if label:
        ax.text((x0 + x1) / 2, y_center, label,
                ha="center", va="center", fontsize=label_fontsize,
                family="monospace", color="black", zorder=2)


def draw_panel(
    ax,
    rows: List[Row],
    pitch: float,
    cells_x0: float,
    font_size: int = 11,
    label_x: float = 0.02,
) -> None:
    n = max(len(rows), 1)
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)

    line_h = 1.0 / n
    cell_w = pitch
    cell_h = line_h * 0.85

    for i, row in enumerate(rows):
        label, cells, colors, merge_runs = row[:4]
        bottom_colors = _field(row, 4)
        text_colors = _field(row, 5)
        label_fontsize = _field(row, 6, font_size)
        y = 1.0 - (i + 0.5) * line_h
        if label:
            ax.text(
                label_x, y, label, ha="left", va="center",
                fontsize=font_size, family="monospace",
            )
        if merge_runs:
            j = 0
            while j < len(cells):
                k = j
                while k + 1 < len(cells) and colors[k + 1] == colors[j] and (
                    bottom_colors is None
                    or bottom_colors[k + 1] == bottom_colors[j]
                ):
                    k += 1
                run_color = colors[j]
                if run_color is not None:
                    x_left = cells_x0 + j * pitch - cell_w / 2
                    x_right = x_left + (k - j) * pitch + cell_w
                    if bottom_colors is None:
                        draw_block(ax, x_left, x_right, y, cell_h,
                                   run_color, label=cells[j],
                                   label_fontsize=label_fontsize)
                    else:
                        # Inset on the boundary-facing side(s) when the
                        # run abuts another run, so the gap straddles
                        # the boundary.
                        half_gap = pitch * 0.10
                        gap_left = half_gap if j > 0 else 0.0
                        gap_right = half_gap if k + 1 < len(cells) else 0.0
                        draw_two_stripe_block(
                            ax, x_left + gap_left, x_right - gap_right, y,
                            cell_h, run_color, bottom_colors[j],
                            label=cells[j], label_fontsize=label_fontsize,
                        )
                j = k + 1
            continue
        for j, (c, color_key) in enumerate(zip(cells, colors)):
            x = cells_x0 + j * pitch
            if color_key is not None:
                ax.add_patch(plt.Rectangle(
                    (x - cell_w / 2, y - cell_h / 2),
                    cell_w, cell_h,
                    facecolor=HAP_PALETTE[color_key],
                    edgecolor="none", zorder=1,
                ))
            tcol = text_colors[j] if text_colors is not None and text_colors[j] else "black"
            ax.text(
                x, y, c, ha="center", va="center",
                fontsize=font_size, family="monospace",
                color=tcol, zorder=2,
            )


def row_y(row_index: int, n_rows: int) -> float:
    """Vertical position (axis coords) of a row inside `draw_panel`'s
    layout — used by callers that want to overlay artists on a row."""
    return 1.0 - (row_index + 0.5) / n_rows


def render_combined_pdf(
    panels: List[List[Row]],
    out_path: Path,
    fig_w: float = 5.5,
    pitch: float = 0.105,
    cells_x0: float = 0.22,
    row_height_in: float = 0.30,
    panel_gap_rows: float = 1.5,
) -> None:
    """Stack multiple panels in one figure with shared x layout (so
    columns align across panels)."""
    counts = [max(len(p), 1) for p in panels]
    total_rows = sum(counts)
    fig_h = (
        row_height_in * total_rows
        + row_height_in * panel_gap_rows * (len(panels) - 1)
        + 0.3
    )
    fig, axes = plt.subplots(
        nrows=len(panels), ncols=1,
        figsize=(fig_w, fig_h),
        gridspec_kw={
            "height_ratios": counts,
            "hspace": panel_gap_rows / (total_rows / len(panels)),
        },
    )
    if len(panels) == 1:
        axes = [axes]
    for ax, rows in zip(axes, panels):
        draw_panel(ax, rows, pitch=pitch, cells_x0=cells_x0)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, bbox_inches="tight", pad_inches=0.1)
    plt.close(fig)
    print(f"wrote {out_path}")
