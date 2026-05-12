"""Render Fig 5 panels for the tapestry manuscript.

Panel A is a stylized cartoon: a stack of 8 founder haplotypes from a
3-generation pedigree. Each haplotype is drawn as in Fig 2A — a
Fig 2A-style colored cell containing the 0/1 allele at the meQTL — and
the methylation level at flanking CpGs is shown as a histogram in the
same style as Fig 4 Panel A/B (the bars rows). REF haps (0) at the
meQTL are unmethylated; ALT haps (1) are methylated, illustrating the
expected covariation that panels B/C/D test in pedigree K1463.

Reuses:
  - `wiki/_panel_grid.HAP_PALETTE` and the cell-rendering idiom from
    `manuscript/fig2.py` (founder rows: 0/1 cells colored by hap key).
  - The bars-row style from
    `wiki/motivation/trio_discovery.render_panel_AB_combined`.

Run:
    .venv/bin/python manuscript/fig5.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt

_REPO_ROOT = Path(__file__).resolve().parent.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from wiki._panel_grid import HAP_PALETTE  # noqa: E402

mpl.rcParams["font.family"] = "Arial"

OUT = Path(__file__).resolve().parent / "fig5"
OUT.mkdir(parents=True, exist_ok=True)


# One palette key per founder haplotype. 8 founder haps = 4 grandparents
# x 2. We reuse the Fig 2 palette by repeating the four keys so the
# colors match the rest of the manuscript.
HAP_KEYS = ["A", "B", "C", "D"]
HAP_LABELS = ["A", "B", "C", "D"]
# Allele at the meQTL on each hap (0 = REF, 1 = ALT). Mixed so the
# covariation with methylation is visible at a glance.
HAP_ALLELES = [0, 1, 0, 1]

# CpG column positions along the bar (data x); meQTL sits in the middle.
N_CPG = 9
MEQTL_COL = 4


def _meth_profile(allele: int) -> list[float]:
    """Methylation level (0..1) per CpG column. REF -> low, ALT -> high,
    with mild per-CpG jitter so the histogram doesn't look identical
    across haps in the same group."""
    if allele == 0:
        return [0.10, 0.05, 0.15, 0.08, 0.0, 0.05, 0.12, 0.07, 0.10]
    return [0.88, 0.92, 0.85, 0.95, 0.0, 0.90, 0.83, 0.92, 0.88]


def panel_A() -> None:
    PREFIX_X = -2.0  # right edge of the row-label column (data x)
    XLIM = (-4.2, N_CPG - 0.5 + 0.3)
    LABEL_FS = 12
    BIT_FS = 12

    # Row heights (inches). Each hap = one bars row (methylation
    # histogram, Fig 4 Panel A/B style) stacked above one hap row
    # (Fig 2A founder-row style: a colored bar spanning all columns
    # with the 0/1 allele shown at the meQTL column).
    BAR_H_IN = 0.70
    HAP_H_IN = 0.22
    SPACER_TIGHT = 0.02
    SPACER_BETWEEN_HAPS = 0.18

    rows = []
    for i, (key, lab, allele) in enumerate(zip(HAP_KEYS, HAP_LABELS, HAP_ALLELES)):
        rows.append(("bars", BAR_H_IN, {
            "color_key": key,
            "bars": _meth_profile(allele),
            "label": "Meth",
        }))
        rows.append(("spacer", SPACER_TIGHT, None))
        rows.append(("hap", HAP_H_IN, {
            "color_key": key,
            "allele": allele,
            "label": lab,
        }))
        if i < len(HAP_LABELS) - 1:
            rows.append(("spacer", SPACER_BETWEEN_HAPS, None))


    heights = [r[1] for r in rows]
    fig_h = sum(heights) + 0.3
    fig_w = 6.0
    fig = plt.figure(figsize=(fig_w, fig_h))
    LEFT = 0.10
    RIGHT = 0.98
    TOP = fig_h - 0.05
    BOTTOM = 0.05  # noqa: F841

    cur_top_in = TOP
    for kind, h_in, payload in rows:
        bottom_in = cur_top_in - h_in
        ax = fig.add_axes([
            LEFT, bottom_in / fig_h,
            RIGHT - LEFT, h_in / fig_h,
        ])
        ax.set_xlim(*XLIM)
        ax.set_ylim(0, 1)
        ax.set_axis_off()

        if kind == "spacer":
            cur_top_in = bottom_in
            continue

        if kind == "footer":
            ax.text(MEQTL_COL, 0.55, "meQTL",
                    ha="center", va="top",
                    fontsize=14, fontstyle="italic")
            ax.annotate(
                "", xy=(MEQTL_COL, 1.0), xytext=(MEQTL_COL, 0.62),
                arrowprops=dict(arrowstyle="-", lw=1.0, color="black"),
            )
            cur_top_in = bottom_in
            continue

        if kind == "bars":
            # Methylation histogram (Fig 4 Panel A/B bars-row style).
            ax.set_axis_on()
            for spine in ("top", "right", "bottom"):
                ax.spines[spine].set_visible(False)
            ax.spines["left"].set_visible(True)
            ax.spines["left"].set_position(("data", -0.8))
            ax.spines["left"].set_bounds(0, 1)
            ax.set_xticks([])
            ax.set_yticks([0, 1])
            ax.set_yticklabels(["0", "1"], fontsize=12)
            ax.tick_params(axis="y", length=2, labelsize=12)
            ax.text(PREFIX_X, 0.5, payload["label"],
                    ha="right", va="center",
                    fontsize=LABEL_FS, family="monospace")
            bar_color = HAP_PALETTE[payload["color_key"]]
            for i, lev in enumerate(payload["bars"]):
                if lev > 0:
                    ax.bar(i, lev, width=0.7,
                           color=bar_color, edgecolor="none")
            ax.set_ylim(0, 1.05)
            cur_top_in = bottom_in
            continue

        # kind == "hap": Fig 2A-style founder row — a single colored bar
        # spanning every column, with the 0/1 allele drawn only at the
        # meQTL column (no other digits).
        ax.set_axis_off()
        x0 = -0.5
        x1 = N_CPG - 0.5
        cell_h_axes = 0.85
        ax.add_patch(plt.Rectangle(
            (x0, 0.5 - cell_h_axes / 2),
            x1 - x0, cell_h_axes,
            facecolor=HAP_PALETTE[payload["color_key"]],
            edgecolor="none", zorder=1,
        ))
        ax.text(PREFIX_X, 0.5, payload["label"],
                ha="right", va="center",
                fontsize=LABEL_FS, family="monospace")
        ax.text(MEQTL_COL, 0.5, str(payload["allele"]),
                ha="center", va="center",
                fontsize=BIT_FS, family="monospace",
                color="black", zorder=3)
        cur_top_in = bottom_in

    out_path = OUT / "panel_A_meqtl_cartoon.pdf"
    fig.savefig(out_path, bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)
    print(f"wrote {out_path}")


PANELS = {
    "A": panel_A,
}


def main() -> None:
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--panel", choices=sorted(PANELS.keys()),
        help="Render only one panel (default: render every panel).",
    )
    args = parser.parse_args()
    targets = [PANELS[args.panel]] if args.panel else list(PANELS.values())
    for fn in targets:
        fn()


if __name__ == "__main__":
    main()
