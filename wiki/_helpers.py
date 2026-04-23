"""
Shared helpers for tapestry's wiki generator.

Two things go here:

1. Permalink utilities (`head_sha`, `permalink`) — pin GitHub source links
   to the current HEAD SHA at generate time.
2. The `draw_methylation_bars` helper, used by both `motivation/`
   figures and the §13 bit-vector cartoon, so methylation-level bars
   render on a directly comparable y-axis everywhere they appear in
   the wiki.

`text_panel` and friends from the upstream
`Platinum-Pedigree-Inheritance/wiki/generate_wiki.py` are deliberately
not yet ported — the motivation cartoons render as visual pileups +
bar plots, not as monospace text. They'll be ported when the
`founder_phased_methylation/` page (Phase 5) needs them.
"""

from __future__ import annotations

import subprocess
from typing import Sequence

import matplotlib.pyplot as plt

TAPESTRY_REPO = "quinlan-lab/tapestry"


def head_sha() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"], text=True
        ).strip()
    except Exception:
        return "main"


def permalink(path: str, line: int, sha: str, repo: str = TAPESTRY_REPO) -> str:
    return f"https://github.com/{repo}/blob/{sha}/{path}#L{line}"


METH_YLIM = (0, 100)
METH_YTICKS = (0, 25, 50, 75, 100)


def draw_methylation_bars(
    ax: plt.Axes,
    labels: Sequence[str],
    fractions: Sequence[float],
    colors: Sequence[str] | None = None,
    title: str | None = None,
) -> None:
    """
    Draw methylation-level bars on a fixed 0–100 % y-axis.

    Both the motivation page (Figs 1 and 2 Panel B) and the §13 hero
    cartoon (Panel D) call this helper so the bars are byte-comparable
    across the wiki.

    `fractions` are in [0, 1]; converted to percent inside the helper.
    """
    n = len(labels)
    assert len(fractions) == n, "labels and fractions must align"
    if colors is None:
        colors = ["#888888"] * n

    xs = list(range(n))
    heights = [100.0 * f for f in fractions]
    ax.bar(xs, heights, color=list(colors), edgecolor="black", linewidth=0.6)
    ax.set_xticks(xs)
    ax.set_xticklabels(list(labels), fontsize=9)
    ax.set_ylim(*METH_YLIM)
    ax.set_yticks(list(METH_YTICKS))
    ax.set_ylabel("methylation (%)", fontsize=9)
    ax.tick_params(axis="y", labelsize=8)
    if title is not None:
        ax.set_title(title, fontsize=10)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
