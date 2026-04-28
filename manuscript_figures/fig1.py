"""Render Fig 1 panels for the tapestry manuscript.

Each panel is produced by reusing a render function from the wiki and
overriding only the parameters that need to change for the manuscript
layout (output path, figure size, etc.). The wiki render functions are
the single source of truth; this script is the only place where Fig 1's
manuscript-specific styling lives.

Run:
    .venv/bin/python manuscript_figures/fig1.py

Outputs land in `manuscript_figures/fig1/` as PNG, ready to drop into
the Illustrator composite via File → Place with Link checked. (PNG is
used in the working draft because Illustrator cannot link SVG; switch
to PDF for the final manuscript.)
"""
from __future__ import annotations

import sys
from pathlib import Path

import cairosvg

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "wiki"))

from motivation import single_indiv_phasing, trio_discovery  # noqa: E402

OUT = Path(__file__).resolve().parent / "fig1"
OUT.mkdir(parents=True, exist_ok=True)

CAIROSVG_DPI = 300


def _trio_svg_to_png(scenario, out_path: Path, label_fontsize: int = 16) -> None:
    svg = trio_discovery.build(scenario, label_fontsize=label_fontsize)
    cairosvg.svg2png(
        bytestring=svg.encode("utf-8"),
        write_to=str(out_path),
        dpi=CAIROSVG_DPI,
    )
    print(f"wrote {out_path}")


# Panel A — single-individual, before phasing (pooled bigwig + read pile-up).
def panel_A() -> None:
    single_indiv_phasing.render_before_phasing(
        out_path=OUT / "panel_A_before_phasing.png",
        glyph_fontsize=18,
        ytick_fontsize=18,
        show_left_labels=False,
        show_title=False,
    )


# Panel B — single-individual, after phasing (two stacked bigwigs + grouped reads).
def panel_B() -> None:
    single_indiv_phasing.render_after_phasing(
        out_path=OUT / "panel_B_after_phasing.png",
        glyph_fontsize=18,
        ytick_fontsize=18,
        show_left_labels=False,
        show_title=False,
    )


# Panel C — trio cartoon for de novo epimutation.
def panel_C() -> None:
    _trio_svg_to_png(
        trio_discovery.SCENARIO_DENOVO,
        out_path=OUT / "panel_C_trio_denovo.png",
        label_fontsize=28,
    )


# Panel D — haplotype-specific methylation table (de novo scenario).
def panel_D() -> None:
    trio_discovery.render_trio_denovo_meth_table(
        out_path=OUT / "panel_D_trio_denovo_meth_table.png",
        show_title=False,
        col_dx=0.105,
        cell_fontsize=9,
        trim_whitespace=True,
    )


# Panel E — polars discovery snippet (de novo scenario).
def panel_E() -> None:
    trio_discovery.render_trio_denovo_polars_code(
        out_path=OUT / "panel_E_trio_denovo_polars_code.png",
        show_title=False,
        trim_whitespace=True,
    )


# Panel F — trio cartoon for compound genetic-epigenetic heterozygote.
def panel_F() -> None:
    _trio_svg_to_png(
        trio_discovery.SCENARIO_COMPOUND,
        out_path=OUT / "panel_F_trio_compound_het.png",
        label_fontsize=28,
    )


PANELS = {
    "A": panel_A,
    "B": panel_B,
    "C": panel_C,
    "D": panel_D,
    "E": panel_E,
    "F": panel_F,
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
