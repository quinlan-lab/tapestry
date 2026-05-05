"""Render Supporting Figure 1 panels for the tapestry manuscript.

Panel A: polars discovery snippet for the de novo scenario.
Panel B: polars discovery snippet for the compound genetic-epigenetic
heterozygote scenario.

Run:
    .venv/bin/python manuscript_figures/supp_fig1.py

Outputs land in `manuscript_figures/supp_fig1/` as PDF, ready to drop
into the Illustrator composite via File → Place with Link checked.
"""
from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "wiki"))

from motivation import trio_discovery  # noqa: E402

OUT = Path(__file__).resolve().parent / "supp_fig1"
OUT.mkdir(parents=True, exist_ok=True)


# Panel A — polars discovery snippet (de novo scenario).
def panel_A() -> None:
    trio_discovery.render_trio_denovo_polars_code(
        out_path=OUT / "panel_A_trio_denovo_polars_code.pdf",
        show_title=False,
        trim_whitespace=True,
    )


# Panel B — polars discovery snippet (compound scenario).
def panel_B() -> None:
    trio_discovery.render_trio_compound_het_polars_code(
        out_path=OUT / "panel_B_trio_compound_het_polars_code.pdf",
        show_title=False,
        trim_whitespace=True,
    )


PANELS = {
    "A": panel_A,
    "B": panel_B,
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
