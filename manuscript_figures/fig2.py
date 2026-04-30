"""Render Fig 2 panels for the tapestry manuscript.

The four panels strip down figures from the upstream nuclear-family wiki
page (`wiki/pedigree_wise_workflow/inheritance_mapping/nuclear_family/`)
to the paternal-haplotype rows only, with no figure titles, explanatory
text, or informative-site markers — the manuscript caption supplies that
context. The simulation/labeling logic is replicated from upstream
`generate_wiki.py` at the vendored SHA pinned in
`wiki/pedigree_wise_workflow/inheritance_mapping/README.md`
(`7448e5e946adbc7969ad5fd5e0730d7cace23a8d`).

Run:
    .venv/bin/python manuscript_figures/fig2.py
    .venv/bin/python manuscript_figures/fig2.py --panel A

Outputs land in `manuscript_figures/fig2/` as PNG, ready for Illustrator
placement.
"""
from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib.pyplot as plt

OUT = Path(__file__).resolve().parent / "fig2"
OUT.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Simulation — copied from upstream generate_wiki.py at the pinned SHA.
# ---------------------------------------------------------------------------
NUM_SITES = 9
HAP_DAD_ALPHA = "100110010"
HAP_DAD_BETA  = "010101011"
HAP_MOM_GAMMA = "101110010"
HAP_MOM_DELTA = "100010100"

KID1_LABELS = [("α", "γ")] * NUM_SITES
KID2_LABELS = [("β", "δ")] * NUM_SITES
KID3_LABELS = [("α", "γ")] * 4 + [("β", "γ")] * 5

KID_MISSING = {"Kid1": {8}, "Kid2": set(), "Kid3": set()}

_GREEK_TO_HAP = {
    "α": HAP_DAD_ALPHA,
    "β": HAP_DAD_BETA,
    "γ": HAP_MOM_GAMMA,
    "δ": HAP_MOM_DELTA,
}


def _hap_to_alleles(hap: str) -> List[int]:
    return [int(c) for c in hap]


def _child_genotype(dad_label: str, mom_label: str, site: int) -> Tuple[int, int]:
    paternal = _hap_to_alleles(_GREEK_TO_HAP[dad_label])[site]
    maternal = _hap_to_alleles(_GREEK_TO_HAP[mom_label])[site]
    return paternal, maternal


def _unphased_gt(a: int, b: int) -> str:
    if a == b:
        return f"{a}/{a}"
    return "0/1"


def _build_simulation() -> Dict:
    dad_alpha = _hap_to_alleles(HAP_DAD_ALPHA)
    dad_beta = _hap_to_alleles(HAP_DAD_BETA)
    mom_gamma = _hap_to_alleles(HAP_MOM_GAMMA)
    mom_delta = _hap_to_alleles(HAP_MOM_DELTA)

    def gts(labels):
        return [_child_genotype(la, lb, i) for i, (la, lb) in enumerate(labels)]

    kid1_phased = gts(KID1_LABELS)
    kid2_phased = gts(KID2_LABELS)
    kid3_phased = gts(KID3_LABELS)

    def kid_unphased(kid_name, phased):
        out = []
        for i, (a, b) in enumerate(phased):
            if i in KID_MISSING[kid_name]:
                out.append("./.")
            else:
                out.append(_unphased_gt(a, b))
        return out

    return {
        "kid1_unphased": kid_unphased("Kid1", kid1_phased),
        "kid2_unphased": kid_unphased("Kid2", kid2_phased),
        "kid3_unphased": kid_unphased("Kid3", kid3_phased),
        "dad_alpha": dad_alpha,
        "dad_beta": dad_beta,
        "mom_gamma": mom_gamma,
        "mom_delta": mom_delta,
        "kid1_phased": kid1_phased,
        "kid2_phased": kid2_phased,
        "kid3_phased": kid3_phased,
    }


def _informative_sites_dad(sim: Dict) -> List[int]:
    return [
        i for i in range(NUM_SITES)
        if sim["dad_alpha"][i] != sim["dad_beta"][i]
        and sim["mom_gamma"][i] == sim["mom_delta"][i]
    ]


def _informative_sites_mom(sim: Dict) -> List[int]:
    return [
        i for i in range(NUM_SITES)
        if sim["mom_gamma"][i] != sim["mom_delta"][i]
        and sim["dad_alpha"][i] == sim["dad_beta"][i]
    ]


def _per_site_parent_labels(
    sim: Dict, info_sites: List[int],
    parent_hap_a_key: str, parent_hap_b_key: str,
    other_hap_a_key: str, other_hap_b_key: str,
    kid_unphased_key: Dict[str, str],
    letter_first: str, letter_second: str,
) -> Tuple[Dict[str, List[str]], Dict[str, List[str]], Dict[str, List[str]]]:
    kids = ["Kid1", "Kid2", "Kid3"]
    stage1 = {k: ["?"] * NUM_SITES for k in kids}
    stage2 = {k: ["?"] * NUM_SITES for k in kids}
    stage3 = {k: ["?"] * NUM_SITES for k in kids}

    for s in info_sites:
        parent_alleles = {sim[parent_hap_a_key][s], sim[parent_hap_b_key][s]}
        other_alleles = {sim[other_hap_a_key][s], sim[other_hap_b_key][s]}
        unique_set = parent_alleles - other_alleles
        if not unique_set:
            continue
        unique = next(iter(unique_set))

        per_kid: Dict[str, str] = {}
        for k in kids:
            gt = sim[kid_unphased_key[k]][s]
            if gt == "./.":
                per_kid[k] = "?"
                continue
            kid_alleles = {int(gt[0]), int(gt[2])}
            per_kid[k] = letter_first if unique in kid_alleles else "?"
        for k in kids:
            stage1[k][s] = per_kid[k]

        if any(per_kid[k] == letter_first for k in kids):
            for k in kids:
                if per_kid[k] == "?":
                    per_kid[k] = letter_second
        for k in kids:
            stage2[k][s] = per_kid[k]

        first_count = sum(1 for k in kids if per_kid[k] == letter_first)
        second_count = sum(1 for k in kids if per_kid[k] == letter_second)
        if first_count < second_count:
            for k in kids:
                if per_kid[k] == letter_first:
                    per_kid[k] = letter_second
                elif per_kid[k] == letter_second:
                    per_kid[k] = letter_first
        for k in kids:
            stage3[k][s] = per_kid[k]

    return stage1, stage2, stage3


def _flip_only(
    per_site: Dict[str, List[str]],
    info_sites: List[int],
    letters: Tuple[str, str],
) -> Dict[str, List[str]]:
    kids = list(per_site.keys())
    out = {k: ["?"] * NUM_SITES for k in kids}
    a, b = letters
    prev: Optional[Dict[str, str]] = None
    for s in info_sites:
        cur = {k: per_site[k][s] for k in kids}
        if prev is not None:
            flipped = {
                k: (b if cur[k] == a else a if cur[k] == b else cur[k])
                for k in kids
            }
            same = sum(1 for k in kids if cur[k] == prev[k])
            flip = sum(1 for k in kids if flipped[k] == prev[k])
            if flip > same:
                cur = flipped
        for k in kids:
            out[k][s] = cur[k]
        prev = cur
    return out


# ---------------------------------------------------------------------------
# Rendering — same monospace text-block layout as upstream
# `_render_panel_image`, with no title or caption around the body.
# ---------------------------------------------------------------------------
def _render_panel_image(body_lines: List[str], out_path: Path) -> None:
    n_lines = max(len(body_lines), 1)
    fig_w = 12.0
    fig_h = 0.30 * n_lines + 0.4
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    fig.subplots_adjust(left=0.02, right=0.98, top=0.96, bottom=0.04)
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)

    top = 0.96
    line_h = 0.92 / n_lines
    font_size = 12
    for i, line in enumerate(body_lines):
        y = top - i * line_h
        ax.text(
            0.01, y, line,
            ha="left", va="top",
            fontsize=font_size, family="monospace",
            transform=ax.transAxes,
        )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
    _trim_whitespace_png(out_path)
    print(f"wrote {out_path}")


def _trim_whitespace_png(path: Path) -> None:
    """Crop pure-white margins off all four sides of a PNG in place."""
    from PIL import Image, ImageChops
    img = Image.open(path).convert("RGB")
    bg = Image.new("RGB", img.size, (255, 255, 255))
    diff = ImageChops.difference(img, bg)
    bbox = diff.getbbox()
    if bbox is not None:
        pad = 10
        left = max(0, bbox[0] - pad)
        top = max(0, bbox[1] - pad)
        right = min(img.size[0], bbox[2] + pad)
        bottom = min(img.size[1], bbox[3] + pad)
        img.crop((left, top, right, bottom)).save(path)


# ---------------------------------------------------------------------------
# Panels
# ---------------------------------------------------------------------------
KIDS = ["Kid1", "Kid2", "Kid3"]
KID_UNPHASED_KEY = {
    "Kid1": "kid1_unphased",
    "Kid2": "kid2_unphased",
    "Kid3": "kid3_unphased",
}


def _kid_pat_row(kid_phased) -> str:
    return " ".join(str(a) for a, _ in kid_phased)


def _paternal_only_rows(state: Dict[str, List[str]], info_sites: List[int]) -> List[str]:
    rows = []
    for k in KIDS:
        cells = [state[k][i] if i in info_sites else "." for i in range(NUM_SITES)]
        rows.append(f"{k} p:  " + " ".join(cells))
    return rows


def panel_A() -> None:
    """Ground-truth founder haplotypes + kid paternal allele rows."""
    sim = _build_simulation()
    body = [
        "Dad  α:  " + " ".join(str(x) for x in sim["dad_alpha"]),
        "Dad  β:  " + " ".join(str(x) for x in sim["dad_beta"]),
        "Mom  γ:  " + " ".join(str(x) for x in sim["mom_gamma"]),
        "Mom  δ:  " + " ".join(str(x) for x in sim["mom_delta"]),
        "",
        "Kid1 p:  " + _kid_pat_row(sim["kid1_phased"]),
        "Kid2 p:  " + _kid_pat_row(sim["kid2_phased"]),
        "Kid3 p:  " + _kid_pat_row(sim["kid3_phased"]),
    ]
    _render_panel_image(body, OUT / "panel_A_founder_haps.png")


def _pat_stage(sim, dad_info, stage_idx: int) -> Dict[str, List[str]]:
    pat_stage1, pat_stage2, pat_stage3 = _per_site_parent_labels(
        sim, dad_info,
        parent_hap_a_key="dad_alpha", parent_hap_b_key="dad_beta",
        other_hap_a_key="mom_gamma", other_hap_b_key="mom_delta",
        kid_unphased_key=KID_UNPHASED_KEY,
        letter_first="A", letter_second="B",
    )
    return [pat_stage1, pat_stage2, pat_stage3][stage_idx]


def panel_B() -> None:
    """fig3_2 — backfill non-carrier fill, paternal-only rows."""
    sim = _build_simulation()
    dad_info = _informative_sites_dad(sim)
    rows = _paternal_only_rows(_pat_stage(sim, dad_info, 1), dad_info)
    _render_panel_image(rows, OUT / "panel_B_backfill_before_swap.png")


def panel_C() -> None:
    """fig3_3 — backfill swap-by-majority, paternal-only rows."""
    sim = _build_simulation()
    dad_info = _informative_sites_dad(sim)
    rows = _paternal_only_rows(_pat_stage(sim, dad_info, 2), dad_info)
    _render_panel_image(rows, OUT / "panel_C_backfill_after_swap.png")


def panel_D() -> None:
    """fig4_1 — perform_flips_in_place #1, paternal-only rows."""
    sim = _build_simulation()
    dad_info = _informative_sites_dad(sim)
    pat_flipped = _flip_only(_pat_stage(sim, dad_info, 2), dad_info, ("A", "B"))
    rows = []
    for k in KIDS:
        cells = [pat_flipped[k][i] if pat_flipped[k][i] != "?" else "." for i in range(NUM_SITES)]
        rows.append(f"{k} p:  " + " ".join(cells))
    _render_panel_image(rows, OUT / "panel_D_after_flip.png")


PANELS = {"A": panel_A, "B": panel_B, "C": panel_C, "D": panel_D}


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
