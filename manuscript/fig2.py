"""Render Fig 2 for the tapestry manuscript.

Emits a single PDF, `manuscript/fig2/fig2.pdf`, that
stacks three panels with aligned columns:

  A. Ground truth — fig1 stripped of title/text: founder haps plus each
     kid's allele on both slots, with each cell filled by the inherited
     parental haplotype (A/B/C/D).
  B. Maternal — three tables: unphased VCF GTs for all 5 individuals at
     mom-informative sites; raw inherited maternal allele (fig5_1); and
     backfill swap-by-majority labels (fig3_3). A bottom track shows
     the ground-truth maternal haplotype as colored runs (C/D).
  C. Paternal — same three tables on the paternal side, plus the
     perform_flips_in_place step (fig4_1), which is a no-op on the
     maternal side and so is omitted there. A bottom track shows the
     ground-truth paternal haplotype (A/B); Kid3's row switches from A
     (blue) to B (orange) at the paternal recombination.

Informative sites for each panel are the columns whose top-row cells
are non-`.`. The trailing simulated site is dropped. PDFs are kept as
vector text so collaborators can edit labels directly inside Illustrator
after placing them.

Simulation/labeling logic is replicated from upstream `generate_wiki.py`
at the SHA pinned in
`wiki/pedigree_wise_workflow/inheritance_mapping/README.md`
(`7448e5e946adbc7969ad5fd5e0730d7cace23a8d`).

Run:
    .venv/bin/python manuscript/fig2.py
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
DISPLAY_SITES = 8  # render the first N sites only; the trailing site is dropped
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
# Rendering — single-page PDF with per-cell text and optional palette fill.
# Vector text survives placement in Illustrator so collaborators can edit
# labels; per-cell fills make haplotype provenance visible at a glance.
# ---------------------------------------------------------------------------
KIDS = ["Kid1", "Kid2", "Kid3"]
KID_PHASED_KEY = {"Kid1": "kid1_phased", "Kid2": "kid2_phased", "Kid3": "kid3_phased"}
KID_UNPHASED_KEY = {
    "Kid1": "kid1_unphased",
    "Kid2": "kid2_unphased",
    "Kid3": "kid3_unphased",
}

HAP_PALETTE = {
    "A": "#a6cee3",  # light blue   — dad hap1
    "B": "#fdbf6f",  # light orange — dad hap2
    "C": "#b2df8a",  # light green  — mom hap1
    "D": "#fb9a99",  # light pink   — mom hap2
}
_GREEK_TO_LETTER = {"α": "A", "β": "B", "γ": "C", "δ": "D"}
KID_LABELS_MAP = {"Kid1": KID1_LABELS, "Kid2": KID2_LABELS, "Kid3": KID3_LABELS}

# Row = (label, cells, colors, merge_runs[, bottom_colors]).
#   colors[i] is a palette key or None.
#   merge_runs=True draws one rectangle per run of identical colors with a
#   single centered letter, instead of one rectangle+glyph per cell.
#   bottom_colors (optional 5th element) — when supplied alongside
#   merge_runs=True, each merged rectangle is split horizontally: top half
#   uses colors[j], bottom half uses bottom_colors[j], and a run extends
#   only while *both* colour streams stay constant. Used by panel D
#   (bilateral IBD blocks).
Row = Tuple
SPACER: Row = ("", [], [], False)


def _plain(label: str, cells: List[str]) -> Row:
    return (label, cells, [None] * len(cells), False)


def _draw_panel(
    ax,
    rows: List[Row],
    pitch: float,
    cells_x0: float,
) -> None:
    n = max(len(rows), 1)
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)

    line_h = 1.0 / n
    font_size = 11
    label_x = 0.02
    cell_w = pitch
    cell_h = line_h * 0.85

    for i, row in enumerate(rows):
        label, cells, colors, merge_runs = row[:4]
        bottom_colors = row[4] if len(row) > 4 else None
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
                    width = (k - j) * pitch + cell_w
                    if bottom_colors is None:
                        ax.add_patch(plt.Rectangle(
                            (x_left, y - cell_h / 2),
                            width, cell_h,
                            facecolor=HAP_PALETTE[run_color],
                            edgecolor="none", zorder=1,
                        ))
                    else:
                        # Split top/bottom: top half = paternal colour,
                        # bottom half = maternal colour. Letter centred on
                        # the rectangle's mid-line. Where another bilateral
                        # block follows (or precedes), inset the rectangle
                        # by half a gap on the boundary-facing side so the
                        # whitespace straddles the true boundary rather
                        # than sitting entirely on one side.
                        half_gap = pitch * 0.10
                        gap_left  = half_gap if j > 0                 else 0.0
                        gap_right = half_gap if k + 1 < len(cells)    else 0.0
                        rect_x = x_left + gap_left
                        rect_w = width - gap_left - gap_right
                        ax.add_patch(plt.Rectangle(
                            (rect_x, y),
                            rect_w, cell_h / 2,
                            facecolor=HAP_PALETTE[run_color],
                            edgecolor="none", zorder=1,
                        ))
                        ax.add_patch(plt.Rectangle(
                            (rect_x, y - cell_h / 2),
                            rect_w, cell_h / 2,
                            facecolor=HAP_PALETTE[bottom_colors[j]],
                            edgecolor="none", zorder=1,
                        ))
                    cx = cells_x0 + (j + k) * pitch / 2
                    ax.text(
                        cx, y, cells[j], ha="center", va="center",
                        fontsize=font_size, family="monospace",
                        color="black", zorder=2,
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
            ax.text(
                x, y, c, ha="center", va="center",
                fontsize=font_size, family="monospace",
                color="black", zorder=2,
            )


def _render_combined_pdf(
    panels: List[List[Row]],
    out_path: Path,
    fig_w: float = 5.5,
    pitch: float = 0.105,
    cells_x0: float = 0.22,
    row_height_in: float = 0.30,
    panel_gap_rows: float = 1.5,
) -> None:
    """Stack multiple panels in one figure, sharing x layout so columns
    align across panels. Each panel's vertical extent is proportional to
    its row count, giving every line the same physical height."""
    counts = [max(len(p), 1) for p in panels]
    total_rows = sum(counts)
    fig_h = row_height_in * total_rows + row_height_in * panel_gap_rows * (len(panels) - 1) + 0.3
    fig, axes = plt.subplots(
        nrows=len(panels), ncols=1,
        figsize=(fig_w, fig_h),
        gridspec_kw={"height_ratios": counts, "hspace": panel_gap_rows / (total_rows / len(panels))},
    )
    if len(panels) == 1:
        axes = [axes]
    for ax, rows in zip(axes, panels):
        _draw_panel(ax, rows, pitch=pitch, cells_x0=cells_x0)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, bbox_inches="tight", pad_inches=0.1)
    plt.close(fig)
    print(f"wrote {out_path}")


# ---------------------------------------------------------------------------
# Sub-table builders
# ---------------------------------------------------------------------------
def _ground_truth_rows(sim: Dict) -> List[Row]:
    rows: List[Row] = []

    def founder(label, hap_key, color_key):
        cells = [str(x) for x in sim[hap_key][:DISPLAY_SITES]]
        rows.append((label, cells, [color_key] * DISPLAY_SITES, False))

    founder("Dad  A:", "dad_alpha", "A")
    founder("Dad  B:", "dad_beta", "B")
    founder("Mom  C:", "mom_gamma", "C")
    founder("Mom  D:", "mom_delta", "D")
    rows.append(SPACER)

    for k in KIDS:
        labels = KID_LABELS_MAP[k]
        phased = sim[KID_PHASED_KEY[k]]
        cells_p = [str(phased[i][0]) for i in range(DISPLAY_SITES)]
        colors_p = [_GREEK_TO_LETTER[labels[i][0]] for i in range(DISPLAY_SITES)]
        rows.append((f"{k} p:", cells_p, colors_p, False))
        cells_m = [str(phased[i][1]) for i in range(DISPLAY_SITES)]
        colors_m = [_GREEK_TO_LETTER[labels[i][1]] for i in range(DISPLAY_SITES)]
        rows.append((f"{k} m:", cells_m, colors_m, False))
    return rows


def _unphased_section(sim: Dict, info_sites: List[int]) -> List[Row]:
    """Top of paternal/maternal panels: unphased GT for all 5 individuals,
    shown only at the parent-specific informative sites."""
    dad_gt = [_unphased_gt(sim["dad_alpha"][i], sim["dad_beta"][i]) for i in range(DISPLAY_SITES)]
    mom_gt = [_unphased_gt(sim["mom_gamma"][i], sim["mom_delta"][i]) for i in range(DISPLAY_SITES)]

    def gated(values):
        return [v if i in info_sites else "." for i, v in enumerate(values)]

    rows: List[Row] = [
        _plain("Dad:", gated(dad_gt)),
        _plain("Mom:", gated(mom_gt)),
    ]
    for k in KIDS:
        unphased = sim[KID_UNPHASED_KEY[k]][:DISPLAY_SITES]
        rows.append(_plain(f"{k}:", gated(unphased)))
    return rows


def _inherited_section(sim: Dict, info_sites: List[int], side: str) -> List[Row]:
    """fig5_1 — raw 0/1 allele on the informative slot, `.` elsewhere."""
    slot = 0 if side == "paternal" else 1
    suffix = "p" if side == "paternal" else "m"
    rows: List[Row] = []
    for k in KIDS:
        phased = sim[KID_PHASED_KEY[k]]
        cells: List[str] = []
        for i in range(DISPLAY_SITES):
            if i not in info_sites:
                cells.append(".")
            elif i in KID_MISSING[k]:
                cells.append("?")
            else:
                cells.append(str(phased[i][slot]))
        rows.append(_plain(f"{k} {suffix}:", cells))
    return rows


def _collapse_to_blocks(
    state: Dict[str, List[str]],
    info_sites: List[int],
) -> Dict[str, List[str]]:
    """Fill non-informative sites by extending each kid's most-recent
    informative-site letter forward (and the first informative letter
    backward), mirroring collapse_identical_iht's `?`-as-wildcard merge.
    The result is a contiguous A/B (or C/D) string per kid, ready to be
    drawn with merge_runs=True so equal-letter runs render as one block."""
    kids = list(state.keys())
    out = {k: ["?"] * DISPLAY_SITES for k in kids}
    for k in kids:
        first = next((state[k][i] for i in info_sites if state[k][i] != "?"), "?")
        cur = first
        for i in range(DISPLAY_SITES):
            if i in info_sites and state[k][i] != "?":
                cur = state[k][i]
            out[k][i] = cur
    return out


def _reconstruction_track_rows(
    state: Dict[str, List[str]],
    info_sites: List[int],
    side: str,
) -> List[Row]:
    """Bottom table: the per-site reconstruction collapsed into blocks
    by run-merging identical letters across (informative + filled-in)
    sites. Cells are colored by their A/B/C/D letter."""
    suffix = "p" if side == "paternal" else "m"
    filled = _collapse_to_blocks(state, info_sites)
    rows: List[Row] = []
    for k in KIDS:
        letters = filled[k]
        rows.append((f"{k} {suffix}:", letters, list(letters), True))
    return rows


def _stage(sim: Dict, info_sites: List[int], side: str, stage_idx: int) -> Dict[str, List[str]]:
    if side == "paternal":
        s1, s2, s3 = _per_site_parent_labels(
            sim, info_sites,
            parent_hap_a_key="dad_alpha", parent_hap_b_key="dad_beta",
            other_hap_a_key="mom_gamma", other_hap_b_key="mom_delta",
            kid_unphased_key=KID_UNPHASED_KEY,
            letter_first="A", letter_second="B",
        )
    else:
        s1, s2, s3 = _per_site_parent_labels(
            sim, info_sites,
            parent_hap_a_key="mom_gamma", parent_hap_b_key="mom_delta",
            other_hap_a_key="dad_alpha", other_hap_b_key="dad_beta",
            kid_unphased_key=KID_UNPHASED_KEY,
            letter_first="C", letter_second="D",
        )
    return [s1, s2, s3][stage_idx]


def _label_rows(state: Dict[str, List[str]], info_sites: List[int], side: str) -> List[Row]:
    suffix = "p" if side == "paternal" else "m"
    rows: List[Row] = []
    for k in KIDS:
        cells = [state[k][i] if i in info_sites else "." for i in range(DISPLAY_SITES)]
        rows.append(_plain(f"{k} {suffix}:", cells))
    return rows


def _flipped_rows(state: Dict[str, List[str]], info_sites: List[int], side: str) -> List[Row]:
    letters = ("A", "B") if side == "paternal" else ("C", "D")
    suffix = "p" if side == "paternal" else "m"
    flipped = _flip_only(state, info_sites, letters)
    rows: List[Row] = []
    for k in KIDS:
        cells = [flipped[k][i] if flipped[k][i] != "?" else "." for i in range(DISPLAY_SITES)]
        rows.append(_plain(f"{k} {suffix}:", cells))
    return rows


def _compose(side: str) -> List[Row]:
    sim = _build_simulation()
    info_full = _informative_sites_dad(sim) if side == "paternal" else _informative_sites_mom(sim)
    info = [i for i in info_full if i < DISPLAY_SITES]
    body: List[Row] = []
    body += _unphased_section(sim, info)
    body.append(SPACER)
    body += _inherited_section(sim, info, side)
    body.append(SPACER)
    stage3 = _stage(sim, info, side, 2)
    body += _label_rows(stage3, info, side)  # fig3_3
    if side == "paternal":
        body.append(SPACER)
        body += _flipped_rows(stage3, info, side)  # fig4_1
        flipped = _flip_only(stage3, info, ("A", "B"))
        body.append(SPACER)
        body += _reconstruction_track_rows(flipped, info, side)
    else:
        body.append(SPACER)
        body += _reconstruction_track_rows(stage3, info, side)
    return body


def _bilateral_block_rows(sim: Dict) -> List[Row]:
    """Panel D: one row per kid, listing the bilateral IBD segments as
    labelled rectangles whose top half is filled with the paternal-letter
    colour and bottom half with the maternal-letter colour. The label is
    the founder-haplotype pair (e.g. ``A|C``); a colour change in either
    stream marks a bilateral-block boundary, so Kid3's paternal
    recombination splits its row into two adjacent rectangles
    (``A|C`` then ``B|C``)."""
    info_pat = [i for i in _informative_sites_dad(sim) if i < DISPLAY_SITES]
    info_mat = [i for i in _informative_sites_mom(sim) if i < DISPLAY_SITES]

    pat_state = _stage(sim, info_pat, "paternal", 2)
    pat_state = _flip_only(pat_state, info_pat, ("A", "B"))
    pat_filled = _collapse_to_blocks(pat_state, info_pat)

    mat_state = _stage(sim, info_mat, "maternal", 2)
    mat_filled = _collapse_to_blocks(mat_state, info_mat)

    rows: List[Row] = []
    for k in KIDS:
        pat = pat_filled[k]
        mat = mat_filled[k]
        cells = [f"{pat[i]}|{mat[i]}" for i in range(DISPLAY_SITES)]
        rows.append((f"{k}:", cells, list(pat), True, list(mat)))
    return rows


def main() -> None:
    sim = _build_simulation()
    panels = [
        _ground_truth_rows(sim),
        _compose("maternal"),
        _compose("paternal"),
        _bilateral_block_rows(sim),
    ]
    _render_combined_pdf(panels, OUT / "fig2.pdf")


if __name__ == "__main__":
    main()
