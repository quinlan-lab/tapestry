"""Render Fig 3 for the tapestry manuscript as a single combined PDF.

Fig 3 walks the reader through exhaustive-enumeration phasing at a clean
non-informative site (N1). Algorithm inputs (observed unphased genotypes
and the inferred founder-haplotype label each kid carries) are coloured
red wherever they appear; match ticks are coloured green. The winning
row (0 mismatches) is boxed.

Ground truth chosen so that some kids show a change of allele order on
phasing (`1|0` → unphased `0/1`) and at least one kid is homozygous, so
the enumeration has a unique winner (otherwise the two-fold
parent-relabeling symmetry at an all-heterozygous non-informative site
ties two orientations at zero mismatches).

Orientation-enumeration logic is replicated from upstream
`Platinum-Pedigree-Inheritance/wiki/generate_wiki.py` at the SHA pinned
in `wiki/pedigree_wise_workflow/inheritance_mapping/README.md`
(`7448e5e946adbc7969ad5fd5e0730d7cace23a8d`); see also
`wiki/pedigree_wise_workflow/inheritance_mapping/concordance/concordance.md`
section 3.

Run:
    .venv/bin/python manuscript/fig3.py
"""
from __future__ import annotations

import re
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

OUT = Path(__file__).resolve().parent / "fig3"
OUT.mkdir(parents=True, exist_ok=True)


# Site N1 ground truth (non-informative: dad=0/1, mom=0/1).
N1_TRUTH = {"A": 1, "B": 0, "C": 0, "D": 1}

BLOCK_LABELS: Dict[str, Tuple[str, str]] = {
    "Kid1": ("A", "C"),
    "Kid2": ("B", "D"),
    "Kid3": ("A", "D"),
}

ORIENTATIONS: List[Tuple[Tuple[str, str], Tuple[str, str]]] = [
    (("A", "B"), ("C", "D")),
    (("B", "A"), ("C", "D")),
    (("A", "B"), ("D", "C")),
    (("B", "A"), ("D", "C")),
]


def _sorted_unphased(pair: Tuple[int, int]) -> Tuple[int, int]:
    return (min(pair), max(pair))


def _site_observed() -> Dict[str, Tuple[int, int]]:
    a, b, c, d = N1_TRUTH["A"], N1_TRUTH["B"], N1_TRUTH["C"], N1_TRUTH["D"]
    truth = {
        "Dad": (a, b),
        "Mom": (c, d),
        "Kid1": (a, c),
        "Kid2": (b, d),
        "Kid3": (a, d),
    }
    return {k: _sorted_unphased(v) for k, v in truth.items()}


def _expected_under_orientation(
    dad_letters: Tuple[str, str],
    mom_letters: Tuple[str, str],
    dad_vcf: Tuple[int, int],
    mom_vcf: Tuple[int, int],
) -> Dict[str, Tuple[int, int]]:
    letter_to_allele: Dict[str, int] = {
        dad_letters[0]: dad_vcf[0],
        dad_letters[1]: dad_vcf[1],
        mom_letters[0]: mom_vcf[0],
        mom_letters[1]: mom_vcf[1],
    }
    expected: Dict[str, Tuple[int, int]] = {
        "Dad": (letter_to_allele[dad_letters[0]], letter_to_allele[dad_letters[1]]),
        "Mom": (letter_to_allele[mom_letters[0]], letter_to_allele[mom_letters[1]]),
    }
    for kid, (lp, lm) in BLOCK_LABELS.items():
        expected[kid] = (letter_to_allele[lp], letter_to_allele[lm])
    return expected


def _count_mismatches(
    observed: Dict[str, Tuple[int, int]],
    expected: Dict[str, Tuple[int, int]],
) -> int:
    return sum(
        1 for k in observed
        if _sorted_unphased(observed[k]) != _sorted_unphased(expected[k])
    )


def _kid_hap_label(kid: str) -> str:
    p, m = BLOCK_LABELS[kid]
    return f"{p}|{m}"


# ---------------------------------------------------------------------------
# Layout constants and colours.
# ---------------------------------------------------------------------------
COL_BASE = 0.85
COL_WIDE = 1.15   # unphased kid columns (need room for tick / right-side gap)
COL_MATCH = 1.40  # "# Matches" header column
ROW_H = 0.42
LEFT_OFFSET = 0.2
RIGHT_PAD = 0.2

# Per-column widths (inches). Cols 5, 7, 9 are kids' "expected unphased"
# columns (need horizontal room between the unphased number, the tick,
# and the right-side separator); col 10 is "# Matches".
COL_WIDTHS: List[float] = [
    COL_BASE, COL_BASE,       # A, B  (Dad)
    COL_BASE, COL_BASE,       # C, D  (Mom)
    COL_BASE, COL_WIDE,       # Kid1 phased, Kid1 unphased
    COL_BASE, COL_WIDE,       # Kid2 phased, Kid2 unphased
    COL_BASE, COL_WIDE,       # Kid3 phased, Kid3 unphased
    COL_MATCH,                # # Matches
]
TICK_RIGHT_INSET = 0.20       # inches from the right edge of the unphased cell

HEADER_FONT = 16     # super, group, and sub all share this size
BODY_FONT = 16
TICK_FONT = 28       # body ticks rendered larger for emphasis
PROSE_FONT = 17

INPUT_COLOR = "#cc0000"   # red: algorithm inputs
TICK_COLOR  = "#1b8a1b"   # green: match tick
BOX_COLOR   = "#1f77b4"   # blue: winning-row box

N_COLS = 11   # 0..3 parents (A,B,C,D); 4..9 kids (phased,unphased) × 3; 10 = # Mis


def _build_enum_rows(obs):
    dad_vcf = obs["Dad"]
    mom_vcf = obs["Mom"]
    rows = []
    for dad_letters, mom_letters in ORIENTATIONS:
        expected = _expected_under_orientation(dad_letters, mom_letters, dad_vcf, mom_vcf)
        letter_to_allele = {
            dad_letters[0]: dad_vcf[0], dad_letters[1]: dad_vcf[1],
            mom_letters[0]: mom_vcf[0], mom_letters[1]: mom_vcf[1],
        }
        a, b, c, d = (letter_to_allele[L] for L in "ABCD")
        n_mis = _count_mismatches(obs, expected)
        kids = []
        for k in ("Kid1", "Kid2", "Kid3"):
            ph = expected[k]
            unph = _sorted_unphased(ph)
            kids.append((ph, unph, unph == obs[k]))
        rows.append(((a, b, c, d), kids, n_mis))
    rows.sort(key=lambda r: r[0])
    return rows


# Segment = (text, colour). Used to render multi-coloured text in one line.
Segment = Tuple[str, str]


def _wrap_segments(segments, max_w_inches, fontsize, family, measure):
    """Greedy word-wrap a list of segments to fit in `max_w_inches`.

    Returns a list of lines, each a list of segments. Continuation lines
    inherit the leading-whitespace indent of the original first segment so
    bulleted / numbered items hang properly.
    """
    if not segments:
        return [[]]

    first_text = segments[0][0]
    # Continuation indent matches the start of the line's content. For a
    # numbered/bulleted item like "    2. ...", the hanging indent extends
    # past the number+period+space so continuation lines align with the
    # body of the item rather than the leading "2.".
    m = re.match(r"^(\s*)((?:\d+\.|[•·–-])\s+)", first_text)
    if m:
        indent = m.group(1) + " " * len(m.group(2))
    else:
        indent = first_text[:len(first_text) - len(first_text.lstrip(" "))]

    # Tokenise: split each segment into runs of spaces and runs of non-spaces,
    # tagging each token with its colour.
    tokens = []   # (text, colour, is_space)
    for seg_text, color in segments:
        if not seg_text:
            continue
        for m in re.finditer(r"\s+|\S+", seg_text):
            tok = m.group()
            tokens.append((tok, color, tok.isspace()))

    lines = [[]]
    cur_w = 0.0
    pending = []  # buffered spaces; dropped on wrap

    for text, color, is_sp in tokens:
        if is_sp:
            pending.append((text, color, measure(text, fontsize, family)))
            continue
        w = measure(text, fontsize, family)
        pend_w = sum(p[2] for p in pending)
        if lines[-1] and cur_w + pend_w + w > max_w_inches:
            # Wrap.
            lines.append([])
            cur_w = 0.0
            pending = []
            if indent:
                lines[-1].append((indent, color))
                cur_w += measure(indent, fontsize, family)
        for p_text, p_color, p_w in pending:
            lines[-1].append((p_text, p_color))
            cur_w += p_w
        pending = []
        lines[-1].append((text, color))
        cur_w += w

    return lines


def _render_figure(out_path: Path) -> None:
    obs = _site_observed()
    enum_rows = _build_enum_rows(obs)

    win_idx = next(i for i, (_, _, m) in enumerate(enum_rows) if m == 0)
    win_kids = enum_rows[win_idx][1]
    deduced = {
        name: f"{ph[0]}|{ph[1]}"
        for name, (ph, _, _) in zip(("K1", "K2", "K3"), win_kids)
    }

    def fmt_unph(p): return f"{p[0]}/{p[1]}"
    def fmt_ph(p):   return f"{p[0]}|{p[1]}"

    obs_strs = {k: fmt_unph(obs[k]) for k in ("Dad", "Mom", "Kid1", "Kid2", "Kid3")}

    # ---- Prose lines (lists of coloured segments) ----
    BLK = "black"
    R = INPUT_COLOR
    G = TICK_COLOR

    inputs_lines: List[List[Segment]] = [
        [("Inputs", BLK)],
        [("  • Observed unphased genotypes:  Dad = ", BLK), (obs_strs["Dad"], R),
         (",  Mom = ", BLK),  (obs_strs["Mom"], R),
         (",  Kid1 = ", BLK), (obs_strs["Kid1"], R),
         (",  Kid2 = ", BLK), (obs_strs["Kid2"], R),
         (",  Kid3 = ", BLK), (obs_strs["Kid3"], R), (".", BLK)],
        [("  • Inferred founder haplotype each kid carries:  Kid1 = ", BLK),
         (_kid_hap_label("Kid1"), R),
         (",  Kid2 = ", BLK), (_kid_hap_label("Kid2"), R),
         (",  Kid3 = ", BLK), (_kid_hap_label("Kid3"), R), (".", BLK)],
    ]

    algo_lines: List[List[Segment]] = [
        [("Algorithm", BLK)],
        [("  Enumerate the 4 ways to assign alleles (0 and 1) to dad's two haplotypes A, B and mom's two haplotypes C, D consistent with the parents' unphased genotypes. For each such assignment:", BLK)],
        [("    1. derive each kid's expected phased genotype by reading off the alleles assigned to the two founder haplotypes it carries (e.g., A|C → dad's A allele assignment paired with mom's C allele assignment);", BLK)],
        [("    2. for each kid, unphase the expected genotype and tick (", BLK), ("✓", G),
         (") if it matches the kid's observed unphased genotype (from Inputs);", BLK)],
        [("    3. count kids whose expected unphased genotype agrees with observed (# Matches).", BLK)],
    ]

    output_lines: List[List[Segment]] = [
        [("Output", BLK)],
        [("  The assignment with the most matches (boxed) gives the deduced phased genotypes:", BLK)],
        [(f"    Kid1 = {deduced['K1']},  Kid2 = {deduced['K2']},  Kid3 = {deduced['K3']}.", BLK)],
    ]

    # ---- Vertical layout ----
    line_h = 0.36
    gap = 0.30
    pad_top = pad_bot = 0.15
    table_h = 7 * ROW_H  # super, group, sub, 4 body
    fig_w = LEFT_OFFSET + sum(COL_WIDTHS) + RIGHT_PAD
    prose_family = ["Arial", "DejaVu Sans"]
    max_prose_w = sum(COL_WIDTHS)

    def is_section_header(line):
        return len(line) == 1 and line[0][0] in ("Inputs", "Algorithm", "Output")

    # ---- Pass 1: measure on a temp figure, wrap each prose line ----
    measure_fig = plt.figure(figsize=(fig_w, 1))
    measure_ax = measure_fig.add_axes([0, 0, 1, 1])
    measure_ax.set_xlim(0, fig_w)
    measure_ax.set_ylim(0, 1)
    measure_ax.set_axis_off()
    measure_fig.canvas.draw()
    m_renderer = measure_fig.canvas.get_renderer()

    def _measure_inches_static(text, fontsize, family):
        t = measure_ax.text(0, 0, text, fontsize=fontsize, family=family)
        bbox = t.get_window_extent(renderer=m_renderer)
        t.remove()
        return bbox.width / measure_fig.dpi

    def wrap_block(lines):
        out = []
        for line in lines:
            sublines = _wrap_segments(line, max_prose_w, PROSE_FONT,
                                      prose_family, _measure_inches_static)
            weight = "bold" if is_section_header(line) else "normal"
            for sub in sublines:
                out.append((sub, weight))
        return out

    inputs_w = wrap_block(inputs_lines)
    algo_w   = wrap_block(algo_lines)
    output_w = wrap_block(output_lines)
    plt.close(measure_fig)

    n_top = len(inputs_w) + 1 + len(algo_w)  # blank line between Inputs and Algorithm
    n_bot = len(output_w)
    top_text_h = n_top * line_h
    bot_text_h = n_bot * line_h
    fig_h = pad_top + top_text_h + gap + table_h + gap + bot_text_h + pad_bot

    fig = plt.figure(figsize=(fig_w, fig_h))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, fig_w)
    ax.set_ylim(0, fig_h)
    ax.set_axis_off()

    def y_from_top(y_top): return fig_h - y_top
    def col_x(j): return LEFT_OFFSET + sum(COL_WIDTHS[:j])

    # ---- Helpers ----
    # Renderer needed to measure proportional-font text widths.
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    _inv = ax.transData.inverted()

    def _measure_inches(text, fontsize, family):
        t = ax.text(0, 0, text, fontsize=fontsize, family=family)
        bbox = t.get_window_extent(renderer=renderer)
        (x0d, _), (x1d, _) = _inv.transform([(0, 0), (bbox.width, 0)])
        t.remove()
        return x1d - x0d

    def write_segments(x0, y_center, segments, fontsize, family="monospace", ha="left", weight="normal"):
        """Render segments left-to-right. If ha='center', the run is centred on x0."""
        widths = [_measure_inches(text, fontsize, family) for text, _ in segments]
        total_w = sum(widths)
        cur_x = x0 - total_w / 2 if ha == "center" else x0
        for (text, color), w in zip(segments, widths):
            ax.text(cur_x, y_center, text,
                    fontsize=fontsize, family=family, color=color,
                    weight=weight, ha="left", va="center")
            cur_x += w

    def write_line(segments, y_top, fontsize=PROSE_FONT,
                   family=["Arial", "DejaVu Sans"], weight="normal"):
        write_segments(LEFT_OFFSET, y_from_top(y_top + line_h * 0.5),
                       segments, fontsize, family=family, ha="left", weight=weight)

    def cell(x_col, y_row, w_cols, h_rows, segments, fontsize, fc="white"):
        """Draw a rectangle and centre `segments` (list of (text,color)) inside it."""
        x0 = col_x(x_col)
        x1 = col_x(x_col + w_cols)
        cell_w = x1 - x0
        y_top = table_top + y_row * ROW_H
        y_bot = y_top + h_rows * ROW_H
        rect_y = y_from_top(y_bot)
        ax.add_patch(Rectangle(
            (x0, rect_y), cell_w, h_rows * ROW_H,
            facecolor=fc, edgecolor="none",
        ))
        if segments:
            write_segments(
                x0 + cell_w / 2,
                rect_y + h_rows * ROW_H / 2,
                segments, fontsize, ha="center",
            )

    # ---- Top prose (wrapped) ----
    y = pad_top
    for sub, weight in inputs_w:
        write_line(sub, y, weight=weight); y += line_h
    y += line_h  # blank between sections
    for sub, weight in algo_w:
        write_line(sub, y, weight=weight); y += line_h

    y += gap
    table_top = y

    # ---- Table ----
    # Reversed shading: group row gets the lighter shade; sub-header row
    # gets the darker shade.
    HEADER_FC = "#f4f4f4"
    SUB_FC    = "#e4e4e4"

    # Super-header (row 0) — no shading.
    cell(0, 0, 4, 1, [("Hap-allele assignment", BLK)], HEADER_FONT)
    cell(4, 0, 6, 1, [("Expected kid genotypes", BLK)], HEADER_FONT)

    # Group row (row 1). Dad/Mom carry their observed unphased genotype as a
    # parenthetical (no tick column gives a redundant cue for them); kid
    # columns rely on the tick instead.
    groups = [
        ([("Dad (", BLK), (obs_strs["Dad"], R), (")", BLK)], 0,  2),
        ([("Mom (", BLK), (obs_strs["Mom"], R), (")", BLK)], 2,  4),
        ([("Kid1", BLK)],    4,  6),
        ([("Kid2", BLK)],    6,  8),
        ([("Kid3", BLK)],    8, 10),
        ([("# Matches", BLK)],10, 11),
    ]
    for segs, x0, x1 in groups:
        cell(x0, 1, x1 - x0, 1, segs, HEADER_FONT, fc=HEADER_FC)

    # Sub-header row (row 2) — founder-hap labels are inputs → red.
    sub_segs = [
        [("A", BLK)], [("B", BLK)], [("C", BLK)], [("D", BLK)],
        [(_kid_hap_label("Kid1"), R)], [],
        [(_kid_hap_label("Kid2"), R)], [],
        [(_kid_hap_label("Kid3"), R)], [],
        [],
    ]
    for j, segs in enumerate(sub_segs):
        cell(j, 2, 1, 1, segs, HEADER_FONT, fc=SUB_FC)

    # Body rows.
    for i, (founders, kids, n_mis) in enumerate(enum_rows):
        ry = 3 + i
        for j, val in enumerate(founders):
            cell(j, ry, 1, 1, [(str(val), BLK)], BODY_FONT)
        for ki, (ph, unph, match) in enumerate(kids):
            ph_x = 4 + ki * 2
            cell(ph_x, ry, 1, 1, [(fmt_ph(ph), BLK)], BODY_FONT)
            # Render the unphased number always centred in its cell (so the
            # numbers line up across rows regardless of whether a tick is
            # present); the tick is drawn separately, larger, to its right.
            cell(ph_x + 1, ry, 1, 1, [(fmt_unph(unph), BLK)], BODY_FONT)
            if match:
                # Tick sits a fixed inset from the right edge of the cell so
                # there is whitespace between it and the vertical separator.
                tick_x = col_x(ph_x + 2) - TICK_RIGHT_INSET
                tick_y = y_from_top(table_top + ry * ROW_H + ROW_H / 2)
                t = ax.text(tick_x, tick_y, "✔",
                            fontsize=TICK_FONT, family="monospace",
                            color=G, weight="bold", ha="center", va="center")
                t.set_path_effects([
                    pe.Stroke(linewidth=0.8, foreground=G),
                    pe.Normal(),
                ])
        n_match = sum(1 for _, _, m in kids if m)
        cell(10, ry, 1, 1, [(str(n_match), BLK)], BODY_FONT)

    # Separators.
    line_kw = dict(color="black", linewidth=0.5)
    x_left = col_x(0)
    x_right = col_x(N_COLS)
    # Horizontal under sub-row (between sub-header and body).
    y_hbar = y_from_top(table_top + 3 * ROW_H)
    ax.plot([x_left, x_right], [y_hbar, y_hbar], **line_kw)
    # Vertical separators at group boundaries. Most start below the
    # super-header row so they do not cross the super-header text. The
    # boundaries between super-header spans (col 4: between "Hap-allele
    # assignment" and "Expected kid genotypes"; col 10: right edge of
    # "Expected kid genotypes") extend up to the top of the super-header.
    sep_top_below_super = y_from_top(table_top + 1 * ROW_H)
    sep_top_full        = y_from_top(table_top + 0 * ROW_H)
    sep_bottom          = y_from_top(table_top + 7 * ROW_H)
    for xc in (2, 4, 6, 8, 10):
        x = col_x(xc)
        top_y = sep_top_full if xc in (4, 10) else sep_top_below_super
        ax.plot([x, x], [top_y, sep_bottom], **line_kw)

    # Box around winning row.
    win_y_top = table_top + (3 + win_idx) * ROW_H
    ax.add_patch(Rectangle(
        (x_left, y_from_top(win_y_top + ROW_H)),
        sum(COL_WIDTHS), ROW_H,
        facecolor="none", edgecolor=BOX_COLOR, linewidth=1.8,
    ))

    # ---- Bottom prose (wrapped) ----
    y = table_top + 7 * ROW_H + gap
    for sub, weight in output_w:
        write_line(sub, y, weight=weight); y += line_h

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)
    print(f"wrote {out_path}")


def main() -> None:
    _render_figure(OUT / "fig3.pdf")


if __name__ == "__main__":
    main()
