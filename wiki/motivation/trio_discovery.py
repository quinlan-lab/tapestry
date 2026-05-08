"""SVG trio with haplotype lollipops + matplotlib companion.

Generates:
  - trio_denovo.svg: de novo methylation scenario
  - trio_compound_het.svg: compound genetic-epigenetic heterozygote
  - trio_denovo_bed.png: BED snippet + polars query for discovering de novo epimutations
  - trio_compound_het_bed.png: BED snippet + polars query for compound genetic-epigenetic heterozygotes

Run: .venv/bin/python wiki/motivation/trio_discovery.py
"""
from __future__ import annotations
from dataclasses import dataclass, field
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle

mpl.rcParams["font.family"] = "Arial"

COLOR_NEUTRAL = "#444444"
OUT = Path(__file__).resolve().parent

X0, Y0, VW, VH = 220, 125, 880, 335
SVG_HEAD = (
    f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{X0} {Y0} {VW} {VH}" '
    f'font-family="Arial, sans-serif" font-size="14">'
    f'<rect x="{X0}" y="{Y0}" width="{VW}" height="{VH}" fill="white"/>'
)
SVG_TAIL = '</svg>'


# ---------- primitives ----------

def rect(x, y, w, h, fill='none', stroke='black', sw=1.5, dash=None, rx=4):
    da = f' stroke-dasharray="{dash}"' if dash else ''
    return (f'<rect x="{x}" y="{y}" width="{w}" height="{h}" fill="{fill}" '
            f'stroke="{stroke}" stroke-width="{sw}"{da} rx="{rx}"/>')


def circle(cx, cy, r, fill='white', stroke='black', sw=1.5):
    return (f'<circle cx="{cx}" cy="{cy}" r="{r}" fill="{fill}" '
            f'stroke="{stroke}" stroke-width="{sw}"/>')


def line(x1, y1, x2, y2, sw=1.5):
    return f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="black" stroke-width="{sw}"/>'


def star(cx, cy, r=7, fill='#d62728', stroke='black', sw=1.0):
    """5-point star marker for a genetic variant."""
    import math
    pts = []
    for i in range(10):
        ang = -math.pi / 2 + i * math.pi / 5
        rr = r if i % 2 == 0 else r * 0.45
        pts.append(f'{cx + rr * math.cos(ang):.2f},{cy + rr * math.sin(ang):.2f}')
    return (f'<polygon points="{" ".join(pts)}" fill="{fill}" '
            f'stroke="{stroke}" stroke-width="{sw}"/>')


# ---------- haplotype ----------

def haplotype(x, y, w, color, cpgs, variants=(), label=None, label_at_start=False,
              label_fontsize=16, stems_below=False, tss=None, right_extent=None,
              bar_length=None):
    """cpgs: list of (rel_x, methylated_bool). variants: list of rel_x drawn as red stars.
    label: optional single-letter hap id (drawn just outside the bar, adjacent
    to the person's shape). label_at_start puts it at bar_x0 (left end).
    stems_below: when True, CpG stems extend downward from the bar instead of upward."""
    rels = [r for r, _ in cpgs] + list(variants)
    if tss is not None:
        rels = rels + [tss]
    pad = min(rels)
    bar_x0 = x
    if bar_length is not None:
        bar_x1 = x + bar_length
    elif right_extent is not None:
        bar_x1 = x + (right_extent + pad) * w
    else:
        bar_x1 = x + (max(rels) + pad) * w
    out = [f'<rect x="{bar_x0}" y="{y - 4}" width="{bar_x1 - bar_x0}" height="8" fill="{color}"/>']
    if label is not None:
        if label_at_start:
            out.append(
                f'<text x="{bar_x0 - 10}" y="{y + 5}" text-anchor="end" '
                f'font-family="Arial, sans-serif" font-size="{label_fontsize}">{label}</text>'
            )
        else:
            out.append(
                f'<text x="{bar_x1 + 10}" y="{y + 5}" text-anchor="start" '
                f'font-family="Arial, sans-serif" font-size="{label_fontsize}">{label}</text>'
            )
    stem_h = 24
    r = 7
    sign = 1 if stems_below else -1
    bar_edge = y + sign * 4  # bottom edge if stems_below else top edge
    stem_end = bar_edge + sign * stem_h
    circle_cy = bar_edge + sign * (stem_h + r - 2)
    for rx, meth in cpgs:
        cx = x + rx * w
        out.append(line(cx, bar_edge, cx, stem_end, sw=1.2))
        out.append(circle(cx, circle_cy, r,
                          fill='black' if meth else 'white', sw=1.2))
    for rx in variants:
        cx = x + rx * w
        out.append(star(cx, y, r=12))
    if tss is not None:
        tx = x + tss * w
        bar_top = y - 4
        h_up = 30
        h_right = 22
        ay = bar_top - h_up
        out.append(line(tx, bar_top, tx, ay, sw=1.8))
        out.append(line(tx, ay, tx + h_right, ay, sw=1.8))
        ah = 5
        out.append(
            f'<polygon points="{tx + h_right + ah},{ay} '
            f'{tx + h_right - 2},{ay - ah} '
            f'{tx + h_right - 2},{ay + ah}" fill="black"/>'
        )
    return '\n'.join(out)


def person(cx, cy, kind, size=100):
    if kind == 'F':
        return circle(cx, cy, size / 2, sw=1.8)
    return rect(cx - size / 2, cy - size / 2, size, size, sw=1.8, rx=0)


# ---------- scenario ----------

@dataclass
class HapData:
    cpgs: list          # [(rel_x, methylated_bool), ...]
    color: str
    variants: tuple = ()
    label: str | None = None

@dataclass
class Scenario:
    mom_top: HapData
    mom_bot: HapData
    dad_top: HapData
    dad_bot: HapData
    kid_top: HapData      # paternal (orange)
    kid_bot: HapData      # maternal (blue)


# Shared CpG positions
RELS = [0.08, 0.18, 0.28, 0.50, 0.60]
METH5 = [(r, True) for r in RELS]
UNMETH5 = [(r, False) for r in RELS]

KID_RELS = [0.10, 0.20, 0.30, 0.55, 0.65]


# ---------- scenarios ----------

SCENARIO_DENOVO = Scenario(
    mom_top=HapData(METH5, '#b2df8a', label='C'),
    mom_bot=HapData(UNMETH5, '#fb9a99', label='D'),
    dad_top=HapData(UNMETH5, '#a6cee3', label='A'),
    dad_bot=HapData(METH5, '#fdbf6f', label='B'),
    # Kid pat = dad A; kid mat = mom C. Two de novo events both sit at
    # the right CpG cluster (positions KID_RELS[3:]). Paternal A inherits
    # dad's unmethylated state at the left cluster but is methylated at
    # the right cluster (de novo gain). Maternal C inherits mom's
    # methylated state at the left cluster but is unmethylated at the
    # right cluster (de novo loss).
    kid_top=HapData([(r, m) for r, m in zip(KID_RELS, [False, False, False, True, True])],
                    '#a6cee3', label='A'),
    kid_bot=HapData([(r, m) for r, m in zip(KID_RELS, [True, True, True, False, False])],
                    '#b2df8a', label='C'),
)

# Compound genetic-epigenetic heterozygote:
#   paternal allele in kid carries a genetic variant (star)
#   maternal allele in kid carries aberrant methylation (silenced promoter-like region)
# Each parent is a silent carrier of one of the two hits.
SCENARIO_COMPOUND = Scenario(
    # Mom: C = hyper-methylated carrier; D = normal unmethylated
    mom_top=HapData([(r, True) for r in RELS[:3]],
                    '#b2df8a', label='C'),
    mom_bot=HapData([(r, False) for r in RELS[:3]], '#fb9a99', label='D'),
    # Dad: A = genetic-variant carrier; B = normal
    dad_top=HapData([(r, False) for r in RELS[:3]], '#a6cee3',
                    variants=(0.80,), label='A'),
    dad_bot=HapData([(r, False) for r in RELS[:3]], '#fdbf6f', label='B'),
    # Kid: paternal = A (variant); maternal = C (aberrant meth) — compound het in trans
    kid_top=HapData([(r, False) for r in KID_RELS[:3]], '#a6cee3',
                    variants=(0.82,), label='A'),
    kid_bot=HapData([(r, True) for r in KID_RELS[:3]],
                   '#b2df8a', label='C'),
)


# ---------- layout ----------

def _vhap_leftmost_box(cx, bar_y, w, hap, n=3,
                       tss_rel=0.50, bar_right_rel=0.90):
    """Dashed rectangle around the n leftmost CpGs on a vertically-laid
    haplotype bar (CpGs above the bar)."""
    rels = [r for r, _ in hap.cpgs] + list(hap.variants) + [tss_rel]
    span = (bar_right_rel + min(rels)) * w
    bar_x0 = cx - span / 2
    cpg_rels = sorted(r for r, _ in hap.cpgs)[:n]
    pad_x = 12
    pad_top = 6
    pad_bot = 4
    x_left = bar_x0 + cpg_rels[0] * w - pad_x
    x_right = bar_x0 + cpg_rels[-1] * w + pad_x
    y_top = bar_y - 42 - pad_top
    y_bot = bar_y - 4 - pad_bot
    return rect(x_left, y_top, x_right - x_left, y_bot - y_top,
                fill='none', stroke='black', sw=2.0, dash='6,4', rx=0)


def _build_vertical(scenario: Scenario, label_fontsize: int,
                    swap_parents: bool, highlight_compound: bool = False) -> str:
    """Variant layout: both of each parent's haplotypes sit above the parent
    shape, and both of the kid's sit below the kid shape. Stems are NOT
    flipped — CpGs always sit above their bar (default haplotype orientation)."""
    X0v, Y0v, VWv, VHv = 335, 80, 725, 495
    head = (
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{X0v} {Y0v} {VWv} {VHv}" '
        f'font-family="Arial, sans-serif" font-size="14">'
        f'<rect x="{X0v}" y="{Y0v}" width="{VWv}" height="{VHv}" fill="white"/>'
    )
    parts = [head]

    if swap_parents:
        left_shape, right_shape = 'M', 'F'
        left_top, left_bot = scenario.dad_top, scenario.dad_bot
        right_top, right_bot = scenario.mom_top, scenario.mom_bot
    else:
        left_shape, right_shape = 'F', 'M'
        left_top, left_bot = scenario.mom_top, scenario.mom_bot
        right_top, right_bot = scenario.dad_top, scenario.dad_bot

    SHAPE = 80
    HAP_W = 280
    LEFT_X, RIGHT_X = 480, 880
    PARENT_Y = 320
    KID_X = (LEFT_X + RIGHT_X) // 2
    KID_Y = 500
    OUTER_OFFSET = 160       # parent outer hap bar offset above shape center
    INNER_OFFSET = 85        # parent inner hap bar offset above shape center
    KID_OUTER_OFFSET = 180   # kid outer hap bar offset below shape center
    KID_INNER_OFFSET = 125   # kid inner hap bar offset below shape center

    parts.append(person(LEFT_X, PARENT_Y, left_shape, size=SHAPE))
    parts.append(person(RIGHT_X, PARENT_Y, right_shape, size=SHAPE))

    TSS_REL = 0.50
    BAR_RIGHT_REL = 0.90

    def vhap(cx, y, hap):
        # Center the bar's drawn extent on cx. With right_extent fixed,
        # the bar spans (BAR_RIGHT_REL + min(rels)) * HAP_W in x.
        rels = [r for r, _ in hap.cpgs] + list(hap.variants) + [TSS_REL]
        span = (BAR_RIGHT_REL + min(rels)) * HAP_W
        return haplotype(cx - span / 2, y, HAP_W, hap.color,
                         hap.cpgs, hap.variants,
                         label=hap.label, label_fontsize=label_fontsize,
                         tss=TSS_REL, right_extent=BAR_RIGHT_REL)

    # Both parent haps above the shape: top hap higher, bot hap closer to shape.
    parts.append(vhap(LEFT_X, PARENT_Y - OUTER_OFFSET, left_top))
    parts.append(vhap(RIGHT_X, PARENT_Y - OUTER_OFFSET, right_top))
    parts.append(vhap(LEFT_X, PARENT_Y - INNER_OFFSET, left_bot))
    parts.append(vhap(RIGHT_X, PARENT_Y - INNER_OFFSET, right_bot))

    # Marriage line + line down to kid (terminating at top of kid figure).
    parts.append(line(LEFT_X + SHAPE / 2, PARENT_Y, RIGHT_X - SHAPE / 2, PARENT_Y))
    parts.append(line(KID_X, PARENT_Y, KID_X, KID_Y - SHAPE / 2))

    parts.append(person(KID_X, KID_Y, 'M', size=SHAPE))
    # Kid hap pair sits to the RIGHT of the kid shape, stacked vertically
    # (kid_top above kid_bot). Compute bar centers so that bar_x0 sits
    # just past the shape's right edge with a small gap.
    KID_HAP_GAP = 20  # gap between shape right edge and bar left edge
    kid_rels = [r for r, _ in scenario.kid_top.cpgs] \
        + list(scenario.kid_top.variants) + [TSS_REL]
    kid_span = (BAR_RIGHT_REL + min(kid_rels)) * HAP_W
    KID_HAP_CX = KID_X + SHAPE / 2 + KID_HAP_GAP + kid_span / 2
    KID_TOP_DY = -38
    KID_BOT_DY = 38
    parts.append(vhap(KID_HAP_CX, KID_Y + KID_TOP_DY, scenario.kid_top))
    parts.append(vhap(KID_HAP_CX, KID_Y + KID_BOT_DY, scenario.kid_bot))

    if highlight_compound:
        # Dashed boxes around the 3 leftmost CpGs on hap C (mom_top, the
        # outer parent bar on the mom side after swap_parents) and the
        # kid's maternal hap C (kid_bot) — marking where hap C is the
        # methylation outlier relative to A, B, D.
        parts.append(_vhap_leftmost_box(
            RIGHT_X, PARENT_Y - OUTER_OFFSET, HAP_W, right_top))
        parts.append(_vhap_leftmost_box(
            KID_HAP_CX, KID_Y + KID_BOT_DY, HAP_W, scenario.kid_bot))

    parts.append(SVG_TAIL)
    return '\n'.join(parts)


def _denovo_highlight_box(bar_x0, bar_y, w, cpgs):
    """Dashed rectangle around the two rightmost CpGs on a haplotype bar."""
    last_two = sorted(r for r, _ in cpgs)[-2:]
    pad_x = 12
    pad_top = 6
    pad_bot = 4
    x_left = bar_x0 + last_two[0] * w - pad_x
    x_right = bar_x0 + last_two[1] * w + pad_x
    # CpG circle top sits at bar_y - 4 - stem_h(24) - r(7) - r(7) ≈ bar_y - 42.
    y_top = bar_y - 42 - pad_top
    y_bot = bar_y - 4 - pad_bot
    return rect(x_left, y_top, x_right - x_left, y_bot - y_top,
                fill='none', stroke='black', sw=2.0, dash='6,4', rx=0)


def build(scenario: Scenario, label_fontsize: int = 16,
          swap_parents: bool = False,
          vertical_haps: bool = False,
          highlight_denovo: bool = False,
          highlight_compound: bool = False,
          kid_y_offset: int = 0,
          uniform_bars: bool = False) -> str:
    if vertical_haps:
        return _build_vertical(scenario, label_fontsize, swap_parents,
                               highlight_compound=highlight_compound)
    vh = VH + kid_y_offset
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{X0} {Y0} {VW} {vh}" '
        f'font-family="Arial, sans-serif" font-size="14">'
        f'<rect x="{X0}" y="{Y0}" width="{VW}" height="{vh}" fill="white"/>'
    ]

    left_pos = (560, 190)
    right_pos = (800, 190)
    if swap_parents:
        left_shape, right_shape = 'M', 'F'
        left_top, left_bot = scenario.dad_top, scenario.dad_bot
        right_top, right_bot = scenario.mom_top, scenario.mom_bot
    else:
        left_shape, right_shape = 'F', 'M'
        left_top, left_bot = scenario.mom_top, scenario.mom_bot
        right_top, right_bot = scenario.dad_top, scenario.dad_bot
    parts.append(person(*left_pos, left_shape))
    parts.append(person(*right_pos, right_shape))

    # Left-parent labels sit on the bar's tail (person side, right end);
    # right-parent labels sit on the head (person side, left end).
    # Pin all 6 haplotype bars to a single pixel length so partner haps
    # (and kid vs. parents) match visually even when only some carry a
    # variant. Length is the largest natural bar extent across all six,
    # computed per-call from each hap's own w/pad/max-rel so no CpG or
    # variant glyph is clipped.
    # In uniform_bars mode (wiki SVGs), use a single narrower w for all
    # six haplotype calls so that bars fit inside the viewBox and don't
    # overlap the parent/kid shapes, even when one hap carries a variant
    # at a large rel offset (~0.80). The unified pixel length is the
    # max natural extent across the six haps under that w.
    w_left, w_right, w_kid = (200, 200, 200) if uniform_bars else (340, 280, 280)
    def _natural_len(h, w):
        rels = [r for r, _ in h.cpgs] + list(h.variants)
        return (max(rels) + min(rels)) * w
    bar_len = max(
        _natural_len(left_top, w_left), _natural_len(left_bot, w_left),
        _natural_len(right_top, w_right), _natural_len(right_bot, w_right),
        _natural_len(scenario.kid_top, w_kid), _natural_len(scenario.kid_bot, w_kid),
    ) if uniform_bars else None
    parts.append(haplotype(240, 184, w_left, left_top.color,
                           left_top.cpgs, left_top.variants,
                           label=left_top.label,
                           label_fontsize=label_fontsize,
                           bar_length=bar_len))
    parts.append(haplotype(240, 234, w_left, left_bot.color,
                           left_bot.cpgs, left_bot.variants,
                           label=left_bot.label,
                           label_fontsize=label_fontsize,
                           bar_length=bar_len))
    parts.append(haplotype(890, 184, w_right, right_top.color,
                           right_top.cpgs, right_top.variants,
                           label=right_top.label, label_at_start=True,
                           label_fontsize=label_fontsize,
                           bar_length=bar_len))
    parts.append(haplotype(890, 234, w_right, right_bot.color,
                           right_bot.cpgs, right_bot.variants,
                           label=right_bot.label, label_at_start=True,
                           label_fontsize=label_fontsize,
                           bar_length=bar_len))

    parts.append(line(left_pos[0] + 50, left_pos[1], right_pos[0] - 50, right_pos[1]))
    midx = (left_pos[0] + right_pos[0]) / 2
    son = (int(midx), 400 + kid_y_offset)
    parts.append(line(midx, left_pos[1], son[0], son[1] - 50))
    parts.append(person(*son, 'M'))

    parts.append(haplotype(380, 394 + kid_y_offset, w_kid, scenario.kid_top.color,
                           scenario.kid_top.cpgs, scenario.kid_top.variants,
                           label=scenario.kid_top.label,
                           label_fontsize=label_fontsize,
                           bar_length=bar_len))
    parts.append(haplotype(380, 444 + kid_y_offset, w_kid, scenario.kid_bot.color,
                           scenario.kid_bot.cpgs, scenario.kid_bot.variants,
                           label=scenario.kid_bot.label,
                           label_fontsize=label_fontsize,
                           bar_length=bar_len))

    if highlight_denovo:
        # Dashed boxes around the two rightmost CpGs on dad's hap A (top
        # parent bar on the dad side) and on the kid's paternal hap A —
        # marking where the de novo gain of methylation arises.
        parts.append(_denovo_highlight_box(240, 184, 340, left_top.cpgs))
        parts.append(_denovo_highlight_box(380, 394 + kid_y_offset, 280, scenario.kid_top.cpgs))

    parts.append(SVG_TAIL)
    return '\n'.join(parts)


# ---------- matplotlib companion: BED snippet + polars query ----------

# Shared data for the de novo BED panel (methylation table + polars query).
# CpG genomic coordinates mirror KID_RELS in the trio_denovo SVG: a left
# cluster of 3 CpGs (concordant unmethylated on dad_A → kid_pat) and a
# right cluster of 2 CpGs (de novo gain — kid_pat methylated where dad_A
# is not). dad_B is the non-transmitted, fully methylated homolog.
# Methylation levels jittered around the SVG values for realism.
DENOVO_TABLE_TITLE = "Haplotype-specific methylation levels"
DENOVO_TABLE_COLS = ("chrom", "start", "kid_pat", "pat_hap", "dad_A")
DENOVO_TABLE_ROWS = [
    ("chr1", "1400", "0.07", "A", "0.05"),
    ("chr1", "1450", "0.06", "A", "0.04"),
    ("chr1", "1550", "0.94", "A", "0.05"),
    ("chr1", "1650", "0.96", "A", "0.04"),
    ("chr1", "1700", "0.07", "A", "0.05"),
    ("chr1", "1800", "0.05", "A", "0.06"),
]
DENOVO_TABLE_HIGHLIGHT = {2, 3}
DENOVO_CODE_TITLE = "Discover de novo gain of methylation on paternal haplotype"
DENOVO_CODE = (
    '# Pick the dad haplotype that was transmitted to the kid.\n'
    'df_meth = df_meth.with_columns(\n'
    '    dad_transmitted = pl.when(pl.col("pat_hap") == "A")\n'
    '        .then(pl.col("dad_A"))\n'
    '        .otherwise(pl.col("dad_B"))\n'
    ')\n'
    '\n'
    '# Find runs of consecutive CpGs where the kid_pat is methylated\n'
    '# but the dad-transmitted haplotype is not — i.e., a de novo\n'
    '# *gain* of methylation on the kid paternal homolog.\n'
    'WINDOW = 2\n'
    'df_denovo_gain_pat = (\n'
    '    df_meth.with_columns(\n'
    '        delta = pl.col("kid_pat") - pl.col("dad_transmitted"),\n'
    '    )\n'
    '    .with_columns(\n'
    '        mean_delta = pl.col("delta").rolling_mean(WINDOW),\n'
    '    )\n'
    '    .filter(pl.col("mean_delta") > 0.5)\n'
    ')\n'
)


def render_trio_denovo_meth_table(out_path: Path | None = None,
                                   show_title: bool = True,
                                   fig_width: float = 8.0,
                                   col_dx: float = 0.085,
                                   cell_fontsize: float = 10.5,
                                   trim_whitespace: bool = False):
    """Methylation table only (no polars-query subpanel)."""
    fig_h = 2.6 if show_title else 2.2
    fig = plt.figure(figsize=(fig_width, fig_h))
    ax = fig.add_axes([0.04, 0.05, 0.94, 0.9])
    _draw_simple_table(
        ax, DENOVO_TABLE_COLS, DENOVO_TABLE_ROWS,
        title=DENOVO_TABLE_TITLE if show_title else "",
        highlight_rows=DENOVO_TABLE_HIGHLIGHT,
        bottom_rule=False,
        col_dx=col_dx,
        cell_fontsize=cell_fontsize,
        highlight_color="#a6cee3",  # haplotype A color (where the de novo arises)
    )
    out = out_path if out_path is not None else OUT / "trio_denovo_meth_table.png"
    _save_fig(fig, out, trim_whitespace=trim_whitespace)
    plt.close(fig)
    print(f"[scratch] Wrote {out}")


def _save_fig(fig, out, trim_whitespace: bool = False) -> None:
    if trim_whitespace and Path(out).suffix.lower() == ".pdf":
        # Manual bbox union over visible children — gives a tighter crop
        # on the right and bottom than bbox_inches='tight', which uses
        # axes-level bboxes that include empty axes area. Pad the bottom
        # by half a line height to compensate for matplotlib's
        # Text.get_window_extent() under-reporting the height of long
        # multi-line code blocks (otherwise the last 1-3 lines clip).
        from matplotlib.transforms import Bbox
        fig.canvas.draw()
        r = fig.canvas.get_renderer()
        from matplotlib.spines import Spine
        from matplotlib.axis import XAxis, YAxis
        from matplotlib.text import Text
        bboxes = []
        max_text_h_px = 0.0
        for ax in fig.axes:
            for child in ax.get_children():
                if child is ax.patch or not child.get_visible():
                    continue
                if isinstance(child, (Spine, XAxis, YAxis)):
                    continue
                if isinstance(child, Text) and not child.get_text():
                    continue
                try:
                    b = child.get_window_extent(r)
                except Exception:
                    continue
                if b.width > 0 and b.height > 0:
                    bboxes.append(b)
                    if isinstance(child, Text) and "\n" in child.get_text():
                        max_text_h_px = max(max_text_h_px, b.height)
        u = Bbox.union(bboxes)
        # Pad the bottom by ~5% of the tallest multi-line text's bbox to
        # absorb matplotlib's under-counted height on long snippets.
        # Scales with text length so short panels keep tight bottoms.
        u = Bbox.from_extents(
            u.x0, u.y0 - 0.05 * max_text_h_px, u.x1, u.y1)
        bbox_in = u.transformed(fig.dpi_scale_trans.inverted())
        fig.savefig(out, dpi=180, bbox_inches=bbox_in, pad_inches=0.02)
    else:
        fig.savefig(out, dpi=180)
        if trim_whitespace:
            _trim_whitespace_png(out)


def _trim_whitespace_png(path) -> None:
    """Crop pure-white margins off all four sides of a PNG in place."""
    from PIL import Image, ImageChops
    img = Image.open(path).convert("RGB")
    bg = Image.new("RGB", img.size, (255, 255, 255))
    diff = ImageChops.difference(img, bg)
    bbox = diff.getbbox()
    if bbox is not None:
        # Add a small breathing-room margin (10 px) on all sides.
        pad = 10
        left = max(0, bbox[0] - pad)
        top = max(0, bbox[1] - pad)
        right = min(img.size[0], bbox[2] + pad)
        bottom = min(img.size[1], bbox[3] + pad)
        img.crop((left, top, right, bottom)).save(path)


def render_trio_denovo_polars_code(out_path: Path | None = None,
                                    show_title: bool = True,
                                    trim_whitespace: bool = False):
    """Polars discovery code snippet only (no table subpanel)."""
    fig_h = 5.5 if show_title else 5.0
    fig = plt.figure(figsize=(12.0, fig_h))
    ax = fig.add_axes([0.04, 0.04, 0.94, 0.92])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_axis_off()
    ax.patch.set_visible(False)

    if show_title:
        ax.text(
            0.0, 0.97, DENOVO_CODE_TITLE,
            ha="left", va="top", fontsize=12, fontweight="bold",
            transform=ax.transAxes,
        )
        code_y = 0.85
    else:
        code_y = 0.97
    ax.text(
        0.0, code_y, DENOVO_CODE,
        ha="left", va="top", fontsize=10.5,
        family="Menlo", color=COLOR_NEUTRAL,
        transform=ax.transAxes,
    )

    out = out_path if out_path is not None else OUT / "trio_denovo_polars_code.png"
    _save_fig(fig, out, trim_whitespace=trim_whitespace)
    plt.close(fig)
    print(f"[scratch] Wrote {out}")


def render_trio_denovo_bed(out_path: Path | None = None):
    """BED snippet + polars query for discovering de novo epimutations."""
    FIG_W, FIG_H = 12.0, 7.5
    fig = plt.figure(figsize=(FIG_W, FIG_H))
    LEFT_IN, RIGHT_IN = 0.04 * FIG_W, 0.98 * FIG_W
    TOP_IN = FIG_H - 0.06 * FIG_H
    BOTTOM_IN = 0.04 * FIG_H

    TABLE_H_IN = 2.04
    TOP_EMPTY_FRAC = 0.03
    META_CONTENT_FRAC = 1.0
    VISIBLE_GAP_IN = 0.4

    meth_empty_in = TABLE_H_IN * (1 - META_CONTENT_FRAC)
    table_top_empty_in = TABLE_H_IN * TOP_EMPTY_FRAC
    gap_table_code_in = VISIBLE_GAP_IN - meth_empty_in - table_top_empty_in
    code_h_in = (TOP_IN - BOTTOM_IN) - TABLE_H_IN - gap_table_code_in
    width_in = RIGHT_IN - LEFT_IN

    def add_axes_in(left_in, bottom_in, w_in, h_in):
        return fig.add_axes([
            left_in / FIG_W, bottom_in / FIG_H,
            w_in / FIG_W,    h_in / FIG_H,
        ])

    ax_meth_top_in = TOP_IN
    ax_meth = add_axes_in(LEFT_IN, ax_meth_top_in - TABLE_H_IN, width_in, TABLE_H_IN)
    ax_code_top_in = ax_meth_top_in - TABLE_H_IN - gap_table_code_in
    ax_code = add_axes_in(LEFT_IN, ax_code_top_in - code_h_in, width_in, code_h_in)

    _draw_simple_table(
        ax_meth, DENOVO_TABLE_COLS, DENOVO_TABLE_ROWS,
        title=DENOVO_TABLE_TITLE,
        highlight_rows=DENOVO_TABLE_HIGHLIGHT,
        bottom_rule=False,
    )

    ax_code.set_xlim(0, 1)
    ax_code.set_ylim(0, 1)
    ax_code.set_axis_off()
    CODE_TITLE_Y = 0.97
    CODE_TITLE_TO_CODE_IN = (0.97 - 0.72) * TABLE_H_IN
    code_y = CODE_TITLE_Y - CODE_TITLE_TO_CODE_IN / code_h_in
    ax_code.text(
        0.0, CODE_TITLE_Y, DENOVO_CODE_TITLE,
        ha="left", va="top", fontsize=12, fontweight="bold",
        transform=ax_code.transAxes,
    )
    ax_code.text(
        0.0, code_y, DENOVO_CODE,
        ha="left", va="top", fontsize=10.5,
        family="Menlo", color=COLOR_NEUTRAL,
        transform=ax_code.transAxes,
    )

    out = out_path if out_path is not None else OUT / "trio_denovo_bed.png"
    fig.savefig(out, dpi=180)
    plt.close(fig)
    print(f"[scratch] Wrote {out}")


def _draw_simple_table(ax, cols, rows, title, highlight_rows=frozenset(),
                        bottom_rule=True, col_dx=0.085, cell_fontsize=10.5,
                        highlight_color="#e5e5e5"):
    """Render a small fixed-width table on a transAxes axis."""
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_axis_off()
    ax.patch.set_visible(False)
    if title:
        ax.text(
            0.0, 0.97, title,
            ha="left", va="top", fontsize=12, fontweight="bold",
            transform=ax.transAxes,
        )
    n_col = len(cols)
    col_xs = [0.02 + col_dx * i for i in range(n_col)]
    table_right = col_xs[-1] + 0.07

    header_y = 0.72
    for x, name in zip(col_xs, cols):
        ax.text(
            x, header_y, name,
            ha="left", va="center", fontsize=cell_fontsize,
            family="monospace", color=COLOR_NEUTRAL,
            transform=ax.transAxes,
        )
    ax.plot(
        [0.005, table_right], [header_y - 0.06, header_y - 0.06],
        color=COLOR_NEUTRAL, linewidth=0.7, transform=ax.transAxes,
    )

    row_dy = 0.09
    for r_idx, row in enumerate(rows):
        y = header_y - 0.12 - r_idx * row_dy
        if r_idx in highlight_rows:
            ax.add_patch(Rectangle(
                (0.005, y - row_dy / 2 + 0.01),
                table_right - 0.005, row_dy - 0.02,
                facecolor=highlight_color, edgecolor="none",
                transform=ax.transAxes,
            ))
        for x, cell in zip(col_xs, row):
            ax.text(
                x, y, cell,
                ha="left", va="center", fontsize=cell_fontsize,
                family="monospace", color=COLOR_NEUTRAL,
                transform=ax.transAxes,
            )
    if bottom_rule:
        rule_y = header_y - 0.12 - (len(rows) - 1) * row_dy - row_dy / 2 + 0.01
        ax.plot(
            [0.005, table_right], [rule_y, rule_y],
            color=COLOR_NEUTRAL, linewidth=0.7, transform=ax.transAxes,
        )


# Shared data for the compound genetic-epigenetic heterozygote BED panel.
COMPOUND_METH_TITLE = "Haplotype-specific methylation levels"
COMPOUND_METH_COLS = ("chrom", "start",
                      "kid_mat", "mat_hap",
                      "mom_C", "mom_D")
COMPOUND_METH_ROWS = [
    ("chr1", "0900", "0.06", "C", "0.07", "0.05"),
    ("chr1", "1000", "0.05", "C", "0.06", "0.04"),
    ("chr1", "1100", "0.94", "C", "0.93", "0.06"),
    ("chr1", "1200", "0.96", "C", "0.95", "0.05"),
    ("chr1", "1300", "0.91", "C", "0.92", "0.04"),
    ("chr1", "1400", "0.04", "C", "0.05", "0.06"),
    ("chr1", "1500", "0.06", "C", "0.04", "0.05"),
]
COMPOUND_METH_HIGHLIGHT = {2, 3, 4}

COMPOUND_GENO_TITLE = "Phased genotypes"
COMPOUND_GENO_COLS = ("chrom", "pos",
                      "kid_pat", "dad_A", "dad_B")
COMPOUND_GENO_ROWS = [
    ("chr1", "0500", "0", "0", "1"),
    ("chr1", "0800", "0", "0", "1"),
    ("chr1", "1420", "1", "1", "0"),
    ("chr1", "2000", "0", "0", "1"),
    ("chr1", "2300", "0", "0", "1"),
]
COMPOUND_GENO_HIGHLIGHT = {2}

COMPOUND_CODE_TITLE = "Discover compound genetic-epigenetic heterozygous locus"
COMPOUND_CODE = (
    '# --- Step 1 -----------------------------------------------------\n'
    '# Find runs of consecutive CpGs where the same parental hap is the\n'
    '# meth outlier, was transmitted to the kid, and the kid hap on that\n'
    '# side has ~the same meth level.\n'
    'TOL = 0.10\n'
    'df_meth_hits = (df_meth\n'
    '    # sites where exactly one parental hap is the meth outlier\n'
    '    .with_columns(\n'
    '        n_hi         = ((pl.col("dad_A") > 0.5)\n'
    '                      + (pl.col("dad_B") > 0.5)\n'
    '                      + (pl.col("mom_C") > 0.5)\n'
    '                      + (pl.col("mom_D") > 0.5)),\n'
    '        outlier_hap  = pl.when(pl.col("dad_A") > 0.5).then(pl.lit("A"))\n'
    '                         .when(pl.col("dad_B") > 0.5).then(pl.lit("B"))\n'
    '                         .when(pl.col("mom_C") > 0.5).then(pl.lit("C"))\n'
    '                         .otherwise(pl.lit("D")),\n'
    '        outlier_meth = pl.max_horizontal("dad_A", "dad_B", "mom_C", "mom_D"),\n'
    '    )\n'
    '    .filter(pl.col("n_hi") == 1)\n'
    '    # group consecutive CpGs sharing the same outlier hap\n'
    '    .with_columns(\n'
    '        run_id = (pl.col("outlier_hap")\n'
    '                  != pl.col("outlier_hap").shift(1)).cum_sum()\n'
    '    )\n'
    '    # collapse each run; outlier_meth = mean over its CpGs\n'
    '    .group_by("run_id", maintain_order=True)\n'
    '    .agg(\n'
    '        chrom        = pl.col("chrom").first(),\n'
    '        start        = pl.col("start").min(),\n'
    '        end          = pl.col("start").max(),\n'
    '        n            = pl.len(),\n'
    '        outlier_hap  = pl.col("outlier_hap").first(),\n'
    '        outlier_meth = pl.col("outlier_meth").mean(),\n'
    '        pat_hap      = pl.col("pat_hap").first(),\n'
    '        mat_hap      = pl.col("mat_hap").first(),\n'
    '        kid_pat      = pl.col("kid_pat").mean(),\n'
    '        kid_mat      = pl.col("kid_mat").mean(),\n'
    '    )\n'
    '    .filter(pl.col("n") >= 3)\n'
    '    # kid inherits the outlier hap with ~unchanged meth\n'
    '    .filter(\n'
    '        ((pl.col("outlier_hap") == pl.col("pat_hap"))\n'
    '         & ((pl.col("kid_pat") - pl.col("outlier_meth")).abs() < TOL))\n'
    '        | ((pl.col("outlier_hap") == pl.col("mat_hap"))\n'
    '           & ((pl.col("kid_mat") - pl.col("outlier_meth")).abs() < TOL))\n'
    '    )\n'
    ')\n'
    '\n'
    '# --- Step 2 -----------------------------------------------------\n'
    '# Exactly one parental hap carries the ALT; kid inherits it on the\n'
    '# matching side.\n'
    'df_geno_hits = (df_geno\n'
    '    # tag each SNV with the ALT count and the (unique) outlier hap\n'
    '    .with_columns(\n'
    '        n_alt       = pl.sum_horizontal("dad_A", "dad_B", "mom_C", "mom_D"),\n'
    '        outlier_hap = pl.when(pl.col("dad_A") == 1).then(pl.lit("A"))\n'
    '                        .when(pl.col("dad_B") == 1).then(pl.lit("B"))\n'
    '                        .when(pl.col("mom_C") == 1).then(pl.lit("C"))\n'
    '                        .otherwise(pl.lit("D")),\n'
    '    )\n'
    '    # keep SNVs with a single ALT-carrying parental hap\n'
    '    .filter(pl.col("n_alt") == 1)\n'
    '    # kid inherits that ALT on the matching side\n'
    '    .filter(\n'
    '        (pl.col("outlier_hap").is_in(["A", "B"]) & (pl.col("kid_pat") == 1))\n'
    '        | (pl.col("outlier_hap").is_in(["C", "D"]) & (pl.col("kid_mat") == 1))\n'
    '    )\n'
    ')\n'
    '\n'
    '# --- Step 3 -----------------------------------------------------\n'
    '# Compound het = meth-outlier and geno-outlier from different\n'
    '# parents, same locus.\n'
    'df_compound_het = (df_meth_hits\n'
    '    .rename({"outlier_hap": "outlier_hap_meth"})\n'
    '    # pair each meth-outlier region with the nearest geno-outlier\n'
    '    # SNV on the same chromosome, within 500 bp\n'
    '    .join_asof(\n'
    '        df_geno_hits.rename({"outlier_hap": "outlier_hap_geno"}),\n'
    '        on="start", by="chrom",\n'
    '        strategy="nearest", tolerance=500,\n'
    '    )\n'
    '    # keep loci where the meth and geno outliers are in trans\n'
    '    # (i.e., come from different parents)\n'
    '    .filter(\n'
    '        (pl.col("outlier_hap_meth").is_in(["A", "B"])\n'
    '         & pl.col("outlier_hap_geno").is_in(["C", "D"]))\n'
    '        | (pl.col("outlier_hap_meth").is_in(["C", "D"])\n'
    '           & pl.col("outlier_hap_geno").is_in(["A", "B"]))\n'
    '    )\n'
    ')\n'
)


def render_trio_compound_het_tables(out_path: Path | None = None,
                                     show_title: bool = True,
                                     fig_width: float = 12.0,
                                     trim_whitespace: bool = False):
    """Methylation + phased-genotype tables only (no polars-query subpanel)."""
    TABLE_H_IN = 2.04
    # _draw_simple_table places the header at axes-y 0.72 and rows below it
    # at row_dy=0.09. The bottom of the last data row therefore sits at
    # axes-y = 0.565 - (n-1)*0.09; convert to a fraction of TABLE_H_IN.
    n_meth = len(COMPOUND_METH_ROWS)
    METH_CONTENT_BOTTOM_FRAC = max(0.025, 0.565 - (n_meth - 1) * 0.09)
    GENO_CONTENT_TOP_FRAC = 0.78
    VISIBLE_GAP_IN = 0.15
    # Grow the figure if the meth content extends below the original 5-row
    # budget, so the geno table still has room beneath it.
    extra_meth_in = max(0.0, (0.205 - METH_CONTENT_BOTTOM_FRAC)) * TABLE_H_IN
    FIG_W = fig_width
    FIG_H = 5.0 + extra_meth_in
    fig = plt.figure(figsize=(FIG_W, FIG_H))
    LEFT_IN, RIGHT_IN = 0.04 * FIG_W, 0.98 * FIG_W
    TOP_IN = FIG_H - 0.06 * FIG_H
    BOTTOM_IN = 0.04 * FIG_H
    width_in = RIGHT_IN - LEFT_IN

    def add_axes_in(left_in, bottom_in, w_in, h_in):
        return fig.add_axes([
            left_in / FIG_W, bottom_in / FIG_H,
            w_in / FIG_W,    h_in / FIG_H,
        ])

    ax_meth_bottom_in = TOP_IN - TABLE_H_IN
    ax_meth = add_axes_in(LEFT_IN, ax_meth_bottom_in, width_in, TABLE_H_IN)
    # Overlap geno axis into meth's empty bottom so the *visible* content gap
    # equals VISIBLE_GAP_IN, not the axis-edge separation.
    meth_content_bottom_in = ax_meth_bottom_in + METH_CONTENT_BOTTOM_FRAC * TABLE_H_IN
    geno_content_top_in = meth_content_bottom_in - VISIBLE_GAP_IN
    ax_geno_bottom_in = geno_content_top_in - GENO_CONTENT_TOP_FRAC * TABLE_H_IN
    ax_geno = add_axes_in(LEFT_IN, ax_geno_bottom_in, width_in, TABLE_H_IN)

    _draw_simple_table(
        ax_meth, COMPOUND_METH_COLS, COMPOUND_METH_ROWS,
        title=COMPOUND_METH_TITLE if show_title else "",
        highlight_rows=COMPOUND_METH_HIGHLIGHT,
        highlight_color="#b2df8a",
        bottom_rule=False,
    )
    _draw_simple_table(
        ax_geno, COMPOUND_GENO_COLS, COMPOUND_GENO_ROWS,
        title=COMPOUND_GENO_TITLE if show_title else "",
        highlight_rows=COMPOUND_GENO_HIGHLIGHT,
        highlight_color="#a6cee3",
        bottom_rule=False,
    )

    out = out_path if out_path is not None else OUT / "trio_compound_het_tables.png"
    _save_fig(fig, out, trim_whitespace=trim_whitespace)
    plt.close(fig)
    print(f"[scratch] Wrote {out}")


# Rebucket panel: shows the per-CpG methylation values of one kid first
# under the arbitrary hap1/hap2 partition produced by hiphase, then under
# the founder-aware (pat / mat / founder-letter) labels produced by Step
# 3's relabel — the literal "relabel, don't recompute" message of the
# wiki section "Relabelling per-CpG methylation". Same numerical values
# in both sub-tables; only the column names change.
REBUCKET_TOP_COLS = ("chrom", "start", "Kid1_hap1", "Kid1_hap2")
REBUCKET_BOTTOM_COLS = ("chrom", "start",
                         "Kid1_pat", "pat_hap",
                         "Kid1_mat", "mat_hap")
REBUCKET_ROWS_TOP = [
    ("chr1", "1100", "0.05", "0.93"),
]
REBUCKET_ROWS_BOTTOM = [
    ("chr1", "1100", "0.05", "A", "0.93", "C"),
]


def render_rebucket_panel(out_path: Path | None = None,
                          show_title: bool = True,
                          fig_width: float = 8.0,
                          col_dx: float = 0.15,
                          cell_fontsize: float = 10.5,
                          trim_whitespace: bool = False) -> None:
    """Two stacked sub-tables joined by a downward arrow: same per-CpG
    methylation values, first under hiphase's arbitrary hap1/hap2 labels,
    then under founder-aware (pat/mat + founder-letter) labels.

    Used as the "rebucketing" panel in Fig 3 of the manuscript; corresponds
    to the wiki section
    `wiki/pedigree_wise_workflow/founder_phased_methylation/founder_phased_methylation.md`
    "Relabelling per-CpG methylation".
    """
    fig_h = 4.0
    fig = plt.figure(figsize=(fig_width, fig_h))
    LEFT = 0.04
    WIDTH = 0.94
    TOP_H = 0.42
    BOTTOM_H = 0.42
    TOP_BOTTOM = 1.0 - TOP_H - 0.02       # top subtable bottom edge
    BOT_BOTTOM = 0.04                      # bottom subtable bottom edge
    ax_top = fig.add_axes([LEFT, TOP_BOTTOM, WIDTH, TOP_H])
    ax_bot = fig.add_axes([LEFT, BOT_BOTTOM, WIDTH, BOTTOM_H])

    title_top = "Per-CpG methylation under hiphase hap1/hap2 partition" if show_title else ""
    title_bot = "After relabelling to founder-aware labels (pat/mat + founder letter)" if show_title else ""

    _draw_simple_table(
        ax_top, REBUCKET_TOP_COLS, REBUCKET_ROWS_TOP,
        title=title_top,
        bottom_rule=False,
        col_dx=col_dx,
        cell_fontsize=cell_fontsize,
    )
    _draw_simple_table(
        ax_bot, REBUCKET_BOTTOM_COLS, REBUCKET_ROWS_BOTTOM,
        title=title_bot,
        bottom_rule=False,
        col_dx=col_dx,
        cell_fontsize=cell_fontsize,
    )

    # Downward arrow between the two sub-tables (figure-level).
    arrow_x = LEFT + 0.10
    arrow_y_tail = TOP_BOTTOM - 0.005
    arrow_y_head = BOT_BOTTOM + BOTTOM_H + 0.005
    fig.patches.append(
        plt.matplotlib.patches.FancyArrow(
            arrow_x, arrow_y_tail,
            0.0, arrow_y_head - arrow_y_tail,
            width=0.002, head_width=0.012, head_length=0.012,
            length_includes_head=True,
            color=COLOR_NEUTRAL,
            transform=fig.transFigure,
            figure=fig,
        )
    )
    fig.text(
        arrow_x + 0.02, (arrow_y_tail + arrow_y_head) / 2,
        "relabel using bit-vector match (Fig 3B)\n"
        "+ founder-letter map (Fig 3A)",
        ha="left", va="center", fontsize=cell_fontsize - 1,
        color=COLOR_NEUTRAL, style="italic",
    )

    out = out_path if out_path is not None else OUT / "trio_rebucket_panel.png"
    _save_fig(fig, out, trim_whitespace=trim_whitespace)
    plt.close(fig)
    print(f"[scratch] Wrote {out}")


def render_trio_compound_het_polars_code(out_path: Path | None = None,
                                          show_title: bool = True,
                                          trim_whitespace: bool = False):
    """Polars discovery code snippet only (no table subpanel)."""
    # ~95 lines of code at fontsize=10.5 need ~19in of vertical space.
    n_lines = COMPOUND_CODE.count("\n") + 1
    code_h_in = n_lines * 0.20
    fig_h = code_h_in + (0.6 if show_title else 0.3)
    fig = plt.figure(figsize=(12.0, fig_h))
    ax = fig.add_axes([0.04, 0.02, 0.94, 0.96])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_axis_off()
    ax.patch.set_visible(False)

    if show_title:
        ax.text(
            0.0, 0.99, COMPOUND_CODE_TITLE,
            ha="left", va="top", fontsize=12, fontweight="bold",
            transform=ax.transAxes,
        )
        code_y = 0.95
    else:
        code_y = 0.99
    ax.text(
        0.0, code_y, COMPOUND_CODE,
        ha="left", va="top", fontsize=10.5,
        family="Menlo", color=COLOR_NEUTRAL,
        transform=ax.transAxes,
    )

    out = out_path if out_path is not None else OUT / "trio_compound_het_polars_code.png"
    _save_fig(fig, out, trim_whitespace=trim_whitespace)
    plt.close(fig)
    print(f"[scratch] Wrote {out}")


def render_trio_compound_het_bed(out_path: Path | None = None):
    """BED snippet (methylation + genotype) + polars query for compound
    genetic-epigenetic heterozygotes."""
    FIG_W, FIG_H = 12.0, 22.0
    fig = plt.figure(figsize=(FIG_W, FIG_H))
    # Manual positions in absolute inches (independent of figure size).
    # Both table axes have the same height so title→header→hline spacings
    # match. Table 2 has fewer rows than table 1, so its axis has empty
    # space at the bottom; we widen the gap above table 2 to compensate so
    # the *visible* content-to-content gaps are equal absolute distances.
    LEFT_IN, RIGHT_IN = 0.04 * FIG_W, 0.98 * FIG_W
    TOP_IN = FIG_H - 0.06 * FIG_H
    BOTTOM_IN = 0.04 * FIG_H

    TABLE_H_IN = 2.04         # absolute height of each table axis
    TOP_EMPTY_FRAC = 0.03     # empty strip above title (axes fraction)
    META_CONTENT_FRAC = 1.0   # table 1 content fills the axis
    GENO_CONTENT_FRAC = 0.53  # table 2 content stops well above the bottom
    VISIBLE_GAP_IN = 0.4      # target visible content-to-content gap (in)

    meth_empty_in = TABLE_H_IN * (1 - META_CONTENT_FRAC)
    geno_empty_in = TABLE_H_IN * (1 - GENO_CONTENT_FRAC)
    table_top_empty_in = TABLE_H_IN * TOP_EMPTY_FRAC

    gap_meth_geno_in = VISIBLE_GAP_IN - meth_empty_in - table_top_empty_in
    # For the gap below table 2, we don't yet know the code axis height;
    # but table 2's top-empty doesn't affect this gap — only its empty
    # bottom and the code axis's top-empty (which we approximate the same).
    gap_geno_code_in = VISIBLE_GAP_IN - geno_empty_in - table_top_empty_in

    code_h_in = (TOP_IN - BOTTOM_IN) - 2 * TABLE_H_IN - gap_meth_geno_in - gap_geno_code_in
    width_in = RIGHT_IN - LEFT_IN

    def add_axes_in(left_in, bottom_in, w_in, h_in):
        return fig.add_axes([
            left_in / FIG_W, bottom_in / FIG_H,
            w_in / FIG_W,    h_in / FIG_H,
        ])

    ax_meth_top_in = TOP_IN
    ax_meth = add_axes_in(LEFT_IN, ax_meth_top_in - TABLE_H_IN, width_in, TABLE_H_IN)
    ax_geno_top_in = ax_meth_top_in - TABLE_H_IN - gap_meth_geno_in
    ax_geno = add_axes_in(LEFT_IN, ax_geno_top_in - TABLE_H_IN, width_in, TABLE_H_IN)
    ax_code_top_in = ax_geno_top_in - TABLE_H_IN - gap_geno_code_in
    ax_code = add_axes_in(LEFT_IN, ax_code_top_in - code_h_in, width_in, code_h_in)

    # Methylation at CpGs in the locus: all 4 parental haps + both kid haps.
    # Outside the region every hap is unmethylated; inside it, exactly one
    # parental hap (mom_C) is hyper-methylated, and the kid inherits that
    # pattern on its maternal side.
    _draw_simple_table(
        ax_meth, COMPOUND_METH_COLS, COMPOUND_METH_ROWS,
        title=COMPOUND_METH_TITLE,
        highlight_rows=COMPOUND_METH_HIGHLIGHT,
        bottom_rule=False,
    )

    # Genotype at a SNV inside the same locus: exactly one parental hap
    # (dad_A) carries the ALT allele, and the kid inherits it on the
    # paternal side — in trans with the meth outlier.
    _draw_simple_table(
        ax_geno, COMPOUND_GENO_COLS, COMPOUND_GENO_ROWS,
        title=COMPOUND_GENO_TITLE,
        highlight_rows=COMPOUND_GENO_HIGHLIGHT,
        bottom_rule=False,
    )

    ax_code.set_xlim(0, 1)
    ax_code.set_ylim(0, 1)
    ax_code.set_axis_off()
    CODE_TITLE_Y = 0.97
    # Match the table's title→header-row gap: title at y=0.97, header at y=0.72.
    CODE_TITLE_TO_CODE_IN = (0.97 - 0.72) * TABLE_H_IN
    code_y = CODE_TITLE_Y - CODE_TITLE_TO_CODE_IN / code_h_in
    ax_code.text(
        0.0, CODE_TITLE_Y,
        COMPOUND_CODE_TITLE,
        ha="left", va="top", fontsize=12, fontweight="bold",
        transform=ax_code.transAxes,
    )
    ax_code.text(
        0.0, code_y, COMPOUND_CODE,
        ha="left", va="top", fontsize=10.5,
        family="Menlo", color=COLOR_NEUTRAL,
        transform=ax_code.transAxes,
    )

    out = out_path if out_path is not None else OUT / "trio_compound_het_bed.png"
    fig.savefig(out, dpi=180)
    plt.close(fig)
    print(f"[scratch] Wrote {out}")


def render_rebucket_visual(out_path: Path | None = None,
                           show_title: bool = False,
                           trim_whitespace: bool = False) -> None:
    """Visual rebucketing: same per-CpG methylation bars shown first
    under hiphase's arbitrary hap1/hap2 labels (neutral grey, two
    read-backed phase blocks inside one IBD segment, with the
    hap1↔hap2 assignment flipping between the two phase blocks), then
    under founder-aware labels (pat_A meth / mat_C meth). Each bar in
    the bottom view is the SAME value as in the top view, just routed
    to the founder-aware track determined by the per-hap-map-block
    matching exercise from Fig 4A.
    """
    N_SITES = 14
    RBP1 = (0, 5)
    RBP2 = (8, 13)

    # Canonical (founder-truth) methylation profiles. Identical biology
    # being measured at every CpG; only the *labels* differ across the
    # transformation. A is hypermethylated and C is hypomethylated at
    # this locus (the founder-haplotype analogue of Fig 1B's
    # discordant-region asymmetry), so the hap1↔hap2 flip across phase
    # blocks is visually obvious in the top "before-rebucketing" view.
    pat_A_meth = [0.92, 0.85, 0.95, 0.88, 0.95, 0.85, 0.0, 0.0,
                  0.88, 0.95, 0.83, 0.92, 0.95, 0.85]
    mat_C_meth = [0.10, 0.15, 0.05, 0.12, 0.05, 0.15, 0.0, 0.0,
                  0.12, 0.05, 0.17, 0.08, 0.05, 0.15]

    def slice_(profile, lo, hi):
        return [(i, profile[i]) for i in range(lo, hi + 1)]

    # In RBP1 hiphase happens to call hap1 = pat_A, hap2 = mat_C.
    # In RBP2 the labels flip: hap1 = mat_C, hap2 = pat_A. Phase blocks
    # are smaller than IBD segments, so the hap1/hap2 assignment can
    # readily flip between neighbouring blocks within the same IBD.
    hap1_bars = slice_(pat_A_meth, *RBP1) + slice_(mat_C_meth, *RBP2)
    hap2_bars = slice_(mat_C_meth, *RBP1) + slice_(pat_A_meth, *RBP2)
    patA_bars = slice_(pat_A_meth, *RBP1) + slice_(pat_A_meth, *RBP2)
    matC_bars = slice_(mat_C_meth, *RBP1) + slice_(mat_C_meth, *RBP2)

    import sys as _sys
    _repo_root = Path(__file__).resolve().parents[2]
    if str(_repo_root) not in _sys.path:
        _sys.path.insert(0, str(_repo_root))
    from wiki._panel_grid import (  # noqa: E402
        HAP_PALETTE, draw_two_stripe_block,
    )

    NEUTRAL_TOP = HAP_PALETTE["neutral_top"]
    NEUTRAL_BOT = HAP_PALETTE["neutral_bot"]
    A_COLOR = HAP_PALETTE["A"]
    C_COLOR = HAP_PALETTE["C"]
    BAR_W = 0.7
    PREFIX_X = -6.5
    # Match Fig 4A: monospace 11 for row labels and the short RBP-style
    # in-block label; 9 for the longer hap-map conclusion that has to
    # fit inside a small rectangle.
    LABEL_FS = 11
    BLOCK_FS_SHORT = 11
    BLOCK_FS_LONG = 9
    LABEL_FAMILY = "monospace"

    fig = plt.figure(figsize=(8.5, 6.0))
    # Layout (top → bottom). Block rows use the same row pitch as Fig 4A
    # (≈0.30 in) so the rectangles render at the same physical size; bar
    # rows are taller because they need vertical room for the bigwig-
    # style methylation bars.
    #   0  RBP rectangles (two)
    #   1  hap1 meth bars
    #   2  hap2 meth bars
    #   3  spacer
    #   4  IBD segment (the founder pair A|C is constant across the
    #      whole region; what differs between RBPs is which hiphase
    #      slot — hap1 vs hap2 — carries the A-bearing reads)
    #   5  hap-map blocks (RBP ∩ IBD; the slot↔founder mapping flips
    #      between them — the rebucketing event)
    #   6  spacer
    #   7  pat_A meth bars
    #   8  mat_C meth bars
    gs = GridSpec(
        nrows=9, ncols=1,
        height_ratios=(0.30, 0.95, 0.95, 0.15, 0.30, 0.30, 0.15, 0.95, 0.95),
        hspace=0.25,
        left=0.04, right=0.98, top=0.97, bottom=0.04,
    )
    ax_rbp = fig.add_subplot(gs[0])
    ax_hap1 = fig.add_subplot(gs[1], sharex=ax_rbp)
    ax_hap2 = fig.add_subplot(gs[2], sharex=ax_rbp)
    ax_ibd = fig.add_subplot(gs[4], sharex=ax_rbp)
    ax_hmb = fig.add_subplot(gs[5], sharex=ax_rbp)
    ax_patA = fig.add_subplot(gs[7], sharex=ax_rbp)
    ax_matC = fig.add_subplot(gs[8], sharex=ax_rbp)

    def setup_bar_axes(ax, color, label):
        del color  # row label rendered in black for consistency
        ax.set_xlim(PREFIX_X - 0.3, N_SITES - 0.5)
        ax.set_ylim(0, 1.05)
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ("top", "right", "bottom", "left"):
            ax.spines[spine].set_visible(False)
        ax.text(PREFIX_X, 0.5, label, ha="left", va="center",
                fontsize=LABEL_FS, family=LABEL_FAMILY, color="black")

    def draw_bars(ax, bars, color, label):
        for i, lev in bars:
            if lev > 0:
                ax.bar(i, lev, width=BAR_W, color=color, edgecolor="none",
                       alpha=1.0)
        setup_bar_axes(ax, color, label)

    # Block visual proportions match Fig 4A: rectangles fill ~85% of
    # the row's vertical extent, drawn via the shared helper from
    # `wiki/_panel_grid.py` so styling is identical.
    BLOCK_H = 0.85
    BLOCK_Y = 0.5

    def block_axes(ax, prefix_label):
        ax.text(PREFIX_X, 0.5, prefix_label, ha="left", va="center",
                fontsize=LABEL_FS, family=LABEL_FAMILY, color="black")
        ax.set_xlim(PREFIX_X - 0.3, N_SITES - 0.5)
        ax.set_ylim(0, 1)
        ax.set_axis_off()

    # ---- RBP row: two neutral two-stripe rectangles (Fig 4A vocabulary).
    for lo, hi in (RBP1, RBP2):
        draw_two_stripe_block(
            ax_rbp, lo - 0.5, hi + 0.5, BLOCK_Y, BLOCK_H,
            "neutral_top", "neutral_bot",
            label="hap1  /  hap2", label_fontsize=BLOCK_FS_SHORT,
        )
    block_axes(ax_rbp, "Read-backed phase blocks:")

    # ---- Top tracks: hap1 / hap2 meth, coloured to match the phase
    # block stripes (hap1 = top stripe shade, hap2 = bottom stripe
    # shade).
    draw_bars(ax_hap1, hap1_bars, NEUTRAL_TOP, "hap1 meth")
    draw_bars(ax_hap2, hap2_bars, NEUTRAL_BOT, "hap2 meth")

    # ---- IBD segment: founder pair A|C across the whole region.
    # Pipe convention: left = pat (A, blue), right = mat (C, green) —
    # consistent with Fig 2D.
    draw_two_stripe_block(
        ax_ibd, RBP1[0] - 0.5, RBP2[1] + 0.5, BLOCK_Y, BLOCK_H,
        "A", "C", label="A  |  C", label_fontsize=BLOCK_FS_SHORT,
    )
    block_axes(ax_ibd, "IBD segment:")

    # ---- Hap-map blocks: one per RBP. Both keep top = pat (A, blue)
    # and bottom = mat (C, green); only the slot↔founder mapping
    # changes (HMB1: hap1 = A; HMB2: hap2 = A — i.e., the A-bearing
    # reads are called hap1 in the first phase block and hap2 in the
    # second). Pipe convention: left = pat, right = mat, so the label
    # ordering reflects which hiphase slot is paternal in each block.
    draw_two_stripe_block(
        ax_hmb, RBP1[0] - 0.5, RBP1[1] + 0.5, BLOCK_Y, BLOCK_H,
        "A", "C", label="hap1 = A  |  hap2 = C",
        label_fontsize=BLOCK_FS_LONG,
    )
    draw_two_stripe_block(
        ax_hmb, RBP2[0] - 0.5, RBP2[1] + 0.5, BLOCK_Y, BLOCK_H,
        "A", "C", label="hap2 = A  |  hap1 = C",
        label_fontsize=BLOCK_FS_LONG,
    )
    block_axes(ax_hmb, "Hap-map blocks:")

    # ---- Bottom tracks: founder-labelled (per hap-map block).
    draw_bars(ax_patA, patA_bars, A_COLOR, "pat_A meth")
    draw_bars(ax_matC, matC_bars, C_COLOR, "mat_C meth")

    if show_title:
        fig.suptitle("Rebucketing methylation bars by founder haplotype",
                     fontsize=11, y=0.99)

    out = out_path if out_path is not None else OUT / "rebucket_visual.png"
    if trim_whitespace:
        fig.savefig(out, dpi=180, bbox_inches="tight", pad_inches=0.05)
    else:
        fig.savefig(out, dpi=180)
    plt.close(fig)
    print(f"[scratch] Wrote {out}")


def render_trio_svg(scenario, out_path: Path, kid_y_offset: int = 0) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(build(scenario, kid_y_offset=kid_y_offset,
                               uniform_bars=True))
    print(f'wrote {out_path}')


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    render_trio_svg(SCENARIO_DENOVO, OUT / "trio_denovo.svg", kid_y_offset=40)
    render_trio_svg(SCENARIO_COMPOUND, OUT / "trio_compound_het.svg")
    render_trio_denovo_bed()
    render_trio_compound_het_bed()


if __name__ == '__main__':
    main()
