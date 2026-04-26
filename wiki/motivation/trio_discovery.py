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
    f'font-family="Helvetica, Arial, sans-serif" font-size="14">'
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

def haplotype(x, y, w, color, cpgs, variants=(), label=None, label_at_start=False):
    """cpgs: list of (rel_x, methylated_bool). variants: list of rel_x drawn as red stars.
    label: optional single-letter hap id (drawn just outside the bar, adjacent
    to the person's shape). label_at_start puts it at bar_x0 (left end)."""
    rels = [r for r, _ in cpgs] + list(variants)
    pad = min(rels)
    bar_x0 = x
    bar_x1 = x + (max(rels) + pad) * w
    out = [f'<rect x="{bar_x0}" y="{y - 4}" width="{bar_x1 - bar_x0}" height="8" fill="{color}"/>']
    if label is not None:
        if label_at_start:
            out.append(
                f'<text x="{bar_x0 - 10}" y="{y + 5}" text-anchor="end" '
                f'font-family="Arial, sans-serif" font-size="16">{label}</text>'
            )
        else:
            out.append(
                f'<text x="{bar_x1 + 10}" y="{y + 5}" text-anchor="start" '
                f'font-family="Arial, sans-serif" font-size="16">{label}</text>'
            )
    stem_h = 24
    r = 7
    for rx, meth in cpgs:
        cx = x + rx * w
        out.append(line(cx, y - 4, cx, y - 4 - stem_h, sw=1.2))
        out.append(circle(cx, y - 4 - stem_h - r + 2, r,
                          fill='black' if meth else 'white', sw=1.2))
    for rx in variants:
        cx = x + rx * w
        out.append(star(cx, y, r=12))
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
    mom_top=HapData(METH5, '#8ecae6', label='C'),
    mom_bot=HapData(UNMETH5, '#1f6ea8', label='D'),
    dad_top=HapData(UNMETH5, '#f4b860', label='A'),
    dad_bot=HapData(METH5, '#e8734a', label='B'),
    # Kid pat = dad A; kid mat = mom C. Two de novo events both sit at
    # the right CpG cluster (positions KID_RELS[3:]). Paternal A inherits
    # dad's unmethylated state at the left cluster but is methylated at
    # the right cluster (de novo gain). Maternal C inherits mom's
    # methylated state at the left cluster but is unmethylated at the
    # right cluster (de novo loss).
    kid_top=HapData([(r, m) for r, m in zip(KID_RELS, [False, False, False, True, True])],
                    '#f4b860', label='A'),
    kid_bot=HapData([(r, m) for r, m in zip(KID_RELS, [True, True, True, False, False])],
                    '#8ecae6', label='C'),
)

# Compound genetic-epigenetic heterozygote:
#   paternal allele in kid carries a genetic variant (star)
#   maternal allele in kid carries aberrant methylation (silenced promoter-like region)
# Each parent is a silent carrier of one of the two hits.
SCENARIO_COMPOUND = Scenario(
    # Mom: C = hyper-methylated carrier; D = normal unmethylated
    mom_top=HapData([(r, True) for r in RELS[:3]] + [(r, False) for r in RELS[3:]],
                    '#8ecae6', label='C'),
    mom_bot=HapData(UNMETH5, '#1f6ea8', label='D'),
    # Dad: A = genetic-variant carrier; B = normal
    dad_top=HapData(UNMETH5, '#f4b860', variants=(0.40,), label='A'),
    dad_bot=HapData(UNMETH5, '#e8734a', label='B'),
    # Kid: paternal = A (variant); maternal = C (aberrant meth) — compound het in trans
    kid_top=HapData([(r, False) for r in KID_RELS], '#f4b860', variants=(0.42,), label='A'),
    kid_bot=HapData([(r, True) for r in KID_RELS[:3]] + [(r, False) for r in KID_RELS[3:]],
                   '#8ecae6', label='C'),
)


# ---------- layout ----------

def build(scenario: Scenario) -> str:
    parts = [SVG_HEAD]

    mom = (560, 190)
    dad = (800, 190)
    parts.append(person(*mom, 'F'))
    parts.append(person(*dad, 'M'))

    # Mom/kid labels sit on the bar's tail (person side); dad labels sit on the head.
    parts.append(haplotype(240, 184, 340, scenario.mom_top.color,
                           scenario.mom_top.cpgs, scenario.mom_top.variants,
                           label=scenario.mom_top.label))
    parts.append(haplotype(240, 234, 340, scenario.mom_bot.color,
                           scenario.mom_bot.cpgs, scenario.mom_bot.variants,
                           label=scenario.mom_bot.label))
    parts.append(haplotype(890, 184, 280, scenario.dad_top.color,
                           scenario.dad_top.cpgs, scenario.dad_top.variants,
                           label=scenario.dad_top.label, label_at_start=True))
    parts.append(haplotype(890, 234, 280, scenario.dad_bot.color,
                           scenario.dad_bot.cpgs, scenario.dad_bot.variants,
                           label=scenario.dad_bot.label, label_at_start=True))

    parts.append(line(mom[0] + 50, mom[1], dad[0] - 50, dad[1]))
    midx = (mom[0] + dad[0]) / 2
    son = (int(midx), 400)
    parts.append(line(midx, mom[1], son[0], son[1] - 50))
    parts.append(person(*son, 'M'))

    parts.append(haplotype(380, 394, 280, scenario.kid_top.color,
                           scenario.kid_top.cpgs, scenario.kid_top.variants,
                           label=scenario.kid_top.label))
    parts.append(haplotype(380, 444, 280, scenario.kid_bot.color,
                           scenario.kid_bot.cpgs, scenario.kid_bot.variants,
                           label=scenario.kid_bot.label))

    parts.append(SVG_TAIL)
    return '\n'.join(parts)


# ---------- matplotlib companion: BED snippet + polars query ----------

def render_trio_denovo_bed():
    """BED snippet + polars query for discovering de novo epimutations."""
    FIG_W, FIG_H = 12.0, 9.0
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

    cols = ("chrom", "start", "kid_pat", "pat_hap", "dad_A", "dad_B")
    # CpG genomic coordinates mirror KID_RELS in the trio_denovo SVG:
    #   left cluster of 3 CpGs at 1100 / 1200 / 1300 (concordant
    #     unmethylated on dad_A → kid_pat),
    #   right cluster of 2 CpGs at 1550 / 1650 (de novo gain — kid_pat
    #     methylated where dad_A is not).
    # dad_B is the non-transmitted, fully methylated homolog.
    # Methylation levels jittered around the SVG values for realism.
    rows = [
        ("chr1", "1100", "0.05", "A", "0.04", "0.93"),
        ("chr1", "1200", "0.04", "A", "0.06", "0.95"),
        ("chr1", "1300", "0.06", "A", "0.05", "0.92"),
        ("chr1", "1550", "0.94", "A", "0.05", "0.96"),
        ("chr1", "1650", "0.96", "A", "0.04", "0.94"),
    ]
    _draw_simple_table(
        ax_meth, cols, rows,
        title="Haplotype-specific methylation levels",
        highlight_rows={3, 4},
        bottom_rule=False,
    )

    ax_code.set_xlim(0, 1)
    ax_code.set_ylim(0, 1)
    ax_code.set_axis_off()
    CODE_TITLE_Y = 0.97
    CODE_TITLE_TO_CODE_IN = (0.97 - 0.72) * TABLE_H_IN
    code_y = CODE_TITLE_Y - CODE_TITLE_TO_CODE_IN / code_h_in
    ax_code.text(
        0.0, CODE_TITLE_Y,
        "Discover de novo gain of methylation on paternal haplotype",
        ha="left", va="top", fontsize=12, fontweight="bold",
        transform=ax_code.transAxes,
    )

    CODE_FONTSIZE = 10.5
    code = (
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
    ax_code.text(
        0.0, code_y, code,
        ha="left", va="top", fontsize=CODE_FONTSIZE,
        family="Menlo", color=COLOR_NEUTRAL,
        transform=ax_code.transAxes,
    )

    # Place the italic note at the same absolute gap below the code block as
    # VISIBLE_GAP_IN (the table↔code gap).
    code_n_lines = code.count('\n')
    line_h_in = CODE_FONTSIZE / 72.0 * 1.2  # ~matplotlib default linespacing
    code_h_used_in = code_n_lines * line_h_in
    note_y = code_y - (code_h_used_in + VISIBLE_GAP_IN) / code_h_in
    ax_code.text(
        0.0, note_y,
        "Note: in the highlighted rows, |kid_pat − dad_A| ≈ 0.9 (the true epimutation signal) "
        "but |kid_pat − dad_B| ≈ 0.5 — picking the wrong dad\n"
        "haplotype roughly halves the delta and would barely clear the threshold, "
        "motivating the need for inheritance-based phasing.",
        ha="left", va="top", fontsize=10,
        family="Arial", color=COLOR_NEUTRAL, style="italic",
        transform=ax_code.transAxes,
    )

    out = OUT / "trio_denovo_bed.png"
    fig.savefig(out, dpi=180)
    plt.close(fig)
    print(f"[scratch] Wrote {out}")


def _draw_simple_table(ax, cols, rows, title, highlight_rows=frozenset(),
                        bottom_rule=True):
    """Render a small fixed-width table on a transAxes axis."""
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_axis_off()
    ax.text(
        0.0, 0.97, title,
        ha="left", va="top", fontsize=12, fontweight="bold",
        transform=ax.transAxes,
    )
    n_col = len(cols)
    col_xs = [0.02 + 0.085 * i for i in range(n_col)]
    table_right = col_xs[-1] + 0.07

    header_y = 0.72
    for x, name in zip(col_xs, cols):
        ax.text(
            x, header_y, name,
            ha="left", va="center", fontsize=10.5,
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
                facecolor="#e5e5e5", edgecolor="none",
                transform=ax.transAxes,
            ))
        for x, cell in zip(col_xs, row):
            ax.text(
                x, y, cell,
                ha="left", va="center", fontsize=10.5,
                family="monospace", color=COLOR_NEUTRAL,
                transform=ax.transAxes,
            )
    if bottom_rule:
        rule_y = header_y - 0.12 - (len(rows) - 1) * row_dy - row_dy / 2 + 0.01
        ax.plot(
            [0.005, table_right], [rule_y, rule_y],
            color=COLOR_NEUTRAL, linewidth=0.7, transform=ax.transAxes,
        )


def render_trio_compound_het_bed():
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
    meth_cols = ("chrom", "start",
                 "kid_pat", "pat_hap", "kid_mat", "mat_hap",
                 "dad_A", "dad_B", "mom_C", "mom_D")
    meth_rows = [
        ("chr1", "1100", "0.05", "A", "0.94", "C", "0.04", "0.05", "0.93", "0.06"),
        ("chr1", "1200", "0.04", "A", "0.96", "C", "0.06", "0.04", "0.95", "0.05"),
        ("chr1", "1300", "0.06", "A", "0.91", "C", "0.05", "0.06", "0.92", "0.04"),
        ("chr1", "1550", "0.05", "A", "0.04", "C", "0.04", "0.05", "0.05", "0.06"),
        ("chr1", "1650", "0.05", "A", "0.04", "C", "0.05", "0.04", "0.06", "0.03"),
    ]
    _draw_simple_table(
        ax_meth, meth_cols, meth_rows,
        title="Haplotype-specific methylation levels",
        highlight_rows={0, 1, 2},
        bottom_rule=False,
    )

    # Genotype at a SNV inside the same locus: exactly one parental hap
    # (dad_A) carries the ALT allele, and the kid inherits it on the
    # paternal side — in trans with the meth outlier.
    geno_cols = ("chrom", "pos",
                 "dad_A", "dad_B", "mom_C", "mom_D",
                 "kid_pat", "kid_mat")
    geno_rows = [
        ("chr1", "1420", "1", "0", "0", "0", "1", "0"),
    ]
    _draw_simple_table(
        ax_geno, geno_cols, geno_rows,
        title="Phased genotypes",
        highlight_rows={0},
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
        "Discover compound genetic-epigenetic heterozygous locus",
        ha="left", va="top", fontsize=12, fontweight="bold",
        transform=ax_code.transAxes,
    )
    code = (
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
    ax_code.text(
        0.0, code_y, code,
        ha="left", va="top", fontsize=10.5,
        family="Menlo", color=COLOR_NEUTRAL,
        transform=ax_code.transAxes,
    )

    out = OUT / "trio_compound_het_bed.png"
    fig.savefig(out, dpi=180)
    plt.close(fig)
    print(f"[scratch] Wrote {out}")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    for name, scen in [('trio_denovo', SCENARIO_DENOVO),
                       ('trio_compound_het', SCENARIO_COMPOUND)]:
        path = OUT / f'{name}.svg'
        path.write_text(build(scen))
        print(f'wrote {path}')
    render_trio_denovo_bed()
    render_trio_compound_het_bed()


if __name__ == '__main__':
    main()
