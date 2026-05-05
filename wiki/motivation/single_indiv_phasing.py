"""Single-individual phasing figures.

Shows how phasing reads exposes haplotype-specific methylation in a
single individual.

Run:
    .venv/bin/python wiki/motivation/single_indiv_phasing.py

Outputs:
    wiki/motivation/single_indiv_before_phasing.png
    wiki/motivation/single_indiv_after_phasing.png

Style requirements:
- Pure monospace ASCII pile, *lots* of 0/1 SNV bits.
- CpG cells occupy the same character cell as a 0 or 1 (single
  monospace column). Rendered as a small filled red/blue rectangle.
- Lots of CpGs (~20).
- Reads staggered: ragged starts/ends; '-' at out-of-window cells.
- Methylation profile above the pile: bigwig-style bars, NO x-axis,
  NO column headers.
- Before phasing: uniform grey, reads in scrambled order, single
  pooled methylation track.
- After phasing: reads grouped pat above mat (with a single
  "Pat hap" / "Mat hap" label per group on the left), two stacked
  methylation tracks.
"""

import random
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle

# Use Arial for all non-monospace text (suptitles, labels, etc.). Texts that
# explicitly request family="monospace" (e.g. the bit string and the polars
# code) are unaffected — they still need a monospace font for column alignment.
mpl.rcParams["font.family"] = "Arial"

random.seed(42)  # deterministic toy

OUT = Path(__file__).resolve().parent
OUT.mkdir(parents=True, exist_ok=True)


# --- Site layout: 13 sites, mixed SNVs (4) + CpGs (9) ---
# 'S' = SNV, 'C' = CpG. 3 flanking sites on each side of the 7-site
# discordant block (originally 5 flanking; trimmed 2 sites from each end).
SITE_KIND = list("CSCCSCCSCCSCC")
N_SITES = len(SITE_KIND)
assert SITE_KIND.count("S") == 4
assert SITE_KIND.count("C") == 9

CPG_INDICES = [i for i, k in enumerate(SITE_KIND) if k == "C"]
SNV_INDICES = [i for i, k in enumerate(SITE_KIND) if k == "S"]


# --- Toy haplotypes ---
# Paternal: SNV bits all '0' (REF). CpG meth pattern picked to give
# visual structure (runs of methylated/unmethylated, plus contrast
# vs the maternal pattern).
# Maternal: SNV bits all '1' (ALT). CpG meth pattern picked to differ
# from pat at multiple positions so the two phased profiles look
# distinct.

PAT = [None] * N_SITES
MAT = [None] * N_SITES
# SNV alleles per haplotype: a mixture of 0s and 1s on each haplotype
# (each SNV is still heterozygous in the individual — pat and mat carry
# opposite alleles — but neither haplotype is purely-REF or purely-ALT).
SNV_PAT_BITS = [1, 0, 1, 1]
SNV_MAT_BITS = [1 - b for b in SNV_PAT_BITS]
for k, i in enumerate(SNV_INDICES):
    PAT[i] = SNV_PAT_BITS[k]
    MAT[i] = SNV_MAT_BITS[k]

# Discordant subregion: pat ~unmethylated, mat ~methylated across CpGs
# whose site index falls in [DISCORDANT_LO, DISCORDANT_HI], with a small
# noise rate so even the discordant region carries some heterogeneity.
DISCORDANT_LO, DISCORDANT_HI = 3, 9
DISCORDANT_NOISE = 0.10  # fraction of reads flipped against the discordant pattern

# Target methylation level per (haplotype, CpG):
#   - inside discordant region: pat near 0, mat near 1.
#   - outside discordant region: drawn from a U-shaped Beta distribution
#     skewed toward 1 (Beta(0.5, 0.2), mean ≈ 0.71) so most CpGs are
#     near-fully-methylated and a few are near-fully-unmethylated. The
#     two haplotypes share the *same* target level at each CpG outside
#     the discordant region (consistent with "all methylated or not,
#     regardless of haplotype, with some heterogeneity").
TARGET_LEVEL = {}  # (source, cpg_i) -> level in [0, 1]
for cpg_i in CPG_INDICES:
    if DISCORDANT_LO <= cpg_i <= DISCORDANT_HI:
        TARGET_LEVEL[("pat", cpg_i)] = DISCORDANT_NOISE
        TARGET_LEVEL[("mat", cpg_i)] = 1.0 - DISCORDANT_NOISE
    else:
        shared = random.betavariate(0.5, 0.5)
        TARGET_LEVEL[("pat", cpg_i)] = shared
        TARGET_LEVEL[("mat", cpg_i)] = shared


# --- Reads: 12 staggered windows, mix of pat/mat sources ---
# Format: (source, start_site, end_site_inclusive)
READS_BY_SOURCE = {
    "pat": [
        (0, 5),
        (0, 6),
        (1, 8),
        (3, 10),
        (5, 12),
        (7, 12),
    ],
    "mat": [
        (0, 4),
        (0, 7),
        (2, 9),
        (4, 11),
        (6, 12),
        (8, 12),
    ],
}
READS_GROUPED = (
    [("pat", s, e) for (s, e) in READS_BY_SOURCE["pat"]]
    + [("mat", s, e) for (s, e) in READS_BY_SOURCE["mat"]]
)
N_READS = len(READS_GROUPED)

# Interleaved order for Fig 1 (no haplotype info).
READS_INTERLEAVED = []
pat_iter = iter(READS_BY_SOURCE["pat"])
mat_iter = iter(READS_BY_SOURCE["mat"])
for _ in range(max(len(READS_BY_SOURCE["pat"]), len(READS_BY_SOURCE["mat"]))):
    p = next(pat_iter, None)
    if p is not None:
        READS_INTERLEAVED.append(("pat", *p))
    m = next(mat_iter, None)
    if m is not None:
        READS_INTERLEAVED.append(("mat", *m))


# --- Colours ---
COLOR_NEUTRAL = "#444444"
COLOR_METH = "#dc2626"    # red
COLOR_UNMETH = "#1d4ed8"  # blue
# Pre-blended (≈55% color over white) so cells can be drawn opaque and
# fully hide the underlying read-tether line.
COLOR_METH_FILL = "#ec8888"
COLOR_UNMETH_FILL = "#82a3e9"
COLOR_OUT = "#cccccc"     # light grey for `-`


# --- Rendering ---
GLYPH_FONTSIZE = 12
PREFIX_FONTSIZE = 12
PREFIX_X = -10.5        # left of column 0; track-label anchors here
CPG_BOX_W = 0.78        # box width in column units (1.0 = one column)
CPG_BOX_H = 0.78        # box height in row units (1.0 = one read row)


# Per-read overrides at specific CpG sites: simulates within-haplotype
# methylation heterogeneity (cell-to-cell variability), so a per-haplotype
# bar can come out at a fraction strictly between 0 and 1.
# Key: (source, read_idx_within_source, site_idx) -> override value.
# Per-read CpG bits sampled from Bernoulli(TARGET_LEVEL).
READ_CPG_VAL = {}  # (source, read_idx, cpg_i) -> 0 or 1
for src in ("pat", "mat"):
    for k, (start, end) in enumerate(READS_BY_SOURCE[src]):
        for cpg_i in CPG_INDICES:
            if start <= cpg_i <= end:
                p = TARGET_LEVEL[(src, cpg_i)]
                READ_CPG_VAL[(src, k, cpg_i)] = 1 if random.random() < p else 0


def _read_data(source: str, read_idx: int, start: int, end: int):
    hap = PAT if source == "pat" else MAT
    out = []
    for i in range(N_SITES):
        in_window = start <= i <= end
        if not in_window:
            out.append((None, False))
        elif SITE_KIND[i] == "S":
            out.append((hap[i], True))
        else:
            out.append((READ_CPG_VAL[(source, read_idx, i)], True))
    return out


def _draw_read_row(ax, row_y, source, read_idx, start, end, hap_color,
                   prefix, prefix_color, glyph_fontsize=GLYPH_FONTSIZE,
                   snv_cell_background=False):
    """Draw one read row.

    snv_cell_background=False (option 1): a single continuous connector
    spans the read; SNV digits get a white text-bbox so they mask the
    connector locally.
    snv_cell_background=True (option 2): the connector is drawn only in
    the gaps between sites, and SNV digits get a CpG-sized neutral cell
    rectangle behind them so they look like cells too."""
    ax.text(
        PREFIX_X, row_y, prefix,
        ha="left", va="center", fontsize=PREFIX_FONTSIZE,
        family="Arial", color=prefix_color,
    )
    if snv_cell_background:
        half_w = CPG_BOX_W / 2
        for i in range(start, end):
            ax.plot([i + half_w, (i + 1) - half_w], [row_y, row_y],
                    color="#bbbbbb", linewidth=0.8, solid_capstyle="butt",
                    zorder=0)
    else:
        ax.plot([start, end], [row_y, row_y],
                color="#bbbbbb", linewidth=0.8, solid_capstyle="butt",
                zorder=0)
    for i, (val, in_window) in enumerate(_read_data(source, read_idx, start, end)):
        if not in_window:
            continue
        if SITE_KIND[i] == "S":
            if snv_cell_background:
                ax.add_patch(
                    Rectangle(
                        (i - CPG_BOX_W / 2, row_y - CPG_BOX_H / 2),
                        CPG_BOX_W, CPG_BOX_H,
                        facecolor="#e6e6e6", edgecolor="none",
                    )
                )
                ax.text(
                    i, row_y, str(val),
                    ha="center", va="center",
                    fontsize=glyph_fontsize, family="Arial",
                    color=hap_color,
                )
            else:
                ax.text(
                    i, row_y, str(val),
                    ha="center", va="center",
                    fontsize=glyph_fontsize, family="Arial",
                    color=hap_color,
                    bbox=dict(facecolor="white", edgecolor="none", pad=2.0),
                )
        else:  # CpG: same-cell-size colored rectangle (opaque, masks tether)
            box_color = COLOR_METH_FILL if val else COLOR_UNMETH_FILL
            ax.add_patch(
                Rectangle(
                    (i - CPG_BOX_W / 2, row_y - CPG_BOX_H / 2),
                    CPG_BOX_W, CPG_BOX_H,
                    facecolor=box_color, edgecolor="none",
                )
            )


DISCORDANT_TRACK_COLOR = "#000000"   # black


def _draw_discordant_track(ax, label="Discordant region", xlim_left=None):
    """An IGV-style feature track: one horizontal bar spanning the
    [DISCORDANT_LO, DISCORDANT_HI] interval, with a side label."""
    ax.add_patch(Rectangle(
        (DISCORDANT_LO - 0.5, 0.25),
        (DISCORDANT_HI - DISCORDANT_LO) + 1.0, 0.5,
        facecolor=DISCORDANT_TRACK_COLOR, edgecolor="none",
    ))
    if label:
        ax.text(
            PREFIX_X, 0.5, label,
            ha="left", va="center", fontsize=PREFIX_FONTSIZE,
            family="Arial", color=DISCORDANT_TRACK_COLOR,
        )
    ax.set_xlim(PREFIX_X - 0.5 if xlim_left is None else xlim_left,
                N_SITES - 0.5)
    ax.set_ylim(0, 1)
    ax.set_axis_off()


LABELED_SITES = (3,)  # one informative discordant CpG (pat all-unmethylated, mat all-methylated)
COUNT_LABEL_FONTSIZE = 28
TALL_BAR_THRESHOLD = 0.5


def _draw_count_labels(ax, levels_dict, counts_dict, color, sites=LABELED_SITES,
                       fontsize=COUNT_LABEL_FONTSIZE, force_outside=False):
    """Annotate selected CpG bars with stacked m-over-n read counts (rendered
    as a mathtext fraction with a horizontal rule; mathtext uses matplotlib's
    own math font, not Arial). Tall bars get a white label inside the bar
    near its top; short bars get a black label just above the bar. With
    force_outside=True, every label is placed above the bar in black
    regardless of height. Bars below MIN_COVERAGE are skipped."""
    for cpg_i in sites:
        if cpg_i not in counts_dict:
            continue
        m, n = counts_dict[cpg_i]
        level = levels_dict.get(cpg_i, 0.0)
        text = rf'$\frac{{{m}}}{{{n}}}$'
        # Mathtext \frac bbox is slightly left-skewed (the fraction-bar tail
        # extends past the numerator/denominator on the left), so nudge the
        # anchor right by a small amount to visually center it on the bar.
        x_nudge = 0.07
        if not force_outside and level >= TALL_BAR_THRESHOLD:
            ax.text(cpg_i + x_nudge, level - 0.05, text,
                    ha="center", va="top",
                    fontsize=fontsize, color="white",
                    fontweight="bold")
        else:
            ax.text(cpg_i + x_nudge, level + 0.05, text,
                    ha="center", va="bottom",
                    fontsize=fontsize, color="black")


def _draw_meth_bars(ax, levels_dict, color, track_label=None,
                    ytick_fontsize=8, xlim_left=None,
                    ylabel=None, ylabel_fontsize=10):
    """Bigwig-style methylation track. No x-axis, no column annotations."""
    BAR_ALPHA = 0.55
    for cpg_i, level in levels_dict.items():
        if level == 0:
            # Bar of height 0 is invisible; draw a short horizontal tick
            # at y=0 over this CpG's column to mark "covered, zero methylation".
            ax.plot(
                [cpg_i - CPG_BOX_W / 2, cpg_i + CPG_BOX_W / 2],
                [0, 0],
                color=color, linewidth=1.4, solid_capstyle="butt",
                alpha=BAR_ALPHA,
            )
        else:
            ax.bar(
                cpg_i, level, width=CPG_BOX_W, color=color,
                edgecolor="none", align="center", alpha=BAR_ALPHA,
            )
    if track_label:
        ax.text(
            PREFIX_X, 0.5, track_label,
            ha="left", va="center", fontsize=PREFIX_FONTSIZE,
            family="Arial", color=color,
        )
    ax.set_xlim(PREFIX_X - 0.5 if xlim_left is None else xlim_left,
                N_SITES - 0.5)
    ax.set_ylim(0, 1.05)
    ax.set_xticks([])
    ax.set_yticks([0, 1])
    ax.set_yticklabels(["0", "1"], fontsize=ytick_fontsize)
    for spine in ("top", "right", "bottom"):
        ax.spines[spine].set_visible(False)
    # Move the y-axis (with its 0/1 ticks) right next to the bars and
    # bound the spine so it stops at the "1" tick.
    ax.spines["left"].set_position(("data", -0.8))
    ax.spines["left"].set_bounds(0, 1)
    ax.tick_params(axis="y", length=2, labelsize=ytick_fontsize)
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=ylabel_fontsize, color=color,
                      family="Arial", rotation=90, labelpad=4)


def _per_source_indexed_reads():
    """Yield (source, read_idx_within_source, start, end) tuples."""
    for src in ("pat", "mat"):
        for k, (start, end) in enumerate(READS_BY_SOURCE[src]):
            yield (src, k, start, end)


MIN_COVERAGE = 2  # don't report a methylation level below this read count


def _pooled_methylation():
    levels = {}
    counts = {}
    for cpg_i in CPG_INDICES:
        meth = 0
        cov = 0
        for src, k, start, end in _per_source_indexed_reads():
            if start <= cpg_i <= end:
                cov += 1
                meth += READ_CPG_VAL[(src, k, cpg_i)]
        if cov >= MIN_COVERAGE:
            levels[cpg_i] = meth / cov
            counts[cpg_i] = (meth, cov)
    return levels, counts


def _phased_methylation():
    levels = {"pat": {}, "mat": {}}
    counts = {"pat": {}, "mat": {}}
    for cpg_i in CPG_INDICES:
        for src in ("pat", "mat"):
            meth = 0
            cov = 0
            for k, (start, end) in enumerate(READS_BY_SOURCE[src]):
                if start <= cpg_i <= end:
                    cov += 1
                    meth += READ_CPG_VAL[(src, k, cpg_i)]
            if cov >= MIN_COVERAGE:
                levels[src][cpg_i] = meth / cov
                counts[src][cpg_i] = (meth, cov)
    return levels, counts


# --- Figures ---

def render_before_phasing(out_path: Path | None = None,
                          glyph_fontsize: int = GLYPH_FONTSIZE,
                          ytick_fontsize: int = 8,
                          show_left_labels: bool = True,
                          show_title: bool = True,
                          ylabel: str | None = None,
                          ylabel_fontsize: int = 10):
    fig = plt.figure(figsize=(5.5, 5.7))
    gs = GridSpec(
        nrows=3, ncols=1,
        height_ratios=(1.0, 0.25, 4.5),
        hspace=0.18,
        left=0.04, right=0.98, top=0.92, bottom=0.04,
    )
    ax_meth = fig.add_subplot(gs[0])
    ax_discordant = fig.add_subplot(gs[1], sharex=ax_meth)
    ax_reads = fig.add_subplot(gs[2], sharex=ax_meth)

    pooled, pooled_counts = _pooled_methylation()
    xlim_left = None if show_left_labels else -1.5
    _draw_meth_bars(
        ax_meth, pooled, COLOR_NEUTRAL,
        track_label="Meth level" if show_left_labels else None,
        ytick_fontsize=ytick_fontsize,
        xlim_left=xlim_left,
        ylabel=ylabel,
        ylabel_fontsize=ylabel_fontsize,
    )
    _draw_count_labels(ax_meth, pooled, pooled_counts, COLOR_NEUTRAL,
                       force_outside=True)
    _draw_discordant_track(
        ax_discordant,
        label="Discordant region" if show_left_labels else None,
        xlim_left=xlim_left,
    )

    src_seen = {"pat": 0, "mat": 0}
    for row, (src, start, end) in enumerate(READS_INTERLEAVED):
        k = src_seen[src]
        src_seen[src] += 1
        _draw_read_row(
            ax_reads, row_y=row, source=src, read_idx=k,
            start=start, end=end,
            hap_color=COLOR_NEUTRAL, prefix="",
            prefix_color=COLOR_NEUTRAL,
            glyph_fontsize=glyph_fontsize,
        )
    ax_reads.set_xlim(PREFIX_X - 0.5 if xlim_left is None else xlim_left,
                      N_SITES - 0.5)
    ax_reads.set_ylim(N_READS - 0.5, -1.0)
    ax_reads.set_axis_off()

    if show_title:
        fig.suptitle(
            'Before phasing: unphased reads, one pooled methylation profile',
            fontsize=12, y=0.97,
        )
    out = out_path if out_path is not None else OUT / "single_indiv_before_phasing.png"
    fig.savefig(out, dpi=180)
    plt.close(fig)
    print(f"[scratch] Wrote {out}")


def render_after_phasing(out_path: Path | None = None,
                         glyph_fontsize: int = GLYPH_FONTSIZE,
                         ytick_fontsize: int = 8,
                         show_left_labels: bool = True,
                         show_title: bool = True,
                         pat_ylabel: str | None = None,
                         mat_ylabel: str | None = None,
                         ylabel_fontsize: int = 10):
    fig = plt.figure(figsize=(5.5, 6.7))
    # gs[1] is an empty spacer row that widens the visible gap between the
    # pat and mat methylation tracks without affecting the other gaps.
    gs = GridSpec(
        nrows=5, ncols=1,
        height_ratios=(0.85, 0.3, 0.85, 0.25, 4.5),
        hspace=0.18,
        left=0.04, right=0.98, top=0.93, bottom=0.04,
    )
    ax_pat = fig.add_subplot(gs[0])
    ax_mat = fig.add_subplot(gs[2], sharex=ax_pat)
    ax_discordant = fig.add_subplot(gs[3], sharex=ax_pat)
    ax_reads = fig.add_subplot(gs[4], sharex=ax_pat)

    xlim_left = None if show_left_labels else -1.5
    phased, phased_counts = _phased_methylation()
    _draw_meth_bars(
        ax_pat, phased["pat"], COLOR_NEUTRAL,
        track_label="Pat meth" if show_left_labels else None,
        ytick_fontsize=ytick_fontsize, xlim_left=xlim_left,
        ylabel=pat_ylabel, ylabel_fontsize=ylabel_fontsize,
    )
    _draw_count_labels(ax_pat, phased["pat"], phased_counts["pat"], COLOR_NEUTRAL)
    _draw_meth_bars(
        ax_mat, phased["mat"], COLOR_NEUTRAL,
        track_label="Mat meth" if show_left_labels else None,
        ytick_fontsize=ytick_fontsize, xlim_left=xlim_left,
        ylabel=mat_ylabel, ylabel_fontsize=ylabel_fontsize,
    )
    _draw_count_labels(ax_mat, phased["mat"], phased_counts["mat"], COLOR_NEUTRAL)
    _draw_discordant_track(
        ax_discordant,
        label="Discordant region" if show_left_labels else None,
        xlim_left=xlim_left,
    )

    n_pat = len(READS_BY_SOURCE["pat"])
    GROUP_GAP = 1.2  # vertical space between the pat and mat groups
    last_y = 0.0
    src_seen = {"pat": 0, "mat": 0}
    for row, (src, start, end) in enumerate(READS_GROUPED):
        y = row + (GROUP_GAP if row >= n_pat else 0.0)
        k = src_seen[src]
        src_seen[src] += 1
        _draw_read_row(
            ax_reads, row_y=y, source=src, read_idx=k,
            start=start, end=end,
            hap_color=COLOR_NEUTRAL, prefix="", prefix_color=COLOR_NEUTRAL,
            glyph_fontsize=glyph_fontsize,
        )
        last_y = y
    n_mat = len(READS_BY_SOURCE["mat"])
    pat_center = (n_pat - 1) / 2.0
    mat_center = n_pat + GROUP_GAP + (n_mat - 1) / 2.0
    if show_left_labels:
        for y, label in ((pat_center, "Pat hap"), (mat_center, "Mat hap")):
            ax_reads.text(
                PREFIX_X, y, label,
                ha="left", va="center", fontsize=PREFIX_FONTSIZE,
                family="Arial", color=COLOR_NEUTRAL,
            )

    ax_reads.set_xlim(PREFIX_X - 0.5 if xlim_left is None else xlim_left,
                      N_SITES - 0.5)
    ax_reads.set_ylim(last_y + 0.6, -1.0)
    ax_reads.set_axis_off()

    if show_title:
        fig.suptitle(
            'After phasing: phased reads, one methylation profile per haplotype',
            fontsize=12, y=0.98,
        )
    out = out_path if out_path is not None else OUT / "single_indiv_after_phasing.png"
    fig.savefig(out, dpi=180)
    plt.close(fig)
    print(f"[scratch] Wrote {out}")



def main() -> None:
    render_before_phasing()
    render_after_phasing()


if __name__ == "__main__":
    main()
