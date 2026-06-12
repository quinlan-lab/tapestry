"""Rank chromosome-spanning hap-map blocks by single-crossover cleanliness.

Each pedMEC phase block spans (most of) a chromosome, so a block's concordance
already tells you *roughly where* a crossover sits along the chromosome (a
crossover at fractional position c down the chromosome leaves the kid matching
the assigned haplotype on a fraction ~c of the het sites). What concordance
CANNOT tell you is whether that crossover is CLEAN -- one sharp transition --
or messy -- many phase switches that happen to average to the same value.

This tool adds the missing spatial test using the bit-vector mismatch sites.
For a block assigned the majority haplotype H, every mismatch site is a het SNV
where the kid matches the OTHER haplotype. Under a single clean crossover those
mismatches are absent on the majority side and densely contiguous on the
minority side, with the transition edge = the breakpoint. Under phase-switch
noise they are scattered across the whole block.

Cleanliness score
-----------------
Let c = block concordance (majority fraction), f = 1 - c (minority fraction),
and map each mismatch position to u = (pos - start) / (end - start) in [0, 1].
ASSUMING roughly uniform het-SNV density along the block, a single crossover
puts the minority segment either on the right (u in (c, 1]) or left (u in
[0, f)). We take whichever side holds more mismatches:

    frac = max( mean(u > c), mean(u < f) )

Random scatter gives frac ~= f; a perfectly clean crossover gives frac ~= 1.
We rescale to a 0..1 score:

    cleanliness = clip( (frac - f) / (1 - f), 0, 1 )

so 0 = mismatches scattered as if random, 1 = mismatches entirely on the
predicted minority side (a single sharp crossover).

The uniform-density assumption is approximate (SNV density dips at centromeres
etc.), so treat this as a RANKING to choose which candidates to confirm in IGV,
not a final breakpoint caller. The reported breakpoint interval brackets the
predicted split position and the innermost minority mismatch; zoom there in IGV
to read the flanking het coordinates.

True crossover vs kid phase-switch
----------------------------------
A real meiotic crossover flips only ONE parent's transmitted haplotype, so a
clean transition should appear in the paternal OR maternal track but not both
at the same kid coordinate. When both tracks have a clean breakpoint within
--same-site-kb of each other on the same chromosome, that is the signature of a
kid phase-switch error, and the candidate is flagged (prefer single-parent
flips for a dev crossover).

IGV confirmation workflow
-------------------------
--igv-bed writes a navigation track of the candidate breakpoint intervals
(cleanliness-ranked); load it in IGV alongside the
bit-vector-sites-mismatches.{paternal,maternal}.vcf.gz tracks and use
Ctrl-F / Ctrl-B to step through candidates. --igv-batch writes a batch script
(`igv.sh -b <path>` or Tools > Run Batch Script) that does a `goto` zoomed on
each estimated breakpoint (+/- --igv-window bp), cleanest-first, optionally
snapshotting a PNG per candidate with --snapshot-dir. At each candidate confirm:
the mismatch density should show a SHARP transition from sparse ticks (kid
agrees) to a dense contiguous run to the block edge (past the crossover), and
that transition should appear in only ONE parent's track (a true meiotic
crossover); a transition at the same kid coordinate in BOTH tracks is a kid
phase-switch artifact.

Usage:
    PYTHONPATH=src:src/util .venv/bin/python src/util/find_clean_crossover.py \\
        --paternal-blocks    <out>/NA12878.hap-map-blocks.paternal.sorted.bed.gz \\
        --maternal-blocks    <out>/NA12878.hap-map-blocks.maternal.sorted.bed.gz \\
        --paternal-mismatch  <out>/NA12878.bit-vector-sites-mismatches.paternal.bed \\
        --maternal-mismatch  <out>/NA12878.bit-vector-sites-mismatches.maternal.bed \\
        --out candidates.tsv --igv-bed candidates.igv.bed --igv-batch candidates.igv.bat
"""

import argparse
import logging

import polars as pl


def _load_blocks(bed):
    """hap-map-blocks BED -> chrom,start,end,haplotype,concordance,num_het_SNVs."""
    return (
        pl.read_csv(bed, separator="\t", has_header=False,
                    new_columns=["chrom", "start", "end", "name"],
                    comment_prefix="#")
        .with_columns(
            haplotype=pl.col("name").str.split(",").list.get(0),
            concordance=pl.col("name").str.split(",").list.get(1).cast(pl.Float64),
            num_het_SNVs=pl.col("name").str.split(",").list.get(2).cast(pl.Int64),
        )
        .select(["chrom", "start", "end", "haplotype", "concordance", "num_het_SNVs"])
    )


def _load_mismatch(bed):
    """bit-vector mismatch BED -> chrom,start (only the first two columns needed)."""
    try:
        df = pl.read_csv(bed, separator="\t", has_header=False,
                         new_columns=["chrom", "start", "end", "REF", "ALT"],
                         comment_prefix="#")
    except pl.exceptions.NoDataError:
        return pl.DataFrame({"chrom": [], "start": []},
                            schema={"chrom": pl.String, "start": pl.Int64})
    return df.select(["chrom", "start"])


def _score_block(block, mm_pos):
    """Cleanliness score + estimated breakpoint for one block.

    block: dict with chrom,start,end,concordance,num_het_SNVs.
    mm_pos: sorted list of mismatch positions inside the block.
    Returns dict or None if too few mismatches to judge.
    """
    start, end = block["start"], block["end"]
    L = end - start
    c = block["concordance"]
    f = 1.0 - c
    m = len(mm_pos)
    if L <= 0 or m < 2 or f <= 0:
        return None

    u = [(p - start) / L for p in mm_pos]
    frac_right = sum(1 for x in u if x > c) / m       # minority on the right
    frac_left = sum(1 for x in u if x < f) / m        # minority on the left

    if frac_right >= frac_left:
        side, frac = "right", frac_right
        break_pos = start + round(c * L)
        inner = min((p for p, x in zip(mm_pos, u) if x > c), default=break_pos)
    else:
        side, frac = "left", frac_left
        break_pos = start + round(f * L)
        inner = max((p for p, x in zip(mm_pos, u) if x < f), default=break_pos)

    cleanliness = max(0.0, min(1.0, (frac - f) / (1.0 - f)))
    lo, hi = sorted((break_pos, inner))
    return {
        "cleanliness": round(cleanliness, 3),
        "side": side,
        "break_pos": break_pos,
        "break_interval": f"{block['chrom']}:{lo}-{hi}",
        "num_mismatch": m,
    }


def write_igv_bed(df, path):
    """Write candidate breakpoint intervals as an IGV navigation track.

    Load in IGV with the mismatch VCF tracks; select the track and step with
    Ctrl-F / Ctrl-B (cleanest-first by file order is not guaranteed -- IGV walks
    genomic order -- but the feature name carries the cleanliness so you can
    prioritize visually). The feature name encodes
    "<parent>,clean=<cleanliness>,conc=<concordance>,<crossover_kind>".
    """
    (
        df.with_columns(
            chrom=pl.col("break_interval").str.split(":").list.get(0),
            _range=pl.col("break_interval").str.split(":").list.get(1),
        )
        .with_columns(
            start=pl.col("_range").str.split("-").list.get(0).cast(pl.Int64),
            end=pl.col("_range").str.split("-").list.get(1).cast(pl.Int64),
            name=pl.col("parent") + pl.lit(",clean=") + pl.col("cleanliness").cast(pl.String)
            + pl.lit(",conc=") + pl.col("concordance").cast(pl.String)
            + pl.lit(",") + pl.col("crossover_kind"),
        )
        .select(["chrom", "start", "end", "name"])
        .sort(["chrom", "start", "end"])
        .write_csv(path, separator="\t", include_header=False)
    )


def write_igv_batch(df, path, window, snapshot_dir=None):
    """Write an IGV batch script that visits each candidate breakpoint zoomed.

    Cleanest-first. Each candidate produces a `goto` to break_pos +/- window so
    you land on the transition (not the whole chromosome). With snapshot_dir a
    PNG is written per candidate. Run via `igv.sh -b <path>` or Tools > Run
    Batch Script.
    """
    lines = []
    if snapshot_dir:
        lines.append(f"snapshotDirectory {snapshot_dir}")
    for r in df.iter_rows(named=True):
        lo = max(0, r["break_pos"] - window)
        hi = r["break_pos"] + window
        lines.append(f"goto {r['chrom']}:{lo}-{hi}")
        if snapshot_dir:
            lines.append(f"snapshot {r['parent']}.{r['chrom']}_{r['break_pos']}.png")
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")


def _analyze_track(blocks_bed, mismatch_bed, parent, min_num_het, max_concordance):
    blocks = _load_blocks(blocks_bed)
    mism = _load_mismatch(mismatch_bed)
    rows = []
    for block in blocks.iter_rows(named=True):
        if block["num_het_SNVs"] < min_num_het or block["concordance"] >= max_concordance:
            continue
        mm_pos = (
            mism.filter(
                (pl.col("chrom") == block["chrom"])
                & (pl.col("start") >= block["start"])
                & (pl.col("start") <= block["end"])
            )
            .get_column("start").sort().to_list()
        )
        scored = _score_block(block, mm_pos)
        if scored is None:
            continue
        rows.append({
            "parent": parent,
            "chrom": block["chrom"],
            "concordance": block["concordance"],
            "num_het_SNVs": block["num_het_SNVs"],
            **scored,
        })
    return rows


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--paternal-blocks", required=True)
    parser.add_argument("--maternal-blocks", required=True)
    parser.add_argument("--paternal-mismatch", required=True)
    parser.add_argument("--maternal-mismatch", required=True)
    parser.add_argument("--min-num-het", type=int, default=500,
                        help="Ignore blocks with fewer het SNVs (default: 500)")
    parser.add_argument("--max-concordance", type=float, default=0.98,
                        help="Only consider blocks below this concordance, i.e. that "
                             "actually contain a switch (default: 0.98)")
    parser.add_argument("--same-site-kb", type=int, default=500,
                        help="Flag a candidate as a likely kid phase-switch (not a "
                             "meiotic crossover) when the other parent has a clean "
                             "breakpoint within this many kb (default: 500)")
    parser.add_argument("--top", type=int, default=10,
                        help="Print this many cleanest candidates (default: 10)")
    parser.add_argument("--out", default=None,
                        help="Optional TSV path for the full ranked candidate table")
    parser.add_argument("--igv-bed", default=None,
                        help="Optional BED path: candidate breakpoint intervals as an "
                             "IGV navigation track")
    parser.add_argument("--igv-batch", default=None,
                        help="Optional IGV batch-script path with a zoomed `goto` per "
                             "candidate (cleanest-first)")
    parser.add_argument("--igv-window", type=int, default=100_000,
                        help="Half-width (bp) of the --igv-batch goto window around each "
                             "estimated breakpoint (default: 100,000)")
    parser.add_argument("--snapshot-dir", default=None,
                        help="If set with --igv-batch, snapshot a PNG per candidate here")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s - %(levelname)s - %(filename)s - %(message)s",
                        datefmt="%Y-%m-%d %H:%M:%S")
    logger = logging.getLogger(__name__)

    rows = (
        _analyze_track(args.paternal_blocks, args.paternal_mismatch, "paternal",
                       args.min_num_het, args.max_concordance)
        + _analyze_track(args.maternal_blocks, args.maternal_mismatch, "maternal",
                         args.min_num_het, args.max_concordance)
    )
    if not rows:
        logger.info("No candidate blocks passed the filters.")
        return

    df = pl.DataFrame(rows)

    # Cross-parent flag: a clean breakpoint in BOTH parents at ~the same kid
    # coordinate is a kid phase-switch, not a meiotic crossover.
    other = {"paternal": "maternal", "maternal": "paternal"}
    flags = []
    for r in df.iter_rows(named=True):
        twin = df.filter(
            (pl.col("parent") == other[r["parent"]])
            & (pl.col("chrom") == r["chrom"])
            & ((pl.col("break_pos") - r["break_pos"]).abs() <= args.same_site_kb * 1000)
        )
        flags.append("BOTH-parents(phase-switch?)" if len(twin) else "single-parent")
    df = df.with_columns(crossover_kind=pl.Series(flags))

    df = df.sort(["cleanliness", "num_het_SNVs"], descending=[True, True])

    logger.info(f"{len(df)} candidate crossover block(s); cleanest first:")
    with pl.Config(tbl_rows=args.top, tbl_cols=-1, fmt_str_lengths=60):
        print(df.head(args.top))

    if args.out:
        df.write_csv(args.out, separator="\t")
        logger.info(f"Wrote ranked candidate table to '{args.out}'")
    if args.igv_bed:
        write_igv_bed(df, args.igv_bed)
        logger.info(f"Wrote IGV navigation track to '{args.igv_bed}'")
    if args.igv_batch:
        write_igv_batch(df, args.igv_batch, args.igv_window, snapshot_dir=args.snapshot_dir)
        logger.info(f"Wrote IGV batch script to '{args.igv_batch}'")

    best = df.filter(pl.col("crossover_kind") == "single-parent").head(1)
    if len(best):
        b = best.row(0, named=True)
        logger.info(
            f"Cleanest single-parent candidate: {b['parent']} {b['break_interval']} "
            f"(cleanliness={b['cleanliness']}, concordance={b['concordance']}, "
            f"num_het={b['num_het_SNVs']}). Confirm the sharp mismatch-density "
            f"transition in IGV before recutting the dev set around it."
        )


if __name__ == "__main__":
    main()
