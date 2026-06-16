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

Sex chromosomes and chrM are excluded by default (--exclude-chroms
chrX,chrY,chrM): for a female proband there is no chrY, and the paternal chrX
is transmitted without meiotic recombination outside the PARs, so candidates
there are artifacts. Pass --exclude-chroms '' to keep all chromosomes.

IGV confirmation
----------------
The printed table is cleanliness-ranked and its break_interval column is a
ready-to-paste chrom:lo-hi locus, so confirm a candidate by pasting that locus
into IGV alongside the bit-vector-sites-mismatches.{paternal,maternal}.vcf.gz
tracks. There the mismatch density should show a SHARP transition from sparse
ticks (kid agrees) to a dense contiguous run to the block edge (past the
crossover), appearing in only ONE parent's track (a true meiotic crossover); a
transition at the same kid coordinate in BOTH tracks is a kid phase-switch
artifact.

Usage:
    PYTHONPATH=src:src/util .venv/bin/python src/util/find_clean_single_crossover.py \\
        --paternal-blocks    <out>/NA12878.hap-map-blocks.paternal.sorted.bed.gz \\
        --maternal-blocks    <out>/NA12878.hap-map-blocks.maternal.sorted.bed.gz \\
        --paternal-mismatch  <out>/NA12878.bit-vector-sites-mismatches.paternal.bed \\
        --maternal-mismatch  <out>/NA12878.bit-vector-sites-mismatches.maternal.bed \\
        --out candidates.tsv
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
    # Fraction of THIS block's mismatches that fall in each predicted minority
    # interval -- i.e. the realised value of "frac" for each candidate side.
    # frac_right: share landing in (c, 1] (minority anchored to the right end).
    # frac_left:  share landing in [0, f) (minority anchored to the left end).
    # Both intervals have length f, so uniformly scattered (phase-switch noise)
    # mismatches give frac ~= f; a clean crossover concentrates ~all of them in
    # the true side, giving frac ~= 1. The larger of the two names the side.
    frac_right = sum(1 for x in u if x > c) / m       # minority on the right
    frac_left = sum(1 for x in u if x < f) / m        # minority on the left

    # break_pos: the crossover position PREDICTED FROM CONCORDANCE ALONE. A single
    #   clean crossover leaves a contiguous minority segment of length f, whose
    #   inner edge sits a distance f from the end it is anchored to -- coordinate
    #   c (= 1 - f) when the minority is on the right, f when it is on the left.
    # inner: the innermost OBSERVED minority mismatch (the one nearest the majority
    #   side). break_pos and inner are then bracketed (see lo, hi below) so the
    #   reported interval spans the disagreement between the concordance-predicted
    #   cut and where the mismatches actually begin; zoom there in IGV for the true
    #   flanking het coordinates. (default=break_pos guards the no-mismatch side.)
    if frac_right >= frac_left:
        side, frac = "right", frac_right
        break_pos = start + round(c * L)
        inner = min((p for p, x in zip(mm_pos, u) if x > c), default=break_pos)
    else:
        side, frac = "left", frac_left
        break_pos = start + round(f * L)
        inner = max((p for p, x in zip(mm_pos, u) if x < f), default=break_pos)

    # frac has two reference values, both following from uniform het-SNV density:
    #   * Random scatter (phase-switch noise) -> frac ~= f. Uniformly scattered
    #     mismatches fill any sub-interval in proportion to its length, and both
    #     predicted sides have length f ((c, 1] and [0, f)), so neither side is
    #     favoured: frac = max(f, f) ~= f.
    #   * Perfectly clean crossover -> frac ~= 1. All mismatches lie in the one
    #     true minority segment, which IS one of those two intervals, so that
    #     side captures ~all of them: frac = max(...) ~= 1.
    # Rescaling maps f -> 0 and 1 -> 1; subtracting f removes the baseline a
    # length-f interval picks up by chance, so the score measures concentration
    # BEYOND random scatter (low-concordance blocks, with larger f, must clear a
    # larger baseline).
    cleanliness = max(0.0, min(1.0, (frac - f) / (1.0 - f)))
    lo, hi = sorted((break_pos, inner))
    return {
        "cleanliness": round(cleanliness, 3),
        "side": side,
        "break_pos": break_pos,
        "break_interval": f"{block['chrom']}:{lo}-{hi}",
        "num_mismatch": m,
    }


def _analyze_track(blocks_bed, mismatch_bed, parent, min_num_het, max_concordance,
                   exclude_chroms):
    blocks = _load_blocks(blocks_bed)
    mism = _load_mismatch(mismatch_bed)
    rows = []
    for block in blocks.iter_rows(named=True):
        if block["chrom"] in exclude_chroms:
            continue
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
    parser.add_argument("--exclude-chroms", default="chrX,chrY,chrM",
                        help="Comma-separated chromosomes to skip. Default 'chrX,chrY,chrM' "
                             "drops sex-chromosome and mitochondrial artifacts, which is "
                             "correct for a female proband (no chrY; the paternal chrX is "
                             "transmitted without meiotic recombination outside the PARs). "
                             "Pass '' to keep all chromosomes.")
    parser.add_argument("--out", default=None,
                        help="Optional TSV path for the full ranked candidate table")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s - %(levelname)s - %(filename)s - %(message)s",
                        datefmt="%Y-%m-%d %H:%M:%S")
    logger = logging.getLogger(__name__)

    exclude_chroms = {c for c in args.exclude_chroms.split(",") if c}
    if exclude_chroms:
        logger.info(f"Excluding chromosomes: {', '.join(sorted(exclude_chroms))}")
    rows = (
        _analyze_track(args.paternal_blocks, args.paternal_mismatch, "paternal",
                       args.min_num_het, args.max_concordance, exclude_chroms)
        + _analyze_track(args.maternal_blocks, args.maternal_mismatch, "maternal",
                         args.min_num_het, args.max_concordance, exclude_chroms)
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
