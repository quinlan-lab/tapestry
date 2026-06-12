"""Regression check on produced hap-map-blocks BEDs: no duplicate blocks.

The triplicate-row bug emitted byte-identical hap-map blocks (e.g. the same
chr14:16089625-16091338 "A,1.0,3" row three times) because a SNV in a nested
phase-set region was fanned out across overlapping phase-block spans. The fix
assigns each SNV to exactly one phase set, so no two output rows can describe
the same interval. This script enforces that invariant on the actual output
BEDs after a (regression) run on CHPC:

    PYTHONPATH=src:src/util .venv/bin/python src/util/check_hap_map_blocks.py \\
        <out>/NA12878.hap-map-blocks.paternal.sorted.bed.gz \\
        <out>/NA12878.hap-map-blocks.maternal.sorted.bed.gz

Exit code is non-zero if any (chrom, start, end) appears more than once.

NOTE on nesting: genuinely nested / overlapping blocks are NOT a bug. A short
phase set in one individual, intersected with a chromosome-spanning block in the
other, legitimately yields a small hap-map block nested inside the big one. So
overlaps are reported for information only and do not fail the check; only
exact-duplicate intervals (the fan-out signature) fail.
"""

import argparse
import logging
import sys

import polars as pl


def load_blocks(bed_gz):
    return pl.read_csv(
        bed_gz, separator="\t", has_header=False,
        new_columns=["chrom", "start", "end", "name"],
    )


def duplicate_blocks(df):
    return (
        df.group_by(["chrom", "start", "end"]).len()
        .filter(pl.col("len") > 1)
        .sort("len", descending=True)
    )


def overlapping_blocks(df):
    """Rows whose interval overlaps the previous interval on the same chrom.

    Informational only. Sorting by (chrom, start) and flagging start < running
    max end catches both nesting and partial overlap.
    """
    d = df.sort(["chrom", "start", "end"])
    d = d.with_columns(
        prev_max_end=pl.col("end").cum_max().shift(1).over("chrom"),
    )
    return d.filter(
        pl.col("prev_max_end").is_not_null() & (pl.col("start") < pl.col("prev_max_end"))
    )


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("beds", nargs="+", help="hap-map-blocks BED files (.bed or .bed.gz)")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s - %(levelname)s - %(filename)s - %(message)s",
                        datefmt="%Y-%m-%d %H:%M:%S")
    logger = logging.getLogger(__name__)

    any_dup = False
    for bed in args.beds:
        df = load_blocks(bed)
        dups = duplicate_blocks(df)
        overlaps = overlapping_blocks(df)
        logger.info(
            f"{bed}: {len(df)} blocks, {len(dups)} duplicated interval(s), "
            f"{len(overlaps)} overlapping/nested interval(s) (informational)"
        )
        if len(dups):
            any_dup = True
            logger.error(f"{bed}: DUPLICATE blocks (fan-out signature):")
            print(dups)

    if any_dup:
        logger.error("FAIL: duplicate hap-map blocks present.")
        sys.exit(1)
    logger.info("PASS: no duplicate hap-map blocks.")


if __name__ == "__main__":
    main()
