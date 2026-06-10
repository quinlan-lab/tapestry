import argparse
import logging
from pathlib import Path

import polars as pl


def scan_blocks(bed_gz, threshold):
    """Scan one hap-map-blocks BED for blocks whose concordance is below threshold.

    The BED is written by write_hap_map_blocks() in src/util/hap_map.py with
    columns: chrom, start, end, name, where name encodes the assigned
    haplotype, the block-wide concordance, and the number of het SNVs in the
    block as "<hap>,<concordance>,<num_het_SNVs>", e.g. "B,0.87,30".
    Concordance is the fraction of het SNVs in the block whose kid/parent
    bit-vector entries agree with the assigned haplotype.

    A concordance well below 1.0 means the kid's allele sequence matches the
    assigned parental haplotype for only part of the block, which is the
    signature of a recombination within the block (the kid switches parental
    haplotype partway through). This is expected when pedMEC phase blocks span
    an entire chromosome (~1 recombination per chromosome). The num_het_SNVs
    field gives the weight behind the concordance: a low concordance over many
    SNVs is far more convincing than the same value over a handful.

    Returns a DataFrame of below-threshold blocks, sorted by concordance
    (most discordant first), with an added 'length' column.
    """
    df = (
        pl.read_csv(
            bed_gz,
            separator="\t",
            has_header=False,
            new_columns=["chrom", "start", "end", "name"],
        )
        .with_columns(
            haplotype=pl.col("name").str.split(",").list.get(0),
            concordance=pl.col("name").str.split(",").list.get(1).cast(pl.Float64),
            num_het_SNVs=pl.col("name").str.split(",").list.get(2).cast(pl.Int64),
        )
        .with_columns(
            length=pl.col("end") - pl.col("start"),
        )
        .filter(pl.col("concordance") < threshold)
        .select(["chrom", "start", "end", "haplotype", "concordance", "num_het_SNVs", "length"])
        .sort("concordance")
    )
    return df


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Scan hap-map-blocks BED(s) for blocks whose concordance is below a "
            "threshold. Low concordance flags a likely recombination within the "
            "block. Input BEDs are produced by phase_meth_to_parent_haps.py, e.g. "
            "NA12878.hap-map-blocks.paternal.sorted.bed.gz"
        )
    )
    parser.add_argument(
        "beds",
        nargs="+",
        help="One or more hap-map-blocks BED files (.bed or .bed.gz)",
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=0.99,
        help="Report blocks with concordance strictly below this value (default: 0.99)",
    )
    parser.add_argument(
        "--out",
        default=None,
        help="Optional TSV path to write the combined results to (in addition to logging)",
    )
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(filename)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    logger = logging.getLogger(__name__)

    results = []
    for bed in args.beds:
        df = scan_blocks(bed, args.threshold)
        logger.info(
            f"{bed}: {len(df)} block(s) with concordance < {args.threshold}"
        )
        if len(df):
            print(df.with_columns(source=pl.lit(Path(bed).name)))
        results.append(df.with_columns(source=pl.lit(Path(bed).name)))

    df_all = pl.concat(results).sort("concordance") if results else pl.DataFrame()
    logger.info(f"Total: {len(df_all)} candidate recombination block(s)")

    if args.out and len(df_all):
        df_all.write_csv(args.out, separator="\t")
        logger.info(f"Wrote candidate blocks to '{args.out}'")


if __name__ == "__main__":
    main()
