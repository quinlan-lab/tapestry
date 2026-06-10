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


def write_igv_bed(df_all, path):
    """Write candidate blocks as a BED navigation track for IGV.

    Load this track in IGV, click to select it, then use Ctrl-F / Ctrl-B
    (next/previous feature) to step through candidate blocks in genomic order
    without copy-pasting coordinates. The feature name carries the same
    "<hap>,<concordance>,<num_het_SNVs>" label as the source hap-map-blocks
    BED, plus the source file, so you can tell paternal from maternal at a
    glance. The track is sorted by position (BED feature navigation walks
    genomic order regardless of input order).
    """
    df_bed = (
        df_all
        .with_columns(
            name=pl.col("haplotype")
            + pl.lit(",")
            + pl.col("concordance").cast(pl.String)
            + pl.lit(",")
            + pl.col("num_het_SNVs").cast(pl.String)
            + pl.lit("[")
            + pl.col("source").str.replace(".hap-map-blocks.", ".").str.replace(".sorted.bed.gz", "")
            + pl.lit("]"),
        )
        .select(["chrom", "start", "end", "name"])
        .sort(["chrom", "start", "end"])
    )
    df_bed.write_csv(path, separator="\t", include_header=False)


def write_igv_batch(df_all, path, snapshot_dir=None):
    """Write an IGV batch script that visits each candidate block in turn.

    Run via IGV's Tools > Run Batch Script, or `igv.sh -b <path>`. Each
    candidate produces a `goto` to the block interval (chromosome-spanning for
    pedMEC blocks, so this positions the chromosome with whatever mismatch/
    block tracks you already have loaded). If snapshot_dir is given, a PNG is
    written per candidate so you can review them outside IGV. Candidates are
    ordered most-discordant-first to match the scan's ranking.
    """
    lines = []
    if snapshot_dir:
        lines.append(f"snapshotDirectory {snapshot_dir}")
    for row in df_all.iter_rows(named=True):
        lines.append(f"goto {row['chrom']}:{row['start']}-{row['end']}")
        if snapshot_dir:
            tag = row["source"].replace(".hap-map-blocks.", ".").replace(".sorted.bed.gz", "")
            lines.append(f"snapshot {tag}.{row['chrom']}_{row['start']}_{row['end']}.png")
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=(
            "Scan hap-map-blocks BED(s) for blocks whose concordance is below a "
            "threshold. Low concordance flags a likely recombination within the "
            "block. Input BEDs are produced by phase_meth_to_parent_haps.py, e.g. "
            "NA12878.hap-map-blocks.paternal.sorted.bed.gz"
        ),
        epilog=(
            "IGV workflow for localizing a crossover:\n"
            "  1. Generate a navigation track and load it in IGV alongside the\n"
            "     bit-vector-sites-mismatches.{paternal,maternal}.vcf.gz tracks:\n"
            "       --igv-bed candidates.igv.bed\n"
            "     Select the track and use Ctrl-F / Ctrl-B to step through the\n"
            "     candidates in genomic order (or --igv-batch for scripted goto/\n"
            "     snapshot stepping).\n"
            "  2. The mismatch VCF holds the het sites where the kid disagrees\n"
            "     with the haplotype assigned to the whole (chromosome-spanning)\n"
            "     block. A single crossover therefore shows up not as scattered\n"
            "     noise but as a SHARP TRANSITION from sparse mismatch ticks\n"
            "     (the majority segment, where the kid agrees) to a dense\n"
            "     contiguous run extending to the block edge (the minority\n"
            "     segment, past the crossover). That transition edge is the\n"
            "     breakpoint; zoom in to read the flanking het-SNV coordinates,\n"
            "     which bound the breakpoint interval.\n"
            "  3. Confirm it is a true meiotic crossover, not a phasing artifact:\n"
            "     compare the paternal vs maternal mismatch tracks. A real\n"
            "     crossover flips only ONE parent's transmitted haplotype, so the\n"
            "     density transition appears in only one track. A kid phase-switch\n"
            "     error flips both at the SAME kid coordinate. Prefer candidates\n"
            "     where only one track transitions.\n"
            "\n"
            "The num_het_SNVs field weights the concordance: a low value over\n"
            "many SNVs is far more convincing than the same value over a handful."
        ),
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
    parser.add_argument(
        "--igv-bed",
        default=None,
        help=(
            "Optional BED path to write candidate blocks as an IGV navigation "
            "track. Load it in IGV and use Ctrl-F/Ctrl-B to step through "
            "candidates in genomic order."
        ),
    )
    parser.add_argument(
        "--igv-batch",
        default=None,
        help=(
            "Optional IGV batch-script path with a `goto` per candidate "
            "(most-discordant-first). Run via Tools > Run Batch Script or "
            "`igv.sh -b <path>`."
        ),
    )
    parser.add_argument(
        "--snapshot-dir",
        default=None,
        help=(
            "If set with --igv-batch, the batch script also snapshots a PNG per "
            "candidate into this directory."
        ),
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

    if args.igv_bed and len(df_all):
        write_igv_bed(df_all, args.igv_bed)
        logger.info(f"Wrote IGV navigation track to '{args.igv_bed}'")

    if args.igv_batch and len(df_all):
        write_igv_batch(df_all, args.igv_batch, snapshot_dir=args.snapshot_dir)
        logger.info(f"Wrote IGV batch script to '{args.igv_batch}'")


if __name__ == "__main__":
    main()
