import numpy as np
import polars as pl
import pysam
from pathlib import Path

from write_data import write_bed


def extract_bit_vector(l):
    return np.array([int(x) for x in l], dtype=np.uint8)


def write_hap_map_blocks(df_hap_map, uid, parental, output_dir):
    """Write hap-map blocks BED for IGV visualization.

    The name column combines the assigned haplotype, the block's concordance,
    and the number of het SNVs in the block, e.g. "B,0.87,30". Concordance is
    the fraction of het SNVs in the block whose kid/parent bit-vector entries
    agree with the assigned haplotype. A block-wide concordance well below 1.0
    flags a likely recombination within the block (the kid switches parental
    haplotype partway through), which is expected when pedMEC phase blocks span
    an entire chromosome (~1 recombination per chromosome). The SNV count gives
    the weight behind that concordance: a low concordance over many SNVs is far
    more convincing than the same value over a handful.
    """
    df_blocks = (
        df_hap_map
        .with_columns(
            name=pl.col(f"{parental}_haplotype")
            + pl.lit(",")
            + pl.col("haplotype_concordance").round(2).cast(pl.String)
            + pl.lit(",")
            + pl.col("num_het_SNVs").cast(pl.String),
        )
        .select([
            pl.col("chrom"),
            pl.col("start"),
            pl.col("end"),
            pl.col("name"),
        ])
    )
    write_bed(output_dir, df_blocks, f"{uid}.hap-map-blocks.{parental}")

    plain = f"{output_dir}/{uid}.hap-map-blocks.{parental}.bed"
    compressed = f"{output_dir}/{uid}.hap-map-blocks.{parental}.bed.gz"
    pysam.tabix_compress(plain, compressed, force=True)
    pysam.tabix_index(compressed, preset="bed", force=True)
    Path(plain).unlink()
