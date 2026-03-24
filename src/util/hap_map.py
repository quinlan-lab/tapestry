import numpy as np
import polars as pl

from shell import shell
from write_data import write_bed


def extract_bit_vector(l):
    return np.array([int(x) for x in l], dtype=np.uint8)


def write_hap_map_blocks(df_hap_map, uid, parental, output_dir):
    """Write hap-map blocks BED for IGV visualization."""
    df_blocks = df_hap_map.select([
        pl.col("chrom"),
        pl.col("start"),
        pl.col("end"),
        pl.col(f"{parental}_haplotype"),
    ])
    write_bed(output_dir, df_blocks, f"{uid}.hap-map-blocks.{parental}")

    cmd = (
        f'cat {output_dir}/{uid}.hap-map-blocks.{parental}.bed'
        f' | src/util/sort-compress-index-bed'
        f' --name {output_dir}/{uid}.hap-map-blocks.{parental}'
    )
    shell(cmd)
    shell(f'rm {output_dir}/{uid}.hap-map-blocks.{parental}.bed')
