import numpy as np
import polars as pl
import bioframe as bf

from util.shell import shell
from util.write_data import write_bed


def get_hap_map_blocks(
    df_blocks_kid: pl.DataFrame,
    df_blocks_dad: pl.DataFrame,
    df_blocks_mom: pl.DataFrame,
) -> pl.DataFrame:
    """
    Compute the intersection of phase blocks across all three individuals.
    Each resulting interval is a hap-map block.
    """
    # Intersect kid and dad blocks
    df_kd = pl.from_pandas(bf.overlap(
        df_blocks_kid.select(["chrom", "start", "end"]).to_pandas(),
        df_blocks_dad.select(["chrom", "start", "end"]).to_pandas(),
        how='inner',
        suffixes=('_kid', '_dad'),
    ))

    # Compute the intersection interval
    df_kd = df_kd.with_columns(
        pl.max_horizontal("start_kid", "start_dad").alias("start"),
        pl.min_horizontal("end_kid", "end_dad").alias("end"),
    ).filter(pl.col("start") < pl.col("end"))

    # Now intersect with mom blocks
    df_kdm = pl.from_pandas(bf.overlap(
        df_kd.select(["chrom_kid", "start", "end"]).rename({"chrom_kid": "chrom"}).to_pandas(),
        df_blocks_mom.select(["chrom", "start", "end"]).to_pandas(),
        how='inner',
        suffixes=('_kd', '_mom'),
    ))

    df_hap_map_blocks = (
        df_kdm
        .with_columns(
            pl.max_horizontal("start_kd", "start_mom").alias("start"),
            pl.min_horizontal("end_kd", "end_mom").alias("end"),
        )
        .filter(pl.col("start") < pl.col("end"))
        .select([
            pl.col("chrom_kd").alias("chrom"),
            "start",
            "end",
        ])
        .sort(["chrom", "start", "end"])
    )

    return df_hap_map_blocks


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


def assign_snvs_to_hap_map_blocks(df_phasing: pl.DataFrame, df_hap_map_blocks: pl.DataFrame) -> pl.DataFrame:
    """
    Overlap het SNVs with hap-map blocks so each SNV is assigned to a block.
    """
    df = pl.from_pandas(bf.overlap(
        df_phasing.to_pandas(),
        df_hap_map_blocks.to_pandas(),
        how='inner',
        suffixes=('', '_block'),
    ))
    df = df.rename({
        "chrom_block": "chrom_hap_map_block",
        "start_block": "start_hap_map_block",
        "end_block": "end_hap_map_block",
    })
    return df


def get_hap_map(
    df_kid: pl.DataFrame,
    df_dad: pl.DataFrame,
    df_mom: pl.DataFrame,
    df_blocks_kid, df_blocks_dad, df_blocks_mom,
):
    """
    Build the haplotype map by comparing bit vectors.

    Labeling convention:
        A = dad's hap1, B = dad's hap2 (fixed)
        C = mom's hap1, D = mom's hap2 (fixed)

    For each hap-map block, we determine which of A or B the kid's paternal
    haplotype (hap1) matches, and which of C or D the kid's maternal
    haplotype (hap2) matches.

    We compare the kid's paternal bit-vector (allele_hap1, since whatshap
    phases hap1=paternal) to the dad's hap1 bit-vector:
        - similarity > 0.5 → kid's hap1 matches dad's hap1 → paternal_haplotype = "A"
        - similarity <= 0.5 → kid's hap1 matches dad's hap2 → paternal_haplotype = "B"
    Similarly for maternal:
        - similarity > 0.5 → kid's hap2 matches mom's hap1 → maternal_haplotype = "C"
        - similarity <= 0.5 → kid's hap2 matches mom's hap2 → maternal_haplotype = "D"

    Returns:
        df_hap_map with columns:
            chrom, start, end,
            paternal_haplotype ("A" or "B"),
            maternal_haplotype ("C" or "D"),
            paternal_concordance, maternal_concordance,
            num_het_SNVs_pat, num_het_SNVs_mat
    """

    df_hap_map_blocks = get_hap_map_blocks(df_blocks_kid, df_blocks_dad, df_blocks_mom)

    # Assign kid, dad, mom SNVs to hap-map blocks
    df_kid_in_blocks = assign_snvs_to_hap_map_blocks(df_kid, df_hap_map_blocks)
    df_dad_in_blocks = assign_snvs_to_hap_map_blocks(df_dad, df_hap_map_blocks)
    df_mom_in_blocks = assign_snvs_to_hap_map_blocks(df_mom, df_hap_map_blocks)

    # Join kid with dad on shared het SNV positions within the same hap-map block
    df_kid_dad = (
        df_kid_in_blocks
        .join(
            df_dad_in_blocks,
            on=["chrom", "start", "end", "REF", "ALT",
                "start_hap_map_block", "end_hap_map_block"],
            how="inner",
            suffix="_dad",
        )
    )

    # Join kid with mom on shared het SNV positions within the same hap-map block
    df_kid_mom = (
        df_kid_in_blocks
        .join(
            df_mom_in_blocks,
            on=["chrom", "start", "end", "REF", "ALT",
                "start_hap_map_block", "end_hap_map_block"],
            how="inner",
            suffix="_mom",
        )
    )

    # Group by hap-map block and build bit vectors for paternal comparison
    df_pat_grouped = (
        df_kid_dad
        .group_by(["chrom", "start_hap_map_block", "end_hap_map_block"])
        .agg([
            pl.col("allele_hap1").implode().alias("kid_pat_seq"),
            pl.col("allele_hap1_dad").implode().alias("dad_hap1_seq"),
        ])
    )

    # Group by hap-map block and build bit vectors for maternal comparison
    df_mat_grouped = (
        df_kid_mom
        .group_by(["chrom", "start_hap_map_block", "end_hap_map_block"])
        .agg([
            pl.col("allele_hap2").implode().alias("kid_mat_seq"),
            pl.col("allele_hap1_mom").implode().alias("mom_hap1_seq"),
        ])
    )

    # Process paternal comparisons
    pat_records = []
    for row in df_pat_grouped.iter_rows(named=True):
        kid_pat = np.array([int(x) for x in row["kid_pat_seq"][0]], dtype=np.uint8)
        dad_hap1 = np.array([int(x) for x in row["dad_hap1_seq"][0]], dtype=np.uint8)

        n = len(kid_pat)
        match_hap1 = np.sum(kid_pat == dad_hap1)
        similarity_to_dad_hap1 = match_hap1 / n

        if similarity_to_dad_hap1 > 0.5:
            pat_haplotype = "A"
            pat_concordance = similarity_to_dad_hap1
        else:
            pat_haplotype = "B"
            pat_concordance = 1.0 - similarity_to_dad_hap1

        pat_records.append({
            "chrom": row["chrom"],
            "start_hap_map_block": row["start_hap_map_block"],
            "end_hap_map_block": row["end_hap_map_block"],
            "paternal_haplotype": pat_haplotype,
            "paternal_concordance": pat_concordance,
            "num_het_SNVs_pat": n,
        })

    # Process maternal comparisons
    mat_records = []
    for row in df_mat_grouped.iter_rows(named=True):
        kid_mat = np.array([int(x) for x in row["kid_mat_seq"][0]], dtype=np.uint8)
        mom_hap1 = np.array([int(x) for x in row["mom_hap1_seq"][0]], dtype=np.uint8)

        n = len(kid_mat)
        match_hap1 = np.sum(kid_mat == mom_hap1)
        similarity_to_mom_hap1 = match_hap1 / n

        if similarity_to_mom_hap1 > 0.5:
            mat_haplotype = "C"
            mat_concordance = similarity_to_mom_hap1
        else:
            mat_haplotype = "D"
            mat_concordance = 1.0 - similarity_to_mom_hap1

        mat_records.append({
            "chrom": row["chrom"],
            "start_hap_map_block": row["start_hap_map_block"],
            "end_hap_map_block": row["end_hap_map_block"],
            "maternal_haplotype": mat_haplotype,
            "maternal_concordance": mat_concordance,
            "num_het_SNVs_mat": n,
        })

    df_pat = pl.DataFrame(pat_records)
    df_mat = pl.DataFrame(mat_records)

    # Join paternal and maternal results
    df_hap_map = (
        df_pat
        .join(
            df_mat,
            on=["chrom", "start_hap_map_block", "end_hap_map_block"],
            how="full",
            coalesce=True,
        )
        .rename({
            "start_hap_map_block": "start",
            "end_hap_map_block": "end",
        })
        .select([
            "chrom", "start", "end",
            "paternal_haplotype", "maternal_haplotype",
            "paternal_concordance", "maternal_concordance",
            "num_het_SNVs_pat", "num_het_SNVs_mat",
        ])
        .sort(["chrom", "start", "end"])
    )

    return df_hap_map


