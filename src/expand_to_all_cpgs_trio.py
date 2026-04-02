import argparse
import logging
from pathlib import Path
from typing import Literal

import polars as pl
import bioframe as bf

from read_data import read_bed_and_header
from get_meth_hap1_hap2 import read_meth_level

REFERENCE_GENOME = "hg38"
CPG_SITE_MISMATCH_SITE_DISTANCE = 50 # bp


def read_meth_unphased_individual(bed_meth_count, bed_meth_model, individual):
    """Read count-based and model-based unphased methylation for one individual.

    Returns a dataframe with columns:
        chrom, start, end,
        total_read_count_{individual},
        methylation_level_{individual}_count,
        methylation_level_{individual}_model
    """
    df_count = (
        read_meth_level(bed_meth_count, pb_cpg_tool_mode='count')
        .rename({
            'chromosome': 'chrom',
            'total_read_count': f'total_read_count_{individual}_count',
            'methylation_level': f'methylation_level_{individual}_count',
        })
    )
    df_model = (
        read_meth_level(bed_meth_model, pb_cpg_tool_mode='model')
        .rename({
            'chromosome': 'chrom',
            'total_read_count': f'total_read_count_{individual}_model',
            'methylation_level': f'methylation_level_{individual}_model',
        })
    )
    df = df_count.join(
        df_model,
        on=['chrom', 'start', 'end'],
        how="full",
        coalesce=True,
    )

    # Assert count and model total_read_counts are identical, then collapse
    count_col = f'total_read_count_{individual}_count'
    model_col = f'total_read_count_{individual}_model'
    if df[count_col].equals(df[model_col]):
        df = (
            df
            .drop([model_col])
            .rename({count_col: f'total_read_count_{individual}'})
        )
    else:
        raise ValueError(
            f'total_read_count for {individual}: '
            f'{count_col} is not identical to {model_col}'
        )

    return df


def read_meth_unphased_trio(
    bed_meth_count_kid, bed_meth_model_kid,
    bed_meth_count_dad, bed_meth_model_dad,
    bed_meth_count_mom, bed_meth_model_mom,
):
    """Read unphased methylation for kid, dad, and mom, and join into one dataframe.

    Returns a dataframe with columns:
        chrom, start, end,
        total_read_count_kid, methylation_level_kid_count, methylation_level_kid_model,
        total_read_count_dad, methylation_level_dad_count, methylation_level_dad_model,
        total_read_count_mom, methylation_level_mom_count, methylation_level_mom_model
    """
    df_kid = read_meth_unphased_individual(bed_meth_count_kid, bed_meth_model_kid, 'kid')
    df_dad = read_meth_unphased_individual(bed_meth_count_dad, bed_meth_model_dad, 'dad')
    df_mom = read_meth_unphased_individual(bed_meth_count_mom, bed_meth_model_mom, 'mom')

    df = df_kid.join(
        df_dad,
        on=['chrom', 'start', 'end'],
        how="full",
        coalesce=True,
    )
    df = df.join(
        df_mom,
        on=['chrom', 'start', 'end'],
        how="full",
        coalesce=True,
    )

    return df


def read_meth_parent_phased(bed_meth_parent_phased):
    return read_bed_and_header(bed_meth_parent_phased)


def expand_meth_to_all_cpgs(df_all_cpgs_in_reference, df_meth_unphased, df_meth_parent_phased):
    df_all_cpgs = (
        df_all_cpgs_in_reference
        .join(
            df_meth_unphased,
            on=['chrom', 'start', 'end'],
            # "how=full" (and "coalesce=True") is REQUIRED
            # because there are records unique to df_all_cpgs_in_reference AND other records unique to df_meth_unphased
            how="full",
            coalesce=True,
        )
    )

    df = (
        df_all_cpgs
        .join(
            df_meth_parent_phased,
            on=['chrom', 'start', 'end'],
            # "how=left" is appropriate
            # because the set of sites in df_meth_parent_phased is a subset of the set of sites in df_all_cpgs
            how="left"
        )
    )

    return df


def _compute_proximity_one_side(df, bed_het_site_mismatches, side):
    """Compute proximity of each CpG site to the nearest mismatch site for one parental side.

    Args:
        df: DataFrame with all CpG sites
        bed_het_site_mismatches: path to mismatch BED file (paternal or maternal)
        side: 'pat' or 'mat'

    Returns:
        DataFrame with an added is_within_{distance}bp_of_mismatch_site_{side} column
    """
    df_sites_mismatch = read_bed_and_header(bed_het_site_mismatches)

    # Integer columns that bioframe's closest() will convert to float (due to NaN)
    int_columns = {
        "start_hap_map_block_pat": pl.Int64,
        "end_hap_map_block_pat": pl.Int64,
        "num_het_SNVs_in_dad": pl.Int64,
        "start_hap_map_block_mat": pl.Int64,
        "end_hap_map_block_mat": pl.Int64,
        "num_het_SNVs_in_mom": pl.Int64,
        "total_read_count_kid_pat": pl.Int64,
        "total_read_count_kid_mat": pl.Int64,
        "total_read_count_dad_A": pl.Int64,
        "total_read_count_dad_B": pl.Int64,
        "total_read_count_mom_C": pl.Int64,
        "total_read_count_mom_D": pl.Int64,
        "total_read_count_kid": pl.Int64,
        "total_read_count_dad": pl.Int64,
        "total_read_count_mom": pl.Int64,
    }

    df = (
        pl
        .from_pandas(
            bf.closest(
                df.to_pandas(),
                df_sites_mismatch.to_pandas()
            )
        )
        .cast(int_columns) # type: ignore 
        .with_columns(
            (pl.col('distance') < CPG_SITE_MISMATCH_SITE_DISTANCE)
            .alias(f'is_within_{CPG_SITE_MISMATCH_SITE_DISTANCE}bp_of_mismatch_site_{side}'),
        )
        .drop([
            "chrom_", "start_", "end_",
            "REF_", "ALT_",
            "distance"
        ])
    )

    return df


def compute_proximity_to_mismatched_heterozygous_sites(
    df_meth_parent_phased_all_cpgs,
    bed_het_site_mismatches_pat,
    bed_het_site_mismatches_mat,
):
    """Compute proximity of each CpG site to nearest paternal and maternal mismatch sites."""
    df = _compute_proximity_one_side(
        df_meth_parent_phased_all_cpgs,
        bed_het_site_mismatches_pat,
        'pat',
    )
    df = _compute_proximity_one_side(
        df,
        bed_het_site_mismatches_mat,
        'mat',
    )
    return df


def reduce_to_phasable_chromosomes(df):
    df = df.filter(
        (pl.col('chrom') != 'chrM') &
        (pl.col('chrom') != 'chrX') &
        (pl.col('chrom') != 'chrY')
    )
    for side in ['pat', 'mat']:
        col = f"is_within_{CPG_SITE_MISMATCH_SITE_DISTANCE}bp_of_mismatch_site_{side}"
        assert df[col].is_not_null().all()
    return df


def compute_fraction_of_cpgs_that_are_close_to_mismatches(df, logger=None):
    df = reduce_to_phasable_chromosomes(df)
    for side in ['pat', 'mat']:
        col = f"is_within_{CPG_SITE_MISMATCH_SITE_DISTANCE}bp_of_mismatch_site_{side}"
        fraction = df[col].mean()
        report = (
            f"Percentage of CpG sites (in reference and sample genome, and on phasable chroms) "
            f"that are within {CPG_SITE_MISMATCH_SITE_DISTANCE}bp of a {side}ernal heterozygous mismatch site: "
            f"{fraction*100:.3f}%"
        )
        if logger:
            logger.info(report)
        else:
            print(report)


def compute_fraction_of_cpgs_at_which_unphased_meth_is_reported(df, individual, mode, logger=None):
    df = reduce_to_phasable_chromosomes(df)
    fraction = (
        df
        .select(pl.col(f"methylation_level_{individual}_{mode}").is_not_null())
        .mean()
        .item()
    )
    report = (
        f"Percentage of CpG sites (in reference and sample genomes, and on phasable chroms) "
        f"at which {mode}-based unphased methylation is reported for {individual}: "
        f"{fraction*100:.2f}%"
    )
    if logger:
        logger.info(report)
    else:
        print(report)


def compute_fraction_of_cpgs_at_which_meth_is_phased_to_given_haplotype(df, haplotype, mode, logger=None):
    df = reduce_to_phasable_chromosomes(df)
    fraction = (
        df
        .select(pl.col(f"methylation_level_{haplotype}_{mode}").is_not_null())
        .mean()
        .item()
    )
    report = (
        f"Percentage of CpG sites (in reference and sample genomes, and on phasable chroms) "
        f"at which {mode}-based methylation is phased to {haplotype} haplotype: "
        f"{fraction*100:.2f}%"
    )
    if logger:
        logger.info(report)
    else:
        print(report)


def create_meth_non_null_expr(
    haplotype_1: str,
    haplotype_2: str,
    mode: str,
    logic: Literal["any", "all"],
    alias: str,
) -> pl.Expr:
    col_1 = pl.col(f"methylation_level_{haplotype_1}_{mode}")
    col_2 = pl.col(f"methylation_level_{haplotype_2}_{mode}")

    if logic == "any":
        expression = col_1.is_not_null() | col_2.is_not_null()
    elif logic == "all":
        expression = col_1.is_not_null() & col_2.is_not_null()
    else:
        raise ValueError("`logic` parameter must be 'any' or 'all'")

    return expression.alias(alias)


def compute_fraction_of_cpgs_at_which_meth_is_phased(df, haplotype_1, haplotype_2, mode, logic, alias, logger=None):
    df = reduce_to_phasable_chromosomes(df)
    expr = create_meth_non_null_expr(haplotype_1, haplotype_2, mode, logic, alias)
    fraction = (
        df
        .with_columns(expr)
        .select(pl.col(alias))
        .mean()
        .item()
    )
    report = (
        f"Percentage of CpG sites (in reference and sample genomes, and on phasable chroms) "
        f"at which {mode}-based methylation is phased to {alias}: "
        f"{fraction*100:.2f}%"
    )
    if logger:
        logger.info(report)
    else:
        print(report)


def compute_methylation_coverage_qc(df, logger=None):
    # Phasing can be partial even within a hap-map block: if coverage on one haplotype
    # is low but high on the other, methylation will be reported on one haplotype but
    # not the other. This is why we distinguish "any" (at least one haplotype phased)
    # from "all" (both haplotypes phased) below.
    for mode in ['count', 'model']:
        for haplotype in ['kid_pat', 'kid_mat']:
            compute_fraction_of_cpgs_at_which_meth_is_phased_to_given_haplotype(
                df, haplotype, mode, logger,
            )
        compute_fraction_of_cpgs_at_which_meth_is_phased(
            df,
            haplotype_1='kid_pat',
            haplotype_2='kid_mat',
            mode=mode,
            logic='any',
            alias='at least one parental haplotype (kid)',
            logger=logger,
        )
        compute_fraction_of_cpgs_at_which_meth_is_phased(
            df,
            haplotype_1='kid_pat',
            haplotype_2='kid_mat',
            mode=mode,
            logic='all',
            alias='both parental haplotypes (kid)',
            logger=logger,
        )
        compute_fraction_of_cpgs_at_which_unphased_meth_is_reported(
            df, 'kid', mode, logger,
        )


def write_combined_bigwig(bed_path, pb_cpg_tool_mode, uid, output_dir, logger=None):
    """
    Write a combined (unphased) bigwig file on a 0-1 scale from a pb-CpG-tools BED file.
    """
    df = read_meth_level(bed_path, pb_cpg_tool_mode)
    df_bed_graph = (
        df
        .filter(pl.col("methylation_level").is_not_null())
        .select([
            pl.col("chromosome").alias("chrom"),
            pl.col("start"),
            pl.col("end"),
            pl.col("methylation_level"),
        ])
        .to_pandas()
    )

    file_path = Path(output_dir) / f"{uid}.dna-methylation.combined.{pb_cpg_tool_mode}.{REFERENCE_GENOME}.bw"

    if file_path.exists() or file_path.is_symlink():
        file_path.unlink()

    bf.to_bigwig(
        df_bed_graph,
        bf.fetch_chromsizes(db=REFERENCE_GENOME),
        outpath=file_path,
        path_to_binary="bedGraphToBigWig",  # type: ignore
    )

    if logger:
        logger.info(f"Wrote combined {pb_cpg_tool_mode}-based methylation bigwig, ASSUMING {REFERENCE_GENOME}, to: '{file_path}'")


def main():
    parser = argparse.ArgumentParser(
        description='Write combined (unphased) methylation bigwig files for a trio'
    )
    parser.add_argument('--kid_id', required=True, help='Child sample ID')
    parser.add_argument('--dad_id', required=True, help='Father sample ID')
    parser.add_argument('--mom_id', required=True, help='Mother sample ID')

    parser.add_argument('--bed_meth_count_combined_kid', required=True)
    parser.add_argument('--bed_meth_model_combined_kid', required=True)
    parser.add_argument('--bed_meth_count_combined_dad', required=True)
    parser.add_argument('--bed_meth_model_combined_dad', required=True)
    parser.add_argument('--bed_meth_count_combined_mom', required=True)
    parser.add_argument('--bed_meth_model_combined_mom', required=True)

    parser.add_argument('--output_dir', required=True, help='Output directory')
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(filename)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    logger = logging.getLogger(__name__)

    logger.info(f"Starting '{__file__}'")
    logger.info("Args: %s", vars(args))

    Path(args.output_dir).mkdir(parents=True, exist_ok=True)

    logger.info("Writing combined (unphased) methylation bigwig files...")
    combined_specs = [
        (args.kid_id, args.bed_meth_count_combined_kid, args.bed_meth_model_combined_kid),
        (args.dad_id, args.bed_meth_count_combined_dad, args.bed_meth_model_combined_dad),
        (args.mom_id, args.bed_meth_count_combined_mom, args.bed_meth_model_combined_mom),
    ]
    for uid, bed_count, bed_model in combined_specs:
        write_combined_bigwig(bed_count, 'count', uid, args.output_dir, logger)
        write_combined_bigwig(bed_model, 'model', uid, args.output_dir, logger)

    logger.info(f"Done running '{__file__}'")


if __name__ == "__main__":
    main()
