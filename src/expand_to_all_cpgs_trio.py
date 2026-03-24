import polars as pl

from read_data import read_bed_and_header
from get_meth_hap1_hap2 import read_meth_level


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
