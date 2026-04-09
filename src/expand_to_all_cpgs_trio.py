import argparse
import logging
from pathlib import Path
from typing import Literal
import os

# polars hogs memory by default:
# https://quinlangroup.slack.com/archives/C02KEHXJ274/p1765314327680539
# but we can force aggressive memory release using the following configuration:
# https://github.com/pola-rs/polars/issues/23128#issuecomment-2978695783
# Note: this configuration must be executed prior to importing polars
JEMALLOC_CONFIG = "background_thread:true,dirty_decay_ms:1000,muzzy_decay_ms:1000,narenas:2"
os.environ["MALLOC_CONF"] = JEMALLOC_CONFIG
# Set this as well to cover specific Polars builds (e.g. PyPI wheels)
os.environ["_RJEM_MALLOC_CONF"] = JEMALLOC_CONFIG
import polars as pl

import bioframe as bf

from memory_profiler import profile

from read_data import read_bed_and_header
from write_data import write_methylation
from get_meth_hap1_hap2 import read_meth_level
from shell import shell
from util import read_all_cpgs_in_reference
from logging_util import report_size

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


def _allele_exprs_for_member(uid, role):
    """Return a list of polars expressions that extract allele_1 and allele_2 for one family member.

    Args:
        uid: Sample UID used to reference the intermediate gt_raw column (from the VCF).
        role: Role label ('kid', 'dad', or 'mom') used for the output column names.
    """
    gt_raw_col = f"gt_raw_{uid}"
    return [
        # Parse GT string: take first field before ':', normalize separator, split into list, extract alleles
        pl.col(gt_raw_col)
            .str.split(":")
            .list.get(0)
            .str.replace("/", "|")
            .str.split("|")
            .list.get(0)
            .alias(f"allele_1_{role}"),
        pl.col(gt_raw_col)
            .str.split(":")
            .list.get(0)
            .str.replace("/", "|")
            .str.split("|")
            .list.get(1)
            .alias(f"allele_2_{role}"),
    ]


def get_joint_called_variants(kid_uid, dad_uid, mom_uid, vcf_path):
    """Parse joint-called VCF once and extract SNV alleles for kid, dad, and mom.

    UIDs are used to select columns from the VCF; output columns use role-based
    names (kid, dad, mom) for consistency with the rest of the pipeline.

    Returns a dataframe with columns:
        chrom, start, end,
        allele_1_kid, allele_2_kid,
        allele_1_dad, allele_2_dad,
        allele_1_mom, allele_2_mom
    """
    return (
        pl.scan_csv(
            vcf_path,
            separator='\t',
            comment_prefix='##',
            has_header=True,
            null_values=['.'],
            ignore_errors=True,
        )
        .rename({"#CHROM": "chrom", "POS": "pos"})
        .select([
            pl.col("chrom"),
            pl.col("pos"),
            pl.col("REF"),
            pl.col("ALT"),
            pl.col(kid_uid).alias(f"gt_raw_{kid_uid}"),
            pl.col(dad_uid).alias(f"gt_raw_{dad_uid}"),
            pl.col(mom_uid).alias(f"gt_raw_{mom_uid}"),
        ])
        # Strict SNV Filter (Multi-allelic safe)
        .filter(
            (pl.col("REF").str.len_chars() == 1) &
            (pl.col("ALT").is_not_null()) &
            (
                pl
                .col("ALT")
                .str.split(",")
                .list.eval(pl.element().str.len_chars() == 1)
                .list.all()
            )
        )
        .with_columns([
            (pl.col("pos") - 1).alias("start"),
            pl.col("pos").alias("end"),
        ])
        # Extract alleles for each family member (UID for VCF access, role for output names)
        .with_columns(
            _allele_exprs_for_member(kid_uid, 'kid')
            + _allele_exprs_for_member(dad_uid, 'dad')
            + _allele_exprs_for_member(mom_uid, 'mom')
        )
        .select([
            "chrom", "start", "end",
            "allele_1_kid", "allele_2_kid",
            "allele_1_dad", "allele_2_dad",
            "allele_1_mom", "allele_2_mom",
        ])
        .collect()
    )


def label_with_variants(df_meth_parent_phased_all_cpgs, df_joint_called_variants):
    """Overlap CpG dinucleotides with joint-called SNVs for the trio.

    Returns the input dataframe augmented with variant overlap information.
    """
    df_meth_dinucleotides = (
        df_meth_parent_phased_all_cpgs
        .with_columns((pl.col("end") + 1)) # enlarge interval from C nucleotide to CG dinucleotide
    )
    df_meth_dinucleotides_pd = df_meth_dinucleotides.to_pandas()

    df_joint_called_variants_pd = df_joint_called_variants.to_pandas()

    df = bf.overlap(
        df_meth_dinucleotides_pd,
        df_joint_called_variants_pd,
        how='left',
        suffixes=('', '_variant')
    )

    # Integer columns that bioframe's overlap() converts to float (due to NaN)
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

    df = pl.from_pandas(df).cast(int_columns) # type: ignore

    # determine the number of overlapping SNVs for each CpG record
    CG_cols = df_meth_dinucleotides.columns
    df = df.with_columns(
        pl
        .when(pl.col("start_variant").is_null())
        .then(pl.lit(0))
        .otherwise(pl.col("start_variant").count().over(CG_cols))
        .alias("num_SNVs_overlapping_CG")
    )

    # Rename allele columns: strip '_variant' suffix added by bioframe overlap
    allele_renames = {
        col + '_variant': col
        for col in df_joint_called_variants.columns
        if col.startswith('allele_')
    }

    df = (
        df
        .rename({
            'start': 'start_cpg',
            'end': 'end_cpg',
            'is_within_50bp_of_mismatch_site_pat': 'cpg_is_within_50bp_of_mismatch_site_pat',
            'is_within_50bp_of_mismatch_site_mat': 'cpg_is_within_50bp_of_mismatch_site_mat',
            **allele_renames,
        })
        .drop(['chrom_variant'])
    )

    return df


def label_cpgs_as_allele_specific(df):
    """Label each unique CpG record with per-member allele-specific flags.

    For each family member (kid, dad, mom), adds columns:
        snv_genotypes_{role}: comma-separated genotype string ("het", "hom", or ".")
        cpg_is_allele_specific_{role}: True if any overlapping SNV is het for that member
    Also adds cpg_overlaps_at_least_one_snv (shared across members).
    """
    roles = ['kid', 'dad', 'mom']

    # Columns that define a unique CpG record (everything except per-variant fields)
    variant_cols = {'start_variant', 'end_variant', 'num_SNVs_overlapping_CG'}
    for role in roles:
        variant_cols.add(f'allele_1_{role}')
        variant_cols.add(f'allele_2_{role}')
    CG_cols = [col for col in df.columns if col not in variant_cols]

    # Compute per-member genotype at each overlapping SNV
    genotype_exprs = []
    for role in roles:
        a1 = f'allele_1_{role}'
        a2 = f'allele_2_{role}'
        genotype_exprs.append(
            pl
            .when(
                pl.col(a1).is_null() &
                pl.col(a2).is_null()
            )
            .then(pl.lit("."))
            .when(
                (pl.col(a1) == '.') |
                (pl.col(a2) == '.')
            )
            .then(pl.lit("."))
            .when(pl.col(a1) == pl.col(a2))
            .then(pl.lit("hom"))
            .otherwise(pl.lit("het"))
            .alias(f"genotype_{role}")
        )

    df = df.with_columns(genotype_exprs)

    # Aggregate: group by CpG record, join genotypes, flag allele-specific
    agg_exprs = []
    for role in roles:
        agg_exprs.extend([
            pl.col(f"genotype_{role}").str.join(",").alias(f"snv_genotypes_{role}"),
            (pl.col(f"genotype_{role}") == "het").any().alias(f"cpg_is_allele_specific_{role}"),
        ])

    df = (
        df
        .group_by(CG_cols)
        .agg(agg_exprs)
        .with_columns(
            pl
            .when(pl.col('snv_genotypes_kid') == '.')
            .then(False)
            .otherwise(True)
            .alias("cpg_overlaps_at_least_one_snv")
        )
        .sort(['chrom', 'start_cpg'])
    )

    # Reorder: move the new columns to the end
    new_cols = ['cpg_overlaps_at_least_one_snv']
    for role in roles:
        new_cols.extend([f'snv_genotypes_{role}', f'cpg_is_allele_specific_{role}'])
    base_cols = [col for col in df.columns if col not in new_cols]
    return df.select(base_cols + new_cols)


@profile # type:ignore
def main():
    parser = argparse.ArgumentParser(
        description='Expand output of tapestry trio workflow to include all CpG sites and unphased DNA methylation levels'
    )
    parser.add_argument('--bed_all_cpgs_in_reference', required=True, help='All CpG sites observed in reference genome')
    parser.add_argument('--bed_meth_count_unphased_kid', required=True, help='Unphased count-based methylation levels for kid')
    parser.add_argument('--bed_meth_model_unphased_kid', required=True, help='Unphased model-based methylation levels for kid')
    parser.add_argument('--bed_meth_count_unphased_dad', required=True, help='Unphased count-based methylation levels for dad')
    parser.add_argument('--bed_meth_model_unphased_dad', required=True, help='Unphased model-based methylation levels for dad')
    parser.add_argument('--bed_meth_count_unphased_mom', required=True, help='Unphased count-based methylation levels for mom')
    parser.add_argument('--bed_meth_model_unphased_mom', required=True, help='Unphased model-based methylation levels for mom')
    parser.add_argument('--bed_meth_parent_phased', required=True, help='Parent-phased methylation levels')
    parser.add_argument('--bed_het_site_mismatches_pat', required=True, help='Heterozygous sites where bit vectors are mismatched (paternal)')
    parser.add_argument('--bed_het_site_mismatches_mat', required=True, help='Heterozygous sites where bit vectors are mismatched (maternal)')
    parser.add_argument('--bed_meth_parent_phased_all_cpgs', required=True, help='Parent-phased methylation levels at all CpG sites, both in reference and sample, including null methylation levels, and unphased methylation levels')
    parser.add_argument('--kid_id', required=True, help='Child sample ID in joint-called multi-sample vcf')
    parser.add_argument('--dad_id', required=True, help='Father sample ID in joint-called multi-sample vcf')
    parser.add_argument('--mom_id', required=True, help='Mother sample ID in joint-called multi-sample vcf')
    parser.add_argument('--vcf_joint_called', required=True, help='Joint-called multi-sample vcf')
    parser.add_argument('--output_dir', required=True, help='Output directory for bigwig files representing unphased methylation levels')
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(filename)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    logger = logging.getLogger(__name__)

    logger.info(f"Starting '{__file__}'")
    logger.info("Script started with the following arguments: %s", vars(args))
    logger.info(f"Polars version: {pl.__version__}")
    logger.info(f"Polars thread pool size: {pl.thread_pool_size()}")

    Path(args.output_dir).mkdir(parents=True, exist_ok=True)

    logger.info("Writing combined (unphased) methylation bigwig files...")
    combined_specs = [
        (args.kid_id, args.bed_meth_count_unphased_kid, args.bed_meth_model_unphased_kid),
        (args.dad_id, args.bed_meth_count_unphased_dad, args.bed_meth_model_unphased_dad),
        (args.mom_id, args.bed_meth_count_unphased_mom, args.bed_meth_model_unphased_mom),
    ]
    for uid, bed_count, bed_model in combined_specs:
        write_combined_bigwig(bed_count, 'count', uid, args.output_dir, logger)
        write_combined_bigwig(bed_model, 'model', uid, args.output_dir, logger)
    logger.info("Done writing combined (unphased) methylation bigwig files")

    df_all_cpgs_in_reference = read_all_cpgs_in_reference(args.bed_all_cpgs_in_reference)
    logger.info(f"Read all CpG sites in reference genome")
    report_size(df_all_cpgs_in_reference, 'All CpGs in reference', logger)

    df_meth_unphased = read_meth_unphased_trio(
        bed_meth_count_kid=args.bed_meth_count_unphased_kid,
        bed_meth_model_kid=args.bed_meth_model_unphased_kid,
        bed_meth_count_dad=args.bed_meth_count_unphased_dad,
        bed_meth_model_dad=args.bed_meth_model_unphased_dad,
        bed_meth_count_mom=args.bed_meth_count_unphased_mom,
        bed_meth_model_mom=args.bed_meth_model_unphased_mom,
    )
    logger.info(f"Read CpG sites at which unphased methylation levels are available for kid, dad, and mom")
    report_size(df_meth_unphased, 'Unphased methylation levels', logger)

    if Path(args.bed_meth_parent_phased).exists():
        df_meth_parent_phased = read_meth_parent_phased(args.bed_meth_parent_phased)
        logger.info(f"Read CpG sites at which count- and model-based methylation levels have been phased to parents")
        report_size(df_meth_parent_phased, 'Phased methylation levels', logger)
    else:
        logger.warning(f"Could not read CpG sites at which count- and model-based methylation levels have been phased to parental haplotypes")
        logger.warning(f"Required file does not exist: '{args.bed_meth_parent_phased}'")
        logger.info(f"Done running '{__file__}'")
        return

    # In the following, to save memory, we update what "df" points to
    # (instead of creating new references for each updated dataframe of CpG sites)

    df = expand_meth_to_all_cpgs(df_all_cpgs_in_reference, df_meth_unphased, df_meth_parent_phased)
    logger.info(f"Expanded DNA methylation dataset to include (1) all CpG sites observed in reference and sample genomes, and (2) unphased methylation levels (where available)")
    report_size(df, 'CpG Methylation', logger)

    # Release memory:
    del df_all_cpgs_in_reference
    del df_meth_unphased
    del df_meth_parent_phased
    import gc
    gc.collect() # type:ignore
    logger.info(f"Removed df_all_cpgs_in_reference, df_meth_unphased, df_meth_parent_phased, which are no longer needed")

    df = compute_proximity_to_mismatched_heterozygous_sites(
        df,
        args.bed_het_site_mismatches_pat,
        args.bed_het_site_mismatches_mat,
    )
    logger.info(f"Computed proximity of all CpG sites to paternal and maternal heterozygous mismatch sites")
    report_size(df, 'CpG Methylation', logger)
    compute_fraction_of_cpgs_that_are_close_to_mismatches(df, logger)
    compute_methylation_coverage_qc(df, logger)

    df_joint_called_variants = get_joint_called_variants(args.kid_id, args.dad_id, args.mom_id, args.vcf_joint_called)
    logger.info(f"Got SNVs from joint-called vcf for kid, dad, and mom")
    report_size(df_joint_called_variants, 'Joint-called SNVs', logger)

    df = label_with_variants(df, df_joint_called_variants)
    logger.info(f"Determined which CpG sites overlap 1 or 2 SNVs")
    report_size(df, 'CpG Methylation', logger)

    # Release memory:
    del df_joint_called_variants
    gc.collect() # type:ignore
    logger.info(f"Removed df_joint_called_variants, which is no longer needed")

    df = label_cpgs_as_allele_specific(df)
    logger.info(f"Flagged CpG sites that are allele-specific by assessing overlap with het SNVs, for each family member")
    report_size(df, 'CpG Methylation', logger)

    methylation_data_root = write_methylation(df, args.bed_meth_parent_phased_all_cpgs, source=f"{__file__} with args {vars(args)}")
    logger.info(f"Wrote expanded and allele-specific-flagged methylation dataframe to: '{methylation_data_root}.sorted.bed.gz'")
    logger.info(f"Index exists at: '{methylation_data_root}.sorted.bed.gz.tbi'")

    logger.info(f"Done running '{__file__}'")


if __name__ == "__main__":
    main()
