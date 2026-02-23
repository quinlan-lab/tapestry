import argparse
import logging
from pathlib import Path
import bioframe as bf
import polars as pl

from get_meth_hap1_hap2 import get_meth_hap1_hap2
from phase_meth_to_founder_haps import (
    write_bigwig,
    combine_count_and_model_based_methylation_levels,
)
from util.write_data import write_bed
from util.version_sort import version_sort
from util.shell import shell

REFERENCE_GENOME = "hg38"


def read_ped(ped_path: str, child_uid: str):
    """
    Parse a PED file (family sample father mother sex) and return (father_id, mother_id)
    for the given child sample.
    """
    with open(ped_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if parts[1] == child_uid:
                return parts[2], parts[3]
    raise ValueError(f"Sample '{child_uid}' not found in PED file '{ped_path}'")


def get_whatshap_blocks(uid: str, vcf_trio_phased: str, output_dir: str, logger) -> pl.DataFrame:
    """
    Run `whatshap stats --block-list` to get phase block intervals and het SNV counts,
    then read the resulting TSV into a Polars DataFrame.

    whatshap stats --block-list produces a TSV with columns:
        #sample  chromosome  phase_set  from  to  variants
    """
    blocks_tsv = f"{output_dir}/{uid}.whatshap.blocks.tsv"
    cmd = (
        f"whatshap stats"
        f" --sample {uid}"
        f" --block-list {blocks_tsv}"
        f" {vcf_trio_phased}"
    )
    logger.info(f"Running: {cmd}")
    shell(cmd)
    logger.info(f"Wrote whatshap blocks to '{blocks_tsv}'")

    df = (
        pl.read_csv(
            blocks_tsv,
            separator='\t',
            has_header=True,
        )
        .rename({
            '#sample': 'sample',
            'chromosome': 'chrom',
            'from': 'start',
            'to': 'end',
            'variants': 'num_het_SNVs_in_hap_map_block',
        })
        .filter(pl.col('sample') == uid)
        .select(['chrom', 'start', 'end', 'num_het_SNVs_in_hap_map_block'])
        .with_columns([
            (pl.col('start') - 1).alias('start'),  # convert 1-based to 0-based
        ])
    )
    return df


def phase_meth_to_parent_haps(df_meth_hap1_hap2: pl.DataFrame, df_blocks: pl.DataFrame, dad_id: str, mom_id: str) -> pl.DataFrame:
    """
    Left-join methylation sites onto whatshap phase blocks.
    Since WhatsHap --ped phases hap1=paternal and hap2=maternal, the mapping is direct.
    """
    df = bf.overlap(
        df_meth_hap1_hap2.to_pandas(),
        df_blocks.to_pandas(),
        how='left',
        suffixes=('', '_'),
    )

    df = (
        pl
        .from_pandas(df)
        .drop(['chrom_'])
        .rename({
            'start_': 'start_hap_map_block',
            'end_': 'end_hap_map_block',
            'num_het_SNVs_in_hap_map_block_': 'num_het_SNVs_in_hap_map_block',
        })
        .cast({'num_het_SNVs_in_hap_map_block': pl.Int64})
        .with_columns([
            pl.col('methylation_level_hap1').alias('methylation_level_pat'),
            pl.col('methylation_level_hap2').alias('methylation_level_mat'),
            pl.col('total_read_count_hap1').alias('total_read_count_pat'),
            pl.col('total_read_count_hap2').alias('total_read_count_mat'),
            pl.lit(dad_id).alias('founder_haplotype_pat'),
            pl.lit(mom_id).alias('founder_haplotype_mat'),
            pl.lit(None).cast(pl.Float64).alias('haplotype_concordance_in_hap_map_block'),
        ])
        .drop([
            'methylation_level_hap1', 'methylation_level_hap2',
            'total_read_count_hap1', 'total_read_count_hap2',
        ])
    )

    return df


def main():
    parser = argparse.ArgumentParser(description='Phase HiFi-based DNA methylation data to parental haplotypes for a trio using WhatsHap pedigree phasing')
    parser.add_argument('--uid', required=True, help='Child sample ID')
    parser.add_argument('--vcf_trio_phased', required=True, help='WhatsHap trio-phased VCF (bgzipped + tabixed)')
    parser.add_argument('--ped', required=True, help='Pedigree file (PED format)')
    parser.add_argument('--bed_meth_count_hap1', required=True, help='Bed file of count-based methylation levels for hap1 (paternal)')
    parser.add_argument('--bed_meth_count_hap2', required=True, help='Bed file of count-based methylation levels for hap2 (maternal)')
    parser.add_argument('--bed_meth_model_hap1', required=True, help='Bed file of model-based methylation levels for hap1 (paternal)')
    parser.add_argument('--bed_meth_model_hap2', required=True, help='Bed file of model-based methylation levels for hap2 (maternal)')
    parser.add_argument('--output_dir', required=True, help='Output directory')
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(filename)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    logger = logging.getLogger(__name__)

    logger.info(f"Starting '{__file__}'")
    logger.info("Script started with the following arguments: %s", vars(args))

    Path(args.output_dir).mkdir(parents=True, exist_ok=True)

    dad_id, mom_id = read_ped(args.ped, args.uid)
    logger.info(f"Read PED file: father='{dad_id}', mother='{mom_id}'")

    df_blocks = get_whatshap_blocks(args.uid, args.vcf_trio_phased, args.output_dir, logger)
    logger.info(f"Got whatshap phase blocks: {len(df_blocks)} rows, {len(df_blocks.columns)} columns")

    df_meth_count_hap1_hap2 = get_meth_hap1_hap2(
        pb_cpg_tool_mode='count',
        bed_hap1=args.bed_meth_count_hap1,
        bed_hap2=args.bed_meth_count_hap2,
    )
    logger.info(f"Got count-based methylation levels: {len(df_meth_count_hap1_hap2)} rows")

    df_meth_model_hap1_hap2 = get_meth_hap1_hap2(
        pb_cpg_tool_mode='model',
        bed_hap1=args.bed_meth_model_hap1,
        bed_hap2=args.bed_meth_model_hap2,
    )
    logger.info(f"Got model-based methylation levels: {len(df_meth_model_hap1_hap2)} rows")

    df_meth_count_parent_phased = phase_meth_to_parent_haps(df_meth_count_hap1_hap2, df_blocks, dad_id, mom_id)
    logger.info(f"Phased count-based methylation levels to parental haplotypes: {len(df_meth_count_parent_phased)} rows")

    df_meth_model_parent_phased = phase_meth_to_parent_haps(df_meth_model_hap1_hap2, df_blocks, dad_id, mom_id)
    logger.info(f"Phased model-based methylation levels to parental haplotypes: {len(df_meth_model_parent_phased)} rows")

    df_meth_parent_phased = combine_count_and_model_based_methylation_levels(df_meth_count_parent_phased, df_meth_model_parent_phased)
    logger.info(f"Combined count- and model-based methylation levels: {len(df_meth_parent_phased)} rows")

    write_bed(args.output_dir, df_meth_parent_phased, filename_stem=f"{args.uid}.dna-methylation.founder-phased")
    logger.info(f"Wrote parent-phased methylation levels to: '{args.output_dir}'")

    for parental in ['pat', 'mat']:
        for pb_cpg_tool_mode in ['count', 'model']:
            write_bigwig(df_meth_parent_phased, args.uid, parental, pb_cpg_tool_mode, args.output_dir, logger)

    logger.info(f"Done running '{__file__}'")


if __name__ == "__main__":
    main()
