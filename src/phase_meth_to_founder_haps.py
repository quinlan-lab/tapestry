import argparse
import json
import logging
from pathlib import Path
from typing import cast
import bioframe as bf # https://bioframe.readthedocs.io/en/latest/index.html
import polars as pl
import pyBigWig

from phasing_pedigree import (
    get_read_phasing, 
    get_read_phase_blocks, 
    get_iht_phasing, 
    get_iht_blocks, 
    get_all_phasing
)
from hap_map_pedigree import get_hap_map
from util.hap_map import write_hap_map_blocks
from get_meth_hap1_hap2 import read_meth_hap1_hap2
from util.write_data import write_bed, write_bit_vector_mismatches_bed, write_bit_vector_mismatches_vcf
from util.version_sort import version_sort

REFERENCE_GENOME = "hg38"

def get_all_phasing_wrapper(uid, vcf_read_phased, tsv_read_phase_blocks, vcf_iht_phased, txt_iht_blocks, logger):
    df_read_phasing = get_read_phasing(vcf_read_phased)
    logger.info(f"Got read-based phasing data: {len(df_read_phasing)} rows, {len(df_read_phasing.columns)} columns")
    
    df_read_phase_blocks = get_read_phase_blocks(tsv_read_phase_blocks)
    logger.info(f"Got read-based phase blocks: {len(df_read_phase_blocks)} rows, {len(df_read_phase_blocks.columns)} columns")

    df_iht_phasing = get_iht_phasing(uid, vcf_iht_phased)
    logger.info(f"Got inheritance-based phasing data: {len(df_iht_phasing)} rows, {len(df_iht_phasing.columns)} columns")

    df_iht_blocks = get_iht_blocks(uid, txt_iht_blocks)
    logger.info(f"Got inheritance-based phase blocks: {len(df_iht_blocks)} rows, {len(df_iht_blocks.columns)} columns")

    df_all_phasing = get_all_phasing(
        df_read_phasing, 
        df_read_phase_blocks, 
        df_iht_phasing, 
        df_iht_blocks
    )
    
    return df_all_phasing

def phase_meth_to_founder_haps(df_meth_hap1_hap2, df_hap_map):
    if df_hap_map.is_empty():
        return (
            df_meth_hap1_hap2
            .with_columns(
                pl.lit(None, dtype=pl.Int64).alias("start_hap_map_block"),
                pl.lit(None, dtype=pl.Int64).alias("end_hap_map_block"),
                pl.lit(None, dtype=pl.Float64).alias(
                    "haplotype_concordance_in_hap_map_block"
                ),
                pl.lit(None, dtype=pl.Int64).alias(
                    "num_het_SNVs_in_hap_map_block"
                ),
                pl.lit(None, dtype=pl.Float64).alias("methylation_level_pat"),
                pl.lit(None, dtype=pl.Float64).alias("methylation_level_mat"),
                pl.lit(None, dtype=pl.Int64).alias("total_read_count_pat"),
                pl.lit(None, dtype=pl.Int64).alias("total_read_count_mat"),
                pl.lit(None, dtype=pl.String).alias("founder_haplotype_pat"),
                pl.lit(None, dtype=pl.String).alias("founder_haplotype_mat"),
            )
            .drop(
                [
                    "total_read_count_hap1",
                    "total_read_count_hap2",
                    "methylation_level_hap1",
                    "methylation_level_hap2",
                ]
            )
        )

    df = bf.overlap(
        df_meth_hap1_hap2.to_pandas(), 
        df_hap_map.to_pandas(), 
        how='left', # we want to retain all CpG sites, including those that we cannot phase to founder haplotypes
        suffixes=('','_'),
    )

    df = (
        pl
        .from_pandas(df)
        .drop([
            "chrom_"
        ])
        .rename({ 
            "start_": "start_hap_map_block",
            "end_": "end_hap_map_block",
            "paternal_haplotype_": "paternal_haplotype_in_hap_map_block", 
            "maternal_haplotype_": "maternal_haplotype_in_hap_map_block",
            "haplotype_concordance_": "haplotype_concordance_in_hap_map_block",
            "num_het_SNVs_": "num_het_SNVs_in_hap_map_block",
        })
        .cast({
            "num_het_SNVs_in_hap_map_block": pl.Int64,
            "total_read_count_hap1": pl.Int64,
            "total_read_count_hap2": pl.Int64,
        })
        .with_columns(
            methylation_level_pat=(
                pl
                .when(pl.col("paternal_haplotype_in_hap_map_block").str.ends_with("_hap1"))
                .then(pl.col("methylation_level_hap1"))
                .when(pl.col("paternal_haplotype_in_hap_map_block").str.ends_with("_hap2"))
                .then(pl.col("methylation_level_hap2"))
                .otherwise(None)
            ),
            methylation_level_mat=(
                pl
                .when(pl.col("maternal_haplotype_in_hap_map_block").str.ends_with("_hap1"))
                .then(pl.col("methylation_level_hap1"))
                .when(pl.col("maternal_haplotype_in_hap_map_block").str.ends_with("_hap2"))
                .then(pl.col("methylation_level_hap2"))
                .otherwise(None)
            ),
            total_read_count_pat=(
                pl
                .when(pl.col("paternal_haplotype_in_hap_map_block").str.ends_with("_hap1"))
                .then(pl.col("total_read_count_hap1"))
                .when(pl.col("paternal_haplotype_in_hap_map_block").str.ends_with("_hap2"))
                .then(pl.col("total_read_count_hap2"))
                .otherwise(None)
            ),
            total_read_count_mat=(
                pl
                .when(pl.col("maternal_haplotype_in_hap_map_block").str.ends_with("_hap1"))
                .then(pl.col("total_read_count_hap1"))
                .when(pl.col("maternal_haplotype_in_hap_map_block").str.ends_with("_hap2"))
                .then(pl.col("total_read_count_hap2"))
                .otherwise(None)
            ),
            founder_haplotype_pat=pl.col("paternal_haplotype_in_hap_map_block").str.split("_").list.get(0),
            founder_haplotype_mat=pl.col("maternal_haplotype_in_hap_map_block").str.split("_").list.get(0),
        )
        .drop([
            "total_read_count_hap1", "total_read_count_hap2",
            "methylation_level_hap1", "methylation_level_hap2",
            "paternal_haplotype_in_hap_map_block", "maternal_haplotype_in_hap_map_block",
        ])
    )

    return df 

def _read_fai_chromsizes(reference_fai):
    chromsizes = []
    with open(reference_fai) as handle:
        for line_number, line in enumerate(handle, 1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                raise ValueError(
                    f"Malformed FASTA index line {line_number} in {reference_fai}"
                )
            chromsizes.append((fields[0], int(fields[1])))
    return chromsizes


def write_bigwig(
    df,
    uid,
    parental,
    pb_cpg_tool_mode,
    output_dir,
    logger=None,
    reference_fai=None,
    reference_name=REFERENCE_GENOME,
):
    """
    Write a bigwig file for a given parental haplotype and given pb-cpg-tools pileup mode
    """
    df_bed_graph = (
        df 
        .filter(pl.col(f"methylation_level_{parental}_{pb_cpg_tool_mode}").is_not_null())        
        .select([
            pl.col("chrom"),
            pl.col("start"),
            pl.col("end"),
            pl.col(f"methylation_level_{parental}_{pb_cpg_tool_mode}")
        ])
        .sort(["chrom", "start", "end"])
    )

    file_path = f"{output_dir}/{uid}.dna-methylation.{parental}.{pb_cpg_tool_mode}.{reference_name}.bw"
    if reference_fai is not None:
        chromsizes = _read_fai_chromsizes(reference_fai)
    else:
        # Legacy CLI compatibility. Generic workflow calls must provide an FAI.
        legacy_sizes = bf.fetch_chromsizes(db=REFERENCE_GENOME)
        chromsizes = list(zip(legacy_sizes.index, legacy_sizes.values))

    size_by_chrom = dict(chromsizes)
    unknown = sorted(set(df_bed_graph["chrom"].to_list()) - set(size_by_chrom))
    if unknown:
        raise ValueError(f"BigWig records use contigs absent from the FASTA index: {unknown}")

    bigwig = pyBigWig.open(file_path, "w")
    try:
        bigwig.addHeader(chromsizes)
        if not df_bed_graph.is_empty():
            bigwig.addEntries(
                df_bed_graph["chrom"].to_list(),
                df_bed_graph["start"].to_list(),
                ends=df_bed_graph["end"].to_list(),
                values=df_bed_graph[
                    f"methylation_level_{parental}_{pb_cpg_tool_mode}"
                ].to_list(),
            )
    finally:
        bigwig.close()

    if logger: 
        logger.info(
            "Wrote bigwig file for %s %s-based methylation levels against %s to: '%s'",
            parental,
            pb_cpg_tool_mode,
            reference_name,
            file_path,
        )

def add_suffix(
    df: pl.DataFrame, suffix: str, cols_to_suffix: list[str]
) -> pl.DataFrame:
    df.columns = [f"{col}{suffix}" if col in cols_to_suffix else col for col in df.columns]
    return df 

def combine_count_and_model_based_methylation_levels(
        df_meth_count: pl.DataFrame | None,
        df_meth_model: pl.DataFrame
    ):
    suffix_count, suffix_model = '_count', '_model'
    unique_cols = ['methylation_level_pat', 'methylation_level_mat']

    df_meth_model = add_suffix(df_meth_model, suffix=suffix_model, cols_to_suffix=unique_cols)

    unique_cols_count = [f'{unique_col}{suffix_count}' for unique_col in unique_cols]
    unique_cols_model = [f'{unique_col}{suffix_model}' for unique_col in unique_cols]

    if df_meth_count is None:
        common_cols = [col for col in df_meth_model.columns if col not in unique_cols_model]
        return version_sort(
            df_meth_model
            .with_columns(
                pl.lit(None, dtype=pl.Float64).alias(unique_cols_count[0]),
                pl.lit(None, dtype=pl.Float64).alias(unique_cols_count[1]),
            )
            .select(common_cols + unique_cols_count + unique_cols_model)
        )

    count_frame = add_suffix(
        cast(pl.DataFrame, df_meth_count),
        suffix=suffix_count,
        cols_to_suffix=unique_cols,
    )
    common_cols = [col for col in count_frame.columns if col not in unique_cols_count]

    df = count_frame.join(
        df_meth_model,
        on=common_cols,
        nulls_equal=True, # join on nulls, e.g., total_read_count_XXX is null in both df_meth_count and df_meth_model
        # "how=full" will capture CpG sites where count-based meth levels are available, but not model-based meth levels, and vice versa,
        # though, for sample 200081, I didn't see any such CpG sites: 
        how="full", 
        coalesce=True, 
    )   

    df = df.select(common_cols + unique_cols_count + unique_cols_model)

    return version_sort(df)

def main():
    parser = argparse.ArgumentParser(description='Phase HiFi-based DNA methylation data to founder haplotypes')
    parser.add_argument('--uid', required=True, help='Sample UID in joint-called multi-sample vcf')
    parser.add_argument('--vcf_read_phased', required=True, help='Single-sample vcf from hiphase')
    parser.add_argument('--tsv_read_phase_blocks', required=True, help='Single-sample tsv from hiphase')
    parser.add_argument('--vcf_iht_phased', required=True, help='Joint-called multi-sample vcf from gtg-ped-map/gtg-concordance')
    parser.add_argument('--txt_iht_blocks', required=True, help='Multi-sample iht blocks file from gtg-ped-map/gtg-concordance')
    parser.add_argument('--bed_meth_count_hap1', help='Optional BED of count-based methylation levels for hap1')
    parser.add_argument('--bed_meth_count_hap2', help='Optional BED of count-based methylation levels for hap2')
    parser.add_argument('--bed_meth_model_hap1', required=True, help='Bed file of model-based methylation levels from aligned_bam_to_cpg_scores for hap1')
    parser.add_argument('--bed_meth_model_hap2', required=True, help='Bed file of model-based methylation levels from aligned_bam_to_cpg_scores for hap2')
    parser.add_argument('--reference_fai', help='FASTA index used for BigWig chromosome lengths')
    parser.add_argument('--reference_name', default=REFERENCE_GENOME, help='Reference label used in BigWig filenames')
    parser.add_argument('--no-bigwig', action='store_true', help='Do not write BigWig outputs')
    parser.add_argument('--output_dir', required=True, help='Output directory')
    args = parser.parse_args()

    if bool(args.bed_meth_count_hap1) != bool(args.bed_meth_count_hap2):
        parser.error("--bed_meth_count_hap1 and --bed_meth_count_hap2 must be provided together")

    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(filename)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    logger = logging.getLogger(__name__)

    logger.info(f"Starting '{__file__}'")
    logger.info("Script started with the following arguments: %s", vars(args))

    Path(args.output_dir).mkdir(parents=True, exist_ok=True)
    
    df_all_phasing = get_all_phasing_wrapper(args.uid, args.vcf_read_phased, args.tsv_read_phase_blocks, args.vcf_iht_phased, args.txt_iht_blocks, logger)
    logger.info(f"Got all phasing data: {len(df_all_phasing)} rows, {len(df_all_phasing.columns)} columns")
    
    df_hap_map, df_sites, df_sites_mismatch = get_hap_map(df_all_phasing)
    logger.info(f"Got hap map: {len(df_hap_map)} rows, {len(df_hap_map.columns)} columns")
    logger.info(f"Got heterozygous sites: {len(df_sites)} rows, {len(df_sites.columns)} columns")
    logger.info(f"Got heterozygous sites where read-based and inheritance-based bit vectors don't match: {len(df_sites_mismatch)} rows, {len(df_sites_mismatch.columns)} columns")
    
    write_hap_map_blocks(df_hap_map, args.uid, "paternal", args.output_dir)
    write_hap_map_blocks(df_hap_map, args.uid, "maternal", args.output_dir)
    logger.info(f"Wrote paternal and maternal hap-map blocks for IGV visualization to '{args.output_dir}'")
    
    write_bed(args.output_dir, df_hap_map, f"{args.uid}.hap-map-blocks")
    logger.info(f"Wrote hap-map blocks to '{args.output_dir}'")
    
    write_bit_vector_mismatches_vcf(args.output_dir, df_sites_mismatch, logger, uid=args.uid)
    write_bit_vector_mismatches_bed(args.output_dir, df_sites_mismatch, logger, uid=args.uid)    

    df_meth_count_hap1_hap2 = None
    if args.bed_meth_count_hap1:
        df_meth_count_hap1_hap2 = read_meth_hap1_hap2(
            pb_cpg_tool_mode='count',
            bed_hap1=args.bed_meth_count_hap1,
            bed_hap2=args.bed_meth_count_hap2
        )
        logger.info(f"Got read-based phasing of count-based methylation levels: {len(df_meth_count_hap1_hap2)} rows, {len(df_meth_count_hap1_hap2.columns)} columns")
    df_meth_model_hap1_hap2 = read_meth_hap1_hap2(
        pb_cpg_tool_mode='model', 
        bed_hap1=args.bed_meth_model_hap1, 
        bed_hap2=args.bed_meth_model_hap2
    )    
    logger.info(f"Got read-based phasing of model-based methylation levels: {len(df_meth_model_hap1_hap2)} rows, {len(df_meth_model_hap1_hap2.columns)} columns")
    
    df_meth_count_founder_phased = None
    if df_meth_count_hap1_hap2 is not None:
        df_meth_count_founder_phased = phase_meth_to_founder_haps(df_meth_count_hap1_hap2, df_hap_map)
        logger.info(f"Phased count-based methylation levels to founder haplotypes: {len(df_meth_count_founder_phased)} rows, {len(df_meth_count_founder_phased.columns)} columns")
    df_meth_model_founder_phased = phase_meth_to_founder_haps(df_meth_model_hap1_hap2, df_hap_map)
    logger.info(f"Phased model-based methylation levels to founder haplotypes: {len(df_meth_model_founder_phased)} rows, {len(df_meth_model_founder_phased.columns)} columns")
    df_meth_founder_phased = combine_count_and_model_based_methylation_levels(df_meth_count_founder_phased, df_meth_model_founder_phased)
    enabled_modes = ['model'] if df_meth_count_founder_phased is None else ['count', 'model']
    logger.info(f"Combined enabled methylation levels {enabled_modes}: {len(df_meth_founder_phased)} rows, {len(df_meth_founder_phased.columns)} columns")

    write_bed(args.output_dir, df_meth_founder_phased, filename_stem=f"{args.uid}.dna-methylation")
    logger.info(f"Wrote enabled methylation levels, phased to founder haplotypes, to: '{args.output_dir}'")

    qc_status = "no_inheritance_phase" if df_hap_map.is_empty() else "complete"
    qc = {
        "sample_id": args.uid,
        "status": qc_status,
        "enabled_modes": enabled_modes,
        "hap_map_blocks": len(df_hap_map),
        "informative_sites": len(df_sites),
        "mismatch_sites": len(df_sites_mismatch),
    }
    (Path(args.output_dir) / f"{args.uid}.phasing-qc.json").write_text(
        json.dumps(qc, indent=2, sort_keys=True) + "\n"
    )
    
    if not args.no_bigwig:
        for parental in ['pat', 'mat']:
            for pb_cpg_tool_mode in enabled_modes:
                write_bigwig(
                    df_meth_founder_phased,
                    args.uid,
                    parental,
                    pb_cpg_tool_mode,
                    args.output_dir,
                    logger,
                    reference_fai=args.reference_fai,
                    reference_name=args.reference_name,
                )

    logger.info(f"Done running '{__file__}'")

if __name__ == "__main__":
    main()
