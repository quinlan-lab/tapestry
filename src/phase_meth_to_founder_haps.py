import argparse
import logging
from pathlib import Path
import bioframe as bf # https://bioframe.readthedocs.io/en/latest/index.html
import polars as pl

from get_all_phasing import (
    get_read_phasing, 
    get_read_phase_blocks, 
    get_iht_phasing, 
    get_iht_blocks, 
    get_all_phasing
)
from get_hap_map import (
    get_hap_map,
    write_hap_map_blocks,
    write_bit_vector_sites_and_mismatches,
)
from get_meth_hap1_hap2 import get_meth_hap1_hap2
from util.write_data import write_bed, write_bed_and_header
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

def write_bit_vector_mismatches(output_dir, uid, df_sites_mismatch, logger):
    bed = f"{output_dir}/{uid}.bit-vector-sites-mismatches.bed"
    write_bed_and_header(bed, df_sites_mismatch)
    logger.info(f"Wrote sites where bit vectors are mismatched, for later computation of proximity of (all) cpg sites to mismatch sites to: '{bed}'")

def phase_meth_to_founder_haps(df_meth_hap1_hap2, df_hap_map):
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

def write_bigwig(df, uid, parental, pb_cpg_tool_mode, output_dir, logger=None):
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
        .to_pandas()  # Convert to pandas DataFrame for bioframe compatibility
    )

    file_path = f"{output_dir}/{uid}.dna-methylation.founder-phased.{parental}.{pb_cpg_tool_mode}.{REFERENCE_GENOME}.bw"

    # https://bioframe.readthedocs.io/en/latest/api-fileops.html#bioframe.io.fileops.to_bigwig
    bf.to_bigwig(
        df_bed_graph, 
        # https://bioframe.readthedocs.io/en/latest/api-resources.html#bioframe.io.resources.fetch_chromsizes 
        bf.fetch_chromsizes(db=REFERENCE_GENOME), 
        outpath=file_path, 
        # we assume that user has bedGraphToBigWig in their PATH: 
        path_to_binary="bedGraphToBigWig" # type: ignore
    )

    if logger: 
        logger.info(f"Wrote bigwig file for {parental} {pb_cpg_tool_mode}-based methylation levels, ASSUMING {REFERENCE_GENOME}, to: '{file_path}'")

def add_suffix(df, suffix, cols_to_suffix): 
    df.columns = [f"{col}{suffix}" if col in cols_to_suffix else col for col in df.columns]
    return df 

def combine_count_and_model_based_methylation_levels(
        df_meth_count: pl.DataFrame, 
        df_meth_model: pl.DataFrame
    ):
    suffix_count, suffix_model = '_count', '_model'
    unique_cols = ['methylation_level_pat', 'methylation_level_mat']

    df_meth_count = add_suffix(df_meth_count, suffix=suffix_count, cols_to_suffix=unique_cols)
    df_meth_model = add_suffix(df_meth_model, suffix=suffix_model, cols_to_suffix=unique_cols)

    unique_cols_count = [f'{unique_col}{suffix_count}' for unique_col in unique_cols]
    common_cols = [col for col in df_meth_count.columns if col not in unique_cols_count]

    df = df_meth_count.join(
        df_meth_model,
        on=common_cols,
        join_nulls=True,
        # "how=full" will capture CpG sites where count-based meth levels are available, but not model-based meth levels, and vice versa,
        # though, for sample 200081, I didn't see any such CpG sites: 
        how="full", 
        coalesce=True, 
    )   

    unique_cols_model = [f'{unique_col}{suffix_model}' for unique_col in unique_cols]
    df = df.select(common_cols + unique_cols_count + unique_cols_model)

    return version_sort(df)

def main():
    parser = argparse.ArgumentParser(description='Phase HiFi-based DNA methylation data to founder haplotypes')
    parser.add_argument('--uid', required=True, help='Sample UID in joint-called multi-sample vcf')
    parser.add_argument('--vcf_read_phased', required=True, help='Single-sample vcf from hiphase')
    parser.add_argument('--tsv_read_phase_blocks', required=True, help='Single-sample tsv from hiphase')
    parser.add_argument('--vcf_iht_phased', required=True, help='Joint-called multi-sample vcf from gtg-ped-map/gtg-concordance')
    parser.add_argument('--txt_iht_blocks', required=True, help='Multi-sample iht blocks file from gtg-ped-map/gtg-concordance')
    parser.add_argument('--bed_meth_count_hap1', required=True, help='Bed file of count-based methylation levels from aligned_bam_to_cpg_scores for hap1')
    parser.add_argument('--bed_meth_count_hap2', required=True, help='Bed file of count-based methylation levels from aligned_bam_to_cpg_scores for hap2')
    parser.add_argument('--bed_meth_model_hap1', required=True, help='Bed file of model-based methylation levels from aligned_bam_to_cpg_scores for hap1')
    parser.add_argument('--bed_meth_model_hap2', required=True, help='Bed file of model-based methylation levels from aligned_bam_to_cpg_scores for hap2')
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
    
    write_bit_vector_sites_and_mismatches(df_sites, df_sites_mismatch, args.uid, args.output_dir, logger)
    write_bit_vector_mismatches(args.output_dir, args.uid, df_sites_mismatch, logger)    

    df_meth_count_hap1_hap2 = get_meth_hap1_hap2(
        pb_cpg_tool_mode='count', 
        bed_hap1=args.bed_meth_count_hap1, 
        bed_hap2=args.bed_meth_count_hap2
    )
    logger.info(f"Got read-based phasing of count-based methylation levels: {len(df_meth_count_hap1_hap2)} rows, {len(df_meth_count_hap1_hap2.columns)} columns")
    df_meth_model_hap1_hap2 = get_meth_hap1_hap2(
        pb_cpg_tool_mode='model', 
        bed_hap1=args.bed_meth_model_hap1, 
        bed_hap2=args.bed_meth_model_hap2
    )    
    logger.info(f"Got read-based phasing of model-based methylation levels: {len(df_meth_model_hap1_hap2)} rows, {len(df_meth_model_hap1_hap2.columns)} columns")
    
    df_meth_count_founder_phased = phase_meth_to_founder_haps(df_meth_count_hap1_hap2, df_hap_map)
    logger.info(f"Phased count-based methylation levels to founder haplotypes: {len(df_meth_count_founder_phased)} rows, {len(df_meth_count_founder_phased.columns)} columns")
    df_meth_model_founder_phased = phase_meth_to_founder_haps(df_meth_model_hap1_hap2, df_hap_map)
    logger.info(f"Phased model-based methylation levels to founder haplotypes: {len(df_meth_model_founder_phased)} rows, {len(df_meth_model_founder_phased.columns)} columns")
    df_meth_founder_phased = combine_count_and_model_based_methylation_levels(df_meth_count_founder_phased, df_meth_model_founder_phased)
    logger.info(f"Combined count- and model-based methylation levels: {len(df_meth_founder_phased)} rows, {len(df_meth_founder_phased.columns)} columns")

    write_bed(args.output_dir, df_meth_founder_phased, filename_stem=f"{args.uid}.dna-methylation.founder-phased")
    logger.info(f"Wrote count- and model-based methylation levels, phased to founder haplotypes, to: '{args.output_dir}'")
    
    for parental in ['pat', 'mat']: 
        for pb_cpg_tool_mode in ['count', 'model']:
            write_bigwig(df_meth_founder_phased, args.uid, parental, pb_cpg_tool_mode, args.output_dir, logger)

    logger.info(f"Done running '{__file__}'")

if __name__ == "__main__":
    main()
