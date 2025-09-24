import argparse
import logging
from pathlib import Path
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
from util.write_data import write_bed
import bioframe as bf # https://bioframe.readthedocs.io/en/latest/index.html
import polars as pl

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)

CPG_SITE_MISMATCH_SITE_DISTANCE = 50 # bp
REFERENCE_GENOME = "hg38"

def get_all_phasing_wrapper(uid, vcf_read_phased, tsv_read_phase_blocks, vcf_iht_phased, txt_iht_blocks):
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

def phase_meth_to_founder_haps(df_meth_hap1_hap2, df_hap_map, df_sites_mismatch):
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

    df = (
        pl
        .from_pandas(
            # Use bf.closest() to find the single nearest mismatch for EACH CpG site.
            # The result includes a 'distance' column indicating distance between CpG site and nearest mismatch site. 
            bf.closest(
                df.to_pandas(), 
                df_sites_mismatch.to_pandas()
            )
        )
        .cast({
            "start_hap_map_block": pl.Int64,
            "end_hap_map_block": pl.Int64,
            "num_het_SNVs_in_hap_map_block": pl.Int64,
            "total_read_count_pat": pl.Int64,
            "total_read_count_mat": pl.Int64,
        })
        .with_columns(
            (pl.col('distance') < CPG_SITE_MISMATCH_SITE_DISTANCE)
            .alias(f'is_within_{CPG_SITE_MISMATCH_SITE_DISTANCE}bp_of_mismatch_site'),
        )
        .drop([
            "chrom_", "start_", "end_",
            "REF_", "ALT_",
            "distance"
        ])
        .sort(["chrom", "start", "end"])
    )

    return df 

def write_bigwig(df, uid, parental, output_dir):
    """
    Write a bigwig file for a given parental haplotype.
    """
    
    df_bed_graph = (
        df
        .filter(pl.col(f"methylation_level_{parental}").is_not_null())        
        .select([
            pl.col("chrom"),
            pl.col("start"),
            pl.col("end"),
            pl.col(f"methylation_level_{parental}")
        ])
        .to_pandas()  # Convert to pandas DataFrame for bioframe compatibility
    )

    # https://bioframe.readthedocs.io/en/latest/api-fileops.html#bioframe.io.fileops.to_bigwig
    bf.to_bigwig(
        df_bed_graph, 
        # https://bioframe.readthedocs.io/en/latest/api-resources.html#bioframe.io.resources.fetch_chromsizes 
        bf.fetch_chromsizes(db=REFERENCE_GENOME), 
        outpath=f"{output_dir}/{uid}.dna-methylation.founder-phased.{parental}.bw", 
        path_to_binary="bedGraphToBigWig" # we assume that user has bedGraphToBigWig in their PATH
    )

def main(args):
    Path(args.output_dir).mkdir(parents=True, exist_ok=True)
    
    df_all_phasing = get_all_phasing_wrapper(args.uid, args.vcf_read_phased, args.tsv_read_phase_blocks, args.vcf_iht_phased, args.txt_iht_blocks)
    logger.info(f"Got all phasing data: {len(df_all_phasing)} rows, {len(df_all_phasing.columns)} columns")
    
    df_hap_map, df_sites, df_sites_mismatch = get_hap_map(df_all_phasing)
    logger.info(f"Got hap map: {len(df_hap_map)} rows, {len(df_hap_map.columns)} columns")
    logger.info(f"Got sites: {len(df_sites)} rows, {len(df_sites.columns)} columns")
    logger.info(f"Got sites where read-based and inheritance-based bit vectors don't match: {len(df_sites_mismatch)} rows, {len(df_sites_mismatch.columns)} columns")
    
    write_hap_map_blocks(df_hap_map, args.uid, "paternal", args.output_dir)
    write_hap_map_blocks(df_hap_map, args.uid, "maternal", args.output_dir)
    logger.info("Wrote paternal and maternal hap-map blocks for IGV visualization")
    
    write_bed(args.output_dir, df_hap_map, f"{args.uid}.hap-map-blocks")
    logger.info("Wrote hap-map blocks")
    
    write_bit_vector_sites_and_mismatches(df_sites, df_sites_mismatch, args.uid, args.output_dir)
    logger.info("Wrote sites of bit-vectors, and sites where bit vectors are mismatched, for IGV visualization")
    
    df_meth_hap1_hap2 = get_meth_hap1_hap2(args.pb_cpg_tool_mode, args.bed_meth_hap1, args.bed_meth_hap2)
    logger.info(f"Got read-based phasing of methylation levels: {len(df_meth_hap1_hap2)} rows, {len(df_meth_hap1_hap2.columns)} columns")
    
    df_meth_founder_phased = phase_meth_to_founder_haps(df_meth_hap1_hap2, df_hap_map, df_sites_mismatch)
    logger.info(f"Phased methylation levels to founder haplotypes: {len(df_meth_founder_phased)} rows, {len(df_meth_founder_phased.columns)} columns")

    fraction_of_CpG_sites_that_are_phased_to_founder_haplotypes = (
        df_meth_founder_phased.select(pl.col("founder_haplotype_pat").is_not_null().mean())
    )
    logger.info(f"Percentage of CpG sites that are phased to founder haplotypes: {fraction_of_CpG_sites_that_are_phased_to_founder_haplotypes[0, 0]*100:.0f}%")

    fraction_of_CpG_sites_that_are_near_mismatch_sites = (
        df_meth_founder_phased.select(pl.col(f'is_within_{CPG_SITE_MISMATCH_SITE_DISTANCE}bp_of_mismatch_site').mean())
    )
    logger.info(f"Percentage of CpG sites that are within {CPG_SITE_MISMATCH_SITE_DISTANCE}bp of a mismatch site: {fraction_of_CpG_sites_that_are_near_mismatch_sites[0, 0]*100:.3f}%")

    write_bed(args.output_dir, df_meth_founder_phased, filename_stem=f"{args.uid}.dna-methylation.founder-phased")
    logger.info("Wrote methylation levels phased to founder haplotypes")
    
    write_bigwig(df_meth_founder_phased, args.uid, "pat", args.output_dir)
    logger.info("Wrote bigwig file for paternal methylation levels")
    
    write_bigwig(df_meth_founder_phased, args.uid, "mat", args.output_dir)
    logger.info("Wrote bigwig file for maternal methylation levels")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Phase HiFi-based DNA methylation data to founder haplotypes')
    parser.add_argument('--pb_cpg_tool_mode', required=True, help='Mode of aligned_bam_to_cpg_scores')
    parser.add_argument('--uid', required=True, help='Sample UID in joint-called multi-sample vcf')
    parser.add_argument('--vcf_read_phased', required=True, help='Single-sample vcf from hiphase')
    parser.add_argument('--tsv_read_phase_blocks', required=True, help='Single-sample tsv from hiphase')
    parser.add_argument('--vcf_iht_phased', required=True, help='Joint-called multi-sample vcf from gtg-ped-map/gtg-concordance')
    parser.add_argument('--txt_iht_blocks', required=True, help='Multi-sample iht blocks file from gtg-ped-map/gtg-concordance')
    parser.add_argument('--bed_meth_hap1', required=True, help='Bed file from aligned_bam_to_cpg_scores for hap1')
    parser.add_argument('--bed_meth_hap2', required=True, help='Bed file from aligned_bam_to_cpg_scores for hap2')
    parser.add_argument('--output_dir', required=True, help='Output directory')
    args = parser.parse_args()
    main(args)
    logger.info(f"Done running {__file__}")