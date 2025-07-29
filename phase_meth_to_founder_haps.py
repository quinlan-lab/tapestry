import argparse
import logging
from pathlib import Path
from get_all_phasing import (
    get_read_based_phasing, 
    get_phase_blocks, 
    get_parental_phasing, 
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

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)

CPG_SITE_MISMATCH_SITE_DISTANCE = 50 # bp
REFERENCE_GENOME = "hg38"

# Constants from the notebook
READ_BACKED_PHASED_DIR = Path('/scratch/ucgd/lustre-labs/quinlan/data-shared/read-backed-phasing')
HAPLOTYPE_MAPS_DIR = Path('/scratch/ucgd/lustre-labs/quinlan/data-shared/haplotype-maps/CEPH1463.GRCh38')
PB_CPG_TOOL_MODE = 'model'
METH_READ_BACKED_PHASED_DIR = Path(f'/scratch/ucgd/lustre-labs/quinlan/data-shared/dna-methylation/CEPH1463.GRCh38.hifi.{PB_CPG_TOOL_MODE}.read-backed-phased')
METH_FOUNDER_PHASED_DIR = Path(f'/scratch/ucgd/lustre-labs/quinlan/data-shared/dna-methylation/CEPH1463.GRCh38.hifi.{PB_CPG_TOOL_MODE}.founder-phased') 

def get_all_phasing_from_uid(uid):
    """Wrapper function that creates all required DataFrames and calls get_all_phasing"""

    # Get read-based phasing data
    df_read_based_phasing = get_read_based_phasing(
        uid=uid, 
        read_backed_phased_dir=READ_BACKED_PHASED_DIR, 
    )
    logger.info(f"Got read-based phasing data: {len(df_read_based_phasing)} rows, {len(df_read_based_phasing.columns)} columns")
    
    # Get read-based phase blocks
    df_phase_blocks = get_phase_blocks(
        uid=uid, 
        read_backed_phased_dir=READ_BACKED_PHASED_DIR
    )
    logger.info(f"Got phase blocks: {len(df_phase_blocks)} rows, {len(df_phase_blocks.columns)} columns")

    # Get inheritance-based phasing data
    df_parental_phasing = get_parental_phasing(
        uid=uid, 
        haplotype_maps_dir=HAPLOTYPE_MAPS_DIR, 
    )
    logger.info(f"Got parental phasing data: {len(df_parental_phasing)} rows, {len(df_parental_phasing.columns)} columns")

    # Get IHT blocks
    df_iht_blocks = get_iht_blocks(
        uid=uid, 
        haplotype_maps_dir=HAPLOTYPE_MAPS_DIR
    )
    logger.info(f"Got IHT blocks: {len(df_iht_blocks)} rows, {len(df_iht_blocks.columns)} columns")

    # Combine all phasing data
    df_all_phasing = get_all_phasing(
        df_read_based_phasing, 
        df_phase_blocks, 
        df_parental_phasing, 
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

def write_bigwig(df, uid, parental, meth_founder_phased_dir):
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
        outpath=f"{meth_founder_phased_dir}/{uid}.dna-methylation.founder-phased.{parental}.bw", 
        path_to_binary="bedGraphToBigWig" # we assume that user has bedGraphToBigWig in their PATH
    )

def main(uid):
    # Create output directory if it doesn't exist
    METH_FOUNDER_PHASED_DIR.mkdir(parents=True, exist_ok=True)
    
    df_all_phasing = get_all_phasing_from_uid(uid)
    logger.info(f"Got all phasing data: {len(df_all_phasing)} rows, {len(df_all_phasing.columns)} columns")
    
    df_hap_map, df_sites, df_sites_mismatch = get_hap_map(df_all_phasing)
    logger.info(f"Got hap map: {len(df_hap_map)} rows, {len(df_hap_map.columns)} columns")
    logger.info(f"Got sites: {len(df_sites)} rows, {len(df_sites.columns)} columns")
    logger.info(f"Got sites mismatch: {len(df_sites_mismatch)} rows, {len(df_sites_mismatch.columns)} columns")
    
    write_hap_map_blocks(df_hap_map, uid, "paternal", METH_FOUNDER_PHASED_DIR)
    write_hap_map_blocks(df_hap_map, uid, "maternal", METH_FOUNDER_PHASED_DIR)
    logger.info("Wrote paternal and maternal hap-map blocks for IGV visualization")
    
    write_bed(METH_FOUNDER_PHASED_DIR, df_hap_map, f"{uid}.hap-map-blocks")
    logger.info("Wrote hap-map blocks")
    
    write_bit_vector_sites_and_mismatches(df_sites, df_sites_mismatch, uid, METH_FOUNDER_PHASED_DIR)
    logger.info("Wrote bit-vector sites and mismatches for IGV visualization")
    
    df_meth_hap1_hap2 = get_meth_hap1_hap2(uid=uid, pb_cpg_tool_mode=PB_CPG_TOOL_MODE, meth_read_backed_phased_dir=METH_READ_BACKED_PHASED_DIR)
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

    write_bed(METH_FOUNDER_PHASED_DIR, df_meth_founder_phased, filename_stem=f"{uid}.dna-methylation.founder-phased")
    logger.info("Wrote methylation levels phased to founder haplotypes")
    
    write_bigwig(df_meth_founder_phased, uid, "pat", METH_FOUNDER_PHASED_DIR)
    logger.info("Wrote bigwig file for paternal methylation levels")
    
    write_bigwig(df_meth_founder_phased, uid, "mat", METH_FOUNDER_PHASED_DIR)
    logger.info("Wrote bigwig file for maternal methylation levels")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Phase DNA methylation data to founder haplotypes')
    parser.add_argument('-uid', required=True, help='Sample UID to process')
    args = parser.parse_args()

    main(args.uid)