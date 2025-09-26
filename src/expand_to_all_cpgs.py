import argparse
import logging
import polars as pl 
import bioframe as bf
from pathlib import Path

from read_data import read_data

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)

def read_all_cpgs(bed_all_cpgs): 
    return pl.read_csv(
        bed_all_cpgs, 
        separator='\t', 
        new_columns=['chrom', 'start', 'end'],
        n_rows=10000 # TESTING 
    ) 

def read_hap_map_blocks(bed_hap_map): 
    bed_hap_map = Path(bed_hap_map)
    parent_dir = bed_hap_map.parent
    file_stem = bed_hap_map.stem
    file_suffix = bed_hap_map.suffix
    assert file_suffix == ".bed"
    return (
        read_data(parent_dir, file_stem)
        .drop([
            'paternal_haplotype',
            'maternal_haplotype',
            'haplotype_concordance',
            'num_het_SNVs'
        ])
    )   

def assign_hap_map_blocks_to_cpgs(df_all_cpgs, df_hap_map_blocks): 
    return (
        pl
        .from_pandas(bf.overlap(
            df_all_cpgs.to_pandas(), 
            df_hap_map_blocks.to_pandas(), 
            how='left',
            suffixes=('', '_hap_map_block')
        ))
        .drop([
            'chrom_hap_map_block'
        ])
    )
    
def expand_meth_to_all_cpgs(): 
    # TODO: [EASY] implement this by merging df_all_cpgs_with_hap_map_blocks (see below)
    # with bed_meth_founder_phased on all cols in df_all_cpgs_with_hap_map_blocks, 
    # and writing to disk 
    pass 

def main(args): 
    df_all_cpgs = read_all_cpgs(args.bed_all_cpgs)
    logger.info(f"Read all CpG sites from {args.bed_all_cpgs}.")

    df_hap_map_blocks = read_hap_map_blocks(args.bed_hap_map)
    logger.info(f"Read all hap-map blocks from {args.bed_hap_map}.")

    df_all_cpgs_with_hap_map_blocks = assign_hap_map_blocks_to_cpgs(df_all_cpgs, df_hap_map_blocks)
    logger.info(f"Each CpG site has been assigned an overlapping hap-map block, if such exists.")

    # TODO: call expand_meth_to_all_cpgs here 

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Expand output of tapestry to include all CpG sites')
    parser.add_argument('--bed_all_cpgs', required=True, help='All CpG sites')
    parser.add_argument('--bed_hap_map', required=True, help='Hap-map blocks')
    parser.add_argument('--bed_meth_founder_phased', required=True, help='Founder-phased methylation levels')
    parser.add_argument('--bed_meth_founder_phased_all_cpgs', required=True, help='Founder-phased methylation levels at all CpG sites (including null methylation levels)')
    args = parser.parse_args()
    main(args)
    logger.info(f"Done running {__file__}")

# TODO: need to also pull in combined meth and report that for cpgs that are not in hap-map blocks 
# TOOO: need to run pb-cpg-tools for both "modes" (model and count)