import argparse
import logging
import polars as pl 
import bioframe as bf
from pathlib import Path

from read_data import read_data
from write_data import write_bed

def read_all_cpgs(bed_all_cpgs): 
    return pl.read_csv(
        bed_all_cpgs, 
        separator='\t', 
        new_columns=['chrom', 'start', 'end'],
        # n_rows=100000 # TESTING 
    ) 

def read_bed_and_header(file_path): 
    file_path = Path(file_path)
    parent_dir = file_path.parent
    file_stem = file_path.stem
    file_suffix = file_path.suffix
    assert file_suffix == ".bed"
    return read_data(parent_dir, file_stem)

def write_bed_and_header(file_path, df): 
    file_path = Path(file_path)
    parent_dir = file_path.parent
    file_stem = file_path.stem
    file_suffix = file_path.suffix
    assert file_suffix == ".bed"
    write_bed(parent_dir, df, file_stem)

def read_hap_map_blocks(bed_hap_map): 
    return (
        read_bed_and_header(bed_hap_map)
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

def read_meth_founder_phased(bed_meth_founder_phased): 
    return read_bed_and_header(bed_meth_founder_phased)

def expand_meth_to_all_cpgs(df_all_cpgs_with_hap_map_blocks, df_meth_founder_phased, bed_meth_founder_phased_all_cpgs): 
    return (
        df_all_cpgs_with_hap_map_blocks
        .join(
            df_meth_founder_phased,
            on=df_all_cpgs_with_hap_map_blocks.columns, 
            how="left"
        )
    )

def main(): 
    parser = argparse.ArgumentParser(description='Expand output of tapestry to include all CpG sites')
    parser.add_argument('--bed_all_cpgs', required=True, help='All CpG sites')
    parser.add_argument('--bed_hap_map', required=True, help='Hap-map blocks')
    parser.add_argument('--bed_meth_founder_phased', required=True, help='Founder-phased methylation levels')
    parser.add_argument('--bed_meth_founder_phased_all_cpgs', required=True, help='Founder-phased methylation levels at all CpG sites (including null methylation levels)')
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    logger = logging.getLogger(__name__)

    df_all_cpgs = read_all_cpgs(args.bed_all_cpgs)
    logger.info(f"Read all CpG sites from {args.bed_all_cpgs}.")

    df_hap_map_blocks = read_hap_map_blocks(args.bed_hap_map)
    logger.info(f"Read all hap-map blocks from {args.bed_hap_map}.")

    df_all_cpgs_with_hap_map_blocks = assign_hap_map_blocks_to_cpgs(df_all_cpgs, df_hap_map_blocks)
    logger.info(f"Each CpG site has been assigned an overlapping hap-map block, if such exists.")

    df_meth_founder_phased = read_meth_founder_phased(args.bed_meth_founder_phased)
    logger.info(f"Read CpG sites at which methylation levels have been phased to founder haplotypes")

    df_meth_founder_phased_all_cpgs = expand_meth_to_all_cpgs(df_all_cpgs_with_hap_map_blocks, df_meth_founder_phased, args.bed_meth_founder_phased_all_cpgs)
    logger.info(f"Expanded DNA methylation to all CpG sites.")

    write_bed_and_header(args.bed_meth_founder_phased_all_cpgs, df_meth_founder_phased_all_cpgs)
    logger.info(f"Wrote expanded methylation dataframe to {args.bed_meth_founder_phased_all_cpgs}")

    logger.info(f"Done running {__file__}")

if __name__ == "__main__":
    main()
