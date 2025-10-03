import argparse
import logging
import polars as pl 

from read_data import read_bed_and_header
from write_data import write_dataframe_to_bed
from get_meth_hap1_hap2 import read_meth_level
from version_sort import version_sort

def read_all_cpgs(bed_all_cpgs): 
    return pl.read_csv(
        bed_all_cpgs, 
        separator='\t', 
        new_columns=['chrom', 'start', 'end'],
        # n_rows=100000 # TESTING 
    ) 

def read_meth_founder_phased(bed_meth_founder_phased): 
    return read_bed_and_header(bed_meth_founder_phased)

def read_meth_unphased(bed_meth_unphased, pb_cpg_tool_mode): 
    return (
        read_meth_level(bed_meth_unphased, pb_cpg_tool_mode)
        .rename({
            'chromosome': 'chrom'
        })
    )

def expand_meth_to_all_cpgs(df_all_cpgs, df_meth_founder_phased, df_meth_unphased): 
    df = (
        df_all_cpgs
        .join(
            df_meth_founder_phased,
            on=['chrom', 'start', 'end'],
            how="left"
        )
        .join(
            df_meth_unphased,
            on=['chrom', 'start', 'end'],
            how="left"
        )
    )
    return version_sort(df)

def main(): 
    parser = argparse.ArgumentParser(description='Expand output of tapestry to include all CpG sites')
    parser.add_argument('--bed_all_cpgs', required=True, help='All CpG sites')
    parser.add_argument('--bed_meth_founder_phased', required=True, help='Founder-phased methylation levels')
    parser.add_argument('--bed_meth_unphased', required=True, help='Unphased methylation levels')
    parser.add_argument('--pb_cpg_tool_mode', required=True, help='--pileup-mode argument of aligned_bam_to_cpg_scores')
    parser.add_argument('--bed_meth_founder_phased_all_cpgs', required=True, help='Founder-phased methylation levels at all CpG sites (including null methylation levels)')
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    logger = logging.getLogger(__name__)

    df_all_cpgs = read_all_cpgs(args.bed_all_cpgs)
    logger.info(f"Read all CpG sites from {args.bed_all_cpgs}")

    df_meth_founder_phased = read_meth_founder_phased(args.bed_meth_founder_phased)
    logger.info(f"Read CpG sites at which methylation levels have been phased to founder haplotypes")

    df_meth_unphased = read_meth_unphased(args.bed_meth_unphased, args.pb_cpg_tool_mode) 
    logger.info(f"Read CpG sites at which unphased methylation levels are available")

    df_meth_founder_phased_all_cpgs = expand_meth_to_all_cpgs(df_all_cpgs, df_meth_founder_phased, df_meth_unphased)
    logger.info(f"Expanded DNA methylation to include all CpG sites, and unphased methylation levels (where available)")

    write_dataframe_to_bed(df_meth_founder_phased_all_cpgs, args.bed_meth_founder_phased_all_cpgs)
    logger.info(f"Wrote expanded methylation dataframe to {args.bed_meth_founder_phased_all_cpgs}")

    logger.info(f"Done running {__file__}")

if __name__ == "__main__":
    main()
