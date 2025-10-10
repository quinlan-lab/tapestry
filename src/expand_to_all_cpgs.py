import argparse
import logging
import polars as pl 
import bioframe as bf # https://bioframe.readthedocs.io/en/latest/index.html

from read_data import read_bed_and_header
from write_data import write_dataframe_to_bed
from get_meth_hap1_hap2 import read_meth_level
from version_sort import version_sort

CPG_SITE_MISMATCH_SITE_DISTANCE = 50 # bp

def read_all_cpgs_in_reference(bed): 
    df = pl.read_csv(
        bed, 
        separator='\t', 
        has_header=False,
        new_columns=['chrom', 'start', 'end'],
        # n_rows=100000 # TESTING 
    ) 
    return version_sort(df)

def read_meth_unphased(bed_meth_count_unphased, bed_meth_model_unphased): 
    df_meth_count_unphased = (
        read_meth_level(bed_meth_count_unphased, pb_cpg_tool_mode='count')
        .rename({
            'chromosome': 'chrom',
            'total_read_count': 'total_read_count_count',
            'methylation_level': 'methylation_level_count',
        })
    )
    df_meth_model_unphased = (
        read_meth_level(bed_meth_model_unphased, pb_cpg_tool_mode='model')
        .rename({
            'chromosome': 'chrom',
            'total_read_count': 'total_read_count_model',
            'methylation_level': 'methylation_level_model',
        })
    )
    df = df_meth_count_unphased.join(
        df_meth_model_unphased,
        on=['chrom', 'start', 'end'],
        join_nulls=True,
        # "how=full" will capture CpG sites where count-based meth levels are available, but not model-based meth levels, and vice versa.
        # This is also why I did not join on "total_read_count", e.g., this could be null for count-based meth and non-null for model-based meth in the same record. 
        how="full", 
        coalesce=True, 
    )   
    return version_sort(df)

def read_meth_founder_phased(bed_meth_founder_phased): 
    return read_bed_and_header(bed_meth_founder_phased)

def expand_meth_to_all_cpgs(df_all_cpgs_in_reference, df_meth_unphased, df_meth_founder_phased): 
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
            df_meth_founder_phased,
            on=['chrom', 'start', 'end'],
            # "how=left" is appropriate 
            # because the set of sites in df_meth_founder_phased is a subset of the set of sites in df_all_cpgs
            how="left"
        )
    )

    if df['total_read_count_count'].equals(df['total_read_count_model']):
        df = (
            df
            .drop(['total_read_count_model'])
            .rename({
                'total_read_count_count': 'total_read_count'
            })
        )
    else: 
        raise ValueError('total_read_count_count is not identical to total_read_count_model')
    
    return version_sort(df)

def compute_proximity_to_mismatched_heterozygous_sites(df_meth_founder_phased_all_cpgs, bed_het_site_mismatches): 
    df_sites_mismatch = read_bed_and_header(bed_het_site_mismatches)

    df_meth_founder_phased_all_cpgs = (
        pl
        .from_pandas(
            # Use bf.closest() to find the single nearest mismatch site for EACH input site.
            # The result includes a 'distance' column indicating distance between input site and nearest mismatch site. 
            bf.closest(
                df_meth_founder_phased_all_cpgs.to_pandas(), 
                df_sites_mismatch.to_pandas()
            )
        )
        .cast({
            "total_read_count": pl.Int64,
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
    )

    return version_sort(df_meth_founder_phased_all_cpgs)

def reduce_to_phasable_chromosomes(df): 
    df = df.filter(
        (pl.col('chrom') != 'chrM') &
        (pl.col('chrom') != 'chrY')
    )
    assert df[f"is_within_{CPG_SITE_MISMATCH_SITE_DISTANCE}bp_of_mismatch_site"].is_not_null().all()
    return df 

def compute_fraction_of_cpgs_that_are_close_to_mismatches(df, logger=None):
    df = reduce_to_phasable_chromosomes(df)
    fraction = df[f"is_within_{CPG_SITE_MISMATCH_SITE_DISTANCE}bp_of_mismatch_site"].mean()
    report = f"Percentage of CpG sites (in reference and sample genome, and on phasable chroms) that are within {CPG_SITE_MISMATCH_SITE_DISTANCE}bp of a heterozygous mismatch site: {fraction*100:.3f}%"
    if logger: 
        logger.info(report)
    else: 
        print(report)

def compute_fraction_of_cpgs_at_which_meth_is_phased(df, parental, mode, logger=None): 
    df = reduce_to_phasable_chromosomes(df)
    fraction = df.select(pl.col(f"methylation_level_{parental}_{mode}").is_not_null().mean())[0,0]
    report = f"Percentage of CpG sites (in reference and sample genomes, and on phasable chroms) at which {mode}-based methylation is phased to {parental} haplotype: {fraction*100:.2f}%"
    if logger: 
        logger.info(report)
    else: 
        print(report)

def main(): 
    parser = argparse.ArgumentParser(description='Expand output of tapestry to include all CpG sites and unphased DNA methylation levels')
    parser.add_argument('--bed_all_cpgs_in_reference', required=True, help='All CpG sites observed in reference genome')
    parser.add_argument('--bed_meth_count_unphased', required=True, help='Unphased count-based methylation levels')
    parser.add_argument('--bed_meth_model_unphased', required=True, help='Unphased model-based methylation levels')
    parser.add_argument('--bed_meth_founder_phased', required=True, help='Founder-phased methylation levels')
    parser.add_argument('--bed_het_site_mismatches', required=True, help='Heterozygous sites where bit vectors are mismatched')
    parser.add_argument('--bed_meth_founder_phased_all_cpgs', required=True, help='Founder-phased methylation levels at all CpG sites, both in reference and sample, including null methylation levels, and unphased methylation levels')
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(filename)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    logger = logging.getLogger(__name__)

    logger.info(f"Starting '{__file__}'")
    logger.info("Script started with the following arguments: %s", vars(args))

    df_all_cpgs_in_reference = read_all_cpgs_in_reference(args.bed_all_cpgs_in_reference)
    logger.info(f"Read all CpG sites in reference genome")

    df_meth_unphased = read_meth_unphased(args.bed_meth_count_unphased, args.bed_meth_model_unphased) 
    logger.info(f"Read CpG sites at which unphased methylation levels are available, both count-based and model-based")

    df_meth_founder_phased = read_meth_founder_phased(args.bed_meth_founder_phased)
    logger.info(f"Read CpG sites at which count- and model-based methylation levels have been phased to founder haplotypes")

    df_meth_founder_phased_all_cpgs = expand_meth_to_all_cpgs(df_all_cpgs_in_reference, df_meth_unphased, df_meth_founder_phased)
    logger.info(f"Expanded DNA methylation dataset to include (1) all CpG sites observed in reference and sample genomes, and (2) unphased methylation levels (where available)")

    df_meth_founder_phased_all_cpgs = compute_proximity_to_mismatched_heterozygous_sites(df_meth_founder_phased_all_cpgs, args.bed_het_site_mismatches)
    logger.info(f"Computed proximity of all CpG sites to heterozygous sites at which bit-vectors are mismatched")

    write_dataframe_to_bed(df_meth_founder_phased_all_cpgs, args.bed_meth_founder_phased_all_cpgs, source=f"{__file__} with args {vars(args)}")
    logger.info(f"Wrote expanded methylation dataframe to: '{args.bed_meth_founder_phased_all_cpgs}'")

    compute_fraction_of_cpgs_that_are_close_to_mismatches(df_meth_founder_phased_all_cpgs, logger)

    for parental in ['pat', 'mat']: 
        for mode in ['count', 'model']:
            compute_fraction_of_cpgs_at_which_meth_is_phased(df_meth_founder_phased_all_cpgs, parental, mode, logger)

    logger.info(f"Done running '{__file__}'")

if __name__ == "__main__":
    main()
