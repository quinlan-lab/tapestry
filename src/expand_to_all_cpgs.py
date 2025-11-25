import argparse
import logging
import polars as pl 
import bioframe as bf # https://bioframe.readthedocs.io/en/latest/index.html
from pathlib import Path
from typing import Literal

from read_data import read_bed_and_header
from write_data import write_dataframe_to_bed
from get_meth_hap1_hap2 import read_meth_level
from version_sort import version_sort

# https://quinlangroup.slack.com/archives/C02KEHXJ274/p1724332317643529
# /scratch/ucgd/lustre-labs/quinlan/u6018199/cyvcf2
# https://quinlangroup.slack.com/archives/C449KJT3J/p1751389842484399
from cyvcf2 import VCF  # type: ignore

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
        (pl.col('chrom') != 'chrX') &
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

def compute_fraction_of_cpgs_at_which_umphased_meth_is_reported(df, mode, logger=None): 
    df = reduce_to_phasable_chromosomes(df)
    fraction = (
        df
        .select(pl.col(f"methylation_level_{mode}").is_not_null())
        .mean()
        .item()
    )
    report = f"Percentage of CpG sites (in reference and sample genomes, and on phasable chroms) at which {mode}-based unphased methylation is reported: {fraction*100:.2f}%"
    if logger: 
        logger.info(report)
    else: 
        print(report)

def compute_fraction_of_cpgs_at_which_meth_is_phased_to_given_parent(df, parental, mode, logger=None): 
    df = reduce_to_phasable_chromosomes(df)
    fraction = (
        df
        .select(pl.col(f"methylation_level_{parental}_{mode}").is_not_null())
        .mean()
        .item()
    )
    report = f"Percentage of CpG sites (in reference and sample genomes, and on phasable chroms) at which {mode}-based methylation is phased to {parental} haplotype: {fraction*100:.2f}%"
    if logger: 
        logger.info(report)
    else: 
        print(report)

def create_meth_non_null_expr(
    mode: str, 
    logic: Literal["any", "all"], 
    alias: str
) -> pl.Expr:
    """
    Creates a Polars expression to check for non-null methylation levels.

    Args:
        mode: A string suffix for the column names ('count' or 'model').
        logic: 'any' for OR (|), 'all' for AND (&).
        alias: The name for the resulting boolean column.
    
    Returns:
        A Polars expression.
    """
    pat_col = pl.col(f"methylation_level_pat_{mode}")
    mat_col = pl.col(f"methylation_level_mat_{mode}")
    
    if logic == "any":
        expression = pat_col.is_not_null() | mat_col.is_not_null()
    elif logic == "all":
        expression = pat_col.is_not_null() & mat_col.is_not_null()
    else:
        raise ValueError("`logic` parameter must be 'any' or 'all'")
        
    return expression.alias(alias)

def compute_fraction_of_cpgs_at_which_meth_is_phased(df, mode, logic, alias, logger=None):
    df = reduce_to_phasable_chromosomes(df)
    expr = create_meth_non_null_expr(mode, logic, alias)
    fraction = (
        df
        .with_columns(expr)
        .select(pl.col(alias))
        .mean()
        .item()
    )
    report = f"Percentage of CpG sites (in reference and sample genomes, and on phasable chroms) at which {mode}-based methylation is phased to {alias}: {fraction*100:.2f}%"
    if logger: 
        logger.info(report)
    else: 
        print(report)

def compute_fraction_of_cpgs_at_which_meth_is_phased_wrapper(df, logger=None):
    for mode in ['count', 'model']:
        for parental in ['pat', 'mat']: 
            compute_fraction_of_cpgs_at_which_meth_is_phased_to_given_parent(df, parental, mode, logger)
        compute_fraction_of_cpgs_at_which_meth_is_phased(df, mode, logic='any', alias='at least one parental haplotype', logger=logger)
        compute_fraction_of_cpgs_at_which_meth_is_phased(df, mode, logic='all', alias='both parental haplotypes', logger=logger)
        compute_fraction_of_cpgs_at_which_umphased_meth_is_reported(df, mode, logger)

def is_snv(variant):
    # is_snp allows multiallelic sites (len(variant.ALT) > 1)
    # https://github.com/brentp/cyvcf2/blob/541ab16a255a5287c331843d8180ed6b9ef10e00/cyvcf2/cyvcf2.pyx#L1903-L1911
    return variant.is_snp

def get_iht_phased_variants(uid, vcf): 
    # assume vcf is phased as: "paternal | maternal" 
    # https://quinlangroup.slack.com/archives/C08U7NLC9PZ/p1748885496941579

    # assume vcf is a joint vcf 

    records = []
    with VCF(vcf, strict_gt=True) as vcf_reader: # cyvcf2 handles .vcf.gz directly
        samples = vcf_reader.samples
        sample_index = samples.index(uid) if uid in samples else None

        # for variant in tqdm(vcf_reader, total=vcf_reader.num_records): # testing 
        for variant in vcf_reader:
        # for i, variant in enumerate(vcf_reader): # testing
        #     if i >= 100: # testing
        #         break # testing

            if not is_snv(variant): 
                continue

            chrom = variant.CHROM
            pos = variant.POS # pos is 1-based
            start = pos - 1 
            end = pos
            REF = variant.REF
            ALT = variant.ALT

            # single sample of 0|1 in vcf becomes [[0, 1, True]]
            # 2 samples of 0/0 and 1|1 would be [[0, 0, False], [1, 1, True]]
            # https://brentp.github.io/cyvcf2/#cyvcf2
            genotype_all_samples = variant.genotypes
            genotype = genotype_all_samples[sample_index]

            # "-1" in "genotype" indicates missing data: 
            # c.f., "test_set_gts" at: https://github.com/brentp/cyvcf2/issues/31#issuecomment-275195917
            allele_pat = str(genotype[0]) if genotype[0] != -1 else '.'
            allele_mat = str(genotype[1]) if genotype[1] != -1 else '.'

            phased = genotype[2]

            if not phased: 
                raise ValueError(f"Expected phased genotype, but found unphased: {genotype}")

            records.append({ 
                "chrom": chrom,
                "start": start,
                "end": end,
                "REF": REF,
                "ALT": ALT,
                "allele_pat": allele_pat,
                "allele_mat": allele_mat,
            })
 
    df = pl.DataFrame(records)
    return df 

def label_with_variants(df_meth_founder_phased_all_cpgs, df_iht_phased_variants):
    # enlarge interval defining CpG sites from C nucleotide to CG dinucleotide
    df_meth_founder_phased_all_cpgs_dinucleotides = df_meth_founder_phased_all_cpgs.clone()
    df_meth_founder_phased_all_cpgs_dinucleotides = df_meth_founder_phased_all_cpgs_dinucleotides.with_columns((pl.col("end") + 1))

    df = bf.overlap(
        df_meth_founder_phased_all_cpgs_dinucleotides.to_pandas(),
        df_iht_phased_variants.to_pandas(),
        how='left',
        suffixes=('', '_variant')
    )
    df = (
        pl
        .from_pandas(df)
        .cast({
            "total_read_count": pl.Int64,
            "start_hap_map_block": pl.Int64,
            "end_hap_map_block": pl.Int64,
            "num_het_SNVs_in_hap_map_block": pl.Int64,
            "total_read_count_pat": pl.Int64,
            "total_read_count_mat": pl.Int64,
        })
    )

    # determine the number of overlapping SNVs for each cpg record 
    CG_cols = df_meth_founder_phased_all_cpgs_dinucleotides.columns 
    df = df.with_columns(
        pl
        .when(pl.col("start_variant").is_null())
        .then(pl.lit(0))
        .otherwise(pl.col("start_variant").count().over(CG_cols)) # The .over() expression is Polars' implementation of Window Functions
        .alias("num_SNVs_overlapping_CG")
    )    

    df = (
        df
        .rename({
            'start': 'start_cpg',
            'end': 'end_cpg',
            'is_within_50bp_of_mismatch_site': 'cpg_is_within_50bp_of_mismatch_site',
            'REF_variant': 'REF',
            'ALT_variant': 'ALT',
            'allele_pat_variant': 'allele_pat',
            'allele_mat_variant': 'allele_mat'
        })
        .drop(['chrom_variant'])
    )

    return df 

def label_with_imprinting_flag(df): 
    CG_cols = [
        col for col in df.columns 
        if col not in [
            'start_variant', 
            'end_variant', 
            'REF', 
            'ALT',	
            'allele_pat', 
            'allele_mat', 
            'num_SNVs_overlapping_CG'
        ]
    ]
    return ( 
        df
        .with_columns(
            pl
            .when(pl.col("allele_pat").is_null())
            .then(None)
            .when(pl.col("allele_pat") == pl.col("allele_mat"))
            .then(pl.lit("hom"))
            .otherwise(pl.lit("het"))
            .alias("genotype")
        )
        .group_by(CG_cols) # there are two records in df for each CG that overlaps 2 SNVs 
        .agg([
            pl.col("genotype").str.join(",").alias("genotypes"),
            (pl.col("genotype") == "het").any().not_().alias("include_for_imprinting"),
        ])
        .sort(['chrom', 'start_cpg'])  
    )

def main(): 
    parser = argparse.ArgumentParser(description='Expand output of tapestry to include all CpG sites and unphased DNA methylation levels')
    parser.add_argument('--bed_all_cpgs_in_reference', required=True, help='All CpG sites observed in reference genome')
    parser.add_argument('--bed_meth_count_unphased', required=True, help='Unphased count-based methylation levels')
    parser.add_argument('--bed_meth_model_unphased', required=True, help='Unphased model-based methylation levels')
    parser.add_argument('--bed_meth_founder_phased', required=True, help='Founder-phased methylation levels')
    parser.add_argument('--bed_het_site_mismatches', required=True, help='Heterozygous sites where bit vectors are mismatched')
    parser.add_argument('--bed_meth_founder_phased_all_cpgs', required=True, help='Founder-phased methylation levels at all CpG sites, both in reference and sample, including null methylation levels, and unphased methylation levels')
    parser.add_argument('--uid', required=True, help='Sample UID in joint-called multi-sample vcf')
    parser.add_argument('--vcf_iht_phased', required=True, help='Joint-called multi-sample vcf from gtg-ped-map/gtg-concordance')
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

    if Path(args.bed_meth_founder_phased).exists():
        df_meth_founder_phased = read_meth_founder_phased(args.bed_meth_founder_phased)
        logger.info(f"Read CpG sites at which count- and model-based methylation levels have been phased to founder haplotypes")
    else: 
        logger.warning(f"Could not read CpG sites at which count- and model-based methylation levels have been phased to founder haplotypes") 
        logger.warning(f"Required file does not exist: '{args.bed_meth_founder_phased}'")
        logger.warning(f"This may be because this sample is a founder and therefore cannot be inheritance-based phased")
        logger.info(f"Done running '{__file__}'")
        return 

    df_meth_founder_phased_all_cpgs = expand_meth_to_all_cpgs(df_all_cpgs_in_reference, df_meth_unphased, df_meth_founder_phased)
    logger.info(f"Expanded DNA methylation dataset to include (1) all CpG sites observed in reference and sample genomes, and (2) unphased methylation levels (where available)")

    df_meth_founder_phased_all_cpgs = compute_proximity_to_mismatched_heterozygous_sites(df_meth_founder_phased_all_cpgs, args.bed_het_site_mismatches)
    logger.info(f"Computed proximity of all CpG sites to heterozygous sites at which bit-vectors are mismatched")

    compute_fraction_of_cpgs_that_are_close_to_mismatches(df_meth_founder_phased_all_cpgs, logger)
    compute_fraction_of_cpgs_at_which_meth_is_phased_wrapper(df_meth_founder_phased_all_cpgs, logger)

    df_iht_phased_variants = get_iht_phased_variants(args.uid, args.vcf_iht_phased)
    logger.info(f"Got inheritance-phased SNVs from joint vcf")

    df_meth_founder_phased_all_cpgs_with_variant_label = label_with_variants(df_meth_founder_phased_all_cpgs, df_iht_phased_variants)
    logger.info(f"Determined which CpG sites overlap 1 or 2 SNVs")

    df_meth_founder_phased_all_cpgs_with_imprinting_flag = label_with_imprinting_flag(df_meth_founder_phased_all_cpgs_with_variant_label) 
    logger.info(f"Flagged CpG sites that have been created or destroyed by assessing overlap with SNVs, e.g., for use in scanning the genome for imprinted loci across a pedigree")

    write_dataframe_to_bed(
        df_meth_founder_phased_all_cpgs_with_imprinting_flag, 
        args.bed_meth_founder_phased_all_cpgs, 
        source=f"{__file__} with args {vars(args)}"
    )
    logger.info(f"Wrote expanded and imprinting-flagged methylation dataframe to: '{args.bed_meth_founder_phased_all_cpgs}'")

    logger.info(f"Done running '{__file__}'")

if __name__ == "__main__":
    main()
