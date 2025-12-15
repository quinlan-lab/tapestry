import bioframe as bf # https://bioframe.readthedocs.io/en/latest/index.html
import polars as pl
from tqdm import tqdm 
from pathlib import Path 

from version_sort import version_sort
from get_palladium_prefixes import get_prefixes_wrapper
from read_data import read_tapestry
from prefix_columns import prefix_columns

def generate_methylation_expressions():
    """
    Generates a generalized list of Polars expressions for aggregating methylation data.

    This function creates expressions to calculate the mean and count for various
    methylation level columns, looping through different alleles ('pat', 'mat')
    and level types ('count', 'model').

    Returns:
        A list of Polars expressions that can be used in an .agg() clause.
    """
    expressions = [
        # This counts total rows in the group (including nulls)
        pl.len().alias("num_cpgs"),
    ]

    level_types = ['count', 'model']
    alleles = ['pat', 'mat']

    for allele in alleles: 
        expressions.append(
            pl.col(f"founder_haplotype_{allele}").unique().alias(f"founder_{allele}")
        )
    
    for level_type in level_types:
        column_name = f"methylation_level_{level_type}"

        # Always calculate mean
        expressions.append(
            # .mean() works by summing all the non-null values and dividing by the count of those non-null values
            pl.col(column_name).mean().alias(f"{level_type}_based_meth")
        )

        # Only calculate non-null count if level_type is 'count'
        if level_type == 'count':
            expressions.append(
                # .count() only counts non-null values
                pl.col(column_name).count().alias(f"num_cpgs_with_non_null_{level_type}_based_meth")
            )        

    for allele in alleles:
        for level_type in level_types:
            expressions.append(
                pl.col(f"methylation_level_{allele}_{level_type}").mean().alias(f"{level_type}_based_meth_{allele}")
            )
            
    return expressions

def compute_methylation(df_intervals, df_meth, aggregation_expressions=generate_methylation_expressions()): 
    assert df_intervals.columns == ['chrom', 'start', 'end']

    df_intersected = bf.overlap(
        df_intervals.to_pandas(),
        df_meth.to_pandas(),
        how='inner', # "inner" is sufficient when df_intervals covers the genome 
        suffixes=('_intervals', ''),
        return_overlap=False
    )
    df_intersected = pl.from_pandas(df_intersected)
    
    return version_sort(
        df_intersected
        .group_by(["chrom_intervals", "start_intervals", "end_intervals"])
        .agg(aggregation_expressions)
        .rename({
            "chrom_intervals": "chrom", 
            "start_intervals": "start", 
            "end_intervals": "end"
        }) 
    )

def compute_methylation_all_samples_at_given_loci(df_loci, meth_read_phased_dir): 
    df_loci = df_loci.select(['chrom', 'start', 'end'])

    prefixes = get_prefixes_wrapper()
    # prefixes = prefixes[:2] # TESTING
    df_all_samples = None
    for prefix in tqdm(prefixes):
        bed_meth = f"{meth_read_phased_dir}/{prefix}.dna-methylation.founder-phased.all_cpgs.bed"
        if Path(bed_meth).exists():
            df_meth = read_tapestry(bed_meth)
        else: 
            print(f"Could not read CpG sites at which count- and model-based methylation levels have been phased to founder haplotypes") 
            print(f"Required file does not exist: '{bed_meth}'")
            print(f"This may be because this sample is a founder and therefore cannot be inheritance-based phased")
            continue

        df_loci_with_meth = compute_methylation(df_loci, df_meth)
        df_loci_with_meth = df_loci_with_meth.drop([
            'num_cpgs', 
            'count_based_meth', 
            'model_based_meth', 
            'founder_pat', 
            'founder_mat'
        ])
        join_keys = ['chrom', 'start', 'end']
        df_loci_with_meth = prefix_columns(df_loci_with_meth, prefix=prefix, join_keys=join_keys)

        if df_all_samples is None:
            df_all_samples = df_loci_with_meth
        else:
            df_all_samples = df_all_samples.join(
                df_loci_with_meth,
                on=join_keys,
                # capture loci in which at least one of the two dfs has a record: 
                how="full", 
                coalesce=True, 
            )   

    return version_sort(df_all_samples)

def test_polars_expressions(): 
    # Generate the list of expressions by calling the function.
    aggregation_expressions = generate_methylation_expressions()

    # Print the generated expressions to see what they look like.
    print("--- Generated Polars Expressions ---")
    for expr in aggregation_expressions:
        print(expr)

if __name__ == '__main__':
    test_polars_expressions()
