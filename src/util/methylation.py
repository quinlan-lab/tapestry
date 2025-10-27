import bioframe as bf # https://bioframe.readthedocs.io/en/latest/index.html
import polars as pl

from version_sort import version_sort

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
        pl.len().alias("num_cpgs"),
    ]

    level_types = ['count', 'model']
    alleles = ['pat', 'mat']

    for allele in alleles: 
        expression = pl.col(f"founder_haplotype_{allele}").unique().alias(f"founder_{allele}")
        expressions.append(expression)
    
    for level_type in level_types:
        column_name = f"methylation_level_{level_type}"
        expressions.extend([
            # .mean() works by summing all the non-null values and dividing by the count of those non-null values
            pl.col(column_name).mean().alias(f"{level_type}_based_meth"), 
            # pl.col(column_name).count().alias(f"num_cpgs_with_non_null_{level_type}_based_meth"),
        ])

    for allele in alleles:
        for level_type in level_types:
            expression = pl.col(f"methylation_level_{allele}_{level_type}").mean().alias(f"{level_type}_based_meth_{allele}")
            expressions.append(expression)
            
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

def compute_delta_methylation(df): 
    # List of the metric types you want to calculate a delta for
    metrics_to_diff = [
        "count_based_meth",
        "model_based_meth"
    ]
    delta_expressions = []
    for metric in metrics_to_diff:
        expr = (pl.col(f"{metric}_pat") - pl.col(f"{metric}_mat")).alias(f"delta_of_{metric}")
        delta_expressions.append(expr)
    df = df.with_columns(delta_expressions)

    df = df.drop(['num_cpgs'])

    return df 

def test_polars_expressions(): 
    # Generate the list of expressions by calling the function.
    aggregation_expressions = generate_methylation_expressions()

    # Print the generated expressions to see what they look like.
    print("--- Generated Polars Expressions ---")
    for expr in aggregation_expressions:
        print(expr)

if __name__ == '__main__':
    test_polars_expressions()
