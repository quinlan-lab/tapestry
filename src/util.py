import polars as pl


def read_all_cpgs_in_reference(bed):
    df = pl.read_csv(
        bed,
        separator='\t',
        has_header=False,
        new_columns=['chrom', 'start', 'end'],
        # n_rows=100000 # TESTING
    )
    return df
