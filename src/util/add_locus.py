import polars as pl 

def add_locus(df, col_name="locus"):
    return df.with_columns(
        pl
        .format(
            "{}:{}-{}", 
            pl.col("chrom"),
            pl.col("start"),
            pl.col("end")
        )
        .alias(col_name)
    )

