import polars as pl 

def version_sort(df):
    return (
        df
        .sort(
            # 1. Sort based on a new expression, extracting only the numbers with `str.extract`, and casting to `int`
            pl.col("chrom").str.extract(r'(\d+)').cast(pl.Int64, strict=False),
            # 2. Sort a second time based on the string values themselves to place "chrX" and "chrY" correctly.
            "chrom", 
            "start",
            nulls_last=True
        )
    )

