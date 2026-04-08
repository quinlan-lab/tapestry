import polars as pl


def report_size(df, df_name, logger):
    """Log the estimated memory size of a polars DataFrame in GB."""
    size_in_gb = df.estimated_size(unit="gb")
    logger.info(f"{df_name}: {size_in_gb:.4f} GB")
