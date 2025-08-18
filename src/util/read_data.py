import polars as pl
from get_header import get_header

def read_data(data_dir, filename_stem): 
    records_filename = data_dir / f"{filename_stem}.bed"      
    header_filename = data_dir / f"{filename_stem}.bed.header"
    df = pl.read_csv(
        records_filename,
        separator='\t',
        has_header=False,
        new_columns=get_header(header_filename),
        infer_schema_length=1000000,
        # n_rows=100000  # TODO: testing
    )
    return df
