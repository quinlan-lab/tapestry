import polars as pl
from pathlib import Path

def get_header(filename):
    with open(filename) as fh: 
        lines = fh.readlines()
        lines = [line.strip() for line in lines]
    return lines

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

def read_bed_and_header(file_path): 
    file_path = Path(file_path)
    parent_dir = file_path.parent
    file_stem = file_path.stem
    file_suffix = file_path.suffix
    assert file_suffix == ".bed"
    return read_data(parent_dir, file_stem)


