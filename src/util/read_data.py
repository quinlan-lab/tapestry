import polars as pl
from pathlib import Path
import gzip 

def get_header(filename):
    with open(filename) as fh: 
        lines = fh.readlines()
        lines = [line.strip() for line in lines]
    return lines

def read_bed_and_header_core(data_dir, filename_stem): 
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
    return read_bed_and_header_core(parent_dir, file_stem)

def read_tapestry(bed) -> pl.DataFrame:
    """
    Reads the final output file of tapestry into a Polars DataFrame.
    
    The function asserts that the file contains header lines starting with '##' 
    and a single header line starting with '#'. 
    """
    is_gzipped = str(bed).endswith('.gz')
    _open = gzip.open if is_gzipped else open

    has_comment_lines = False
    has_header_line = False

    with _open(bed, 'rt') as f:
        for line in f:            
            if line.startswith('##'):
                has_comment_lines = True
            elif line.startswith('#'):
                has_header_line = True
                break
            else:
                # Data line reached before header
                break
    
    assert has_comment_lines, f"File {bed} is missing comment lines starting with '##'"
    assert has_header_line, f"File {bed} is missing a header line starting with '#'"

    df = (
        pl
        .read_csv(
            bed,
            separator='\t',
            comment_prefix='##',
            has_header=True,
            # The header line starts with '#', which polars handles automatically.
        )
        .rename({
            '#chrom': 'chrom',
        })
    )

    return df

