from pathlib import Path
import polars as pl 

# https://samtools.github.io/hts-specs/BEDv1.pdf
# "We recommend that only a single tab (\t) be used as field separator."
# "Comment lines start with # with no horizontal whitespace beforehand. A # appearing anywhere else in a data line is treated as feature data, not a comment."
def write_dataframe_to_bed(df: pl.DataFrame, file_path: str):
    """
    Writes a Polars DataFrame to a file in BED format.

    The BED format is defined here as:
    1. Fields are separated by tabs.
    2. The header line is prefixed with a '#' character.

    Args:
        df: The Polars DataFrame to be written to disk.
        file_path: The path of the output file.
    """
    with open(file_path, "w") as f:
        # Create the custom header string and write it to the file
        header_string = "#" + "\t".join(df.columns) + "\n"
        f.write(header_string)

        # Write the DataFrame content directly to the file handler,
        # without a header.
        df.write_csv(f, separator='\t', include_header=False)

def write_header(filename, header):
    with open(filename, 'w') as f:
        f.write('\n'.join(header) + '\n')

def write_data(output_dir, df, filename_stem, suffix='bed'):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    df.write_csv(
        output_dir / f"{filename_stem}.{suffix}",         
        separator='\t', 
        include_header=False
    )
    write_header(
        output_dir / f"{filename_stem}.{suffix}.header", 
        header=df.columns
    )

def write_bed(output_dir, df, filename_stem):
    write_data(output_dir, df, filename_stem, suffix='bed')

def write_bedgraph(output_dir, df, filename_stem):
    write_data(output_dir, df, filename_stem, suffix='bedgraph')

def write_bed_and_header(file_path, df): 
    file_path = Path(file_path)
    parent_dir = file_path.parent
    file_stem = file_path.stem
    file_suffix = file_path.suffix
    assert file_suffix == ".bed"
    write_bed(parent_dir, df, file_stem)

