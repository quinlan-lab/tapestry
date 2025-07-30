from util.write_header import write_header
from pathlib import Path

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