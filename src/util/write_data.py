from pathlib import Path

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

