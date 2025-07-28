from write_header import write_header

def write_data(data_dir, df, filename_stem, suffix='bed'):
    data_dir.mkdir(parents=True, exist_ok=True)
    df.write_csv(
        data_dir / f"{filename_stem}.{suffix}",         
        separator='\t', 
        include_header=False
    )
    write_header(
        data_dir / f"{filename_stem}.{suffix}.header", 
        header=df.columns
    )

def write_bed(data_dir, df, filename_stem):
    write_data(data_dir, df, filename_stem, suffix='bed')

def write_bedgraph(data_dir, df, filename_stem):
    write_data(data_dir, df, filename_stem, suffix='bedgraph')