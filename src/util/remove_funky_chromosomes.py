import sys
import polars as pl

def funky(chrom): 
    return "_random" in chrom or "chrEBV" in chrom or "chrUn_" in chrom

def remove_funky_chromosomes(df, chrom_column):
    """
    Removes rows where the chromosome name indicates a non-standard contig 
    e.g., chr1_KI270706v1_random, chrEBV, chrUn_GL000195v1  
    """
    # using .filter() with pl.col() expressions allows us to handle DataFrames and LazyFrames:
    return df.filter(
        ~pl.col(chrom_column).str.contains("_random") &
        ~pl.col(chrom_column).str.contains("chrEBV") &
        ~pl.col(chrom_column).str.contains("chrUn_") 
    )

def main(): 
    df = pl.read_csv(
        source=sys.argv[1],
        separator='\t',
        has_header=False, # assume no header
        new_columns=['chrom', 'start', 'end'], # rename columns
        infer_schema_length=1000000,
    )
    old_number_records = len(df)
    df = remove_funky_chromosomes(df, chrom_column='chrom')
    new_number_records = len(df)
    number_records_removed = old_number_records - new_number_records
    print(f'Removed {number_records_removed} records from {sys.argv[1]} with funky chromosomes')
    df.write_csv(
        file=sys.argv[2], 
        separator='\t', 
        include_header=False
    )

if __name__ == "__main__":
    main()