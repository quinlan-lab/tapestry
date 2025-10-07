import argparse
import logging
import pysam
import re
import os
import polars as pl 

from remove_funky_chromosomes import remove_funky_chromosomes

def find_all_cpgs(reference, logger):
    """
    Scans a reference genome FASTA file for all CpG sites and returns them
    as a Polars DataFrame.

    Args:
        reference (str): The path to the reference FASTA file.

    Returns:
        pl.DataFrame: A DataFrame with columns ['chrom', 'start', 'end']
                      representing the CpG sites in BED format.
    """
    cpg_sites = []
    with pysam.FastaFile(reference) as fa:
        for chrom in fa.references:
            logger.info(f"Scanning {chrom} from: {os.path.abspath(reference)}")
            sequence = fa.fetch(chrom).upper()
            # Use a list comprehension for efficiency
            cpg_sites.extend(
                [
                    (chrom, match.start(), match.start() + 1)
                    for match in re.finditer('CG', sequence)
                ]
            )

    logger.info(f"Found {len(cpg_sites)} CpG sites. Creating DataFrame.")

    # Create a Polars DataFrame from the list of tuples
    df = pl.DataFrame(
        cpg_sites,
        schema={'chrom': pl.Utf8, 'start': pl.Int64, 'end': pl.Int64}
    )

    return df

def main(args): 
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    logger = logging.getLogger(__name__)

    df = find_all_cpgs(args.reference, logger)
    df = remove_funky_chromosomes(df, chrom_column='chrom')
    df.write_csv(args.bed_all_cpgs, separator='\t', include_header=False)
    logger.info(f"Wrote {args.bed_all_cpgs}.")
    logger.info(f"Done running {__file__}.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Get all CpG sites and write them to disk')
    parser.add_argument('--reference', required=True, help='Genome reference sequence to find locations of CpG sites')
    parser.add_argument('--bed_all_cpgs', required=True, help='Bed file to store CpG sites in')
    args = parser.parse_args()
    main(args)
