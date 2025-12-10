import argparse
import logging
import pysam
import re

from remove_funky_chromosomes import funky

def scan_and_write_cpgs(reference, output_file, logger):
    """
    Scans a reference genome for CpG sites, filters out 'funky' chromosomes,
    and writes valid sites directly to a BED file record-by-record.
    
    Args:
        reference (str): Path to the reference FASTA.
        output_file (str): Path to the output BED/TSV file.
        logger: Logger instance.
    """
    count = 0
    
    # Open both the FASTA and the output file simultaneously
    with pysam.FastaFile(reference) as fa, open(output_file, "w") as out_f:
        
        for chrom in fa.references:
            # 1. Filter Logic:
            # Check the chromosome name *before* fetching the sequence. 
            if funky(chrom):
                logger.debug(f"Skipping funky chromosome: {chrom}")
                continue

            logger.info(f"Scanning '{chrom}'")
            
            # 2. Fetch Sequence:
            sequence = fa.fetch(chrom).upper()

            # 3. Stream & Write:
            # re.finditer yields matches one by one (generator), avoiding big lists.
            for match in re.finditer('CG', sequence):
                start = match.start()
                # Write immediately to file: chrom, start, end (start+1)
                out_f.write(f"{chrom}\t{start}\t{start + 1}\n")
                count += 1

    logger.info(f"Finished. Wrote {count} CpG sites to {output_file}")

@profile # type:ignore 
def main(): 
    parser = argparse.ArgumentParser(description='Get all CpG sites in reference genome and write them to disk')
    parser.add_argument('--reference', required=True, help='Genome reference sequence to find locations of CpG sites')
    parser.add_argument('--bed_all_cpgs_in_reference', required=True, help='Bed file to store CpG sites observed in reference genome')
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG,  # <--- This allows debug messages through
        format='%(asctime)s - %(levelname)s - %(filename)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    logger = logging.getLogger(__name__)

    logger.info(f"Starting '{__file__}'")
    logger.info("Script started with the following arguments: %s", vars(args))

    scan_and_write_cpgs(
        reference=args.reference, 
        output_file=args.bed_all_cpgs_in_reference, 
        logger=logger
    )
    logger.info(f"Done running '{__file__}'")

if __name__ == "__main__":
    main()
