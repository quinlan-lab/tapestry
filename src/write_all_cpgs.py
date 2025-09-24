import argparse
import logging
from pathlib import Path
import pysam
import re
import os

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)

def write_all_cpgs(reference, output_dir):
    bed_filename = f"{output_dir}/all_cpg_sites.bed"
    with (
        pysam.FastaFile(reference) as fa,
        open(bed_filename, 'w') as bed
    ):
        for chrom in fa.references:
            logger.info(f"Scanning {chrom} from: {os.path.abspath(reference)}")
            sequence = fa.fetch(chrom).upper()
            for match in re.finditer('CG', sequence):
                start = match.start()
                end = start + 1 # The C in the CpG is at pos start
                bed.write(f"{chrom}\t{start}\t{end}\n")
    logger.info(f"CpG sites saved to: {os.path.abspath(bed_filename)}")

def main(args): 
    Path(args.output_dir).mkdir(parents=True, exist_ok=True)
    write_all_cpgs(args.reference, args.output_dir)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Get all CpG sites and write them to disk')
    parser.add_argument('--reference', required=True, help='Genome reference sequence to find locations of CpG sites')
    parser.add_argument('--output_dir', required=True, help='Directory to store CpG sites in')
    args = parser.parse_args()
    main(args)
    logger.info(f"Done running {__file__}")