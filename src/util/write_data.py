import os
from pathlib import Path
import polars as pl

from shell import shell

# https://samtools.github.io/hts-specs/BEDv1.pdf
# "We recommend that only a single tab (\t) be used as field separator."
# "Comment lines start with # with no horizontal whitespace beforehand. A # appearing anywhere else in a data line is treated as feature data, not a comment."
def write_dataframe_to_bed(df: pl.DataFrame, file_path: str, source: str):
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
        header_string = f"##source='{source}'\n"
        header_string += "#" + "\t".join(df.columns) + "\n"
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


def write_bit_vector_mismatches_bed(output_dir, df_sites_mismatch, logger, uid=None, parental=None):
    """Write bitvector mismatch sites to BED for later proximity computation.

    Used by both pedigree-wise and trio-wise workflows.
    """
    prefix = f"{uid}." if uid else ""
    suffix = f".{parental}" if parental else ""
    bed = f"{output_dir}/{prefix}bit-vector-sites-mismatches{suffix}.bed"
    write_bed_and_header(bed, df_sites_mismatch)
    logger.info(f"Wrote sites where bit vectors are mismatched to: '{bed}'")
    logger.info("These will be used for later computation of proximity of (all) cpg sites to mismatch sites")


def write_bit_vector_mismatches_vcf(output_dir, df_sites_mismatch, logger, uid=None, parental=None):
    """Write bitvector mismatch sites to VCF for IGV visualization, then compress and index.

    Used by both pedigree-wise and trio-wise workflows.
    """
    prefix = f"{uid}." if uid else ""
    suffix = f".{parental}" if parental else ""
    stem = f"{output_dir}/{prefix}bit-vector-sites-mismatches{suffix}"
    vcf = f"{stem}.vcf"
    write_df_to_vcf(df_sites_mismatch, vcf, uid=uid)

    cmd = f'src/util/compress-index-vcf --name {stem}'
    shell(cmd)

    logger.info(f"Wrote bit-vector-sites-mismatches (for IGV) to: '{vcf}.gz'")


def write_df_to_vcf(df, vcf, uid=None):
    df = df.sort(["chrom", "start", "end"])

    sample = uid if uid else "SAMPLE"
    header = [
        '##fileformat=VCFv4.2',
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + sample
    ]

    with open(vcf, 'w') as f:
        for line in header:
            f.write(line + '\n')

        for row in df.iter_rows(named=True):
            chrom = row['chrom']
            pos = row['start'] + 1  # VCF is 1-based
            id_ = '.'
            ref = row['REF']
            alt = row['ALT']
            qual = '.'
            flt = '.'
            info = '.'
            fmt_keys = '.'
            fmt_vals = '.'
            f.write(f"{chrom}\t{pos}\t{id_}\t{ref}\t{alt}\t{qual}\t{flt}\t{info}\t{fmt_keys}\t{fmt_vals}\n")


def write_methylation(df, file_path, source):
    """Write a methylation dataframe to a sorted, compressed, and indexed BED file.

    Writes the dataframe to a BED file, then sorts, compresses (bgzip), and indexes (tabix).
    The original uncompressed file is removed.

    Returns the root path (without .sorted.bed.gz suffix).
    """
    write_dataframe_to_bed(df, file_path, source)
    root, suffix = os.path.splitext(file_path)

    cmd = (
        f'cat {file_path}'
        f' | src/util/sort-compress-index-bed'
        f' --name {root}'
    )
    shell(cmd)
    shell(f'rm {file_path}')

    return root



