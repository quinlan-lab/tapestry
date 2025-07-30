import gzip
from pathlib import Path
import polars as pl
from util.remove_funky_chromosomes import remove_funky_chromosomes

def read_meth_level(bed: Path, pb_cpg_tool_mode: str) -> pl.DataFrame:
    """
    Reads a methylation level bed file from pb-cpg-tools into a Polars DataFrame.
    
    The function asserts that the file contains header lines starting with '##' 
    and a single header line starting with '#'. It also asserts that the
    pileup-mode in the file matches the expected mode.
    """
    is_gzipped = str(bed).endswith('.gz')
    _open = gzip.open if is_gzipped else open

    has_comment_lines = False
    has_header_line = False
    found_pileup_mode = False

    with _open(bed, 'rt') as f:
        for line in f:
            if line.startswith('##pileup-mode='):
                pileup_mode_in_file = line.strip().split('=')[1]
                assert pileup_mode_in_file == pb_cpg_tool_mode, \
                    f"Expected pileup-mode '{pb_cpg_tool_mode}' but found '{pileup_mode_in_file}' in {bed}"
                found_pileup_mode = True
            
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
    assert found_pileup_mode, f"File {bed} is missing '##pileup-mode' comment line."

    df = (
        pl
        .read_csv(
            bed,
            separator='\t',
            comment_prefix='##',
            has_header=True,
            # The header line starts with '#', which polars handles automatically.
        )
        # "The bed file columns will differ between the model and count pileup methods, but both share the first six columns"
        # https://github.com/PacificBiosciences/pb-CpG-tools?tab=readme-ov-file#bed-file-format
        .rename({
            '#chrom': 'chromosome',
            'begin': 'start',
            'end': 'end',
            'mod_score': 'methylation_level_percent',
            'type': 'type',
            'cov': 'total_read_count',
        })
        .select([
            'chromosome', 
            'start', 
            'end', 
            'methylation_level_percent', 
            'total_read_count'
        ]) 
        .with_columns(
            (pl.col('methylation_level_percent') / 100)
            .alias('methylation_level')
        )
        .drop('methylation_level_percent')
    )

    df = remove_funky_chromosomes(df, chrom_column='chromosome')

    return df

def get_meth_hap1_hap2(pb_cpg_tool_mode, bed_hap1, bed_hap2): 
    df_meth_hap1 = read_meth_level(
        bed_hap1,
        pb_cpg_tool_mode
    ).select(pl.all().name.suffix("_hap1"))

    df_meth_hap2 = read_meth_level(
        bed_hap2,
        pb_cpg_tool_mode
    ).select(pl.all().name.suffix("_hap2"))

    df_meth = (
        df_meth_hap1        
        .join(
            df_meth_hap2,
            left_on=["chromosome_hap1", "start_hap1", "end_hap1"],
            right_on=["chromosome_hap2", "start_hap2", "end_hap2"],
            how="inner", # consider only CpG sites for which methylation levels are present in both haplotypes
        )
        .rename({
            "chromosome_hap1": "chrom", 
            "start_hap1": "start", 
            "end_hap1": "end",
        })
    )
    return df_meth
