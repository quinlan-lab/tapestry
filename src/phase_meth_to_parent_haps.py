import argparse
import logging
from pathlib import Path
import bioframe as bf
import polars as pl

from phasing_trio import (
    get_phase_blocks,
    get_all_phasing,
)
from hap_map_trio import get_hap_map
from util.hap_map import write_hap_map_blocks
from get_meth_hap1_hap2 import read_meth_hap1_hap2
from util.shell import shell
from util.write_data import (
    write_bit_vector_mismatches_bed, 
    write_bit_vector_mismatches_vcf,
    write_bed, 
)

REFERENCE_GENOME = "hg38"


def _combine_count_model(
    df_count: pl.DataFrame,
    df_model: pl.DataFrame,
) -> pl.DataFrame:
    """Combine count-based and model-based methylation into one dataframe.

    The two dataframes share chrom, start, end, total_read_count_hap1/hap2.
    The methylation_level columns are suffixed with _count and _model.
    """
    df_count = df_count.rename({
        'methylation_level_hap1': 'methylation_level_hap1_count',
        'methylation_level_hap2': 'methylation_level_hap2_count',
    })
    df_model = df_model.select([
        'chrom', 'start', 'end',
        pl.col('methylation_level_hap1').alias('methylation_level_hap1_model'),
        pl.col('methylation_level_hap2').alias('methylation_level_hap2_model'),
    ])
    return df_count.join(df_model, on=['chrom', 'start', 'end'], how='full', coalesce=True)


def phase_meth_to_parent_haplotypes(
    df_meth_count_hap1_hap2_kid: pl.DataFrame,
    df_meth_model_hap1_hap2_kid: pl.DataFrame,
    df_meth_count_hap1_hap2_dad: pl.DataFrame,
    df_meth_model_hap1_hap2_dad: pl.DataFrame,
    df_meth_count_hap1_hap2_mom: pl.DataFrame,
    df_meth_model_hap1_hap2_mom: pl.DataFrame,
    df_hap_map_pat: pl.DataFrame,
    df_hap_map_mat: pl.DataFrame,
) -> pl.DataFrame:
    """
    Phase count-based and model-based methylation to parent haplotypes for a trio.

    Combines count and model methylation for each individual, overlaps kid's CpG
    methylation with the paternal and maternal hap maps, then joins dad's and
    mom's methylation at the same CpG positions.

    Haplotype labeling:
    - Kid: hap1=paternal (kid_pat), hap2=maternal (kid_mat)
    - Dad: hap1=A, hap2=B (fixed by definition)
    - Mom: hap1=C, hap2=D (fixed by definition)

    Columns are nulled out when the corresponding haplotype is unphased:
    - kid_pat, dad_A, dad_B columns are null when paternal_haplotype is null
    - kid_mat, mom_C, mom_D columns are null when maternal_haplotype is null
    """
    # Combine count and model for each individual
    df_kid = _combine_count_model(df_meth_count_hap1_hap2_kid, df_meth_model_hap1_hap2_kid)
    df_dad = _combine_count_model(df_meth_count_hap1_hap2_dad, df_meth_model_hap1_hap2_dad)
    df_mom = _combine_count_model(df_meth_count_hap1_hap2_mom, df_meth_model_hap1_hap2_mom)

    # Overlap kid meth with paternal hap map
    df = bf.overlap(
        df_kid.to_pandas(),
        df_hap_map_pat.to_pandas(),
        how='left', # we want to retain all CpG sites, including those that we cannot phase to haplotypes in the dad
        suffixes=('', '_'),
    )
    df = (
        pl.from_pandas(df)
        .drop(['chrom_'])
        .rename({
            'start_': 'start_hap_map_block_pat',
            'end_': 'end_hap_map_block_pat',
            'paternal_haplotype_': 'paternal_haplotype',
            'haplotype_concordance_': 'paternal_concordance',
            'num_het_SNVs_in_parent_': 'num_het_SNVs_in_dad',
        })
    )

    # Overlap with maternal hap map
    df = bf.overlap(
        df.to_pandas(),
        df_hap_map_mat.to_pandas(),
        how='left', # we want to retain all CpG sites, including those that we cannot phase to haplotypes in the mom
        suffixes=('', '_'),
    )
    df = (
        pl.from_pandas(df)
        .drop(['chrom_'])
        .rename({
            'start_': 'start_hap_map_block_mat',
            'end_': 'end_hap_map_block_mat',
            'maternal_haplotype_': 'maternal_haplotype',
            'haplotype_concordance_': 'maternal_concordance',
            'num_het_SNVs_in_parent_': 'num_het_SNVs_in_mom',
        })
    )

    # Kid: hap1=paternal, hap2=maternal; null when haplotype is unphased
    df = (
        df
        .with_columns(
            methylation_level_kid_pat_count=(
                pl.when(pl.col("paternal_haplotype").is_not_null())
                .then(pl.col("methylation_level_hap1_count"))
                .otherwise(None)
            ),
            methylation_level_kid_mat_count=(
                pl.when(pl.col("maternal_haplotype").is_not_null())
                .then(pl.col("methylation_level_hap2_count"))
                .otherwise(None)
            ),
            methylation_level_kid_pat_model=(
                pl.when(pl.col("paternal_haplotype").is_not_null())
                .then(pl.col("methylation_level_hap1_model"))
                .otherwise(None)
            ),
            methylation_level_kid_mat_model=(
                pl.when(pl.col("maternal_haplotype").is_not_null())
                .then(pl.col("methylation_level_hap2_model"))
                .otherwise(None)
            ),
            total_read_count_kid_pat=(
                pl.when(pl.col("paternal_haplotype").is_not_null())
                .then(pl.col("total_read_count_hap1"))
                .otherwise(None)
            ),
            total_read_count_kid_mat=(
                pl.when(pl.col("maternal_haplotype").is_not_null())
                .then(pl.col("total_read_count_hap2"))
                .otherwise(None)
            ),
        )
        .drop(['methylation_level_hap1_count', 'methylation_level_hap2_count',
               'methylation_level_hap1_model', 'methylation_level_hap2_model',
               'total_read_count_hap1', 'total_read_count_hap2'])
    )

    # Join with dad: hap1=A, hap2=B
    df_dad = df_dad.rename({
        'methylation_level_hap1_count': 'methylation_level_dad_A_count',
        'methylation_level_hap2_count': 'methylation_level_dad_B_count',
        'methylation_level_hap1_model': 'methylation_level_dad_A_model',
        'methylation_level_hap2_model': 'methylation_level_dad_B_model',
        'total_read_count_hap1': 'total_read_count_dad_A',
        'total_read_count_hap2': 'total_read_count_dad_B',
    })
    # Capture CpG sites present in dad's meth data but not in kid's
    df = df.join(df_dad, on=['chrom', 'start', 'end'], how='full', coalesce=True)

    # Join with mom: hap1=C, hap2=D
    df_mom = df_mom.rename({
        'methylation_level_hap1_count': 'methylation_level_mom_C_count',
        'methylation_level_hap2_count': 'methylation_level_mom_D_count',
        'methylation_level_hap1_model': 'methylation_level_mom_C_model',
        'methylation_level_hap2_model': 'methylation_level_mom_D_model',
        'total_read_count_hap1': 'total_read_count_mom_C',
        'total_read_count_hap2': 'total_read_count_mom_D',
    })
    # Capture CpG sites present in mom's meth data but not in kid's
    df = df.join(df_mom, on=['chrom', 'start', 'end'], how='full', coalesce=True)

    # Null out parent columns where the corresponding haplotype is unphased
    pat_not_null = pl.col("paternal_haplotype").is_not_null()
    mat_not_null = pl.col("maternal_haplotype").is_not_null()
    df = df.with_columns(
        # Dad columns: null when paternal_haplotype is null
        methylation_level_dad_A_count=(
            pl.when(pat_not_null)
            .then(pl.col("methylation_level_dad_A_count"))
            .otherwise(None)
        ),
        methylation_level_dad_B_count=(
            pl.when(pat_not_null)
            .then(pl.col("methylation_level_dad_B_count"))
            .otherwise(None)
        ),
        methylation_level_dad_A_model=(
            pl.when(pat_not_null)
            .then(pl.col("methylation_level_dad_A_model"))
            .otherwise(None)
        ),
        methylation_level_dad_B_model=(
            pl.when(pat_not_null)
            .then(pl.col("methylation_level_dad_B_model"))
            .otherwise(None)
        ),
        total_read_count_dad_A=(
            pl.when(pat_not_null)
            .then(pl.col("total_read_count_dad_A"))
            .otherwise(None)
        ),
        total_read_count_dad_B=(
            pl.when(pat_not_null)
            .then(pl.col("total_read_count_dad_B"))
            .otherwise(None)
        ),
        # Mom columns: null when maternal_haplotype is null
        methylation_level_mom_C_count=(
            pl.when(mat_not_null)
            .then(pl.col("methylation_level_mom_C_count"))
            .otherwise(None)
        ),
        methylation_level_mom_D_count=(
            pl.when(mat_not_null)
            .then(pl.col("methylation_level_mom_D_count"))
            .otherwise(None)
        ),
        methylation_level_mom_C_model=(
            pl.when(mat_not_null)
            .then(pl.col("methylation_level_mom_C_model"))
            .otherwise(None)
        ),
        methylation_level_mom_D_model=(
            pl.when(mat_not_null)
            .then(pl.col("methylation_level_mom_D_model"))
            .otherwise(None)
        ),
        total_read_count_mom_C=(
            pl.when(mat_not_null)
            .then(pl.col("total_read_count_mom_C"))
            .otherwise(None)
        ),
        total_read_count_mom_D=(
            pl.when(mat_not_null)
            .then(pl.col("total_read_count_mom_D"))
            .otherwise(None)
        ),
    )

    # Cast columns that bf.overlap converts to f64 back to Int64
    int_cols = (
        [c for c in df.columns if c.startswith('total_read_count_')]
        + ['start_hap_map_block_pat', 'end_hap_map_block_pat',
           'start_hap_map_block_mat', 'end_hap_map_block_mat',
           'num_het_SNVs_in_dad', 'num_het_SNVs_in_mom']
    )
    df = df.with_columns(pl.col(c).cast(pl.Int64) for c in int_cols)

    return df


def write_bigwig(df, uid, hap_label, person_hap, pb_cpg_tool_mode, output_dir, logger=None):
    """
    Write a bigwig file for a given haplotype and pb-cpg-tools pileup mode,
    recording only non-null methylation values.
    """
    meth_col = f"methylation_level_{person_hap}_{pb_cpg_tool_mode}"
    df_bed_graph = (
        df
        .filter(pl.col(meth_col).is_not_null())
        .select([
            pl.col("chrom"),
            pl.col("start"),
            pl.col("end"),
            pl.col(meth_col),
        ])
        .to_pandas()
    )

    file_path = Path(output_dir) / f"{uid}.dna-methylation.{hap_label}.{pb_cpg_tool_mode}.{REFERENCE_GENOME}.bw"

    # Remove any existing file or symlink (bedGraphToBigWig refuses to overwrite)
    if file_path.exists() or file_path.is_symlink():
        file_path.unlink()

    bf.to_bigwig(
        df_bed_graph,
        bf.fetch_chromsizes(db=REFERENCE_GENOME),
        outpath=file_path,
        # we assume that user has bedGraphToBigWig in their PATH: 
        path_to_binary="bedGraphToBigWig",  # type: ignore
    )

    if logger:
        logger.info(f"Wrote bigwig file for {person_hap} {pb_cpg_tool_mode}-based methylation levels, ASSUMING {REFERENCE_GENOME}, to: '{file_path}'")


def main():
    parser = argparse.ArgumentParser(
        description='Phase HiFi DNA methylation in a trio to parental haplotypes (A/B in dad and C/D in mom)'
    )
    parser.add_argument('--kid_id', required=True, help='Child sample ID')
    parser.add_argument('--dad_id', required=True, help='Father sample ID')
    parser.add_argument('--mom_id', required=True, help='Mother sample ID')
    parser.add_argument('--vcf_pedmec_phased', required=True, help='WhatsHap pedmec-phased trio VCF')
    parser.add_argument('--blocks_tsv_kid', required=True, help='WhatsHap blocks TSV for child')
    parser.add_argument('--blocks_tsv_dad', required=True, help='WhatsHap blocks TSV for father')
    parser.add_argument('--blocks_tsv_mom', required=True, help='WhatsHap blocks TSV for mother')

    # Methylation BEDs for kid
    parser.add_argument('--bed_meth_count_hap1_kid', required=True)
    parser.add_argument('--bed_meth_count_hap2_kid', required=True)
    parser.add_argument('--bed_meth_model_hap1_kid', required=True)
    parser.add_argument('--bed_meth_model_hap2_kid', required=True)

    # Methylation BEDs for dad
    parser.add_argument('--bed_meth_count_hap1_dad', required=True)
    parser.add_argument('--bed_meth_count_hap2_dad', required=True)
    parser.add_argument('--bed_meth_model_hap1_dad', required=True)
    parser.add_argument('--bed_meth_model_hap2_dad', required=True)

    # Methylation BEDs for mom
    parser.add_argument('--bed_meth_count_hap1_mom', required=True)
    parser.add_argument('--bed_meth_count_hap2_mom', required=True)
    parser.add_argument('--bed_meth_model_hap1_mom', required=True)
    parser.add_argument('--bed_meth_model_hap2_mom', required=True)

    parser.add_argument('--output_dir', required=True, help='Output directory')
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(filename)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    logger = logging.getLogger(__name__)

    logger.info(f"Starting '{__file__}'")
    logger.info("Args: %s", vars(args))

    Path(args.output_dir).mkdir(parents=True, exist_ok=True)

    kid_id, dad_id, mom_id = args.kid_id, args.dad_id, args.mom_id

    # Step 1: Get phase blocks for each individual (pre-computed by run-whatshap.sh)
    logger.info("Getting phase blocks for each individual...")
    df_blocks_kid = get_phase_blocks(args.blocks_tsv_kid, kid_id)
    df_blocks_dad = get_phase_blocks(args.blocks_tsv_dad, dad_id)
    df_blocks_mom = get_phase_blocks(args.blocks_tsv_mom, mom_id)
    logger.info(f"Phase blocks: kid={len(df_blocks_kid)}, dad={len(df_blocks_dad)}, mom={len(df_blocks_mom)}")

    # Step 2: Get paired parent-child alleles from the WhatsHap pedmec-phased VCF
    logger.info("Getting paired kid-parent allele sequences from the WhatsHap pedmec-phased VCF...")
    df_kid_dad, df_kid_mom = get_all_phasing(
        args.vcf_pedmec_phased, kid_id, dad_id, mom_id,
        df_blocks_kid, df_blocks_dad, df_blocks_mom,
    )
    logger.info(f"Number of het sites in dad: {len(df_kid_dad)}, Number of het sites in mom: {len(df_kid_mom)}")

    # Step 3: Construct a hap map, consisting of intervals ("hap-map blocks") in which
    # haplotypes in the kid are mapped to haplotypes in the parents
    logger.info("Constructing hap map (hap-map blocks mapping kid haplotypes to parent haplotypes)...")
    df_hap_map_pat, df_hap_map_mat, df_mismatch_pat, df_mismatch_mat = get_hap_map(df_kid_dad, df_kid_mom)
    logger.info(f"Paternal hap map: {len(df_hap_map_pat)} blocks, Maternal hap map: {len(df_hap_map_mat)} blocks")
    logger.info(f"Paternal bitvector mismatches: {len(df_mismatch_pat)}, Maternal bitvector mismatches: {len(df_mismatch_mat)}")

    # Write hap-map blocks for IGV (separately for paternal and maternal)
    write_hap_map_blocks(df_hap_map_pat, kid_id, "paternal", args.output_dir)
    write_hap_map_blocks(df_hap_map_mat, kid_id, "maternal", args.output_dir)
    logger.info(f"Wrote hap-map blocks to '{args.output_dir}'")

    # Write bitvector mismatch sites to BED (for later proximity computation) and VCF (for IGV visualization)
    write_bit_vector_mismatches_bed(args.output_dir, df_mismatch_pat, logger, uid=kid_id, parental="paternal")
    write_bit_vector_mismatches_bed(args.output_dir, df_mismatch_mat, logger, uid=kid_id, parental="maternal")
    write_bit_vector_mismatches_vcf(args.output_dir, df_mismatch_pat, logger, uid=kid_id, parental="paternal")
    write_bit_vector_mismatches_vcf(args.output_dir, df_mismatch_mat, logger, uid=kid_id, parental="maternal")

    # Step 4: Get HiFi DNA methylation levels (both count-based and model-based)
    # in kid, dad, and mom at CpG sites phased to hap1/hap2
    logger.info("Getting HiFi DNA methylation levels (count and model) for kid, dad, and mom...")
    df_meth_count_kid = read_meth_hap1_hap2('count', args.bed_meth_count_hap1_kid, args.bed_meth_count_hap2_kid)
    df_meth_model_kid = read_meth_hap1_hap2('model', args.bed_meth_model_hap1_kid, args.bed_meth_model_hap2_kid)
    logger.info(f"Kid methylation: {len(df_meth_count_kid)} count rows, {len(df_meth_model_kid)} model rows")

    df_meth_count_dad = read_meth_hap1_hap2('count', args.bed_meth_count_hap1_dad, args.bed_meth_count_hap2_dad)
    df_meth_model_dad = read_meth_hap1_hap2('model', args.bed_meth_model_hap1_dad, args.bed_meth_model_hap2_dad)
    logger.info(f"Dad methylation: {len(df_meth_count_dad)} count rows, {len(df_meth_model_dad)} model rows")

    df_meth_count_mom = read_meth_hap1_hap2('count', args.bed_meth_count_hap1_mom, args.bed_meth_count_hap2_mom)
    df_meth_model_mom = read_meth_hap1_hap2('model', args.bed_meth_model_hap1_mom, args.bed_meth_model_hap2_mom)
    logger.info(f"Mom methylation: {len(df_meth_count_mom)} count rows, {len(df_meth_model_mom)} model rows")

    # Step 5: Phase DNA methylation levels in trio to parent haplotypes
    logger.info("Phasing DNA methylation levels in trio to parent haplotypes...")
    df_meth_phased = phase_meth_to_parent_haplotypes(
        df_meth_count_kid, df_meth_model_kid,
        df_meth_count_dad, df_meth_model_dad,
        df_meth_count_mom, df_meth_model_mom,
        df_hap_map_pat, df_hap_map_mat,
    )
    logger.info(f"Phased methylation: {len(df_meth_phased)} rows")

    # Step 6: Write phased methylation BED, then sort, compress and index
    bed_stem = f"trio.dna-methylation"
    write_bed(args.output_dir, df_meth_phased, filename_stem=bed_stem)
    logger.info(f"Wrote '{args.output_dir}/{bed_stem}.bed'")
    cmd = (
        f'cat {args.output_dir}/{bed_stem}.bed'
        f' | src/util/sort-compress-index-bed'
        f' --name {args.output_dir}/{bed_stem}'
    )
    shell(cmd)
    logger.info(f"Wrote '{args.output_dir}/{bed_stem}.sorted.bed.gz'")

    # Step 7: Write bigwig files with only non-null methylation values
    logger.info("Writing bigwig files for phased methylation...")
    bigwig_specs = [
        # (uid, hap_label for filename, person_hap for column name)
        # Kid: hap1=paternal, hap2=maternal
        (kid_id, 'pat', 'kid_pat'),
        (kid_id, 'mat', 'kid_mat'),
        # Dad: hap1=A, hap2=B
        (dad_id, 'A', 'dad_A'),
        (dad_id, 'B', 'dad_B'),
        # Mom: hap1=C, hap2=D
        (mom_id, 'C', 'mom_C'),
        (mom_id, 'D', 'mom_D'),
    ]
    for uid, hap_label, person_hap in bigwig_specs:
        for pb_cpg_tool_mode in ['count', 'model']:
            write_bigwig(df_meth_phased, uid, hap_label, person_hap, pb_cpg_tool_mode, args.output_dir, logger)

    logger.info(f"Done running '{__file__}'")


if __name__ == "__main__":
    main()
