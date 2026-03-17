import argparse
import logging
from pathlib import Path
import bioframe as bf
import polars as pl

from phasing_trio import (
    get_pedmec_phasing,
    get_phase_blocks,
)
from hap_map_trio import (
    get_hap_map,
    write_hap_map_blocks,
    get_hap_map_blocks,
)
from get_meth_hap1_hap2 import read_meth_hap1_hap2
from phase_meth_to_founder_haps import (
    write_bigwig,
    add_suffix,
)
from util.write_data import write_bed
from util.version_sort import version_sort

REFERENCE_GENOME = "hg38"


def read_ped(ped_path: str, child_uid: str):
    """Parse a PED file and return (father_id, mother_id) for the given child."""
    with open(ped_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if parts[1] == child_uid:
                return parts[2], parts[3]
    raise ValueError(f"Sample '{child_uid}' not found in PED file '{ped_path}'")


def phase_meth_to_haplotypes(
    df_meth_hap1_hap2: pl.DataFrame,
    df_hap_map: pl.DataFrame,
    individual: str,
) -> pl.DataFrame:
    """
    Phase methylation to A/B/C/D haplotypes for a given individual.

    For the kid: hap1=paternal (A or B), hap2=maternal (C or D).
    For the dad: hap1=A, hap2=B (by definition).
    For the mom: hap1=C, hap2=D (by definition).

    The hap_map tells us, for each block, which of A/B is the kid's paternal
    haplotype and which of C/D is the kid's maternal haplotype. For the parents
    themselves, A/B and C/D are fixed as hap1/hap2 respectively.

    Args:
        individual: one of "kid", "dad", "mom"
    """
    df = bf.overlap(
        df_meth_hap1_hap2.to_pandas(),
        df_hap_map.to_pandas(),
        how='left',
        suffixes=('', '_'),
    )

    df = (
        pl.from_pandas(df)
        .drop(['chrom_'])
        .rename({
            'start_': 'start_hap_map_block',
            'end_': 'end_hap_map_block',
            'paternal_haplotype_': 'paternal_haplotype',
            'maternal_haplotype_': 'maternal_haplotype',
            'paternal_concordance_': 'paternal_concordance',
            'maternal_concordance_': 'maternal_concordance',
            'num_het_SNVs_pat_': 'num_het_SNVs_pat',
            'num_het_SNVs_mat_': 'num_het_SNVs_mat',
        })
    )

    if individual == "kid":
        # Kid: hap1=paternal (from dad), hap2=maternal (from mom)
        # paternal_haplotype (A or B) and maternal_haplotype (C or D) indicate
        # which founder haplotype was transmitted per block.
        df = df.with_columns(
            pl.col('methylation_level_hap1').alias('methylation_level_pat'),
            pl.col('methylation_level_hap2').alias('methylation_level_mat'),
            pl.col('total_read_count_hap1').alias('total_read_count_pat'),
            pl.col('total_read_count_hap2').alias('total_read_count_mat'),
        )
    elif individual == "dad":
        # Dad: hap1=A, hap2=B (fixed by definition)
        df = df.with_columns(
            pl.col('methylation_level_hap1').alias('methylation_level_A'),
            pl.col('methylation_level_hap2').alias('methylation_level_B'),
            pl.col('total_read_count_hap1').alias('total_read_count_A'),
            pl.col('total_read_count_hap2').alias('total_read_count_B'),
        )
    elif individual == "mom":
        # Mom: hap1=C, hap2=D (fixed by definition)
        df = df.with_columns(
            pl.col('methylation_level_hap1').alias('methylation_level_C'),
            pl.col('methylation_level_hap2').alias('methylation_level_D'),
            pl.col('total_read_count_hap1').alias('total_read_count_C'),
            pl.col('total_read_count_hap2').alias('total_read_count_D'),
        )

    cols_to_drop = [
        'methylation_level_hap1', 'methylation_level_hap2',
        'total_read_count_hap1', 'total_read_count_hap2',
        'paternal_concordance', 'maternal_concordance',
        'num_het_SNVs_pat', 'num_het_SNVs_mat',
    ]
    if individual != "kid":
        # Parents don't need the hap map indicator columns
        cols_to_drop += ['paternal_haplotype', 'maternal_haplotype']
    df = df.drop(cols_to_drop)

    return df


METH_LEVEL_SUFFIXES = {
    "kid": ["pat", "mat"],
    "dad": ["A", "B"],
    "mom": ["C", "D"],
}


def combine_count_and_model_based_methylation_levels(
    df_meth_count: pl.DataFrame,
    df_meth_model: pl.DataFrame,
    individual: str,
) -> pl.DataFrame:
    suffix_count, suffix_model = '_count', '_model'
    haps = METH_LEVEL_SUFFIXES[individual]
    unique_cols = [f'methylation_level_{h}' for h in haps]

    df_meth_count = add_suffix(df_meth_count, suffix=suffix_count, cols_to_suffix=unique_cols)
    df_meth_model = add_suffix(df_meth_model, suffix=suffix_model, cols_to_suffix=unique_cols)

    unique_cols_count = [f'{col}{suffix_count}' for col in unique_cols]
    common_cols = [col for col in df_meth_count.columns if col not in unique_cols_count]

    df = df_meth_count.join(
        df_meth_model,
        on=common_cols,
        nulls_equal=True,
        how="full",
        coalesce=True,
    )

    unique_cols_model = [f'{col}{suffix_model}' for col in unique_cols]
    df = df.select(common_cols + unique_cols_count + unique_cols_model)

    return version_sort(df)


def process_individual(
    uid: str,
    individual: str,
    bed_meth_count_hap1: str,
    bed_meth_count_hap2: str,
    bed_meth_model_hap1: str,
    bed_meth_model_hap2: str,
    df_hap_map: pl.DataFrame,
    output_dir: str,
    logger,
):
    """Phase methylation to A/B/C/D for one individual and write outputs."""
    df_meth_count_hap1_hap2 = read_meth_hap1_hap2(
        pb_cpg_tool_mode='count',
        bed_hap1=bed_meth_count_hap1,
        bed_hap2=bed_meth_count_hap2,
    )
    logger.info(f"[{uid}] Got count-based methylation: {len(df_meth_count_hap1_hap2)} rows")

    df_meth_model_hap1_hap2 = read_meth_hap1_hap2(
        pb_cpg_tool_mode='model',
        bed_hap1=bed_meth_model_hap1,
        bed_hap2=bed_meth_model_hap2,
    )
    logger.info(f"[{uid}] Got model-based methylation: {len(df_meth_model_hap1_hap2)} rows")

    df_meth_count_phased = phase_meth_to_haplotypes(df_meth_count_hap1_hap2, df_hap_map, individual)
    logger.info(f"[{uid}] Phased count-based methylation: {len(df_meth_count_phased)} rows")

    df_meth_model_phased = phase_meth_to_haplotypes(df_meth_model_hap1_hap2, df_hap_map, individual)
    logger.info(f"[{uid}] Phased model-based methylation: {len(df_meth_model_phased)} rows")

    df_meth_phased = combine_count_and_model_based_methylation_levels(
        df_meth_count_phased, df_meth_model_phased, individual
    )
    logger.info(f"[{uid}] Combined count+model methylation: {len(df_meth_phased)} rows")

    write_bed(output_dir, df_meth_phased, filename_stem=f"{uid}.dna-methylation.parent-phased")
    logger.info(f"[{uid}] Wrote phased methylation to '{output_dir}'")

    haps = METH_LEVEL_SUFFIXES[individual]
    for hap in haps:
        for pb_cpg_tool_mode in ['count', 'model']:
            write_bigwig(df_meth_phased, uid, hap, pb_cpg_tool_mode, output_dir, logger)


def main():
    parser = argparse.ArgumentParser(
        description='Phase HiFi DNA methylation to parental founder haplotypes (A,B,C,D) for a trio'
    )
    parser.add_argument('--kid_id', required=True, help='Child sample ID')
    parser.add_argument('--dad_id', required=True, help='Father sample ID')
    parser.add_argument('--mom_id', required=True, help='Mother sample ID')
    parser.add_argument('--ped', required=True, help='Pedigree file (PED format)')
    parser.add_argument('--vcf_trio_phased', required=True, help='WhatsHap trio-phased VCF (multi-sample)')
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

    # Step 1: Read phasing data from the multi-sample WhatsHap VCF
    logger.info("Reading trio phasing data from VCF...")
    phasing = get_pedmec_phasing(args.vcf_trio_phased, kid_id, dad_id, mom_id)
    df_kid = phasing[kid_id]
    df_dad = phasing[dad_id]
    df_mom = phasing[mom_id]
    logger.info(f"Het SNVs: kid={len(df_kid)}, dad={len(df_dad)}, mom={len(df_mom)}")

    # Step 2: Get phase blocks for each individual (pre-computed by run-whatshap.sh)
    logger.info("Reading phase blocks for each individual...")
    df_blocks_kid = get_phase_blocks(args.blocks_tsv_kid, kid_id)
    df_blocks_dad = get_phase_blocks(args.blocks_tsv_dad, dad_id)
    df_blocks_mom = get_phase_blocks(args.blocks_tsv_mom, mom_id)
    logger.info(f"Phase blocks: kid={len(df_blocks_kid)}, dad={len(df_blocks_dad)}, mom={len(df_blocks_mom)}")

    # Step 3: Compute hap-map blocks (intersection of all three)
    logger.info("Computing hap-map blocks (intersection of phase blocks)...")
    df_hap_map_blocks = get_hap_map_blocks(df_blocks_kid, df_blocks_dad, df_blocks_mom)
    logger.info(f"Hap-map blocks: {len(df_hap_map_blocks)}")

    # Step 4: Build the haplotype map via bit-vector comparison
    logger.info("Building trio haplotype map...")
    df_hap_map = get_hap_map(
        df_kid, df_dad, df_mom,
        df_blocks_kid, df_blocks_dad, df_blocks_mom
    )
    logger.info(f"Hap map: {len(df_hap_map)} blocks")

    # Write hap-map blocks for IGV
    write_hap_map_blocks(df_hap_map, kid_id, "paternal", args.output_dir)
    write_hap_map_blocks(df_hap_map, kid_id, "maternal", args.output_dir)
    write_bed(args.output_dir, df_hap_map, f"{kid_id}.hap-map-blocks")
    logger.info(f"Wrote hap-map blocks to '{args.output_dir}'")

    # Step 5: Phase methylation to A, B, C, D for all three individuals
    individuals = {
        "kid": (kid_id, args.bed_meth_count_hap1_kid, args.bed_meth_count_hap2_kid,
                args.bed_meth_model_hap1_kid, args.bed_meth_model_hap2_kid),
        "dad": (dad_id, args.bed_meth_count_hap1_dad, args.bed_meth_count_hap2_dad,
                args.bed_meth_model_hap1_dad, args.bed_meth_model_hap2_dad),
        "mom": (mom_id, args.bed_meth_count_hap1_mom, args.bed_meth_count_hap2_mom,
                args.bed_meth_model_hap1_mom, args.bed_meth_model_hap2_mom),
    }

    for individual, (uid, bc1, bc2, bm1, bm2) in individuals.items():
        logger.info(f"Processing {individual} ({uid})...")
        process_individual(
            uid=uid,
            individual=individual,
            bed_meth_count_hap1=bc1,
            bed_meth_count_hap2=bc2,
            bed_meth_model_hap1=bm1,
            bed_meth_model_hap2=bm2,
            df_hap_map=df_hap_map,
            output_dir=args.output_dir,
            logger=logger,
        )

    logger.info(f"Done running '{__file__}'")


if __name__ == "__main__":
    main()
