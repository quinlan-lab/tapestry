from cyvcf2 import VCF  # type: ignore
import polars as pl
import bioframe as bf


def is_snv_het(variant, sample_index):
    gt_types = variant.gt_types
    is_het = gt_types[sample_index] == 1
    is_snp = variant.is_snp
    has_single_ALT_allele = len(variant.ALT) == 1
    return is_het and is_snp and has_single_ALT_allele


def get_trio_read_phasing(vcf_path: str, kid_id: str, dad_id: str, mom_id: str):
    """
    Read the multi-sample WhatsHap trio-phased VCF and extract per-sample
    het SNV alleles and phase block IDs for all three individuals.

    Returns one DataFrame per individual with columns:
        chrom, start, end, REF, ALT, phase_block_id, allele_hap1, allele_hap2
    """
    dfs = {uid: [] for uid in [kid_id, dad_id, mom_id]}

    # Don't use context manager — cyvcf2 supports it on Linux but not macOS
    vcf_reader = VCF(vcf_path, strict_gt=True)
    samples = vcf_reader.samples
    sample_indices = {uid: samples.index(uid) for uid in [kid_id, dad_id, mom_id]}

    for variant in vcf_reader:
        chrom = variant.CHROM
        pos = variant.POS
        start = pos - 1
        end = pos
        REF = variant.REF
        ALT = variant.ALT[0] if len(variant.ALT) == 1 else None

        if ALT is None:
            continue
        if not variant.is_snp:
            continue

        genotype_all = variant.genotypes
        phase_block_all = variant.format('PS')

        for uid in [kid_id, dad_id, mom_id]:
            si = sample_indices[uid]

            if not is_snv_het(variant, si):
                continue

            genotype = genotype_all[si]
            phased = genotype[2]
            if not phased:
                continue

            allele_hap1 = str(genotype[0]) if genotype[0] != -1 else '.'
            allele_hap2 = str(genotype[1]) if genotype[1] != -1 else '.'

            if phase_block_all is None:
                phase_block_id = '.'
            else:
                pb = phase_block_all[si, 0]
                phase_block_id = '.' if pb == -2147483648 else str(pb)

            dfs[uid].append({
                "chrom": chrom,
                "start": start,
                "end": end,
                "REF": REF,
                "ALT": ALT,
                "phase_block_id": phase_block_id,
                "allele_hap1": allele_hap1,
                "allele_hap2": allele_hap2,
            })

    vcf_reader.close()
    return {uid: pl.DataFrame(records) for uid, records in dfs.items()}


def get_trio_phase_blocks(blocks_tsv: str, uid: str) -> pl.DataFrame:
    """
    Read the whatshap stats --block-list TSV (produced by run-whatshap.sh) for a sample.
    Returns DataFrame with columns: chrom, start, end, phase_block_id, num_variants
    """
    df = (
        pl.read_csv(blocks_tsv, separator='\t', has_header=True)
        .rename({
            '#sample': 'sample',
            'chromosome': 'chrom',
            'from': 'start',
            'to': 'end',
            'variants': 'num_variants',
        })
        .filter(pl.col('sample') == uid)
        .with_columns([
            (pl.col('start') - 1).alias('start'),  # 1-based to 0-based
            pl.col('phase_set').cast(pl.String).alias('phase_block_id'),
        ])
        .select(['chrom', 'start', 'end', 'phase_block_id', 'num_variants'])
    )
    return df


def get_trio_hap_map_blocks(
    df_blocks_kid: pl.DataFrame,
    df_blocks_dad: pl.DataFrame,
    df_blocks_mom: pl.DataFrame,
) -> pl.DataFrame:
    """
    Compute the intersection of phase blocks across all three individuals.
    Each resulting interval is a hap-map block.
    """
    # Intersect kid and dad blocks
    df_kd = pl.from_pandas(bf.overlap(
        df_blocks_kid.select(["chrom", "start", "end"]).to_pandas(),
        df_blocks_dad.select(["chrom", "start", "end"]).to_pandas(),
        how='inner',
        suffixes=('_kid', '_dad'),
    ))

    # Compute the intersection interval
    df_kd = df_kd.with_columns(
        pl.max_horizontal("start_kid", "start_dad").alias("start"),
        pl.min_horizontal("end_kid", "end_dad").alias("end"),
    ).filter(pl.col("start") < pl.col("end"))

    # Now intersect with mom blocks
    df_kdm = pl.from_pandas(bf.overlap(
        df_kd.select(["chrom_kid", "start", "end"]).rename({"chrom_kid": "chrom"}).to_pandas(),
        df_blocks_mom.select(["chrom", "start", "end"]).to_pandas(),
        how='inner',
        suffixes=('_kd', '_mom'),
    ))

    df_hap_map_blocks = (
        df_kdm
        .with_columns(
            pl.max_horizontal("start_kd", "start_mom").alias("start"),
            pl.min_horizontal("end_kd", "end_mom").alias("end"),
        )
        .filter(pl.col("start") < pl.col("end"))
        .select([
            pl.col("chrom_kd").alias("chrom"),
            "start",
            "end",
        ])
        .sort(["chrom", "start", "end"])
    )

    return df_hap_map_blocks
