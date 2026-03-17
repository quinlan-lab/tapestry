from cyvcf2 import VCF  # type: ignore
import bioframe as bf
import polars as pl

from phasing import is_snv_het

def get_phase_blocks(blocks_tsv: str, uid: str) -> pl.DataFrame:
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


def _annotate_phase_blocks(df_snvs, df_blocks, label):
    """Annotate SNVs with the phase block they fall in for a given individual."""
    df = pl.from_pandas(bf.overlap(
        df_snvs.to_pandas(),
        df_blocks.select(["chrom", "start", "end"]).to_pandas(),
        how='left',
        suffixes=('', f'_{label}'),
    ))
    return (
        df
        .rename({
            f"start_{label}": f"start_phase_block_{label}",
            f"end_{label}": f"end_phase_block_{label}",
        })
        .drop(f"chrom_{label}")
    )

def get_all_phasing(
    vcf_path: str,
    kid_id: str, dad_id: str, mom_id: str,
    df_blocks_kid: pl.DataFrame,
    df_blocks_dad: pl.DataFrame,
    df_blocks_mom: pl.DataFrame,
):
    """
    Read the multi-sample WhatsHap trio-phased VCF and construct paired
    parent-child allele DataFrames at sites heterozygous in each parent.

    At dad-het SNV sites (where dad is phased), report the kid's paternal
    allele (hap1) alongside dad's hap1 and hap2 alleles. The kid must be
    phased or homozygous to determine hap1.

    At mom-het SNV sites (where mom is phased), report the kid's maternal
    allele (hap2) alongside mom's hap1 and hap2 alleles. The kid must be
    phased or homozygous to determine hap2.

    Naming convention (matching get_hap_map):
        kid: pat = hap1, mat = hap2
        dad: A = hap1, B = hap2
        mom: C = hap1, D = hap2

    Returns:
        df_kid_dad: DataFrame at dad-het sites with columns:
            chrom, start, end, REF, ALT,
            kid_allele_pat, dad_allele_A, dad_allele_B,
            start_phase_block_kid, end_phase_block_kid,
            start_phase_block_dad, end_phase_block_dad
        df_kid_mom: DataFrame at mom-het sites with columns:
            chrom, start, end, REF, ALT,
            kid_allele_mat, mom_allele_C, mom_allele_D,
            start_phase_block_kid, end_phase_block_kid,
            start_phase_block_mom, end_phase_block_mom
    """
    records_kid_dad = []
    records_kid_mom = []

    # Don't use context manager — cyvcf2 supports it on Linux but not macOS
    vcf_reader = VCF(vcf_path, strict_gt=True)
    samples = vcf_reader.samples
    si_kid = samples.index(kid_id)
    si_dad = samples.index(dad_id)
    si_mom = samples.index(mom_id)

    for variant in vcf_reader:
        if not variant.is_snp:
            continue
        if len(variant.ALT) != 1:
            continue

        chrom = variant.CHROM
        pos = variant.POS
        start = pos - 1
        end = pos
        REF = variant.REF
        ALT = variant.ALT[0]

        genotype_all = variant.genotypes
        gt_kid = genotype_all[si_kid]
        kid_phased = gt_kid[2]
        kid_is_hom = (gt_kid[0] == gt_kid[1])

        # Dad-het sites → paternal allele in kid
        if is_snv_het(variant, si_dad):
            gt_dad = genotype_all[si_dad]
            if gt_dad[2]:  # dad must be phased
                if kid_phased or kid_is_hom:
                    records_kid_dad.append({
                        "chrom": chrom,
                        "start": start,
                        "end": end,
                        "REF": REF,
                        "ALT": ALT,
                        "kid_allele_pat": str(gt_kid[0]) if gt_kid[0] != -1 else '.',
                        "dad_allele_A": str(gt_dad[0]),
                        "dad_allele_B": str(gt_dad[1]),
                    })

        # Mom-het sites → maternal allele in kid
        if is_snv_het(variant, si_mom):
            gt_mom = genotype_all[si_mom]
            if gt_mom[2]:  # mom must be phased
                if kid_phased or kid_is_hom:
                    records_kid_mom.append({
                        "chrom": chrom,
                        "start": start,
                        "end": end,
                        "REF": REF,
                        "ALT": ALT,
                        "kid_allele_mat": str(gt_kid[1]) if gt_kid[1] != -1 else '.',
                        "mom_allele_C": str(gt_mom[0]),
                        "mom_allele_D": str(gt_mom[1]),
                    })

    vcf_reader.close()

    df_kid_dad = pl.DataFrame(records_kid_dad)
    df_kid_mom = pl.DataFrame(records_kid_mom)

    # Annotate with phase block boundaries
    df_kid_dad = _annotate_phase_blocks(df_kid_dad, df_blocks_kid, "kid")
    df_kid_dad = _annotate_phase_blocks(df_kid_dad, df_blocks_dad, "dad")
    df_kid_mom = _annotate_phase_blocks(df_kid_mom, df_blocks_kid, "kid")
    df_kid_mom = _annotate_phase_blocks(df_kid_mom, df_blocks_mom, "mom")

    return df_kid_dad, df_kid_mom


