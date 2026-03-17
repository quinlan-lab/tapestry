from cyvcf2 import VCF  # type: ignore
import polars as pl

from phasing import stringify, is_snv_het

def get_pedmec_phasing(vcf_path: str, kid_id: str, dad_id: str, mom_id: str):
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
                phase_block_id = stringify(pb)

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


