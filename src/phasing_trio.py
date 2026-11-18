import logging

from cyvcf2 import VCF  # type: ignore
import bioframe as bf
import polars as pl

from phasing import is_snv_het
# Shared with crossover flank annotation; keep the historical helper name.
from phase_sets import read_ps as _read_ps

logger = logging.getLogger(__name__)


def get_phase_blocks(blocks_tsv: str, uid: str) -> pl.DataFrame:
    """
    Read the whatshap stats --block-list TSV (produced by run-whatshap.sh) for a sample.
    Returns DataFrame with columns: chrom, start, end, phase_block_id, num_variants
    """
    df = (
        pl.read_csv(blocks_tsv, separator='\t', schema_overrides={
            '#sample': pl.String, 'chromosome': pl.String, 'phase_set': pl.String,
            'from': pl.Int64, 'to': pl.Int64, 'variants': pl.Int64,
        })
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


def _annotate_phase_blocks(df_snvs, df_blocks, label, require_phase_sets):
    """Attach each SNV's phase block for one individual, resolving overlaps by
    phase-set identity.

    bf.overlap matches a SNV to EVERY phase block whose [start, end) span
    contains it. WhatsHap blocks can nest — a short read-backed phase set sits
    inside the genomic span of a chromosome-spanning pedMEC set — so a single
    SNV can overlap several blocks. Assigning it to all of them fans the row
    out, and once the kid and parent annotations are combined a SNV multiplies
    into duplicate (and nested) hap-map blocks. We therefore keep exactly one
    block per SNV:

      * het + phased sites carry a PS (phase-set) tag in ``ps_{label}``; we keep
        the overlapping block whose id equals that PS — the phase set the SNV
        was actually phased in.
      * homozygous sites (kid only; the parent is always het at these sites)
        can lack PS. Their transmitted allele is phase-independent, so we attach
        them to the largest-span overlapping block (the chromosome-spanning
        pedMEC block) rather than to an incidental nested block.

    The ``ps_{label}`` column on ``df_snvs`` is dropped on the way out; the
    chosen block's id is kept as ``phase_block_id_{label}`` for grouping in
    _build_hap_map.
    """
    ps_col = f"ps_{label}"
    if df_snvs.is_empty() or df_blocks.is_empty():
        return df_snvs.clear().drop(ps_col).with_columns(
            pl.lit(None, dtype=pl.Int64).alias(f"start_phase_block_{label}"),
            pl.lit(None, dtype=pl.Int64).alias(f"end_phase_block_{label}"),
            pl.lit(None, dtype=pl.String).alias(f"phase_block_id_{label}"),
        )
    df = pl.from_pandas(bf.overlap(
        df_snvs.with_row_index("_marker_row").to_pandas(),
        df_blocks.select(["chrom", "start", "end", "phase_block_id"]).to_pandas(),
        how='inner',  # restrict to phased SNVs (those overlapping a block)
        suffixes=('', f'_{label}'),
    ))
    df = (
        df
        .rename({
            f"start_{label}": f"start_phase_block_{label}",
            f"end_{label}": f"end_phase_block_{label}",
        })
        .drop(f"chrom_{label}")
        .with_columns(
            _is_ps_match=(
                pl.col(ps_col).is_not_null()
                & (pl.col(f"phase_block_id_{label}").cast(pl.String) == pl.col(ps_col))
            ),
            _span=pl.col(f"end_phase_block_{label}") - pl.col(f"start_phase_block_{label}"),
        )
    )
    if require_phase_sets:
        # At this point a missing PS is allowed only for a homozygous child.
        df = df.filter(pl.col("_is_ps_match") | pl.col(ps_col).is_null())
    return (
        df
        # Keep one block per VCF row, with a deterministic choice for homozygotes.
        .sort(["_is_ps_match", "_span", f"phase_block_id_{label}"], descending=[True, True, False])
        .unique(subset="_marker_row", keep="first", maintain_order=True)
        .drop(["_is_ps_match", "_span", "_marker_row", ps_col])
        .sort(["chrom", "start", "end"])
    )

def get_all_phasing(
    vcf_path: str,
    kid_id: str, dad_id: str, mom_id: str,
    df_blocks_kid: pl.DataFrame,
    df_blocks_dad: pl.DataFrame,
    df_blocks_mom: pl.DataFrame,
    *,
    require_phase_sets: bool = False,
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

    With require_phase_sets, phased heterozygotes must have a PS that matches
    a reported block. A PS-free homozygous child uses the widest overlapping
    child block because its transmitted allele is phase-independent.

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

    # The nested-block disambiguation relies on the per-sample PS (phase-set)
    # FORMAT tag. The legacy path falls back to span-based assignment if the
    # VCF does not declare PS, which can assign markers to the wrong block.
    try:
        vcf_reader.get_header_type('PS')
        ps_in_header = True
    except KeyError:
        if require_phase_sets:
            vcf_reader.close()
            raise ValueError(f"No PS FORMAT tag declared in {vcf_path}") from None
        ps_in_header = False
        logger.warning(
            "No 'PS' FORMAT tag declared in %s; phase-set disambiguation will "
            "fall back to widest-span assignment and may produce duplicate or "
            "nested hap-map blocks. Verify the pedmec-phased VCF carries PS.",
            vcf_path,
        )

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
        ps_all = variant.format('PS') if ps_in_header else None
        gt_kid = genotype_all[si_kid]
        if len(gt_kid) != 3 or any(allele not in (0, 1) for allele in gt_kid[:2]):
            continue
        kid_phased = gt_kid[2]
        kid_is_hom = (gt_kid[0] == gt_kid[1])
        ps_kid = _read_ps(ps_all, si_kid)
        if not (kid_is_hom or (kid_phased and (ps_kid is not None or not require_phase_sets))):
            continue

        # Dad-het sites → paternal allele in kid
        if is_snv_het(variant, si_dad):
            gt_dad = genotype_all[si_dad]
            ps_dad = _read_ps(ps_all, si_dad)
            if len(gt_dad) == 3 and gt_dad[2] and (ps_dad is not None or not require_phase_sets):
                records_kid_dad.append({
                    "chrom": chrom,
                    "start": start,
                    "end": end,
                    "REF": REF,
                    "ALT": ALT,
                    "kid_allele_pat": str(gt_kid[0]),
                    "dad_allele_A": str(gt_dad[0]),
                    "dad_allele_B": str(gt_dad[1]),
                    "ps_kid": ps_kid,
                    "ps_dad": ps_dad,
                })

        # Mom-het sites → maternal allele in kid
        if is_snv_het(variant, si_mom):
            gt_mom = genotype_all[si_mom]
            ps_mom = _read_ps(ps_all, si_mom)
            if len(gt_mom) == 3 and gt_mom[2] and (ps_mom is not None or not require_phase_sets):
                records_kid_mom.append({
                    "chrom": chrom,
                    "start": start,
                    "end": end,
                    "REF": REF,
                    "ALT": ALT,
                    "kid_allele_mat": str(gt_kid[1]),
                    "mom_allele_C": str(gt_mom[0]),
                    "mom_allele_D": str(gt_mom[1]),
                    "ps_kid": ps_kid,
                    "ps_mom": ps_mom,
                })

    vcf_reader.close()

    marker_schema = {"chrom": pl.String, "start": pl.Int64, "end": pl.Int64,
                     "REF": pl.String, "ALT": pl.String, "ps_kid": pl.String}
    df_kid_dad = pl.DataFrame(records_kid_dad, schema={
        **marker_schema, "kid_allele_pat": pl.String, "dad_allele_A": pl.String,
        "dad_allele_B": pl.String, "ps_dad": pl.String,
    })
    df_kid_mom = pl.DataFrame(records_kid_mom, schema={
        **marker_schema, "kid_allele_mat": pl.String, "mom_allele_C": pl.String,
        "mom_allele_D": pl.String, "ps_mom": pl.String,
    })

    # Annotate with phase block boundaries
    df_kid_dad = _annotate_phase_blocks(df_kid_dad, df_blocks_kid, "kid", require_phase_sets)
    df_kid_dad = _annotate_phase_blocks(df_kid_dad, df_blocks_dad, "dad", require_phase_sets)
    df_kid_mom = _annotate_phase_blocks(df_kid_mom, df_blocks_kid, "kid", require_phase_sets)
    df_kid_mom = _annotate_phase_blocks(df_kid_mom, df_blocks_mom, "mom", require_phase_sets)

    # Cast phase block columns from f64 (pandas round-trip artifact) to Int64
    phase_block_cols_dad = {
        col: pl.Int64 for col in df_kid_dad.columns if col.startswith(("start_phase_block", "end_phase_block"))
    }
    phase_block_cols_mom = {
        col: pl.Int64 for col in df_kid_mom.columns if col.startswith(("start_phase_block", "end_phase_block"))
    }
    df_kid_dad = df_kid_dad.cast(phase_block_cols_dad)
    df_kid_mom = df_kid_mom.cast(phase_block_cols_mom)

    return df_kid_dad, df_kid_mom


