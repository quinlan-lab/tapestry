import argparse
from contextlib import ExitStack
import json
import logging
from pathlib import Path
from typing import cast
import bioframe as bf # https://bioframe.readthedocs.io/en/latest/index.html
import polars as pl
import pyBigWig

from phasing_pedigree import (
    get_read_phasing, 
    get_read_phase_blocks, 
    get_iht_phasing, 
    get_iht_blocks, 
    get_all_phasing
)
from hap_map_pedigree import get_hap_map
from util.hap_map import write_hap_map_blocks
from get_meth_hap1_hap2 import (
    IndexedMethylationBed,
    read_meth_hap1_hap2_chromosome,
)
from util.write_data import (
    write_bed,
    write_bit_vector_mismatches_bed,
    write_bit_vector_mismatches_vcf,
    write_header,
)
from util.version_sort import version_sort

REFERENCE_GENOME = "hg38"

def phase_meth_to_founder_haps(df_meth_hap1_hap2, df_hap_map):
    if df_hap_map.is_empty() or df_meth_hap1_hap2.is_empty():
        return (
            df_meth_hap1_hap2
            .with_columns(
                pl.lit(None, dtype=pl.Int64).alias("start_hap_map_block"),
                pl.lit(None, dtype=pl.Int64).alias("end_hap_map_block"),
                pl.lit(None, dtype=pl.Float64).alias(
                    "haplotype_concordance_in_hap_map_block"
                ),
                pl.lit(None, dtype=pl.Int64).alias(
                    "num_het_SNVs_in_hap_map_block"
                ),
                pl.lit(None, dtype=pl.Float64).alias("methylation_level_pat"),
                pl.lit(None, dtype=pl.Float64).alias("methylation_level_mat"),
                pl.lit(None, dtype=pl.Int64).alias("total_read_count_pat"),
                pl.lit(None, dtype=pl.Int64).alias("total_read_count_mat"),
                pl.lit(None, dtype=pl.String).alias("founder_haplotype_pat"),
                pl.lit(None, dtype=pl.String).alias("founder_haplotype_mat"),
            )
            .drop(
                [
                    "total_read_count_hap1",
                    "total_read_count_hap2",
                    "methylation_level_hap1",
                    "methylation_level_hap2",
                ]
            )
        )

    df = bf.overlap(
        df_meth_hap1_hap2.to_pandas(), 
        df_hap_map.to_pandas(), 
        how='left', # we want to retain all CpG sites, including those that we cannot phase to founder haplotypes
        suffixes=('','_'),
    )

    df = (
        pl
        .from_pandas(df)
        .drop([
            "chrom_"
        ])
        .rename({ 
            "start_": "start_hap_map_block",
            "end_": "end_hap_map_block",
            "paternal_haplotype_": "paternal_haplotype_in_hap_map_block", 
            "maternal_haplotype_": "maternal_haplotype_in_hap_map_block",
            "haplotype_concordance_": "haplotype_concordance_in_hap_map_block",
            "num_het_SNVs_": "num_het_SNVs_in_hap_map_block",
        })
        .cast({
            "num_het_SNVs_in_hap_map_block": pl.Int64,
            "total_read_count_hap1": pl.Int64,
            "total_read_count_hap2": pl.Int64,
        })
        .with_columns(
            methylation_level_pat=(
                pl
                .when(pl.col("paternal_haplotype_in_hap_map_block").str.ends_with("_hap1"))
                .then(pl.col("methylation_level_hap1"))
                .when(pl.col("paternal_haplotype_in_hap_map_block").str.ends_with("_hap2"))
                .then(pl.col("methylation_level_hap2"))
                .otherwise(None)
            ),
            methylation_level_mat=(
                pl
                .when(pl.col("maternal_haplotype_in_hap_map_block").str.ends_with("_hap1"))
                .then(pl.col("methylation_level_hap1"))
                .when(pl.col("maternal_haplotype_in_hap_map_block").str.ends_with("_hap2"))
                .then(pl.col("methylation_level_hap2"))
                .otherwise(None)
            ),
            total_read_count_pat=(
                pl
                .when(pl.col("paternal_haplotype_in_hap_map_block").str.ends_with("_hap1"))
                .then(pl.col("total_read_count_hap1"))
                .when(pl.col("paternal_haplotype_in_hap_map_block").str.ends_with("_hap2"))
                .then(pl.col("total_read_count_hap2"))
                .otherwise(None)
            ),
            total_read_count_mat=(
                pl
                .when(pl.col("maternal_haplotype_in_hap_map_block").str.ends_with("_hap1"))
                .then(pl.col("total_read_count_hap1"))
                .when(pl.col("maternal_haplotype_in_hap_map_block").str.ends_with("_hap2"))
                .then(pl.col("total_read_count_hap2"))
                .otherwise(None)
            ),
            founder_haplotype_pat=pl.col("paternal_haplotype_in_hap_map_block").str.split("_").list.get(0),
            founder_haplotype_mat=pl.col("maternal_haplotype_in_hap_map_block").str.split("_").list.get(0),
        )
        .drop([
            "total_read_count_hap1", "total_read_count_hap2",
            "methylation_level_hap1", "methylation_level_hap2",
            "paternal_haplotype_in_hap_map_block", "maternal_haplotype_in_hap_map_block",
        ])
    )

    return df 

def _read_fai_chromsizes(reference_fai):
    chromsizes = []
    with open(reference_fai) as handle:
        for line_number, line in enumerate(handle, 1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                raise ValueError(
                    f"Malformed FASTA index line {line_number} in {reference_fai}"
                )
            chromsizes.append((fields[0], int(fields[1])))
    return chromsizes


def add_bigwig_entries(bigwig, df, parental, pb_cpg_tool_mode, chunk_size=100_000):
    """Append one unambiguous value per interval in bounded Python lists."""
    value_column = f"methylation_level_{parental}_{pb_cpg_tool_mode}"
    interval_columns = ["chrom", "start", "end"]
    df_bed_graph = (
        df
        .filter(pl.col(value_column).is_not_null())
        .select(interval_columns + [value_column])
        .group_by(interval_columns)
        .agg(
            pl.col(value_column).first().alias(value_column),
            pl.col(value_column).n_unique().alias("_distinct_values"),
        )
        .sort(["chrom", "start", "end"])
    )

    conflicts = df_bed_graph.filter(pl.col("_distinct_values") > 1)
    if not conflicts.is_empty():
        examples = ", ".join(
            f"{chrom}:{start}-{end}"
            for chrom, start, end in conflicts.select(interval_columns)
            .head(5)
            .iter_rows()
        )
        raise ValueError(
            f"Conflicting {parental} {pb_cpg_tool_mode} BigWig values at "
            f"{len(conflicts)} interval(s); examples: {examples}"
        )

    df_bed_graph = df_bed_graph.drop("_distinct_values")
    for offset in range(0, len(df_bed_graph), chunk_size):
        chunk = df_bed_graph.slice(offset, chunk_size)
        bigwig.addEntries(
            chunk["chrom"].to_list(),
            chunk["start"].to_list(),
            ends=chunk["end"].to_list(),
            values=chunk[value_column].to_list(),
        )


def write_bigwig(
    df,
    uid,
    parental,
    pb_cpg_tool_mode,
    output_dir,
    logger=None,
    reference_fai=None,
    reference_name=REFERENCE_GENOME,
):
    """Write one complete BigWig, adding entries in bounded chunks."""

    file_path = f"{output_dir}/{uid}.dna-methylation.{parental}.{pb_cpg_tool_mode}.{reference_name}.bw"
    if reference_fai is not None:
        chromsizes = _read_fai_chromsizes(reference_fai)
    else:
        # Legacy CLI compatibility. Generic workflow calls must provide an FAI.
        legacy_sizes = bf.fetch_chromsizes(db=REFERENCE_GENOME)
        chromsizes = list(zip(legacy_sizes.index, legacy_sizes.values))

    size_by_chrom = dict(chromsizes)
    unknown = sorted(set(df["chrom"].unique().to_list()) - set(size_by_chrom))
    if unknown:
        raise ValueError(f"BigWig records use contigs absent from the FASTA index: {unknown}")

    bigwig = pyBigWig.open(file_path, "w")
    try:
        bigwig.addHeader(chromsizes)
        add_bigwig_entries(bigwig, df, parental, pb_cpg_tool_mode)
    finally:
        bigwig.close()

    if logger: 
        logger.info(
            "Wrote bigwig file for %s %s-based methylation levels against %s to: '%s'",
            parental,
            pb_cpg_tool_mode,
            reference_name,
            file_path,
        )

def add_suffix(
    df: pl.DataFrame, suffix: str, cols_to_suffix: list[str]
) -> pl.DataFrame:
    df.columns = [f"{col}{suffix}" if col in cols_to_suffix else col for col in df.columns]
    return df 

def combine_count_and_model_based_methylation_levels(
        df_meth_count: pl.DataFrame | None,
        df_meth_model: pl.DataFrame
    ):
    suffix_count, suffix_model = '_count', '_model'
    unique_cols = ['methylation_level_pat', 'methylation_level_mat']

    df_meth_model = add_suffix(df_meth_model, suffix=suffix_model, cols_to_suffix=unique_cols)

    unique_cols_count = [f'{unique_col}{suffix_count}' for unique_col in unique_cols]
    unique_cols_model = [f'{unique_col}{suffix_model}' for unique_col in unique_cols]

    if df_meth_count is None:
        common_cols = [col for col in df_meth_model.columns if col not in unique_cols_model]
        return version_sort(
            df_meth_model
            .with_columns(
                pl.lit(None, dtype=pl.Float64).alias(unique_cols_count[0]),
                pl.lit(None, dtype=pl.Float64).alias(unique_cols_count[1]),
            )
            .select(common_cols + unique_cols_count + unique_cols_model)
        )

    count_frame = add_suffix(
        cast(pl.DataFrame, df_meth_count),
        suffix=suffix_count,
        cols_to_suffix=unique_cols,
    )
    common_cols = [col for col in count_frame.columns if col not in unique_cols_count]

    df = count_frame.join(
        df_meth_model,
        on=common_cols,
        nulls_equal=True, # join on nulls, e.g., total_read_count_XXX is null in both df_meth_count and df_meth_model
        # "how=full" will capture CpG sites where count-based meth levels are available, but not model-based meth levels, and vice versa,
        # though, for sample 200081, I didn't see any such CpG sites: 
        how="full", 
        coalesce=True, 
    )   

    df = df.select(common_cols + unique_cols_count + unique_cols_model)

    return version_sort(df)


def _autosome_sort_key(chrom: str) -> int:
    if not chrom.startswith("chr") or not chrom[3:].isdigit():
        raise ValueError(f"Expected an autosome named chrN, found {chrom!r}")
    number = int(chrom[3:])
    if number < 1 or number > 22:
        raise ValueError(f"Expected chr1-chr22, found {chrom!r}")
    return number


def _ordered_autosomes(regions: list[str]) -> list[str]:
    ordered = sorted(regions, key=_autosome_sort_key)
    if len(ordered) != len(set(ordered)):
        raise ValueError("Regions contain duplicates")
    return ordered


def _phase_one_chromosome(
    *,
    uid: str,
    chrom: str,
    vcf_read_phased: str,
    vcf_iht_phased: str,
    read_phase_blocks: pl.DataFrame,
    iht_blocks: pl.DataFrame,
) -> tuple[pl.DataFrame, pl.DataFrame, int, int, int]:
    """Build compact hap-map outputs while large chromosome tables are local."""
    read_phasing = get_read_phasing(vcf_read_phased, chrom)
    iht_phasing = get_iht_phasing(uid, vcf_iht_phased, chrom)
    all_phasing = get_all_phasing(
        read_phasing,
        read_phase_blocks.filter(pl.col("chrom") == chrom),
        iht_phasing,
        iht_blocks.filter(pl.col("chrom") == chrom),
    )
    hap_map, informative_sites, mismatches = get_hap_map(all_phasing)
    return (
        hap_map,
        mismatches,
        len(informative_sites),
        len(read_phasing),
        len(iht_phasing),
    )


def build_hap_map_by_chromosome(
    *,
    uid: str,
    regions: list[str],
    vcf_read_phased: str,
    tsv_read_phase_blocks: str,
    vcf_iht_phased: str,
    txt_iht_blocks: str,
    logger: logging.Logger,
) -> tuple[pl.DataFrame, pl.DataFrame, int]:
    """Reconcile VCF phasing in chromosome-sized batches."""
    ordered_regions = _ordered_autosomes(regions)
    read_phase_blocks = get_read_phase_blocks(tsv_read_phase_blocks)
    iht_blocks = get_iht_blocks(uid, txt_iht_blocks)
    logger.info(
        "Loaded %s HiPhase blocks and %s inheritance blocks",
        len(read_phase_blocks),
        len(iht_blocks),
    )

    hap_maps: list[pl.DataFrame] = []
    mismatch_frames: list[pl.DataFrame] = []
    informative_site_count = 0
    for chrom in ordered_regions:
        hap_map, mismatches, informative, read_hets, inheritance_hets = (
            _phase_one_chromosome(
                uid=uid,
                chrom=chrom,
                vcf_read_phased=vcf_read_phased,
                vcf_iht_phased=vcf_iht_phased,
                read_phase_blocks=read_phase_blocks,
                iht_blocks=iht_blocks,
            )
        )
        hap_maps.append(hap_map)
        mismatch_frames.append(mismatches)
        informative_site_count += informative
        logger.info(
            "%s: read-phased hets=%s, inheritance hets=%s, informative "
            "sites=%s, hap-map blocks=%s, mismatches=%s",
            chrom,
            read_hets,
            inheritance_hets,
            informative,
            len(hap_map),
            len(mismatches),
        )

    return (
        version_sort(pl.concat(hap_maps, rechunk=True)),
        version_sort(pl.concat(mismatch_frames, rechunk=True)),
        informative_site_count,
    )


def phase_methylation_by_chromosome(
    *,
    uid: str,
    regions: list[str],
    df_hap_map: pl.DataFrame,
    model_hap1: str,
    model_hap2: str,
    output_dir: Path,
    reference_fai: str | None,
    reference_name: str,
    write_bigwigs: bool,
    logger: logging.Logger,
    count_hap1: str | None = None,
    count_hap2: str | None = None,
) -> tuple[list[str], int]:
    """Phase and publish methylation one chromosome at a time."""
    ordered_regions = _ordered_autosomes(regions)
    if (count_hap1 is None) != (count_hap2 is None):
        raise ValueError("Count hap1 and hap2 BEDs must be provided together")
    enabled_modes = ["model"] if count_hap1 is None else ["count", "model"]
    output_path = output_dir / f"{uid}.dna-methylation.bed"
    header_path = output_dir / f"{uid}.dna-methylation.bed.header"
    output_columns: list[str] | None = None
    total_rows = 0

    chromsizes: list[tuple[str, int]] = []
    if write_bigwigs:
        if reference_fai is not None:
            chromsizes = _read_fai_chromsizes(reference_fai)
        else:
            legacy_sizes = bf.fetch_chromsizes(db=REFERENCE_GENOME)
            chromsizes = list(zip(legacy_sizes.index, legacy_sizes.values))
        size_by_chrom = dict(chromsizes)
        missing_regions = [
            region for region in ordered_regions if region not in size_by_chrom
        ]
        if missing_regions:
            raise ValueError(f"Regions are absent from the FASTA index: {missing_regions}")

    with ExitStack() as stack:
        model_hap1_reader = stack.enter_context(
            IndexedMethylationBed(model_hap1, "model")
        )
        model_hap2_reader = stack.enter_context(
            IndexedMethylationBed(model_hap2, "model")
        )
        count_hap1_reader = None
        count_hap2_reader = None
        if count_hap1 is not None and count_hap2 is not None:
            count_hap1_reader = stack.enter_context(
                IndexedMethylationBed(count_hap1, "count")
            )
            count_hap2_reader = stack.enter_context(
                IndexedMethylationBed(count_hap2, "count")
            )

        output_handle = stack.enter_context(output_path.open("w", encoding="utf-8"))
        bigwigs: dict[tuple[str, str], object] = {}
        if write_bigwigs:
            for parental in ("pat", "mat"):
                for mode in enabled_modes:
                    path = output_dir / (
                        f"{uid}.dna-methylation.{parental}.{mode}."
                        f"{reference_name}.bw"
                    )
                    handle = pyBigWig.open(str(path), "w")
                    handle.addHeader(chromsizes)
                    stack.callback(handle.close)
                    bigwigs[(parental, mode)] = handle

        for chrom in ordered_regions:
            model = read_meth_hap1_hap2_chromosome(
                model_hap1_reader, model_hap2_reader, chrom
            )
            chrom_hap_map = df_hap_map.filter(pl.col("chrom") == chrom)
            model_founder = phase_meth_to_founder_haps(model, chrom_hap_map)

            count_founder = None
            if count_hap1_reader is not None and count_hap2_reader is not None:
                count = read_meth_hap1_hap2_chromosome(
                    count_hap1_reader, count_hap2_reader, chrom
                )
                count_founder = phase_meth_to_founder_haps(count, chrom_hap_map)

            combined = combine_count_and_model_based_methylation_levels(
                count_founder, model_founder
            )
            if output_columns is None:
                output_columns = combined.columns
                write_header(header_path, output_columns)
            elif set(combined.columns) != set(output_columns):
                raise ValueError(
                    f"Output schema changed while processing {chrom}: {combined.columns}"
                )
            else:
                combined = combined.select(output_columns)

            combined.write_csv(output_handle, separator="\t", include_header=False)
            total_rows += len(combined)
            for (parental, mode), bigwig in bigwigs.items():
                add_bigwig_entries(bigwig, combined, parental, mode)
            logger.info(
                "Phased %s model CpGs to founder haplotypes on %s",
                len(combined),
                chrom,
            )

    if output_columns is None:
        raise ValueError("At least one autosome is required")
    return enabled_modes, total_rows

def main():
    parser = argparse.ArgumentParser(description='Phase HiFi-based DNA methylation data to founder haplotypes')
    parser.add_argument('--uid', required=True, help='Sample UID in joint-called multi-sample vcf')
    parser.add_argument('--vcf_read_phased', required=True, help='Single-sample vcf from hiphase')
    parser.add_argument('--tsv_read_phase_blocks', required=True, help='Single-sample tsv from hiphase')
    parser.add_argument('--vcf_iht_phased', required=True, help='Joint-called multi-sample vcf from gtg-ped-map/gtg-concordance')
    parser.add_argument('--txt_iht_blocks', required=True, help='Multi-sample iht blocks file from gtg-ped-map/gtg-concordance')
    parser.add_argument('--bed_meth_count_hap1', help='Optional BED of count-based methylation levels for hap1')
    parser.add_argument('--bed_meth_count_hap2', help='Optional BED of count-based methylation levels for hap2')
    parser.add_argument('--bed_meth_model_hap1', required=True, help='Bed file of model-based methylation levels from aligned_bam_to_cpg_scores for hap1')
    parser.add_argument('--bed_meth_model_hap2', required=True, help='Bed file of model-based methylation levels from aligned_bam_to_cpg_scores for hap2')
    parser.add_argument('--reference_fai', help='FASTA index used for BigWig chromosome lengths')
    parser.add_argument('--reference_name', default=REFERENCE_GENOME, help='Reference label used in BigWig filenames')
    parser.add_argument('--regions', required=True, help='Comma-separated autosomes to process')
    parser.add_argument('--no-bigwig', action='store_true', help='Do not write BigWig outputs')
    parser.add_argument('--output_dir', required=True, help='Output directory')
    args = parser.parse_args()

    if bool(args.bed_meth_count_hap1) != bool(args.bed_meth_count_hap2):
        parser.error("--bed_meth_count_hap1 and --bed_meth_count_hap2 must be provided together")

    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(filename)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    logger = logging.getLogger(__name__)

    logger.info(f"Starting '{__file__}'")
    logger.info("Script started with the following arguments: %s", vars(args))

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    regions = [region.strip() for region in args.regions.split(",") if region.strip()]
    if not regions:
        parser.error("--regions must contain at least one autosome")
    
    df_hap_map, df_sites_mismatch, informative_site_count = (
        build_hap_map_by_chromosome(
            uid=args.uid,
            regions=regions,
            vcf_read_phased=args.vcf_read_phased,
            tsv_read_phase_blocks=args.tsv_read_phase_blocks,
            vcf_iht_phased=args.vcf_iht_phased,
            txt_iht_blocks=args.txt_iht_blocks,
            logger=logger,
        )
    )
    logger.info(f"Got hap map: {len(df_hap_map)} rows, {len(df_hap_map.columns)} columns")
    logger.info("Got %s informative heterozygous sites", informative_site_count)
    logger.info(f"Got heterozygous sites where read-based and inheritance-based bit vectors don't match: {len(df_sites_mismatch)} rows, {len(df_sites_mismatch.columns)} columns")
    
    write_hap_map_blocks(df_hap_map, args.uid, "paternal", args.output_dir)
    write_hap_map_blocks(df_hap_map, args.uid, "maternal", args.output_dir)
    logger.info(f"Wrote paternal and maternal hap-map blocks for IGV visualization to '{args.output_dir}'")
    
    write_bed(args.output_dir, df_hap_map, f"{args.uid}.hap-map-blocks")
    logger.info(f"Wrote hap-map blocks to '{args.output_dir}'")
    
    write_bit_vector_mismatches_vcf(args.output_dir, df_sites_mismatch, logger, uid=args.uid)
    write_bit_vector_mismatches_bed(args.output_dir, df_sites_mismatch, logger, uid=args.uid)    

    mismatch_site_count = len(df_sites_mismatch)

    enabled_modes, methylation_rows = phase_methylation_by_chromosome(
        uid=args.uid,
        regions=regions,
        df_hap_map=df_hap_map,
        model_hap1=args.bed_meth_model_hap1,
        model_hap2=args.bed_meth_model_hap2,
        count_hap1=args.bed_meth_count_hap1,
        count_hap2=args.bed_meth_count_hap2,
        output_dir=output_dir,
        reference_fai=args.reference_fai,
        reference_name=args.reference_name,
        write_bigwigs=not args.no_bigwig,
        logger=logger,
    )
    logger.info(
        "Wrote %s enabled methylation rows in chromosome-sized batches to %s",
        methylation_rows,
        output_dir,
    )

    qc_status = "no_inheritance_phase" if df_hap_map.is_empty() else "complete"
    qc = {
        "sample_id": args.uid,
        "status": qc_status,
        "enabled_modes": enabled_modes,
        "hap_map_blocks": len(df_hap_map),
        "informative_sites": informative_site_count,
        "mismatch_sites": mismatch_site_count,
    }
    (Path(args.output_dir) / f"{args.uid}.phasing-qc.json").write_text(
        json.dumps(qc, indent=2, sort_keys=True) + "\n"
    )
    
    logger.info(f"Done running '{__file__}'")

if __name__ == "__main__":
    main()
