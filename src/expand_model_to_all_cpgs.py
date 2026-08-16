#!/usr/bin/env python3
"""Expand model-only founder-phased methylation to reference and sample CpGs."""

from __future__ import annotations

import argparse
import json
import sys
import tempfile
from pathlib import Path
from typing import Any

import numpy as np
import polars as pl
import pysam


class ExpansionError(RuntimeError):
    pass


PIPELINE_VERSION = "0.1.0-dev"


FOUNDER_COLUMNS = {
    "start_hap_map_block": pl.Int64,
    "end_hap_map_block": pl.Int64,
    "haplotype_concordance_in_hap_map_block": pl.Float64,
    "num_het_SNVs_in_hap_map_block": pl.Int64,
    "total_read_count_pat": pl.Int64,
    "total_read_count_mat": pl.Int64,
    "founder_haplotype_pat": pl.String,
    "founder_haplotype_mat": pl.String,
    "methylation_level_pat_count": pl.Float64,
    "methylation_level_mat_count": pl.Float64,
    "methylation_level_pat_model": pl.Float64,
    "methylation_level_mat_model": pl.Float64,
}


OUTPUT_COLUMNS = [
    "chrom",
    "start_cpg",
    "end_cpg",
    "total_read_count",
    "methylation_level_count",
    "methylation_level_model",
    "start_hap_map_block",
    "end_hap_map_block",
    "haplotype_concordance_in_hap_map_block",
    "num_het_SNVs_in_hap_map_block",
    "total_read_count_pat",
    "total_read_count_mat",
    "founder_haplotype_pat",
    "founder_haplotype_mat",
    "methylation_level_pat_count",
    "methylation_level_mat_count",
    "methylation_level_pat_model",
    "methylation_level_mat_model",
    "cpg_is_within_mismatch_window",
    "cpg_overlaps_at_least_one_snv",
    "snv_genotypes",
    "cpg_is_allele_specific",
]


def _read_model_bed(path: Path) -> pl.DataFrame:
    frame = pl.read_csv(
        path,
        separator="\t",
        comment_prefix="##",
        has_header=True,
        infer_schema_length=1_000_000,
    ).rename({"#chrom": "chrom", "begin": "start"})
    required = {"chrom", "start", "end", "mod_score", "cov"}
    if not required.issubset(frame.columns):
        raise ExpansionError(f"model BED is missing columns {sorted(required - set(frame.columns))}")
    return frame.select(
        "chrom",
        pl.col("start").cast(pl.Int64),
        pl.col("end").cast(pl.Int64),
        pl.col("cov").cast(pl.Int64).alias("total_read_count"),
        (pl.col("mod_score").cast(pl.Float64) / 100).alias("methylation_level_model"),
    )


def _read_reference_cpgs(path: Path) -> pl.DataFrame:
    return (
        pl.read_csv(
            path,
            separator="\t",
            comment_prefix="##",
            has_header=True,
            infer_schema_length=0,
        )
        .rename({"#chrom": "chrom"})
        .cast({"chrom": pl.String, "start": pl.Int64, "end": pl.Int64})
    )


def _read_legacy_table(data: Path, header: Path) -> pl.DataFrame:
    columns = [line.strip() for line in header.read_text(encoding="utf-8").splitlines()]
    if not columns:
        raise ExpansionError(f"empty table header: {header}")
    if data.stat().st_size == 0:
        schema = {column: FOUNDER_COLUMNS.get(column, pl.String) for column in columns}
        return pl.DataFrame(schema=schema)
    return pl.read_csv(
        data,
        separator="\t",
        has_header=False,
        new_columns=columns,
        infer_schema_length=1_000_000,
        null_values=[".", "null"],
    )


def _mismatch_proximity(
    cpgs: pl.DataFrame, mismatches: pl.DataFrame, window_bp: int
) -> pl.Series:
    result = np.zeros(len(cpgs), dtype=bool)
    if mismatches.is_empty() or cpgs.is_empty():
        return pl.Series("cpg_is_within_mismatch_window", result)
    chrom_values = cpgs["chrom"].to_numpy()
    starts = cpgs["start"].to_numpy()
    for chrom in cpgs["chrom"].unique(maintain_order=True).to_list():
        cpg_indexes = np.flatnonzero(chrom_values == chrom)
        mismatch_starts = (
            mismatches.filter(pl.col("chrom") == chrom)["start"].to_numpy()
        )
        if len(mismatch_starts) == 0:
            continue
        mismatch_starts = np.sort(mismatch_starts)
        query = starts[cpg_indexes]
        insertion = np.searchsorted(mismatch_starts, query)
        right = mismatch_starts[np.minimum(insertion, len(mismatch_starts) - 1)]
        left = mismatch_starts[np.maximum(insertion - 1, 0)]
        distance = np.minimum(np.abs(query - left), np.abs(query - right))
        result[cpg_indexes] = distance < window_bp
    return pl.Series("cpg_is_within_mismatch_window", result)


def _variant_annotations(vcf_path: Path, sample_id: str) -> pl.DataFrame:
    rows: list[dict[str, Any]] = []
    with pysam.VariantFile(str(vcf_path)) as vcf:
        if sample_id not in vcf.header.samples:
            raise ExpansionError(f"sample {sample_id!r} is absent from {vcf_path}")
        for record in vcf:
            reference = record.ref
            if (
                reference is None
                or len(reference) != 1
                or not record.alts
                or any(len(alt) != 1 for alt in record.alts)
            ):
                continue
            call = record.samples[sample_id]
            alleles = call.get("GT")
            if not alleles or any(allele is None for allele in alleles):
                genotype = "."
            elif len(set(alleles)) == 1:
                genotype = "hom"
            else:
                genotype = "het"
            variant_start = record.pos - 1
            for cpg_start in (variant_start - 1, variant_start):
                if cpg_start >= 0:
                    rows.append(
                        {
                            "chrom": record.contig,
                            "start": cpg_start,
                            "snv_genotype": genotype,
                        }
                    )
    if not rows:
        return pl.DataFrame(
            schema={
                "chrom": pl.String,
                "start": pl.Int64,
                "snv_genotypes": pl.String,
                "cpg_is_allele_specific": pl.Boolean,
            }
        )
    return (
        pl.DataFrame(rows)
        .group_by("chrom", "start", maintain_order=True)
        .agg(
            pl.col("snv_genotype").str.join(",").alias("snv_genotypes"),
            (pl.col("snv_genotype") == "het")
            .any()
            .alias("cpg_is_allele_specific"),
        )
    )


def expand_model_to_all_cpgs(
    reference_cpgs: Path,
    combined_model_bed: Path,
    founder_bed: Path,
    founder_header: Path,
    mismatch_bed: Path,
    mismatch_header: Path,
    joint_vcf: Path,
    sample_id: str,
    regions: list[str],
    mismatch_window_bp: int,
) -> pl.DataFrame:
    reference = _read_reference_cpgs(reference_cpgs)
    model = _read_model_bed(combined_model_bed)
    founder = _read_legacy_table(founder_bed, founder_header)
    mismatches = _read_legacy_table(mismatch_bed, mismatch_header)

    result = reference.join(model, on=["chrom", "start", "end"], how="full", coalesce=True)
    founder_keep = ["chrom", "start", "end", *FOUNDER_COLUMNS.keys()]
    founder_keep = list(dict.fromkeys(column for column in founder_keep if column in founder.columns))
    result = result.join(
        founder.select(founder_keep), on=["chrom", "start", "end"], how="left"
    )
    for column, dtype in FOUNDER_COLUMNS.items():
        if column not in result.columns:
            result = result.with_columns(pl.lit(None, dtype=dtype).alias(column))
        else:
            result = result.with_columns(pl.col(column).cast(dtype, strict=False))
    result = result.with_columns(
        pl.lit(None, dtype=pl.Float64).alias("methylation_level_count"),
        _mismatch_proximity(result, mismatches, mismatch_window_bp),
    )
    annotations = _variant_annotations(joint_vcf, sample_id)
    result = result.join(annotations, on=["chrom", "start"], how="left")
    result = (
        result.with_columns(
            pl.col("snv_genotypes")
            .is_not_null()
            .alias("cpg_overlaps_at_least_one_snv"),
            pl.col("snv_genotypes").fill_null("."),
            pl.col("cpg_is_allele_specific").fill_null(False),
        )
        .rename({"start": "start_cpg", "end": "end_cpg"})
        .filter(pl.col("chrom").is_in(regions))
        .sort(
            pl.col("chrom").str.extract(r"(\d+)").cast(pl.Int64),
            "start_cpg",
            "end_cpg",
        )
        .select(OUTPUT_COLUMNS)
    )
    return result


def write_output(
    frame: pl.DataFrame,
    output_bed: Path,
    metadata: dict[str, Any],
) -> None:
    output_bed.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".all-cpg-", dir=output_bed.parent) as temporary:
        temporary_dir = Path(temporary)
        plain = temporary_dir / "all-cpg.bed"
        compressed = temporary_dir / "all-cpg.bed.gz"
        with plain.open("w", encoding="utf-8") as destination:
            for key, value in metadata.items():
                destination.write(f"##{key}={json.dumps(value, sort_keys=True)}\n")
            destination.write("#" + "\t".join(frame.columns) + "\n")
            frame.write_csv(destination, separator="\t", include_header=False, null_value=".")
        pysam.tabix_compress(str(plain), str(compressed), force=True)
        pysam.tabix_index(str(compressed), preset="bed", force=True)
        compressed.replace(output_bed)
        Path(str(compressed) + ".tbi").replace(Path(str(output_bed) + ".tbi"))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-cpgs", required=True, type=Path)
    parser.add_argument("--combined-model-bed", required=True, type=Path)
    parser.add_argument("--founder-bed", required=True, type=Path)
    parser.add_argument("--founder-header", required=True, type=Path)
    parser.add_argument("--mismatch-bed", required=True, type=Path)
    parser.add_argument("--mismatch-header", required=True, type=Path)
    parser.add_argument("--phasing-qc", required=True, type=Path)
    parser.add_argument("--joint-vcf", required=True, type=Path)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--regions", required=True)
    parser.add_argument("--mismatch-window-bp", required=True, type=int)
    parser.add_argument("--reference-name", required=True)
    parser.add_argument("--min-coverage", required=True, type=int)
    parser.add_argument("--config-fingerprint", required=True)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    try:
        phasing_qc = json.loads(args.phasing_qc.read_text(encoding="utf-8"))
        if phasing_qc.get("sample_id") != args.sample_id:
            raise ExpansionError("phasing QC sample_id does not match --sample-id")
        phasing_status = phasing_qc.get("status")
        if phasing_status not in {"complete", "no_inheritance_phase"}:
            raise ExpansionError(f"unsupported phasing QC status: {phasing_status!r}")
        frame = expand_model_to_all_cpgs(
            args.reference_cpgs,
            args.combined_model_bed,
            args.founder_bed,
            args.founder_header,
            args.mismatch_bed,
            args.mismatch_header,
            args.joint_vcf,
            args.sample_id,
            args.regions.split(","),
            args.mismatch_window_bp,
        )
        output = args.output_dir / f"{args.sample_id}.dna-methylation.all-cpgs.bed.gz"
        write_output(
            frame,
            output,
            {
                "pipeline": "tapestry",
                "pipeline-version": PIPELINE_VERSION,
                "sample-id": args.sample_id,
                "reference": args.reference_name,
                "regions": args.regions.split(","),
                "min-coverage": args.min_coverage,
                "enabled-modes": ["model"],
                "config-fingerprint": args.config_fingerprint,
            },
        )
    except (ExpansionError, OSError, ValueError, pl.exceptions.PolarsError) as exc:
        print(f"all-CpG expansion failed: {exc}", file=sys.stderr)
        return 2
    qc = {
        "sample_id": args.sample_id,
        "status": "complete",
        "founder_phasing_status": phasing_status,
        "enabled_modes": ["model"],
        "records": len(frame),
        "model_records": frame["methylation_level_model"].is_not_null().sum(),
        "founder_model_records": (
            frame["methylation_level_pat_model"].is_not_null()
            | frame["methylation_level_mat_model"].is_not_null()
        ).sum(),
        "mismatch_proximal_records": frame["cpg_is_within_mismatch_window"].sum(),
        "allele_specific_records": frame["cpg_is_allele_specific"].sum(),
    }
    (args.output_dir / f"{args.sample_id}.all-cpgs-qc.json").write_text(
        json.dumps(qc, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"expanded {args.sample_id} to {len(frame)} CpGs")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
