#!/usr/bin/env python3
"""Filter upstream pb-CpG-tools model BEDs to the Tapestry coverage contract."""

from __future__ import annotations

import argparse
import gzip
import json
import sys
import tempfile
from pathlib import Path
from typing import Any, TextIO

import pysam


class ModelBedError(RuntimeError):
    pass


def _open_text(path: Path) -> TextIO:
    return gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else path.open(encoding="utf-8")


def filter_model_bed(
    input_bed: Path,
    input_index: Path,
    output_bed: Path,
    min_coverage: int,
    regions: list[str],
) -> dict[str, Any]:
    if min_coverage < 1:
        raise ModelBedError("minimum coverage must be at least one")
    for label, path in (("model BED", input_bed), ("model BED index", input_index)):
        if not path.is_file():
            raise ModelBedError(f"{label} does not exist: {path}")
    try:
        with pysam.TabixFile(str(input_bed), index=str(input_index)):
            pass
    except (OSError, ValueError) as exc:
        raise ModelBedError(f"cannot open model BED/index pair: {exc}") from exc

    output_bed.parent.mkdir(parents=True, exist_ok=True)
    upstream_min_coverage: str | None = None
    input_records = 0
    retained_records = 0
    excluded_region_records = 0
    low_coverage_records = 0
    region_set = set(regions)

    with tempfile.TemporaryDirectory(prefix=".model-bed-", dir=output_bed.parent) as temporary:
        temporary_dir = Path(temporary)
        plain = temporary_dir / "filtered.bed"
        compressed = temporary_dir / "filtered.bed.gz"
        header: list[str] | None = None
        with _open_text(input_bed) as source, plain.open("w", encoding="utf-8") as destination:
            for line_number, line in enumerate(source, 1):
                if line.startswith("##"):
                    if line.startswith("##min-coverage="):
                        upstream_min_coverage = line.rstrip("\n").split("=", 1)[1]
                    destination.write(line)
                    continue
                if line.startswith("#"):
                    header = line.rstrip("\n")[1:].split("\t")
                    if "chrom" not in header or "cov" not in header:
                        raise ModelBedError(
                            f"model BED header must contain chrom and cov columns: {input_bed}"
                        )
                    destination.write(f"##tapestry-min-coverage={min_coverage}\n")
                    destination.write("##tapestry-enabled-modes=model\n")
                    destination.write(line)
                    continue
                if not line.strip():
                    continue
                if header is None:
                    raise ModelBedError(
                        f"data appears before the model BED header at line {line_number}"
                    )
                values = line.rstrip("\n").split("\t")
                if len(values) != len(header):
                    raise ModelBedError(
                        f"line {line_number}: expected {len(header)} fields, found {len(values)}"
                    )
                row = dict(zip(header, values))
                input_records += 1
                if row["chrom"] not in region_set:
                    excluded_region_records += 1
                    continue
                try:
                    coverage = int(row["cov"])
                except ValueError as exc:
                    raise ModelBedError(
                        f"line {line_number}: invalid cov value {row['cov']!r}"
                    ) from exc
                if coverage < min_coverage:
                    low_coverage_records += 1
                    continue
                destination.write(line)
                retained_records += 1
        if header is None:
            raise ModelBedError(f"model BED has no # header: {input_bed}")

        pysam.tabix_compress(str(plain), str(compressed), force=True)
        pysam.tabix_index(str(compressed), preset="bed", force=True)
        compressed.replace(output_bed)
        Path(str(compressed) + ".tbi").replace(Path(str(output_bed) + ".tbi"))

    return {
        "input": input_bed.name,
        "output": output_bed.name,
        "input_records": input_records,
        "retained_records": retained_records,
        "excluded_region_records": excluded_region_records,
        "low_coverage_records": low_coverage_records,
        "upstream_min_coverage": upstream_min_coverage,
        "tapestry_min_coverage": min_coverage,
        "regions": regions,
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--combined-bed", required=True, type=Path)
    parser.add_argument("--combined-index", required=True, type=Path)
    parser.add_argument("--hap1-bed", required=True, type=Path)
    parser.add_argument("--hap1-index", required=True, type=Path)
    parser.add_argument("--hap2-bed", required=True, type=Path)
    parser.add_argument("--hap2-index", required=True, type=Path)
    parser.add_argument("--min-coverage", required=True, type=int)
    parser.add_argument("--regions", required=True, help="Comma-separated autosomes")
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    regions = args.regions.split(",")
    inputs = {
        "combined": (args.combined_bed, args.combined_index),
        "hap1": (args.hap1_bed, args.hap1_index),
        "hap2": (args.hap2_bed, args.hap2_index),
    }
    report: dict[str, Any] = {
        "sample_id": args.sample_id,
        "enabled_modes": ["model"],
        "beds": {},
    }
    try:
        for label, (bed, index) in inputs.items():
            output = args.output_dir / f"{args.sample_id}.model.{label}.bed.gz"
            report["beds"][label] = filter_model_bed(
                bed, index, output, args.min_coverage, regions
            )
    except ModelBedError as exc:
        print(f"model BED filtering failed: {exc}", file=sys.stderr)
        return 2
    qc_path = args.output_dir / f"{args.sample_id}.model-filter-qc.json"
    qc_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(
        f"filtered model BEDs for {args.sample_id} at coverage {args.min_coverage}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
