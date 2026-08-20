#!/usr/bin/env python3
"""Summarize founder-phased methylation over raw inheritance-map rows."""

from __future__ import annotations

import argparse
import gzip
import json
import sys
from pathlib import Path
from typing import IO, Any


class MethylationSummaryError(RuntimeError):
    """Raised when visualization methylation inputs violate their contract."""


REQUIRED_COLUMNS = {
    "#chrom",
    "start_cpg",
    "methylation_level_model",
    "founder_haplotype_pat",
    "founder_haplotype_mat",
    "methylation_level_pat_model",
    "methylation_level_mat_model",
    "cpg_is_within_mismatch_window",
    "cpg_is_allele_specific",
}


def _open_text(path: Path) -> IO[str]:
    if path.name.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open(encoding="utf-8")


def _number(value: str) -> float | None:
    if value in {"", ".", "null", "NA", "nan"}:
        return None
    return float(value)


def _truth(value: str) -> bool:
    return value.lower() in {"1", "true", "yes"}


def read_inheritance_rows(path: Path, sample_id: str) -> dict[str, list[dict[str, Any]]]:
    """Read raw, shared-coordinate IHT rows for one sample."""
    header: list[str] | None = None
    rows: dict[str, list[dict[str, Any]]] = {}
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            fields = line.strip().split()
            if not fields:
                continue
            if header is None:
                header = fields
                if sample_id not in header:
                    raise MethylationSummaryError(
                        f"{path}: sample {sample_id!r} is absent from the inheritance map"
                    )
                continue
            if len(fields) != len(header):
                raise MethylationSummaryError(
                    f"{path} line {line_number}: expected {len(header)} fields, "
                    f"found {len(fields)}"
                )
            values = dict(zip(header, fields, strict=True))
            chromosome = values.get("#chrom", values.get("chrom"))
            if chromosome is None:
                raise MethylationSummaryError(f"{path}: missing chromosome column")
            start = int(values["start"])
            end = int(values["end"])
            if end < start:
                raise MethylationSummaryError(
                    f"{path}: invalid interval {chromosome}:{start}-{end}"
                )
            if end == start:
                end = start + 1
            cell = values[sample_id]
            labels = cell.split("|") if "|" in cell else []
            pat_label, mat_label = labels if len(labels) == 2 else (None, None)
            chromosome_rows = rows.setdefault(chromosome, [])
            if chromosome_rows and start < int(chromosome_rows[-1]["end"]):
                raise MethylationSummaryError(
                    f"{path}: overlapping inheritance rows on {chromosome}"
                )
            chromosome_rows.append(
                {
                    "start": start,
                    "end": end,
                    "pat_label": pat_label,
                    "mat_label": mat_label,
                    "pat_sum": 0.0,
                    "pat_n": 0,
                    "mat_sum": 0.0,
                    "mat_n": 0,
                    "unphased_sum": 0.0,
                    "unphased_n": 0,
                    "delta_sum": 0.0,
                    "delta_n": 0,
                    "mismatch_n": 0,
                    "allele_specific_n": 0,
                    "cpg_n": 0,
                }
            )
    if header is None or not rows:
        raise MethylationSummaryError(f"{path}: inheritance map has no data rows")
    return rows


def summarize(
    iht_path: Path, all_cpg_path: Path, sample_id: str
) -> dict[str, Any]:
    rows = read_inheritance_rows(iht_path, sample_id)
    cursors = {chromosome: 0 for chromosome in rows}
    columns: dict[str, int] | None = None
    last_position: dict[str, int] = {}

    with _open_text(all_cpg_path) as handle:
        for line_number, line in enumerate(handle, 1):
            if line.startswith("##") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if columns is None:
                if not fields[0].startswith("#"):
                    raise MethylationSummaryError(
                        f"{all_cpg_path} line {line_number}: missing column header"
                    )
                columns = {name: index for index, name in enumerate(fields)}
                missing = REQUIRED_COLUMNS - set(columns)
                if missing:
                    raise MethylationSummaryError(
                        f"{all_cpg_path}: missing columns {sorted(missing)}"
                    )
                continue
            if len(fields) != len(columns):
                raise MethylationSummaryError(
                    f"{all_cpg_path} line {line_number}: expected {len(columns)} "
                    f"fields, found {len(fields)}"
                )
            chromosome = fields[columns["#chrom"]]
            chromosome_rows = rows.get(chromosome)
            if not chromosome_rows:
                continue
            position = int(fields[columns["start_cpg"]])
            if position < last_position.get(chromosome, -1):
                raise MethylationSummaryError(
                    f"{all_cpg_path}: CpGs are not sorted on {chromosome}"
                )
            last_position[chromosome] = position
            cursor = cursors[chromosome]
            while cursor < len(chromosome_rows) and int(
                chromosome_rows[cursor]["end"]
            ) <= position:
                cursor += 1
            cursors[chromosome] = cursor
            if cursor == len(chromosome_rows):
                continue
            row = chromosome_rows[cursor]
            if position < int(row["start"]):
                continue

            row["cpg_n"] += 1
            unphased = _number(fields[columns["methylation_level_model"]])
            if unphased is not None:
                row["unphased_sum"] += unphased
                row["unphased_n"] += 1

            pat = _number(fields[columns["methylation_level_pat_model"]])
            mat = _number(fields[columns["methylation_level_mat_model"]])
            pat_matches = (
                row["pat_label"] is not None
                and fields[columns["founder_haplotype_pat"]] == row["pat_label"]
            )
            mat_matches = (
                row["mat_label"] is not None
                and fields[columns["founder_haplotype_mat"]] == row["mat_label"]
            )
            if pat_matches and pat is not None:
                row["pat_sum"] += pat
                row["pat_n"] += 1
            if mat_matches and mat is not None:
                row["mat_sum"] += mat
                row["mat_n"] += 1
            if pat_matches and mat_matches and pat is not None and mat is not None:
                row["delta_sum"] += pat - mat
                row["delta_n"] += 1
            if _truth(fields[columns["cpg_is_within_mismatch_window"]]):
                row["mismatch_n"] += 1
            if _truth(fields[columns["cpg_is_allele_specific"]]):
                row["allele_specific_n"] += 1

    if columns is None:
        raise MethylationSummaryError(f"{all_cpg_path}: no column header")

    compact: dict[str, list[list[Any]]] = {}
    for chromosome, chromosome_rows in rows.items():
        compact[chromosome] = [
            [
                row["start"],
                row["end"],
                row["pat_label"],
                row["mat_label"],
                round(float(row["pat_sum"]), 8),
                row["pat_n"],
                round(float(row["mat_sum"]), 8),
                row["mat_n"],
                round(float(row["unphased_sum"]), 8),
                row["unphased_n"],
                round(float(row["delta_sum"]), 8),
                row["delta_n"],
                row["mismatch_n"],
                row["allele_specific_n"],
                row["cpg_n"],
            ]
            for row in chromosome_rows
        ]
    return {
        "schema_version": 1,
        "sample_id": sample_id,
        "sources": {"inheritance_map": iht_path.name, "all_cpgs": all_cpg_path.name},
        "chromosomes": compact,
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--iht", required=True, type=Path)
    parser.add_argument("--all-cpgs", required=True, type=Path)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--output", required=True, type=Path)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        result = summarize(args.iht, args.all_cpgs, args.sample_id)
        args.output.write_text(
            json.dumps(result, separators=(",", ":")) + "\n", encoding="utf-8"
        )
    except (MethylationSummaryError, OSError, ValueError) as exc:
        print(f"haplotype methylation summary failed: {exc}", file=sys.stderr)
        return 2
    row_count = sum(len(value) for value in result["chromosomes"].values())
    print(f"summarized {args.sample_id}: {row_count} inheritance rows")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
