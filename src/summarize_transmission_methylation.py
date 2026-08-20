#!/usr/bin/env python3
"""Summarize methylation concordance on transmitted haplotypes by chromosome."""

from __future__ import annotations

import argparse
import gzip
import json
import re
import sys
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path
from typing import IO, Any


class TransmissionSummaryError(RuntimeError):
    """Raised when parent-child methylation inputs violate their contract."""


MISSING_PARENTS = {"0", "NA", "N/A", ".", "-"}
BED_SUFFIX = ".dna-methylation.all-cpgs.bed.gz"
REQUIRED_COLUMNS = {
    "#chrom",
    "start_cpg",
    "end_cpg",
    "founder_haplotype_pat",
    "founder_haplotype_mat",
    "methylation_level_pat_model",
    "methylation_level_mat_model",
    "cpg_is_within_mismatch_window",
}


@dataclass(frozen=True)
class Cpg:
    chromosome: str
    position: int
    end: int
    pat_label: str | None
    mat_label: str | None
    pat_value: float | None
    mat_value: float | None
    mismatch: bool

    def value_for_label(self, label: str | None) -> tuple[float | None, float | None]:
        """Return methylation on ``label`` and on the other haplotype."""
        if label is None:
            return None, None
        if self.pat_label == label:
            return self.pat_value, self.mat_value
        if self.mat_label == label:
            return self.mat_value, self.pat_value
        return None, None


def _open_text(path: Path) -> IO[str]:
    if path.name.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open(encoding="utf-8")


def _number(value: str) -> float | None:
    if value in {"", ".", "null", "NA", "nan"}:
        return None
    result = float(value)
    if not 0.0 <= result <= 1.0:
        raise TransmissionSummaryError(f"methylation value outside [0, 1]: {value}")
    return result


def _label(value: str) -> str | None:
    return None if value in {"", ".", "null", "NA", "?"} else value


def _truth(value: str) -> bool:
    normalized = value.lower()
    if normalized in {"1", "true", "yes"}:
        return True
    if normalized in {"0", "false", "no"}:
        return False
    raise TransmissionSummaryError(f"invalid boolean value: {value!r}")


def _chromosome_key(chromosome: str) -> tuple[int, int | str]:
    match = re.fullmatch(r"chr(\d+)", chromosome, flags=re.IGNORECASE)
    if match:
        return (0, int(match.group(1)))
    return (1, chromosome)


def read_pedigree(path: Path) -> dict[str, dict[str, str | None]]:
    people: dict[str, dict[str, str | None]] = {}
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            fields = stripped.split()
            if len(fields) != 6:
                raise TransmissionSummaryError(
                    f"{path} line {line_number}: expected six PED columns"
                )
            family, sample, father, mother, _sex, _phenotype = fields
            if sample in people:
                raise TransmissionSummaryError(f"{path}: duplicate sample {sample!r}")
            people[sample] = {
                "family": family,
                "father": None if father.upper() in MISSING_PARENTS else father,
                "mother": None if mother.upper() in MISSING_PARENTS else mother,
            }
    if not people:
        raise TransmissionSummaryError(f"{path}: no pedigree records")
    return people


def sample_id_from_path(path: Path) -> str:
    if not path.name.endswith(BED_SUFFIX):
        raise TransmissionSummaryError(
            f"{path}: expected a filename ending in {BED_SUFFIX!r}"
        )
    sample = path.name[: -len(BED_SUFFIX)]
    if not sample:
        raise TransmissionSummaryError(f"{path}: could not determine sample ID")
    return sample


def _same_number(left: float | None, right: float | None) -> bool:
    return left is None or right is None or abs(left - right) <= 1e-12


def _merge_duplicate(left: Cpg, right: Cpg, path: Path) -> Cpg:
    if (left.chromosome, left.position, left.end) != (
        right.chromosome,
        right.position,
        right.end,
    ):
        raise AssertionError("duplicate merge received different CpGs")
    if left.pat_label != right.pat_label or left.mat_label != right.mat_label:
        raise TransmissionSummaryError(
            f"{path}: conflicting founder labels at {left.chromosome}:{left.position}"
        )
    if not _same_number(left.pat_value, right.pat_value) or not _same_number(
        left.mat_value, right.mat_value
    ):
        raise TransmissionSummaryError(
            f"{path}: conflicting methylation values at "
            f"{left.chromosome}:{left.position}"
        )
    return Cpg(
        chromosome=left.chromosome,
        position=left.position,
        end=left.end,
        pat_label=left.pat_label,
        mat_label=left.mat_label,
        pat_value=left.pat_value if left.pat_value is not None else right.pat_value,
        mat_value=left.mat_value if left.mat_value is not None else right.mat_value,
        mismatch=left.mismatch or right.mismatch,
    )


def iter_cpgs(path: Path) -> Iterator[Cpg]:
    """Stream sorted CpGs, collapsing only identical duplicate observations."""
    columns: dict[str, int] | None = None
    pending: Cpg | None = None
    last_key: tuple[tuple[int, int | str], int, int] | None = None
    with _open_text(path) as handle:
        for line_number, line in enumerate(handle, 1):
            if line.startswith("##") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if columns is None:
                if not fields[0].startswith("#"):
                    raise TransmissionSummaryError(
                        f"{path} line {line_number}: missing column header"
                    )
                if len(fields) != len(set(fields)):
                    raise TransmissionSummaryError(f"{path}: duplicate column name")
                columns = {name: index for index, name in enumerate(fields)}
                missing = REQUIRED_COLUMNS - set(columns)
                if missing:
                    raise TransmissionSummaryError(
                        f"{path}: missing columns {sorted(missing)}"
                    )
                continue
            if len(fields) != len(columns):
                raise TransmissionSummaryError(
                    f"{path} line {line_number}: expected {len(columns)} fields, "
                    f"found {len(fields)}"
                )
            current = Cpg(
                chromosome=fields[columns["#chrom"]],
                position=int(fields[columns["start_cpg"]]),
                end=int(fields[columns["end_cpg"]]),
                pat_label=_label(fields[columns["founder_haplotype_pat"]]),
                mat_label=_label(fields[columns["founder_haplotype_mat"]]),
                pat_value=_number(fields[columns["methylation_level_pat_model"]]),
                mat_value=_number(fields[columns["methylation_level_mat_model"]]),
                mismatch=_truth(fields[columns["cpg_is_within_mismatch_window"]]),
            )
            if current.end <= current.position:
                raise TransmissionSummaryError(
                    f"{path}: invalid CpG interval "
                    f"{current.chromosome}:{current.position}-{current.end}"
                )
            current_key = (
                _chromosome_key(current.chromosome),
                current.position,
                current.end,
            )
            if last_key is not None and current_key < last_key:
                raise TransmissionSummaryError(f"{path}: CpGs are not sorted")
            last_key = current_key
            if pending is None:
                pending = current
            elif (pending.chromosome, pending.position, pending.end) == (
                current.chromosome,
                current.position,
                current.end,
            ):
                pending = _merge_duplicate(pending, current, path)
            else:
                yield pending
                pending = current
    if columns is None:
        raise TransmissionSummaryError(f"{path}: no column header")
    if pending is not None:
        yield pending


def _paired_rows(parent_path: Path, child_path: Path) -> Iterator[tuple[Cpg, Cpg]]:
    parent_rows = iter(iter_cpgs(parent_path))
    child_rows = iter(iter_cpgs(child_path))
    parent = next(parent_rows, None)
    child = next(child_rows, None)
    while parent is not None and child is not None:
        parent_key = (
            _chromosome_key(parent.chromosome),
            parent.position,
            parent.end,
        )
        child_key = (
            _chromosome_key(child.chromosome),
            child.position,
            child.end,
        )
        if parent_key < child_key:
            parent = next(parent_rows, None)
        elif child_key < parent_key:
            child = next(child_rows, None)
        else:
            yield parent, child
            parent = next(parent_rows, None)
            child = next(child_rows, None)


def _empty_metrics() -> dict[str, int | float]:
    return {
        "shared_cpgs": 0,
        "eligible_cpgs": 0,
        "mismatch_excluded_cpgs": 0,
        "paired_cpgs": 0,
        "abs_difference_sum": 0.0,
        "difference_sum": 0.0,
        "discordant_cpgs": 0,
        "specificity_cpgs": 0,
        "transmitted_common_abs_sum": 0.0,
        "nontransmitted_abs_sum": 0.0,
    }


def _finalize_metrics(
    values: dict[str, int | float],
    discordance_threshold: float,
    minimum_paired_cpgs: int,
) -> dict[str, int | float | None]:
    eligible = int(values["eligible_cpgs"])
    excluded = int(values["mismatch_excluded_cpgs"])
    evaluated = eligible - excluded
    paired = int(values["paired_cpgs"])
    specificity_n = int(values["specificity_cpgs"])
    mean_abs = float(values["abs_difference_sum"]) / paired if paired else None
    return {
        "shared_cpgs": int(values["shared_cpgs"]),
        "eligible_cpgs": eligible,
        "mismatch_excluded_cpgs": excluded,
        "evaluated_cpgs": evaluated,
        "paired_cpgs": paired,
        "sufficient_cpgs": paired >= minimum_paired_cpgs,
        "callable_fraction": paired / evaluated if evaluated else None,
        "agreement": 1.0 - mean_abs if mean_abs is not None else None,
        "mean_difference": (
            float(values["difference_sum"]) / paired if paired else None
        ),
        "discordant_fraction": (
            int(values["discordant_cpgs"]) / paired if paired else None
        ),
        "discordance_threshold": discordance_threshold,
        "specificity_cpgs": specificity_n,
        "sufficient_specificity_cpgs": specificity_n >= minimum_paired_cpgs,
        "inherited_specificity": (
            (
                float(values["nontransmitted_abs_sum"])
                - float(values["transmitted_common_abs_sum"])
            )
            / specificity_n
            if specificity_n
            else None
        ),
    }


def compare_transmission(
    parent_path: Path,
    child_path: Path,
    relationship: str,
    discordance_threshold: float,
    minimum_paired_cpgs: int,
) -> dict[str, dict[str, int | float | None]]:
    if relationship not in {"paternal", "maternal"}:
        raise ValueError(f"invalid relationship: {relationship}")
    by_chromosome: dict[str, dict[str, int | float]] = {}
    for parent, child in _paired_rows(parent_path, child_path):
        values = by_chromosome.setdefault(child.chromosome, _empty_metrics())
        values["shared_cpgs"] += 1
        transmitted_label = (
            child.pat_label if relationship == "paternal" else child.mat_label
        )
        child_value = (
            child.pat_value if relationship == "paternal" else child.mat_value
        )
        parent_value, nontransmitted_value = parent.value_for_label(transmitted_label)
        if transmitted_label is None or (
            parent.pat_label != transmitted_label
            and parent.mat_label != transmitted_label
        ):
            continue
        values["eligible_cpgs"] += 1
        if parent.mismatch or child.mismatch:
            values["mismatch_excluded_cpgs"] += 1
            continue
        if parent_value is None or child_value is None:
            continue
        difference = child_value - parent_value
        absolute = abs(difference)
        values["paired_cpgs"] += 1
        values["abs_difference_sum"] += absolute
        values["difference_sum"] += difference
        if absolute >= discordance_threshold:
            values["discordant_cpgs"] += 1
        if nontransmitted_value is not None:
            values["specificity_cpgs"] += 1
            values["transmitted_common_abs_sum"] += absolute
            values["nontransmitted_abs_sum"] += abs(
                child_value - nontransmitted_value
            )
    return {
        chromosome: _finalize_metrics(
            values, discordance_threshold, minimum_paired_cpgs
        )
        for chromosome, values in sorted(
            by_chromosome.items(), key=lambda item: _chromosome_key(item[0])
        )
    }


def summarize(
    ped_path: Path,
    all_cpg_paths: list[Path],
    discordance_threshold: float = 0.4,
    minimum_paired_cpgs: int = 100,
) -> dict[str, Any]:
    if not 0.0 <= discordance_threshold <= 1.0:
        raise TransmissionSummaryError("discordance threshold must be in [0, 1]")
    if minimum_paired_cpgs < 1:
        raise TransmissionSummaryError("minimum paired CpGs must be at least 1")
    people = read_pedigree(ped_path)
    sample_paths: dict[str, Path] = {}
    for path in sorted(all_cpg_paths, key=lambda value: value.name):
        sample = sample_id_from_path(path)
        if sample in sample_paths:
            raise TransmissionSummaryError(f"duplicate all-CpG input for {sample!r}")
        if sample not in people:
            raise TransmissionSummaryError(f"all-CpG sample {sample!r} is absent from PED")
        sample_paths[sample] = path

    comparisons: list[dict[str, Any]] = []
    edges: list[dict[str, Any]] = []
    unavailable: list[dict[str, str]] = []
    for child, person in people.items():
        # A QC edge is expected only for a downstream-selected child. PED
        # branches with neither endpoint processed are outside this run, not
        # missing QC data.
        if child not in sample_paths:
            continue
        for relationship, parent_key in (("paternal", "father"), ("maternal", "mother")):
            parent = person[parent_key]
            if parent is None:
                continue
            available = parent in sample_paths
            edge = {
                "family_id": str(person["family"]),
                "parent_id": parent,
                "child_id": child,
                "relationship": relationship,
                "has_methylation_outputs": available,
            }
            edges.append(edge)
            if not available:
                unavailable.append(
                    {
                        **{key: str(value) for key, value in edge.items() if key != "has_methylation_outputs"},
                        "reason": "methylation output is required for both samples",
                    }
                )
                continue
            chromosome_metrics = compare_transmission(
                sample_paths[parent],
                sample_paths[child],
                relationship,
                discordance_threshold,
                minimum_paired_cpgs,
            )
            for chromosome, metrics in chromosome_metrics.items():
                comparisons.append(
                    {
                        "family_id": str(person["family"]),
                        "parent_id": parent,
                        "child_id": child,
                        "relationship": relationship,
                        "chromosome": chromosome,
                        **metrics,
                    }
                )
    return {
        "schema_version": 1,
        "measurement": "model-based methylation",
        "discordance_threshold": discordance_threshold,
        "minimum_paired_cpgs": minimum_paired_cpgs,
        "samples": sorted(sample_paths),
        "edges": edges,
        "comparisons": comparisons,
        "unavailable_edges": unavailable,
        "sources": {
            "ped": ped_path.name,
            "all_cpgs": sorted(path.name for path in all_cpg_paths),
        },
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ped", required=True, type=Path)
    parser.add_argument("--all-cpgs", nargs="*", type=Path, default=[])
    parser.add_argument("--discordance-threshold", type=float, default=0.4)
    parser.add_argument("--minimum-paired-cpgs", type=int, default=100)
    parser.add_argument("--output", required=True, type=Path)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        result = summarize(
            args.ped,
            args.all_cpgs,
            args.discordance_threshold,
            args.minimum_paired_cpgs,
        )
        args.output.write_text(
            json.dumps(result, separators=(",", ":")) + "\n", encoding="utf-8"
        )
    except (TransmissionSummaryError, OSError, ValueError) as exc:
        print(f"transmission methylation summary failed: {exc}", file=sys.stderr)
        return 2
    print(
        f"summarized {len(result['comparisons'])} chromosome-level "
        f"parent-child comparisons"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
