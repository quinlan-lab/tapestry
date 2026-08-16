#!/usr/bin/env python3
"""Create all-site and complete-family VCF views for pedigree inheritance."""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any

import pysam


class NormalizationError(RuntimeError):
    pass


def _usable_depth(sample: Any) -> tuple[str, float] | None:
    """Return the first usable gtg-compatible depth, preferring DP."""
    for field in ("DP", "AD", "SD"):
        value = sample.get(field)
        values = value if isinstance(value, (tuple, list)) else (value,)
        numeric = [float(item) for item in values if isinstance(item, (int, float))]
        if numeric:
            return field, sum(numeric)
    return None


def _run(command: list[str]) -> None:
    try:
        completed = subprocess.run(command, check=False, text=True, capture_output=True)
    except OSError as exc:
        raise NormalizationError(f"cannot execute {command[0]!r}: {exc}") from exc
    if completed.returncode != 0:
        stderr = completed.stderr.strip()
        raise NormalizationError(
            f"command failed with exit {completed.returncode}: {command!r}"
            + (f"\n{stderr}" if stderr else "")
        )


def _ped_samples(path: Path) -> list[str]:
    samples: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            fields = stripped.split()
            if len(fields) != 6:
                raise NormalizationError(
                    f"{path} line {line_number}: expected six PED fields"
                )
            samples.append(fields[1])
    if not samples:
        raise NormalizationError(f"{path}: PED contains no samples")
    if len(samples) != len(set(samples)):
        raise NormalizationError(f"{path}: PED contains duplicate sample IDs")
    return samples


def _vcf_stats(path: Path) -> dict[str, Any]:
    stats: dict[str, Any] = {
        "records": 0,
        "non_pass_records": 0,
        "multiallelic_records": 0,
        "complete_family_records": 0,
        "samples": [],
        "format_fields": [],
    }
    try:
        vcf = pysam.VariantFile(str(path))
    except (OSError, ValueError) as exc:
        raise NormalizationError(f"cannot open normalized VCF {path}: {exc}") from exc
    with vcf:
        samples = list(vcf.header.samples)
        stats["samples"] = samples
        stats["format_fields"] = sorted(vcf.header.formats)
        sample_accumulators = {
            sample: {
                "called_genotypes": 0,
                "missing_genotypes": 0,
                "depth_observations": 0,
                "depth_sum": 0.0,
                "depth_min": None,
                "depth_max": None,
                "depth_source_records": {"DP": 0, "AD": 0, "SD": 0},
            }
            for sample in samples
        }
        for record in vcf:
            stats["records"] += 1
            if record.filter.keys() and "PASS" not in record.filter.keys():
                stats["non_pass_records"] += 1
            if len(record.alts or ()) > 1:
                stats["multiallelic_records"] += 1
            complete = True
            for sample in samples:
                genotype = record.samples[sample].get("GT")
                if genotype is None or any(allele is None for allele in genotype):
                    sample_accumulators[sample]["missing_genotypes"] += 1
                    complete = False
                    continue
                accumulator = sample_accumulators[sample]
                accumulator["called_genotypes"] += 1
                depth = _usable_depth(record.samples[sample])
                if depth is not None:
                    source, value = depth
                    accumulator["depth_observations"] += 1
                    accumulator["depth_sum"] += value
                    accumulator["depth_source_records"][source] += 1
                    current_min = accumulator["depth_min"]
                    current_max = accumulator["depth_max"]
                    accumulator["depth_min"] = (
                        value if current_min is None else min(current_min, value)
                    )
                    accumulator["depth_max"] = (
                        value if current_max is None else max(current_max, value)
                    )
            if complete:
                stats["complete_family_records"] += 1
        stats["sample_stats"] = {}
        for sample, accumulator in sample_accumulators.items():
            observations = accumulator["depth_observations"]
            stats["sample_stats"][sample] = {
                "called_genotypes": accumulator["called_genotypes"],
                "missing_genotypes": accumulator["missing_genotypes"],
                "called_genotype_depth": {
                    "observations": observations,
                    "missing": accumulator["called_genotypes"] - observations,
                    "min": accumulator["depth_min"],
                    "max": accumulator["depth_max"],
                    "mean": (
                        round(accumulator["depth_sum"] / observations, 6)
                        if observations
                        else None
                    ),
                    "source_records": accumulator["depth_source_records"],
                },
            }
    return stats


def _atomic_move_indexed_vcf(temporary_vcf: Path, destination_vcf: Path) -> None:
    temporary_index = Path(str(temporary_vcf) + ".tbi")
    destination_index = Path(str(destination_vcf) + ".tbi")
    destination_vcf.parent.mkdir(parents=True, exist_ok=True)
    temporary_vcf.replace(destination_vcf)
    temporary_index.replace(destination_index)


def normalize_joint_vcf(
    input_vcf: Path,
    input_index: Path,
    ped: Path,
    regions: list[str],
    output_dir: Path,
    prefix: str,
    bcftools: str = "bcftools",
) -> dict[str, Any]:
    for label, path in (
        ("input VCF", input_vcf),
        ("input VCF index", input_index),
        ("PED", ped),
    ):
        if not path.is_file():
            raise NormalizationError(f"{label} does not exist: {path}")
    executable = shutil.which(bcftools)
    if executable is None:
        raise NormalizationError(f"bcftools executable not found: {bcftools!r}")
    if not regions:
        raise NormalizationError("at least one autosome is required")

    samples = _ped_samples(ped)
    try:
        with pysam.VariantFile(str(input_vcf)) as source:
            input_samples = list(source.header.samples)
            if set(input_samples) != set(samples) or len(input_samples) != len(samples):
                raise NormalizationError(
                    f"input VCF sample set {input_samples} does not equal PED sample set {samples}"
                )
            removable = [
                field for field in ("PS", "PF") if field in source.header.formats
            ]
    except (OSError, ValueError) as exc:
        if isinstance(exc, NormalizationError):
            raise
        raise NormalizationError(f"cannot inspect input VCF {input_vcf}: {exc}") from exc

    output_dir.mkdir(parents=True, exist_ok=True)
    all_site_output = output_dir / f"{prefix}.all-sites.vcf.gz"
    map_output = output_dir / f"{prefix}.complete-family-map.vcf.gz"
    report_output = output_dir / f"{prefix}.vcf-normalization.json"

    with tempfile.TemporaryDirectory(prefix=".normalize-vcf-", dir=output_dir) as temporary:
        temporary_dir = Path(temporary)
        stripped_bcf = temporary_dir / "stripped.bcf"
        all_site_temporary = temporary_dir / "all-sites.vcf.gz"
        map_temporary = temporary_dir / "complete-family-map.vcf.gz"

        if removable:
            _run(
                [
                    executable,
                    "annotate",
                    "--remove",
                    ",".join(f"FORMAT/{field}" for field in removable),
                    "--output-type",
                    "b",
                    "--output",
                    str(stripped_bcf),
                    str(input_vcf),
                ]
            )
        else:
            _run(
                [
                    executable,
                    "view",
                    "--output-type",
                    "b",
                    "--output",
                    str(stripped_bcf),
                    str(input_vcf),
                ]
            )
        _run([executable, "index", "--force", str(stripped_bcf)])

        _run(
            [
                executable,
                "view",
                "--regions",
                ",".join(regions),
                "--samples",
                ",".join(samples),
                "--no-update",
                "--output-type",
                "z",
                "--output",
                str(all_site_temporary),
                str(stripped_bcf),
            ]
        )
        _run([executable, "index", "--tbi", "--force", str(all_site_temporary)])

        _run(
            [
                executable,
                "view",
                "--include",
                "N_MISSING=0",
                "--no-update",
                "--output-type",
                "z",
                "--output",
                str(map_temporary),
                str(all_site_temporary),
            ]
        )
        _run([executable, "index", "--tbi", "--force", str(map_temporary)])

        all_site_stats = _vcf_stats(all_site_temporary)
        map_stats = _vcf_stats(map_temporary)
        if all_site_stats["samples"] != samples:
            raise NormalizationError(
                "normalized all-site VCF sample order does not match normalized PED order"
            )
        if "PS" in all_site_stats["format_fields"] or "PF" in all_site_stats["format_fields"]:
            raise NormalizationError("normalized all-site VCF still contains FORMAT/PS or FORMAT/PF")
        if map_stats["records"] != map_stats["complete_family_records"]:
            raise NormalizationError("map-only VCF contains an incomplete family genotype")

        _atomic_move_indexed_vcf(all_site_temporary, all_site_output)
        _atomic_move_indexed_vcf(map_temporary, map_output)

    report = {
        "input_vcf": str(input_vcf.resolve()),
        "input_index": str(input_index.resolve()),
        "ped": ped.name,
        "regions": regions,
        "samples": samples,
        "removed_format_fields": removable,
        "all_sites": {
            "vcf": all_site_output.name,
            "index": f"{all_site_output.name}.tbi",
            **all_site_stats,
        },
        "complete_family_map": {
            "vcf": map_output.name,
            "index": f"{map_output.name}.tbi",
            **map_stats,
        },
        "bcftools_version": _bcftools_version(executable),
    }
    temporary_report = report_output.with_name(f".{report_output.name}.tmp")
    temporary_report.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    temporary_report.replace(report_output)
    return report


def _bcftools_version(executable: str) -> str:
    try:
        completed = subprocess.run(
            [executable, "--version"], check=True, text=True, capture_output=True
        )
    except (OSError, subprocess.CalledProcessError) as exc:
        raise NormalizationError(f"cannot obtain bcftools version: {exc}") from exc
    return completed.stdout.splitlines()[0]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Normalize the WDL joint VCF for concordance/annotation and create a "
            "complete-family map-only VCF."
        )
    )
    parser.add_argument("--resolved-run", required=True, type=Path)
    parser.add_argument("--resolved-manifest", required=True, type=Path)
    parser.add_argument("--ped", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--bcftools", default="bcftools")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        run_config = json.loads(args.resolved_run.read_text(encoding="utf-8"))
        manifest = json.loads(args.resolved_manifest.read_text(encoding="utf-8"))
        joint = manifest["joint_small_variants"]
        regions = run_config["_tapestry"]["selected_autosomes"]
        prefix = run_config["project"]["id"]
        report = normalize_joint_vcf(
            Path(joint["vcf"]),
            Path(joint["index"]),
            args.ped,
            regions,
            args.output_dir,
            prefix,
            args.bcftools,
        )
    except (OSError, KeyError, json.JSONDecodeError, NormalizationError) as exc:
        print(f"joint VCF normalization failed: {exc}", file=sys.stderr)
        return 2
    print(
        f"normalized {report['all_sites']['records']} all-site records; "
        f"{report['complete_family_map']['records']} complete-family records"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
