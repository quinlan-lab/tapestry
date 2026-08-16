#!/usr/bin/env python3
"""Run pinned gtg map/concordance and publish deterministic inheritance outputs."""

from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any

import pysam


class InheritanceError(RuntimeError):
    pass


def _run_logged(command: list[str], log_path: Path) -> None:
    try:
        with log_path.open("w", encoding="utf-8") as log:
            completed = subprocess.run(
                command,
                check=False,
                text=True,
                stdout=log,
                stderr=subprocess.STDOUT,
            )
    except OSError as exc:
        raise InheritanceError(f"cannot execute {command[0]!r}: {exc}") from exc
    if completed.returncode != 0:
        raise InheritanceError(
            f"command failed with exit {completed.returncode}: {command!r}; see {log_path}"
        )


def _run(command: list[str]) -> None:
    try:
        completed = subprocess.run(command, check=False, text=True, capture_output=True)
    except OSError as exc:
        raise InheritanceError(f"cannot execute {command[0]!r}: {exc}") from exc
    if completed.returncode != 0:
        raise InheritanceError(
            f"command failed with exit {completed.returncode}: {command!r}\n"
            f"{completed.stderr.strip()}"
        )


def _version(executable: str) -> str:
    try:
        completed = subprocess.run(
            [executable, "--version"], check=True, text=True, capture_output=True
        )
    except (OSError, subprocess.CalledProcessError) as exc:
        raise InheritanceError(f"cannot obtain version from {executable!r}: {exc}") from exc
    return completed.stdout.strip().splitlines()[0]


def _natural(value: str) -> tuple[Any, ...]:
    return tuple(
        int(part) if part.isdigit() else part.lower()
        for part in re.split(r"(\d+)", value)
    )


def _sort_gtg_table(source: Path, destination: Path) -> int:
    lines = source.read_text(encoding="utf-8").splitlines()
    if not lines:
        raise InheritanceError(f"gtg output is empty: {source}")
    header, records = lines[0], [line for line in lines[1:] if line.strip()]

    def key(line: str) -> tuple[Any, ...]:
        fields = line.split()
        if len(fields) < 2:
            raise InheritanceError(f"malformed gtg output line in {source}: {line!r}")
        try:
            position = int(fields[1])
        except ValueError:
            position = 0
        return (_natural(fields[0]), position, line)

    records.sort(key=key)
    destination.write_text("\n".join([header, *records]) + "\n", encoding="utf-8")
    return len(records)


def _sort_and_index_vcf(
    bcftools: str, source: Path, destination: Path
) -> dict[str, int]:
    _run(
        [
            bcftools,
            "sort",
            "--output-type",
            "z",
            "--output",
            str(destination),
            str(source),
        ]
    )
    _run([bcftools, "index", "--tbi", "--force", str(destination)])
    records = 0
    with pysam.VariantFile(str(destination)) as vcf:
        for _ in vcf:
            records += 1
    return {"records": records}


def _replace(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    source.replace(destination)


def run_inheritance(
    resolved_run: Path,
    normalized_ped: Path,
    map_vcf: Path,
    all_sites_vcf: Path,
    normalization_report: Path,
    output_dir: Path,
    gtg_ped_map: str = "gtg-ped-map",
    gtg_concordance: str = "gtg-concordance",
    bcftools: str = "bcftools",
) -> dict[str, Any]:
    for label, path in (
        ("resolved run", resolved_run),
        ("normalized PED", normalized_ped),
        ("complete-family map VCF", map_vcf),
        ("all-site VCF", all_sites_vcf),
        ("normalization report", normalization_report),
    ):
        if not path.is_file():
            raise InheritanceError(f"{label} does not exist: {path}")

    executables: dict[str, str] = {}
    for label, requested in (
        ("gtg-ped-map", gtg_ped_map),
        ("gtg-concordance", gtg_concordance),
        ("bcftools", bcftools),
    ):
        executable = shutil.which(requested)
        if executable is None:
            raise InheritanceError(f"required executable not found: {requested!r}")
        executables[label] = executable

    try:
        config = json.loads(resolved_run.read_text(encoding="utf-8"))
        normalization = json.loads(normalization_report.read_text(encoding="utf-8"))
        prefix = config["project"]["id"]
        map_parameters = config["inheritance"]["map"]
        concordance_parameters = config["inheritance"]["concordance"]
    except (OSError, KeyError, json.JSONDecodeError) as exc:
        raise InheritanceError(f"cannot read normalized contracts: {exc}") from exc

    output_dir.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".gtg-", dir=output_dir) as temporary:
        temporary_dir = Path(temporary)
        gtg_prefix = temporary_dir / prefix
        map_log = temporary_dir / f"{prefix}.gtg-ped-map.log"
        concordance_log = temporary_dir / f"{prefix}.gtg-concordance.log"

        _run_logged(
            [
                executables["gtg-ped-map"],
                "--ped",
                str(normalized_ped),
                "--vcf",
                str(map_vcf),
                "--prefix",
                str(gtg_prefix),
                "--qual",
                str(map_parameters["min_qual"]),
                "--depth",
                str(map_parameters["min_depth"]),
                "--run",
                str(map_parameters["min_run_markers"]),
                "--verbose",
            ],
            map_log,
        )

        iht_raw = Path(str(gtg_prefix) + ".iht.txt")
        markers_raw = Path(str(gtg_prefix) + ".markers.txt")
        recombinants_raw = Path(str(gtg_prefix) + ".recombinants.txt")
        for path in (iht_raw, markers_raw, recombinants_raw):
            if not path.is_file():
                raise InheritanceError(f"gtg-ped-map did not create expected output: {path}")

        iht_sorted = temporary_dir / f"{prefix}.iht.sorted.txt"
        markers_sorted = temporary_dir / f"{prefix}.markers.sorted.txt"
        recombinants_sorted = temporary_dir / f"{prefix}.recombinants.sorted.txt"
        iht_records = _sort_gtg_table(iht_raw, iht_sorted)
        marker_records = _sort_gtg_table(markers_raw, markers_sorted)
        recombinant_records = _sort_gtg_table(recombinants_raw, recombinants_sorted)
        if iht_records == 0:
            raise InheritanceError("gtg-ped-map produced no inheritance blocks")

        _run_logged(
            [
                executables["gtg-concordance"],
                "--ped",
                str(normalized_ped),
                "--inheritance",
                str(iht_sorted),
                "--vcf",
                str(all_sites_vcf),
                "--prefix",
                str(gtg_prefix),
                "--qual",
                str(concordance_parameters["min_qual"]),
                "--depth",
                str(concordance_parameters["min_depth"]),
                "--verbose",
            ],
            concordance_log,
        )

        pass_raw = Path(str(gtg_prefix) + ".pass.vcf")
        fail_raw = Path(str(gtg_prefix) + ".fail.vcf")
        filtering_stats_raw = Path(str(gtg_prefix) + ".filtering_stats.txt")
        for path in (pass_raw, fail_raw, filtering_stats_raw):
            if not path.is_file():
                raise InheritanceError(
                    f"gtg-concordance did not create expected output: {path}"
                )

        pass_sorted = temporary_dir / f"{prefix}.pass.vcf.gz"
        fail_sorted = temporary_dir / f"{prefix}.fail.vcf.gz"
        pass_stats = _sort_and_index_vcf(executables["bcftools"], pass_raw, pass_sorted)
        fail_stats = _sort_and_index_vcf(executables["bcftools"], fail_raw, fail_sorted)

        qc = {
            "family_id": normalization.get("family_id", prefix),
            "parameters": {
                "map": map_parameters,
                "concordance": concordance_parameters,
            },
            "normalization": {
                "all_site_records": normalization["all_sites"]["records"],
                "complete_family_map_records": normalization["complete_family_map"]["records"],
            },
            "sample_stats": normalization["all_sites"]["sample_stats"],
            "inheritance": {
                "blocks": iht_records,
                "markers": marker_records,
                "recombinant_rows": recombinant_records,
                "pass_records": pass_stats["records"],
                "fail_records": fail_stats["records"],
            },
            "versions": {
                "gtg_ped_map": _version(executables["gtg-ped-map"]),
                "gtg_concordance": _version(executables["gtg-concordance"]),
                "bcftools": _version(executables["bcftools"]),
            },
        }
        qc_path = temporary_dir / f"{prefix}.inheritance-qc.json"
        qc_path.write_text(json.dumps(qc, indent=2, sort_keys=True) + "\n", encoding="utf-8")

        destinations = {
            iht_sorted: output_dir / iht_sorted.name,
            markers_sorted: output_dir / markers_sorted.name,
            recombinants_sorted: output_dir / recombinants_sorted.name,
            pass_sorted: output_dir / pass_sorted.name,
            Path(str(pass_sorted) + ".tbi"): output_dir / f"{pass_sorted.name}.tbi",
            fail_sorted: output_dir / fail_sorted.name,
            Path(str(fail_sorted) + ".tbi"): output_dir / f"{fail_sorted.name}.tbi",
            filtering_stats_raw: output_dir / f"{prefix}.filtering-stats.txt",
            map_log: output_dir / map_log.name,
            concordance_log: output_dir / concordance_log.name,
            qc_path: output_dir / qc_path.name,
        }
        for source, destination in destinations.items():
            _replace(source, destination)
    return qc


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--resolved-run", required=True, type=Path)
    parser.add_argument("--normalized-ped", required=True, type=Path)
    parser.add_argument("--map-vcf", required=True, type=Path)
    parser.add_argument("--all-sites-vcf", required=True, type=Path)
    parser.add_argument("--normalization-report", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--gtg-ped-map", default="gtg-ped-map")
    parser.add_argument("--gtg-concordance", default="gtg-concordance")
    parser.add_argument("--bcftools", default="bcftools")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        qc = run_inheritance(
            args.resolved_run,
            args.normalized_ped,
            args.map_vcf,
            args.all_sites_vcf,
            args.normalization_report,
            args.output_dir,
            args.gtg_ped_map,
            args.gtg_concordance,
            args.bcftools,
        )
    except InheritanceError as exc:
        print(f"inheritance failed: {exc}", file=sys.stderr)
        return 2
    print(
        f"created {qc['inheritance']['blocks']} inheritance blocks and "
        f"{qc['inheritance']['pass_records']} concordant variants"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
