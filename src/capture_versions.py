#!/usr/bin/env python3
"""Capture runtime, upstream, and container provenance for a Tapestry run."""

from __future__ import annotations

import argparse
import importlib.metadata
import json
import platform
import subprocess
from pathlib import Path
from typing import Any


GTG_COMMIT = "e12aca6b49ee7208952467db4a2a9e2f79b98efb"
BEDGRAPH_TO_BIGWIG_SHA256 = (
    "43c1f8f2aecf2647dc02ec15a8b8c43af9c2639ab850ffbb717e0c5cb18634da"
)


def _command_version(command: list[str], *, check: bool = True) -> str:
    completed = subprocess.run(command, check=check, text=True, capture_output=True)
    output = completed.stdout or completed.stderr
    return output.strip().splitlines()[0]


def capture_versions(
    resolved_manifest: Path, container: str, nextflow_version: str
) -> dict[str, Any]:
    manifest = json.loads(resolved_manifest.read_text(encoding="utf-8"))
    packages = {}
    for package in (
        "bioframe",
        "cyvcf2",
        "numpy",
        "pandas",
        "polars",
        "pyBigWig",
        "pyarrow",
        "pyfaidx",
        "pysam",
        "PyYAML",
    ):
        packages[package] = importlib.metadata.version(package)
    digest = container.split("@", 1)[1] if "@" in container else None
    return {
        "tapestry": "0.1.0-dev",
        "python": platform.python_version(),
        "nextflow": nextflow_version,
        "container": {"reference": container, "digest": digest},
        "executables": {
            "bcftools": _command_version(["bcftools", "--version"]),
            "gtg_ped_map": _command_version(["gtg-ped-map", "--version"]),
            "gtg_concordance": _command_version(["gtg-concordance", "--version"]),
            "gtg_commit": GTG_COMMIT,
            "bedgraph_to_bigwig": _command_version(
                ["bedGraphToBigWig"], check=False
            ),
            "bedgraph_to_bigwig_sha256": BEDGRAPH_TO_BIGWIG_SHA256,
        },
        "python_packages": packages,
        "upstream_workflow": manifest["workflow"],
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--resolved-manifest", required=True, type=Path)
    parser.add_argument("--container", required=True)
    parser.add_argument("--nextflow-version", required=True)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args(argv)
    result = capture_versions(
        args.resolved_manifest, args.container, args.nextflow_version
    )
    args.output.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"captured runtime versions for {result['container']['reference']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
