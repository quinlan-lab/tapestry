#!/usr/bin/env python3
"""Generate one indexed reference-CpG BED for the configured autosomes."""

from __future__ import annotations

import argparse
import json
import re
import sys
import tempfile
from pathlib import Path
from typing import Any

import pysam


class ReferenceCpgError(RuntimeError):
    pass


def generate_reference_cpgs(
    fasta: Path,
    fasta_index: Path,
    regions: list[str],
    output_bed: Path,
) -> dict[str, Any]:
    for label, path in (("FASTA", fasta), ("FASTA index", fasta_index)):
        if not path.is_file():
            raise ReferenceCpgError(f"{label} does not exist: {path}")
    if Path(str(fasta) + ".fai").resolve() != fasta_index.resolve():
        raise ReferenceCpgError(
            "FASTA index must use the conventional <fasta>.fai name for pysam"
        )

    output_bed.parent.mkdir(parents=True, exist_ok=True)
    counts: dict[str, int] = {}
    try:
        with pysam.FastaFile(str(fasta)) as reference:
            missing = [region for region in regions if region not in reference.references]
            if missing:
                raise ReferenceCpgError(f"FASTA is missing configured regions: {missing}")
            with tempfile.TemporaryDirectory(
                prefix=".reference-cpgs-", dir=output_bed.parent
            ) as temporary:
                temporary_dir = Path(temporary)
                plain = temporary_dir / "reference-cpgs.bed"
                compressed = temporary_dir / "reference-cpgs.bed.gz"
                with plain.open("w", encoding="utf-8") as destination:
                    destination.write("##reference-cpg-definition=cytosine-base-of-CG\n")
                    destination.write("#chrom\tstart\tend\n")
                    for region in regions:
                        count = 0
                        sequence = reference.fetch(region).upper()
                        for match in re.finditer("CG", sequence):
                            start = match.start()
                            destination.write(f"{region}\t{start}\t{start + 1}\n")
                            count += 1
                        counts[region] = count
                pysam.tabix_compress(str(plain), str(compressed), force=True)
                pysam.tabix_index(str(compressed), preset="bed", force=True)
                compressed.replace(output_bed)
                Path(str(compressed) + ".tbi").replace(
                    Path(str(output_bed) + ".tbi")
                )
    except (OSError, ValueError) as exc:
        raise ReferenceCpgError(f"cannot generate reference CpGs: {exc}") from exc

    return {
        "fasta": fasta.name,
        "regions": regions,
        "records_by_region": counts,
        "records": sum(counts.values()),
        "output": output_bed.name,
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", required=True, type=Path)
    parser.add_argument("--fasta-index", required=True, type=Path)
    parser.add_argument("--regions", required=True)
    parser.add_argument("--reference-name", required=True)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    output = args.output_dir / f"{args.reference_name}.autosomes.cpgs.bed.gz"
    try:
        report = generate_reference_cpgs(
            args.fasta,
            args.fasta_index,
            args.regions.split(","),
            output,
        )
    except ReferenceCpgError as exc:
        print(f"reference CpG generation failed: {exc}", file=sys.stderr)
        return 2
    (args.output_dir / "reference-cpgs-qc.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"generated {report['records']} reference CpGs")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
