#!/usr/bin/env python3
"""Write the stable, work-directory-free manifest of published Tapestry results."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any


def _artifact(path: str, index: str | None = None) -> dict[str, Any]:
    result: dict[str, Any] = {"path": path, "status": "complete"}
    if index is not None:
        result["index"] = index
    return result


def build_results_manifest(
    resolved_run: Path,
    selected_samples: Path,
    sample_qc_paths: list[Path],
    output_path: Path,
) -> dict[str, Any]:
    config = json.loads(resolved_run.read_text(encoding="utf-8"))
    sample_ids = [
        line.split("\t", 1)[0]
        for line in selected_samples.read_text(encoding="utf-8").splitlines()[1:]
        if line.strip()
    ]
    sample_statuses: dict[str, str] = {}
    for qc_path in sample_qc_paths:
        qc = json.loads(qc_path.read_text(encoding="utf-8"))
        sample_id = qc.get("sample_id")
        status = qc.get("founder_phasing_status")
        if sample_id not in sample_ids:
            raise ValueError(f"unexpected sample QC for {sample_id!r}: {qc_path}")
        if sample_id in sample_statuses:
            raise ValueError(f"duplicate sample QC for {sample_id!r}")
        if status not in {"complete", "no_inheritance_phase"}:
            raise ValueError(f"invalid founder phasing status for {sample_id!r}: {status!r}")
        sample_statuses[sample_id] = status
    missing_qc = sorted(set(sample_ids) - set(sample_statuses))
    if missing_qc:
        raise ValueError(f"missing sample QC for: {', '.join(missing_qc)}")
    project_id = config["project"]["id"]
    bigwig_enabled = config["outputs"]["bigwig"]
    reference_name = config["reference"]["name"]
    reference_bed = f"reference/{reference_name}.autosomes.cpgs.bed.gz"
    inheritance_prefix = f"inheritance/{project_id}"
    samples: dict[str, Any] = {}
    for sample_id in sample_ids:
        root = f"samples/{sample_id}"
        model_root = f"{root}/model_inputs/{sample_id}.model"
        founder_phasing: dict[str, Any] = {
            "table": _artifact(f"{root}/{sample_id}.dna-methylation.bed"),
            "table_header": _artifact(
                f"{root}/{sample_id}.dna-methylation.bed.header"
            ),
            "hap_map": _artifact(f"{root}/{sample_id}.hap-map-blocks.bed"),
            "hap_map_header": _artifact(
                f"{root}/{sample_id}.hap-map-blocks.bed.header"
            ),
            "paternal_hap_map": _artifact(
                f"{root}/{sample_id}.hap-map-blocks.paternal.bed.gz",
                f"{root}/{sample_id}.hap-map-blocks.paternal.bed.gz.tbi",
            ),
            "paternal_hap_map_header": _artifact(
                f"{root}/{sample_id}.hap-map-blocks.paternal.bed.header"
            ),
            "maternal_hap_map": _artifact(
                f"{root}/{sample_id}.hap-map-blocks.maternal.bed.gz",
                f"{root}/{sample_id}.hap-map-blocks.maternal.bed.gz.tbi",
            ),
            "maternal_hap_map_header": _artifact(
                f"{root}/{sample_id}.hap-map-blocks.maternal.bed.header"
            ),
            "mismatch_bed": _artifact(
                f"{root}/{sample_id}.bit-vector-sites-mismatches.bed"
            ),
            "mismatch_bed_header": _artifact(
                f"{root}/{sample_id}.bit-vector-sites-mismatches.bed.header"
            ),
            "mismatch_vcf": _artifact(
                f"{root}/{sample_id}.bit-vector-sites-mismatches.vcf.gz",
                f"{root}/{sample_id}.bit-vector-sites-mismatches.vcf.gz.tbi",
            ),
            "qc": _artifact(f"{root}/{sample_id}.phasing-qc.json"),
        }
        if bigwig_enabled:
            founder_phasing.update(
                {
                    "paternal_bigwig": _artifact(
                        f"{root}/{sample_id}.dna-methylation.pat.model.{reference_name}.bw"
                    ),
                    "maternal_bigwig": _artifact(
                        f"{root}/{sample_id}.dna-methylation.mat.model.{reference_name}.bw"
                    ),
                }
            )
        else:
            founder_phasing["bigwig"] = {"status": "disabled"}

        samples[sample_id] = {
            "status": sample_statuses[sample_id],
            "enabled_modes": ["model"],
            "disabled_modes": ["count"],
            "model_inputs": {
                label: _artifact(
                    f"{model_root}.{label}.bed.gz",
                    f"{model_root}.{label}.bed.gz.tbi",
                )
                for label in ("combined", "hap1", "hap2")
            },
            "model_filter_qc": _artifact(
                f"{root}/model_inputs/{sample_id}.model-filter-qc.json"
            ),
            "founder_phasing": founder_phasing,
            "all_cpgs": {
                "table": _artifact(
                    f"{root}/{sample_id}.dna-methylation.all-cpgs.bed.gz",
                    f"{root}/{sample_id}.dna-methylation.all-cpgs.bed.gz.tbi",
                ),
                "qc": _artifact(f"{root}/{sample_id}.all-cpgs-qc.json"),
            },
        }

    result = {
        "schema_version": 1,
        "project_id": project_id,
        "status": "complete",
        "config_fingerprint": config["_tapestry"]["config_fingerprint"],
        "pipeline_info": {
            filename: _artifact(f"pipeline_info/{filename}")
            for filename in (
                "resolved-run.json",
                "resolved-manifest.json",
                "normalized.ped",
                "selected-samples.tsv",
                "selected-artifacts.json",
                "validation-report.json",
                "config-fingerprint.txt",
                "validation.success",
                "versions.json",
            )
        },
        "reference": {
            "name": reference_name,
            "regions": config["_tapestry"]["selected_autosomes"],
            "cpgs": _artifact(reference_bed, f"{reference_bed}.tbi"),
            "qc": _artifact("reference/reference-cpgs-qc.json"),
        },
        "inheritance": {
            "all_sites_vcf": _artifact(
                f"{inheritance_prefix}.all-sites.vcf.gz",
                f"{inheritance_prefix}.all-sites.vcf.gz.tbi",
            ),
            "complete_family_map_vcf": _artifact(
                f"{inheritance_prefix}.complete-family-map.vcf.gz",
                f"{inheritance_prefix}.complete-family-map.vcf.gz.tbi",
            ),
            "iht": _artifact(f"{inheritance_prefix}.iht.sorted.txt"),
            "markers": _artifact(f"{inheritance_prefix}.markers.sorted.txt"),
            "recombinants": _artifact(f"{inheritance_prefix}.recombinants.sorted.txt"),
            "pass_vcf": _artifact(
                f"{inheritance_prefix}.pass.vcf.gz",
                f"{inheritance_prefix}.pass.vcf.gz.tbi",
            ),
            "fail_vcf": _artifact(
                f"{inheritance_prefix}.fail.vcf.gz",
                f"{inheritance_prefix}.fail.vcf.gz.tbi",
            ),
            "qc": _artifact(f"{inheritance_prefix}.inheritance-qc.json"),
            "normalization_qc": _artifact(
                f"{inheritance_prefix}.vcf-normalization.json"
            ),
            "filtering_stats": _artifact(
                f"{inheritance_prefix}.filtering-stats.txt"
            ),
            "map_log": _artifact(f"{inheritance_prefix}.gtg-ped-map.log"),
            "concordance_log": _artifact(
                f"{inheritance_prefix}.gtg-concordance.log"
            ),
        },
        "selected_samples": sample_ids,
        "samples": samples,
        "disabled_output_modes": ["count"],
    }
    output_path.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return result


def format_completion_summary(result: dict[str, Any], output_dir: Path) -> str:
    """Return a concise human-readable summary of a completed run."""
    statuses = [sample["status"] for sample in result["samples"].values()]
    complete = statuses.count("complete")
    no_phase = statuses.count("no_inheritance_phase")
    return "\n".join(
        [
            "Tapestry completed",
            f"  Results manifest: {output_dir / 'results-manifest.json'}",
            f"  Samples: {complete} complete, {no_phase} without inheritance phase",
        ]
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--resolved-run", required=True, type=Path)
    parser.add_argument("--selected-samples", required=True, type=Path)
    parser.add_argument("--sample-qc", required=True, action="append", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args(argv)
    result = build_results_manifest(
        args.resolved_run, args.selected_samples, args.sample_qc, args.output
    )
    config = json.loads(args.resolved_run.read_text(encoding="utf-8"))
    print(format_completion_summary(result, Path(config["project"]["outdir"])))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
