#!/usr/bin/env python3
"""Validate and normalize a generic Tapestry pedigree/model-only run."""

from __future__ import annotations

import argparse
import copy
import gzip
import hashlib
import json
import re
import sys
from pathlib import Path
from typing import Any, TextIO

import pysam


VALIDATOR_VERSION = "0.2.0"
PIPELINE_VERSION = "0.1.0-dev"
GTG_COMMIT = "e12aca6b49ee7208952467db4a2a9e2f79b98efb"
RUNTIME_LOCK_VERSION = 1
SUPPORTED_WDL_RELEASES = {
    "v3.3.0": "db06f0af2354d847b971b0548eaade9ff145c912",
    "v3.3.1": "477ef39ad69e86e90897ea7e313b86bfc12a2a96",
}
SAMPLE_OUTPUTS = (
    "phased_small_variant_vcf",
    "phased_small_variant_vcf_index",
    "phase_blocks",
    "cpg_combined_bed",
    "cpg_combined_bed_index",
    "cpg_hap1_bed",
    "cpg_hap1_bed_index",
    "cpg_hap2_bed",
    "cpg_hap2_bed_index",
)
AUTOSOMES = tuple(f"chr{i}" for i in range(1, 23))
MISSING_PARENTS = frozenset({"0", ".", "NA"})
PHASE_BLOCK_COLUMNS = frozenset(
    {
        "source_block_index",
        "sample_name",
        "phase_block_id",
        "chrom",
        "start",
        "end",
        "num_variants",
    }
)
MODEL_BED_COLUMNS = frozenset(
    {"chrom", "begin", "end", "mod_score", "type", "cov"}
)


class InputValidationError(ValueError):
    """An actionable error in a user-provided run contract or artifact."""


def _unique_json_object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise InputValidationError(f"duplicate JSON key: {key!r}")
        result[key] = value
    return result


def _load_json(path: Path, label: str) -> dict[str, Any]:
    try:
        text = path.read_text(encoding="utf-8")
    except OSError as exc:
        raise InputValidationError(f"cannot read {label} {path}: {exc}") from exc

    try:
        value = json.loads(text, object_pairs_hook=_unique_json_object)
    except InputValidationError:
        raise
    except json.JSONDecodeError as exc:
        raise InputValidationError(f"cannot parse {label} {path}: {exc}") from exc

    if not isinstance(value, dict):
        raise InputValidationError(f"{label} must contain a JSON object")
    return value


def _require_file(path_value: str, label: str) -> Path:
    path = Path(path_value)
    if not path.is_file():
        raise InputValidationError(f"{label}: file does not exist: {path}")
    try:
        with path.open("rb") as handle:
            handle.read(1)
    except OSError as exc:
        raise InputValidationError(f"{label}: file is not readable: {path}: {exc}") from exc
    return path


def _normalize_parent(value: str) -> str | None:
    return None if value.upper() in MISSING_PARENTS else value


def parse_ped(path: Path) -> list[dict[str, str | None]]:
    records: list[dict[str, str | None]] = []
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            fields = stripped.split()
            if len(fields) != 6:
                raise InputValidationError(
                    f"pedigree.ped line {line_number}: expected 6 fields, found {len(fields)}"
                )
            family, sample, father, mother, sex, phenotype = fields
            if sex not in {"0", "1", "2"}:
                raise InputValidationError(
                    f"pedigree.ped line {line_number}: sex must be 0, 1, or 2; found {sex!r}"
                )
            father_value = _normalize_parent(father)
            mother_value = _normalize_parent(mother)
            if father_value is not None and father_value == mother_value:
                raise InputValidationError(
                    f"pedigree.ped line {line_number}: father and mother are both {father_value!r}"
                )
            records.append(
                {
                    "family_id": family,
                    "sample_id": sample,
                    "father_id": father_value,
                    "mother_id": mother_value,
                    "sex": sex,
                    "phenotype": phenotype,
                }
            )
    if not records:
        raise InputValidationError("pedigree.ped: no records found")
    return records


def validate_pedigree(records: list[dict[str, str | None]]) -> str:
    family_ids = {str(record["family_id"]) for record in records}
    if len(family_ids) != 1:
        raise InputValidationError(
            f"pedigree.ped: expected one family_id, found {sorted(family_ids)}"
        )
    by_id: dict[str, dict[str, str | None]] = {}
    for record in records:
        sample_id = str(record["sample_id"])
        if sample_id in by_id:
            raise InputValidationError(f"pedigree.ped: duplicate sample_id {sample_id!r}")
        by_id[sample_id] = record

    for sample_id, record in by_id.items():
        father = record["father_id"]
        mother = record["mother_id"]
        if father is not None:
            if father not in by_id:
                raise InputValidationError(
                    f"pedigree.ped: father {father!r} of {sample_id!r} is not in the PED"
                )
            if by_id[father]["sex"] == "2":
                raise InputValidationError(
                    f"pedigree.ped: father {father!r} of {sample_id!r} is coded female"
                )
        if mother is not None:
            if mother not in by_id:
                raise InputValidationError(
                    f"pedigree.ped: mother {mother!r} of {sample_id!r} is not in the PED"
                )
            if by_id[mother]["sex"] == "1":
                raise InputValidationError(
                    f"pedigree.ped: mother {mother!r} of {sample_id!r} is coded male"
                )

    visiting: set[str] = set()
    visited: set[str] = set()

    def visit(sample_id: str) -> None:
        if sample_id in visiting:
            raise InputValidationError(
                f"pedigree.ped: ancestry cycle detected at {sample_id!r}"
            )
        if sample_id in visited:
            return
        visiting.add(sample_id)
        record = by_id[sample_id]
        for parent in (record["father_id"], record["mother_id"]):
            if parent is not None:
                visit(parent)
        visiting.remove(sample_id)
        visited.add(sample_id)

    for sample_id in by_id:
        visit(sample_id)
    return next(iter(family_ids))


def select_samples(
    config: dict[str, Any], records: list[dict[str, str | None]]
) -> list[str]:
    eligible = {
        str(record["sample_id"])
        for record in records
        if record["father_id"] is not None and record["mother_id"] is not None
    }
    if not eligible:
        raise InputValidationError(
            "pedigree.ped: no eligible samples have both parents present"
        )
    requested = config.get("samples", {}).get("include")
    if requested is None:
        return sorted(eligible)
    ineligible = sorted(set(requested) - eligible)
    if ineligible:
        raise InputValidationError(
            "samples.include contains ineligible or unknown samples: "
            + ", ".join(ineligible)
        )
    return list(requested)


def parse_fai(path: Path) -> dict[str, int]:
    contigs: dict[str, int] = {}
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                raise InputValidationError(
                    f"reference.fasta_index line {line_number}: expected at least 2 fields"
                )
            contig = fields[0]
            if contig in contigs:
                raise InputValidationError(
                    f"reference.fasta_index: duplicate contig {contig!r}"
                )
            try:
                length = int(fields[1])
            except ValueError as exc:
                raise InputValidationError(
                    f"reference.fasta_index line {line_number}: invalid length {fields[1]!r}"
                ) from exc
            if length <= 0:
                raise InputValidationError(
                    f"reference.fasta_index line {line_number}: length must be positive"
                )
            contigs[contig] = length
    return contigs


def _check_indexed_file(data: Path, index: Path, label: str) -> None:
    _require_file(str(data), label)
    _require_file(str(index), f"{label}.index")
    expected_index = Path(str(data) + ".tbi")
    if index.resolve() != expected_index.resolve():
        raise InputValidationError(
            f"{label}.index must be the conventional tabix path {expected_index}; "
            f"found {index}"
        )
    try:
        with pysam.TabixFile(str(data), index=str(index)):
            pass
    except (OSError, ValueError) as exc:
        raise InputValidationError(
            f"{label}: cannot open data/index pair {data}, {index}: {exc}"
        ) from exc


def inspect_vcf_header(
    vcf_path: Path,
    index_path: Path,
    label: str,
    expected_samples: list[str],
    regions: list[str],
    reference_lengths: dict[str, int],
    *,
    require_depth: bool,
    require_phase_set: bool,
) -> dict[str, Any]:
    """Validate an indexed VCF contract without scanning its records."""
    _check_indexed_file(vcf_path, index_path, label)
    try:
        vcf = pysam.VariantFile(str(vcf_path), index_filename=str(index_path))
    except (OSError, ValueError) as exc:
        raise InputValidationError(f"{label}: cannot open VCF: {exc}") from exc

    with vcf:
        samples = list(vcf.header.samples)
        if set(samples) != set(expected_samples) or len(samples) != len(expected_samples):
            raise InputValidationError(
                f"{label}: VCF sample set {samples} does not equal expected set "
                f"{expected_samples}"
            )
        if "GT" not in vcf.header.formats:
            raise InputValidationError(f"{label}: VCF header is missing FORMAT/GT")
        if require_phase_set and "PS" not in vcf.header.formats:
            raise InputValidationError(f"{label}: VCF header is missing FORMAT/PS")
        depth_fields = sorted(set(vcf.header.formats) & {"DP", "AD", "SD"})
        if require_depth and not depth_fields:
            raise InputValidationError(
                f"{label}: VCF header needs at least one of FORMAT/DP, FORMAT/AD, FORMAT/SD"
            )
        for contig in regions:
            if contig not in vcf.header.contigs:
                raise InputValidationError(f"{label}: VCF header is missing {contig}")
            header_length = vcf.header.contigs[contig].length
            if header_length is None:
                raise InputValidationError(
                    f"{label}: VCF header contig {contig} has no length"
                )
            if header_length != reference_lengths[contig]:
                raise InputValidationError(
                    f"{label}: {contig} length {header_length} differs from reference "
                    f"length {reference_lengths[contig]}"
                )

        return {
            "inspection": "header-and-index",
            "depth_fields": depth_fields,
            "samples": samples,
            "autosomes": regions,
        }


def _open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open(encoding="utf-8")


def inspect_phase_blocks_header(
    path: Path,
    label: str,
) -> dict[str, Any]:
    """Validate the HiPhase table header without scanning genome-wide blocks."""
    _require_file(str(path), label)
    with _open_text(path) as handle:
        header_line = handle.readline().rstrip("\n")
        if header_line.startswith("#"):
            header_line = header_line[1:]
        fieldnames = header_line.split("\t")
        missing = sorted(PHASE_BLOCK_COLUMNS - set(fieldnames))
        if missing:
            raise InputValidationError(
                f"{label}: HiPhase block header is missing columns {missing}"
            )
    return {"inspection": "header", "columns": fieldnames}


def inspect_model_bed_header(
    bed_path: Path,
    index_path: Path,
    label: str,
) -> dict[str, Any]:
    """Validate an indexed model BED contract without scanning CpG records."""
    _check_indexed_file(bed_path, index_path, label)
    metadata: dict[str, str] = {}
    header: list[str] | None = None
    with _open_text(bed_path) as handle:
        for line in handle:
            stripped = line.rstrip("\n")
            if stripped.startswith("##"):
                key, separator, value = stripped[2:].partition("=")
                if separator:
                    metadata[key] = value
                continue
            if stripped.startswith("#"):
                header = stripped[1:].split("\t")
                missing = sorted(MODEL_BED_COLUMNS - set(header))
                if missing:
                    raise InputValidationError(
                        f"{label}: model BED header is missing columns {missing}"
                    )
                break
            if not stripped:
                continue
            raise InputValidationError(f"{label}: data appears before the header")
    if header is None:
        raise InputValidationError(f"{label}: model BED is missing a # header line")
    if metadata.get("pileup-mode") != "model":
        raise InputValidationError(
            f"{label}: expected ##pileup-mode=model, found {metadata.get('pileup-mode')!r}"
        )
    return {
        "inspection": "header-and-index",
        "columns": header,
        "pileup_mode": metadata["pileup-mode"],
        "upstream_min_coverage": metadata.get("min-coverage"),
    }


def _json_fingerprint(*values: Any) -> str:
    payload = json.dumps(values, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _check_publish_collision(outdir: Path, fingerprint: str) -> None:
    if outdir.exists() and not outdir.is_dir():
        raise InputValidationError(f"project.outdir exists but is not a directory: {outdir}")
    if outdir == Path(outdir.anchor):
        raise InputValidationError(f"project.outdir may not be a filesystem root: {outdir}")
    if not outdir.exists() or not any(outdir.iterdir()):
        return
    fingerprint_path = outdir / "pipeline_info" / "config-fingerprint.txt"
    if not fingerprint_path.is_file():
        raise InputValidationError(
            f"project.outdir is nonempty and has no Tapestry fingerprint: {outdir}"
        )
    existing = fingerprint_path.read_text(encoding="utf-8").strip()
    if existing != fingerprint:
        raise InputValidationError(
            f"project.outdir belongs to a different configuration fingerprint: {outdir}"
        )


def _atomic_write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp")
    temporary.write_text(text, encoding="utf-8")
    temporary.replace(path)


def _write_json(path: Path, value: Any) -> None:
    _atomic_write_text(path, json.dumps(value, indent=2, sort_keys=True) + "\n")


def _write_normalized_ped(
    path: Path, records: list[dict[str, str | None]]
) -> None:
    lines = []
    for record in records:
        lines.append(
            "\t".join(
                [
                    str(record["family_id"]),
                    str(record["sample_id"]),
                    str(record["father_id"] or "NA"),
                    str(record["mother_id"] or "NA"),
                    str(record["sex"]),
                    str(record["phenotype"]),
                ]
            )
        )
    _atomic_write_text(path, "\n".join(lines) + "\n")


def _write_selected_samples(
    path: Path,
    selected: list[str],
    records: list[dict[str, str | None]],
    manifest: dict[str, Any],
) -> None:
    ped_by_id = {str(record["sample_id"]): record for record in records}
    manifest_index = {
        sample["id"]: index for index, sample in enumerate(manifest["samples"])
    }
    lines = ["sample_id\tfamily_id\tfather_id\tmother_id\tsex\tmanifest_index"]
    for sample_id in selected:
        record = ped_by_id[sample_id]
        lines.append(
            "\t".join(
                [
                    sample_id,
                    str(record["family_id"]),
                    str(record["father_id"]),
                    str(record["mother_id"]),
                    str(record["sex"]),
                    str(manifest_index[sample_id]),
                ]
            )
        )
    _atomic_write_text(path, "\n".join(lines) + "\n")


def _write_selected_artifacts(
    path: Path,
    selected: list[str],
    manifest: dict[str, Any],
    min_coverage: int,
    mismatch_window_bp: int,
    bigwig: bool,
    regions: list[str],
    reference_fasta: str,
    reference_fai: str,
    reference_name: str,
    config_fingerprint: str,
) -> None:
    manifest_by_id = {sample["id"]: sample for sample in manifest["samples"]}
    samples = []
    for sample_id in selected:
        sample = copy.deepcopy(manifest_by_id[sample_id])
        sample["min_coverage"] = min_coverage
        sample["mismatch_window_bp"] = mismatch_window_bp
        sample["bigwig"] = bigwig
        sample["regions"] = regions
        sample["reference_fasta"] = reference_fasta
        sample["reference_fai"] = reference_fai
        sample["reference_name"] = reference_name
        sample["config_fingerprint"] = config_fingerprint
        samples.append(sample)
    _write_json(path, {"samples": samples})


def _miniwdl_output(outputs: dict[str, Any], name: str) -> Any:
    matches = [
        value for key, value in outputs.items() if key.rsplit(".", 1)[-1] == name
    ]
    if len(matches) != 1:
        raise InputValidationError(
            f"expected exactly one miniwdl output named {name!r}; found {len(matches)}"
        )
    return matches[0]


def _load_miniwdl_outputs(path: Path) -> dict[str, Any]:
    raw = _load_json(path, "miniwdl outputs")
    outputs = raw.get("outputs", raw)
    if not isinstance(outputs, dict):
        raise InputValidationError("miniwdl 'outputs' value must contain an object")
    return outputs


def _detect_wdl_release(outputs: dict[str, Any]) -> str:
    matches = [
        value
        for key, value in outputs.items()
        if key.rsplit(".", 1)[-1] == "workflow_version"
    ]
    if len(matches) > 1:
        raise InputValidationError(
            "found more than one miniwdl workflow_version output"
        )
    detected = matches[0] if matches else None
    if detected is not None and not isinstance(detected, str):
        raise InputValidationError("miniwdl workflow_version must be a string")
    if detected is None:
        raise InputValidationError(
            "miniwdl outputs do not contain workflow_version"
        )
    release = detected
    if release not in SUPPORTED_WDL_RELEASES:
        raise InputValidationError(
            f"miniwdl workflow_version {release!r} is not audited; supported "
            f"releases: {', '.join(sorted(SUPPORTED_WDL_RELEASES))}"
        )
    return release


def _miniwdl_artifact(value: Any, outputs_dir: Path, label: str) -> str | None:
    if value is None:
        return None
    if not isinstance(value, str) or not value:
        raise InputValidationError(
            f"miniwdl output {label!r} must be a path or null"
        )
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = outputs_dir / path
    # Preserve miniwdl's stable out/ symlink rather than recording an
    # engine-private call/work target.
    return str(path.absolute())


def _paired_miniwdl_artifact(
    data_path: str | None,
    index_path: str | None,
    *,
    data_key: str,
    label: str,
) -> dict[str, str] | None:
    if data_path is None and index_path is None:
        return None
    if data_path is None or index_path is None:
        raise InputValidationError(
            f"incomplete {label}: file and index must both be present"
        )
    return {data_key: data_path, "index": index_path}


def _build_miniwdl_manifest(
    outputs: dict[str, Any],
    outputs_path: Path,
    family_id: str,
    ped_ids: list[str],
    release: str,
) -> dict[str, Any]:
    sample_ids = _miniwdl_output(outputs, "sample_ids")
    if not isinstance(sample_ids, list) or not sample_ids or not all(
        isinstance(sample, str) and sample for sample in sample_ids
    ):
        raise InputValidationError(
            "miniwdl sample_ids must be a non-empty string array"
        )
    if len(sample_ids) != len(set(sample_ids)):
        raise InputValidationError("miniwdl sample_ids contains duplicates")
    if set(sample_ids) != set(ped_ids):
        raise InputValidationError(
            f"miniwdl sample set {sorted(sample_ids)} does not equal PED sample set "
            f"{sorted(ped_ids)}"
        )

    arrays = {name: _miniwdl_output(outputs, name) for name in SAMPLE_OUTPUTS}
    for name, values in arrays.items():
        if not isinstance(values, list) or len(values) != len(sample_ids):
            raise InputValidationError(
                f"miniwdl output {name!r} must align with sample_ids "
                f"({len(sample_ids)} entries)"
            )

    samples: list[dict[str, Any]] = []
    for index, sample_id in enumerate(sample_ids):
        artifacts = {
            name: _miniwdl_artifact(values[index], outputs_path.parent, name)
            for name, values in arrays.items()
        }
        cpg_values = [
            artifacts[name]
            for name in (
                "cpg_combined_bed",
                "cpg_combined_bed_index",
                "cpg_hap1_bed",
                "cpg_hap1_bed_index",
                "cpg_hap2_bed",
                "cpg_hap2_bed_index",
            )
        ]
        if any(value is None for value in cpg_values) and not all(
            value is None for value in cpg_values
        ):
            raise InputValidationError(
                f"incomplete CpG model output set for {sample_id!r}"
            )
        cpg_model = None
        if cpg_values[0] is not None:
            cpg_model = {
                "combined": {"bed": cpg_values[0], "index": cpg_values[1]},
                "hap1": {"bed": cpg_values[2], "index": cpg_values[3]},
                "hap2": {"bed": cpg_values[4], "index": cpg_values[5]},
            }
        samples.append(
            {
                "id": sample_id,
                "phased_small_variants": _paired_miniwdl_artifact(
                    artifacts["phased_small_variant_vcf"],
                    artifacts["phased_small_variant_vcf_index"],
                    data_key="vcf",
                    label=f"phased small-variant VCF for {sample_id!r}",
                ),
                "phase_blocks": artifacts["phase_blocks"],
                "cpg_model": cpg_model,
            }
        )

    joint = _paired_miniwdl_artifact(
        _miniwdl_artifact(
            _miniwdl_output(outputs, "joint_small_variants_vcf"),
            outputs_path.parent,
            "joint_small_variants_vcf",
        ),
        _miniwdl_artifact(
            _miniwdl_output(outputs, "joint_small_variants_vcf_index"),
            outputs_path.parent,
            "joint_small_variants_vcf_index",
        ),
        data_key="vcf",
        label="joint small-variant VCF",
    )
    if joint is None:
        raise InputValidationError("miniwdl run has no joint small-variant VCF")
    return {
        "schema_version": 1,
        "provider": "pacbio_hifi_human_wgs_wdl",
        "workflow": {
            "name": "humanwgs_family",
            "release": release,
            "commit": SUPPORTED_WDL_RELEASES[release],
        },
        "family_id": family_id,
        "joint_small_variants": joint,
        "samples": samples,
        "provenance": {
            "engine": "miniwdl",
            "outputs_manifest": str(outputs_path),
        },
    }


def validate_miniwdl_run(
    *,
    outputs_json: Path,
    ped: Path,
    reference_fasta: Path,
    project_outdir: Path,
    output_dir: Path,
    reference_index: Path | None = None,
    reference_gzi: Path | None = None,
    project_id: str | None = None,
    samples: list[str] | None = None,
    regions: str | None = None,
    map_min_qual: float = 20,
    map_min_depth: int = 10,
    min_run_markers: int = 10,
    concordance_min_qual: float = 20,
    concordance_min_depth: int = 5,
    min_coverage: int = 10,
    mismatch_window_bp: int = 50,
    qc_discordance_threshold: float = 0.4,
    qc_min_paired_cpgs: int = 100,
    bigwig: bool = True,
) -> dict[str, Any]:
    """Validate direct miniwdl inputs and publish normalized pipeline records."""
    outputs_path = outputs_json.expanduser().absolute()
    ped_path = ped.expanduser().absolute()
    fasta_path = reference_fasta.expanduser().absolute()
    outdir = project_outdir.expanduser().absolute()
    fai_path = (
        reference_index.expanduser().absolute()
        if reference_index is not None
        else Path(f"{fasta_path}.fai")
    )
    gzi_path = reference_gzi.expanduser().absolute() if reference_gzi else None

    selected_regions = list(AUTOSOMES)
    if regions is not None:
        selected_regions = [region.strip() for region in regions.split(",")]
        if not selected_regions or any(not region for region in selected_regions):
            raise InputValidationError("--regions must be a comma-separated list")
        invalid_regions = sorted(set(selected_regions) - set(AUTOSOMES))
        if invalid_regions:
            raise InputValidationError(
                "--regions contains non-autosomes: " + ", ".join(invalid_regions)
            )
        if len(selected_regions) != len(set(selected_regions)):
            raise InputValidationError("--regions contains duplicates")
    if samples and len(samples) != len(set(samples)):
        raise InputValidationError("--sample selects a sample more than once")
    nonnegative = {
        "--map-min-qual": map_min_qual,
        "--map-min-depth": map_min_depth,
        "--concordance-min-qual": concordance_min_qual,
        "--concordance-min-depth": concordance_min_depth,
        "--mismatch-window-bp": mismatch_window_bp,
    }
    for option, value in nonnegative.items():
        if value < 0:
            raise InputValidationError(f"{option} must be nonnegative")
    if min_run_markers < 1:
        raise InputValidationError("--min-run-markers must be at least 1")
    if min_coverage < 1:
        raise InputValidationError("--min-coverage must be at least 1")
    if not 0.0 <= qc_discordance_threshold <= 1.0:
        raise InputValidationError("--qc-discordance-threshold must be in [0, 1]")
    if qc_min_paired_cpgs < 1:
        raise InputValidationError("--qc-min-paired-cpgs must be at least 1")

    _require_file(str(outputs_path), "outputs-json")
    _require_file(str(ped_path), "pedigree.ped")
    records = parse_ped(ped_path)
    family_id = validate_pedigree(records)
    effective_project_id = project_id or family_id
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", effective_project_id):
        raise InputValidationError(
            "project ID must contain only letters, numbers, '.', '_' or '-'"
        )
    outputs = _load_miniwdl_outputs(outputs_path)
    release = _detect_wdl_release(outputs)
    manifest = _build_miniwdl_manifest(
        outputs,
        outputs_path,
        family_id,
        [str(record["sample_id"]) for record in records],
        release,
    )

    pipeline_info = outdir / "pipeline_info"
    published_run_path = pipeline_info / "resolved-run.json"
    published_manifest_path = pipeline_info / "resolved-manifest.json"
    config: dict[str, Any] = {
        "schema_version": 1,
        "mode": "pedigree",
        "project": {
            "id": effective_project_id,
            "outdir": str(outdir),
        },
        "pedigree": {"ped": str(ped_path)},
        "reference": {
            "name": "GRCh38",
            "fasta": str(fasta_path),
            "fasta_index": str(fai_path),
        },
        "upstream": {
            "manifest": str(published_manifest_path),
        },
        "inheritance": {
            "method": "gtg",
            "map": {
                "min_qual": map_min_qual,
                "min_depth": map_min_depth,
                "min_run_markers": min_run_markers,
            },
            "concordance": {
                "min_qual": concordance_min_qual,
                "min_depth": concordance_min_depth,
            },
        },
        "methylation": {
            "modes": ["model"],
            "min_coverage": min_coverage,
            "mismatch_window_bp": mismatch_window_bp,
            "transmission_qc": {
                "discordance_threshold": qc_discordance_threshold,
                "minimum_paired_cpgs": qc_min_paired_cpgs,
            },
        },
        "regions": {"include": selected_regions},
        "outputs": {"bigwig": bigwig},
    }
    if gzi_path is not None:
        config["reference"]["fasta_gzi"] = str(gzi_path)
    if samples:
        config["samples"] = {"include": samples}
    ped_ids = [str(record["sample_id"]) for record in records]
    selected = select_samples(config, records)
    config.setdefault("samples", {})["include"] = selected

    manifest_by_id = {sample["id"]: sample for sample in manifest["samples"]}
    for sample_id in selected:
        sample = manifest_by_id[sample_id]
        missing = [
            key
            for key in ("phased_small_variants", "phase_blocks", "cpg_model")
            if sample[key] is None
        ]
        if missing:
            raise InputValidationError(
                f"manifest sample {sample_id!r} is selected but has null artifacts: {missing}"
            )

    fasta = _require_file(config["reference"]["fasta"], "reference.fasta")
    fai = _require_file(config["reference"]["fasta_index"], "reference.fasta_index")
    if config["reference"].get("fasta_gzi"):
        _require_file(config["reference"]["fasta_gzi"], "reference.fasta_gzi")
    reference_lengths = parse_fai(fai)
    missing_regions = [
        region for region in selected_regions if region not in reference_lengths
    ]
    if missing_regions:
        raise InputValidationError(
            f"reference.fasta_index is missing configured autosomes: {missing_regions}"
        )

    joint = manifest["joint_small_variants"]
    joint_stats = inspect_vcf_header(
        Path(joint["vcf"]),
        Path(joint["index"]),
        "manifest.joint_small_variants",
        ped_ids,
        selected_regions,
        reference_lengths,
        require_depth=True,
        require_phase_set=False,
    )

    sample_stats: dict[str, Any] = {}
    for sample_id in selected:
        sample = manifest_by_id[sample_id]
        phased = sample["phased_small_variants"]
        model = sample["cpg_model"]
        stats: dict[str, Any] = {
            "phased_small_variants": inspect_vcf_header(
                Path(phased["vcf"]),
                Path(phased["index"]),
                f"manifest.samples[{sample_id}].phased_small_variants",
                [sample_id],
                selected_regions,
                reference_lengths,
                require_depth=False,
                require_phase_set=True,
            ),
            "phase_blocks": inspect_phase_blocks_header(
                Path(sample["phase_blocks"]),
                f"manifest.samples[{sample_id}].phase_blocks",
            ),
            "cpg_model": {},
        }
        for model_label in ("combined", "hap1", "hap2"):
            artifact = model[model_label]
            stats["cpg_model"][model_label] = inspect_model_bed_header(
                Path(artifact["bed"]),
                Path(artifact["index"]),
                f"manifest.samples[{sample_id}].cpg_model.{model_label}",
            )
        sample_stats[sample_id] = stats

    fingerprint = _json_fingerprint(
        {
            "pipeline_version": PIPELINE_VERSION,
            "validator_version": VALIDATOR_VERSION,
            "gtg_commit": GTG_COMMIT,
            "runtime_lock_version": RUNTIME_LOCK_VERSION,
        },
        config,
        manifest,
        {region: reference_lengths[region] for region in selected_regions},
    )
    _check_publish_collision(Path(config["project"]["outdir"]), fingerprint)

    config["_tapestry"] = {
        "config_path": str(published_run_path),
        "manifest_path": str(published_manifest_path),
        "selected_autosomes": selected_regions,
        "config_fingerprint": fingerprint,
        "validator_version": VALIDATOR_VERSION,
        "pipeline_version": PIPELINE_VERSION,
        "gtg_commit": GTG_COMMIT,
        "runtime_lock_version": RUNTIME_LOCK_VERSION,
    }
    manifest["_tapestry"] = {"manifest_path": str(published_manifest_path)}
    report = {
        "status": "valid",
        "validator_version": VALIDATOR_VERSION,
        "config_fingerprint": fingerprint,
        "project_id": config["project"]["id"],
        "output_dir": config["project"]["outdir"],
        "family_id": family_id,
        "selected_samples": selected,
        "selected_autosomes": selected_regions,
        "reference": {
            "name": config["reference"]["name"],
            "fasta": str(fasta),
            "contig_lengths": {
                region: reference_lengths[region] for region in selected_regions
            },
        },
        "settings": {
            "inheritance_map": config["inheritance"]["map"],
            "inheritance_concordance": config["inheritance"]["concordance"],
            "model_min_coverage": config["methylation"]["min_coverage"],
            "mismatch_window_bp": config["methylation"]["mismatch_window_bp"],
            "transmission_qc": config["methylation"]["transmission_qc"],
            "bigwig": config["outputs"]["bigwig"],
        },
        "joint_small_variants": joint_stats,
        "samples": sample_stats,
        "versions": {
            "wdl_release": manifest["workflow"]["release"],
            "wdl_commit": manifest["workflow"]["commit"],
            "pipeline": PIPELINE_VERSION,
            "gtg_commit": GTG_COMMIT,
            "pysam": pysam.__version__,
        },
    }

    output_dir.mkdir(parents=True, exist_ok=True)
    _write_json(output_dir / "resolved-run.json", config)
    _write_json(output_dir / "resolved-manifest.json", manifest)
    _write_normalized_ped(output_dir / "normalized.ped", records)
    _write_selected_samples(
        output_dir / "selected-samples.tsv", selected, records, manifest
    )
    _write_selected_artifacts(
        output_dir / "selected-artifacts.json",
        selected,
        manifest,
        config["methylation"]["min_coverage"],
        config["methylation"]["mismatch_window_bp"],
        config["outputs"]["bigwig"],
        selected_regions,
        config["reference"]["fasta"],
        config["reference"]["fasta_index"],
        config["reference"]["name"],
        fingerprint,
    )
    _write_json(output_dir / "validation-report.json", report)
    _atomic_write_text(output_dir / "config-fingerprint.txt", fingerprint + "\n")
    _atomic_write_text(output_dir / "validation.success", "valid\n")
    return report


def format_validation_summary(report: dict[str, Any]) -> str:
    """Return a concise human-readable summary of a validated run."""
    regions = report["selected_autosomes"]
    region_text = "chr1-chr22" if regions == list(AUTOSOMES) else ", ".join(regions)
    inheritance_map = report["settings"]["inheritance_map"]
    concordance = report["settings"]["inheritance_concordance"]
    lines = [
        "Tapestry validation succeeded",
        f"  Family: {report['family_id']}",
        f"  Samples: {', '.join(report['selected_samples'])}",
        f"  Upstream: PacBio WDL {report['versions']['wdl_release']}",
        f"  Reference: {report['reference']['name']} ({region_text})",
        "  Inheritance map: "
        f"QUAL>={inheritance_map['min_qual']:g}, "
        f"depth>={inheritance_map['min_depth']}, "
        f"run markers>={inheritance_map['min_run_markers']}",
        "  Concordance: "
        f"QUAL>={concordance['min_qual']:g}, "
        f"depth>={concordance['min_depth']}",
        f"  Model minimum coverage: {report['settings']['model_min_coverage']}",
        f"  Mismatch window: {report['settings']['mismatch_window_bp']} bp",
        "  Transmission QC: "
        f"difference>={report['settings']['transmission_qc']['discordance_threshold']:g}, "
        f"paired CpGs>={report['settings']['transmission_qc']['minimum_paired_cpgs']}",
        f"  BigWig: {'enabled' if report['settings']['bigwig'] else 'disabled'}",
        f"  Output: {report['output_dir']}",
    ]
    return "\n".join(lines)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Validate PacBio family-WDL miniwdl outputs and normalize the "
            "Tapestry run contract."
        )
    )
    parser.add_argument("--outputs-json", required=True, type=Path)
    parser.add_argument("--ped", required=True, type=Path)
    parser.add_argument("--reference-fasta", required=True, type=Path)
    parser.add_argument("--reference-index", type=Path)
    parser.add_argument("--reference-gzi", type=Path)
    parser.add_argument("--project-outdir", required=True, type=Path)
    parser.add_argument("--project-id")
    parser.add_argument("--sample", action="append")
    parser.add_argument("--regions")
    parser.add_argument("--map-min-qual", type=float, default=20)
    parser.add_argument("--map-min-depth", type=int, default=10)
    parser.add_argument("--min-run-markers", type=int, default=10)
    parser.add_argument("--concordance-min-qual", type=float, default=20)
    parser.add_argument("--concordance-min-depth", type=int, default=5)
    parser.add_argument("--min-coverage", type=int, default=10)
    parser.add_argument("--mismatch-window-bp", type=int, default=50)
    parser.add_argument("--qc-discordance-threshold", type=float, default=0.4)
    parser.add_argument("--qc-min-paired-cpgs", type=int, default=100)
    parser.add_argument("--no-bigwig", action="store_true")
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        report = validate_miniwdl_run(
            outputs_json=args.outputs_json,
            ped=args.ped,
            reference_fasta=args.reference_fasta,
            reference_index=args.reference_index,
            reference_gzi=args.reference_gzi,
            project_outdir=args.project_outdir,
            output_dir=args.output_dir,
            project_id=args.project_id,
            samples=args.sample,
            regions=args.regions,
            map_min_qual=args.map_min_qual,
            map_min_depth=args.map_min_depth,
            min_run_markers=args.min_run_markers,
            concordance_min_qual=args.concordance_min_qual,
            concordance_min_depth=args.concordance_min_depth,
            min_coverage=args.min_coverage,
            mismatch_window_bp=args.mismatch_window_bp,
            qc_discordance_threshold=args.qc_discordance_threshold,
            qc_min_paired_cpgs=args.qc_min_paired_cpgs,
            bigwig=not args.no_bigwig,
        )
    except InputValidationError as exc:
        args.output_dir.mkdir(parents=True, exist_ok=True)
        _write_json(
            args.output_dir / "validation-report.json",
            {"status": "invalid", "errors": [str(exc)]},
        )
        print(f"tapestry validation failed: {exc}", file=sys.stderr)
        return 2
    print(format_validation_summary(report))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
