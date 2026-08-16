#!/usr/bin/env python3
"""Validate and normalize a generic Tapestry pedigree/model-only run."""

from __future__ import annotations

import argparse
import copy
import csv
import gzip
import hashlib
import importlib.metadata
import json
import re
import sys
from pathlib import Path
from typing import Any, Iterable, TextIO

import jsonschema
import pysam
import yaml


VALIDATOR_VERSION = "0.1.0"
PIPELINE_VERSION = "0.1.0-dev"
GTG_COMMIT = "e12aca6b49ee7208952467db4a2a9e2f79b98efb"
RUNTIME_LOCK_VERSION = 1
SUPPORTED_WDL_RELEASES = {
    "v3.3.0": "db06f0af2354d847b971b0548eaade9ff145c912",
    "v3.3.1": "477ef39ad69e86e90897ea7e313b86bfc12a2a96",
}
STABLE_V3_RELEASE = re.compile(r"^v3\.\d+\.\d+$")
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


class _UniqueKeyLoader(yaml.SafeLoader):
    pass


def _construct_unique_mapping(
    loader: _UniqueKeyLoader, node: yaml.nodes.MappingNode, deep: bool = False
) -> dict[Any, Any]:
    mapping: dict[Any, Any] = {}
    for key_node, value_node in node.value:
        key = loader.construct_object(key_node, deep=deep)
        if key in mapping:
            raise InputValidationError(f"duplicate YAML key: {key!r}")
        mapping[key] = loader.construct_object(value_node, deep=deep)
    return mapping


_UniqueKeyLoader.add_constructor(
    yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, _construct_unique_mapping
)


def _unique_json_object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise InputValidationError(f"duplicate JSON key: {key!r}")
        result[key] = value
    return result


def _format_json_path(parts: Iterable[Any]) -> str:
    result = "$"
    for part in parts:
        result += f"[{part}]" if isinstance(part, int) else f".{part}"
    return result


def load_document(path: Path) -> dict[str, Any]:
    """Load strict YAML/JSON without custom tags or duplicate keys."""
    try:
        text = path.read_text(encoding="utf-8")
    except OSError as exc:
        raise InputValidationError(f"cannot read {path}: {exc}") from exc

    try:
        if path.suffix.lower() == ".json":
            value = json.loads(text, object_pairs_hook=_unique_json_object)
        else:
            value = yaml.load(text, Loader=_UniqueKeyLoader)
    except InputValidationError:
        raise
    except (json.JSONDecodeError, yaml.YAMLError) as exc:
        raise InputValidationError(f"cannot parse {path}: {exc}") from exc

    if not isinstance(value, dict):
        raise InputValidationError(f"{path}: document root must be an object")
    return value


def load_schema(path: Path) -> dict[str, Any]:
    try:
        schema = json.loads(path.read_text(encoding="utf-8"))
        jsonschema.Draft202012Validator.check_schema(schema)
    except (OSError, json.JSONDecodeError, jsonschema.SchemaError) as exc:
        raise RuntimeError(f"invalid bundled schema {path}: {exc}") from exc
    return schema


def validate_schema(
    document: dict[str, Any], schema: dict[str, Any], label: str
) -> None:
    validator = jsonschema.Draft202012Validator(schema)
    errors = sorted(validator.iter_errors(document), key=lambda error: list(error.path))
    if not errors:
        return
    details = "; ".join(
        f"{_format_json_path(error.absolute_path)}: {error.message}"
        for error in errors
    )
    raise InputValidationError(f"{label} schema validation failed: {details}")


def _resolve_path(value: str, base_dir: Path) -> str:
    path = Path(value)
    if not path.is_absolute():
        path = base_dir / path
    return str(path.resolve(strict=False))


def resolve_run_config(config: dict[str, Any], config_path: Path) -> dict[str, Any]:
    resolved = copy.deepcopy(config)
    base_dir = config_path.parent
    resolved["project"]["outdir"] = _resolve_path(
        resolved["project"]["outdir"], base_dir
    )
    resolved["pedigree"]["ped"] = _resolve_path(
        resolved["pedigree"]["ped"], base_dir
    )
    for key in ("fasta", "fasta_index", "fasta_gzi"):
        if key in resolved["reference"]:
            resolved["reference"][key] = _resolve_path(
                resolved["reference"][key], base_dir
            )
    resolved["upstream"]["manifest"] = _resolve_path(
        resolved["upstream"]["manifest"], base_dir
    )
    resolved["upstream"].setdefault("allow_unaudited_release", False)
    return resolved


def _resolve_indexed_artifact(
    artifact: dict[str, Any] | None, base_dir: Path, data_key: str
) -> None:
    if artifact is None:
        return
    artifact[data_key] = _resolve_path(artifact[data_key], base_dir)
    artifact["index"] = _resolve_path(artifact["index"], base_dir)


def resolve_manifest(
    manifest: dict[str, Any], manifest_path: Path
) -> dict[str, Any]:
    resolved = copy.deepcopy(manifest)
    base_dir = manifest_path.parent
    _resolve_indexed_artifact(resolved["joint_small_variants"], base_dir, "vcf")
    for sample in resolved["samples"]:
        _resolve_indexed_artifact(sample["phased_small_variants"], base_dir, "vcf")
        if sample["phase_blocks"] is not None:
            sample["phase_blocks"] = _resolve_path(sample["phase_blocks"], base_dir)
        if sample["cpg_model"] is not None:
            for label in ("combined", "hap1", "hap2"):
                _resolve_indexed_artifact(sample["cpg_model"][label], base_dir, "bed")
        for label in ("phase_stats", "phase_haplotags"):
            if sample.get(label) is not None:
                sample[label] = _resolve_path(sample[label], base_dir)
        _resolve_indexed_artifact(sample.get("haplotagged_bam"), base_dir, "bam")
    provenance = resolved.get("provenance", {})
    for label in ("outputs_manifest", "metadata"):
        if provenance.get(label) is not None:
            provenance[label] = _resolve_path(provenance[label], base_dir)
    return resolved


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


def validate_release(config: dict[str, Any], manifest: dict[str, Any]) -> list[str]:
    release = manifest["workflow"]["release"]
    if release in SUPPORTED_WDL_RELEASES:
        commit = manifest["workflow"]["commit"].lower()
        expected_commit = SUPPORTED_WDL_RELEASES[release]
        if commit != expected_commit:
            raise InputValidationError(
                f"upstream release {release!r} must use commit {expected_commit}; "
                f"manifest names {commit}"
            )
        return []
    allow_unaudited = config["upstream"].get("allow_unaudited_release", False)
    if not STABLE_V3_RELEASE.fullmatch(release):
        raise InputValidationError(
            f"upstream release {release!r} is not a stable PacBio WDL v3.x release"
        )
    if not allow_unaudited:
        raise InputValidationError(
            f"upstream release {release!r} is unaudited; supported releases are "
            f"{', '.join(sorted(SUPPORTED_WDL_RELEASES))}. Set "
            "upstream.allow_unaudited_release: true only for a deliberate evaluation run"
        )
    return [
        f"UNAUDITED WDL RELEASE: {release} was accepted by explicit "
        "upstream.allow_unaudited_release opt-in"
    ]


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


def inspect_vcf(
    vcf_path: Path,
    index_path: Path,
    label: str,
    expected_samples: list[str],
    regions: list[str],
    reference_lengths: dict[str, int],
    *,
    require_depth: bool,
    require_phase_set: bool,
    require_numeric_qual: bool,
) -> dict[str, Any]:
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

        region_rank = {contig: rank for rank, contig in enumerate(regions)}
        previous: tuple[int, int] | None = None
        stats: dict[str, Any] = {
            "records": 0,
            "non_pass_records": 0,
            "multiallelic_records": 0,
            "records_missing_qual": 0,
            "depth_fields": depth_fields,
            "samples": samples,
            "sample_called_genotypes": {sample: 0 for sample in samples},
            "sample_missing_genotypes": {sample: 0 for sample in samples},
        }
        try:
            for record in vcf:
                if record.contig not in region_rank:
                    continue
                current = (region_rank[record.contig], record.pos)
                if previous is not None and current < previous:
                    raise InputValidationError(
                        f"{label}: records are not sorted at {record.contig}:{record.pos}"
                    )
                previous = current
                if record.pos > reference_lengths[record.contig]:
                    raise InputValidationError(
                        f"{label}: record {record.contig}:{record.pos} is outside the reference"
                    )
                stats["records"] += 1
                if record.qual is None:
                    stats["records_missing_qual"] += 1
                    if require_numeric_qual:
                        raise InputValidationError(
                            f"{label}: record {record.contig}:{record.pos} has missing QUAL"
                        )
                if record.filter.keys() and "PASS" not in record.filter.keys():
                    stats["non_pass_records"] += 1
                if len(record.alts or ()) > 1:
                    stats["multiallelic_records"] += 1
                for sample in samples:
                    genotype = record.samples[sample].get("GT")
                    if genotype is None or any(allele is None for allele in genotype):
                        stats["sample_missing_genotypes"][sample] += 1
                    else:
                        stats["sample_called_genotypes"][sample] += 1
        except (OSError, ValueError) as exc:
            if isinstance(exc, InputValidationError):
                raise
            raise InputValidationError(f"{label}: failed while reading VCF: {exc}") from exc
    return stats


def _open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open(encoding="utf-8")


def inspect_phase_blocks(
    path: Path,
    label: str,
    expected_sample: str,
    regions: list[str],
    reference_lengths: dict[str, int],
) -> dict[str, Any]:
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
        reader = csv.DictReader(handle, fieldnames=fieldnames, delimiter="\t")
        region_rank = {contig: rank for rank, contig in enumerate(regions)}
        previous: tuple[int, int] | None = None
        records = 0
        for line_number, row in enumerate(reader, 2):
            if row["sample_name"] != expected_sample:
                raise InputValidationError(
                    f"{label} line {line_number}: sample_name {row['sample_name']!r} "
                    f"does not equal {expected_sample!r}"
                )
            chrom = str(row["chrom"])
            if chrom not in region_rank:
                continue
            try:
                start = int(str(row["start"]))
                end = int(str(row["end"]))
                int(str(row["num_variants"]))
            except ValueError as exc:
                raise InputValidationError(
                    f"{label} line {line_number}: start/end/num_variants must be integers"
                ) from exc
            if start < 1 or end < start or end > reference_lengths[chrom]:
                raise InputValidationError(
                    f"{label} line {line_number}: invalid 1-based block {chrom}:{start}-{end}"
                )
            current = (region_rank[chrom], start)
            if previous is not None and current < previous:
                raise InputValidationError(
                    f"{label}: blocks are not sorted at line {line_number}"
                )
            previous = current
            records += 1
    return {"records": records, "sample": expected_sample}


def inspect_model_bed(
    bed_path: Path,
    index_path: Path,
    label: str,
    regions: list[str],
    reference_lengths: dict[str, int],
) -> dict[str, Any]:
    _check_indexed_file(bed_path, index_path, label)
    metadata: dict[str, str] = {}
    header: list[str] | None = None
    records = 0
    region_rank = {contig: rank for rank, contig in enumerate(regions)}
    previous: tuple[int, int, int] | None = None
    with _open_text(bed_path) as handle:
        for line_number, line in enumerate(handle, 1):
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
                continue
            if not stripped:
                continue
            if header is None:
                raise InputValidationError(
                    f"{label} line {line_number}: data appears before the header"
                )
            values = stripped.split("\t")
            if len(values) != len(header):
                raise InputValidationError(
                    f"{label} line {line_number}: expected {len(header)} columns, "
                    f"found {len(values)}"
                )
            row = dict(zip(header, values))
            chrom = row["chrom"]
            if chrom not in region_rank:
                continue
            try:
                start = int(row["begin"])
                end = int(row["end"])
                float(row["mod_score"])
                coverage = int(row["cov"])
            except ValueError as exc:
                raise InputValidationError(
                    f"{label} line {line_number}: invalid coordinate, mod_score, or cov"
                ) from exc
            if start < 0 or end <= start or end > reference_lengths[chrom]:
                raise InputValidationError(
                    f"{label} line {line_number}: invalid BED interval {chrom}:{start}-{end}"
                )
            if coverage < 0:
                raise InputValidationError(
                    f"{label} line {line_number}: cov must be non-negative"
                )
            current = (region_rank[chrom], start, end)
            if previous is not None and current < previous:
                raise InputValidationError(
                    f"{label}: records are not sorted at line {line_number}"
                )
            previous = current
            records += 1
    if header is None:
        raise InputValidationError(f"{label}: model BED is missing a # header line")
    if metadata.get("pileup-mode") != "model":
        raise InputValidationError(
            f"{label}: expected ##pileup-mode=model, found {metadata.get('pileup-mode')!r}"
        )
    return {
        "records": records,
        "pileup_mode": metadata["pileup-mode"],
        "upstream_min_coverage": metadata.get("min-coverage"),
    }


def _iter_manifest_paths(manifest: dict[str, Any]) -> Iterable[tuple[str, str]]:
    provenance = manifest.get("provenance", {})
    for key in ("outputs_manifest", "metadata"):
        if provenance.get(key):
            yield f"manifest.provenance.{key}", provenance[key]
    for sample in manifest["samples"]:
        sample_id = sample["id"]
        for key in ("phase_stats", "phase_haplotags"):
            if sample.get(key):
                yield f"manifest.samples[{sample_id}].{key}", sample[key]
        bam = sample.get("haplotagged_bam")
        if bam:
            yield f"manifest.samples[{sample_id}].haplotagged_bam.bam", bam["bam"]
            yield f"manifest.samples[{sample_id}].haplotagged_bam.index", bam["index"]


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


def validate_run(
    run_config_path: Path,
    output_dir: Path,
    schema_dir: Path | None = None,
) -> dict[str, Any]:
    run_config_path = run_config_path.resolve(strict=False)
    if schema_dir is None:
        schema_dir = Path(__file__).resolve().parents[1] / "schemas"
    run_schema = load_schema(schema_dir / "run.schema.json")
    manifest_schema = load_schema(schema_dir / "upstream-manifest.schema.json")

    raw_config = load_document(run_config_path)
    validate_schema(raw_config, run_schema, "run config")
    config = resolve_run_config(raw_config, run_config_path)

    manifest_path = _require_file(config["upstream"]["manifest"], "upstream.manifest")
    raw_manifest = load_document(manifest_path)
    validate_schema(raw_manifest, manifest_schema, "upstream manifest")
    manifest = resolve_manifest(raw_manifest, manifest_path)

    warnings = validate_release(config, manifest)
    ped_path = _require_file(config["pedigree"]["ped"], "pedigree.ped")
    records = parse_ped(ped_path)
    family_id = validate_pedigree(records)
    if manifest["family_id"] != family_id:
        raise InputValidationError(
            f"manifest.family_id {manifest['family_id']!r} does not equal PED family_id "
            f"{family_id!r}"
        )

    manifest_ids = [sample["id"] for sample in manifest["samples"]]
    if len(manifest_ids) != len(set(manifest_ids)):
        raise InputValidationError("upstream manifest contains duplicate sample IDs")
    ped_ids = [str(record["sample_id"]) for record in records]
    if set(manifest_ids) != set(ped_ids):
        raise InputValidationError(
            f"manifest sample set {sorted(manifest_ids)} does not equal PED sample set "
            f"{sorted(ped_ids)}"
        )
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
    regions = config.get("regions", {}).get("include", list(AUTOSOMES))
    missing_regions = [region for region in regions if region not in reference_lengths]
    if missing_regions:
        raise InputValidationError(
            f"reference.fasta_index is missing configured autosomes: {missing_regions}"
        )

    for label, path in _iter_manifest_paths(manifest):
        _require_file(path, label)

    joint = manifest["joint_small_variants"]
    joint_stats = inspect_vcf(
        Path(joint["vcf"]),
        Path(joint["index"]),
        "manifest.joint_small_variants",
        ped_ids,
        regions,
        reference_lengths,
        require_depth=True,
        require_phase_set=False,
        require_numeric_qual=True,
    )

    sample_stats: dict[str, Any] = {}
    for sample_id in selected:
        sample = manifest_by_id[sample_id]
        phased = sample["phased_small_variants"]
        model = sample["cpg_model"]
        stats: dict[str, Any] = {
            "phased_small_variants": inspect_vcf(
                Path(phased["vcf"]),
                Path(phased["index"]),
                f"manifest.samples[{sample_id}].phased_small_variants",
                [sample_id],
                regions,
                reference_lengths,
                require_depth=False,
                require_phase_set=True,
                require_numeric_qual=False,
            ),
            "phase_blocks": inspect_phase_blocks(
                Path(sample["phase_blocks"]),
                f"manifest.samples[{sample_id}].phase_blocks",
                sample_id,
                regions,
                reference_lengths,
            ),
            "cpg_model": {},
        }
        for model_label in ("combined", "hap1", "hap2"):
            artifact = model[model_label]
            stats["cpg_model"][model_label] = inspect_model_bed(
                Path(artifact["bed"]),
                Path(artifact["index"]),
                f"manifest.samples[{sample_id}].cpg_model.{model_label}",
                regions,
                reference_lengths,
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
        {region: reference_lengths[region] for region in regions},
    )
    _check_publish_collision(Path(config["project"]["outdir"]), fingerprint)

    config["_tapestry"] = {
        "config_path": str(run_config_path),
        "manifest_path": str(manifest_path),
        "selected_autosomes": regions,
        "config_fingerprint": fingerprint,
        "validator_version": VALIDATOR_VERSION,
        "pipeline_version": PIPELINE_VERSION,
        "gtg_commit": GTG_COMMIT,
        "runtime_lock_version": RUNTIME_LOCK_VERSION,
    }
    manifest["_tapestry"] = {"manifest_path": str(manifest_path)}
    report = {
        "status": "valid",
        "warnings": warnings,
        "validator_version": VALIDATOR_VERSION,
        "config_fingerprint": fingerprint,
        "family_id": family_id,
        "selected_samples": selected,
        "selected_autosomes": regions,
        "reference": {"fasta": str(fasta), "contig_lengths": {
            region: reference_lengths[region] for region in regions
        }},
        "joint_small_variants": joint_stats,
        "samples": sample_stats,
        "versions": {
            "wdl_release": manifest["workflow"]["release"],
            "wdl_commit": manifest["workflow"]["commit"],
            "pipeline": PIPELINE_VERSION,
            "gtg_commit": GTG_COMMIT,
            "pysam": pysam.__version__,
            "jsonschema": importlib.metadata.version("jsonschema"),
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
        regions,
        config["reference"]["fasta"],
        config["reference"]["fasta_index"],
        config["reference"]["name"],
        fingerprint,
    )
    _write_json(output_dir / "validation-report.json", report)
    _atomic_write_text(output_dir / "config-fingerprint.txt", fingerprint + "\n")
    _atomic_write_text(output_dir / "validation.success", "valid\n")
    return report


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Validate and normalize a Tapestry schema-v1 run configuration."
    )
    parser.add_argument("--run-config", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--schema-dir", type=Path)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        report = validate_run(args.run_config, args.output_dir, args.schema_dir)
    except (InputValidationError, RuntimeError) as exc:
        args.output_dir.mkdir(parents=True, exist_ok=True)
        _write_json(
            args.output_dir / "validation-report.json",
            {"status": "invalid", "errors": [str(exc)]},
        )
        print(f"tapestry validation failed: {exc}", file=sys.stderr)
        return 2
    for warning in report["warnings"]:
        print(f"WARNING: {warning}", file=sys.stderr)
    print(
        f"validated family {report['family_id']} with samples "
        f"{', '.join(report['selected_samples'])}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
