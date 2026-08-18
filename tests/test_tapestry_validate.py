"""Deterministic tests for the generic run-contract validator."""

from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path
from typing import Any

import pysam

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from tapestry_validate import (  # noqa: E402
    InputValidationError,
    format_validation_summary,
    validate_miniwdl_run,
)


COMMIT_331 = "477ef39ad69e86e90897ea7e313b86bfc12a2a96"
COMMIT_330 = "db06f0af2354d847b971b0548eaade9ff145c912"


def _write_bgzip(path: Path, text: str, preset: str) -> tuple[Path, Path]:
    plain = path.with_suffix("")
    plain.write_text(text, encoding="utf-8")
    pysam.tabix_compress(str(plain), str(path), force=True)
    pysam.tabix_index(str(path), preset=preset, force=True)
    plain.unlink()
    return path, Path(str(path) + ".tbi")


def _write_vcf(path: Path, samples: list[str], phased: bool) -> tuple[Path, Path]:
    formats = [
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">',
    ]
    if phased:
        formats.append(
            '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">'
        )
    fmt = "GT:DP:PS" if phased else "GT:DP"
    calls = ["0|1:20:2" if phased else "0/1:20" for _ in samples]
    text = "\n".join(
        [
            "##fileformat=VCFv4.2",
            "##contig=<ID=chr1,length=20>",
            *formats,
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
            + "\t".join(samples),
            "chr1\t2\t.\tC\tT\t50\tPASS\t.\t" + fmt + "\t" + "\t".join(calls),
            "",
        ]
    )
    return _write_bgzip(path, text, "vcf")


def _write_e2e_vcf(
    path: Path,
    samples: list[str],
    calls: list[str] | list[list[str]],
    phased: bool,
) -> tuple[Path, Path]:
    formats = [
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">',
    ]
    if phased:
        formats.append('##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">')
    fmt = "GT:DP:PS" if phased else "GT:DP"
    rows = []
    for index, position in enumerate((2, 6, 10, 14, 18)):
        selected_calls = calls[index] if calls and isinstance(calls[0], list) else calls
        assert isinstance(selected_calls, list)
        assert not selected_calls or isinstance(selected_calls[0], str)
        rows.append(
            f"chr1\t{position}\t.\tC\tT\t50\tPASS\t.\t{fmt}\t"
            + "\t".join(selected_calls)
        )
    text = "\n".join(
        [
            "##fileformat=VCFv4.2",
            "##contig=<ID=chr1,length=20>",
            *formats,
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
            + "\t".join(samples),
            *rows,
            "",
        ]
    )
    return _write_bgzip(path, text, "vcf")


def _write_model_bed(path: Path) -> tuple[Path, Path]:
    text = "\n".join(
        [
            "##pileup-mode=model",
            "##min-coverage=4",
            "#chrom\tbegin\tend\tmod_score\ttype\tcov",
            "chr1\t1\t2\t75\tCG\t12",
            "",
        ]
    )
    return _write_bgzip(path, text, "bed")


def make_fixture(
    root: Path,
    *,
    release: str = "v3.3.1",
    joint_samples: list[str] | None = None,
    selected_artifacts: bool = True,
    inheritance_informative: bool = False,
) -> None:
    data = root / "data"
    data.mkdir()
    fasta = data / "reference.fa"
    fasta.write_text(">chr1\nACGTCGACGTCGACGTCGAC\n", encoding="utf-8")
    pysam.faidx(str(fasta))

    ped = root / "family.ped"
    ped.write_text(
        "FAM FATHER NA NA 1 0\n"
        "FAM MOTHER . . 2 0\n"
        "FAM CHILD FATHER MOTHER 0 0\n",
        encoding="utf-8",
    )

    all_samples = ["FATHER", "MOTHER", "CHILD"]
    if inheritance_informative and joint_samples is None:
        joint_vcf, joint_index = _write_e2e_vcf(
            data / "family.vcf.gz",
            all_samples,
            [
                ["0/1:20", "0/0:20", "1/0:20"],
                ["0/0:20", "0/1:20", "0/1:20"],
                ["0/1:20", "0/0:20", "1/0:20"],
                ["0/0:20", "0/1:20", "0/1:20"],
                ["0/1:20", "0/0:20", "1/0:20"],
            ],
            phased=False,
        )
        child_vcf, child_index = _write_e2e_vcf(
            data / "CHILD.vcf.gz",
            ["CHILD"],
            [
                ["1|0:20:2"],
                ["0|1:20:2"],
                ["1|0:20:2"],
                ["0|1:20:2"],
                ["1|0:20:2"],
            ],
            phased=True,
        )
    else:
        joint_vcf, joint_index = _write_vcf(
            data / "family.vcf.gz", joint_samples or all_samples, phased=False
        )
        child_vcf, child_index = _write_vcf(
            data / "CHILD.vcf.gz", ["CHILD"], phased=True
        )
    blocks = data / "CHILD.blocks.tsv"
    blocks.write_text(
        "source_block_index\tsample_name\tphase_block_id\tchrom\tstart\tend\tnum_variants\n"
        + (
            "0\tCHILD\t2\tchr1\t1\t20\t5\n"
            if inheritance_informative
            else "0\tCHILD\t2\tchr1\t2\t2\t1\n"
        ),
        encoding="utf-8",
    )
    beds: dict[str, tuple[Path, Path]] = {}
    for label in ("combined", "hap1", "hap2"):
        beds[label] = _write_model_bed(data / f"CHILD.{label}.bed.gz")

    if selected_artifacts:
        child: dict[str, Any] = {
            "id": "CHILD",
            "phased_small_variants": {
                "vcf": str(child_vcf.relative_to(root)),
                "index": str(child_index.relative_to(root)),
            },
            "phase_blocks": str(blocks.relative_to(root)),
            "cpg_model": {
                label: {
                    "bed": str(paths[0].relative_to(root)),
                    "index": str(paths[1].relative_to(root)),
                }
                for label, paths in beds.items()
            },
        }
    else:
        child = {
            "id": "CHILD",
            "phased_small_variants": None,
            "phase_blocks": None,
            "cpg_model": None,
        }

    commit = COMMIT_330 if release == "v3.3.0" else COMMIT_331
    manifest = {
        "schema_version": 1,
        "provider": "pacbio_hifi_human_wgs_wdl",
        "workflow": {
            "name": "humanwgs_family",
            "release": release,
            "commit": commit,
        },
        "family_id": "FAM",
        "joint_small_variants": {
            "vcf": str(joint_vcf.relative_to(root)),
            "index": str(joint_index.relative_to(root)),
        },
        "samples": [
            {
                "id": "FATHER",
                "phased_small_variants": None,
                "phase_blocks": None,
                "cpg_model": None,
            },
            {
                "id": "MOTHER",
                "phased_small_variants": None,
                "phase_blocks": None,
                "cpg_model": None,
            },
            child,
        ],
    }
    miniwdl_outputs = {
        "humanwgs_family.sample_ids": [
            sample["id"] for sample in manifest["samples"]
        ],
        "humanwgs_family.workflow_version": release,
        "humanwgs_family.joint_small_variants_vcf": manifest[
            "joint_small_variants"
        ]["vcf"],
        "humanwgs_family.joint_small_variants_vcf_index": manifest[
            "joint_small_variants"
        ]["index"],
        "humanwgs_family.phased_small_variant_vcf": [
            sample["phased_small_variants"]["vcf"]
            if sample["phased_small_variants"]
            else None
            for sample in manifest["samples"]
        ],
        "humanwgs_family.phased_small_variant_vcf_index": [
            sample["phased_small_variants"]["index"]
            if sample["phased_small_variants"]
            else None
            for sample in manifest["samples"]
        ],
        "humanwgs_family.phase_blocks": [
            sample["phase_blocks"] for sample in manifest["samples"]
        ],
    }
    for label in ("combined", "hap1", "hap2"):
        miniwdl_outputs[f"humanwgs_family.cpg_{label}_bed"] = [
            sample["cpg_model"][label]["bed"] if sample["cpg_model"] else None
            for sample in manifest["samples"]
        ]
        miniwdl_outputs[f"humanwgs_family.cpg_{label}_bed_index"] = [
            sample["cpg_model"][label]["index"] if sample["cpg_model"] else None
            for sample in manifest["samples"]
        ]
    (root / "miniwdl-outputs.json").write_text(
        json.dumps(miniwdl_outputs), encoding="utf-8"
    )


def validate_fixture(root: Path, **overrides: Any) -> dict[str, Any]:
    arguments: dict[str, Any] = {
        "outputs_json": root / "miniwdl-outputs.json",
        "ped": root / "family.ped",
        "reference_fasta": root / "data/reference.fa",
        "project_outdir": root / "results/fixture",
        "output_dir": root / "validation",
        "samples": ["CHILD"],
        "regions": "chr1",
    }
    arguments.update(overrides)
    return validate_miniwdl_run(**arguments)


class ValidateRunTests(unittest.TestCase):
    def test_direct_miniwdl_inputs_are_validated_and_normalized(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            make_fixture(root)
            output = root / "validation"
            report = validate_fixture(root, output_dir=output)

            self.assertEqual(report["status"], "valid")
            self.assertEqual(report["versions"]["wdl_release"], "v3.3.1")
            self.assertEqual(report["selected_samples"], ["CHILD"])
            self.assertEqual(report["selected_autosomes"], ["chr1"])
            self.assertEqual(
                report["joint_small_variants"]["inspection"], "header-and-index"
            )
            self.assertNotIn("records", report["joint_small_variants"])
            self.assertEqual(
                report["samples"]["CHILD"]["phase_blocks"]["inspection"],
                "header",
            )
            self.assertEqual(
                report["samples"]["CHILD"]["cpg_model"]["combined"][
                    "inspection"
                ],
                "header-and-index",
            )
            self.assertTrue((output / "validation.success").is_file())
            self.assertEqual(
                (output / "normalized.ped").read_text(encoding="utf-8").splitlines()[0],
                "FAM\tFATHER\tNA\tNA\t1\t0",
            )
            resolved = json.loads(
                (output / "resolved-run.json").read_text(encoding="utf-8")
            )
            self.assertEqual(resolved["mode"], "pedigree")
            self.assertEqual(resolved["reference"]["name"], "GRCh38")
            self.assertEqual(resolved["inheritance"]["method"], "gtg")
            self.assertEqual(resolved["methylation"]["modes"], ["model"])
            self.assertEqual(resolved["regions"]["include"], ["chr1"])
            self.assertEqual(
                resolved["upstream"]["manifest"],
                str(root / "results/fixture/pipeline_info/resolved-manifest.json"),
            )
            summary = format_validation_summary(report)
            self.assertIn("Tapestry validation succeeded", summary)
            self.assertIn("Reference: GRCh38 (chr1)", summary)
            self.assertIn("Inheritance map: QUAL>=20", summary)
            self.assertIn("Mismatch window: 50 bp", summary)

    def test_direct_miniwdl_inputs_reject_partial_cpg_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            make_fixture(root)
            outputs_path = root / "miniwdl-outputs.json"
            outputs = json.loads(outputs_path.read_text(encoding="utf-8"))
            outputs["humanwgs_family.cpg_hap2_bed"][-1] = None
            outputs_path.write_text(json.dumps(outputs), encoding="utf-8")

            with self.assertRaisesRegex(
                InputValidationError, "incomplete CpG model output set"
            ):
                validate_fixture(root)

    def test_miniwdl_outputs_require_workflow_version(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            make_fixture(root)
            outputs_path = root / "miniwdl-outputs.json"
            outputs = json.loads(outputs_path.read_text(encoding="utf-8"))
            del outputs["humanwgs_family.workflow_version"]
            outputs_path.write_text(json.dumps(outputs), encoding="utf-8")
            with self.assertRaisesRegex(InputValidationError, "workflow_version"):
                validate_fixture(root)

    def test_missing_index_and_output_collision_are_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            make_fixture(root)
            (root / "data" / "family.vcf.gz.tbi").unlink()
            with self.assertRaisesRegex(InputValidationError, "does not exist"):
                validate_fixture(root)

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            make_fixture(root)
            fingerprint = root / "results" / "fixture" / "pipeline_info" / "config-fingerprint.txt"
            fingerprint.parent.mkdir(parents=True)
            fingerprint.write_text("different-run\n", encoding="utf-8")
            with self.assertRaisesRegex(InputValidationError, "different configuration fingerprint"):
                validate_fixture(root)

    def test_pedigree_cycle_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            make_fixture(root)
            (root / "family.ped").write_text(
                "FAM FATHER CHILD NA 1 0\n"
                "FAM MOTHER NA NA 2 0\n"
                "FAM CHILD FATHER MOTHER 0 0\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(InputValidationError, "cycle"):
                validate_fixture(root)

    def test_v330_is_supported_without_opt_in(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            make_fixture(root, release="v3.3.0")
            report = validate_fixture(root)
            self.assertEqual(report["versions"]["wdl_release"], "v3.3.0")

    def test_other_wdl_release_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            make_fixture(root, release="v3.4.0")
            with self.assertRaisesRegex(InputValidationError, "not audited"):
                validate_fixture(root)

    def test_non_autosomal_and_duplicate_regions_are_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            make_fixture(root)
            with self.assertRaisesRegex(InputValidationError, "non-autosomes"):
                validate_fixture(root, regions="chrX")
            with self.assertRaisesRegex(InputValidationError, "duplicates"):
                validate_fixture(root, regions="chr1,chr1")

    def test_joint_vcf_sample_set_must_equal_ped(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            make_fixture(root, joint_samples=["FATHER", "CHILD"])
            with self.assertRaisesRegex(InputValidationError, "VCF sample set"):
                validate_fixture(root)

    def test_non_vcf_artifact_headers_are_checked(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            make_fixture(root)
            (root / "data" / "CHILD.blocks.tsv").write_text(
                "sample_name\tchrom\tstart\tend\n", encoding="utf-8"
            )
            with self.assertRaisesRegex(InputValidationError, "missing columns"):
                validate_fixture(root)

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            make_fixture(root)
            _write_bgzip(
                root / "data" / "CHILD.combined.bed.gz",
                "##pileup-mode=count\n"
                "#chrom\tbegin\tend\tmod_score\ttype\tcov\n"
                "chr1\t1\t2\t75\tCG\t12\n",
                "bed",
            )
            with self.assertRaisesRegex(InputValidationError, "pileup-mode=model"):
                validate_fixture(root)

    def test_selected_sample_artifacts_may_not_be_null(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            make_fixture(root, selected_artifacts=False)
            with self.assertRaisesRegex(InputValidationError, "selected but has null"):
                validate_fixture(root)

    def test_invalid_direct_settings_are_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            make_fixture(root)
            with self.assertRaisesRegex(InputValidationError, "project ID"):
                validate_fixture(root, project_id="bad id")
            with self.assertRaisesRegex(InputValidationError, "at least 1"):
                validate_fixture(root, min_coverage=-1)
            with self.assertRaisesRegex(InputValidationError, "more than once"):
                validate_fixture(root, samples=["CHILD", "CHILD"])


if __name__ == "__main__":
    unittest.main()
