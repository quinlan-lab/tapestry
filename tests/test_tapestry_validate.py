"""Deterministic tests for the generic run-contract validator."""

from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path
from typing import Any

import pysam
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from tapestry_validate import InputValidationError, validate_run  # noqa: E402


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
    allow_unaudited: bool = False,
    joint_samples: list[str] | None = None,
    selected_artifacts: bool = True,
    region: str = "chr1",
    inheritance_informative: bool = False,
) -> Path:
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
    manifest_path = root / "manifest.json"
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")

    config = {
        "schema_version": 1,
        "mode": "pedigree",
        "project": {"id": "fixture", "outdir": "results/fixture"},
        "pedigree": {"ped": "family.ped"},
        "reference": {
            "name": "GRCh38",
            "fasta": "data/reference.fa",
            "fasta_index": "data/reference.fa.fai",
        },
        "upstream": {
            "manifest": "manifest.json",
            "allow_unaudited_release": allow_unaudited,
        },
        "samples": {"include": ["CHILD"]},
        "inheritance": {
            "method": "gtg",
            "map": {
                "min_qual": 20,
                "min_depth": 10,
                "min_run_markers": 1 if inheritance_informative else 10,
            },
            "concordance": {"min_qual": 20, "min_depth": 5},
        },
        "methylation": {
            "modes": ["model"],
            "min_coverage": 10,
            "mismatch_window_bp": 50,
        },
        "regions": {"include": [region]},
        "outputs": {"bigwig": True},
    }
    config_path = root / "run.yaml"
    config_path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    return config_path


class ValidateRunTests(unittest.TestCase):
    def test_valid_run_emits_normalized_contracts(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config = make_fixture(root)
            output = root / "validation"
            report = validate_run(config, output)

            self.assertEqual(report["status"], "valid")
            self.assertEqual(report["selected_samples"], ["CHILD"])
            self.assertEqual(report["selected_autosomes"], ["chr1"])
            self.assertEqual(report["joint_small_variants"]["records"], 1)
            self.assertTrue((output / "validation.success").is_file())
            self.assertEqual(
                (output / "normalized.ped").read_text(encoding="utf-8").splitlines()[0],
                "FAM\tFATHER\tNA\tNA\t1\t0",
            )
            resolved = json.loads(
                (output / "resolved-run.json").read_text(encoding="utf-8")
            )
            self.assertTrue(Path(resolved["reference"]["fasta"]).is_absolute())

    def test_json_and_yaml_have_the_same_fingerprint(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            yaml_config = make_fixture(root)
            yaml_report = validate_run(yaml_config, root / "yaml-validation")
            json_config = root / "run.json"
            json_config.write_text(
                json.dumps(yaml.safe_load(yaml_config.read_text(encoding="utf-8"))),
                encoding="utf-8",
            )
            json_report = validate_run(json_config, root / "json-validation")
            self.assertEqual(
                yaml_report["config_fingerprint"],
                json_report["config_fingerprint"],
            )

    def test_unknown_and_duplicate_config_keys_are_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config = make_fixture(root)
            parsed = yaml.safe_load(config.read_text(encoding="utf-8"))
            parsed["unexpected"] = True
            config.write_text(yaml.safe_dump(parsed), encoding="utf-8")
            with self.assertRaisesRegex(InputValidationError, "run config schema"):
                validate_run(config, root / "validation")

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config = make_fixture(root)
            config.write_text(
                config.read_text(encoding="utf-8") + "mode: pedigree\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(InputValidationError, "duplicate YAML key"):
                validate_run(config, root / "validation")

    def test_missing_index_and_output_collision_are_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config = make_fixture(root)
            (root / "data" / "family.vcf.gz.tbi").unlink()
            with self.assertRaisesRegex(InputValidationError, "does not exist"):
                validate_run(config, root / "validation")

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config = make_fixture(root)
            fingerprint = root / "results" / "fixture" / "pipeline_info" / "config-fingerprint.txt"
            fingerprint.parent.mkdir(parents=True)
            fingerprint.write_text("different-run\n", encoding="utf-8")
            with self.assertRaisesRegex(InputValidationError, "different configuration fingerprint"):
                validate_run(config, root / "validation")

    def test_pedigree_cycle_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config = make_fixture(root)
            (root / "family.ped").write_text(
                "FAM FATHER CHILD NA 1 0\n"
                "FAM MOTHER NA NA 2 0\n"
                "FAM CHILD FATHER MOTHER 0 0\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(InputValidationError, "cycle"):
                validate_run(config, root / "validation")

    def test_v330_is_supported_without_opt_in(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            report = validate_run(
                make_fixture(root, release="v3.3.0"), root / "validation"
            )
            self.assertEqual(report["warnings"], [])

    def test_supported_release_commit_must_match_tag(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config = make_fixture(root)
            manifest_path = root / "manifest.json"
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            manifest["workflow"]["commit"] = "0" * 40
            manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
            with self.assertRaisesRegex(InputValidationError, "must use commit"):
                validate_run(config, root / "validation")

    def test_other_v3_release_requires_explicit_opt_in(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config = make_fixture(root, release="v3.4.0")
            with self.assertRaisesRegex(InputValidationError, "unaudited"):
                validate_run(config, root / "validation")

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config = make_fixture(
                root, release="v3.4.0", allow_unaudited=True
            )
            report = validate_run(config, root / "validation")
            self.assertEqual(len(report["warnings"]), 1)
            self.assertIn("UNAUDITED", report["warnings"][0])

    def test_non_autosomal_region_is_rejected_by_schema(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config = make_fixture(root, region="chrX")
            with self.assertRaisesRegex(InputValidationError, "run config schema"):
                validate_run(config, root / "validation")

    def test_joint_vcf_sample_set_must_equal_ped(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config = make_fixture(root, joint_samples=["FATHER", "CHILD"])
            with self.assertRaisesRegex(InputValidationError, "VCF sample set"):
                validate_run(config, root / "validation")

    def test_selected_sample_artifacts_may_not_be_null(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config = make_fixture(root, selected_artifacts=False)
            with self.assertRaisesRegex(InputValidationError, "selected but has null"):
                validate_run(config, root / "validation")


if __name__ == "__main__":
    unittest.main()
