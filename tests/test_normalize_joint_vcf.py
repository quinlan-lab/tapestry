"""Regression tests for WDL joint-VCF normalization."""

from __future__ import annotations

import json
import shutil
import sys
import tempfile
import unittest
from pathlib import Path

import pysam

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from normalize_joint_vcf import normalize_joint_vcf  # noqa: E402


BCFTOOLS = shutil.which("bcftools")


def _fixture_vcf(root: Path) -> tuple[Path, Path]:
    plain = root / "joint.vcf"
    compressed = root / "joint.vcf.gz"
    samples = ["MOTHER", "CHILD", "FATHER"]
    header = [
        "##fileformat=VCFv4.2",
        "##contig=<ID=chr1,length=50>",
        '##FILTER=<ID=LowQual,Description="Low quality">',
        '##INFO=<ID=TAG,Number=1,Type=String,Description="Fixture tag">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">',
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">',
        '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">',
        '##FORMAT=<ID=PF,Number=1,Type=String,Description="Phase flag">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
        + "\t".join(samples),
    ]
    records = [
        "chr1\t2\trs1\tC\tT\t50\tPASS\tTAG=keep1\tGT:DP:AD:PS:PF"
        "\t0|1:20:10,10:2:PASS\t0|1:22:11,11:2:PASS\t0|0:18:18,0:2:PASS",
        "chr1\t4\trs2\tT\tA\t21\tLowQual\tTAG=keep2\tGT:DP:AD:PS:PF"
        "\t0|1:20:10,10:4:PASS\t./.:.:.:.:.\t0|0:18:18,0:4:PASS",
        "chr1\t6\trs3\tC\tG,T\t35\tLowQual\tTAG=keep3\tGT:DP:AD:PS:PF"
        "\t0|2:20:8,0,12:6:PASS\t1|2:.:0,11,11:6:PASS\t0|1:18:9,9,0:6:PASS",
    ]
    plain.write_text("\n".join(header + records + [""]), encoding="utf-8")
    pysam.tabix_compress(str(plain), str(compressed), force=True)
    pysam.tabix_index(str(compressed), preset="vcf", force=True)
    plain.unlink()
    return compressed, Path(str(compressed) + ".tbi")


@unittest.skipUnless(BCFTOOLS, "bcftools is required")
class NormalizeJointVcfTests(unittest.TestCase):
    def test_preserves_all_sites_and_builds_complete_family_map(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            input_vcf, input_index = _fixture_vcf(root)
            ped = root / "family.ped"
            ped.write_text(
                "FAM FATHER 0 0 1 0\n"
                "FAM MOTHER 0 0 2 0\n"
                "FAM CHILD FATHER MOTHER 0 0\n",
                encoding="utf-8",
            )
            output = root / "out"
            report = normalize_joint_vcf(
                input_vcf,
                input_index,
                ped,
                ["chr1"],
                output,
                "FAM",
                BCFTOOLS or "bcftools",
            )

            self.assertEqual(report["removed_format_fields"], ["PS", "PF"])
            self.assertEqual(report["all_sites"]["records"], 3)
            self.assertEqual(report["all_sites"]["non_pass_records"], 2)
            self.assertEqual(report["all_sites"]["multiallelic_records"], 1)
            self.assertEqual(report["complete_family_map"]["records"], 2)
            self.assertEqual(report["complete_family_map"]["non_pass_records"], 1)
            self.assertEqual(
                report["all_sites"]["samples"], ["FATHER", "MOTHER", "CHILD"]
            )
            child_stats = report["all_sites"]["sample_stats"]["CHILD"]
            self.assertEqual(child_stats["called_genotypes"], 2)
            self.assertEqual(child_stats["missing_genotypes"], 1)
            self.assertEqual(child_stats["called_genotype_depth"]["mean"], 22.0)
            self.assertEqual(
                child_stats["called_genotype_depth"]["source_records"],
                {"DP": 1, "AD": 1, "SD": 0},
            )

            all_sites = output / "FAM.all-sites.vcf.gz"
            map_vcf = output / "FAM.complete-family-map.vcf.gz"
            for path in (all_sites, map_vcf):
                self.assertTrue(path.is_file())
                self.assertTrue(Path(str(path) + ".tbi").is_file())

            with pysam.VariantFile(str(all_sites)) as vcf:
                self.assertNotIn("PS", vcf.header.formats)
                self.assertNotIn("PF", vcf.header.formats)
                rows = list(vcf)
                self.assertEqual([row.id for row in rows], ["rs1", "rs2", "rs3"])
                self.assertEqual([row.info["TAG"] for row in rows], ["keep1", "keep2", "keep3"])
                self.assertEqual(rows[0].samples["CHILD"]["DP"], 22)

            with pysam.VariantFile(str(map_vcf)) as vcf:
                rows = list(vcf)
                self.assertEqual([row.id for row in rows], ["rs1", "rs3"])

            on_disk_report = json.loads(
                (output / "FAM.vcf-normalization.json").read_text(encoding="utf-8")
            )
            self.assertEqual(on_disk_report["complete_family_map"]["records"], 2)


if __name__ == "__main__":
    unittest.main()
