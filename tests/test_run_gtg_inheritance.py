"""Tests for deterministic gtg orchestration and publication."""

from __future__ import annotations

import json
import shutil
import sys
import tempfile
import textwrap
import unittest
from pathlib import Path

import pysam

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from run_gtg_inheritance import run_inheritance  # noqa: E402


BCFTOOLS = shutil.which("bcftools")


def _vcf(root: Path, name: str) -> Path:
    plain = root / f"{name}.vcf"
    compressed = root / f"{name}.vcf.gz"
    plain.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr1,length=20>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tFATHER\tMOTHER\tCHILD\n"
        "chr1\t2\t.\tC\tT\t50\tPASS\t.\tGT:DP\t0/0:20\t0/1:20\t0/1:20\n",
        encoding="utf-8",
    )
    pysam.tabix_compress(str(plain), str(compressed), force=True)
    pysam.tabix_index(str(compressed), preset="vcf", force=True)
    plain.unlink()
    return compressed


def _fake_gtg(root: Path) -> tuple[Path, Path]:
    mapper = root / "gtg-ped-map"
    mapper.write_text(
        textwrap.dedent(
            """\
            #!/usr/bin/env python3
            import sys
            from pathlib import Path
            if '--version' in sys.argv:
                print('gtg-ped-map 0.1.0-test')
                raise SystemExit(0)
            prefix = Path(sys.argv[sys.argv.index('--prefix') + 1])
            Path(str(prefix) + '.iht.txt').write_text(
                '#chrom start end FATHER MOTHER CHILD marker_count len markers\\n'
                'chr10 5 8 0/1 2/3 0|2 2 3 b,a\\n'
                'chr2 1 4 0/1 2/3 1|3 2 3 d,c\\n'
            )
            Path(str(prefix) + '.markers.txt').write_text(
                '#chrom pos marker\\nchr10 7 z\\nchr2 3 y\\n'
            )
            Path(str(prefix) + '.recombinants.txt').write_text(
                '#chrom pos sample\\nchr10 8 CHILD\\nchr2 4 CHILD\\n'
            )
            """
        ),
        encoding="utf-8",
    )
    concordance = root / "gtg-concordance"
    concordance.write_text(
        textwrap.dedent(
            """\
            #!/usr/bin/env python3
            import sys
            from pathlib import Path
            import pysam
            if '--version' in sys.argv:
                print('gtg-concordance 0.1.0-test')
                raise SystemExit(0)
            vcf_path = Path(sys.argv[sys.argv.index('--vcf') + 1])
            prefix = Path(sys.argv[sys.argv.index('--prefix') + 1])
            with pysam.VariantFile(str(vcf_path)) as source:
                with pysam.VariantFile(str(prefix) + '.pass.vcf', 'w', header=source.header) as passed:
                    for record in source:
                        passed.write(record)
                with pysam.VariantFile(str(prefix) + '.fail.vcf', 'w', header=source.header):
                    pass
            Path(str(prefix) + '.filtering_stats.txt').write_text(
                '#chrom start end passing failing nocall\\nchr1 1 20 1 0 0\\n'
            )
            """
        ),
        encoding="utf-8",
    )
    mapper.chmod(0o755)
    concordance.chmod(0o755)
    return mapper, concordance


@unittest.skipUnless(BCFTOOLS, "bcftools is required")
class RunGtgInheritanceTests(unittest.TestCase):
    def test_runs_tools_sorts_outputs_and_writes_qc(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            run = root / "resolved-run.json"
            run.write_text(
                json.dumps(
                    {
                        "project": {"id": "FAM"},
                        "inheritance": {
                            "map": {
                                "min_qual": 20,
                                "min_depth": 10,
                                "min_run_markers": 10,
                            },
                            "concordance": {"min_qual": 20, "min_depth": 5},
                        },
                    }
                ),
                encoding="utf-8",
            )
            ped = root / "normalized.ped"
            ped.write_text(
                "FAM FATHER NA NA 1 0\n"
                "FAM MOTHER NA NA 2 0\n"
                "FAM CHILD FATHER MOTHER 0 0\n",
                encoding="utf-8",
            )
            all_sites = _vcf(root, "all-sites")
            map_vcf = _vcf(root, "map")
            normalization = root / "normalization.json"
            normalization.write_text(
                json.dumps(
                    {
                        "all_sites": {
                            "records": 1,
                            "sample_stats": {
                                sample: {
                                    "called_genotypes": 1,
                                    "missing_genotypes": 0,
                                    "called_genotype_depth": {
                                        "observations": 1,
                                        "missing": 0,
                                        "min": 20.0,
                                        "max": 20.0,
                                        "mean": 20.0,
                                        "source_records": {"DP": 1, "AD": 0, "SD": 0},
                                    },
                                }
                                for sample in ("FATHER", "MOTHER", "CHILD")
                            },
                        },
                        "complete_family_map": {"records": 1},
                    }
                ),
                encoding="utf-8",
            )
            mapper, concordance = _fake_gtg(root)
            output = root / "inheritance"

            qc = run_inheritance(
                run,
                ped,
                map_vcf,
                all_sites,
                normalization,
                output,
                str(mapper),
                str(concordance),
                BCFTOOLS or "bcftools",
            )

            self.assertEqual(qc["inheritance"]["blocks"], 2)
            self.assertEqual(qc["inheritance"]["pass_records"], 1)
            self.assertEqual(qc["sample_stats"]["CHILD"]["called_genotypes"], 1)
            self.assertEqual(
                qc["sample_stats"]["CHILD"]["called_genotype_depth"]["mean"],
                20.0,
            )
            iht_lines = (output / "FAM.iht.sorted.txt").read_text(
                encoding="utf-8"
            ).splitlines()
            self.assertTrue(iht_lines[1].startswith("chr2 "))
            self.assertTrue(iht_lines[2].startswith("chr10 "))
            for suffix in ("pass.vcf.gz", "fail.vcf.gz"):
                self.assertTrue((output / f"FAM.{suffix}").is_file())
                self.assertTrue((output / f"FAM.{suffix}.tbi").is_file())
            self.assertTrue((output / "FAM.inheritance-qc.json").is_file())


if __name__ == "__main__":
    unittest.main()
