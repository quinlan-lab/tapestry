"""End-to-end test for the model-only founder-phasing command."""

from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import pyBigWig
import pysam


def _bgzip(path: Path, text: str, preset: str) -> Path:
    plain = path.with_suffix("")
    plain.write_text(text, encoding="utf-8")
    pysam.tabix_compress(str(plain), str(path), force=True)
    pysam.tabix_index(str(path), preset=preset, force=True)
    plain.unlink()
    return path


class FounderPhasingCliTests(unittest.TestCase):
    def test_model_only_cli_publishes_founder_and_mismatch_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            read_vcf = _bgzip(
                root / "child.vcf.gz",
                "##fileformat=VCFv4.2\n"
                "##contig=<ID=chr1,length=20>\n"
                '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
                '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">\n'
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tCHILD\n"
                "chr1\t2\t.\tC\tT\t50\tPASS\t.\tGT:PS\t0|1:2\n",
                "vcf",
            )
            inheritance_vcf = _bgzip(
                root / "inheritance.vcf.gz",
                "##fileformat=VCFv4.2\n"
                "##contig=<ID=chr1,length=20>\n"
                '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tFATHER\tMOTHER\tCHILD\n"
                "chr1\t2\t.\tC\tT\t50\tPASS\t.\tGT\t0/0\t1/1\t0|1\n",
                "vcf",
            )
            blocks = root / "blocks.tsv"
            blocks.write_text(
                "source_block_index\tsample_name\tphase_block_id\tchrom\tstart\tend\tnum_variants\n"
                "0\tCHILD\t2\tchr1\t1\t20\t1\n",
                encoding="utf-8",
            )
            iht = root / "family.iht.txt"
            iht.write_text(
                "#chrom start end FATHER MOTHER CHILD marker_count len markers\n"
                "chr1 0 20 0/1 2/3 0|2 1 20 marker\n",
                encoding="utf-8",
            )
            beds = []
            for haplotype, score in (("hap1", 80), ("hap2", 20)):
                beds.append(
                    _bgzip(
                        root / f"{haplotype}.bed.gz",
                        "##pileup-mode=model\n"
                        "##min-coverage=10\n"
                        "#chrom\tbegin\tend\tmod_score\ttype\tcov\n"
                        f"chr1\t1\t2\t{score}\tCG\t12\n",
                        "bed",
                    )
                )
            fai = root / "reference.fa.fai"
            fai.write_text("chr1\t20\t6\t20\t21\n", encoding="utf-8")
            output = root / "output"
            command = [
                sys.executable,
                str(Path(__file__).resolve().parents[1] / "src" / "phase_meth_to_founder_haps.py"),
                "--uid",
                "CHILD",
                "--vcf_read_phased",
                str(read_vcf),
                "--tsv_read_phase_blocks",
                str(blocks),
                "--vcf_iht_phased",
                str(inheritance_vcf),
                "--txt_iht_blocks",
                str(iht),
                "--bed_meth_model_hap1",
                str(beds[0]),
                "--bed_meth_model_hap2",
                str(beds[1]),
                "--reference_fai",
                str(fai),
                "--reference_name",
                "GRCh38",
                "--regions",
                "chr1",
                "--output_dir",
                str(output),
            ]
            environment = os.environ.copy()
            source = str(Path(__file__).resolve().parents[1] / "src")
            environment["PYTHONPATH"] = f"{source}:{source}/util"
            completed = subprocess.run(
                command, check=False, env=environment, capture_output=True, text=True
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)

            qc = json.loads((output / "CHILD.phasing-qc.json").read_text(encoding="utf-8"))
            self.assertEqual(qc["status"], "complete")
            self.assertEqual(qc["enabled_modes"], ["model"])
            self.assertTrue((output / "CHILD.bit-vector-sites-mismatches.vcf.gz.tbi").is_file())
            self.assertTrue((output / "CHILD.dna-methylation.bed").is_file())
            with pyBigWig.open(
                str(output / "CHILD.dna-methylation.pat.model.GRCh38.bw")
            ) as bigwig:
                self.assertAlmostEqual(bigwig.values("chr1", 1, 2)[0], 0.8)

            no_bigwig_output = root / "output-no-bigwig"
            no_bigwig_command = command[:-1] + [str(no_bigwig_output), "--no-bigwig"]
            completed = subprocess.run(
                no_bigwig_command,
                check=False,
                env=environment,
                capture_output=True,
                text=True,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            self.assertFalse(list(no_bigwig_output.glob("*.bw")))
            self.assertTrue((no_bigwig_output / "CHILD.phasing-qc.json").is_file())


if __name__ == "__main__":
    unittest.main()
