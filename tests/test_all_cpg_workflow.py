"""Tests for reusable reference CpGs, model-only expansion, and publication."""

from __future__ import annotations

import gzip
import json
import sys
import tempfile
import unittest
from pathlib import Path

import polars as pl
import pysam

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from expand_model_to_all_cpgs import (  # noqa: E402
    expand_model_to_all_cpgs,
    write_output,
)
from generate_reference_cpgs import generate_reference_cpgs  # noqa: E402
from write_results_manifest import (  # noqa: E402
    build_results_manifest,
    format_completion_summary,
)


def _bgzip(path: Path, text: str, preset: str) -> Path:
    plain = path.with_suffix("")
    plain.write_text(text, encoding="utf-8")
    pysam.tabix_compress(str(plain), str(path), force=True)
    pysam.tabix_index(str(path), preset=preset, force=True)
    plain.unlink()
    return path


def _sample_qc(root: Path, status: str = "complete") -> Path:
    path = root / "CHILD.all-cpgs-qc.json"
    path.write_text(
        json.dumps(
            {
                "sample_id": "CHILD",
                "status": "complete",
                "founder_phasing_status": status,
            }
        ),
        encoding="utf-8",
    )
    return path


class AllCpgWorkflowTests(unittest.TestCase):
    def test_reference_and_sample_created_cpgs_are_retained(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fasta = root / "reference.fa"
            fasta.write_text(">chr1\nACGTAACGTT\n", encoding="utf-8")
            pysam.faidx(str(fasta))
            reference_cpgs = root / "reference-cpgs.bed.gz"
            report = generate_reference_cpgs(
                fasta, Path(str(fasta) + ".fai"), ["chr1"], reference_cpgs
            )
            self.assertEqual(report["records"], 2)

            model = _bgzip(
                root / "combined.bed.gz",
                "##pileup-mode=model\n"
                "#chrom\tbegin\tend\tmod_score\ttype\tcov\n"
                "chr1\t1\t2\t75\tCG\t12\n"
                "chr1\t8\t9\t25\tCG\t11\n",
                "bed",
            )
            founder = root / "founder.bed"
            founder_header = root / "founder.bed.header"
            founder_columns = [
                "chrom",
                "start",
                "end",
                "methylation_level_pat_count",
                "methylation_level_mat_count",
                "methylation_level_pat_model",
                "methylation_level_mat_model",
            ]
            founder_header.write_text("\n".join(founder_columns) + "\n", encoding="utf-8")
            founder.write_text("chr1\t1\t2\t.\t.\t0.75\t0.25\n", encoding="utf-8")
            mismatch = root / "mismatch.bed"
            mismatch_header = root / "mismatch.bed.header"
            mismatch_header.write_text("chrom\nstart\nend\nREF\nALT\n", encoding="utf-8")
            mismatch.write_text("chr1\t1\t2\tC\tT\n", encoding="utf-8")
            joint_vcf = _bgzip(
                root / "joint.vcf.gz",
                "##fileformat=VCFv4.2\n"
                "##contig=<ID=chr1,length=10>\n"
                '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tCHILD\n"
                "chr1\t2\t.\tC\tT\t50\tPASS\t.\tGT\t0|1\n",
                "vcf",
            )

            frame = expand_model_to_all_cpgs(
                reference_cpgs,
                model,
                founder,
                founder_header,
                mismatch,
                mismatch_header,
                joint_vcf,
                "CHILD",
                ["chr1"],
                50,
            )
            self.assertEqual(frame["start_cpg"].to_list(), [1, 6, 8])
            self.assertEqual(frame["methylation_level_count"].dtype, pl.Float64)
            sample_created = frame.filter(pl.col("start_cpg") == 8)
            self.assertEqual(sample_created["methylation_level_model"].item(), 0.25)
            self.assertTrue(
                frame.filter(pl.col("start_cpg") == 1)[
                    "cpg_is_allele_specific"
                ].item()
            )

            output = root / "all-cpgs.bed.gz"
            write_output(frame, output, {"enabled-modes": ["model"]})
            self.assertTrue(Path(str(output) + ".tbi").is_file())
            with gzip.open(output, "rt", encoding="utf-8") as handle:
                text = handle.read()
            self.assertIn('##enabled-modes=["model"]', text)
            self.assertIn("#chrom\tstart_cpg\tend_cpg", text)

    def test_results_manifest_contains_no_work_paths(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            resolved = root / "resolved-run.json"
            resolved.write_text(
                json.dumps(
                    {
                        "project": {"id": "FAM"},
                        "reference": {"name": "TESTREF"},
                        "outputs": {"bigwig": True},
                        "_tapestry": {
                            "config_fingerprint": "abc123",
                            "selected_autosomes": ["chr1"],
                        },
                    }
                ),
                encoding="utf-8",
            )
            selected = root / "selected-samples.tsv"
            selected.write_text(
                "sample_id\tfamily_id\nCHILD\tFAM\n", encoding="utf-8"
            )
            output = root / "results-manifest.json"
            manifest = build_results_manifest(
                resolved, selected, [_sample_qc(root)], output
            )
            self.assertEqual(manifest["selected_samples"], ["CHILD"])
            self.assertEqual(
                manifest["reference"]["cpgs"]["path"],
                "reference/TESTREF.autosomes.cpgs.bed.gz",
            )
            self.assertEqual(
                manifest["samples"]["CHILD"]["founder_phasing"]["paternal_bigwig"]["path"],
                "samples/CHILD/CHILD.dna-methylation.pat.model.TESTREF.bw",
            )
            self.assertEqual(
                manifest["visualizations"]["haplotype_ancestry"]["path"],
                "visualizations/haplotype-ancestry/index.html",
            )
            self.assertEqual(
                manifest["visualizations"]["haplotype_ancestry"]["bundle"],
                "visualizations/haplotype-ancestry",
            )
            self.assertNotIn(str(root), output.read_text(encoding="utf-8"))
            summary = format_completion_summary(manifest, root / "published")
            self.assertIn("Tapestry completed", summary)
            self.assertIn(
                f"Results manifest: {root / 'published' / 'results-manifest.json'}",
                summary,
            )
            self.assertIn("Samples: 1 complete, 0 without inheritance phase", summary)

    def test_results_manifest_marks_bigwig_disabled(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            resolved = root / "resolved-run.json"
            resolved.write_text(
                json.dumps(
                    {
                        "project": {"id": "FAM"},
                        "reference": {"name": "GRCh38"},
                        "outputs": {"bigwig": False},
                        "_tapestry": {
                            "config_fingerprint": "abc123",
                            "selected_autosomes": ["chr1"],
                        },
                    }
                ),
                encoding="utf-8",
            )
            selected = root / "selected-samples.tsv"
            selected.write_text("sample_id\nCHILD\n", encoding="utf-8")
            manifest = build_results_manifest(
                resolved,
                selected,
                [_sample_qc(root, "no_inheritance_phase")],
                root / "results-manifest.json",
            )
            self.assertEqual(
                manifest["samples"]["CHILD"]["status"],
                "no_inheritance_phase",
            )
            founder = manifest["samples"]["CHILD"]["founder_phasing"]
            self.assertEqual(founder["bigwig"]["status"], "disabled")
            self.assertNotIn("paternal_bigwig", founder)


if __name__ == "__main__":
    unittest.main()
