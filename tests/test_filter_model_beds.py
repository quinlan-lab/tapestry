"""Coverage and metadata tests for model BED filtering."""

from __future__ import annotations

import gzip
import sys
import tempfile
import unittest
from pathlib import Path

import pysam

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from filter_model_beds import filter_model_bed  # noqa: E402


class FilterModelBedsTests(unittest.TestCase):
    def test_filters_coverage_and_regions_while_preserving_metadata(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            plain = root / "input.bed"
            compressed = root / "input.bed.gz"
            plain.write_text(
                "##pileup-mode=model\n"
                "##min-coverage=4\n"
                "#chrom\tbegin\tend\tmod_score\ttype\tcov\n"
                "chr1\t1\t2\t70\tCG\t9\n"
                "chr1\t3\t4\t80\tCG\t10\n"
                "chrX\t5\t6\t90\tCG\t20\n",
                encoding="utf-8",
            )
            pysam.tabix_compress(str(plain), str(compressed), force=True)
            pysam.tabix_index(str(compressed), preset="bed", force=True)
            plain.unlink()
            output = root / "filtered.bed.gz"

            report = filter_model_bed(
                compressed,
                Path(str(compressed) + ".tbi"),
                output,
                10,
                ["chr1"],
            )

            self.assertEqual(report["input_records"], 3)
            self.assertEqual(report["retained_records"], 1)
            self.assertEqual(report["low_coverage_records"], 1)
            self.assertEqual(report["excluded_region_records"], 1)
            self.assertEqual(report["upstream_min_coverage"], "4")
            self.assertTrue(Path(str(output) + ".tbi").is_file())
            with gzip.open(output, "rt", encoding="utf-8") as handle:
                text = handle.read()
            self.assertIn("##pileup-mode=model", text)
            self.assertIn("##tapestry-min-coverage=10", text)
            self.assertIn("chr1\t3\t4", text)
            self.assertNotIn("chr1\t1\t2", text)
            self.assertNotIn("chrX", text)


if __name__ == "__main__":
    unittest.main()
