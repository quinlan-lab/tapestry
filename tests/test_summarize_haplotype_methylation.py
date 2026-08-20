"""Tests for compact inheritance-row methylation summaries."""

from __future__ import annotations

import gzip
import tempfile
import unittest
from pathlib import Path

from summarize_haplotype_methylation import summarize


class HaplotypeMethylationSummaryTests(unittest.TestCase):
    def test_streams_phased_values_and_preserves_no_data_rows(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            iht = root / "FAM.iht.sorted.txt"
            iht.write_text(
                "#chrom start end FATHER MOTHER CHILD marker_count len markers\n"
                "chr1 0 10 A/B C/D A|C 2 10 a\n"
                "chr1 10 20 A/B C/D A|D 2 10 b\n"
                "chr1 20 30 A/B C/D A/C 2 10 c\n",
                encoding="utf-8",
            )
            all_cpgs = root / "CHILD.all-cpgs.bed.gz"
            header = (
                "#chrom\tstart_cpg\tmethylation_level_model\t"
                "founder_haplotype_pat\tfounder_haplotype_mat\t"
                "methylation_level_pat_model\tmethylation_level_mat_model\t"
                "cpg_is_within_mismatch_window\tcpg_is_allele_specific\n"
            )
            records = (
                "chr1\t1\t0.5\tA\tC\t0.8\t0.2\ttrue\tfalse\n"
                "chr1\t5\t0.5\tA\tC\t0.6\t0.4\tfalse\tfalse\n"
                "chr1\t12\t0.7\tA\tD\t0.7\t.\tfalse\ttrue\n"
                "chr1\t25\t0.9\tA\tC\t0.9\t0.1\tfalse\tfalse\n"
            )
            with gzip.open(all_cpgs, "wt", encoding="utf-8") as handle:
                handle.write(header + records)

            result = summarize(iht, all_cpgs, "CHILD")
            first, second, unphased = result["chromosomes"]["chr1"]

            self.assertEqual(first[:4], [0, 10, "A", "C"])
            self.assertEqual(first[4:8], [1.4, 2, 0.6, 2])
            self.assertEqual(first[10:15], [0.8, 2, 1, 0, 2])
            self.assertEqual(second[4:8], [0.7, 1, 0.0, 0])
            self.assertEqual(second[13:15], [1, 1])
            self.assertEqual(unphased[2:4], [None, None])
            self.assertEqual(unphased[4:8], [0.0, 0, 0.0, 0])
            self.assertEqual(unphased[8:10], [0.9, 1])


if __name__ == "__main__":
    unittest.main()
