"""Tests for chromosome-level parent-child methylation concordance."""

from __future__ import annotations

import gzip
import tempfile
import unittest
from pathlib import Path

from summarize_transmission_methylation import TransmissionSummaryError, main, summarize


HEADER = (
    "#chrom\tstart_cpg\tend_cpg\tfounder_haplotype_pat\tfounder_haplotype_mat\t"
    "methylation_level_pat_model\tmethylation_level_mat_model\t"
    "cpg_is_within_mismatch_window\n"
)


def _write_bed(path: Path, rows: list[str]) -> None:
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write(HEADER)
        for row in rows:
            fields = row.split("\t")
            fields.insert(2, str(int(fields[1]) + 2))
            handle.write("\t".join(fields) + "\n")


class TransmissionMethylationSummaryTests(unittest.TestCase):
    def test_transmitted_haplotype_metrics_are_paired_by_cpg(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            ped = root / "family.ped"
            ped.write_text(
                "FAM\tFATHER\t0\t0\t1\t0\n"
                "FAM\tMOTHER\t0\t0\t2\t0\n"
                "FAM\tCHILD\tFATHER\tMOTHER\t1\t0\n",
                encoding="utf-8",
            )
            father = root / "FATHER.dna-methylation.all-cpgs.bed.gz"
            child = root / "CHILD.dna-methylation.all-cpgs.bed.gz"
            _write_bed(
                father,
                [
                    "chr1\t0\tA\tB\t0.20\t0.80\tfalse",
                    "chr1\t1\tA\tB\t0.80\t0.10\tfalse",
                    "chr1\t2\tA\tB\t0.30\t0.70\tfalse",
                    "chr1\t3\tA\tB\t0.40\t0.60\tfalse",
                    "chr2\t0\tA\tB\t0.50\t0.50\tfalse",
                ],
            )
            _write_bed(
                child,
                [
                    "chr1\t0\tA\tC\t0.25\t0.30\tfalse",
                    "chr1\t1\tA\tC\t0.50\t0.30\tfalse",
                    "chr1\t2\tA\tC\t0.90\t0.30\ttrue",
                    "chr1\t3\tA\tC\t.\t0.30\tfalse",
                    "chr2\t0\tA\tC\t0.50\t0.30\tfalse",
                ],
            )

            result = summarize(ped, [father, child])

            paternal = next(
                row
                for row in result["comparisons"]
                if row["relationship"] == "paternal" and row["chromosome"] == "chr1"
            )
            self.assertEqual(paternal["eligible_cpgs"], 4)
            self.assertEqual(paternal["mismatch_excluded_cpgs"], 1)
            self.assertEqual(paternal["paired_cpgs"], 2)
            self.assertAlmostEqual(paternal["callable_fraction"], 2 / 3)
            self.assertAlmostEqual(paternal["agreement"], 0.825)
            self.assertAlmostEqual(paternal["mean_difference"], -0.125)
            self.assertAlmostEqual(paternal["discordant_fraction"], 0.0)
            self.assertEqual(paternal["discordance_threshold"], 0.4)
            self.assertEqual(paternal["specificity_cpgs"], 2)
            self.assertAlmostEqual(paternal["inherited_specificity"], 0.3)
            self.assertFalse(paternal["sufficient_cpgs"])
            self.assertFalse(paternal["sufficient_specificity_cpgs"])
            self.assertEqual(len(result["edges"]), 2)
            self.assertTrue(result["edges"][0]["has_methylation_outputs"])
            self.assertFalse(result["edges"][1]["has_methylation_outputs"])
            self.assertTrue(
                any(
                    edge["parent_id"] == "MOTHER"
                    for edge in result["unavailable_edges"]
                )
            )

    def test_identical_duplicate_rows_are_collapsed(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            ped = root / "family.ped"
            ped.write_text(
                "FAM\tP\t0\t0\t1\t0\nFAM\tC\tP\t0\t1\t0\n",
                encoding="utf-8",
            )
            parent = root / "P.dna-methylation.all-cpgs.bed.gz"
            child = root / "C.dna-methylation.all-cpgs.bed.gz"
            duplicate = "chr1\t0\tA\tB\t0.2\t0.8\tfalse"
            _write_bed(parent, [duplicate, duplicate])
            _write_bed(child, ["chr1\t0\tA\tC\t0.2\t0.4\tfalse"])

            result = summarize(ped, [parent, child])

            self.assertEqual(result["comparisons"][0]["paired_cpgs"], 1)
            self.assertEqual(result["comparisons"][0]["agreement"], 1.0)

    def test_maternal_transmission_follows_founder_label_across_recombination(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            ped = root / "family.ped"
            ped.write_text(
                "FAM\tMOTHER\t0\t0\t2\t0\nFAM\tCHILD\t0\tMOTHER\t1\t0\n",
                encoding="utf-8",
            )
            mother = root / "MOTHER.dna-methylation.all-cpgs.bed.gz"
            child = root / "CHILD.dna-methylation.all-cpgs.bed.gz"
            _write_bed(
                mother,
                [
                    "chr1\t0\tX\tC\t0.90\t0.40\tfalse",
                    "chr1\t1\tD\tY\t0.20\t0.80\tfalse",
                    "chr1\t2\tD\tY\t0.20\t0.80\ttrue",
                ],
            )
            _write_bed(
                child,
                [
                    "chr1\t0\tA\tC\t0.10\t0.30\tfalse",
                    "chr1\t1\tA\tD\t0.10\t0.30\tfalse",
                    "chr1\t2\tA\tD\t0.10\t0.30\tfalse",
                ],
            )

            result = summarize(
                ped, [mother, child], minimum_paired_cpgs=2
            )

            maternal = result["comparisons"][0]
            self.assertEqual(maternal["relationship"], "maternal")
            self.assertEqual(maternal["paired_cpgs"], 2)
            self.assertEqual(maternal["mismatch_excluded_cpgs"], 1)
            self.assertTrue(maternal["sufficient_cpgs"])
            self.assertTrue(maternal["sufficient_specificity_cpgs"])
            self.assertAlmostEqual(maternal["agreement"], 0.9)
            self.assertAlmostEqual(maternal["mean_difference"], 0.0)
            self.assertAlmostEqual(maternal["inherited_specificity"], 0.45)
            self.assertEqual(
                result["sources"]["all_cpgs"],
                sorted([mother.name, child.name]),
            )

    def test_conflicting_duplicate_values_are_excluded_as_ambiguous(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            ped = root / "family.ped"
            ped.write_text(
                "FAM\tP\t0\t0\t1\t0\nFAM\tC\tP\t0\t1\t0\n",
                encoding="utf-8",
            )
            parent = root / "P.dna-methylation.all-cpgs.bed.gz"
            child = root / "C.dna-methylation.all-cpgs.bed.gz"
            _write_bed(
                parent,
                [
                    "chr1\t0\tA\tB\t0.95\t0.634\tfalse",
                    "chr1\t0\tA\tB\t0.634\t0.95\tfalse",
                    "chr1\t1\tA\tB\t0.4\t0.6\tfalse",
                ],
            )
            _write_bed(
                child,
                [
                    "chr1\t0\tA\tC\t0.2\t0.4\tfalse",
                    "chr1\t1\tA\tC\t0.4\t0.3\tfalse",
                ],
            )

            result = summarize(ped, [parent, child], minimum_paired_cpgs=1)

            paternal = result["comparisons"][0]
            self.assertEqual(paternal["shared_cpgs"], 2)
            self.assertEqual(paternal["ambiguous_cpgs"], 1)
            self.assertEqual(paternal["eligible_cpgs"], 1)
            self.assertEqual(paternal["paired_cpgs"], 1)
            self.assertEqual(paternal["agreement"], 1.0)

    def test_conflicting_duplicate_labels_are_excluded_as_ambiguous(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            ped = root / "family.ped"
            ped.write_text(
                "FAM\tP\t0\t0\t1\t0\nFAM\tC\tP\t0\t1\t0\n",
                encoding="utf-8",
            )
            parent = root / "P.dna-methylation.all-cpgs.bed.gz"
            child = root / "C.dna-methylation.all-cpgs.bed.gz"
            _write_bed(
                parent,
                [
                    "chr1\t0\tA\tB\t0.2\t0.8\tfalse",
                    "chr1\t0\tB\tA\t0.8\t0.2\tfalse",
                ],
            )
            _write_bed(child, ["chr1\t0\tA\tC\t0.2\t0.4\tfalse"])

            result = summarize(ped, [parent, child])

            paternal = result["comparisons"][0]
            self.assertEqual(paternal["ambiguous_cpgs"], 1)
            self.assertEqual(paternal["eligible_cpgs"], 0)
            self.assertEqual(paternal["paired_cpgs"], 0)

    def test_summarize_emits_every_parent_child_edge_with_both_beds(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            ped = root / "family.ped"
            ped.write_text(
                "FAM\tGF\t0\t0\t1\t0\n"
                "FAM\tGM\t0\t0\t2\t0\n"
                "FAM\tPARENT\tGF\tGM\t1\t0\n"
                "FAM\tCHILD\tPARENT\t0\t1\t0\n",
                encoding="utf-8",
            )
            beds = {}
            for sample, pat, mat in (
                ("GF", "A", "B"),
                ("GM", "C", "D"),
                ("PARENT", "A", "C"),
                ("CHILD", "A", "X"),
            ):
                path = root / f"{sample}.dna-methylation.all-cpgs.bed.gz"
                _write_bed(path, [f"chr1\t0\t{pat}\t{mat}\t0.20\t0.80\tfalse"])
                beds[sample] = path

            result = summarize(ped, list(beds.values()), minimum_paired_cpgs=1)

            pairs = {
                (row["parent_id"], row["child_id"], row["relationship"])
                for row in result["comparisons"]
            }
            self.assertEqual(
                pairs,
                {
                    ("GF", "PARENT", "paternal"),
                    ("GM", "PARENT", "maternal"),
                    ("PARENT", "CHILD", "paternal"),
                },
            )
            self.assertEqual(len(result["edges"]), 3)
            self.assertEqual(result["unavailable_edges"], [])
            self.assertEqual(
                result["samples"], ["CHILD", "GF", "GM", "PARENT"]
            )

    def test_cli_keeps_repeated_all_cpgs_flags(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            ped = root / "family.ped"
            ped.write_text(
                "FAM\tP\t0\t0\t1\t0\nFAM\tC\tP\t0\t1\t0\n",
                encoding="utf-8",
            )
            parent = root / "P.dna-methylation.all-cpgs.bed.gz"
            child = root / "C.dna-methylation.all-cpgs.bed.gz"
            _write_bed(parent, ["chr1\t0\tA\tB\t0.20\t0.80\tfalse"])
            _write_bed(child, ["chr1\t0\tA\tC\t0.20\t0.40\tfalse"])
            output = root / "transmission-qc-summary.json"

            status = main(
                [
                    "--ped",
                    str(ped),
                    "--all-cpgs",
                    str(parent),
                    "--all-cpgs",
                    str(child),
                    "--minimum-paired-cpgs",
                    "1",
                    "--output",
                    str(output),
                ]
            )

            self.assertEqual(status, 0)
            document = output.read_text(encoding="utf-8")
            self.assertIn('"parent_id":"P"', document)
            self.assertIn('"child_id":"C"', document)
            self.assertIn("P.dna-methylation.all-cpgs.bed.gz", document)
            self.assertIn("C.dna-methylation.all-cpgs.bed.gz", document)


if __name__ == "__main__":
    unittest.main()
