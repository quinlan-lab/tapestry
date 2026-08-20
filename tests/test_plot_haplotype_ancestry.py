"""Tests for the standalone interactive pedigree haplotype visualization."""

from __future__ import annotations

import tempfile
import unittest
import json
from pathlib import Path

from plot_haplotype_ancestry import build_payload, write_visualization


class HaplotypeAncestryPlotTests(unittest.TestCase):
    def test_payload_merges_runs_and_preserves_pedigree_ancestry(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            ped = root / "family.ped"
            ped.write_text(
                "FAM\tGPF\tNA\tNA\t1\t0\n"
                "FAM\tGPM\tNA\tNA\t2\t0\n"
                "FAM\tMGF\tNA\tNA\t1\t0\n"
                "FAM\tMGM\tNA\tNA\t2\t0\n"
                "FAM\tFATHER\tGPF\tGPM\t1\t0\n"
                "FAM\tMOTHER\tMGF\tMGM\t2\t0\n"
                "FAM\tCHILD\tFATHER\tMOTHER\t1\t0\n",
                encoding="utf-8",
            )
            iht = root / "FAM.iht.sorted.txt"
            iht.write_text(
                "#chrom start end GPF GPM MGF MGM FATHER MOTHER CHILD marker_count len markers\n"
                "chr2 0 10 A/B C/D E/F G/H A|C E|G A|E 4 10 a\n"
                "chr2 10 20 A/B C/D E/F G/H A|C E|G A|E 4 10 b\n"
                "chr2 20 30 A/B C/D E/F G/H B|D F|H B|F 4 10 c\n"
                "chr10 0 5 A/B C/D E/F G/H A|C E|G A|E 4 5 d\n"
                "chr10 8 8 A/B C/D E/F G/H B|D F|H B|F 1 0 e\n",
                encoding="utf-8",
            )
            selected = root / "selected-samples.tsv"
            selected.write_text("sample_id\tfamily_id\nCHILD\tFAM\n", encoding="utf-8")

            payload = build_payload(ped, iht, selected)

            self.assertEqual(payload["initialSample"], "CHILD")
            self.assertEqual(payload["chromosomes"], ["chr2", "chr10"])
            self.assertEqual(payload["people"]["CHILD"]["generation"], 2)
            self.assertEqual(payload["labelNames"]["A"], "GPF hap1")
            self.assertEqual(payload["labelNames"]["H"], "MGM hap2")
            child_paternal = [
                block
                for block in payload["blocks"]
                if block["chrom"] == "chr2"
                and block["sample"] == "CHILD"
                and block["hap"] == 1
            ]
            self.assertEqual(
                [(block["start"], block["end"], block["label"]) for block in child_paternal],
                [(0, 20, "A"), (20, 30, "B")],
            )
            point_block = next(
                block
                for block in payload["blocks"]
                if block["chrom"] == "chr10"
                and block["sample"] == "CHILD"
                and block["hap"] == 1
                and block["start"] == 8
            )
            self.assertEqual((point_block["start"], point_block["end"]), (8, 9))

    def test_generation_labels_align_married_in_founders_with_their_mates(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            ped = root / "family.ped"
            ped.write_text(
                "FAM\tFATHER\tNA\tNA\t1\t0\n"
                "FAM\tMOTHER\tNA\tNA\t2\t0\n"
                "FAM\tCHILD\tFATHER\tMOTHER\t1\t0\n"
                "FAM\tSPOUSE\tNA\tNA\t2\t0\n"
                "FAM\tGRANDCHILD\tCHILD\tSPOUSE\t1\t0\n",
                encoding="utf-8",
            )
            iht = root / "FAM.iht.sorted.txt"
            iht.write_text(
                "#chrom start end FATHER MOTHER CHILD SPOUSE GRANDCHILD marker_count len markers\n"
                "chr1 0 20 A/B C/D A|C E/F A|E 4 20 a\n",
                encoding="utf-8",
            )

            payload = build_payload(ped, iht)

            self.assertEqual(payload["people"]["FATHER"]["generation"], 0)
            self.assertEqual(payload["people"]["CHILD"]["generation"], 1)
            self.assertEqual(payload["people"]["SPOUSE"]["generation"], 1)
            self.assertEqual(payload["people"]["GRANDCHILD"]["generation"], 2)

    def test_bundle_is_offline_lazy_and_exposes_interactions(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            ped = root / "family.ped"
            ped.write_text(
                "FAM\tFATHER\tNA\tNA\t1\t0\n"
                "FAM\tMOTHER\tNA\tNA\t2\t0\n"
                "FAM\tCHILD\tFATHER\tMOTHER\t1\t0\n",
                encoding="utf-8",
            )
            iht = root / "FAM.iht.sorted.txt"
            iht.write_text(
                "#chrom start end FATHER MOTHER CHILD marker_count len markers\n"
                "chr1 0 20 A/B C/D A|C 4 20 a\n",
                encoding="utf-8",
            )
            output = root / "haplotype-ancestry"
            methylation = root / "CHILD.methylation-summary.json"
            methylation.write_text(
                json.dumps(
                    {
                        "schema_version": 1,
                        "sample_id": "CHILD",
                        "chromosomes": {
                            "chr1": [
                                [0, 20, "A", "C", 1.5, 2, 0.5, 2, 2.0, 3, 1.0, 2, 1, 0, 3]
                            ]
                        },
                    }
                ),
                encoding="utf-8",
            )
            transmission_qc = root / "transmission-qc-summary.json"
            transmission_qc.write_text(
                json.dumps(
                    {
                        "schema_version": 1,
                        "measurement": "model-based methylation",
                        "discordance_threshold": 0.4,
                        "minimum_paired_cpgs": 1,
                        "samples": ["FATHER", "CHILD"],
                        "edges": [
                            {
                                "family_id": "FAM",
                                "parent_id": "FATHER",
                                "child_id": "CHILD",
                                "relationship": "paternal",
                                "has_methylation_outputs": True,
                            }
                        ],
                        "comparisons": [
                            {
                                "family_id": "FAM",
                                "parent_id": "FATHER",
                                "child_id": "CHILD",
                                "relationship": "paternal",
                                "chromosome": "chr1",
                                "eligible_cpgs": 3,
                                "mismatch_excluded_cpgs": 1,
                                "evaluated_cpgs": 2,
                                "paired_cpgs": 2,
                                "callable_fraction": 1.0,
                                "agreement": 0.875,
                                "mean_difference": 0.025,
                                "discordant_fraction": 0.0,
                                "specificity_cpgs": 2,
                                "inherited_specificity": 0.3,
                            }
                        ],
                        "unavailable_edges": [],
                    }
                ),
                encoding="utf-8",
            )

            write_visualization(
                build_payload(ped, iht), output, [methylation], transmission_qc
            )
            document = (output / "index.html").read_text(encoding="utf-8")
            metadata = (output / "metadata.js").read_text(encoding="utf-8")
            qc_document = (output / "transmission-qc.html").read_text(
                encoding="utf-8"
            )
            qc_data = (output / "transmission-qc.js").read_text(encoding="utf-8")

            self.assertIn("Plotly.newPlot", document)
            self.assertIn("plotly_click", document)
            self.assertIn("loadShard", document)
            self.assertIn("? `${generation} ${sample}`", document)
            self.assertIn("`${generation} ${sample} (inheritance map only)`", document)
            self.assertNotIn("methylation output)", document)
            self.assertNotIn("(pipeline sample)", document)
            self.assertIn("F${DATA.people[sample].generation}", document)
            self.assertIn("Blank spans have no inheritance-map block", document)
            self.assertIn("loadMethylationShard", document)
            self.assertIn("Mean model methylation", document)
            self.assertIn("Transmission QC", document)
            self.assertIn("Inherited-haplotype agreement", qc_document)
            self.assertIn("Parent–child pair", qc_document)
            self.assertIn("All pairs", qc_document)
            self.assertIn("Missing transmission comparisons", qc_document)
            self.assertIn("insufficient contributing CpGs", qc_document)
            self.assertIn("plotly_click", qc_document)
            self.assertIn('"agreement":0.875', qc_data)
            self.assertIn('"labelNames":{"A":"FATHER hap1"', metadata)
            self.assertIn('"methylationSamples":["CHILD"]', metadata)
            self.assertIn('src="plotly.min.js"', document)
            self.assertIn('"initialSample":"CHILD"', metadata)
            self.assertEqual(len(list((output / "data" / "samples").glob("*.js"))), 3)
            self.assertEqual(
                len(list((output / "data" / "chromosomes").glob("*.js"))), 1
            )
            self.assertEqual(
                len(
                    list(
                        (output / "data" / "methylation" / "chromosomes").glob(
                            "*.js"
                        )
                    )
                ),
                1,
            )
            self.assertGreater((output / "plotly.min.js").stat().st_size, 1_000_000)
            self.assertLess((output / "index.html").stat().st_size, 100_000)


if __name__ == "__main__":
    unittest.main()
