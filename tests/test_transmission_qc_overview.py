"""Tests for genome-wide transmission QC scoring."""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from plot_transmission_qc import (
    compute_overview,
    driver_chromosomes,
    edge_key,
    global_metric,
    robust_outliers,
    write_qc_page,
)


def _row(
    parent: str,
    child: str,
    relationship: str,
    chromosome: str,
    *,
    agreement: float | None = 0.99,
    mean_difference: float = 0.0,
    discordant_fraction: float = 0.0,
    callable_fraction: float = 1.0,
    inherited_specificity: float = 0.1,
    paired_cpgs: int = 1000,
    evaluated_cpgs: int = 1000,
    eligible_cpgs: int = 1000,
    mismatch_excluded_cpgs: int = 0,
    specificity_cpgs: int = 1000,
) -> dict[str, object]:
    return {
        "parent_id": parent,
        "child_id": child,
        "relationship": relationship,
        "chromosome": chromosome,
        "agreement": agreement,
        "mean_difference": mean_difference,
        "discordant_fraction": discordant_fraction,
        "callable_fraction": callable_fraction,
        "inherited_specificity": inherited_specificity,
        "paired_cpgs": paired_cpgs,
        "evaluated_cpgs": evaluated_cpgs,
        "eligible_cpgs": eligible_cpgs,
        "mismatch_excluded_cpgs": mismatch_excluded_cpgs,
        "specificity_cpgs": specificity_cpgs,
    }


class TransmissionQcOverviewTests(unittest.TestCase):
    def test_scores_are_blank_below_four_pairs(self) -> None:
        stats = robust_outliers(
            [("a", 0.99), ("b", 0.98), ("c", 0.50)], "agreement"
        )
        self.assertIsNone(stats["c"]["score"])
        self.assertAlmostEqual(stats["c"]["center"], 0.98)

    def test_low_agreement_is_concerning_and_high_agreement_is_not(self) -> None:
        stats = robust_outliers(
            [
                ("a", 0.990),
                ("b", 0.991),
                ("c", 0.989),
                ("d", 0.992),
                ("bad", 0.80),
            ],
            "agreement",
        )
        self.assertGreater(stats["bad"]["score"], 4)
        self.assertEqual(stats["a"]["score"], 0.0)

    def test_yield_scores_fold_change_not_raw_count_gaps(self) -> None:
        typical = [
            ("a", 1_000_000.0),
            ("b", 1_100_000.0),
            ("c", 900_000.0),
            ("d", 1_050_000.0),
        ]
        slight = robust_outliers([*typical, ("low", 800_000.0)], "paired_cpgs")
        severe = robust_outliers([*typical, ("tiny", 1_000.0)], "paired_cpgs")
        self.assertLess(slight["low"]["score"], severe["tiny"]["score"])
        self.assertGreater(severe["tiny"]["score"], 4)

    def test_concordance_ignores_chromosomes_below_minimum_paired_cpgs(self) -> None:
        rows = [
            _row("P", "C", "paternal", "chr1", agreement=0.99, paired_cpgs=1000),
            _row("P", "C", "paternal", "chr2", agreement=0.10, paired_cpgs=10),
        ]
        self.assertAlmostEqual(
            global_metric(rows, ["chr1", "chr2"], "agreement", 100), 0.99
        )
        self.assertAlmostEqual(
            global_metric(rows, ["chr1", "chr2"], "callable_fraction", 100),
            1010 / 2000,
        )

    def test_drivers_rank_weighted_pull_not_local_mad(self) -> None:
        rows = [
            _row("P", "C", "paternal", "chr1", agreement=0.90, paired_cpgs=10_000),
            _row("P", "C", "paternal", "chr2", agreement=0.10, paired_cpgs=100),
        ]
        drivers = driver_chromosomes(
            rows,
            ["chr1", "chr2"],
            "agreement",
            100,
            cohort_centers={"chr1": 0.99, "chr2": 0.99},
        )
        self.assertEqual([item["chromosome"] for item in drivers], ["chr1", "chr2"])
        self.assertGreater(drivers[0]["contribution"], drivers[1]["contribution"])

    def test_completeness_drivers_keep_low_yield_chromosomes(self) -> None:
        rows = [
            _row(
                "P",
                "C",
                "paternal",
                "chr1",
                callable_fraction=1.0,
                paired_cpgs=1000,
                evaluated_cpgs=1000,
            ),
            _row(
                "P",
                "C",
                "paternal",
                "chr2",
                callable_fraction=0.01,
                paired_cpgs=10,
                evaluated_cpgs=1000,
            ),
        ]
        drivers = driver_chromosomes(
            rows,
            ["chr1", "chr2"],
            "callable_fraction",
            100,
            cohort_centers={"chr1": 0.9, "chr2": 0.9},
        )
        self.assertEqual(drivers[0]["chromosome"], "chr2")
        self.assertEqual(len(drivers), 2)

    def test_paired_cpg_drivers_use_chromosome_count_deficits(self) -> None:
        rows = [
            _row("P", "C", "paternal", "chr1", paired_cpgs=1_000_000),
            _row("P", "C", "paternal", "chr2", paired_cpgs=1_000),
        ]
        drivers = driver_chromosomes(
            rows,
            ["chr1", "chr2"],
            "paired_cpgs",
            100,
            cohort_centers={"chr1": 1_000_000, "chr2": 100_000},
        )
        self.assertEqual(drivers[0]["chromosome"], "chr2")
        self.assertEqual(drivers[0]["contribution"], 99_000)

    def test_overview_splits_completeness_from_concordance(self) -> None:
        summary = {
            "minimum_paired_cpgs": 100,
            "edges": [
                {
                    "parent_id": f"P{index}",
                    "child_id": f"C{index}",
                    "relationship": "paternal",
                    "has_methylation_outputs": True,
                }
                for index in range(4)
            ],
            "comparisons": [
                _row(
                    f"P{index}",
                    f"C{index}",
                    "paternal",
                    "chr1",
                    agreement=0.99 if index < 3 else 0.70,
                    paired_cpgs=1000 if index < 3 else 200,
                    mismatch_excluded_cpgs=0 if index < 3 else 400,
                    eligible_cpgs=1000 if index < 3 else 600,
                    evaluated_cpgs=1000 if index < 3 else 200,
                    callable_fraction=1.0 if index < 3 else 1.0,
                )
                for index in range(4)
            ],
        }
        overview = compute_overview(summary, ["chr1", "chr2"])
        bad = edge_key("P3", "C3", "paternal")
        typical = edge_key("P0", "C0", "paternal")
        self.assertGreater(
            overview["metrics"]["agreement"]["edges"][bad]["score"], 1
        )
        self.assertEqual(
            overview["metrics"]["agreement"]["edges"][typical]["score"], 0.0
        )
        self.assertEqual(
            overview["metrics"]["missing_chromosomes"]["edges"][bad]["value"], 1
        )
        self.assertAlmostEqual(
            overview["metrics"]["mismatch_excluded_fraction"]["edges"][bad]["value"],
            400 / 600,
        )
        self.assertEqual(
            overview["metrics"]["agreement"]["edges"][bad]["drivers"][0]["chromosome"],
            "chr1",
        )

    def test_write_qc_page_embeds_overview_scores(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            summary_path = root / "transmission-qc-summary.json"
            summary = {
                "schema_version": 1,
                "measurement": "model-based methylation",
                "discordance_threshold": 0.4,
                "minimum_paired_cpgs": 100,
                "samples": ["P0", "C0", "P1", "C1", "P2", "C2", "P3", "C3"],
                "edges": [
                    {
                        "family_id": "FAM",
                        "parent_id": f"P{index}",
                        "child_id": f"C{index}",
                        "relationship": "paternal",
                        "has_methylation_outputs": True,
                    }
                    for index in range(4)
                ],
                "comparisons": [
                    _row(
                        f"P{index}",
                        f"C{index}",
                        "paternal",
                        "chr1",
                        agreement=0.99 if index < 3 else 0.70,
                    )
                    for index in range(4)
                ],
                "unavailable_edges": [],
            }
            summary_path.write_text(json.dumps(summary), encoding="utf-8")
            payload = {
                "samples": summary["samples"],
                "chromosomes": ["chr1"],
            }
            result = write_qc_page(root, payload, summary_path)
            bad = edge_key("P3", "C3", "paternal")
            self.assertIn("overview", result)
            self.assertGreater(
                result["overview"]["metrics"]["agreement"]["edges"][bad]["score"], 1
            )
            published = (root / "transmission-qc.js").read_text(encoding="utf-8")
            self.assertIn('"overview"', published)
