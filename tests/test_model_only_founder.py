"""Model-only founder-phasing and FAI-backed BigWig tests."""

from __future__ import annotations

import sys
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

import polars as pl
import pyBigWig

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src" / "util"))

from hap_map_pedigree import get_hap_map  # noqa: E402
from phase_meth_to_founder_haps import (  # noqa: E402
    combine_count_and_model_based_methylation_levels,
    phase_meth_to_founder_haps,
    write_bigwig,
)


class ModelOnlyFounderTests(unittest.TestCase):
    def test_model_only_keeps_typed_null_count_columns(self) -> None:
        model = pl.DataFrame(
            {
                "chrom": ["chr1"],
                "start": [1],
                "end": [2],
                "methylation_level_pat": [0.75],
                "methylation_level_mat": [0.25],
            }
        )
        combined = combine_count_and_model_based_methylation_levels(None, model)
        self.assertEqual(combined["methylation_level_pat_model"].to_list(), [0.75])
        self.assertEqual(combined["methylation_level_mat_model"].to_list(), [0.25])
        self.assertEqual(combined["methylation_level_pat_count"].dtype, pl.Float64)
        self.assertEqual(combined["methylation_level_mat_count"].dtype, pl.Float64)
        self.assertEqual(combined["methylation_level_pat_count"].null_count(), 1)

    def test_no_hap_map_retains_cpg_with_null_founder_fields(self) -> None:
        methylation = pl.DataFrame(
            {
                "chrom": ["chr1"],
                "start": [1],
                "end": [2],
                "total_read_count_hap1": [12],
                "methylation_level_hap1": [0.8],
                "total_read_count_hap2": [11],
                "methylation_level_hap2": [0.2],
            }
        )
        hap_map = pl.DataFrame(
            schema={
                "chrom": pl.String,
                "start": pl.Int64,
                "end": pl.Int64,
                "paternal_haplotype": pl.String,
                "maternal_haplotype": pl.String,
                "haplotype_concordance": pl.Float64,
                "num_het_SNVs": pl.Int64,
            }
        )
        phased = phase_meth_to_founder_haps(methylation, hap_map)
        self.assertEqual(len(phased), 1)
        self.assertEqual(phased["methylation_level_pat"].null_count(), 1)
        self.assertEqual(phased["founder_haplotype_mat"].null_count(), 1)

    def test_empty_phasing_returns_stable_empty_outputs(self) -> None:
        phasing = pl.DataFrame(
            schema={
                "chrom": pl.String,
                "start": pl.Int64,
                "end": pl.Int64,
                "REF": pl.String,
                "ALT": pl.String,
            }
        )
        hap_map, sites, mismatches = get_hap_map(phasing)
        self.assertTrue(hap_map.is_empty())
        self.assertTrue(sites.is_empty())
        self.assertTrue(mismatches.is_empty())
        self.assertEqual(
            hap_map.columns,
            [
                "chrom",
                "start",
                "end",
                "paternal_haplotype",
                "maternal_haplotype",
                "haplotype_concordance",
                "num_het_SNVs",
            ],
        )

    def test_bigwig_uses_fai_lengths_and_exact_values(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fai = root / "reference.fa.fai"
            fai.write_text("chr1\t20\t6\t20\t21\n", encoding="utf-8")
            data = pl.DataFrame(
                {
                    "chrom": ["chr1", "chr1"],
                    "start": [1, 5],
                    "end": [2, 6],
                    "methylation_level_pat_model": [0.75, None],
                },
                schema_overrides={"methylation_level_pat_model": pl.Float64},
            )
            write_bigwig(
                data,
                "CHILD",
                "pat",
                "model",
                root,
                reference_fai=fai,
                reference_name="GRCh38",
            )
            path = root / "CHILD.dna-methylation.pat.model.GRCh38.bw"
            with pyBigWig.open(str(path)) as bigwig:
                self.assertEqual(bigwig.chroms(), {"chr1": 20})
                self.assertAlmostEqual(bigwig.values("chr1", 1, 2)[0], 0.75)
                self.assertTrue(bigwig.values("chr1", 5, 6)[0] != bigwig.values("chr1", 5, 6)[0])

    @unittest.skipUnless(
        shutil.which("bedGraphToBigWig"),
        "legacy bedGraphToBigWig is required for the migration parity gate",
    )
    def test_pybigwig_is_value_equivalent_to_legacy_generator(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fai = root / "reference.fa.fai"
            fai.write_text("chr1\t20\t6\t20\t21\nchr2\t10\t32\t10\t11\n", encoding="utf-8")
            chromsizes = root / "chrom.sizes"
            chromsizes.write_text("chr1\t20\nchr2\t10\n", encoding="utf-8")
            bedgraph = root / "legacy.bedGraph"
            bedgraph.write_text("chr1\t1\t2\t0.75\nchr2\t4\t5\t0.125\n", encoding="utf-8")
            legacy = root / "legacy.bw"
            subprocess.run(
                ["bedGraphToBigWig", str(bedgraph), str(chromsizes), str(legacy)],
                check=True,
            )
            data = pl.DataFrame(
                {
                    "chrom": ["chr1", "chr2", "chr2"],
                    "start": [1, 4, 7],
                    "end": [2, 5, 8],
                    "methylation_level_pat_model": [0.75, 0.125, None],
                },
                schema_overrides={"methylation_level_pat_model": pl.Float64},
            )
            write_bigwig(
                data,
                "SAMPLE",
                "pat",
                "model",
                root,
                reference_fai=fai,
                reference_name="GRCh38",
            )
            generated = root / "SAMPLE.dna-methylation.pat.model.GRCh38.bw"
            with pyBigWig.open(str(legacy)) as expected, pyBigWig.open(str(generated)) as observed:
                self.assertEqual(observed.chroms(), expected.chroms())
                for chrom, length in expected.chroms().items():
                    expected_values = expected.values(chrom, 0, length)
                    observed_values = observed.values(chrom, 0, length)
                    for expected_value, observed_value in zip(expected_values, observed_values):
                        if expected_value != expected_value:
                            self.assertTrue(observed_value != observed_value)
                        else:
                            self.assertAlmostEqual(observed_value, expected_value, places=6)


if __name__ == "__main__":
    unittest.main()
