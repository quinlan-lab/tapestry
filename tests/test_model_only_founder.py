"""Model-only founder-phasing and FAI-backed BigWig tests."""

from __future__ import annotations

import logging
import sys
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import polars as pl
from polars.testing import assert_frame_equal
import pyBigWig
import pysam

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src" / "util"))

from hap_map_pedigree import get_hap_map  # noqa: E402
from get_meth_hap1_hap2 import read_meth_hap1_hap2  # noqa: E402
from phase_meth_to_founder_haps import (  # noqa: E402
    combine_count_and_model_based_methylation_levels,
    phase_methylation_by_chromosome,
    phase_meth_to_founder_haps,
    write_bigwig,
)
from util.read_data import read_bed_and_header  # noqa: E402


def _methylation_bed(path: Path, mode: str, rows: list[str]) -> Path:
    plain = path.with_suffix("")
    plain.write_text(
        f"##pileup-mode={mode}\n"
        "##min-coverage=10\n"
        "#chrom\tbegin\tend\tmod_score\ttype\tcov\n"
        + "\n".join(rows)
        + "\n",
        encoding="utf-8",
    )
    pysam.tabix_compress(str(plain), str(path), force=True)
    pysam.tabix_index(str(path), preset="bed", force=True)
    plain.unlink()
    return path


class ModelOnlyFounderTests(unittest.TestCase):
    def test_chromosome_batches_match_whole_genome_in_count_and_model_modes(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            model_hap1 = _methylation_bed(
                root / "model-hap1.bed.gz",
                "model",
                [
                    "chr1\t1\t2\t80\tCG\t12",
                    "chr2\t3\t4\t70\tCG\t13",
                    "chr2\t6\t7\t60\tCG\t14",
                ],
            )
            model_hap2 = _methylation_bed(
                root / "model-hap2.bed.gz",
                "model",
                [
                    "chr1\t1\t2\t20\tCG\t11",
                    "chr1\t5\t6\t40\tCG\t10",
                    "chr2\t3\t4\t30\tCG\t15",
                ],
            )
            count_hap1 = _methylation_bed(
                root / "count-hap1.bed.gz",
                "count",
                [
                    "chr1\t1\t2\t75\tCG\t12",
                    "chr2\t3\t4\t65\tCG\t13",
                    "chr2\t6\t7\t55\tCG\t14",
                ],
            )
            count_hap2 = _methylation_bed(
                root / "count-hap2.bed.gz",
                "count",
                [
                    "chr1\t1\t2\t25\tCG\t11",
                    "chr1\t5\t6\t45\tCG\t10",
                    "chr2\t3\t4\t35\tCG\t15",
                ],
            )
            hap_map = pl.DataFrame(
                {
                    "chrom": ["chr1", "chr2"],
                    "start": [0, 0],
                    "end": [10, 10],
                    "paternal_haplotype": ["A_hap1", "B_hap2"],
                    "maternal_haplotype": ["C_hap2", "D_hap1"],
                    "haplotype_concordance": [1.0, 0.9],
                    "num_het_SNVs": [2, 3],
                }
            )
            expected = combine_count_and_model_based_methylation_levels(
                phase_meth_to_founder_haps(
                    read_meth_hap1_hap2("count", count_hap1, count_hap2), hap_map
                ),
                phase_meth_to_founder_haps(
                    read_meth_hap1_hap2("model", model_hap1, model_hap2), hap_map
                ),
            )
            fai = root / "reference.fa.fai"
            fai.write_text(
                "chr1\t20\t6\t20\t21\n"
                "chr2\t20\t32\t20\t21\n"
                "chr3\t20\t58\t20\t21\n",
                encoding="utf-8",
            )
            output = root / "output"
            output.mkdir()
            modes, rows = phase_methylation_by_chromosome(
                uid="CHILD",
                regions=["chr2", "chr3", "chr1"],
                df_hap_map=hap_map,
                model_hap1=str(model_hap1),
                model_hap2=str(model_hap2),
                count_hap1=str(count_hap1),
                count_hap2=str(count_hap2),
                output_dir=output,
                reference_fai=str(fai),
                reference_name="GRCh38",
                write_bigwigs=True,
                logger=logging.getLogger(__name__),
            )
            observed = read_bed_and_header(output / "CHILD.dna-methylation.bed")
            observed = observed.cast(expected.schema).select(expected.columns)

            self.assertEqual(modes, ["count", "model"])
            self.assertEqual(rows, len(expected))
            assert_frame_equal(
                observed,
                expected,
                check_exact=False,
                rel_tol=1e-12,
                abs_tol=1e-12,
            )
            self.assertEqual(
                observed["chrom"].to_list(), ["chr1", "chr1", "chr2", "chr2"]
            )
            with pyBigWig.open(
                str(output / "CHILD.dna-methylation.pat.model.GRCh38.bw")
            ) as bigwig:
                self.assertAlmostEqual(bigwig.values("chr1", 1, 2)[0], 0.8)
                self.assertAlmostEqual(bigwig.values("chr2", 3, 4)[0], 0.3)
            with pyBigWig.open(
                str(output / "CHILD.dna-methylation.pat.count.GRCh38.bw")
            ) as bigwig:
                self.assertAlmostEqual(bigwig.values("chr1", 1, 2)[0], 0.75)
                self.assertAlmostEqual(bigwig.values("chr2", 3, 4)[0], 0.35)

            no_bigwig_output = root / "no-bigwig-output"
            no_bigwig_output.mkdir()
            with patch(
                "phase_meth_to_founder_haps.bf.fetch_chromsizes",
                side_effect=AssertionError("chromsizes must not be fetched"),
            ):
                no_bigwig_modes, no_bigwig_rows = phase_methylation_by_chromosome(
                    uid="CHILD",
                    regions=["chr2", "chr3", "chr1"],
                    df_hap_map=hap_map,
                    model_hap1=str(model_hap1),
                    model_hap2=str(model_hap2),
                    count_hap1=str(count_hap1),
                    count_hap2=str(count_hap2),
                    output_dir=no_bigwig_output,
                    reference_fai=None,
                    reference_name="GRCh38",
                    write_bigwigs=False,
                    logger=logging.getLogger(__name__),
                )
            no_bigwig_observed = read_bed_and_header(
                no_bigwig_output / "CHILD.dna-methylation.bed"
            )
            no_bigwig_observed = no_bigwig_observed.cast(expected.schema).select(
                expected.columns
            )
            self.assertEqual(no_bigwig_modes, ["count", "model"])
            self.assertEqual(no_bigwig_rows, len(expected))
            self.assertEqual(list(no_bigwig_output.glob("*.bw")), [])
            assert_frame_equal(
                no_bigwig_observed,
                expected,
                check_exact=False,
                rel_tol=1e-12,
                abs_tol=1e-12,
            )

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
                    "chrom": ["chr1", "chr1", "chr1"],
                    "start": [1, 1, 5],
                    "end": [2, 2, 6],
                    "methylation_level_pat_model": [0.75, 0.75, None],
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

    def test_bigwig_rejects_conflicting_values_at_one_interval(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fai = root / "reference.fa.fai"
            fai.write_text("chr1\t20\t6\t20\t21\n", encoding="utf-8")
            data = pl.DataFrame(
                {
                    "chrom": ["chr1", "chr1"],
                    "start": [1, 1],
                    "end": [2, 2],
                    "methylation_level_pat_model": [0.75, 0.25],
                }
            )
            with self.assertRaisesRegex(
                ValueError,
                r"Conflicting pat model BigWig values.*chr1:1-2",
            ):
                write_bigwig(
                    data,
                    "CHILD",
                    "pat",
                    "model",
                    root,
                    reference_fai=fai,
                    reference_name="GRCh38",
                )

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
