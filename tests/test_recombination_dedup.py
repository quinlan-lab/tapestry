"""Regression test for the triplicate / nested hap-map-block bug.

Background
----------
WhatsHap phase-block spans can nest: a short read-backed phase set sits inside
the genomic span of a chromosome-spanning pedMEC phase set. The old code
assigned each SNV to a phase block by genomic-range overlap (bf.overlap), so a
SNV in a nested region was assigned to BOTH the small and the big block. Done
independently for kid and parent, this fanned a single SNV out into several
rows and collapsed them into duplicate (triplicate) hap-map blocks, while also
letting the chromosome-spanning block absorb SNVs that belong to a nested phase
set.

The fix assigns each SNV to the phase set it was actually phased in, using the
per-sample PS (phase-set) FORMAT tag, so every SNV lands in exactly one
(kid_phase_set, parent_phase_set) group.

This test synthesises that nested scenario and asserts the invariant that
directly distinguishes the fixed code from the buggy code:

    sum(num_het_SNVs over hap-map blocks) == number of input phased het SNVs

Under the old fan-out this sum was inflated (each nested SNV counted 2-3x) and
the nested block appeared as duplicate rows. It also checks there are no
duplicate (chrom, start, end) rows.

It also runs the trio hap map through the shared
src/util/hap_map.py:write_hap_map_blocks writer exactly as
src/phase_meth_to_parent_haps.py does, so a change to that writer which breaks
the trio workflow (e.g. hardcoding the pedigree "num_het_SNVs" column) fails CI.

Synthetic genotypes are Mendelian-valid: both parents are het (0|1) at every
site, so any kid genotype/phasing is consistent with inheritance.

Run with the CI discovery command, or directly as a script:
    PYTHONPATH=src:src/util .venv/bin/python -m unittest discover -s tests -p 'test_*.py'
    .venv/bin/python tests/test_recombination_dedup.py
"""

import gzip
import sys
import tempfile
import unittest
from pathlib import Path

import polars as pl

_REPO = Path(__file__).resolve().parents[1]
for _path in (_REPO / "src" / "util", _REPO / "src"):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))

from phasing_trio import get_phase_blocks, get_all_phasing  # noqa: E402
from hap_map_trio import get_hap_map  # noqa: E402
from hap_map import write_hap_map_blocks  # noqa: E402
from report_whatshap_crossovers import COLUMNS, read_unresolved_spans  # noqa: E402

# Trio hap-map SNV-count column, as passed by src/phase_meth_to_parent_haps.py.
TRIO_NUM_HET_SNVS_COLUMN = "num_het_SNVs_in_parent"

KID, DAD, MOM = "NA_KID", "NA_DAD", "NA_MOM"

# Big chromosome-spanning phase set (PS=100) covers pos 100..900 for everyone.
# A small nested phase set (PS=500) covers pos 500..502 in the KID and DAD only;
# the MOM stays in her big block there. So sites 500-502 fall inside the span of
# BOTH the big and the small block for kid and dad -- the trigger for the old
# fan-out.
VCF_ROWS = [
    # pos, kid_PS, dad_PS, mom_PS
    (100, 100, 100, 100),
    (200, 100, 100, 100),
    (300, 100, 100, 100),
    (500, 500, 500, 100),
    (501, 500, 500, 100),
    (502, 500, 500, 100),
    (900, 100, 100, 100),
]

VCF_HEADER = """\
##fileformat=VCFv4.2
##contig=<ID=chr1>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set identifier">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{kid}\t{dad}\t{mom}
""".format(kid=KID, dad=DAD, mom=MOM)

# All samples het and phased everywhere (0|1), so every site is informative for
# both the paternal (dad-het) and maternal (mom-het) comparisons.
BLOCKS_TSV = """\
#sample\tchromosome\tphase_set\tfrom\tto\tvariants
{kid}\tchr1\t100\t100\t900\t5
{kid}\tchr1\t500\t500\t502\t3
{dad}\tchr1\t100\t100\t900\t5
{dad}\tchr1\t500\t500\t502\t3
{mom}\tchr1\t100\t100\t900\t7
""".format(kid=KID, dad=DAD, mom=MOM)


def _write_inputs(d: Path):
    vcf = VCF_HEADER
    for pos, kid_ps, dad_ps, mom_ps in VCF_ROWS:
        vcf += (
            f"chr1\t{pos}\t.\tA\tG\t.\t.\t.\tGT:PS"
            f"\t0|1:{kid_ps}\t0|1:{dad_ps}\t0|1:{mom_ps}\n"
        )
    vcf_path = d / "trio.vcf"
    vcf_path.write_text(vcf)
    blocks_path = d / "blocks.tsv"
    blocks_path.write_text(BLOCKS_TSV)
    return str(vcf_path), str(blocks_path)


def _build_trio_hap_map(d: Path):
    """Run the trio phasing + hap-map stages on the synthetic fixture."""
    d.mkdir(parents=True, exist_ok=True)
    vcf_path, blocks_path = _write_inputs(d)

    df_blocks_kid = get_phase_blocks(blocks_path, KID)
    df_blocks_dad = get_phase_blocks(blocks_path, DAD)
    df_blocks_mom = get_phase_blocks(blocks_path, MOM)

    df_kid_dad, df_kid_mom = get_all_phasing(
        vcf_path, KID, DAD, MOM,
        df_blocks_kid, df_blocks_dad, df_blocks_mom,
    )

    df_pat, df_mat, _, _ = get_hap_map(df_kid_dad, df_kid_mom)
    return df_kid_dad, df_kid_mom, df_pat, df_mat


def _read_bgzipped_bed(path: Path):
    with gzip.open(path, "rt") as f:
        return [
            line.rstrip("\n").split("\t")
            for line in f
            if line.strip() and not line.startswith("#")
        ]


class TrioHapMapNestedPhaseSetTests(unittest.TestCase):
    """Each phased het SNV lands in exactly one (kid PS, parent PS) block."""

    @classmethod
    def setUpClass(cls):
        with tempfile.TemporaryDirectory() as tmp:
            (
                cls.df_kid_dad,
                cls.df_kid_mom,
                cls.df_pat,
                cls.df_mat,
            ) = _build_trio_hap_map(Path(tmp))
        cls.hap_maps = {
            "paternal": (cls.df_pat, len(cls.df_kid_dad)),
            "maternal": (cls.df_mat, len(cls.df_kid_mom)),
        }

    def test_no_duplicate_hap_map_blocks(self):
        for parental, (df, _) in self.hap_maps.items():
            with self.subTest(parental=parental):
                dups = (
                    df.group_by(["chrom", "start", "end"])
                    .len()
                    .filter(pl.col("len") > 1)
                )
                self.assertEqual(
                    len(dups), 0,
                    f"{parental}: duplicate hap-map blocks:\n{dups}",
                )

    def test_num_het_SNVs_sum_equals_input_phased_het_SNVs(self):
        for parental, (df, n_input) in self.hap_maps.items():
            with self.subTest(parental=parental):
                self.assertGreater(n_input, 0)
                self.assertEqual(
                    df[TRIO_NUM_HET_SNVS_COLUMN].sum(), n_input,
                    f"{parental}: {TRIO_NUM_HET_SNVS_COLUMN} over blocks must "
                    "equal the number of input phased het SNVs "
                    "(fan-out double-counting)",
                )

    def test_nested_phase_set_is_its_own_paternal_block(self):
        # The nested set yields its OWN single block of 3 SNVs, and the big
        # block holds the remaining 4 -- not 7-in-the-big + a triplicate.
        self.assertEqual(
            sorted(self.df_pat[TRIO_NUM_HET_SNVS_COLUMN].to_list()), [3, 4],
        )


def _write_crossover_tsv(path, calls):
    with path.open("w") as handle:
        print(*COLUMNS, sep="\t", file=handle)
        for start, end, parent, nearest, clustered in calls:
            event = f"{KID}:chr1:{start}-{end}:{parent}"
            print(event, KID, "chr1", parent, start, end, 1, 0, 45, nearest, clustered,
                  sep="\t", file=handle)


class TrioCrossoverMatchingTests(unittest.TestCase):
    def test_independent_comparisons_split_nested_phase_pairs(self):
        positions = [100, 150, 200, 230, 240, 260, 280, 300, 340, 350, 360, 450,
                     500, 501, 502, 510, 599, 600, 630, 640, 650, 670, 680, 681, 800, 850, 900]
        vcf = VCF_HEADER
        for pos in positions:
            nested = 500 <= pos <= 502
            kid_ps = 500 if nested else 100
            dad_ps = 999 if pos == 510 else kid_ps  # No matching block in the father.
            pat = int(pos == 150 or (260 < pos < 600 and not nested))
            mat = int((pos < 340 and pos != 300) or nested)
            kid = f"{pat}|{mat}:{kid_ps}" if pat != mat else f"{pat}/{mat}:."
            vcf += f"chr1\t{pos}\t.\tA\tG\t.\t.\t.\tGT:PS\t{kid}\t0|1:{dad_ps}\t0|1:100\n"

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            vcf_path, blocks, calls = root / "trio.vcf", root / "blocks.tsv", root / "calls.tsv"
            vcf_path.write_text(vcf)
            blocks.write_text(BLOCKS_TSV)
            _write_crossover_tsv(calls, [
                (230, 260, "father", 80, 0), (340, 360, "mother", 80, 0),
                (600, 630, "father", 20, 1), (650, 680, "mother", 20, 1),
            ])
            spans = read_unresolved_spans(calls, KID, 20)
            self.assertEqual(spans, {
                "father": [("chr1", 229, 260), ("chr1", 599, 680)],
                "mother": [("chr1", 339, 360), ("chr1", 599, 680)],
            })
            kid_dad, kid_mom = get_all_phasing(
                str(vcf_path), KID, DAD, MOM,
                *(get_phase_blocks(str(blocks), sample) for sample in (KID, DAD, MOM)),
                require_phase_sets=True,
            )

        for paired, parent in [(kid_dad, "dad"), (kid_mom, "mom")]:
            expected_positions = [pos for pos in positions if parent == "mom" or pos != 510]
            self.assertEqual(paired["end"].to_list(), expected_positions)
            nested = paired.filter(pl.col("end").is_between(500, 502))
            self.assertEqual(nested["phase_block_id_kid"].to_list(), ["500"] * 3)
            self.assertEqual(nested[f"phase_block_id_{parent}"].to_list(),
                             ["500" if parent == "dad" else "100"] * 3)

        pat, mat, mismatch_pat, mismatch_mat = get_hap_map(kid_dad, kid_mom, unresolved_spans=spans)
        self.assertEqual(pat.select("start", "end", "paternal_haplotype", TRIO_NUM_HET_SNVS_COLUMN).rows(), [
            (99, 229, "A", 3), (260, 599, "B", 7), (499, 502, "A", 3), (680, 900, "A", 4),
        ])
        self.assertEqual(mat.select("start", "end", "maternal_haplotype", TRIO_NUM_HET_SNVS_COLUMN).rows(), [
            (99, 339, "D", 8), (360, 599, "C", 3), (499, 502, "D", 3), (680, 900, "C", 4),
        ])
        self.assertEqual(mismatch_pat["end"].to_list(), [150])
        self.assertEqual(mismatch_mat["end"].to_list(), [300])
        for parent, paired, hap_map in [("father", kid_dad, pat), ("mother", kid_mom, mat)]:
            eligible = [pos for pos in paired["end"]
                        if all(not start < pos <= end for _, start, end in spans[parent])]
            self.assertEqual(hap_map[TRIO_NUM_HET_SNVS_COLUMN].sum(), len(eligible))

        empty, _, mismatches, _ = get_hap_map(kid_dad, kid_mom, unresolved_spans={
            "father": [("chr1", 0, 900)],
        })
        self.assertTrue(empty.is_empty())
        self.assertEqual(empty.schema, pat.schema)
        self.assertTrue(mismatches.is_empty())
        self.assertEqual(mismatches.schema, mismatch_pat.schema)

    def test_transitive_cluster_uses_full_span_for_affected_parent(self):
        with tempfile.TemporaryDirectory() as tmp:
            calls = Path(tmp) / "calls.tsv"
            _write_crossover_tsv(calls, [
                (36, 40, "father", 10, 1), (1, 10, "father", 10, 1),
                (20, 26, "father", 10, 1),
            ])
            self.assertEqual(read_unresolved_spans(calls, KID, 10), {
                "father": [("chr1", 0, 40)], "mother": [],
            })
            with self.assertRaisesRegex(ValueError, "clustered flags disagree"):
                read_unresolved_spans(calls, KID, 9)
            _write_crossover_tsv(calls, [])
            self.assertEqual(read_unresolved_spans(calls, KID, 10), {"father": [], "mother": []})


class TrioHapMapBlocksWriterTests(unittest.TestCase):
    """The shared hap-map-blocks writer must keep working for the trio."""

    def test_writes_indexed_hap_map_blocks_with_trio_snv_counts(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            _, _, df_pat, df_mat = _build_trio_hap_map(d / "inputs")
            output_dir = d / "outputs"
            hap_maps = {"paternal": df_pat, "maternal": df_mat}

            for parental, df in hap_maps.items():
                write_hap_map_blocks(
                    df, KID, parental, str(output_dir),
                    num_het_SNVs_column=TRIO_NUM_HET_SNVS_COLUMN,
                )

            for parental, df in hap_maps.items():
                with self.subTest(parental=parental):
                    bed = output_dir / f"{KID}.hap-map-blocks.{parental}.bed.gz"
                    self.assertTrue(bed.is_file(), f"missing {bed}")
                    self.assertTrue(
                        Path(f"{bed}.tbi").is_file(), f"missing {bed}.tbi",
                    )
                    self.assertFalse(
                        (output_dir / f"{KID}.hap-map-blocks.{parental}.bed").exists(),
                        "uncompressed BED should be removed after bgzip",
                    )

                    expected = [
                        [
                            row["chrom"],
                            str(row["start"]),
                            str(row["end"]),
                            ",".join([
                                row[f"{parental}_haplotype"],
                                str(round(row["haplotype_concordance"], 2)),
                                str(row[TRIO_NUM_HET_SNVS_COLUMN]),
                            ]),
                        ]
                        for row in df.sort(["chrom", "start", "end"]).iter_rows(named=True)
                    ]
                    rows = _read_bgzipped_bed(bed)
                    self.assertEqual(
                        sorted(rows, key=lambda r: (r[0], int(r[1]), int(r[2]))),
                        expected,
                    )

                    # name = "<hap>,<concordance>,<count>"; counts are the
                    # trio's parent-het SNV counts.
                    labels = {"paternal": {"A", "B"}, "maternal": {"C", "D"}}[parental]
                    counts = []
                    for row in rows:
                        hap, concordance, count = row[3].split(",")
                        self.assertIn(hap, labels)
                        self.assertTrue(0.0 <= float(concordance) <= 1.0)
                        counts.append(int(count))
                    self.assertEqual(
                        sorted(counts),
                        sorted(df[TRIO_NUM_HET_SNVS_COLUMN].to_list()),
                    )


if __name__ == "__main__":
    unittest.main()
