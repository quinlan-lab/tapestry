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

Run locally (no CHPC, no real data needed):
    PYTHONPATH=src:src/util .venv/bin/python tests/test_recombination_dedup.py
"""

import sys
import tempfile
from pathlib import Path

import polars as pl

from phasing_trio import get_phase_blocks, get_all_phasing
from hap_map_trio import get_hap_map

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


def _check(df, parental, n_input):
    """Assert no fan-out: every input het SNV counted in exactly one block."""
    failures = []

    dups = (
        df.group_by(["chrom", "start", "end"]).len().filter(pl.col("len") > 1)
    )
    if len(dups):
        failures.append(f"{parental}: duplicate hap-map blocks:\n{dups}")

    total = df["num_het_SNVs_in_parent"].sum()
    if total != n_input:
        failures.append(
            f"{parental}: num_het_SNVs sums to {total} but there were "
            f"{n_input} input phased het SNVs (fan-out double-counting)."
        )
    return failures


def main():
    with tempfile.TemporaryDirectory() as tmp:
        d = Path(tmp)
        vcf_path, blocks_path = _write_inputs(d)

        df_blocks_kid = get_phase_blocks(blocks_path, KID)
        df_blocks_dad = get_phase_blocks(blocks_path, DAD)
        df_blocks_mom = get_phase_blocks(blocks_path, MOM)

        df_kid_dad, df_kid_mom = get_all_phasing(
            vcf_path, KID, DAD, MOM,
            df_blocks_kid, df_blocks_dad, df_blocks_mom,
        )

        df_pat, df_mat, _, _ = get_hap_map(df_kid_dad, df_kid_mom)

    failures = []
    failures += _check(df_pat, "paternal", len(df_kid_dad))
    failures += _check(df_mat, "maternal", len(df_kid_mom))

    # Concrete expectation: the nested set yields its OWN single block of 3 SNVs,
    # and the big block holds the remaining 4 -- not 7-in-the-big + a triplicate.
    pat_sizes = sorted(df_pat["num_het_SNVs_in_parent"].to_list())
    if pat_sizes != [3, 4]:
        failures.append(
            f"paternal block sizes {pat_sizes}, expected [3, 4] "
            "(4 big-block + 3 nested-block SNVs, no double counting)."
        )

    print("paternal hap-map blocks:")
    print(df_pat)
    print("maternal hap-map blocks:")
    print(df_mat)

    if failures:
        print("\nFAIL:")
        for f in failures:
            print(" -", f)
        sys.exit(1)
    print("\nPASS: each SNV assigned to exactly one phase set; no fan-out, no duplicate blocks.")


if __name__ == "__main__":
    main()
