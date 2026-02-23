"""
Generate mock input data for testing src/trio_phase_meth_to_parent_haps.py
and src/expand_to_all_cpgs.py (trio mode).

Run from the repo root:
    python tests/create_mock_data.py
"""
import gzip
import os
from pathlib import Path

OUTPUT_DIR = Path(__file__).parent / "mock_data"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# trio.ped
# ---------------------------------------------------------------------------
def write_ped():
    path = OUTPUT_DIR / "trio.ped"
    lines = [
        "TRIO\tkid\tdad\tmom\t1",
        # "TRIO\tdad\t0\t0\t1",
        # "TRIO\tmom\t0\t0\t2",
    ]
    path.write_text("\n".join(lines) + "\n")
    print(f"Wrote {path}")


# ---------------------------------------------------------------------------
# whatshap_blocks.tsv  (pre-computed output from `whatshap stats --block-list`)
# ---------------------------------------------------------------------------
def write_whatshap_blocks():
    path = OUTPUT_DIR / "whatshap_blocks.tsv"
    lines = [
        "#sample\tchromosome\tphase_set\tfrom\tto\tvariants",
        "kid\tchr1\t1001\t1001\t100000\t12",
        "kid\tchr1\t120001\t120001\t190000\t8",
    ]
    path.write_text("\n".join(lines) + "\n")
    print(f"Wrote {path}")


# ---------------------------------------------------------------------------
# Methylation BED files  (gzipped, with tapestry-style header)
# CpG sites on chr1 in 4 regions:
#   block 1: 1000–100000   (0-based start of block)
#   block 2: 120000–190000
#   outside blocks: a few sites
# ---------------------------------------------------------------------------
CHR1_CPG_SITES = (
    # (start, end)  — 0-based half-open
    # Inside block 1
    [(1000 + i * 1800, 1001 + i * 1800) for i in range(20)]
    # Inside block 2
    + [(120000 + i * 2000, 120001 + i * 2000) for i in range(15)]
    # Outside both blocks
    + [(200000 + i * 5000, 200001 + i * 5000) for i in range(5)]
)


def meth_level(i: int, seed: float) -> float:
    return round((seed + i * 0.03) % 1.0, 4)


def write_meth_bed(filename: str, mode: str, hap: str):
    """Write a synthetic pb-CpG-tools methylation BED file (gzipped)."""
    path = OUTPUT_DIR / filename
    seed = {"hap1": 0.2, "hap2": 0.6}[hap]
    lines = [
        f"##pileup-mode={mode}",
        "#chrom\tbegin\tend\tmod_score\ttype\tcov",
    ]
    for i, (start, end) in enumerate(CHR1_CPG_SITES):
        ml_pct = round(meth_level(i, seed) * 100, 2)
        cov = 15 + (i % 10)
        lines.append(f"chr1\t{start}\t{end}\t{ml_pct}\tCG\t{cov}")
    with gzip.open(path, "wt") as f:
        f.write("\n".join(lines) + "\n")
    print(f"Wrote {path}")


def write_meth_combined_bed(filename: str, mode: str):
    """Write a synthetic combined (unphased) pb-CpG-tools BED (gzipped)."""
    path = OUTPUT_DIR / filename
    seed = 0.4
    lines = [
        f"##pileup-mode={mode}",
        "#chrom\tbegin\tend\tmod_score\ttype\tcov",
    ]
    for i, (start, end) in enumerate(CHR1_CPG_SITES):
        ml_pct = round(meth_level(i, seed) * 100, 2)
        cov = 20 + (i % 8)
        lines.append(f"chr1\t{start}\t{end}\t{ml_pct}\tCG\t{cov}")
    with gzip.open(path, "wt") as f:
        f.write("\n".join(lines) + "\n")
    print(f"Wrote {path}")


# ---------------------------------------------------------------------------
# all_cpg_sites_in_reference.bed  (plain TSV, no header)
# Superset of CHR1_CPG_SITES plus a few extra reference-only CpGs
# ---------------------------------------------------------------------------
def write_all_cpg_sites_in_reference():
    path = OUTPUT_DIR / "all_cpg_sites_in_reference.bed"
    extra = [(50000 + i * 3000, 50001 + i * 3000) for i in range(10)]
    all_sites = sorted(CHR1_CPG_SITES + extra)
    lines = [f"chr1\t{s}\t{e}" for s, e in all_sites]
    path.write_text("\n".join(lines) + "\n")
    print(f"Wrote {path}")


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    write_ped()
    write_whatshap_blocks()
    for mode in ["count", "model"]:
        for hap in ["hap1", "hap2"]:
            write_meth_bed(f"kid.GRCh38.haplotagged.{hap}.{mode}.bed.gz", mode=mode, hap=hap)
        write_meth_combined_bed(f"kid.GRCh38.haplotagged.combined.{mode}.bed.gz", mode=mode)
    write_all_cpg_sites_in_reference()
    print("Done creating mock data.")
