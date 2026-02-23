"""
Verification test for the trio workflow.

Run from the repo root on a cluster with all dependencies available:
    PYTHONPATH=src:src/util python tests/test_trio_workflow.py

Requires: whatshap, cyvcf2, bgzip, tabix, bedGraphToBigWig in PATH,
          and all Python packages from requirements.txt.

Generate mock input data first (one-time):
    python tests/create_mock_data.py
"""

import subprocess
import sys
import os
from pathlib import Path

REPO   = Path(__file__).parent.parent
MOCK   = REPO / "tests" / "mock_data"
OUT    = REPO / "tests" / "test_output"
PYTHON = sys.executable

sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "src" / "util"))


def header(msg):
    print(f"\n{'='*60}")
    print(f"  {msg}")
    print('='*60)


def ok(msg):   print(f"  [PASS] {msg}")
def fail(msg): print(f"  [FAIL] {msg}"); sys.exit(1)


def run(cmd, **kwargs):
    env = {**os.environ, "PYTHONPATH": f"{REPO/'src'}:{REPO/'src/util'}"}
    result = subprocess.run(cmd, capture_output=True, text=True, env=env, **kwargs)
    if result.returncode != 0:
        print(result.stdout[-2000:])
        print(result.stderr[-2000:])
        fail(f"Command failed: {' '.join(str(c) for c in cmd)}")
    return result


UID = "kid"


# ---------------------------------------------------------------------------
# Ensure mock data exists
# ---------------------------------------------------------------------------
header("Generating mock data")
run([PYTHON, str(REPO / "tests" / "create_mock_data.py")])
ok("Mock data ready")


# ---------------------------------------------------------------------------
# Step 1: trio_phase_meth_to_parent_haps.py
# ---------------------------------------------------------------------------
header("Step 1: whatshap unphase + whatshap phase + trio_phase_meth_to_parent_haps.py")

OUT.mkdir(parents=True, exist_ok=True)

# Create a minimal joint-called VCF (no variants) as input to run-whatshap-phase.sh.
# In production this would be the real joint-called trio VCF.
mock_joint_vcf_raw = OUT / "trio_joint_called_mock.vcf"
mock_joint_vcf = OUT / "trio_joint_called_mock.vcf.gz"
mock_joint_vcf_raw.write_text(
    "##fileformat=VCFv4.2\n"
    f"##contig=<ID=chr1,length=248956422>\n"
    f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{UID}\tdad\tmom\n"
)
run(["bgzip", "-f", str(mock_joint_vcf_raw)])
run(["tabix", "-p", "vcf", str(mock_joint_vcf)])

# Step 1a: unphase (mirrors run-whatshap-phase.sh step 1a)
mock_vcf_unphased_raw = OUT / "trio_unphased_mock.vcf"
result = run(["whatshap", "unphase", str(mock_joint_vcf)])
mock_vcf_unphased_raw.write_text(result.stdout)
ok("whatshap unphase completed")

# Step 1b: phase (mirrors run-whatshap-phase.sh step 1b)
# With no variants, whatshap phase produces a valid but empty phased VCF.
mock_vcf_raw = OUT / "trio_phased_mock.vcf"
mock_vcf = OUT / "trio_phased_mock.vcf.gz"
run([
    "whatshap", "phase",
    "--ped",          str(MOCK / "trio.ped"),
    "--sample",       UID,
    "--use-ped-samples",
    "--output",       str(mock_vcf_raw),
    str(mock_vcf_unphased_raw),
])
run(["bgzip", "-f", str(mock_vcf_raw)])
run(["tabix", "-p", "vcf", str(mock_vcf)])

run([
    PYTHON, str(REPO / "src" / "trio_phase_meth_to_parent_haps.py"),
    "--uid",               UID,
    "--vcf_trio_phased",   str(mock_vcf),
    "--ped",               str(MOCK / "trio.ped"),
    "--bed_meth_count_hap1", str(MOCK / f"{UID}.GRCh38.haplotagged.hap1.count.bed.gz"),
    "--bed_meth_count_hap2", str(MOCK / f"{UID}.GRCh38.haplotagged.hap2.count.bed.gz"),
    "--bed_meth_model_hap1", str(MOCK / f"{UID}.GRCh38.haplotagged.hap1.model.bed.gz"),
    "--bed_meth_model_hap2", str(MOCK / f"{UID}.GRCh38.haplotagged.hap2.model.bed.gz"),
    "--output_dir",        str(OUT),
])
ok("trio_phase_meth_to_parent_haps.py completed")

# Verify output BED
import polars as pl
from util.read_data import read_bed_and_header

bed_path = OUT / f"{UID}.dna-methylation.founder-phased.bed"
assert bed_path.exists(), f"Missing output BED: {bed_path}"

df = read_bed_and_header(str(bed_path))
expected_cols = {
    'chrom', 'start', 'end',
    'start_hap_map_block', 'end_hap_map_block',
    'haplotype_concordance_in_hap_map_block', 'num_het_SNVs_in_hap_map_block',
    'total_read_count_pat', 'total_read_count_mat',
    'founder_haplotype_pat', 'founder_haplotype_mat',
    'methylation_level_pat_count', 'methylation_level_mat_count',
    'methylation_level_pat_model', 'methylation_level_mat_model',
}
missing = expected_cols - set(df.columns)
assert not missing, f"Output BED missing columns: {missing}"
assert len(df.columns) == 15, f"Expected 15 columns, got {len(df.columns)}: {df.columns}"
ok(f"Output BED has correct 15 columns")

assert df['haplotype_concordance_in_hap_map_block'].is_null().all(), \
    "haplotype_concordance_in_hap_map_block should be null for trio"
ok("haplotype_concordance_in_hap_map_block is null for all rows")

assert (df['founder_haplotype_pat'] == 'dad').all(), "founder_haplotype_pat should be 'dad'"
assert (df['founder_haplotype_mat'] == 'mom').all(), "founder_haplotype_mat should be 'mom'"
ok("founder_haplotype_pat='dad', founder_haplotype_mat='mom' for all rows")

for parental in ['pat', 'mat']:
    for mode in ['count', 'model']:
        bw = OUT / f"{UID}.dna-methylation.founder-phased.{parental}.{mode}.hg38.bw"
        assert bw.exists(), f"Missing BigWig: {bw}"
ok("BigWig files produced for pat/mat x count/model")


# ---------------------------------------------------------------------------
# Step 2: expand_to_all_cpgs.py (no --bed_het_site_mismatches)
# ---------------------------------------------------------------------------
header("Step 2: expand_to_all_cpgs.py (trio mode — no --bed_het_site_mismatches)")

bed_all_cpgs_out = OUT / f"{UID}.dna-methylation.founder-phased.all_cpgs.bed"

run([
    PYTHON, str(REPO / "src" / "expand_to_all_cpgs.py"),
    "--bed_all_cpgs_in_reference",       str(MOCK / "all_cpg_sites_in_reference.bed"),
    "--bed_meth_founder_phased",         str(bed_path),
    "--bed_meth_count_unphased",         str(MOCK / f"{UID}.GRCh38.haplotagged.combined.count.bed.gz"),
    "--bed_meth_model_unphased",         str(MOCK / f"{UID}.GRCh38.haplotagged.combined.model.bed.gz"),
    "--bed_meth_founder_phased_all_cpgs", str(bed_all_cpgs_out),
    "--uid",            UID,
    "--vcf_joint_called", str(mock_vcf),
])
ok("expand_to_all_cpgs.py completed")

# Verify output
from util.read_data import read_dataframe_from_bed

sorted_bed = bed_all_cpgs_out.with_name(
    bed_all_cpgs_out.stem.replace(".bed", "") + ".sorted.bed.gz"
)
# expand_to_all_cpgs.py writes <stem>.sorted.bed.gz
sorted_gz = OUT / f"{UID}.dna-methylation.founder-phased.all_cpgs.sorted.bed.gz"
assert sorted_gz.exists(), f"Missing output: {sorted_gz}"

df2 = read_dataframe_from_bed(str(sorted_gz))
expected_22 = {
    'chrom', 'start_cpg', 'end_cpg',
    'total_read_count', 'methylation_level_count', 'methylation_level_model',
    'start_hap_map_block', 'end_hap_map_block',
    'haplotype_concordance_in_hap_map_block', 'num_het_SNVs_in_hap_map_block',
    'total_read_count_pat', 'total_read_count_mat',
    'founder_haplotype_pat', 'founder_haplotype_mat',
    'methylation_level_pat_count', 'methylation_level_mat_count',
    'methylation_level_pat_model', 'methylation_level_mat_model',
    'cpg_is_within_50bp_of_mismatch_site',
    'cpg_overlaps_at_least_one_snv', 'snv_genotypes', 'cpg_is_allele_specific',
}
missing = expected_22 - set(df2.columns)
assert not missing, f"Output missing columns: {missing}"
assert len(df2.columns) == 22, f"Expected 22 columns, got {len(df2.columns)}"
ok("Output has correct 22 columns")

assert (df2['cpg_is_within_50bp_of_mismatch_site'] == False).all(), \
    "cpg_is_within_50bp_of_mismatch_site should be False for all rows (trio mode)"
ok("cpg_is_within_50bp_of_mismatch_site is False for all rows")


header("ALL CHECKS PASSED")