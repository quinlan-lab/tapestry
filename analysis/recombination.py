import polars as pl
import bioframe as bf
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys
import importlib

REPO_DIR = Path('/scratch/ucgd/lustre-labs/quinlan/u6018199/tapestry')
sys.path.append(str(REPO_DIR / 'src')) 
sys.path.append(str(REPO_DIR / 'src/util')) 

from get_meth_hap1_hap2 import read_meth_level
from get_all_phasing import get_iht_blocks

import version_sort
importlib.reload(version_sort)
from version_sort import version_sort

PB_CPG_TOOL_MODE = 'model'

IHT_PHASED_DIR = Path('/scratch/ucgd/lustre-labs/quinlan/data-shared/haplotype-maps/CEPH1463.GRCh38')
METH_READ_PHASED_DIR = Path(f'/scratch/ucgd/lustre-labs/quinlan/data-shared/dna-methylation/CEPH1463.GRCh38.hifi.{PB_CPG_TOOL_MODE}.read-backed-phased')

TXT_IHT_BLOCKS = f"{IHT_PHASED_DIR}/CEPH1463.GRCh38.iht.sorted.txt" # multi-sample iht blocks file from gtg-ped-map/gtg-concordance

REFERENCE_GENOME = "hg38"
RNG_SEED = None

SAMPLES = ['200081', 'NA12877'] # TODO: extend to all of CEPH1463

def read_meth_level_wrapper(sample, meth_read_phased_dir, pb_cpg_tool_mode): 
    bed_meth_combined = f"{meth_read_phased_dir}/{sample}.GRCh38.haplotagged.combined.bed.gz" # bed file from aligned_bam_to_cpg_scores for both haplotypes combined    

    return read_meth_level(
        bed_meth_combined, # type: ignore
        pb_cpg_tool_mode
    ).rename({ "chromosome": "chrom" })

def get_recombination_blocks(sample, haplotype, txt_iht_blocks):
    df_iht = get_iht_blocks(sample, txt_iht_blocks)

    recomb_blocks = []
    for _, df_iht_chrom in df_iht.group_by("chrom", maintain_order=True):
        df_iht_chrom = df_iht_chrom.sort("start")

        # Find where founder_label changes
        df_recomb_chrom = (
            df_iht_chrom
            .with_columns(
                pl.col("start").shift(1).alias("start_prev"),
                pl.col("end").shift(1).alias("end_prev"),
                pl.col(f"founder_label_{haplotype}").shift(1).alias(f"founder_label_{haplotype}_prev"),
            )
            .with_columns(
                (pl.col(f"founder_label_{haplotype}") != pl.col(f"founder_label_{haplotype}_prev")).alias(f"founder_label_{haplotype}_change")
            )
            .filter(pl.col(f"founder_label_{haplotype}_change").is_not_null())
            .filter(pl.col(f"founder_label_{haplotype}_change"))
            .with_columns(
                pl.col("end_prev").alias(f"{haplotype}_recomb_block_start"), 
                pl.col("start").alias(f"{haplotype}_recomb_block_end")
            )
            .select(["chrom", f"{haplotype}_recomb_block_start", f"{haplotype}_recomb_block_end"])
            .rename({
                f"{haplotype}_recomb_block_start": "start",
                f"{haplotype}_recomb_block_end": "end"
            })
            .with_columns(
                (pl.col("end") - pl.col("start")).alias("length"),
                pl.lit(haplotype).alias("haplotype")
            )
        )
        recomb_blocks.extend(list(df_recomb_chrom.iter_rows(named=True)))

    return version_sort(pl.DataFrame(recomb_blocks))

def get_recombination_blocks_both_haplotypes(sample, txt_iht_blocks):
    df_recomb_pat = get_recombination_blocks(sample, "pat", txt_iht_blocks)
    df_recomb_mat = get_recombination_blocks(sample, "mat", txt_iht_blocks)
    return version_sort(pl.concat([df_recomb_pat, df_recomb_mat]))

def sample_intervals_from_complement(
    df: pl.DataFrame,
    reference_genome: str,
    rng_seed: int
) -> pl.DataFrame:
    """
    For each interval in df, sample a random interval of the same length
    from the complement of df in the genome, using polars DataFrames where possible.

    Args:
        df: pl.DataFrame with columns ["chrom", "start", "end"]
        reference_genome: str, genome name for bioframe.fetch_chromsizes
        rng_seed: int, random seed for reproducibility

    Returns:
        pl.DataFrame with columns ["chrom", "start", "end"]
    """
    # Compute complement intervals 
    complement_intervals = bf.complement(
        df=df.select(["chrom", "start", "end"]).to_pandas(),
        view_df=bf.fetch_chromsizes(db=reference_genome),
    )
    complement_intervals["length"] = complement_intervals["end"] - complement_intervals["start"]
    complement_intervals = pl.from_pandas(complement_intervals)
    complement_intervals = complement_intervals.drop("view_region")
    
    # For reproducibility
    rng = np.random.default_rng(rng_seed)

    # Build a flat list of all complement intervals for sampling
    all_comps = list(complement_intervals.iter_rows(named=True))

    # Get original lengths
    original_lengths = df["end"] - df["start"]

    # Sample intervals
    sampled_intervals = []
    for length in original_lengths:
        # Find all complement intervals that are at least as long as "length"
        candidates = [iv for iv in all_comps if iv["length"] >= length]
        if not candidates:
            continue  # skip if no suitable region
        chosen_iv = candidates[rng.integers(len(candidates))]
        chrom, start, end = chosen_iv["chrom"], chosen_iv["start"], chosen_iv["end"]
        max_start = end - length
        if max_start == start:
            chosen_start = start
        else:
            chosen_start = rng.integers(start, max_start + 1)
        sampled_intervals.append({
            "chrom": chrom,
            "start": chosen_start,
            "end": chosen_start + length,
        })

    return version_sort(
        pl
        .DataFrame(sampled_intervals)
        .with_columns(
            (pl.col("end") - pl.col("start")).alias("length")
        )
    )

def compute_methylation(df, df_meth): 
    df_intersected = bf.overlap(
        df.to_pandas(),
        df_meth.to_pandas(),
        how='inner',
        suffixes=('_df', '_meth')
    )
    df_intersected = pl.from_pandas(df_intersected)
    return version_sort(
        df_intersected
        .group_by(["chrom_df", "start_df", "end_df", "length_df"])
        .agg([
            pl.col("total_read_count_meth").mean().alias("mean_total_read_count"),
            pl.col("methylation_level_meth").mean().alias("mean_methylation_level"),
            pl.len().alias("num_cpg_sites"),
        ])
        .rename({
            "chrom_df": "chrom",
            "start_df": "start",
            "end_df": "end",
            "length_df": "length"
        })
        .with_columns(
            (pl.col("num_cpg_sites") / pl.col("length")).alias("num_cpg_sites_per_bp")
        )
    )

def downsample_to_equal_size(df1, df2, seed=42):
    """
    Downsample the larger of two dataframes to match the size of the smaller one.
    Returns the two dataframes, both with the same number of rows.
    """
    n1 = len(df1)
    n2 = len(df2)
    if n1 < n2:
        df2 = df2.sample(n=n1, shuffle=True, seed=seed)
    else:
        df1 = df1.sample(n=n2, shuffle=True, seed=seed)
    return df1, df2

def plot_histograms(
    df_recomb, 
    df_control, 
    column="mean_methylation_level", 
    bins=np.arange(0, 1.05, 0.01), 
    xlabel="Methylation Level", 
    ylabel="Counts", 
    legend_labels=("Recombination Blocks", "Control Intervals"),
    figsize=(8, 5),
    df_meth=None
):
    """
    Plot histograms for a specified column from two dataframes.

    Args:
        df_recomb: DataFrame for recombination blocks.
        df_control: DataFrame for control intervals.
        column: The column to plot histograms for.
        bins: Number of bins for the histogram.
        xlabel: Label for the x-axis. 
        ylabel: Label for the y-axis.
        legend_labels: Tuple of labels for the legend.
        figsize: Figure size.
    """
    df_recomb, df_control = downsample_to_equal_size(df_recomb, df_control)

    plt.figure(figsize=figsize)
    # Superimpose the distribution of methylation from the base methylation dataframe, if plotting methylation
    if column == "mean_methylation_level":
        # Try to use df_meth if available
        try:
            df_meth = df_meth.sample(n=len(df_recomb), shuffle=True) # type: ignore
            plt.hist(
                df_meth["methylation_level"], 
                bins=bins, # type: ignore
                alpha=0.4, 
                label="Random CpG Sites",
                color="gray"
            )
        except Exception:
            pass
    plt.hist(
        df_recomb[column], 
        bins=bins, # type: ignore
        alpha=0.6, 
        label=legend_labels[0]
    )
    plt.hist(
        df_control[column], 
        bins=bins, # type: ignore
        alpha=0.6, 
        label=legend_labels[1]
    )
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend(frameon=False)


def plot_cdf(
    df_recomb, 
    df_control, 
    column, 
    xlabel=None, 
    ylabel="Cumulative Fraction", 
    legend_labels=("Recombination Blocks", "Control Intervals"),
    figsize=(8, 5),
    df_meth=None
):
    """
    Plot CDFs for a specified column from two dataframes.

    Args:
        df_recomb: DataFrame for recombination blocks.
        df_control: DataFrame for control intervals.
        column: The column to plot CDFs for.
        xlabel: Label for the x-axis. If None, uses the column name.
        ylabel: Label for the y-axis.
        legend_labels: Tuple of labels for the legend.
        figsize: Figure size.
    """
    df_recomb, df_control = downsample_to_equal_size(df_recomb, df_control)

    plt.figure(figsize=figsize)
    # Plot CDF for random CpG sites if methylation
    if column == "mean_methylation_level":
        try:
            df_meth = df_meth.sample(n=len(df_recomb), shuffle=True) # type: ignore
            values = np.sort(df_meth["methylation_level"])
            cdf = np.arange(1, len(values)+1) / len(values)
            plt.plot(values, cdf, label="Random CpG Sites", color="gray", alpha=0.7)
        except Exception:
            pass

    for df, label in zip([df_recomb, df_control], legend_labels):
        values = np.sort(df[column])
        cdf = np.arange(1, len(values)+1) / len(values)
        plt.plot(values, cdf, label=label, alpha=0.8)

    plt.xlabel(xlabel if xlabel is not None else column)
    plt.ylabel(ylabel)
    plt.legend(frameon=False)
    plt.tight_layout()

def run_statistical_tests(df_recomb, df_control):
    """
    Perform Kolmogorov-Smirnov (K-S) test and Mann-Whitney U test
    on mean methylation levels between recombination blocks and control intervals.
    Prints the results of both tests.
    """
    from scipy.stats import ks_2samp, mannwhitneyu

    data1 = df_recomb["mean_methylation_level"] # type: ignore
    data2 = df_control["mean_methylation_level"] # type: ignore

    # Kolmogorov-Smirnov test
    ks_stat, ks_pvalue = ks_2samp(data1, data2) 
    print('Kolmogorov-Smirnov test:')
    print(f"Maximum absolute difference between CDFs: {ks_stat:.4f}")
    print(f"p-value: {ks_pvalue:.4e}")
    print('')

    # Mann-Whitney U test
    mw_stat, mw_pvalue = mannwhitneyu(data1, data2) 
    print('Mann-Whitney U test:')
    print(f"U statistic: {mw_stat:.4f}")
    print(f"p-value: {mw_pvalue:.4e}")

def compute_methylation_in_blocks(sample, txt_iht_blocks, reference_genome, rng_seed, meth_read_phased_dir, pb_cpg_tool_mode):
    df_recomb = get_recombination_blocks_both_haplotypes(sample, txt_iht_blocks)
    df_control = sample_intervals_from_complement(df_recomb, reference_genome, rng_seed)
    df_meth = read_meth_level_wrapper(sample, meth_read_phased_dir, pb_cpg_tool_mode)
    df_recomb_meth = compute_methylation(df_recomb, df_meth)
    df_control_meth = compute_methylation(df_control, df_meth)
    return df_recomb_meth, df_control_meth

def main(samples, txt_iht_blocks, reference_genome, rng_seed, meth_read_phased_dir, pb_cpg_tool_mode):
    df_recomb_meth_all = []
    df_control_meth_all = []
    for i, sample in enumerate(samples):
        if i == 2: # TODO: remove
            break
        print(f"Computing methylation in blocks for sample {i+1} of {len(samples)}...")
        df_recomb_meth, df_control_meth = compute_methylation_in_blocks(
            sample, txt_iht_blocks, reference_genome, rng_seed, meth_read_phased_dir, pb_cpg_tool_mode
        )
        df_recomb_meth = df_recomb_meth.with_columns(pl.lit(sample).alias("sample"))
        df_control_meth = df_control_meth.with_columns(pl.lit(sample).alias("sample"))
        print(f"...done")
        df_recomb_meth_all.append(df_recomb_meth)
        df_control_meth_all.append(df_control_meth)

    df_recomb_meth_all = pl.concat(df_recomb_meth_all)
    df_control_meth_all = pl.concat(df_control_meth_all)
    
    script_dir = Path(__file__).parent
    df_recomb_meth_all.write_csv(script_dir / "recomb_blocks_methylation.csv")
    df_control_meth_all.write_csv(script_dir / "control_intervals_methylation.csv")

    print("Wrote methylation data to csv files")
    print("Done")

if __name__ == "__main__":
    main(SAMPLES, TXT_IHT_BLOCKS, REFERENCE_GENOME, RNG_SEED, METH_READ_PHASED_DIR, PB_CPG_TOOL_MODE)