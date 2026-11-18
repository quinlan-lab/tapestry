from itertools import groupby
from operator import itemgetter

import numpy as np
import polars as pl

from hap_map import extract_bit_vector
from report_whatshap_crossovers import UnresolvedSpan


def _assign_resolved_segments(df: pl.DataFrame, spans: list[UnresolvedSpan]) -> pl.DataFrame:
    """Assign each marker to the gap-free segment between excluded spans."""
    if not spans:
        return df.with_columns(
            segment_start=pl.lit(0, dtype=pl.Int64),
            segment_end=pl.lit(None, dtype=pl.Int64),
        )

    segments = []
    for chrom, exclusions in groupby(sorted(spans), key=itemgetter(0)):
        segment_start = 0
        for _, start, end in exclusions:
            segments.append((chrom, segment_start, start))
            segment_start = end
        segments.append((chrom, segment_start, None))
    bounds = pl.DataFrame(
        segments, schema={"chrom": pl.String, "segment_start": pl.Int64, "segment_end": pl.Int64},
        orient="row",
    )

    # Within an excluded span, the preceding segment is the last as-of match,
    # so its end filters the marker out. Unscreened chromosomes have null bounds.
    return (
        df.sort("start").join_asof(
            bounds.sort("segment_start"), left_on="start", right_on="segment_start", by="chrom",
            check_sortedness=False,
        )
        .filter(pl.col("segment_end").is_null() | (pl.col("end") <= pl.col("segment_end")))
        .with_columns(pl.col("segment_start").fill_null(0))
    )


def _build_hap_map(df, kid_allele_col, parent_allele_col,
                   phase_block_kid_id_col, phase_block_parent_id_col,
                   phase_block_kid_cols, phase_block_parent_cols,
                   hap_labels, hap_col_name, unresolved_spans):
    """
    Build hap-map for one parental comparison.

    Groups SNVs by the intersection of kid's and parent's phase blocks,
    compares bit vectors, and assigns haplotype labels.

    Blocks are keyed by phase-set identity (the PS-derived phase_block_id),
    not by their genomic start/end span: distinct phase sets can share a span
    after nesting, so grouping on (kid_phase_set, parent_phase_set) is the only
    key that maps each SNV to exactly one hap-map block. The block's start/end
    are carried along (identical within a phase-set group) to compute the
    intersection interval for output.

    Args:
        df: DataFrame from get_all_phasing (df_kid_dad or df_kid_mom)
        kid_allele_col: column name for kid's allele (e.g. "kid_allele_pat")
        parent_allele_col: column name for parent's hap1 allele (e.g. "dad_allele_A")
        phase_block_kid_id_col: phase-set id column for the kid (e.g. "phase_block_id_kid")
        phase_block_parent_id_col: phase-set id column for the parent (e.g. "phase_block_id_dad")
        phase_block_kid_cols: [start_col, end_col] for kid's phase block
        phase_block_parent_cols: [start_col, end_col] for parent's phase block
        hap_labels: (match_label, mismatch_label) e.g. ("A", "B")
        hap_col_name: output column name e.g. "paternal_haplotype"
        unresolved_spans: Sorted or unsorted disjoint BED intervals excluded for this parent.

    Returns:
        df_hap_map: DataFrame with chrom, start, end, haplotype, concordance, num_het_SNVs
        df_mismatch: DataFrame of mismatch sites with chrom, start, end, REF, ALT
    """
    group_cols = ["chrom", phase_block_kid_id_col, phase_block_parent_id_col,
                  "segment_start", "segment_end"]
    span_cols = phase_block_kid_cols + phase_block_parent_cols

    # Compare the two bit vectors within each phase-set pair and segment.
    # The phase-block bounds are constant within each group.
    # Assumes SNVs are phased in kid and parent.
    df_grouped = (
        _assign_resolved_segments(df, unresolved_spans)
        .group_by(group_cols)
        .agg([
            *[pl.col(c).first().alias(c) for c in span_cols],
            pl.col("start").alias("start_seq"),
            pl.col("end").alias("end_seq"),
            pl.col(kid_allele_col).alias("kid_allele_seq"),
            pl.col(parent_allele_col).alias("parent_allele_seq"),
            pl.col("REF").alias("REF_seq"),
            pl.col("ALT").alias("ALT_seq"),
        ])
        .sort(group_cols)  # Stable output if distinct PS pairs share a span.
    )

    records = []
    mismatches = []
    for row in df_grouped.iter_rows(named=True):
        kid_vec = extract_bit_vector(row["kid_allele_seq"])
        parent_vec = extract_bit_vector(row["parent_allele_seq"])

        n = len(kid_vec)

        mismatch = kid_vec != parent_vec
        edit_distance = np.sum(mismatch)
        similarity = 1 - (edit_distance / n)

        # Similarity of kid to parent's hap1 vs hap2 sum to 1.0
        # (the parent is het, so hap1 and hap2 are complementary)
        if similarity > 0.5:
            haplotype = hap_labels[0]
            concordance = similarity
            mismatch_mask = mismatch
        else:
            haplotype = hap_labels[1]
            concordance = 1.0 - similarity
            mismatch_mask = ~mismatch

        start = max(row[phase_block_kid_cols[0]], row[phase_block_parent_cols[0]], row["segment_start"])
        end = min(row[phase_block_kid_cols[1]], row[phase_block_parent_cols[1]])
        if row["segment_end"] is not None:
            end = min(end, row["segment_end"])
        records.append((row["chrom"], start, end, haplotype, concordance, n))

        for i in np.flatnonzero(mismatch_mask):
            mismatches.append((row["chrom"], row["start_seq"][i], row["end_seq"][i],
                               row["REF_seq"][i], row["ALT_seq"][i]))

    coords = {"chrom": pl.String, "start": pl.Int64, "end": pl.Int64}
    df_hap_map = pl.DataFrame(records, schema={
        **coords, hap_col_name: pl.String,
        "haplotype_concordance": pl.Float64, "num_het_SNVs_in_parent": pl.Int64,
    }, orient="row").sort(list(coords), maintain_order=True)
    df_mismatch = pl.DataFrame(mismatches, schema={
        **coords, "REF": pl.String, "ALT": pl.String,
    }, orient="row").sort([*coords, "REF", "ALT"])

    return df_hap_map, df_mismatch


def get_hap_map(
    df_kid_dad: pl.DataFrame, df_kid_mom: pl.DataFrame, *,
    unresolved_spans: dict[str, list[UnresolvedSpan]] | None = None,
):
    """
    Build the haplotype map by comparing bit vectors.

    Labeling convention:
        A = dad's hap1, B = dad's hap2 (fixed)
        C = mom's hap1, D = mom's hap2 (fixed)

    For each block (intersection of kid's and parent's phase blocks),
    compare the kid's allele sequence to the parent's hap1 allele sequence:
        - Paternal: kid_allele_pat vs dad_allele_A
          similarity > 0.5 → paternal_haplotype = "A", else "B"
        - Maternal: kid_allele_mat vs mom_allele_C
          similarity > 0.5 → maternal_haplotype = "C", else "D"

    With unresolved spans, compare each resolved portion of a phase-set pair
    independently. Excluded spans are keyed by "father" and "mother".

    Returns:
        df_pat: DataFrame with chrom, start, end, paternal_haplotype,
                haplotype_concordance, num_het_SNVs
        df_mat: DataFrame with chrom, start, end, maternal_haplotype,
                haplotype_concordance, num_het_SNVs
        df_mismatch_pat: mismatch sites for paternal comparison
        df_mismatch_mat: mismatch sites for maternal comparison
    """
    unresolved_spans = unresolved_spans or {}
    df_pat, df_mismatch_pat = _build_hap_map(
        df_kid_dad, "kid_allele_pat", "dad_allele_A",
        "phase_block_id_kid", "phase_block_id_dad",
        ["start_phase_block_kid", "end_phase_block_kid"],
        ["start_phase_block_dad", "end_phase_block_dad"],
        ("A", "B"), "paternal_haplotype", unresolved_spans.get("father", []),
    )

    df_mat, df_mismatch_mat = _build_hap_map(
        df_kid_mom, "kid_allele_mat", "mom_allele_C",
        "phase_block_id_kid", "phase_block_id_mom",
        ["start_phase_block_kid", "end_phase_block_kid"],
        ["start_phase_block_mom", "end_phase_block_mom"],
        ("C", "D"), "maternal_haplotype", unresolved_spans.get("mother", []),
    )

    return df_pat, df_mat, df_mismatch_pat, df_mismatch_mat
