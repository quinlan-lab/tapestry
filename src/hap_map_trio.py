import numpy as np
import polars as pl

from util.shell import shell
from util.write_data import write_bed


def extract_bit_vector(l):
    return np.array([int(x) for x in l[0]], dtype=np.uint8)


def write_df_to_vcf(df, vcf, uid):
    df = df.sort(["chrom", "start", "end"])

    header = [
        '##fileformat=VCFv4.2',
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + uid
    ]

    with open(vcf, 'w') as f:
        for line in header:
            f.write(line + '\n')

        for row in df.iter_rows(named=True):
            chrom = row['chrom']
            pos = row['start'] + 1  # VCF is 1-based
            id_ = '.'
            ref = row['REF']
            alt = row['ALT']
            qual = '.'
            flt = '.'
            info = '.'
            fmt_keys = '.'
            fmt_vals = '.'
            f.write(f"{chrom}\t{pos}\t{id_}\t{ref}\t{alt}\t{qual}\t{flt}\t{info}\t{fmt_keys}\t{fmt_vals}\n")


def write_bit_vector_sites_and_mismatches(df_sites_mismatch, uid, parental, output_dir, logger):
    vcf = f"{output_dir}/{uid}.bit-vector-sites-mismatches.{parental}.vcf"
    write_df_to_vcf(df_sites_mismatch, vcf, uid)
    logger.info(f"Wrote {parental} bit-vector-sites-mismatches (for IGV) to: '{vcf}'")

    cmd = (
        f'src/util/compress-index-vcf'
        f' --name {output_dir}/{uid}.bit-vector-sites-mismatches.{parental}'
    )
    shell(cmd)


def write_hap_map_blocks(df_hap_map, uid, parental, output_dir):
    """Write hap-map blocks BED for IGV visualization."""
    df_blocks = df_hap_map.select([
        pl.col("chrom"),
        pl.col("start"),
        pl.col("end"),
        pl.col(f"{parental}_haplotype"),
    ])
    write_bed(output_dir, df_blocks, f"{uid}.hap-map-blocks.{parental}")

    cmd = (
        f'cat {output_dir}/{uid}.hap-map-blocks.{parental}.bed'
        f' | src/util/sort-compress-index-bed'
        f' --name {output_dir}/{uid}.hap-map-blocks.{parental}'
    )
    shell(cmd)
    shell(f'rm {output_dir}/{uid}.hap-map-blocks.{parental}.bed')


def _build_hap_map(df, kid_allele_col, parent_allele_col,
                   phase_block_kid_cols, phase_block_parent_cols,
                   hap_labels, hap_col_name):
    """
    Build hap-map for one parental comparison.

    Groups SNVs by the intersection of kid's and parent's phase blocks,
    compares bit vectors, and assigns haplotype labels.

    Args:
        df: DataFrame from get_all_phasing (df_kid_dad or df_kid_mom)
        kid_allele_col: column name for kid's allele (e.g. "kid_allele_pat")
        parent_allele_col: column name for parent's hap1 allele (e.g. "dad_allele_A")
        phase_block_kid_cols: [start_col, end_col] for kid's phase block
        phase_block_parent_cols: [start_col, end_col] for parent's phase block
        hap_labels: (match_label, mismatch_label) e.g. ("A", "B")
        hap_col_name: output column name e.g. "paternal_haplotype"

    Returns:
        df_hap_map: DataFrame with chrom, start, end, haplotype, concordance, num_het_SNVs
        df_mismatch: DataFrame of mismatch sites with chrom, start, end, REF, ALT
    """
    group_cols = ["chrom"] + phase_block_kid_cols + phase_block_parent_cols

    # Group SNVs by (chrom, kid_phase_block, parent_phase_block),
    # which implicitly finds the intersection of these two types of blocks,
    # and compute the kid's and parent's allele bit vectors in those intersections
    df_grouped = (
        df
        .group_by(group_cols)
        .agg([
            pl.col("start").implode().alias("start_seq"),
            pl.col("end").implode().alias("end_seq"),
            pl.col(kid_allele_col).implode().alias("kid_allele_seq"),
            pl.col(parent_allele_col).implode().alias("parent_allele_seq"),
            pl.col("REF").implode().alias("REF_seq"),
            pl.col("ALT").implode().alias("ALT_seq"),
        ])
        .sort(phase_block_kid_cols[0])  # Sort for reproducibility
    )

    records = []
    data_mismatch = []
    for row in df_grouped.iter_rows(named=True):
        kid_vec = extract_bit_vector(row["kid_allele_seq"])
        parent_vec = extract_bit_vector(row["parent_allele_seq"])

        n = len(kid_vec)

        mismatch = kid_vec != parent_vec
        edit_distance = np.sum(mismatch)
        # Similarity of kid to parent's hap1 vs hap2 sum to 1.0
        # (the parent is het, so hap1 and hap2 are complementary)
        similarity = 1 - (edit_distance / n)

        starts = np.array(row["start_seq"][0])
        ends = np.array(row["end_seq"][0])
        REFs = np.array(row["REF_seq"][0])
        ALTs = np.array(row["ALT_seq"][0])

        if similarity > 0.5:
            haplotype = hap_labels[0]
            concordance = similarity
            mismatch_mask = mismatch
        else:
            haplotype = hap_labels[1]
            concordance = 1.0 - similarity
            mismatch_mask = ~mismatch

        record = {col: row[col] for col in group_cols}
        record[hap_col_name] = haplotype
        record["haplotype_concordance"] = concordance
        record["num_het_SNVs"] = n
        records.append(record)

        data_mismatch.append({
            "chrom": row["chrom"],
            "start": starts[mismatch_mask],
            "end": ends[mismatch_mask],
            "REF": REFs[mismatch_mask],
            "ALT": ALTs[mismatch_mask],
        })

    df_hap_map = (
        pl.DataFrame(records)
        .with_columns(
            pl.max_horizontal(phase_block_kid_cols[0], phase_block_parent_cols[0]).alias("start"),  # start of intersection
            pl.min_horizontal(phase_block_kid_cols[1], phase_block_parent_cols[1]).alias("end"),  # end of intersection
        )
        .select([
            "chrom", "start", "end",
            hap_col_name, "haplotype_concordance", "num_het_SNVs",
        ])
        .sort(["chrom", "start", "end"])
    )

    schema = {"chrom": pl.String, "start": pl.Int64, "end": pl.Int64, "REF": pl.String, "ALT": pl.String}
    df_mismatch = pl.concat([pl.DataFrame(dm, schema=schema) for dm in data_mismatch])

    return df_hap_map, df_mismatch


def get_hap_map(df_kid_dad, df_kid_mom):
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

    Returns:
        df_pat: DataFrame with chrom, start, end, paternal_haplotype,
                haplotype_concordance, num_het_SNVs
        df_mat: DataFrame with chrom, start, end, maternal_haplotype,
                haplotype_concordance, num_het_SNVs
        df_mismatch_pat: mismatch sites for paternal comparison
        df_mismatch_mat: mismatch sites for maternal comparison
    """
    df_pat, df_mismatch_pat = _build_hap_map(
        df_kid_dad, "kid_allele_pat", "dad_allele_A",
        ["start_phase_block_kid", "end_phase_block_kid"],
        ["start_phase_block_dad", "end_phase_block_dad"],
        ("A", "B"), "paternal_haplotype",
    )

    df_mat, df_mismatch_mat = _build_hap_map(
        df_kid_mom, "kid_allele_mat", "mom_allele_C",
        ["start_phase_block_kid", "end_phase_block_kid"],
        ["start_phase_block_mom", "end_phase_block_mom"],
        ("C", "D"), "maternal_haplotype",
    )

    return df_pat, df_mat, df_mismatch_pat, df_mismatch_mat
