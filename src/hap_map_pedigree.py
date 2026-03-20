import numpy as np
import polars as pl

from util.hap_map import extract_bit_vector


def extract_bit_vectors(record):
    bit_vector_hap1 = extract_bit_vector(record["allele_seq_hap1"])
    bit_vector_pat = extract_bit_vector(record["allele_seq_pat"])
    return bit_vector_hap1, bit_vector_pat


def get_hap_map(df_all_phasing):
    df_sites = df_all_phasing.select(["chrom", "start", "end", "REF", "ALT"])

    # Group het SNVs by (chrom, phase_block, iht_block), 
    # which implicitly finds the intersection of these two types of blocks, 
    # and compute the read-backed and inheritance-backed bit vectors in those intersections.
    # Assumes SNVs are phased in kid and parent. 
    df_grouped = (
        df_all_phasing
        .group_by([
            "chrom", 
            "start_phase_block", "end_phase_block", 
            "start_iht_block", "end_iht_block", 
            "founder_label_pat_iht_block", "founder_label_mat_iht_block"
        ])
        .agg([
            pl.col("start").alias("start_seq"),
            pl.col("end").alias("end_seq"),
            pl.col("allele_hap1").alias("allele_seq_hap1"),
            pl.col("allele_pat").alias("allele_seq_pat"),
            pl.col("REF").alias("REF_seq"),
            pl.col("ALT").alias("ALT_seq"),
        ])
        .sort("start_phase_block")  # Sort by phase_block_id for reproducibility 
    )

    records_hap_map = []
    data_mismatch = []
    for record_hap_map in df_grouped.iter_rows(named=True):
        (
            bit_vector_hap1, 
            bit_vector_pat, 
        ) = extract_bit_vectors(record_hap_map)

        record_hap_map["num_het_SNVs"] = len(bit_vector_hap1)

        starts = np.array(record_hap_map["start_seq"])
        ends = np.array(record_hap_map["end_seq"])
        REFs = np.array(record_hap_map["REF_seq"])
        ALTs = np.array(record_hap_map["ALT_seq"])

        hap1_pat_mismatch = bit_vector_hap1 != bit_vector_pat 
        edit_distance_hap1_pat = np.sum(hap1_pat_mismatch)
        similarity_hap1_pat = 1 - (edit_distance_hap1_pat/len(bit_vector_hap1)) 

        if (similarity_hap1_pat > 0.5): # similarity_hap1_pat + similarity_hap1_mat = 1.0 (SNVs are hets)
            record_hap_map["paternal_haplotype"] = record_hap_map["founder_label_pat_iht_block"] + "_hap1"
            record_hap_map["maternal_haplotype"] = record_hap_map["founder_label_mat_iht_block"] + "_hap2"
            record_hap_map["haplotype_concordance"] = similarity_hap1_pat

            starts_mismatch = starts[hap1_pat_mismatch]
            ends_mismatch = ends[hap1_pat_mismatch]
            REFs_mismatch = REFs[hap1_pat_mismatch]
            ALTs_mismatch = ALTs[hap1_pat_mismatch]
        else: 
            record_hap_map["maternal_haplotype"] = record_hap_map["founder_label_mat_iht_block"] + "_hap1"
            record_hap_map["paternal_haplotype"] = record_hap_map["founder_label_pat_iht_block"] + "_hap2"
            record_hap_map["haplotype_concordance"] = 1.0 - similarity_hap1_pat

            starts_mismatch = starts[~hap1_pat_mismatch] 
            ends_mismatch = ends[~hap1_pat_mismatch]
            REFs_mismatch = REFs[~hap1_pat_mismatch]
            ALTs_mismatch = ALTs[~hap1_pat_mismatch]

        records_hap_map.append(record_hap_map)

        data_mismatch.append({
            "chrom": record_hap_map["chrom"],
            "start": starts_mismatch,
            "end": ends_mismatch,
            "REF": REFs_mismatch,
            "ALT": ALTs_mismatch
        })

    df_hap_map = (
        pl
        .DataFrame(records_hap_map)
        .drop([
            "allele_seq_hap1", 
            "allele_seq_pat", 
            "founder_label_pat_iht_block",
            "founder_label_mat_iht_block",
        ])
        .with_columns(
            pl.max_horizontal("start_phase_block", "start_iht_block").alias("start"), # start of intersection
            pl.min_horizontal("end_phase_block", "end_iht_block").alias("end"), # end of intersection
        )
        .select([
            "chrom", "start", "end", 
            "paternal_haplotype", "maternal_haplotype", 
            "haplotype_concordance", "num_het_SNVs",
        ])
        .sort(["chrom", "start", "end"])
    )

    schema = {"chrom": pl.String, "start": pl.Int64, "end": pl.Int64, "REF": pl.String, "ALT": pl.String}
    df_sites_mismatch = pl.concat([pl.DataFrame(dm, schema=schema) for dm in data_mismatch])

    return df_hap_map, df_sites, df_sites_mismatch
