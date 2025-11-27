import polars as pl 
import bioframe as bf
import matplotlib.pyplot as plt 
from matplotlib_venn import venn2

from version_sort import version_sort
from add_locus import add_locus

def find_unique_and_common_records(df1, df2, min_overlap):
    # track records using an index:  
    df1 = df1.with_row_index("idx") # do not use "index"
    df2 = df2.with_row_index("idx") # do not use "index"

    # find overlaps
    matches = bf.overlap(
        df1.to_pandas(), 
        df2.to_pandas(), 
        how='inner', 
        return_overlap=True,
        suffixes=('', '_2') 
    )

    # find significant overlaps
    matches['overlap_len'] = matches['overlap_end'] - matches['overlap_start']
    significant_matches = matches[matches['overlap_len'] >= min_overlap]

    # find tracking indices in df1 and df2 of significant overlaps
    matched_indices_1 = significant_matches["idx"].unique() 
    matched_indices_2 = significant_matches["idx_2"].unique()

    # construct a dataframe containing records unique to df1 (i.e., no significant match in df2), 
    # another one with records unique to df2, 
    # and a final one with records from df1 and df2 that significantly match 
    df1_unique = df1.filter(~pl.col("idx").is_in(matched_indices_1)).drop("idx")
    df2_unique = df2.filter(~pl.col("idx").is_in(matched_indices_2)).drop("idx")
    df_common = (
        pl
        .from_pandas(significant_matches)
        .drop(['idx', 'idx_2', 'overlap_start', 'overlap_end'])
    )

    return version_sort(df1_unique), version_sort(df2_unique), version_sort(df_common)

def find_unique_and_common_records_with_venn_diagram(df1, df2, min_overlap, labels):
    df1_unique, df2_unique, df_common = find_unique_and_common_records(df1, df2, min_overlap)

    plt.figure(figsize=(9, 9))
    venn2(
        subsets=(len(df1_unique), len(df2_unique), len(df_common)),
        set_labels=labels
    )
    plt.show()

    return (
        add_locus(df1_unique), 
        add_locus(df2_unique), 
        add_locus(df_common, col_name=f"locus_{labels[0]}")
    )
