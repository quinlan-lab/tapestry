"""Shared Tapestry utilities."""

import polars as pl


def read_all_cpgs_in_reference(bed):
    return pl.read_csv(
        bed,
        separator="\t",
        has_header=False,
        new_columns=["chrom", "start", "end"],
    )


__all__ = ["read_all_cpgs_in_reference"]
