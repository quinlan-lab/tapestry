"""Shared Tapestry utilities.

The repository historically had both ``src/util.py`` and a ``src/util/``
directory. Making the directory an explicit package keeps imports such as
``util.hap_map`` deterministic while retaining the legacy
``from util import read_all_cpgs_in_reference`` API.
"""

import polars as pl


def read_all_cpgs_in_reference(bed):
    return pl.read_csv(
        bed,
        separator="\t",
        has_header=False,
        new_columns=["chrom", "start", "end"],
    )


__all__ = ["read_all_cpgs_in_reference"]
