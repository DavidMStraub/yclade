"""Main functions for the yclade package."""

from __future__ import annotations

from pathlib import Path

from . import find, snps, tree, types


def find_clade(
    snp_string: str,
    version: str | None = None,
    data_dir: Path | str | None = None,
    cache: bool = False,
) -> list[types.CladeInfo]:
    """Find the clade that matches the given SNP string.

    Set `cache` to True to keep the parsed tree in memory, which speeds up
    repeated calls at the cost of a large amount of memory. See
    `yclade.tree.get_yfull_tree_data`.
    """
    snp_results = snps.parse_snp_results(snp_string)
    tree_data = tree.get_yfull_tree_data(
        version=version, data_dir=data_dir, cache=cache
    )
    snp_results = snps.normalize_snp_results(
        snp_results=snp_results,
        snp_aliases=tree_data.snp_aliases,
    )
    return find.get_ordered_clade_details(tree=tree_data, snps=snp_results)
