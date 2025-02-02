"""Main functions for the yclade package."""

from __future__ import annotations

from pathlib import Path

from . import find, snps, tree, types


def find_clade(
    snp_string: str, version: str | None = None, data_dir: Path | None = None
) -> list[types.CladeInfo]:
    """Find the clade that matches the given SNP string."""
    snp_results = snps.parse_snp_results(snp_string)
    tree_data = tree.get_yfull_tree_data(version=version, data_dir=data_dir)
    snp_results = snps.normalize_snp_results(
        snp_results=snp_results,
        snp_aliases=tree_data.snp_aliases,
    )
    return find.get_ordered_clade_details(tree=tree_data, snps=snp_results)
