"""Tools for finding the best sublade for a given set of SNPs."""

from yclade.types import SnpResults, YTreeData, CladeName

def find_nodes_with_positive_matches(tree: YTreeData, snps: SnpResults) -> set[CladeName]:
    """Find the nodes in the tree that have at least one matching positive SNP."""
    nodes = set()
    for clade, clade_snps in tree.snps.items():
        if snps.positive & clade_snps:
            nodes.add(clade)
    return nodes