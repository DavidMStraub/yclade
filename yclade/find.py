"""Tools for finding the best sublade for a given set of SNPs."""

from yclade.types import CladeMatchInfo, CladeName, SnpResults, YTreeData


def find_nodes_with_positive_matches(
    tree: YTreeData, snps: SnpResults
) -> set[CladeName]:
    """Find the nodes in the tree that have at least one matching positive SNP."""
    nodes = set()
    for clade, clade_snps in tree.clade_snps.items():
        if snps.positive & clade_snps:
            nodes.add(clade)
    return nodes


def get_nodes_with_match_info(
    tree: YTreeData, snps: SnpResults
) -> dict[CladeName, CladeMatchInfo]:
    """Find the nodes in the tree that have overlap with postive or negative SNPs."""
    node_info = {}
    for clade, clade_snps in tree.clade_snps.items():
        if len(clade_snps) == 0:
            continue
        positive = len(snps.positive & clade_snps)
        negative = len(snps.negative & clade_snps)
        if positive or negative:
            node_info[clade] = CladeMatchInfo(
                positive=positive,
                negative=negative,
                length=len(clade_snps),
            )
    return node_info
