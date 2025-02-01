"""Tools for finding the best sublade for a given set of SNPs."""

import networkx as nx

from yclade.types import CladeMatchInfo, CladeName, SnpResults, YTreeData, Snp


def find_nodes_with_positive_matches(
    tree: YTreeData, snps: SnpResults
) -> set[CladeName]:
    """Find the nodes in the tree that have at least one matching positive SNP."""
    nodes = set()
    positives = set(tree.snp_aliases.get(snp, snp) for snp in snps.positive)
    for clade, clade_snps in tree.clade_snps.items():
        if positives & clade_snps:
            nodes.add(clade)
    return nodes


def get_node_match_info(tree: YTreeData, node: CladeName, snps: SnpResults) -> CladeMatchInfo:
    """Get the match info for a single node."""
    clade_snps = tree.clade_snps[node]
    positives = set(tree.snp_aliases.get(snp, snp) for snp in snps.positive)
    negatives = set(tree.snp_aliases.get(snp, snp) for snp in snps.negative)
    clade_positives = len(positives & clade_snps)
    clade_negatives = len(negatives & clade_snps)
    return CladeMatchInfo(
        positive=clade_positives,
        negative=clade_negatives,
        length=len(clade_snps),
    )

def get_all_nodes_match_info(
    tree: YTreeData, snps: SnpResults
) -> dict[CladeName, CladeMatchInfo]:
    """Find the nodes in the tree that have overlap with postive or negative SNPs."""
    node_info = {}
    for clade, clade_snps in tree.clade_snps.items():
        if len(clade_snps) == 0:
            continue
        node_info[clade] = get_node_match_info(tree, clade, snps)
    return node_info
