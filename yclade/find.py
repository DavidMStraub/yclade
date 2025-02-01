"""Tools for finding the best sublade for a given set of SNPs."""

from yclade.types import CladeMatchInfo, CladeName, SnpResults, YTreeData


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


def get_nodes_with_match_info(
    tree: YTreeData, snps: SnpResults
) -> dict[CladeName, CladeMatchInfo]:
    """Find the nodes in the tree that have overlap with postive or negative SNPs."""
    node_info = {}
    positives = set(tree.snp_aliases.get(snp, snp) for snp in snps.positive)
    negatives = set(tree.snp_aliases.get(snp, snp) for snp in snps.negative)
    for clade, clade_snps in tree.clade_snps.items():
        if len(clade_snps) == 0:
            continue
        clade_positives = len(positives & clade_snps)
        clade_negatives = len(negatives & clade_snps)
        if clade_positives or clade_negatives:
            node_info[clade] = CladeMatchInfo(
                positive=clade_positives,
                negative=clade_negatives,
                length=len(clade_snps),
            )
    return node_info
