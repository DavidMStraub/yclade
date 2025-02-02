"""Tools for finding the best sublade for a given set of SNPs."""

import logging
from collections.abc import Callable

import networkx as nx

from yclade.types import CladeInfo, CladeMatchInfo, CladeName, SnpResults, YTreeData


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


def get_node_match_info(
    tree: YTreeData, node: CladeName, snps: SnpResults
) -> CladeMatchInfo:
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
    """Get the match info for all nodes in the tree."""
    node_info = {}
    for clade, clade_snps in tree.clade_snps.items():
        if len(clade_snps) == 0:
            continue
        node_info[clade] = get_node_match_info(tree, clade, snps)
    return node_info


def get_node_path_scores(
    tree: YTreeData,
    node: CladeName,
    snps: SnpResults,
    scoring_function: Callable[[CladeMatchInfo], float],
) -> dict[CladeName, float]:
    """Get the score for a single node's path (from the root to the node).

    Args:
        tree: The YTreeData object.
        node: The ID of the node to score.
        snps: The SnpResults to use for scoring.
        scoring_function: The scoring function to use.

    The scoring function should take a CladeMatchInfo object and return a number.
    """
    scores = {}
    for ancestor_node in nx.ancestors(tree.graph, node):
        match_info = get_node_match_info(tree=tree, node=ancestor_node, snps=snps)
        scores[ancestor_node] = scoring_function(match_info)
    match_info = get_node_match_info(tree=tree, node=node, snps=snps)
    scores[node] = scoring_function(match_info)
    return scores


def simple_scoring_function(match_info: CladeMatchInfo) -> float:
    """Score a node based on the number of positive and negative SNPs."""
    return match_info.positive - match_info.negative


def get_ordered_clade_details(tree: YTreeData, snps: SnpResults) -> list[CladeInfo]:
    """Get an ordered list of clades with their match info.

    The first clade in the list is the best match.
    """
    warn_unknown_snps(tree=tree, snps=snps)
    candidates = find_nodes_with_positive_matches(tree=tree, snps=snps)
    candidate_scores = {}
    for node in candidates:
        scores = get_node_path_scores(
            tree=tree,
            node=node,
            snps=snps,
            scoring_function=simple_scoring_function,
        )
        total_score = sum(scores.values())
        candidate_scores[node] = total_score
    return [
        CladeInfo(
            name=node,
            age_info=tree.clade_age_infos[node],
            score=candidate_scores[node],
        )
        for node in sorted(candidates, key=lambda x: candidate_scores[x], reverse=True)
    ]


def warn_unknown_snps(tree: YTreeData, snps: SnpResults) -> None:
    """Log a warning if the query contains SNPs not contained in the Y tree."""
    known_snps = set(tree.snp_aliases.keys()) | set(
        snp for clade_snps in tree.clade_snps.values() for snp in clade_snps
    )
    unknown_snps = (snps.positive | snps.negative) - known_snps
    if unknown_snps:
        logger = logging.getLogger("yclade")
        logger.warning(
            "The SNP query contains %s SNP%s "
            "not contained in the Y tree: %s",
            len(unknown_snps),
            "s" if len(unknown_snps) > 1 else "",
            unknown_snps,
        )
