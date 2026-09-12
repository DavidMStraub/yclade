"""Tools for finding the best sublade for a given set of SNPs."""

import logging
from collections.abc import Callable

import networkx as nx

from yclade.types import (
    CladeInfo,
    CladeMatchInfo,
    CladeName,
    Snp,
    SnpResults,
    YTreeData,
)


def _canonical_snps(tree: YTreeData, snps: SnpResults) -> tuple[set[Snp], set[Snp]]:
    """Cast the tested SNPs to their canonical form."""
    return (
        {tree.snp_aliases.get(snp, snp) for snp in snps.positive},
        {tree.snp_aliases.get(snp, snp) for snp in snps.negative},
    )


def _match_info(
    clade_snps: set[Snp], positives: set[Snp], negatives: set[Snp]
) -> CladeMatchInfo:
    """Match a clade's SNPs against already canonicalized test results."""
    return CladeMatchInfo(
        positive=len(positives & clade_snps),
        negative=len(negatives & clade_snps),
        length=len(clade_snps),
    )


def find_nodes_with_positive_matches(
    tree: YTreeData, snps: SnpResults
) -> set[CladeName]:
    """Find the nodes in the tree that have at least one matching positive SNP."""
    positives, _ = _canonical_snps(tree, snps)
    return {
        clade
        for clade, clade_snps in tree.clade_snps.items()
        if positives & clade_snps
    }


def get_node_match_info(
    tree: YTreeData, node: CladeName, snps: SnpResults
) -> CladeMatchInfo:
    """Get the match info for a single node."""
    return _match_info(tree.clade_snps[node], *_canonical_snps(tree, snps))


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
    positives, negatives = _canonical_snps(tree, snps)
    scores: dict[CladeName, float] = {}

    def score(node: CladeName) -> float:
        """Score a single clade, reusing the scores of already seen clades."""
        if node not in scores:
            scores[node] = simple_scoring_function(
                _match_info(tree.clade_snps[node], positives, negatives)
            )
        return scores[node]

    candidates = find_nodes_with_positive_matches(tree=tree, snps=snps)
    # Candidates share most of their ancestors, so scoring each clade once and
    # reusing the result is much faster than scoring every path independently.
    candidate_scores = {
        node: score(node) + sum(score(a) for a in nx.ancestors(tree.graph, node))
        for node in candidates
    }
    return [
        CladeInfo(
            name=node,
            age_info=tree.clade_age_infos.get(node),
            score=candidate_scores[node],
        )
        for node in sorted(candidates, key=lambda x: candidate_scores[x], reverse=True)
    ]


def warn_unknown_snps(tree: YTreeData, snps: SnpResults) -> None:
    """Log a warning if the query contains SNPs not contained in the Y tree."""
    # Whittle down the tested SNPs instead of collecting every SNP in the tree,
    # which needs a lot of time and memory for a large query.
    unknown_snps = (snps.positive | snps.negative) - tree.snp_aliases.keys()
    for clade_snps in tree.clade_snps.values():
        if not unknown_snps:
            break
        unknown_snps -= clade_snps
    if unknown_snps:
        logger = logging.getLogger("yclade")
        logger.warning(
            "The SNP query contains %s SNP%s not contained in the Y tree: %s",
            len(unknown_snps),
            "s" if len(unknown_snps) > 1 else "",
            unknown_snps,
        )


def _get_ancestors_ordered(graph: nx.DiGraph, node: str) -> list[str]:
    """Get the ancestors of a node in order from the root to the node."""
    ancestors = []
    visited = {node}  # Prevent the start node appearing as its own ancestor if the graph has a cycle.

    def depth_first_search(n):
        for parent in graph.predecessors(n):
            if parent not in visited:
                visited.add(parent)
                depth_first_search(parent)
                ancestors.append(parent)

    depth_first_search(node)
    return ancestors


def get_clade_lineage(tree: YTreeData, node: str) -> list[CladeInfo]:
    """Get the ancestors of a clade, including the node itself, as a list of
    `CladeInfo`.

    The clades are ordered from root to leaf."""
    # The root of the YFull tree has an empty ID and is not a real clade.
    names = [name for name in _get_ancestors_ordered(tree.graph, node) if name]
    names.append(node)
    clade_age_infos = [tree.clade_age_infos.get(clade) for clade in names]
    return [
        CladeInfo(name=clade, age_info=age_info)
        for clade, age_info in zip(names, clade_age_infos)
    ]
