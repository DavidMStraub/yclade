"""Types for yclade."""

from dataclasses import dataclass

import networkx as nx


Snp = str
CladeName = str
CladeSnps = dict[CladeName, set[Snp]]


@dataclass
class YTreeData:
    """Y tree data structure."""

    graph: nx.DiGraph
    snps: CladeSnps


@dataclass
class SnpResults:
    """A set of positive and negative Y SNP test results."""
    positive: set[Snp]
    negative: set[Snp]
