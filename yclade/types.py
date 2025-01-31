"""Types for yclade."""

from dataclasses import dataclass

import networkx as nx


SNP = str
CLADE_NAME = str
CladeSnps = dict[CLADE_NAME, set[SNP]]


@dataclass
class YTreeData:
    """Y tree data structure."""

    graph: nx.DiGraph
    snps: CladeSnps
