"""Types for yclade."""

from __future__ import annotations

from dataclasses import dataclass

import networkx as nx

Snp = str
CladeName = str
CladeSnps = dict[CladeName, set[Snp]]


@dataclass
class CladeAgeInfo:
    """A data type containing estimated age information."""

    formed: float | None
    """How many years ago the clade was formed."""

    formed_confidence_interval: tuple[float, float] | None
    """95 % confidence interval for the formed age."""

    most_recent_common_ancestor: float | None
    """How many years ago the most recent common ancestor was born."""

    most_recent_common_ancestor_confidence_interval: tuple[float, float] | None
    """95 % confidence interval for the most recent common ancestor age."""


CladeAgeInfos = dict[CladeName, CladeAgeInfo]
SnpAliases = dict[Snp, Snp]


@dataclass
class YTreeData:
    """A data type representing the information in the Y tree."""

    graph: nx.DiGraph
    """A networkx DiGraph representing the tree's nodes as clade IDs."""

    clade_snps: CladeSnps
    """A dictionary of clade IDs to sets of SNPs."""

    clade_age_infos: CladeAgeInfos
    """A dictionary of clade IDs to estimated age information."""

    snp_aliases: SnpAliases
    """A dictionary of SNP aliases to their canonical form."""


@dataclass
class SnpResults:
    """A set of positive and negative Y SNP test results."""

    positive: set[Snp]
    """The set of positively tested SNPs."""

    negative: set[Snp]
    """The set of negatively tested SNPs."""


@dataclass
class CladeMatchInfo:
    """A data type containing the number of positive and negative SNPs matched."""

    positive: int
    """Number of positively tested SNPs that are present in the clade."""

    negative: int
    """Number of negatively tested SNPs that are present in the clade."""

    length: int
    """Total number of SNPs present in the clade."""


@dataclass
class CladeInfo:
    """A data type containing info about a specific clade."""

    name: str
    """The ID of the clade."""

    age_info: CladeAgeInfo
    """The age information for the clade."""

    score: float | None = None
    """The score for the clade (optional)."""
