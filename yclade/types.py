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

    formed_confidence_interval: tuple[float, float]  | None
    """95 % confidence interval for the formed age."""

    most_recent_common_ancestor: float | None
    """How many years ago the most recent common ancestor was born."""

    most_recent_common_ancestor_confidence_interval: tuple[float, float] | None
    """95 % confidence interval for the most recent common ancestor age."""


CladeAgeInfos = dict[CladeName, CladeAgeInfo]


@dataclass
class YTreeData:
    """Y tree data structure."""

    graph: nx.DiGraph
    clade_snps: CladeSnps
    clade_age_infos: CladeAgeInfos
    snp_aliases: dict[Snp, Snp]


@dataclass
class SnpResults:
    """A set of positive and negative Y SNP test results."""

    positive: set[Snp]
    negative: set[Snp]


@dataclass
class CladeMatchInfo:
    """A data type containing the number of positive and negative SNPs matched."""

    positive: int
    negative: int
    length: int


@dataclass
class CladeInfo:
    """A data type containing info about a specific clade."""

    name: str
    """The ID of the clade."""

    age_info: CladeAgeInfo
    """The age information for the clade."""

    score: float | None = None
    """The score for the clade (optional)."""
