"""Utilities to handle the Y tree data."""

from __future__ import annotations

import json
import logging
import urllib.request
import zipfile
from pathlib import Path

import networkx as nx
from platformdirs import user_data_dir

from yclade.const import (
    YTREE_DEFAULT_VERSION,
    YTREE_JSON_FILENAME,
    YTREE_URL,
    YTREE_ZIP_FILENAME,
)
from yclade.types import CladeSnps, YTreeData, Snp, CladeAgeInfos, CladeAgeInfo


def _get_actual_data_dir(data_dir: Path | None = None) -> Path:
    """Get the data directory, creating it if it doesn't exist.

    Args:
        data_dir: The data directory path. If None, the default user data diretory
            is used, e.g. '~/.local/share/yclade' on Linux.
    """
    data_dir = data_dir or Path(user_data_dir("yclade"))
    data_dir.mkdir(parents=True, exist_ok=True)
    return data_dir


def download_yfull_tree(
    version: str | None = None, data_dir: Path | None = None, force: bool = False
) -> None:
    """Download the YFull tree.

    Args:
        version: The version of the YFull tree to download. If None, the default
            version is used.
        data_dir: The directory to store the YFull data. If None, the default user
            data directory is used.
        force: If True, the tree is downloaded even if it already exists.
    """
    version = version or YTREE_DEFAULT_VERSION
    data_dir = _get_actual_data_dir(data_dir)
    file_path = data_dir / YTREE_ZIP_FILENAME.format(version=version)
    if file_path.exists() and not force:
        logging.info("YFull tree already downloaded to %s", file_path)
        return
    url = YTREE_URL.format(version=version)
    urllib.request.urlretrieve(url, file_path)
    logging.info("Downloaded YFull tree to %s", file_path)
    with zipfile.ZipFile(file_path, "r") as zip_ref:
        zip_ref.extractall(data_dir)


def _build_graph(
    node, graph: nx.DiGraph | None = None, parent: str | None = None
) -> nx.DiGraph:
    """Recursively build a DiGraph from the tree data."""
    if graph is None:
        graph = nx.DiGraph()
    graph.add_node(node["id"])
    if parent:
        graph.add_edge(parent, node["id"])
    for child in node.get("children", []):
        _build_graph(node=child, graph=graph, parent=node["id"])
    return graph


def _get_clade_snps(tree_data, snps: CladeSnps | None = None) -> CladeSnps:
    """Recursively get the dictionary of clade SNPs from the tree data."""
    if not snps:
        snps = {}
    if "snps" in tree_data:
        if tree_data["snps"]:
            snps[tree_data["id"]] = set(tree_data["snps"].split(", "))
        else:
            snps[tree_data["id"]] = set()
    else:
        snps[tree_data["id"]] = set()
    for child in tree_data.get("children", []):
        _get_clade_snps(child, snps)
    return snps


def _clade_snps_to_snp_aliases(clade_snps: CladeSnps) -> dict[Snp, Snp]:
    """Create a dictionary of SNP aliases from the clade SNPs."""
    snp_aliases = {}
    for snps in clade_snps.values():
        for snp in snps:
            if "/" in snp:
                alias1, alias2 = snp.split("/", 1)
                snp_aliases[alias1] = snp
                snp_aliases[alias2] = snp
    return snp_aliases


def _get_clade_age_infos(tree_data, age_infos=None) -> CladeAgeInfos:
    """Recursively get the clade age information from the tree data."""
    if age_infos is None:
        age_infos = {}
    if "formed" not in tree_data:
        return age_infos
    age_infos[tree_data["id"]] = CladeAgeInfo(
        formed=None if tree_data["formed"] == "-" else tree_data["formed"],
        formed_confidence_interval=(
            None
            if tree_data["formedlowage"] == "-" or tree_data["formedhighage"] == "-"
            else (
                tree_data["formedlowage"],
                tree_data["formedhighage"],
            )
        ),
        most_recent_common_ancestor=tree_data["tmrca"],
        most_recent_common_ancestor_confidence_interval=(
            None
            if tree_data["tmrcalowage"] == "-" or tree_data["tmrcahighage"] == "-"
            else (
                tree_data["tmrcalowage"],
                tree_data["tmrcahighage"],
            )
        ),
    )
    for child in tree_data.get("children", []):
        _get_clade_age_infos(child, age_infos)
    return age_infos


def yfull_tree_to_tree_data(file_path: Path) -> YTreeData:
    """Convert a YFull tree file to a tree data dictionary.

    Args:
        file_path: The path to the YFull tree JSON file.
    """
    with open(file_path) as f:
        tree_data = json.load(f)
    graph = _build_graph(tree_data)
    clade_snps = _get_clade_snps(tree_data)
    snp_aliases = _clade_snps_to_snp_aliases(clade_snps)
    clade_age_infos = _get_clade_age_infos(tree_data)
    return YTreeData(
        graph=graph,
        clade_snps=clade_snps,
        clade_age_infos=clade_age_infos,
        snp_aliases=snp_aliases,
    )


def get_yfull_tree_data(
    version: str | None = None, data_dir: Path | None = None
) -> YTreeData:
    """Get the YFull tree data from cache, download if needed.

    Args:
        version: The version of the YFull tree to download. If None, the default
            version is used.
        data_dir: The directory to store the YFull data. If None, the default user
            data directory is used.
    """
    download_yfull_tree(version=version, data_dir=data_dir, force=False)
    data_dir = _get_actual_data_dir(data_dir)
    version = version or YTREE_DEFAULT_VERSION
    file_path = data_dir / YTREE_JSON_FILENAME.format(version=version)
    return yfull_tree_to_tree_data(file_path)
