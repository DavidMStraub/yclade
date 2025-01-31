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
    YTREE_URL,
    YTREE_ZIP_FILENAME,
)
from yclade.types import CladeSnps, YTreeData


def download_yfull_tree(
    version: str | None = None, data_dir: Path | None = None, force: bool = False
) -> None:
    """Download the YFull tree."""
    version = version or YTREE_DEFAULT_VERSION
    data_dir = data_dir or Path(user_data_dir("yclade"))
    data_dir.mkdir(parents=True, exist_ok=True)
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
    print(tree_data["id"])
    if not snps:
        snps = {}
    if "snps" in tree_data:
        if tree_data["snps"]:
            snps[tree_data["id"]] = set(tree_data["snps"].split(", "))
            print(snps)
        else:
            snps[tree_data["id"]] = set()
            print(snps)
    else:
        snps[tree_data["id"]] = set()
        print(snps)
    for child in tree_data.get("children", []):
        _get_clade_snps(child, snps)
    return snps


def yfull_tree_to_tree_data(file_path: Path) -> YTreeData:
    """Convert a YFull tree file to a tree data dictionary."""
    with open(file_path) as f:
        tree_data = json.load(f)
    graph = _build_graph(tree_data)
    snps = _get_clade_snps(tree_data)
    return YTreeData(graph=graph, snps=snps)
