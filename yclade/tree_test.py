"""Unit tests for the yclade.tree module."""

import json
import tempfile
from pathlib import Path

import networkx as nx
import pytest

import yclade.tree


@pytest.fixture
def tree_data():
    return {
        "id": "root",
        "children": [
            {
                "id": "A",
                "snps": "",
                "children": [
                    {
                        "id": "B",
                        "snps": "s1, s2",
                        "children": [
                            {"id": "C", "children": []},
                            {"id": "D", "children": []},
                        ],
                    }
                ],
            },
            {
                "id": "E",
                "snps": "s3",
                "children": [
                    {"id": "F", "children": [{"id": "G", "children": []}]},
                    {"id": "H", "children": []},
                ],
            },
        ],
    }


def test_build_graph(tree_data):
    graph = yclade.tree._build_graph(tree_data)
    assert set(graph.nodes) == {
        "root",
        "A",
        "B",
        "C",
        "D",
        "E",
        "F",
        "G",
        "H",
    }
    assert list(nx.dfs_preorder_nodes(graph, source="root")) == [
        "root",
        "A",
        "B",
        "C",
        "D",
        "E",
        "F",
        "G",
        "H",
    ]


def test_yfull_tree_to_tree_data(tree_data):
    content = json.dumps(tree_data)
    with tempfile.NamedTemporaryFile("w") as f:
        f.write(content)
        f.seek(0)
        tree_data_object = yclade.tree.yfull_tree_to_tree_data(Path(f.name))
    assert set(tree_data_object.graph.nodes) == {
        "root",
        "A",
        "B",
        "C",
        "D",
        "E",
        "F",
        "G",
        "H",
    }
    assert tree_data_object.clade_snps == {
        "root": set(),
        "A": set(),
        "B": {"s1", "s2"},
        "C": set(),
        "D": set(),
        "E": {"s3"},
        "F": set(),
        "G": set(),
        "H": set(),
    }
