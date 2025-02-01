"""Unit tests for the yclade.tree module."""

import json
import tempfile
from pathlib import Path

import networkx as nx
import pytest

import yclade.tree


def _age_info(formed):
    """Create dummy age info."""
    return {
        "formed": formed,
        "formedhighage": formed + 100,
        "formedlowage": formed - 100,
        "tmrca": formed + 200,
        "tmrcahighage": formed + 300,
        "tmrcalowage": formed + 100,
    }


@pytest.fixture
def tree_data():
    return {
        "id": "root",
        **_age_info(22000),
        "children": [
            {
                "id": "A",
                **_age_info(21000),
                "snps": "",
                "children": [
                    {
                        "id": "B",
                        **_age_info(20000),
                        "snps": "s1/t1, s2",
                        "children": [
                            {"id": "C", **_age_info(18000), "children": []},
                            {"id": "D", **_age_info(18000), "children": []},
                        ],
                    }
                ],
            },
            {
                "id": "E",
                **_age_info(20000),
                "snps": "s3",
                "children": [
                    {
                        "id": "F",
                        **_age_info(12000),
                        "children": [{"id": "G", "children": []}],
                    },
                    {"id": "H", **_age_info(13000), "children": []},
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
        "B": {"s1/t1", "s2"},
        "C": set(),
        "D": set(),
        "E": {"s3"},
        "F": set(),
        "G": set(),
        "H": set(),
    }
    assert tree_data_object.snp_aliases == {"s1": "s1/t1", "t1": "s1/t1"}
    assert tree_data_object.clade_age_infos["root"].formed == 22000
    assert tree_data_object.clade_age_infos["root"].formed_confidence_interval == (
        21900,
        22100,
    )
    assert tree_data_object.clade_age_infos["root"].most_recent_common_ancestor == 22200
    assert tree_data_object.clade_age_infos[
        "root"
    ].most_recent_common_ancestor_confidence_interval == (22100, 22300)
