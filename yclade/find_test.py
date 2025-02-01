"""Unit tests for the yclade.find module."""

import pytest

import yclade.find
import yclade.snps
from yclade.tree import _build_graph, _get_clade_snps
from yclade.types import CladeMatchInfo, YTreeData


@pytest.fixture
def tree_data():
    raw_data = {
        "id": "root",
        "children": [
            {
                "id": "A",
                "snps": "a, b, c",
                "children": [
                    {
                        "id": "B",
                        "snps": "a, b, d, e, f",
                        "children": [
                            {"id": "C", "snps": "a, b, d, e, f, g", "children": []},
                        ],
                    }
                ],
            }
        ],
    }
    graph = _build_graph(raw_data)
    clade_snps = _get_clade_snps(raw_data)
    return YTreeData(graph=graph, clade_snps=clade_snps)


def test_find_nodes_with_positive_matches(tree_data):
    snps = yclade.snps.parse_snp_results("a+,b-,d+")
    nodes = yclade.find.find_nodes_with_positive_matches(tree_data, snps)
    assert nodes == {"A", "B", "C"}

def test_find_nodes_with_positive_matches_unknown(tree_data):
    snps = yclade.snps.parse_snp_results("x+")
    nodes = yclade.find.find_nodes_with_positive_matches(tree_data, snps)
    assert nodes == set()

def test_find_nodes_with_positive_matches_single(tree_data):
    snps = yclade.snps.parse_snp_results("g+")
    nodes = yclade.find.find_nodes_with_positive_matches(tree_data, snps)
    assert nodes == {"C"}

def test_get_nodes_with_match_info(tree_data):
    snps = yclade.snps.parse_snp_results("a+,b-,d+")
    node_info = yclade.find.get_nodes_with_match_info(tree_data, snps)
    assert node_info == {
        "A": CladeMatchInfo(positive=1, negative=1, length=3),
        "B": CladeMatchInfo(positive=2, negative=1, length=5),
        "C": CladeMatchInfo(positive=2, negative=1, length=6),
    }