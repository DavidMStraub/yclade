"""Unit tests for the yclade.find module."""

import pytest

from yclade import tree
import yclade.find
import yclade.snps
from yclade.tree import _build_graph, _get_clade_snps, _clade_snps_to_snp_aliases
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
                        "snps": "a, b, d, e, f/z",
                        "children": [
                            {"id": "C", "snps": "a, b, d, e, f/z, g", "children": []},
                        ],
                    }
                ],
            }
        ],
    }
    graph = _build_graph(raw_data)
    clade_snps = _get_clade_snps(raw_data)
    snp_aliases = _clade_snps_to_snp_aliases(clade_snps)
    return YTreeData(graph=graph, clade_snps=clade_snps, snp_aliases=snp_aliases)


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

def test_find_nodes_with_positive_matches_alias(tree_data):
    snps = yclade.snps.parse_snp_results("f+")
    nodes = yclade.find.find_nodes_with_positive_matches(tree_data, snps)
    assert nodes == {"B", "C"}
    snps = yclade.snps.parse_snp_results("z+")
    nodes = yclade.find.find_nodes_with_positive_matches(tree_data, snps)
    assert nodes == {"B", "C"}
    snps = yclade.snps.parse_snp_results("f/z+")
    nodes = yclade.find.find_nodes_with_positive_matches(tree_data, snps)
    assert nodes == {"B", "C"}