"""Tests for the way yclade is used by downstream consumers.

These mirror the call pattern of the Gramps Web API, the main consumer of
yclade, to catch changes that keep the signatures intact but break callers.
"""

import json
from dataclasses import asdict

import pytest

import yclade.find
import yclade.snps
import yclade.tree


@pytest.fixture
def tree_file(tmp_path):
    """A small YFull style tree file, with an empty root ID as in the real tree."""
    raw_data = {
        "id": "",
        "formed": "-",
        "formedlowage": "-",
        "formedhighage": "-",
        "tmrca": 235900,
        "tmrcalowage": 200000,
        "tmrcahighage": 300000,
        "snps": "",
        "children": [
            {
                "id": "A-1",
                "formed": 20000,
                "formedlowage": 19000,
                "formedhighage": 21000,
                "tmrca": 18000,
                "tmrcalowage": 17000,
                "tmrcahighage": 19000,
                "snps": "a1, a2",
                "children": [
                    {
                        "id": "A-2",
                        "formed": 18000,
                        "formedlowage": 17000,
                        "formedhighage": 19000,
                        "tmrca": 16000,
                        "tmrcalowage": 15000,
                        "tmrcahighage": 17000,
                        "snps": "b1, b2/b3",
                        "children": [],
                    }
                ],
            }
        ],
    }
    file_path = tmp_path / "tree_1.00.0.json"
    file_path.write_text(json.dumps(raw_data))
    return file_path


def test_gramps_web_call_pattern(tree_file):
    """The sequence of calls the Gramps Web API makes for a Y-DNA request."""
    snp_results = yclade.snps.parse_snp_results("a1+, a2+, b1+, b3+")
    tree_data = yclade.tree.yfull_tree_to_tree_data(tree_file, version="1.00.0")
    snp_results = yclade.snps.normalize_snp_results(
        snp_results=snp_results,
        snp_aliases=tree_data.snp_aliases,
    )
    ordered_clade_details = yclade.find.get_ordered_clade_details(
        tree=tree_data, snps=snp_results
    )
    assert len(ordered_clade_details) > 0
    most_likely_clade = ordered_clade_details[0].name
    assert most_likely_clade == "A-2"
    clade_lineage = yclade.find.get_clade_lineage(
        tree=tree_data, node=most_likely_clade
    )
    # The lineage is returned to API clients, so the empty root must stay out of it.
    assert [clade_info.name for clade_info in clade_lineage] == ["A-1", "A-2"]
    assert tree_data.version == "1.00.0"


def test_clade_info_is_serializable(tree_file):
    """Consumers serialize the results with dataclasses.asdict."""
    tree_data = yclade.tree.yfull_tree_to_tree_data(tree_file, version="1.00.0")
    clade_lineage = yclade.find.get_clade_lineage(tree=tree_data, node="A-2")
    serialized = [asdict(clade_info) for clade_info in clade_lineage]
    assert serialized[-1] == {
        "name": "A-2",
        "age_info": {
            "formed": 18000,
            "formed_confidence_interval": (17000, 19000),
            "most_recent_common_ancestor": 16000,
            "most_recent_common_ancestor_confidence_interval": (15000, 17000),
        },
        "score": None,
    }
    assert json.dumps(serialized)


def test_no_matching_clades_returns_an_empty_list(tree_file):
    """Consumers check the result for emptiness before taking the first clade."""
    snp_results = yclade.snps.parse_snp_results("unknown1+, unknown2-")
    tree_data = yclade.tree.yfull_tree_to_tree_data(tree_file, version="1.00.0")
    ordered_clade_details = yclade.find.get_ordered_clade_details(
        tree=tree_data, snps=snp_results
    )
    assert ordered_clade_details == []
