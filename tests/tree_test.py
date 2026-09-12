"""Unit tests for the yclade.tree module."""

import json
import shutil
import tempfile
import zipfile
from pathlib import Path

import networkx as nx
import pytest

import yclade.tree
from yclade.const import YTREE_JSON_FILENAME, YTREE_ZIP_FILENAME


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


def test_build_graph_with_empty_root_id():
    """The YFull tree's root has an empty ID but must still be connected."""
    raw_data = {"id": "", "children": [{"id": "A00", "children": []}]}
    graph = yclade.tree._build_graph(raw_data)
    assert set(graph.edges) == {("", "A00")}
    assert [node for node, degree in graph.in_degree() if degree == 0] == [""]


def test_get_actual_data_dir_accepts_a_string(tmp_path):
    """A data directory given as a string is usable, e.g. from the command line."""
    data_dir = yclade.tree._get_actual_data_dir(str(tmp_path / "data"))
    assert data_dir == tmp_path / "data"
    assert data_dir.is_dir()


@pytest.fixture
def data_dir_with_zip(tmp_path, tree_data):
    """A data directory holding only a tree ZIP file, as if never extracted."""
    version = "0.00.0"
    zip_path = tmp_path / YTREE_ZIP_FILENAME.format(version=version)
    with zipfile.ZipFile(zip_path, "w") as zip_ref:
        zip_ref.writestr(
            YTREE_JSON_FILENAME.format(version=version), json.dumps(tree_data)
        )
    return version, tmp_path


def test_get_yfull_tree_data_extracts_cached_zip(data_dir_with_zip):
    """A cached ZIP file is extracted even though nothing is downloaded."""
    version, data_dir = data_dir_with_zip
    tree_data_object = yclade.tree.get_yfull_tree_data(
        version=version, data_dir=data_dir
    )
    assert tree_data_object.version == version
    assert "B" in tree_data_object.graph
    assert (data_dir / YTREE_JSON_FILENAME.format(version=version)).exists()


def test_get_yfull_tree_data_is_cached_on_request(data_dir_with_zip):
    """Parsing is expensive, so the same tree can be reused."""
    version, data_dir = data_dir_with_zip
    first = yclade.tree.get_yfull_tree_data(
        version=version, data_dir=data_dir, cache=True
    )
    second = yclade.tree.get_yfull_tree_data(
        version=version, data_dir=data_dir, cache=True
    )
    assert first is second


def test_get_yfull_tree_data_is_not_cached_by_default(data_dir_with_zip):
    """A parsed tree needs a lot of memory, so it is not kept around uninvited."""
    version, data_dir = data_dir_with_zip
    first = yclade.tree.get_yfull_tree_data(version=version, data_dir=data_dir)
    second = yclade.tree.get_yfull_tree_data(version=version, data_dir=data_dir)
    assert first is not second


def test_download_yfull_tree_replaces_a_leftover_json(tmp_path, tree_data, monkeypatch):
    """A newly downloaded archive wins over a JSON file left from an earlier one."""
    version = "0.00.0"
    served_zip = tmp_path / "served.zip"
    with zipfile.ZipFile(served_zip, "w") as zip_ref:
        zip_ref.writestr(
            YTREE_JSON_FILENAME.format(version=version), json.dumps(tree_data)
        )
    data_dir = tmp_path / "data"
    data_dir.mkdir()
    json_path = data_dir / YTREE_JSON_FILENAME.format(version=version)
    json_path.write_text(json.dumps({"id": "stale", "children": []}))
    monkeypatch.setattr(
        yclade.tree.urllib.request,
        "urlretrieve",
        lambda url, path: shutil.copy(served_zip, path),
    )

    tree_data_object = yclade.tree.get_yfull_tree_data(
        version=version, data_dir=data_dir
    )

    assert "stale" not in tree_data_object.graph
    assert "B" in tree_data_object.graph
