"""Bump yclade to the latest YFull tree release.

Updates YTREE_DEFAULT_VERSION in yclade/const.py and the package version in
pyproject.toml if a newer tree is available upstream. Does nothing otherwise.
"""

from __future__ import annotations

import json
import os
import re
import urllib.request
from pathlib import Path

CONTENTS_URL = "https://api.github.com/repos/YFullTeam/YTree/contents/ytree"
FILENAME_RE = re.compile(r"^tree_(\d+\.\d+\.\d+)\.zip$")

ROOT = Path(__file__).resolve().parents[2]
CONST_PATH = ROOT / "yclade" / "const.py"
PYPROJECT_PATH = ROOT / "pyproject.toml"
CONST_RE = re.compile(r'(?m)^(YTREE_DEFAULT_VERSION = ")([^"]+)(")$')
PYPROJECT_RE = re.compile(r'(?m)^(version = ")(\d+)\.(\d+)\.(\d+)(")$')


def version_key(version: str) -> tuple[int, ...]:
    """Turn a version string into a tuple of ints for comparison."""
    return tuple(int(part) for part in version.split("."))


def latest_version() -> str:
    """Return the newest YFull tree version available as a ZIP file."""
    request = urllib.request.Request(CONTENTS_URL)
    token = os.environ.get("GH_TOKEN")
    if token:  # unauthenticated API requests are rate limited
        request.add_header("Authorization", f"Bearer {token}")
    with urllib.request.urlopen(request, timeout=60) as response:
        entries = json.load(response)
    versions = [m.group(1) for e in entries if (m := FILENAME_RE.match(e["name"]))]
    # Compare numerically, but keep the version as spelled by the file name: the
    # minor version is zero padded, e.g. 'tree_14.05.0.zip'.
    return max(versions, key=version_key)


def main() -> None:
    const = CONST_PATH.read_text()
    const_match = CONST_RE.search(const)
    assert const_match, f"no YTREE_DEFAULT_VERSION in {CONST_PATH}"
    current, latest = const_match.group(2), latest_version()
    print(f"Current YFull tree version: {current}, latest: {latest}")
    if version_key(latest) <= version_key(current):
        return
    CONST_PATH.write_text(CONST_RE.sub(rf"\g<1>{latest}\g<3>", const, count=1))
    pyproject = PYPROJECT_PATH.read_text()
    version_match = PYPROJECT_RE.search(pyproject)
    assert version_match, f"no package version in {PYPROJECT_PATH}"
    package_version = f"{version_match.group(2)}.{int(version_match.group(3)) + 1}.0"
    PYPROJECT_PATH.write_text(
        PYPROJECT_RE.sub(rf"\g<1>{package_version}\g<5>", pyproject, count=1)
    )
    print(f"Bumped to YFull tree {latest}, package version {package_version}.")


if __name__ == "__main__":
    main()
