"""Bump yclade to the latest YFull tree release.

Updates YTREE_DEFAULT_VERSION in yclade/const.py if a newer tree is available
upstream. Does nothing otherwise.
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
CONST_RE = re.compile(r'(?m)^(YTREE_DEFAULT_VERSION = ")([^"]+)(")$')


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
    # The package version is derived from the git tag by setuptools_scm, so
    # releasing the new tree is a matter of tagging the merged pull request.
    print(f"Bumped to YFull tree {latest}.")


if __name__ == "__main__":
    main()
