"""Find the position on the human Y chromosome tree for a set of SNPs."""

from importlib.metadata import PackageNotFoundError
from importlib.metadata import version as _get_version

from . import const, find, snps, tree, types
from ._yclade import find_clade
from .find import get_clade_lineage

try:
    __version__ = _get_version("yclade")
except PackageNotFoundError:  # not installed, e.g. running from a source checkout
    __version__ = "unknown"

__all__ = [
    "__version__",
    "const",
    "find",
    "find_clade",
    "get_clade_lineage",
    "snps",
    "tree",
    "types",
]
