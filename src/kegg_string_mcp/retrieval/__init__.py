"""Vector retrieval arm, measured against the keyword tools already in the repo.

Imports lazily: this subpackage depends on the optional `[vector]` extra (torch,
chromadb, sentence-transformers), and importing the rest of the library must not
require any of it.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:  # pragma: no cover
    from kegg_string_mcp.retrieval.corpus import Corpus as Corpus
    from kegg_string_mcp.retrieval.corpus import Passage as Passage
    from kegg_string_mcp.retrieval.corpus import build as build

_LAZY = {"Corpus": "kegg_string_mcp.retrieval.corpus",
         "Passage": "kegg_string_mcp.retrieval.corpus",
         "build": "kegg_string_mcp.retrieval.corpus"}

__all__ = ["Corpus", "Passage", "build"]


def __getattr__(name: str) -> Any:
    if name not in _LAZY:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    from importlib import import_module

    return getattr(import_module(_LAZY[name]), name)
