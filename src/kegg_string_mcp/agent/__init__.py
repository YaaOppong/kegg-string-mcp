"""Agent layer.

Submodules are resolved lazily (PEP 562). `store` and `validate` depend on nothing
but the standard library, while `pipeline` pulls in anthropic, httpx and mcp -- and
an eager re-export here meant that importing the validator dragged the whole HTTP
and model stack in with it. That cost the demo its ability to run anywhere those
packages are awkward (a browser under Pyodide, most obviously), for no benefit.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:  # pragma: no cover - re-exported for type checkers only
    from kegg_string_mcp.agent.pipeline import Tools as Tools
    from kegg_string_mcp.agent.pipeline import annotate_epistasis as annotate_epistasis
    from kegg_string_mcp.agent.pipeline import annotate_gene as annotate_gene
    from kegg_string_mcp.agent.store import RunStore as RunStore
    from kegg_string_mcp.agent.validate import ValidationReport as ValidationReport
    from kegg_string_mcp.agent.validate import validate as validate

_LAZY = {
    "Tools": "kegg_string_mcp.agent.pipeline",
    "annotate_epistasis": "kegg_string_mcp.agent.pipeline",
    "annotate_gene": "kegg_string_mcp.agent.pipeline",
    "RunStore": "kegg_string_mcp.agent.store",
    "ValidationReport": "kegg_string_mcp.agent.validate",
    "validate": "kegg_string_mcp.agent.validate",
}

__all__ = [
    "RunStore",
    "Tools",
    "ValidationReport",
    "annotate_epistasis",
    "annotate_gene",
    "validate",
]


def __getattr__(name: str) -> Any:
    if name not in _LAZY:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    from importlib import import_module

    return getattr(import_module(_LAZY[name]), name)


def __dir__() -> list[str]:
    return __all__
