"""Thin JSON boundary between the browser UI and the real replay layer.

The static build runs on bare Pyodide with no packages installed at all. An
earlier attempt used Gradio-Lite, which installs gradio through micropip and so
drags in huggingface-hub -- a package with no pure-Python wheel, which fails at
load. Nothing here imports anything outside the standard library, so there is no
dependency resolution to go wrong.

The UI is hand-written HTML and JavaScript, but every verdict on the page comes
from this module calling the same `validate()` the library uses. The browser
formats; it never decides.
"""

from __future__ import annotations

import json

from replay import LABELS, ORDERED, available, citation_rows, load, quote_rows, tool_call_rows


def index() -> str:
    """The gene picker. `fails` drives the warning marker, so which runs catch
    something is computed here rather than hard-coded in the page."""
    entries = []
    for name in ORDERED:
        if name not in available():
            continue
        entries.append({"id": name, "label": LABELS.get(name, name),
                        "fails": not load(name).clean})
    return json.dumps(entries)


def run(name: str) -> str:
    replay = load(name)
    model_calls = sum(1 for c in replay.calls if c.get("origin", "").startswith("turn"))
    return json.dumps({
        "target": replay.target,
        "mode": replay.mode,
        "clean": replay.clean,
        "model_calls": model_calls,
        "pipeline_calls": len(replay.calls) - model_calls,
        "turns": len(replay.turns),
        "calls": tool_call_rows(replay),
        "summary": replay.summary,
        "citations": citation_rows(replay),
        "quotes": quote_rows(replay),
        "n_failures": len(replay.failures),
    })
