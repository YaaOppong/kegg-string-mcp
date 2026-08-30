"""Generate a static, serverless build of the demo for GitHub Pages.

Gradio-Lite runs Gradio in the browser under Pyodide, so the page needs no server
at all -- which suits a portfolio link better than a hosted Space that sleeps and
makes a visitor wait thirty seconds for a cold start.

This is only possible because the replay layer is standard-library only: nothing
to `micropip install`, so the page loads without fetching wheels. `agent/__init__`
imports lazily for exactly this reason, and a test enforces it.

Two transformations happen here, and both are mechanical so they cannot drift from
the real code:

* **Flattening.** The browser filesystem gets `store.py` and `validate.py` at the
  top level, and the one import line in the replay layer is rewritten to match.
  The module bodies are copied verbatim -- the page runs the same validator the
  library does, not a reimplementation of it.
* **Trimming.** `abstract` and `abstract_sections` duplicate `quotable_text`,
  which is the only one the checking reads. Dropping them halves the payload.
"""

from __future__ import annotations

import html
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RUNS = ROOT / "demo" / "runs"
OUT = ROOT / "docs" / "index.html"

# Fields the page never reads. quotable_text is kept: it is what a quoted claim is
# checked against.
DROP_DETAIL = ("abstract", "abstract_sections")

MODULES = {
    "store.py": ROOT / "src/kegg_string_mcp/agent/store.py",
    "validate.py": ROOT / "src/kegg_string_mcp/agent/validate.py",
}


def trimmed_run(path: Path) -> dict:
    record = json.loads(path.read_text(encoding="utf-8"))
    for call in record["calls"]:
        for item in call.get("result", {}).get("records", []):
            for key in DROP_DETAIL:
                item.get("detail", {}).pop(key, None)
    return record


def flatten(source: str) -> str:
    """Rewrite package imports for the browser's flat filesystem.

    Only import lines change; the module bodies are copied verbatim, so the page
    runs the same validator and the same layout the library does rather than a
    reimplementation of either.
    """
    for old, new in (
        ("from kegg_string_mcp.agent.store import RunStore", "from store import RunStore"),
        ("from kegg_string_mcp.agent.validate import ValidationReport, validate",
         "from validate import ValidationReport, validate"),
        ("from app.replay import (", "from replay import ("),
        ('RUNS_DIR = Path(__file__).resolve().parent.parent / "demo" / "runs"',
         'RUNS_DIR = Path("runs")'),
    ):
        source = source.replace(old, new)
    return source


def gradio_file(name: str, body: str) -> str:
    return f'<gradio-file name="{name}">\n{html.escape(body)}\n</gradio-file>'


def build() -> Path:
    parts = [gradio_file(name, path.read_text(encoding="utf-8"))
             for name, path in MODULES.items()]
    parts.append(gradio_file("replay.py",
                             flatten((ROOT / "app" / "replay.py").read_text(encoding="utf-8"))))
    parts.append(gradio_file("ui.py",
                             flatten((ROOT / "app" / "ui.py").read_text(encoding="utf-8"))))
    for run in sorted(RUNS.glob("*.json")):
        parts.append(gradio_file(f"runs/{run.name}",
                                 json.dumps(trimmed_run(run), separators=(",", ":"))))
    # Gradio-Lite runs app.py and looks for `demo`.
    parts.append(gradio_file("app.py", "from ui import build\n\ndemo = build()\n"))

    page = TEMPLATE.replace("__FILES__", "\n".join(parts))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(page, encoding="utf-8")
    return OUT


TEMPLATE = """<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8" />
<meta name="viewport" content="width=device-width, initial-scale=1" />
<title>Gene annotation with checked citations</title>
<meta name="description"
      content="Watch a language model annotate a tuberculosis gene, then watch every
               citation it wrote get checked against what the tools actually returned." />
<script type="module" crossorigin
        src="https://cdn.jsdelivr.net/npm/@gradio/lite/dist/lite.js"></script>
<link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/@gradio/lite/dist/lite.css" />
<style>
  body { margin: 0; font-family: system-ui, sans-serif; }
  .boot { padding: 2rem; max-width: 46rem; margin: 0 auto; color: #444; }
</style>
</head>
<body>
<gradio-lite>
__FILES__
</gradio-lite>
<noscript>
  <div class="boot">
    <h1>Gene annotation with checked citations</h1>
    <p>This page runs Python in your browser, so it needs JavaScript enabled.
       The code and the same runs are on
       <a href="https://github.com/YaaOppong/kegg-string-mcp">GitHub</a>.</p>
  </div>
</noscript>
</body>
</html>
"""


if __name__ == "__main__":
    written = build()
    print(f"  wrote {written.relative_to(ROOT)}  ({written.stat().st_size // 1024} KB)")
