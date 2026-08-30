"""Generate the static, serverless build of the demo for GitHub Pages.

Runs on **bare Pyodide with no packages installed**. An earlier attempt used
Gradio-Lite, which installs gradio through micropip and so pulls in
huggingface-hub -- a package with no pure-Python wheel, which fails at load with
`Can't find a pure Python 3 wheel`. Since the replay layer is standard-library
only, dropping the framework removes dependency resolution from the page
entirely: there is nothing left to fail.

The UI is therefore hand-written HTML and JavaScript. That trade is deliberate and
narrow: the browser formats, but every verdict comes from `browser_api` calling
the same `validate()` the library uses. Nothing about the checking is
reimplemented in JavaScript.

Two mechanical transformations, so the page cannot drift from the code:

* **Flattening** -- `store.py` and `validate.py` go in at the top level and the
  import lines in `replay.py` are rewritten. Module bodies are copied verbatim.
* **Trimming** -- `abstract` and `abstract_sections` duplicate `quotable_text`,
  the only one the checking reads. Dropping them halves the payload.
"""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RUNS = ROOT / "demo" / "runs"
OUT = ROOT / "docs" / "index.html"
PYODIDE = "https://cdn.jsdelivr.net/pyodide/v0.27.3/full/"

DROP_DETAIL = ("abstract", "abstract_sections")

MODULES = {
    "store.py": ROOT / "src/kegg_string_mcp/agent/store.py",
    "validate.py": ROOT / "src/kegg_string_mcp/agent/validate.py",
    "replay.py": ROOT / "app/replay.py",
    "browser_api.py": ROOT / "app/browser_api.py",
}

REWRITES = (
    ("from kegg_string_mcp.agent.store import RunStore", "from store import RunStore"),
    ("from kegg_string_mcp.agent.validate import ValidationReport, validate",
     "from validate import ValidationReport, validate"),
    # Absolute, so nothing depends on the interpreter's working directory.
    ('RUNS_DIR = Path(__file__).resolve().parent.parent / "demo" / "runs"',
     'RUNS_DIR = Path("/demo/runs")'),
)


def flatten(source: str) -> str:
    for old, new in REWRITES:
        source = source.replace(old, new)
    return source


def trimmed_run(path: Path) -> dict:
    record = json.loads(path.read_text(encoding="utf-8"))
    for call in record["calls"]:
        for item in call.get("result", {}).get("records", []):
            for key in DROP_DETAIL:
                item.get("detail", {}).pop(key, None)
    return record


def payload() -> dict[str, str]:
    files = {name: flatten(path.read_text(encoding="utf-8"))
             for name, path in MODULES.items()}
    for run in sorted(RUNS.glob("*.json")):
        files[f"runs/{run.name}"] = json.dumps(trimmed_run(run), separators=(",", ":"))
    return files


def build() -> Path:
    blob = json.dumps(payload()).replace("</", "<\\/")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(TEMPLATE.replace("__PAYLOAD__", blob).replace("__PYODIDE__", PYODIDE),
                   encoding="utf-8")
    return OUT


TEMPLATE = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8" />
<meta name="viewport" content="width=device-width, initial-scale=1" />
<title>Gene annotation with checked citations</title>
<meta name="description" content="Watch a language model annotate a tuberculosis gene, then watch every citation it wrote get checked against what the tools returned." />
<style>
 :root { --fail:#b3261e; --fail-bg:#fdecea; --ok:#1e8e3e; --ok-bg:#e9f7ef; --line:#d8dce2; }
 * { box-sizing: border-box; }
 body { margin:0; font:16px/1.6 system-ui,-apple-system,Segoe UI,Roboto,sans-serif; color:#1b1f24; }
 main { max-width: 60rem; margin: 0 auto; padding: 2rem 1.25rem 5rem; }
 h1 { font-size: 1.9rem; line-height:1.25; margin:.2em 0 .6em; }
 h2 { font-size: 1.15rem; margin-top: 2.4rem; border-top:1px solid var(--line); padding-top:1.4rem; }
 a { color:#0b57d0; }
 .lede { font-size:1.03rem; }
 .banner { border-left:5px solid var(--ok); background:var(--ok-bg); padding:1rem 1.25rem;
           border-radius:6px; margin:1.4rem 0; }
 .banner.fail { border-left-color:var(--fail); background:var(--fail-bg); }
 select { font:inherit; padding:.55rem .7rem; width:100%; max-width:44rem;
          border:1px solid var(--line); border-radius:6px; background:#fff; }
 table { border-collapse:collapse; width:100%; font-size:.86rem; margin-top:.8rem; display:block;
         overflow-x:auto; }
 th,td { border:1px solid var(--line); padding:.45rem .6rem; text-align:left; vertical-align:top; }
 th { background:#f4f6f8; font-weight:600; }
 td.bad { color:var(--fail); font-weight:600; white-space:nowrap; }
 td.good { color:var(--ok); white-space:nowrap; }
 code, .mono { font-family: ui-monospace,SFMono-Regular,Menlo,monospace; font-size:.86em; }
 #summary { background:#fbfcfd; border:1px solid var(--line); border-radius:6px; padding:1rem 1.25rem; }
 #boot { color:#555; }
 .note { color:#555; font-size:.92rem; }
</style>
</head>
<body>
<main>
<h1>Does the annotation say what its sources say?</h1>

<p class="lede">This annotates <em>Mycobacterium tuberculosis</em> genes. It hands a
language model a set of lookup tools &mdash; KEGG for pathways, STRING for protein
interactions, PubMed for literature &mdash; lets it decide which to call and when it has
enough, then writes up what it found.</p>

<p class="lede">Then the last step. <strong>Every identifier in the write-up is checked
against what the tools returned.</strong> Not by another
model judging it, but by looking: was this exact record retrieved, for this exact gene,
and where text is quoted, do those words really appear in the source? A model can produce
a fluent, plausible, correctly-formatted citation for something it was never shown. This
catches that.</p>

<p class="note">Nothing here is live &mdash; these are real runs captured earlier and
replayed. The checking itself runs fresh in your browser each time, using the same code
as the library. <a href="https://github.com/YaaOppong/kegg-string-mcp">Source code, and
the reasoning behind every design decision</a>.</p>

<p id="boot">Starting Python in your browser&hellip; (a few seconds, first time only)</p>

<div id="app" hidden>
  <label for="gene"><strong>Gene</strong></label>
  <select id="gene"></select>
  <div id="verdict" class="banner"></div>

  <h2>1. What the model asked for</h2>
  <p class="note">The model chooses which tools to call and when it has enough. Calls
  marked <em>pipeline</em> were made deterministically by the code, not chosen by the
  model.</p>
  <p id="turns" class="mono"></p>
  <table id="calls"></table>

  <h2>2. What it wrote</h2>
  <div id="summary"></div>

  <h2>3. What survived checking</h2>
  <table id="citations"></table>
  <div id="quotes-wrap" hidden>
    <h2>Quoted passages</h2>
    <table id="quotes"></table>
  </div>

  <p class="note" style="margin-top:2rem"><strong>CROSS-TARGET</strong> means the record
  was retrieved during this run, but for a different gene than the sentence attributes it
  to &mdash; a real identifier in the wrong place, which a check of &ldquo;does this ID
  exist?&rdquo; would miss. <strong>UNSUPPORTED</strong> means no tool returned it at all.
  <br/><br/>Research use only. Not for clinical decisions.
  <a href="https://github.com/YaaOppong/kegg-string-mcp">Repository</a>.</p>
</div>
</main>

<script type="application/json" id="payload">__PAYLOAD__</script>
<script src="__PYODIDE__pyodide.js"></script>
<script>
const FILES = JSON.parse(document.getElementById("payload").textContent);
const esc = s => String(s ?? "").replace(/[&<>"]/g, c =>
  ({ "&":"&amp;", "<":"&lt;", ">":"&gt;", '"':"&quot;" }[c]));

// Minimal Markdown: headings, bold, italics, inline code, paragraphs. Deliberately
// not a library -- another CDN dependency is another thing that can fail to load,
// which is the mistake this rewrite exists to undo.
function markdown(src) {
  return esc(src).split(/\n{2,}/).map(block => {
    const h = block.match(/^(#{1,6})\s+(.*)$/s);
    let body = (h ? h[2] : block)
      .replace(/\*\*(.+?)\*\*/gs, "<strong>$1</strong>")
      .replace(/(^|\W)\*([^*\n]+)\*/gs, "$1<em>$2</em>")
      .replace(/`([^`]+)`/g, "<code>$1</code>")
      .replace(/\n/g, "<br/>");
    if (h) return `<h${Math.min(h[1].length + 2, 6)}>${body}</h${Math.min(h[1].length + 2, 6)}>`;
    return `<p>${body}</p>`;
  }).join("\n");
}

function table(el, headers, rows, classify) {
  el.innerHTML =
    "<thead><tr>" + headers.map(h => `<th>${esc(h)}</th>`).join("") + "</tr></thead><tbody>" +
    rows.map(r => "<tr>" + r.map((cell, i) => {
      const cls = classify ? classify(i, cell) : "";
      return `<td${cls ? ` class="${cls}"` : ""}>${esc(cell)}</td>`;
    }).join("") + "</tr>").join("") + "</tbody>";
}

const FAIL_BANNER = `<strong>&#9888;&#65039; On this run, the checking caught a bad citation.</strong>
<br/>The write-up below cites records that were <strong>not</strong>
retrieved for this gene. The model produced real, correctly-formatted identifiers for
things it was never shown for this gene &mdash; that is the failure this project exists
to catch, and it is caught here on a real run. Details in section 3.`;

const CLEAN_BANNER = `<strong>&#9989; On this run, every citation checked out.</strong>
<br/>Each identifier in the write-up names a record a tool returned for this
gene, and every quoted span appears verbatim in its source. Pick one of the genes marked
&#9888;&#65039; to watch the checking catch something.`;

async function main() {
  const step = label => { window.__step = label; };

  step("loading Pyodide");
  const pyodide = await loadPyodide({ indexURL: "__PYODIDE__" });

  // An absolute directory, not the cwd: Pyodide's working directory is not on
  // sys.path, and a relative write lands somewhere the import machinery will not
  // look. Bytes, not a JS string: Emscripten's FS.writeFile wants a Uint8Array,
  // and handing it a string is silently wrong or throws depending on the build.
  step("writing files into the virtual filesystem");
  const encoder = new TextEncoder();
  pyodide.FS.mkdirTree("/demo/runs");
  for (const [name, body] of Object.entries(FILES)) {
    pyodide.FS.writeFile("/demo/" + name, encoder.encode(body));
  }

  step("importing the checker");
  // Nothing is installed here: every module in the payload is standard library.
  // runPython only, deliberately: it is the most basic entry point in the JS API,
  // so the page depends on as little of Pyodide's surface as possible.
  pyodide.runPython("import sys; sys.path.insert(0, '/demo')\nimport browser_api");
  const api = {
    index: () => pyodide.runPython("import browser_api; browser_api.index()"),
    run: name => pyodide.runPython(
      `import browser_api; browser_api.run(${JSON.stringify(name)})`),
  };

  step("building the gene list");
  const picker = document.getElementById("gene");
  for (const entry of JSON.parse(api.index())) {
    const option = document.createElement("option");
    option.value = entry.id;
    // Mark failing runs in the picker itself: if someone has to select the right
    // gene to find the point, most will not.
    option.textContent = (entry.fails ? "\u26A0\uFE0F  " : "") + entry.label +
                         (entry.fails ? "  \u2014 citation check FAILS" : "");
    picker.append(option);
  }

  function show(name) {
    const d = JSON.parse(api.run(name));
    const verdict = document.getElementById("verdict");
    verdict.className = "banner" + (d.clean ? "" : " fail");
    verdict.innerHTML = d.clean ? CLEAN_BANNER : FAIL_BANNER;

    document.getElementById("turns").textContent =
      `The model made ${d.model_calls} tool call(s) across ${d.turns} turn(s)` +
      (d.pipeline_calls ? `, after the pipeline fetched ${d.pipeline_calls} deterministically` : "");

    table(document.getElementById("calls"),
          ["requested by", "tool", "arguments", "returned", "note from the tool"], d.calls);
    document.getElementById("summary").innerHTML = markdown(d.summary);
    table(document.getElementById("citations"),
          ["identifier", "kind", "status", "detail"], d.citations,
          (i, cell) => i === 2 ? (cell === "verified" ? "good" : "bad") : "");

    const wrap = document.getElementById("quotes-wrap");
    wrap.hidden = d.quotes.length === 0;
    if (d.quotes.length) {
      table(document.getElementById("quotes"),
            ["record", "status", "quoted in the write-up", "closest text in the source"],
            d.quotes, (i, cell) => i === 1 ? (cell === "verified" ? "good" : "bad") : "");
    }
  }

  step("rendering the first run");
  picker.addEventListener("change", () => show(picker.value));
  show(picker.value);
  document.getElementById("boot").hidden = true;
  document.getElementById("app").hidden = false;
}

main().catch(err => {
  // esc(err) on an Error prints "[object Object]", which says nothing. Report the
  // step that failed and the actual message, so a bug report is usable.
  const detail = (err && (err.message || err.toString())) || String(err);
  const where = window.__step ? `while ${window.__step}` : "during startup";
  console.error(err);
  document.getElementById("boot").innerHTML =
    `<strong>The demo failed to start</strong> ${esc(where)}.` +
    `<pre style="white-space:pre-wrap;background:#fdecea;padding:1rem;border-radius:6px">` +
    `${esc(detail)}</pre>` +
    `<p>Please <a href="https://github.com/YaaOppong/kegg-string-mcp/issues">open an issue</a> ` +
    `with this message. The same runs are in ` +
    `<a href="https://github.com/YaaOppong/kegg-string-mcp">the repository</a>.</p>`;
});
</script>
</body>
</html>
"""


if __name__ == "__main__":
    written = build()
    print(f"  wrote {written.relative_to(ROOT)}  ({written.stat().st_size // 1024} KB)")
