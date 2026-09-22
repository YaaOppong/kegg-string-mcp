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

# The rule classification goes in as a package tree rather than flattened. Its
# modules import each other by their real names, and both `__init__.py` files are
# import-safe -- the fetching modules are not in the tree, so nothing reaches for
# httpx. That means no import rewriting at all: every file is the file.
RULE_PACKAGE = (
    "__init__.py", "rules/__init__.py", "rules/parse.py", "rules/annotation.py",
    "rules/features.py", "rules/catalogue.py", "rules/enrichment.py", "rules/sets.py",
    "rules/nesting.py", "rules/signature.py", "rules/evidence.py", "rules/questions.py",
)
RULE_FIXTURE = ROOT / "demo" / "rule_fixture.json"

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

    # The rule tab, when a fixture has been built. Absent, the page still works
    # and the tab is hidden rather than showing an error.
    if RULE_FIXTURE.exists():
        for name in RULE_PACKAGE:
            source = ROOT / "src" / "kegg_string_mcp" / name
            files[f"kegg_string_mcp/{name}"] = source.read_text(encoding="utf-8")
        files["rule_replay.py"] = (ROOT / "app" / "rule_replay.py").read_text(encoding="utf-8")
        files["rule_fixture.json"] = RULE_FIXTURE.read_text(encoding="utf-8")
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
<title>LLM Gene Annotation Demo</title>
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
 nav#tabs { margin:1.5rem 0 .5rem; border-bottom:1px solid var(--line); }
 .tab { font:inherit; background:none; border:0; border-bottom:3px solid transparent;
        padding:.6rem .9rem; cursor:pointer; color:#4a5158; }
 .tab.on { border-bottom-color:#0b57d0; color:#1b1f24; font-weight:600; }
 tr.pick { cursor:pointer; }
 tr.pick:hover td { background:#f3f6fb; }
 tr.chosen td { background:#e8f0fe; }
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
<h1>LLM Gene Annotation Demo</h1>
<p class="lede" style="font-size:1.12rem;color:#444;margin-top:-.4rem">
Does the annotation say what its sources say?</p>

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

<nav id="tabs" hidden>
  <button id="tab-genes" class="tab on" type="button">Gene annotation</button>
  <button id="tab-rules" class="tab" type="button" hidden>Rule classification</button>
</nav>

<div id="app" hidden>
  <label for="gene"><strong>Gene</strong></label>
  <select id="gene"></select>
  <div id="verdict" class="banner"></div>

  <h2>1. What the model asked for</h2>
  <p class="note">The model chooses which tools to call and when it has enough.</p>
  <p class="note" id="pipeline-note" hidden>This is a two-gene run, so the calls marked
  <em>pipeline</em> were made deterministically by the code before the model started, not
  chosen by it. Single-gene runs have none.</p>
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
  was retrieved during this run, but not for the gene being annotated &mdash; a real
  identifier that a check of &ldquo;does this ID exist?&rdquo; would wave through. It is
  a mechanical comparison against the run&rsquo;s own record of which gene returned what,
  not a reading of the sentence, so a neighbour gene discussed as context is flagged as
  well. <strong>UNSUPPORTED</strong> means no tool returned it at all. Two-gene runs have
  no single target and so are not checked this way.
  <br/><br/>Research use only. Not for clinical decisions.
  <a href="https://github.com/YaaOppong/kegg-string-mcp">Repository</a>.</p>
</div>

<div id="rules-app" hidden>
  <p class="lede">An upstream classifier emits rules over per-locus variant states.
  This page reads each one against the WHO catalogue, the annotation and the
  interaction data, and says what is already accounted for.</p>

  <div class="banner" id="rules-banner"></div>

  <p class="note"><strong>Nothing here was stored.</strong> The gene tab above replays a
  model&rsquo;s summary, because a summary cannot be regenerated without an API key. A
  rule classification has no model in it &mdash; given the catalogue rows, the annotation
  and the rules, every verdict is arithmetic. So the page ships those inputs and runs the
  same classifier the library does. What you see is computed here, now, and cannot drift
  from the code.</p>

  <h2>Rules</h2>
  <p class="note">Pick one to see the reasoning. <em>minimal</em> means no smaller rule in
  the set is contained in it.</p>
  <table id="rule-list"></table>

  <div id="rule-detail" hidden>
    <h2>Why</h2>
    <div id="rule-verdict"></div>
    <table id="rule-roles"></table>
  </div>

  <p class="note" style="margin-top:2rem">A locus being catalogued means the
  <em>gene</em> is known to resistance, never that the variant driving the rule is one of
  the graded ones &mdash; katG holds 1,771 catalogued variants of which 139 are graded
  associated. Enrichment is measured against the full annotation, not the subset shipped
  here, and a p measures surprise against a uniform draw rather than evidence that the
  loci are related.
  <br/><br/>Research use only. Not for clinical decisions.</p>
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

const FAIL_BANNER = `<strong>&#9888;&#65039; On this run, the checking flagged a citation.</strong>
<br/>The write-up below cites records that were <strong>not</strong>
retrieved for this gene. The check compares each identifier against the gene being
annotated; it does not read the sentence, so a neighbour gene cited as context is
flagged too. Read section 3 against the write-up and judge it yourself &mdash; that is
the point of showing the working rather than a verdict.`;

const CLEAN_BANNER = `<strong>&#9989; On this run, every citation checked out.</strong>
<br/>Each identifier in the write-up names a record a tool returned for this
gene, and every quoted span appears verbatim in its source. Pick one of the genes marked
&#9888;&#65039; to watch the checking fire.`;

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
  pyodide.FS.mkdirTree("/demo/kegg_string_mcp/rules");
  for (const [name, body] of Object.entries(FILES)) {
    pyodide.FS.writeFile("/demo/" + name, encoder.encode(body));
  }

  step("importing the checker");
  // Nothing is installed here: every module in the payload is standard library.
  // runPython only, deliberately: it is the most basic entry point in the JS API,
  // so the page depends on as little of Pyodide's surface as possible.
  pyodide.runPython("import sys; sys.path.insert(0, '/demo')\nimport browser_api");
  // The rule tab computes rather than replays: same classifier, same inputs, no
  // stored verdict anywhere on the page.
  const hasRules = "rule_fixture.json" in FILES;
  if (hasRules) {
    pyodide.runPython(
      "import json, rule_replay\n"
      + "_RULES = json.dumps(rule_replay.classified("
      + "rule_replay.load('/demo/rule_fixture.json'), workdir='/demo'))");
  }
  const api = {
    index: () => pyodide.runPython("import browser_api; browser_api.index()"),
    rules: () => hasRules ? JSON.parse(pyodide.runPython("_RULES")) : null,
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
    // Shown only where there is something to see: the note used to sit above every
    // run, including the nine with no pipeline row in the table below it.
    document.getElementById("pipeline-note").hidden = !d.pipeline_calls;

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

  step("classifying the rules");
  const rules = api.rules();
  if (rules) renderRules(rules);

  document.getElementById("boot").hidden = true;
  document.getElementById("tabs").hidden = false;
  document.getElementById("app").hidden = false;
}

// --- the rule tab ---------------------------------------------------------
// Rendering only. Every verdict, role and q value on this tab was computed by
// the library a moment ago; nothing here recomputes or reinterprets any of it.

function renderRules(data) {
  document.getElementById("tab-rules").hidden = false;
  const counts = Object.entries(data.by_signature)
    .map(([k, n]) => `${esc(k)} &times;${n}`).join(" &middot; ");
  document.getElementById("rules-banner").innerHTML =
    `<strong>${data.rules.length} rules</strong> over ${data.loci} loci, `
    + `${data.nesting.minimal} minimal. ${counts}`;

  const list = document.getElementById("rule-list");
  list.innerHTML = "<tr><th>Rule</th><th>Predicts</th><th>Classification</th>"
    + "<th>Minimal</th></tr>"
    + data.rules.map((r, i) =>
        `<tr class="pick" data-i="${i}"><td class="mono">${esc(r.conditions)}</td>`
        + `<td>${esc(r.predicted_class)}</td><td>${esc(r.primary_signature)}</td>`
        + `<td>${r.is_minimal === "True" ? "yes" : ""}</td></tr>`).join("");

  const detail = document.getElementById("rule-detail");
  list.querySelectorAll("tr.pick").forEach(row => row.addEventListener("click", () => {
    list.querySelectorAll("tr").forEach(r => r.classList.remove("chosen"));
    row.classList.add("chosen");
    const r = data.rules[Number(row.dataset.i)];
    document.getElementById("rule-verdict").innerHTML = `<p>${esc(r.verdict)}</p>`;
    const rows = [["role of each condition", r.roles],
                  ["signatures", r.signatures],
                  ["enriched category", r.enriched_term === "NA" ? "NA"
                     : `${r.enriched_term} (${r.enriched_m_of_k}, expected `
                       + `${r.enriched_expected}, q=${r.enriched_q}, `
                       + `${r.terms_tested} term(s) tested)`],
                  ["drugs spanned", r.drugs],
                  ["links", r.links],
                  ["contained in", r.contained_in],
                  ["contradicted by", r.contradicted_by]];
    document.getElementById("rule-roles").innerHTML =
      rows.filter(([, v]) => v && v !== "NA")
          .map(([k, v]) => `<tr><th>${esc(k)}</th><td class="mono">${esc(v)}</td></tr>`)
          .join("");
    detail.hidden = false;
  }));
}

function showTab(which) {
  for (const [name, pane] of [["genes", "app"], ["rules", "rules-app"]]) {
    document.getElementById(pane).hidden = name !== which;
    document.getElementById("tab-" + name).classList.toggle("on", name === which);
  }
}
document.addEventListener("click", event => {
  if (event.target.id === "tab-genes") showTab("genes");
  if (event.target.id === "tab-rules") showTab("rules");
});

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
