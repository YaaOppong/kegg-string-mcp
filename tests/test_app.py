"""The demo replays committed runs and re-validates them. No gradio needed here:
the replay layer is deliberately separate from the UI so it can be tested without
installing a web framework, and so the UI cannot quietly change what is checked.
"""


from pathlib import Path

import pytest

from app.replay import RUNS_DIR, available, citation_rows, load, quote_rows, tool_call_rows

RUNS = available()


def test_demo_runs_are_committed():
    assert RUNS_DIR.is_dir()
    assert len(RUNS) >= 6, "the demo should show a range of genes, not one"


@pytest.mark.parametrize("name", RUNS)
def test_every_run_replays(name):
    replay = load(name)
    assert replay.summary.strip(), f"{name} has no write-up to show"
    assert replay.calls, f"{name} made no tool calls"
    assert replay.turns, f"{name} has no decision log"


@pytest.mark.parametrize("name", RUNS)
def test_every_tool_call_is_shown_with_its_origin(name):
    """The multi-turn loop must be visible. Calls the pipeline makes itself are
    labelled as such rather than being attributed to a turn the model never took."""
    replay = load(name)
    rows = tool_call_rows(replay)
    assert len(rows) == len(replay.calls)
    assert all(row[0] == "pipeline" or row[0].startswith("turn ") for row in rows)


def test_model_turns_are_visible_in_single_gene_mode():
    rows = tool_call_rows(load("furA"))
    assert {row[0] for row in rows} >= {"turn 1", "turn 2"}, "the loop should be legible"


def test_pipeline_precomputation_is_distinguished_from_model_choices():
    """Epistasis mode pre-fetches deterministically before the model runs."""
    rows = tool_call_rows(load("katG-ahpC-epistasis"))
    assert any(row[0] == "pipeline" for row in rows)
    assert any(row[0].startswith("turn") for row in rows)


def test_at_least_one_run_shows_a_caught_failure():
    """Non-negotiable. A demo where everything passes proves nothing -- the point of
    the project is that the checking catches things, so a visitor must be able to
    see it catch something real."""
    failing = [name for name in RUNS if not load(name).clean]
    assert failing, "no demo run exhibits a validation failure"


def test_the_default_run_is_one_that_fails():
    """A visitor who changes nothing should still see the point. Asserted against
    the gradio-free module so the whole matrix checks it, not just the demo job."""
    from app.replay import ORDERED

    assert not load(ORDERED[0]).clean, "the first gene offered should be a caught failure"


def test_every_offered_run_exists_and_is_labelled():
    from app.replay import LABELS, ORDERED

    assert set(ORDERED) <= set(RUNS), "the picker offers a run that is not committed"
    assert set(ORDERED) <= set(LABELS), "an offered run has no human-readable label"


def test_failures_name_their_class_and_reason():
    for name in RUNS:
        replay = load(name)
        for row in citation_rows(replay):
            if row[2] != "verified":
                assert row[2] in {"UNSUPPORTED", "CROSS-TARGET"}
                assert row[3], "a failure must say why"


def test_no_run_shows_a_stale_verdict():
    """Verdicts are recomputed, not stored, so the demo cannot show conclusions the
    current validator no longer reaches."""
    import json

    raw = json.loads((RUNS_DIR / f"{RUNS[0]}.json").read_text())
    assert "validation" not in raw, "a stored verdict would drift from the library"


@pytest.mark.parametrize("name", RUNS)
def test_quote_rows_render_for_runs_with_quotes(name):
    for row in quote_rows(load(name)):
        assert row[1] in {"verified", "NOT IN SOURCE"}


def test_app_module_builds_without_a_network_or_key(monkeypatch):
    """Importing and building the UI must not need credentials or a live service."""
    pytest.importorskip("gradio")
    monkeypatch.delenv("ANTHROPIC_API_KEY", raising=False)
    from app.app import build

    assert build() is not None


def test_the_replay_layer_needs_only_the_standard_library():
    """The demo must be able to run where anthropic, httpx and mcp are awkward --
    a browser under Pyodide, most obviously. store.py and validate.py are stdlib
    only; an eager re-export in agent/__init__.py used to drag the whole HTTP and
    model stack in behind them."""
    import json
    import subprocess
    import sys

    probe = (
        "import sys, json\n"
        "before = set(sys.modules)\n"
        "from app.replay import load, available\n"
        "heavy = {'anthropic', 'mcp', 'httpx', 'pydantic', 'pydantic_core', 'gradio'}\n"
        "print(json.dumps(sorted({m.split('.')[0] for m in set(sys.modules) - before} & heavy)))"
    )
    result = subprocess.run([sys.executable, "-c", probe], capture_output=True, text=True,
                            cwd=Path(__file__).resolve().parent.parent, check=False)
    assert result.returncode == 0, result.stderr
    pulled = json.loads(result.stdout.strip().splitlines()[-1])
    assert pulled == [], f"replay pulled in {pulled}"


def test_lazy_agent_exports_still_resolve():
    from kegg_string_mcp import agent

    assert agent.validate is not None
    assert agent.RunStore is not None
    assert "annotate_gene" in dir(agent)


# --- the failure must be impossible to miss ---------------------------------

def test_failing_runs_are_marked_in_the_picker():
    """If someone has to select the right gene to find the point, most will not."""
    pytest.importorskip("gradio")
    from app.replay import ORDERED
    from app.ui import _label

    for name in ORDERED:
        marked = "⚠️" in _label(name)
        assert marked is (not load(name).clean), f"{name} is mislabelled in the picker"


def test_the_verdict_is_the_first_thing_rendered():
    """The verdict used to sit below the tool-call table, so a visitor had to scroll
    to find out whether anything was caught."""
    pytest.importorskip("gradio")
    from app.ui import render

    verdict = render("furA")[0]
    assert "caught a bad citation" in verdict


def test_a_clean_run_points_at_the_failing_ones():
    pytest.importorskip("gradio")
    from app.ui import render

    verdict = render("katG")[0]
    assert "every citation checked out" in verdict.lower()
    assert "⚠️" in verdict, "a clean run should say where to find a caught failure"


def test_the_static_build_rewrites_only_import_lines():
    """A divergent copy of the validator would let the page show a verdict the
    library does not produce, so the module bodies must be copied verbatim."""
    from demo.build_pages import REWRITES, flatten

    source = (Path(__file__).resolve().parent.parent / "app" / "replay.py").read_text()
    flat = flatten(source)
    assert "from store import RunStore" in flat
    assert "from kegg_string_mcp" not in flat

    restored = flat
    for old, new in REWRITES:
        restored = restored.replace(new, old)
    assert restored == source, "the build changed something other than an import line"


def test_the_static_build_installs_no_packages():
    """Gradio-Lite pulled gradio through micropip, which pulls huggingface-hub --
    a package with no pure-Python wheel, so the page failed to start. Bare Pyodide
    with a standard-library-only payload has no dependency resolution to fail."""
    import ast
    import sys

    from demo.build_pages import payload

    files = payload()
    assert "browser_api.py" in files and "validate.py" in files

    # Parse rather than grep: a docstring explaining why micropip is avoided is
    # not an install, and a substring check cannot tell the difference.
    stdlib = set(sys.stdlib_module_names)
    local = {name[:-3] for name in files if name.endswith(".py")}
    for name, body in files.items():
        if not name.endswith(".py"):
            continue
        for node in ast.walk(ast.parse(body)):
            if isinstance(node, ast.Import):
                roots = [alias.name.split(".")[0] for alias in node.names]
            elif isinstance(node, ast.ImportFrom):
                roots = [(node.module or "").split(".")[0]] if node.level == 0 else []
            else:
                continue
            for root in roots:
                assert root in stdlib or root in local, (
                    f"{name} imports {root!r}, which the browser would have to install")


def test_the_static_payload_carries_every_offered_run():
    from app.replay import ORDERED
    from demo.build_pages import payload

    files = payload()
    for name in ORDERED:
        assert f"runs/{name}.json" in files, f"{name} is offered but not embedded"
