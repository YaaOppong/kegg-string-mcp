"""The demo replays committed runs and re-validates them. No gradio needed here:
the replay layer is deliberately separate from the UI so it can be tested without
installing a web framework, and so the UI cannot quietly change what is checked.
"""


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
