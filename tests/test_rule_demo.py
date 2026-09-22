"""The browser path: recompute a classification from a captured fixture.

The gene demo replays a model's summary and re-validates it, because the summary
cannot be regenerated without an API key. A rule classification has no model in
it, so the page ships the inputs and runs the real classifier -- which means a
verdict it shows is one the library produces rather than one that was stored.

These tests pin the two properties that make that trustworthy: the page's
arithmetic uses the real background, and nothing it needs has to be installed.
"""

import sys
from pathlib import Path

import pytest

FIXTURE = Path(__file__).resolve().parents[1] / "demo" / "rule_fixture.json"

pytestmark = pytest.mark.skipif(not FIXTURE.exists(),
                                reason="demo/rule_fixture.json has not been built")


@pytest.fixture
def replayed(tmp_path):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from app.rule_replay import classified, load

    return classified(load(FIXTURE), workdir=tmp_path)


def test_every_rule_comes_back_classified(replayed):
    assert replayed["rules"]
    for row in replayed["rules"]:
        assert row["primary_signature"]
        assert row["verdict"]
        assert row["roles"] != "NA"


def test_the_enrichment_uses_the_real_background_not_the_shipped_subset(replayed):
    """The page holds a few dozen loci; the pipeline's universe is thousands. A
    universe rebuilt from the subset would answer a different question under the
    same name, which is worse than showing nothing."""
    import json

    fixture = json.loads(FIXTURE.read_text())
    captured = fixture["universe"]
    shipped = fixture["annotation_gff"].count("\n")

    assert captured["total"] > shipped * 10, "the background should dwarf what is shipped"
    for row in replayed["rules"]:
        if row["enriched_term"] != "NA":
            # n is the universe the q was computed against.
            assert int(row["enriched_m_of_k"].split("/")[1]) <= shipped
            assert float(row["enriched_q"]) <= 1.0
    assert replayed["terms"] == len(captured["sizes"])


def test_the_fixture_carries_interaction_evidence(replayed):
    """Without it a compensation verdict claims no link was found, which is a
    different and weaker statement than the pipeline's -- so the page would be
    misleading rather than merely incomplete."""
    import json

    fixture = json.loads(FIXTURE.read_text())
    assert fixture["links"], "no links captured"
    assert fixture["partners"], "no partner lists captured"

    compensation = [r for r in replayed["rules"]
                    if r["primary_signature"] == "known:compensation"]
    if compensation:
        assert any("beyond literature co-mention" in r["verdict"]
                   or "textmining channel" in r["verdict"] for r in compensation)


def test_the_browser_path_runs_with_httpx_and_pydantic_blocked(tmp_path):
    """Bare Pyodide installs nothing, so anything the page imports must be in the
    standard library."""
    import importlib

    blocked = ("httpx", "pydantic", "pydantic_core")

    class _Block:
        def find_spec(self, name, path=None, target=None):
            if name.split(".")[0] in blocked:
                raise ImportError(f"{name} is blocked by this test")

    saved = {name: sys.modules.pop(name) for name in list(sys.modules)
             if name.startswith(("kegg_string_mcp", "app.rule_replay", *blocked))}
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    sys.meta_path.insert(0, _Block())
    try:
        module = importlib.import_module("app.rule_replay")
        out = module.classified(module.load(FIXTURE), workdir=tmp_path)
        assert out["rules"]
        for name in blocked:
            assert name not in sys.modules, f"{name} was imported anyway"
    finally:
        sys.meta_path.pop(0)
        for name in list(sys.modules):
            if name.startswith(("kegg_string_mcp", "app.rule_replay")):
                del sys.modules[name]
        sys.modules.update(saved)


def test_the_pages_payload_runs_the_classification(tmp_path, monkeypatch):
    """Simulate the page's bootstrap: write the payload files where Pyodide would
    and run the same Python the page runs.

    The page cannot be opened here, so this is the closest check there is -- and
    it covers the failure that would otherwise only appear in a browser: a module
    left out of the payload, or one that imports something Pyodide has not got.
    """
    import importlib

    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from demo.build_pages import payload

    files = payload()
    assert "rule_fixture.json" in files, "the rule tab was not included"

    for name, body in files.items():
        target = tmp_path / name
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(body, encoding="utf-8")

    blocked = ("httpx", "pydantic", "pydantic_core")

    class _Block:
        def find_spec(self, name, path=None, target=None):
            if name.split(".")[0] in blocked:
                raise ImportError(f"{name} is blocked by this test")

    saved = {name: sys.modules.pop(name) for name in list(sys.modules)
             if name.startswith(("kegg_string_mcp", "rule_replay", *blocked))}
    monkeypatch.syspath_prepend(str(tmp_path))
    sys.meta_path.insert(0, _Block())
    try:
        replay = importlib.import_module("rule_replay")
        out = replay.classified(replay.load(tmp_path / "rule_fixture.json"),
                                workdir=tmp_path)
        assert out["rules"] and out["by_signature"]
        assert all(row["verdict"] for row in out["rules"])
    finally:
        sys.meta_path.pop(0)
        for name in list(sys.modules):
            if name.startswith(("kegg_string_mcp", "rule_replay")):
                del sys.modules[name]
        sys.modules.update(saved)
