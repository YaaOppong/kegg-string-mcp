"""Store, evidence assembly and citation validation. No model, no network."""

from pathlib import Path

from kegg_string_mcp.agent.evidence import all_pairs, classify_pathway, pair_evidence
from kegg_string_mcp.agent.store import RunStore
from kegg_string_mcp.agent.validate import extract_citations, validate

KEGG_RESULT = {
    "resolved": {"kegg_gene_id": "mtu:Rv1908c", "matched_by": "symbol"},
    "records": [{"record_id": "mtu00360", "name": "Phenylalanine metabolism"},
                {"record_id": "mtu01100", "name": "Metabolic pathways"}],
    "record_ids": ["mtu00360", "mtu01100"],
}
STRING_RESULT = {
    "resolved": {"string_id": "83332.Rv1908c", "preferred_name": "katG"},
    "records": [{"record_id": "83332.Rv1909c", "name": "furA",
                 "detail": {"combined_score": 0.979, "textmining_score": 0.878,
                            "max_non_textmining_score": 0.829,
                            "evidence_beyond_textmining": True}}],
    "record_ids": ["83332.Rv1909c"],
}


def _store(tmp_path: Path) -> RunStore:
    store = RunStore(path=tmp_path / "run.jsonl", run_id="test")
    store.tool_result("kegg_pathways", {"gene": "katG"}, KEGG_RESULT)
    store.tool_result("string_partners", {"gene": "katG"}, STRING_RESULT)
    return store


# --- store -----------------------------------------------------------------


def test_resolved_identifiers_are_citable(tmp_path):
    """Citing the subject of your own query is correct, not a hallucination.
    A validator that flags it teaches the reader to ignore every flag."""
    store = _store(tmp_path)
    assert "mtu:Rv1908c" in store.citable_ids
    assert "83332.Rv1908c" in store.citable_ids


def test_records_are_citable_and_attributed_to_their_gene(tmp_path):
    store = _store(tmp_path)
    assert {"mtu00360", "mtu01100", "83332.Rv1909c"} <= store.citable_ids
    assert "mtu00360" in store.per_target["KATG"]


def test_store_is_append_only_and_replayable(tmp_path):
    store = _store(tmp_path)
    store.decision({"turn": 1, "stop_reason": "tool_use"})
    kinds = [entry["kind"] for entry in store.replay()]
    assert kinds == ["tool_result", "tool_result", "decision"]


# --- validation ------------------------------------------------------------


def test_unsupported_identifier_is_flagged(tmp_path):
    store = _store(tmp_path)
    report = validate("katG is in mtu00360 and also mtu00627.", store.citable_ids)
    assert [c.identifier for c in report.unsupported] == ["mtu00627"]
    assert not report.passed


def test_clean_summary_passes(tmp_path):
    store = _store(tmp_path)
    report = validate("katG (mtu:Rv1908c) is in mtu00360 and partners with 83332.Rv1909c.",
                      store.citable_ids)
    assert report.passed, report.summary_line()


def test_cross_target_citation_is_distinguished_from_a_hallucination():
    """An ID that was retrieved, but for a different gene. Plausible-looking, and
    invisible to a validator that only checks global membership."""
    report = validate("katG is in mtu99999 and mtu00777.",
                      citable_ids={"mtu00777", "mtu00360"},
                      per_target={"KATG": {"mtu00360"}, "FURA": {"mtu00777"}},
                      claimed_target="KATG")
    statuses = {c.identifier: c.status for c in report.citations}
    assert statuses["mtu99999"] == "unsupported"
    assert statuses["mtu00777"] == "cross_target"


def test_extracts_both_identifier_shapes():
    found = extract_citations("see mtu00360 and 83332.Rv1909c, plus noise like 42 and ABC")
    assert found == ["mtu00360", "83332.Rv1909c"]


# --- evidence --------------------------------------------------------------


def test_broad_pathway_is_not_treated_as_a_link():
    """mtu01100 holds ~17% of the M. tuberculosis genome."""
    assert classify_pathway(698, 4008)[0] == "broad"
    assert classify_pathway(11, 4008)[0] == "specific"


def test_pair_sharing_only_a_broad_pathway_gets_a_negative_verdict():
    ev = pair_evidence(
        "geneA", "geneB",
        pathways={"geneA": [{"record_id": "mtu01100", "name": "Metabolic pathways"}],
                  "geneB": [{"record_id": "mtu01100", "name": "Metabolic pathways"}]},
        partners={"geneA": [], "geneB": []},
        pathway_sizes={"mtu01100": 698}, genome_size=4008,
    )
    assert "not evidence of a mechanistic link" in ev.verdict


def test_pair_sharing_a_specific_pathway_is_reported_as_a_link():
    ev = pair_evidence(
        "geneA", "geneB",
        pathways={"geneA": [{"record_id": "mtu00983", "name": "Drug metabolism"}],
                  "geneB": [{"record_id": "mtu00983", "name": "Drug metabolism"}]},
        partners={"geneA": [], "geneB": []},
        pathway_sizes={"mtu00983": 11}, genome_size=4008,
    )
    assert "Share specific pathway" in ev.verdict


def test_direct_interaction_is_detected_and_its_support_characterised():
    ev = pair_evidence(
        "katG", "furA",
        pathways={}, pathway_sizes={}, genome_size=4008,
        partners={"katG": STRING_RESULT["records"], "furA": []},
    )
    assert ev.direct_interaction is not None
    assert "supported beyond literature co-mention" in ev.verdict


def test_no_evidence_is_an_explicit_verdict_not_silence():
    ev = pair_evidence("geneA", "geneB", pathways={}, partners={},
                       pathway_sizes={}, genome_size=4008)
    assert ev.verdict.startswith("No known link")


def test_all_pairs_covers_every_combination():
    assert len(all_pairs(["a", "b", "c"], {}, {}, {}, 4008)) == 3
