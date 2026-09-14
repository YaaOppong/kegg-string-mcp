"""The residue gate decides what reaches hypothesis generation."""

from __future__ import annotations

import pytest

from kegg_string_mcp.hypothesis.residue import (
    CO_MENTIONED,
    DEFAULT_EXPLAINING,
    SHARED_PATHWAY,
    STRING_EXPERIMENTAL,
    STRING_TEXTMINING,
    assess,
    residue,
    summarise,
    undetermined,
)

PAIRS = [("katG", "ahpC"), ("rpoB", "rpoC"), ("pks13", "zur")]


def test_unexplained_pair_reaches_the_residue():
    out = assess([("pks13", "zur")])
    assert out[0].reasons == []
    assert out[0].is_residue()


def test_non_textmining_string_edge_explains_a_pair():
    status = {("rpoB", "rpoC"): {"status": "corroborating", "max_non_textmining": 0.999}}
    out = assess([("rpoB", "rpoC")], string_status=status)
    assert out[0].codes() == {STRING_EXPERIMENTAL}
    assert not out[0].is_residue()


def test_textmining_alone_does_not_explain_by_default():
    """katG-pncA scores 0.965 textmining because both appear in resistance review
    tables. Discarding it as 'already known' would lose a candidate to a shared
    table row, so the reason is recorded but does not count."""
    status = {("katG", "pncA"): {"status": "textmining_only", "textmining": 0.965}}
    out = assess([("katG", "pncA")], string_status=status)
    assert out[0].codes() == {STRING_TEXTMINING}
    assert out[0].is_residue()
    # ...but it is available to anyone who wants it to count.
    assert not out[0].is_residue(frozenset({STRING_TEXTMINING}))


def test_shared_pathway_is_recorded_but_does_not_explain():
    """Two genes in one pathway can still be an undiscovered compensatory pair."""
    out = assess([("embA", "embB")], pathways={"embA": {"mtu01501"}, "embB": {"mtu01501"}})
    assert out[0].codes() == {SHARED_PATHWAY}
    assert out[0].is_residue()
    assert SHARED_PATHWAY not in DEFAULT_EXPLAINING


def test_co_mention_explains_a_pair():
    out = assess([("katG", "ahpC")], co_mentions={("katG", "ahpC"): 7})
    assert out[0].codes() == {CO_MENTIONED}
    assert not out[0].is_residue()
    assert out[0].reasons[0].value == 7.0


def test_zero_co_mentions_is_not_a_reason():
    out = assess([("katG", "ahpC")], co_mentions={("katG", "ahpC"): 0})
    assert out[0].reasons == []


def test_pair_order_and_case_do_not_matter():
    """Callers build pair keys from different sources; none should have to agree
    on ordering or capitalisation for a lookup to hit."""
    out = assess([("ahpC", "KATG")], co_mentions={("katg", "ahpc"): 3})
    assert out[0].codes() == {CO_MENTIONED}


def test_all_applicable_reasons_are_attached_not_just_the_first():
    status = {("katG", "ahpC"): {"status": "corroborating", "max_non_textmining": 0.8}}
    out = assess([("katG", "ahpC")], string_status=status,
                 pathways={"katG": {"mtu01501"}, "ahpC": {"mtu01501"}},
                 co_mentions={("katG", "ahpC"): 4})
    assert out[0].codes() == {STRING_EXPERIMENTAL, SHARED_PATHWAY, CO_MENTIONED}


def test_residue_is_recomputable_without_refetching():
    """The whole point of recording every reason: change the definition of
    'explained' and re-filter the same assessments."""
    out = assess(PAIRS, string_status={
        ("katG", "ahpC"): {"status": "textmining_only", "textmining": 0.965},
        ("rpoB", "rpoC"): {"status": "corroborating", "max_non_textmining": 0.999}})
    assert [(a.gene_a, a.gene_b) for a in residue(out)] == [("katG", "ahpC"), ("pks13", "zur")]
    strict = frozenset(DEFAULT_EXPLAINING | {STRING_TEXTMINING})
    assert [(a.gene_a, a.gene_b) for a in residue(out, strict)] == [("pks13", "zur")]


def test_summary_key_order_does_not_depend_on_set_iteration():
    """Counting over a set of codes made the written JSON differ between runs by
    key order alone -- same values, different order, which defeats diffing an
    artefact against itself."""
    out = assess(PAIRS, string_status={
        ("katG", "ahpC"): {"status": "corroborating", "max_non_textmining": 0.8}},
        pathways={"katG": {"mtu01501"}, "ahpC": {"mtu01501"}},
        co_mentions={("katG", "ahpC"): 2})
    counts = summarise(out)["reason_counts"]
    assert list(counts) == sorted(counts)


def test_summary_reports_the_configuration_it_used():
    out = assess(PAIRS, co_mentions={("katG", "ahpC"): 2})
    summary = summarise(out)
    assert summary["pairs"] == 3
    assert summary["residue"] == 2
    assert summary["residue_fraction"] == 0.667
    assert summary["reason_counts"] == {CO_MENTIONED: 1}
    assert summary["explaining"] == sorted(DEFAULT_EXPLAINING)


def test_summary_of_nothing_does_not_divide_by_zero():
    assert summarise([])["residue_fraction"] == 0.0


@pytest.mark.parametrize("status", ["silent", "unknown", None])
def test_string_statuses_that_assert_nothing_add_no_reason(status):
    out = assess([("a", "b")], string_status={("a", "b"): {"status": status}})
    assert out[0].reasons == []


def test_serialised_shape_is_what_downstream_reads():
    """to_dict is the artefact format scripts/residue.py writes; a renamed or
    dropped field would break a consumer with nothing failing here."""
    out = assess([("katG", "ahpC")], co_mentions={("katG", "ahpC"): 2})
    payload = out[0].to_dict()
    assert set(payload) == {"gene_a", "gene_b", "reasons", "undetermined"}
    assert payload["reasons"] == [
        {"code": CO_MENTIONED, "detail": "2 corpus paper(s) name both genes", "value": 2.0}]


# -- the third answer: pairs nothing could speak to ---------------------------


def test_an_unresolved_pair_is_undetermined_not_a_candidate():
    """`classify()` used to call a resolution failure 'silent', which produced no
    reasons here, which promoted the pair to novel candidate -- the strongest
    claim in the pipeline, made from a failed lookup."""
    status = {("katG", "fakeGene1"): {"status": "unresolved",
                                      "note": "fakeGene1 did not resolve"}}
    out = assess([("katG", "fakeGene1")], string_status=status)
    assert out[0].is_undetermined()
    assert not out[0].is_residue()
    assert undetermined(out) == out


def test_a_truncated_lookup_is_undetermined_too():
    status = {("hubA", "hubB"): {"status": "truncated", "note": "both lists were full"}}
    assert assess([("hubA", "hubB")], string_status=status)[0].is_undetermined()


def test_a_source_that_never_answered_makes_the_pair_undetermined():
    out = assess([("katG", "ahpC")], unanswered={"ahpC": ["kegg", "uniprot"]})
    assert out[0].is_undetermined()
    assert "never answered for ahpC" in out[0].undetermined[0].detail


def test_an_explained_pair_stays_explained_even_if_another_source_failed():
    """Undetermined is about what is unknown; it does not erase what is known."""
    out = assess([("katG", "ahpC")], co_mentions={("katG", "ahpC"): 7},
                 unanswered={"ahpC": ["kegg"]})
    assert out[0].codes() == {CO_MENTIONED}
    assert out[0].is_undetermined()
    assert not out[0].is_residue()


def test_undetermined_pairs_are_kept_out_of_the_residue_fraction():
    """A residue fraction computed over pairs nobody could assess measures the
    outage, not the biology."""
    out = assess([("a", "b"), ("c", "d"), ("e", "f")],
                 string_status={("a", "b"): {"status": "unresolved"}},
                 co_mentions={("c", "d"): 3})
    summary = summarise(out)
    assert summary["pairs"] == 3
    assert summary["undetermined"] == 1
    assert summary["assessable"] == 2
    assert summary["residue"] == 1
    assert summary["residue_fraction"] == 0.5


def test_summary_of_nothing_assessable_does_not_divide_by_zero():
    out = assess([("a", "b")], string_status={("a", "b"): {"status": "unresolved"}})
    assert summarise(out)["residue_fraction"] == 0.0
