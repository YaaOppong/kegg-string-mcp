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
