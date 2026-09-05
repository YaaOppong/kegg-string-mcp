"""The classifier decides which pairs the literature can speak to independently."""

from __future__ import annotations

import pathlib
from types import SimpleNamespace

import pytest

from kegg_string_mcp.retrieval.corpus import Corpus
from kegg_string_mcp.retrieval.independence import (
    STRING_RELEASE_YEAR,
    IndependenceReport,
    classify,
    gene_partner_map,
    post_release_fraction,
)

FIXTURE = pathlib.Path(__file__).parent / "fixtures" / "corpus_small.json"


def _partners(**edges):
    """{'katG': {'ahpC': (combined, textmining, max_non_tm)}} -> partner map."""
    return {gene: {name.lower(): {"combined": c, "textmining": t, "max_non_textmining": n}
                   for name, (c, t, n) in partners.items()}
            for gene, partners in edges.items()}


def test_no_edge_is_silent():
    verdicts = classify(["katG", "ahpC"], _partners(katG={}, ahpC={}))
    assert [v.status for v in verdicts] == ["silent"]


def test_strong_non_textmining_channel_is_corroboration():
    partners = _partners(katG={"ahpC": (0.97, 0.96, 0.62)}, ahpC={})
    assert classify(["katG", "ahpC"], partners)[0].status == "corroborating"


def test_textmining_only_edge_is_not_independent_evidence():
    """A high combined score built from textmining is the literature, restated.

    Counting the retrieved abstracts as evidence on top of STRING's score would be
    one line of evidence twice -- the failure the repo already guards against at
    the tool level.
    """
    partners = _partners(katG={"ahpC": (0.968, 0.965, 0.096)}, ahpC={})
    verdict = classify(["katG", "ahpC"], partners)[0]
    assert verdict.status == "textmining_only"
    assert verdict.combined > 0.9


def test_edge_found_from_either_direction():
    """STRING returns partner lists, not a matrix; a limit can truncate one side."""
    partners = _partners(katG={}, ahpC={"katG": (0.97, 0.10, 0.90)})
    assert classify(["katG", "ahpC"], partners)[0].status == "corroborating"


def test_every_pair_is_classified_exactly_once():
    genes = ["a", "b", "c", "d"]
    verdicts = classify(genes, _partners(**{g: {} for g in genes}))
    assert len(verdicts) == 6
    assert len({(v.gene_a, v.gene_b) for v in verdicts}) == 6


def test_partner_map_makes_one_call_per_gene_not_per_pair():
    """820 pairs would be 820 requests against a rate-limited service."""
    calls = []

    class FakeString:
        def partners(self, gene, **kwargs):
            calls.append(gene)
            return SimpleNamespace(records=[
                SimpleNamespace(name="ahpC", detail={"combined_score": 0.9,
                                                     "textmining_score": 0.9,
                                                     "max_non_textmining_score": 0.05})])

    genes = [f"g{i}" for i in range(10)]
    out = gene_partner_map(genes, FakeString())
    assert calls == genes
    assert out["g0"]["ahpc"]["combined"] == 0.9


def test_report_counts_and_filters_by_status():
    partners = _partners(a={"b": (0.9, 0.89, 0.05)}, b={}, c={})
    report = IndependenceReport(verdicts=classify(["a", "b", "c"], partners))
    assert report.by_status() == {"textmining_only": 1, "silent": 2}
    assert ("a", "c") in report.pairs_with_status("silent")
    assert report.pairs_with_status("textmining_only") == [("a", "b")]


def test_post_release_fraction_counts_papers_not_chunks():
    corpus = Corpus.read(FIXTURE)
    recent, total = post_release_fraction(corpus)
    assert total == len({p.pmid for p in corpus.passages})
    assert 0 <= recent <= total


def test_post_release_uses_the_string_release_year():
    corpus = SimpleNamespace(passages=[
        SimpleNamespace(pmid="1", year=str(STRING_RELEASE_YEAR - 1)),
        SimpleNamespace(pmid="2", year=str(STRING_RELEASE_YEAR)),
        SimpleNamespace(pmid="3", year=str(STRING_RELEASE_YEAR + 1)),
        SimpleNamespace(pmid="4", year=""),          # missing dates are not recent
    ])
    assert post_release_fraction(corpus) == (2, 4)


@pytest.mark.parametrize("status,edge", [
    ("silent", None),
    ("textmining_only", (0.9, 0.89, 0.05)),
    ("corroborating", (0.9, 0.10, 0.80)),
])
def test_every_status_carries_a_note_explaining_it(status, edge):
    partners = _partners(a={} if edge is None else {"b": edge}, b={})
    verdict = classify(["a", "b"], partners)[0]
    assert verdict.status == status
    assert verdict.note
