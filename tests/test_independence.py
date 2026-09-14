"""The classifier decides which pairs the literature can speak to independently."""

from __future__ import annotations

import pathlib
from types import SimpleNamespace

import pytest

from kegg_string_mcp.identity import GeneIdentity, IdentitySet
from kegg_string_mcp.retrieval.corpus import Corpus
from kegg_string_mcp.retrieval.independence import (
    STRING_RELEASE_YEAR,
    Edge,
    IndependenceReport,
    StringEdges,
    classify,
    network_edges,
    partner_edges,
    post_release_fraction,
)

FIXTURE = pathlib.Path(__file__).parent / "fixtures" / "corpus_small.json"


def _identities(**genes) -> IdentitySet:
    """{'katG': ('83332.Rv1908c', ['katG', 'Rv1908c'])} -> an IdentitySet."""
    return IdentitySet(identities={
        query: GeneIdentity(query=query, string_id=string_id, aliases=list(aliases),
                            matched_by={"string": "get_string_ids" if string_id else "none"})
        for query, (string_id, aliases) in genes.items()})


def _edges(identities: IdentitySet, source: str = "network", truncated=(), **pairs) -> StringEdges:
    """{'katG__ahpC': (combined, textmining, max_non_tm)} -> a StringEdges."""
    out = StringEdges(identities=identities, source=source, truncated=list(truncated))
    for key, (combined, textmining, max_non_tm) in pairs.items():
        a, b = key.split("__")
        ids = tuple(sorted((identities.string_id(a), identities.string_id(b))))
        out.edges[ids] = Edge(combined=combined, textmining=textmining,
                              max_non_textmining=max_non_tm)
    return out


KATG_AHPC = {"katG": ("83332.Rv1908c", ["katG", "Rv1908c"]),
             "ahpC": ("83332.Rv2428", ["ahpC", "Rv2428"])}


def test_no_edge_is_silent():
    identities = _identities(**KATG_AHPC)
    verdicts = classify(["katG", "ahpC"], _edges(identities))
    assert [v.status for v in verdicts] == ["silent"]


def test_strong_non_textmining_channel_is_corroboration():
    identities = _identities(**KATG_AHPC)
    edges = _edges(identities, katG__ahpC=(0.97, 0.96, 0.62))
    assert classify(["katG", "ahpC"], edges)[0].status == "corroborating"


def test_textmining_only_edge_is_not_independent_evidence():
    """A high combined score built from textmining is the literature, restated.

    Counting the retrieved abstracts as evidence on top of STRING's score would be
    one line of evidence twice -- the failure the repo already guards against at
    the tool level.
    """
    identities = _identities(**KATG_AHPC)
    edges = _edges(identities, katG__ahpC=(0.968, 0.965, 0.096))
    verdict = classify(["katG", "ahpC"], edges)[0]
    assert verdict.status == "textmining_only"
    assert verdict.combined > 0.9


def test_verdict_is_the_same_whichever_name_the_caller_used():
    """The bug this module was rebuilt for. `Rv0678` and `mmpR5` are one protein;
    keying partners by preferred name made the locus-tag spelling read as silent
    (0.0) while the symbol spelling read as corroborating (0.989)."""
    identities = IdentitySet(identities={
        "Rv0678": GeneIdentity(query="Rv0678", string_id="83332.Rv0678",
                               aliases=["Rv0678", "mmpR5"], matched_by={"string": "x"}),
        "Rv0676c": GeneIdentity(query="Rv0676c", string_id="83332.Rv0676c",
                                aliases=["Rv0676c", "mmpL5"], matched_by={"string": "x"}),
    })
    edges = StringEdges(identities=identities)
    edges.edges[("83332.Rv0676c", "83332.Rv0678")] = Edge(0.989, 0.30, 0.81)

    by_locus = classify(["Rv0678", "Rv0676c"], edges)[0]
    by_symbol = classify(["mmpR5", "mmpL5"], edges)[0]
    assert by_locus.status == by_symbol.status == "corroborating"
    assert by_locus.combined == by_symbol.combined == 0.989


def test_unresolved_gene_is_not_silent():
    """`string_partners` already says a resolution failure is not evidence of no
    partners. Calling it silent handed the residue a failed lookup as a negative
    result, and `hypothesis/residue.py` then counted the pair as a candidate."""
    identities = _identities(katG=("83332.Rv1908c", ["katG"]), fakeGene1=("", ["fakeGene1"]))
    verdict = classify(["katG", "fakeGene1"], _edges(identities))[0]
    assert verdict.status == "unresolved"
    assert verdict.undetermined
    assert "not evidence of no interaction" in verdict.note


def test_two_full_partner_lists_report_truncated_not_silent():
    """At the module defaults every gene tested returned exactly `limit` partners,
    so an edge below both cuts is invisible -- which is not the same as absent."""
    identities = _identities(hubA=("83332.A", ["hubA"]), hubB=("83332.B", ["hubB"]))
    edges = _edges(identities, source="partners", truncated=("hubA", "hubB"))
    assert classify(["hubA", "hubB"], edges)[0].status == "truncated"


def test_one_truncated_list_is_still_a_real_answer():
    """Only both lists being full hides an edge; one full list cannot."""
    identities = _identities(hubA=("83332.A", ["hubA"]), small=("83332.B", ["small"]))
    edges = _edges(identities, source="partners", truncated=("hubA",))
    assert classify(["hubA", "small"], edges)[0].status == "silent"


def test_every_pair_is_classified_exactly_once():
    identities = _identities(**{g: (f"83332.{g}", [g]) for g in "abcd"})
    verdicts = classify(list("abcd"), _edges(identities))
    assert len(verdicts) == 6
    assert len({(v.gene_a, v.gene_b) for v in verdicts}) == 6


class FakeString:
    """Enough of StringClient to drive both edge builders offline."""

    def __init__(self, ids: dict[str, str], rows: list[dict] | None = None,
                 partners: dict[str, list] | None = None):
        self.ids = ids
        self.rows = rows or []
        self.partners_by_gene = partners or {}
        self.network_calls: list[list[str]] = []
        self.partner_calls: list[str] = []

    def resolve(self, gene, species):
        string_id = self.ids.get(gene, "")
        hit = {"stringId": string_id, "preferredName": gene} if string_id else None
        return hit, SimpleNamespace(url="", retrieved_at="", cached=True, status=200,
                                    content_sha256="")

    def network(self, identifiers, species=83332, required_score=150):
        self.network_calls.append(list(identifiers))
        return SimpleNamespace(records=[SimpleNamespace(detail=row) for row in self.rows],
                               notes=[])

    def partners(self, gene, species=83332, limit=200, required_score=150):
        self.partner_calls.append(gene)
        return SimpleNamespace(records=self.partners_by_gene.get(gene, []), notes=[])


def test_network_edges_asks_once_for_the_whole_set():
    """820 pairs would be 820 requests against a service that asks for roughly one
    request a second -- and /network answers every pair in a single call."""
    genes = [f"g{i}" for i in range(10)]
    string = FakeString(ids={g: f"83332.{g}" for g in genes}, rows=[
        {"string_id_a": "83332.g0", "string_id_b": "83332.g1", "combined_score": 0.9,
         "textmining_score": 0.1, "max_non_textmining_score": 0.8}])
    edges = network_edges(genes, string)
    assert len(string.network_calls) == 1
    assert len(string.network_calls[0]) == 10
    assert classify(genes, edges)[0].status == "corroborating"


def test_network_edges_records_genes_string_could_not_resolve():
    string = FakeString(ids={"katG": "83332.Rv1908c"})
    edges = network_edges(["katG", "nosuchgene"], string)
    assert edges.identities.unresolved("string") == ["nosuchgene"]
    assert classify(["katG", "nosuchgene"], edges)[0].status == "unresolved"


def test_partner_edges_are_keyed_by_string_id_not_by_name():
    """The partner path is the fallback, so it must not reintroduce name keying."""
    partner = SimpleNamespace(record_id="83332.Rv0676c", name="mmpL5",
                              detail={"combined_score": 0.989, "textmining_score": 0.3,
                                      "max_non_textmining_score": 0.81})
    string = FakeString(ids={"Rv0678": "83332.Rv0678", "Rv0676c": "83332.Rv0676c"},
                        partners={"Rv0678": [partner]})
    edges = partner_edges(["Rv0678", "Rv0676c"], string, limit=200)
    assert classify(["Rv0678", "Rv0676c"], edges)[0].status == "corroborating"


def test_partner_edges_skip_the_call_for_an_unresolved_gene():
    """A gene STRING cannot resolve has no partner list to fetch, and asking
    anyway spends a request to produce a result that means nothing."""
    string = FakeString(ids={"katG": "83332.Rv1908c"})
    partner_edges(["katG", "nosuchgene"], string)
    assert string.partner_calls == ["katG"]


def test_report_counts_and_filters_by_status():
    identities = _identities(**{g: (f"83332.{g}", [g]) for g in "abc"})
    edges = _edges(identities, a__b=(0.9, 0.89, 0.05))
    report = IndependenceReport(verdicts=classify(list("abc"), edges))
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
    identities = _identities(a=("83332.a", ["a"]), b=("83332.b", ["b"]))
    edges = _edges(identities) if edge is None else _edges(identities, a__b=edge)
    verdict = classify(["a", "b"], edges)[0]
    assert verdict.status == status
    assert verdict.note


def test_serialised_verdict_shape():
    """The independence artefact is read by scripts/residue.py."""
    identities = _identities(a=("83332.a", ["a"]), b=("83332.b", ["b"]))
    payload = classify(["a", "b"], _edges(identities, a__b=(0.9, 0.89, 0.05)))[0].to_dict()
    assert set(payload) == {"gene_a", "gene_b", "combined", "textmining",
                            "max_non_textmining", "status", "note", "string_ids"}
    assert payload["status"] == "textmining_only"
    assert payload["string_ids"] == ["83332.a", "83332.b"]
