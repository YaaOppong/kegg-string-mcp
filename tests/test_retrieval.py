"""Retrieval arm tests.

The corpus builder, rank fusion and the comparison logic depend on nothing beyond
the standard library and run in the core matrix. Everything else needs part of the
optional `[vector]` extra and skips cleanly without it -- BM25 needs `rank_bm25`,
the graph needs `langgraph`, the dense index needs `chromadb`.

The fakes matter: the graph's control flow is tested against a scripted retriever
rather than a real index, so the cycle is under test rather than the retriever.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from kegg_string_mcp.retrieval.corpus import Corpus, Passage, build
from kegg_string_mcp.retrieval.index import Hit, reciprocal_rank_fusion

# Each arm needs a different piece of the optional extra, so they are gated
# separately rather than skipping the whole file.
bm25 = pytest.importorskip("rank_bm25", reason="needs the [vector] extra")
requires_graph = pytest.mark.skipif(
    __import__("importlib.util", fromlist=["util"]).find_spec("langgraph") is None,
    reason="needs the [vector] extra")

CORPUS = Path(__file__).resolve().parent.parent / "data" / "corpus_tb.json"


from kegg_string_mcp.retrieval.index import KeywordIndex


def passage(pid, text, mentions=(), title="t"):
    return Passage(passage_id=pid, pmid=pid, text=text, title=title, year="2020",
                   journal="J", doi="", pmcid="", in_pmc=False, mentions=list(mentions))


def tiny_corpus():
    return Corpus(genes=["katG", "ahpC"], passages=[
        passage("1", "katG encodes a catalase-peroxidase activating isoniazid", ["katG"]),
        passage("2", "ahpC is an alkyl hydroperoxide reductase under oxyR control", ["ahpC"]),
        passage("3", "katG loss is compensated by ahpC promoter mutations", ["katG", "ahpC"]),
        passage("4", "rifampicin resistance arises from rpoB mutations", ["rpoB"]),
    ])


# --- corpus ----------------------------------------------------------------


class FakePubMed:
    """Stands in for PubMedClient so the corpus builder needs no network."""

    def __init__(self, by_gene):
        self.by_gene = by_gene
        self.calls = []

    def abstracts(self, gene, organism="Mycobacterium tuberculosis", limit=20):
        from kegg_string_mcp.provenance import Record, ToolResult

        self.calls.append((gene, limit))
        records = [
            Record(record_id=pmid, type="article", name=f"paper {pmid}",
                   url=f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/", source="pubmed",
                   retrieved_at="2026-09-03T00:00:00+00:00", cached=False,
                   detail={"quotable_text": text, "title": f"paper {pmid}", "year": "2020",
                           "journal": "J", "doi": "", "pmcid": "", "in_pmc": False,
                           "mentions": mentions})
            for pmid, text, mentions in self.by_gene.get(gene, [])
        ]
        return ToolResult.build({"gene": gene}, records)


def test_corpus_deduplicates_papers_across_genes():
    """The same paper is routinely returned for several genes. Counting it once per
    gene would inflate any metric that depends on corpus size."""
    fake = FakePubMed({"katG": [("1", "katG text", ["katG"]), ("9", "shared", ["katG"])],
                       "ahpC": [("2", "ahpC text", ["ahpC"]), ("9", "shared", ["ahpC"])]})
    corpus = build(["katG", "ahpC"], client=fake)
    assert len(corpus.passages) == 3
    shared = next(p for p in corpus.passages if p.pmid == "9")
    assert sorted(shared.queried_for) == ["ahpC", "katG"]
    assert sorted(shared.mentions) == ["ahpC", "katG"]


def test_corpus_records_a_gene_with_no_literature():
    fake = FakePubMed({"katG": [("1", "katG text", ["katG"])]})
    corpus = build(["katG", "ghost"], client=fake)
    assert any("ghost" in n for n in corpus.notes)


def test_corpus_round_trips(tmp_path):
    original = tiny_corpus()
    path = original.write(tmp_path / "c.json")
    assert Corpus.read(path).passages[0].text == original.passages[0].text


def test_every_passage_keeps_its_provenance():
    """A passage surfaced by vector search must be as checkable as one from the
    keyword tools -- same PMID, same resolvable URL."""
    for p in tiny_corpus().passages:
        assert p.pmid and p.url.endswith(f"/{p.pmid}/")


# --- keyword arm -----------------------------------------------------------


def test_bm25_matches_exact_gene_symbols():
    hits = KeywordIndex(tiny_corpus()).search("ahpC hydroperoxide", k=3)
    assert hits[0].pmid in {"2", "3"}


def test_tokeniser_keeps_identifiers_intact():
    """Gene symbols and locus tags must survive tokenisation as single tokens."""
    tokens = KeywordIndex._tokenise("Rv1908c and beta-lactamase in M. tuberculosis")
    assert "rv1908c" in tokens
    assert "beta-lactamase" in tokens


def test_keyword_hits_carry_mentions_for_scoring():
    hits = KeywordIndex(tiny_corpus()).search("katG ahpC compensation", k=4)
    assert any(set(h.mentions) >= {"katG", "ahpC"} for h in hits)


# --- fusion ----------------------------------------------------------------


def test_rank_fusion_merges_on_rank_not_score():
    """BM25 scores and cosine similarities are not comparable -- different scales,
    different distributions -- so fusing on rank is the standard answer."""
    a = [Hit("x", "x", score=12.5, rank=1, title="", text=""),
         Hit("y", "y", score=11.0, rank=2, title="", text="")]
    b = [Hit("y", "y", score=0.44, rank=1, title="", text=""),
         Hit("z", "z", score=0.41, rank=2, title="", text="")]
    fused = reciprocal_rank_fusion([a, b], k=3)
    assert next(h.passage_id for h in fused) == "y", "ranked highly in both runs"
    assert {h.passage_id for h in fused} == {"x", "y", "z"}


def test_rank_fusion_is_stable_with_one_run():
    a = [Hit(str(i), str(i), score=1.0, rank=i, title="", text="") for i in range(1, 4)]
    assert [h.passage_id for h in reciprocal_rank_fusion([a], k=3)] == ["1", "2", "3"]


# --- the LangGraph cycle ---------------------------------------------------


class ScriptedRetriever:
    """Returns a fixed list per call, so the graph's control flow is under test
    rather than the retriever's behaviour."""

    def __init__(self, rounds):
        self.rounds = list(rounds)
        self.queries = []

    def search(self, query, k=10):
        self.queries.append(query)
        return self.rounds.pop(0) if self.rounds else []


_ids = iter(range(1, 10_000))


def hits_naming(*mention_sets):
    """Fresh passage ids each call: reusing them makes a later round look like it
    returned nothing new, which trips the graph's no-progress guard."""
    return [Hit(str(n), str(n), score=1.0, rank=i + 1, title="", text="", mentions=list(m))
            for i, (n, m) in enumerate((next(_ids), m) for m in mention_sets)]


@requires_graph
def test_graph_stops_when_one_paper_names_both_genes():
    from kegg_string_mcp.retrieval.graph import RetrievalGraph

    retriever = ScriptedRetriever([hits_naming(["katG", "ahpC"])])
    out = RetrievalGraph(retriever=retriever).run(["katG", "ahpC"])
    assert out["sufficient"] and out["rounds"] == 1


@requires_graph
def test_graph_rewrites_when_no_single_paper_covers_the_pair():
    """The union across hits is the wrong test: a corpus built one gene at a time
    almost always mentions both genes between its hits while no paper discusses the
    pair. Judging on the union made the cycle unreachable."""
    from kegg_string_mcp.retrieval.graph import RetrievalGraph

    retriever = ScriptedRetriever([hits_naming(["katG"], ["ahpC"]),
                                   hits_naming(["katG", "ahpC"])])
    out = RetrievalGraph(retriever=retriever).run(["katG", "ahpC"])
    assert out["rounds"] == 2 and out["sufficient"]
    assert len(retriever.queries) == 2 and retriever.queries[0] != retriever.queries[1]


@requires_graph
def test_graph_gives_up_rather_than_inventing_a_link():
    """'No joint evidence' is a result. The loop is bounded so a corpus that does
    not contain the answer terminates instead of spinning."""
    from kegg_string_mcp.retrieval.graph import RetrievalGraph

    retriever = ScriptedRetriever([hits_naming(["katG"]), hits_naming(["ahpC"]),
                                   hits_naming(["rpoB"])])
    out = RetrievalGraph(retriever=retriever, max_rounds=3).run(["katG", "ahpC"])
    assert not out["sufficient"] and out["rounds"] <= 3
    assert "never in the same paper" in out["reason"] or "nothing retrieved" in out["reason"]


@requires_graph
def test_graph_stops_when_a_rewrite_finds_nothing_new():
    """A rewrite that returns only passages already seen has exhausted its
    phrasing; spinning further rounds costs time for no evidence."""
    from kegg_string_mcp.retrieval.graph import RetrievalGraph

    same = hits_naming(["katG"])          # identical ids on every round, on purpose
    retriever = ScriptedRetriever([same, list(same), list(same)])
    out = RetrievalGraph(retriever=retriever, max_rounds=5).run(["katG", "ahpC"])
    assert out["rounds"] < 5


@requires_graph
def test_rewrites_vary_between_rounds():
    from kegg_string_mcp.retrieval.graph import RetrievalGraph

    first = RetrievalGraph._default_rewrite(["a", "b"], [], 1)
    second = RetrievalGraph._default_rewrite(["a", "b"], [], 2)
    assert first != second, "a repeated rewrite spins the loop for nothing"


# --- comparison ------------------------------------------------------------


def test_on_target_uses_what_the_text_names_not_what_was_queried():
    """PubMed matches on metadata the reader never sees, so a paper can be returned
    for a gene it never names. Scoring on the query would call that a hit."""
    from kegg_string_mcp.retrieval.compare import on_target

    named, unnamed = hits_naming(["katG"], [])
    assert on_target([named, unnamed], ["katG"]) == [named]


def test_comparison_reports_overlap_and_precision():
    from kegg_string_mcp.retrieval.compare import compare

    corpus = tiny_corpus()
    arms = {"keyword": KeywordIndex(corpus)}
    result = compare(arms, [("katG ahpC compensation", ["katG", "ahpC"])], k=3)
    assert result.results[0].precision["keyword"] > 0
    assert result.mean_overlap("keyword", "keyword") == 1.0


# --- dense index, needs the optional extra ---------------------------------


@pytest.mark.skipif(not CORPUS.exists(), reason="corpus not built")
def test_dense_and_keyword_arms_disagree_on_the_real_corpus():
    """If the arms returned the same thing there would be nothing to measure."""
    pytest.importorskip("chromadb")
    from kegg_string_mcp.retrieval.index import VectorIndex

    corpus = Corpus.read(CORPUS)
    query = "how does katG relate to ahpC in isoniazid resistance?"
    kw = {h.pmid for h in KeywordIndex(corpus).search(query, k=10)}
    vec = {h.pmid for h in VectorIndex(corpus).search(query, k=10)}
    assert kw and vec
    assert kw != vec, "the arms are interchangeable, so the comparison is vacuous"


@pytest.mark.skipif(not CORPUS.exists(), reason="corpus not built")
def test_recorded_comparison_has_the_headline_numbers():
    path = CORPUS.parent / "comparison.json"
    if not path.exists():
        pytest.skip("comparison not run")
    data = json.loads(path.read_text())
    assert set(data["arms"]) == {"keyword", "vector", "hybrid"}
    assert data["queries"] >= 10
    for arm, value in data["mean_on_target_precision"].items():
        assert 0.0 <= value <= 1.0, arm


# --- index identity --------------------------------------------------------


def test_corpus_fingerprint_tracks_content_not_metadata():
    """Re-running the corpus build reorders `queried_for` without changing what is
    searchable, so the fingerprint must ignore it or every rebuild re-embeds."""
    from kegg_string_mcp.retrieval.index import corpus_fingerprint

    a = tiny_corpus()
    b = tiny_corpus()
    b.passages[0].queried_for = ["ahpC", "katG"]
    b.passages[0].mentions = ["katG", "extra"]
    assert corpus_fingerprint(a) == corpus_fingerprint(b)

    c = tiny_corpus()
    c.passages[0].text += " and more"
    assert corpus_fingerprint(a) != corpus_fingerprint(c)


@pytest.mark.skipif(
    __import__("importlib.util", fromlist=["util"]).find_spec("chromadb") is None,
    reason="needs the [vector] extra")
def test_a_second_corpus_does_not_inherit_the_first_ones_index(tmp_path):
    """Keyed on collection name alone, a second corpus silently reused the first
    one's vectors: the load was skipped because the collection was non-empty, and
    retrieval then answered from the wrong text. Nothing failed -- the comparison
    numbers would simply have been wrong."""
    from kegg_string_mcp.retrieval.index import VectorIndex

    first = Corpus(genes=["katG"], passages=[passage("1", "katG catalase peroxidase")])
    second = Corpus(genes=["gyrA"], passages=[passage("99", "gyrA gyrase supercoiling")])

    VectorIndex(first, path=str(tmp_path))
    hits = VectorIndex(second, path=str(tmp_path)).search("gyrase", k=2)
    assert [h.passage_id for h in hits] == ["99"]

    back = VectorIndex(first, path=str(tmp_path)).search("catalase", k=2)
    assert [h.passage_id for h in back] == ["1"]


@pytest.mark.skipif(
    __import__("importlib.util", fromlist=["util"]).find_spec("chromadb") is None,
    reason="needs the [vector] extra")
def test_empty_corpus_returns_nothing_rather_than_raising():
    from kegg_string_mcp.retrieval.index import VectorIndex

    assert VectorIndex(Corpus(genes=[], passages=[])).search("anything", k=5) == []
