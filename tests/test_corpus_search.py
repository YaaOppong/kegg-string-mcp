"""Corpus search: the measured retrieval arm, exposed as a tool."""

from __future__ import annotations

import json
import pathlib
import tempfile

import pytest

from kegg_string_mcp.corpus_search import CORPUS_ENV, CorpusSearchClient

pytest.importorskip("chromadb", reason="needs the vector extra")
pytest.importorskip("rank_bm25", reason="needs the vector extra")

FIXTURE = "tests/fixtures/corpus_small.json"


@pytest.fixture
def client():
    return CorpusSearchClient(FIXTURE)


def test_no_corpus_configured_is_an_explicit_answer(monkeypatch):
    """Inert, not broken: the tool exists whether or not a corpus does, and says
    which it is rather than returning a bare empty list."""
    monkeypatch.delenv(CORPUS_ENV, raising=False)
    result = CorpusSearchClient().search("katG and ahpC")
    assert result.records == []
    assert "No local corpus is configured" in result.notes[0]
    assert "pubmed_abstracts" in result.notes[0]
    assert "not evidence" in result.notes[0]


def test_a_missing_corpus_file_names_itself(monkeypatch):
    monkeypatch.delenv(CORPUS_ENV, raising=False)
    result = CorpusSearchClient("data/nope.json").search("katG")
    assert result.records == []
    assert "does not exist" in result.notes[0]


def test_an_empty_query_does_not_build_an_index(client):
    result = client.search("   ")
    assert result.records == []
    assert client._index is None          # nothing was embedded


def test_search_returns_citable_pubmed_records(client):
    result = client.search("katG and ahpC in oxidative stress", limit=3)
    assert result.records
    for record in result.records:
        assert record.record_id.isdigit()          # a PMID, citable like any other
        assert record.type == "article"
        assert record.source == "pubmed"
        assert record.url.endswith(f"/{record.record_id}/")


def test_quotable_text_is_the_whole_abstract_not_the_matched_chunk():
    """The trap this avoids: the run store keeps the FIRST record for an ID and
    merges only `mentions`. A chunk arriving before the same PMID's full abstract
    would become the text every later quote is checked against, and a correct
    quote from elsewhere in that abstract would read as fabricated."""
    raw = json.loads(pathlib.Path(FIXTURE).read_text())
    # Split one paper into two overlapping chunks, as scripts/build_corpus.py does.
    original = raw["passages"][0]
    words = original["text"].split()
    half = len(words) // 2
    a = dict(original, passage_id=original["passage_id"] + "-0", chunk_index=0, chunk_of=2,
             text=" ".join(words[:half + 8]))
    b = dict(original, passage_id=original["passage_id"] + "-1", chunk_index=1, chunk_of=2,
             text=" ".join(words[half:]))
    raw["passages"] = [a, b] + raw["passages"][1:]

    with tempfile.TemporaryDirectory() as tmp:
        path = pathlib.Path(tmp) / "chunked.json"
        path.write_text(json.dumps(raw))
        result = CorpusSearchClient(path).search(original["title"], limit=5)

    record = next(r for r in result.records if r.record_id == original["pmid"])
    quotable = record.detail["quotable_text"]
    assert record.detail["passages_in_corpus"] == 2
    # Every span of the original abstract still contains-matches, which is the
    # property the validator relies on.
    assert " ".join(words[:6]) in quotable
    assert " ".join(words[-6:]) in quotable
    assert len(quotable) > len(record.detail["matched_passage"])


def test_records_are_one_per_paper_not_one_per_passage(client):
    result = client.search("katG", limit=4)
    ids = [r.record_id for r in result.records]
    assert len(ids) == len(set(ids))
    assert len(ids) <= 4


def test_notes_name_the_genes_the_corpus_covers(client):
    """A gene outside the corpus retrieves nothing relevant however it is queried,
    so the caller has to be able to see the boundary."""
    result = client.search("katG", limit=2)
    covered = " ".join(result.notes)
    assert "katG" in covered
    assert "pubmed_abstracts" in covered


def test_provenance_points_at_the_corpus_by_content(client):
    """A corpus is a file, not a fetch: its provenance is a content hash and the
    file's own mtime, never an invented retrieval time."""
    result = client.search("katG", limit=1)
    trace = result.requests[0]
    assert trace.url.startswith("file://")
    assert len(trace.content_sha256) >= 16
    assert trace.retrieved_at.startswith("20")


def test_the_index_is_built_once_and_reused(client):
    client.search("katG", limit=1)
    first = client._index
    client.search("ahpC", limit=1)
    assert client._index is first
