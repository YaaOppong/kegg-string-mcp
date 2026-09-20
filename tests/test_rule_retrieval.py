"""Stage B: retrieval for the questions the structured sources could not answer.

No network here -- the retrievers are stand-ins. What is tested is the query, the
order of the two sources, and that nothing in this stage decides anything.
"""

from pathlib import Path

import pytest

from kegg_string_mcp.rules import parse_annotation
from kegg_string_mcp.rules.questions import COMPENSATION, PRIOR_ART, RESISTANCE_ROLE, Question
from kegg_string_mcp.rules.retrieve import CORPUS, PUBMED, Clients, build_query, retrieve

BED = [
    ("AL123456", 99, 3617, "Rv0667", ".", "+", "rpoB"),
    ("AL123456", 3699, 7700, "Rv0668", ".", "+", "rpoC"),
    ("AL123456", 7999, 10221, "Rv1908c", ".", "-", "katG"),
    ("AL123456", 10399, 10800, "Rv1909c", ".", "-", "furA"),
    ("AL123456", 19999, 20743, "Rv1483", ".", "+", "fabG1"),
    ("AL123456", 20799, 21608, "Rv1484", ".", "+", "inhA"),
]


@pytest.fixture
def annotation(tmp_path: Path):
    path = tmp_path / "genes.bed"
    path.write_text("\n".join("\t".join(str(c) for c in r) for r in BED) + "\n")
    return parse_annotation(path)


class _Source:
    def __init__(self, records=None, raises=None):
        self.records, self.raises, self.calls = records or [], raises, []

    def __call__(self, query, organism=None, limit=10):
        self.calls.append(query)
        if self.raises:
            raise self.raises
        return {"records": self.records}


def _record(pmid, title="T", year="2020", score=0.9):
    return {"record_id": pmid, "detail": {"title": title, "year": year, "score": score}}


# --- the query -------------------------------------------------------------


def test_a_query_names_every_spelling_of_a_locus(annotation):
    """A paper writes katG, not Rv1908c, and a vocabulary may carry either."""
    query = build_query(Question(kind=RESISTANCE_ROLE, locus="Rv1908c",
                                 drugs=("isoniazid",)), annotation)
    assert '("katG" OR "Rv1908c")' in query
    assert '"isoniazid"' in query
    assert '"Mycobacterium tuberculosis"' in query


def test_a_query_never_names_the_conclusion(annotation):
    """A query containing `compensatory` retrieves papers that use the word, so
    finding a passage saying so becomes near-certain and the verdict circular.
    The query asks for papers about the loci; a quote settles the rest."""
    query = build_query(Question(kind=COMPENSATION, locus="Rv0668", partner="Rv0667",
                                 drugs=("rifampicin",)), annotation)
    for word in ("compensat", "fitness", "resistance", "evidence"):
        assert word not in query.lower()
    assert '("rpoC" OR "Rv0668")' in query and '("rpoB" OR "Rv0667")' in query


def test_a_pair_question_is_not_narrowed_by_the_drug(annotation):
    """Two loci are specific enough. Adding the drug would restrict a pair query
    to papers already framing them as a resistance story, which is the reading
    being tested."""
    paired = build_query(Question(kind=PRIOR_ART, locus="Rv0667", partner="Rv0668",
                                  drugs=("rifampicin",)), annotation)
    assert "rifampicin" not in paired

    single = build_query(Question(kind=RESISTANCE_ROLE, locus="Rv0667",
                                  drugs=("rifampicin",)), annotation)
    assert "rifampicin" in single


def test_an_intergenic_feature_is_queried_by_its_flanking_genes(annotation):
    """No paper names a gap; a gene's promoter is discussed under the gene."""
    query = build_query(Question(kind=RESISTANCE_ROLE, locus="Rv1483-Rv1484",
                                 drugs=("isoniazid",)), annotation)
    assert "fabG1" in query and "inhA" in query


# --- the retrieval ---------------------------------------------------------


def test_pubmed_is_only_asked_where_the_corpus_reached_nothing(annotation):
    corpus = _Source([_record("111")])
    pubmed = _Source([_record("222")])
    question = Question(kind=RESISTANCE_ROLE, locus="Rv1908c", drugs=("isoniazid",), text="q")

    report = retrieve([question], annotation, Clients(corpus=corpus, pubmed=pubmed))
    assert [c.source for c in report.candidates] == [CORPUS]
    assert pubmed.calls == []

    empty = _Source([])
    report = retrieve([question], annotation, Clients(corpus=empty, pubmed=pubmed))
    assert [c.source for c in report.candidates] == [PUBMED]
    assert pubmed.calls


def test_a_retriever_that_fails_does_not_take_the_batch_down(annotation):
    """A few hundred questions must not be all-or-nothing."""
    question = Question(kind=RESISTANCE_ROLE, locus="Rv1908c", text="q")
    report = retrieve([question], annotation,
                      Clients(corpus=_Source(raises=ValueError("boom")),
                              pubmed=_Source([_record("222")])))
    assert [c.record_id for c in report.candidates] == ["222"]
    assert any("corpus unavailable" in n for n in report.notes)


def test_the_headline_number_is_how_many_questions_got_nothing(annotation):
    """The number that sizes the corpus-widening job."""
    questions = [Question(kind=RESISTANCE_ROLE, locus="Rv1908c", text="a"),
                 Question(kind=RESISTANCE_ROLE, locus="Rv0667", text="b")]
    report = retrieve(questions, annotation, Clients(pubmed=_Source([])))
    assert report.summary()["with_a_candidate"] == 0
    assert report.summary()["without"] == 2


def test_every_retrieval_is_written_to_the_store(annotation, tmp_path):
    """Stage D validates quotes against the pipeline's copy of the text, never
    the model's account of it. Without the store there is nothing to check."""
    from kegg_string_mcp.agent.store import RunStore

    store = RunStore(path=tmp_path / "r.jsonl", run_id="r")
    question = Question(kind=RESISTANCE_ROLE, locus="Rv1908c", text="q")
    retrieve([question], annotation, Clients(pubmed=_Source([_record("222")])), store=store)
    assert "222" in store.citable_ids


def test_an_unset_ncbi_email_is_said_out_loud(annotation, monkeypatch):
    """One unidentified lookup is unremarkable; a few hundred is where NCBI
    notices, and the failure mode is a throttle with no local signal."""
    monkeypatch.delenv("NCBI_EMAIL", raising=False)
    report = retrieve([], annotation, Clients(pubmed=_Source([])))
    assert any("NCBI_EMAIL is not set" in n for n in report.notes)

    monkeypatch.setenv("NCBI_EMAIL", "you@example.org")
    report = retrieve([], annotation, Clients(pubmed=_Source([])))
    assert not any("NCBI_EMAIL" in n for n in report.notes)


def test_the_model_facing_tool_still_refuses_query_syntax(http):
    """`abstracts()` phrase-quotes its argument and refuses `"`, `[`, `]`, because
    a model must not be able to change the meaning of a search. Composing a
    boolean query is a different caller, and gets a different door -- but the
    guard on the model-facing one is unchanged."""
    from kegg_string_mcp.pubmed import PubMedClient

    result = PubMedClient(http).abstracts('("katG" OR "Rv1908c") AND "isoniazid"')
    assert result.records == []
    assert any("PubMed query syntax" in n for n in result.notes)


def test_an_empty_term_is_refused_rather_than_searched(http):
    from kegg_string_mcp.pubmed import PubMedClient

    result = PubMedClient(http).search_abstracts("   ")
    assert result.records == []
    assert any("no search term" in n for n in result.notes)
