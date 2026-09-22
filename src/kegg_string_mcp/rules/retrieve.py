"""Stage B: fetch candidate passages for each question. No model, no judgement.

Retrieval only. Nothing here decides whether a question is answered -- that needs
a verbatim quote, which is stage C, checked by stage D. This stage's job is to
put the candidate text in front of those stages and to record it where a quote
can be validated against it.

**Queries are gene-centric, never mechanism-centric.** A query containing
`compensatory` retrieves papers that use the word, so finding a passage saying
so becomes near-certain and the verdict circular -- the same circularity
`docs/RETRIEVAL.md` identifies when relevance is scored by whether a passage
names the queried genes. The query asks for papers about the loci; whether a
relationship is *stated* is settled later, by a quote.

**Every alias, not one spelling.** A paper writes `ahpC`, not `Rv2428`, and the
two halves of a vocabulary disagree about which it uses. The annotation supplies
both free, so both are OR'd into the query. This is `NEXT_STEPS.md` workstream 1
applied to questions rather than to a corpus build.

Corpus first, PubMed second. The corpus ranks better over what it holds -- 0.917
precision@10 against 0.844 for lexical alone -- but it covers a prebuilt gene
set, and a rule set names loci far outside it. An empty corpus result means the
corpus does not reach the question, never that no literature exists.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from typing import Any

from kegg_string_mcp.http import FetchError
from kegg_string_mcp.rules.annotation import Annotation
from kegg_string_mcp.rules.questions import COMPENSATION, PRIOR_ART, Question

ORGANISM = "Mycobacterium tuberculosis"
DEFAULT_LIMIT = 10

CORPUS = "corpus"
PUBMED = "pubmed"


@dataclass
class Candidate:
    """One retrieved record, with enough to trace it and nothing interpreted."""

    question_id: str
    kind: str
    source: str                 # corpus | pubmed
    record_id: str              # PMID
    title: str = ""
    year: str = ""
    rank: int = 0
    score: float = 0.0
    query: str = ""

    def to_dict(self) -> dict[str, Any]:
        return {"question_id": self.question_id, "kind": self.kind, "source": self.source,
                "record_id": self.record_id, "title": self.title, "year": self.year,
                "rank": self.rank, "score": f"{self.score:.4f}" if self.score else "NA",
                "query": self.query}


def aliases_for(locus: str, annotation: Annotation) -> list[str]:
    """The spellings a paper might use, cheapest first.

    Symbol and locus tag come from the annotation at no cost. A richer alias set
    is available from `identity.resolve`, at one cached lookup per locus; it is
    deliberately not used here so stage B stays offline-capable and free.
    """
    gene = annotation.gene(locus)
    out = [locus]
    if gene is not None and gene.symbol and gene.symbol != locus:
        out.insert(0, gene.symbol)
    return out


def _any_of(names: list[str]) -> str:
    quoted = " OR ".join(f'"{n}"' for n in names)
    return f"({quoted})" if len(names) > 1 else quoted


def build_query(question: Question, annotation: Annotation) -> str:
    """A PubMed query for one question. Loci and drugs only.

    An intergenic feature is queried by its flanking genes: no paper names a gap,
    and the promoter of a gene is discussed under that gene's name.
    """
    def names(locus: str) -> list[str]:
        interval = annotation.intergenic(locus)
        if interval is not None:
            return (aliases_for(interval.left.locus, annotation)
                    + aliases_for(interval.right.locus, annotation))
        return aliases_for(locus, annotation)

    parts = [_any_of(names(question.locus))]
    if question.partner:
        parts.append(_any_of(names(question.partner)))
    # Two loci are specific enough on their own; adding a drug would narrow a
    # pair query to the papers that already frame them as a resistance story,
    # which is the reading being tested. A single-locus question has no such
    # constraint, so the drug supplies the context instead.
    if question.drugs and question.kind not in (COMPENSATION, PRIOR_ART):
        parts.append(_any_of(list(question.drugs)))
    parts.append(f'"{ORGANISM}"')
    return " AND ".join(parts)


@dataclass
class Report:
    candidates: list[Candidate] = field(default_factory=list)
    answered: set[str] = field(default_factory=set)      # question ids with >=1 candidate
    asked: list[str] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)

    def summary(self) -> dict[str, Any]:
        by_source: dict[str, int] = {}
        for candidate in self.candidates:
            by_source[candidate.source] = by_source.get(candidate.source, 0) + 1
        return {"questions": len(self.asked),
                "with_a_candidate": len(self.answered),
                "without": len(self.asked) - len(self.answered),
                "records": len(self.candidates),
                "distinct_records": len({c.record_id for c in self.candidates}),
                "by_source": by_source}


def _records(result: dict[str, Any]) -> list[dict[str, Any]]:
    return result.get("records", []) if isinstance(result, dict) else []


def retrieve(questions: list[Question], annotation: Annotation, clients: Any,
             store: Any = None, limit: int = DEFAULT_LIMIT) -> Report:
    """Fetch candidates for every question. Fail-soft, one question at a time."""
    report = Report()
    if not os.environ.get("NCBI_EMAIL", "").strip():
        # One unidentified lookup is unremarkable; a few hundred in a batch is
        # where NCBI starts noticing, and the failure mode is a throttle with no
        # local signal.
        report.notes.append(
            "NCBI_EMAIL is not set. NCBI asks callers to identify themselves and contacts "
            "you before blocking, so a batch of this size should set it.")

    corpus = getattr(clients, "corpus", None)
    pubmed = getattr(clients, "pubmed", None)

    for question in questions:
        report.asked.append(question.id)
        query = build_query(question, annotation)
        found: list[Candidate] = []

        if corpus is not None:
            try:
                result = corpus(question.text, limit=limit)
                if store is not None:
                    store.tool_result("corpus_search", {"query": question.text}, result)
                for rank, record in enumerate(_records(result), start=1):
                    detail = record.get("detail", {})
                    found.append(Candidate(
                        question.id, question.kind, CORPUS, record["record_id"],
                        title=detail.get("title", ""), year=str(detail.get("year", "")),
                        rank=rank, score=float(detail.get("score", 0.0)), query=question.text))
            except (FetchError, ValueError, KeyError) as exc:
                report.notes.append(f"{question.id}: corpus unavailable ({exc})")

        # PubMed only where the corpus reached nothing. An empty corpus result is
        # a statement about the corpus, not about the literature.
        if not found and pubmed is not None:
            try:
                arguments = {"term": query, "limit": limit}
                result = pubmed(query, organism=ORGANISM, limit=limit)
                if store is not None:
                    store.tool_result("pubmed_abstracts", arguments, result)
                for rank, record in enumerate(_records(result), start=1):
                    detail = record.get("detail", {})
                    found.append(Candidate(
                        question.id, question.kind, PUBMED, record["record_id"],
                        title=detail.get("title", ""), year=str(detail.get("year", "")),
                        rank=rank, query=query))
            except (FetchError, ValueError, KeyError) as exc:
                report.notes.append(f"{question.id}: PubMed unavailable ({exc})")

        report.candidates.extend(found)
        if found:
            report.answered.add(question.id)
    return report


@dataclass
class Clients:
    """The two retrievers, as plain callables returning a ToolResult payload.

    `corpus` is None when no corpus is configured, which is the ordinary case for
    a rule set naming loci outside a prebuilt one -- and is reported as a note
    rather than left to look like an empty corpus.
    """

    corpus: Any = None
    pubmed: Any = None

    @classmethod
    def live(cls, http: Any, corpus_path: str | None = None) -> Clients:
        from kegg_string_mcp.corpus_search import CorpusSearchClient
        from kegg_string_mcp.pubmed import PubMedClient

        search = CorpusSearchClient(corpus_path)
        corpus = ((lambda q, limit=DEFAULT_LIMIT: search.search(q, limit=limit).model_dump())
                  if search.configured else None)
        client = PubMedClient(http)
        # `search_abstracts`, not `abstracts`: the latter refuses query syntax and
        # phrase-quotes its argument, so a composed boolean query never runs. That
        # guard is right for the model-facing tool and wrong for pipeline code.
        return cls(corpus=corpus,
                   pubmed=lambda q, organism=ORGANISM, limit=DEFAULT_LIMIT:
                   client.search_abstracts(q, limit=limit).model_dump())
