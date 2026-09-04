"""Build a literature corpus from the existing PubMed tool.

The corpus is assembled with the tool the repo already has, not a new fetcher.
That matters for the comparison: both retrieval arms then search *the same text*,
retrieved the same way, so any difference between them is the retrieval method
rather than a difference in what was fetched.

Every chunk keeps the PMID it came from, so a passage surfaced by vector search is
as checkable as one surfaced by keyword search -- the existing validator works on
either without modification.
"""

from __future__ import annotations

import json
import re
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

from kegg_string_mcp.cache import DiskCache
from kegg_string_mcp.http import PoliteClient
from kegg_string_mcp.pubmed import PubMedClient


@dataclass
class Passage:
    """One retrievable unit, carrying enough provenance to be validated.

    `text` is a span of the record's `quotable_text`, so a quote drawn from a
    passage still contains-matches the source the validator checks against.
    """

    passage_id: str
    pmid: str
    text: str
    title: str
    year: str
    journal: str
    doi: str
    pmcid: str
    in_pmc: bool
    mentions: list[str] = field(default_factory=list)
    queried_for: list[str] = field(default_factory=list)
    # Genes of the WHOLE corpus that this passage's text names. Distinct from
    # `mentions`, which PubMed's tool computes against the terms of the query that
    # fetched the paper -- so a paper found by searching katG can never record
    # ahpC there, however prominently it discusses it. Relevance scoring needs
    # this field; using `mentions` measured corpus construction instead of
    # retrieval, and undercounted pair evidence 25 vs 88 on the TB corpus.
    genes_named: list[str] = field(default_factory=list)
    # Set when a passage is a chunk of a longer abstract. Paper-level metrics
    # dedupe on `pmid`, so a paper split into three chunks still counts once.
    chunk_index: int = 0
    chunk_of: int = 1

    @property
    def url(self) -> str:
        return f"https://pubmed.ncbi.nlm.nih.gov/{self.pmid}/"

    def to_dict(self) -> dict[str, Any]:
        return asdict(self) | {"url": self.url}


@dataclass
class Corpus:
    genes: list[str]
    passages: list[Passage] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)

    @property
    def pmids(self) -> set[str]:
        return {p.pmid for p in self.passages}

    def write(self, path: Path) -> Path:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(
            {"genes": self.genes, "notes": self.notes,
             "passages": [p.to_dict() for p in self.passages]}, indent=1), encoding="utf-8")
        return path

    @classmethod
    def read(cls, path: Path) -> Corpus:
        raw = json.loads(path.read_text(encoding="utf-8"))
        return cls(genes=raw["genes"], notes=raw.get("notes", []),
                   passages=[Passage(**{k: v for k, v in p.items() if k != "url"})
                             for p in raw["passages"]])


def build(genes: list[str], organism: str = "Mycobacterium tuberculosis",
          limit: int = 20, client: PubMedClient | None = None) -> Corpus:
    """Fetch abstracts for each gene and flatten them into deduplicated passages.

    Deduplication is by PMID: the same paper is routinely returned for several
    genes, and counting it once per gene would inflate any retrieval metric that
    depends on corpus size.
    """
    client = client or PubMedClient(PoliteClient(DiskCache()))
    corpus = Corpus(genes=list(genes))
    by_pmid: dict[str, Passage] = {}

    for gene in genes:
        result = client.abstracts(gene, organism=organism, limit=limit)
        if not result.records:
            corpus.notes.append(f"{gene}: no abstracts retrieved. "
                                + (result.notes[0][:160] if result.notes else ""))
            continue
        for record in result.records:
            detail = record.detail
            existing = by_pmid.get(record.record_id)
            if existing is not None:
                # Same paper, different query. Keep one copy and record both.
                if gene not in existing.queried_for:
                    existing.queried_for.append(gene)
                for mention in detail.get("mentions", []):
                    if mention not in existing.mentions:
                        existing.mentions.append(mention)
                continue
            by_pmid[record.record_id] = Passage(
                passage_id=record.record_id, pmid=record.record_id,
                text=detail.get("quotable_text", ""), title=detail.get("title", ""),
                year=detail.get("year", ""), journal=detail.get("journal", ""),
                doi=detail.get("doi", ""), pmcid=detail.get("pmcid", ""),
                in_pmc=bool(detail.get("in_pmc")),
                mentions=list(detail.get("mentions", [])), queried_for=[gene],
            )

    corpus.passages = [by_pmid[k] for k in sorted(by_pmid)]
    annotate_genes_named(corpus)
    return corpus


# MiniLM truncates at 256 wordpieces. ~180 English words sits inside that with
# room for subword splitting, which is heavy on gene symbols and Latin binomials.
CHUNK_WORDS = 180
CHUNK_OVERLAP = 40


def chunk(corpus: Corpus, max_words: int = CHUNK_WORDS,
          overlap: int = CHUNK_OVERLAP) -> Corpus:
    """Split long passages so both arms index the same units.

    Without this the comparison is confounded rather than measured: the embedding
    model truncates at its context window while BM25 indexes every word, so on
    this corpus 80% of abstracts were only partly visible to the dense arm and in
    18 of them a gene appeared solely past the cut. Any dense-vs-lexical gap then
    conflates "worse retrieval" with "read less of the document".

    Overlapping windows because a relationship stated across a boundary would
    otherwise be split in half and lost to both arms.
    """
    out = Corpus(genes=list(corpus.genes), notes=list(corpus.notes))
    for passage in corpus.passages:
        words = passage.text.split()
        if len(words) <= max_words:
            out.passages.append(passage)
            continue
        step = max(max_words - overlap, 1)
        starts = list(range(0, len(words), step))
        # Drop a trailing window that the previous one already covers entirely.
        starts = [i for i in starts if i == 0 or i + max_words > starts[-1] or len(words) - i > overlap]
        pieces = [" ".join(words[i:i + max_words]) for i in starts]
        pieces = [p for p in pieces if p.strip()]
        for index, text in enumerate(pieces):
            out.passages.append(Passage(
                passage_id=f"{passage.passage_id}#{index}", pmid=passage.pmid, text=text,
                title=passage.title, year=passage.year, journal=passage.journal,
                doi=passage.doi, pmcid=passage.pmcid, in_pmc=passage.in_pmc,
                mentions=list(passage.mentions), queried_for=list(passage.queried_for),
                chunk_index=index, chunk_of=len(pieces)))
    return annotate_genes_named(out)


def annotate_genes_named(corpus: Corpus) -> Corpus:
    """Record, per passage, which of the corpus's genes its text actually names.

    Word-boundary matching so `embB` does not match inside another token. Done
    once over the finished corpus rather than per query, which is the whole point:
    the answer must not depend on which gene's search happened to return the paper.
    """
    patterns = {gene: re.compile(rf"\b{re.escape(gene)}\b", re.IGNORECASE)
                for gene in corpus.genes}
    for passage in corpus.passages:
        passage.genes_named = [g for g, pattern in patterns.items()
                               if pattern.search(passage.text)]
    return corpus
