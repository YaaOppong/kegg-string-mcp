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
    return corpus
