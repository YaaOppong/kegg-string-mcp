"""Measure the retrieval arms against each other on one corpus and one query set.

The comparison is the contribution. Adding vector search to a repo is a tutorial;
measuring where it beats keyword search, where it loses, and by how much, is a
result.

**Relevance is judged deterministically**, without hand-labelling and without a
model. Each retrieved paper carries `mentions`: the genes its own text names, as
opposed to the query that happened to find it. So "did this retrieval return
papers that actually discuss the gene asked about?" has an exact answer. That is
weaker than human relevance judgements -- a paper can name a gene in passing --
but it is reproducible, costs nothing, and cannot be tuned after the fact.

Three things are measured:

* **overlap** -- how much the arms agree, which says whether the second arm is
  earning its place at all;
* **on-target precision@k** -- of what each returned, how much names the gene;
* **exact-term retrieval** -- querying a bare gene symbol or locus tag. This is
  where dense retrieval is expected to fail, and finding it in your own data is
  worth more than citing that it happens.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from itertools import combinations
from pathlib import Path
from typing import Any

from kegg_string_mcp.retrieval.index import DEFAULT_K, Hit


def on_target(hits: list[Hit], genes: list[str]) -> list[Hit]:
    """Hits whose own text names at least one of the queried genes."""
    wanted = {g.lower() for g in genes}
    return [h for h in hits if wanted & {m.lower() for m in h.genes_named}]


def naming_all(hits: list[Hit], genes: list[str]) -> list[Hit]:
    """Passages naming every gene of the pair, deduplicated to one per paper.

    Chunking splits an abstract into overlapping windows, so a paper can appear
    several times in one result list. Counting each chunk would inflate the
    pair-evidence metric by the chunking parameters rather than by retrieval.
    """
    wanted = {g.lower() for g in genes}
    matched = [h for h in hits if wanted <= {m.lower() for m in h.genes_named}]
    seen: set[str] = set()
    unique: list[Hit] = []
    for hit in matched:
        if hit.pmid not in seen:
            seen.add(hit.pmid)
            unique.append(hit)
    return unique


@dataclass
class QueryResult:
    query: str
    genes: list[str]
    per_arm: dict[str, list[str]] = field(default_factory=dict)      # arm -> pmids
    precision: dict[str, float] = field(default_factory=dict)
    joint: dict[str, int] = field(default_factory=dict)              # papers naming ALL genes

    def overlap(self, a: str, b: str) -> float:
        """Jaccard over returned PMIDs. 1.0 means the arms are interchangeable here."""
        sa, sb = set(self.per_arm.get(a, [])), set(self.per_arm.get(b, []))
        union = sa | sb
        return round(len(sa & sb) / len(union), 3) if union else 0.0

    def only_in(self, a: str, b: str) -> list[str]:
        return sorted(set(self.per_arm.get(a, [])) - set(self.per_arm.get(b, [])))


@dataclass
class Comparison:
    k: int
    results: list[QueryResult] = field(default_factory=list)
    exact_term: list[dict[str, Any]] = field(default_factory=list)

    def mean_overlap(self, a: str, b: str) -> float:
        vals = [r.overlap(a, b) for r in self.results]
        return round(sum(vals) / len(vals), 3) if vals else 0.0

    def mean_precision(self, arm: str) -> float:
        vals = [r.precision[arm] for r in self.results if arm in r.precision]
        return round(sum(vals) / len(vals), 3) if vals else 0.0

    def mean_joint(self, arm: str) -> float:
        vals = [r.joint[arm] for r in self.results if arm in r.joint]
        return round(sum(vals) / len(vals), 2) if vals else 0.0

    def to_dict(self) -> dict[str, Any]:
        arms = sorted({a for r in self.results for a in r.per_arm})
        return {
            "k": self.k,
            "queries": len(self.results),
            "arms": arms,
            "mean_on_target_precision": {a: self.mean_precision(a) for a in arms},
            "mean_papers_naming_both": {a: self.mean_joint(a) for a in arms},
            "mean_overlap": {f"{a}|{b}": self.mean_overlap(a, b)
                             for a, b in combinations(arms, 2)},
            "exact_term": self.exact_term,
            "per_query": [asdict(r) for r in self.results],
        }

    def write(self, path: Path) -> Path:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(self.to_dict(), indent=1), encoding="utf-8")
        return path


def pair_queries(genes: list[str]) -> list[tuple[str, list[str]]]:
    return [(f"What is the relationship between {a} and {b} in Mycobacterium tuberculosis?",
             [a, b]) for a, b in combinations(genes, 2)]


def compare(arms: dict[str, Any], queries: list[tuple[str, list[str]]],
            k: int = DEFAULT_K) -> Comparison:
    out = Comparison(k=k)
    for query, genes in queries:
        result = QueryResult(query=query, genes=genes)
        for name, arm in arms.items():
            hits = arm.search(query, k=k)
            result.per_arm[name] = [h.pmid for h in hits]
            result.precision[name] = round(len(on_target(hits, genes)) / max(len(hits), 1), 3)
            result.joint[name] = len(naming_all(hits, genes))
        out.results.append(result)
    return out


def exact_term_probe(arms: dict[str, Any], terms: list[str], k: int = DEFAULT_K,
                     corpus: Any | None = None) -> list[dict]:
    """Query a bare identifier and ask whether the top-k passages contain it.

    The known weakness of dense retrieval: a gene symbol or locus tag has no
    useful neighbourhood in embedding space, while BM25 matches it exactly. This
    is the measurement that argues for hybrid search rather than a swap.
    """
    import re

    rows = []
    for term in terms:
        # Scored against the passage TEXT, not the corpus gene list: a locus tag
        # like Rv1908c is never a corpus gene, so scoring it the same way as a
        # gene-pair query guaranteed a zero for every arm and made the probe read
        # as "BM25 fails at exact terms too" -- the opposite of what it measures.
        pattern = re.compile(rf"\b{re.escape(term)}\b", re.IGNORECASE)
        row: dict[str, Any] = {"term": term}
        for name, arm in arms.items():
            hits = arm.search(term, k=k)
            containing = [h for h in hits if pattern.search(h.text)]
            row[name] = {"hits": len(hits), "containing_the_term": len(containing),
                         "top_pmid": hits[0].pmid if hits else None}
        # How many passages in the whole corpus contain it at all. Without this a
        # zero is ambiguous between "retrieval missed it" and "it is not there" --
        # and for a locus tag it is usually the latter, which is a finding about
        # the corpus rather than about any retriever.
        row["in_corpus"] = (sum(1 for psg in corpus.passages if pattern.search(psg.text))
                            if corpus is not None else None)
        rows.append(row)
    return rows
