"""Two retrieval arms over one corpus: BM25 keyword, and dense vector search.

They deliberately share an interface and a corpus. The comparison is only
meaningful if the *only* thing that differs between them is how a passage is
scored -- same text, same fetch, same provenance.

Embeddings come from Chroma's bundled ONNX MiniLM rather than
sentence-transformers. That is not a compromise: it removes torch from the
dependency set entirely, which turns a ~2 GB install that will not build on every
platform into an 80 MB model download that runs anywhere, CI included.

Every hit carries its PMID, so a passage surfaced by vector search is exactly as
checkable as one surfaced by keyword search. The validator already in the repo
works on either without modification -- which is the point of adding the second
arm here rather than in a separate project.
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass, field
from typing import Any, Protocol

from kegg_string_mcp.retrieval.corpus import Corpus, Passage

DEFAULT_K = 10


def corpus_fingerprint(corpus: Corpus) -> str:
    """Short content hash of the corpus: its passage ids and their text.

    Two corpora with the same passages share an index; any difference gets a fresh
    one. Metadata is excluded deliberately -- re-running the corpus build changes
    `queried_for` ordering without changing what is searchable.
    """
    digest = hashlib.sha256()
    for passage in sorted(corpus.passages, key=lambda p: p.passage_id):
        digest.update(passage.passage_id.encode())
        digest.update(b"\x00")
        digest.update(passage.text.encode())
        digest.update(b"\x00")
    return digest.hexdigest()[:12]


@dataclass
class Hit:
    """One retrieved passage, with the provenance needed to validate a claim from it."""

    passage_id: str
    pmid: str
    score: float
    rank: int
    title: str
    text: str
    mentions: list[str] = field(default_factory=list)

    @property
    def url(self) -> str:
        return f"https://pubmed.ncbi.nlm.nih.gov/{self.pmid}/"

    def to_dict(self) -> dict[str, Any]:
        return {"passage_id": self.passage_id, "pmid": self.pmid, "score": round(self.score, 4),
                "rank": self.rank, "title": self.title, "mentions": self.mentions,
                "url": self.url}


class Retriever(Protocol):
    name: str

    def search(self, query: str, k: int = DEFAULT_K) -> list[Hit]: ...


class KeywordIndex:
    """BM25 over the same passages. The baseline, and the arm already implicit in
    the repo: PubMed's own search is lexical, so this makes that comparable."""

    name = "keyword"

    def __init__(self, corpus: Corpus):
        from rank_bm25 import BM25Okapi

        self.passages: list[Passage] = list(corpus.passages)
        self._tokens = [self._tokenise(p.text) for p in self.passages]
        self._bm25 = BM25Okapi(self._tokens)

    @staticmethod
    def _tokenise(text: str) -> list[str]:
        # Lowercase word characters, keeping hyphens and digits so gene symbols and
        # identifiers survive: "Rv1908c" and "beta-lactamase" must stay single tokens.
        import re

        return re.findall(r"[a-z0-9][a-z0-9\-]*", text.lower())

    def search(self, query: str, k: int = DEFAULT_K) -> list[Hit]:
        scores = self._bm25.get_scores(self._tokenise(query))
        order = sorted(range(len(scores)), key=lambda i: -scores[i])[:k]
        return [
            Hit(passage_id=self.passages[i].passage_id, pmid=self.passages[i].pmid,
                score=float(scores[i]), rank=rank, title=self.passages[i].title,
                text=self.passages[i].text, mentions=list(self.passages[i].mentions))
            for rank, i in enumerate(order, start=1)
        ]


class VectorIndex:
    """Dense retrieval over the same passages, via Chroma."""

    name = "vector"

    def __init__(self, corpus: Corpus, collection: str = "tb-literature",
                 path: str | None = None):
        import chromadb
        from chromadb.utils import embedding_functions

        self.passages = {p.passage_id: p for p in corpus.passages}
        client = (chromadb.PersistentClient(path=path) if path
                  else chromadb.EphemeralClient())
        self._embed = embedding_functions.ONNXMiniLM_L6_V2()
        # The collection name carries a fingerprint of the corpus. Keyed on name
        # alone, a second corpus silently reused the first one's vectors -- the
        # load was skipped because count() was non-zero -- and every retrieval
        # then answered from the wrong text. That is fatal for a comparison
        # harness: the numbers would be wrong without anything failing.
        # Same instinct as the HTTP cache elsewhere in the repo: address by content.
        name = f"{collection}-{corpus_fingerprint(corpus)}"
        # Cosine, not the default L2: these embeddings are direction-normalised and
        # cosine is what the model was trained against.
        self._collection = client.get_or_create_collection(
            name, metadata={"hnsw:space": "cosine"})

        if self._collection.count() != len(corpus.passages) and corpus.passages:
            self._collection.add(
                ids=[p.passage_id for p in corpus.passages],
                documents=[p.text for p in corpus.passages],
                embeddings=self._embed([p.text for p in corpus.passages]),
                metadatas=[{"pmid": p.pmid, "title": p.title[:300]} for p in corpus.passages],
            )

    def search(self, query: str, k: int = DEFAULT_K) -> list[Hit]:
        if not self.passages:
            return []
        result = self._collection.query(
            query_embeddings=self._embed([query]),
            n_results=min(k, max(self._collection.count(), 1)))
        hits: list[Hit] = []
        for rank, (pid, dist) in enumerate(
                zip(result["ids"][0], result["distances"][0]), start=1):
            passage = self.passages.get(pid)
            if passage is None:
                # Belt and braces behind the fingerprint: an id the corpus does not
                # contain means the collection is not this corpus, and returning it
                # would attribute someone else's text to this run.
                raise RuntimeError(
                    f"vector collection returned id {pid!r}, which is not in this corpus -- "
                    f"the index does not match the corpus it was opened with")
            # Chroma returns cosine *distance*; report similarity so both arms
            # score in the same direction (higher is better).
            hits.append(Hit(passage_id=pid, pmid=passage.pmid, score=1.0 - float(dist),
                            rank=rank, title=passage.title, text=passage.text,
                            mentions=list(passage.mentions)))
        return hits


def reciprocal_rank_fusion(runs: list[list[Hit]], k: int = DEFAULT_K,
                           smoothing: int = 60) -> list[Hit]:
    """Merge ranked lists by rank rather than score.

    Scores from BM25 and cosine similarity are not comparable -- different scales,
    different distributions -- so fusing on rank is the standard answer. The
    smoothing constant damps the influence of the top position; 60 is the value
    from the original RRF paper and is what most implementations use.
    """
    fused: dict[str, float] = {}
    best: dict[str, Hit] = {}
    for run in runs:
        for hit in run:
            fused[hit.passage_id] = fused.get(hit.passage_id, 0.0) + 1.0 / (smoothing + hit.rank)
            if hit.passage_id not in best or hit.rank < best[hit.passage_id].rank:
                best[hit.passage_id] = hit
    order = sorted(fused, key=lambda pid: -fused[pid])[:k]
    return [
        Hit(passage_id=pid, pmid=best[pid].pmid, score=fused[pid], rank=rank,
            title=best[pid].title, text=best[pid].text, mentions=best[pid].mentions)
        for rank, pid in enumerate(order, start=1)
    ]


class HybridIndex:
    """Both arms, fused by reciprocal rank. What production systems actually use."""

    name = "hybrid"

    def __init__(self, keyword: KeywordIndex, vector: VectorIndex):
        self.keyword, self.vector = keyword, vector

    def search(self, query: str, k: int = DEFAULT_K) -> list[Hit]:
        return reciprocal_rank_fusion(
            [self.keyword.search(query, k * 2), self.vector.search(query, k * 2)], k=k)
