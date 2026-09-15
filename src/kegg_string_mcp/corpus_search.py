"""Search a locally built literature corpus, using the measured retrieval arm.

`pubmed_abstracts` asks NCBI a keyword query and works for any gene. This asks a
corpus already assembled for a gene set, and ranks with the hybrid arm that
measured best on it: 0.917 precision@10 against 0.844 for lexical alone and
0.739 for dense alone (`docs/RETRIEVAL.md`). The two are complements rather than
alternatives -- live search reaches anything, corpus search ranks better over
what has been gathered -- so this tool says plainly when a query falls outside
the corpus, and the agent falls back.

**Records are keyed by PMID, and `quotable_text` is the whole abstract.** What
retrieval ranks is a *passage*: a 180-word chunk with 40 words of overlap, so a
long abstract occupies several. Returning the matched chunk as the quotable text
would be a trap. The run store keeps the first record it sees for an ID and
merges only `mentions`, so a chunk arriving before the same PMID's full abstract
would become the text every later quote is checked against -- and a correct quote
from elsewhere in that abstract would be reported as fabricated. The chunks of a
PMID are therefore rejoined in order. They overlap, so the join repeats a span;
that is harmless, because the check is containment and every span of the original
abstract still appears contiguously within it.

The corpus is not committed, so this tool is inert until one is built. Set
`KEGG_STRING_MCP_CORPUS` to a corpus JSON written by `scripts/build_corpus.py`.
"""

from __future__ import annotations

import hashlib
import os
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from kegg_string_mcp.provenance import Record, RequestTrace, ToolResult

CORPUS_ENV = "KEGG_STRING_MCP_CORPUS"
DEFAULT_LIMIT = 5
PUBMED_URL = "https://pubmed.ncbi.nlm.nih.gov/"

_NO_CORPUS = (
    "No local corpus is configured, so this tool has nothing to search. Set "
    f"{CORPUS_ENV} to a corpus built by scripts/build_corpus.py, or use "
    "pubmed_abstracts, which queries NCBI live and covers any gene. This is not "
    "evidence that no literature exists."
)


class CorpusSearchClient:
    """Hybrid retrieval over a prebuilt corpus.

    The index is built on first use and reused. Embedding 1,350 passages takes
    seconds and is pure overhead for a session that never searches, so nothing
    happens until a query arrives -- which also keeps the optional `vector`
    dependencies out of the import path of a server that has no corpus.
    """

    def __init__(self, corpus_path: str | Path | None = None):
        self._configured = corpus_path or os.environ.get(CORPUS_ENV) or None
        # expanduser: the value arrives from an environment variable, where a
        # leading ~ is a literal character and not the shell's home directory.
        self._path = (Path(self._configured).expanduser()
                      if self._configured else None)
        self._index: Any = None
        self._corpus: Any = None
        self._fingerprint = ""
        self._built_at = ""
        self._error: str | None = None

    @property
    def configured(self) -> bool:
        return self._path is not None

    def _build(self) -> bool:
        """True once an index is ready. Failures are recorded, never raised."""
        if self._index is not None or self._error:
            return self._index is not None
        if self._path is None:
            self._error = _NO_CORPUS
            return False
        if not self._path.exists():
            self._error = (f"The configured corpus {str(self._path)!r} does not exist. "
                           f"Build one with scripts/build_corpus.py, or use "
                           f"pubmed_abstracts. No search was performed.")
            return False
        try:
            from kegg_string_mcp.retrieval.corpus import Corpus
            from kegg_string_mcp.retrieval.index import HybridIndex, KeywordIndex, VectorIndex
        except ImportError as exc:
            self._error = (f"Corpus search needs the optional retrieval extra "
                           f"(pip install -e '.[vector]'): {exc}. Use pubmed_abstracts "
                           f"instead. No search was performed.")
            return False
        self._corpus = Corpus.read(self._path)
        # The corpus is a file, not a fetch. Its provenance is the content hash
        # and the file's own mtime -- there is no retrieval timestamp to borrow,
        # and inventing one would put a fabricated time in the run store.
        # A real digest of the bytes read, because the field it fills is named
        # content_sha256 and every other tool puts one there. retrieval's
        # corpus_fingerprint is a truncated hash for naming a vector collection;
        # putting that here would quietly break what the field asserts.
        self._fingerprint = hashlib.sha256(self._path.read_bytes()).hexdigest()
        self._built_at = (datetime.fromtimestamp(self._path.stat().st_mtime, tz=timezone.utc)
                          .isoformat())
        keyword = KeywordIndex(self._corpus)
        self._index = HybridIndex(keyword, VectorIndex(self._corpus))
        return True

    def search(self, query: str, limit: int = DEFAULT_LIMIT) -> ToolResult:
        request: dict[str, Any] = {"query": query, "limit": limit}
        query = (query or "").strip()
        if not query:
            return ToolResult.build(request, [], resolved={"matched_by": "none"},
                                    notes=["No query was supplied."])
        if not self._build():
            return ToolResult.build(request, [], resolved={"matched_by": "none"},
                                    notes=[self._error or _NO_CORPUS])

        # Over-fetch: several passages can belong to one paper, and the caller
        # asked for papers.
        hits = self._index.search(query, k=max(limit * 3, limit))
        by_pmid: dict[str, list[Any]] = {}
        for hit in hits:
            by_pmid.setdefault(hit.pmid, []).append(hit)

        records = []
        for pmid, group in list(by_pmid.items())[:limit]:
            passages = [p for p in self._corpus.passages if p.pmid == pmid]
            passages.sort(key=lambda p: p.chunk_index)
            best = min(group, key=lambda h: h.rank)
            records.append(Record(
                record_id=pmid,
                type="article",
                name=best.title,
                url=f"{PUBMED_URL}{pmid}/",
                source="pubmed",
                retrieved_at=self._built_at,
                cached=True,
                detail={
                    "quotable_text": " ".join(p.text for p in passages),
                    "matched_passage": best.text,
                    "rank": best.rank,
                    "score": round(best.score, 4),
                    "genes_named": list(best.genes_named),
                    "mentions": list(best.mentions),
                    "passages_in_corpus": len(passages),
                    # The corpus carries what the manifest needs for a downstream
                    # full-text fetch, and omitting it here did not leave those
                    # fields empty -- it left `in_pmc: false` on every paper of a
                    # corpus-only run, which reads as "not in PMC" rather than
                    # "not recorded" and is the field that decides whether full
                    # text can be fetched at all.
                    "title": best.title,
                    "year": passages[0].year,
                    "journal": passages[0].journal,
                    "doi": passages[0].doi,
                    "pmcid": passages[0].pmcid,
                    "in_pmc": bool(passages[0].in_pmc),
                },
            ))

        notes = [
            ("Ranked by the hybrid arm (BM25 fused with dense embeddings), which measured "
             "0.917 precision@10 on this corpus against 0.844 lexical and 0.739 dense. "
             "`quotable_text` is the whole abstract; `matched_passage` is the span that "
             "ranked."),
            (f"This corpus holds {len(self._corpus.passages)} passages from "
             f"{len(self._corpus.pmids)} papers, gathered for: "
             f"{', '.join(self._corpus.genes[:12])}"
             f"{' ...' if len(self._corpus.genes) > 12 else ''}. A gene outside that set "
             f"will retrieve nothing relevant however it is queried -- use "
             f"pubmed_abstracts for those."),
        ]
        if not records:
            notes.append("Nothing in the corpus matched. This says the corpus does not "
                         "cover the query, not that no literature exists.")
        traces = [RequestTrace(url=f"file://{self._path}", retrieved_at=self._built_at,
                               cached=True, status=200,
                               content_sha256=self._fingerprint)]
        return ToolResult.build(request, records, resolved={"matched_by": "corpus_search"},
                                requests=traces, notes=notes)
