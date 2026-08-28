"""PubMed lookups via NCBI E-utilities: gene -> abstracts.

This tool is epistemically weaker than the other two, and most of the code below
exists to say so precisely rather than to hide it.

KEGG and STRING **resolve** an identifier: there is a right answer, the lookup
either finds it or does not, and a miss is reportable as a miss. PubMed does not.
`esearch` is a relevance-ranked free-text search that essentially always returns
something plausible, so there is no "did not resolve" signal to report -- a query
for a gene symbol will happily return papers that merely contain the string. The
tool therefore returns PubMed's own `querytranslation` (what it actually searched
for, after Automatic Term Mapping) and the total hit count alongside the
truncated result set, so a reader can see how much was left behind.

The second difference matters more. A KEGG pathway ID and a STRING score are
*structured*: the record means one thing, and the model's only job is to relay
it. An abstract is prose, and any claim drawn from it is the model deciding what
a record says -- which is precisely what this server's other tools are built to
prevent. Set membership on `record_ids` cannot catch that: a real PMID attached
to a fabricated finding passes it.

So every record carries `quotable_text`: the exact, whitespace-normalised text
that was retrieved for it. A claim sourced from an abstract should quote a span
of that string, and the quote can then be checked by containment -- deterministic,
no embeddings, and it upgrades validation from "this record was retrieved" to
"this record was retrieved AND the quoted span really is in it". Building that
check is the agent layer's job; storing the exact text to check against is this
module's.

Text is flattened from PubMed's XML, which marks up titles and abstracts with
`<i>` and `<sup>`. Flattening drops the markup and normalises whitespace, so
"10<sup>5</sup>" becomes "105". `quotable_text` is the flattened form, and a
validator must normalise a candidate quote the same way before comparing.

NCBI asks callers to identify themselves with `tool` and `email` on every
request, and rate-limits by IP. Set NCBI_EMAIL; set NCBI_API_KEY if you have one.
"""

from __future__ import annotations

import json
import os
import re
import xml.etree.ElementTree as ET
from typing import Any

from kegg_string_mcp.http import FetchError, PoliteClient
from kegg_string_mcp.provenance import Record, RequestTrace, ToolResult

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
ARTICLE = "https://pubmed.ncbi.nlm.nih.gov/"

MTB_H37RV_NAME = "Mycobacterium tuberculosis"

DEFAULT_LIMIT = 10
MAX_LIMIT = 100

# PubMed's term grammar. The gene and organism are interpolated into a query
# string, so a value containing these would change what is searched for rather
# than what is searched -- `katG" AND cancer OR "` is a different question than
# the caller asked. Quoting alone does not contain it, so reject instead: the
# codebase's rule is that a bad argument is an explicit error, never a silently
# empty result that reads as "no literature exists".
_QUERY_SYNTAX = re.compile(r'["\[\]]')
_PMID = re.compile(r"^\d+$")

METADATA_MATCH_CAVEAT = (
    "PubMed searched [All Fields], which spans MeSH terms, author keywords and substance "
    "lists as well as the title and abstract -- and it does NOT search full text. So a "
    "record can match without the gene appearing in the text retrieved here, and a paper "
    "that discusses the gene only in its Results is not findable this way."
)

SEARCH_CAVEAT = (
    "PubMed search is relevance-ranked text retrieval, NOT identifier resolution. These "
    "articles matched the query string; matching is not the same as being about this gene, "
    "and unlike KEGG or STRING there is no 'did not resolve' signal to distinguish the two. "
    "Read each record's title and abstract before citing it for a claim about the gene."
)

TEXTMINING_CAVEAT = (
    "If STRING scored an interaction for this gene via its textmining channel, these "
    "abstracts may be the very literature that produced that score. Citing a STRING "
    "textmining-supported interaction and a PubMed abstract about the same pair is one "
    "line of evidence counted twice, not independent corroboration."
)


def ncbi_identity() -> dict[str, str]:
    """NCBI's equivalent of STRING's caller_identity, plus an optional key.

    `email` and `api_key` are stripped from the cache key by IDENTITY_PARAMS, and
    `api_key` is additionally scrubbed from the audit URL by SECRET_PARAMS -- a
    credential must not reach the run store on disk.
    """
    params = {"tool": os.environ.get("NCBI_TOOL", "kegg-string-mcp")}
    for variable, param in (("NCBI_EMAIL", "email"), ("NCBI_API_KEY", "api_key")):
        value = os.environ.get(variable, "").strip()
        if value:
            params[param] = value
    return params


def _trace(resp) -> RequestTrace:
    return RequestTrace(url=resp.audit_url, retrieved_at=resp.fetched_at, cached=resp.cached,
                        status=resp.status, content_sha256=resp.content_sha256)


def _text(element: ET.Element | None) -> str:
    """Flatten an element's text, markup children included, to normalised text.

    `itertext()` rather than `.text`: PubMed puts `<i>`/`<sup>` inside titles and
    abstracts, and `.text` stops at the first child -- which silently truncates
    an abstract at its first italicised species name, storing a fragment as if it
    were the whole retrieved text.
    """
    if element is None:
        return ""
    return re.sub(r"\s+", " ", "".join(element.itertext())).strip()


def _year(pub_date: ET.Element | None) -> str:
    """PubDate is either a <Year>, or a <MedlineDate> holding a range ('1998 Nov-Dec')."""
    if pub_date is None:
        return ""
    year = (pub_date.findtext("Year") or "").strip()
    if year:
        return year
    medline = _text(pub_date.find("MedlineDate"))
    found = re.search(r"\d{4}", medline)
    return found.group(0) if found else medline


class PubMedClient:
    def __init__(self, http: PoliteClient):
        self.http = http

    def search(self, term: str, limit: int) -> tuple[list[str] | None, dict[str, Any], RequestTrace]:
        """Run esearch. Returns (pmids, meta, trace); pmids is None if unreadable."""
        resp = self.http.get(
            f"{EUTILS}/esearch.fcgi",
            {"db": "pubmed", "term": term, "retmode": "json", "retmax": limit,
             "sort": "relevance", **ncbi_identity()},
        )
        try:
            payload = json.loads(resp.body) if resp.body.strip() else None
        except json.JSONDecodeError:
            payload = None  # NCBI serves an HTML error page with HTTP 200 under load
        if not isinstance(payload, dict):
            return None, {}, _trace(resp)

        result = payload.get("esearchresult")
        if not isinstance(result, dict):
            return None, {}, _trace(resp)
        # E-utilities reports bad terms inside a 200 body, not as a status code.
        error = result.get("ERROR") or result.get("error")
        if error:
            return None, {"error": str(error)}, _trace(resp)

        idlist = result.get("idlist")
        if not isinstance(idlist, list):
            return None, {}, _trace(resp)
        # Filter to digits before these reach a URL path/param: the same reasoning
        # as kegg._ORGANISM. An upstream that returned something else must not get
        # to shape the next request.
        pmids = [str(p) for p in idlist if _PMID.match(str(p).strip())]

        total = result.get("count")
        meta = {
            "total": int(total) if isinstance(total, str) and total.isdigit() else None,
            "query_translation": str(result.get("querytranslation", "")),
        }
        return pmids, meta, _trace(resp)

    def fetch(self, pmids: list[str]) -> tuple[list[Record] | None, list[str], RequestTrace]:
        """Run efetch and parse the articles. Returns (records, missing, trace).

        `records` is None if the response could not be parsed at all. `missing`
        holds PMIDs that were requested but absent from the response -- PubMed
        returns book chapters as <PubmedBookArticle>, which has no MedlineCitation
        and is deliberately not parsed here rather than guessed at.
        """
        resp = self.http.get(
            f"{EUTILS}/efetch.fcgi",
            {"db": "pubmed", "id": ",".join(pmids), "retmode": "xml", **ncbi_identity()},
        )
        try:
            # ElementTree does not fetch the DTD the response declares, and does not
            # expand external entities.
            root = ET.fromstring(resp.body)
        except ET.ParseError:
            return None, [], _trace(resp)

        found: dict[str, Record] = {}
        for article in root.iter("PubmedArticle"):
            record = self._record(article, resp)
            if record is not None:
                found[record.record_id] = record

        # Ordered by the PMID list, not by document order: esearch ranked these by
        # relevance and efetch does not preserve that ranking.
        records = [found[pmid] for pmid in pmids if pmid in found]
        missing = [pmid for pmid in pmids if pmid not in found]
        return records, missing, _trace(resp)

    def _record(self, article: ET.Element, resp) -> Record | None:
        citation = article.find("MedlineCitation")
        if citation is None:
            return None
        pmid = (citation.findtext("PMID") or "").strip()
        entry = citation.find("Article")
        if not _PMID.match(pmid) or entry is None:
            return None

        title = _text(entry.find("ArticleTitle"))
        sections: list[dict[str, str]] = []
        parts: list[str] = []
        for node in entry.findall("Abstract/AbstractText"):
            body = _text(node)
            if not body:
                continue
            label = (node.get("Label") or "").strip()
            sections.append({"label": label, "text": body})
            parts.append(f"{label}: {body}" if label else body)
        abstract = "\n\n".join(parts)

        # Scoped to PubmedData/ArticleIdList, NOT './/ArticleIdList/ArticleId'. Every
        # entry in <ReferenceList> carries its own ArticleIdList -- 35 of them in one
        # fixture article, of which exactly one is the article's. A descendant search
        # is saved only by document order (ArticleIdList precedes ReferenceList inside
        # PubmedData), and that accident runs out for an article with no DOI of its
        # own, where it would attribute a cited paper's DOI to this one.
        # PMCID matters as much as the DOI for anything downstream: a DOI resolves
        # to the publisher, which is usually paywalled and whose terms forbid bulk
        # retrieval. A PMCID means the article is in PubMed Central, which is the
        # licit route to full text. Same ArticleIdList, same scoping reasoning.
        ids = {"doi": "", "pmc": ""}
        for identifier in article.findall("PubmedData/ArticleIdList/ArticleId"):
            id_type = identifier.get("IdType", "")
            if id_type in ids and not ids[id_type]:
                ids[id_type] = (identifier.text or "").strip()
        doi, pmcid = ids["doi"], ids["pmc"]

        return Record(
            record_id=pmid,
            type="article",
            name=title or f"PMID {pmid}",
            url=f"{ARTICLE}{pmid}/",
            source="pubmed",
            retrieved_at=resp.fetched_at,
            cached=resp.cached,
            detail={
                "title": title,
                "abstract": abstract,
                "has_abstract": bool(abstract),
                "abstract_sections": sections,
                # The exact retrieved text a quote must be checked against. Kept as
                # one field so the check is a containment test on a single string
                # rather than a walk over title/section fragments.
                "quotable_text": "\n\n".join(p for p in (title, abstract) if p),
                "journal": _text(entry.find("Journal/Title")),
                "year": _year(entry.find("Journal/JournalIssue/PubDate")),
                "doi": doi,
                "pmcid": pmcid,
                # Presence in PMC is necessary but NOT sufficient for full-text
                # reuse: the PMC Open Access subset is a subset of PMC, and
                # per-article licences vary within it. Checking that is the
                # downstream consumer's job; recording the handle is this one's.
                "in_pmc": bool(pmcid),
            },
        )

    def abstracts(
        self, gene: str, organism: str = MTB_H37RV_NAME, limit: int = DEFAULT_LIMIT
    ) -> ToolResult:
        query: dict[str, Any] = {"gene": gene, "organism": organism, "limit": limit}
        gene, organism = gene.strip(), organism.strip()

        problems = []
        if not gene:
            problems.append("no gene identifier was supplied")
        if _QUERY_SYNTAX.search(gene):
            problems.append('gene contains PubMed query syntax (" [ ]), which would change '
                            "the meaning of the search rather than the terms searched for")
        if _QUERY_SYNTAX.search(organism):
            problems.append('organism contains PubMed query syntax (" [ ])')
        if not 1 <= limit <= MAX_LIMIT:
            problems.append(f"limit={limit} is outside 1-{MAX_LIMIT}")
        if problems:
            return ToolResult.build(
                query, [], resolved={"matched_by": "none"},
                notes=[f"Invalid argument(s), so no search was performed: {'; '.join(problems)}. "
                       f"An empty result here does NOT mean there is no literature on this gene."],
            )

        # Phrase-quoted so a multi-word organism stays one concept and a gene symbol
        # is not split by Automatic Term Mapping. `resolved.query_translation` below
        # reports what PubMed actually ran, which is the auditable part.
        term = f'"{gene}"'
        if organism:
            term = f'{term} AND "{organism}"'
        resolved: dict[str, Any] = {"term": term, "matched_by": "pubmed_search"}
        traces: list[RequestTrace] = []

        try:
            pmids, meta, trace = self.search(term, limit)
        except FetchError as exc:
            return ToolResult.build(
                query, [], resolved=resolved,
                notes=[f"PubMed esearch failed: HTTP {exc.status}. No articles were retrieved."],
            )
        traces.append(trace)
        resolved["query_translation"] = meta.get("query_translation", "")

        if pmids is None:
            detail = f" NCBI reported: {meta['error']}." if meta.get("error") else ""
            return ToolResult.build(
                query, [], resolved=resolved, requests=traces,
                notes=[f"PubMed esearch returned an unreadable or error response for {term!r}.{detail} "
                       f"No articles were retrieved. This is a retrieval failure, not evidence that "
                       f"no literature exists."],
            )

        total = meta.get("total")
        if not pmids:
            return ToolResult.build(
                query, [], resolved=resolved, requests=traces,
                notes=[f"PubMed returned no articles for {term!r}. Try a different gene name or drop "
                       f"the organism term; this is a search miss, not evidence that the gene is "
                       f"unstudied.", SEARCH_CAVEAT],
            )

        try:
            records, missing, fetch_trace = self.fetch(pmids)
        except FetchError as exc:
            return ToolResult.build(
                query, [], resolved=resolved, requests=traces,
                notes=[f"PubMed esearch found {len(pmids)} article(s) but efetch failed: HTTP "
                       f"{exc.status}, so none of their content was retrieved. The PMIDs are listed "
                       f"here for reference and are NOT citable from this result: "
                       f"{', '.join(pmids)}."],
            )
        traces.append(fetch_trace)

        if records is None:
            return ToolResult.build(
                query, [], resolved=resolved, requests=traces,
                notes=[f"PubMed efetch returned a response that could not be parsed as PubMed XML, so "
                       f"no article content was retrieved. The PMIDs are listed here for reference and "
                       f"are NOT citable from this result: {', '.join(pmids)}."],
            )

        notes = [SEARCH_CAVEAT, TEXTMINING_CAVEAT]
        if total is not None and total > len(records):
            notes.insert(0, f"{total} articles matched this query; the {len(records)} most relevant "
                            f"were retrieved. A finding absent from these abstracts may still be "
                            f"reported in the {total - len(records)} not retrieved.")
        # A record can match on metadata the model never sees -- MeSH terms,
        # keywords, substance lists -- and arrive with an abstract that never
        # mentions the gene. Observed: a query for "ahpC katG" returned a general
        # thioredoxin review whose abstract contains neither symbol. Its PMID is
        # still citable, so without this note a model could cite it for a claim
        # about the gene with nothing to quote. Name them.
        terms = [t for t in re.split(r"\s+", gene) if len(t) > 2]
        # Which query terms are actually present in the retrieved text. Recorded
        # per record so a downstream corpus can be filtered to papers that really
        # discuss a gene, rather than ones that merely matched its metadata.
        for record in records:
            haystack = record.detail["quotable_text"].lower()
            record.detail["mentions"] = [t for t in terms if t.lower() in haystack]
        unquotable = [r.record_id for r in records if terms and not r.detail["mentions"]]
        if unquotable:
            notes.append(
                f"PMID(s) {', '.join(unquotable)} matched the search but do NOT mention "
                f"{' or '.join(terms)} anywhere in the retrieved title or abstract -- they matched "
                f"on record metadata not shown here. There is nothing in them to quote for a claim "
                f"about this gene. {METADATA_MATCH_CAVEAT}"
            )

        # Named, not counted: a validator checking a quote against these records
        # needs to know which ones have no abstract to quote from.
        without = [r.record_id for r in records if not r.detail["has_abstract"]]
        if without:
            notes.append(f"No abstract is available for PMID(s) {', '.join(without)} -- PubMed holds "
                         f"only the title and metadata. Their `quotable_text` is the title alone.")
        if missing:
            notes.append(f"PMID(s) {', '.join(missing)} were returned by the search but no article "
                         f"record came back for them (book chapters and withdrawn records do this), "
                         f"so they are NOT citable from this result.")

        return ToolResult.build(query, records, resolved=resolved, requests=traces, notes=notes)
