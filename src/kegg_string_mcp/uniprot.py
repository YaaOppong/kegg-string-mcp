"""UniProt lookups via the REST API: gene -> curated protein function.

This is the tool that covers the hole the others leave. KEGG assigns a pathway to
only 1,171 of 4,008 M. tuberculosis genes -- 29% -- and gyrA, one of the most
studied genes in TB, is in the 71% with none. UniProt annotates the whole
proteome, so "KEGG has no pathway for this gene" stops being the end of the
answer.

Like a PubMed abstract and unlike a KEGG pathway ID, a UniProt function statement
is **prose**: a claim drawn from it is the model deciding what a record says. So
each record carries `quotable_text`, and the agent layer's span validator applies
unchanged.

What UniProt adds that neither other source does is **evidence codes**. Every
statement is tagged with how it is known, and experimental statements carry the
PMIDs that support them:

    ECO:0000269  experimental evidence, manual assertion   -> PubMed IDs attached
    ECO:0000255  match to a sequence model (a HAMAP rule)  -> inferred, no experiment

That distinction is the same shape as STRING's textmining channel: a statement
inferred from sequence similarity is not independent evidence about *this*
protein, and presenting it as though it were is the failure this codebase keeps
guarding against. Statements are therefore returned grouped by evidence tier, and
the supporting PMIDs are surfaced so a claim can be traced to the paper behind it
-- and those PMIDs are already citable record IDs elsewhere in this pipeline.
"""

from __future__ import annotations

import json
import re
from typing import Any

from kegg_string_mcp.http import FetchError, PoliteClient
from kegg_string_mcp.provenance import Record, RequestTrace, ToolResult

REST = "https://rest.uniprot.org/uniprotkb/search"
ENTRY = "https://www.uniprot.org/uniprotkb/"

MTB_H37RV = 83332
DEFAULT_LIMIT = 3
MAX_LIMIT = 10

FIELDS = ("accession,id,protein_name,gene_names,organism_name,cc_function,"
          "cc_catalytic_activity,cc_subunit,xref_pdb,ec")

# UniProt's query grammar. A value containing these would change what is asked
# rather than what is asked about -- the same reasoning as pubmed._QUERY_SYNTAX.
_QUERY_SYNTAX = re.compile(r'[:()\[\]"\s]')
_ACCESSION = re.compile(r"^[A-Z0-9]{6,10}$")

# ECO code -> how the statement is known. The tier is what matters downstream;
# the code is kept verbatim in the record so nothing is lost in translation.
EVIDENCE_TIERS = {
    "ECO:0000269": "experimental",     # experimental evidence, manual assertion
    "ECO:0000305": "curator_inference",
    "ECO:0000250": "sequence_similarity",
    "ECO:0000255": "sequence_model",   # e.g. a HAMAP rule
    "ECO:0000256": "automatic",
    "ECO:0000312": "imported",
    "ECO:0000313": "imported",
    "ECO:0007829": "automatic",
}
INFERRED_TIERS = {"sequence_similarity", "sequence_model", "automatic", "imported"}

SIMILARITY_CAVEAT = (
    "Statements tiered as sequence_similarity, sequence_model, automatic or imported are "
    "INFERRED -- from a rule or from homology -- not measured on this protein. They are not "
    "independent evidence about this gene, and should not be reported as though an "
    "experiment established them."
)


def _trace(resp) -> RequestTrace:
    return RequestTrace(url=resp.audit_url, retrieved_at=resp.fetched_at, cached=resp.cached,
                        status=resp.status, content_sha256=resp.content_sha256)


def _statements(entry: dict[str, Any], kind: str) -> list[dict[str, Any]]:
    """Comment texts of one type, each with its evidence tier and supporting PMIDs."""
    out: list[dict[str, Any]] = []
    for comment in entry.get("comments", []):
        if comment.get("commentType") != kind:
            continue
        for text in comment.get("texts", []):
            value = (text.get("value") or "").strip()
            if not value:
                continue
            codes, pmids = [], []
            for evidence in text.get("evidences", []):
                code = evidence.get("evidenceCode", "")
                if code:
                    codes.append(code)
                if evidence.get("source") == "PubMed" and evidence.get("id"):
                    pmids.append(str(evidence["id"]))
            tiers = sorted({EVIDENCE_TIERS.get(c, "unknown") for c in codes}) or ["unstated"]
            out.append({
                "text": value,
                "evidence_codes": sorted(set(codes)),
                "tiers": tiers,
                # An experimental tier anywhere means at least one statement-level
                # experiment; mixed statements are common and are reported as mixed.
                "experimental": "experimental" in tiers,
                "supporting_pmids": sorted(set(pmids)),
            })
    return out


class UniProtClient:
    def __init__(self, http: PoliteClient):
        self.http = http

    def search(self, gene: str, organism_id: int,
               limit: int) -> tuple[list[dict] | None, str, RequestTrace]:
        """Returns (entries, release, trace). `release` is UniProt's own statement of
        which release answered -- the same class of fact as a KEGG release date, and
        the thing that makes an annotation reproducible."""
        resp = self.http.get(REST, {"query": f"organism_id:{organism_id} AND gene:{gene}",
                                    "fields": FIELDS, "size": limit, "format": "json"})
        release = resp.headers.get("x-uniprot-release", "")
        try:
            payload = json.loads(resp.body) if resp.body.strip() else None
        except json.JSONDecodeError:
            payload = None
        if not isinstance(payload, dict) or not isinstance(payload.get("results"), list):
            return None, release, _trace(resp)
        return [r for r in payload["results"] if isinstance(r, dict)], release, _trace(resp)

    def _record(self, entry: dict[str, Any], resp_meta: tuple[str, bool]) -> Record | None:
        accession = str(entry.get("primaryAccession") or "").strip()
        if not _ACCESSION.match(accession):
            return None

        description = entry.get("proteinDescription", {})
        name = ((description.get("recommendedName") or {}).get("fullName") or {}).get("value", "")
        if not name:
            submitted = description.get("submissionNames") or []
            name = (submitted[0].get("fullName", {}).get("value", "") if submitted else "")

        functions = _statements(entry, "FUNCTION")
        catalytic = [c.get("reaction", {}).get("name", "")
                     for c in entry.get("comments", []) if c.get("commentType") == "CATALYTIC ACTIVITY"]
        subunit = _statements(entry, "SUBUNIT")

        genes = [g.get("geneName", {}).get("value", "") for g in entry.get("genes", [])]
        ordered = [ol.get("value", "") for g in entry.get("genes", [])
                   for ol in g.get("orderedLocusNames", [])]
        pdb = [x.get("id", "") for x in entry.get("uniProtKBCrossReferences", [])
               if x.get("database") == "PDB"]

        fetched_at, cached = resp_meta
        return Record(
            record_id=accession,
            type="protein",
            name=name or accession,
            url=f"{ENTRY}{accession}",
            source="uniprot",
            retrieved_at=fetched_at,
            cached=cached,
            detail={
                "entry_name": entry.get("uniProtkbId", ""),
                "gene_names": [g for g in genes if g],
                "locus_tags": [o for o in ordered if o],
                "reviewed": entry.get("entryType", "").startswith("UniProtKB reviewed"),
                "function_statements": functions,
                "catalytic_activity": [c for c in catalytic if c],
                "subunit": subunit,
                "pdb": pdb,
                "has_experimental_function": any(f["experimental"] for f in functions),
                # Everything a claim may be quoted against, in one string, matching
                # how pubmed.py stores it so the span validator needs no special case.
                "quotable_text": "\n\n".join(
                    [name] + [f["text"] for f in functions]
                    + [s["text"] for s in subunit] + [c for c in catalytic if c]
                ).strip(),
            },
        )

    def protein(self, gene: str, organism_id: int = MTB_H37RV,
                limit: int = DEFAULT_LIMIT) -> ToolResult:
        query: dict[str, Any] = {"gene": gene, "organism_id": organism_id, "limit": limit}
        gene = gene.strip()

        problems = []
        if not gene:
            problems.append("no gene identifier was supplied")
        elif _QUERY_SYNTAX.search(gene):
            problems.append('gene contains UniProt query syntax (: ( ) [ ] " or whitespace), '
                            "which would change the meaning of the search")
        if organism_id <= 0:
            problems.append(f"organism_id={organism_id} is not a valid NCBI taxon ID")
        if not 1 <= limit <= MAX_LIMIT:
            problems.append(f"limit={limit} is outside 1-{MAX_LIMIT}")
        if problems:
            return ToolResult.build(
                query, [], resolved={"matched_by": "none"},
                notes=[(f"Invalid argument(s), so no lookup was performed: {'; '.join(problems)}. "
                       f"An empty result here does NOT mean the protein is unannotated.")],
            )

        try:
            entries, release, trace = self.search(gene, organism_id, limit)
        except FetchError as exc:
            return ToolResult.build(query, [], resolved={"matched_by": "none"},
                                    notes=[(f"UniProt search failed: HTTP {exc.status}. "
                                           f"No records were retrieved.")])
        traces = [trace]

        resolved_base = {"matched_by": "none", "uniprot_release": release}
        if entries is None:
            return ToolResult.build(
                query, [], resolved=resolved_base, requests=traces,
                notes=[("UniProt returned an unreadable response (expected a JSON object with a "
                       "results list). No records were retrieved. This is a retrieval failure, "
                       "not evidence that the protein is unannotated.")])
        if not entries:
            return ToolResult.build(
                query, [], resolved=resolved_base, requests=traces,
                notes=[(f"UniProt returned no entry for gene '{gene}' in organism {organism_id}. "
                       f"Try a locus tag (e.g. Rv1908c) or a different symbol; this is a "
                       f"resolution failure, not evidence that the protein is unannotated.")])

        records = [r for r in (self._record(e, (trace.retrieved_at, trace.cached))
                               for e in entries) if r is not None]
        notes: list[str] = []

        if not records:
            # UniProt answered with entries but none of them parsed -- a malformed or
            # changed payload. Returning an empty record list with no note would read
            # exactly like "this protein is unannotated", which is the silence this
            # codebase exists to avoid.
            return ToolResult.build(
                query, [], resolved=resolved_base, requests=traces,
                notes=[(f"UniProt returned {len(entries)} entr(y/ies) for '{gene}' but none could "
                        f"be parsed into a record (no valid accession). This is a parsing "
                        f"failure, not evidence that the protein is unannotated.")])

        distinct = sorted({r.record_id for r in records})
        if len(distinct) > 1:
            # Several proteins matched the gene name -- paralogues, or a symbol shared
            # across entries. Reporting only the first accession as "the" answer would
            # hide that a choice was made.
            notes.append(f"{len(distinct)} UniProt entries matched gene '{gene}': "
                         f"{', '.join(distinct)}. They are all returned; "
                         f"resolved.accession names only the first.")
        if any(any(t in INFERRED_TIERS for f in r.detail["function_statements"] for t in f["tiers"])
               for r in records):
            notes.append(SIMILARITY_CAVEAT)
        unevidenced = [r.record_id for r in records
                       if r.detail["function_statements"] and not r.detail["has_experimental_function"]]
        if unevidenced:
            notes.append(f"No function statement for {', '.join(unevidenced)} carries experimental "
                         f"evidence (ECO:0000269); all are inferred.")
        no_function = [r.record_id for r in records if not r.detail["function_statements"]]
        if no_function:
            notes.append(f"UniProt holds no FUNCTION statement for {', '.join(no_function)} -- "
                         f"the entry exists but its function is not described. Nothing to quote.")

        return ToolResult.build(
            query, records,
            resolved={"accession": records[0].record_id if records else "",
                      "matched_by": "uniprot_search", "n_entries": len(records),
                      "uniprot_release": release},
            requests=traces, notes=notes,
        )
