"""KEGG REST lookups: gene -> pathways.

KEGG's REST interface returns bare TSV, so parsing is the bulk of this module.

Identifier resolution is deliberately **exact or nothing**. KEGG's `/find`
endpoint returns loose matches and would make the tool non-deterministic, so
symbols are resolved against the organism's full gene list, fetched once and
cached, and the match type is always reported.

Getting that list parsed correctly is subtler than it looks -- see `gene_index`.

Licence note: KEGG is free for academic use; commercial use requires a licence
from Pathway Solutions. This tool does not redistribute KEGG content -- it fetches
per query and caches locally for the caller -- but a commercial deployment is the
caller's responsibility, and the README says so.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Any

from kegg_string_mcp.http import FetchError, PoliteClient
from kegg_string_mcp.provenance import Record, RequestTrace, ToolResult

REST = "https://rest.kegg.jp"
ENTRY = "https://www.kegg.jp/entry/"

# Organism codes are 3-4 letters. Accept any case (KEGG's own IDs are lowercase,
# but callers paste mixed case) and normalise on use.
_KEGG_ID = re.compile(r"^[A-Za-z]{3,4}:\S+$")

# Prefixes that look like an organism code but namespace something other than a
# gene. Without this, `path:mtu00360` parses as organism "path".
_NON_GENE_PREFIXES = {"path", "map", "ko", "ec", "rn", "rc", "cpd", "gl", "dr", "ds", "br", "kot"}


def _trace(resp) -> RequestTrace:
    return RequestTrace(url=resp.audit_url, retrieved_at=resp.fetched_at, cached=resp.cached,
                        status=resp.status, content_sha256=resp.content_sha256)


def _rows(body: str) -> list[list[str]]:
    return [line.split("\t") for line in body.strip().splitlines() if line.strip()]


def _strip_prefix(value: str) -> str:
    """KEGG returns pathway IDs as 'path:mtu00360' in /link but 'mtu00360' in /list."""
    return value.split(":", 1)[1] if ":" in value else value


@dataclass
class GeneIndex:
    """Lookup table for one organism, plus the symbols that are not unique."""

    entries: dict[str, str] = field(default_factory=dict)
    ambiguous: set[str] = field(default_factory=set)
    locus_tags: set[str] = field(default_factory=set)

    def __contains__(self, key: str) -> bool:
        return key in self.entries

    def __getitem__(self, key: str) -> str:
        return self.entries[key]

    def get(self, key: str, default: str | None = None) -> str | None:
        return self.entries.get(key, default)


class KeggClient:
    def __init__(self, http: PoliteClient):
        self.http = http

    def gene_index(self, organism: str) -> tuple[GeneIndex, RequestTrace]:
        """Build symbol/locus-tag -> KEGG gene ID for one organism.

        KEGG writes the description column as ``"symbolA, symbolB; product name"``.
        When a gene has no assigned symbol the column is *only* the product name,
        with no semicolon -- and in M. tuberculosis that is the majority of genes
        (221 of 403 rows in the test fixture).

        Splitting on ';' unconditionally therefore indexed whole product names as
        if they were gene symbols, so `toxin`, `hydrolase` and `pseudogene` all
        resolved to arbitrary genes and were reported as exact symbol matches. The
        comma split compounded it: `beta-1,3-glucanase` became the keys `BETA-1`
        and `3-GLUCANASE`, so the real name failed while two fragments succeeded.

        Hence: only treat the leading field as symbols when a ';' is actually
        present, and index locus tags in a second pass so a real identifier can
        never be shadowed by a symbol from an earlier row.
        """
        resp = self.http.get(f"{REST}/list/{organism}")
        rows = _rows(resp.body)

        symbols: dict[str, str] = {}
        ambiguous: set[str] = set()
        for row in rows:
            if len(row) < 2:
                continue
            kegg_id, description = row[0], row[-1]
            if ";" not in description:
                continue  # no symbol field at all -- product name only
            for symbol in description.split(";")[0].split(","):
                symbol = symbol.strip().upper()
                if not symbol:
                    continue
                if symbol in symbols and symbols[symbol] != kegg_id:
                    ambiguous.add(symbol)
                symbols.setdefault(symbol, kegg_id)

        locus_tags = {}
        for row in rows:
            if len(row) >= 2:
                locus_tags.setdefault(_strip_prefix(row[0]).upper(), row[0])

        # Locus tags win: a real identifier must never lose to a symbol collision.
        return (
            GeneIndex(entries={**symbols, **locus_tags}, ambiguous=ambiguous,
                      locus_tags=set(locus_tags)),
            _trace(resp),
        )

    def pathway_names(self, organism: str) -> tuple[dict[str, str], RequestTrace]:
        resp = self.http.get(f"{REST}/list/pathway/{organism}")
        names = {_strip_prefix(row[0]): row[1] for row in _rows(resp.body) if len(row) >= 2}
        return names, _trace(resp)

    def pathways(self, gene: str, organism: str = "mtu") -> ToolResult:
        query: dict[str, Any] = {"gene": gene, "organism": organism}
        traces: list[RequestTrace] = []
        notes: list[str] = []
        gene = gene.strip()

        prefix = gene.split(":", 1)[0].lower() if ":" in gene else ""
        if prefix in _NON_GENE_PREFIXES:
            return ToolResult.build(
                query, [], resolved={"matched_by": "none"},
                notes=[f"'{gene}' is a KEGG '{prefix}:' identifier, which namespaces "
                       f"{'pathways' if prefix in {'path', 'map'} else 'a non-gene entity'}, "
                       f"not a gene. This tool takes a gene ID, locus tag, or symbol."],
            )

        if _KEGG_ID.match(gene):
            id_organism, locus = gene.split(":", 1)
            id_organism = id_organism.lower()          # KEGG organism codes are lowercase
            kegg_id, matched_by = f"{id_organism}:{locus}", "kegg_id"
            # A fully-qualified ID carries its own organism. Trust it over the
            # `organism` argument, otherwise pathway names get looked up in the
            # wrong organism and every record comes back "(name unavailable)"
            # with correct IDs and no explanation.
            if id_organism != organism:
                notes.append(
                    f"'{gene}' is a KEGG ID for organism '{id_organism}', but organism="
                    f"'{organism}' was requested. Used '{id_organism}', from the identifier."
                )
                organism = id_organism
                query["organism_used"] = id_organism
        else:
            try:
                index, trace = self.gene_index(organism)
            except FetchError as exc:
                return ToolResult.build(
                    query, [], resolved={"matched_by": "none"},
                    notes=[f"Could not fetch the gene list for organism '{organism}': HTTP "
                           f"{exc.status}. '{organism}' may not be a valid KEGG organism code. "
                           f"No lookup was performed."],
                )
            traces.append(trace)
            key = gene.upper()
            kegg_id = index.get(key)
            matched_by = "locus_tag_or_symbol" if kegg_id else "none"
            # A locus tag that happens to collide with an ambiguous symbol resolved
            # unambiguously -- locus tags win the index, so do not warn about it.
            matched_symbol = key not in index.locus_tags
            if kegg_id and matched_symbol and key in index.ambiguous:
                notes.append(
                    f"'{gene}' is not a unique symbol in organism '{organism}'; it matches more "
                    f"than one gene. Used {kegg_id}. Pass a locus tag to disambiguate."
                )
                query["ambiguous_symbol"] = True

        if not kegg_id:
            return ToolResult.build(
                query, [], resolved={"matched_by": "none"}, requests=traces,
                notes=[f"'{gene}' did not match any locus tag or gene symbol in KEGG organism "
                       f"'{organism}'. No pathways were looked up. This is a resolution failure, "
                       f"not evidence that the gene has no pathways."],
            )

        try:
            link = self.http.get(f"{REST}/link/pathway/{kegg_id}")
        except FetchError as exc:
            return ToolResult.build(
                query, [], resolved={"kegg_gene_id": kegg_id, "matched_by": matched_by},
                requests=traces, notes=notes + [f"KEGG /link failed: HTTP {exc.status}. "
                                                f"No pathway data retrieved."],
            )
        traces.append(_trace(link))

        pathway_ids = [_strip_prefix(row[1]) for row in _rows(link.body) if len(row) >= 2]
        resolved = {"kegg_gene_id": kegg_id, "matched_by": matched_by}
        if not pathway_ids:
            if matched_by == "kegg_id":
                # This branch skips the gene index, so nothing has verified the ID.
                # KEGG answers /link for an unknown gene with HTTP 200 and an empty
                # body, which is indistinguishable from a real gene with no pathways
                # -- so `mtu:NOTAGENE`, and the very natural `mtu:katG`, both drew a
                # confident "the gene exists in KEGG" claim.
                notes.append(
                    f"KEGG returned no pathway assignments for '{kegg_id}'. This does NOT confirm "
                    f"the identifier exists: KEGG answers identically for an unknown gene. Note "
                    f"that '{organism}:SYMBOL' is not a valid KEGG gene ID -- pass a bare symbol, "
                    f"or the locus tag form '{organism}:LOCUS'."
                )
            else:
                notes.append(f"KEGG returned no pathway assignments for {kegg_id}. The gene exists "
                             f"in KEGG but is not mapped to any pathway in this organism.")
            return ToolResult.build(query, [], resolved=resolved, requests=traces, notes=notes)

        # Names are a nicety; the IDs are the citable part. A failure here degrades
        # rather than losing the successfully retrieved pathway IDs.
        try:
            names, name_trace = self.pathway_names(organism)
            traces.append(name_trace)
        except FetchError as exc:
            names = {}
            notes.append(f"Could not fetch pathway names for '{organism}': HTTP {exc.status}. "
                         f"The pathway IDs below are still valid and citable.")

        # No `and names` guard: an empty-but-successful /list/pathway response is
        # precisely the case that needs the note, and guarding on `names` produced
        # records reading "(name unavailable)" with an empty notes list.
        missing = [pid for pid in pathway_ids if pid not in names]
        if missing:
            notes.append(f"KEGG's pathway list for '{organism}' had no name for: "
                         f"{', '.join(missing)}. The IDs are still valid and citable.")

        records = [
            Record(
                record_id=pid, type="pathway", name=names.get(pid, "(name unavailable)"),
                url=f"{ENTRY}{pid}", source="kegg",
                retrieved_at=link.fetched_at, cached=link.cached,
                detail={"kegg_gene_id": kegg_id},
            )
            for pid in pathway_ids
        ]
        return ToolResult.build(query, records, resolved=resolved, requests=traces, notes=notes)
