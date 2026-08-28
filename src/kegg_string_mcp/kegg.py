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

# KEGG organism codes are 3-4 lowercase letters. The value goes into a URL *path*,
# and a model can pass anything: unvalidated, "../../etc" normalises to
# rest.kegg.jp/etc, so the tool issues a request the caller never intended and
# reports the resulting 404 as if it described the gene.
_ORGANISM = re.compile(r"^[a-z]{3,4}$")


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

    @staticmethod
    def _require_organism(organism: str) -> str:
        """`organism` goes into a URL PATH. pathways() guards this, but the agent
        pipeline calls gene_index and pathway_sizes directly and bypassed it."""
        organism = (organism or "").strip().lower()
        if not _ORGANISM.match(organism):
            raise ValueError(f"{organism!r} is not a valid KEGG organism code "
                             f"(3-4 lowercase letters, e.g. 'mtu')")
        return organism

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
        organism = self._require_organism(organism)
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
        # _KEGG_ID deliberately accepts mixed case and line ~165 lowercases the
        # organism from a qualified ID for that reason. Validating before
        # normalising made pathways("mtu:Rv1908c", "MTU") fail outright.
        organism = organism.strip().lower()

        if not _ORGANISM.match(organism):
            return ToolResult.build(
                query, [], resolved={"matched_by": "none"},
                notes=[(f"'{organism}' is not a valid KEGG organism code (3-4 lowercase letters, "
                       f"e.g. 'mtu' for M. tuberculosis H37Rv). No lookup was performed.")],
            )

        if not gene:
            return ToolResult.build(
                query, [], resolved={"matched_by": "none"},
                notes=["No gene identifier was supplied. No lookup was performed."],
            )

        prefix = gene.split(":", 1)[0].lower() if ":" in gene else ""
        if prefix in _NON_GENE_PREFIXES:
            return ToolResult.build(
                query, [], resolved={"matched_by": "none"},
                notes=[(f"'{gene}' is a KEGG '{prefix}:' identifier, which namespaces "
                       f"{'pathways' if prefix in {'path', 'map'} else 'a non-gene entity'}, "
                       f"not a gene. This tool takes a gene ID, locus tag, or symbol.")],
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
                    notes=[(f"Could not fetch the gene list for organism '{organism}': HTTP "
                           f"{exc.status}. '{organism}' may not be a valid KEGG organism code. "
                           f"No lookup was performed.")],
                )
            traces.append(trace)
            key = gene.upper()
            kegg_id = index.get(key)
            if not kegg_id:
                matched_by = "none"
            else:
                # Report WHICH path matched. Collapsing both into one value left the
                # caller unable to tell a locus-tag hit from a symbol hit, which
                # matters when the two interpretations disagree.
                matched_by = "locus_tag" if key in index.locus_tags else "symbol"
            if kegg_id and key in index.ambiguous:
                if matched_by == "locus_tag":
                    # Resolved unambiguously (locus tags win the index), but the
                    # caller still deserves to know a second reading existed.
                    notes.append(
                        f"'{gene}' resolved as a locus tag ({kegg_id}). It is also an ambiguous "
                        f"gene symbol in organism '{organism}'; the locus tag was used."
                    )
                else:
                    # Only a genuine symbol hit is ambiguous. Setting this for a
                    # locus-tag hit made agents branching on the flag hedge on an
                    # exact match.
                    query["ambiguous_symbol"] = True
                    notes.append(
                        f"'{gene}' is not a unique symbol in organism '{organism}'; it matches more "
                        f"than one gene. Used {kegg_id}. Pass a locus tag to disambiguate."
                    )

        if not kegg_id:
            return ToolResult.build(
                query, [], resolved={"matched_by": "none"}, requests=traces,
                notes=[(f"'{gene}' did not match any locus tag or gene symbol in KEGG organism "
                       f"'{organism}'. No pathways were looked up. This is a resolution failure, "
                       f"not evidence that the gene has no pathways.")],
            )

        try:
            link = self.http.get(f"{REST}/link/pathway/{kegg_id}")
        except FetchError as exc:
            return ToolResult.build(
                query, [], resolved={"kegg_gene_id": kegg_id, "matched_by": matched_by},
                requests=traces, notes=notes + [(f"KEGG /link failed: HTTP {exc.status}. "
                                                f"No pathway data retrieved.")],
            )
        traces.append(_trace(link))

        pathway_ids = [_strip_prefix(row[1]) for row in _rows(link.body) if len(row) >= 2]
        resolved = {"kegg_gene_id": kegg_id, "matched_by": matched_by}
        if not pathway_ids:
            existence = True  # the index branch already proved the gene exists
            if matched_by == "kegg_id":
                # The qualified-ID branch skips the index, and KEGG answers /link for
                # an unknown gene with HTTP 200 and an empty body -- indistinguishable
                # from a real gene with no pathways. Rather than hedge in prose, spend
                # one (cached) request and report which case it actually is.
                index = None
                try:
                    index, index_trace = self.gene_index(organism)
                    traces.append(index_trace)
                    if not index.entries:
                        # Empty-but-successful /list/{organism} -- proxy truncation,
                        # a maintenance 200 -- is NOT proof the gene is absent, and
                        # that 200 is cached for the full TTL so the false claim
                        # would repeat. Same fetched-vs-empty distinction as
                        # names_fetched below; getting it wrong here fabricates.
                        existence = None
                    else:
                        # Locus tags ONLY: the index also holds gene symbols, and the
                        # locus part of a qualified KEGG ID must be a locus tag.
                        # Checking the whole index reported 'mtu:katG' as existing.
                        existence = _strip_prefix(kegg_id).upper() in index.locus_tags
                except FetchError:
                    existence = None

            if existence is True:
                notes.append(f"KEGG returned no pathway assignments for {kegg_id}. The gene exists "
                             f"in KEGG but is not mapped to any pathway in this organism.")
            elif existence is False:
                # Drop kegg_gene_id: leaving a concrete, citable-looking ID beside
                # matched_by="none" invites a model to cite the thing we just said
                # did not resolve. The sibling failure return omits it too.
                resolved = {"matched_by": "none"}
                # The old `if ":" in kegg_id` guard was always true on this branch
                # (it is reachable only via _KEGG_ID, which requires a colon), so the
                # "that's a symbol, not a gene ID" advice was given even for
                # mtu:NOTAGENE, where it is simply wrong. Test what it claims.
                locus_key = _strip_prefix(kegg_id).upper()
                looks_like_symbol = (index is not None and locus_key in index.entries
                                     and locus_key not in index.locus_tags)
                hint = (f" Note that '{organism}:SYMBOL' is not a valid KEGG gene ID -- pass a bare "
                        f"symbol, or the locus tag form.") if looks_like_symbol else ""
                notes.append(f"'{kegg_id}' was not found in KEGG organism '{organism}', so no "
                             f"pathways were looked up. This is a resolution failure, not evidence "
                             f"that the gene has no pathways.{hint}")
            else:
                notes.append(f"KEGG returned no pathway assignments for {kegg_id}, and the gene "
                             f"list could not be fetched to confirm the identifier exists.")
            return ToolResult.build(query, [], resolved=resolved, requests=traces, notes=notes)

        # Names are a nicety; the IDs are the citable part. A failure here degrades
        # rather than losing the successfully retrieved pathway IDs.
        names_fetched = True
        try:
            names, name_trace = self.pathway_names(organism)
            traces.append(name_trace)
        except FetchError as exc:
            names, names_fetched = {}, False
            notes.append(f"Could not fetch pathway names for '{organism}': HTTP {exc.status}. "
                         f"The pathway IDs below are still valid and citable.")

        # Gate on whether the list was actually FETCHED, not on whether it is
        # non-empty. Guarding on truthiness suppressed the note for an empty-but-
        # successful response; removing the guard entirely made the tool assert
        # "KEGG's pathway list had no name for X" about a response it never
        # received, contradicting the fetch-failure note directly above.
        missing = [pid for pid in pathway_ids if pid not in names]
        if missing and names_fetched:
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

    def pathway_sizes(self, organism: str) -> tuple[dict[str, int], RequestTrace]:
        """Gene count per pathway, from one cached request.

        Essential for interpreting *shared* pathways. In M. tuberculosis,
        `mtu01100` ("Metabolic pathways") holds 698 of ~4000 genes, so two genes
        sharing it says almost nothing; `mtu00983` holds 11, and sharing that is a
        real signal. Reporting "these genes share a pathway" without the size is
        the pathway-shaped version of the hub-gene trap.
        """
        organism = self._require_organism(organism)
        resp = self.http.get(f"{REST}/link/pathway/{organism}")
        sizes: dict[str, int] = {}
        for row in _rows(resp.body):
            if len(row) >= 2:
                sizes[_strip_prefix(row[1])] = sizes.get(_strip_prefix(row[1]), 0) + 1
        return sizes, _trace(resp)
