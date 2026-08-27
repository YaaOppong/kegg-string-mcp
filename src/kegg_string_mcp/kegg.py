"""KEGG REST lookups: gene -> pathways.

KEGG's REST interface returns bare TSV, so parsing is the bulk of this module.
Two identifier subtleties matter:

* KEGG gene IDs are `{organism}:{locus}` (e.g. `mtu:Rv1908c`). For M. tuberculosis
  the locus part is the Rv number, so a locus tag resolves directly.
* Gene *symbols* (`katG`) do not. Rather than use KEGG's fuzzy `/find` endpoint --
  which can return several loose matches and would make the tool non-deterministic
  -- this resolves symbols against the organism's full gene list, fetched once and
  cached. Exact match or nothing, and the match type is reported.

Licence note: KEGG is free for academic use; commercial use requires a licence
from Pathway Solutions. This tool does not redistribute KEGG content -- it fetches
per query and caches locally for the caller -- but a commercial deployment is the
caller's responsibility, and the README says so.
"""

from __future__ import annotations

import re
from typing import Any

from kegg_string_mcp.http import FetchError, PoliteClient
from kegg_string_mcp.provenance import Record, RequestTrace, ToolResult

REST = "https://rest.kegg.jp"
ENTRY = "https://www.kegg.jp/entry/"
_KEGG_ID = re.compile(r"^[a-z]{3,4}:\S+$")


def _trace(resp) -> RequestTrace:
    return RequestTrace(url=resp.url, retrieved_at=resp.fetched_at, cached=resp.cached,
                        status=resp.status, content_sha256=resp.content_sha256)


def _rows(body: str) -> list[list[str]]:
    return [line.split("\t") for line in body.strip().splitlines() if line.strip()]


def _strip_prefix(value: str) -> str:
    """KEGG returns pathway IDs as 'path:mtu00360' in /link but 'mtu00360' in /list."""
    return value.split(":", 1)[1] if ":" in value else value


class KeggClient:
    def __init__(self, http: PoliteClient):
        self.http = http

    def gene_index(self, organism: str) -> tuple[dict[str, str], RequestTrace]:
        """symbol (upper) -> KEGG gene ID, from the organism's full gene list."""
        resp = self.http.get(f"{REST}/list/{organism}")
        index: dict[str, str] = {}
        for row in _rows(resp.body):
            if len(row) < 2:
                continue
            # /list/{org} is 4 columns: id, type, location, description -- but the
            # column count has varied across KEGG versions, so take the description
            # as the last field rather than a fixed index.
            kegg_id, description = row[0], row[-1]
            index.setdefault(_strip_prefix(kegg_id).upper(), kegg_id)  # locus tag
            # Description is "symbolA, symbolB; long product name".
            for symbol in description.split(";")[0].split(","):
                symbol = symbol.strip()
                if symbol:
                    index.setdefault(symbol.upper(), kegg_id)
        return index, _trace(resp)

    def pathway_names(self, organism: str) -> tuple[dict[str, str], RequestTrace]:
        resp = self.http.get(f"{REST}/list/pathway/{organism}")
        names = {_strip_prefix(row[0]): row[1] for row in _rows(resp.body) if len(row) >= 2}
        return names, _trace(resp)

    def pathways(self, gene: str, organism: str = "mtu") -> ToolResult:
        query: dict[str, Any] = {"gene": gene, "organism": organism}
        traces: list[RequestTrace] = []
        notes: list[str] = []

        if _KEGG_ID.match(gene):
            kegg_id, matched_by = gene, "kegg_id"
        else:
            index, trace = self.gene_index(organism)
            traces.append(trace)
            kegg_id = index.get(gene.strip().upper())
            matched_by = "locus_tag_or_symbol" if kegg_id else "none"

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
                requests=traces, notes=[f"KEGG /link failed: HTTP {exc.status}. No pathway data retrieved."],
            )
        traces.append(_trace(link))

        pathway_ids = [_strip_prefix(row[1]) for row in _rows(link.body) if len(row) >= 2]
        if not pathway_ids:
            notes.append(f"KEGG returned no pathway assignments for {kegg_id}. The gene exists in "
                         f"KEGG but is not mapped to any pathway in this organism.")
            return ToolResult.build(query, [], resolved={"kegg_gene_id": kegg_id, "matched_by": matched_by},
                                    requests=traces, notes=notes)

        names, name_trace = self.pathway_names(organism)
        traces.append(name_trace)

        records = [
            Record(
                record_id=pid, type="pathway", name=names.get(pid, "(name unavailable)"),
                url=f"{ENTRY}{pid}", source="kegg",
                retrieved_at=link.fetched_at, cached=link.cached,
                detail={"kegg_gene_id": kegg_id},
            )
            for pid in pathway_ids
        ]
        return ToolResult.build(query, records,
                                resolved={"kegg_gene_id": kegg_id, "matched_by": matched_by},
                                requests=traces, notes=notes)
