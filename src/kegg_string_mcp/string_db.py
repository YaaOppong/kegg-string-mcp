"""STRING API lookups: gene -> interaction partners.

STRING returns a combined confidence score plus a per-channel breakdown. The
combined score **includes a textmining channel** -- co-mention in the literature.
That matters for any downstream summary: an interaction supported only by
textmining is not independent evidence from a literature citation about the same
pair, and treating it as such double-counts.

So this tool returns every channel score verbatim and adds a derived
`evidence_beyond_textmining` flag. It deliberately does not attempt to recompute
a textmining-free combined score: STRING combines channels probabilistically with
a prior correction, and a hand-rolled approximation would look authoritative while
being wrong.

STRING asks callers to identify themselves via `caller_identity`; set
STRING_CALLER_IDENTITY to something traceable to you.
"""

from __future__ import annotations

import json
import os
from typing import Any

from kegg_string_mcp.http import FetchError, PoliteClient
from kegg_string_mcp.provenance import Record, RequestTrace, ToolResult

API = "https://string-db.org/api"
NETWORK = "https://string-db.org/network/"
MTB_H37RV = 83332

# STRING's per-channel score fields. tscore is textmining and is held apart.
CHANNELS = {
    "nscore": "neighborhood",
    "fscore": "fusion",
    "pscore": "cooccurrence",
    "ascore": "coexpression",
    "escore": "experimental",
    "dscore": "database",
}
TEXTMINING = "tscore"

# STRING's own confidence bands: 0.15 low, 0.4 medium, 0.7 high, 0.9 highest.
# `evidence_beyond_textmining` thresholds at medium rather than testing for
# non-zero, because near-zero channel values are noise. Real example: katG-embB
# scores 0.964 combined, of which tscore is 0.963 and the only other non-zero
# channel is ascore 0.044. A "> 0" test calls that non-textmining evidence; it
# plainly is not.
MEDIUM_CONFIDENCE = 0.4

# STRING's confidence scores are 0-1000. Out-of-range values are not rejected by
# the API -- required_score=99999 simply returns zero partners, which this tool
# would otherwise report as a genuine empty result rather than a bad argument.
MIN_SCORE, MAX_SCORE = 0, 1000
MAX_LIMIT = 1000

# Identifiers per /network request. STRING accepts a list in one parameter, but a
# gene set of a few hundred would build a URL long enough for an intermediary to
# truncate -- and a silently shortened identifier list returns fewer edges, which
# reads exactly like "those pairs do not interact".
NETWORK_BATCH = 100


def _batched(items: list[str], size: int) -> list[list[str]]:
    return [items[i:i + size] for i in range(0, len(items), size)]


def _scores(row: dict[str, Any]) -> dict[str, Any]:
    """The score block both STRING endpoints return, parsed once.

    /network and /interaction_partners share a row schema, so the textmining
    bookkeeping this module exists for must not be re-implemented per endpoint.
    """
    channels = {label: float(row.get(field, 0) or 0) for field, label in CHANNELS.items()}
    best_other = max(channels.values(), default=0.0)
    return {"combined_score": float(row.get("score", 0) or 0),
            "channels": channels,
            "textmining_score": float(row.get(TEXTMINING, 0) or 0),
            "max_non_textmining_score": best_other,
            "evidence_beyond_textmining": best_other >= MEDIUM_CONFIDENCE}


def _edge_record(row: dict[str, Any], species: int, resp) -> Record:
    """One undirected edge, citable by the pair of STRING IDs it connects."""
    a, b = row.get("stringId_A", ""), row.get("stringId_B", "")
    name_a, name_b = row.get("preferredName_A", a), row.get("preferredName_B", b)
    first, second = sorted((a, b))
    return Record(
        record_id=f"{first}|{second}",
        type="interaction",
        name=f"{name_a}--{name_b}",
        url=f"https://string-db.org/cgi/network?identifiers={first}%0d{second}&species={species}",
        source="string",
        retrieved_at=resp.fetched_at,
        cached=resp.cached,
        detail=_scores(row) | {"string_id_a": a, "string_id_b": b,
                               "preferred_name_a": name_a, "preferred_name_b": name_b},
    )


def caller_identity() -> str:
    return os.environ.get("STRING_CALLER_IDENTITY", "kegg-string-mcp")


def _trace(resp) -> RequestTrace:
    return RequestTrace(url=resp.audit_url, retrieved_at=resp.fetched_at, cached=resp.cached,
                        status=resp.status, content_sha256=resp.content_sha256)


class StringClient:
    def __init__(self, http: PoliteClient):
        self.http = http

    def resolve(self, gene: str, species: int) -> tuple[dict[str, Any] | None, RequestTrace]:
        """Map a free-text identifier to a STRING protein ID via get_string_ids."""
        resp = self.http.get(
            f"{API}/json/get_string_ids",
            {"identifiers": gene, "species": species, "limit": 1, "echo_query": 1,
             "caller_identity": caller_identity()},
        )
        try:
            hits = json.loads(resp.body) if resp.body.strip() else []
        except json.JSONDecodeError:
            # STRING occasionally serves an HTML maintenance page with HTTP 200.
            hits = None
        # STRING's documented error shape is a JSON *object*
        # ({"Error": ..., "ErrorMessage": ...}), frequently with HTTP 200. That
        # decodes cleanly and is truthy, so checking only for a decode failure
        # left `hits[0]` raising KeyError straight out of the tool.
        # Checking only the container is not enough: a 200 decoding to
        # ["no such identifier"] passes an isinstance(list) test, and hit.get(...)
        # below then raises AttributeError out of the tool.
        if not isinstance(hits, list) or not hits or not isinstance(hits[0], dict):
            return None, _trace(resp)
        return hits[0], _trace(resp)

    def network(self, identifiers: list[str], species: int = MTB_H37RV,
                required_score: int = 0) -> ToolResult:
        """Every edge STRING holds *among* the identifiers given.

        `partners()` answers "who does this protein interact with, best first",
        which a `limit` truncates. Asking whether A and B interact by looking for B
        in A's top 20 is therefore only reliable when neither list is full: katG
        and rpoC score 0.823 with rpoC at rank 24 of katG's 31 partners, and the
        pair verdict said "No direct interaction" -- a confident negative produced
        by a pagination parameter.

        /network is not ranked and not truncated: it returns the edges between the
        proteins asked about, so a missing pair means "below required_score",
        which is a statement the caller can act on. One call answers every pair in
        a gene set, which is also why it is cheaper than the per-pair alternative.
        """
        query: dict[str, Any] = {"identifiers": list(identifiers), "species": species,
                                 "required_score": required_score}
        ids = [i.strip() for i in identifiers if i and i.strip()]
        problems = []
        if len(ids) < 2:
            problems.append("at least two identifiers are needed to ask for edges between them")
        if not MIN_SCORE <= required_score <= MAX_SCORE:
            problems.append(f"required_score={required_score} is outside STRING's range "
                            f"{MIN_SCORE}-{MAX_SCORE}")
        if species <= 0:
            problems.append(f"species={species} is not a valid NCBI taxon ID")
        if problems:
            return ToolResult.build(
                query, [], resolved={"matched_by": "none"},
                notes=[(f"Invalid argument(s), so no lookup was performed: {'; '.join(problems)}. "
                       f"An empty result here does NOT mean the proteins do not interact.")])

        records: list[Record] = []
        traces: list[RequestTrace] = []
        notes: list[str] = []
        seen: set[tuple[str, str]] = set()
        for batch in _batched(ids, NETWORK_BATCH):
            try:
                # STRING takes a list in one parameter, carriage-return separated.
                resp = self.http.get(
                    f"{API}/json/network",
                    {"identifiers": "\r".join(batch), "species": species,
                     "required_score": required_score, "caller_identity": caller_identity()},
                )
            except FetchError as exc:
                notes.append(f"STRING /network failed for {len(batch)} identifier(s): "
                             f"HTTP {exc.status}. Those pairs were NOT checked; this is a "
                             f"retrieval failure, not evidence of no interaction.")
                continue
            traces.append(_trace(resp))
            try:
                rows = json.loads(resp.body) if resp.body.strip() else []
            except json.JSONDecodeError:
                rows = None
            if not isinstance(rows, list) or any(not isinstance(r, dict) for r in rows):
                notes.append("STRING returned an unreadable or error response for /network "
                             "(expected a JSON list of edge objects). Those pairs were NOT checked.")
                continue
            for row in rows:
                record = _edge_record(row, species, resp)
                # A batch reports each undirected edge once, but overlapping
                # batches would repeat it; dedupe on the unordered ID pair so a
                # caller counting edges is not counting batches.
                key = tuple(sorted((row.get("stringId_A", ""), row.get("stringId_B", ""))))
                if key in seen:
                    continue
                seen.add(key)
                records.append(record)

        notes.append(
            f"STRING /network returns only edges scoring at or above required_score "
            f"({required_score}) among the {len(ids)} identifier(s) given. A pair that is "
            f"absent is below that threshold -- not unrelated, and not untested.")
        notes.append(
            "combined_score includes the textmining channel (literature co-mention); "
            f"evidence_beyond_textmining is True only when another channel reaches "
            f"STRING's medium-confidence threshold ({MEDIUM_CONFIDENCE}).")
        return ToolResult.build(
            query, records,
            resolved={"matched_by": "identifiers_as_given", "n_identifiers": len(ids)},
            requests=traces, notes=notes)

    def partners(
        self, gene: str, species: int = MTB_H37RV, limit: int = 20, required_score: int = 700
    ) -> ToolResult:
        query: dict[str, Any] = {"gene": gene, "species": species, "limit": limit,
                                 "required_score": required_score}
        traces: list[RequestTrace] = []
        notes: list[str] = []

        problems = []
        if not gene.strip():
            problems.append("no gene identifier was supplied")
        if not MIN_SCORE <= required_score <= MAX_SCORE:
            problems.append(f"required_score={required_score} is outside STRING's range "
                            f"{MIN_SCORE}-{MAX_SCORE} (700 = high confidence)")
        if not 1 <= limit <= MAX_LIMIT:
            problems.append(f"limit={limit} is outside 1-{MAX_LIMIT}")
        if species <= 0:
            problems.append(f"species={species} is not a valid NCBI taxon ID")
        if problems:
            return ToolResult.build(
                query, [], resolved={"matched_by": "none"},
                notes=[(f"Invalid argument(s), so no lookup was performed: {'; '.join(problems)}. "
                       f"An empty result here does NOT mean the gene has no partners.")],
            )

        try:
            hit, trace = self.resolve(gene, species)
        except FetchError as exc:
            return ToolResult.build(query, [], notes=[f"STRING identifier lookup failed: HTTP {exc.status}."])
        traces.append(trace)

        if hit is None:
            return ToolResult.build(
                query, [], resolved={"matched_by": "none"}, requests=traces,
                notes=[(f"'{gene}' did not resolve to a STRING protein in species {species} (or STRING "
                       f"returned an unreadable response). No partners were looked up. This is a "
                       f"resolution failure, not evidence of no partners.")],
            )

        string_id = hit.get("stringId", "")
        # STRING resolves fuzzily and synonym-matches. Say so when the protein it
        # picked is not literally what was asked for, rather than presenting a
        # best-guess match as if it were exact.
        preferred = hit.get("preferredName", "")
        # STRING protein IDs are "{taxon}.{locus}", so a locus-tag query is an exact
        # match even though preferredName is the gene symbol. Comparing only against
        # preferredName warned on correct input, which teaches the reader to ignore
        # the warning entirely.
        locus = string_id.split(".", 1)[-1]
        asked = gene.strip().upper()
        exact = asked in {preferred.upper(), locus.upper(), string_id.upper()}
        if preferred and not exact:
            notes.append(f"STRING resolved '{gene}' to '{preferred}' ({string_id}) by its own "
                         f"synonym matching, not by exact match. Verify this is the intended protein.")
        try:
            resp = self.http.get(
                f"{API}/json/interaction_partners",
                {"identifiers": string_id, "species": species, "limit": limit,
                 "required_score": required_score, "caller_identity": caller_identity()},
            )
        except FetchError as exc:
            return ToolResult.build(
                query, [], resolved={"string_id": string_id}, requests=traces,
                # notes + [...]: the sibling branch below was fixed for this and
                # this one was left behind, so a caller learned a request failed
                # without learning it was for a protein they never asked for.
                notes=notes + [f"STRING interaction_partners failed: HTTP {exc.status}."],
            )
        traces.append(_trace(resp))

        try:
            rows = json.loads(resp.body) if resp.body.strip() else []
        except json.JSONDecodeError:
            rows = None
        # Same reasoning as in resolve(): an error object decodes fine, and
        # iterating it would yield dict *keys* -- strings -- so `row.get(...)`
        # below would raise AttributeError out of the tool.
        if not isinstance(rows, list) or any(not isinstance(r, dict) for r in rows):
            return ToolResult.build(
                query, [], resolved={"string_id": string_id}, requests=traces,
                # notes + [...] rather than a fresh list: dropping the accumulated
                # notes would discard the synonym-match warning, so the caller would
                # not learn that string_id is not what they asked for.
                notes=notes + [(f"STRING returned an unreadable or error response for {string_id} "
                               f"(expected a JSON list of partner objects). No partner data retrieved.")],
            )
        if not rows:
            notes.append(f"STRING returned no partners for {string_id} at required_score>={required_score}. "
                         f"Lowering the threshold may return partners; this is not evidence of isolation.")

        records = []
        for row in rows:
            partner_id = row.get("stringId_B", "")
            records.append(
                Record(
                    record_id=partner_id,
                    type="partner",
                    name=row.get("preferredName_B", partner_id),
                    url=f"{NETWORK}{partner_id}",
                    source="string",
                    retrieved_at=resp.fetched_at,
                    cached=resp.cached,
                    detail=_scores(row) | {"partner_of": string_id},
                )
            )

        # Only name a partner as textmining-driven when textmining ACTUALLY carries it.
        # Testing `not evidence_beyond_textmining` alone also caught partners whose
        # tscore is 0 and whose support is spread across several sub-medium channels --
        # asserting literature support that the data does not show, which is exactly
        # the fabrication this module exists to prevent.
        textmining_only = [
            r.name for r in records
            if r.detail["textmining_score"] >= MEDIUM_CONFIDENCE
            and r.detail["max_non_textmining_score"] < MEDIUM_CONFIDENCE
        ]
        notes.append(
            "STRING's combined_score includes the textmining channel (literature co-mention). "
            f"evidence_beyond_textmining is True only when some other channel reaches "
            f"STRING's medium-confidence threshold ({MEDIUM_CONFIDENCE})."
        )
        if textmining_only:
            notes.append(
                "Supported essentially only by textmining, so NOT independent of literature "
                f"evidence about the same pair: {', '.join(textmining_only)}."
            )
        return ToolResult.build(
            query, records,
            resolved={"string_id": string_id, "preferred_name": hit.get("preferredName", ""),
                      "matched_by": "get_string_ids"},
            requests=traces, notes=notes,
        )
