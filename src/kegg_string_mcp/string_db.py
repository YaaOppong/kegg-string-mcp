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


def caller_identity() -> str:
    return os.environ.get("STRING_CALLER_IDENTITY", "kegg-string-mcp")


def _trace(resp) -> RequestTrace:
    return RequestTrace(url=resp.url, retrieved_at=resp.fetched_at, cached=resp.cached,
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
            # Letting that raise would surface as a hard tool error.
            hits = None
        return (hits[0] if hits else None), _trace(resp)

    def partners(
        self, gene: str, species: int = MTB_H37RV, limit: int = 20, required_score: int = 700
    ) -> ToolResult:
        query: dict[str, Any] = {"gene": gene, "species": species, "limit": limit,
                                 "required_score": required_score}
        traces: list[RequestTrace] = []
        notes: list[str] = []

        try:
            hit, trace = self.resolve(gene, species)
        except FetchError as exc:
            return ToolResult.build(query, [], notes=[f"STRING identifier lookup failed: HTTP {exc.status}."])
        traces.append(trace)

        if hit is None:
            return ToolResult.build(
                query, [], resolved={"matched_by": "none"}, requests=traces,
                notes=[f"'{gene}' did not resolve to a STRING protein in species {species} (or STRING "
                       f"returned an unreadable response). No partners were looked up. This is a "
                       f"resolution failure, not evidence of no partners."],
            )

        string_id = hit.get("stringId", "")
        # STRING resolves fuzzily and synonym-matches. Say so when the protein it
        # picked is not literally what was asked for, rather than presenting a
        # best-guess match as if it were exact.
        preferred = hit.get("preferredName", "")
        if preferred and preferred.upper() != gene.strip().upper() and gene.strip().upper() != string_id.upper():
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
                notes=[f"STRING interaction_partners failed: HTTP {exc.status}."],
            )
        traces.append(_trace(resp))

        try:
            rows = json.loads(resp.body) if resp.body.strip() else []
        except json.JSONDecodeError:
            return ToolResult.build(
                query, [], resolved={"string_id": string_id}, requests=traces,
                notes=[f"STRING returned an unreadable (non-JSON) response for {string_id}; it may be "
                       f"serving an error page. No partner data retrieved."],
            )
        if not rows:
            notes.append(f"STRING returned no partners for {string_id} at required_score>={required_score}. "
                         f"Lowering the threshold may return partners; this is not evidence of isolation.")

        records = []
        for row in rows:
            channels = {label: float(row.get(field, 0) or 0) for field, label in CHANNELS.items()}
            textmining = float(row.get(TEXTMINING, 0) or 0)
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
                    detail={
                        "combined_score": float(row.get("score", 0) or 0),
                        "channels": channels,
                        "textmining_score": textmining,
                        "max_non_textmining_score": max(channels.values(), default=0.0),
                        "evidence_beyond_textmining": max(channels.values(), default=0.0)
                        >= MEDIUM_CONFIDENCE,
                        "partner_of": string_id,
                    },
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
