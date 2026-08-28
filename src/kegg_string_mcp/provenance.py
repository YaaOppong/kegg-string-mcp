"""Provenance types shared by both tools.

Every record a tool returns carries the four things the validation layer needs to
check a citation: a stable `record_id` the model can cite, a resolvable `url`, the
`source` it came from, and when it was actually retrieved.

`record_ids` on the envelope is the flat list the citation validator checks
against. Keeping it as a top-level field rather than making the validator walk
the record list means the check is a set membership test, which is hard to get
subtly wrong.
"""

from __future__ import annotations

import hashlib
from datetime import datetime, timezone
from typing import Any, Literal

from pydantic import BaseModel, Field

Source = Literal["kegg", "string", "pubmed", "uniprot"]


def sha256(data: str | bytes) -> str:
    if isinstance(data, str):
        data = data.encode("utf-8")
    return hashlib.sha256(data).hexdigest()


def utcnow() -> str:
    return datetime.now(timezone.utc).isoformat()


class Record(BaseModel):
    """One retrieved fact, citable by `record_id`."""

    record_id: str = Field(description="Stable ID within the source, e.g. 'mtu00360', '83332.Rv1908c', '35038342'.")
    type: str = Field(description="'pathway' | 'partner' | 'article'.")
    name: str
    url: str = Field(description="Resolvable URL for a human to check.")
    source: Source
    retrieved_at: str = Field(
        description="When the underlying HTTP response was FETCHED -- on a cache hit this is the "
        "original fetch time, not now. Otherwise every cached record would claim to be fresh."
    )
    cached: bool
    detail: dict[str, Any] = Field(
        default_factory=dict, description="Source-specific fields (scores, channel breakdown, ...)."
    )


class RequestTrace(BaseModel):
    """One HTTP request the tool made. The audit trail for how a result was obtained."""

    url: str
    retrieved_at: str
    cached: bool
    status: int
    content_sha256: str


class ToolResult(BaseModel):
    """The envelope every tool returns.

    Empty `records` is a legitimate, explicit answer -- see `notes` for why. A tool
    that returns nothing without saying so invites the model to fill the gap.
    """

    query: dict[str, Any]
    resolved: dict[str, Any] = Field(
        default_factory=dict, description="How the input identifier was resolved, and by what match."
    )
    records: list[Record] = Field(default_factory=list)
    record_ids: list[str] = Field(
        default_factory=list, description="Flat citable ID list -- what the citation validator checks against."
    )
    notes: list[str] = Field(
        default_factory=list, description="Anything the caller must not infer from silence."
    )
    requests: list[RequestTrace] = Field(default_factory=list)

    @classmethod
    def build(cls, query: dict[str, Any], records: list[Record], **kw: Any) -> ToolResult:
        return cls(query=query, records=records, record_ids=[r.record_id for r in records], **kw)
