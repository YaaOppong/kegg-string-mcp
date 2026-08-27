"""Content-addressed disk cache for HTTP responses.

Two jobs. The obvious one is staying inside KEGG's and STRING's rate limits. The
less obvious one is reproducibility: a cached response replays byte-for-byte, so
a pipeline run over the same genes returns the same records, and the agent's
decisions can be re-examined later against exactly what it saw.

The cached `fetched_at` is returned unchanged on a hit. A record that says it was
retrieved now, when it actually came off disk from three weeks ago, is a false
provenance claim -- small, silent, and exactly the kind of thing this project
exists to not do.
"""

from __future__ import annotations

import json
import os
import time
from dataclasses import dataclass
from pathlib import Path

from kegg_string_mcp.provenance import sha256, utcnow

DEFAULT_TTL_SECONDS = 30 * 24 * 3600  # KEGG and STRING release infrequently.


def default_cache_dir() -> Path:
    return Path(os.environ.get("KEGG_STRING_MCP_CACHE", Path.home() / ".cache" / "kegg-string-mcp"))


@dataclass(frozen=True)
class CachedResponse:
    url: str
    status: int
    body: str
    fetched_at: str
    content_sha256: str
    cached: bool  # True if served from disk


class DiskCache:
    def __init__(self, root: Path | None = None, ttl_seconds: int | None = DEFAULT_TTL_SECONDS):
        self.root = Path(root) if root else default_cache_dir()
        self.ttl_seconds = ttl_seconds

    def _path(self, url: str) -> Path:
        digest = sha256(url)
        # Two-level fan-out keeps directory listings usable at tens of thousands of entries.
        return self.root / digest[:2] / f"{digest}.json"

    def get(self, url: str) -> CachedResponse | None:
        path = self._path(url)
        if not path.exists():
            return None
        if self.ttl_seconds is not None and (time.time() - path.stat().st_mtime) > self.ttl_seconds:
            return None
        payload = json.loads(path.read_text(encoding="utf-8"))
        return CachedResponse(
            url=payload["url"],
            status=payload["status"],
            body=payload["body"],
            fetched_at=payload["fetched_at"],
            content_sha256=payload["content_sha256"],
            cached=True,
        )

    def put(self, url: str, status: int, body: str) -> CachedResponse:
        response = CachedResponse(
            url=url,
            status=status,
            body=body,
            fetched_at=utcnow(),
            content_sha256=sha256(body),
            cached=False,
        )
        path = self._path(url)
        path.parent.mkdir(parents=True, exist_ok=True)
        # Write-then-rename: a crash mid-write must not leave a truncated cache entry
        # that later reads as a valid response.
        tmp = path.with_suffix(".tmp")
        tmp.write_text(
            json.dumps(
                {
                    "url": url,
                    "status": status,
                    "body": body,
                    "fetched_at": response.fetched_at,
                    "content_sha256": response.content_sha256,
                },
                indent=0,
            ),
            encoding="utf-8",
        )
        tmp.replace(path)
        return response
