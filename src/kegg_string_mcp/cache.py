"""Content-addressed disk cache for HTTP responses.

Two jobs. The obvious one is staying inside KEGG's and STRING's rate limits. The
less obvious one is reproducibility: a cached response replays byte-for-byte, so
a pipeline run over the same genes returns the same records, and an agent's
decisions can be re-examined later against exactly what it saw.

The cached `fetched_at` is returned unchanged on a hit. A record that says it was
retrieved now, when it actually came off disk from three weeks ago, is a false
provenance claim -- small, silent, and exactly the kind of thing this project
exists to not do.
"""

from __future__ import annotations

import json
import os
import tempfile
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
        """A damaged entry is a miss, never an exception.

        Without this, one corrupt file makes every future call for that URL raise
        for the whole 30-day TTL, with no recovery short of hand-deleting the cache.
        A miss just refetches.
        """
        path = self._path(url)
        try:
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
        except (OSError, ValueError, KeyError, TypeError):
            return None

    def put(self, url: str, status: int, body: str) -> CachedResponse:
        response = CachedResponse(
            url=url, status=status, body=body, fetched_at=utcnow(),
            content_sha256=sha256(body), cached=False,
        )
        path = self._path(url)
        path.parent.mkdir(parents=True, exist_ok=True)
        # Unique temp file per writer, then atomic rename. A fixed `.tmp` name would
        # let two concurrent writers (worker threads, or two processes sharing the
        # cache dir) interleave into the same file and rename spliced JSON into
        # place as a "valid" entry.
        fd, tmp_name = tempfile.mkstemp(dir=path.parent, suffix=".tmp")
        try:
            with os.fdopen(fd, "w", encoding="utf-8") as fh:
                json.dump(
                    {"url": url, "status": status, "body": body,
                     "fetched_at": response.fetched_at, "content_sha256": response.content_sha256},
                    fh,
                )
            os.replace(tmp_name, path)
        except BaseException:
            Path(tmp_name).unlink(missing_ok=True)
            raise
        return response
