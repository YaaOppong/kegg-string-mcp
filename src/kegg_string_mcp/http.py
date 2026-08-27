"""Polite HTTP with disk caching, per-host rate limiting and bounded retries.

Both upstreams are free public services run by academic groups. KEGG asks for
modest request rates; STRING asks for roughly one request per second and for a
`caller_identity` on every call so they can see who is using the API. Honouring
both is not politeness theatre -- it is the difference between a tool that can be
run over a whole gene list and one that gets the caller's IP blocked.

Only 2xx responses are cached. Caching a 500 would freeze a transient upstream
failure into a permanent wrong answer.
"""

from __future__ import annotations

import os
import time
from urllib.parse import urlencode

import httpx

from kegg_string_mcp.cache import CachedResponse, DiskCache

USER_AGENT = os.environ.get("KEGG_STRING_MCP_USER_AGENT", "kegg-string-mcp/0.1 (+MCP tool server)")

# Conservative floors, in seconds between requests to the same host.
MIN_INTERVAL = {
    "rest.kegg.jp": 0.34,      # ~3 req/s
    "string-db.org": 1.05,     # STRING asks for ~1 req/s
}
DEFAULT_MIN_INTERVAL = 1.0
RETRY_STATUSES = {429, 500, 502, 503, 504}
MAX_ATTEMPTS = 4


class FetchError(RuntimeError):
    def __init__(self, url: str, status: int, body: str = ""):
        self.url, self.status, self.body = url, status, body
        super().__init__(f"GET {url} -> HTTP {status}: {body[:200]}")


class PoliteClient:
    def __init__(self, cache: DiskCache | None = None, timeout: float = 30.0, sleep=time.sleep):
        self.cache = cache or DiskCache()
        self.timeout = timeout
        self._sleep = sleep            # injectable so tests do not actually wait
        self._last_request: dict[str, float] = {}

    def _throttle(self, host: str) -> None:
        interval = MIN_INTERVAL.get(host, DEFAULT_MIN_INTERVAL)
        last = self._last_request.get(host)
        if last is not None:
            wait = interval - (time.monotonic() - last)
            if wait > 0:
                self._sleep(wait)
        self._last_request[host] = time.monotonic()

    def get(self, url: str, params: dict[str, object] | None = None) -> CachedResponse:
        """Fetch, or replay from cache. Params are sorted so the cache key is stable
        regardless of the order the caller happened to build the dict in."""
        if params:
            query = urlencode(sorted((k, str(v)) for k, v in params.items() if v is not None))
            url = f"{url}?{query}"

        hit = self.cache.get(url)
        if hit is not None:
            return hit

        host = httpx.URL(url).host
        last_error: Exception | None = None

        for attempt in range(MAX_ATTEMPTS):
            self._throttle(host)
            try:
                response = httpx.get(url, timeout=self.timeout, headers={"User-Agent": USER_AGENT},
                                     follow_redirects=True)
            except httpx.HTTPError as exc:
                last_error = exc
                self._sleep(2**attempt)
                continue

            if response.status_code in RETRY_STATUSES and attempt < MAX_ATTEMPTS - 1:
                retry_after = response.headers.get("Retry-After")
                delay = float(retry_after) if retry_after and retry_after.isdigit() else 2**attempt
                self._sleep(delay)
                continue

            if response.status_code >= 400:
                raise FetchError(url, response.status_code, response.text)

            return self.cache.put(url, response.status_code, response.text)

        raise FetchError(url, 0, f"exhausted {MAX_ATTEMPTS} attempts: {last_error}")
