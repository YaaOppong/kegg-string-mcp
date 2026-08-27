"""Polite HTTP with disk caching, per-host rate limiting and bounded retries.

Both upstreams are free public services run by academic groups. KEGG asks for
modest request rates; STRING asks for roughly one request per second and for a
`caller_identity` on every call so they can see who is using the API. Honouring
both is not politeness theatre -- it is the difference between a tool that can be
run over a whole gene list and one that gets the caller's IP blocked.

Only successful responses are cached. Caching an error -- or an unfollowable
redirect, whose empty body reads downstream as a confident "no results" -- would
freeze a transient upstream problem into a wrong answer for the whole TTL.
"""

from __future__ import annotations

import os
import threading
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
MAX_RETRY_AFTER = 60.0

# Params that identify the caller rather than the resource. They belong in the
# request but not in the cache key -- otherwise changing STRING_CALLER_IDENTITY
# silently invalidates the entire STRING cache and forces a refetch at 1 req/s.
IDENTITY_PARAMS = {"caller_identity"}


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
        # MCP dispatches sync tool functions onto worker threads, so concurrent tool
        # calls genuinely race here. Without a lock, two callers read the same
        # `last`, sleep the same interval and fire together -- doubling the request
        # rate, which is precisely what this class exists to prevent.
        #
        # One lock PER HOST: a single shared lock is held across the sleep, so a
        # caller waiting on KEGG would also block an unrelated STRING call that had
        # no reason to wait.
        self._locks: dict[str, threading.Lock] = {}
        self._locks_guard = threading.Lock()

    def _host_lock(self, host: str) -> threading.Lock:
        with self._locks_guard:
            return self._locks.setdefault(host, threading.Lock())

    def _throttle(self, host: str) -> None:
        interval = MIN_INTERVAL.get(host, DEFAULT_MIN_INTERVAL)
        # Held across read, sleep and write: releasing before the sleep would let a
        # second thread through immediately.
        with self._host_lock(host):
            last = self._last_request.get(host)
            if last is not None:
                wait = interval - (time.monotonic() - last)
                if wait > 0:
                    self._sleep(wait)
            self._last_request[host] = time.monotonic()

    @staticmethod
    def _urls(url: str, params: dict[str, object] | None) -> tuple[str, str]:
        """Return (request_url, cache_key). Params are sorted so the key is stable
        regardless of the order the caller happened to build the dict in."""
        if not params:
            return url, url
        items = sorted((k, str(v)) for k, v in params.items() if v is not None)
        request_url = f"{url}?{urlencode(items)}"
        keyed = [(k, v) for k, v in items if k not in IDENTITY_PARAMS]
        cache_key = f"{url}?{urlencode(keyed)}" if keyed else url
        return request_url, cache_key

    def get(self, url: str, params: dict[str, object] | None = None) -> CachedResponse:
        """Fetch, or replay from cache."""
        request_url, cache_key = self._urls(url, params)

        hit = self.cache.get(cache_key)
        if hit is not None:
            return hit

        host = httpx.URL(request_url).host
        last_error: Exception | None = None

        for attempt in range(MAX_ATTEMPTS):
            final_attempt = attempt == MAX_ATTEMPTS - 1
            self._throttle(host)
            try:
                response = httpx.get(request_url, timeout=self.timeout,
                                     headers={"User-Agent": USER_AGENT}, follow_redirects=True)
            except httpx.HTTPError as exc:
                last_error = exc
                if not final_attempt:      # no point sleeping before giving up
                    self._sleep(2**attempt)
                continue

            if response.status_code in RETRY_STATUSES and not final_attempt:
                retry_after = response.headers.get("Retry-After")
                delay = float(retry_after) if retry_after and retry_after.isdigit() else 2**attempt
                # An uncapped Retry-After lets an upstream park a worker thread for
                # an hour, well past any client timeout.
                self._sleep(min(delay, MAX_RETRY_AFTER))
                continue

            if not response.is_success:
                # Covers 4xx, 5xx, and any 3xx that could not be followed. An
                # unfollowable redirect has an empty body, which downstream would
                # otherwise read as a genuine "no results".
                raise FetchError(request_url, response.status_code, response.text)

            return self.cache.put(cache_key, response.status_code, response.text,
                                  request_url=request_url)

        raise FetchError(request_url, 0, f"exhausted {MAX_ATTEMPTS} attempts: {last_error}")
