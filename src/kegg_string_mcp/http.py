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
    # NCBI caps unauthenticated callers at 3 req/s and blocks over it. An API key
    # raises that to 10 req/s, but the floor stays put: the extra rate is worth
    # far less than never being the client that gets the user's IP banned.
    "eutils.ncbi.nlm.nih.gov": 0.35,
    # UniProt does not publish a hard rate limit; this is deliberately polite.
    "rest.uniprot.org": 0.4,
}
DEFAULT_MIN_INTERVAL = 1.0
RETRY_STATUSES = {429, 500, 502, 503, 504}
MAX_ATTEMPTS = 4
MAX_RETRY_AFTER = 60.0

# Params that identify the caller rather than the resource. They belong in the
# request but not in the cache key -- otherwise changing STRING_CALLER_IDENTITY
# silently invalidates the entire STRING cache and forces a refetch at 1 req/s.
# NCBI's tool/email/api_key are the same shape of thing: who is asking, not what
# is being asked for.
IDENTITY_PARAMS = {"caller_identity", "tool", "email", "api_key"}

# Credentials. Stripped from the cache key like every identity param, but ALSO
# redacted from the audit URL, which the others are not: a live fetch records the
# URL actually sent, and that URL travels into ToolResult.requests and from there
# into the run store on disk. An api_key is the one identity param that must not
# be written down anywhere, so it is scrubbed before provenance ever sees it.
SECRET_PARAMS = {"api_key"}

# Identifies a person rather than authenticating one. NCBI asks callers to send
# an email so it can contact whoever is hammering the service; that address then
# travelled into the audit URL on every live fetch, into ToolResult.requests,
# into the run store, into the demo runs committed to a public repository, and
# into the published demo page built from them. Redacted from the audit URL for
# the same reason as a credential: a provenance trail is written down and shared,
# and a URL a reader can check does not need the operator's address in it. The
# request still carries it -- NCBI's terms ask for that -- and it is already
# absent from the cache key via IDENTITY_PARAMS.
PERSONAL_PARAMS = {"email"}

# Everything redacted before provenance sees it, as opposed to merely dropped
# from the cache key.
REDACTED_IN_AUDIT = SECRET_PARAMS | PERSONAL_PARAMS

# Response headers kept as provenance. Deliberately a whitelist: response headers
# can carry cookies and tokens, and everything here lands in the run store on disk.
CAPTURED_HEADERS = {"x-uniprot-release", "x-uniprot-release-date"}
REDACTED = "REDACTED"


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
    def _urls(url: str, params: dict[str, object] | None) -> tuple[str, str, str]:
        """Return (request_url, cache_key, audit_url). Params are sorted so the key
        is stable regardless of the order the caller happened to build the dict in.

        `request_url` is what goes on the wire; `audit_url` is the same thing with
        credentials scrubbed, and is the only one of the two that is allowed to be
        recorded as provenance.
        """
        if not params:
            return url, url, url
        items = sorted((k, str(v)) for k, v in params.items() if v is not None)
        request_url = f"{url}?{urlencode(items)}"
        keyed = [(k, v) for k, v in items if k not in IDENTITY_PARAMS]
        cache_key = f"{url}?{urlencode(keyed)}" if keyed else url
        scrubbed = [(k, REDACTED if k in REDACTED_IN_AUDIT else v) for k, v in items]
        audit_url = f"{url}?{urlencode(scrubbed)}"
        return request_url, cache_key, audit_url

    def get(self, url: str, params: dict[str, object] | None = None) -> CachedResponse:
        """Fetch, or replay from cache."""
        request_url, cache_key, audit_url = self._urls(url, params)

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
                # audit_url, not request_url: FetchError's message is printed on an
                # uncaught traceback, and request_url still carries the api_key.
                raise FetchError(audit_url, response.status_code, response.text)

            # audit_url, not request_url: CachedResponse.request_url exists solely to
            # be reported as provenance, so a credential must never reach it.
            kept = {k: v for k, v in response.headers.items() if k.lower() in CAPTURED_HEADERS}
            return self.cache.put(cache_key, response.status_code, response.text,
                                  request_url=audit_url, headers=kept)

        raise FetchError(audit_url, 0, f"exhausted {MAX_ATTEMPTS} attempts: {last_error}")
