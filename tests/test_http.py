"""PoliteClient behaviour. No real sockets: httpx.get is monkeypatched."""

import threading

import httpx
import pytest

from kegg_string_mcp.cache import DiskCache
import kegg_string_mcp.http as h
from kegg_string_mcp.http import MAX_RETRY_AFTER, FetchError, PoliteClient


class Resp:
    def __init__(self, status=200, text="ok", headers=None):
        self.status_code, self.text, self.headers = status, text, headers or {}

    @property
    def is_success(self):
        return 200 <= self.status_code < 300


def client(tmp_path, **kw):
    slept = []
    c = PoliteClient(DiskCache(tmp_path), sleep=slept.append, **kw)
    return c, slept


def test_identity_params_are_excluded_from_the_cache_key(tmp_path, monkeypatch):
    """Changing STRING_CALLER_IDENTITY used to invalidate the whole STRING cache."""
    calls = []
    monkeypatch.setattr(httpx, "get", lambda url, **kw: (calls.append(url), Resp())[1])
    c, _ = client(tmp_path)
    base = "https://string-db.org/api/json/x"
    c.get(base, {"identifiers": "katG", "caller_identity": "alice"})
    c.get(base, {"identifiers": "katG", "caller_identity": "bob"})
    assert len(calls) == 1, "second call should have hit the cache"
    assert "caller_identity=alice" in calls[0], "identity must still be sent upstream"


def test_param_order_does_not_change_the_cache_key(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(httpx, "get", lambda url, **kw: (calls.append(url), Resp())[1])
    c, _ = client(tmp_path)
    c.get("https://x.test/a", {"b": 2, "a": 1})
    c.get("https://x.test/a", {"a": 1, "b": 2})
    assert len(calls) == 1


def test_unfollowable_redirect_is_not_cached_as_success(tmp_path, monkeypatch):
    """A 3xx body is empty, and downstream reads an empty body as a genuine
    'no results' -- freezing a false negative on disk for the whole TTL."""
    monkeypatch.setattr(httpx, "get", lambda url, **kw: Resp(status=300, text=""))
    c, _ = client(tmp_path)
    with pytest.raises(FetchError):
        c.get("https://x.test/a")
    assert c.cache.get("https://x.test/a") is None


def test_error_responses_are_not_cached(tmp_path, monkeypatch):
    monkeypatch.setattr(httpx, "get", lambda url, **kw: Resp(status=404, text="nope"))
    c, _ = client(tmp_path)
    with pytest.raises(FetchError):
        c.get("https://x.test/a")
    assert c.cache.get("https://x.test/a") is None


def test_retry_after_is_capped(tmp_path, monkeypatch):
    """An uncapped Retry-After parks a worker thread far past any client timeout."""
    monkeypatch.setattr(httpx, "get",
                        lambda url, **kw: Resp(status=429, headers={"Retry-After": "3600"}))
    c, slept = client(tmp_path)
    with pytest.raises(FetchError):
        c.get("https://x.test/a")
    assert slept and max(slept) <= MAX_RETRY_AFTER


def test_no_backoff_sleep_after_the_final_attempt(tmp_path, monkeypatch):
    """The caller used to wait 8s for a failure that was already decided."""
    def boom(url, **kw):
        raise httpx.ConnectError("down")

    monkeypatch.setattr(httpx, "get", boom)
    c, slept = client(tmp_path)
    with pytest.raises(FetchError):
        c.get("https://x.test/a")
    # Throttle waits (fractional) interleave with backoffs (exact powers of two).
    backoffs = [s for s in slept if float(s).is_integer()]
    assert backoffs == [1, 2, 4], f"expected 3 backoffs and none after the last attempt, got {backoffs}"


def test_throttle_is_serialised_across_threads(tmp_path, monkeypatch):
    """MCP dispatches sync tools onto worker threads, so this genuinely races.
    Unlocked, both threads read the same timestamp and fire together."""
    monkeypatch.setattr(httpx, "get", lambda url, **kw: Resp())
    overlaps = []
    inside = threading.Event()

    def watching_sleep(seconds):
        if inside.is_set():
            overlaps.append(seconds)
        inside.set()
        inside.clear()

    c = PoliteClient(DiskCache(tmp_path), sleep=watching_sleep)
    c._last_request["x.test"] = -1e9  # force a throttle decision

    barrier = threading.Barrier(4)

    def hit(i):
        barrier.wait()
        c.get(f"https://x.test/{i}")

    threads = [threading.Thread(target=hit, args=(i,)) for i in range(4)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    assert not overlaps, "two threads were inside the throttle window at once"


def test_throttle_locks_are_per_host(tmp_path, monkeypatch):
    """A single shared lock is held across the sleep, so a caller waiting on KEGG
    would block an unrelated STRING call that had no reason to wait."""
    import time

    monkeypatch.setattr(httpx, "get", lambda url, **kw: Resp())
    c = PoliteClient(DiskCache(tmp_path), sleep=time.sleep)
    assert c._host_lock("rest.kegg.jp") is not c._host_lock("string-db.org")
    assert c._host_lock("rest.kegg.jp") is c._host_lock("rest.kegg.jp")


def test_concurrent_same_host_calls_are_serialised(tmp_path, monkeypatch):
    """Must distinguish locked from unlocked. Counting sleeps cannot: unlocked,
    every thread also reads a stale timestamp and sleeps. What differs is whether
    the sleeps overlap -- so measure elapsed time. Locked: N x interval.
    Unlocked: ~1 x interval, because they all sleep at once.
    """
    import time

    monkeypatch.setattr(httpx, "get", lambda url, **kw: Resp())
    monkeypatch.setitem(h.MIN_INTERVAL, "slow.test", 0.05)

    c = PoliteClient(DiskCache(tmp_path), sleep=time.sleep)
    c._last_request["slow.test"] = time.monotonic()   # force every caller to wait

    threads = [threading.Thread(target=c.get, args=(f"https://slow.test/{i}",)) for i in range(4)]
    start_time = time.monotonic()
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    elapsed = time.monotonic() - start_time

    # 4 waiters at 0.05s each: ~0.20s serialised, ~0.05s if they overlap.
    assert elapsed >= 0.15, f"throttle waits overlapped ({elapsed:.3f}s) -- lock not held across sleep"


def test_cache_hit_does_not_replay_another_callers_identity(tmp_path, monkeypatch):
    """The cache is shared. An entry written under caller_identity=alice must not
    report alice's URL as bob's audit trail, nor hand alice's identity to a model."""
    monkeypatch.setattr(httpx, "get", lambda url, **kw: Resp())
    c, _ = client(tmp_path)
    base = "https://string-db.org/api/json/x"

    live = c.get(base, {"identifiers": "katG", "caller_identity": "alice"})
    assert "alice" in live.audit_url, "a live fetch should record the URL actually sent"

    replay = c.get(base, {"identifiers": "katG", "caller_identity": "bob"})
    assert replay.cached
    assert "alice" not in replay.audit_url
    assert "caller_identity" not in replay.audit_url


def test_caller_identity_is_never_written_to_disk(tmp_path, monkeypatch):
    monkeypatch.setattr(httpx, "get", lambda url, **kw: Resp())
    c, _ = client(tmp_path)
    c.get("https://string-db.org/api/json/x", {"q": "katG", "caller_identity": "alice@example.com"})
    on_disk = "".join(p.read_text() for p in tmp_path.rglob("*.json"))
    assert "alice@example.com" not in on_disk
