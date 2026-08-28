import json
import time

from kegg_string_mcp.cache import DiskCache
from kegg_string_mcp.provenance import sha256

URL = "https://rest.kegg.jp/link/pathway/mtu:Rv1908c"


def test_roundtrip_preserves_body_and_hash(tmp_path):
    cache = DiskCache(tmp_path)
    written = cache.put(URL, 200, "mtu:Rv1908c\tpath:mtu00360\n")
    read = cache.get(URL)
    assert read.body == written.body
    assert read.content_sha256 == sha256(written.body)


def test_cache_hit_reports_the_original_fetch_time(tmp_path):
    """A record served from disk must not claim it was retrieved just now --
    that would be a false provenance claim on every cached result."""
    cache = DiskCache(tmp_path)
    written = cache.put(URL, 200, "body")
    time.sleep(0.01)
    read = cache.get(URL)
    assert read.fetched_at == written.fetched_at
    assert written.cached is False and read.cached is True


def test_miss_returns_none(tmp_path):
    assert DiskCache(tmp_path).get("https://example.org/never-fetched") is None


def test_expired_entry_is_a_miss(tmp_path):
    cache = DiskCache(tmp_path, ttl_seconds=0)
    cache.put(URL, 200, "body")
    time.sleep(0.01)
    assert cache.get(URL) is None


def test_ttl_none_never_expires(tmp_path):
    cache = DiskCache(tmp_path, ttl_seconds=None)
    cache.put(URL, 200, "body")
    assert cache.get(URL) is not None


def test_distinct_urls_do_not_collide(tmp_path):
    cache = DiskCache(tmp_path)
    cache.put(URL, 200, "first")
    cache.put(URL + "?x=1", 200, "second")
    assert cache.get(URL).body == "first"
    assert cache.get(URL + "?x=1").body == "second"


def test_partial_write_is_not_left_readable(tmp_path):
    """Write-then-rename: a .tmp file must never be picked up as a cache entry."""
    cache = DiskCache(tmp_path)
    cache.put(URL, 200, json.dumps({"ok": True}))
    assert not list(tmp_path.rglob("*.tmp"))


def test_corrupt_entry_is_a_miss_not_an_exception(tmp_path):
    """One corrupt file used to make every call for that URL raise for the whole
    30-day TTL, with no recovery short of deleting the cache by hand."""
    cache = DiskCache(tmp_path)
    cache.put(URL, 200, "body")
    entry = next(tmp_path.rglob("*.json"))
    entry.write_text("{ truncated json")
    assert cache.get(URL) is None


def test_entry_missing_a_required_key_is_a_miss(tmp_path):
    cache = DiskCache(tmp_path)
    cache.put(URL, 200, "body")
    entry = next(tmp_path.rglob("*.json"))
    entry.write_text(json.dumps({"url": URL, "status": 200}))
    assert cache.get(URL) is None


def test_concurrent_writers_do_not_corrupt_an_entry(tmp_path):
    """A fixed '.tmp' filename let two writers interleave into the same file and
    rename spliced JSON into place as a valid-looking entry."""
    import threading

    cache = DiskCache(tmp_path)
    bodies = [f"payload-{i}" * 500 for i in range(8)]
    threads = [threading.Thread(target=cache.put, args=(URL, 200, b)) for b in bodies]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    result = cache.get(URL)
    assert result is not None and result.body in bodies


def test_audit_url_falls_back_to_the_key_when_request_url_is_unset(tmp_path):
    from kegg_string_mcp.cache import CachedResponse

    direct = CachedResponse(url=URL, status=200, body="b", fetched_at="t",
                            content_sha256="h", cached=False)
    assert direct.audit_url == URL, "a live response must never report an empty URL"


def test_audit_url_prefers_the_request_url_on_a_live_fetch(tmp_path):
    from kegg_string_mcp.cache import CachedResponse

    live = CachedResponse(url=URL, status=200, body="b", fetched_at="t", content_sha256="h",
                          cached=False, request_url=URL + "?caller_identity=alice")
    assert "alice" in live.audit_url


def test_audit_url_uses_the_key_on_a_cache_hit(tmp_path):
    from kegg_string_mcp.cache import CachedResponse

    hit = CachedResponse(url=URL, status=200, body="b", fetched_at="t", content_sha256="h",
                         cached=True, request_url=URL + "?caller_identity=alice")
    assert hit.audit_url == URL and "alice" not in hit.audit_url
