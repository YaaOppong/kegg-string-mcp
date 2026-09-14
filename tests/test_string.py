def test_resolves_symbol_and_returns_partners(string):
    result = string.partners("katG")
    assert result.resolved["string_id"] == "83332.Rv1908c"
    assert "83332.Rv1909c" in result.record_ids


def test_channel_scores_are_returned_verbatim(string):
    result = string.partners("katG")
    fur_a = next(r for r in result.records if r.name == "furA")
    assert fur_a.detail["channels"]["neighborhood"] == 0.829
    assert fur_a.detail["textmining_score"] == 0.878
    assert fur_a.detail["combined_score"] == 0.979


def test_textmining_only_partners_are_flagged(string):
    """katG-embB scores 0.964, of which 0.963 is textmining and the only other
    non-zero channel is coexpression 0.044. That is not independent evidence."""
    by_name = {r.name: r for r in string.partners("katG").records}
    assert by_name["embB"].detail["evidence_beyond_textmining"] is False
    assert by_name["pncA"].detail["evidence_beyond_textmining"] is False
    # furA has neighborhood 0.829, sodA has database 0.5 -- both real support.
    assert by_name["furA"].detail["evidence_beyond_textmining"] is True
    assert by_name["sodA"].detail["evidence_beyond_textmining"] is True


def test_notes_name_the_textmining_only_partners(string):
    notes = " ".join(string.partners("katG").notes)
    assert "combined_score includes the textmining channel" in notes
    assert "embB" in notes and "pncA" in notes


def test_partner_urls_resolve_to_string_network_pages(string):
    for record in string.partners("katG").records:
        assert record.url == f"https://string-db.org/network/{record.record_id}"


def test_partners_with_no_textmining_are_not_called_textmining_driven(string, monkeypatch):
    """A partner whose support is spread across several sub-medium channels with
    tscore 0 was being named in the 'supported essentially only by textmining' note
    -- asserting literature support the data does not show."""
    import json

    rows = [{
        "stringId_A": "83332.Rv1908c", "stringId_B": "83332.Rv9999", "preferredName_A": "katG",
        "preferredName_B": "spreadEvidence", "ncbiTaxonId": 83332, "score": 0.72,
        "nscore": 0.25, "fscore": 0, "pscore": 0, "ascore": 0.3, "escore": 0.35,
        "dscore": 0, "tscore": 0.0,
    }]
    original = string.http.get

    def fake(url, params=None):
        resp = original(url, params)
        if "interaction_partners" in url:
            return type(resp)(url=resp.url, status=200, body=json.dumps(rows),
                              fetched_at=resp.fetched_at, content_sha256=resp.content_sha256,
                              cached=False)
        return resp

    string.http.get = fake
    result = string.partners("katG")
    record = result.records[0]
    assert record.detail["textmining_score"] == 0.0
    assert record.detail["evidence_beyond_textmining"] is False  # no channel reaches 0.4
    assert not any("spreadEvidence" in n for n in result.notes), \
        "a partner with tscore 0 must never be named as textmining-supported"


def test_unreadable_upstream_response_is_a_note_not_an_exception(string):
    """STRING serving an HTML error page with HTTP 200 raised JSONDecodeError."""

    def html(url, params=None):
        from kegg_string_mcp.cache import CachedResponse
        return CachedResponse(url=url, status=200, body="<html>maintenance</html>",
                              fetched_at="2026-08-27T00:00:00+00:00", content_sha256="x", cached=False)

    string.http.get = html
    result = string.partners("katG")
    assert result.records == []
    assert "unreadable" in " ".join(result.notes)


def test_locus_tag_query_is_not_flagged_as_a_synonym_match(string):
    """STRING IDs are '{taxon}.{locus}', so querying by locus tag is exact even
    though preferredName is the symbol. Warning on correct input teaches the
    reader to ignore the warning."""
    result = string.partners("Rv1908c")
    assert not any("synonym matching" in n for n in result.notes)


def test_symbol_query_is_not_flagged_either(string):
    assert not any("synonym matching" in n for n in string.partners("katG").notes)


def test_genuine_synonym_match_is_still_flagged(string):
    """An input that matches neither the symbol nor the locus must still warn."""
    result = string.partners("catalase-peroxidase")
    assert any("synonym matching" in n for n in result.notes)


def _json_http(payload):
    """STRING serving a JSON *object* — its documented error shape — with HTTP 200."""
    import json

    from kegg_string_mcp.cache import CachedResponse

    def get(url, params=None):
        return CachedResponse(url=url, status=200, body=json.dumps(payload),
                              fetched_at="2026-08-27T00:00:00+00:00", content_sha256="x",
                              cached=False, request_url=url)

    return get


def test_json_error_object_from_resolve_is_a_note_not_a_keyerror(string):
    """An error object decodes cleanly and is truthy, so hits[0] raised KeyError: 0."""
    string.http.get = _json_http({"Error": "not found", "ErrorMessage": "no such identifier"})
    result = string.partners("katG")
    assert result.records == []
    assert "resolution failure" in " ".join(result.notes)


def test_json_error_object_from_partners_is_a_note_not_an_attributeerror(string):
    """Iterating a dict yields its keys, so row.get(...) raised AttributeError."""
    original = string.http.get

    def get(url, params=None):
        if "interaction_partners" in url:
            return _json_http({"Error": "bad request"})(url, params)
        return original(url, params)

    string.http.get = get
    result = string.partners("katG")
    assert result.records == []
    assert "unreadable or error response" in " ".join(result.notes)


def test_json_list_of_non_objects_is_a_note_not_a_crash(string):
    """isinstance(x, list) admits ["no such identifier"], and .get() then raises."""
    original = string.http.get

    def get(url, params=None):
        if "interaction_partners" in url:
            return _json_http(["no such identifier"])(url, params)
        return original(url, params)

    string.http.get = get
    result = string.partners("katG")
    assert result.records == []
    assert "unreadable or error response" in " ".join(result.notes)


def test_synonym_warning_survives_a_later_error_response(string):
    """Building a fresh notes list discarded the warning that string_id is not
    what the caller asked for."""
    original = string.http.get

    def get(url, params=None):
        if "interaction_partners" in url:
            return _json_http({"Error": "boom"})(url, params)
        return original(url, params)

    string.http.get = get
    notes = " ".join(string.partners("catalase-peroxidase").notes)
    assert "synonym matching" in notes
    assert "unreadable or error response" in notes


def test_every_request_trace_carries_a_resolvable_url(string):
    result = string.partners("katG")
    assert result.requests
    for trace in result.requests:
        assert trace.url.startswith("http"), f"empty/invalid provenance URL: {trace.url!r}"
        assert trace.content_sha256 and trace.retrieved_at


def test_out_of_range_arguments_are_rejected_not_silently_empty(string):
    """STRING accepts required_score=99999 and returns zero partners, which would
    otherwise be reported as a genuine 'no partners' result."""
    for kwargs in ({"required_score": 99999}, {"required_score": -5},
                   {"limit": 0}, {"limit": -5}, {"species": -1}):
        result = string.partners("katG", **kwargs)
        assert result.records == [], kwargs
        joined = " ".join(result.notes)
        assert "Invalid argument" in joined, kwargs
        assert "does NOT mean the gene has no partners" in joined, kwargs


def test_valid_boundary_arguments_are_accepted(string):
    assert string.partners("katG", required_score=0).records
    assert string.partners("katG", required_score=1000).records
    assert string.partners("katG", limit=1).records


def test_synonym_warning_survives_a_fetch_error_too(string):
    """The JSON-error branch was fixed; the adjacent FetchError branch was not."""
    from kegg_string_mcp.http import FetchError

    original = string.http.get

    def get(url, params=None):
        if "interaction_partners" in url:
            raise FetchError(url, 500, "")
        return original(url, params)

    string.http.get = get
    notes = " ".join(string.partners("catalase-peroxidase").notes)
    assert "synonym matching" in notes
    assert "HTTP 500" in notes


# -- /network: the pair query -------------------------------------------------
#
# Served from the captured interaction_partners fixture rather than a new file:
# the two endpoints return the same row schema, and these tests are about what
# this module does with the rows, not about parsing a second format.

import json as _json
import pathlib as _pathlib

_ROWS = _json.loads((_pathlib.Path(__file__).parent / "fixtures"
                     / "string_interaction_partners_Rv1908c.json").read_text())


def _serving_network(string, rows):
    """Point /network at `rows`, leaving every other route on its fixture."""
    from kegg_string_mcp.cache import CachedResponse
    from kegg_string_mcp.provenance import sha256

    original = string.http.get

    def get(url, params=None):
        if "/json/network" in url:
            body = _json.dumps(rows)
            string.http.calls.append(url)
            return CachedResponse(url=url, status=200, body=body,
                                  fetched_at="2026-08-27T09:00:00+00:00",
                                  content_sha256=sha256(body), cached=False)
        return original(url, params)

    string.http.get = get
    return string


def test_network_returns_one_record_per_edge(string):
    result = _serving_network(string, _ROWS).network(["83332.Rv1908c", "83332.Rv1909c"])
    assert len(result.records) == len(_ROWS)
    edge = next(r for r in result.records if r.detail["preferred_name_b"] == "furA")
    assert edge.detail["combined_score"] == 0.979
    assert edge.detail["channels"]["neighborhood"] == 0.829
    assert edge.detail["evidence_beyond_textmining"] is True


def test_network_edge_ids_are_order_independent(string):
    """An undirected edge must have one ID whichever way STRING reports it, or a
    caller keying on it counts the same interaction twice."""
    flipped = [dict(row, stringId_A=row["stringId_B"], stringId_B=row["stringId_A"],
                    preferredName_A=row["preferredName_B"],
                    preferredName_B=row["preferredName_A"]) for row in _ROWS]
    forward = _serving_network(string, _ROWS).network(["83332.Rv1908c", "83332.Rv1909c"])
    backward = _serving_network(string, flipped).network(["83332.Rv1908c", "83332.Rv1909c"])
    assert sorted(r.record_id for r in forward.records) == \
           sorted(r.record_id for r in backward.records)


def test_network_says_that_an_absent_pair_is_below_threshold(string):
    """The whole point of the endpoint: a missing pair is a threshold statement,
    not the absence of a relationship."""
    notes = " ".join(_serving_network(string, []).network(["a", "b"], required_score=700).notes)
    assert "below that threshold" in notes
    assert "not untested" in notes


def test_network_needs_two_identifiers(string):
    result = string.network(["83332.Rv1908c"])
    assert result.records == []
    joined = " ".join(result.notes)
    assert "Invalid argument" in joined
    assert "does NOT mean the proteins do not interact" in joined


def test_network_failure_is_a_note_not_an_exception(string):
    from kegg_string_mcp.http import FetchError

    original = string.http.get

    def get(url, params=None):
        if "/json/network" in url:
            raise FetchError(url, 503, "")
        return original(url, params)

    string.http.get = get
    result = string.network(["83332.Rv1908c", "83332.Rv1909c"])
    joined = " ".join(result.notes)
    assert result.records == []
    assert "HTTP 503" in joined
    assert "not evidence of no interaction" in joined


def test_network_batches_long_identifier_lists(string):
    """A URL long enough for an intermediary to truncate returns fewer edges,
    which reads exactly like 'those pairs do not interact'."""
    from kegg_string_mcp.string_db import NETWORK_BATCH

    identifiers = [f"83332.G{i}" for i in range(NETWORK_BATCH * 2 + 1)]
    _serving_network(string, []).network(identifiers)
    calls = [c for c in string.http.calls if "/json/network" in c]
    assert len(calls) == 3
