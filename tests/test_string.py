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
    original_url = "https://string-db.org/api/json/get_string_ids"

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
