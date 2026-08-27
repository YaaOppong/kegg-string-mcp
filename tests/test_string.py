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
