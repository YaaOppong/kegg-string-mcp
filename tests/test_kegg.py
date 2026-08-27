def test_gene_index_reads_description_from_last_column(kegg):
    """/list/{org} returns 4 columns: id, type, location, description. Taking column 1
    would index every gene under the symbol 'CDS'."""
    index, _ = kegg.gene_index("mtu")
    assert "CDS" not in index
    assert index["KATG"] == "mtu:Rv1908c"


def test_symbol_and_locus_tag_both_resolve(kegg):
    index, _ = kegg.gene_index("mtu")
    assert index["KATG"] == index["RV1908C"] == "mtu:Rv1908c"


def test_pathways_returns_citable_records(kegg):
    result = kegg.pathways("katG", "mtu")
    assert result.resolved["kegg_gene_id"] == "mtu:Rv1908c"
    assert "mtu00360" in result.record_ids
    for record in result.records:
        assert record.url.startswith("https://www.kegg.jp/entry/")
        assert record.name and record.source == "kegg"


def test_record_ids_mirrors_records_exactly(kegg):
    """The citation validator checks record_ids; it must not drift from records."""
    result = kegg.pathways("katG", "mtu")
    assert result.record_ids == [r.record_id for r in result.records]


def test_pathway_names_are_resolved_not_left_as_ids(kegg):
    result = kegg.pathways("katG", "mtu")
    names = {r.record_id: r.name for r in result.records}
    assert "Metabolic pathways" in names["mtu01100"]


def test_unresolved_gene_is_reported_as_resolution_failure(kegg):
    """An empty result must not be readable as 'this gene has no pathways'."""
    result = kegg.pathways("NOT_A_REAL_GENE", "mtu")
    assert result.records == [] and result.resolved["matched_by"] == "none"
    assert "did not match" in result.notes[0]
    assert "not evidence" in result.notes[0]


def test_kegg_id_input_skips_the_gene_index(kegg, http):
    """Passing a full KEGG ID should not trigger the 4000-line gene list fetch."""
    kegg.pathways("mtu:Rv1908c", "mtu")
    assert not any(url.endswith("/list/mtu") for url in http.calls)


def test_qualified_kegg_id_overrides_the_organism_argument(kegg, http):
    """A fully-qualified ID carries its own organism. Ignoring it looked up pathway
    names in the wrong organism, producing correct IDs with no names and no note."""
    result = kegg.pathways("mtu:Rv1908c", organism="hsa")
    assert result.query["organism_used"] == "mtu"
    assert any("organism='hsa' was requested" in n for n in result.notes)
    # Names resolve because the name lookup followed the identifier, not the argument.
    assert all(r.name != "(name unavailable)" for r in result.records)
    assert any("/list/pathway/mtu" in url for url in http.calls)


def test_unnamed_pathways_are_never_silent(kegg):
    """If a name genuinely cannot be resolved, say so rather than emitting
    '(name unavailable)' with no explanation."""
    result = kegg.pathways("katG", "mtu")
    unnamed = [r.record_id for r in result.records if r.name == "(name unavailable)"]
    assert not unnamed or any("no name for" in n for n in result.notes)


def test_product_names_are_not_indexed_as_gene_symbols(kegg):
    """KEGG writes 'symbolA, symbolB; product name'. With no ';' there is no symbol
    field -- the cell is only a product name. Indexing it anyway made `toxin`,
    `hydrolase` and `pseudogene` resolve to arbitrary genes as exact matches."""
    index, _ = kegg.gene_index("mtu")
    for bogus in ("TOXIN", "ANTITOXIN", "HYDROLASE", "PSEUDOGENE", "HYPOTHETICAL PROTEIN"):
        assert bogus not in index, f"{bogus!r} should not be a symbol key"


def test_comma_in_a_product_name_does_not_create_fragment_keys(kegg):
    """'beta-1,3-glucanase' must not become the keys 'BETA-1' and '3-GLUCANASE'."""
    index, _ = kegg.gene_index("mtu")
    assert "BETA-1" not in index and "3-GLUCANASE" not in index


def test_bogus_symbol_is_reported_as_unresolved(kegg):
    result = kegg.pathways("toxin", "mtu")
    assert result.records == [] and result.resolved["matched_by"] == "none"
    assert "did not match" in result.notes[0]


def test_real_symbols_still_resolve(kegg):
    """The fix must not throw away genuine symbols along with the bogus ones."""
    index, _ = kegg.gene_index("mtu")
    assert index["KATG"] == "mtu:Rv1908c"
    assert index["DNAA"] == "mtu:Rv0001"


def test_locus_tag_is_never_shadowed_by_a_symbol(kegg):
    index, _ = kegg.gene_index("mtu")
    for locus in ("RV0001", "RV1908C", "RV0315"):
        assert index[locus].upper().endswith(locus)


def test_mixed_case_kegg_id_resolves(kegg):
    """'MTU:Rv1908c' is a valid identifier; it used to be reported as no match."""
    result = kegg.pathways("MTU:Rv1908c", "mtu")
    assert result.resolved["kegg_gene_id"] == "mtu:Rv1908c"
    assert result.record_ids


def test_pathway_identifier_is_rejected_not_treated_as_an_organism(kegg):
    """'path:mtu00360' matched the organism-code pattern and set organism='path'."""
    result = kegg.pathways("path:mtu00360", "mtu")
    assert result.records == [] and result.resolved["matched_by"] == "none"
    assert "not a gene" in result.notes[0]


def test_invalid_organism_returns_a_note_not_an_exception(kegg, http):
    """A bad organism code raised FetchError out of the tool."""
    import pytest
    from kegg_string_mcp.http import FetchError

    def boom(url, params=None):
        raise FetchError(url, 400, "")

    http.get = boom
    result = kegg.pathways("katG", "notanorg")
    assert result.records == []
    assert "may not be a valid KEGG organism code" in result.notes[0]


def test_qualified_id_with_no_pathways_does_not_assert_the_gene_exists(kegg, http):
    """KEGG answers /link identically (200, empty body) for an unknown gene, so
    'the gene exists in KEGG' was a fabricated claim on this path. 'mtu:katG' is
    the natural way to trigger it -- and it is not a valid KEGG gene ID."""
    from kegg_string_mcp.cache import CachedResponse

    original = http.get

    def get(url, params=None):
        if "/link/pathway/" in url:
            return CachedResponse(url=url, status=200, body="", fetched_at="2026-08-27T00:00:00+00:00",
                                  content_sha256="x", cached=False, request_url=url)
        return original(url, params)

    http.get = get
    result = kegg.pathways("mtu:NOTAGENE", "mtu")
    joined = " ".join(result.notes)
    assert "does NOT confirm" in joined
    assert "the gene exists in KEGG" not in joined.lower()


def test_missing_names_are_noted_even_when_the_name_list_is_empty(kegg, http):
    """The `and names` guard suppressed the note in exactly the case needing it:
    a successful but empty /list/pathway response."""
    from kegg_string_mcp.cache import CachedResponse

    original = http.get

    def get(url, params=None):
        if "/list/pathway/" in url:
            return CachedResponse(url=url, status=200, body="", fetched_at="2026-08-27T00:00:00+00:00",
                                  content_sha256="x", cached=False, request_url=url)
        return original(url, params)

    http.get = get
    result = kegg.pathways("katG", "mtu")
    assert result.records, "pathway IDs should still be returned"
    assert any("no name for" in n for n in result.notes), "silent (name unavailable) again"
