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
