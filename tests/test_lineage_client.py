"""Lineage-marker lookup: a stage-1 annotation alongside KEGG, STRING and UniProt."""

from __future__ import annotations

from types import SimpleNamespace

from kegg_string_mcp.lineage import (
    BARCODE_ORGANISM,
    GeneSpan,
    LineageClient,
    LineageSnp,
    annotate,
    is_marker,
    lineages,
    marks_position,
    parse_barcode,
    parse_gene_spans,
)

# Real rows from tbdb/barcode.bed, plus a header line that must not become a SNP.
BARCODE = """\
Chromosome\t1130\t1131\tlineage4.2.2.1\tA\tEuro-American\tLAM7-TUR\tNone
Chromosome\t4205\t4206\tlineage4.1.3\tT\tEuro-American\tT;X;H\tNone
Chromosome\t10726\t10727\tLa1.8\tG\tM.bovis\tNone\tNone
track name=barcode"""

# Rv0004 has no symbol; Rv0009 is on the reverse strand.
KEGG = """\
mtu:Rv0001\tCDS\t1..1524\tdnaA; chromosomal replication initiator protein DnaA
mtu:Rv0004\tCDS\t4000..5000\thypothetical protein
mtu:Rv0009\tCDS\tcomplement(10000..11000)\tppiA; peptidyl-prolyl cis-trans isomerase"""


def test_barcode_parses_and_skips_non_data_lines():
    snps = parse_barcode(BARCODE)
    assert [s.position for s in snps] == [1131, 4206, 10727]
    assert snps[0].lineage == "lineage4.2.2.1"
    assert snps[0].lineage_name == "Euro-American"
    assert snps[0].allele == "A"


def test_gene_spans_handle_complement_and_missing_symbols():
    spans = parse_gene_spans(KEGG)
    assert [(s.locus, s.symbol, s.start, s.end) for s in spans] == [
        ("Rv0001", "dnaA", 1, 1524),
        ("Rv0004", None, 4000, 5000),
        ("Rv0009", "ppiA", 10000, 11000),
    ]


def test_annotation_is_keyed_by_both_locus_tag_and_symbol():
    """Keying on the symbol alone returned "no marker" for every locus tag -- an
    informative-looking negative that was an indexing artefact."""
    index = annotate(parse_gene_spans(KEGG), parse_barcode(BARCODE))
    assert set(index) == {"Rv0001", "dnaA", "Rv0004", "Rv0009", "ppiA"}
    assert index["dnaA"] == index["Rv0001"]
    assert index["dnaA"][0].position == 1131


def test_snp_outside_every_gene_is_dropped_not_misassigned():
    spans = [GeneSpan("Rv0001", "dnaA", 1, 100)]
    snps = [LineageSnp(500, "lineage4", "Euro-American")]
    assert annotate(spans, snps) == {}


def test_marker_flag_is_yes_or_no():
    index = annotate(parse_gene_spans(KEGG), parse_barcode(BARCODE))
    assert is_marker("dnaA", index)
    assert is_marker("Rv0001", index)
    assert not is_marker("katG", index)          # absent from the index entirely


def test_lineages_are_deduplicated_and_sorted():
    index = {"a": [LineageSnp(1, "lineage4.6", "Euro-American"),
                   LineageSnp(2, "lineage2.2.2", "East-Asian"),
                   LineageSnp(3, "lineage4.6", "Euro-American")]}
    assert lineages("a", index) == ["lineage2.2.2", "lineage4.6"]
    assert lineages("absent", index) == []


def test_position_level_check_is_separate_from_the_gene_flag():
    """855 of 4,008 genes carry a barcode position, so the gene flag is coarse.
    Whether a specific hit is a lineage marker is a position question."""
    snps = [LineageSnp(1131, "lineage4.2.2.1", "Euro-American")]
    assert marks_position(snps, 1131).lineage == "lineage4.2.2.1"
    assert marks_position(snps, 1132) is None


class FakeHttp:
    """Stands in for PoliteClient, recording what was fetched."""

    def __init__(self):
        self.seen = []

    def get(self, url):
        self.seen.append(url)
        return SimpleNamespace(
            body=BARCODE if "barcode" in url else KEGG,
            audit_url=url, fetched_at="2026-01-01T00:00:00+00:00", cached=False,
            status=200, content_sha256="0" * 64)


def test_client_returns_citable_records_with_traces():
    result = LineageClient(FakeHttp()).markers("dnaA")
    assert [r.record_id for r in result.records] == ["tbdb:1131"]
    record = result.records[0]
    assert record.type == "lineage_marker"
    assert record.source == "tbdb"
    assert record.detail["lineage"] == "lineage4.2.2.1"
    assert record.detail["position"] == 1131
    # Structured, like a KEGG pathway ID: nothing to quote from.
    assert not getattr(record, "quotable_text", None)
    assert len(result.requests) == 2


def test_client_resolves_a_locus_tag_and_is_case_insensitive():
    client = LineageClient(FakeHttp())
    assert client.markers("Rv0001").resolved["matched_by"] == "locus_tag"
    assert client.markers("DNAA").resolved["matched_by"] == "symbol"


def test_absence_is_reported_as_a_finding_not_a_blank():
    result = LineageClient(FakeHttp()).markers("katG")
    assert result.records == []
    assert result.resolved["matched_by"] == "none"
    assert "No lineage-defining position" in result.notes[0]


def test_a_fetch_failure_is_not_reported_as_absence():
    """The dangerous failure: a network error rendered as 'not a lineage marker'."""
    from kegg_string_mcp.http import FetchError

    class Broken:
        def get(self, url):
            raise FetchError(url, 503)

    result = LineageClient(Broken()).markers("phoP")
    assert result.records == []
    assert "not evidence" in result.notes[0]
    assert "503" in result.notes[0]


def test_sources_are_fetched_once_per_client_not_once_per_gene():
    http = FakeHttp()
    client = LineageClient(http)
    for gene in ("dnaA", "ppiA", "katG", "Rv0004"):
        client.markers(gene)
    assert len(http.seen) == 2


def test_another_organism_is_refused_not_answered():
    """The barcode is a list of H37Rv coordinates and the file does not say so.
    Applied to E. coli spans it matches 789 of 4,639 genes by coordinate
    collision -- confident nonsense with citable record IDs. A refusal is the
    only safe answer."""
    http = FakeHttp()
    result = LineageClient(http).markers("recA", organism="eco")
    assert result.records == []
    assert http.seen == []                      # nothing was even fetched
    assert "cannot be applied" in result.notes[0]
    assert "not evidence" in result.notes[0]


def test_the_barcode_organism_is_accepted_in_any_case():
    client = LineageClient(FakeHttp())
    assert client.markers("dnaA", organism=BARCODE_ORGANISM.upper()).records


def test_a_position_in_overlapping_genes_is_recorded_against_both():
    """Overlapping reading frames are ordinary in a compact bacterial genome, and
    3 of the 1,111 barcode positions sit in two genes. Stopping at the first
    match made the second look unmarked, depending on KEGG's output order."""
    spans = [GeneSpan("Rv1009", "a", 100, 200), GeneSpan("Rv1010", "b", 150, 250)]
    index = annotate(spans, [LineageSnp(175, "lineage4", "Euro-American")])
    assert is_marker("Rv1009", index)
    assert is_marker("Rv1010", index)
    assert index["a"] == index["b"]
