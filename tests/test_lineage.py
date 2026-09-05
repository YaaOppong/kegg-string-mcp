"""Lineage annotation, and the nested/sibling distinction that decides its meaning."""

from __future__ import annotations

import pytest

from kegg_string_mcp.hypothesis.lineage import (
    GeneSpan,
    LineageSnp,
    annotate,
    flag_pair,
    marks_position,
    parse_barcode,
    parse_gene_spans,
    relate,
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


def test_annotation_keys_on_symbol_and_falls_back_to_locus():
    index = annotate(parse_gene_spans(KEGG), parse_barcode(BARCODE))
    assert set(index) == {"dnaA", "Rv0004", "ppiA"}
    assert index["dnaA"][0].position == 1131


def test_snp_outside_every_gene_is_dropped_not_misassigned():
    spans = [GeneSpan("Rv0001", "dnaA", 1, 100)]
    snps = [LineageSnp(500, "lineage4", "Euro-American")]
    assert annotate(spans, snps) == {}


def test_position_level_check_is_separate_from_the_gene_flag():
    """853 of 4,008 genes carry a barcode position, so the gene-level flag is a
    prior. Whether a specific hit is a lineage marker is a position question."""
    snps = [LineageSnp(1131, "lineage4.2.2.1", "Euro-American")]
    assert marks_position(snps, 1131).lineage == "lineage4.2.2.1"
    assert marks_position(snps, 1132) is None


@pytest.mark.parametrize("a,b,expected", [
    ("lineage4", "lineage4.6", "nested"),           # every 4.6 isolate is also lineage 4
    ("lineage4.6.1", "lineage4", "nested"),         # order must not matter
    ("lineage4.6.1", "lineage4.6.1", "nested"),     # same clade
    ("lineage4.6.1", "lineage4.6.3", "sibling"),    # sister clades; no isolate has both
    ("lineage4.1.3", "lineage4.2.2.1", "sibling"),
    ("lineage4.6.1", "lineage2.2.2", "unrelated"),
    ("lineage4", "La1.8", "unrelated"),             # M. bovis labels are not lineageN
])
def test_lineage_relations(a, b, expected):
    assert relate(a, b) == expected


def test_prefix_match_is_on_levels_not_characters():
    """String prefixing would call lineage4 an ancestor of lineage41."""
    assert relate("lineage4", "lineage41") == "unrelated"


def test_nested_markers_flag_a_positive_confound():
    index = {"a": [LineageSnp(1, "lineage4", "Euro-American")],
             "b": [LineageSnp(2, "lineage4.6", "Euro-American")]}
    flag = flag_pair("a", "b", index)
    assert flag.risk == "confounding_positive"
    assert "condition on lineage" in flag.note


def test_sibling_markers_flag_a_negative_confound():
    index = {"a": [LineageSnp(1, "lineage4.6.1", "Euro-American")],
             "b": [LineageSnp(2, "lineage4.6.3", "Euro-American")]}
    assert flag_pair("a", "b", index).risk == "confounding_negative"


def test_one_marked_gene_is_not_a_confound():
    index = {"a": [LineageSnp(1, "lineage4", "Euro-American")], "b": []}
    flag = flag_pair("a", "b", index)
    assert flag.risk == "none"
    assert flag.relations == []


def test_unrelated_clades_are_marked_but_not_flagged_as_confounding():
    index = {"a": [LineageSnp(1, "lineage4.6.1", "Euro-American")],
             "b": [LineageSnp(2, "lineage2.2.2", "East-Asian")]}
    assert flag_pair("a", "b", index).risk == "both_marked"


def test_nesting_wins_when_a_gene_carries_several_markers():
    """A gene with SNPs for two clades confounds if ANY pairing is nested."""
    index = {"a": [LineageSnp(1, "lineage2.2.2", "East-Asian"),
                   LineageSnp(2, "lineage4", "Euro-American")],
             "b": [LineageSnp(3, "lineage4.6", "Euro-American")]}
    flag = flag_pair("a", "b", index)
    assert flag.risk == "confounding_positive"
    assert ("lineage4", "lineage4.6", "nested") in flag.relations


def test_load_uses_the_caching_client_for_provenance():
    """A lineage call that cannot be traced to a barcode version is not evidence."""
    from types import SimpleNamespace

    from kegg_string_mcp.hypothesis import lineage as mod

    seen = []

    class FakeHttp:
        def get(self, url):
            seen.append(url)
            return SimpleNamespace(body=BARCODE if "barcode" in url else KEGG)

    index = mod.load(FakeHttp())
    assert seen == [mod.BARCODE_URL, mod.KEGG_GENE_LIST]
    assert "dnaA" in index
