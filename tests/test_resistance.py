"""WHO-graded resistance variants: a stage 1 annotation like KEGG or UniProt."""

from __future__ import annotations

from types import SimpleNamespace

import pytest

from kegg_string_mcp.resistance import (
    UNGRADED,
    ResistanceClient,
    grading_counts,
    parse_catalogue,
)

# Real rows from tbdb/mutations.csv, including a quoted comment containing a
# comma -- splitting on commas shifts every field after it.
CSV = '''Gene,Mutation,type,drug,original_mutation,confidence,source,comment
katG,p.Ser315Thr,who_confidence,isoniazid,p.Ser315Thr,Assoc w R,WHO catalogue v2,
katG,p.Arg463Leu,who_confidence,isoniazid,p.Arg463Leu,Not assoc w R,WHO catalogue v2,
katG,p.Gly279Asp,who_confidence,isoniazid,p.Gly279Asp,Uncertain significance,WHO catalogue v2,
katG,p.Trp90Arg,who_confidence,isoniazid,p.Trp90Arg,,tbdb,
atpE,p.Ala63Pro,drug_resistance,bedaquiline,p.Ala63Pro,Assoc w R - Interim,WHO catalogue v2,"One site submitted only resistant strains, which may have inflated the PPV"
bacA,c.102G>A,who_confidence,amikacin,c.102G>A,Not assoc w R,WHO catalogue v2,
'''


class FakeHttp:
    def __init__(self, body=CSV):
        self.body, self.seen = body, []

    def get(self, url):
        self.seen.append(url)
        return SimpleNamespace(body=self.body, audit_url=url, cached=False, status=200,
                               fetched_at="2026-01-01T00:00:00+00:00", content_sha256="0" * 64)


def test_quoted_commas_in_the_comment_do_not_shift_fields():
    catalogue = parse_catalogue(CSV)
    variant = catalogue["atpE"][0]
    assert variant.confidence == "Assoc w R - Interim"
    assert variant.drug == "bedaquiline"
    assert "inflated the PPV" in variant.comment


def test_an_empty_confidence_is_ungraded_not_absent():
    """A third state: not association, not non-association, never graded."""
    ungraded = next(v for v in parse_catalogue(CSV)["katG"]
                    if v.mutation == "p.Trp90Arg")
    assert ungraded.confidence == UNGRADED
    assert not ungraded.associated


def test_both_association_tiers_count_as_associated():
    catalogue = parse_catalogue(CSV)
    assert catalogue["katG"][0].associated          # Assoc w R
    assert catalogue["atpE"][0].associated          # Assoc w R - Interim


def test_grading_counts_cover_every_tier():
    assert grading_counts(parse_catalogue(CSV)["katG"]) == {
        UNGRADED: 1, "Assoc w R": 1, "Not assoc w R": 1, "Uncertain significance": 1}


def test_one_associated_variant_marks_the_gene():
    """The asymmetric rule: three of katG's four listed variants are not
    associated, and the gene is still resistance-associated."""
    result = ResistanceClient(FakeHttp()).variants("katG")
    assert result.resolved["resistance_associated"] is True
    assert result.resolved["drugs"] == ["isoniazid"]
    assert result.resolved["variants_in_catalogue"] == 4
    assert "RESISTANCE-ASSOCIATED" in result.notes[0]


def test_only_associated_variants_are_returned_by_default():
    """Returning all of katG's real 1,771 rows would bury the answer in the
    1,254 graded 'Uncertain significance'."""
    result = ResistanceClient(FakeHttp()).variants("katG")
    assert [r.detail["mutation"] for r in result.records] == ["p.Ser315Thr"]


def test_filtering_by_drug_never_unmarks_the_gene():
    """The flag is computed over every variant, not the filtered subset."""
    result = ResistanceClient(FakeHttp()).variants("katG", drug="bedaquiline")
    assert result.records == []
    assert result.resolved["resistance_associated"] is True


def test_a_gene_assessed_and_negative_is_not_a_gene_never_assessed():
    assessed = ResistanceClient(FakeHttp()).variants("bacA")
    assert assessed.resolved["resistance_associated"] is False
    assert assessed.resolved["variants_in_catalogue"] == 1
    assert "was assessed and no variant met" in assessed.notes[0]

    absent = ResistanceClient(FakeHttp()).variants("sodA")
    assert absent.resolved["resistance_associated"] is False
    assert "not in the WHO catalogue" in absent.notes[0]
    assert "not a finding that the gene is unrelated" in absent.notes[0]


def test_record_ids_are_citable_and_name_gene_and_variant():
    record = ResistanceClient(FakeHttp()).variants("katG").records[0]
    assert record.record_id == "tbdb:katG:p.Ser315Thr"
    assert record.type == "resistance_variant"
    assert record.source == "tbdb"
    assert not getattr(record, "quotable_text", None)


def test_lookup_is_case_insensitive_on_the_gene():
    assert ResistanceClient(FakeHttp()).variants("KATG").records


def test_the_catalogue_is_fetched_once_per_client():
    """4.7 MB and 49,330 rows; re-parsing per gene would dominate every lookup."""
    http = FakeHttp()
    client = ResistanceClient(http)
    for gene in ("katG", "atpE", "bacA", "sodA"):
        client.variants(gene)
    assert len(http.seen) == 1


def test_a_fetch_failure_is_not_reported_as_no_resistance():
    from kegg_string_mcp.http import FetchError

    class Broken:
        def get(self, url):
            raise FetchError(url, 503)

    result = ResistanceClient(Broken()).variants("katG")
    assert result.records == []
    assert "not evidence" in result.notes[0]


@pytest.mark.parametrize("drug,expected", [("isoniazid", 1), ("bedaquiline", 0)])
def test_filtering_by_drug_narrows_the_associated_set(drug, expected):
    """Not the whole gene: every katG row is an isoniazid row, so filtering by
    drug alone must not return the 1,254 'Uncertain significance' variants that
    the default view exists to keep out."""
    result = ResistanceClient(FakeHttp()).variants("katG", drug=drug)
    assert len(result.records) == expected
    assert all(r.detail["associated"] for r in result.records)


def test_only_associated_variants_are_ever_returned():
    """The tool takes loci, not variants. No query shape surfaces a row graded
    'Not assoc w R' or 'Uncertain significance' as though it were a finding
    about a variant the caller asked about."""
    client = ResistanceClient(FakeHttp())
    for kwargs in ({}, {"drug": "isoniazid"}):
        result = client.variants("katG", **kwargs)
        assert all(r.detail["associated"] for r in result.records)
