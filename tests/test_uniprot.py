"""UniProt client. Fixtures are verbatim API responses captured 2026-08-28."""

from pathlib import Path

import pytest

from kegg_string_mcp.cache import CachedResponse
from kegg_string_mcp.uniprot import EVIDENCE_TIERS, UniProtClient

FIXTURES = Path(__file__).parent / "fixtures"


class FixtureHttp:
    def __init__(self, name: str, headers: dict | None = None):
        self.name, self.headers = name, headers or {"x-uniprot-release": "2026_02"}
        self.calls: list[str] = []

    def get(self, url, params=None):
        self.calls.append(url)
        body = (FIXTURES / f"{self.name}.json").read_text(encoding="utf-8")
        return CachedResponse(url=url, status=200, body=body,
                              fetched_at="2026-08-28T00:00:00+00:00", content_sha256="x",
                              cached=False, request_url=url, headers=self.headers)


def client(name: str) -> UniProtClient:
    return UniProtClient(FixtureHttp(name))


# --- records ---------------------------------------------------------------


def test_returns_a_citable_accession_and_resolvable_url():
    result = client("uniprot_search_Rv1908c").protein("Rv1908c")
    record = result.records[0]
    assert record.record_id == "P9WIE5"
    assert record.url == "https://www.uniprot.org/uniprotkb/P9WIE5"
    assert result.record_ids == [r.record_id for r in result.records]


def test_covers_a_gene_kegg_has_no_pathway_for():
    """gyrA is a KEGG negative control -- zero pathways. UniProt annotates it fully,
    which is the entire reason this tool exists."""
    record = client("uniprot_search_Rv0006").protein("Rv0006").records[0]
    assert record.record_id == "P9WG47"
    assert "gyrase" in record.name.lower()
    assert record.detail["has_experimental_function"]
    assert record.detail["pdb"]


def test_function_text_is_stored_for_span_checking():
    """Function statements are prose, so a claim drawn from one needs the same
    quote check as a PubMed abstract."""
    record = client("uniprot_search_Rv1908c").protein("Rv1908c").records[0]
    quotable = record.detail["quotable_text"]
    assert len(quotable) > 200
    assert record.detail["function_statements"][0]["text"] in quotable


# --- evidence tiers --------------------------------------------------------


def test_experimental_statements_carry_their_supporting_pmids():
    """The link that makes a function claim traceable to a paper -- and those PMIDs
    are already citable record IDs elsewhere in this pipeline."""
    record = client("uniprot_search_Rv1908c").protein("Rv1908c").records[0]
    experimental = [f for f in record.detail["function_statements"] if f["experimental"]]
    assert experimental
    assert any(f["supporting_pmids"] for f in experimental)


def test_inferred_statements_are_tiered_and_caveated():
    """A statement inferred from a HAMAP rule is not evidence about this protein."""
    result = client("uniprot_search_Rv1908c").protein("Rv1908c")
    tiers = {t for f in result.records[0].detail["function_statements"] for t in f["tiers"]}
    assert "sequence_model" in tiers
    assert any("INFERRED" in n for n in result.notes)


def test_evidence_code_map_covers_the_experimental_code():
    assert EVIDENCE_TIERS["ECO:0000269"] == "experimental"
    assert EVIDENCE_TIERS["ECO:0000255"] == "sequence_model"


def test_release_is_recorded_as_provenance():
    result = client("uniprot_search_Rv1908c").protein("Rv1908c")
    assert result.resolved["uniprot_release"] == "2026_02"


def test_release_is_recorded_even_when_nothing_matched():
    result = client("uniprot_search_nohit").protein("NOTAGENE")
    assert result.resolved["uniprot_release"] == "2026_02"


# --- failure modes ---------------------------------------------------------


def test_no_match_is_a_resolution_failure_not_an_absence_of_annotation():
    result = client("uniprot_search_nohit").protein("NOTAGENE")
    assert result.records == []
    assert "resolution failure" in " ".join(result.notes)
    assert "not evidence that the protein is unannotated" in " ".join(result.notes)


@pytest.mark.parametrize("gene", ["", "  ", "kat G", "gene:katG", 'katG"', "kat(G)"])
def test_query_syntax_and_empty_genes_are_rejected_before_any_request(gene):
    """A value containing UniProt query syntax would change what is asked."""
    http = FixtureHttp("uniprot_search_nohit")
    result = UniProtClient(http).protein(gene)
    assert result.records == [] and "Invalid argument" in result.notes[0]
    assert http.calls == []


@pytest.mark.parametrize("kwargs", [{"organism_id": 0}, {"organism_id": -1},
                                    {"limit": 0}, {"limit": 99}])
def test_out_of_range_arguments_are_rejected(kwargs):
    result = client("uniprot_search_nohit").protein("katG", **kwargs)
    assert result.records == []
    assert "does NOT mean the protein is unannotated" in " ".join(result.notes)


def test_unreadable_response_is_a_note_not_a_crash():
    class Html(FixtureHttp):
        def get(self, url, params=None):
            return CachedResponse(url=url, status=200, body="<html>down</html>",
                                  fetched_at="t", content_sha256="x", cached=False,
                                  request_url=url, headers={})

    result = UniProtClient(Html("x")).protein("katG")
    assert result.records == [] and "unreadable" in " ".join(result.notes)


def test_unparseable_entries_are_not_a_silent_empty_success():
    """Entries returned but none with a valid accession produced zero records and an
    empty notes list -- indistinguishable from 'this protein is unannotated'."""
    import json

    class Malformed(FixtureHttp):
        def get(self, url, params=None):
            body = json.dumps({"results": [{"primaryAccession": "not-an-accession"}]})
            return CachedResponse(url=url, status=200, body=body, fetched_at="t",
                                  content_sha256="x", cached=False, request_url=url,
                                  headers={"x-uniprot-release": "2026_02"})

    result = UniProtClient(Malformed("x")).protein("katG")
    assert result.records == []
    joined = " ".join(result.notes)
    assert "none could be parsed" in joined
    assert "not evidence that the protein is unannotated" in joined
    assert result.resolved["matched_by"] == "none"


def test_multiple_matching_entries_are_disclosed():
    """resolved.accession names only the first; reporting it as 'the' answer would
    hide that a choice was made between paralogues."""
    import json

    body = json.loads((FIXTURES / "uniprot_search_Rv0006.json").read_text())
    entry = body["results"][0]
    second = dict(entry, primaryAccession="A0A123XYZ9")

    class Two(FixtureHttp):
        def get(self, url, params=None):
            return CachedResponse(url=url, status=200,
                                  body=json.dumps({"results": [entry, second]}),
                                  fetched_at="t", content_sha256="x", cached=False,
                                  request_url=url, headers={"x-uniprot-release": "2026_02"})

    result = UniProtClient(Two("x")).protein("gyrA")
    assert len(result.records) == 2
    assert "2 UniProt entries matched" in " ".join(result.notes)


def test_single_entry_gets_no_ambiguity_note():
    result = client("uniprot_search_Rv1908c").protein("Rv1908c")
    assert not any("entries matched" in n for n in result.notes)
