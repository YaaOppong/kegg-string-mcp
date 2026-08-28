"""Fixtures replay real KEGG, STRING and PubMed responses captured 2026-08-27/28.

Tests never hit the network. The fixtures are verbatim API output, so the parsers
are tested against the formats the services actually return -- which is how the
four-column `/list/{org}` layout was caught, and how the PubMed reference-list
DOIs were caught (one fixture article carries 35 ArticleIdList elements, of which
exactly one is its own).

The two efetch fixtures are the real responses for exactly the PMIDs their
paired esearch returned, so the search -> fetch path replays end to end.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from kegg_string_mcp.cache import CachedResponse
from kegg_string_mcp.kegg import KeggClient
from kegg_string_mcp.provenance import sha256
from kegg_string_mcp.pubmed import PubMedClient
from kegg_string_mcp.string_db import StringClient

FIXTURES = Path(__file__).parent / "fixtures"

# Longest match wins, so /list/pathway/ is not shadowed by /list/, and the
# single-PMID efetch is not shadowed by the multi-PMID one.
ROUTES = [
    ("rest.kegg.jp/list/pathway/mtu", "kegg_list_pathway_mtu.tsv"),
    ("rest.kegg.jp/link/pathway/mtu:Rv1908c", "kegg_link_pathway_Rv1908c.tsv"),
    ("rest.kegg.jp/list/mtu", "kegg_list_mtu_subset.tsv"),
    ("string-db.org/api/json/get_string_ids", "string_get_string_ids_katG.json"),
    ("string-db.org/api/json/interaction_partners", "string_interaction_partners_Rv1908c.json"),
    ("esearch.fcgi", "pubmed_esearch_katG_mtb.json"),
    ("efetch.fcgi", "pubmed_efetch_katG_mtb.xml"),
    ("id=8950806&retmode", "pubmed_efetch_no_abstract.xml"),
]


class FixtureHttp:
    """Stands in for PoliteClient. Records the URLs requested so tests can assert
    on call patterns (e.g. that the gene index is fetched once, not per gene)."""

    def __init__(self, fetched_at: str = "2026-08-27T09:00:00+00:00"):
        self.fetched_at = fetched_at
        self.calls: list[str] = []

    def get(self, url: str, params: dict | None = None) -> CachedResponse:
        from urllib.parse import urlencode

        if params:
            url = f"{url}?{urlencode(sorted((k, str(v)) for k, v in params.items() if v is not None))}"
        self.calls.append(url)
        for pattern, filename in sorted(ROUTES, key=lambda r: -len(r[0])):
            if pattern in url:
                body = (FIXTURES / filename).read_text(encoding="utf-8")
                return CachedResponse(url=url, status=200, body=body, fetched_at=self.fetched_at,
                                      content_sha256=sha256(body), cached=False)
        raise AssertionError(f"no fixture registered for {url}")


@pytest.fixture
def http() -> FixtureHttp:
    return FixtureHttp()


@pytest.fixture
def kegg(http: FixtureHttp) -> KeggClient:
    return KeggClient(http)


@pytest.fixture
def string(http: FixtureHttp) -> StringClient:
    return StringClient(http)


@pytest.fixture
def pubmed(http: FixtureHttp, monkeypatch) -> PubMedClient:
    # A developer's real NCBI_EMAIL/NCBI_API_KEY would otherwise leak into the
    # fixture URLs and change which route matches.
    for variable in ("NCBI_EMAIL", "NCBI_API_KEY", "NCBI_TOOL"):
        monkeypatch.delenv(variable, raising=False)
    return PubMedClient(http)
