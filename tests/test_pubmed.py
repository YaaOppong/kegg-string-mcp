"""PubMed client behaviour, replayed against verbatim E-utilities responses.

The fixture search is `"katG" AND "Mycobacterium tuberculosis"`, retmax 4, and the
efetch fixture is the real response for exactly the PMIDs it returned.
"""

import xml.etree.ElementTree as ET

import pytest

from kegg_string_mcp.cache import CachedResponse
from kegg_string_mcp.http import FetchError
from kegg_string_mcp.pubmed import MAX_LIMIT

FIXTURE_PMIDS = ["35919400", "30808448", "42239534", "35038342"]


def _response(body, url="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"):
    def get(_url, params=None):
        return CachedResponse(url=url, status=200, body=body,
                              fetched_at="2026-08-28T00:00:00+00:00", content_sha256="x",
                              cached=False, request_url=url)

    return get


def _only(pubmed, endpoint, body):
    """Replace one endpoint's response, leaving the other on its fixture."""
    original = pubmed.http.get

    def get(url, params=None):
        if endpoint in url:
            return _response(body, url)(url, params)
        return original(url, params)

    return get


# -- the happy path -------------------------------------------------------


def test_search_and_fetch_return_pmids_as_record_ids(pubmed):
    result = pubmed.abstracts("katG")
    assert result.record_ids == FIXTURE_PMIDS
    assert all(r.source == "pubmed" and r.type == "article" for r in result.records)


def test_records_keep_the_relevance_order_of_the_search(pubmed):
    """efetch does not preserve the ranking esearch returned, so the client
    reorders by the PMID list. Serving records in document order would silently
    change which article the model reads first."""
    assert [r.record_id for r in pubmed.abstracts("katG").records] == FIXTURE_PMIDS


def test_the_searched_term_and_pubmeds_translation_are_both_reported(pubmed):
    """The term this tool built and the term PubMed actually ran are different
    things, and only the second one explains the result set."""
    resolved = pubmed.abstracts("katG").resolved
    assert resolved["term"] == '"katG" AND "Mycobacterium tuberculosis"'
    assert resolved["query_translation"] == (
        '"katG"[All Fields] AND "Mycobacterium tuberculosis"[All Fields]'
    )


def test_article_urls_resolve_to_pubmed(pubmed):
    for record in pubmed.abstracts("katG").records:
        assert record.url == f"https://pubmed.ncbi.nlm.nih.gov/{record.record_id}/"


def test_every_request_trace_carries_a_resolvable_url(pubmed):
    result = pubmed.abstracts("katG")
    assert len(result.requests) == 2  # esearch + efetch
    for trace in result.requests:
        assert trace.url.startswith("http"), f"empty/invalid provenance URL: {trace.url!r}"
        assert trace.content_sha256 and trace.retrieved_at


# -- quotable_text, the span-grounding contract ---------------------------


def test_quotable_text_is_the_title_and_abstract_actually_retrieved(pubmed):
    record = next(r for r in pubmed.abstracts("katG").records if r.record_id == "35038342")
    assert record.detail["has_abstract"] is True
    assert record.detail["quotable_text"].startswith(record.detail["title"])
    assert record.detail["abstract"] in record.detail["quotable_text"]
    # The check a validator will actually run.
    quote = "Mycobacterium tuberculosis utilizes several mechanisms"
    assert quote in record.detail["quotable_text"]


def test_markup_inside_titles_and_abstracts_is_flattened(pubmed):
    """PubMed marks up titles and abstracts with <i> and <sup>. Reading .text
    instead of itertext() stops at the first child, which truncates an abstract at
    its first italicised species name and stores a fragment as the retrieved text."""
    record = next(r for r in pubmed.abstracts("katG").records if r.record_id == "35919400")
    text = record.detail["quotable_text"]
    assert "<" not in text and ">" not in text
    # The title's <i>Mycobacterium tuberculosis</i> survives as plain text...
    assert "Mycobacterium tuberculosis" in record.detail["title"]
    # ...and the abstract continues past its first markup child rather than stopping.
    assert len(record.detail["abstract"]) > 500


def test_title_markup_is_not_silently_truncated(pubmed):
    """Guards the specific .text-vs-itertext failure: the fixture title has an <i>
    child, so .text alone would end at 'Drug resistant '."""
    record = next(r for r in pubmed.abstracts("katG").records if r.record_id == "35919400")
    assert record.detail["title"] != "Drug resistant"
    assert record.detail["title"].endswith(".") or len(record.detail["title"]) > 40


def test_structured_abstract_labels_are_kept(pubmed):
    record = next(r for r in pubmed.abstracts("katG").records if r.record_id == "35919400")
    labels = [s["label"] for s in record.detail["abstract_sections"]]
    assert labels == ["BACKGROUND", "METHODS", "RESULTS", "CONCLUSIONS"]
    assert "BACKGROUND:" in record.detail["abstract"]


def test_unstructured_abstracts_get_no_invented_label(pubmed):
    record = next(r for r in pubmed.abstracts("katG").records if r.record_id == "35038342")
    assert [s["label"] for s in record.detail["abstract_sections"]] == [""]
    assert not record.detail["abstract"].startswith(":")


# -- articles with no abstract --------------------------------------------


def test_article_without_an_abstract_is_flagged_not_dropped(pubmed):
    """PMID 8950806 is a 1996 Trends piece PubMed holds without an abstract. It is
    still a real, citable record -- but nothing can be quoted from it beyond the
    title, and a validator has to be told that."""
    records, missing, _ = pubmed.fetch(["8950806"])
    assert missing == []
    record = records[0]
    assert record.detail["has_abstract"] is False
    assert record.detail["abstract"] == ""
    assert record.detail["quotable_text"] == record.detail["title"] == "Life without KatG."


def test_missing_abstracts_are_named_in_the_notes(pubmed):
    """Naming them, not counting them: a validator needs to know WHICH records
    have no abstract to check a quote against."""
    original = pubmed.http.get
    pubmed.http.get = _only(pubmed, "esearch.fcgi",
                            '{"esearchresult":{"count":"1","idlist":["8950806"],'
                            '"querytranslation":"x"}}')
    try:
        result = pubmed.abstracts("katG")
    finally:
        pubmed.http.get = original
    assert result.record_ids == ["8950806"]
    assert any("8950806" in n and "No abstract" in n for n in result.notes)


# -- the DOI scoping bug ---------------------------------------------------


def test_doi_is_the_articles_own_not_a_cited_papers(pubmed):
    by_pmid = {r.record_id: r for r in pubmed.abstracts("katG").records}
    assert by_pmid["35919400"].detail["doi"] == "10.7717/peerj.13645"
    assert by_pmid["42239534"].detail["doi"] == "10.3389/fcimb.2026.1844184"


def test_an_article_with_no_doi_does_not_borrow_one_from_its_references(pubmed):
    """The scoping bug this guards is invisible in the fixtures, because
    PubmedData/ArticleIdList precedes ReferenceList in document order, so taking
    the first DOI found anywhere happens to be right whenever the article HAS one.
    It stops being right here."""
    xml = """<PubmedArticleSet><PubmedArticle>
      <MedlineCitation><PMID>999</PMID><Article>
        <ArticleTitle>No DOI of its own</ArticleTitle>
        <Journal><Title>J Test</Title></Journal>
      </Article></MedlineCitation>
      <PubmedData><ArticleIdList>
        <ArticleId IdType="pubmed">999</ArticleId>
      </ArticleIdList><ReferenceList><Reference><ArticleIdList>
        <ArticleId IdType="doi">10.9999/not-this-articles-doi</ArticleId>
      </ArticleIdList></Reference></ReferenceList></PubmedData>
    </PubmedArticle></PubmedArticleSet>"""
    record = pubmed._record(ET.fromstring(xml).find("PubmedArticle"),
                            CachedResponse(url="u", status=200, body="", fetched_at="t",
                                           content_sha256="x", cached=False))
    assert record.detail["doi"] == "", "a cited paper's DOI was attributed to the article"


# -- the caveats the model has to see -------------------------------------


def test_search_is_always_labelled_as_search_not_resolution(pubmed):
    notes = " ".join(pubmed.abstracts("katG").notes)
    assert "NOT identifier resolution" in notes


def test_textmining_double_counting_is_flagged(pubmed):
    """STRING's combined_score already includes literature co-mention, so an
    abstract about the same pair may be the source of that score rather than
    independent support for it."""
    assert any("independent corroboration" in n for n in pubmed.abstracts("katG").notes)


def test_truncated_result_sets_say_how_much_was_left_behind(pubmed):
    """1382 articles matched and 4 were returned. Silence here reads as
    'this is the literature', which is the more damaging error."""
    notes = " ".join(pubmed.abstracts("katG").notes)
    assert "1382 articles matched" in notes
    assert "1378 not retrieved" in notes


def test_pmids_returned_by_search_but_not_by_fetch_are_reported(pubmed):
    """Book chapters come back as <PubmedBookArticle>, which has no MedlineCitation.
    Dropping them silently would leave a PMID the search found but that nothing
    can be cited from."""
    original = pubmed.http.get
    pubmed.http.get = _only(
        pubmed, "esearch.fcgi",
        '{"esearchresult":{"count":"5","idlist":["35919400","1"],"querytranslation":"x"}}')
    try:
        result = pubmed.abstracts("katG")
    finally:
        pubmed.http.get = original
    assert result.record_ids == ["35919400"]
    assert any("1" in n and "NOT citable" in n for n in result.notes)


# -- arguments -------------------------------------------------------------


@pytest.mark.parametrize("kwargs", [
    {"gene": ""},
    {"gene": "   "},
    {"gene": 'katG" AND cancer OR "'},
    {"gene": "katG[Title]"},
    {"gene": "katG", "organism": 'x" OR "y'},
    {"gene": "katG", "limit": 0},
    {"gene": "katG", "limit": -1},
    {"gene": "katG", "limit": MAX_LIMIT + 1},
])
def test_invalid_arguments_are_rejected_not_silently_empty(pubmed, kwargs):
    """An empty result from a malformed query reads as 'no literature exists',
    which is the fabrication this module is built to avoid."""
    result = pubmed.abstracts(**kwargs)
    assert result.records == []
    joined = " ".join(result.notes)
    assert "Invalid argument" in joined
    assert "does NOT mean there is no literature" in joined
    assert not pubmed.http.calls, "no request should be made for a rejected argument"


def test_query_syntax_injection_cannot_reach_pubmed(pubmed):
    """`katG" AND cancer OR "` would close the phrase quote and change the question
    asked, while still returning a plausible-looking result set."""
    assert pubmed.abstracts('katG" AND cancer OR "').records == []
    assert not pubmed.http.calls


def test_boundary_limits_are_accepted(pubmed):
    assert pubmed.abstracts("katG", limit=1).records
    assert pubmed.abstracts("katG", limit=MAX_LIMIT).records


def test_empty_organism_searches_without_organism_context(pubmed):
    assert pubmed.abstracts("katG", organism="").resolved["term"] == '"katG"'


# -- upstream failures -----------------------------------------------------


def test_html_error_page_from_esearch_is_a_note_not_an_exception(pubmed):
    """NCBI serves an HTML holding page with HTTP 200 when overloaded."""
    pubmed.http.get = _response("<html>service unavailable</html>")
    result = pubmed.abstracts("katG")
    assert result.records == []
    assert "unreadable or error response" in " ".join(result.notes)


def test_esearch_error_object_is_reported_with_ncbis_message(pubmed):
    pubmed.http.get = _response('{"esearchresult":{"ERROR":"Invalid db name"}}')
    result = pubmed.abstracts("katG")
    assert result.records == []
    joined = " ".join(result.notes)
    assert "Invalid db name" in joined
    assert "not evidence that no literature exists" in joined


def test_zero_results_is_an_explicit_note(pubmed):
    pubmed.http.get = _response('{"esearchresult":{"count":"0","idlist":[],'
                                '"querytranslation":"nope"}}')
    result = pubmed.abstracts("katG")
    assert result.records == []
    joined = " ".join(result.notes)
    assert "no articles" in joined
    assert "search miss" in joined


def test_non_numeric_ids_from_esearch_are_dropped_before_reaching_a_url(pubmed):
    """The same reasoning as kegg's organism validation: an upstream that returns
    something unexpected must not get to shape the next request."""
    pmids, _, _ = pubmed.search('"katG"', 4)
    assert pmids == FIXTURE_PMIDS
    pubmed.http.get = _response('{"esearchresult":{"count":"1",'
                                '"idlist":["../../etc/passwd"],"querytranslation":"x"}}')
    pmids, _, _ = pubmed.search('"katG"', 4)
    assert pmids == []


def test_unparseable_efetch_lists_the_pmids_as_not_citable(pubmed):
    """The search succeeded, so the PMIDs are real -- but nothing was retrieved for
    them, and a record_id with no retrieved content is exactly what must not become
    citable."""
    original = pubmed.http.get
    pubmed.http.get = _only(pubmed, "efetch.fcgi", "<html>nope</html>")
    try:
        result = pubmed.abstracts("katG")
    finally:
        pubmed.http.get = original
    assert result.records == [] and result.record_ids == []
    joined = " ".join(result.notes)
    assert "NOT citable" in joined
    assert "35919400" in joined


def test_efetch_fetch_error_keeps_the_pmids_out_of_the_citable_set(pubmed):
    original = pubmed.http.get

    def get(url, params=None):
        if "efetch.fcgi" in url:
            raise FetchError(url, 500, "")
        return original(url, params)

    pubmed.http.get = get
    result = pubmed.abstracts("katG")
    assert result.record_ids == []
    joined = " ".join(result.notes)
    assert "HTTP 500" in joined and "NOT citable" in joined


def test_esearch_fetch_error_is_a_note_not_an_exception(pubmed):
    def get(url, params=None):
        raise FetchError(url, 503, "")

    pubmed.http.get = get
    result = pubmed.abstracts("katG")
    assert result.records == []
    assert "HTTP 503" in " ".join(result.notes)
