"""Store, evidence assembly and citation validation. No model, no network."""

from pathlib import Path

from kegg_string_mcp.agent.evidence import all_pairs, classify_pathway, pair_evidence
from kegg_string_mcp.agent.store import RunStore
from kegg_string_mcp.agent.validate import extract_citations, validate

KEGG_RESULT = {
    "resolved": {"kegg_gene_id": "mtu:Rv1908c", "matched_by": "symbol"},
    "records": [{"record_id": "mtu00360", "name": "Phenylalanine metabolism"},
                {"record_id": "mtu01100", "name": "Metabolic pathways"}],
    "record_ids": ["mtu00360", "mtu01100"],
}
STRING_RESULT = {
    "resolved": {"string_id": "83332.Rv1908c", "preferred_name": "katG"},
    "records": [{"record_id": "83332.Rv1909c", "name": "furA",
                 "detail": {"combined_score": 0.979, "textmining_score": 0.878,
                            "max_non_textmining_score": 0.829,
                            "evidence_beyond_textmining": True}}],
    "record_ids": ["83332.Rv1909c"],
}


def _store(tmp_path: Path) -> RunStore:
    store = RunStore(path=tmp_path / "run.jsonl", run_id="test")
    store.tool_result("kegg_pathways", {"gene": "katG"}, KEGG_RESULT)
    store.tool_result("string_partners", {"gene": "katG"}, STRING_RESULT)
    return store


# --- store -----------------------------------------------------------------


def test_resolved_identifiers_are_citable(tmp_path):
    """Citing the subject of your own query is correct, not a hallucination.
    A validator that flags it teaches the reader to ignore every flag."""
    store = _store(tmp_path)
    assert "mtu:Rv1908c" in store.citable_ids
    assert "83332.Rv1908c" in store.citable_ids


def test_records_are_citable_and_attributed_to_their_gene(tmp_path):
    store = _store(tmp_path)
    assert {"mtu00360", "mtu01100", "83332.Rv1909c"} <= store.citable_ids
    assert "mtu00360" in store.per_target["KATG"]


def test_store_is_append_only_and_replayable(tmp_path):
    store = _store(tmp_path)
    store.decision({"turn": 1, "stop_reason": "tool_use"})
    kinds = [entry["kind"] for entry in store.replay()]
    assert kinds == ["tool_result", "tool_result", "decision"]


# --- validation ------------------------------------------------------------


def test_unsupported_identifier_is_flagged(tmp_path):
    store = _store(tmp_path)
    report = validate("katG is in mtu00360 and also mtu00627.", store.citable_ids)
    assert [c.identifier for c in report.unsupported] == ["mtu00627"]
    assert not report.passed


def test_clean_summary_passes(tmp_path):
    store = _store(tmp_path)
    report = validate("katG (mtu:Rv1908c) is in mtu00360 and partners with 83332.Rv1909c.",
                      store.citable_ids)
    assert report.passed, report.summary_line()


def test_cross_target_citation_is_distinguished_from_a_hallucination():
    """An ID that was retrieved, but for a different gene. Plausible-looking, and
    invisible to a validator that only checks global membership."""
    report = validate("katG is in mtu99999 and mtu00777.",
                      citable_ids={"mtu00777", "mtu00360"},
                      per_target={"KATG": {"mtu00360"}, "FURA": {"mtu00777"}},
                      claimed_target="KATG")
    statuses = {c.identifier: c.status for c in report.citations}
    assert statuses["mtu99999"] == "unsupported"
    assert statuses["mtu00777"] == "cross_target"


def test_extracts_both_identifier_shapes():
    found = extract_citations("see mtu00360 and 83332.Rv1909c, plus noise like 42 and ABC")
    assert found == ["mtu00360", "83332.Rv1909c"]


# --- evidence --------------------------------------------------------------


def test_broad_pathway_is_not_treated_as_a_link():
    """mtu01100 holds ~17% of the M. tuberculosis genome."""
    assert classify_pathway(698, 4008)[0] == "broad"
    assert classify_pathway(11, 4008)[0] == "specific"


def test_pair_sharing_only_a_broad_pathway_gets_a_negative_verdict():
    ev = pair_evidence(
        "geneA", "geneB",
        pathways={"geneA": [{"record_id": "mtu01100", "name": "Metabolic pathways"}],
                  "geneB": [{"record_id": "mtu01100", "name": "Metabolic pathways"}]},
        partners={"geneA": [], "geneB": []},
        pathway_sizes={"mtu01100": 698}, genome_size=4008,
    )
    assert "not evidence of a mechanistic link" in ev.verdict


def test_pair_sharing_a_specific_pathway_is_reported_as_a_link():
    ev = pair_evidence(
        "geneA", "geneB",
        pathways={"geneA": [{"record_id": "mtu00983", "name": "Drug metabolism"}],
                  "geneB": [{"record_id": "mtu00983", "name": "Drug metabolism"}]},
        partners={"geneA": [], "geneB": []},
        pathway_sizes={"mtu00983": 11}, genome_size=4008,
    )
    assert "Share specific pathway" in ev.verdict


def test_direct_interaction_is_detected_and_its_support_characterised():
    ev = pair_evidence(
        "katG", "furA",
        pathways={}, pathway_sizes={}, genome_size=4008,
        partners={"katG": STRING_RESULT["records"], "furA": []},
    )
    assert ev.direct_interaction is not None
    assert "supported beyond literature co-mention" in ev.verdict


def test_no_evidence_is_an_explicit_verdict_not_silence():
    ev = pair_evidence("geneA", "geneB", pathways={}, partners={},
                       pathway_sizes={}, genome_size=4008)
    assert ev.verdict.startswith("No known link")


def test_all_pairs_covers_every_combination():
    assert len(all_pairs(["a", "b", "c"], {}, {}, {}, 4008)) == 3


# --- quote validation -------------------------------------------------------

ABSTRACT = {
    "record_id": "10609885", "type": "article", "name": "AhpC and KatG",
    "detail": {"quotable_text": "All the isoniazid-resistant, AhpC-overexpressing strains\n"
                                "were also deficient in activity of the mycobacterial\n"
                                "catalase-peroxidase KatG."},
}


def test_verbatim_quote_is_verified():
    from kegg_string_mcp.agent.validate import validate

    text = 'PMID:10609885 "were also deficient in activity of the mycobacterial catalase-peroxidase KatG"'
    report = validate(text, {"10609885"}, records={"10609885": ABSTRACT})
    assert [q.status for q in report.quotes] == ["verified"]
    assert report.passed


def test_fabricated_quote_on_a_real_pmid_is_caught():
    """The failure set membership cannot see: a genuinely retrieved PMID carrying
    words that are not in it."""
    from kegg_string_mcp.agent.validate import validate

    text = 'PMID:10609885 "KatG binds directly to AhpC in a stable complex"'
    report = validate(text, {"10609885"}, records={"10609885": ABSTRACT})
    assert [q.status for q in report.quotes] == ["not_in_source"]
    assert not report.passed


def test_quote_survives_line_breaks_and_recapitalisation():
    """The stored text is flattened from XML; a model legitimately recapitalises a
    span at the start of a sentence. Neither is a fabrication."""
    from kegg_string_mcp.agent.validate import validate

    text = 'PMID:10609885 "Were also deficient in activity of the mycobacterial catalase-peroxidase"'
    report = validate(text, {"10609885"}, records={"10609885": ABSTRACT})
    assert [q.status for q in report.quotes] == ["verified"]


def test_quote_is_read_in_either_written_order():
    from kegg_string_mcp.agent.validate import extract_quotes

    after = extract_quotes('PMID:10609885 "were also deficient in activity"')
    before = extract_quotes('"were also deficient in activity" (PMID:10609885)')
    assert after == before == [("10609885", "were also deficient in activity")]


def test_quote_against_a_record_with_no_text_is_not_silently_passed():
    from kegg_string_mcp.agent.validate import validate

    empty = {"record_id": "999", "detail": {"quotable_text": ""}}
    report = validate('PMID:999 "some claimed finding here"', {"999"}, records={"999": empty})
    assert [q.status for q in report.quotes] == ["no_source_text"]
    assert not report.passed


def test_bare_numbers_are_not_treated_as_pmid_citations():
    """A bare 8-digit number could be a coordinate or a score. Only PMID: form
    counts, or ordinary prose would be flagged as unsupported."""
    from kegg_string_mcp.agent.validate import extract_citations

    assert extract_citations("the region spans 1673400 to 12345678 bases") == []
    assert extract_citations("see PMID:12345678") == ["12345678"]


def test_short_quoted_fragments_are_ignored():
    """A three-character 'quote' would pass containment against almost anything."""
    from kegg_string_mcp.agent.validate import extract_quotes

    assert extract_quotes('PMID:10609885 "KatG"') == []


# --- corpus manifest --------------------------------------------------------

PAPER = {
    "resolved": {"term": '"katG"'},
    "records": [
        {"record_id": "34086544", "type": "article", "name": "Review",
         "url": "https://pubmed.ncbi.nlm.nih.gov/34086544/",
         "detail": {"doi": "10.1080/x", "pmcid": "PMC8812758", "in_pmc": True,
                    "title": "Review", "journal": "J", "year": "2021",
                    "mentions": ["katG"], "quotable_text": "katG stuff"}},
        {"record_id": "23899494", "type": "article", "name": "Thioredoxin review",
         "url": "https://pubmed.ncbi.nlm.nih.gov/23899494/",
         "detail": {"doi": "10.1016/y", "pmcid": "", "in_pmc": False,
                    "title": "Thioredoxin", "journal": "J", "year": "2013",
                    "mentions": [], "quotable_text": "thioredoxin stuff"}},
    ],
    "record_ids": ["34086544", "23899494"],
}


def test_manifest_records_pmcid_not_just_doi(tmp_path):
    """A DOI resolves to a paywalled publisher; a PMCID is the licit full-text route."""
    store = RunStore(path=tmp_path / "r.jsonl", run_id="t")
    store.tool_result("pubmed_abstracts", {"gene": "katG"}, PAPER)
    by_pmid = {p["pmid"]: p for p in store.corpus_manifest()}
    assert by_pmid["34086544"]["pmcid"] == "PMC8812758" and by_pmid["34086544"]["in_pmc"]
    assert by_pmid["23899494"]["pmcid"] == "" and not by_pmid["23899494"]["in_pmc"]


def test_manifest_mentions_reflect_the_text_not_the_query(tmp_path):
    """PubMed matches on metadata the model never sees, so a paper can come back
    for a gene it never names. Filtering on the query term would admit those."""
    store = RunStore(path=tmp_path / "r.jsonl", run_id="t")
    store.tool_result("pubmed_abstracts", {"gene": "katG"}, PAPER)
    by_pmid = {p["pmid"]: p for p in store.corpus_manifest()}
    assert by_pmid["34086544"]["mentions"] == ["katG"]
    assert by_pmid["23899494"]["mentions"] == []


def test_manifest_excludes_non_article_records(tmp_path):
    store = RunStore(path=tmp_path / "r.jsonl", run_id="t")
    store.tool_result("kegg_pathways", {"gene": "katG"}, KEGG_RESULT)
    assert store.corpus_manifest() == []


def test_manifest_dedupes_papers_across_calls(tmp_path):
    store = RunStore(path=tmp_path / "r.jsonl", run_id="t")
    store.tool_result("pubmed_abstracts", {"gene": "katG"}, PAPER)
    store.tool_result("pubmed_abstracts", {"gene": "ahpC"}, PAPER)
    assert len(store.corpus_manifest()) == 2


def test_trailing_full_stop_does_not_fail_an_honest_quote():
    """Observed live: a model quoted a span verbatim and closed it with a full stop
    where the source sentence continues. Strict containment called that fabrication."""
    from kegg_string_mcp.agent.validate import validate

    source = {"detail": {"quotable_text":
              "There was also a significant inverse association between katG315 mutations and "
              "mutations in ahpC or inhA and between mutations in kasA and mutations in ahpC."}}
    text = ('PMID:16870753 "a significant inverse association between katG315 mutations and '
            'mutations in ahpC or inhA."')
    report = validate(text, {"16870753"}, records={"16870753": source})
    assert [q.status for q in report.quotes] == ["verified"]


def test_elided_quote_is_checked_fragment_by_fragment_in_order():
    from kegg_string_mcp.agent.validate import quote_in_source

    source = "The first finding was clear. Much intervening text. The second finding was not."
    assert quote_in_source("The first finding was clear ... The second finding was not", source)
    # Out of order must still fail: an elision cannot reorder the source.
    assert not quote_in_source("The second finding was not ... The first finding was clear", source)


def test_fabricated_quote_still_fails_after_the_punctuation_fix():
    """The tolerance must not swallow real fabrication."""
    from kegg_string_mcp.agent.validate import quote_in_source

    source = "KatG is a catalase-peroxidase required for isoniazid activation."
    assert not quote_in_source("KatG binds directly to AhpC in a stable complex.", source)
