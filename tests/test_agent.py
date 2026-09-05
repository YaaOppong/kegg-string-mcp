"""Store, evidence assembly and citation validation. No model, no network."""

from pathlib import Path

from kegg_string_mcp.agent.evidence import all_pairs, classify_pathway, pair_evidence
from kegg_string_mcp.agent.store import RunStore
from kegg_string_mcp.agent.validate import (
    check_quotes,
    extract_citations,
    extract_quotes,
    validate,
)

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


# --- failure triage ---------------------------------------------------------

SOURCE = {"detail": {"quotable_text":
          "There was also a significant inverse association between katG315 mutations and "
          "mutations in ahpC or inhA and between mutations in kasA and mutations in ahpC."}}


def test_quoting_artefact_and_fabrication_are_separated_by_similarity():
    """Ranks failures for human attention. Deliberately does NOT change the verdict:
    a model must never be able to argue its way past a deterministic check."""
    from kegg_string_mcp.agent.validate import check_quotes

    artefact = check_quotes(
        'PMID:1 "a significant inverse assoc between katG315 mutations and mutations in ahpC or inhA"',
        {"1": SOURCE})[0]
    fabricated = check_quotes(
        'PMID:1 "KatG binds directly to AhpC forming a stable heterodimeric complex"',
        {"1": SOURCE})[0]

    assert artefact.triage == "likely_quoting_artefact" and artefact.similarity >= 0.9
    assert fabricated.triage == "likely_fabricated" and fabricated.similarity < 0.6
    # Both remain failures. Triage informs; it does not absolve.
    assert artefact.status == fabricated.status == "not_in_source"


def test_triage_reports_the_closest_span_so_a_human_can_adjudicate():
    from kegg_string_mcp.agent.validate import check_quotes

    check = check_quotes('PMID:1 "a significant inverse assoc between katG315 mutations"',
                         {"1": SOURCE})[0]
    assert "inverse association between katg315" in check.closest_span


def test_a_near_miss_still_fails_validation_overall():
    from kegg_string_mcp.agent.validate import validate

    report = validate('PMID:1 "a significant inverse assoc between katG315 mutations"',
                      {"1"}, records={"1": SOURCE})
    assert not report.passed, "a high-similarity near miss must not silently pass"


def test_nearest_span_is_empty_when_there_is_no_source():
    from kegg_string_mcp.agent.validate import nearest_span

    assert nearest_span("anything", "") == (0.0, "")


# --- regressions from probing the agent layer -------------------------------

def _partners(n, prefix):
    return [{"record_id": f"83332.{prefix}{i}", "name": f"{prefix}{i}", "detail": {}} for i in range(n)]


def test_truncated_partner_lists_are_flagged_not_reported_as_degree():
    """`limit` caps retrieval, so counting retrieved partners made a hub with 500
    indistinguishable from a gene with exactly 20 -- disabling the hub check the
    field exists for."""
    from kegg_string_mcp.agent.evidence import pair_evidence

    ev = pair_evidence("hub", "small", pathways={}, pathway_sizes={}, genome_size=4008,
                       partners={"hub": _partners(20, "X"), "small": _partners(3, "Y")},
                       partner_limit=20)
    assert ev.truncated == ["hub"]
    assert ev.partners_retrieved == {"hub": 20, "small": 3}


def test_shared_partners_verdict_is_caveated_when_degree_is_unknown():
    from kegg_string_mcp.agent.evidence import pair_evidence

    shared = _partners(20, "S")
    ev = pair_evidence("a", "b", pathways={}, pathway_sizes={}, genome_size=4008,
                       partners={"a": shared, "b": list(shared)}, partner_limit=20)
    assert "true network degree is unknown" in ev.verdict
    assert "cannot be distinguished" in ev.verdict


def test_untruncated_lists_get_no_spurious_caveat():
    from kegg_string_mcp.agent.evidence import pair_evidence

    shared = _partners(3, "S")
    ev = pair_evidence("a", "b", pathways={}, pathway_sizes={}, genome_size=4008,
                       partners={"a": shared, "b": list(shared)}, partner_limit=20)
    assert ev.truncated == [] and "degree is unknown" not in ev.verdict


def test_empty_summary_does_not_pass_validation():
    """A run that produced nothing had zero citations, so every check trivially
    passed and a silent failure was reported as a success."""
    from kegg_string_mcp.agent.validate import validate

    report = validate("", citable_ids={"mtu00360"}, records={})
    assert not report.passed
    assert "nothing could be validated" in report.summary_line()


def test_whitespace_only_summary_does_not_pass_either():
    from kegg_string_mcp.agent.validate import validate

    assert not validate("   \n  ", citable_ids=set(), records={}).passed


def test_a_real_summary_still_passes():
    from kegg_string_mcp.agent.validate import validate

    assert validate("katG is in mtu00360.", citable_ids={"mtu00360"}, records={}).passed


# --- regressions from the fifth review pass ---------------------------------

def test_consecutive_quotes_are_not_cross_attributed():
    """Two quotes in a row, each closing with its own citation, is the natural way
    to write this. The leading-citation pattern reached across the sentence
    boundary and bound the second quote to the first PMID as well."""
    from kegg_string_mcp.agent.validate import extract_quotes

    pairs = extract_quotes('a "encodes catalase-peroxidase" (PMID: 111). '
                           'b "ahpC compensates for katG loss" (PMID: 222).')
    assert pairs == [("111", "encodes catalase-peroxidase"),
                     ("222", "ahpC compensates for katG loss")]


def test_decimals_and_dois_are_not_read_as_string_ids():
    """Epistasis mode is handed pathway sizes and asked to weigh base rates, so
    percentages and fold-changes are routine prose. Each was being reported as an
    unsupported citation and failing the run."""
    from kegg_string_mcp.agent.validate import extract_citations

    found = extract_citations("mtu01100 holds 698 of ~4000 genes (17.5%). "
                              "doi 10.1038/nature12345, a 12.3-fold change, see 83332.Rv1908c")
    assert found == ["mtu01100", "83332.Rv1908c"]


def test_alias_spellings_of_the_same_gene_do_not_fail_cross_target(tmp_path):
    """The tool schemas invite a different spelling, so a model annotating katG may
    call string_partners(gene='Rv1908c'). Keying only on the literal argument made
    every ID from that call report as cross_target."""
    from kegg_string_mcp.agent.validate import validate

    store = RunStore(path=tmp_path / "r.jsonl", run_id="t")
    store.tool_result("string_partners", {"gene": "Rv1908c"}, STRING_RESULT)
    report = validate("katG partners with 83332.Rv1909c.", store.citable_ids,
                      store.per_target, "KATG")
    assert report.passed, report.summary_line()


def test_unknown_parameter_returns_an_envelope_not_a_typeerror():
    """Model-supplied arguments are untrusted input; a schema deviation used to
    kill the run mid-flight instead of giving the model something to correct."""
    from kegg_string_mcp.agent.pipeline import Tools

    result = Tools.__call__(object.__new__(Tools), "string_partners",
                            {"gene": "katG", "organism": "mtu"})
    assert result["records"] == []
    assert "not a parameter of string_partners" in " ".join(result["notes"])


def test_string_typed_limit_is_coerced_not_crashed(monkeypatch):
    from kegg_string_mcp.agent.pipeline import _coerce

    clean, problems = _coerce("string_partners", {"gene": "katG", "limit": "20"})
    assert clean == {"gene": "katG", "limit": 20} and not problems

    _, bad = _coerce("string_partners", {"gene": "katG", "limit": "twenty"})
    assert bad and "must be int" in bad[0]


def test_unknown_genome_size_does_not_downgrade_a_container_pathway():
    """Falling through a zero denominator labelled mtu01100 (698 genes) 'moderate',
    reopening the base-rate trap the function exists to close."""
    from kegg_string_mcp.agent.evidence import classify_pathway

    label, note = classify_pathway(698, 0)
    assert label == "unknown" and "could not be determined" in note


def test_mentions_merge_across_calls_for_the_same_paper(tmp_path):
    """The same PMID returned for two genes kept only the last record, losing the
    first gene's mention -- the field a downstream corpus filters on."""
    store = RunStore(path=tmp_path / "r.jsonl", run_id="t")
    paper = lambda gene: {
        "resolved": {}, "record_ids": ["999"],
        "records": [{"record_id": "999", "type": "article", "name": "p", "url": "u",
                     "detail": {"mentions": [gene], "quotable_text": "t", "doi": "",
                                "pmcid": "", "in_pmc": False, "title": "p",
                                "journal": "", "year": ""}}],
    }
    store.tool_result("pubmed_abstracts", {"gene": "katG"}, paper("katG"))
    store.tool_result("pubmed_abstracts", {"gene": "ahpC"}, paper("ahpC"))
    assert store.corpus_manifest()[0]["mentions"] == ["ahpC", "katG"]


def test_invalid_organism_is_rejected_by_the_direct_call_sites(kegg=None):
    """pathways() guarded this; gene_index and pathway_sizes are called directly by
    the agent pipeline and bypassed it."""
    import pytest

    from kegg_string_mcp.kegg import KeggClient

    client = KeggClient(http=None)
    for method in (client.gene_index, client.pathway_sizes):
        with pytest.raises(ValueError, match="not a valid KEGG organism code"):
            method("../../info/mtu")


def test_free_text_literature_queries_do_not_fail_cross_target(tmp_path):
    """pubmed_abstracts takes free text by design, so keying the target on the raw
    argument flagged every PMID from a multi-word query as cross-target."""
    from kegg_string_mcp.agent.validate import validate

    store = RunStore(path=tmp_path / "r.jsonl", run_id="t")
    store.tool_result("pubmed_abstracts",
                      {"gene": "furA katG regulation oxidative stress"},
                      {"resolved": {}, "record_ids": ["27328747"],
                       "records": [{"record_id": "27328747", "type": "article", "name": "p",
                                    "url": "u", "detail": {"mentions": ["furA"],
                                                           "quotable_text": "furA text"}}]})
    report = validate("See PMID:27328747.", store.citable_ids, store.per_target, "FURA")
    assert report.passed, report.summary_line()


def test_literature_is_attributed_by_what_the_paper_mentions(tmp_path):
    store = RunStore(path=tmp_path / "r.jsonl", run_id="t")
    store.tool_result("pubmed_abstracts", {"gene": "oxidative stress"},
                      {"resolved": {}, "record_ids": ["1"],
                       "records": [{"record_id": "1", "type": "article", "name": "p", "url": "u",
                                    "detail": {"mentions": ["katG", "ahpC"], "quotable_text": "t"}}]})
    assert "1" in store.per_target["KATG"] and "1" in store.per_target["AHPC"]


def test_uniprot_is_dispatchable_and_argument_checked():
    """The fourth tool goes through the same untrusted-input guard as the others."""
    from kegg_string_mcp.agent.pipeline import TOOL_PARAMS, _coerce

    assert "uniprot_protein" in TOOL_PARAMS
    clean, problems = _coerce("uniprot_protein", {"gene": "gyrA", "limit": "2"})
    assert clean == {"gene": "gyrA", "limit": 2} and not problems

    _, bad = _coerce("uniprot_protein", {"gene": "gyrA", "organism": "mtu"})
    assert bad and "not a parameter of uniprot_protein" in bad[0]


def test_tool_schemas_come_from_the_server_not_a_copy():
    """The agent used to hold its own copy of the schemas, and they drifted -- the
    model driving the pipeline saw shorter descriptions than an external MCP client
    did. There is now one definition, so the copies cannot diverge."""
    import asyncio

    from kegg_string_mcp.agent import loop
    from kegg_string_mcp.agent.pipeline import server_tool_schemas
    from kegg_string_mcp.server import mcp

    assert not hasattr(loop, "TOOL_SCHEMAS"), "a second schema definition reappeared"

    schemas = asyncio.run(server_tool_schemas())
    served = {t.name: (t.description or "") for t in asyncio.run(mcp.list_tools())}
    assert {s["name"] for s in schemas} == set(served)
    for schema in schemas:
        assert schema["description"] == served[schema["name"]]


def test_dispatch_may_be_sync_or_async():
    """The direct dispatch is sync and the MCP one is async; the loop must take both."""
    import asyncio

    from kegg_string_mcp.agent.pipeline import _call

    def sync_tool(name, arguments):
        return {"records": [], "record_ids": [], "notes": ["sync"]}

    async def async_tool(name, arguments):
        return {"records": [], "record_ids": [], "notes": ["async"]}

    assert asyncio.run(_call(sync_tool, "x", {}))["notes"] == ["sync"]
    assert asyncio.run(_call(async_tool, "x", {}))["notes"] == ["async"]


# --- regressions from the develop review ------------------------------------

def test_uniprot_accessions_are_recognised_as_citations():
    """UniProt was added as a fourth tool without teaching the validator about it,
    so a fabricated accession was not flagged -- it was not even seen."""
    from kegg_string_mcp.agent.validate import extract_citations

    found = extract_citations("P9WIE5 and A0A123XYZ9 alongside mtu00360, 83332.Rv1909c, PMID:123456")
    assert set(found) == {"P9WIE5", "A0A123XYZ9", "mtu00360", "83332.Rv1909c", "123456"}


def test_ordinary_prose_is_not_read_as_an_accession():
    from kegg_string_mcp.agent.validate import extract_citations

    assert extract_citations("KatG and AhpC bind DNA and RNA in vitro.") == []


def test_fabricated_quote_on_a_uniprot_record_is_caught():
    """The docstring claimed the span validator applied unchanged. It did not apply
    at all: a wholly invented quote on a real accession returned passed=True."""
    from kegg_string_mcp.agent.validate import validate

    report = validate('P9WIE5 "a wholly invented finding about this protein"', {"P9WIE5"},
                      records={"P9WIE5": {"detail": {"quotable_text": "nothing like that"}}})
    assert [q.status for q in report.quotes] == ["not_in_source"]
    assert not report.passed


def test_verbatim_quote_on_a_uniprot_record_passes():
    from kegg_string_mcp.agent.validate import validate

    source = {"detail": {"quotable_text": "Bifunctional enzyme with both catalase and "
                                          "broad-spectrum peroxidase activity."}}
    report = validate('P9WIE5 "Bifunctional enzyme with both catalase"', {"P9WIE5"},
                      records={"P9WIE5": source})
    assert [q.status for q in report.quotes] == ["verified"] and report.passed


def test_a_citation_inside_a_quotation_is_not_bound_to_the_next_one():
    """`"essential (PMID: 1)" and "it also binds NADH"` — the PMID belongs to the
    first quotation. An earlier review argued the gap class was too narrow and I
    widened it; that let this pattern reach across a finished quotation and
    attribute the second span to the first span's citation."""
    from kegg_string_mcp.agent.validate import extract_quotes

    pairs = extract_quotes('The protein is “essential (PMID: 12345678)” and separately '
                           '“it also binds NADH in vitro” per UniProt.')
    assert ("12345678", "it also binds NADH in vitro") not in pairs


def test_pmids_from_a_uniprot_statement_are_citable(tmp_path):
    """UniProt advertises the PMIDs evidencing a function statement, so citing one
    is correct and traceable -- it was being scored as a hallucination."""
    from kegg_string_mcp.agent.validate import validate

    store = RunStore(path=tmp_path / "r.jsonl", run_id="t")
    store.tool_result("uniprot_protein", {"gene": "katG"}, {
        "resolved": {"accession": "P9WIE5"}, "record_ids": ["P9WIE5"],
        "records": [{"record_id": "P9WIE5", "type": "protein", "name": "Catalase", "url": "u",
                     "detail": {"quotable_text": "t", "function_statements": [
                         {"text": "t", "supporting_pmids": ["18178143"], "experimental": True,
                          "tiers": ["experimental"], "evidence_codes": ["ECO:0000269"]}]}}],
    })
    assert "18178143" in store.citable_ids
    assert validate("See PMID:18178143.", store.citable_ids, store.per_target, "KATG").passed


def _fake_tool(name: str = "kegg_pathways"):
    """Stands in for an mcp.types.Tool as advertised by list_tools()."""
    from types import SimpleNamespace

    return SimpleNamespace(
        name=name, description="d",
        input_schema={"type": "object",
                      "properties": {"gene": {"type": "string"},
                                     "organism": {"type": "string"}},
                      "required": ["gene"]},
    )


def test_mcp_dispatch_validates_against_the_servers_advertised_schema():
    """MCP silently drops undeclared arguments, so a model asking for a parameter a
    tool does not have got results as though the constraint applied. The direct
    dispatch refuses that; both paths must behave the same."""
    import asyncio

    from kegg_string_mcp.agent.mcp_tools import McpTools

    class Session:
        called = False

        async def call_tool(self, name, arguments):
            Session.called = True

    tools = McpTools(Session(), [_fake_tool()])
    envelope = asyncio.run(tools("kegg_pathways", {"gene": "katG", "species": 83332}))
    assert envelope["records"] == []
    assert "'species' is not a parameter of kegg_pathways" in envelope["notes"][0]
    assert not Session.called, "an invalid call must not reach the server"

    missing = asyncio.run(tools("kegg_pathways", {"organism": "mtu"}))
    assert "'gene' is required" in missing["notes"][0]


def test_both_dispatch_paths_coerce_types_identically():
    """MCP is the default path, so refusing a type the direct dispatch coerces would
    turn a survivable model deviation into an empty envelope."""
    import asyncio

    from kegg_string_mcp.agent.mcp_tools import McpTools
    from kegg_string_mcp.agent.pipeline import _coerce

    sent = {}

    class Session:
        async def call_tool(self, name, arguments):
            sent.update(arguments)

    tools = McpTools(Session(), [_fake_tool()])
    asyncio.run(tools("kegg_pathways", {"gene": 42}))
    assert sent == {"gene": "42"}, "MCP should coerce, as the direct path does"

    direct, problems = _coerce("kegg_pathways", {"gene": 42})
    assert direct == {"gene": "42"} and not problems


def test_mcp_dispatch_surfaces_a_tool_error_as_an_envelope():
    """A tool error is data the model can act on, not a crash -- same contract as
    the direct dispatch, which returns an envelope for bad arguments."""
    import asyncio

    from kegg_string_mcp.agent.mcp_tools import McpTools

    class Block:
        type, text = "text", "boom"

    class Result:
        is_error, content, structured_content = True, [Block()], None

    class Session:
        async def call_tool(self, name, arguments):
            return Result()

    envelope = asyncio.run(McpTools(Session(), [_fake_tool()])("kegg_pathways", {"gene": "katG"}))
    assert envelope["records"] == [] and "returned an error" in envelope["notes"][0]


def test_mcp_dispatch_falls_back_to_text_when_structured_output_is_absent():
    import asyncio
    import json

    from kegg_string_mcp.agent.mcp_tools import McpTools

    class Block:
        type = "text"
        text = json.dumps({"records": [], "record_ids": ["mtu00360"], "notes": []})

    class Result:
        is_error, content, structured_content = False, [Block()], None

    class Session:
        async def call_tool(self, name, arguments):
            return Result()

    envelope = asyncio.run(McpTools(Session(), [_fake_tool()])("kegg_pathways", {"gene": "katG"}))
    assert envelope["record_ids"] == ["mtu00360"]


def test_child_env_forwards_credentials_but_not_the_whole_environment(monkeypatch):
    from kegg_string_mcp.agent.mcp_tools import child_env

    monkeypatch.setenv("NCBI_EMAIL", "you@example.org")
    monkeypatch.setenv("SOME_UNRELATED_SECRET", "nope")
    env = child_env()
    assert env["NCBI_EMAIL"] == "you@example.org"
    assert "SOME_UNRELATED_SECRET" not in env


def test_quotes_are_not_attributed_to_structured_records():
    """A KEGG pathway ID has no text to quote from. Attributing a nearby quote to
    one made a model quoting a tool's own note fail as 'no source text' -- a false
    positive on correct output."""
    from kegg_string_mcp.agent.validate import extract_quotes

    text = ('katG is in mtu00360, mtu01110. ahpC returned an empty record set with the note '
            'that "the gene exists in KEGG but is not mapped to any pathway".')
    assert extract_quotes(text) == []


def test_quotes_are_still_attributed_to_prose_records():
    from kegg_string_mcp.agent.validate import extract_quotes

    assert extract_quotes('PMID:10609885 "were also deficient in activity"')
    assert extract_quotes('P9WIE5 "Bifunctional enzyme with both catalase"')


def test_markdown_emphasis_around_a_quote_is_not_fabrication():
    """Models routinely italicise a quotation, and the markers land inside the
    captured span. Observed live on the gyrA run: five quotes failed purely because
    the model wrote *"..."* rather than "...".
    """
    from kegg_string_mcp.agent.validate import quote_in_source

    source = "The gyrA mutations occurring most frequently in fluoroquinolone-resistant isolates."
    assert quote_in_source("*gyrA mutations occurring most frequently*", source)
    assert quote_in_source("**gyrA mutations occurring most frequently**", source)
    assert quote_in_source("`gyrA mutations occurring most frequently`", source)


def test_emphasis_tolerance_does_not_swallow_fabrication():
    from kegg_string_mcp.agent.validate import quote_in_source

    source = "The gyrA mutations occurring most frequently in resistant isolates."
    assert not quote_in_source("*gyrA binds directly to the ribosome*", source)


def test_a_pmid_cited_by_uniprot_has_something_to_check_a_quote_against(tmp_path):
    """Making the PMID citable without registering a record left the quote check
    with nothing to compare against, so a quote attached to it failed as 'no source
    text' — on exactly the behaviour the prompt asks for."""
    from kegg_string_mcp.agent.validate import validate

    store = RunStore(path=tmp_path / "r.jsonl", run_id="t")
    store.tool_result("uniprot_protein", {"gene": "katG"}, {
        "resolved": {"accession": "P9WIE5"}, "record_ids": ["P9WIE5"],
        "records": [{"record_id": "P9WIE5", "type": "protein", "name": "Catalase", "url": "u",
                     "detail": {"quotable_text": "a bifunctional catalase-peroxidase enzyme",
                                "function_statements": [
                                    {"text": "a bifunctional catalase-peroxidase enzyme",
                                     "supporting_pmids": ["18178143"], "experimental": True,
                                     "tiers": ["experimental"], "evidence_codes": ["ECO:0000269"]}]}}],
    })
    report = validate('UniProt states "a bifunctional catalase-peroxidase enzyme" (PMID:18178143).',
                      store.citable_ids, store.per_target, "KATG", records=store.records)
    assert report.passed, report.summary_line()
    assert [q.status for q in report.quotes] == ["verified"]


def test_a_fabricated_quote_on_such_a_pmid_still_fails(tmp_path):
    from kegg_string_mcp.agent.validate import validate

    store = RunStore(path=tmp_path / "r.jsonl", run_id="t")
    store.tool_result("uniprot_protein", {"gene": "katG"}, {
        "resolved": {}, "record_ids": ["P9WIE5"],
        "records": [{"record_id": "P9WIE5", "type": "protein", "name": "n", "url": "u",
                     "detail": {"quotable_text": "a catalase-peroxidase",
                                "function_statements": [
                                    {"text": "a catalase-peroxidase",
                                     "supporting_pmids": ["18178143"]}]}}],
    })
    report = validate('PMID:18178143 "binds directly to the ribosome in vitro"',
                      store.citable_ids, store.per_target, "KATG", records=store.records)
    assert [q.status for q in report.quotes] == ["not_in_source"]


# --- lineage-marker record IDs -------------------------------------------------
# tbdb IDs became citable when lineage_markers joined the stage 1 tools. These
# pin the properties that would fail silently: that they are recognised at all,
# that a bare genome coordinate is not mistaken for one, and that they are
# structured rather than quotable.

def test_a_lineage_marker_id_is_a_citation():
    assert extract_citations("phoR carries a lineage marker (tbdb:852641).") == ["tbdb:852641"]


def test_several_lineage_markers_are_all_extracted():
    text = "Two markers, tbdb:852641 and tbdb:853469, fall in phoR."
    assert extract_citations(text) == ["tbdb:852641", "tbdb:853469"]


def test_a_bare_genome_coordinate_is_not_a_citation():
    """An unprefixed 7-digit H37Rv position is indistinguishable from a PMID.
    The prefix is what stops the validator confusing a coordinate with a paper --
    without it, "position 852641" would cite a PubMed article."""
    assert extract_citations("Position 852641 falls inside phoR.") == []


def test_a_lineage_marker_is_not_quotable():
    """A marker is structured, like a KEGG pathway ID: the record means one thing
    and there is no text to quote from it."""
    assert extract_quotes('The gene is "described as a marker" (tbdb:852641).') == []


def test_lineage_markers_mix_with_the_other_sources():
    text = "See tbdb:852641, PMID:35919400 and mtu00360."
    assert set(extract_citations(text)) == {"tbdb:852641", "35919400", "mtu00360"}


def test_an_uncited_lineage_marker_is_caught_as_unsupported():
    """The point of the whole layer: a marker the tools never returned must fail."""
    report = validate("phoR is a lineage marker (tbdb:999999).",
                      citable_ids={"tbdb:852641"})
    assert not report.passed
    assert [c.identifier for c in report.unsupported] == ["tbdb:999999"]


def test_a_retrieved_lineage_marker_passes():
    report = validate("phoR carries a lineage marker (tbdb:852641).",
                      citable_ids={"tbdb:852641"})
    assert report.passed
    assert report.unsupported == []


def test_a_citation_inside_a_quotation_does_not_bind_the_next_span():
    """From a real run. Quoting a tool note that names its own record put an
    accession inside the quotation with no citation after the closing mark. The
    leading-citation pattern ran from that accession, across a closing mark it
    could not tell from an opening one, and bound the model's own prose as a
    quotation -- which then failed as "not in source". A false accusation of
    fabrication against correct output is the one failure a validator must not
    have."""
    text = ('The note is explicit: "UniProt holds no FUNCTION statement for P71814 - '
            'the entry exists but its function is not described." '
            '`has_experimental_function` is false, and there is no subunit annotation '
            'to quote. The name is "Possible two component system regulator" (P71814).')
    quotes = extract_quotes(text)
    assert quotes == [("P71814", "Possible two component system regulator")]


def test_a_leading_citation_outside_any_quotation_still_binds():
    """The guard must not cost the ordinary leading-citation form."""
    text = 'PMID:35919400 reported that "the isolates carried katG mutations".'
    assert extract_quotes(text) == [("35919400", "the isolates carried katG mutations")]


# --- quoting a tool note ----------------------------------------------------
# A note is part of what the tool said in this run, so quoting one verbatim is
# quoting a real source. Searching only record quotable_text reported it as
# likely fabricated.

NOTE = ("No function statement for A0A0N7EHL5, A0A0N9DRZ6 carries experimental "
        "evidence (ECO:0000269); all are inferred.")
_SUMMARY = f'Two further entries matched (A0A0N7EHL5, A0A0N9DRZ6) but the tool notes "{NOTE}"'


def test_a_verbatim_tool_note_is_not_a_fabrication():
    """From a real katG run. The quote is note 3 of the UniProt result, word for
    word, and was reported NOT_IN_SOURCE [likely_fabricated] because notes were
    never searched -- the record has its own text, so the quote was compared
    against the wrong source rather than found missing."""
    records = {"A0A0N7EHL5": {"record_id": "A0A0N7EHL5",
                              "detail": {"quotable_text": "Catalase-peroxidase."}}}
    assert check_quotes(_SUMMARY, records)[0].status == "not_in_source"
    checks = check_quotes(_SUMMARY, records, {"A0A0N7EHL5": [NOTE]})
    assert [c.status for c in checks] == ["verified"]


def test_a_note_belonging_to_another_record_does_not_verify_a_quote():
    """The loophole to avoid: notes are searched per record, so a quote cannot
    be verified against a note the tool returned about something else."""
    records = {"A0A0N7EHL5": {"record_id": "A0A0N7EHL5",
                              "detail": {"quotable_text": "Catalase-peroxidase."}}}
    checks = check_quotes(_SUMMARY, records, {"P9WIE5": [NOTE]})
    assert [c.status for c in checks] == ["not_in_source"]


def test_a_genuinely_fabricated_quote_still_fails_with_notes_present():
    """The fix must not turn the check into one that passes everything."""
    records = {"A0A0N7EHL5": {"record_id": "A0A0N7EHL5",
                              "detail": {"quotable_text": "Catalase-peroxidase activity."}}}
    text = 'The entry states "this protein activates isoniazid" (A0A0N7EHL5).'
    checks = check_quotes(text, records, {"A0A0N7EHL5": [NOTE]})
    assert [c.status for c in checks] == ["not_in_source"]


def test_store_keys_notes_to_the_records_of_the_same_result(tmp_path: Path):
    store = _store(tmp_path)
    store.tool_result("uniprot_protein", {"gene": "katG"},
              {"records": [{"record_id": "P9WIE5", "detail": {}}],
               "record_ids": ["P9WIE5"], "notes": ["a note about katG"]})
    store.tool_result("kegg_pathways", {"gene": "katG"},
              {"records": [{"record_id": "mtu00360"}], "record_ids": ["mtu00360"],
               "notes": ["a note about pathways"]})
    assert store.notes["P9WIE5"] == ["a note about katG"]
    assert store.notes["mtu00360"] == ["a note about pathways"]
