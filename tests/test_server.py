import asyncio

from kegg_string_mcp.server import mcp


def _tools():
    return {t.name: t for t in asyncio.run(mcp.list_tools())}


def test_all_tools_are_registered():
    assert set(_tools()) == {"kegg_pathways", "string_partners", "pubmed_abstracts",
                             "uniprot_protein", "lineage_markers"}


def test_tools_declare_read_only():
    """These tools only fetch. Declaring it lets a client skip approval prompts."""
    for tool in _tools().values():
        assert tool.annotations.read_only_hint is True


def test_tools_expose_structured_output_schemas():
    for tool in _tools().values():
        assert tool.output_schema, f"{tool.name} has no structured output schema"


def test_defaults_target_mtb_h37rv():
    props = _tools()["string_partners"].input_schema["properties"]
    assert props["species"]["default"] == 83332
    assert _tools()["kegg_pathways"].input_schema["properties"]["organism"]["default"] == "mtu"
    assert _tools()["pubmed_abstracts"].input_schema["properties"]["organism"]["default"] == (
        "Mycobacterium tuberculosis"
    )


def test_descriptions_warn_the_model_about_textmining():
    """The model reads these. The textmining caveat has to be in the description,
    not only in the code comments."""
    assert "textmining" in _tools()["string_partners"].description.lower()


def test_descriptions_distinguish_no_data_from_no_match():
    assert "resolve" in _tools()["kegg_pathways"].description.lower()


def test_pubmed_description_warns_that_a_hit_is_not_a_resolved_identifier():
    """The weakness of this tool relative to the other two is the thing the model
    most needs to know, so it belongs in the description it always reads."""
    description = _tools()["pubmed_abstracts"].description
    assert "not identifier resolution" in description.lower()


def test_pubmed_description_tells_the_model_to_quote_verbatim():
    """Span grounding only works if the model knows a paraphrase will fail it."""
    description = _tools()["pubmed_abstracts"].description
    assert "quotable_text" in description
    assert "verbatim" in description.lower()


def test_server_instructions_state_research_use_only():
    """The scope boundary should reach any model that connects, not just a human
    reading the README."""
    assert "RESEARCH USE ONLY" in (mcp.instructions or "")
    assert "not a clinical decision support system" in (mcp.instructions or "")


# What the prompts must say about each tool for the model to have a reason to
# reach for it. The prompts name SOURCES rather than tool names -- the model gets
# tool names and descriptions from the server -- so this maps each registered
# tool to the word that has to appear. A checklist, and deliberately one: adding
# a tool without telling the model anything about it is the failure being
# guarded against, and only a list knows what a new tool should have said.
TOOL_CUES = {
    "kegg_pathways": "KEGG",
    "string_partners": "STRING",
    "pubmed_abstracts": "abstract",
    "uniprot_protein": "UniProt",
    "lineage_markers": "lineage_markers",
}


def test_every_registered_tool_has_a_reason_to_be_called_in_the_prompts():
    """A tool the prompts give no reason to call is invisible in practice.

    lineage_markers was registered, citable and correct, and the agent never
    called it: the shared prompt still opened "using only the KEGG and STRING
    tools", written before UniProt, PubMed and lineage existed. Nothing failed --
    the annotation quietly lacked a source, which is the worst way for this to go
    wrong. This fails the moment a tool is registered without the prompts being
    updated to match.
    """
    from kegg_string_mcp.agent.modes import EPISTASIS, SINGLE_GENE

    prompts = SINGLE_GENE + EPISTASIS
    unlisted = sorted(set(_tools()) - set(TOOL_CUES))
    assert not unlisted, f"registered tools with no prompt cue declared: {unlisted}"
    missing = sorted(name for name, cue in TOOL_CUES.items()
                     if name in _tools() and cue not in prompts)
    assert not missing, f"registered but the prompts never motivate calling: {missing}"


def test_the_prompts_do_not_claim_a_stale_tool_set():
    """The opening line enumerated two sources and went stale three times."""
    from kegg_string_mcp.agent.modes import SINGLE_GENE

    assert "only the KEGG and STRING tools" not in SINGLE_GENE


def test_the_lineage_check_is_unconditional():
    """Asked for every locus, not only when the question mentions a scan --
    a negative is as informative as a positive."""
    from kegg_string_mcp.agent.modes import EPISTASIS, SINGLE_GENE

    for prompt in (SINGLE_GENE, EPISTASIS):
        assert "EVERY gene or locus" in prompt
