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
