import asyncio

from kegg_string_mcp.server import mcp


def _tools():
    return {t.name: t for t in asyncio.run(mcp.list_tools())}


def test_both_tools_are_registered():
    assert set(_tools()) == {"kegg_pathways", "string_partners"}


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


def test_descriptions_warn_the_model_about_textmining():
    """The model reads these. The textmining caveat has to be in the description,
    not only in the code comments."""
    assert "textmining" in _tools()["string_partners"].description.lower()


def test_descriptions_distinguish_no_data_from_no_match():
    assert "resolve" in _tools()["kegg_pathways"].description.lower()
