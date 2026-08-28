"""MCP server exposing KEGG, STRING and PubMed as model-callable tools.

Every tool is a deterministic lookup: same input, same cache, same output. The
model chooses *which* to call and *when it has enough*; it never decides what a
record says.

`pubmed_abstracts` is the one that strains that last clause, and it is worth
being exact about how. It does not weaken it at the retrieval layer -- the tool
returns the text verbatim and decides nothing -- but an abstract is prose, so a
*claim* drawn from one is the model interpreting content, which a KEGG pathway ID
never required. The answer is not to trust it: every article record carries
`quotable_text`, the exact retrieved string, so a claim's quoted span can be
checked by containment the way a record_id is checked by set membership. Both
tiers are the pipeline's job, not the model's.

Tool descriptions below are written for the model, not for a human reader. They
state the things it would otherwise get wrong: that an empty result can mean
"identifier did not resolve" rather than "no data exists", that STRING's headline
score already includes literature co-mention, and that a PubMed hit is a text
match rather than a resolved identifier.
"""

from __future__ import annotations

from mcp.server import MCPServer
from mcp.types import ToolAnnotations

from kegg_string_mcp.cache import DiskCache
from kegg_string_mcp.http import PoliteClient
from kegg_string_mcp.kegg import KeggClient
from kegg_string_mcp.provenance import ToolResult
from kegg_string_mcp.pubmed import DEFAULT_LIMIT, MTB_H37RV_NAME, PubMedClient
from kegg_string_mcp.string_db import MTB_H37RV, StringClient

mcp = MCPServer(
    name="kegg-string",
    version="0.1.0",
    instructions=(
        "Deterministic lookups against KEGG, STRING and PubMed for gene annotation. Every record "
        "returned carries a stable record_id and a resolvable URL. Cite only record_ids that "
        "appear in a tool result's record_ids list -- citations are checked programmatically "
        "against what was actually retrieved. Read the `notes` field: an empty `records` list "
        "may mean the identifier failed to resolve rather than that no data exists. When a claim "
        "rests on a PubMed article, quote a verbatim span of that record's `quotable_text` and "
        "attach it to the PMID -- quotes are checked against the retrieved text, so a paraphrase "
        "or a remembered fact will not pass even under a correct PMID."
    ),
)

_http = PoliteClient(DiskCache())
_kegg = KeggClient(_http)
_string = StringClient(_http)
_pubmed = PubMedClient(_http)

READ_ONLY = ToolAnnotations(readOnlyHint=True, openWorldHint=True)


@mcp.tool(
    annotations=READ_ONLY,
    description=(
        "Look up the KEGG pathways a gene belongs to. Accepts a KEGG gene ID (mtu:Rv1908c), a "
        "locus tag (Rv1908c), or a gene symbol (katG); symbols are resolved by exact match "
        "against the organism's gene list, never fuzzily. Returns one record per pathway with "
        "its KEGG pathway ID as record_id. If `records` is empty, check `notes`: it "
        "distinguishes 'identifier did not resolve' from 'gene exists but maps to no pathway'. "
        "KEGG is free for academic use; commercial use requires a licence from Pathway Solutions."
    ),
)
def kegg_pathways(gene: str, organism: str = "mtu") -> ToolResult:
    """KEGG pathway membership for one gene.

    Args:
        gene: KEGG gene ID, locus tag, or gene symbol.
        organism: KEGG organism code. Defaults to 'mtu' (M. tuberculosis H37Rv).
    """
    return _kegg.pathways(gene=gene, organism=organism)


@mcp.tool(
    annotations=READ_ONLY,
    description=(
        "Look up STRING interaction partners for a gene. Returns one record per partner with "
        "the STRING protein ID as record_id, the combined score, and the full per-channel "
        "score breakdown. IMPORTANT: STRING's combined_score includes a textmining channel "
        "(literature co-mention), so a high score is not independent of literature evidence "
        "about the same pair. Each record carries evidence_beyond_textmining -- use it before "
        "describing an interaction as experimentally or database supported. An empty result at "
        "a high required_score is not evidence that a protein has no partners."
    ),
)
def string_partners(
    gene: str, species: int = MTB_H37RV, limit: int = 20, required_score: int = 700
) -> ToolResult:
    """STRING interaction partners for one gene.

    Args:
        gene: Gene symbol, locus tag, or STRING protein ID.
        species: NCBI taxon ID. Defaults to 83332 (M. tuberculosis H37Rv).
        limit: Maximum partners to return.
        required_score: STRING confidence threshold, 0-1000. 700 = high confidence.
    """
    return _string.partners(gene=gene, species=species, limit=limit, required_score=required_score)


@mcp.tool(
    annotations=READ_ONLY,
    description=(
        "Search PubMed for articles about a gene and return their titles and abstracts. Returns "
        "one record per article with the PMID as record_id. IMPORTANT: unlike the KEGG and STRING "
        "tools, this is a relevance-ranked TEXT SEARCH, not identifier resolution -- an article is "
        "returned because it matched the query string, which is not the same as being about this "
        "gene, and there is no 'did not resolve' signal to tell the two apart. Read each title and "
        "abstract before relying on it. Every record carries `quotable_text`: the exact retrieved "
        "text. Any claim you draw from an abstract must quote a verbatim span of that field, "
        "because quotes are checked against it programmatically -- a paraphrase will not pass. "
        "`has_abstract` is false for records where PubMed holds only a title. Check `notes` for how "
        "many articles matched beyond the ones returned."
    ),
)
def pubmed_abstracts(
    gene: str, organism: str = MTB_H37RV_NAME, limit: int = DEFAULT_LIMIT
) -> ToolResult:
    """PubMed titles and abstracts for one gene.

    Args:
        gene: Gene symbol or locus tag, e.g. 'katG' or 'Rv1908c'.
        organism: Organism name added to the query. Defaults to 'Mycobacterium
            tuberculosis'. Pass an empty string to search without organism context.
        limit: Maximum articles to return, 1-100, ranked by PubMed relevance.
    """
    return _pubmed.abstracts(gene=gene, organism=organism, limit=limit)


def main() -> None:
    mcp.run(transport="stdio")


if __name__ == "__main__":
    main()
