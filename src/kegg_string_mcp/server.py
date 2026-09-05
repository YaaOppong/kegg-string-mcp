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
from kegg_string_mcp.lineage import LineageClient
from kegg_string_mcp.provenance import ToolResult
from kegg_string_mcp.pubmed import DEFAULT_LIMIT, MTB_H37RV_NAME, PubMedClient
from kegg_string_mcp.resistance import ResistanceClient
from kegg_string_mcp.string_db import MTB_H37RV, StringClient
from kegg_string_mcp.uniprot import UniProtClient

mcp = MCPServer(
    name="kegg-string",
    version="0.1.0",
    instructions=(
        "Deterministic lookups against KEGG, STRING and PubMed for gene annotation. Every record "
        "returned carries a stable record_id and a resolvable URL. Cite only record_ids that "
        "appear in a tool result's record_ids list -- citations are checked programmatically "
        "against what was retrieved. Read the `notes` field: an empty `records` list "
        "may mean the identifier failed to resolve rather than that no data exists. When a claim "
        "rests on a PubMed article, quote a verbatim span of that record's `quotable_text` and "
        "attach it to the PMID -- quotes are checked against the retrieved text, so a paraphrase "
        "or a remembered fact will not pass even under a correct PMID. "
        "RESEARCH USE ONLY: this server is not a clinical decision support system and must "
        "not be used to guide patient care. It returns drug-resistance-associated genes; "
        "resistance interpretation for clinical purposes requires validated methods and "
        "expert review."
    ),
)

_http = PoliteClient(DiskCache())
_kegg = KeggClient(_http)
_string = StringClient(_http)
_uniprot = UniProtClient(_http)
_lineage = LineageClient(_http)
_resistance = ResistanceClient(_http)
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




@mcp.tool(
    annotations=READ_ONLY,
    description=(
        "Curated protein annotation from UniProt: function, catalytic activity, subunit "
        "structure, PDB cross-references. Use this when KEGG has no pathway for a gene -- KEGG "
        "assigns one to only 29% of M. tuberculosis genes, so 'no KEGG pathway' is usually an "
        "annotation gap rather than a fact about the protein. IMPORTANT: each function statement "
        "carries an evidence tier. 'experimental' statements were measured on this protein and "
        "list the PMIDs that support them; 'sequence_similarity', 'sequence_model', 'automatic' "
        "and 'imported' statements are INFERRED from a rule or a homologue and are not evidence "
        "about this gene. Function text is prose, so quote a verbatim span of quotable_text for "
        "any claim you draw from it."
    ),
)
def uniprot_protein(gene: str, organism_id: int = 83332, limit: int = 3) -> ToolResult:
    """Curated protein annotation for one gene.

    Args:
        gene: Gene symbol or locus tag (e.g. katG, Rv1908c).
        organism_id: NCBI taxon ID. Defaults to 83332 (M. tuberculosis H37Rv).
        limit: Maximum UniProt entries to return.
    """
    return _uniprot.protein(gene=gene, organism_id=organism_id, limit=limit)


@mcp.tool(
    annotations=READ_ONLY,
    description=(
        "Whether a gene contains a lineage-defining SNP from the M. tuberculosis SNP barcode "
        "(TB-Profiler tbdb, after Coll 2014 and Napier 2020). Use this when a gene comes out of "
        "a scan over clinical isolates, because an association between two lineage-marked genes "
        "can be population structure rather than biology. IMPORTANT: this is a lookup, not a "
        "verdict. A positive result says the gene CONTAINS a lineage-defining position -- 855 of "
        "4,008 M. tuberculosis genes do -- not that the variant you are asking about is one. "
        "Compare your variant's H37Rv coordinate against the `position` field to answer that. "
        "Whether an association survives conditioning on lineage is a question for genotype "
        "data and this tool cannot answer it. Records are structured, so cite the record_id; "
        "there is no text to quote."
    ),
)
def lineage_markers(gene: str, organism: str = "mtu") -> ToolResult:
    """Lineage-defining SNPs contained in one gene.

    Args:
        gene: Gene symbol or locus tag (e.g. phoP, Rv0757).
        organism: KEGG organism code, used for the gene coordinates. Defaults to 'mtu'.
    """
    return _lineage.markers(gene=gene, organism=organism)


@mcp.tool(
    annotations=READ_ONLY,
    description=(
        "WHO-graded drug-resistance variants for a gene, from the TB-Profiler tbdb catalogue "
        "(WHO catalogue of mutations v2). Returns whether the gene is resistance-associated, "
        "which drugs, and the graded variants themselves. A gene counts as "
        "resistance-associated if ANY of its variants is graded associated, however many are "
        "not. IMPORTANT: the gene-level flag does not grade a variant. Within katG, "
        "p.Ser315Thr is 'Assoc w R' while p.Arg463Leu, a common polymorphism, is explicitly "
        "'Not assoc w R' -- pass `mutation` to grade a specific one. Distinguish the three "
        "negatives: a gene absent from the catalogue was never assessed (it covers 74 genes "
        "chosen for resistance surveillance), a gene present with no associated variant WAS "
        "assessed, and a variant graded 'Uncertain significance' -- 70% of all rows -- is "
        "neither. Records are structured, so cite the record_id; there is no text to quote."
    ),
)
def resistance_variants(gene: str, mutation: str | None = None,
                        drug: str | None = None) -> ToolResult:
    """WHO-graded resistance variants for one gene.

    Args:
        gene: Gene symbol or locus tag as the catalogue names it (e.g. katG, Rv0678, gid).
        mutation: Optional HGVS variant to grade specifically (e.g. p.Ser315Thr, c.-15C>T).
        drug: Optional drug name to restrict to (e.g. isoniazid, bedaquiline).
    """
    return _resistance.variants(gene=gene, mutation=mutation, drug=drug)


def main() -> None:
    mcp.run(transport="stdio")


if __name__ == "__main__":
    main()
