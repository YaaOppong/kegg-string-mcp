"""Ties the pieces together for both modes.

Sequence, and the order matters:

  1. pipeline fetches (single mode: on demand by the model; epistasis: up front)
  2. pipeline computes the deterministic evidence and writes the store
  3. model summarises, seeing only what the store already holds
  4. pipeline validates every citation against the store

Validation runs on the pipeline's copy of what the tools returned, never on the
model's account of it -- otherwise a model that misreported its own retrieval
would validate against its own mistake.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from kegg_string_mcp.agent.evidence import all_pairs
from kegg_string_mcp.agent.loop import new_store, run_loop
from kegg_string_mcp.agent.validate import validate
from kegg_string_mcp.cache import DiskCache
from kegg_string_mcp.http import PoliteClient
from kegg_string_mcp.kegg import KeggClient
from kegg_string_mcp.pubmed import PubMedClient
from kegg_string_mcp.string_db import StringClient


class Tools:
    """Dispatch table shared by the loop and by the pre-fetch step."""

    def __init__(self, http: PoliteClient | None = None):
        self.http = http or PoliteClient(DiskCache())
        self.kegg = KeggClient(self.http)
        self.string = StringClient(self.http)
        self.pubmed = PubMedClient(self.http)

    def __call__(self, name: str, arguments: dict[str, Any]) -> dict[str, Any]:
        if name == "kegg_pathways":
            return self.kegg.pathways(**arguments).model_dump()
        if name == "string_partners":
            return self.string.partners(**arguments).model_dump()
        if name == "pubmed_abstracts":
            return self.pubmed.abstracts(**arguments).model_dump()
        return {"records": [], "record_ids": [], "notes": [f"unknown tool '{name}'"]}


def annotate_gene(gene: str, organism: str = "mtu", runs: Path = Path("runs"),
                  tools: Tools | None = None, client: Any | None = None) -> dict[str, Any]:
    tools = tools or Tools()
    store = new_store(runs, f"single-{gene}")

    task = (f"Annotate the gene '{gene}' in KEGG organism '{organism}' "
            f"(NCBI taxon 83332 for STRING). Describe its function, pathways and "
            f"notable interaction partners, citing record IDs from the tool results.")
    result = run_loop("single", task, tools, store, client=client)

    report = validate(result.text, store.citable_ids, store.per_target, gene.strip().upper(),
                      records=store.records)
    payload = {"mode": "single", "gene": gene, "organism": organism,
               "summary": result.text, "turns": result.turns,
               "stop_reason": result.stop_reason, "validation": report.to_dict()}
    store.output(payload)
    return payload | {"store": str(store.path)}


def annotate_epistasis(genes: list[str], organism: str = "mtu", runs: Path = Path("runs"),
                       tools: Tools | None = None, client: Any | None = None) -> dict[str, Any]:
    """Pre-compute every pairwise relationship, then let the model interpret it.

    The set intersections are done here rather than by the model: a model doing
    arithmetic over a hundred identifiers will get some of it wrong, and the
    error will not be visible in the prose.
    """
    tools = tools or Tools()
    store = new_store(runs, "epistasis")

    pathways: dict[str, list[dict[str, Any]]] = {}
    partners: dict[str, list[dict[str, Any]]] = {}
    for gene in genes:
        kegg_result = tools("kegg_pathways", {"gene": gene, "organism": organism})
        store.tool_result("kegg_pathways", {"gene": gene, "organism": organism}, kegg_result)
        pathways[gene] = kegg_result.get("records", [])

        string_result = tools("string_partners", {"gene": gene})
        store.tool_result("string_partners", {"gene": gene}, string_result)
        partners[gene] = string_result.get("records", [])

    sizes, _ = tools.kegg.pathway_sizes(organism)
    genome_size = _annotated_gene_count(tools, organism)
    store.derived("pathway_sizes", {"organism": organism, "n_pathways": len(sizes),
                                    "genome_size": genome_size})

    pairs = all_pairs(genes, pathways, partners, sizes, genome_size)
    store.derived("pair_evidence", {"pairs": [p.to_dict() for p in pairs]})

    table = "\n\n".join(
        f"PAIR {p.gene_a} / {p.gene_b}\n"
        f"  deterministic verdict: {p.verdict}\n"
        f"  degrees: {p.degrees}\n"
        f"  direct interaction: {p.direct_interaction}\n"
        f"  shared pathways: " + (", ".join(
            f"{sp.pathway_id} '{sp.name}' [{sp.specificity}, {sp.size} genes]"
            for sp in p.shared_pathways) or "none") + "\n"
        f"  shared partners: " + (", ".join(
            f"{sp['record_id']} ({sp['name']})" for sp in p.shared_partners) or "none")
        for p in pairs
    )

    task = (f"Genes flagged as interacting by an upstream analysis: {', '.join(genes)} "
            f"(KEGG organism '{organism}').\n\n"
            f"Pre-computed pairwise evidence:\n\n{table}\n\n"
            f"Interpret these relationships. Do not contradict the deterministic verdicts.")
    result = run_loop("epistasis", task, tools, store, client=client)

    report = validate(result.text, store.citable_ids, records=store.records)
    payload = {"mode": "epistasis", "genes": genes, "organism": organism,
               "summary": result.text, "turns": result.turns,
               "stop_reason": result.stop_reason,
               "pairs": [p.to_dict() for p in pairs],
               "validation": report.to_dict()}
    store.output(payload)
    return payload | {"store": str(store.path)}


def _annotated_gene_count(tools: Tools, organism: str) -> int:
    """Denominator for the 'is this pathway a container?' judgement."""
    index, _ = tools.kegg.gene_index(organism)
    return len(index.locus_tags)
