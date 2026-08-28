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

import json
from pathlib import Path
from typing import Any

from kegg_string_mcp.agent.evidence import all_pairs
from kegg_string_mcp.agent.loop import new_store, run_loop
from kegg_string_mcp.agent.store import RunStore
from kegg_string_mcp.agent.validate import validate
from kegg_string_mcp.cache import DiskCache
from kegg_string_mcp.http import FetchError, PoliteClient
from kegg_string_mcp.kegg import KeggClient
from kegg_string_mcp.pubmed import PubMedClient
from kegg_string_mcp.string_db import StringClient
from kegg_string_mcp.uniprot import UniProtClient

# Model-supplied arguments are untrusted input. Splatting them into typed clients
# turned a schema deviation -- limit="20", or organism= passed to string_partners --
# into a TypeError that killed the run mid-flight, instead of an error envelope the
# model could correct from, which is this codebase's rule for bad arguments.
TOOL_PARAMS: dict[str, dict[str, type]] = {
    "kegg_pathways": {"gene": str, "organism": str},
    "string_partners": {"gene": str, "species": int, "limit": int, "required_score": int},
    "pubmed_abstracts": {"gene": str, "organism": str, "limit": int},
    "uniprot_protein": {"gene": str, "organism_id": int, "limit": int},
}


def _coerce(name: str, arguments: dict[str, Any]) -> tuple[dict[str, Any], list[str]]:
    spec = TOOL_PARAMS[name]
    clean: dict[str, Any] = {}
    problems: list[str] = []
    for key, value in arguments.items():
        if key not in spec:
            problems.append(f"'{key}' is not a parameter of {name} "
                            f"(accepts: {', '.join(sorted(spec))})")
            continue
        try:
            clean[key] = spec[key](value)
        except (TypeError, ValueError):
            problems.append(f"'{key}' must be {spec[key].__name__}, got {value!r}")
    return clean, problems


class Tools:
    """Dispatch table shared by the loop and by the pre-fetch step."""

    def __init__(self, http: PoliteClient | None = None):
        self.http = http or PoliteClient(DiskCache())
        self.kegg = KeggClient(self.http)
        self.string = StringClient(self.http)
        self.pubmed = PubMedClient(self.http)
        self.uniprot = UniProtClient(self.http)

    def __call__(self, name: str, arguments: dict[str, Any]) -> dict[str, Any]:
        if name not in TOOL_PARAMS:
            return {"query": {"tool": name}, "records": [], "record_ids": [],
                    "notes": [f"unknown tool '{name}'"]}
        clean, problems = _coerce(name, arguments)
        if problems:
            return {"query": dict(arguments), "records": [], "record_ids": [],
                    "notes": [(f"Invalid argument(s), so no lookup was performed: "
                              f"{'; '.join(problems)}. An empty result here does NOT mean "
                              f"there is no data.")]}
        method = {"kegg_pathways": self.kegg.pathways,
                  "string_partners": self.string.partners,
                  "pubmed_abstracts": self.pubmed.abstracts,
                  "uniprot_protein": self.uniprot.protein}[name]
        return method(**clean).model_dump()


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
    return _finish(store, payload)


PARTNER_LIMIT = 20


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

        string_args = {"gene": gene, "limit": PARTNER_LIMIT}
        string_result = tools("string_partners", string_args)
        store.tool_result("string_partners", string_args, string_result)
        partners[gene] = string_result.get("records", [])

    # The only upstream calls in the pipeline that are not already fail-soft.
    # One KEGG 5xx here discarded every per-gene fetch already made and crashed.
    try:
        sizes, _ = tools.kegg.pathway_sizes(organism)
        genome_size = _annotated_gene_count(tools, organism)
        size_note = ""
    except (FetchError, ValueError) as exc:
        sizes, genome_size = {}, 0
        size_note = (f"Pathway sizes could not be fetched for '{organism}' ({exc}), so shared "
                     f"pathways cannot be judged as specific or broad. A shared pathway below "
                     f"may be a container category covering much of the genome.")
    store.derived("pathway_sizes", {"organism": organism, "n_pathways": len(sizes),
                                    "genome_size": genome_size})

    pairs = all_pairs(genes, pathways, partners, sizes, genome_size, PARTNER_LIMIT)
    store.derived("pair_evidence", {"pairs": [p.to_dict() for p in pairs]})

    table = "\n\n".join(
        f"PAIR {p.gene_a} / {p.gene_b}\n"
        f"  deterministic verdict: {p.verdict}\n"
        f"  partners retrieved: {p.partners_retrieved}"
        + (f" (LIMIT REACHED for {', '.join(p.truncated)}; true degree unknown)"
           if p.truncated else "") + "\n"
        f"  direct interaction: {p.direct_interaction}\n"
        f"  shared pathways: " + (", ".join(
            f"{sp.pathway_id} '{sp.name}' [{sp.specificity}, {sp.size} genes]"
            for sp in p.shared_pathways) or "none") + "\n"
        "  shared partners: " + (", ".join(
            f"{sp['record_id']} ({sp['name']})" for sp in p.shared_partners) or "none")
        for p in pairs
    )

    if size_note:
        table = f"WARNING: {size_note}\n\n{table}"

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
    return _finish(store, payload)


def _write_manifest(store: RunStore) -> Path | None:
    """Emit the corpus manifest beside the run store, for a downstream full-text
    pipeline to consume. Written only when the run actually found papers."""
    papers = store.corpus_manifest()
    if not papers:
        return None
    path = store.path.with_suffix(".corpus.jsonl")
    with path.open("w", encoding="utf-8") as fh:
        for paper in papers:
            fh.write(json.dumps(paper) + "\n")
    return path


def _finish(store: RunStore, payload: dict[str, Any]) -> dict[str, Any]:
    manifest = _write_manifest(store)
    payload = payload | {"corpus": [p for p in store.corpus_manifest()]}
    store.output(payload)
    return payload | {"store": str(store.path),
                      "corpus_manifest": str(manifest) if manifest else None}


def _annotated_gene_count(tools: Tools, organism: str) -> int:
    """Denominator for the 'is this pathway a container?' judgement."""
    index, _ = tools.kegg.gene_index(organism)
    return len(index.locus_tags)
