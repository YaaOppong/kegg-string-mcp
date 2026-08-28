"""Run the pipeline over the gold set and score it."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from kegg_string_mcp.agent.pipeline import Tools, annotate_gene
from kegg_string_mcp.evaluate.gold import GoldSet, load
from kegg_string_mcp.evaluate.score import EvalReport, score_gene


def evaluate(gold: GoldSet | None = None, runs: Path = Path("runs/eval"),
             tools: Tools | None = None, client: Any | None = None,
             annotate=annotate_gene) -> EvalReport:
    gold = gold or load()
    report = EvalReport(organism=gold.organism, reference=gold.reference,
                        retrieved_on=gold.retrieved_on, coverage=gold.coverage)
    tools = tools or Tools()

    for entry in gold.genes:
        try:
            payload = annotate(entry.query, gold.organism, runs, tools, client)
            score = score_gene(entry, payload.get("summary") or "",
                               payload.get("validation", {}), gold.organism)
        except Exception as exc:                      # noqa: BLE001
            # One gene failing must not discard the whole evaluation; a crash is
            # itself a result and is reported rather than swallowed.
            score = score_gene(entry, "", {}, gold.organism)
            score.error = f"{type(exc).__name__}: {exc}"
        report.scores.append(score)

    return report


def render(report: EvalReport) -> str:
    """Human-readable table. The three blocks are kept apart on purpose: only two
    of them are ground truth, and citation integrity is the only one that is exact."""
    data = report.to_dict()
    lines = [
        f"Reference: {data['reference']} (retrieved {data['reference_retrieved_on']})",
        (f"KEGG assigns a pathway to {data['coverage']['genes_with_any_pathway']} of "
        f"{data['coverage']['genes_total']} {data['organism']} genes "
        f"({data['coverage']['genes_with_any_pathway'] / data['coverage']['genes_total']:.0%})."),
        "",
        f"{'gene':10} {'kind':4} {'expected':9} {'reported':9} {'hits':5} {'missed':7} {'cites':9} {'quotes':7}",
    ]
    for score in data["per_gene"]:
        kind = "neg" if score["negative_control"] else "pos"
        cites = f"{score['citations_verified']}/{score['citations_total']}"
        quotes = f"{score['quotes_verified']}/{score['quotes_total']}"
        flag = "  FABRICATED" if score["negative_control"] and score["abstained"] is False else ""
        err = f"  ERROR {score['error']}" if score["error"] else ""
        lines.append(
            f"{score['gene']:10} {kind:4} {len(score['expected']):<9} {len(score['reported']):<9} "
            f"{len(score['hits']):<5} {len(score['missed']):<7} {cites:9} {quotes:7}{flag}{err}"
        )

    fidelity, abstention, integrity = (data["retrieval_fidelity"], data["abstention"],
                                       data["citation_integrity"])
    lines += [
        "",
        "Retrieval fidelity (positive controls -- did it report what KEGG assigns?)",
        (f"  recall    {fidelity['recall']}   precision {fidelity['precision']}   "
        f"n={fidelity['n_genes']}"),
        "",
        "Abstention (negative controls -- KEGG assigns nothing; correct answer is 'none')",
        f"  abstained {abstention['rate']}   n={abstention['n_genes']}"
        + (f"   fabricated on: {', '.join(abstention['fabricated_on'])}"
           if abstention["fabricated_on"] else ""),
        "",
        "Citation integrity (exact, computed not judged)",
        (f"  citation precision {integrity['citation_precision']}   "
        f"quote precision {integrity['quote_precision']}"),
    ]
    if data["errors"]:
        lines += ["", "Errors:"] + [f"  {e['gene']}: {e['error']}" for e in data["errors"]]
    return "\n".join(lines)


def write(report: EvalReport, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(report.to_dict(), indent=2), encoding="utf-8")
    return path
