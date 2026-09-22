"""Run the deterministic half of the workflow, end to end.

    rules.tsv + annotation.bed  ->  loci.tsv + rules.tsv

No model. Everything here is a lookup or a count, so the same inputs give the
same tables and a reader can check any classification against the locus row it
rests on.

Order matters for cost. The free work -- resolution, the catalogue, the lineage
barcode, the population counts -- happens first and is enough to classify most
rules. The per-locus fetches follow, once per distinct locus rather than once per
condition. The single `network()` call over the resolved loci comes last and is
what upgrades a compensation-shaped rule to a linked one.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

from kegg_string_mcp.rules.annotate import (
    LinkedLocus,
    annotate_locus,
    co_occurring,
    compute_links,
    population_counts,
)
from kegg_string_mcp.rules.annotation import parse as parse_annotation
from kegg_string_mcp.rules.catalogue import ANCHOR
from kegg_string_mcp.rules.enrichment import universe_for
from kegg_string_mcp.rules.evidence import evidence_for, load_sources, rename_map, summarise
from kegg_string_mcp.rules.nesting import nest
from kegg_string_mcp.rules.nesting import summarise as summarise_nesting
from kegg_string_mcp.rules.parse import parse as parse_rules
from kegg_string_mcp.rules.questions import generate as generate_questions
from kegg_string_mcp.rules.questions import summarise as summarise_questions
from kegg_string_mcp.rules.report import write_loci, write_questions, write_rules
from kegg_string_mcp.rules.sets import set_evidence
from kegg_string_mcp.rules.signature import classify


@dataclass
class RunResult:
    rules: list[Any]
    signatures: list[Any]
    annotations: dict[str, Any]
    written: dict[str, Path]
    summary: dict[str, Any]
    questions: list[Any]
    notes: list[str]


def _edges(annotations: dict[str, Any], clients: Any, notes: list[str]
           ) -> dict[tuple[str, str], dict[str, Any]]:
    """Every STRING edge among the loci, in one call, keyed by locus.

    Asking whether A and B interact by looking for B in A's top-20 partner list
    is only reliable when neither list is full; the pair-level query has no such
    cut. One call covers every pair.
    """
    string = getattr(clients, "string", None)
    by_string_id = {a.string_id: locus for locus, a in annotations.items() if a.string_id}
    if string is None or len(by_string_id) < 2:
        notes.append("no STRING edges computed: fewer than two loci resolved in STRING")
        return {}
    try:
        network = string.network(sorted(by_string_id), required_score=150).model_dump()
    except Exception as exc:                          # noqa: BLE001
        notes.append(f"STRING network unavailable ({type(exc).__name__}: {exc}); links fall "
                     f"back to partner lists and absence of an edge is not a negative")
        return {}

    out: dict[tuple[str, str], dict[str, Any]] = {}
    for record in network.get("records", []):
        detail = record.get("detail", {})
        left = by_string_id.get(str(detail.get("string_id_a", "")))
        right = by_string_id.get(str(detail.get("string_id_b", "")))
        if left and right:
            out[(left, right)] = detail
    notes.append(f"{len(out)} STRING edge(s) among {len(by_string_id)} resolved loci")
    return out


def run(rules_path: str | Path, annotation_path: str | Path, out_dir: str | Path,
        clients: Any, resistance: Any, lineage: Any = None, organism: str = "mtu",
        genome_size: int = 0, pathway_sizes: dict[str, int] | None = None) -> RunResult:
    annotation = parse_annotation(annotation_path)
    rules = parse_rules(rules_path)
    sources = load_sources(annotation, resistance, lineage, organism)
    notes = list(sources.notes)
    notes.append(f"annotation {annotation.sha256[:16]}… {len(annotation.genes):,} genes")

    rename = rename_map(rules, sources)

    # One annotation per distinct locus, with the rule-file spellings that
    # reached it recorded so the table can be compared back to the rules.
    labels_for: dict[str, list[str]] = {}
    for label, locus in rename.items():
        labels_for.setdefault(locus, []).append(label)

    annotations: dict[str, Any] = {}
    for locus, labels in sorted(labels_for.items()):
        feature, _, _ = sources.feature_evidence(labels[0])
        symbol = feature.gene.symbol if feature.gene is not None else ""
        record = annotate_locus(locus, symbol, sources, clients)
        record.labels = tuple(sorted(labels))
        annotations[locus] = record

    anchors = {locus for locus, a in annotations.items()
               if a.catalogue is not None and a.catalogue.status == ANCHOR}
    population_counts(rules, rename, anchors, annotations)

    pairs = co_occurring(rules, rename)
    links = compute_links(pairs, annotations, sources,
                          pathway_sizes=pathway_sizes, genome_size=genome_size,
                          edges=_edges(annotations, clients, notes))

    # The same link belongs to both loci, rendered once so the two rows cannot
    # disagree about it.
    by_locus: dict[str, list[str]] = {}
    for (left, right), found in links.items():
        for link in found:
            by_locus.setdefault(left, []).append(LinkedLocus(right, link.kind, link.detail).render())
            by_locus.setdefault(right, []).append(LinkedLocus(left, link.kind, link.detail).render())

    best = {pair: found[0] for pair, found in links.items() if found}

    # Built once for the population rather than per rule: the enrichment universe
    # is the whole annotation, and the partner lists are already fetched.
    universe = universe_for(annotation)
    if not universe.size:
        notes.append("the supplied annotation carries no functional categories, so no "
                     "enrichment was computed. Not a negative result.")
    else:
        notes.append(f"enrichment universe: {universe.size:,} annotated loci, "
                     f"{universe.terms} category terms")
    partner_sets = {locus: {name or pid for pid, name, _ in a.partners}
                    for locus, a in annotations.items()}
    # Pairs STRING links, for the subgraph shape. Co-mention only still counts as
    # an edge here -- the shape is a description, not a claim about support.
    linked_pairs = {pair for pair, found in links.items()
                    if any(link.kind.startswith("string") for link in found)}

    # Containment across the population. One pass, no lookups.
    nesting = nest(rules, rename)
    notes.append(
        f"{summarise_nesting(nesting)['minimal']:,} of {len(rules):,} rules are minimal; "
        f"the rest contain a smaller rule in the population")

    signatures = []
    for rule in rules:
        conditions = evidence_for(rule, sources)
        signatures.append(classify(
            rule, conditions, best,
            sets=set_evidence(conditions, annotation, universe=universe,
                              partners=partner_sets, linked_pairs=linked_pairs)))

    questions = generate_questions(signatures)

    out_dir = Path(out_dir)
    written = {"loci": write_loci(out_dir / "loci.tsv", list(annotations.values()), by_locus),
               "rules": write_rules(out_dir / "rules.tsv", signatures, rename, nesting),
               "questions": write_questions(out_dir / "questions.tsv", questions)}
    summary = summarise(signatures) | summarise_questions(questions) | summarise_nesting(nesting)
    return RunResult(rules=rules, signatures=signatures, annotations=annotations,
                     written=written, summary=summary, questions=questions, notes=notes)
