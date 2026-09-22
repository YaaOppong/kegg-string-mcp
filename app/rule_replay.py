"""Recompute a rule classification from a captured fixture. No network, no model.

Standard-library only, like `replay.py`, so it runs on bare Pyodide with nothing
installed -- which is what lets the page compute rather than replay.

That difference is the point. The gene demo replays a model's summary and
re-validates it, because the summary cannot be regenerated without an API key.
A rule classification has no model in it: given the same catalogue rows, the same
annotation and the same rules, the verdicts are arithmetic. So the page ships the
inputs and runs the real classifier over them, and a verdict it shows is one the
library produces rather than one that was stored.

The only thing defined here rather than imported is a stand-in for a catalogue
row. The real one lives in `resistance.py`, which imports pydantic for its tool
envelope; the classification never constructs one and only reads four fields, so
a plain record satisfies it. The grades that count as associated come from the
fixture, captured from the real module at build time, so the two cannot drift.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from kegg_string_mcp.rules.annotation import parse as parse_annotation
from kegg_string_mcp.rules.catalogue import Catalogue
from kegg_string_mcp.rules.enrichment import universe_for
from kegg_string_mcp.rules.evidence import Sources, evidence_for, rename_map
from kegg_string_mcp.rules.nesting import nest
from kegg_string_mcp.rules.nesting import summarise as summarise_nesting
from kegg_string_mcp.rules.parse import parse as parse_rules
from kegg_string_mcp.rules.questions import generate as generate_questions
from kegg_string_mcp.rules.sets import set_evidence
from kegg_string_mcp.rules.signature import Link, classify

FIXTURE = Path(__file__).resolve().parent.parent / "demo" / "rule_fixture.json"


@dataclass(frozen=True)
class Row:
    """A catalogue row, as the classification reads one."""

    gene: str
    mutation: str
    drug: str
    confidence: str
    source: str = ""
    comment: str = ""
    associated_grades: tuple[str, ...] = ()

    @property
    def associated(self) -> bool:
        return self.confidence in self.associated_grades


@dataclass(frozen=True)
class Position:
    """A lineage-defining position, as `Sources.lineages_in` reads one."""

    position: int
    lineage: str
    lineage_name: str = ""
    allele: str = ""


def load(path: str | Path | None = None) -> dict[str, Any]:
    return json.loads(Path(path or FIXTURE).read_text(encoding="utf-8"))


def sources_from(fixture: dict[str, Any], workdir: str | Path = "/tmp") -> Sources:
    """Build the real `Sources` from captured inputs.

    The annotation goes through the real parser rather than being reconstructed,
    so the page exercises the same coordinate handling, the same suffix rule and
    the same intergenic construction the pipeline uses.
    """
    grades = tuple(fixture.get("associated_grades", ()))
    path = Path(workdir) / "fixture_annotation.gff"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(fixture["annotation_gff"], encoding="utf-8")
    annotation = parse_annotation(path)

    # `gene` repeats the key it is stored under, so it is restored rather than
    # shipped 5,560 times.
    catalogue = {gene: [Row(gene=gene, associated_grades=grades, **row) for row in rows]
                 for gene, rows in fixture.get("catalogue", {}).items()}
    barcode = [Position(**row) for row in fixture.get("barcode", [])]
    return Sources(annotation=annotation, catalogue=Catalogue(catalogue, annotation),
                   barcode=barcode)


def classified(fixture: dict[str, Any], workdir: str | Path = "/tmp") -> dict[str, Any]:
    """Every rule in the fixture, classified. The page renders this."""
    sources = sources_from(fixture, workdir)
    path = Path(workdir) / "fixture_rules.tsv"
    path.write_text(fixture["rules_tsv"], encoding="utf-8")
    rules = parse_rules(path)

    rename = rename_map(rules, sources)
    universe = universe_for(sources.annotation)
    captured = fixture.get("universe")
    if captured:
        # Keep the real background: the membership map covers the loci held, the
        # counts come from the annotation the pipeline used.
        universe.sizes = dict(captured["sizes"])
        universe.total = int(captured["total"])
    nesting = nest(rules, rename)

    # Interaction evidence is captured, not recomputed: it is the one input that
    # cannot be derived from a file. Without it a compensation verdict would
    # claim no link was found, which is a different statement from the pipeline's.
    partners = {locus: set(names) for locus, names in fixture.get("partners", {}).items()}
    best: dict[tuple[str, str], Link] = {}
    for entry in fixture.get("links", []):
        best.setdefault((entry["a"], entry["b"]), Link(entry["kind"], entry.get("detail", "")))
    linked = {pair for pair, link in best.items() if link.kind.startswith("string")}

    signatures = []
    for rule in rules:
        conditions = evidence_for(rule, sources)
        signatures.append(classify(
            rule, conditions, best,
            sets=set_evidence(conditions, sources.annotation, universe=universe,
                              partners=partners, linked_pairs=linked)))

    counts: dict[str, int] = {}
    for signature in signatures:
        counts[signature.primary] = counts.get(signature.primary, 0) + 1
    return {
        "rules": [signature.to_dict(rename)
                  | nesting[rule.rule_id(rename)].to_dict()
                  # `roles` is added by the TSV writer rather than by to_dict, so
                  # it has to be added here too or the two paths disagree.
                  | {"roles": "|".join(f"{c.locus}:{c.role}" for c in signature.conditions)
                     or "NA"}
                  for rule, signature in zip(rules, signatures)],
        "by_signature": dict(sorted(counts.items())),
        "questions": [q.to_dict() for q in generate_questions(signatures)],
        "nesting": summarise_nesting(nesting),
        "annotation_sha256": fixture.get("annotation_sha256", ""),
        "loci": len(sources.annotation.genes),
        "terms": universe.terms,
    }
