"""Gather the structured evidence each rule condition rests on.

Everything here is a lookup against sources fetched once per run and shared
across every rule. A classifier population reuses its loci heavily -- the same
locus appears in many rules -- so per-condition work is cached by feature, and a
run over a few hundred rules costs a few hundred lookups rather than a few
thousand.

Lineage is taken by interval rather than by gene. `lineage_markers` answers for a
gene, which is what a reader wants; a feature here may be an intergenic interval,
and a lineage-defining position either falls inside a span or it does not. One
rule covers both kinds and needs no special case.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from kegg_string_mcp.lineage import LineageSnp
from kegg_string_mcp.rules.annotation import Annotation
from kegg_string_mcp.rules.catalogue import Catalogue, CatalogueStatus
from kegg_string_mcp.rules.features import CODING, INTERGENIC, Feature, Resolver
from kegg_string_mcp.rules.parse import Rule
from kegg_string_mcp.rules.signature import ConditionEvidence, Link, RuleSignature, classify


@dataclass
class Sources:
    """The three structured sources a condition is read against.

    `barcode` is the tbdb lineage list; `catalogue` the WHO grading. Both are
    H37Rv-specific and neither generalises to another organism, which is why this
    stays a TB tool while the annotation is an input.
    """

    annotation: Annotation
    catalogue: Catalogue
    barcode: list[LineageSnp] = field(default_factory=list)
    organism: str = "mtu"
    notes: list[str] = field(default_factory=list)

    _resolver: Resolver | None = None
    _cache: dict[str, tuple[Feature, CatalogueStatus, tuple[str, ...]]] = \
        field(default_factory=dict, repr=False)

    @property
    def resolver(self) -> Resolver:
        if self._resolver is None:
            self._resolver = Resolver(self.annotation)
        return self._resolver

    def lineages_in(self, start: int, end: int) -> tuple[str, ...]:
        return tuple(sorted({s.lineage for s in self.barcode
                             if start <= s.position <= end and s.lineage}))

    def feature_evidence(self, label: str) -> tuple[Feature, CatalogueStatus, tuple[str, ...]]:
        """Resolve a label and look it up, once per distinct label per run."""
        if label in self._cache:
            return self._cache[label]

        feature = self.resolver.resolve(label)
        if feature.kind == CODING and feature.gene is not None:
            status = self.catalogue.for_locus(feature.gene.locus, feature.gene.symbol)
            lineages = self.lineages_in(feature.gene.start, feature.gene.end)
        elif feature.kind == INTERGENIC and feature.interval is not None:
            status = self.catalogue.for_interval(feature.interval)
            lineages = self.lineages_in(feature.interval.start, feature.interval.end)
        else:
            # An unresolved label has no coordinates, so it has no catalogue or
            # lineage answer either. Reporting "absent from the catalogue" here
            # would state a fact about a gene nothing has identified.
            status = CatalogueStatus(
                "absent", note="label did not resolve, so no catalogue lookup was made")
            lineages = ()
        self._cache[label] = (feature, status, lineages)
        return self._cache[label]


def load_sources(annotation: Annotation, resistance: Any, lineage: Any = None,
                 organism: str = "mtu") -> Sources:
    """Fetch the two TB-specific sources and build a `Sources`.

    Both are single files, fetched once and cached on disk. The lineage barcode
    was previously never loaded -- `lineages_in` worked but every real run saw an
    empty list, so the population-structure confound could not fire.

    Fail-soft, separately: a missing barcode must not cost the catalogue, and a
    run without one reports that it could not check the confound rather than
    reporting that there is none.
    """
    from kegg_string_mcp.lineage import BARCODE_URL, parse_barcode
    from kegg_string_mcp.rules.catalogue import Catalogue

    notes: list[str] = []
    catalogue_rows, _ = resistance.catalogue()
    sources = Sources(annotation=annotation, catalogue=Catalogue(catalogue_rows, annotation),
                      organism=organism, notes=notes)
    if lineage is None:
        notes.append("no lineage client supplied, so population structure was not checked")
        return sources
    try:
        sources.barcode = parse_barcode(lineage.http.get(BARCODE_URL).body)
        notes.append(f"lineage barcode: {len(sources.barcode)} defining positions")
    except Exception as exc:                       # noqa: BLE001
        notes.append(f"lineage barcode unavailable ({type(exc).__name__}: {exc}), so the "
                     f"population-structure confound could not be checked. This is not "
                     f"evidence that no rule is confounded.")
    return sources


def evidence_for(rule: Rule, sources: Sources) -> list[ConditionEvidence]:
    out = []
    for condition in rule.conditions:
        feature, status, lineages = sources.feature_evidence(condition.label)
        out.append(ConditionEvidence(condition=condition, feature=feature,
                                     catalogue=status, lineages=lineages))
    return out


def classify_all(rules: list[Rule], sources: Sources,
                 links: dict[tuple[str, str], Link] | None = None) -> list[RuleSignature]:
    return [classify(rule, evidence_for(rule, sources), links) for rule in rules]


def rename_map(rules: list[Rule], sources: Sources) -> dict[str, str]:
    """label -> resolved locus, for rule ids that survive a respelling.

    A label that does not resolve maps to itself: giving two unresolved labels
    the same id would merge two rules on the strength of a failed lookup.
    """
    out: dict[str, str] = {}
    for rule in rules:
        for condition in rule.conditions:
            if condition.label in out:
                continue
            feature, _, _ = sources.feature_evidence(condition.label)
            out[condition.label] = feature.locus or condition.label
    return out


def summarise(signatures: list[RuleSignature]) -> dict[str, Any]:
    counts: dict[str, int] = {}
    for signature in signatures:
        counts[signature.primary] = counts.get(signature.primary, 0) + 1
    return {"rules": len(signatures), "by_primary_signature": dict(sorted(counts.items()))}
