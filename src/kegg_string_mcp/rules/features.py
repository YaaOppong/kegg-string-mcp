"""Resolve a rule's feature labels against the annotation that named them.

A feature is either a coding locus (`Rv1908c`, or a symbol like `rpoC` where the
annotation carried one) or an intergenic interval named for its flanking genes
(`Rv1482c-Rv1483`). Everything downstream -- catalogue lookup, lineage, pathway
membership, the coordinate hits inside an interval -- keys off what this module
returns, so an unresolved label must stay unresolved rather than being guessed
into the nearest neighbour.

Two rules, both learnt from a real vocabulary of 2,903 labels:

* **Match the strand suffix loosely, the rest exactly.** 1,346 of those labels
  disagreed with a second annotation only in the trailing `c`, in both
  directions. A locus tag that differs only there is the same gene; one that
  differs anywhere else is not.
* **Refuse ambiguity, record it.** If a normalised tag matches more than one
  gene, no answer is given. This is `identity.py`'s contract applied to a
  vocabulary rather than to a query, and for the same reason: a silently wrong
  locus is worse than a missing one, because nothing downstream can detect it.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from kegg_string_mcp.rules.annotation import (
    Annotation,
    Gene,
    Intergenic,
    is_locus_tag,
    normalise_locus,
)

CODING = "coding"
INTERGENIC = "intergenic"
UNRESOLVED = "unresolved"


@dataclass(frozen=True)
class Feature:
    """One resolved feature label.

    `matched_by` is kept for the same reason `identity.py` keeps it: a label that
    resolved exactly and one that resolved only after ignoring a strand suffix
    are different levels of confidence, and collapsing them hides which
    assumption a downstream number rests on.
    """

    label: str
    kind: str                       # CODING | INTERGENIC | UNRESOLVED
    matched_by: str                 # exact | strand_suffix | symbol | flanking | none
    gene: Gene | None = None
    interval: Intergenic | None = None
    candidates: tuple[str, ...] = ()   # populated when a match was refused
    note: str = ""

    @property
    def resolved(self) -> bool:
        return self.kind != UNRESOLVED

    @property
    def locus(self) -> str:
        if self.gene is not None:
            return self.gene.locus
        if self.interval is not None:
            return self.interval.name
        return ""

    def to_dict(self) -> dict[str, Any]:
        out: dict[str, Any] = {"label": self.label, "kind": self.kind,
                               "matched_by": self.matched_by, "note": self.note}
        if self.gene is not None:
            out |= {"locus": self.gene.locus, "symbol": self.gene.symbol,
                    "start": self.gene.start, "end": self.gene.end,
                    "strand": self.gene.strand, "length": self.gene.length}
        if self.interval is not None:
            out |= {"locus": self.interval.name,
                    "left": self.interval.left.locus, "right": self.interval.right.locus,
                    "start": self.interval.start, "end": self.interval.end,
                    "width": self.interval.width,
                    "promoter_of": "|".join(g.locus for g in self.interval.promoter_of())}
        if self.candidates:
            out["candidates"] = "|".join(self.candidates)
        return out


def _unresolved(label: str, note: str, candidates: tuple[str, ...] = ()) -> Feature:
    return Feature(label=label, kind=UNRESOLVED, matched_by="none",
                   candidates=candidates, note=note)


class Resolver:
    """Resolves labels against one annotation, and remembers what it decided."""

    def __init__(self, annotation: Annotation):
        self.annotation = annotation
        # Intergenic intervals addressable by their flanking pair with strand
        # suffixes ignored, so a gap named `Rv1482c-Rv1483` in one vocabulary
        # still resolves when the annotation spells a flank the other way.
        self._by_flanks: dict[tuple[str, str], list[Intergenic]] = {}
        for name in annotation.intergenic_names:
            interval = annotation.intergenic(name)
            assert interval is not None
            key = (normalise_locus(interval.left.locus), normalise_locus(interval.right.locus))
            self._by_flanks.setdefault(key, []).append(interval)

    def resolve(self, label: str) -> Feature:
        label = label.strip()
        if not label:
            return _unresolved(label, "empty label")
        if "-" in label:
            return self._intergenic(label)
        return self._coding(label)

    # -- coding --------------------------------------------------------------

    def _coding(self, label: str) -> Feature:
        gene = self.annotation.gene(label)
        if gene is not None:
            return Feature(label, CODING, "exact", gene=gene)

        if is_locus_tag(label):
            hits = self.annotation.by_normalised(label)
            if len(hits) == 1:
                return Feature(
                    label, CODING, "strand_suffix", gene=hits[0],
                    note=(f"matched {hits[0].locus} once the complement-strand suffix is "
                          f"ignored; the two vocabularies disagree on the suffix only"))
            if len(hits) > 1:
                return _unresolved(
                    label, "more than one gene matches once the suffix is ignored, so no "
                           "match is made", tuple(g.locus for g in hits))

        by_symbol = self.annotation.by_symbol(label)
        if len(by_symbol) == 1:
            return Feature(label, CODING, "symbol", gene=by_symbol[0])
        if len(by_symbol) > 1:
            return _unresolved(label, "the symbol names more than one gene in this annotation",
                               tuple(g.locus for g in by_symbol))

        return _unresolved(label, "no locus tag or symbol in the supplied annotation matches")

    # -- intergenic ----------------------------------------------------------

    def _intergenic(self, label: str) -> Feature:
        interval = self.annotation.intergenic(label)
        if interval is not None:
            return Feature(label, INTERGENIC, "exact", interval=interval)

        left, _, right = label.partition("-")
        if not left or not right:
            return _unresolved(label, "not a well-formed flanking pair")

        hits = self._by_flanks.get((normalise_locus(left), normalise_locus(right)), [])
        if len(hits) == 1:
            return Feature(
                label, INTERGENIC, "flanking", interval=hits[0],
                note=(f"matched {hits[0].name} on the flanking pair once complement-strand "
                      f"suffixes are ignored"))
        if len(hits) > 1:
            return _unresolved(label, "the flanking pair is not unique in this annotation",
                               tuple(i.name for i in hits))

        # Both flanks may exist while the gap does not: genes that overlap have
        # no interval between them, which is a fact about the genome rather than
        # a failure to resolve, and is worth saying so the caller can tell the
        # two apart.
        lf, rf = self.annotation.by_normalised(left), self.annotation.by_normalised(right)
        if lf and rf:
            return _unresolved(
                label, "both flanking genes exist but they are not consecutive with a gap "
                       "between them in this annotation")
        missing = ", ".join(n for n, f in ((left, lf), (right, rf)) if not f)
        return _unresolved(label, f"flanking gene(s) not in this annotation: {missing}")


@dataclass
class Coverage:
    """What a whole vocabulary resolved to. The pre-flight report."""

    features: dict[str, Feature] = field(default_factory=dict)

    def add(self, feature: Feature) -> None:
        self.features[feature.label] = feature

    def counts(self) -> dict[str, int]:
        out = {"total": len(self.features), CODING: 0, INTERGENIC: 0, UNRESOLVED: 0,
               "exact": 0, "strand_suffix": 0, "symbol": 0, "flanking": 0, "ambiguous": 0}
        for feature in self.features.values():
            out[feature.kind] += 1
            if feature.matched_by in out:
                out[feature.matched_by] += 1
            if feature.candidates:
                out["ambiguous"] += 1
        return out

    def unresolved(self) -> list[Feature]:
        return [f for f in self.features.values() if not f.resolved]


def resolve_all(labels: list[str], annotation: Annotation) -> Coverage:
    resolver = Resolver(annotation)
    coverage = Coverage()
    for label in labels:
        coverage.add(resolver.resolve(label))
    return coverage
