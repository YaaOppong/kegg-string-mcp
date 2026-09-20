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
    strip_version,
)

CODING = "coding"
INTERGENIC = "intergenic"
UNRESOLVED = "unresolved"

# How a producer names a non-coding region, mapped to what the name means.
# SnpEff's effect vocabulary supplies these three, so they are the default rather
# than a guess -- but they are a default, not a law: another variant caller will
# use other words, and a caller can pass its own. When nothing matches, the
# refusal names what was tried, so the next vocabulary is a configuration change
# rather than someone else's bug.
UPSTREAM, DOWNSTREAM, PAIR = "upstream", "downstream", "pair"
PREFIXES: dict[str, str] = {
    "upstream_": UPSTREAM,
    "downstream_": DOWNSTREAM,
    "intergenic_": PAIR,
}

# SnpEff calls an upstream_gene_variant within a window of the gene -- 5,000 bp
# by default -- and 3,048 of H37Rv's 3,049 intergenic intervals are narrower than
# that (median 81 bp). So `upstream_X` names a region that usually extends past
# the gap and into the neighbouring gene. Resolving it to the gap is right about
# the location and understates the extent, which the note says rather than the
# coordinates pretending otherwise.
_WINDOW_CAVEAT = ("the producer's upstream/downstream window is typically wider than this "
                  "interval, so these coordinates are the intergenic part of that region "
                  "rather than its full extent")


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

    def __init__(self, annotation: Annotation, prefixes: dict[str, str] | None = None):
        self.annotation = annotation
        self.prefixes = dict(PREFIXES if prefixes is None else prefixes)
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
        # Longest prefix first, so `intergenic_` is not shadowed by a shorter one
        # a caller happens to configure.
        for prefix in sorted(self.prefixes, key=len, reverse=True):
            if label.startswith(prefix):
                return self._region(label, label[len(prefix):], self.prefixes[prefix])
        if "-" in label:
            return self._intergenic(label)
        return self._coding(label)

    def _region(self, label: str, rest: str, kind: str) -> Feature:
        """A prefixed region name: `upstream_X`, `downstream_X`, `intergenic_a-b`."""
        if kind == PAIR:
            return self._intergenic(label, rest)

        target = self._coding(strip_version(rest))
        if target.gene is None:
            return _unresolved(
                label, f"the gene named by this region does not resolve: {target.note}")

        interval = (self.annotation.upstream_of(target.gene) if kind == UPSTREAM
                    else self.annotation.downstream_of(target.gene))
        if interval is None:
            side = "5'" if kind == UPSTREAM else "3'"
            return _unresolved(
                label,
                f"{target.gene.locus} has no intergenic interval {side} of it: its neighbour "
                f"abuts or overlaps it, so a variant called {kind} of it lies inside that "
                f"neighbour. This is a fact about the genome, not a failed lookup -- 830 of "
                f"H37Rv's 4,008 genes are like this.")
        return Feature(label, INTERGENIC, kind, interval=interval,
                       note=(f"the interval {'5' if kind == UPSTREAM else '3'}' of "
                             f"{target.gene.locus} ({target.gene.strand} strand); "
                             f"{_WINDOW_CAVEAT}"))

    # -- coding --------------------------------------------------------------

    def _coding(self, label: str) -> Feature:
        gene = self.annotation.gene(label)
        if gene is not None:
            return Feature(label, CODING, "exact", gene=gene)

        bare = strip_version(label)
        if bare != label and self.annotation.gene(bare) is not None:
            return Feature(label, CODING, "exact", gene=self.annotation.gene(bare),
                           note="matched once the transcript version was stripped")
        label = bare

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

        tried = ", ".join(sorted(self.prefixes)) or "none"
        return _unresolved(
            label, f"no locus tag or symbol in the supplied annotation matches. If this names "
                   f"a non-coding region, its prefix is not one this run recognises (tried: "
                   f"{tried})")

    # -- intergenic ----------------------------------------------------------

    def _intergenic(self, label: str, body: str | None = None) -> Feature:
        body = label if body is None else body
        interval = self.annotation.intergenic(body)
        if interval is not None:
            return Feature(label, INTERGENIC, "exact", interval=interval)

        left, _, right = body.partition("-")
        left, right = strip_version(left), strip_version(right)
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
               "exact": 0, "strand_suffix": 0, "symbol": 0, "flanking": 0,
               UPSTREAM: 0, DOWNSTREAM: 0, "ambiguous": 0}
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
