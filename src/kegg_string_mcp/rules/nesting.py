"""Which rules contain which, across the whole population.

A learning classifier system produces rules that nest: a two-condition rule sits
inside a seven-condition one that adds five more. Read as 1,594 independent
findings, a population like that overstates what it found by a large factor. Read
as a few minimal rules plus their elaborations, it becomes something a person can
work through.

**Nesting is on (locus, state) pairs, not on loci.** `katG=1` is contained in
`katG=1 AND inhA=1`. `katG=0 AND inhA=1` is not, because the states conflict --
those two rules describe disjoint sets of isolates, which is the opposite of one
containing the other.

**A superset predicting the opposite class is the interesting case.** If `A=1`
predicts resistance and `A=1 AND B=1` predicts susceptibility, then B reverses
the outcome in A's presence. That is epistasis in the strict sense -- the effect
of one condition depending on another -- and it is invisible to any per-rule
reading. It is reported separately from ordinary nesting because it means
something different.

Costs one pass. For each rule, only the rules containing its first condition can
possibly contain it, so the candidate set is an index lookup rather than a scan
over the population.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field
from typing import Any

from kegg_string_mcp.rules.parse import Rule


@dataclass
class Nesting:
    """One rule's place in the population."""

    rule_id: str
    supersets: list[str] = field(default_factory=list)      # rules that contain this one
    subsets: list[str] = field(default_factory=list)        # rules this one contains
    contradicted_by: list[str] = field(default_factory=list)  # supersets predicting otherwise

    @property
    def is_minimal(self) -> bool:
        """No rule in the population is strictly contained in this one."""
        return not self.subsets

    @property
    def is_maximal(self) -> bool:
        return not self.supersets

    def to_dict(self) -> dict[str, Any]:
        return {"n_supersets": len(self.supersets), "n_subsets": len(self.subsets),
                "is_minimal": str(self.is_minimal), "is_maximal": str(self.is_maximal),
                "contained_in": "|".join(self.supersets[:5]) or "NA",
                "contradicted_by": "|".join(self.contradicted_by[:5]) or "NA"}


def _terms(rule: Rule, rename: dict[str, str] | None) -> frozenset[tuple[str, int]]:
    rename = rename or {}
    return frozenset((rename.get(c.label, c.label), c.state) for c in rule.conditions)


def nest(rules: list[Rule], rename: dict[str, str] | None = None) -> dict[str, Nesting]:
    """Containment for every rule, keyed by rule id.

    Two rules with identical condition sets contain each other, which would make
    both non-minimal and hide them both. They are treated as one: identical sets
    are not counted as containing one another.
    """
    keyed = [(rule.rule_id(rename), _terms(rule, rename), rule.predicted_class)
             for rule in rules]
    holding: dict[tuple[str, int], set[int]] = defaultdict(set)
    for index, (_, terms, _class) in enumerate(keyed):
        for term in terms:
            holding[term].add(index)

    out = {rule_id: Nesting(rule_id=rule_id) for rule_id, _, _ in keyed}
    for index, (rule_id, terms, predicted) in enumerate(keyed):
        if not terms:
            continue
        # Only rules holding every one of this rule's conditions can contain it.
        candidates = set.intersection(*(holding[term] for term in terms))
        for other in candidates:
            other_id, other_terms, other_class = keyed[other]
            if other == index or other_terms == terms or other_id == rule_id:
                continue
            out[rule_id].supersets.append(other_id)
            out[other_id].subsets.append(rule_id)
            if other_class != predicted:
                out[rule_id].contradicted_by.append(other_id)
    for nesting in out.values():
        nesting.supersets = sorted(dict.fromkeys(nesting.supersets))
        nesting.subsets = sorted(dict.fromkeys(nesting.subsets))
        nesting.contradicted_by = sorted(dict.fromkeys(nesting.contradicted_by))
    return out


def summarise(nestings: dict[str, Nesting]) -> dict[str, Any]:
    total = len(nestings)
    minimal = sum(1 for n in nestings.values() if n.is_minimal)
    contradicted = sum(1 for n in nestings.values() if n.contradicted_by)
    return {"rules": total,
            "minimal": minimal,
            "elaborations": total - minimal,
            "contradicted_by_a_superset": contradicted}
