"""Over-representation of annotation terms among a rule's loci.

Fisher's exact test, one-tailed, with Benjamini-Hochberg correction across the
terms tested for one rule -- the method the *M. tuberculosis* literature uses for
gene-set enrichment over TubercuList functional categories and COG classes.

**The annotation is both the map and the background.** Term membership and the
universe drawn from come from one hashed file, so they cannot disagree about
which loci exist. A KEGG-based enrichment cannot promise that: membership is
fetched, covers 29% of H37Rv, and the background then has to be chosen --
4,008 loci or the 1,172 annotated ones, a 40-fold difference in p for the same
observation. Here the universe is exactly the loci the annotation carries, which
is exactly the set the feature matrix was drawn from.

**Which axis matters more than the test.** Measured on 1,500 random k=7 draws,
the best p across terms falls under 0.05 in 52.4% of draws for KEGG pathways and
21.7% for functional categories -- before any real signal. That is the max
statistic, not the background, and it is why the number of terms tested is
reported beside every result and why BH is applied rather than a raw p being
quoted.

**A p here measures surprise, not significance.** The null is a uniform draw from
the annotated genome. The loci in a rule were not drawn uniformly: a classifier
selected them together because they jointly predict a phenotype. So a small p
says the set is unlike a random draw, which it certainly is -- not that the loci
are functionally related beyond that.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field
from math import comb
from typing import Any

from kegg_string_mcp.rules.annotation import Annotation

DEFAULT_AXIS = "category"
MIN_MEMBERS = 2


@dataclass
class Universe:
    """term -> loci carrying it, plus the loci that carry any term at all.

    `size` is the annotated universe, not the whole annotation. A locus with no
    term on this axis could not have been a hit, so including it in the
    denominator would treat "unknown" as "not in the term" -- the same absent
    versus unknown distinction the rest of this layer turns on.
    """

    axis: str
    members: dict[str, frozenset[str]] = field(default_factory=dict)
    annotated: frozenset[str] = frozenset()

    @property
    def size(self) -> int:
        return len(self.annotated)

    @property
    def terms(self) -> int:
        return len(self.members)


def universe_for(annotation: Annotation, axis: str = DEFAULT_AXIS) -> Universe:
    by_term: dict[str, set[str]] = defaultdict(set)
    annotated: set[str] = set()
    for gene in annotation.genes:
        for term in gene.terms.get(axis, ()):
            by_term[term].add(gene.locus)
            annotated.add(gene.locus)
    return Universe(axis=axis,
                    members={t: frozenset(v) for t, v in by_term.items()},
                    annotated=frozenset(annotated))


def p_at_least(m: int, k: int, s: int, n: int) -> float:
    """Fisher's exact, one-tailed: P(X >= m) drawing k of n, s of which carry the
    term. Computed exactly rather than approximated -- k is small and the counts
    at these sizes make a normal approximation unreliable."""
    if m <= 0 or k <= 0 or s <= 0 or n <= 0 or k > n:
        return 1.0
    top = min(k, s)
    if m > top:
        return 0.0
    total = comb(n, k)
    return sum(comb(s, i) * comb(n - s, k - i) for i in range(m, top + 1)) / total


def benjamini_hochberg(pvalues: list[float]) -> list[float]:
    """Adjusted values, monotone, in the order given.

    Applied across the terms tested for ONE rule. Correcting across every rule in
    a population instead would be a different and much harsher question, and is
    left to the caller who knows how many rules they intend to read.
    """
    if not pvalues:
        return []
    order = sorted(range(len(pvalues)), key=lambda i: pvalues[i])
    total = len(pvalues)
    adjusted = [1.0] * total
    running = 1.0
    for rank, index in enumerate(reversed(order), start=1):
        position = total - rank + 1
        running = min(running, pvalues[index] * total / position)
        adjusted[index] = min(1.0, running)
    return adjusted


@dataclass
class Enriched:
    term: str
    members: tuple[str, ...]        # loci in the set carrying it
    m: int                          # how many of the set
    k: int                          # set size, restricted to annotated loci
    s: int                          # loci in the universe carrying it
    n: int                          # universe size
    expected: float
    p: float
    q: float = 1.0                  # Benjamini-Hochberg adjusted

    @property
    def ratio(self) -> float:
        return self.m / self.expected if self.expected else 0.0

    def render(self) -> str:
        return (f"{self.term}:{self.m}/{self.k}:exp{self.expected:.2f}:q{self.q:.2g}")

    def to_dict(self) -> dict[str, Any]:
        return {"term": self.term, "members": "|".join(self.members), "m": self.m, "k": self.k,
                "s": self.s, "n": self.n, "expected": f"{self.expected:.3f}",
                "ratio": f"{self.ratio:.2f}", "p": f"{self.p:.3g}", "q": f"{self.q:.3g}"}


@dataclass
class Result:
    axis: str
    terms: list[Enriched] = field(default_factory=list)
    tested: int = 0                 # how many terms the p-values came from
    k: int = 0                      # loci in the set carrying any term
    dropped: tuple[str, ...] = ()   # loci with no term on this axis

    def significant(self, q: float = 0.05) -> list[Enriched]:
        return [t for t in self.terms if t.q <= q]

    def note(self) -> str:
        if not self.k:
            return (f"no locus in this rule carries a {self.axis} term, so no enrichment was "
                    f"computed. Not a negative result.")
        part = (f"; {len(self.dropped)} locus/loci carry no {self.axis} term and are excluded "
                f"from k" if self.dropped else "")
        return (f"{self.tested} {self.axis} term(s) tested over {self.k} annotated locus/loci"
                f"{part}. A p here measures surprise against a uniform draw from the annotated "
                f"genome, not evidence that the loci are related.")


def enrich(loci: list[str], universe: Universe, min_members: int = MIN_MEMBERS) -> Result:
    """Terms carried by at least `min_members` of the set, BH-adjusted.

    `min_members` defaults to 2 because a term carried by one locus says nothing
    about the set, and because strict intersection over all of them almost never
    fires: at the median rule size of seven loci, no term is usually shared by
    every one.
    """
    present = [locus for locus in dict.fromkeys(loci) if locus in universe.annotated]
    dropped = tuple(locus for locus in dict.fromkeys(loci) if locus not in universe.annotated)
    result = Result(axis=universe.axis, k=len(present), dropped=dropped)
    if not present:
        return result

    members = set(present)
    found: list[Enriched] = []
    for term, carriers in universe.members.items():
        hit = members & carriers
        if len(hit) < min_members:
            continue
        found.append(Enriched(
            term=term, members=tuple(sorted(hit)), m=len(hit), k=len(present),
            s=len(carriers), n=universe.size,
            expected=len(present) * len(carriers) / universe.size,
            p=p_at_least(len(hit), len(present), len(carriers), universe.size)))

    result.tested = len(found)
    for item, q in zip(found, benjamini_hochberg([f.p for f in found])):
        item.q = q
    result.terms = sorted(found, key=lambda f: (f.q, f.p, -f.m))
    return result
