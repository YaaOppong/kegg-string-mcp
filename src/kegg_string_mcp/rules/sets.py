"""What is true of a rule's loci as a set, rather than of its pairs.

The median rule in a real population has seven conditions, and 96.6% have three
or more. At that size the pairwise reading stops working: a k=14 rule decomposes
to 91 pairs, and a verdict that lists them says nothing a reader can hold. Worse,
some facts are only visible at set level -- three loci can pairwise share three
different pathways and share none in common.

So this is a second layer, not a replacement. Pairwise claims stay pairwise:
compensation is a claim about *this anchor* and *its* compensator, aliasing about
*these two* conditions. Set-level claims come from here.

**The base rate falls with k, and the fixed threshold does not.** A pathway
holding a sixth of the genome is a base rate for two loci and a 1-in-200,000
coincidence for seven. The enrichment test carries k, so a container category
stops being dismissed the moment enough loci share it.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from kegg_string_mcp.rules.annotation import Annotation
from kegg_string_mcp.rules.enrichment import DEFAULT_AXIS, Result, Universe, enrich

# Loci closer than this are reported as a contiguous run. Not an operon call --
# that needs expression evidence -- but a run of neighbours in a rule is worth
# seeing, because co-located genes are frequently co-transcribed.
RUN_GAP = 200

CLIQUE, STAR, CHAIN, SPARSE, DISCONNECTED = "clique", "star", "chain", "sparse", "disconnected"


@dataclass
class Subgraph:
    """The shape of the interaction graph among a rule's loci.

    A star and a clique are different findings: a star means one hub the others
    each touch, a clique means everything touches everything. Reporting only the
    edge count loses that, and at k=7 the edge count alone is unreadable anyway.
    """

    n: int = 0
    edges: int = 0
    shape: str = DISCONNECTED
    hub: str = ""
    components: int = 0

    @property
    def possible(self) -> int:
        return self.n * (self.n - 1) // 2

    @property
    def density(self) -> float:
        return self.edges / self.possible if self.possible else 0.0


def _shape(loci: list[str], pairs: set[tuple[str, str]]) -> Subgraph:
    n = len(loci)
    out = Subgraph(n=n, edges=len(pairs))
    if n < 2:
        out.shape = DISCONNECTED
        out.components = n
        return out

    degree: dict[str, int] = dict.fromkeys(loci, 0)
    adjacency: dict[str, set[str]] = {locus: set() for locus in loci}
    for left, right in pairs:
        if left in degree and right in degree:
            degree[left] += 1
            degree[right] += 1
            adjacency[left].add(right)
            adjacency[right].add(left)

    seen: set[str] = set()
    for locus in loci:
        if locus in seen:
            continue
        out.components += 1
        stack = [locus]
        while stack:
            current = stack.pop()
            if current in seen:
                continue
            seen.add(current)
            stack.extend(adjacency[current] - seen)

    if not pairs:
        out.shape = DISCONNECTED
    elif out.edges == out.possible:
        out.shape = CLIQUE
    else:
        top = max(degree, key=lambda locus: degree[locus])
        # A star: one locus touches everything and the edges are only its own.
        if degree[top] == n - 1 and out.edges == n - 1:
            out.shape, out.hub = STAR, top
        elif (out.components == 1 and out.edges == n - 1
              and max(degree.values()) <= 2):
            # A tree with n-1 edges is only a chain if no vertex branches. A
            # T-shape has the same edge count and is not one.
            out.shape = CHAIN
        else:
            out.shape = SPARSE
    return out


@dataclass
class SetEvidence:
    loci: list[str] = field(default_factory=list)
    enrichment: Result | None = None
    shared_by_all: list[str] = field(default_factory=list)
    common_partners: list[str] = field(default_factory=list)
    # Loci with no partner list at all. A failed STRING lookup must not read as
    # "these loci share no partner"; it means the question was not answered.
    partners_unknown: list[str] = field(default_factory=list)
    subgraph: Subgraph = field(default_factory=Subgraph)
    runs: list[list[str]] = field(default_factory=list)
    drugs: list[str] = field(default_factory=list)

    @property
    def k(self) -> int:
        return len(self.loci)

    def to_dict(self) -> dict[str, Any]:
        top = self.enrichment.terms[0] if self.enrichment and self.enrichment.terms else None
        return {
            "set_k": self.k,
            "enriched_term": top.term if top else "NA",
            "enriched_m_of_k": f"{top.m}/{top.k}" if top else "NA",
            "enriched_expected": f"{top.expected:.2f}" if top else "NA",
            "enriched_q": f"{top.q:.3g}" if top else "NA",
            "terms_tested": self.enrichment.tested if self.enrichment else 0,
            "shared_by_all": "|".join(self.shared_by_all) or "NA",
            "common_partners": "|".join(self.common_partners[:10]) or "NA",
            "partners_unknown": "|".join(self.partners_unknown) or "NA",
            "subgraph_shape": self.subgraph.shape,
            "subgraph_edges": f"{self.subgraph.edges}/{self.subgraph.possible}",
            "subgraph_hub": self.subgraph.hub or "NA",
            "contiguous_runs": "|".join("+".join(r) for r in self.runs) or "NA",
            "n_drugs": len(self.drugs),
            "drugs": "|".join(self.drugs) or "NA",
        }


def contiguous_runs(loci: list[str], annotation: Annotation) -> list[list[str]]:
    """Runs of loci that sit next to each other on the genome.

    Reported as adjacency, never as an operon: that needs strand agreement and
    expression evidence this does not have.
    """
    genes = [g for g in (annotation.gene(locus) for locus in loci) if g is not None]
    genes.sort(key=lambda g: g.start)
    runs: list[list[str]] = []
    current: list[Any] = []
    for gene in genes:
        if current and gene.start - current[-1].end - 1 <= RUN_GAP:
            current.append(gene)
        else:
            if len(current) > 1:
                runs.append([g.locus for g in current])
            current = [gene]
    if len(current) > 1:
        runs.append([g.locus for g in current])
    return runs


def set_evidence(conditions: list[Any], annotation: Annotation,
                 universe: Universe | None = None,
                 partners: dict[str, set[str]] | None = None,
                 linked_pairs: set[tuple[str, str]] | None = None,
                 axis: str = DEFAULT_AXIS) -> SetEvidence:
    """Everything true of the set the rule's present conditions name.

    Only conditions asserting state 1 count. A `=0` condition says the locus
    matches the reference, so including it in a shared-pathway or common-partner
    statement would attribute to the set a locus the rule says is not varying.
    """
    loci = [c.locus for c in conditions
            if c.condition.state == 1 and c.feature.resolved]
    out = SetEvidence(loci=sorted(dict.fromkeys(loci)))
    if not out.loci:
        return out

    if universe is not None:
        out.enrichment = enrich(out.loci, universe)
        terms = [set(annotation.gene(locus).terms.get(axis, ()))
                 for locus in out.loci if annotation.gene(locus) is not None]
        if terms and all(terms):
            out.shared_by_all = sorted(set.intersection(*terms))

    if partners:
        out.partners_unknown = [locus for locus in out.loci if not partners.get(locus)]
        if not out.partners_unknown:
            out.common_partners = sorted(
                set.intersection(*(partners[locus] for locus in out.loci)))

    pairs = {tuple(sorted(p)) for p in (linked_pairs or set())
             if set(p) <= set(out.loci)}
    out.subgraph = _shape(out.loci, pairs)          # type: ignore[arg-type]
    out.runs = contiguous_runs(out.loci, annotation)
    out.drugs = sorted({d for c in conditions if c.condition.state == 1
                        for d in c.catalogue.drugs})
    return out


def describe(evidence: SetEvidence) -> str:
    """A sentence about the set, for a verdict that would otherwise list pairs."""
    if evidence.k < 2:
        return ""
    parts = [f"The rule names {evidence.k} loci present"]
    if evidence.drugs:
        parts.append(f"spanning {len(evidence.drugs)} drug(s) ({', '.join(evidence.drugs)})")

    top = evidence.enrichment.terms[0] if evidence.enrichment and evidence.enrichment.terms \
        else None
    if top is not None and top.q <= 0.05:
        parts.append(
            f"and over-represented for {top.term} ({top.m} of {top.k}, expected "
            f"{top.expected:.2f}, q={top.q:.2g} over {evidence.enrichment.tested} term(s) "
            f"tested)")
    elif evidence.enrichment is not None and evidence.enrichment.k:
        # The test count belongs here MORE than beside a hit: a null result is
        # only interpretable if you know how many chances it had to fire.
        parts.append(f"with no term over-represented beyond its base rate "
                     f"({evidence.enrichment.tested} term(s) tested)")

    sentence = ", ".join(parts) + "."
    extra: list[str] = []
    if evidence.shared_by_all:
        extra.append(f"Every locus carries {', '.join(evidence.shared_by_all)}.")
    if evidence.common_partners:
        extra.append(
            f"All of them partner with {', '.join(evidence.common_partners[:3])}"
            + (f" and {len(evidence.common_partners) - 3} other(s)"
               if len(evidence.common_partners) > 3 else "")
            + " -- a shared neighbour of every member is the shape a complex or a regulon has, "
              "though a partner list capped at 20 cannot distinguish that from a hub.")
    if evidence.partners_unknown:
        extra.append(
            f"No partner list was retrieved for {', '.join(evidence.partners_unknown)}, so "
            f"whether the set shares a neighbour is undetermined rather than answered.")
    if evidence.subgraph.shape == STAR:
        extra.append(f"Their interaction graph is a star around {evidence.subgraph.hub}: the "
                     f"others touch it and not each other.")
    elif evidence.subgraph.shape == CLIQUE and evidence.k > 2:
        extra.append("Every pair among them interacts.")
    elif evidence.subgraph.shape == DISCONNECTED and evidence.k > 2:
        extra.append("No interaction is recorded between any of them.")
    if evidence.runs:
        extra.append(
            "Genomically contiguous: "
            + "; ".join("+".join(run) for run in evidence.runs)
            + " -- adjacency, not an operon call, which needs expression evidence.")
    return " ".join([sentence, *extra])
