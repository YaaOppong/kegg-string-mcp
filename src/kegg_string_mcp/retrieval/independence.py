"""Classify what the literature adds over the structured sources.

Retrieving abstracts about a pair STRING already scores is not new evidence -- and
where STRING's score is textmining, those abstracts may be the very papers that
produced it. The repo already refuses to double-count that at the tool level; this
applies the same rule one level up, to the retrieval arms.

Four cases where the literature is doing real work rather than restating STRING:

* **post-release** -- STRING is a fixed release, so a paper published after it
  cannot be in any channel, textmining included;
* **below threshold** -- a pair under the score cutoff is not returned at all,
  which is not the same as no relationship;
* **not an interaction** -- STRING models association between proteins.
  Compensatory mutation and co-occurrence in clinical isolates are genotype-level
  population phenomena; katG/ahpC is exactly this;
* **poorly covered genes** -- few edges at any threshold.

**Silence has to be earned.** Three things used to arrive here as "STRING says
nothing": a gene STRING never resolved, an edge that sat below the rank cut of
both partner lists, and a pair STRING genuinely scores below the threshold. Only
the third is a statement about the genes. The first two are now `unresolved` and
`truncated`, so a lookup failure cannot be counted as evidence of independence --
the same rule `retrieval/coverage.py` applies to annotation.

This also removes a circularity from the comparison. Scoring relevance by whether
a passage names the genes is correlated with BM25's own ranking function, which
tilts the measurement toward the lexical arm. STRING's assertion is independent of
both retrievers, so restricting to pairs STRING is silent on tests the arms where
neither has a built-in advantage.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from itertools import combinations
from typing import Any

# STRING's own medium-confidence band, matching string_db.MEDIUM_CONFIDENCE.
INDEPENDENT_MIN = 0.4
# STRING v12.0. A paper published after this cannot be in any channel.
STRING_RELEASE_YEAR = 2023


@dataclass
class PairVerdict:
    """What STRING says about a pair, and therefore what literature would add."""

    gene_a: str
    gene_b: str
    combined: float = 0.0
    textmining: float = 0.0
    max_non_textmining: float = 0.0
    # silent | textmining_only | corroborating | unresolved | truncated.
    # The last two are not verdicts about the genes: they say the question was
    # never put to STRING, which is a different fact from a negative answer and
    # must not reach the residue as one.
    status: str = "silent"
    note: str = ""
    string_ids: list[str] = field(default_factory=list)

    @property
    def undetermined(self) -> bool:
        return self.status in UNDETERMINED

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


UNDETERMINED = ("unresolved", "truncated")


@dataclass
class Edge:
    combined: float = 0.0
    textmining: float = 0.0
    max_non_textmining: float = 0.0


@dataclass
class StringEdges:
    """STRING's answer for a whole gene set, plus what it could not answer.

    Keyed by unordered pairs of STRING IDs rather than by name: a gene called
    `Rv0678` in one artefact and `mmpR5` in another is one protein, and keying by
    whichever spelling arrived first is what made the same edge read as 0.989 or
    as silent depending on the caller (see `kegg_string_mcp.identity`).
    """

    edges: dict[tuple[str, str], Edge] = field(default_factory=dict)
    identities: Any = None                                   # IdentitySet, or None
    required_score: int = 150
    source: str = "network"                                  # network | partners
    truncated: list[str] = field(default_factory=list)       # partner lists that hit the limit
    notes: list[str] = field(default_factory=list)

    def string_id(self, gene: str) -> str:
        if self.identities is None:
            return ""
        return self.identities.string_id(gene)

    def lookup(self, a: str, b: str) -> tuple[Edge | None, str, list[str]]:
        """(edge, state, ids). `state` is 'ok', 'unresolved' or 'truncated'."""
        id_a, id_b = self.string_id(a), self.string_id(b)
        ids = [i for i in (id_a, id_b) if i]
        missing = [g for g, i in ((a, id_a), (b, id_b)) if not i]
        if missing:
            return None, "unresolved", ids
        edge = self.edges.get(_pair(id_a, id_b))
        if edge is None and self.source == "partners" and {a, b} <= set(self.truncated):
            # Both partner lists were full, so an edge could sit past the cut of
            # both. /network answers this outright; the partner path cannot.
            return None, "truncated", ids
        return edge, "ok", ids


def _pair(a: str, b: str) -> tuple[str, str]:
    return tuple(sorted((a, b)))  # type: ignore[return-value]


def network_edges(genes: list[str], string_client: Any, identities: Any = None,
                  species: int = 83332, required_score: int = 150) -> StringEdges:
    """One /network call for the whole set: every pair answered, none truncated.

    Preferred over `partner_edges`. A low threshold is deliberate, so that "below
    threshold" stays distinguishable from "STRING holds nothing".
    """
    if identities is None:
        from kegg_string_mcp import identity as identity_module

        identities = identity_module.resolve(genes, string=string_client, species=species)

    ids = [i for i in (identities.string_id(g) for g in genes) if i]
    out = StringEdges(identities=identities, required_score=required_score, source="network")
    unresolved = [g for g in genes if not identities.string_id(g)]
    if unresolved:
        out.notes.append(f"STRING did not resolve {', '.join(unresolved)}; their pairs are "
                         f"reported as unresolved, not as silent.")
    if len(ids) < 2:
        out.notes.append("Fewer than two genes resolved, so no pair could be queried.")
        return out

    result = string_client.network(ids, species=species, required_score=required_score)
    out.notes.extend(result.notes)
    for record in result.records:
        detail = record.detail
        key = _pair(detail.get("string_id_a", ""), detail.get("string_id_b", ""))
        out.edges[key] = Edge(combined=detail.get("combined_score", 0.0),
                              textmining=detail.get("textmining_score", 0.0),
                              max_non_textmining=detail.get("max_non_textmining_score", 0.0))
    return out


def partner_edges(genes: list[str], string_client: Any, identities: Any = None,
                  species: int = 83332, required_score: int = 150,
                  limit: int = 200) -> StringEdges:
    """Fallback: one partner call per gene, keyed by STRING ID.

    Kept because it also measures degree, and because it works when /network is
    unavailable. It cannot see an edge below the rank cut of both genes -- at the
    module defaults every M. tuberculosis gene tested returned exactly `limit`
    partners -- so pairs whose lists are both full are reported as `truncated`
    rather than as `silent`.
    """
    if identities is None:
        from kegg_string_mcp import identity as identity_module

        identities = identity_module.resolve(genes, string=string_client, species=species)

    out = StringEdges(identities=identities, required_score=required_score, source="partners")
    for gene in genes:
        gene_id = identities.string_id(gene)
        if not gene_id:
            continue
        result = string_client.partners(gene, species=species, limit=limit,
                                        required_score=required_score)
        if len(result.records) >= limit:
            out.truncated.append(gene)
        for record in result.records:
            detail = record.detail
            out.edges[_pair(gene_id, record.record_id)] = Edge(
                combined=detail.get("combined_score", 0.0),
                textmining=detail.get("textmining_score", 0.0),
                max_non_textmining=detail.get("max_non_textmining_score", 0.0))
    if out.truncated:
        out.notes.append(f"Partner lists hit the limit of {limit} for "
                         f"{', '.join(out.truncated)}; pairs between two of them are reported "
                         f"as truncated, because an edge can sit below both cuts.")
    return out


def classify(genes: list[str], edges: StringEdges) -> list[PairVerdict]:
    verdicts: list[PairVerdict] = []
    for a, b in combinations(genes, 2):
        edge, state, ids = edges.lookup(a, b)
        verdict = PairVerdict(gene_a=a, gene_b=b, string_ids=ids)
        if state == "unresolved":
            # The pair was never asked about. Calling it silent would hand the
            # residue a resolution failure dressed as a negative result.
            verdict.status = "unresolved"
            verdict.note = ("at least one gene did not resolve to a STRING protein, so this pair "
                            "was never queried -- this is a resolution failure, not evidence of "
                            "no interaction")
        elif state == "truncated":
            verdict.status = "truncated"
            verdict.note = ("both partner lists hit the retrieval limit, so an edge between them "
                            "could sit below both cuts and was not checked")
        elif edge is None:
            verdict.status = "silent"
            verdict.note = (f"STRING returns no edge at or above {edges.required_score}, so "
                            f"literature is the only evidence and is not restating a "
                            f"structured source")
        else:
            verdict.combined = edge.combined
            verdict.textmining = edge.textmining
            verdict.max_non_textmining = edge.max_non_textmining
            if edge.max_non_textmining >= INDEPENDENT_MIN:
                verdict.status = "corroborating"
                verdict.note = ("STRING already asserts this pair on a non-textmining channel; "
                                "literature corroborates rather than adds")
            else:
                verdict.status = "textmining_only"
                verdict.note = ("STRING's score for this pair is essentially textmining, so the "
                                "retrieved abstracts may be the very papers that produced it -- "
                                "counting both is one line of evidence twice")
        verdicts.append(verdict)
    return verdicts


@dataclass
class IndependenceReport:
    verdicts: list[PairVerdict] = field(default_factory=list)
    post_release_papers: int = 0
    total_papers: int = 0
    # Genes no source could resolve, after the cross-source retries. A count of
    # `unresolved` verdicts buried in `counts` is easy to read past; the gene
    # names are what a reader can act on, so they are a field of their own.
    unresolved: list[str] = field(default_factory=list)

    def by_status(self) -> dict[str, int]:
        counts: dict[str, int] = {}
        for v in self.verdicts:
            counts[v.status] = counts.get(v.status, 0) + 1
        return counts

    def pairs_with_status(self, status: str) -> list[tuple[str, str]]:
        return [(v.gene_a, v.gene_b) for v in self.verdicts if v.status == status]

    def to_dict(self) -> dict[str, Any]:
        return {"counts": self.by_status(),
                "unresolved": self.unresolved,
                "post_release_papers": self.post_release_papers,
                "total_papers": self.total_papers,
                "string_release_year": STRING_RELEASE_YEAR,
                "verdicts": [v.to_dict() for v in self.verdicts]}


def post_release_fraction(corpus: Any) -> tuple[int, int]:
    """Papers that cannot be in any STRING channel because they postdate the release."""
    years = {p.pmid: p.year for p in corpus.passages}
    recent = sum(1 for y in years.values()
                 if y.isdigit() and int(y) >= STRING_RELEASE_YEAR)
    return recent, len(years)
