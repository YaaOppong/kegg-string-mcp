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
    status: str = "silent"      # silent | textmining_only | corroborating
    note: str = ""

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def gene_partner_map(genes: list[str], string_client: Any, species: int = 83332,
                     required_score: int = 150, limit: int = 200) -> dict[str, dict[str, dict]]:
    """One STRING call per gene, not one per pair.

    820 pairs would be 820 calls against a service that asks for roughly one
    request a second. Partner lists are symmetric enough to derive every pair from
    41 calls, and a low score threshold is used deliberately so that "below
    threshold" stays distinguishable from "silent".
    """
    out: dict[str, dict[str, dict]] = {}
    for gene in genes:
        result = string_client.partners(gene, species=species, limit=limit,
                                        required_score=required_score)
        by_name: dict[str, dict] = {}
        for record in result.records:
            detail = record.detail
            by_name[record.name.lower()] = {
                "combined": detail.get("combined_score", 0.0),
                "textmining": detail.get("textmining_score", 0.0),
                "max_non_textmining": detail.get("max_non_textmining_score", 0.0),
            }
        out[gene] = by_name
    return out


def classify(genes: list[str], partners: dict[str, dict[str, dict]]) -> list[PairVerdict]:
    verdicts: list[PairVerdict] = []
    for a, b in combinations(genes, 2):
        # Either direction: STRING returns partner lists, not a symmetric matrix,
        # and a limit can truncate one side without truncating the other.
        edge = partners.get(a, {}).get(b.lower()) or partners.get(b, {}).get(a.lower())
        verdict = PairVerdict(gene_a=a, gene_b=b)
        if edge is None:
            verdict.status = "silent"
            verdict.note = ("STRING returns no edge at this threshold, so literature is the "
                            "only evidence and is not restating a structured source")
        else:
            verdict.combined = edge["combined"]
            verdict.textmining = edge["textmining"]
            verdict.max_non_textmining = edge["max_non_textmining"]
            if edge["max_non_textmining"] >= INDEPENDENT_MIN:
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

    def by_status(self) -> dict[str, int]:
        counts: dict[str, int] = {}
        for v in self.verdicts:
            counts[v.status] = counts.get(v.status, 0) + 1
        return counts

    def pairs_with_status(self, status: str) -> list[tuple[str, str]]:
        return [(v.gene_a, v.gene_b) for v in self.verdicts if v.status == status]

    def to_dict(self) -> dict[str, Any]:
        return {"counts": self.by_status(),
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
