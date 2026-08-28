"""Scoring. Pure functions over a run's output -- no network, no model.

Three things are measured, and they are deliberately kept apart because they mean
different things and only two of them are ground truth:

1. **Pathway recall/precision** on positive controls -- retrieval fidelity
   against KEGG's own assignments.
2. **Abstention** on negative controls -- did the pipeline correctly report
   nothing where KEGG knows nothing, or did it invent something?
3. **Citation and quote precision** -- exact, deterministic, and computed by the
   validator rather than judged. This is the number that is genuinely hard to
   argue with, and the one most RAG evaluations do not have.
"""

from __future__ import annotations

import re
from dataclasses import asdict, dataclass, field
from typing import Any

from kegg_string_mcp.evaluate.gold import GoldGene

PATHWAY_ID = re.compile(r"\b[a-z]{3,4}\d{5}\b")


def reported_pathways(summary: str, organism: str,
                      validation: dict[str, Any] | None = None) -> set[str]:
    """Pathway IDs the summary claims FOR THIS GENE.

    Restricted to the organism asked about, and -- importantly -- excluding IDs the
    validator marked `cross_target`. A model annotating furA may legitimately look
    up its neighbour katG and discuss katG's pathways as context; counting those as
    claims about furA scored a correct annotation as a fabrication. The validator
    already distinguishes the two, so reuse its judgement rather than re-deriving
    it from prose.
    """
    found = {p for p in PATHWAY_ID.findall(summary or "") if p.startswith(organism)}
    for citation in (validation or {}).get("citations", []):
        if citation.get("status") == "cross_target":
            found.discard(citation.get("identifier"))
    return found


@dataclass
class GeneScore:
    gene: str
    kegg_gene_id: str
    negative_control: bool
    expected: list[str] = field(default_factory=list)
    reported: list[str] = field(default_factory=list)
    hits: list[str] = field(default_factory=list)
    missed: list[str] = field(default_factory=list)
    extra: list[str] = field(default_factory=list)
    abstained: bool | None = None      # negatives only
    citations_total: int = 0
    citations_verified: int = 0
    quotes_total: int = 0
    quotes_verified: int = 0
    validation_passed: bool = False
    error: str = ""

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def score_gene(gold: GoldGene, summary: str, validation: dict[str, Any],
               organism: str) -> GeneScore:
    expected = set(gold.pathways)
    reported = reported_pathways(summary, organism, validation)
    citations = validation.get("citations", [])
    quotes = validation.get("quotes", [])

    score = GeneScore(
        gene=gold.query,
        kegg_gene_id=gold.kegg_gene_id,
        negative_control=gold.is_negative_control,
        expected=sorted(expected),
        reported=sorted(reported),
        hits=sorted(expected & reported),
        missed=sorted(expected - reported),
        extra=sorted(reported - expected),
        citations_total=len(citations),
        citations_verified=sum(1 for c in citations if c["status"] == "verified"),
        quotes_total=len(quotes),
        quotes_verified=sum(1 for q in quotes if q["status"] == "verified"),
        validation_passed=bool(validation.get("passed")),
    )
    if gold.is_negative_control:
        # The whole point: KEGG assigns nothing, so reporting nothing is correct.
        score.abstained = not reported
    return score


@dataclass
class EvalReport:
    organism: str
    reference: str
    retrieved_on: str
    coverage: dict
    scores: list[GeneScore] = field(default_factory=list)

    # -- positive controls: retrieval fidelity ------------------------------

    @property
    def positives(self) -> list[GeneScore]:
        return [s for s in self.scores if not s.negative_control and not s.error]

    @property
    def recall(self) -> float | None:
        expected = sum(len(s.expected) for s in self.positives)
        return round(sum(len(s.hits) for s in self.positives) / expected, 3) if expected else None

    @property
    def precision(self) -> float | None:
        reported = sum(len(s.reported) for s in self.positives)
        return round(sum(len(s.hits) for s in self.positives) / reported, 3) if reported else None

    # -- negative controls: abstention --------------------------------------

    @property
    def negatives(self) -> list[GeneScore]:
        return [s for s in self.scores if s.negative_control and not s.error]

    @property
    def abstention_rate(self) -> float | None:
        if not self.negatives:
            return None
        return round(sum(1 for s in self.negatives if s.abstained) / len(self.negatives), 3)

    @property
    def fabricated_on(self) -> list[str]:
        return [s.gene for s in self.negatives if not s.abstained]

    # -- citation integrity: exact ------------------------------------------

    @property
    def citation_precision(self) -> float | None:
        total = sum(s.citations_total for s in self.scores)
        return round(sum(s.citations_verified for s in self.scores) / total, 3) if total else None

    @property
    def quote_precision(self) -> float | None:
        total = sum(s.quotes_total for s in self.scores)
        return round(sum(s.quotes_verified for s in self.scores) / total, 3) if total else None

    def to_dict(self) -> dict[str, Any]:
        return {
            "organism": self.organism,
            "reference": self.reference,
            "reference_retrieved_on": self.retrieved_on,
            "coverage": self.coverage,
            "retrieval_fidelity": {"recall": self.recall, "precision": self.precision,
                                   "n_genes": len(self.positives)},
            "abstention": {"rate": self.abstention_rate, "n_genes": len(self.negatives),
                           "fabricated_on": self.fabricated_on},
            "citation_integrity": {"citation_precision": self.citation_precision,
                                   "quote_precision": self.quote_precision},
            "errors": [{"gene": s.gene, "error": s.error} for s in self.scores if s.error],
            "per_gene": [s.to_dict() for s in self.scores],
        }
