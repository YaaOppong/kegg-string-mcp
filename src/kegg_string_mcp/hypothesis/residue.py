"""Which pairs are still unexplained after the structured sources have spoken.

This is the interface between what the repo already does and hypothesis
generation. Stages 1 and 2 explain what they can; whatever survives is the input
to stage 3, and is worth an LLM's attention precisely because nothing cheaper
accounted for it.

Every pair carries the reasons it was explained, and every reason is recorded
whether or not it counts. Nothing is silently dropped: the residue is a filter
you can re-run with a different definition of "explained" without refetching
anything, and a pair excluded in one configuration can be inspected to see why.

**On using STRING's textmining channel as a novelty filter.** It is tempting --
it is computed over all of PubMed rather than over one corpus of a few hundred
abstracts, so it should catch prior art that local co-mention counting misses. It
does not survive contact with this domain. STRING's textmining channel scores
co-occurrence in an abstract, not a reported relationship, and M. tuberculosis
resistance genes co-occur constantly in review tables listing genes that confer
resistance. Measured on the 41-gene set: katG-pncA scores 0.965 textmining with
every other channel below 0.05, as do katG-gyrA, embB-rpsL and 40-odd more pairs
whose only published connection is appearing in the same list. Treating that as
"already known" would discard candidates on the strength of a shared table row,
so `string_textmining` is recorded and off by default.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any

# Reasons a pair might already be accounted for. Recorded independently of
# whether they are configured to count, so the residue can be recomputed.
STRING_EXPERIMENTAL = "string_experimental"   # non-textmining channel, i.e. real evidence
STRING_TEXTMINING = "string_textmining"       # co-occurrence; see module docstring
SHARED_PATHWAY = "shared_pathway"             # both in one KEGG pathway
CO_MENTIONED = "co_mentioned"                 # a corpus paper names both

# Shared pathway is deliberately absent. A compensatory relationship between two
# genes in one pathway is a target of this work, not an explanation of it -- the
# question is whether the pair jointly confers resistance, and "both are in
# tuberculosis drug resistance (mtu01501)" does not answer it.
DEFAULT_EXPLAINING = frozenset({STRING_EXPERIMENTAL, CO_MENTIONED})


@dataclass
class Reason:
    code: str
    detail: str
    value: float = 0.0

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass
class PairAssessment:
    gene_a: str
    gene_b: str
    reasons: list[Reason] = field(default_factory=list)

    def codes(self) -> set[str]:
        return {r.code for r in self.reasons}

    def is_residue(self, explaining: frozenset[str] = DEFAULT_EXPLAINING) -> bool:
        return not (self.codes() & explaining)

    def to_dict(self) -> dict[str, Any]:
        return {"gene_a": self.gene_a, "gene_b": self.gene_b,
                "reasons": [r.to_dict() for r in self.reasons]}


def assess(pairs: list[tuple[str, str]],
           string_status: dict[tuple[str, str], dict] | None = None,
           pathways: dict[str, set[str]] | None = None,
           co_mentions: dict[tuple[str, str], int] | None = None) -> list[PairAssessment]:
    """Attach every applicable reason to every pair. Pure -- no network.

    Inputs are keyed by unordered pair, so callers need not agree on gene order:
    `_key` sorts before lookup.
    """
    string_status = _normalise(string_status or {})
    co_mentions = _normalise(co_mentions or {})
    pathways = pathways or {}

    out: list[PairAssessment] = []
    for a, b in pairs:
        assessment = PairAssessment(gene_a=a, gene_b=b)
        key = _key(a, b)

        edge = string_status.get(key)
        if edge:
            status = edge.get("status")
            if status == "corroborating":
                assessment.reasons.append(Reason(
                    STRING_EXPERIMENTAL,
                    "STRING asserts this pair on a non-textmining channel",
                    edge.get("max_non_textmining", 0.0)))
            elif status == "textmining_only":
                assessment.reasons.append(Reason(
                    STRING_TEXTMINING,
                    "STRING links this pair by co-occurrence in text only, which in this "
                    "domain is often a shared review table rather than a reported relationship",
                    edge.get("textmining", 0.0)))

        shared = pathways.get(a, set()) & pathways.get(b, set())
        if shared:
            assessment.reasons.append(Reason(
                SHARED_PATHWAY,
                f"both genes are in KEGG {', '.join(sorted(shared))}",
                float(len(shared))))

        count = co_mentions.get(key, 0)
        if count:
            assessment.reasons.append(Reason(
                CO_MENTIONED, f"{count} corpus paper(s) name both genes", float(count)))

        out.append(assessment)
    return out


def residue(assessments: list[PairAssessment],
            explaining: frozenset[str] = DEFAULT_EXPLAINING) -> list[PairAssessment]:
    return [a for a in assessments if a.is_residue(explaining)]


def summarise(assessments: list[PairAssessment],
              explaining: frozenset[str] = DEFAULT_EXPLAINING) -> dict[str, Any]:
    # Counted over sorted codes, not raw set iteration: a set of strings iterates
    # in an order that depends on PYTHONHASHSEED, so the written JSON differed
    # between runs by key order alone. Identical values in a different order
    # still breaks diffing an artefact against itself.
    counts: dict[str, int] = {}
    for assessment in assessments:
        for code in sorted(assessment.codes()):
            counts[code] = counts.get(code, 0) + 1
    remaining = residue(assessments, explaining)
    return {"pairs": len(assessments),
            "reason_counts": dict(sorted(counts.items())),
            "explaining": sorted(explaining),
            "residue": len(remaining),
            "residue_fraction": round(len(remaining) / len(assessments), 3)
            if assessments else 0.0}


def _key(a: str, b: str) -> tuple[str, str]:
    return tuple(sorted((a.lower(), b.lower())))  # type: ignore[return-value]


def _normalise(mapping: dict) -> dict:
    return {_key(*k): v for k, v in mapping.items()}
