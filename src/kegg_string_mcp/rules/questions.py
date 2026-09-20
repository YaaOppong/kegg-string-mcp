"""What the structured sources could not answer, phrased as questions.

Stage A of the literature layer, and it costs nothing: the question list falls
out of the classification, so it can be counted and priced before a single
retrieval happens.

Questions are about **loci and pairs, never about rules**. A classifier
population reuses its loci heavily, so the same gap is raised by many rules; ask
once, record which rules raised it, and the answer serves all of them. Asking per
rule would multiply the cost by the redundancy in the rule set and produce
answers that could not be reused.

Each question is narrow enough to be answered with a verbatim quote. "Explain
this rule" cannot be checked; "is there published evidence that variants in X
compensate for the fitness cost of Y resistance mutations" either has a passage
behind it or does not.

The wording matters for a second reason. A question that names a conclusion
invites the model to confirm it, so these ask for the evidence and leave the
verdict to the code that checks the quote.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from kegg_string_mcp.rules.signature import (
    ROLE_ANCHOR,
    ROLE_NEGATED_ANCHOR,
    ROLE_UNKNOWN,
    RuleSignature,
)

# What is being asked, in the order a reviewer would care about it.
COMPENSATION = "compensation"        # does this locus compensate for that one?
RESISTANCE_ROLE = "resistance_role"  # is this locus implicated in resistance at all?
ALT_ROUTE = "alt_route"              # what confers resistance without the canonical locus?
PRIOR_ART = "prior_art"              # is any relationship between these two reported?

ORDER = (COMPENSATION, ALT_ROUTE, RESISTANCE_ROLE, PRIOR_ART)

ORGANISM = "Mycobacterium tuberculosis"


@dataclass
class Question:
    kind: str
    locus: str
    partner: str = ""
    drugs: tuple[str, ...] = ()
    text: str = ""
    names: dict[str, str] = field(default_factory=dict)     # locus -> name in papers
    raised_by: list[str] = field(default_factory=list)      # rule ids
    raised_by_signature: set[str] = field(default_factory=set)

    @property
    def key(self) -> tuple[str, str, str]:
        return (self.kind, self.locus, self.partner)

    @property
    def id(self) -> str:
        parts = [self.kind, self.locus] + ([self.partner] if self.partner else [])
        return "-".join(parts + (["-".join(self.drugs)] if self.drugs else []))

    def to_dict(self) -> dict[str, Any]:
        return {"question_id": self.id, "kind": self.kind, "locus": self.locus,
                "partner": self.partner or "NA",
                "drugs": "|".join(self.drugs) or "NA",
                "n_rules": len(self.raised_by),
                "raised_by": "|".join(self.raised_by[:10]) or "NA",
                "signatures": "|".join(sorted(self.raised_by_signature)) or "NA",
                "question": self.text}


def _drugs_of(signature: RuleSignature) -> tuple[str, ...]:
    """The drugs this rule is about.

    Context for a question about an unexplained locus: "involved in resistance"
    is unanswerable, "involved in isoniazid resistance" is a search.

    Negated anchors count. `katG=0 AND <locus>=1 -> R` is still a rule about
    isoniazid: the canonical locus being at reference is what makes the question
    interesting, not a reason to drop the drug from it.
    """
    return tuple(sorted({d for c in signature.conditions
                         if c.role in (ROLE_ANCHOR, ROLE_NEGATED_ANCHOR)
                         for d in c.catalogue.drugs}))


def _names(signature: RuleSignature) -> dict[str, str]:
    """locus -> the name a paper would use for it.

    Questions are search strings before they are anything else, and the
    literature says "ahpC", not "Rv2428". Both are given -- the symbol to find
    the papers, the locus tag so the answer joins back to the tables.
    """
    out: dict[str, str] = {}
    for condition in signature.conditions:
        gene = condition.feature.gene
        symbol = (gene.symbol if gene is not None else "") or ""
        out[condition.locus] = (f"{symbol} ({condition.locus})" if symbol
                                else condition.locus)
    return out


def _drug_phrase(drugs: tuple[str, ...]) -> str:
    return " or ".join(drugs) if drugs else "drug"


def _render(question: Question) -> str:
    """The question text, built after every rule that raises it has been merged."""
    def named(locus: str) -> str:
        return question.names.get(locus, locus)

    drugs = _drug_phrase(question.drugs)
    if question.kind == COMPENSATION:
        return (f"Is there published evidence that variants in {named(question.locus)} "
                f"compensate for the fitness cost of {named(question.partner)} {drugs} "
                f"resistance mutations in {ORGANISM}? Quote any passage that states such a "
                f"relationship.")
    if question.kind == ALT_ROUTE:
        return (f"What mechanisms are reported to confer {drugs} resistance in {ORGANISM} "
                f"isolates that carry no {named(question.locus)} mutation? Quote any passage "
                f"naming such a mechanism.")
    if question.kind == PRIOR_ART:
        return (f"Is there published evidence of a functional relationship between "
                f"{named(question.locus)} and {named(question.partner)} in {ORGANISM} -- "
                f"shared pathway, regulation, physical interaction, or joint involvement in "
                f"resistance? Quote any passage that states one. Co-occurrence in a list of "
                f"resistance genes is not such a statement.")
    return (f"Is there published evidence that {named(question.locus)} is involved in {drugs} "
            f"resistance in {ORGANISM}, whether as a resistance determinant, a compensatory "
            f"locus, or a regulator? Quote any passage that states a role.")


def generate(signatures: list[RuleSignature]) -> list[Question]:
    """One question per distinct gap, with every rule that raised it attached."""
    seen: dict[tuple, Question] = {}

    def add(question: Question, signature: RuleSignature) -> None:
        existing = seen.get(question.key)
        if existing is None:
            seen[question.key] = question
            existing = question
        else:
            # Same gap, different rule context. Merge rather than ask twice.
            existing.drugs = tuple(sorted(set(existing.drugs) | set(question.drugs)))
            existing.names |= question.names
        rule_id = signature.rule.rule_id()
        if rule_id not in existing.raised_by:
            existing.raised_by.append(rule_id)
        existing.raised_by_signature.add(signature.primary)

    for signature in signatures:
        drugs = _drugs_of(signature)
        names = _names(signature)

        # A pair already linked by evidence beyond co-mention needs no literature
        # to be plausible, so it is not asked about.
        linked = {(a, b) for a, b, link in signature.links
                  if link.kind != "string_textmining_only"}
        for anchor, candidate, pair_drugs in signature.compensation_pairs:
            if (anchor, candidate) in linked:
                continue
            add(Question(kind=COMPENSATION, locus=candidate, partner=anchor,
                         drugs=tuple(pair_drugs), names=dict(names)), signature)

        for condition in signature.conditions:
            if condition.role == ROLE_UNKNOWN:
                add(Question(kind=RESISTANCE_ROLE, locus=condition.locus, drugs=drugs,
                             names=dict(names)), signature)
            elif condition.role == ROLE_NEGATED_ANCHOR:
                add(Question(kind=ALT_ROUTE, locus=condition.locus,
                             drugs=tuple(condition.catalogue.drugs),
                             names=dict(names)), signature)

        # Prior art for pairs the structured sources say nothing about. Cheap
        # co-mention is NOT the test: M. tuberculosis resistance genes co-occur
        # constantly in review tables -- katG/pncA scores 0.965 textmining with
        # every other STRING channel below 0.05 -- so the question asks for a
        # passage stating a relationship, which a shared table row cannot supply.
        unexplained = sorted({c.locus for c in signature.conditions
                              if c.role == ROLE_UNKNOWN and c.condition.state == 1})
        for index, left in enumerate(unexplained):
            for right in unexplained[index + 1:]:
                add(Question(kind=PRIOR_ART, locus=left, partner=right, drugs=drugs,
                             names=dict(names)), signature)

    for question in seen.values():
        question.text = _render(question)
    return sorted(seen.values(),
                  key=lambda q: (ORDER.index(q.kind), -len(q.raised_by), q.locus))


def summarise(questions: list[Question]) -> dict[str, Any]:
    counts: dict[str, int] = {}
    for question in questions:
        counts[question.kind] = counts.get(question.kind, 0) + 1
    return {"questions": len(questions),
            "by_kind": {kind: counts.get(kind, 0) for kind in ORDER if counts.get(kind)},
            "distinct_loci": len({q.locus for q in questions} | {q.partner for q in questions
                                                                 if q.partner})}


def load(path) -> list[Question]:
    """Read a questions.tsv back.

    Stage B runs off the written artefact rather than re-deriving from the rules,
    so what is retrieved for is exactly what was reported -- and a hand-edited
    question file works without special handling.
    """
    import csv
    from pathlib import Path

    out: list[Question] = []
    with Path(path).open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            drugs = (row.get("drugs") or "").strip()
            raised = (row.get("raised_by") or "").strip()
            out.append(Question(
                kind=(row.get("kind") or "").strip(),
                locus=(row.get("locus") or "").strip(),
                partner="" if (row.get("partner") or "NA").strip() == "NA"
                        else row["partner"].strip(),
                drugs=() if drugs in ("", "NA") else tuple(drugs.split("|")),
                text=(row.get("question") or "").strip(),
                raised_by=[] if raised in ("", "NA") else raised.split("|")))
    return out
