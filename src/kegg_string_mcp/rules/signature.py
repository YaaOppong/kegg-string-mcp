"""Classify a rule by what the structured sources already account for.

Pure: it takes the evidence gathered for each condition and returns roles, a
signature and a sentence. Nothing here fetches, and nothing here is decided by a
model -- a rule's classification has to be reproducible from the sources, or the
count of "novel" rules means whatever the model felt that day.

The roles come from two things about each condition: the state the rule asserts,
and what the catalogue holds for the locus.

    state 1  anchor            a graded-associated locus: known pharmacology
    state 1  compensator       assessed, none associated: the rpoC/rpoA/ahpC shape
    state 1  unknown           absent from the catalogue: never assessed
    state 0  negated anchor    resistance WITHOUT the canonical locus
    state 0  negated           reference at a locus nothing is known about

`=0` is an assertion, not silence. A rule that does not mention katG says nothing
about it; a rule that says `katG=0` says every isolate it covers carries no
qualifying katG variant, and that is where alternative-route patterns live.

Class matters as much as state. For an R-predicting rule an anchor is
explanatory; for an S-predicting rule the same anchor at state 1 is *discordant*
-- a known resistance locus varying in susceptible isolates, which is either a
non-causal variant, a suppressor, or a calling artefact, and all three are worth
looking at.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from kegg_string_mcp.rules.catalogue import ANCHOR, ASSESSED_NEGATIVE, CatalogueStatus
from kegg_string_mcp.rules.features import UNRESOLVED, Feature
from kegg_string_mcp.rules.parse import Condition, Rule

# roles
ROLE_ANCHOR = "anchor"
ROLE_COMPENSATOR = "compensator_candidate"
ROLE_UNKNOWN = "unknown"
ROLE_NEGATED_ANCHOR = "negated_anchor"
ROLE_NEGATED = "negated"
ROLE_UNRESOLVED = "unresolved"

# signatures
SIG_RESISTANCE = "known:resistance"
SIG_COMPENSATION = "known:compensation"
SIG_ALT_ROUTE = "known:alt_route"
SIG_CONFOUNDED = "confounded:lineage"
SIG_DISCORDANT = "discordant"
SIG_SUSCEPTIBLE = "known:susceptible_consistent"
SIG_UNKNOWN = "unknown"

# Most actionable first. A lineage confound leads because it undermines whatever
# else the rule looks like: if both loci mark the same lineage, the isolates
# carry both alleles by descent and the pattern needs no mechanism at all.
PRECEDENCE = (SIG_CONFOUNDED, SIG_COMPENSATION, SIG_ALT_ROUTE, SIG_DISCORDANT,
              SIG_UNKNOWN, SIG_RESISTANCE, SIG_SUSCEPTIBLE)

RESISTANT = "R"
SUSCEPTIBLE = "S"


def phenotype(predicted_class: str) -> str:
    """R, S, or "" when the label is neither -- in which case no R/S reading is
    asserted rather than one being assumed."""
    first = (predicted_class or "").strip()[:1].upper()
    return first if first in (RESISTANT, SUSCEPTIBLE) else ""


@dataclass(frozen=True)
class Link:
    """A functional relationship between two loci, from the structured sources.

    Co-occurrence alone does not make a compensator. What makes rpoB + rpoC
    convincing is that they are subunits of the same complex, which STRING scores
    high with evidence beyond literature co-mention.
    """

    kind: str        # string_beyond_textmining | shared_specific_pathway | adjacent
    detail: str = ""


@dataclass
class ConditionEvidence:
    condition: Condition
    feature: Feature
    catalogue: CatalogueStatus
    lineages: tuple[str, ...] = ()
    role: str = ""

    @property
    def label(self) -> str:
        return self.condition.label

    @property
    def locus(self) -> str:
        return self.feature.locus or self.condition.label

    def to_dict(self) -> dict[str, Any]:
        return ({"label": self.label, "state": self.condition.state, "role": self.role,
                 "locus": self.locus, "kind": self.feature.kind,
                 "matched_by": self.feature.matched_by,
                 "lineages": "|".join(self.lineages) or "NA"}
                | self.catalogue.to_dict())


def assign_role(evidence: ConditionEvidence) -> str:
    if evidence.feature.kind == UNRESOLVED:
        return ROLE_UNRESOLVED
    anchor = evidence.catalogue.status == ANCHOR
    if evidence.condition.state == 1:
        if anchor:
            return ROLE_ANCHOR
        if evidence.catalogue.status == ASSESSED_NEGATIVE:
            return ROLE_COMPENSATOR
        return ROLE_UNKNOWN                     # ABSENT: never assessed
    return ROLE_NEGATED_ANCHOR if anchor else ROLE_NEGATED


@dataclass
class RuleSignature:
    rule: Rule
    conditions: list[ConditionEvidence] = field(default_factory=list)
    signatures: list[str] = field(default_factory=list)
    primary: str = ""
    links: list[tuple[str, str, Link]] = field(default_factory=list)
    shared_lineages: list[str] = field(default_factory=list)
    shared_drugs: list[str] = field(default_factory=list)
    verdict: str = ""

    def of_role(self, role: str) -> list[ConditionEvidence]:
        return [c for c in self.conditions if c.role == role]

    def to_dict(self, rename: dict[str, str] | None = None) -> dict[str, Any]:
        return self.rule.to_dict(rename) | {
            "primary_signature": self.primary,
            "signatures": "|".join(self.signatures) or "NA",
            "anchors": "|".join(c.locus for c in self.of_role(ROLE_ANCHOR)) or "NA",
            "compensator_candidates":
                "|".join(c.locus for c in self.of_role(ROLE_COMPENSATOR)) or "NA",
            "unknown_loci": "|".join(c.locus for c in self.of_role(ROLE_UNKNOWN)) or "NA",
            "negated_anchors":
                "|".join(c.locus for c in self.of_role(ROLE_NEGATED_ANCHOR)) or "NA",
            "unresolved": "|".join(c.label for c in self.of_role(ROLE_UNRESOLVED)) or "NA",
            "links": "|".join(f"{a}~{b}:{link.kind}" for a, b, link in self.links) or "NA",
            "shared_lineages": "|".join(self.shared_lineages) or "NA",
            "shared_drugs": "|".join(self.shared_drugs) or "NA",
            "verdict": self.verdict,
        }


def _shared(values: list[tuple[str, ...]]) -> list[str]:
    if len(values) < 2 or not all(values):
        return []
    common = set(values[0])
    for other in values[1:]:
        common &= set(other)
    return sorted(common)


def classify(rule: Rule, conditions: list[ConditionEvidence],
             links: dict[tuple[str, str], Link] | None = None) -> RuleSignature:
    """Roles, signatures and a verdict for one rule."""
    links = links or {}
    for evidence in conditions:
        evidence.role = assign_role(evidence)
    result = RuleSignature(rule=rule, conditions=conditions)

    present = [c for c in conditions if c.condition.state == 1]
    anchors = result.of_role(ROLE_ANCHOR)
    candidates = result.of_role(ROLE_COMPENSATOR)
    unknowns = result.of_role(ROLE_UNKNOWN)
    negated_anchors = result.of_role(ROLE_NEGATED_ANCHOR)

    # Confounds, over the conditions the rule asserts are present. Two loci
    # marking one lineage are inherited together; two anchors for one drug are
    # co-selected by treating with it.
    result.shared_lineages = _shared([c.lineages for c in present])
    result.shared_drugs = _shared([tuple(c.catalogue.drugs) for c in anchors])

    for anchor in anchors:
        for candidate in candidates:
            link = links.get((anchor.locus, candidate.locus)) or \
                links.get((candidate.locus, anchor.locus))
            if link is not None:
                result.links.append((anchor.locus, candidate.locus, link))

    signatures: list[str] = []
    kind = phenotype(rule.predicted_class)

    if result.shared_lineages and len(present) > 1:
        signatures.append(SIG_CONFOUNDED)

    if kind == RESISTANT:
        if anchors and candidates:
            signatures.append(SIG_COMPENSATION)
        if negated_anchors and (anchors or candidates or unknowns):
            signatures.append(SIG_ALT_ROUTE)
        if anchors and not candidates and not unknowns:
            signatures.append(SIG_RESISTANCE)
        if not anchors and (unknowns or candidates):
            signatures.append(SIG_UNKNOWN)
    elif kind == SUSCEPTIBLE:
        if anchors:
            signatures.append(SIG_DISCORDANT)
        elif negated_anchors:
            signatures.append(SIG_SUSCEPTIBLE)
        else:
            signatures.append(SIG_UNKNOWN)
    else:
        signatures.append(SIG_UNKNOWN)

    result.signatures = [s for s in PRECEDENCE if s in signatures]
    result.primary = result.signatures[0] if result.signatures else SIG_UNKNOWN
    result.verdict = _verdict(result, kind)
    return result


def _names(items: list[ConditionEvidence]) -> str:
    return ", ".join(c.locus for c in items)


def _verdict(result: RuleSignature, kind: str) -> str:
    """A defensible sentence, stating what was used to reach it.

    Never causal: a rule is an association in one cohort, and the wording must
    not let that drift into a mechanism the sources do not support.
    """
    anchors = result.of_role(ROLE_ANCHOR)
    candidates = result.of_role(ROLE_COMPENSATOR)
    unknowns = result.of_role(ROLE_UNKNOWN)
    negated = result.of_role(ROLE_NEGATED_ANCHOR)
    unresolved = result.of_role(ROLE_UNRESOLVED)

    parts: list[str] = []
    if unresolved:
        parts.append(
            f"{len(unresolved)} condition(s) did not resolve against the supplied annotation "
            f"({', '.join(c.label for c in unresolved)}), so this rule is classified on the "
            f"rest and the classification may be incomplete.")

    if result.primary == SIG_CONFOUNDED:
        parts.append(
            f"CONFOUND: every locus the rule requires present contains positions defining "
            f"{', '.join(result.shared_lineages)}, so isolates of that lineage carry them "
            f"together by descent. Population structure explains the co-occurrence without "
            f"any mechanism, and must be excluded first.")
    elif result.primary == SIG_COMPENSATION:
        link = (f" {_names(anchors)} and {_names(candidates)} are functionally linked ("
                f"{', '.join(sorted({lk.kind for _, _, lk in result.links}))}), which is what "
                f"distinguishes compensation from coincidence."
                if result.links else
                " No functional link between them was found in the structured sources, so this "
                "is co-occurrence rather than evidence of compensation.")
        parts.append(
            f"Compensation-shaped: {_names(anchors)} carries graded-associated variants for "
            f"{', '.join(result.shared_drugs or sorted({d for c in anchors for d in c.catalogue.drugs}))}, "
            f"while {_names(candidates)} is catalogued with none graded associated -- the shape "
            f"the known compensatory loci have.{link}")
    elif result.primary == SIG_ALT_ROUTE:
        parts.append(
            f"Alternative route: the rule requires {_names(negated)} to match the reference "
            f"while predicting resistance, so these isolates are resistant without a "
            f"qualifying variant at a canonical locus. Note the condition asserts no "
            f"QUALIFYING variant, not that the locus is wild-type.")
    elif result.primary == SIG_DISCORDANT:
        parts.append(
            f"Discordant: {_names(anchors)} carries graded-associated variants yet the rule "
            f"predicts susceptibility. The variant driving the condition may not be one of the "
            f"graded ones, may be suppressed, or may be a calling artefact.")
    elif result.primary == SIG_RESISTANCE:
        drugs = sorted({d for c in anchors for d in c.catalogue.drugs})
        extra = (f" Both loci are associated with {', '.join(result.shared_drugs)}, so treating "
                 f"with it selects them together: the rule may be recovering co-selection "
                 f"rather than a relationship between the loci."
                 if len(anchors) > 1 and result.shared_drugs else "")
        parts.append(
            f"Recapitulates the catalogue: {_names(anchors)} is graded resistance-associated "
            f"for {', '.join(drugs)}. The condition says the locus carries some qualifying "
            f"variant, not that it carries a graded one.{extra}")
    elif result.primary == SIG_SUSCEPTIBLE:
        parts.append(
            f"Consistent with susceptibility: the rule requires {_names(negated)} to match the "
            f"reference. Expected rather than informative.")
    else:
        loci = _names(unknowns + candidates) or "its loci"
        parts.append(
            f"Nothing in the catalogue accounts for this rule: {loci} carries no graded "
            f"resistance association. Either a mechanism the catalogue does not cover, or a "
            f"confound the structured sources cannot see. The catalogue covers 74 genes, so "
            f"absence here is unassessed rather than negative.")

    if SIG_CONFOUNDED in result.signatures and result.primary != SIG_CONFOUNDED:
        parts.append(f"Also lineage-confounded on {', '.join(result.shared_lineages)}.")
    return " ".join(parts)
