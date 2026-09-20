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
from itertools import combinations
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
# Two more ways a rule can describe treatment history rather than biology. Both
# are about anchors only: loci already known to confer resistance.
#   co_selection  anchors sharing a drug -- treating with it selects them together
#   multidrug     anchors for DIFFERENT drugs -- an MDR isolate, which in TB is
#                 defined as resistance to at least isoniazid and rifampicin, so
#                 katG + rpoB says "this isolate is MDR", not "these loci interact"
SIG_CO_SELECTION = "confounded:co_selection"
SIG_MULTIDRUG = "confounded:multidrug"
# Two conditions that one variant can satisfy. Not a confounded relationship --
# no relationship at all, the same observation entered twice, which reads as
# epistasis because the rule names two features.
SIG_ALIASED = "confounded:feature_overlap"
SIG_DISCORDANT = "discordant"
SIG_SUSCEPTIBLE = "known:susceptible_consistent"
SIG_UNKNOWN = "unknown"

# Most actionable first. A lineage confound leads because it undermines whatever
# else the rule looks like: if both loci mark the same lineage, the isolates
# carry both alleles by descent and the pattern needs no mechanism at all.
# Aliasing leads everything. A lineage confound says the co-occurrence has a
# non-biological cause; aliasing says there may be no co-occurrence to explain.
PRECEDENCE = (SIG_ALIASED, SIG_CONFOUNDED, SIG_COMPENSATION, SIG_ALT_ROUTE, SIG_DISCORDANT,
              SIG_UNKNOWN, SIG_MULTIDRUG, SIG_CO_SELECTION, SIG_RESISTANCE,
              SIG_SUSCEPTIBLE)

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
    # (anchor, candidate, drugs) pairs where the candidate was assessed for a
    # drug the anchor confers resistance to, and the pairs where it was not.
    compensation_pairs: list[tuple[str, str, list[str]]] = field(default_factory=list)
    drug_mismatched: list[tuple[str, str]] = field(default_factory=list)
    # Anchor pairs sharing a drug, and anchor pairs for disjoint drugs.
    co_selected: list[tuple[str, str, list[str]]] = field(default_factory=list)
    multidrug: list[tuple[str, str, list[str]]] = field(default_factory=list)
    # A single locus graded for more than one drug: one mechanism, several drugs.
    cross_resistant: list[tuple[str, list[str]]] = field(default_factory=list)
    # (a, b, kind, detail) for condition pairs one variant could satisfy.
    aliased: list[tuple[str, str, str, str]] = field(default_factory=list)
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
            "aliased": "|".join(f"{a}~{b}:{kind}" for a, b, kind, _ in self.aliased) or "NA",
            "verdict": self.verdict,
        }


SAME_FEATURE = "same_feature"
OVERLAPPING = "overlapping_spans"
ADJACENT = "adjacent_undetermined"


def aliasing(conditions: list[ConditionEvidence]) -> list[tuple[str, str, str, str]]:
    """Condition pairs that one variant could satisfy.

    Three kinds, and they are not equally certain:

    * `same_feature` -- both resolve to the same locus or the same interval, so
      any variant satisfying one satisfies the other. Certain.
    * `overlapping_spans` -- the resolved spans intersect. H37Rv has 917
      overlapping consecutive gene pairs, mostly 4 bp start/stop junctions but
      56 of at least 50 bp, and a non-synonymous variant in the intersection is
      annotated to both genes. Certain if such a variant exists.
    * `adjacent_undetermined` -- one condition is a region 5' or 3' of a gene and
      the other IS one of that region's flanking genes. A caller that reports
      upstream variants within a window -- SnpEff's default is 5,000 bp against a
      median H37Rv gap of 81 bp -- annotates a variant inside the flank to both
      features. Whether any did is in the VCF's distance field, which this run
      does not have, so the answer is undetermined rather than either.

    Keyed on adjacency in the rule rather than on a window size, so it fires on
    the specific pairing that can alias and not on every gene with a near
    neighbour.
    """
    out: list[tuple[str, str, str, str]] = []
    for index, first in enumerate(conditions):
        for second in conditions[index + 1:]:
            a, b = first.feature, second.feature
            if not (a.resolved and b.resolved):
                continue
            if a.locus == b.locus:
                out.append((first.label, second.label, SAME_FEATURE,
                            f"both name {a.locus}"))
                continue
            left, right = a.span, b.span
            if left and right and left[0] <= right[1] and right[0] <= left[1]:
                width = min(left[1], right[1]) - max(left[0], right[0]) + 1
                out.append((first.locus, second.locus, OVERLAPPING, f"{width}bp shared"))
                continue
            for region, gene in ((a, b), (b, a)):
                flanks = region.flanks
                if (flanks and gene.gene is not None and gene.gene.locus in flanks
                        and region.matched_by in ("upstream", "downstream")):
                    out.append((region.locus, gene.locus, ADJACENT,
                                f"{gene.gene.locus} bounds this region"))
                    break
    return out


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

    # Pairwise rather than over all anchors at once: three anchors where two share
    # isoniazid and the third is rifampicin is both co-selection and multidrug,
    # and an intersection over all three would report neither.
    for first, second in combinations(anchors, 2):
        left, right = set(first.catalogue.drugs), set(second.catalogue.drugs)
        both = sorted(left & right)
        # Symmetric, so ordered by name rather than by the order the rule listed
        # them: the same pair must produce the same row from either spelling.
        one, two = sorted((first.locus, second.locus))
        if both:
            result.co_selected.append((one, two, both))
        else:
            result.multidrug.append((one, two, sorted(left | right)))
    result.cross_resistant = [(c.locus, sorted(c.catalogue.drugs))
                              for c in anchors if len(c.catalogue.drugs) > 1]
    result.aliased = aliasing(conditions)

    # A compensator is only plausible for the drug its locus was catalogued
    # under. rpoA and rpoC are assessed for rifampicin only; pairing either with
    # an isoniazid anchor and calling it compensation is a false positive the
    # catalogue can rule out, and STRING will happily supply a weak edge between
    # any two well-studied genes to support it.
    for anchor in anchors:
        for candidate in candidates:
            shared = sorted(set(anchor.catalogue.drugs) & set(candidate.catalogue.assessed_drugs))
            if shared:
                result.compensation_pairs.append((anchor.locus, candidate.locus, shared))
            else:
                result.drug_mismatched.append((anchor.locus, candidate.locus))
            link = links.get((anchor.locus, candidate.locus)) or \
                links.get((candidate.locus, anchor.locus))
            if link is not None and shared:
                result.links.append((anchor.locus, candidate.locus, link))

    signatures: list[str] = []
    kind = phenotype(rule.predicted_class)

    if result.aliased:
        signatures.append(SIG_ALIASED)
    if result.shared_lineages and len(present) > 1:
        signatures.append(SIG_CONFOUNDED)

    if kind == RESISTANT:
        # Each test states one fact about the rule and they are NOT exclusive: a
        # rule can carry a known anchor AND a locus nothing accounts for, which is
        # a known mechanism with something unexplained riding along. Written as
        # mutually exclusive branches, that combination matched none of them and
        # fell through to an empty signature list -- losing the category the whole
        # exercise is aimed at.
        if anchors:
            signatures.append(SIG_RESISTANCE)
        if result.co_selected:
            signatures.append(SIG_CO_SELECTION)
        if result.multidrug:
            signatures.append(SIG_MULTIDRUG)
        if result.compensation_pairs:
            signatures.append(SIG_COMPENSATION)
        if negated_anchors and (anchors or candidates or unknowns):
            signatures.append(SIG_ALT_ROUTE)
        if unknowns or (candidates and not anchors):
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
    # Nothing may fall through silently: an empty list means a rule shape the
    # tests above do not describe, and it should be visible as unknown rather
    # than arrived at by default.
    if not result.signatures:
        result.signatures = [SIG_UNKNOWN]
    result.primary = result.signatures[0]
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

    if result.primary == SIG_ALIASED:
        certain = [x for x in result.aliased if x[2] != ADJACENT]
        possible = [x for x in result.aliased if x[2] == ADJACENT]
        if certain:
            parts.append(
                "ALIASED: " + "; ".join(f"{a} and {b} ({detail})" for a, b, _, detail in certain)
                + ". One variant can satisfy both conditions, so this rule may be one "
                  "observation entered twice rather than two pieces of evidence -- which is "
                  "what makes it read as epistasis.")
        if possible:
            parts.append(
                "POSSIBLY ALIASED: " + "; ".join(f"{a} and {b} ({detail})"
                                                 for a, b, _, detail in possible)
                + ". A caller reporting upstream variants within a window annotates a variant "
                  "inside the flanking gene to both features. Whether any did is in the "
                  "caller's distance field, which this run does not have, so this is "
                  "undetermined rather than a finding either way.")
    elif result.primary == SIG_CONFOUNDED:
        parts.append(
            f"CONFOUND: every locus the rule requires present contains positions defining "
            f"{', '.join(result.shared_lineages)}, so isolates of that lineage carry them "
            f"together by descent. Population structure explains the co-occurrence without "
            f"any mechanism, and must be excluded first.")
    elif result.primary == SIG_COMPENSATION:
        kinds = {lk.kind for _, _, lk in result.links}
        if "string_beyond_textmining" in kinds or "kegg_shared_pathway" in kinds \
                or "adjacent" in kinds:
            link = (f" They are linked by evidence beyond literature co-mention "
                    f"({', '.join(sorted(kinds))}), which is what distinguishes compensation "
                    f"from coincidence.")
        elif kinds:
            # STRING's textmining channel IS co-mention in papers, so an edge
            # supported only by it is the same evidence a literature search would
            # return -- not a second, independent line of it.
            link = (" The only link between them is STRING's textmining channel, which is "
                    "co-mention in papers rather than independent support.")
        else:
            link = (" No link between them was found in the structured sources, so this is "
                    "co-occurrence rather than evidence of compensation.")
        pairs = "; ".join(f"{a} + {b} ({', '.join(drugs)})"
                          for a, b, drugs in result.compensation_pairs)
        parts.append(
            f"Compensation-shaped: {pairs}. The first carries graded-associated variants for "
            f"that drug; the second was assessed against it and graded as conferring none -- "
            f"the shape the known compensatory loci have.{link}")
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
    elif result.primary == SIG_MULTIDRUG:
        pairs = "; ".join(f"{a} + {b} ({', '.join(drugs)})" for a, b, drugs in result.multidrug)
        parts.append(
            f"Multi-drug: {pairs}. These loci confer resistance to different drugs, so an "
            f"isolate carrying both is one that acquired resistance to each under combination "
            f"therapy. In M. tuberculosis, MDR is DEFINED as resistance to at least isoniazid "
            f"and rifampicin, so their co-occurrence is the diagnosis rather than a "
            f"relationship between the loci.")
    elif result.primary == SIG_CO_SELECTION:
        pairs = "; ".join(f"{a} + {b} ({', '.join(drugs)})" for a, b, drugs in result.co_selected)
        parts.append(
            f"Co-selection: {pairs}. Treating with that drug selects every locus conferring "
            f"resistance to it at once, so these co-occur across isolates under treatment "
            f"rather than through any link between them.")
    elif result.primary == SIG_RESISTANCE:
        drugs = sorted({d for c in anchors for d in c.catalogue.drugs})
        parts.append(
            f"Recapitulates the catalogue: {_names(anchors)} is graded resistance-associated "
            f"for {', '.join(drugs)}. The condition says the locus carries some qualifying "
            f"variant, not that it carries a graded one.")
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

    if SIG_ALIASED in result.signatures and result.primary != SIG_ALIASED:
        parts.append("Also possibly aliased: "
                     + "; ".join(f"{a}/{b}" for a, b, _, _ in result.aliased) + ".")
    if result.cross_resistant:
        parts.append(
            "Cross-resistance: "
            + "; ".join(f"{locus} is graded for {', '.join(drugs)}" for locus, drugs in
                        result.cross_resistant)
            + " -- one mechanism covering several drugs, so resistance to all of them can "
              "follow from this locus alone.")
    if result.drug_mismatched:
        parts.append(
            "Not counted as compensation: "
            + "; ".join(f"{a} + {b}" for a, b in result.drug_mismatched)
            + " -- the second locus was never assessed against the drug the first confers "
              "resistance to, so it cannot be compensating for it.")
    if SIG_CONFOUNDED in result.signatures and result.primary != SIG_CONFOUNDED:
        parts.append(f"Also lineage-confounded on {', '.join(result.shared_lineages)}.")
    return " ".join(parts)
