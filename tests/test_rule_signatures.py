"""Catalogue lookups and rule classification. No network, no model.

The shapes tested here are the ones the real catalogue shows: rpoC, rpoA and
ahpC are heavily catalogued with zero graded-associated variants, which is what
separates a compensatory locus from a resistance one, and the catalogue's
promoter variants are named for the gene downstream of them rather than for
where they physically sit.
"""

from pathlib import Path

import pytest

from kegg_string_mcp.resistance import Variant
from kegg_string_mcp.rules import parse_annotation, parse_rules
from kegg_string_mcp.rules.catalogue import ABSENT, ANCHOR, ASSESSED_NEGATIVE, Catalogue
from kegg_string_mcp.rules.evidence import Sources, classify_all, evidence_for, summarise
from kegg_string_mcp.rules.signature import (
    ROLE_ANCHOR,
    ROLE_COMPENSATOR,
    ROLE_NEGATED_ANCHOR,
    ROLE_UNKNOWN,
    SIG_ALT_ROUTE,
    SIG_COMPENSATION,
    SIG_CONFOUNDED,
    SIG_DISCORDANT,
    SIG_RESISTANCE,
    SIG_UNKNOWN,
    Link,
    classify,
)

# name, 1-based start, end, strand, symbol
GENES = [
    ("Rv0667", 100, 3617, "+", "rpoB"),
    ("Rv0668", 3700, 7700, "+", "rpoC"),
    ("Rv1908c", 8000, 10221, "-", "katG"),
    ("Rv1909c", 10400, 10800, "-", "furA"),
    ("Rv1483", 20000, 20743, "+", "fabG1"),
    ("Rv1484", 20800, 21608, "+", "inhA"),
    ("Rv9000", 30000, 30500, "+", ""),
]


def _variant(gene, mutation, drug, confidence):
    return Variant(gene=gene, mutation=mutation, drug=drug, confidence=confidence,
                   source="WHO catalogue v2", comment="")


CATALOGUE = {
    # An anchor: graded-associated variants exist.
    "rpoB": [_variant("rpoB", "p.Ser450Leu", "rifampicin", "Assoc w R"),
             _variant("rpoB", "p.His445Tyr", "rifampicin", "Uncertain significance")],
    # Assessed and negative: the compensatory shape.
    "rpoC": [_variant("rpoC", "p.Val483Gly", "rifampicin", "Uncertain significance"),
             _variant("rpoC", "p.Gly332Arg", "rifampicin", "Uncertain significance")],
    "katG": [_variant("katG", "p.Ser315Thr", "isoniazid", "Assoc w R"),
             _variant("katG", "c.-20A>G", "isoniazid", "Uncertain significance")],
    # The promoter case: named for inhA, physically upstream of it.
    "inhA": [_variant("inhA", "c.-15C>T", "isoniazid", "Assoc w R"),
             _variant("inhA", "c.-15C>T", "ethionamide", "Assoc w R"),
             _variant("inhA", "p.Ser94Ala", "isoniazid", "Uncertain significance")],
}


@pytest.fixture
def annotation(tmp_path: Path):
    path = tmp_path / "genes.bed"
    path.write_text("\n".join(
        f"AL123456\t{start - 1}\t{end}\t{name}\t.\t{strand}\t{symbol}"
        for name, start, end, strand, symbol in GENES) + "\n")
    return parse_annotation(path)


@pytest.fixture
def catalogue(annotation):
    return Catalogue(CATALOGUE, annotation)


@pytest.fixture
def sources(annotation, catalogue):
    return Sources(annotation=annotation, catalogue=catalogue)


def _rule(conditions: str, predicted: str, tmp_path: Path, name="r.tsv"):
    path = tmp_path / name
    path.write_text(f"conditions\tpredicted_class\n{conditions}\t{predicted}\n")
    return parse_rules(path)[0]


# --- the catalogue ---------------------------------------------------------


def test_the_catalogue_keeps_three_states_apart(catalogue):
    """An anchor, a locus assessed and found negative, and a locus never
    assessed. Collapsing the last two would call every unstudied gene 'not
    resistance-associated', which the catalogue does not say."""
    assert catalogue.for_locus("Rv0667", "rpoB").status == ANCHOR
    assert catalogue.for_locus("Rv0668", "rpoC").status == ASSESSED_NEGATIVE
    assert catalogue.for_locus("Rv9000", "").status == ABSENT


def test_an_anchor_reports_its_drugs_and_keeps_the_granularity_caveat(catalogue):
    status = catalogue.for_locus("Rv0667", "rpoB")
    assert status.drugs == ["rifampicin"]
    assert (status.associated, status.catalogued) == (1, 2)
    assert "about the GENE" in status.note


def test_an_upstream_variant_is_placed_in_front_of_a_forward_strand_gene(catalogue, annotation):
    """`inhA c.-15` is 15 bases 5' of inhA, which on the forward strand is a lower
    coordinate -- and lands in the gap behind it, not inside inhA."""
    interval = annotation.intergenic("Rv1483-Rv1484")
    status = catalogue.for_interval(interval)
    assert status.status == ANCHOR
    # One catalogue row per (variant, drug), so a variant graded for two drugs
    # places twice at the same base.
    assert {p.position for p in status.placed} == {20785}
    assert status.drugs == ["ethionamide", "isoniazid"]
    # The catalogue named it for inhA; it physically sits between fabG1 and inhA.
    assert status.placed[0].named_for == "inhA"
    assert "downstream of it" in status.note


def test_an_upstream_variant_on_the_reverse_strand_goes_the_other_way(catalogue, annotation):
    """katG is on the reverse strand, so 5' of it is a HIGHER coordinate. Getting
    this backwards puts every promoter variant inside the gene instead."""
    interval = annotation.intergenic("Rv1908c-Rv1909c")
    status = catalogue.for_interval(interval)
    assert [p.position for p in status.placed] == [10241]     # katG end 10221 + 20
    assert status.status == ASSESSED_NEGATIVE                  # graded uncertain


def test_an_interval_with_no_catalogued_variant_says_so_without_claiming_a_negative(
        catalogue, annotation):
    status = catalogue.for_interval(annotation.intergenic("Rv0667-Rv0668"))
    assert status.status == ABSENT
    assert "not evidence that variation here is unimportant" in status.note


# --- roles -----------------------------------------------------------------


def test_state_and_catalogue_status_together_give_the_role(sources, tmp_path):
    rule = _rule("Rv0667=1 AND Rv0668=1 AND Rv9000=1 AND Rv1908c=0", "R", tmp_path)
    roles = {e.condition.label: e.role
             for e in classify(rule, evidence_for(rule, sources)).conditions}
    assert roles == {"Rv0667": ROLE_ANCHOR,            # anchor, asserted present
                     "Rv0668": ROLE_COMPENSATOR,       # assessed negative, present
                     "Rv9000": ROLE_UNKNOWN,           # never assessed, present
                     "Rv1908c": ROLE_NEGATED_ANCHOR}   # anchor, asserted ABSENT


# --- signatures ------------------------------------------------------------


def test_anchor_plus_assessed_negative_is_compensation_shaped(sources, tmp_path):
    """rpoB + rpoC: the canonical rifampicin compensation pair."""
    rule = _rule("Rv0667=1 AND Rv0668=1", "R", tmp_path)
    result = classify(rule, evidence_for(rule, sources))
    assert result.primary == SIG_COMPENSATION
    assert "Compensation-shaped" in result.verdict
    # Without a functional link it is co-occurrence, and the verdict says so.
    assert "co-occurrence rather than evidence" in result.verdict


def test_a_functional_link_is_what_separates_compensation_from_coincidence(sources, tmp_path):
    rule = _rule("Rv0667=1 AND Rv0668=1", "R", tmp_path)
    links = {("Rv0667", "Rv0668"): Link("string_beyond_textmining", "combined 0.999")}
    result = classify(rule, evidence_for(rule, sources), links)
    assert result.primary == SIG_COMPENSATION
    assert result.links and result.links[0][2].kind == "string_beyond_textmining"
    assert "beyond literature co-mention" in result.verdict


def test_an_anchor_alone_recapitulates_the_catalogue(sources, tmp_path):
    rule = _rule("Rv1908c=1", "R", tmp_path)
    result = classify(rule, evidence_for(rule, sources))
    assert result.primary == SIG_RESISTANCE
    assert "Recapitulates the catalogue" in result.verdict
    assert "not that it carries a graded one" in result.verdict


def test_two_anchors_for_one_drug_read_as_co_selection(sources, tmp_path):
    """Treating with isoniazid selects every isoniazid locus at once, so their
    co-occurrence needs no relationship between them."""
    from kegg_string_mcp.rules.signature import SIG_CO_SELECTION

    rule = _rule("Rv1908c=1 AND Rv1484=1", "R", tmp_path)
    result = classify(rule, evidence_for(rule, sources))
    assert result.primary == SIG_CO_SELECTION
    assert result.co_selected == [("Rv1484", "Rv1908c", ["isoniazid"])]
    assert "Co-selection" in result.verdict
    assert "selects every locus conferring resistance to it" in result.verdict


def test_anchors_for_different_drugs_describe_a_multidrug_isolate(sources, tmp_path):
    """katG + rpoB is MDR-TB by definition. The co-occurrence is the diagnosis,
    not a relationship between the loci -- and those rules will be everywhere."""
    from kegg_string_mcp.rules.signature import SIG_MULTIDRUG

    rule = _rule("Rv1908c=1 AND Rv0667=1", "R", tmp_path)
    result = classify(rule, evidence_for(rule, sources))
    assert result.primary == SIG_MULTIDRUG
    assert result.multidrug == [("Rv0667", "Rv1908c", ["isoniazid", "rifampicin"])]
    assert result.co_selected == []
    assert "MDR is DEFINED as resistance to at least isoniazid and rifampicin" in result.verdict


def test_one_locus_graded_for_several_drugs_is_cross_resistance(sources, tmp_path):
    """inhA is graded for isoniazid AND ethionamide: one mechanism, two drugs.
    Different from two loci for two drugs, and worth saying separately."""
    rule = _rule("Rv1484=1", "R", tmp_path)
    result = classify(rule, evidence_for(rule, sources))
    assert result.cross_resistant == [("Rv1484", ["ethionamide", "isoniazid"])]
    assert "Cross-resistance" in result.verdict
    assert "one mechanism covering several drugs" in result.verdict


def test_a_negated_anchor_is_an_alternative_route(sources, tmp_path):
    """Resistance while the canonical locus is at reference: the inhA-promoter
    route taken instead of katG."""
    rule = _rule("Rv1908c=0 AND Rv1483-Rv1484=1", "R", tmp_path)
    result = classify(rule, evidence_for(rule, sources))
    assert result.primary == SIG_ALT_ROUTE
    assert "without a qualifying variant at a canonical locus" in result.verdict
    # And the condition's real meaning is preserved, not overstated.
    assert "not that the locus is wild-type" in result.verdict


def test_an_anchor_predicting_susceptibility_is_discordant(sources, tmp_path):
    rule = _rule("Rv0667=1", "S", tmp_path)
    result = classify(rule, evidence_for(rule, sources))
    assert result.primary == SIG_DISCORDANT
    assert "may not be one of the graded ones" in result.verdict


def test_a_rule_the_catalogue_cannot_account_for_is_unknown(sources, tmp_path):
    rule = _rule("Rv9000=1", "R", tmp_path)
    result = classify(rule, evidence_for(rule, sources))
    assert result.primary == SIG_UNKNOWN
    assert "covers 74 genes" in result.verdict


def test_a_lineage_confound_outranks_whatever_else_the_rule_looks_like(annotation,
                                                                       catalogue, tmp_path):
    """If both loci mark one lineage the isolates carry them together by descent,
    so the pattern needs no mechanism -- and that has to lead, not trail."""
    from kegg_string_mcp.lineage import LineageSnp

    barcode = [LineageSnp(position=200, lineage="lineage2.2.1", lineage_name="Beijing",
                          allele="A"),
               LineageSnp(position=4000, lineage="lineage2.2.1", lineage_name="Beijing",
                          allele="T")]
    sources = Sources(annotation=annotation, catalogue=catalogue, barcode=barcode)
    rule = _rule("Rv0667=1 AND Rv0668=1", "R", tmp_path)
    result = classify(rule, evidence_for(rule, sources))

    assert result.primary == SIG_CONFOUNDED
    assert result.shared_lineages == ["lineage2.2.1"]
    assert result.verdict.startswith("CONFOUND")
    # The compensation reading is still recorded, just not leading.
    assert SIG_COMPENSATION in result.signatures


def test_an_unresolved_condition_is_said_out_loud(sources, tmp_path):
    rule = _rule("Rv0667=1 AND Rv7777=1", "R", tmp_path)
    result = classify(rule, evidence_for(rule, sources))
    assert "did not resolve" in result.verdict
    assert "may be incomplete" in result.verdict


def test_the_run_summarises_by_primary_signature(sources, tmp_path):
    path = tmp_path / "many.tsv"
    path.write_text("conditions\tpredicted_class\n"
                    "Rv1908c=1\tR\n"
                    "Rv0667=1 AND Rv0668=1\tR\n"
                    "Rv9000=1\tR\n")
    results = classify_all(parse_rules(path), sources)
    assert summarise(results)["by_primary_signature"] == {
        SIG_COMPENSATION: 1, SIG_RESISTANCE: 1, SIG_UNKNOWN: 1}


def test_a_known_anchor_beside_an_unaccounted_locus_is_not_lost(sources, tmp_path):
    """The signatures are facts about a rule, not exclusive branches. Written as
    exclusive, `anchor + unaccounted locus` matched none of them and fell through
    to an empty list -- losing the very category this exists to surface: a known
    mechanism with something unexplained riding along."""
    rule = _rule("Rv1908c=1 AND Rv9000=1", "R", tmp_path)
    result = classify(rule, evidence_for(rule, sources))

    assert set(result.signatures) == {SIG_UNKNOWN, SIG_RESISTANCE}
    # Unknown leads, so the rule lands in the bucket a reviewer would open.
    assert result.primary == SIG_UNKNOWN


def test_no_rule_shape_ends_with_an_empty_signature_list(sources, tmp_path):
    """Falling through to a default hides a shape the tests do not describe."""
    shapes = ["Rv1908c=1", "Rv1908c=0", "Rv9000=0", "Rv0668=1", "Rv9000=1",
              "Rv1908c=1 AND Rv9000=1", "Rv0667=1 AND Rv0668=1 AND Rv1908c=0"]
    for index, conditions in enumerate(shapes):
        for predicted in ("R", "S", ""):
            rule = _rule(conditions, predicted, tmp_path, name=f"s{index}{predicted}.tsv")
            result = classify(rule, evidence_for(rule, sources))
            assert result.signatures, f"{conditions} -> {predicted} produced no signature"
            assert result.primary in result.signatures
            assert result.verdict


# --- per-locus annotation and links ----------------------------------------


class _Result:
    def __init__(self, payload): self._payload = payload
    def model_dump(self): return self._payload


class _Clients:
    """Stand-ins. Any client may be absent or may raise; both must be survivable."""

    def __init__(self, kegg=None, uniprot=None, string=None):
        if kegg is not None: self.kegg = kegg
        if uniprot is not None: self.uniprot = uniprot
        if string is not None: self.string = string


class _Kegg:
    def __init__(self, records): self._records = records
    def pathways(self, gene, organism="mtu"): return _Result({"records": self._records})


class _Raises:
    def pathways(self, *a, **k): raise ValueError("KEGG said no")
    def protein(self, *a, **k): raise ValueError("UniProt said no")
    def partners(self, *a, **k): raise ValueError("STRING said no")


def test_a_locus_is_annotated_from_whichever_sources_answer(sources):
    from kegg_string_mcp.rules.annotate import annotate_locus

    clients = _Clients(kegg=_Kegg([{"record_id": "mtu00360", "name": "Phenylalanine"}]))
    out = annotate_locus("Rv1908c", "katG", sources, clients)
    assert out.pathways == [("mtu00360", "Phenylalanine")]
    assert out.catalogue.status == ANCHOR
    assert out.gene.length == 2222
    # No UniProt or STRING client was supplied, so those fields stay empty and
    # nothing pretends the lookup happened.
    assert out.product == "" and out.partners == []


def test_a_source_that_fails_does_not_take_the_locus_down(sources):
    """Annotating a hundred loci must not be all-or-nothing."""
    from kegg_string_mcp.rules.annotate import annotate_locus

    out = annotate_locus("Rv1908c", "katG", sources, _Clients(kegg=_Raises()))
    assert out.pathways == []
    assert any("KEGG pathways unavailable" in n for n in out.notes)
    assert out.catalogue.status == ANCHOR       # the free lookups still ran


def test_only_pairs_that_share_a_rule_are_linked(sources, tmp_path):
    """Asking about every pair of a hundred loci is thousands of lookups for
    links no rule would use."""
    from kegg_string_mcp.rules.annotate import co_occurring

    path = tmp_path / "many.tsv"
    path.write_text("conditions\tpredicted_class\n"
                    "Rv0667=1 AND Rv0668=1\tR\n"
                    "Rv1908c=1 AND Rv9000=1\tR\n")
    rules = parse_rules(path)
    rename = {label: label for rule in rules for label in rule.labels}
    assert co_occurring(rules, rename) == {("Rv0667", "Rv0668"), ("Rv1908c", "Rv9000")}


def test_link_kinds_keep_the_string_channel_apart(sources):
    """An edge supported only by textmining is co-mention in papers -- the same
    evidence a literature search finds, not a second line of it."""
    from kegg_string_mcp.rules.annotate import LocusAnnotation, compute_links

    annotations = {"Rv0667": LocusAnnotation("Rv0667"), "Rv0668": LocusAnnotation("Rv0668")}
    edges = {("Rv0667", "Rv0668"): {"combined_score": 0.999,
                                    "evidence_beyond_textmining": True}}
    links = compute_links({("Rv0667", "Rv0668")}, annotations, sources, edges=edges)
    assert [link.kind for link in links[("Rv0667", "Rv0668")]] == ["string_beyond_textmining"]

    edges[("Rv0667", "Rv0668")]["evidence_beyond_textmining"] = False
    links = compute_links({("Rv0667", "Rv0668")}, annotations, sources, edges=edges)
    assert [link.kind for link in links[("Rv0667", "Rv0668")]] == ["string_textmining_only"]


def test_a_container_pathway_is_not_a_link(sources):
    """mtu01100 holds a sixth of the genome. Sharing it is a base rate."""
    from kegg_string_mcp.rules.annotate import LocusAnnotation, compute_links

    annotations = {
        "Rv0667": LocusAnnotation("Rv0667", pathways=[("mtu01100", "Metabolic"),
                                                      ("mtu00983", "Drug metabolism")]),
        "Rv0668": LocusAnnotation("Rv0668", pathways=[("mtu01100", "Metabolic"),
                                                      ("mtu00983", "Drug metabolism")])}
    links = compute_links({("Rv0667", "Rv0668")}, annotations, sources,
                          pathway_sizes={"mtu01100": 698, "mtu00983": 11}, genome_size=4008)
    kinds = [link.detail for link in links[("Rv0667", "Rv0668")]]
    assert any("mtu00983" in d for d in kinds)
    assert not any("mtu01100" in d for d in kinds)


def test_neighbouring_loci_are_reported_with_their_distance(sources, annotation):
    """katG and furA are 6 bp apart. Reported as a distance, not as an operon
    call -- that needs strand and expression evidence this does not have."""
    from kegg_string_mcp.rules.annotate import LocusAnnotation, compute_links

    annotations = {locus: LocusAnnotation(locus, gene=annotation.gene(locus))
                   for locus in ("Rv1908c", "Rv1909c")}
    links = compute_links({("Rv1908c", "Rv1909c")}, annotations, sources)
    adjacent = [link for link in links[("Rv1908c", "Rv1909c")] if link.kind == "adjacent"]
    assert adjacent and adjacent[0].detail == "178bp"


def test_a_compensator_appears_beside_an_anchor_and_rarely_alone(sources, tmp_path):
    """The strongest compensation evidence available comes from the rule set
    itself, not from any external source, and costs one pass."""
    from kegg_string_mcp.rules.annotate import LocusAnnotation, population_counts

    path = tmp_path / "pop.tsv"
    path.write_text("conditions\tpredicted_class\n"
                    "Rv1908c=1 AND Rv0668=1\tR\n"       # anchor + candidate
                    "Rv0667=1 AND Rv0668=1\tR\n"        # anchor + candidate
                    "Rv0668=1 AND Rv9000=1\tR\n")       # candidate alone
    rules = parse_rules(path)
    rename = {label: label for rule in rules for label in rule.labels}
    annotations = {locus: LocusAnnotation(locus)
                   for locus in ("Rv1908c", "Rv0667", "Rv0668", "Rv9000")}
    population_counts(rules, rename, {"Rv1908c", "Rv0667"}, annotations)

    candidate = annotations["Rv0668"]
    assert candidate.n_present == 3
    assert candidate.with_anchor == 2
    assert candidate.without_anchor == 1


def test_a_compensator_for_another_drug_is_not_compensation(sources, tmp_path):
    """katG confers isoniazid resistance; rpoC was only ever assessed against
    rifampicin. Pairing them and calling it compensation is a false positive the
    catalogue can rule out -- and STRING will supply a weak edge between any two
    well-studied genes to support it if nothing checks."""
    rule = _rule("Rv1908c=1 AND Rv0668=1", "R", tmp_path)
    links = {("Rv1908c", "Rv0668"): Link("string_textmining_only", "0.512")}
    result = classify(rule, evidence_for(rule, sources), links)

    assert result.primary != SIG_COMPENSATION
    assert result.compensation_pairs == []
    assert result.drug_mismatched == [("Rv1908c", "Rv0668")]
    assert "never assessed against the drug" in result.verdict


def test_a_textmining_only_edge_is_not_called_a_functional_link(sources, tmp_path):
    """STRING's textmining channel IS co-mention in papers, so an edge supported
    only by it is the same evidence a literature search returns -- not a second,
    independent line of it."""
    rule = _rule("Rv0667=1 AND Rv0668=1", "R", tmp_path)
    links = {("Rv0667", "Rv0668"): Link("string_textmining_only", "0.910")}
    result = classify(rule, evidence_for(rule, sources), links)

    assert result.primary == SIG_COMPENSATION
    assert "co-mention in papers rather than independent support" in result.verdict
    assert "beyond literature co-mention" not in result.verdict


def test_the_drugs_a_locus_was_assessed_against_are_recorded(catalogue):
    """Not the same as the drugs it is associated with: an assessed-negative
    locus has the second empty and the first populated, and that is what makes
    the concordance check possible."""
    rpoc = catalogue.for_locus("Rv0668", "rpoC")
    assert rpoc.drugs == []
    assert rpoc.assessed_drugs == ["rifampicin"]
    assert "assessed against rifampicin" in rpoc.note


# --- the gold set ----------------------------------------------------------


def test_every_gold_rule_parses_and_states_an_expectation():
    """An entry with no expectation is a rule that can never fail, which is worse
    than no entry at all."""
    from kegg_string_mcp.rules.gold import load

    gold = load()
    assert len(gold.rules) >= 10
    for entry in gold.rules:
        rule = entry.as_rule(1)
        assert rule.conditions, f"{entry.id} parsed to no conditions"
        assert not rule.problems, f"{entry.id}: {rule.problems}"
        assert entry.expect_primary or entry.expect_contains or entry.expect_absent \
            or entry.expect_cross_resistant, f"{entry.id} asserts nothing"
        assert entry.why, f"{entry.id} has no rationale"


def test_the_set_holds_a_negative_control():
    """Positives test that the classifier finds what is there. Only a negative
    tests that it declines what is not -- and the drug-concordance check has
    nothing else guarding it."""
    from kegg_string_mcp.rules.gold import load

    negatives = [r for r in load().rules if r.kind == "negative"]
    assert negatives
    assert all(r.expect_absent for r in negatives)


def test_the_scorer_reports_a_wrong_expectation(sources, tmp_path):
    """A check that cannot fail is not a check."""
    from kegg_string_mcp.rules.gold import GoldRule, score_one

    wrong = GoldRule(id="deliberately-wrong", kind="positive",
                     conditions="Rv1908c=1", predicted_class="R",
                     expect_primary="known:compensation", why="katG alone is not compensation")
    result = score_one(wrong, 1, sources)
    assert not result.passed
    assert "expected 'known:compensation'" in result.failures[0]


def test_a_locus_the_annotation_lacks_is_skipped_not_failed(sources):
    """The classifier was never given the chance to be right, so scoring it as a
    miss blames the wrong component."""
    from kegg_string_mcp.rules.gold import GoldRule, score_one

    entry = GoldRule(id="absent-locus", kind="positive", conditions="Rv7777=1",
                     predicted_class="R", expect_primary="known:resistance", why="n/a")
    result = score_one(entry, 1, sources)
    assert result.unresolved == ["Rv7777"]
    assert not result.failures
    assert not result.passed


def test_the_gold_set_ships_with_the_wheel():
    """`gold.py` reads a JSON file beside it, so a wheel without that file has
    code that cannot run. Same trap as the evaluation gold set and the skill."""
    pyproject = (Path(__file__).resolve().parents[1] / "pyproject.toml").read_text()
    assert '"src/kegg_string_mcp/rules/gold_rules.json"' in pyproject
    assert '"kegg_string_mcp/rules/gold_rules.json"' in pyproject


# --- stage A: the question list --------------------------------------------


def test_a_pair_already_linked_beyond_co_mention_is_not_asked_about(sources, tmp_path):
    """The structured sources already support it, so spending a retrieval on it
    buys nothing. Only the textmining-only and unlinked pairs need literature."""
    from kegg_string_mcp.rules.questions import COMPENSATION, generate

    rule = _rule("Rv0667=1 AND Rv0668=1", "R", tmp_path)
    linked = {("Rv0667", "Rv0668"): Link("string_beyond_textmining", "0.999")}
    assert not [q for q in generate([classify(rule, evidence_for(rule, sources), linked)])
                if q.kind == COMPENSATION]

    weak = {("Rv0667", "Rv0668"): Link("string_textmining_only", "0.910")}
    asked = [q for q in generate([classify(rule, evidence_for(rule, sources), weak)])
             if q.kind == COMPENSATION]
    assert len(asked) == 1


def test_a_question_uses_the_name_a_paper_would_use(sources, tmp_path):
    """Questions are search strings before they are anything else, and the
    literature says ahpC, not Rv2428. Both appear: the symbol finds the papers,
    the locus tag joins the answer back to the tables."""
    from kegg_string_mcp.rules.questions import generate

    rule = _rule("Rv1908c=1 AND Rv9000=1", "R", tmp_path)
    question = generate([classify(rule, evidence_for(rule, sources))])[0]
    assert "katG (Rv1908c)" in question.text or "Rv9000" in question.text
    assert "Mycobacterium tuberculosis" in question.text
    # Asks for the evidence, never for the verdict -- a question naming a
    # conclusion invites the model to confirm it.
    assert "Quote any passage" in question.text


def test_a_negated_anchor_still_supplies_the_drug(sources, tmp_path):
    """`katG=0 AND X=1 -> R` is still a rule about isoniazid. Dropping the drug
    because the anchor is negated leaves the question unanswerable."""
    from kegg_string_mcp.rules.questions import RESISTANCE_ROLE, generate

    rule = _rule("Rv1908c=0 AND Rv9000=1", "R", tmp_path)
    questions = generate([classify(rule, evidence_for(rule, sources))])
    role = next(q for q in questions if q.kind == RESISTANCE_ROLE)
    assert role.drugs == ("isoniazid",)
    assert "isoniazid" in role.text


def test_the_same_gap_raised_by_many_rules_is_asked_once(sources, tmp_path):
    """A classifier population reuses its loci, so asking per rule multiplies the
    cost by the redundancy in the rule set."""
    from kegg_string_mcp.rules.questions import generate

    path = tmp_path / "rep.tsv"
    path.write_text("conditions\tpredicted_class\n"
                    "Rv1908c=1 AND Rv9000=1\tR\n"
                    "Rv0667=1 AND Rv9000=1\tR\n"
                    "Rv9000=1\tR\n")
    signatures = [classify(r, evidence_for(r, sources)) for r in parse_rules(path)]
    asked = [q for q in generate(signatures) if q.locus == "Rv9000"]
    assert len(asked) == 1
    assert len(asked[0].raised_by) == 3


# --- one variant, two conditions -------------------------------------------


def test_two_conditions_on_the_same_locus_are_aliased(sources, tmp_path):
    """`Rv0667=1 AND rpoB=1` names one gene twice. A rule reads as two pieces of
    evidence when it is one."""
    from kegg_string_mcp.rules.signature import SAME_FEATURE, SIG_ALIASED

    rule = _rule("Rv0667=1 AND rpoB=1", "R", tmp_path)
    result = classify(rule, evidence_for(rule, sources))
    assert result.primary == SIG_ALIASED
    assert result.aliased[0][2] == SAME_FEATURE
    assert "one observation entered twice" in result.verdict


def test_overlapping_gene_spans_are_aliased(sources, tmp_path):
    """H37Rv has 917 overlapping consecutive gene pairs -- mostly 4 bp start/stop
    junctions, but 56 of at least 50 bp. A non-synonymous variant in the shared
    span is annotated to both genes."""
    from kegg_string_mcp.rules.signature import OVERLAPPING, SIG_ALIASED

    # Rv1908c 8000-10221 and a gene overlapping its 3' end.
    path = tmp_path / "ov.bed"
    path.write_text("AL123456\t7999\t10221\tRv1908c\t.\t-\tkatG\n"
                    "AL123456\t10100\t10500\tRv1909c\t.\t+\tfurA\n")
    from kegg_string_mcp.rules import parse_annotation
    from kegg_string_mcp.rules.catalogue import Catalogue
    from kegg_string_mcp.rules.evidence import Sources

    annotation = parse_annotation(path)
    local = Sources(annotation=annotation, catalogue=Catalogue(CATALOGUE, annotation))
    rule = _rule("Rv1908c=1 AND Rv1909c=1", "R", tmp_path, name="ov.tsv")
    result = classify(rule, evidence_for(rule, local))

    assert result.primary == SIG_ALIASED
    kind, detail = result.aliased[0][2], result.aliased[0][3]
    assert kind == OVERLAPPING and detail == "121bp shared"  # 10221 - 10101 + 1


def test_a_region_and_its_own_flanking_gene_are_undetermined(sources, tmp_path):
    """The case that matters: `upstream_X` and the gene bounding that region can
    be driven by one variant if the caller used a window. Whether any did is in
    the caller's distance field, which this run does not have -- so undetermined,
    not a finding either way."""
    from kegg_string_mcp.rules.signature import ADJACENT, SIG_ALIASED

    rule = _rule("upstream_Rv1484=1 AND Rv1483=1", "R", tmp_path)
    result = classify(rule, evidence_for(rule, sources))
    assert result.primary == SIG_ALIASED
    assert result.aliased[0][2] == ADJACENT
    assert "POSSIBLY ALIASED" in result.verdict
    assert "does not have, so this is undetermined" in result.verdict


def test_a_region_and_an_unrelated_gene_are_not_aliased(sources, tmp_path):
    """Keyed on adjacency in the rule, not on a window size, so it fires on the
    pairing that can alias and not on every gene with a near neighbour."""
    from kegg_string_mcp.rules.signature import SIG_ALIASED

    rule = _rule("upstream_Rv1484=1 AND Rv0667=1", "R", tmp_path)
    result = classify(rule, evidence_for(rule, sources))
    assert result.aliased == []
    assert SIG_ALIASED not in result.signatures


def test_aliasing_leads_whatever_else_the_rule_looks_like(sources, tmp_path):
    """A lineage confound says the co-occurrence has a non-biological cause.
    Aliasing says there may be no co-occurrence to explain."""
    from kegg_string_mcp.rules.signature import SIG_ALIASED, SIG_RESISTANCE

    rule = _rule("Rv0667=1 AND rpoB=1", "R", tmp_path)
    result = classify(rule, evidence_for(rule, sources))
    assert result.signatures[0] == SIG_ALIASED
    assert SIG_RESISTANCE in result.signatures


def _overlapping(tmp_path: Path):
    """An annotation where Rv1909c's 5' neighbour overlaps it, so the region 5'
    of Rv1909c has no interval and lies inside katG. Every unresolved label in a
    real five-cohort run was this shape, and all seven sampled overlapped."""
    from kegg_string_mcp.rules import parse_annotation
    from kegg_string_mcp.rules.catalogue import Catalogue
    from kegg_string_mcp.rules.evidence import Sources

    path = tmp_path / "overlap.bed"
    path.write_text("AL123456\t99\t3617\tRv0667\t.\t+\trpoB\n"
                    "AL123456\t7999\t10221\tRv1908c\t.\t-\tkatG\n"
                    "AL123456\t10100\t10500\tRv1909c\t.\t+\tfurA\n")
    annotation = parse_annotation(path)
    return Sources(annotation=annotation, catalogue=Catalogue(CATALOGUE, annotation))


def test_a_region_with_no_gap_names_the_gene_it_lies_inside(sources, tmp_path):
    """Every unresolved label in a real five-cohort run was this case, and the
    note said 'its neighbour' without saying which. Naming it matters twice: the
    reader learns which gene, and the aliasing check can see the relation."""
    local = _overlapping(tmp_path)
    rule = _rule("upstream_Rv1909c=1", "R", tmp_path, name="nb.tsv")
    feature = evidence_for(rule, local)[0].feature
    assert not feature.resolved
    assert feature.bounded_by == "Rv1908c"
    assert "katG (Rv1908c)" in feature.note
    # The coordinates stay unresolved on purpose.
    assert "belongs to that gene, not to this region" in feature.note


def test_a_region_inside_a_gene_aliases_with_that_gene(sources, tmp_path):
    """The most certain aliasing case, and the one that escaped every check: the
    region does not resolve, so a pairwise test over resolved spans skipped it.
    A non-synonymous variant inside the bounding gene sets both conditions."""
    from kegg_string_mcp.rules.signature import INSIDE, SIG_ALIASED

    local = _overlapping(tmp_path)
    rule = _rule("upstream_Rv1909c=1 AND Rv1908c=1", "R", tmp_path, name="al.tsv")
    result = classify(rule, evidence_for(rule, local))
    assert result.primary == SIG_ALIASED
    assert result.aliased[0][2] == INSIDE
    assert "the region lies inside Rv1908c" in result.verdict


def test_a_region_inside_one_gene_does_not_alias_with_another(sources, tmp_path):
    from kegg_string_mcp.rules.signature import SIG_ALIASED

    local = _overlapping(tmp_path)
    rule = _rule("upstream_Rv1909c=1 AND Rv0667=1", "R", tmp_path, name="na.tsv")
    result = classify(rule, evidence_for(rule, local))
    assert SIG_ALIASED not in result.signatures
