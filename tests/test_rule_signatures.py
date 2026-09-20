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
    assert [p.position for p in status.placed] == [20785]
    assert status.drugs == ["isoniazid"]
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
    assert "functionally linked" in result.verdict


def test_an_anchor_alone_recapitulates_the_catalogue(sources, tmp_path):
    rule = _rule("Rv1908c=1", "R", tmp_path)
    result = classify(rule, evidence_for(rule, sources))
    assert result.primary == SIG_RESISTANCE
    assert "Recapitulates the catalogue" in result.verdict
    assert "not that it carries a graded one" in result.verdict


def test_two_anchors_for_one_drug_read_as_co_selection(sources, tmp_path):
    """Treating with isoniazid selects every isoniazid locus at once, so their
    co-occurrence needs no relationship between them."""
    rule = _rule("Rv1908c=1 AND Rv1484=1", "R", tmp_path)
    result = classify(rule, evidence_for(rule, sources))
    assert result.shared_drugs == ["isoniazid"]
    assert "co-selection" in result.verdict


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
