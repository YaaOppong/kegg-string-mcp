"""What is true of a rule's loci as a set. No network.

The median real rule has seven conditions and 96.6% have three or more, so the
pairwise reading is the exception rather than the rule. These tests cover the
facts that only exist at set level and the ones that must NOT leak into it.
"""

from pathlib import Path

import pytest

from kegg_string_mcp.rules import parse_annotation
from kegg_string_mcp.rules.enrichment import universe_for
from kegg_string_mcp.rules.parse import Condition
from kegg_string_mcp.rules.sets import (
    CLIQUE,
    DISCONNECTED,
    STAR,
    contiguous_runs,
    describe,
    set_evidence,
)
from kegg_string_mcp.rules.signature import ConditionEvidence

ROWS = [
    ("Rv0001", 100, 200, "lipid metabolism"),
    ("Rv0002", 250, 350, "lipid metabolism"),
    ("Rv0003", 360, 460, "lipid metabolism"),          # 10 bp after Rv0002
    ("Rv0004", 5000, 5100, "information pathways"),
    ("Rv0005", 9000, 9100, "cell wall and cell processes"),
    ("Rv0006", 12000, 12100, ""),                      # no category
]


@pytest.fixture
def annotation(tmp_path: Path):
    path = tmp_path / "sets.gff"
    path.write_text("\n".join(
        f"NC_000962.3\tMB\tCDS\t{s}\t{e}\t.\t+\t\tLocus={locus}"
        + (f";Functional_Category={cat}" if cat else "")
        for locus, s, e, cat in ROWS) + "\n")
    return parse_annotation(path)


def _conditions(annotation, spec):
    """spec: [(locus, state)] -> ConditionEvidence list."""
    from kegg_string_mcp.rules.catalogue import ABSENT, CatalogueStatus
    from kegg_string_mcp.rules.features import Resolver

    resolver = Resolver(annotation)
    return [ConditionEvidence(condition=Condition(locus, state),
                              feature=resolver.resolve(locus),
                              catalogue=CatalogueStatus(ABSENT))
            for locus, state in spec]


def test_only_the_loci_the_rule_says_are_present_count(annotation):
    """A `=0` condition asserts the locus matches the reference. Including it in
    a shared-term or common-partner statement would attribute to the set a locus
    the rule says is not varying."""
    conditions = _conditions(annotation, [("Rv0001", 1), ("Rv0002", 1), ("Rv0004", 0)])
    evidence = set_evidence(conditions, annotation, universe=universe_for(annotation))
    assert evidence.loci == ["Rv0001", "Rv0002"]
    assert evidence.k == 2


def test_a_term_every_member_carries_is_reported_separately_from_enrichment(annotation):
    """Enrichment asks whether a term is over-represented; shared_by_all asks
    whether every member has it. A three-locus set can be enriched for a term
    only two of them carry."""
    conditions = _conditions(annotation, [("Rv0001", 1), ("Rv0002", 1), ("Rv0003", 1)])
    evidence = set_evidence(conditions, annotation, universe=universe_for(annotation))
    assert evidence.shared_by_all == ["lipid metabolism"]

    mixed = set_evidence(_conditions(annotation, [("Rv0001", 1), ("Rv0004", 1)]),
                         annotation, universe=universe_for(annotation))
    assert mixed.shared_by_all == []


def test_a_member_with_no_term_empties_the_shared_set(annotation):
    """Rv0006 carries no category, so nothing is shared by ALL -- which is not
    the same as the others sharing nothing."""
    conditions = _conditions(annotation, [("Rv0001", 1), ("Rv0002", 1), ("Rv0006", 1)])
    evidence = set_evidence(conditions, annotation, universe=universe_for(annotation))
    assert evidence.shared_by_all == []
    assert evidence.enrichment.terms[0].term == "lipid metabolism"
    assert evidence.enrichment.dropped == ("Rv0006",)


def test_a_partner_every_member_shares_is_the_complex_signature(annotation):
    conditions = _conditions(annotation, [("Rv0001", 1), ("Rv0002", 1), ("Rv0003", 1)])
    partners = {"Rv0001": {"hubA", "x"}, "Rv0002": {"hubA", "y"}, "Rv0003": {"hubA", "z"}}
    evidence = set_evidence(conditions, annotation, partners=partners)
    assert evidence.common_partners == ["hubA"]
    assert "a shared neighbour of every member" in describe(evidence)
    # ...with the caveat that a capped partner list cannot rule out a hub.
    assert "cannot distinguish that from a hub" in describe(evidence)


def test_the_subgraph_shape_distinguishes_a_star_from_a_clique(annotation):
    conditions = _conditions(annotation, [(f"Rv000{i}", 1) for i in (1, 2, 3)])
    star = {("Rv0001", "Rv0002"), ("Rv0001", "Rv0003")}
    assert set_evidence(conditions, annotation, linked_pairs=star).subgraph.shape == STAR
    assert set_evidence(conditions, annotation, linked_pairs=star).subgraph.hub == "Rv0001"

    clique = star | {("Rv0002", "Rv0003")}
    assert set_evidence(conditions, annotation, linked_pairs=clique).subgraph.shape == CLIQUE
    assert set_evidence(conditions, annotation, linked_pairs=set()).subgraph.shape == DISCONNECTED


def test_contiguous_loci_are_reported_as_adjacency_not_as_an_operon(annotation):
    """Rv0002 ends at 350 and Rv0003 begins at 360. An operon call needs strand
    agreement and expression evidence this does not have."""
    runs = contiguous_runs(["Rv0001", "Rv0002", "Rv0003", "Rv0005"], annotation)
    assert runs == [["Rv0001", "Rv0002", "Rv0003"]]

    conditions = _conditions(annotation, [("Rv0001", 1), ("Rv0002", 1), ("Rv0003", 1)])
    assert "not an operon call" in describe(set_evidence(conditions, annotation))


def test_the_description_says_how_many_terms_were_tested(annotation):
    """The max statistic is the larger effect, so a q without the test count
    hides it."""
    conditions = _conditions(annotation, [("Rv0001", 1), ("Rv0002", 1), ("Rv0003", 1)])
    described = describe(set_evidence(conditions, annotation,
                                      universe=universe_for(annotation)))
    assert "term(s) tested" in described


def test_a_set_with_nothing_over_represented_says_so(annotation):
    conditions = _conditions(annotation, [("Rv0001", 1), ("Rv0004", 1), ("Rv0005", 1)])
    described = describe(set_evidence(conditions, annotation,
                                      universe=universe_for(annotation)))
    assert "no term over-represented beyond its base rate" in described


def test_a_single_locus_rule_gets_no_set_description(annotation):
    """There is no set to describe, and a sentence about one locus would repeat
    what the condition already says."""
    conditions = _conditions(annotation, [("Rv0001", 1)])
    assert describe(set_evidence(conditions, annotation)) == ""
