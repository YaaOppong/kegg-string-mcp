"""Over-representation of annotation terms among a rule's loci. No network.

The method the M. tuberculosis literature uses for gene-set enrichment: Fisher's
exact, one-tailed, Benjamini-Hochberg across the terms tested. What is tested
here is mostly the bookkeeping around it -- which loci count, which universe, and
what the numbers are allowed to claim.
"""

from pathlib import Path

import pytest

from kegg_string_mcp.rules import parse_annotation
from kegg_string_mcp.rules.enrichment import benjamini_hochberg, enrich, p_at_least, universe_for


def _gff(tmp_path: Path, rows) -> Path:
    path = tmp_path / "terms.gff"
    path.write_text("\n".join(
        f"NC_000962.3\tMB\tCDS\t{i*100+1}\t{i*100+50}\t.\t+\t\tLocus={locus};{attrs}"
        for i, (locus, attrs) in enumerate(rows)) + "\n")
    return path


@pytest.fixture
def annotation(tmp_path: Path):
    # 10 loci: 4 lipid metabolism, 3 information pathways, 2 uncategorised, 1 other.
    rows = [(f"Rv000{i}", a) for i, a in enumerate([
        "Functional_Category=lipid metabolism;PFAM=PF00109,PF02801",
        "Functional_Category=lipid metabolism;PFAM=PF00109",
        "Functional_Category=lipid metabolism;PFAM=PF00109",
        "Functional_Category=lipid metabolism",
        "Functional_Category=information pathways",
        "Functional_Category=information pathways",
        "Functional_Category=information pathways",
        "Functional_Category=cell wall and cell processes",
        "Product=hypothetical",
        "Product=hypothetical",
    ])]
    return parse_annotation(_gff(tmp_path, rows))


def test_a_spaced_attribute_key_is_read_whole(tmp_path):
    """`Gene Ontology` and `Protein Data Bank` have spaces. A `\\w+` key pattern
    matched `Ontology` and `Bank`, so their values were unreachable."""
    path = _gff(tmp_path, [("Rv0001", "Gene Ontology=GO:0006629,GO:0008152;PFAM=PF00109")])
    gene = parse_annotation(path).gene("Rv0001")
    assert gene.terms["go"] == ("GO:0006629", "GO:0008152")
    assert gene.terms["pfam"] == ("PF00109",)


def test_a_category_containing_a_comma_is_not_split(tmp_path):
    """One of the eleven is `virulence, detoxification, adaptation`. Splitting it
    invents three categories that do not exist."""
    path = _gff(tmp_path, [("Rv0001",
                            "Functional_Category=virulence, detoxification, adaptation")])
    assert parse_annotation(path).gene("Rv0001").terms["category"] == (
        "virulence, detoxification, adaptation",)


def test_the_universe_is_the_annotated_loci_not_every_locus(annotation):
    """A locus carrying no term could not have been a hit, so counting it in the
    denominator treats `unknown` as `not in the term`."""
    universe = universe_for(annotation, "category")
    assert universe.size == 8            # the two uncategorised are out
    assert len(annotation.genes) == 10
    assert universe.terms == 3


def test_a_coherent_set_is_enriched_and_an_arbitrary_one_is_not(annotation):
    universe = universe_for(annotation, "category")

    coherent = enrich(["Rv0000", "Rv0001", "Rv0002"], universe)
    top = coherent.terms[0]
    assert top.term == "lipid metabolism"
    assert (top.m, top.k, top.s, top.n) == (3, 3, 4, 8)
    # Exact rather than thresholded: all 3 drawn land in a 4-of-8 term, which is
    # C(4,3)/C(8,3) = 4/56. A small universe cannot produce a small p, and a test
    # asserting one would be asserting something about the fixture, not the maths.
    assert top.p == pytest.approx(4 / 56)
    assert top.expected == pytest.approx(1.5)
    assert top.ratio == pytest.approx(2.0)

    mixed = enrich(["Rv0000", "Rv0004", "Rv0007"], universe)
    assert mixed.terms == []             # no term carried by two of them


def test_a_locus_with_no_term_is_excluded_from_k_and_said_so(annotation):
    """Keeping it in k would dilute every ratio by a locus that could not have
    contributed to any term."""
    result = enrich(["Rv0000", "Rv0001", "Rv0008", "Rv0009"], universe_for(annotation))
    assert result.k == 2
    assert set(result.dropped) == {"Rv0008", "Rv0009"}
    assert "carry no category term and are excluded from k" in result.note()


def test_a_set_with_no_annotated_locus_is_not_a_negative(annotation):
    result = enrich(["Rv0008", "Rv0009"], universe_for(annotation))
    assert result.terms == [] and result.k == 0
    assert "Not a negative result" in result.note()


def test_the_number_of_terms_tested_travels_with_the_result(annotation):
    """The max statistic is the larger effect: over random k=7 draws the best p
    falls under 0.05 in 52.4% of draws for KEGG pathways and 21.7% for functional
    categories, before any signal. A p without the test count hides that."""
    result = enrich(["Rv0000", "Rv0001", "Rv0004", "Rv0005"], universe_for(annotation))
    assert result.tested == len(result.terms) == 2
    assert "2 category term(s) tested" in result.note()


def test_the_note_never_claims_the_loci_are_related(annotation):
    """The null is a uniform draw. A classifier selected these loci together, so
    they are not a uniform draw and never were."""
    note = enrich(["Rv0000", "Rv0001"], universe_for(annotation)).note()
    assert "measures surprise" in note
    assert "not evidence that the loci are related" in note


# --- the statistics --------------------------------------------------------


def test_fisher_matches_the_hand_computation():
    # 2 of 3 drawn from 10, of which 4 carry the term.
    assert p_at_least(2, 3, 4, 10) == pytest.approx((36 + 4) / 120)
    assert p_at_least(0, 3, 4, 10) == 1.0
    assert p_at_least(5, 3, 4, 10) == 0.0        # more hits than the draw allows


def test_benjamini_hochberg_is_monotone_and_ordered():
    raw = [0.001, 0.008, 0.039, 0.041, 0.9]
    adjusted = benjamini_hochberg(raw)
    assert adjusted == sorted(adjusted)          # monotone in p order
    assert all(a >= p for a, p in zip(adjusted, raw))
    assert adjusted[0] == pytest.approx(0.005)
    assert benjamini_hochberg([]) == []


def test_correction_is_within_a_rule_not_across_the_population():
    """Correcting across every rule in a population is a different and much
    harsher question, and belongs to whoever decides how many they will read."""
    single = benjamini_hochberg([0.01])
    assert single == [pytest.approx(0.01)]


# --- the enrichment gold set -----------------------------------------------


def test_every_gold_set_states_an_expectation_and_a_reason():
    from kegg_string_mcp.rules.gold import load_sets

    axis, note, entries = load_sets()
    assert axis and note
    assert len(entries) >= 5
    for entry in entries:
        assert entry.loci, f"{entry.id} names no loci"
        assert entry.why, f"{entry.id} has no rationale"
        assert (entry.expect_term or entry.expect_no_term_below_q), f"{entry.id} asserts nothing"


def test_the_set_holds_a_base_rate_negative_control():
    """A term present at exactly its expected count must not be called enriched.
    That entry fails if the universe, the expected value or the correction is
    wrong, and no positive control would notice."""
    from kegg_string_mcp.rules.gold import load_sets

    negatives = [e for e in load_sets()[2] if e.kind == "negative"]
    assert negatives
    assert all(e.expect_no_term_below_q for e in negatives)


def test_an_annotation_without_the_axis_skips_rather_than_fails(tmp_path):
    """These are properties of the supplied annotation, not of this code. An
    NCBI GFF carries no functional categories, so nothing was tested and nothing
    failed."""
    from kegg_string_mcp.rules.gold import score_sets

    path = tmp_path / "plain.gff"
    path.write_text("NC_000962.3\tRefSeq\tgene\t1\t99\t.\t+\t.\tlocus_tag=Rv0001\n")
    scores = score_sets(parse_annotation(path))
    assert scores and all(s.skipped for s in scores)
    assert not any(s.failures for s in scores)


def test_a_wrong_expectation_is_reported(annotation, tmp_path):
    """A check that cannot fail is not a check."""
    import json

    from kegg_string_mcp.rules.gold import score_sets

    path = tmp_path / "wrong.json"
    path.write_text(json.dumps({"axis": "category", "sets": [
        {"id": "deliberately-wrong", "kind": "positive", "loci": ["Rv0000", "Rv0001"],
         "expect_term": "information pathways", "expect_q_below": 1.0,
         "why": "these two are lipid metabolism"}]}))
    score = score_sets(annotation, path)[0]
    assert not score.passed
    assert "expected 'information pathways'" in score.failures[0]


def test_the_gold_sets_file_ships_with_the_wheel():
    pyproject = (Path(__file__).resolve().parents[1] / "pyproject.toml").read_text()
    assert '"src/kegg_string_mcp/rules/gold_sets.json"' in pyproject
    assert '"kegg_string_mcp/rules/gold_sets.json"' in pyproject
