"""Which rules contain which. No network, no annotation.

A learning classifier system nests its rules, so a population read as N
independent findings overstates what it found. These tests cover the two things
that make the containment meaningful rather than superficial: it is on states
rather than loci, and a superset predicting the opposite class is not nesting at
all but an interaction.
"""

from pathlib import Path

from kegg_string_mcp.rules.nesting import nest, summarise
from kegg_string_mcp.rules.parse import parse as parse_rules


def _rules(tmp_path: Path, rows: list[tuple[str, str]]):
    path = tmp_path / "n.tsv"
    path.write_text("conditions\tpredicted_class\n"
                    + "".join(f"{c}\t{k}\n" for c, k in rows))
    return parse_rules(path)


def _by_conditions(rules, nesting, conditions):
    rule = next(r for r in rules if " AND ".join(str(c) for c in r.conditions) == conditions)
    return nesting[rule.rule_id()]


def test_a_smaller_rule_is_contained_in_a_larger_one(tmp_path):
    rules = _rules(tmp_path, [("A=1", "R"), ("A=1 AND B=1", "R"),
                              ("A=1 AND B=1 AND C=1", "R")])
    nesting = nest(rules)
    smallest = _by_conditions(rules, nesting, "A=1")
    largest = _by_conditions(rules, nesting, "A=1 AND B=1 AND C=1")

    assert smallest.is_minimal and not smallest.is_maximal
    assert len(smallest.supersets) == 2
    assert largest.is_maximal and not largest.is_minimal
    assert len(largest.subsets) == 2


def test_containment_is_on_states_not_on_loci(tmp_path):
    """`A=0 AND B=1` is not inside `A=1 AND B=1`: the states conflict, so the two
    describe disjoint sets of isolates, which is the opposite of containment."""
    rules = _rules(tmp_path, [("A=0 AND B=1", "R"), ("A=1 AND B=1", "R"), ("B=1", "R")])
    nesting = nest(rules)

    conflicting = _by_conditions(rules, nesting, "A=0 AND B=1")
    assert conflicting.supersets == []
    # ...while the one that agrees on state is contained in both.
    shared = _by_conditions(rules, nesting, "B=1")
    assert len(shared.supersets) == 2


def test_a_superset_predicting_the_opposite_class_is_reported_apart(tmp_path):
    """If A=1 predicts resistance and A=1 AND B=1 predicts susceptibility, B
    reverses the outcome in A's presence. That is an interaction, not an
    elaboration, and no per-rule reading can see it."""
    rules = _rules(tmp_path, [("A=1", "R"), ("A=1 AND B=1", "S")])
    nesting = nest(rules)
    smaller = _by_conditions(rules, nesting, "A=1")

    assert len(smaller.contradicted_by) == 1
    assert smaller.contradicted_by == smaller.supersets
    assert summarise(nesting)["contradicted_by_a_superset"] == 1


def test_identical_rules_do_not_contain_each_other(tmp_path):
    """Two rules with the same conditions would otherwise each be non-minimal,
    hiding both from a minimal-rules reading."""
    rules = _rules(tmp_path, [("A=1 AND B=1", "R"), ("B=1 AND A=1", "R")])
    nesting = nest(rules)
    assert all(n.is_minimal and n.is_maximal for n in nesting.values())


def test_two_spellings_of_one_locus_are_one_condition(tmp_path):
    """`Rv0006=1` and `Rv0006c=1` are the same locus, so a rule naming one is
    contained in a rule naming the other -- once the rename map says so."""
    rules = _rules(tmp_path, [("Rv0006=1", "R"), ("Rv0006c=1 AND B=1", "R")])
    rename = {"Rv0006": "Rv0006", "Rv0006c": "Rv0006", "B": "B"}

    assert nest(rules)[rules[0].rule_id()].supersets == []
    nested = nest(rules, rename)
    assert len(nested[rules[0].rule_id(rename)].supersets) == 1


def test_the_summary_counts_what_a_reader_would_skip(tmp_path):
    rules = _rules(tmp_path, [("A=1", "R"), ("A=1 AND B=1", "R"),
                              ("A=1 AND B=1 AND C=1", "R"), ("D=1", "R")])
    counts = summarise(nest(rules))
    assert counts == {"rules": 4, "minimal": 2, "elaborations": 2,
                      "contradicted_by_a_superset": 0}


def test_a_rule_with_no_conditions_is_skipped_not_crashed(tmp_path):
    rules = _rules(tmp_path, [("", "R"), ("A=1", "R")])
    nesting = nest(rules)
    assert len(nesting) == 2


def test_the_population_count_is_distinct_rules_not_input_rows(tmp_path):
    """Two rows spelling one rule are one rule, and every figure reported beside
    `minimal` has to use that same denominator. `N of M minimal` with M counting
    rows and N counting identities is a ratio of two different things."""
    path = tmp_path / "dupes.tsv"
    path.write_text(
        "conditions\tpredicted_class\n"
        "Rv0001=1\tR\n"
        "Rv0001=1\tR\n"                       # the same rule again
        "Rv0001=1 AND Rv0002=1\tR\n")
    rules = parse_rules(path)
    nesting = nest(rules)

    assert len(rules) == 3
    assert len({r.rule_id() for r in rules}) == 2
    assert summarise(nesting)["rules"] == 2
    assert summarise(nesting)["minimal"] <= summarise(nesting)["rules"]
