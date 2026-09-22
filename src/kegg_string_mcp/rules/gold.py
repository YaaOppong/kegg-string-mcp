"""Score the classifier against rules whose answer is already known.

The classification is deterministic, which makes it checkable -- and something
checkable that is never checked is just an assertion. The gold set is small on
purpose: a dozen rules drawn from canonical M. tuberculosis pharmacology, each
with the reason it must classify the way it does written next to it.

Two kinds of entry, and the second is what keeps the first honest:

* **positive** -- must reach a named primary signature. Tests that the classifier
  finds what is there.
* **negative** -- must NOT reach a named signature. Tests that it declines what
  is not. `katG + rpoA` is the one that matters: it has the compensation shape
  and STRING will supply a co-mention edge for it, so only the drug-concordance
  check stops it, and nothing else in the suite would notice if that check broke.

A failure is not automatically a bug. The WHO catalogue is versioned, and a locus
regraded upstream changes what the classifier is able to say. The report prints
what it got and why the expectation existed, so the reader can tell a broken rule
from a moved source.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from kegg_string_mcp.rules.evidence import Sources, evidence_for
from kegg_string_mcp.rules.parse import Rule, parse_conditions
from kegg_string_mcp.rules.signature import Link, RuleSignature, classify

GOLD = Path(__file__).with_name("gold_rules.json")


@dataclass
class GoldRule:
    id: str
    kind: str
    conditions: str
    predicted_class: str
    why: str = ""
    expect_primary: str = ""
    expect_contains: list[str] = field(default_factory=list)
    expect_absent: list[str] = field(default_factory=list)
    expect_cross_resistant: list[str] = field(default_factory=list)

    def as_rule(self, row: int) -> Rule:
        conditions, problems = parse_conditions(self.conditions)
        return Rule(row=row, conditions=conditions,
                    predicted_class=self.predicted_class, problems=problems)


@dataclass
class GoldSet:
    organism: str
    reference: str
    validated_on: str
    note: str
    rules: list[GoldRule]


def load(path: str | Path | None = None) -> GoldSet:
    data = json.loads(Path(path or GOLD).read_text(encoding="utf-8"))
    return GoldSet(organism=data.get("organism", "mtu"), reference=data.get("reference", ""),
                   validated_on=data.get("validated_on", ""), note=data.get("note", ""),
                   rules=[GoldRule(**entry) for entry in data["rules"]])


@dataclass
class Score:
    gold: GoldRule
    signature: RuleSignature | None
    failures: list[str] = field(default_factory=list)
    unresolved: list[str] = field(default_factory=list)

    @property
    def passed(self) -> bool:
        return not self.failures and not self.unresolved

    @property
    def got(self) -> str:
        return self.signature.primary if self.signature else "-"


def score_one(gold: GoldRule, row: int, sources: Sources,
              links: dict[tuple[str, str], Link] | None = None) -> Score:
    rule = gold.as_rule(row)
    conditions = evidence_for(rule, sources)
    signature = classify(rule, conditions, links)
    result = Score(gold=gold, signature=signature)

    # A label the supplied annotation cannot resolve is reported separately from
    # a wrong answer: the classifier was never given the chance to be right, and
    # scoring it as a miss would blame the wrong component.
    result.unresolved = [c.label for c in conditions if not c.feature.resolved]
    if result.unresolved:
        return result

    if gold.expect_primary and signature.primary != gold.expect_primary:
        result.failures.append(f"primary {signature.primary!r}, expected {gold.expect_primary!r}")
    for expected in gold.expect_contains:
        if expected not in signature.signatures:
            result.failures.append(f"missing signature {expected!r}")
    for forbidden in gold.expect_absent:
        if forbidden in signature.signatures:
            result.failures.append(f"must not carry {forbidden!r}")
    if gold.expect_cross_resistant:
        got = {locus for locus, _ in signature.cross_resistant}
        for locus in gold.expect_cross_resistant:
            if locus not in got:
                result.failures.append(f"{locus} not reported cross-resistant")
    return result


def score(sources: Sources, gold: GoldSet | None = None,
          links: dict[tuple[str, str], Link] | None = None) -> list[Score]:
    gold = gold or load()
    return [score_one(entry, row, sources, links)
            for row, entry in enumerate(gold.rules, start=1)]


def render(scores: list[Score], gold: GoldSet) -> str:
    passed = sum(1 for s in scores if s.passed)
    skipped = [s for s in scores if s.unresolved]
    lines = [
        f"Reference: {gold.reference}",
        f"Validated: {gold.validated_on}",
        "",
        f"{'rule':38} {'expected':26} {'got':26} ",
    ]
    for result in scores:
        if result.unresolved:
            state = f"SKIPPED unresolved: {', '.join(result.unresolved)}"
        elif result.passed:
            state = "ok"
        else:
            state = "FAIL " + "; ".join(result.failures)
        lines.append(f"{result.gold.id:38} {result.gold.expect_primary or '-':26} "
                     f"{result.got:26} {state}")

    scored = len(scores) - len(skipped)
    lines += ["", f"passed {passed}/{scored} scored"
                  + (f", {len(skipped)} skipped for unresolved loci" if skipped else "")]
    failures = [s for s in scores if s.failures]
    if failures:
        lines += ["", "why each expectation exists:"]
        for result in failures:
            lines.append(f"  {result.gold.id}: {result.gold.why}")
        lines += ["", gold.note]
    return "\n".join(lines)


def summary(scores: list[Score]) -> dict[str, Any]:
    scored = [s for s in scores if not s.unresolved]
    return {"total": len(scores), "scored": len(scored),
            "passed": sum(1 for s in scored if s.passed),
            "failed": sum(1 for s in scored if s.failures),
            "skipped": len(scores) - len(scored)}


# --- the enrichment gold set -----------------------------------------------
#
# The rule gold set scores the classifier. This one scores the enrichment, and
# it is a different kind of check: these are properties of the supplied
# annotation, not of this code, so an annotation without the axis skips rather
# than fails, and a changed annotation may legitimately change an expected term.

GOLD_SETS = Path(__file__).with_name("gold_sets.json")


@dataclass
class GoldLocusSet:
    id: str
    kind: str
    loci: list[str]
    why: str = ""
    expect_term: str = ""
    expect_q_below: float = 1.0
    expect_no_term_below_q: float = 0.0


@dataclass
class SetScore:
    gold: GoldLocusSet
    top: Any = None
    failures: list[str] = field(default_factory=list)
    skipped: str = ""

    @property
    def passed(self) -> bool:
        return not self.failures and not self.skipped


def load_sets(path: str | Path | None = None) -> tuple[str, str, list[GoldLocusSet]]:
    data = json.loads(Path(path or GOLD_SETS).read_text(encoding="utf-8"))
    return (data.get("axis", "category"), data.get("note", ""),
            [GoldLocusSet(**entry) for entry in data["sets"]])


def score_sets(annotation: Any, path: str | Path | None = None) -> list[SetScore]:
    from kegg_string_mcp.rules.enrichment import enrich, universe_for

    axis, _, entries = load_sets(path)
    universe = universe_for(annotation, axis)
    out: list[SetScore] = []
    for entry in entries:
        score = SetScore(gold=entry)
        if not universe.size:
            # The annotation carries no terms on this axis -- an NCBI GFF, say.
            # Nothing was tested, so nothing failed.
            score.skipped = f"the supplied annotation carries no {axis} terms"
            out.append(score)
            continue

        result = enrich(entry.loci, universe)
        missing = [locus for locus in entry.loci if universe.members and locus in entry.loci
                   and locus not in universe.annotated]
        score.top = result.terms[0] if result.terms else None

        if entry.kind == "negative":
            called = result.significant(entry.expect_no_term_below_q)
            if called:
                score.failures.append(
                    "called " + ", ".join(f"{t.term} (q={t.q:.3g})" for t in called)
                    + f" at q<={entry.expect_no_term_below_q}")
        else:
            if score.top is None:
                score.failures.append("no term reached two members")
            else:
                if score.top.term != entry.expect_term:
                    score.failures.append(
                        f"top term {score.top.term!r}, expected {entry.expect_term!r}")
                if score.top.q > entry.expect_q_below:
                    score.failures.append(
                        f"q={score.top.q:.3g}, expected below {entry.expect_q_below:g}")
        if missing and not score.failures:
            score.failures.append(f"{len(missing)} locus/loci absent from the annotation")
        out.append(score)
    return out


def render_sets(scores: list[SetScore], note: str = "") -> str:
    lines = [f"{'set':22} {'expected term':42} {'got':42} "]
    for score in scores:
        if score.skipped:
            state, got = f"SKIPPED {score.skipped}", "-"
        elif score.passed:
            state = f"ok  q={score.top.q:.2g}" if score.top else "ok"
            got = score.top.term if score.top else "-"
        else:
            state, got = "FAIL " + "; ".join(score.failures), (
                score.top.term if score.top else "-")
        expected = score.gold.expect_term or f"nothing below q={score.gold.expect_no_term_below_q}"
        lines.append(f"{score.gold.id:22} {expected:42} {got:42} {state}")
    scored = [s for s in scores if not s.skipped]
    lines += ["", f"passed {sum(1 for s in scored if s.passed)}/{len(scored)} scored"
                  + (f", {len(scores) - len(scored)} skipped" if len(scored) != len(scores) else "")]
    failures = [s for s in scores if s.failures]
    if failures:
        lines += ["", "why each expectation exists:"]
        lines += [f"  {s.gold.id}: {s.gold.why}" for s in failures]
        if note:
            lines += ["", note]
    return "\n".join(lines)
