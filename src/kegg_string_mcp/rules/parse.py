"""Read a rule-learner's output without deciding anything about it.

A rule is a conjunction of conditions over feature states, predicting a class:

    Rv9012c=1 AND Rv9107=0 AND Rv9233=1    ->  R

Three things this module refuses to throw away, because each is a claim the rule
is making and each disappears under a naive "set of genes" reading:

* **State.** `Rv9107=0` says those isolates carry no qualifying variant at that
  locus. That is an assertion about every isolate the rule covers, not silence
  about the gene, and it is where alternative-route patterns live.
* **Predicted class.** A rule predicting susceptibility is a different claim from
  one predicting resistance and cannot share a verdict template with it.
* **The learner's own statistics.** Numerosity, accuracy, coverage and precision
  are the scan's judgement of the rule. They pass through untouched and prefixed,
  so nothing here can be mistaken for something this package computed.
"""

from __future__ import annotations

import csv
import hashlib
import io
import re
from collections.abc import Iterable
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

# Column names seen in the wild, lowercased. The first match wins.
CONDITION_COLUMNS = ("conditions", "condition", "rule", "antecedent", "if")
CLASS_COLUMNS = ("predicted_class", "class", "consequent", "prediction", "then")
COUNT_COLUMNS = ("n_conditions", "n_cond", "length", "specificity")

_SPLIT = re.compile(r"\s+AND\s+", re.IGNORECASE)
_CONDITION = re.compile(r"^\s*(?P<label>[^=<>!\s]+)\s*=+\s*(?P<state>[01])\s*$")

PASSTHROUGH_PREFIX = "scan_"


@dataclass(frozen=True)
class Condition:
    label: str
    state: int           # 1 = carries a qualifying variant, 0 = matches reference

    def __str__(self) -> str:
        return f"{self.label}={self.state}"


@dataclass(frozen=True)
class Rule:
    """One rule, exactly as the learner emitted it."""

    row: int
    conditions: tuple[Condition, ...]
    predicted_class: str
    extras: dict[str, str] = field(default_factory=dict)
    problems: tuple[str, ...] = ()

    @property
    def k(self) -> int:
        return len(self.conditions)

    @property
    def labels(self) -> tuple[str, ...]:
        return tuple(c.label for c in self.conditions)

    def canonical(self, rename: dict[str, str] | None = None) -> str:
        """Order-independent signature of the rule.

        Conditions are unordered in a conjunction, so two files listing the same
        rule differently must produce the same key. `rename` maps each label to
        its resolved identifier, which is what stops `Rv0006=1` and `Rv0006c=1`
        -- the same locus under two spellings -- counting as two rules.
        """
        rename = rename or {}
        parts = sorted(f"{rename.get(c.label, c.label)}={c.state}" for c in self.conditions)
        return f"{self.predicted_class}|" + "&".join(parts)

    def rule_id(self, rename: dict[str, str] | None = None) -> str:
        return hashlib.sha256(self.canonical(rename).encode()).hexdigest()[:12]

    def to_dict(self, rename: dict[str, str] | None = None) -> dict[str, Any]:
        return {"rule_id": self.rule_id(rename), "row": self.row, "k": self.k,
                "predicted_class": self.predicted_class,
                "conditions": " AND ".join(str(c) for c in self.conditions),
                "problems": "|".join(self.problems)} | self.extras


def _pick(header: list[str], candidates: Iterable[str]) -> str | None:
    lower = {h.strip().lower(): h for h in header}
    for name in candidates:
        if name in lower:
            return lower[name]
    return None


def parse_conditions(text: str) -> tuple[tuple[Condition, ...], tuple[str, ...]]:
    """`A=1 AND B=0` -> conditions, plus whatever could not be read.

    A condition that does not parse is reported rather than dropped: a rule
    silently one condition shorter is still a valid-looking rule, and every
    downstream count would be wrong with nothing to show for it.
    """
    conditions: list[Condition] = []
    problems: list[str] = []
    for chunk in _SPLIT.split(text.strip()):
        if not chunk.strip():
            continue
        match = _CONDITION.match(chunk)
        if match is None:
            problems.append(f"unparsed condition: {chunk.strip()!r}")
            continue
        conditions.append(Condition(match.group("label"), int(match.group("state"))))
    seen: dict[str, int] = {}
    for condition in conditions:
        if condition.label in seen and seen[condition.label] != condition.state:
            problems.append(f"{condition.label} is required to be both 0 and 1")
        seen[condition.label] = condition.state
    return tuple(conditions), tuple(problems)


def parse(path: str | Path, delimiter: str | None = None) -> list[Rule]:
    """Read a rule table. Tab-separated by default; any dialect csv can sniff."""
    path = Path(path)
    text = path.read_text(encoding="utf-8")
    if delimiter is None:
        try:
            delimiter = csv.Sniffer().sniff(text[:4096], delimiters="\t,;").delimiter
        except csv.Error:
            delimiter = "\t"
    reader = csv.DictReader(io.StringIO(text), delimiter=delimiter)
    header = [h for h in (reader.fieldnames or []) if h]
    if not header:
        raise ValueError(f"{path} has no header row")

    condition_column = _pick(header, CONDITION_COLUMNS)
    class_column = _pick(header, CLASS_COLUMNS)
    count_column = _pick(header, COUNT_COLUMNS)
    if condition_column is None:
        raise ValueError(
            f"{path} has no conditions column; looked for {', '.join(CONDITION_COLUMNS)}")

    rules: list[Rule] = []
    for row_number, row in enumerate(reader, start=2):
        conditions, problems = parse_conditions(row.get(condition_column) or "")
        problems = list(problems)

        if count_column:
            # The learner states how many conditions it wrote. Disagreement means
            # the conjunction was split wrongly, and the rule must not be scored
            # as though it had been read correctly.
            declared = (row.get(count_column) or "").strip()
            if declared.isdigit() and int(declared) != len(conditions):
                problems.append(
                    f"{count_column} says {declared} but {len(conditions)} condition(s) parsed")

        extras = {f"{PASSTHROUGH_PREFIX}{k.strip().lower()}": (v or "").strip()
                  for k, v in row.items()
                  if k and k not in {condition_column, class_column}}
        rules.append(Rule(
            row=row_number, conditions=conditions,
            predicted_class=(row.get(class_column) or "").strip() if class_column else "",
            extras=extras, problems=tuple(problems)))
    return rules


def vocabulary(rules: list[Rule]) -> list[str]:
    """Every distinct label the rule set uses, in first-appearance order.

    This, not a separate labels file, is what has to resolve for a run to mean
    anything -- a vocabulary file may list features no rule ever used.
    """
    seen: dict[str, None] = {}
    for rule in rules:
        for condition in rule.conditions:
            seen.setdefault(condition.label, None)
    return list(seen)
