"""Citation validation: does the summary cite only what the tools actually returned?

This is the part most RAG pipelines skip. Retrieval and generation are easy to
wire together; checking that the generated text corresponds to the retrieval is
what makes the output trustworthy rather than merely fluent.

The check is a set-membership test against `RunStore.citable_ids`, which is
populated from tool results *before* the model sees them. Deliberately not a
similarity score: a citation either names a record a tool returned, or it does
not. Fuzzy matching here would reintroduce exactly the ambiguity the store exists
to remove.

Two failure classes are distinguished, because they mean different things:

* **unsupported** -- an identifier that is well-formed and looks authoritative
  but was never returned for this run. This is the hallucination case.
* **cross-target** -- an identifier that WAS retrieved, but for a different gene
  than the sentence attributes it to. Subtler, more plausible-looking, and
  invisible to a validator that only checks global membership.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Any

# KEGG pathway IDs: 3-4 letter organism code + 5 digits (mtu00360).
# STRING protein IDs: NCBI taxon + '.' + locus (83332.Rv1909c).
CITATION_PATTERNS = [
    re.compile(r"\b[a-z]{3,4}\d{5}\b"),
    re.compile(r"\b\d{2,7}\.[A-Za-z0-9_]+\b"),
]


@dataclass
class Citation:
    identifier: str
    status: str            # "verified" | "unsupported" | "cross_target"
    detail: str = ""


@dataclass
class ValidationReport:
    citations: list[Citation] = field(default_factory=list)
    uncited_records: list[str] = field(default_factory=list)

    @property
    def unsupported(self) -> list[Citation]:
        return [c for c in self.citations if c.status == "unsupported"]

    @property
    def cross_target(self) -> list[Citation]:
        return [c for c in self.citations if c.status == "cross_target"]

    @property
    def passed(self) -> bool:
        return not self.unsupported and not self.cross_target

    def summary_line(self) -> str:
        verified = sum(1 for c in self.citations if c.status == "verified")
        return (f"{verified}/{len(self.citations)} citations verified, "
                f"{len(self.unsupported)} unsupported, {len(self.cross_target)} cross-target")

    def to_dict(self) -> dict[str, Any]:
        return {
            "passed": self.passed,
            "summary": self.summary_line(),
            "citations": [vars(c) for c in self.citations],
            "uncited_records": self.uncited_records,
        }


def extract_citations(text: str) -> list[str]:
    """Identifier-shaped tokens, in order of first appearance."""
    found: list[str] = []
    for pattern in CITATION_PATTERNS:
        for match in pattern.findall(text or ""):
            if match not in found:
                found.append(match)
    return found


def validate(
    text: str,
    citable_ids: set[str],
    per_target: dict[str, set[str]] | None = None,
    claimed_target: str | None = None,
) -> ValidationReport:
    """Check every identifier in `text` against what the tools returned.

    `per_target` maps a gene to the IDs retrieved *for that gene*, enabling the
    cross-target check: an ID that exists in the run but was never returned for
    the gene the sentence is about.
    """
    report = ValidationReport()
    cited = extract_citations(text)

    for identifier in cited:
        if identifier not in citable_ids:
            report.citations.append(Citation(
                identifier, "unsupported",
                "not returned by any tool call in this run",
            ))
        elif (per_target and claimed_target
              and identifier not in per_target.get(claimed_target, set())):
            owners = sorted(g for g, ids in per_target.items() if identifier in ids)
            report.citations.append(Citation(
                identifier, "cross_target",
                f"retrieved for {', '.join(owners) or 'another target'}, not for {claimed_target}",
            ))
        else:
            report.citations.append(Citation(identifier, "verified"))

    report.uncited_records = sorted(citable_ids - set(cited))
    return report
