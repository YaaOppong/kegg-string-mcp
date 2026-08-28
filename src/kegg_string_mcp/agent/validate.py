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
    # PMIDs only in explicit PMID: form. A bare 8-digit number is ambiguous --
    # it could be a coordinate, a score, a year range -- and treating every one
    # as a citation would flag ordinary prose as unsupported.
    re.compile(r"\bPMID:?\s*(\d{1,8})\b", re.IGNORECASE),
]

# A quoted span attributed to a PMID, in either order. The span must be long
# enough that matching it means something; a three-character "quote" would pass
# containment against almost any abstract.
MIN_QUOTE_CHARS = 12
QUOTE_THEN_PMID = re.compile(
    r"[\"\u201c]([^\"\u201c\u201d]{%d,})[\"\u201d]\s*[\(\[]?\s*PMID:?\s*(\d{1,8})" % MIN_QUOTE_CHARS,
    re.IGNORECASE,
)
PMID_THEN_QUOTE = re.compile(
    r"PMID:?\s*(\d{1,8})[^\"\u201c]{0,60}?[\"\u201c]([^\"\u201c\u201d]{%d,})[\"\u201d]" % MIN_QUOTE_CHARS,
    re.IGNORECASE,
)


def normalise(text: str) -> str:
    """Whitespace-normalise and casefold for containment comparison.

    Must match how `pubmed.py` flattened the source text, or a quote that really
    is in the abstract fails on a line break. Case is folded because a model
    legitimately recapitalises a span at the start of a sentence; that is not a
    fabrication, and flagging it would train the reader to ignore quote failures.
    """
    return re.sub(r"\s+", " ", text or "").strip().casefold()


@dataclass
class QuoteCheck:
    record_id: str
    quote: str
    status: str        # "verified" | "not_in_source" | "no_source_text"
    detail: str = ""


@dataclass
class Citation:
    identifier: str
    status: str            # "verified" | "unsupported" | "cross_target"
    detail: str = ""


@dataclass
class ValidationReport:
    citations: list[Citation] = field(default_factory=list)
    quotes: list[QuoteCheck] = field(default_factory=list)
    uncited_records: list[str] = field(default_factory=list)

    @property
    def unsupported(self) -> list[Citation]:
        return [c for c in self.citations if c.status == "unsupported"]

    @property
    def cross_target(self) -> list[Citation]:
        return [c for c in self.citations if c.status == "cross_target"]

    @property
    def bad_quotes(self) -> list[QuoteCheck]:
        return [q for q in self.quotes if q.status != "verified"]

    @property
    def passed(self) -> bool:
        return not self.unsupported and not self.cross_target and not self.bad_quotes

    def summary_line(self) -> str:
        verified = sum(1 for c in self.citations if c.status == "verified")
        line = (f"{verified}/{len(self.citations)} citations verified, "
                f"{len(self.unsupported)} unsupported, {len(self.cross_target)} cross-target")
        if self.quotes:
            good = sum(1 for q in self.quotes if q.status == "verified")
            line += f"; {good}/{len(self.quotes)} quotes found in source"
        return line

    def to_dict(self) -> dict[str, Any]:
        return {
            "passed": self.passed,
            "summary": self.summary_line(),
            "citations": [vars(c) for c in self.citations],
            "quotes": [vars(q) for q in self.quotes],
            "uncited_records": self.uncited_records,
        }


def extract_citations(text: str) -> list[str]:
    """Identifier-shaped tokens, in order of first appearance.

    PMIDs are returned bare (the form `record_id` uses) even though they are only
    matched in `PMID:` form, so membership is a plain set test.
    """
    found: list[str] = []
    for pattern in CITATION_PATTERNS:
        for match in pattern.findall(text or ""):
            token = match if isinstance(match, str) else match[0]
            if token not in found:
                found.append(token)
    return found


def extract_quotes(text: str) -> list[tuple[str, str]]:
    """(record_id, quoted span) pairs, in either written order."""
    pairs: list[tuple[str, str]] = []
    for span, pmid in QUOTE_THEN_PMID.findall(text or ""):
        pairs.append((pmid, span))
    for pmid, span in PMID_THEN_QUOTE.findall(text or ""):
        if (pmid, span) not in pairs:
            pairs.append((pmid, span))
    return pairs


def check_quotes(text: str, records: dict[str, dict[str, Any]]) -> list[QuoteCheck]:
    """Verify each quoted span really appears in the record it is attributed to.

    This is the step that upgrades validation from "the model cited a record that
    was retrieved" to "the model cited a record that was retrieved AND the words
    it put in quotation marks are actually in it". Containment on normalised
    text -- deterministic, no embeddings, no similarity threshold.

    It does NOT verify that the surrounding sentence characterises the quote
    fairly. A correctly quoted span can still be framed misleadingly; that is a
    harder problem and this check does not claim to solve it.
    """
    checks: list[QuoteCheck] = []
    for record_id, quote in extract_quotes(text):
        record = records.get(record_id)
        source = (record or {}).get("detail", {}).get("quotable_text", "")
        if not source:
            checks.append(QuoteCheck(record_id, quote, "no_source_text",
                                     "no retrieved text for this record to check the quote against"))
        elif normalise(quote) in normalise(source):
            checks.append(QuoteCheck(record_id, quote, "verified"))
        else:
            checks.append(QuoteCheck(record_id, quote, "not_in_source",
                                     "this span does not appear in the retrieved title or abstract"))
    return checks


def validate(
    text: str,
    citable_ids: set[str],
    per_target: dict[str, set[str]] | None = None,
    claimed_target: str | None = None,
    records: dict[str, dict[str, Any]] | None = None,
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

    if records:
        report.quotes = check_quotes(text, records)

    report.uncited_records = sorted(citable_ids - set(cited))
    return report
