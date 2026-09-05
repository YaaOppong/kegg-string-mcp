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

import difflib
import re
from dataclasses import dataclass, field
from typing import Any

# KEGG pathway IDs: 3-4 letter organism code + 5 digits (mtu00360).
# STRING protein IDs: NCBI taxon + '.' + locus (83332.Rv1909c).
# UniProt's own accession grammar. Precise rather than "an uppercase token", so
# ordinary words and gene symbols in prose are not read as citations.
UNIPROT_ACCESSION = r"[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2}"
KEGG_PATHWAY = r"[a-z]{3,4}\d{5}"
# A lineage-marker record ID. Prefixed rather than a bare H37Rv coordinate: an
# unprefixed 7-digit position is indistinguishable from a PMID, and a citation
# validator that cannot tell a genome coordinate from a paper is worse than none.
TBDB_MARKER = r"tbdb:\d{1,7}"
# A WHO-catalogue variant, e.g. tbdb:katG:p.Ser315Thr. The gene segment starts
# with a letter and the lineage form is all digits, so the two never collide.
# Mutation strings carry HGVS punctuation (p.Ser315Thr, c.-15C>T, n.1401A>G) and
# consequence terms (frameshift_variant), so the tail is deliberately permissive.
TBDB_VARIANT = r"tbdb:[A-Za-z]\w*:[A-Za-z0-9_.>*-]+"
STRING_PROTEIN = r"\d{2,7}\.(?=[A-Za-z0-9_]*[A-Za-z])[A-Za-z0-9_]+"
# PMIDs only in explicit PMID: form. A bare 8-digit number is ambiguous -- it
# could be a coordinate, a score, a year range -- and treating every one as a
# citation would flag ordinary prose as unsupported.
PMID = r"(?i:PMID:?\s*\d{1,8})"

CITATION_PATTERNS = [
    re.compile(r"\b(?:" + KEGG_PATHWAY + r")\b"),
    re.compile(r"\b(?:" + STRING_PROTEIN + r")\b"),
    re.compile(r"\b(?:" + PMID + r")\b"),
    re.compile(r"\b(?:" + UNIPROT_ACCESSION + r")\b"),
    # Variant before marker: both start "tbdb:", and the shorter pattern would
    # otherwise claim the prefix of a variant ID and leave the rest dangling.
    re.compile(r"\b(?:" + TBDB_VARIANT + r")"),
    re.compile(r"\b(?:" + TBDB_MARKER + r")"),
]

# Any identifier that can appear as a citation.
CITE_TOKEN = ("(?:" + PMID + "|" + KEGG_PATHWAY + "|" + STRING_PROTEIN
              + "|" + UNIPROT_ACCESSION + "|" + TBDB_VARIANT
              + "|" + TBDB_MARKER + ")")

# Only the PROSE sources may carry a quote. A KEGG pathway ID and a STRING score
# are structured -- the record means one thing and there is no text to quote from.
# Including them attributed a quote of a tool's own note to a nearby pathway ID,
# which then failed as "no source text": a false positive on correct output, which
# is the one failure a validator must not have.
QUOTABLE_CITE_TOKEN = "(?:" + PMID + "|" + UNIPROT_ACCESSION + ")"

# A quoted span attributed to a record, in either order. The span must be long
# enough that matching it means something; a three-character "quote" would pass
# containment against almost any source text.
MIN_QUOTE_CHARS = 12

# Built by concatenation rather than %-format or f-string: the pattern contains
# {n,} quantifiers, so brace-based formatting would need every brace doubled.
_QUOTE_CHARS = r"\"\u201c\u201d"
_SPAN = r"([^" + _QUOTE_CHARS + r"]{" + str(MIN_QUOTE_CHARS) + r",})"
_OPEN = r"[\"\u201c]"
_CLOSE = r"[\"\u201d]"
# The gap excludes every quote character, opening and closing. Allowing a closing
# curly quote through let this pattern reach across a FINISHED quotation and bind a
# later, unrelated span to a citation that sat inside the earlier one -- in
# `"essential (PMID: 1)" and "it also binds NADH"`, the PMID belongs to the first
# quotation, not the second. An earlier review argued the narrow class suppressed a
# valid match and I widened it; that was wrong, and the test written at the time
# encoded the misattribution as expected behaviour.
_GAP = r"[^\"\u201c\u201d]{0,60}?"

QUOTE_THEN_CITE = re.compile(
    # No global IGNORECASE: UniProt accessions are case-sensitive, and matching
    # them loosely made the quote extractor recognise lowercase accession-shaped
    # tokens the citation extractor never sees -- so a run failed on a token that
    # was never a citation. Only the PMID alternative is case-insensitive.
    _OPEN + _SPAN + _CLOSE + r"\s*[\(\[]?\s*(" + QUOTABLE_CITE_TOKEN + r")",
)
CITE_THEN_QUOTE = re.compile(
    r"(" + QUOTABLE_CITE_TOKEN + r")" + _GAP + _OPEN + _SPAN + _CLOSE,
)


def as_record_id(token: str) -> str:
    """Normalise a cited token to the form `record_id` uses (PMIDs are bare digits)."""
    stripped = token.strip()
    if stripped.upper().startswith("PMID"):
        return stripped.split(":", 1)[-1].strip() if ":" in stripped else stripped[4:].strip()
    return stripped


def normalise(text: str) -> str:
    """Whitespace-normalise and casefold for containment comparison.

    Must match how `pubmed.py` flattened the source text, or a quote that really
    is in the abstract fails on a line break. Case is folded because a model
    legitimately recapitalises a span at the start of a sentence; that is not a
    fabrication, and flagging it would train the reader to ignore quote failures.
    """
    return re.sub(r"\s+", " ", text or "").strip().casefold()


# Punctuation a quoter legitimately adds or drops at the edges of a span. Markdown
# emphasis is included: models routinely italicise or bold a quotation, and those
# markers land inside the captured span. Treating "*quoted text*" as absent from a
# source containing "quoted text" is a false positive on correct output.
_EDGE_PUNCT = " \t\n.,;:!?\"'\u201c\u201d\u2018\u2019()[]*_`"
_ELLIPSIS = re.compile(r"\s*(?:\.\.\.|\u2026)\s*")


def _fragments(quote: str) -> list[str]:
    """Split a quote on ellipses. A model eliding a passage writes "A ... B", and
    requiring that to appear contiguously would fail an honest quotation."""
    return [f for f in (part.strip(_EDGE_PUNCT) for part in _ELLIPSIS.split(quote)) if f]


def quote_in_source(quote: str, source: str) -> bool:
    """Containment, tolerant of the edits honest quoting actually involves.

    Terminal punctuation is stripped: truncating mid-sentence and closing with a
    full stop is normal quoting, not fabrication. Observed live -- a model quoted
    "...mutations in ahpC or inhA." from a sentence that continues "and between
    mutations in kasA", and strict containment called that a fabricated quote.

    Fragments split on an ellipsis must each appear, and in order, so an elision
    cannot silently join two unrelated passages.
    """
    haystack = normalise(source)
    position = 0
    for fragment in _fragments(quote):
        found = haystack.find(normalise(fragment), position)
        if found < 0:
            return False
        position = found + len(normalise(fragment))
    return True


# Above this similarity, a failed quote is near-certainly a quoting artefact --
# a dropped word, altered punctuation, a British/American spelling -- rather than
# an invented claim. It does NOT change the verdict: the check stays deterministic
# and binary, and a model must never be able to argue its way past it. It ranks
# failures so a human looks at the fabrications first.
LIKELY_ARTEFACT = 0.90


@dataclass
class QuoteCheck:
    record_id: str
    quote: str
    status: str        # "verified" | "not_in_source" | "no_source_text"
    detail: str = ""
    similarity: float | None = None     # best match found in the source, 0-1
    closest_span: str = ""              # the passage it most resembles
    triage: str = ""                    # "likely_quoting_artefact" | "likely_fabricated"


@dataclass
class Citation:
    identifier: str
    status: str            # "verified" | "unsupported" | "cross_target"
    detail: str = ""


@dataclass
class ValidationReport:
    produced_output: bool = True
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
        # An empty summary has nothing to contradict, so every check trivially
        # passed -- and a run that produced nothing was reported as a success.
        # Absence of failure is not success.
        return (self.produced_output and not self.unsupported
                and not self.cross_target and not self.bad_quotes)

    def summary_line(self) -> str:
        if not self.produced_output:
            return "no summary was produced, so nothing could be validated"
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
            token = as_record_id(match if isinstance(match, str) else match[0])
            if token and token not in found:
                found.append(token)
    return found


def _quotation_ranges(text: str) -> list[tuple[int, int]]:
    """Character ranges enclosed by quotation marks, paired left to right.

    Straight quotes carry no direction, so a regex cannot tell an opening mark
    from a closing one; pairing them in order can. This is what stops a citation
    that sits INSIDE a quotation from being read as introducing the NEXT one.
    """
    marks = [i for i, c in enumerate(text or "") if c in "\"\u201c\u201d"]
    return [(a, b) for a, b in zip(marks[::2], marks[1::2], strict=False)]


def extract_quotes(text: str) -> list[tuple[str, str]]:
    """(record_id, quoted span) pairs, in either written order.

    A span already claimed by a *trailing* citation is never re-attributed to a
    different PMID. Two quotes in a row, each closing with its own citation, is
    the natural way to write this -- and the leading-citation pattern would reach
    forward across the sentence boundary and bind the second quote to the first
    PMID as well, failing a correctly-cited summary.

    A leading citation that sits INSIDE a quotation is ignored for the same
    reason. Quoting a tool note that happens to name its own record --
    `"UniProt holds no FUNCTION statement for P71814 -- ..."` -- put an accession
    inside the quotation with no citation after the closing mark. The
    leading-citation pattern then ran from that accession, over the closing mark
    it could not recognise as closing, and bound the model's own following prose
    as a quotation attributed to the record. It failed as "not in source": a
    false accusation of fabrication against correct output, which is the one
    failure a validator must not have.
    """
    pairs: list[tuple[str, str]] = []
    claimed: set[str] = set()
    for span, token in QUOTE_THEN_CITE.findall(text or ""):
        pairs.append((as_record_id(token), span))
        claimed.add(span)
    inside = _quotation_ranges(text)
    for match in CITE_THEN_QUOTE.finditer(text or ""):
        token, span = match.group(1), match.group(2)
        start = match.start(1)
        if any(open_at < start < close_at for open_at, close_at in inside):
            continue
        record_id = as_record_id(token)
        if span not in claimed and (record_id, span) not in pairs:
            pairs.append((record_id, span))
            claimed.add(span)
    return pairs


def nearest_span(quote: str, source: str) -> tuple[float, str]:
    """Best-matching window in the source, and how close it is.

    Purely diagnostic. Adjudicating a failed quote otherwise means reading the
    abstract by hand, and both false positives this validator has produced would
    have been obvious in one line of output: a 0.98 match differing only in a
    trailing full stop is a quoting artefact, and a 0.31 match is a fabrication.
    """
    needle, haystack = normalise(quote), normalise(source)
    if not needle or not haystack:
        return 0.0, ""
    window = len(needle)
    best_ratio, best_span = 0.0, ""
    # Step by a fraction of the window so a match straddling a boundary is still
    # found, without comparing every offset in a long abstract.
    step = max(1, window // 4)
    for start in range(0, max(1, len(haystack) - window + 1) + step, step):
        candidate = haystack[start:start + window]
        if not candidate:
            break
        ratio = difflib.SequenceMatcher(None, needle, candidate).ratio()
        if ratio > best_ratio:
            best_ratio, best_span = ratio, candidate
    return round(best_ratio, 3), best_span


def check_quotes(text: str, records: dict[str, dict[str, Any]],
                 notes: dict[str, list[str]] | None = None) -> list[QuoteCheck]:
    """Verify each quoted span really appears in the record it is attributed to.

    This is the step that upgrades validation from "the model cited a record that
    was retrieved" to "the model cited a record that was retrieved AND the words
    it put in quotation marks are actually in it". Containment on normalised
    text -- deterministic, no embeddings, no similarity threshold.

    It does NOT verify that the surrounding sentence characterises the quote
    fairly. A correctly quoted span can still be framed misleadingly; that is a
    harder problem and this check does not claim to solve it.
    """
    notes = notes or {}
    checks: list[QuoteCheck] = []
    for record_id, quote in extract_quotes(text):
        record = records.get(record_id)
        source = (record or {}).get("detail", {}).get("quotable_text", "")
        # A tool NOTE is source text too. UniProt's "No function statement for
        # A0A0N7EHL5 ..." was quoted verbatim and reported as likely fabricated,
        # because only record quotable_text was searched -- a false accusation
        # against correct output. Only the notes from the result that returned
        # THIS record are searched, so a quote cannot be verified against a note
        # about something else.
        if any(quote_in_source(quote, note) for note in notes.get(record_id, [])):
            checks.append(QuoteCheck(record_id, quote, "verified"))
        elif not source:
            checks.append(QuoteCheck(record_id, quote, "no_source_text",
                                     "no retrieved text for this record to check the quote against"))
        elif quote_in_source(quote, source):
            checks.append(QuoteCheck(record_id, quote, "verified"))
        else:
            ratio, span = nearest_span(quote, source)
            artefact = ratio >= LIKELY_ARTEFACT
            checks.append(QuoteCheck(
                record_id, quote, "not_in_source",
                "this span does not appear in the retrieved title or abstract",
                similarity=ratio, closest_span=span,
                triage="likely_quoting_artefact" if artefact else "likely_fabricated",
            ))
    return checks


def validate(
    text: str,
    citable_ids: set[str],
    per_target: dict[str, set[str]] | None = None,
    claimed_target: str | None = None,
    records: dict[str, dict[str, Any]] | None = None,
    notes: dict[str, list[str]] | None = None,
) -> ValidationReport:
    """Check every identifier in `text` against what the tools returned.

    `per_target` maps a gene to the IDs retrieved *for that gene*, enabling the
    cross-target check: an ID that exists in the run but was never returned for
    the gene the sentence is about.
    """
    report = ValidationReport(produced_output=bool((text or "").strip()))
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
        report.quotes = check_quotes(text, records, notes)

    report.uncited_records = sorted(citable_ids - set(cited))
    return report
