"""Drug-resistance variant lookups: is this gene resistance-associated, and how?

A stage 1 annotation alongside KEGG, STRING, UniProt and the lineage barcode.
The source is TB-Profiler's `tbdb/mutations.csv`, the machine-readable form of
the WHO catalogue of mutations: 49,330 graded variants across 74 genes and 18
drugs, public, versioned in git, no credentials.

**The gene-level rule is deliberately asymmetric.** A gene is marked
resistance-associated if ANY of its variants is graded as associated with
resistance, however many are not. katG lists 1,771 variants of which only 139
carry an association grade -- it is still a resistance gene, and the 1,254
graded "Uncertain significance" do not dilute that. Absence of association is
not evidence of absence at the gene level.

**The grading is the annotation, not a detail of it.** Collapsing this file to
"is the gene present" would report katG and bacA identically, because presence
mostly means the gene was examined. 70% of all rows are "Uncertain
significance". Every record therefore carries its WHO grade, and the two
association tiers are kept apart:

    Assoc w R            15 genes   the confident grade
    Assoc w R - Interim  +9 genes   associated, lower confidence

**This tool takes loci, never variants.** The flag says the gene has known
resistance-associated variants, not that any particular variant in it does.
Within katG, p.Ser315Thr is "Assoc w R" while p.Arg463Leu -- a common
polymorphism -- is explicitly "Not assoc w R", and the gene-level answer cannot
tell them apart. Reading the flag as a verdict on a variant would mark a known
benign polymorphism as resistance-conferring, well-sourced and wrong. The graded
variants are returned as records so the finer detail is visible and citable, but
grading a specific variant is a question this tool does not answer.

An explicit "Not assoc w R" is a finding rather than silence, and rarer than
either: the catalogue graded that variant and found against it. Counts for every
grade are returned so a caller can tell the two apart.
"""

from __future__ import annotations

import csv
import io
from dataclasses import asdict, dataclass
from typing import Any

from kegg_string_mcp.http import FetchError, PoliteClient
from kegg_string_mcp.provenance import Record, RequestTrace, ToolResult

MUTATIONS_URL = "https://raw.githubusercontent.com/jodyphelan/tbdb/master/mutations.csv"
# No per-variant landing page exists, so records resolve to the catalogue itself.
CATALOGUE_SOURCE = "https://github.com/jodyphelan/tbdb"

# WHO grades that assert an association. Interim is a weaker grade, not a
# different kind of claim, so both mark the gene -- and both are reported
# separately so a caller can require the confident one.
ASSOCIATED = ("Assoc w R", "Assoc w R - Interim")
UNGRADED = "(ungraded)"


@dataclass(frozen=True)
class Variant:
    gene: str
    mutation: str
    drug: str
    confidence: str
    source: str
    comment: str = ""

    @property
    def associated(self) -> bool:
        return self.confidence in ASSOCIATED

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def parse_catalogue(text: str) -> dict[str, list[Variant]]:
    """Gene -> its graded variants.

    Parsed with the csv module rather than by splitting on commas: the comment
    column is free text and quotes embedded commas, so a naive split silently
    shifts every field after it.
    """
    out: dict[str, list[Variant]] = {}
    for row in csv.DictReader(io.StringIO(text)):
        gene = (row.get("Gene") or "").strip()
        if not gene:
            continue
        out.setdefault(gene, []).append(Variant(
            gene=gene,
            mutation=(row.get("Mutation") or "").strip(),
            drug=(row.get("drug") or "").strip(),
            # An empty confidence means the row was never graded, which is a
            # third state -- not association and not non-association.
            confidence=(row.get("confidence") or "").strip() or UNGRADED,
            source=(row.get("source") or "").strip(),
            comment=(row.get("comment") or "").strip()))
    return out


def grading_counts(variants: list[Variant]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for variant in variants:
        counts[variant.confidence] = counts.get(variant.confidence, 0) + 1
    return dict(sorted(counts.items()))


def _trace(resp: Any) -> RequestTrace:
    return RequestTrace(url=resp.audit_url, retrieved_at=resp.fetched_at, cached=resp.cached,
                        status=resp.status, content_sha256=resp.content_sha256)


class ResistanceClient:
    """Gene -> WHO-graded resistance variants.

    The catalogue is fetched once per client and reused: it is 4.7 MB, and
    re-parsing 49,330 rows per gene would dominate every lookup.
    """

    def __init__(self, http: PoliteClient):
        self.http = http
        self._catalogue: dict[str, list[Variant]] | None = None
        self._traces: list[RequestTrace] = []

    def _load(self) -> tuple[dict[str, list[Variant]], list[RequestTrace]]:
        if self._catalogue is None:
            response = self.http.get(MUTATIONS_URL)
            self._catalogue = parse_catalogue(response.body)
            self._traces = [_trace(response)]
        return self._catalogue, self._traces

    def variants(self, gene: str, drug: str | None = None) -> ToolResult:
        query: dict[str, Any] = {"gene": gene, "drug": drug}
        gene = gene.strip()
        if not gene:
            return ToolResult.build(query, [], resolved={"matched_by": "none"},
                                    notes=["No gene identifier was supplied."])
        try:
            catalogue, traces = self._load()
        except FetchError as exc:
            return ToolResult.build(
                query, [], resolved={"matched_by": "none"},
                notes=[(f"Could not fetch the resistance catalogue: HTTP {exc.status}. No "
                        f"lookup was made -- this is not evidence that {gene} lacks "
                        f"resistance-associated variants.")])

        found = catalogue.get(gene)
        matched_by = "symbol"
        if found is None:
            for key, value in catalogue.items():
                if key.lower() == gene.lower():
                    found, gene = value, key
                    break
        if found is None:
            return ToolResult.build(
                query, [], resolved={"matched_by": "none", "resistance_associated": False},
                requests=traces,
                notes=[(f"'{gene}' is not in the WHO catalogue. The catalogue covers 74 genes "
                        f"selected for drug-resistance surveillance, so absence means the gene "
                        f"was not assessed -- it is not a finding that the gene is unrelated "
                        f"to resistance.")])

        # Records are the association-graded variants, optionally narrowed by
        # drug. Returning every row would bury the answer: katG lists 1,771, of
        # which 1,254 are graded "Uncertain significance".
        selected = [v for v in found if v.associated]
        if drug:
            selected = [v for v in selected if v.drug.lower() == drug.lower()]

        # The flag is computed over EVERY variant of the gene, never over the
        # filtered subset: one associated variant marks the gene however many are
        # not, so narrowing by drug must not unmark it.
        associated = [v for v in found if v.associated]
        resolved = {
            "matched_by": matched_by,
            "resistance_associated": bool(associated),
            "drugs": sorted({v.drug for v in associated if v.drug}),
            "variants_in_catalogue": len(found),
            "grading_counts": grading_counts(found),
        }

        records_for = selected
        records = [
            Record(
                record_id=f"tbdb:{v.gene}:{v.mutation}",
                type="resistance_variant",
                name=f"{v.mutation} -- {v.confidence} ({v.drug})",
                url=CATALOGUE_SOURCE,
                source="tbdb",
                retrieved_at=traces[0].retrieved_at,
                cached=traces[0].cached,
                detail=v.to_dict() | {"associated": v.associated},
            )
            for v in records_for
        ]

        return ToolResult.build(query, records, resolved=resolved, requests=traces,
                                notes=[_gene_note(gene, found, associated)])


def _gene_note(gene: str, found: list[Variant], associated: list[Variant]) -> str:
    counts = grading_counts(found)
    if associated:
        tiers = ", ".join(f"{c} {t}" for t, c in counts.items() if t in ASSOCIATED)
        return (f"{gene} is RESISTANCE-ASSOCIATED: {len(associated)} of its {len(found)} "
                f"catalogued variants are graded as associated with resistance ({tiers}). "
                f"One associated variant marks the gene however many are not. This is a "
                f"statement about the GENE: it does not mean a particular variant in {gene} "
                f"confers resistance, and the catalogue explicitly grades some variants in "
                f"resistance-associated genes as not associated.")
    return (f"{gene} is in the WHO catalogue with {len(found)} graded variants, none of them "
            f"associated with resistance ({', '.join(f'{c} {t}' for t, c in counts.items())}). "
            f"The gene was assessed and no variant met the association grade.")
