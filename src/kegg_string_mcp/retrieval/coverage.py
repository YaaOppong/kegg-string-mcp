"""Which genes need literature retrieval, and why.

Stage 2 exists to cover what stage 1 cannot. Running it over every gene wastes
the effort on genes KEGG and UniProt already describe well, and -- worse -- hides
the cases it was built for among them. This module makes the routing decision
explicit and auditable: each gene gets a coverage record naming the gaps that
qualify it, so a later reader can see why an annotation rested on literature
rather than on a curated source.

Three gaps qualify, in descending severity:

* **no_function** -- UniProt returns no function statement at all. Nothing
  structured describes what the protein does.
* **no_experimental_function** -- UniProt describes it, but every statement is
  inferred from a rule or a homologue. That is not evidence about *this*
  protein, which is the distinction `uniprot.py` already draws and the same one
  STRING's textmining channel forces elsewhere in this codebase.
* **no_pathway** -- KEGG assigns no pathway. KEGG covers roughly 1,170 of 4,008
  M. tuberculosis genes, so this is common and is the weakest of the three on its
  own: gyrA has no KEGG pathway and is one of the best-characterised genes in TB.

`no_pathway` alone is therefore a weak signal, and a gene qualifying on it alone
is reported separately from one with a genuine functional gap.

**An absent answer is not an absent annotation.** A lookup can fail to describe a
gene in three ways, and only one of them is a finding:

* the request failed -- no request trace was recorded;
* the identifier did not resolve -- the request succeeded but matched nothing,
  which is a fact about the symbol, not the gene. UniProt has no entry for
  `icl1`; the same protein under `icl` or `Rv0467` is annotated;
* the source resolved the gene and holds no such annotation -- the finding.

Only the third routes anything. The first two are reported as an unknown, so a
gene never reaches literature retrieval because its symbol was wrong, and never
records "nothing is known about this gene" on the strength of a failed match.

**The two sources disagree about symbols**, which is the ordinary state of gene
nomenclature rather than a defect. KEGG calls Rv3133c `devR` and Rv3919c `gid`;
UniProt calls the first `devR` too but has no entry for `icl1`, which it holds as
`icl`. A single symbol therefore resolves in one source and not the other. Both
carry the locus tag on records they do return -- UniProt as `locus_tags`, KEGG as
`kegg_gene_id` -- so whichever source answered supplies the identifier to retry
the other with, and `resolved_via` records that it happened.

**The retry never guesses which record to take the tag from.** A symbol search
can match several entries, and taking the first would let an unreviewed fragment
or a paralogue redirect the lookup to a different gene. Where the candidates
disagree, the reviewed entry settles it -- SwissProt versus TrEMBL is exactly
that distinction, and it is already parsed. Where they still disagree, the gene
is flagged `needs_review` with the competing tags in `locus_candidates` and
routed nowhere, because a coverage verdict about the wrong gene is worse than no
verdict. In 17 genes checked, including paralogue-prone families, no symbol
returned two tagged records; the flag exists so that the day one does, it is a
line in the output rather than a silent substitution.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any

NO_PATHWAY = "no_pathway"
NO_FUNCTION = "no_function"
NO_EXPERIMENTAL = "no_experimental_function"

# Gaps that mean nothing structured describes the protein's function. A missing
# KEGG pathway on its own does not: KEGG's coverage is a property of KEGG.
FUNCTIONAL_GAPS = frozenset({NO_FUNCTION, NO_EXPERIMENTAL})


@dataclass
class Coverage:
    gene: str
    kegg_pathways: int = 0
    uniprot_statements: int = 0
    has_experimental_function: bool = False
    reasons: list[str] = field(default_factory=list)
    unknown: list[str] = field(default_factory=list)
    resolved_via: str = ""              # locus tag used to retry, when a symbol failed
    locus_candidates: list[str] = field(default_factory=list)

    @property
    def needs_review(self) -> bool:
        """A human has to pick the locus tag; the code refused to."""
        return bool(self.locus_candidates)

    @property
    def lookup_failed(self) -> bool:
        return bool(self.unknown)

    @property
    def thin(self) -> bool:
        # A gene with an unanswered lookup is not thin, and not well covered
        # either -- it is unknown. Calling it thin would route it to literature
        # retrieval because a symbol was wrong or the network hiccuped; calling
        # it covered would silently drop it. An ambiguous locus is the same
        # situation: the gene the verdict describes is not established.
        return bool(self.reasons) and not self.unknown and not self.locus_candidates

    @property
    def functional_gap(self) -> bool:
        """A gap in what is known, rather than a gap in one database's coverage."""
        return (bool(set(self.reasons) & FUNCTIONAL_GAPS)
                and not self.unknown and not self.locus_candidates)

    def to_dict(self) -> dict[str, Any]:
        return dict(asdict(self), thin=self.thin, functional_gap=self.functional_gap,
                    needs_review=self.needs_review)


def classify(gene: str, kegg_pathways: int, uniprot_statements: int,
             has_experimental_function: bool, kegg_known: bool = True,
             uniprot_known: bool = True) -> Coverage:
    """A source that could not answer contributes no reason, only an unknown."""
    reasons, unknown = [], []
    if uniprot_known:
        if uniprot_statements == 0:
            reasons.append(NO_FUNCTION)
        elif not has_experimental_function:
            reasons.append(NO_EXPERIMENTAL)
    else:
        unknown.append("uniprot")
    if kegg_known:
        if kegg_pathways == 0:
            reasons.append(NO_PATHWAY)
    else:
        unknown.append("kegg")
    return Coverage(gene=gene, kegg_pathways=kegg_pathways,
                    uniprot_statements=uniprot_statements,
                    has_experimental_function=has_experimental_function,
                    reasons=reasons, unknown=unknown)


def assess(genes: list[str], kegg: Any, uniprot: Any,
           organism: str = "mtu", organism_id: int = 83332) -> list[Coverage]:
    """One KEGG and one UniProt call per gene.

    A failed lookup must not read as an absent annotation. These clients return a
    ToolResult with a note rather than raising, so an empty record list means
    either "no annotation" or "the request failed" -- and routing on the second
    would send a gene to literature retrieval because the network hiccuped.

    Two discriminators, one per failure mode. A fetch failure leaves no request
    trace, because those paths return before any trace is appended. A failed
    identifier resolution leaves `resolved.matched_by == "none"`, which is how
    both clients report matching nothing. Either way the source is recorded as
    unknown and the gene is routed nowhere, so it surfaces as a list to fix
    rather than disappearing into a bucket.
    """
    out = []
    for gene in genes:
        pathways = kegg.pathways(gene, organism=organism)
        protein = uniprot.protein(gene, organism_id=organism_id)

        # A symbol that resolves in one source and not the other is ordinary
        # nomenclature drift, not an annotation gap. Retry the source that could
        # not answer with the locus tag the other one returned.
        via = ""
        locus, candidates = _locus(protein)
        if not locus and not candidates:
            locus, candidates = _locus(pathways)
        if locus and locus.lower() != gene.lower():
            if not _answered(pathways):
                pathways, via = kegg.pathways(locus, organism=organism), locus
            if not _answered(protein):
                protein, via = uniprot.protein(locus, organism_id=organism_id), locus

        statements = sum(len(r.detail.get("function_statements", []))
                         for r in protein.records)
        experimental = any(r.detail.get("has_experimental_function")
                           for r in protein.records)
        coverage = classify(
            gene, len(pathways.records), statements, experimental,
            kegg_known=_answered(pathways), uniprot_known=_answered(protein))
        coverage.resolved_via = via
        coverage.locus_candidates = candidates
        out.append(coverage)
    return out


def _locus(result: Any) -> tuple[str, list[str]]:
    """The locus tag a source attached to what it returned, or the candidates.

    Returns (tag, []) when one tag is established and ([], candidates) when
    several compete. Taking the first record's tag would let an unreviewed
    fragment or a paralogue redirect the retry to a different gene, so the
    reviewed entry breaks the tie and an unbroken tie is refused.
    """
    tagged: list[tuple[bool, str]] = []
    for record in result.records:
        detail = getattr(record, "detail", {}) or {}
        tags = detail.get("locus_tags") or []
        if tags:
            tagged.append((bool(detail.get("reviewed")), str(tags[0])))
            continue
        kegg_id = detail.get("kegg_gene_id") or ""
        if ":" in kegg_id:
            # KEGG returns one gene per lookup, so its tag is never ambiguous.
            return kegg_id.split(":", 1)[1], []

    distinct = sorted({tag for _, tag in tagged})
    if len(distinct) <= 1:
        return (distinct[0] if distinct else ""), []
    reviewed = sorted({tag for is_reviewed, tag in tagged if is_reviewed})
    if len(reviewed) == 1:
        return reviewed[0], []
    return "", distinct


def _answered(result: Any) -> bool:
    """Did the source actually speak about this gene?

    Empty records alone do not settle it: "resolved, and holds nothing" is a
    finding, while "never resolved" is not.
    """
    if not result.requests:
        return False
    return bool(result.records) or result.resolved.get("matched_by") not in (None, "none")


def route(coverages: list[Coverage], functional_only: bool = False) -> list[str]:
    """The genes stage 2 should run on.

    Never includes a failed lookup or a gene whose locus tag is ambiguous: both
    mean the verdict may not be about the gene that was asked for.
    """
    return [c.gene for c in coverages
            if (c.functional_gap if functional_only else c.thin)]


def summarise(coverages: list[Coverage]) -> dict[str, Any]:
    counts: dict[str, int] = {}
    for coverage in coverages:
        for reason in coverage.reasons:
            counts[reason] = counts.get(reason, 0) + 1
    failed = {c.gene: c.unknown for c in coverages if c.unknown}
    review = {c.gene: c.locus_candidates for c in coverages if c.needs_review}
    return {"genes": len(coverages),
            "needs_review": review,
            "reason_counts": counts,
            "thin": sum(1 for c in coverages if c.thin),
            "functional_gap": sum(1 for c in coverages if c.functional_gap),
            "well_covered": sum(1 for c in coverages if not c.thin
                                and not c.lookup_failed and not c.needs_review),
            "lookup_failed": failed}
