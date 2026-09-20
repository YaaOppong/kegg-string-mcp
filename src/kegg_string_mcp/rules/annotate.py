"""Annotate every locus a rule set uses, once, and link the ones that co-occur.

The expensive half of the workflow, and the half whose output is the artefact:
one row per distinct locus, which is the supplementary table a reader checks an
interpretation against.

Three shapes of work, in increasing cost:

* **Free.** Coordinates, strand and length from the supplied annotation; the
  catalogue and lineage lookups, which are two files fetched once for the run;
  and the population counts, which are a pass over the rules themselves.
* **Per locus.** KEGG pathways, UniProt product, STRING partners. A classifier
  population reuses its loci heavily, so this is done once per distinct locus
  rather than once per condition -- a few hundred rules over a hundred loci
  costs a hundred lookups.
* **Per co-occurring pair.** Whether two loci are functionally linked. Restricted
  to pairs that actually appear together in some rule: asking about every pair of
  a hundred loci is thousands of lookups for links no rule would ever use.

Every lookup is fail-soft. A source that errors leaves its field empty and says
so in `notes`, because annotating 100 loci must not be an all-or-nothing act.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from itertools import combinations
from typing import Any

from kegg_string_mcp.http import FetchError
from kegg_string_mcp.rules.annotation import Gene
from kegg_string_mcp.rules.catalogue import CatalogueStatus
from kegg_string_mcp.rules.evidence import Sources
from kegg_string_mcp.rules.parse import Rule
from kegg_string_mcp.rules.signature import Link

# Two loci closer than this are reported as adjacent. Operon structure is not
# being called here -- that needs strand and expression evidence -- but genes a
# few bases apart are co-transcribed often enough that the distance is worth
# carrying. katG and furA are 6 bp apart; fabG1 and inhA, 19.
ADJACENT_BP = 200

# STRING's own band for "this edge is worth reporting at all". Below it the score
# is noise, and a link column full of 0.2 edges is worse than an empty one.
LINK_MIN_SCORE = 400


@dataclass
class LinkedLocus:
    """One entry in a locus's `linked_loci` column."""

    locus: str
    source: str            # string | kegg | adjacent
    detail: str

    def render(self) -> str:
        return f"{self.locus}:{self.source}:{self.detail}"


@dataclass
class LocusAnnotation:
    """One row of the supplementary table."""

    locus: str
    symbol: str = ""
    gene: Gene | None = None
    product: str = ""
    string_id: str = ""
    pathways: list[tuple[str, str]] = field(default_factory=list)
    partners: list[tuple[str, str, float]] = field(default_factory=list)
    lineages: tuple[str, ...] = ()
    catalogue: CatalogueStatus | None = None
    linked: list[LinkedLocus] = field(default_factory=list)
    labels: tuple[str, ...] = ()          # the rule-file spellings that reached here
    n_rules: int = 0
    n_present: int = 0                    # rules asserting state 1
    n_absent: int = 0                     # rules asserting state 0
    with_anchor: int = 0                  # present-rules that also contain an anchor
    notes: list[str] = field(default_factory=list)

    @property
    def without_anchor(self) -> int:
        """Present-rules with no known resistance locus in them.

        The number that matters for compensation. A compensator turns up beside a
        resistance locus and rarely alone, so a locus with a high `with_anchor`
        and a `without_anchor` near zero is behaving like one -- computed from
        your own rule set, not from any external source.
        """
        return self.n_present - self.with_anchor


def _fail_soft(note_target: list[str], what: str, call):
    try:
        return call()
    except (FetchError, ValueError, KeyError) as exc:
        note_target.append(f"{what} unavailable: {type(exc).__name__}: {exc}")
        return None


def annotate_locus(locus: str, symbol: str, sources: Sources, clients: Any) -> LocusAnnotation:
    """Everything known about one locus, from whichever sources answer."""
    gene = sources.annotation.gene(locus)
    out = LocusAnnotation(locus=locus, symbol=symbol or (gene.symbol if gene else ""), gene=gene)

    if gene is not None:
        out.lineages = sources.lineages_in(gene.start, gene.end)
    out.catalogue = sources.catalogue.for_locus(locus, out.symbol)

    query = out.symbol or locus
    kegg = getattr(clients, "kegg", None)
    if kegg is not None:
        result = _fail_soft(out.notes, "KEGG pathways",
                            lambda: kegg.pathways(query, sources.organism).model_dump())
        if result:
            out.pathways = [(r["record_id"], r.get("name", "")) for r in result.get("records", [])]

    uniprot = getattr(clients, "uniprot", None)
    if uniprot is not None:
        result = _fail_soft(out.notes, "UniProt", lambda: uniprot.protein(query).model_dump())
        records = (result or {}).get("records", [])
        if records:
            out.product = records[0].get("name", "")

    string = getattr(clients, "string", None)
    if string is not None:
        result = _fail_soft(out.notes, "STRING partners",
                            lambda: string.partners(query, limit=20).model_dump())
        # The locus's own STRING id comes free with its partner list, and one
        # `network()` call over all of them then gives every edge among the loci
        # without a lookup per pair.
        out.string_id = str((result or {}).get("resolved", {}).get("string_id", ""))
        for record in (result or {}).get("records", []):
            detail = record.get("detail", {})
            out.partners.append((record["record_id"], record.get("name", ""),
                                 float(detail.get("combined_score", 0.0))))
    return out


def co_occurring(rules: list[Rule], rename: dict[str, str]) -> set[tuple[str, str]]:
    """Unordered locus pairs that appear together in at least one rule.

    The bound on how much linking work there is. Pairs are keyed on the resolved
    locus so two spellings of the same gene do not produce two pairs.
    """
    pairs: set[tuple[str, str]] = set()
    for rule in rules:
        loci = sorted({rename.get(c.label, c.label) for c in rule.conditions})
        pairs.update(combinations(loci, 2))
    return pairs


def compute_links(pairs: set[tuple[str, str]], annotations: dict[str, LocusAnnotation],
                  sources: Sources, pathway_sizes: dict[str, int] | None = None,
                  genome_size: int = 0,
                  edges: dict[tuple[str, str], dict[str, Any]] | None = None,
                  ) -> dict[tuple[str, str], list[Link]]:
    """Why two co-occurring loci might be related, from the structured sources.

    Three kinds, and they are not equivalent. A STRING edge whose support is only
    textmining is co-mention in papers, which is the same evidence a literature
    search would find rather than a second line of it -- so the channel is
    carried, never flattened into one score.
    """
    edges = edges or {}
    sizes = pathway_sizes or {}
    out: dict[tuple[str, str], list[Link]] = {}

    for left, right in sorted(pairs):
        found: list[Link] = []
        a, b = annotations.get(left), annotations.get(right)
        if a is None or b is None:
            continue

        edge = edges.get((left, right)) or edges.get((right, left))
        if edge and float(edge.get("combined_score", 0)) * 1000 >= LINK_MIN_SCORE:
            beyond = bool(edge.get("evidence_beyond_textmining"))
            found.append(Link(
                kind="string_beyond_textmining" if beyond else "string_textmining_only",
                detail=f"{float(edge.get('combined_score', 0)):.3f}"))

        shared = {p for p, _ in a.pathways} & {p for p, _ in b.pathways}
        for pathway in sorted(shared):
            size = sizes.get(pathway, 0)
            # A pathway holding a sixth of the genome is a base rate, not a link.
            # mtu01100 has 698 of ~4,000 genes; mtu00983 has 11.
            if genome_size and size and size / genome_size >= 0.05:
                continue
            found.append(Link(kind="kegg_shared_pathway",
                              detail=f"{pathway}" + (f":{size}genes" if size else "")))

        if a.gene is not None and b.gene is not None:
            gap = max(a.gene.start, b.gene.start) - min(a.gene.end, b.gene.end) - 1
            if 0 <= gap <= ADJACENT_BP:
                found.append(Link(kind="adjacent", detail=f"{gap}bp"))

        if found:
            out[(left, right)] = found
    return out


def population_counts(rules: list[Rule], rename: dict[str, str],
                      anchors: set[str], annotations: dict[str, LocusAnnotation]) -> None:
    """Fill the appears-with / appears-without counts, in place.

    One pass over the rules, no lookups. `anchors` is the set of loci the
    catalogue grades resistance-associated, so "appeared beside an anchor" is a
    question about this rule set rather than about the literature.
    """
    for rule in rules:
        present = {rename.get(c.label, c.label) for c in rule.conditions if c.state == 1}
        absent = {rename.get(c.label, c.label) for c in rule.conditions if c.state == 0}
        has_anchor = bool(present & anchors)
        for locus in present | absent:
            annotation = annotations.get(locus)
            if annotation is None:
                continue
            annotation.n_rules += 1
            if locus in present:
                annotation.n_present += 1
                if has_anchor and locus not in anchors:
                    annotation.with_anchor += 1
            else:
                annotation.n_absent += 1
