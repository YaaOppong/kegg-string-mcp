"""Deterministic evidence assembly. Code computes; the model interprets.

The organising principle of the whole project applied to epistasis: set
intersections, pathway sizes and network degrees are arithmetic, and a language
model doing arithmetic over a hundred IDs will get some of it wrong in a way
nobody can audit. So the pipeline computes the relationships and hands the model
a finished, checkable evidence table. The model's job is to say what it *means*.

Two traps this module exists to defuse:

* **Pathway size.** In M. tuberculosis, `mtu01100` ("Metabolic pathways") holds
  698 of ~4000 genes. Two genes sharing it is a base rate, not a finding.
  `mtu00983` holds 11, and sharing that is real. Every shared pathway is
  therefore reported with its size and an explicit specificity judgement.
* **Hub proteins.** A STRING partner list is not evidence of a *specific*
  relationship if the protein in question partners with everything. Shared
  partners are reported with the degree that produced them.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from itertools import combinations
from typing import Any

# A pathway holding more than this share of the annotated genome is a container,
# not a mechanism. mtu01100 (698 genes) and mtu01110 (373) are the obvious cases.
BROAD_PATHWAY_FRACTION = 0.05
SPECIFIC_PATHWAY_MAX = 50


@dataclass
class SharedPathway:
    pathway_id: str
    name: str
    size: int
    specificity: str  # "specific" | "moderate" | "broad"
    note: str


@dataclass
class DirectInteraction:
    partner_id: str
    partner_name: str
    combined_score: float
    textmining_score: float
    max_non_textmining_score: float
    evidence_beyond_textmining: bool


@dataclass
class PairEvidence:
    gene_a: str
    gene_b: str
    direct_interaction: DirectInteraction | None = None
    # Was the pair itself put to STRING (via /network), or only inferred from two
    # ranked partner lists? Only the first can support "no direct interaction".
    checked_directly: bool = False
    unresolved: list[str] = field(default_factory=list)   # genes STRING never resolved
    shared_pathways: list[SharedPathway] = field(default_factory=list)
    shared_partners: list[dict[str, Any]] = field(default_factory=list)
    # NOT network degree: the number of partners RETRIEVED, which `limit` caps.
    # Reporting it as degree made a hub with 500 partners indistinguishable from a
    # gene with exactly 20 -- disabling the very hub check it exists for.
    partners_retrieved: dict[str, int] = field(default_factory=dict)
    truncated: list[str] = field(default_factory=list)   # genes whose list hit the limit
    verdict: str = ""

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def classify_pathway(size: int, genome_size: int) -> tuple[str, str]:
    """Judge whether sharing this pathway means anything. Returned as data so the
    reasoning is in the record, not only in a prompt the model may ignore."""
    if not genome_size:
        # Without a denominator the broad/specific judgement cannot be made, and
        # falling through silently labelled mtu01100 (698 genes) "moderate" --
        # reopening the exact base-rate trap this function exists to close.
        return "unknown", (f"{size} genes, but the genome size could not be determined, so whether "
                           f"this is a container category is unknown. Treat co-membership with caution.")
    if size / genome_size >= BROAD_PATHWAY_FRACTION:
        return "broad", (
            f"{size} genes ({size / genome_size:.0%} of annotated genes) -- a container "
            f"category. Co-membership here is a base rate, not evidence of a link."
        )
    if size <= SPECIFIC_PATHWAY_MAX:
        return "specific", f"{size} genes -- narrow enough that co-membership is informative."
    return "moderate", f"{size} genes -- co-membership is weak evidence on its own."


def pair_evidence(
    gene_a: str,
    gene_b: str,
    pathways: dict[str, list[dict[str, Any]]],
    partners: dict[str, list[dict[str, Any]]],
    pathway_sizes: dict[str, int],
    genome_size: int,
    partner_limit: int | None = None,
    edges: dict[tuple[str, str], dict[str, Any]] | None = None,
    unresolved: set[str] | frozenset[str] | None = None,
) -> PairEvidence:
    """Assemble everything known about ONE pair. Pure function of tool output.

    `edges` holds the pair-level STRING answer (from `StringClient.network`),
    keyed by the unordered pair of gene names as given. It is authoritative:
    partner lists are ranked and truncated, so the absence of B from A's top 20 is
    not the absence of an edge. katG/rpoC scores 0.823 at rank 24 of 31 and was
    reported as "No direct interaction" until the pair itself was queried.

    `unresolved` names the genes STRING could not resolve. Their pairs get no
    verdict about interaction at all -- the tool already says a resolution failure
    is not evidence of no partners, and this is that rule applied downstream.
    """
    ev = PairEvidence(gene_a=gene_a, gene_b=gene_b)
    unresolved = set(unresolved or ())
    ev.unresolved = [g for g in (gene_a, gene_b) if g in unresolved]

    a_partners = partners.get(gene_a, [])
    b_partners = partners.get(gene_b, [])
    ev.partners_retrieved = {gene_a: len(a_partners), gene_b: len(b_partners)}
    if partner_limit:
        ev.truncated = [g for g, p in ((gene_a, a_partners), (gene_b, b_partners))
                        if len(p) >= partner_limit]

    edge = (edges or {}).get(_pair_key(gene_a, gene_b))
    if edges is not None and not ev.unresolved:
        # The pair was asked about directly, so "no edge" is now a real answer.
        ev.checked_directly = True
    if edge:
        ev.direct_interaction = DirectInteraction(
            partner_id=edge.get("record_id", ""),
            partner_name=edge.get("name", ""),
            combined_score=edge.get("combined_score", 0.0),
            textmining_score=edge.get("textmining_score", 0.0),
            max_non_textmining_score=edge.get("max_non_textmining_score", 0.0),
            evidence_beyond_textmining=edge.get("evidence_beyond_textmining", False))

    # Direct interaction: is B in A's partner list (or vice versa)?
    for record in [] if ev.direct_interaction else a_partners + b_partners:
        name = (record.get("name") or "").upper()
        rid = (record.get("record_id") or "").upper()
        other = gene_b.upper() if record in a_partners else gene_a.upper()
        if name == other or rid.endswith("." + other):
            detail = record.get("detail", {})
            ev.direct_interaction = DirectInteraction(
                partner_id=record["record_id"],
                partner_name=record.get("name", ""),
                combined_score=detail.get("combined_score", 0.0),
                textmining_score=detail.get("textmining_score", 0.0),
                max_non_textmining_score=detail.get("max_non_textmining_score", 0.0),
                evidence_beyond_textmining=detail.get("evidence_beyond_textmining", False),
            )
            break

    # Shared pathways, each judged by size.
    a_paths = {r["record_id"]: r.get("name", "") for r in pathways.get(gene_a, [])}
    b_paths = {r["record_id"] for r in pathways.get(gene_b, [])}
    for pid in sorted(a_paths.keys() & b_paths):
        size = pathway_sizes.get(pid, 0)
        specificity, note = classify_pathway(size, genome_size)
        ev.shared_pathways.append(SharedPathway(pid, a_paths[pid], size, specificity, note))

    # Shared network neighbours, with the degrees that produced them.
    a_ids = {r["record_id"]: r.get("name", "") for r in a_partners}
    b_ids = {r["record_id"] for r in b_partners}
    ev.shared_partners = [{"record_id": rid, "name": a_ids[rid]} for rid in sorted(a_ids.keys() & b_ids)]

    ev.verdict = _verdict(ev)
    return ev


def _verdict(ev: PairEvidence) -> str:
    """A deterministic, defensible one-liner. The model may elaborate on it but
    must not contradict it -- and 'no known link' is a real answer, not a gap to
    be filled with plausible prose."""
    specific = [p for p in ev.shared_pathways if p.specificity == "specific"]
    broad_only = ev.shared_pathways and not specific

    if ev.unresolved and not ev.direct_interaction:
        # Say what is missing rather than what is absent. The epistasis prompt
        # tells the model not to contradict this line, so a verdict of "no known
        # link" built on a failed lookup propagates into the summary as fact.
        return (f"STRING could not resolve {', '.join(ev.unresolved)}, so no interaction "
                f"evidence was retrieved for this pair. This is a resolution failure, not "
                f"evidence of no link; the KEGG evidence below stands on its own.")
    if ev.direct_interaction:
        di = ev.direct_interaction
        support = ("supported beyond literature co-mention"
                   if di.evidence_beyond_textmining else
                   "supported essentially only by literature co-mention")
        return f"Direct STRING interaction (combined {di.combined_score}), {support}."
    direct = _no_direct_phrase(ev)
    if specific:
        names = ", ".join(f"{p.pathway_id} ({p.size} genes)" for p in specific)
        return f"{direct}. Share specific pathway(s): {names}."
    if ev.shared_partners:
        caveat = ""
        if ev.truncated:
            # Without a true degree there is no base rate to compare against, and
            # shared partners between two hubs is exactly that: a base rate.
            caveat = (f" Partner lists for {', '.join(ev.truncated)} hit the retrieval limit, so "
                      f"true network degree is unknown and this overlap cannot be distinguished "
                      f"from what any two well-connected proteins would share.")
        return (f"{direct} and no specific shared pathway, but "
                f"{len(ev.shared_partners)} shared network partner(s).{caveat}")
    if broad_only:
        return (f"{direct}. Shared pathways are broad container categories "
                f"only, which is not evidence of a mechanistic link.")
    if len(ev.truncated) == 2:
        # The only case where "no known link" would be a claim about the genes
        # rather than about the query: both lists were full, so an edge between
        # them could sit below both cuts and was never seen.
        return (f"{direct}, and no shared pathway. Absence of an edge here is a limit of the "
                f"retrieval, not a negative result.")
    return "No known link in KEGG or STRING at the thresholds queried."


def _no_direct_phrase(ev: PairEvidence) -> str:
    """How far "no direct interaction" can honestly be pushed.

    Three different facts wore the same words. Only a pair queried directly
    supports "no interaction"; two truncated partner lists support nothing at all,
    because the edge can sit below both cuts -- which is how rpoB/rpsC, STRING
    0.991, was reported as having none.
    """
    if ev.checked_directly:
        return "No direct STRING interaction at the threshold queried"
    if len(ev.truncated) == 2:
        limit = min(ev.partners_retrieved.values(), default=0)
        return (f"Direct interaction NOT CHECKED beyond the top {limit} partners of each gene, "
                f"and both lists were full")
    return "No direct interaction in the partner lists retrieved"


def all_pairs(
    genes: list[str],
    pathways: dict[str, list[dict[str, Any]]],
    partners: dict[str, list[dict[str, Any]]],
    pathway_sizes: dict[str, int],
    genome_size: int,
    partner_limit: int | None = None,
    edges: dict[tuple[str, str], dict[str, Any]] | None = None,
    unresolved: set[str] | frozenset[str] | None = None,
) -> list[PairEvidence]:
    return [
        pair_evidence(a, b, pathways, partners, pathway_sizes, genome_size, partner_limit,
                      edges=edges, unresolved=unresolved)
        for a, b in combinations(genes, 2)
    ]


def _pair_key(a: str, b: str) -> tuple[str, str]:
    """Unordered, case-insensitive. Callers build pair keys from different sources."""
    return tuple(sorted((a.strip().lower(), b.strip().lower())))  # type: ignore[return-value]


def edge_index(network_result: dict[str, Any], identities: Any,
               genes: list[str]) -> dict[tuple[str, str], dict[str, Any]]:
    """Map a `string_partners`-style network result onto the caller's gene names.

    The model and the upstream analysis speak in whatever spelling they were given;
    STRING answers in protein IDs. Translating once, here, is what keeps
    `pair_evidence` a pure function of names.
    """
    by_id = {identities.string_id(g): g for g in genes if identities.string_id(g)}
    out: dict[tuple[str, str], dict[str, Any]] = {}
    for record in network_result.get("records", []):
        detail = record.get("detail", {})
        gene_a = by_id.get(detail.get("string_id_a", ""))
        gene_b = by_id.get(detail.get("string_id_b", ""))
        if not gene_a or not gene_b:
            continue
        out[_pair_key(gene_a, gene_b)] = dict(detail) | {
            "record_id": record.get("record_id", ""), "name": record.get("name", "")}
    return out
