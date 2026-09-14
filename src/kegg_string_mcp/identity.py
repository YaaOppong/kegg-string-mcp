"""One place where a gene's names are settled, and one place that says when they
were not.

Before this module, every consumer re-derived gene identity for itself:
`retrieval/independence.py` keyed STRING partners by preferred name,
`agent/evidence.py` matched a partner by name or by the suffix of its STRING ID,
`retrieval/corpus.py` compiled a regex from the literal query string,
`retrieval/coverage.py` had its own locus-tag retry, and `agent/store.py` built a
private alias set for citation checks. Five reasonable rules, no two the same, so
the answer depended on which spelling the caller happened to type:

    classify(["Rv0678", "Rv0676c"])  -> silent, 0.0
    classify(["mmpR5",  "mmpL5"])    -> corroborating, 0.989   # the same edge

That is the failure this module exists to remove. A gene resolves once, against
every source that can speak to it, and the result carries both the identifiers
(`string_id`, `kegg_gene_id`, `locus_tag`, `accession`) and the honest record of
which sources answered (`matched_by`). Downstream code asks the identity rather
than the string.

**Aliases are filtered, not accumulated.** Matching a gene's synonyms in text
fixes an undercount (a paper saying "mmpR5" is about Rv0678) but can create an
overcount if a synonym is a name other genes also carry. KEGG's own gene index
already knows which symbols are not unique, so those are refused, as are names
claimed by more than one gene of the set being resolved. Every refusal is
recorded in `rejected` rather than dropped, because a silently missing alias is
the same class of bug as a silently missing edge.

**An unresolved gene is not a negative result.** `matched_by[source] == "none"`
is the same convention the tools already use in `ToolResult.resolved`, and it is
what lets a consumer distinguish "STRING says these proteins do not interact"
from "STRING never found one of them" -- a distinction the pair verdicts, the
independence classifier and the residue all previously collapsed into silence.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import asdict, dataclass, field
from typing import Any

from kegg_string_mcp.http import FetchError

MTB_H37RV = 83332
DEFAULT_ORGANISM = "mtu"

# An alias shorter than this is not usable for text matching: two-letter tokens
# collide with ordinary words and with parts of other identifiers, and the cost of
# a false co-mention (a pair wrongly called "already known") is higher than the
# cost of missing a synonym.
MIN_ALIAS_LENGTH = 3

SOURCES = ("string", "kegg", "uniprot")

# How many suggested spellings to put to the structured sources before giving up.
# Each one costs a request per unresolved source, and a proposer that needs more
# than a few tries is guessing rather than reading.
MAX_PROPOSALS = 3


@dataclass
class GeneIdentity:
    """Everything known about which gene a caller meant, and what answered."""

    query: str
    string_id: str = ""          # 83332.Rv0678
    preferred_name: str = ""     # STRING's display symbol, e.g. mmpR5
    kegg_gene_id: str = ""       # mtu:Rv0678
    locus_tag: str = ""          # Rv0678
    accession: str = ""          # UniProt, reviewed entry preferred
    aliases: list[str] = field(default_factory=list)
    rejected: dict[str, str] = field(default_factory=dict)   # alias -> why refused
    matched_by: dict[str, str] = field(default_factory=dict)  # source -> how, or "none"
    # source -> the spelling that worked, when it was not the caller's. A silent
    # cross-source retry would make the run unreproducible by hand: the reader
    # could not tell which identifier the answer is actually about.
    resolved_via: dict[str, str] = field(default_factory=dict)
    tried: dict[str, list[str]] = field(default_factory=dict)  # source -> spellings attempted
    # Spellings suggested by a proposer (e.g. a literature search), with the
    # evidence offered for each. Kept whether or not a source confirmed them.
    proposals: dict[str, str] = field(default_factory=dict)
    notes: list[str] = field(default_factory=list)

    def answered(self, source: str) -> bool:
        """Did this source resolve the identifier at all?

        Deliberately not `bool(records)`: a source can resolve a gene and hold
        nothing about it, which is a finding, while failing to resolve it is not.
        """
        return self.matched_by.get(source, "none") not in ("none", "")

    @property
    def in_string(self) -> bool:
        return bool(self.string_id)

    @property
    def unresolved_sources(self) -> list[str]:
        return [s for s in SOURCES if s in self.matched_by and not self.answered(s)]

    def names(self) -> list[str]:
        """Every spelling safe to match this gene by, query first."""
        return list(self.aliases)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass
class IdentitySet:
    """The resolved identities for one gene list, addressable by any alias."""

    identities: dict[str, GeneIdentity] = field(default_factory=dict)
    ambiguous: dict[str, list[str]] = field(default_factory=dict)  # alias -> claimants
    notes: list[str] = field(default_factory=list)

    def __iter__(self):
        return iter(self.identities.values())

    def __len__(self) -> int:
        return len(self.identities)

    def __contains__(self, name: str) -> bool:
        return self.get(name) is not None

    def get(self, name: str) -> GeneIdentity | None:
        """Look a gene up by the query string, or by any alias it kept."""
        if name in self.identities:
            return self.identities[name]
        return self._index().get(name.strip().lower())

    def string_id(self, name: str) -> str:
        identity = self.get(name)
        return identity.string_id if identity else ""

    def alias_map(self) -> dict[str, list[str]]:
        """{query: every accepted spelling}. What `Corpus.aliases` holds."""
        return {q: list(i.aliases) for q, i in self.identities.items()}

    def unanswered(self) -> dict[str, list[str]]:
        """{query: sources that never resolved it}. What the residue calls undetermined."""
        return {q: i.unresolved_sources for q, i in self.identities.items()
                if i.unresolved_sources}

    def unresolved(self, source: str = "string") -> list[str]:
        """Queries this source could not resolve, in input order."""
        return [i.query for i in self.identities.values() if not i.answered(source)]

    def _index(self) -> dict[str, GeneIdentity]:
        index: dict[str, GeneIdentity] = {}
        for identity in self.identities.values():
            for alias in identity.aliases:
                index.setdefault(alias.lower(), identity)
            if identity.string_id:
                index.setdefault(identity.string_id.lower(), identity)
            if identity.kegg_gene_id:
                index.setdefault(identity.kegg_gene_id.lower(), identity)
        return index

    def to_dict(self) -> dict[str, Any]:
        return {"identities": {k: v.to_dict() for k, v in self.identities.items()},
                "ambiguous": self.ambiguous, "notes": self.notes}


def resolve(genes: list[str], string: Any | None = None, kegg: Any | None = None,
            uniprot: Any | None = None, organism: str = DEFAULT_ORGANISM,
            species: int = MTB_H37RV,
            proposer: Callable[[GeneIdentity], list[tuple[str, str]]] | None = None,
            max_proposals: int = MAX_PROPOSALS) -> IdentitySet:
    """Resolve each gene against whichever sources are supplied.

    Every source is optional and every source is fail-soft: a client that raises,
    or is not passed at all, leaves its `matched_by` entry absent or "none" rather
    than aborting the run. Resolving 41 genes must not be an all-or-nothing act,
    because the caller's alternative is the ad-hoc string matching this module
    replaces.

    **Sources rescue each other.** A symbol that resolves in one source and not
    another is ordinary nomenclature drift, not an absent gene, so a source that
    failed on the caller's spelling is retried with what the others learnt -- the
    locus tag above all, which is the one identifier every source agrees on. This
    is `retrieval/coverage.py`'s cross-source retry, generalised: asking STRING
    for `Rv0678` after UniProt supplied that tag resolves a gene that
    `get_string_ids` refused by symbol, and an unresolved gene is a real absence
    rather than a spelling the caller happened to pick.

    Cost is one `gene_index` call for the organism plus at most two calls per gene
    (STRING `get_string_ids`, UniProt search), plus one retry per source that
    failed. All go through the disk cache, so a re-run is free.
    """
    index, kegg_symbols = _kegg_index(kegg, organism)
    identities: dict[str, GeneIdentity] = {}

    for gene in genes:
        identity = GeneIdentity(query=gene.strip())
        candidates = [identity.query]
        opening = [identity.query]

        # KEGG first: its index is already in memory, so a locus tag learnt here
        # costs nothing and may save the other two a failed request.
        if index is not None:
            candidates += _from_kegg(identity, index, kegg_symbols, organism, opening)
        if uniprot is not None:
            candidates += _from_uniprot(identity, uniprot, species, opening)
        if string is not None:
            # The caller's spelling first, then anything KEGG or UniProt already
            # supplied -- so a symbol STRING refuses is tried as its locus tag in
            # the same pass rather than waiting for the retry below.
            candidates += _from_string(identity, string, species,
                                       [identity.query, *_retry_names(identity, opening)])

        # Second pass: whatever is still unresolved, retried with everything the
        # successful sources turned up.
        retry = _retry_names(identity, opening)
        if retry:
            if index is not None and not identity.answered("kegg"):
                candidates += _from_kegg(identity, index, kegg_symbols, organism, retry)
            if uniprot is not None and not identity.answered("uniprot"):
                candidates += _from_uniprot(identity, uniprot, species, retry)
            if string is not None and not identity.answered("string"):
                candidates += _from_string(identity, string, species, retry)

        # Third pass, only if something is still unresolved: a proposer may
        # suggest spellings the structured sources never offered -- a gene renamed
        # since the release the query was written against. Every suggestion is
        # then put to the structured sources, and only a source resolving it makes
        # it true here.
        if proposer is not None and identity.unresolved_sources:
            candidates += _from_proposals(identity, proposer, max_proposals, index=index,
                                          kegg_symbols=kegg_symbols, organism=organism,
                                          uniprot=uniprot, string=string, species=species)

        identity.aliases = _accept(candidates, identity,
                                   ambiguous=index.ambiguous if index else set())
        identities[identity.query] = identity

    return _resolve_collisions(IdentitySet(identities=identities))


def _retry_names(identity: GeneIdentity, already: list[str]) -> list[str]:
    """Spellings learnt from other sources, best first.

    The locus tag leads because it is the identifier KEGG, STRING and UniProt all
    carry; a display symbol is the one most likely to differ between them.
    """
    seen = {n.lower() for n in already}
    out = []
    for name in (identity.locus_tag, identity.kegg_gene_id.split(":")[-1],
                 identity.preferred_name):
        if name and name.lower() not in seen:
            seen.add(name.lower())
            out.append(name)
    return out


def _kegg_index(kegg: Any | None, organism: str) -> tuple[Any | None, dict[str, list[str]]]:
    """KEGG's organism-wide symbol table, plus its inverse.

    `gene_index` maps SYMBOL -> KEGG ID; the inverse gives every symbol KEGG
    records for a gene, which is the alias source STRING's single preferred name
    cannot provide.
    """
    if kegg is None:
        return None, {}
    try:
        index, _ = kegg.gene_index(organism)
    except (FetchError, ValueError):
        return None, {}
    symbols: dict[str, list[str]] = {}
    for symbol, kegg_id in index.entries.items():
        if symbol not in index.locus_tags:
            symbols.setdefault(kegg_id, []).append(symbol)
    return index, symbols


def _from_string(identity: GeneIdentity, string: Any, species: int,
                 names: list[str]) -> list[str]:
    for name in names:
        identity.tried.setdefault("string", []).append(name)
        try:
            hit, _ = string.resolve(name, species)
        except (FetchError, ValueError) as exc:
            identity.matched_by["string"] = "none"
            identity.notes.append(f"STRING lookup failed for '{name}' ({exc}); "
                                  f"this is a retrieval failure, not evidence about the gene.")
            continue
        if not hit:
            identity.matched_by["string"] = "none"
            continue
        identity.string_id = hit.get("stringId", "")
        identity.preferred_name = hit.get("preferredName", "")
        identity.matched_by["string"] = "get_string_ids"
        _record_retry(identity, "string", name)
        # STRING protein IDs are "{taxon}.{locus}", so the suffix is the locus tag
        # whenever the caller asked by symbol.
        suffix = identity.string_id.split(".", 1)[-1] if "." in identity.string_id else ""
        if suffix and not identity.locus_tag:
            identity.locus_tag = suffix
        return [identity.preferred_name, suffix]
    identity.notes.append(f"STRING resolved none of {', '.join(names)} in species {species}.")
    return []


def _from_proposals(identity: GeneIdentity, proposer: Callable[[GeneIdentity], list[tuple[str, str]]],
                    limit: int, *, index: Any, kegg_symbols: dict[str, list[str]], organism: str,
                    uniprot: Any, string: Any, species: int) -> list[str]:
    """Last resort: a suggested spelling, accepted only if a source resolves it.

    The suggestion may come from anywhere -- a literature search read by a model
    is the case this exists for, since a gene renamed after a release is exactly
    what the structured sources cannot tell you about the old name. What the
    proposer may NOT do is settle the question: it returns candidate spellings
    with the evidence for each, and identity here re-asks KEGG, UniProt and
    STRING. A name no source confirms is refused and recorded, and the gene stays
    unresolved. That keeps the rule this codebase is built on -- the model reads,
    the structured sources assert -- and means a hallucinated identifier fails
    closed rather than propagating as a resolved gene.
    """
    try:
        proposals = list(proposer(identity))[:limit]
    except Exception as exc:  # noqa: BLE001 -- external code (a model client, a
        # network call); a last-resort rescue failing must not end a 41-gene run.
        identity.notes.append(f"identifier proposer failed for '{identity.query}' ({exc}).")
        return []

    found: list[str] = []
    for name, evidence in proposals:
        name = (name or "").strip()
        if not name or name.lower() == identity.query.lower():
            continue
        identity.proposals[name] = evidence
        before = set(identity.unresolved_sources)
        if index is not None and "kegg" in before:
            found += _from_kegg(identity, index, kegg_symbols, organism, [name])
        if uniprot is not None and "uniprot" in before:
            found += _from_uniprot(identity, uniprot, species, [name])
        if string is not None and "string" in before:
            found += _from_string(identity, string, species, [name])
        confirmed = before - set(identity.unresolved_sources)
        if confirmed:
            identity.notes.append(
                f"'{identity.query}' was unresolved; the proposed identifier '{name}' "
                f"({evidence}) was confirmed by {', '.join(sorted(confirmed))}.")
            return found
        identity.rejected[name] = ("proposed for this gene but no structured source "
                                   "resolved it, so it was not accepted")
    if proposals:
        identity.notes.append(
            f"No proposed identifier for '{identity.query}' was confirmed by any source; "
            f"it remains unresolved rather than being resolved on the proposal alone.")
    return found


def _record_retry(identity: GeneIdentity, source: str, name: str) -> None:
    if name.lower() != identity.query.lower():
        identity.resolved_via[source] = name
        identity.notes.append(f"{source} did not resolve '{identity.query}' but did resolve "
                              f"'{name}', which another source gave for the same gene.")


def _from_kegg(identity: GeneIdentity, index: Any, symbols: dict[str, list[str]],
               organism: str, names: list[str]) -> list[str]:
    for name in names:
        identity.tried.setdefault("kegg", []).append(name)
        kegg_id = index.get((name or "").strip().upper())
        if kegg_id:
            identity.kegg_gene_id = kegg_id
            identity.matched_by["kegg"] = ("locus_tag" if name.strip().upper() in index.locus_tags
                                           else "symbol")
            _record_retry(identity, "kegg", name)
            tag = kegg_id.split(":", 1)[1] if ":" in kegg_id else kegg_id
            identity.locus_tag = identity.locus_tag or tag
            return [tag, *symbols.get(kegg_id, [])]
    identity.matched_by["kegg"] = "none"
    identity.notes.append(f"None of {', '.join(names)} is in the KEGG gene index for {organism}.")
    return []


def _from_uniprot(identity: GeneIdentity, uniprot: Any, species: int,
                  names: list[str]) -> list[str]:
    result = None
    for name in names:
        identity.tried.setdefault("uniprot", []).append(name)
        try:
            candidate = uniprot.protein(name, organism_id=species)
        except (FetchError, ValueError) as exc:
            identity.matched_by["uniprot"] = "none"
            identity.notes.append(f"UniProt lookup failed for '{name}' ({exc}).")
            continue
        if candidate.records:
            result = candidate
            _record_retry(identity, "uniprot", name)
            break
    if result is None:
        identity.matched_by["uniprot"] = "none"
        return []
    # Reviewed first. A symbol search returns TrEMBL fragments alongside the
    # Swiss-Prot entry, and their gene-name lists are the ones most likely to be
    # wrong -- the same tie-break `retrieval/coverage.py` already makes.
    records = sorted(result.records, key=lambda r: not r.detail.get("reviewed"))
    identity.accession = records[0].record_id
    identity.matched_by["uniprot"] = "uniprot_search"
    names: list[str] = []
    for record in records:
        names += [n for n in record.detail.get("gene_names", []) if n]
        tags = [t for t in record.detail.get("locus_tags", []) if t]
        names += tags
        if tags and not identity.locus_tag:
            identity.locus_tag = tags[0]
    return names


def _accept(candidates: list[str], identity: GeneIdentity,
            ambiguous: set[str]) -> list[str]:
    """Keep the spellings that can only mean this gene. Query string always first."""
    kept: list[str] = []
    seen: set[str] = {identity.query.lower()}
    for raw in candidates:
        name = (raw or "").strip()
        if not name or name.lower() in seen:
            continue
        if len(name) < MIN_ALIAS_LENGTH:
            identity.rejected[name] = f"shorter than {MIN_ALIAS_LENGTH} characters"
            continue
        if name.isdigit():
            identity.rejected[name] = "numeric, not a usable name"
            continue
        if name.upper() in ambiguous:
            identity.rejected[name] = "KEGG records this symbol for more than one gene"
            continue
        seen.add(name.lower())
        kept.append(name)
    return [identity.query, *kept]


def _resolve_collisions(identities: IdentitySet) -> IdentitySet:
    """Refuse any alias claimed by two genes of the same set.

    Within one run the set is small and known, so a collision is decidable rather
    than hypothetical: if `eis` were an alias of two resolved genes, matching it in
    text could not say which was meant, and a wrong co-mention makes a pair look
    explained -- the exact direction of error this work is removing.
    """
    claims: dict[str, list[str]] = {}
    for identity in identities.identities.values():
        for alias in identity.aliases:
            claims.setdefault(alias.lower(), []).append(identity.query)

    for alias, claimants in sorted(claims.items()):
        if len(set(claimants)) < 2:
            continue
        identities.ambiguous[alias] = sorted(set(claimants))
        for query in claimants:
            identity = identities.identities[query]
            # The query string itself is never refused: the caller asked for it,
            # and dropping it would leave a gene with no name at all.
            identity.aliases = [a for a in identity.aliases
                                if a.lower() != alias or a.lower() == identity.query.lower()]
            others = [q for q in sorted(set(claimants)) if q != query]
            if alias != query.lower():
                identity.rejected[alias] = f"also names {', '.join(others)} in this gene set"
    return identities
