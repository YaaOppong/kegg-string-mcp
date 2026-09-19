"""Tabular views of runs already on disk: one row per gene, one row per pair.

Built from the run stores, never from the summaries. Every column here has a
deterministic source in a `tool_result` or a `derived` line that was written
before the model saw it, so the tables carry what the tools returned rather than
what the prose said they returned -- the same rule `validate.py` enforces on
citations, applied to the numbers.

That distinction is invisible at one gene and load-bearing at a hundred. Four
situations read identically in prose and must not read identically in a cell:

    kegg_n = 0      KEGG was asked and assigns no pathway (true of 71% of genes)
    kegg_n = NA     KEGG never resolved the gene -- says nothing about the gene
    string_n = 20, string_truncated = True     the list hit `limit`; absence
                                               below the cut is a query artefact
    resistance_associated = NA                 not in the catalogue: never
                                               assessed, which is not "negative"

`NA` is written for all of them rather than an empty cell, so R and pandas both
read it as missing rather than as zero or as an empty string.

The unit of the gene table is the *resolved* gene, not the queried spelling. A
run may call tools for a neighbour -- annotating furA and looking up katG for
context is legitimate, and the validator already reports it as cross-target -- so
counts are taken only from calls whose gene argument matches the target's alias
set. Without that, a neighbour's pathways would land in the target's row, which
is the tabular form of the error the cross-target check exists to catch.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

NA = "NA"
LIST_SEPARATOR = "|"

GENE_COLUMNS = [
    "query", "resolved_id", "locus_tag", "preferred_name", "accession",
    "aliases", "n_aliases", "aliases_rejected",
    "matched_by_string", "matched_by_kegg", "matched_by_uniprot",
    "kegg_n", "string_n", "string_truncated", "uniprot_n",
    "pubmed_n", "corpus_n", "lineage_n",
    "resistance_associated", "resistance_drugs",
    "variants_associated_n", "variants_in_catalogue",
    "citations_verified", "citations_total", "quotes_verified", "quotes_total",
    "validation_passed", "validation_flags",
    "turns", "stop_reason", "run_id", "store",
]

# The specific variants and markers a locus contains. They are in every run --
# one record each, fully detailed -- but a gene row can only count them, and a
# count is the one thing you cannot act on: "katG contains 138 associated
# variants" does not tell you whether YOUR variant is one of them. These two
# tables are the long form, one row per variant and per marker, joinable to
# genes.tsv on `resolved_id`.
VARIANT_COLUMNS = [
    "run_id", "resolved_id", "gene_query", "variant_id", "gene", "mutation",
    "drug", "confidence", "associated", "source", "comment", "store",
]

MARKER_COLUMNS = [
    "run_id", "resolved_id", "gene_query", "marker_id", "position", "lineage",
    "lineage_name", "allele", "store",
]

PAIR_COLUMNS = [
    "run_id", "organism", "gene_a", "gene_b", "verdict",
    "checked_directly", "direct_partner_id", "direct_score",
    "textmining_score", "max_non_textmining_score", "evidence_beyond_textmining",
    "n_shared_pathways", "shared_specific", "shared_broad",
    "n_shared_partners", "shared_partners",
    "partners_retrieved_a", "partners_retrieved_b", "truncated", "unresolved",
    # Why the pair might co-occur without being related. `confounds` is the token
    # form for filtering; the sentence is in `verdict`.
    "lineage_a_n", "lineage_b_n", "shared_lineages",
    "resistance_a", "resistance_b", "shared_drugs", "confounds",
    "store",
]

# Which resolved-identity field is the join key, in order of preference. The KEGG
# gene ID is organism-qualified (`mtu:Rv1908c`) and so is unique without further
# context; the locus tag is the identifier every source agrees on when KEGG is
# silent, which for 71% of genes it is.
KEY_FIELDS = ("kegg_gene_id", "locus_tag", "string_id", "accession")


@dataclass
class Run:
    """One store, parsed. `derived` keeps the last value per label, which is what
    the pipeline wrote last and therefore what the run finished with."""

    path: Path
    run_id: str = ""
    mode: str = ""
    output: dict[str, Any] = field(default_factory=dict)
    derived: dict[str, Any] = field(default_factory=dict)
    calls: list[dict[str, Any]] = field(default_factory=list)


def load_run(path: Path) -> Run | None:
    """Parse one store. A truncated final line is expected rather than fatal: the
    store is append-only and a run killed mid-write leaves one, and refusing to
    table a whole batch because one gene was interrupted is the wrong trade."""
    run = Run(path=path)
    with path.open(encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            try:
                entry = json.loads(line)
            except json.JSONDecodeError:
                continue
            run.run_id = run.run_id or entry.get("run_id", "")
            kind = entry.get("kind")
            if kind == "tool_result":
                run.calls.append(entry)
            elif kind == "derived":
                run.derived[entry.get("label", "")] = entry.get("payload", {})
            elif kind == "output":
                run.output = entry
    if not run.output:
        # No `output` line means the run never reached `_finish`: it crashed, or is
        # still going. Reporting its partial tool calls as a finished row would put
        # an unvalidated annotation in the table indistinguishable from a good one.
        return None
    run.mode = run.output.get("mode", "")
    return run


def _identity(run: Run) -> dict[str, Any]:
    """The single-mode gene's resolved identity, or an empty one."""
    identities = (run.derived.get("identities") or {}).get("identities") or {}
    if not identities:
        return {}
    target = str(run.output.get("gene", "")).strip()
    return identities.get(target) or next(iter(identities.values()))


def _alias_keys(identity: dict[str, Any], query: str) -> set[str]:
    """Every spelling that means "this gene", lowercased, for matching a call's
    `gene` argument. The query itself is included because a run whose identity
    resolution found nothing still annotated the gene the caller asked for."""
    keys = {query.strip().lower()}
    keys.update(str(a).strip().lower() for a in identity.get("aliases", []))
    for field_name in KEY_FIELDS:
        value = str(identity.get(field_name) or "").strip().lower()
        if not value:
            continue
        keys.add(value)
        # `mtu:Rv1908c` and `83332.Rv1908c` both carry the bare tag the model is
        # just as likely to have typed.
        for separator in (":", "."):
            if separator in value:
                keys.add(value.split(separator, 1)[1])
    keys.discard("")
    return keys


def _calls_for(run: Run, tool: str, keys: set[str]) -> list[dict[str, Any]]:
    """Calls of `tool` made about this gene. `corpus_search` takes free text
    rather than a gene, so its calls belong to the run's target by construction."""
    out = []
    for entry in run.calls:
        if entry.get("tool") != tool:
            continue
        if tool == "corpus_search":
            out.append(entry)
            continue
        gene = str(entry.get("arguments", {}).get("gene", "")).strip().lower()
        if gene in keys:
            out.append(entry)
    return out


# Whether a tool's empty result is a zero or a missing value is NOT uniform, so it
# is declared per tool rather than inferred from one shared convention:
#
#   identity   `resolved.matched_by == "none"` means the gene was never resolved
#              here, so the count is missing rather than zero. KEGG, STRING,
#              UniProt and the WHO catalogue all use it that way.
#   performed  `lineage_markers` returns "none" for a gene with no marker too --
#              its index holds only genes that have one, so it cannot tell the
#              two apart by identity. There "none" is the informative negative
#              the tool's own note calls it ("855 of 4,008 H37Rv genes carry
#              one"), and what separates it from a lookup that never ran is
#              whether the barcode was fetched: the three early returns that skip
#              the fetch carry no `requests`.
#   free       `pubmed_abstracts` and `corpus_search` take a query, not a gene.
#              Nothing is resolved, and empty means empty.
GATES = {
    "kegg_pathways": "identity", "string_partners": "identity",
    "uniprot_protein": "identity", "resistance_variants": "identity",
    "lineage_markers": "performed",
    "pubmed_abstracts": "free", "corpus_search": "free",
}


def _answered(calls: list[dict[str, Any]], tool: str) -> list[dict[str, Any]]:
    """The calls of `tool` that answered the question, by that tool's own rule."""
    gate = GATES.get(tool, "identity")
    if gate == "free":
        return list(calls)
    if gate == "performed":
        return [c for c in calls if c.get("result", {}).get("requests")]
    return [c for c in calls
            if c.get("result", {}).get("resolved", {}).get("matched_by") != "none"]


def _count(calls: list[dict[str, Any]], tool: str) -> int | str:
    """Distinct record IDs, or NA when the question was never answered.

    NA and 0 are different findings and the whole table turns on keeping them
    apart: NA is "this source never answered for this gene", 0 is "it answered
    and holds nothing", which for KEGG is the ordinary case at 71% of genes.
    """
    if not calls:
        return NA
    answered = _answered(calls, tool)
    if not answered:
        return NA
    ids: set[str] = set()
    for call in answered:
        ids.update(call.get("result", {}).get("record_ids", []))
    return len(ids)


def _truncated(calls: list[dict[str, Any]], default_limit: int = 20) -> bool:
    """Did any partner list come back full? A full list means the true degree is
    unknown, so an absent partner is a property of `limit`, not of the network."""
    for call in _answered(calls, "string_partners"):
        limit = call.get("arguments", {}).get("limit") or default_limit
        if len(call.get("result", {}).get("record_ids", [])) >= int(limit):
            return True
    return False


def _join(values: Any) -> str:
    items = [str(v).strip() for v in (values or []) if str(v).strip()]
    return LIST_SEPARATOR.join(items) if items else NA


def _flag(value: Any) -> str:
    return NA if value is None else str(bool(value))


def gene_row(run: Run) -> dict[str, Any] | None:
    """One row for a single-gene run."""
    if run.mode != "single":
        return None
    query = str(run.output.get("gene", "")).strip()
    identity = _identity(run)
    keys = _alias_keys(identity, query)

    resolved_id = next((identity[f] for f in KEY_FIELDS if identity.get(f)), NA)

    matched_by = identity.get("matched_by") or {}
    rejected = identity.get("rejected") or {}

    kegg = _calls_for(run, "kegg_pathways", keys)
    string = _calls_for(run, "string_partners", keys)
    uniprot = _calls_for(run, "uniprot_protein", keys)
    lineage = _calls_for(run, "lineage_markers", keys)
    resistance = _calls_for(run, "resistance_variants", keys)

    # The resistance flag is about the GENE and lives in `resolved`, not in the
    # variant records: a gene with catalogued variants none of which are graded
    # associated is a real negative, and counting records would call it positive.
    associated: Any = None
    drugs: list[str] = []
    in_catalogue: Any = NA
    for call in _answered(resistance, "resistance_variants"):
        block = call.get("result", {}).get("resolved", {})
        if "resistance_associated" in block:
            associated = bool(block["resistance_associated"]) or bool(associated)
            drugs = sorted(set(drugs) | set(block.get("drugs") or []))
        if "variants_in_catalogue" in block:
            in_catalogue = block["variants_in_catalogue"]

    validation = run.output.get("validation") or {}
    citations = validation.get("citations") or []
    quotes = validation.get("quotes") or []
    flags = sorted({c["status"] for c in citations if c.get("status") != "verified"}
                   | {q["status"] for q in quotes if q.get("status") != "verified"})

    return {
        "query": query,
        "resolved_id": resolved_id,
        "locus_tag": identity.get("locus_tag") or NA,
        "preferred_name": identity.get("preferred_name") or NA,
        "accession": identity.get("accession") or NA,
        "aliases": _join(identity.get("aliases")),
        "n_aliases": len(identity.get("aliases") or []),
        # alias:reason pairs, so a refusal is auditable from the table alone. A
        # silently missing alias is the same class of bug as a silently missing
        # edge, and at a hundred rows nobody opens the store to find out.
        "aliases_rejected": _join(f"{alias}:{reason}" for alias, reason in sorted(rejected.items())),
        "matched_by_string": matched_by.get("string") or NA,
        "matched_by_kegg": matched_by.get("kegg") or NA,
        "matched_by_uniprot": matched_by.get("uniprot") or NA,
        "kegg_n": _count(kegg, "kegg_pathways"),
        "string_n": _count(string, "string_partners"),
        "string_truncated": (str(_truncated(string))
                            if _answered(string, "string_partners") else NA),
        "uniprot_n": _count(uniprot, "uniprot_protein"),
        "pubmed_n": _count(_calls_for(run, "pubmed_abstracts", keys), "pubmed_abstracts"),
        "corpus_n": _count(_calls_for(run, "corpus_search", keys), "corpus_search"),
        "lineage_n": _count(lineage, "lineage_markers"),
        "resistance_associated": _flag(associated),
        "resistance_drugs": _join(drugs),
        # Records are the ASSOCIATED variants only -- katG lists 1,771 in the
        # catalogue of which 1,254 are graded "Uncertain significance", so one
        # count would read as the gene's total and overstate it tenfold.
        "variants_associated_n": _count(resistance, "resistance_variants"),
        "variants_in_catalogue": in_catalogue,
        "citations_verified": sum(1 for c in citations if c.get("status") == "verified"),
        "citations_total": len(citations),
        "quotes_verified": sum(1 for q in quotes if q.get("status") == "verified"),
        "quotes_total": len(quotes),
        "validation_passed": _flag(validation.get("passed")),
        "validation_flags": _join(flags),
        "turns": run.output.get("turns", NA),
        "stop_reason": run.output.get("stop_reason") or NA,
        "run_id": run.run_id or NA,
        "store": str(run.path),
    }


def _record_rows(run: Run, tool: str, columns: list[str],
                 detail_map: dict[str, str]) -> list[dict[str, Any]]:
    """One row per record a per-gene tool returned, across every gene in the run.

    Epistasis runs call these tools for each of their genes, so rows are keyed on
    the gene the call was made for rather than on a single run target -- the same
    reason `store.per_target` exists.
    """
    identities = (run.derived.get("identities") or {}).get("identities") or {}
    rows: list[dict[str, Any]] = []
    for entry in run.calls:
        if entry.get("tool") != tool:
            continue
        query = str(entry.get("arguments", {}).get("gene", "")).strip()
        identity = identities.get(query) or {}
        if not identity:
            # Matched by alias rather than by the exact spelling the call used.
            for candidate in identities.values():
                if query.lower() in {str(a).lower() for a in candidate.get("aliases", [])}:
                    identity = candidate
                    break
        resolved_id = next((identity[f] for f in KEY_FIELDS if identity.get(f)), NA)
        for record in entry.get("result", {}).get("records", []):
            detail = record.get("detail", {})
            row = {"run_id": run.run_id or NA, "resolved_id": resolved_id,
                   "gene_query": query or NA,
                   columns[3]: record.get("record_id", NA), "store": str(run.path)}
            for column, key in detail_map.items():
                value = detail.get(key)
                row[column] = NA if value in (None, "") else value
            rows.append(row)
    return rows


def variant_rows(run: Run) -> list[dict[str, Any]]:
    """One row per resistance-associated variant found in the loci of this run.

    Only graded-associated variants are records, which is the tool's deliberate
    narrowing -- katG has 1,771 catalogued rows of which 1,254 are "Uncertain
    significance" -- so `variants_in_catalogue` in genes.tsv is the denominator
    for these rows, and their absence is not the absence of catalogued variants.
    """
    return _record_rows(run, "resistance_variants", VARIANT_COLUMNS,
                        {"gene": "gene", "mutation": "mutation", "drug": "drug",
                         "confidence": "confidence", "associated": "associated",
                         "source": "source", "comment": "comment"})


def marker_rows(run: Run) -> list[dict[str, Any]]:
    """One row per lineage-defining position falling within the loci of this run.

    `position` is the H37Rv coordinate to compare your own variant against. A row
    here says the GENE can carry a marker, never that a particular variant in it
    is one -- 855 of 4,008 genes contain at least one.
    """
    return _record_rows(run, "lineage_markers", MARKER_COLUMNS,
                        {"position": "position", "lineage": "lineage",
                         "lineage_name": "lineage_name", "allele": "allele"})


def _marker_count(pair: dict[str, Any], which: str) -> Any:
    """NA means the barcode was never consulted for that gene; 0 means it was and
    the gene contains no lineage-defining position -- an informative negative,
    since 855 of 4,008 genes contain one."""
    value = (pair.get("lineage_markers") or {}).get(pair.get(which, ""))
    return NA if value is None else value


def _resistance_cell(pair: dict[str, Any], which: str) -> str:
    """Three states, not two. `NA` is absent from the WHO catalogue and therefore
    never assessed; `not_associated` is assessed and negative; anything else is
    the drugs the gene has an associated variant for."""
    drugs = (pair.get("resistance_drugs") or {}).get(pair.get(which, ""), None)
    if drugs is None:
        return NA
    return LIST_SEPARATOR.join(drugs) if drugs else "not_associated"


def pair_rows(run: Run) -> list[dict[str, Any]]:
    """One row per unordered pair of an epistasis run: N genes give N(N-1)/2.

    Read from `pair_evidence`, which the pipeline computed before the model saw
    anything -- the verdict in this column is the deterministic one the model was
    forbidden to contradict, not the model's reading of it.
    """
    if run.mode != "epistasis":
        return []
    pairs = (run.derived.get("pair_evidence") or {}).get("pairs")
    if pairs is None:
        pairs = run.output.get("pairs") or []
    organism = run.output.get("organism") or NA
    unresolved_run = set(run.output.get("unresolved") or [])

    rows = []
    for pair in pairs:
        direct = pair.get("direct_interaction") or {}
        shared = pair.get("shared_pathways") or []
        retrieved = pair.get("partners_retrieved") or {}
        # An unresolved gene is not a zero: its absence from the pair's evidence
        # is a failed lookup, and the verdict string says so. Carry it as a column
        # so a consumer can filter without parsing that sentence.
        unresolved = sorted(set(pair.get("unresolved") or []) | (
            unresolved_run & {pair.get("gene_a", ""), pair.get("gene_b", "")}))
        rows.append({
            "run_id": run.run_id or NA,
            "organism": organism,
            "gene_a": pair.get("gene_a", NA),
            "gene_b": pair.get("gene_b", NA),
            "verdict": pair.get("verdict", NA),
            # Only a direct /network query can support "no interaction"; a pair
            # inferred from two truncated partner lists cannot.
            "checked_directly": _flag(pair.get("checked_directly")),
            "direct_partner_id": direct.get("partner_id") or NA,
            "direct_score": direct.get("combined_score", NA),
            "textmining_score": direct.get("textmining_score", NA),
            "max_non_textmining_score": direct.get("max_non_textmining_score", NA),
            "evidence_beyond_textmining": _flag(direct.get("evidence_beyond_textmining")
                                                if direct else None),
            "n_shared_pathways": len(shared),
            # Split by specificity rather than counted together: mtu01100 holds
            # 698 genes and sharing it is a base rate, so one column mixing it
            # with mtu00983 (11 genes) is the base-rate trap in tabular form.
            "shared_specific": _join(p["pathway_id"] for p in shared
                                     if p.get("specificity") == "specific"),
            "shared_broad": _join(p["pathway_id"] for p in shared
                                  if p.get("specificity") == "broad"),
            "n_shared_partners": len(pair.get("shared_partners") or []),
            "shared_partners": _join(p.get("record_id") for p in
                                     (pair.get("shared_partners") or [])),
            "partners_retrieved_a": retrieved.get(pair.get("gene_a", ""), NA),
            "partners_retrieved_b": retrieved.get(pair.get("gene_b", ""), NA),
            "truncated": _join(pair.get("truncated")),
            "unresolved": _join(unresolved),
            "lineage_a_n": _marker_count(pair, "gene_a"),
            "lineage_b_n": _marker_count(pair, "gene_b"),
            "shared_lineages": _join(pair.get("shared_lineages")),
            "resistance_a": _resistance_cell(pair, "gene_a"),
            "resistance_b": _resistance_cell(pair, "gene_b"),
            "shared_drugs": _join(pair.get("shared_drugs")),
            "confounds": _join(pair.get("confound_kinds")),
            "store": str(run.path),
        })
    return rows


def write_tsv(path: Path, columns: list[str], rows: list[dict[str, Any]]) -> Path:
    """TSV rather than CSV: gene names, pathway names and verdict sentences all
    contain commas, and a verdict is a whole sentence. Tabs need no quoting here
    because nothing upstream can contain one -- newlines are stripped for the
    same reason."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as fh:
        fh.write("\t".join(columns) + "\n")
        for row in rows:
            fh.write("\t".join(
                str(row.get(column, NA)).replace("\t", " ").replace("\n", " ")
                for column in columns) + "\n")
    return path


def build(runs_dir: Path, out_dir: Path) -> dict[str, Path]:
    """Table every finished run under `runs_dir`.

    All four files are always written, empty but for their header when a mode or a
    tool is absent: a missing file is ambiguous between "no epistasis runs" and
    "the step did not run", and a downstream rule should not have to tell them
    apart.
    """
    genes: list[dict[str, Any]] = []
    pairs: list[dict[str, Any]] = []
    variants: list[dict[str, Any]] = []
    markers: list[dict[str, Any]] = []
    for path in sorted(runs_dir.rglob("*.jsonl")):
        if path.name.endswith(".corpus.jsonl"):
            continue
        run = load_run(path)
        if run is None:
            continue
        row = gene_row(run)
        if row:
            genes.append(row)
        pairs.extend(pair_rows(run))
        variants.extend(variant_rows(run))
        markers.extend(marker_rows(run))

    genes.sort(key=lambda r: (str(r["resolved_id"]), str(r["query"])))
    pairs.sort(key=lambda r: (str(r["run_id"]), str(r["gene_a"]), str(r["gene_b"])))
    variants.sort(key=lambda r: (str(r["resolved_id"]), str(r["drug"]), str(r["mutation"])))
    markers.sort(key=lambda r: (str(r["resolved_id"]), r["position"]
                                if isinstance(r["position"], int) else 0))
    return {
        "genes": write_tsv(out_dir / "genes.tsv", GENE_COLUMNS, genes),
        "pairs": write_tsv(out_dir / "pairs.tsv", PAIR_COLUMNS, pairs),
        "variants": write_tsv(out_dir / "resistance_variants.tsv", VARIANT_COLUMNS, variants),
        "markers": write_tsv(out_dir / "lineage_markers.tsv", MARKER_COLUMNS, markers),
    }
