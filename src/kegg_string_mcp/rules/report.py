"""The pre-flight report: does this vocabulary resolve against this annotation?

Run before anything costs money. Every later stage assumes a feature means one
locus or one interval, so a vocabulary that does not resolve cleanly is a finding
about the inputs rather than a step to push past.
"""

from __future__ import annotations

from pathlib import Path

from kegg_string_mcp.rules.annotation import Annotation
from kegg_string_mcp.rules.features import CODING, INTERGENIC, Coverage

FEATURE_COLUMNS = ["label", "kind", "matched_by", "locus", "symbol", "start", "end",
                   "strand", "length", "left", "right", "width", "promoter_of",
                   "candidates", "note"]

NA = "NA"


def render(coverage: Coverage, annotation: Annotation, limit: int = 15) -> str:
    counts = coverage.counts()
    total = max(counts["total"], 1)
    resolved = counts[CODING] + counts[INTERGENIC]
    lines = [
        f"annotation: {annotation.path}",
        (f"  sha256 {annotation.sha256[:16]}…  format {annotation.source_format.upper()}  "
         f"{len(annotation.genes):,} genes  {len(annotation.intergenic_names):,} intergenic"),
        *(f"  note: {n}" for n in annotation.notes),
        "",
        f"labels: {counts['total']:,}",
        f"  resolved                {resolved:,}  ({resolved / total:.1%})",
        f"    coding                {counts[CODING]:,}",
        f"      exact               {counts['exact']:,}",
        f"      via strand suffix   {counts['strand_suffix']:,}",
        f"      via symbol          {counts['symbol']:,}",
        f"    intergenic            {counts[INTERGENIC]:,}",
        f"      via flanking pair   {counts['flanking']:,}",
        f"      5' of a gene        {counts['upstream']:,}",
        f"      3' of a gene        {counts['downstream']:,}",
        f"  unresolved              {counts['unresolved']:,}",
        f"    of which ambiguous    {counts['ambiguous']:,}  (refused, not guessed)",
    ]
    unresolved = coverage.unresolved()
    if unresolved:
        lines += ["", "unresolved labels:"]
        for feature in unresolved[:limit]:
            extra = f"  candidates: {', '.join(feature.candidates)}" if feature.candidates else ""
            lines.append(f"  {feature.label:22} {feature.note}{extra}")
        if len(unresolved) > limit:
            lines.append(f"  … and {len(unresolved) - limit:,} more")
    return "\n".join(lines)


def write_tsv(path: Path, coverage: Coverage) -> Path:
    """One row per label. Tabs, because a note is a sentence and contains commas."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write("\t".join(FEATURE_COLUMNS) + "\n")
        for label in coverage.features:
            row = coverage.features[label].to_dict()
            handle.write("\t".join(
                str(row.get(column, NA)).replace("\t", " ").replace("\n", " ")
                for column in FEATURE_COLUMNS) + "\n")
    return path


# --- the two output tables -------------------------------------------------
#
# One row per locus and one row per rule, and nothing else. A supplementary
# table is read by people, so relationships live in a structured string column
# rather than in a third file: `ahpC:string_textmining_only:0.968|furA:adjacent:6bp`
# stays readable and still parses.

LOCUS_COLUMNS = [
    "locus", "symbol", "labels", "product", "start", "end", "strand", "length",
    "catalogue_status", "catalogue_drugs", "catalogued_variants", "associated_variants",
    "lineage_markers", "kegg_pathways", "string_partners", "linked_loci",
    "n_rules", "n_present", "n_absent", "with_anchor", "without_anchor", "notes",
]

RULE_COLUMNS = [
    "rule_id", "row", "k", "predicted_class", "conditions", "roles",
    "primary_signature", "signatures", "anchors", "compensator_candidates",
    "unknown_loci", "negated_anchors", "unresolved", "links",
    "shared_lineages", "shared_drugs", "verdict", "problems",
]


def _join(values, limit: int | None = None) -> str:
    # Materialise first: callers pass generators, and the truncation branch needs
    # to count what it is truncating.
    items = [str(v) for v in values if str(v).strip()]
    if limit is not None and len(items) > limit:
        items = items[:limit] + [f"…+{len(items) - limit}"]
    return "|".join(items) if items else NA


def locus_row(annotation, links) -> dict:
    """`links` is this locus's entries, already rendered."""
    gene = annotation.gene
    catalogue = annotation.catalogue
    return {
        "locus": annotation.locus,
        "symbol": annotation.symbol or NA,
        # Which spellings in the rule file reached this locus. A reader comparing
        # the supplementary table to the rules needs to see that `Rv0006c` and
        # `Rv0006` are one row, not two.
        "labels": _join(annotation.labels),
        "product": annotation.product or NA,
        "start": gene.start if gene else NA,
        "end": gene.end if gene else NA,
        "strand": gene.strand if gene else NA,
        "length": gene.length if gene else NA,
        "catalogue_status": catalogue.status if catalogue else NA,
        "catalogue_drugs": _join(catalogue.drugs) if catalogue else NA,
        "catalogued_variants": catalogue.catalogued if catalogue else NA,
        "associated_variants": catalogue.associated if catalogue else NA,
        "lineage_markers": _join(annotation.lineages),
        "kegg_pathways": _join(f"{pid}:{name}" for pid, name in annotation.pathways),
        "string_partners": _join((f"{name or pid}:{score:.3f}"
                                  for pid, name, score in annotation.partners), limit=10),
        "linked_loci": _join(links),
        "n_rules": annotation.n_rules,
        "n_present": annotation.n_present,
        "n_absent": annotation.n_absent,
        # Meaningless for an anchor: the counts ask "did this locus appear beside a
        # known resistance locus", and a known resistance locus trivially did.
        # Left blank rather than printed as zero, which would read as "never".
        "with_anchor": (NA if catalogue and catalogue.is_anchor else annotation.with_anchor),
        "without_anchor": (NA if catalogue and catalogue.is_anchor
                           else annotation.without_anchor),
        "notes": _join(annotation.notes),
    }


def write_loci(path: Path, annotations, links_by_locus) -> Path:
    rows = [locus_row(a, links_by_locus.get(a.locus, []))
            for a in sorted(annotations, key=lambda a: (a.gene.start if a.gene else 0, a.locus))]
    return _write(path, LOCUS_COLUMNS, rows)


def write_rules(path: Path, signatures, rename: dict[str, str] | None = None) -> Path:
    """Rules, with the learner's own columns kept and prefixed.

    The `roles` column is what stops a rule row needing a third table: it says
    which locus played which part, in the same conjunction order as `conditions`.
    """
    rows = []
    extra_columns: list[str] = []
    for signature in signatures:
        row = signature.to_dict(rename)
        row["roles"] = _join(f"{c.locus}:{c.role}" for c in signature.conditions)
        rows.append(row)
        for key in row:
            if key.startswith("scan_") and key not in extra_columns:
                extra_columns.append(key)
    return _write(path, RULE_COLUMNS + sorted(extra_columns), rows)


def _write(path: Path, columns: list[str], rows: list[dict]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write("\t".join(columns) + "\n")
        for row in rows:
            handle.write("\t".join(
                str(row.get(column, NA)).replace("\t", " ").replace("\n", " ")
                for column in columns) + "\n")
    return path


QUESTION_COLUMNS = ["question_id", "kind", "locus", "partner", "drugs", "n_rules",
                    "raised_by", "signatures", "question"]


def write_questions(path: Path, questions) -> Path:
    """Stage A of the literature layer, and the only free part of it.

    Written whatever the run found, empty header included: no questions is a
    result -- every rule accounted for by the structured sources -- and a missing
    file is ambiguous between that and a step that did not run.
    """
    return _write(path, QUESTION_COLUMNS, [q.to_dict() for q in questions])
