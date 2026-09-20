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
