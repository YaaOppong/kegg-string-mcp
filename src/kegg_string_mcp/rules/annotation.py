"""The reference annotation a feature vocabulary was named from.

Supplied by the caller, never fetched. The features in a rule set are whatever
the variant caller called them, and that naming comes from one specific
annotation build -- SnpEff's `Mycobacterium_tuberculosis_h37rv`, say. Resolving
those names against a *different* gene list is how a locus silently becomes the
wrong locus, so the annotation is an input and its identity is recorded.

KEGG is deliberately not that input. Its curated `mtu` list holds 4,008 entries
and is authoritative for what a gene *participates in*; it is not authoritative
for what the caller's feature is called or where it sits. Measured against one
real vocabulary, 76 of 2,879 locus tags were absent from KEGG entirely and 1,346
more differed in their strand suffix. Two sources, two jobs.

**Coordinates are normalised to 1-based inclusive on the way in.** BED is 0-based
half-open and GFF/GTF is 1-based inclusive; mixing them shifts every start by one
and every interval width by one, which is invisible until a variant sits on a
boundary and lands in the wrong feature. The conversion happens once, here.
"""

from __future__ import annotations

import gzip
import hashlib
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

# A trailing `c` is the complement-strand marker in H37Rv nomenclature and may or
# may not be present depending on which annotation named the feature. A trailing
# `a`, `A` or `B` is part of the gene's NAME -- Rv0063a is a different gene from
# Rv0063 -- so stripping any trailing letter manufactures ambiguity that does not
# exist. Across 4,008 loci, stripping only `c` collides zero times; stripping any
# letter collides 96 times.
_STRAND_SUFFIX = re.compile(r"^(Rv\d{4})[cC]$")
_LOCUS_TAG = re.compile(r"^Rv\d{4}[A-Za-z]?$", re.IGNORECASE)


def normalise_locus(tag: str) -> str:
    """Drop the complement-strand suffix, and nothing else."""
    match = _STRAND_SUFFIX.match(tag.strip())
    return match.group(1) if match else tag.strip()


def is_locus_tag(label: str) -> bool:
    return bool(_LOCUS_TAG.match(label.strip()))


@dataclass(frozen=True)
class Gene:
    """One annotated locus, 1-based inclusive."""

    locus: str
    symbol: str
    start: int
    end: int
    strand: str          # "+" | "-" | "." when the annotation does not say

    @property
    def length(self) -> int:
        return self.end - self.start + 1


@dataclass(frozen=True)
class Intergenic:
    """The gap between two consecutive genes, named as the caller names it.

    `name` is `{left}-{right}` in genome order, which is the convention the
    upstream matrix uses. Locus tags contain no hyphen, so the name splits back
    into its two flanking genes unambiguously.
    """

    name: str
    left: Gene
    right: Gene
    start: int
    end: int

    @property
    def width(self) -> int:
        return self.end - self.start + 1

    def promoter_of(self) -> list[Gene]:
        """Which flanking gene(s) this gap sits 5' of, by strand.

        An inference from orientation, not an annotation: a gap upstream of a
        forward-strand gene is where that gene's promoter would be. Both flanks
        qualify when they diverge (`<-- -->`), neither when they converge. The
        caller should report this as orientation, never as "the promoter of X".
        """
        out = []
        if self.right.strand == "+":
            out.append(self.right)
        if self.left.strand == "-":
            out.append(self.left)
        return out


@dataclass
class Annotation:
    """A parsed reference annotation, plus the provenance to compare two runs.

    `sha256` is the point of this class as much as the genes are. Every HTTP
    response in a run store already carries a content hash; a local file that
    decides what every feature means deserves the same, or a reference update
    changes the results with nothing in the record to say so.
    """

    path: str
    sha256: str
    source_format: str                  # "bed" | "gff" | "gtf"
    genes: list[Gene] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)

    _by_locus: dict[str, Gene] = field(default_factory=dict, repr=False)
    _by_norm: dict[str, list[Gene]] = field(default_factory=dict, repr=False)
    _by_symbol: dict[str, list[Gene]] = field(default_factory=dict, repr=False)
    _intergenic: dict[str, Intergenic] = field(default_factory=dict, repr=False)

    def __post_init__(self) -> None:
        self.genes = sorted(self.genes, key=lambda g: (g.start, g.end, g.locus))
        for gene in self.genes:
            self._by_locus.setdefault(gene.locus, gene)
            self._by_norm.setdefault(normalise_locus(gene.locus), []).append(gene)
            if gene.symbol:
                self._by_symbol.setdefault(gene.symbol, []).append(gene)
        self._build_intergenic()

    def _build_intergenic(self) -> None:
        """Gaps between consecutive genes. Overlapping neighbours produce none.

        Genes that overlap have no gap between them, so that junction simply has
        no feature -- which is a fact about the genome, not a parsing failure.
        Around 900 of H37Rv's consecutive pairs overlap.
        """
        overlapping = 0
        for left, right in zip(self.genes, self.genes[1:]):
            if right.start > left.end + 1:
                name = f"{left.locus}-{right.locus}"
                self._intergenic[name] = Intergenic(
                    name=name, left=left, right=right,
                    start=left.end + 1, end=right.start - 1)
            else:
                overlapping += 1
        self.notes.append(
            f"{len(self._intergenic)} intergenic intervals; {overlapping} consecutive "
            f"gene pairs overlap and so have none.")

    # -- lookups -------------------------------------------------------------

    def gene(self, locus: str) -> Gene | None:
        return self._by_locus.get(locus.strip())

    def by_normalised(self, tag: str) -> list[Gene]:
        """Genes matching a locus tag once the strand suffix is ignored."""
        return list(self._by_norm.get(normalise_locus(tag), []))

    def by_symbol(self, symbol: str) -> list[Gene]:
        return list(self._by_symbol.get(symbol.strip(), []))

    def intergenic(self, name: str) -> Intergenic | None:
        return self._intergenic.get(name.strip())

    @property
    def intergenic_names(self) -> list[str]:
        return list(self._intergenic)

    def covering(self, position: int) -> list[Gene]:
        """Genes whose span contains a 1-based position. Overlaps mean a list."""
        return [g for g in self.genes if g.start <= position <= g.end]

    def summary(self) -> dict[str, Any]:
        return {"path": self.path, "sha256": self.sha256, "format": self.source_format,
                "genes": len(self.genes), "intergenic": len(self._intergenic),
                "notes": list(self.notes)}


# -- parsing ----------------------------------------------------------------

def _open(path: Path):
    return gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" \
        else path.open(encoding="utf-8")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


_ATTR = re.compile(r'(\w+)[=\s]"?([^";]+)"?')


def _attributes(field_text: str) -> dict[str, str]:
    return {k: v.strip() for k, v in _ATTR.findall(field_text)}


def _detect(first_row: list[str]) -> str:
    """BED or GFF/GTF, from the shape of a row rather than the file extension.

    A file named `.bed` that is really a GTF would otherwise be read with the
    wrong coordinate convention, which shifts every start by one silently.
    """
    if len(first_row) >= 8 and first_row[3].isdigit() and first_row[4].isdigit():
        return "gff"
    return "bed"


def parse(path: str | Path, feature_types: tuple[str, ...] = ("gene", "CDS")) -> Annotation:
    """Read a BED or GFF/GTF into 1-based inclusive `Gene` records.

    Name resolution order is the caller's vocabulary first: `locus_tag`, then
    `gene_id`, then `ID`, then `Name`/`gene_name`. SnpEff falls back to a gene's
    NAME where the annotation has one and the locus tag otherwise, which is why a
    real vocabulary carries a mix of `Rv0001` and `rpoC` -- both must resolve.
    """
    path = Path(path)
    digest = _sha256(path)
    rows: list[list[str]] = []
    with _open(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            rows.append(line.rstrip("\n").split("\t"))
    if not rows:
        raise ValueError(f"{path} has no usable rows")

    fmt = _detect(rows[0])
    genes: list[Gene] = []
    skipped = 0
    for row in rows:
        if fmt == "gff":
            if len(row) < 9:
                skipped += 1
                continue
            if feature_types and row[2] not in feature_types:
                continue
            attrs = _attributes(row[8])
            name = (attrs.get("locus_tag") or attrs.get("gene_id")
                    or attrs.get("ID") or attrs.get("Name") or attrs.get("gene_name"))
            symbol = attrs.get("gene_name") or attrs.get("Name") or attrs.get("gene") or ""
            start, end, strand = int(row[3]), int(row[4]), row[6] or "."
        else:
            if len(row) < 4:
                skipped += 1
                continue
            name = row[3].strip()
            symbol = row[6].strip() if len(row) > 6 else ""
            # BED is 0-based half-open: [start, end) -> [start+1, end].
            start, end = int(row[1]) + 1, int(row[2])
            strand = row[5].strip() if len(row) > 5 and row[5].strip() in "+-." else "."
        if not name:
            skipped += 1
            continue
        if symbol == name:
            symbol = ""
        genes.append(Gene(locus=name, symbol=symbol, start=start, end=end, strand=strand))

    notes = [f"parsed {len(genes)} genes from {fmt.upper()} at {path}"]
    if skipped:
        notes.append(f"{skipped} row(s) skipped for having no usable name or too few columns")
    if not any(g.strand in "+-" for g in genes):
        # Strand drives the promoter-orientation reading, and BED without a
        # strand column silently makes every gap ambiguous rather than wrong.
        notes.append("no strand column found: promoter orientation cannot be inferred")
    return Annotation(path=str(path), sha256=digest, source_format=fmt,
                      genes=genes, notes=notes)
