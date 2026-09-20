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

**One row per locus, whatever the file says.** A GFF3 carries a parent `gene`
row and a child `CDS` row for the same locus tag, so reading both as genes puts
every tag in the index twice. Nothing looks broken -- until a label like
`Rv0003c` has to reach `Rv0003` through the strand-suffix rule, finds two
candidates that are the same gene, and is refused as ambiguous. Measured on
NCBI's H37Rv GFF3: 7,884 "genes" parsed, `via strand suffix` matches fall from
1,339 to zero, and resolution halves from 97% to 51%. The overlap count inflates
from 940 to 4,846, so intergenic naming goes with it.

Gene-level rows (`gene`, `pseudogene`) are the locus set. On that file they cover
all 4,008 distinct locus tags -- 3,978 `gene` plus 30 `pseudogene` -- and carry
every biotype, so rRNA loci like `rrs` and `rrl` come in without special-casing.
A file with no gene-level rows at all (Prokka output, some Ensembl bacterial
dumps) falls back to `CDS` and says so.

**The CDS start is kept separately, and it is not decoration.** HGVS `c.-N`
numbering is relative to the translation start, not to the gene's 5' end. Those
differ on three H37Rv loci -- Rv0614's gene row begins 243 bp before its CDS --
so placing a promoter variant from the gene start would put it 243 bp from where
it is. The gene span is what a locus occupies; the CDS start is where `c.`
counts from.

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


# SnpEff appends a transcript version to the gene part of a region name --
# `upstream_Rv1482c.1`. Locus tags and gene symbols carry no dot, so a trailing
# `.<digits>` is a version and nothing else.
_VERSION = re.compile(r"\.\d+$")


def strip_version(name: str) -> str:
    return _VERSION.sub("", name.strip())


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
    # Translation start, where the annotation distinguishes it from the gene's
    # 5' end. `c.` coordinates count from here, not from `start`.
    cds_start: int | None = None

    @property
    def length(self) -> int:
        return self.end - self.start + 1

    @property
    def coding_start(self) -> int:
        """Where HGVS `c.` numbering counts from, on this gene's own strand."""
        if self.cds_start is not None:
            return self.cds_start
        return self.end if self.strand == "-" else self.start


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
    # Intervals by the locus on each side, for resolving `upstream_X`.
    _by_left: dict[str, Intergenic] = field(default_factory=dict, repr=False)
    _by_right: dict[str, Intergenic] = field(default_factory=dict, repr=False)

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
                interval = Intergenic(name=name, left=left, right=right,
                                      start=left.end + 1, end=right.start - 1)
                self._intergenic[name] = interval
                self._by_left[left.locus] = interval
                self._by_right[right.locus] = interval
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

    def upstream_of(self, gene: Gene) -> Intergenic | None:
        """The interval immediately 5' of a gene, on that gene's own strand.

        None when the neighbour abuts or overlaps it, which is not rare: 830 of
        H37Rv's 4,008 genes have no gap 5' of them. That is a fact about the
        genome and a different answer from "the gene was not found".
        """
        return (self._by_left.get(gene.locus) if gene.strand == "-"
                else self._by_right.get(gene.locus))

    def downstream_of(self, gene: Gene) -> Intergenic | None:
        return (self._by_right.get(gene.locus) if gene.strand == "-"
                else self._by_left.get(gene.locus))

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


# Parent features: one row per locus, covering every biotype. Child rows (CDS,
# exon, tRNA, rRNA) repeat their parent's locus tag and must not be counted again.
GENE_LEVEL = ("gene", "pseudogene")
CODING = ("CDS",)


def parse(path: str | Path, feature_types: tuple[str, ...] = GENE_LEVEL) -> Annotation:
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
    notes_extra: list[str] = []

    cds_starts: dict[str, tuple[int, int]] = {}
    if fmt == "gff":
        present = {row[2] for row in rows if len(row) >= 9}
        if not (set(feature_types) & present):
            # Prokka and some Ensembl bacterial dumps emit CDS rows and no gene
            # rows. Falling back keeps those files usable; saying so keeps the
            # coordinate meaning honest, since a CDS span is the translated
            # region rather than the locus.
            feature_types = CODING
            notes_extra.append(
                f"no {', '.join(GENE_LEVEL)} rows found; fell back to CDS rows, so spans are "
                f"translated regions rather than gene extents")
        # Translation starts, by locus tag, for HGVS `c.` conversion.
        for row in rows:
            if len(row) < 9 or row[2] not in CODING:
                continue
            tag = _attributes(row[8]).get("locus_tag")
            if tag:
                cds_starts.setdefault(tag, (int(row[3]), int(row[4])))

    seen_tags: set[str] = set()
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
            # One row per locus. A second row for a tag already seen is a child
            # feature or a duplicate, and admitting it makes every strand-suffix
            # lookup ambiguous against the gene's own other row.
            if name and name in seen_tags:
                continue
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
        seen_tags.add(name)
        cds = cds_starts.get(name)
        coding = None
        if cds is not None and (cds[0], cds[1]) != (start, end):
            coding = cds[1] if strand == "-" else cds[0]
        genes.append(Gene(locus=name, symbol=symbol, start=start, end=end, strand=strand,
                          cds_start=coding))

    notes = [f"parsed {len(genes)} genes from {fmt.upper()} at {path}", *notes_extra]
    distinct = sum(1 for g in genes if g.cds_start is not None)
    if distinct:
        notes.append(f"{distinct} locus/loci whose CDS start differs from the gene start; "
                     f"`c.` coordinates count from the CDS")
    if skipped:
        notes.append(f"{skipped} row(s) skipped for having no usable name or too few columns")
    if not any(g.strand in "+-" for g in genes):
        # Strand drives the promoter-orientation reading, and BED without a
        # strand column silently makes every gap ambiguous rather than wrong.
        notes.append("no strand column found: promoter orientation cannot be inferred")
    return Annotation(path=str(path), sha256=digest, source_format=fmt,
                      genes=genes, notes=notes)
