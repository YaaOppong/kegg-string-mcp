"""Lineage-marker annotation: the confound that has to be ruled out first.

An epistasis scan over clinical isolates is confounded by population structure as
the rule, not the exception. If two genes each carry a SNP that defines the same
part of the M. tuberculosis phylogeny, every isolate in that clade has both
variants, and the scan reports an association with no biological interaction
behind it. A hypothesis generator that only proposes mechanisms will dress that
up as a discovery.

The lineage-defining positions come from TB-Profiler's `tbdb/barcode.bed`, the
machine-readable form of the Coll (2014) and Napier (2020) SNP barcodes: 1,111
positions in H37Rv coordinates, each labelled with the lineage it defines. It is
public, versioned in git, and needs no credentials. Gene spans come from KEGG's
organism list, which the repo already fetches and caches.

**The direction follows from set containment, not from statistics.** A barcode
SNP is derived inside its clade and ancestral everywhere else, so for two markers
defining clades L_a and L_b there are only two cases -- two nodes of a phylogeny
are either nested or disjoint, never anything else:

* **nested** (L_b subset of L_a, e.g. lineage4 and lineage4.6) -- every isolate
  carrying B also carries A, so the A-B+ cell of the 2x2 is structurally empty
  and the pair appears positively associated;
* **disjoint** (e.g. lineage4.6.1 and lineage4.6.3, or lineage2 and lineage4) --
  no isolate carries both, so the A+B+ cell is structurally empty and the pair
  appears negatively associated.

Phylogenetic distance is not the axis. Sister sublineages and entirely separate
lineages are equally disjoint and confound identically; only nesting reverses the
sign. Both are population structure rather than biology.

**A gene-level flag is a prior, not a verdict.** 853 of 4,008 H37Rv genes contain
at least one barcode position, so a flag here says the gene is capable of
carrying a lineage marker, not that the variant in your hit is one. Position-level
matching is available via `marks_position`, and the verdict comes from your own
genotype matrix: does the association survive conditioning on lineage? No
database answers that.
"""

from __future__ import annotations

import re
from dataclasses import asdict, dataclass, field
from typing import Any

BARCODE_URL = "https://raw.githubusercontent.com/jodyphelan/tbdb/master/barcode.bed"
KEGG_GENE_LIST = "https://rest.kegg.jp/list/mtu"

# KEGG location column: "1..1524", "complement(2052..3260)", "join(a..b,c..d)".
# The outermost span is what matters for containment; strand and exon structure
# do not change which positions fall in the gene.
_SPAN = re.compile(r"(\d+)\.\.(\d+)")


@dataclass(frozen=True)
class LineageSnp:
    position: int
    lineage: str          # e.g. "lineage4.6.1", "La1.8", "M.canetti"
    lineage_name: str     # e.g. "Euro-American"
    allele: str = ""

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass
class GeneSpan:
    locus: str
    symbol: str | None
    start: int
    end: int


def parse_barcode(text: str) -> list[LineageSnp]:
    """BED: chrom, start, end, lineage, allele, lineage_name, spoligotype, RD."""
    out = []
    for line in text.splitlines():
        row = line.rstrip("\n").split("\t")
        if len(row) < 6 or not row[2].isdigit():
            continue
        out.append(LineageSnp(position=int(row[2]), lineage=row[3],
                              allele=row[4], lineage_name=row[5]))
    return out


def parse_gene_spans(text: str) -> list[GeneSpan]:
    out = []
    for line in text.splitlines():
        row = line.rstrip("\n").split("\t")
        if len(row) < 4:
            continue
        spans = _SPAN.findall(row[2])
        if not spans:
            continue
        locus = row[0].split(":", 1)[1] if ":" in row[0] else row[0]
        symbol = row[3].split(";")[0].strip() if ";" in row[3] else None
        out.append(GeneSpan(locus=locus, symbol=symbol,
                            start=int(spans[0][0]), end=int(spans[-1][1])))
    return out


def annotate(spans: list[GeneSpan], snps: list[LineageSnp]) -> dict[str, list[LineageSnp]]:
    """Gene symbol (or locus, when unnamed) -> the lineage SNPs it contains.

    Sorted-span bisection would be faster; at 4,000 genes by 1,111 SNPs the
    quadratic form runs in well under a second and is easier to be sure of.
    """
    index: dict[str, list[LineageSnp]] = {}
    for snp in snps:
        for span in spans:
            if span.start <= snp.position <= span.end:
                index.setdefault(span.symbol or span.locus, []).append(snp)
                break
    return index


def marks_position(snps: list[LineageSnp], position: int) -> LineageSnp | None:
    """Position-level check: is this exact variant a lineage marker?

    The gene-level flag is a prior; this is the question that decides whether a
    specific hit is population structure.
    """
    for snp in snps:
        if snp.position == position:
            return snp
    return None


def _levels(lineage: str) -> list[str]:
    return lineage.split(".")


def relate(a: str, b: str) -> str:
    """How two lineage labels co-segregate: nested | disjoint.

    Two cases only, because two nodes of a phylogeny are either nested or
    disjoint. An earlier version returned a third value for markers on different
    top-level lineages, which conflated phylogenetic distance with
    co-segregation: lineage2 and lineage4 are disjoint clades, so their markers
    exclude each other exactly as sister sublineages do, and confound in the same
    negative direction.
    """
    la, lb = _levels(a), _levels(b)
    shorter, longer = (la, lb) if len(la) <= len(lb) else (lb, la)
    # Compared level by level, so lineage4 is not read as an ancestor of lineage41.
    return "nested" if longer[:len(shorter)] == shorter else "disjoint"


@dataclass
class PairLineageFlag:
    gene_a: str
    gene_b: str
    snps_a: list[LineageSnp] = field(default_factory=list)
    snps_b: list[LineageSnp] = field(default_factory=list)
    relations: list[tuple[str, str, str]] = field(default_factory=list)  # lin_a, lin_b, relation

    @property
    def risk(self) -> str:
        """One risk per pair, but a gene can carry markers for several clades.

        The scan tests variant pairs, so each (marker_a, marker_b) pairing has its
        own direction. Reporting only the strongest would hide that a gene pair
        can be confounded both ways depending on which variants were called, so a
        mixed set is named rather than collapsed.
        """
        if not self.snps_a or not self.snps_b:
            return "none"
        kinds = {rel for _, _, rel in self.relations}
        if kinds == {"nested"}:
            return "confounding_positive"
        if kinds == {"disjoint"}:
            return "confounding_negative"
        return "confounding_mixed"

    @property
    def note(self) -> str:
        return {
            "none": "at most one gene carries a lineage-defining SNP",
            "confounding_mixed": ("the two genes carry markers for several clades, some nested "
                                  "and some disjoint, so the direction depends on which variants "
                                  "the scan called; resolve at position level before interpreting"),
            "confounding_positive": ("both genes mark nested clades, so isolates in the inner "
                                     "clade carry both variants and will appear associated "
                                     "regardless of any interaction -- condition on lineage "
                                     "before treating this pair as a finding"),
            "confounding_negative": ("the genes mark disjoint clades, which no isolate shares, so "
                                     "a negative association here is population structure"),
        }[self.risk]

    def to_dict(self) -> dict[str, Any]:
        return {"gene_a": self.gene_a, "gene_b": self.gene_b, "risk": self.risk,
                "note": self.note,
                "snps_a": [s.to_dict() for s in self.snps_a],
                "snps_b": [s.to_dict() for s in self.snps_b],
                "relations": [list(r) for r in self.relations]}


def flag_pair(a: str, b: str, index: dict[str, list[LineageSnp]]) -> PairLineageFlag:
    snps_a, snps_b = index.get(a, []), index.get(b, [])
    relations = []
    for lin_a in sorted({s.lineage for s in snps_a}):
        for lin_b in sorted({s.lineage for s in snps_b}):
            relations.append((lin_a, lin_b, relate(lin_a, lin_b)))
    return PairLineageFlag(gene_a=a, gene_b=b, snps_a=snps_a, snps_b=snps_b,
                           relations=relations)


def load(http: Any) -> dict[str, list[LineageSnp]]:
    """Fetch both sources through the repo's polite, caching client.

    Routed through PoliteClient rather than a bare request so the barcode lands in
    the content-addressed cache with an audit URL, like every other retrieval the
    repo makes. A lineage call that cannot be traced to a specific version of the
    barcode is not evidence.
    """
    barcode = parse_barcode(http.get(BARCODE_URL).body)
    spans = parse_gene_spans(http.get(KEGG_GENE_LIST).body)
    return annotate(spans, barcode)
