"""Lineage-marker lookups: is this gene a known lineage-defining locus?

A standard annotation alongside KEGG pathways, STRING partners and UniProt
function -- not an interpretation. An epistasis scan over clinical isolates is
confounded by population structure, so whether a gene sits on a lineage-defining
site belongs with the other things you look up about it before anything is
proposed. This tool reports whether it does and which lineages the markers
define, and stops there.

Positions come from TB-Profiler's `tbdb/barcode.bed`, the machine-readable form
of the Coll (2014) and Napier (2020) SNP barcodes: 1,111 positions in H37Rv
coordinates, public, versioned in git, no credentials. Gene spans come from
KEGG's organism list, which the repo already fetches and caches. Both go through
the polite client, so a lineage call traces to a specific barcode version -- an
annotation that cannot name the barcode it came from is not evidence.

Like a KEGG pathway ID and unlike a PubMed abstract, a marker is **structured**:
the record means one thing and there is nothing to quote from it. Records are
therefore citable but not quotable, and carry no `quotable_text`.

A gene-level flag is coarse: 855 of 4,008 H37Rv genes contain at least one
barcode position, so a flag says the gene can carry a marker, not that the
variant in a given hit is one. `marks_position` answers the variant-level
question, and every record carries its exact H37Rv position so a caller can.
"""

from __future__ import annotations

import re
from dataclasses import asdict, dataclass
from typing import Any

from kegg_string_mcp.http import FetchError, PoliteClient
from kegg_string_mcp.provenance import Record, RequestTrace, ToolResult

BARCODE_URL = "https://raw.githubusercontent.com/jodyphelan/tbdb/master/barcode.bed"
KEGG_GENE_LIST = "https://rest.kegg.jp/list/{organism}"
# The barcode has no per-position landing page, so records resolve to the file
# itself. It is a few hundred kB of TSV and a human can find the row by position.
BARCODE_SOURCE = "https://github.com/jodyphelan/tbdb"
# The barcode is a list of H37Rv coordinates and nothing in the file says so.
# Matched against another organism's gene spans it produces confident nonsense:
# 789 of 4,639 E. coli genes "carry" a TB lineage marker by coordinate collision
# alone. Only the organism whose coordinates the barcode is written in is
# accepted; every other one is refused rather than answered.
BARCODE_ORGANISM = "mtu"

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
    """Gene identifier -> the lineage SNPs it contains.

    Keyed under BOTH the locus tag and the symbol, so Rv0757 and phoP resolve to
    the same markers. Callers arrive with whichever identifier their pipeline
    uses, and keying on the symbol alone silently returned "no marker" for every
    locus tag -- an informative-looking negative that was an indexing artefact.
    A gene therefore occupies two keys: len(index) is not a gene count.

    A position inside overlapping CDSs is recorded against every gene containing
    it, not the first one in file order. Overlapping reading frames are ordinary
    in a compact bacterial genome -- 3 of the 1,111 barcode positions sit in two
    genes each -- and stopping at the first match made the second gene look
    unmarked depending on how KEGG happened to order its output.

    Sorted-span bisection would be faster; at 4,000 genes by 1,111 SNPs the
    quadratic form runs in well under a second and is easier to be sure of.
    """
    index: dict[str, list[LineageSnp]] = {}
    for snp in snps:
        for span in spans:
            if span.start <= snp.position <= span.end:
                for key in filter(None, (span.locus, span.symbol)):
                    index.setdefault(key, []).append(snp)
    return index


def is_marker(gene: str, index: dict[str, list[LineageSnp]]) -> bool:
    return bool(index.get(gene))


def lineages(gene: str, index: dict[str, list[LineageSnp]]) -> list[str]:
    """Which lineages the gene's markers define, if any."""
    return sorted({snp.lineage for snp in index.get(gene, [])})


def marks_position(snps: list[LineageSnp], position: int) -> LineageSnp | None:
    """Variant-level check: is this exact position a lineage-defining site?

    The gene-level flag is coarse; this is the question that applies to a
    specific hit.
    """
    for snp in snps:
        if snp.position == position:
            return snp
    return None


def _trace(resp: Any) -> RequestTrace:
    return RequestTrace(url=resp.audit_url, retrieved_at=resp.fetched_at, cached=resp.cached,
                        status=resp.status, content_sha256=resp.content_sha256)


class LineageClient:
    """Gene -> lineage-defining SNPs it contains.

    The barcode and the gene spans are fetched once per client and reused. A
    per-gene fetch would re-parse 1,111 positions and 4,008 spans for every
    lookup, and the underlying responses are identical every time.
    """

    def __init__(self, http: PoliteClient):
        self.http = http
        self._index: dict[str, dict[str, list[LineageSnp]]] = {}
        self._traces: dict[str, list[RequestTrace]] = {}

    def _load(self, organism: str) -> tuple[dict[str, list[LineageSnp]], list[RequestTrace]]:
        if organism not in self._index:
            barcode = self.http.get(BARCODE_URL)
            genes = self.http.get(KEGG_GENE_LIST.format(organism=organism))
            self._index[organism] = annotate(parse_gene_spans(genes.body),
                                             parse_barcode(barcode.body))
            self._traces[organism] = [_trace(barcode), _trace(genes)]
        return self._index[organism], self._traces[organism]

    def markers(self, gene: str, organism: str = BARCODE_ORGANISM) -> ToolResult:
        query: dict[str, Any] = {"gene": gene, "organism": organism}
        gene = gene.strip()
        organism = organism.strip().lower()

        if organism != BARCODE_ORGANISM:
            return ToolResult.build(
                query, [], resolved={"matched_by": "none"},
                notes=[(f"The lineage barcode is a list of positions in M. tuberculosis H37Rv "
                        f"coordinates ({BARCODE_ORGANISM}), so it cannot be applied to "
                        f"'{organism}'. No lookup was performed -- this is not evidence that "
                        f"{gene} lacks a lineage marker. Comparing these positions against "
                        f"another genome's gene spans matches by coordinate collision alone.")])

        if not gene:
            return ToolResult.build(
                query, [], resolved={"matched_by": "none"},
                notes=["No gene identifier was supplied. No lookup was performed."])

        try:
            index, traces = self._load(organism)
        except FetchError as exc:
            return ToolResult.build(
                query, [], resolved={"matched_by": "none"},
                notes=[(f"Could not fetch the lineage barcode or gene coordinates: HTTP "
                        f"{exc.status}. No marker call was made -- this is not evidence "
                        f"that {gene} lacks a lineage marker.")])

        # The index carries both locus tags and symbols; the fallback scan makes
        # the lookup case-insensitive, since KEGG symbols are mixed case (phoP,
        # katG) and callers do not reliably match them.
        snps = index.get(gene, [])
        if not snps:
            for key, value in index.items():
                if key.lower() == gene.lower():
                    snps = value
                    break
        matched_by = ("locus_tag" if snps and gene.lower().startswith("rv")
                      else "symbol" if snps else "none")

        records = [
            Record(
                record_id=f"tbdb:{snp.position}",
                type="lineage_marker",
                name=f"{snp.lineage} ({snp.lineage_name})",
                url=BARCODE_SOURCE,
                source="tbdb",
                retrieved_at=traces[0].retrieved_at,
                cached=traces[0].cached,
                detail={"position": snp.position, "lineage": snp.lineage,
                        "lineage_name": snp.lineage_name, "allele": snp.allele,
                        "gene": gene},
            )
            for snp in sorted(snps, key=lambda s: s.position)
        ]

        notes = []
        if records:
            notes.append(
                f"{gene} contains {len(records)} lineage-defining position(s). This says the gene "
                f"can carry a lineage marker, not that a particular variant in it is one -- "
                f"compare your variant's H37Rv position against the `position` field.")
        else:
            notes.append(
                f"No lineage-defining position falls within {gene} in the tbdb barcode. 855 of "
                f"4,008 H37Rv genes carry one, so this is an informative negative.")
        return ToolResult.build(query, records, resolved={"matched_by": matched_by},
                                requests=traces, notes=notes)
