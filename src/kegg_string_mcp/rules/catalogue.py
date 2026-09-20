"""What the WHO catalogue says about a locus, and about an interval.

Two lookups, and they are not the same question.

**By locus.** Three states that must stay apart. A gene with a graded-associated
variant is an *anchor* -- known pharmacology. A gene present in the catalogue
with none is *assessed and negative*, which is the signature the canonical
compensatory loci carry: rpoC (1,401 catalogued, 0 associated), rpoA (507, 0),
ahpC (250, 0). They are in the catalogue because they recur in resistant
isolates, and the committee graded their variants as not conferring resistance.
A gene absent from the catalogue was never assessed, which is neither.

**By interval.** The catalogue attributes a promoter variant to the gene
downstream of it by convention, so its `c.-N` coordinates are relative to a gene
that the variant may not physically sit in. `inhA c.-770` is 1,673,432, which is
in the Rv1482c-fabG1 gap; `inhA c.-154` is 1,674,048, which is inside fabG1's
coding sequence. Name matching therefore cannot connect a feature to a
catalogued variant -- only coordinates can.

The granularity caveat runs through both. The catalogue grades *specific
mutations*; a feature that says "this locus carries a protein-changing variant"
is a weaker statement. katG holds 1,205 non-synonymous catalogued variants of
which 133 are graded associated and 960 are "Uncertain significance". So a locus
being an anchor means the GENE is known to resistance, never that the variant
driving the rule is one of the graded ones.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Any

from kegg_string_mcp.resistance import Variant
from kegg_string_mcp.rules.annotation import Annotation, Intergenic, normalise_locus

ANCHOR = "anchor"                       # >=1 graded-associated variant
ASSESSED_NEGATIVE = "assessed_negative"  # in the catalogue, none associated
ABSENT = "absent"                        # never assessed

_UPSTREAM = re.compile(r"^c\.-(\d+)")


@dataclass(frozen=True)
class Placed:
    """A catalogue variant put on the genome, with the gene it was named for."""

    position: int
    named_for: str
    variant: Variant


@dataclass
class CatalogueStatus:
    """What the catalogue holds for one locus or interval."""

    status: str                              # ANCHOR | ASSESSED_NEGATIVE | ABSENT
    drugs: list[str] = field(default_factory=list)          # associated FOR these
    # The drugs the gene was assessed AGAINST, associated or not. A compensator
    # is only plausible for the drug its locus was catalogued under: rpoA and
    # rpoC are assessed for rifampicin only, ahpC for isoniazid only. Without
    # this, any assessed-negative locus pairs with any anchor.
    assessed_drugs: list[str] = field(default_factory=list)
    associated: int = 0
    catalogued: int = 0
    placed: list[Placed] = field(default_factory=list)   # intervals only
    note: str = ""

    @property
    def is_anchor(self) -> bool:
        return self.status == ANCHOR

    def to_dict(self) -> dict[str, Any]:
        return {"catalogue_status": self.status, "drugs": "|".join(self.drugs) or "NA",
                "assessed_drugs": "|".join(self.assessed_drugs) or "NA",
                "associated_variants": self.associated,
                "catalogued_variants": self.catalogued,
                "note": self.note}


class Catalogue:
    """Locus and interval lookups over one parsed catalogue and one annotation.

    Built once per run. The placement pass is over the whole catalogue, so doing
    it per interval would re-scan 49,330 rows for every feature.
    """

    def __init__(self, catalogue: dict[str, list[Variant]], annotation: Annotation):
        self.annotation = annotation
        self._by_gene: dict[str, list[Variant]] = {}
        for gene, variants in catalogue.items():
            self._by_gene[gene] = list(variants)
            # The catalogue names genes by symbol or by locus tag depending on the
            # gene, so both spellings have to reach the same entry.
            self._by_gene.setdefault(normalise_locus(gene), list(variants))
        self._placed: list[Placed] = sorted(
            self._place(catalogue), key=lambda p: p.position)

    def _place(self, catalogue: dict[str, list[Variant]]) -> list[Placed]:
        """Upstream variants onto genome coordinates, honouring strand.

        Only `c.-N` forms are placed: they are the ones that can fall outside the
        gene they are named for. A coding or protein-level variant is inside its
        own gene by definition and is found by the locus lookup instead.
        """
        out: list[Placed] = []
        for gene, variants in catalogue.items():
            span = self._span(gene)
            if span is None:
                continue
            for variant in variants:
                match = _UPSTREAM.match(variant.mutation)
                if match is None:
                    continue
                offset = int(match.group(1))
                # Upstream is lower coordinate on the forward strand and higher on
                # the reverse. Getting this backwards puts every promoter variant
                # inside the gene instead of in front of it.
                #
                # Counted from the TRANSLATION start: HGVS `c.` numbering is
                # relative to the ATG, and three H37Rv gene rows begin before
                # their CDS -- Rv0614 by 243 bp -- so the gene's 5' end would put
                # the variant that far from where it is.
                origin = span.coding_start
                position = origin - offset if span.strand != "-" else origin + offset
                out.append(Placed(position=position, named_for=gene, variant=variant))
        return out

    def _span(self, gene: str):
        found = self.annotation.gene(gene)
        if found is not None:
            return found
        by_norm = self.annotation.by_normalised(gene)
        if len(by_norm) == 1:
            return by_norm[0]
        by_symbol = self.annotation.by_symbol(gene)
        return by_symbol[0] if len(by_symbol) == 1 else None

    # -- lookups -------------------------------------------------------------

    def for_locus(self, locus: str, symbol: str = "") -> CatalogueStatus:
        variants = (self._by_gene.get(locus) or self._by_gene.get(normalise_locus(locus))
                    or (self._by_gene.get(symbol) if symbol else None))
        if not variants:
            return CatalogueStatus(
                ABSENT, note=(f"'{locus}' is not in the WHO catalogue, which covers 74 genes "
                              f"selected for resistance surveillance. Absence means the gene "
                              f"was not assessed, not that it is unrelated to resistance."))
        assessed = sorted({v.drug for v in variants if v.drug})
        associated = [v for v in variants if v.associated]
        if associated:
            return CatalogueStatus(
                ANCHOR, drugs=sorted({v.drug for v in associated if v.drug}),
                assessed_drugs=assessed,
                associated=len(associated), catalogued=len(variants),
                note=(f"{len(associated)} of {len(variants)} catalogued variants are graded "
                      f"associated. The flag is about the GENE: a rule condition says the "
                      f"locus carries some qualifying variant, not that it carries one of "
                      f"these."))
        return CatalogueStatus(
            ASSESSED_NEGATIVE, catalogued=len(variants), assessed_drugs=assessed,
            note=(f"assessed against {', '.join(assessed) or 'no drug'}: {len(variants)} "
                  f"catalogued variants, none graded associated. "
                  f"This is the signature the known compensatory loci carry -- present "
                  f"because they recur in resistant isolates, graded as not conferring "
                  f"resistance themselves."))

    def for_interval(self, interval: Intergenic) -> CatalogueStatus:
        inside = [p for p in self._placed if interval.start <= p.position <= interval.end]
        if not inside:
            return CatalogueStatus(
                ABSENT, note=(f"no catalogued variant falls within {interval.name} "
                              f"({interval.width} bp). The catalogue covers 74 genes, so this "
                              f"is not evidence that variation here is unimportant."))
        associated = [p for p in inside if p.variant.associated]
        named = sorted({p.named_for for p in inside})
        if associated:
            drugs = sorted({p.variant.drug for p in associated if p.variant.drug})
            return CatalogueStatus(
                ANCHOR, drugs=drugs, associated=len(associated), catalogued=len(inside),
                placed=inside,
                note=(f"{len(associated)} graded-associated catalogued variant(s) fall inside "
                      f"{interval.name}, attributed by the catalogue to {', '.join(named)}. "
                      f"The catalogue names a promoter variant for the gene downstream of it, "
                      f"so the name and the interval need not agree."))
        return CatalogueStatus(
            ASSESSED_NEGATIVE, catalogued=len(inside), placed=inside,
            note=(f"{len(inside)} catalogued variant(s) fall inside {interval.name}, none "
                  f"graded associated."))
