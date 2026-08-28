"""The gold-standard set, and what it can and cannot measure.

Frozen from KEGG's own `/link/pathway/mtu` on a recorded date, so the reference
is reproducible and its provenance is explicit -- the same discipline the pipeline
applies to its own output.

The reference is deliberately KEGG itself. That makes recall a measure of
**retrieval fidelity** (did the pipeline faithfully report what its tools
returned?) and not of biological truth. Conflating the two would be the mistake:
KEGG assigns no pathway to 71% of M. tuberculosis genes, gyrA among them, so
scoring biological correctness against it would mark a correct "KEGG has nothing"
as a miss.

Hence two classes of gene:

* **positive controls** -- KEGG assigns pathways. Measures recall and precision.
* **negative controls** -- KEGG assigns none, so the correct answer is "nothing
  found". Measures whether the pipeline abstains or fabricates, which for this
  project is the more important number.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

GOLD_PATH = Path(__file__).with_name("gold_mtu.json")


@dataclass(frozen=True)
class GoldGene:
    kegg_gene_id: str
    symbol: str
    pathways: tuple[str, ...]

    @property
    def query(self) -> str:
        """What to ask the pipeline. The symbol where one exists -- that is how a
        user would ask, and it exercises symbol resolution rather than bypassing it."""
        return self.symbol or self.kegg_gene_id

    @property
    def is_negative_control(self) -> bool:
        return not self.pathways


@dataclass(frozen=True)
class GoldSet:
    organism: str
    reference: str
    retrieved_on: str
    coverage: dict
    genes: tuple[GoldGene, ...]

    @property
    def positives(self) -> tuple[GoldGene, ...]:
        return tuple(g for g in self.genes if not g.is_negative_control)

    @property
    def negatives(self) -> tuple[GoldGene, ...]:
        return tuple(g for g in self.genes if g.is_negative_control)


def load(path: Path | None = None) -> GoldSet:
    payload = json.loads((path or GOLD_PATH).read_text(encoding="utf-8"))
    return GoldSet(
        organism=payload["organism"],
        reference=payload["reference"],
        retrieved_on=payload["retrieved_on"],
        coverage=payload["coverage"],
        genes=tuple(GoldGene(g["kegg_gene_id"], g["symbol"], tuple(g["pathways"]))
                    for g in payload["genes"]),
    )
