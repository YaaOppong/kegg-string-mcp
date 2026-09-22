"""Capture what the rule classification needs, so a browser can recompute it.

    python demo/build_rule_fixture.py RULES.tsv --annotation h37rv.gff

The existing demo replays a model's output and re-runs the validator over it.
This one has nothing to replay: the classification is arithmetic over structured
data, so given the data the page can compute every verdict itself. What is
captured is therefore the *inputs*, never a result.

Three trims keep it small enough to inline in a page:

* the annotation is cut to the loci the rules name and their neighbours, since
  an intergenic feature needs the genes on both sides;
* the catalogue is cut to those loci's genes, and each row to the four fields
  the classification reads. `gene` repeats the key, `source` repeats one string
  5,560 times, and `comment` is free text nothing reads -- together three
  quarters of the payload. `mutation` stays, because a promoter variant's `c.-N`
  offset is what places it on the genome;
* the lineage barcode is cut to positions inside the loci kept.

The interaction data is captured rather than recomputed, because it is the one
input that cannot be derived from a file: STRING partner lists per locus, and the
links between co-occurring loci with their channel. Without them a compensation
verdict reads "no link was found, so this is co-occurrence" -- a weaker and
different claim from the one the pipeline makes, and the page would be
misleading rather than merely incomplete.

Everything else is the real thing. The page parses the real GFF text with the
real parser and classifies with the real classifier, so a verdict it shows is one
the library produces.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "demo" / "rule_fixture.json"

# Neighbouring loci to keep on each side, so an intergenic feature and the
# aliasing checks that depend on adjacency still resolve.
FLANK = 2


def _gff_lines(path: Path) -> list[str]:
    return [line.rstrip("\n") for line in path.read_text(encoding="utf-8").splitlines()
            if line.strip() and not line.startswith("#")]


def build(rules_path: Path, annotation_path: Path, out: Path = OUT) -> Path:
    import sys

    sys.path.insert(0, str(ROOT / "src"))
    from kegg_string_mcp.cache import DiskCache
    from kegg_string_mcp.http import PoliteClient
    from kegg_string_mcp.lineage import BARCODE_URL, parse_barcode
    from kegg_string_mcp.resistance import ASSOCIATED, ResistanceClient
    from kegg_string_mcp.rules.annotation import parse as parse_annotation
    from kegg_string_mcp.rules.parse import parse as parse_rules
    from kegg_string_mcp.rules.parse import vocabulary

    annotation = parse_annotation(annotation_path)
    rules = parse_rules(rules_path)

    # Which loci the rules reach, plus their neighbours.
    wanted: set[str] = set()
    order = {gene.locus: index for index, gene in enumerate(annotation.genes)}
    for label in vocabulary(rules):
        for part in label.replace("upstream_", "").replace("downstream_", "") \
                         .replace("intergenic_", "").split("-"):
            bare = part.split(".")[0]
            gene = annotation.gene(bare) or (annotation.by_normalised(bare) or [None])[0] \
                or (annotation.by_symbol(bare) or [None])[0]
            if gene is None:
                continue
            index = order[gene.locus]
            for offset in range(-FLANK, FLANK + 1):
                if 0 <= index + offset < len(annotation.genes):
                    wanted.add(annotation.genes[index + offset].locus)

    keep = {gene.locus for gene in annotation.genes if gene.locus in wanted}
    lines = [line for line in _gff_lines(annotation_path)
             if any(f"Locus={locus};" in line or f"locus_tag={locus};" in line
                    for locus in keep)]

    http = PoliteClient(DiskCache())
    catalogue, _ = ResistanceClient(http).catalogue()
    symbols = {annotation.gene(locus).symbol for locus in keep if annotation.gene(locus)}
    rows: dict[str, list[dict[str, str]]] = {}
    for gene, variants in catalogue.items():
        if gene in keep or gene in symbols:
            # The comment is free text and is never read by the classification.
            rows[gene] = [{"mutation": v.mutation, "drug": v.drug,
                           "confidence": v.confidence} for v in variants]

    spans = [(annotation.gene(locus).start, annotation.gene(locus).end)
             for locus in keep if annotation.gene(locus)]
    barcode = [{"position": s.position, "lineage": s.lineage, "lineage_name": s.lineage_name,
                "allele": s.allele}
               for s in parse_barcode(http.get(BARCODE_URL).body)
               if any(lo <= s.position <= hi for lo, hi in spans)]

    # Interaction evidence, from the same path `gar rules` uses.
    partners, links = _interaction(rules, annotation, keep, http)

    fixture = {
        "annotation_sha256": annotation.sha256,
        "partners": partners,
        "links": links,
        "annotation_gff": "\n".join(lines) + "\n",
        "rules_tsv": rules_path.read_text(encoding="utf-8"),
        "catalogue": rows,
        "associated_grades": list(ASSOCIATED),
        "barcode": barcode,
        "loci_kept": sorted(keep),
        # Enrichment is a statement about a background. The page holds twenty
        # loci; without the real counts it would report a different q under the
        # same name, so the counts travel and the membership map covers only
        # what is tested.
        "universe": _universe_counts(annotation),
    }
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(fixture, separators=(",", ":")) + "\n", encoding="utf-8")
    return out


def _universe_counts(annotation) -> dict[str, Any]:
    from kegg_string_mcp.rules.enrichment import universe_for

    universe = universe_for(annotation)
    return {"axis": universe.axis, "total": universe.size, "sizes": dict(universe.sizes)}


def _interaction(rules, annotation, keep, http):
    """Partner lists per locus, and the links between co-occurring loci."""
    from kegg_string_mcp.kegg import KeggClient
    from kegg_string_mcp.rules.annotate import (
        LocusAnnotation,
        annotate_locus,
        co_occurring,
        compute_links,
    )
    from kegg_string_mcp.rules.catalogue import Catalogue
    from kegg_string_mcp.rules.evidence import Sources, rename_map
    from kegg_string_mcp.string_db import StringClient
    from kegg_string_mcp.uniprot import UniProtClient

    class _Clients:
        pass

    clients = _Clients()
    clients.string, clients.kegg = StringClient(http), KeggClient(http)
    clients.uniprot = UniProtClient(http)
    sources = Sources(annotation=annotation, catalogue=Catalogue({}, annotation))
    rename = rename_map(rules, sources)

    annotations: dict[str, LocusAnnotation] = {}
    for locus in sorted(set(rename.values()) & keep):
        gene = annotation.gene(locus)
        annotations[locus] = annotate_locus(locus, gene.symbol if gene else "",
                                            sources, clients)
    try:
        sizes, _ = clients.kegg.pathway_sizes("mtu")
        index, _ = clients.kegg.gene_index("mtu")
        genome = len(index.locus_tags)
    except Exception:                                    # noqa: BLE001
        sizes, genome = {}, 0

    edges: dict[tuple[str, str], dict] = {}
    by_string_id = {a.string_id: locus for locus, a in annotations.items() if a.string_id}
    if len(by_string_id) >= 2:
        try:
            network = clients.string.network(sorted(by_string_id),
                                             required_score=150).model_dump()
            for record in network.get("records", []):
                detail = record.get("detail", {})
                left = by_string_id.get(str(detail.get("string_id_a", "")))
                right = by_string_id.get(str(detail.get("string_id_b", "")))
                if left and right:
                    edges[(left, right)] = detail
        except Exception as exc:                         # noqa: BLE001
            # A fixture built without edges would silently produce weaker
            # compensation verdicts than the pipeline gives, so say it happened.
            print(f"warning: STRING network unavailable ({exc}); the fixture will "
                  f"carry no interaction links")

    found = compute_links(co_occurring(rules, rename), annotations, sources,
                          pathway_sizes=sizes, genome_size=genome, edges=edges)
    return (
        {locus: sorted({name or pid for pid, name, _ in a.partners})
         for locus, a in annotations.items()},
        [{"a": a, "b": b, "kind": link.kind, "detail": link.detail}
         for (a, b), links in found.items() for link in links],
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("rules", type=Path)
    parser.add_argument("--annotation", type=Path, required=True)
    parser.add_argument("--out", type=Path, default=OUT)
    args = parser.parse_args()
    path = build(args.rules, args.annotation, args.out)
    size = path.stat().st_size
    print(f"{path}  {size / 1024:.0f} KB")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
