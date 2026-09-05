#!/usr/bin/env python
"""Compute the unexplained residue: what stages 1 and 2 could not account for.

Reads the artefacts the earlier steps wrote, adds KEGG pathway membership, and
reports how many pairs survive -- the input to hypothesis generation.

    python scripts/residue.py --tag tb41
"""

from __future__ import annotations

import argparse
import json
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "src"))

from kegg_string_mcp.cache import DiskCache
from kegg_string_mcp.http import PoliteClient
from kegg_string_mcp.hypothesis.residue import assess, residue, summarise
from kegg_string_mcp.kegg import KeggClient
from kegg_string_mcp.lineage import LineageClient
from kegg_string_mcp.retrieval.corpus import Corpus


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--tag", default="tb41")
    ap.add_argument("--data", type=pathlib.Path, default=pathlib.Path("data"))
    args = ap.parse_args()

    independence = json.loads((args.data / f"independence_{args.tag}.json").read_text())
    comparison = json.loads((args.data / f"comparison_{args.tag}.json").read_text())
    corpus = Corpus.read(args.data / f"corpus_{args.tag}.json")

    string_status = {(v["gene_a"], v["gene_b"]): v for v in independence["verdicts"]}
    # Best co-mention count any arm found. Taking the max rather than one arm's
    # figure keeps the gate from calling a pair novel because a single retriever
    # missed the paper connecting it.
    co_mentions = {tuple(r["genes"]): max(r["joint"].values())
                   for r in comparison["per_query"]}

    http = PoliteClient(DiskCache())
    kegg = KeggClient(http)
    pathways: dict[str, set[str]] = {}
    for gene in corpus.genes:
        pathways[gene] = {r.record_id for r in kegg.pathways(gene).records}

    # Lineage markers do not gate the residue. They are recorded against each
    # surviving pair so a known marker is visible before any mechanism is
    # proposed for it.
    lineage = LineageClient(http)

    pairs = [(v["gene_a"], v["gene_b"]) for v in independence["verdicts"]]
    assessments = assess(pairs, string_status=string_status,
                         pathways=pathways, co_mentions=co_mentions)
    summary = summarise(assessments)

    remaining = residue(assessments)
    marked = {g: sorted({r.detail['lineage'] for r in lineage.markers(g).records})
              for g in corpus.genes}
    out = args.data / f"residue_{args.tag}.json"
    out.write_text(json.dumps(
        {"summary": summary,
         "residue": [dict(a.to_dict(),
                          lineage_marker={a.gene_a: marked[a.gene_a],
                                          a.gene_b: marked[a.gene_b]})
                     for a in remaining],
         "assessments": [a.to_dict() for a in assessments]}, indent=1), encoding="utf-8")

    print(f"{summary['pairs']} pairs assessed")
    for code, n in sorted(summary["reason_counts"].items(), key=lambda kv: -kv[1]):
        counts = " (counts as explained)" if code in summary["explaining"] else ""
        print(f"  {code:22} {n:4}{counts}")
    print(f"\nresidue: {summary['residue']} pairs "
          f"({summary['residue_fraction']:.0%}) -> {out}")
    both = sum(1 for a in remaining if marked[a.gene_a] and marked[a.gene_b])
    either = sum(1 for a in remaining if marked[a.gene_a] or marked[a.gene_b])
    flagged = sorted(g for g, lins in marked.items() if lins)
    print(f"\nlineage markers: {len(flagged)}/{len(marked)} genes -- "
          f"{', '.join(flagged)}")
    print(f"  residue pairs with one marked gene : {either - both}")
    print(f"  residue pairs with both marked     : {both}")

    genes_with_no_pathway = [g for g, p in pathways.items() if not p]
    if genes_with_no_pathway:
        print(f"\nno KEGG pathway ({len(genes_with_no_pathway)}): "
              f"{', '.join(sorted(genes_with_no_pathway))}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
