#!/usr/bin/env python
"""Measure the retrieval arms against each other, on all pairs and on the pairs
STRING cannot answer.

The all-pairs numbers have a circularity: relevance is scored by whether a
passage names the queried genes, which is close to what BM25 ranks on. The
restricted run removes it. STRING's verdict comes from neither retriever, so
"pairs STRING is silent on" is a query set chosen independently of both arms --
and it is also the set where literature retrieval is the only evidence rather
than a restatement of a structured source.

    python scripts/run_comparison.py data/corpus_tb41.json --tag tb41
"""

from __future__ import annotations

import argparse
import json
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "src"))

from kegg_string_mcp.cache import DiskCache
from kegg_string_mcp.http import PoliteClient
from kegg_string_mcp.identity import resolve
from kegg_string_mcp.kegg import KeggClient
from kegg_string_mcp.retrieval.compare import (
    compare,
    exact_term_probe,
    pair_queries,
    queries_for_pairs,
)
from kegg_string_mcp.retrieval.corpus import Corpus, annotate_genes_named
from kegg_string_mcp.retrieval.independence import (
    IndependenceReport,
    classify,
    network_edges,
    post_release_fraction,
)
from kegg_string_mcp.retrieval.index import HybridIndex, KeywordIndex, VectorIndex
from kegg_string_mcp.string_db import StringClient
from kegg_string_mcp.uniprot import UniProtClient

PROBE_TERMS = ["katG", "Rv1908c", "ahpC", "Rv2428", "whiB7", "Rv3197A"]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("corpus", type=pathlib.Path)
    ap.add_argument("--tag", default="tb")
    ap.add_argument("--k", type=int, default=10)
    ap.add_argument("--out", type=pathlib.Path, default=pathlib.Path("data"))
    args = ap.parse_args()

    corpus = Corpus.read(args.corpus)
    print(f"corpus: {len(corpus.passages)} passages, {len(corpus.genes)} genes")

    http = PoliteClient(DiskCache())
    string = StringClient(http)
    print(f"resolving {len(corpus.genes)} genes ...")
    identities = resolve(corpus.genes, string=string, kegg=KeggClient(http),
                         uniprot=UniProtClient(http))
    unresolved = identities.unresolved("string")
    if unresolved:
        print(f"  STRING did not resolve {len(unresolved)}: {', '.join(unresolved)}")

    # A corpus built before the alias map existed counts only the literal query
    # string, so both arms and the judge undercount. Backfill and re-annotate
    # before indexing rather than silently measuring the old behaviour.
    if not corpus.aliases:
        corpus.aliases = identities.alias_map()
        annotate_genes_named(corpus)
        extra = sum(len(v) - 1 for v in corpus.aliases.values())
        print(f"  corpus had no alias map: added {extra} synonym(s) and re-annotated")

    keyword = KeywordIndex(corpus)
    vector = VectorIndex(corpus)
    arms = {"lexical": keyword, "dense": vector,
            "hybrid": HybridIndex(keyword, vector)}

    print(f"classifying {len(corpus.genes)} genes against STRING ...")
    edges = network_edges(corpus.genes, string, identities=identities)
    report = IndependenceReport(verdicts=classify(corpus.genes, edges),
                                unresolved=unresolved)
    report.post_release_papers, report.total_papers = post_release_fraction(corpus)
    (args.out / f"independence_{args.tag}.json").write_text(
        json.dumps(report.to_dict() | {"identities": identities.to_dict()}, indent=1),
        encoding="utf-8")
    print("  " + ", ".join(f"{k}={v}" for k, v in sorted(report.by_status().items())))

    def run(queries, label, path):
        result = compare(arms, queries, k=args.k)
        result.exact_term = exact_term_probe(arms, PROBE_TERMS, k=args.k, corpus=corpus)
        result.write(path)
        summary = result.to_dict()
        print(f"\n{label}  ({len(queries)} queries)")
        for arm in summary["arms"]:
            print(f"  {arm:9} precision@{args.k} {summary['mean_on_target_precision'][arm]:.3f}"
                  f"   papers naming both {summary['mean_papers_naming_both'][arm]:.2f}")
        for pair, val in summary["mean_overlap"].items():
            print(f"  overlap {pair:20} {val:.3f}")
        return result

    run(pair_queries(corpus.genes), "ALL PAIRS",
        args.out / f"comparison_{args.tag}.json")

    if unresolved:
        # Last thing printed, after the tables, because it is the finding a reader
        # must act on: those genes contributed no STRING evidence at all.
        print(f"\nUNRESOLVED ({len(unresolved)}): {', '.join(unresolved)} -- "
              f"their pairs are reported as unresolved, not silent, and are excluded "
              f"from the residue rather than counted as novel.")

    silent = report.pairs_with_status("silent")
    run(queries_for_pairs(silent), "PAIRS STRING IS SILENT ON",
        args.out / f"comparison_{args.tag}_silent.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
