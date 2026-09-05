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
from kegg_string_mcp.retrieval.compare import (
    compare,
    exact_term_probe,
    pair_queries,
    queries_for_pairs,
)
from kegg_string_mcp.retrieval.corpus import Corpus
from kegg_string_mcp.retrieval.independence import (
    IndependenceReport,
    classify,
    gene_partner_map,
    post_release_fraction,
)
from kegg_string_mcp.retrieval.index import HybridIndex, KeywordIndex, VectorIndex
from kegg_string_mcp.string_db import StringClient

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

    keyword = KeywordIndex(corpus)
    vector = VectorIndex(corpus)
    arms = {"lexical": keyword, "dense": vector,
            "hybrid": HybridIndex(keyword, vector)}

    string = StringClient(PoliteClient(DiskCache()))
    print(f"classifying {len(corpus.genes)} genes against STRING ...")
    report = IndependenceReport(verdicts=classify(corpus.genes,
                                                 gene_partner_map(corpus.genes, string)))
    report.post_release_papers, report.total_papers = post_release_fraction(corpus)
    (args.out / f"independence_{args.tag}.json").write_text(
        json.dumps(report.to_dict(), indent=1), encoding="utf-8")
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

    silent = report.pairs_with_status("silent")
    run(queries_for_pairs(silent), "PAIRS STRING IS SILENT ON",
        args.out / f"comparison_{args.tag}_silent.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
