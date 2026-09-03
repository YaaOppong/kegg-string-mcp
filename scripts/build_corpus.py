"""Rebuild the retrieval corpus and run the head-to-head comparison.

The corpus itself is not committed -- a few hundred PubMed abstracts under
publisher copyright is bulk redistribution, and it is reproducible from here in
one command:

    NCBI_EMAIL=you@example.org python scripts/build_corpus.py

Both arms then search exactly the same text, retrieved the same way, so any
difference between them is the retrieval method rather than a difference in what
was fetched.
"""

from __future__ import annotations

import argparse
from pathlib import Path

DEFAULT_GENES = ["katG", "inhA", "rpoB", "embB", "pncA", "fabG1",
                 "gyrA", "gyrB", "ahpC", "furA", "pknB", "sodA"]

EXACT_TERMS = ["katG", "Rv1908c", "ahpC", "gyrA", "pncA", "embB"]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genes", nargs="+", default=DEFAULT_GENES)
    parser.add_argument("--limit", type=int, default=20, help="abstracts per gene")
    parser.add_argument("--k", type=int, default=10)
    parser.add_argument("--out", type=Path, default=Path("data"))
    args = parser.parse_args()

    from kegg_string_mcp.retrieval.compare import compare, exact_term_probe, pair_queries
    from kegg_string_mcp.retrieval.corpus import build
    from kegg_string_mcp.retrieval.index import HybridIndex, KeywordIndex, VectorIndex

    corpus = build(args.genes, limit=args.limit)
    corpus_path = corpus.write(args.out / "corpus_tb.json")
    chars = sum(len(p.text) for p in corpus.passages)
    print(f"  corpus: {len(corpus.passages)} passages, {chars:,} chars -> {corpus_path}")
    for note in corpus.notes:
        print(f"    note: {note[:110]}")

    keyword, vector = KeywordIndex(corpus), VectorIndex(corpus)
    arms = {"keyword": keyword, "vector": vector,
            "hybrid": HybridIndex(keyword, vector)}

    result = compare(arms, pair_queries(corpus.genes), k=args.k)
    result.exact_term = exact_term_probe(arms, EXACT_TERMS, k=args.k)
    path = result.write(args.out / "comparison.json")

    data = result.to_dict()
    print(f"\n  {data['queries']} gene-pair queries, k={data['k']}")
    print("  mean on-target precision@k :", data["mean_on_target_precision"])
    print("  mean papers naming both    :", data["mean_papers_naming_both"])
    print("  mean overlap               :", data["mean_overlap"])
    print("\n  exact-term probe (naming the queried term / k):")
    for row in data["exact_term"]:
        cells = "  ".join(f"{a}={row[a]['naming_the_term']}/{args.k}" for a in arms)
        print(f"    {row['term']:10} {cells}")
    print(f"\n  written: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
