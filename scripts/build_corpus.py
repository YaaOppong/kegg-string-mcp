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

# The set the KEGG gold standard uses: canonical TB resistance and oxidative-stress
# loci, half annotated in KEGG and half not.
DEFAULT_GENES = ["katG", "inhA", "rpoB", "embB", "pncA", "fabG1",
                 "gyrA", "gyrB", "ahpC", "furA", "pknB", "sodA"]

# A larger set spanning resistance, virulence, regulation, metabolism and cell
# wall. Not a different question -- a robustness check on the same one. A
# 198-passage corpus is small enough that the arms' ranking could be an artefact
# of it; if the conclusions hold across roughly three times the text and six times
# the pair queries, they are about the retrieval methods rather than this corpus.
EXTENDED_GENES = DEFAULT_GENES + [
    # further resistance loci
    "rpsL", "eis", "ethA", "gidB", "tlyA", "kasA", "ndh", "embA", "embC", "rpoC",
    # virulence and secretion
    "esxA", "esxB", "phoP", "phoR", "whiB3", "whiB7", "sigH", "sigE",
    # dormancy and regulation
    "dosR", "dosS", "mprA", "relA", "zur",
    # metabolism and cell wall
    "icl1", "glpK", "pckA", "pks13", "mmpL3", "embR",
]

EXACT_TERMS = ["katG", "Rv1908c", "ahpC", "gyrA", "pncA", "embB"]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genes", nargs="+", default=None)
    parser.add_argument("--extended", action="store_true",
                        help="use the larger gene set (robustness check)")
    parser.add_argument("--limit", type=int, default=20, help="abstracts per gene")
    parser.add_argument("--k", type=int, default=10)
    parser.add_argument("--out", type=Path, default=Path("data"))
    parser.add_argument("--tag", default="tb", help="suffix for the output files")
    args = parser.parse_args()
    genes = args.genes or (EXTENDED_GENES if args.extended else DEFAULT_GENES)

    from kegg_string_mcp.retrieval.compare import compare, exact_term_probe, pair_queries
    from kegg_string_mcp.retrieval.corpus import build, chunk
    from kegg_string_mcp.retrieval.index import HybridIndex, KeywordIndex, VectorIndex

    corpus = chunk(build(genes, limit=args.limit))
    corpus_path = corpus.write(args.out / f"corpus_{args.tag}.json")
    chars = sum(len(p.text) for p in corpus.passages)
    print(f"  corpus: {len(corpus.passages)} passages, {chars:,} chars -> {corpus_path}")
    for note in corpus.notes:
        print(f"    note: {note[:110]}")

    keyword, vector = KeywordIndex(corpus), VectorIndex(corpus)
    arms = {"keyword": keyword, "vector": vector,
            "hybrid": HybridIndex(keyword, vector)}

    queries = pair_queries(corpus.genes)
    print(f"  {len(queries)} pair queries across {len(corpus.genes)} genes")
    result = compare(arms, queries, k=args.k)
    result.exact_term = exact_term_probe(arms, EXACT_TERMS, k=args.k, corpus=corpus)
    path = result.write(args.out / f"comparison_{args.tag}.json")

    data = result.to_dict()
    print(f"\n  {data['queries']} gene-pair queries, k={data['k']}")
    print("  mean on-target precision@k :", data["mean_on_target_precision"])
    print("  mean papers naming both    :", data["mean_papers_naming_both"])
    print("  mean overlap               :", data["mean_overlap"])
    print("\n  exact-term probe (top-k passages containing the queried term):")
    for row in data["exact_term"]:
        cells = "  ".join(f"{a}={row[a]['containing_the_term']}/{args.k}" for a in arms)
        print(f"    {row['term']:10} in_corpus={row['in_corpus']:<4} {cells}")
    print(f"\n  written: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
