"""Build the stage 2 retrieval corpus for the genes that need it.

Stage 2 covers what stage 1 cannot, so by default it runs on the genes whose
structured annotation is thin -- no UniProt function, only inferred function, or
no KEGG pathway -- rather than on every gene given. The routing decision is
written alongside the corpus so a later reader can see which genes were included
and why.

    NCBI_EMAIL=you@example.org python scripts/build_corpus.py

Pass --all-genes to bypass routing. The retrieval-arm comparison needs it: a
head-to-head on gene pairs requires every gene in the corpus regardless of how
well annotated it is, and that is a measurement rather than a pipeline run.

The corpus itself is not committed -- a few hundred PubMed abstracts under
publisher copyright is bulk redistribution, and it is reproducible from here.
Run scripts/run_comparison.py afterwards to measure the arms against each other.
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

def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genes", nargs="+", default=None)
    parser.add_argument("--extended", action="store_true",
                        help="use the larger gene set (robustness check)")
    parser.add_argument("--limit", type=int, default=20, help="abstracts per gene")
    parser.add_argument("--out", type=Path, default=Path("data"))
    parser.add_argument("--tag", default="tb", help="suffix for the output files")
    parser.add_argument("--all-genes", action="store_true",
                        help="skip coverage routing (needed for the arm comparison)")
    parser.add_argument("--functional-only", action="store_true",
                        help="route only genes with a functional gap, ignoring missing "
                             "KEGG pathways, which say more about KEGG than about the gene")
    args = parser.parse_args()
    genes = args.genes or (EXTENDED_GENES if args.extended else DEFAULT_GENES)

    import json

    from kegg_string_mcp.cache import DiskCache
    from kegg_string_mcp.http import PoliteClient
    from kegg_string_mcp.kegg import KeggClient
    from kegg_string_mcp.retrieval.corpus import build, chunk
    from kegg_string_mcp.retrieval.coverage import assess, route, summarise
    from kegg_string_mcp.uniprot import UniProtClient

    args.out.mkdir(parents=True, exist_ok=True)
    if args.all_genes:
        print(f"  routing bypassed: all {len(genes)} genes")
    else:
        http = PoliteClient(DiskCache())
        coverages = assess(genes, KeggClient(http), UniProtClient(http))
        summary = summarise(coverages)
        (args.out / f"coverage_{args.tag}.json").write_text(
            json.dumps({"summary": summary,
                        "coverage": [c.to_dict() for c in coverages]}, indent=1),
            encoding="utf-8")
        print(f"  coverage: {summary['well_covered']} well covered, "
              f"{summary['thin']} thin ({summary['functional_gap']} with a functional gap)")
        for reason, n in sorted(summary["reason_counts"].items(), key=lambda kv: -kv[1]):
            print(f"    {reason:26} {n}")
        if summary["lookup_failed"]:
            print("    UNANSWERED, not routed (fix the identifier or retry):")
            for gene, sources in sorted(summary["lookup_failed"].items()):
                print(f"      {gene:10} no answer from {', '.join(sources)}")
        genes = route(coverages, functional_only=args.functional_only)
        if not genes:
            print("  nothing to retrieve: every gene is already well annotated.")
            return 0
        print(f"  routed to stage 2 ({len(genes)}): {', '.join(genes)}")

    corpus = chunk(build(genes, limit=args.limit))
    corpus_path = corpus.write(args.out / f"corpus_{args.tag}.json")
    chars = sum(len(p.text) for p in corpus.passages)
    print(f"  corpus: {len(corpus.passages)} passages, {chars:,} chars -> {corpus_path}")
    for note in corpus.notes:
        print(f"    note: {note[:110]}")

    print("\n  next: python scripts/run_comparison.py "
          f"{corpus_path} --tag {args.tag}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
