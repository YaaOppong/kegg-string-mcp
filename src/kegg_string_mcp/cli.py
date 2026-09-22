"""Command line entry point for the annotation pipeline.

    gar single katG
    gar epistasis katG furA ahpC
    gar eval
    gar table                       # tabulate every finished run under runs/
    gar features LABELS --annotation H37RV.bed    # do the feature labels resolve?
    gar rules RULES.tsv --annotation H37RV.bed    # annotate loci, classify rules
    gar gold --annotation H37RV.bed               # score the classifier on known rules
    gar retrieve QUESTIONS.tsv --annotation H37RV.bed   # fetch candidate papers
"""

from __future__ import annotations

import argparse
import asyncio
import json
import os
import sys
from contextlib import asynccontextmanager
from pathlib import Path

from kegg_string_mcp.agent import annotate_epistasis, annotate_gene
from kegg_string_mcp.agent.mcp_tools import mcp_tools


@asynccontextmanager
async def _tools(args):
    """A live MCP session by default; the in-process dispatch under --direct.

    Going over MCP is the point: the pipeline is then a client of its own server,
    and the tool schemas come from the one place they are defined.
    """
    if args.direct:
        yield None
    else:
        async with mcp_tools() as tools:
            yield tools


async def _run(annotate, args, genes):
    async with _tools(args) as tools:
        return await annotate(genes, args.organism, args.runs, tools=tools)


async def _run_eval(evaluate, args):
    async with _tools(args) as tools:
        return await evaluate(runs=args.runs / "eval", tools=tools)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="gar", description=__doc__)
    parser.add_argument("mode", choices=["single", "epistasis", "eval", "table", "features", "rules", "gold", "retrieve"])
    parser.add_argument("genes", nargs="*")
    parser.add_argument("--organism", default="mtu")
    parser.add_argument("--runs", type=Path, default=Path("runs"))
    parser.add_argument("--json", action="store_true", help="emit the full payload")
    # A DIRECTORY for every mode that writes. It previously meant a file for
    # `features` and a directory for `rules`, so running one after the other with
    # the same --out made the second fail on a path that was already a file.
    parser.add_argument("--out", type=Path, default=Path("tables"),
                        help="directory the TSVs are written to")
    # The reference annotation is an input, never fetched: the features in a rule
    # set are named by whichever annotation the variant caller used, and resolving
    # them against a different gene list is how a locus becomes the wrong locus.
    parser.add_argument("--annotation", type=Path,
                        default=(Path(os.environ["KEGG_STRING_MCP_ANNOTATION"])
                                 if os.environ.get("KEGG_STRING_MCP_ANNOTATION") else None),
                        help="gene annotation (BED/GFF/GTF) the feature labels were named "
                             "from; defaults to $KEGG_STRING_MCP_ANNOTATION")
    parser.add_argument("--limit", type=int, default=10,
                        help="records to retrieve per question (`retrieve`)")
    parser.add_argument("--strict", action="store_true",
                        help="`features` exits non-zero if any label fails to resolve")
    parser.add_argument("--direct", action="store_true",
                        help="dispatch tools in-process instead of over MCP (debugging)")
    args = parser.parse_args(argv)

    if args.mode == "table":
        # A pure post-processor over stores already on disk: no model, no network,
        # and re-runnable. That is what lets the fan-out over many genes be someone
        # else's job -- a Snakemake rule per gene, then one aggregate rule here.
        from kegg_string_mcp.agent.tables import build

        written = build(args.runs, args.out)
        for label, path in written.items():
            rows = max(0, sum(1 for _ in path.open(encoding="utf-8")) - 1)
            print(f"{label:9} {rows:5} rows  {path}")
        return 0

    if args.mode == "retrieve":
        # Stage B of the literature layer: retrieval only, no model. Queries are
        # gene-centric so that finding a passage is not made near-certain by the
        # query having named the conclusion.
        from kegg_string_mcp.agent.loop import new_store
        from kegg_string_mcp.cache import DiskCache
        from kegg_string_mcp.http import PoliteClient
        from kegg_string_mcp.rules.annotation import parse as parse_annotation
        from kegg_string_mcp.rules.questions import load as load_questions
        from kegg_string_mcp.rules.report import write_retrieval
        from kegg_string_mcp.rules.retrieve import Clients, retrieve

        if len(args.genes) != 1:
            parser.error("retrieve mode takes exactly one questions file")
        if args.annotation is None:
            parser.error("retrieve mode needs --annotation (or $KEGG_STRING_MCP_ANNOTATION)")

        questions = load_questions(Path(args.genes[0]))
        annotation = parse_annotation(args.annotation)
        clients = Clients.live(PoliteClient(DiskCache()))
        if clients.corpus is None:
            print("note: no corpus configured, so every question goes to PubMed. Set "
                  "KEGG_STRING_MCP_CORPUS to a corpus from scripts/build_corpus.py.")
        store = new_store(args.out / "retrieval", "retrieve")
        report = retrieve(questions, annotation, clients, store=store, limit=args.limit)

        for note in report.notes:
            print(f"note: {note}")
        summary = report.summary()
        print(f"\n{summary['questions']:,} question(s)")
        print(f"  with at least one candidate   {summary['with_a_candidate']:,}")
        print(f"  with none                     {summary['without']:,}")
        print(f"  records retrieved             {summary['records']:,} "
              f"({summary['distinct_records']:,} distinct)")
        for source, count in sorted(summary["by_source"].items()):
            print(f"    from {source:22} {count:,}")
        print(f"\nretrieval: {write_retrieval(args.out / 'retrieval.tsv', report.candidates)}")
        print(f"store:     {store.path}")
        return 0

    if args.mode == "gold":
        # Scores the classifier against rules whose answer is known. Needs the
        # catalogue and the annotation, but no model and no STRING: the expected
        # signatures turn on the catalogue, not on link evidence.
        from kegg_string_mcp.cache import DiskCache
        from kegg_string_mcp.http import PoliteClient
        from kegg_string_mcp.lineage import LineageClient
        from kegg_string_mcp.resistance import ResistanceClient
        from kegg_string_mcp.rules.annotation import parse as parse_annotation
        from kegg_string_mcp.rules.evidence import load_sources
        from kegg_string_mcp.rules.gold import (
            load,
            load_sets,
            render,
            render_sets,
            score,
            score_sets,
            summary,
        )

        if args.annotation is None:
            parser.error("gold mode needs --annotation (or $KEGG_STRING_MCP_ANNOTATION)")
        http = PoliteClient(DiskCache())
        sources = load_sources(parse_annotation(args.annotation), ResistanceClient(http),
                               LineageClient(http), args.organism)
        gold = load()
        annotation = parse_annotation(args.annotation)
        scores = score(sources, gold)
        print(render(scores, gold))

        # The enrichment gold set: properties of the supplied annotation rather
        # than of this code, so it is scored and reported separately.
        set_scores = score_sets(annotation)
        print("\n" + "-" * 70 + "\nenrichment\n")
        print(render_sets(set_scores, load_sets()[1]))
        # Non-zero on a real failure, but not on a skip: a locus the supplied
        # annotation lacks is a fact about the annotation, not a wrong answer.
        failed = summary(scores)["failed"] + sum(1 for s in set_scores if s.failures)
        return 1 if failed else 0

    if args.mode == "rules":
        # The deterministic half: no model, no API key. Clients are the direct
        # ones rather than the MCP dispatch, because nothing here is a tool call
        # the model makes -- it is the pipeline fetching for itself.
        from kegg_string_mcp.cache import DiskCache
        from kegg_string_mcp.http import PoliteClient
        from kegg_string_mcp.kegg import KeggClient
        from kegg_string_mcp.lineage import LineageClient
        from kegg_string_mcp.resistance import ResistanceClient
        from kegg_string_mcp.rules.run import run
        from kegg_string_mcp.string_db import StringClient
        from kegg_string_mcp.uniprot import UniProtClient

        if len(args.genes) != 1:
            parser.error("rules mode takes exactly one rule file")
        if args.annotation is None:
            parser.error("rules mode needs --annotation (or $KEGG_STRING_MCP_ANNOTATION)")

        http = PoliteClient(DiskCache())
        kegg = KeggClient(http)

        class _Clients:
            pass

        clients = _Clients()
        clients.kegg, clients.string = kegg, StringClient(http)
        clients.uniprot = UniProtClient(http)

        # Pathway sizes decide whether a shared pathway is a link or a base rate.
        # A failure here must not take the run down: without them, shared
        # pathways are simply not reported as links.
        try:
            sizes, _ = kegg.pathway_sizes(args.organism)
            index, _ = kegg.gene_index(args.organism)
            genome_size = len(index.locus_tags)
        except Exception as exc:                      # noqa: BLE001
            sizes, genome_size = {}, 0
            print(f"note: pathway sizes unavailable ({exc}); shared pathways not linked")

        result = run(Path(args.genes[0]), args.annotation, args.out,
                     clients=clients, resistance=ResistanceClient(http),
                     lineage=LineageClient(http), organism=args.organism,
                     genome_size=genome_size, pathway_sizes=sizes)

        for note in result.notes:
            print(f"note: {note}")
        # `minimal` and `contradicted` are counts over distinct rules, so the
        # rules figure beside them has to be too, not the input row count.
        print(f"\n{len(result.annotations):,} loci, "
              f"{result.summary.get('rules', len(result.rules)):,} rules "
              f"({result.summary.get('minimal', 0):,} minimal, "
              f"{result.summary.get('contradicted_by_a_superset', 0):,} contradicted by a "
              f"superset)")
        for name, count in result.summary["by_primary_signature"].items():
            print(f"  {name:28} {count:5,}")
        if result.questions:
            print(f"\n{len(result.questions):,} question(s) the sources could not answer, "
                  f"over {result.summary['distinct_loci']:,} loci")
            for kind, count in result.summary["by_kind"].items():
                print(f"  {kind:28} {count:5,}")
        for label, path in result.written.items():
            print(f"\n{label}: {path}")
        return 0

    if args.mode == "features":
        from kegg_string_mcp.rules import parse_annotation, parse_rules, resolve_all, vocabulary
        from kegg_string_mcp.rules.report import render, write_tsv

        if len(args.genes) != 1:
            parser.error("features mode takes exactly one labels or rules file")
        if args.annotation is None:
            parser.error("features mode needs --annotation (or $KEGG_STRING_MCP_ANNOTATION)")

        source = Path(args.genes[0])
        try:
            # A rule file states its own vocabulary, and that is the one that has
            # to resolve: a separate labels list may name features no rule uses.
            labels = vocabulary(parse_rules(source))
            origin = f"{len(labels):,} distinct labels used by the rules in {source}"
        except (ValueError, UnicodeDecodeError):
            labels = [line.strip() for line in source.read_text().splitlines() if line.strip()]
            origin = f"{len(labels):,} labels listed in {source}"

        annotation = parse_annotation(args.annotation)
        coverage = resolve_all(labels, annotation)
        print(origin + "\n")
        print(render(coverage, annotation))
        print(f"\nfeatures: {write_tsv(args.out / 'features.tsv', coverage)}")
        # A report succeeds by reporting. Some labels not resolving is a finding
        # about the inputs, not a failure of this step, and failing the process
        # would stop a pipeline on something the report is there to show you.
        # `--strict` is for the caller who wants the run gated on full coverage.
        return 1 if (args.strict and coverage.unresolved()) else 0

    if args.mode == "eval":
        from kegg_string_mcp.evaluate import evaluate, render, write

        report = asyncio.run(_run_eval(evaluate, args))
        print(render(report))
        path = write(report, args.runs / "eval" / "report.json")
        print(f"\nreport: {path}")
        return 0

    if args.mode == "single":
        if len(args.genes) != 1:
            parser.error("single mode takes exactly one gene")
        payload = asyncio.run(_run(annotate_gene, args, args.genes[0]))
    else:
        if len(args.genes) < 2:
            parser.error("epistasis mode needs at least two genes")
        payload = asyncio.run(_run(annotate_epistasis, args, args.genes))

    if args.json:
        print(json.dumps(payload, indent=2, default=str))
    else:
        print(payload["summary"])
        print("\n--- validation ---")
        print(payload["validation"]["summary"])
        for citation in payload["validation"]["citations"]:
            if citation["status"] != "verified":
                print(f"  {citation['status'].upper():13} {citation['identifier']}  {citation['detail']}")
        for quote in payload["validation"].get("quotes", []):
            if quote["status"] != "verified":
                triage = f" [{quote['triage']}, similarity {quote['similarity']}]" if quote.get("triage") else ""
                # Not every quotable record is a PMID: UniProt accessions are
                # quotable too, and the hard-coded prefix printed
                # "PMID:A0A0N7EHL5" for one.
                print(f"  {quote['status'].upper():13} {quote['record_id']}{triage}")
                print(f"                  quoted:  {quote['quote'][:88]!r}")
                if quote.get("closest_span"):
                    print(f"                  closest: {quote['closest_span'][:88]!r}")
        print(f"\nrun store: {payload['store']}")

    return 0 if payload["validation"]["passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
