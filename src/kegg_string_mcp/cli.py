"""Command line entry point for the annotation pipeline.

    gar single katG
    gar epistasis katG furA ahpC
    gar eval
"""

from __future__ import annotations

import argparse
import asyncio
import json
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
    parser.add_argument("mode", choices=["single", "epistasis", "eval"])
    parser.add_argument("genes", nargs="*")
    parser.add_argument("--organism", default="mtu")
    parser.add_argument("--runs", type=Path, default=Path("runs"))
    parser.add_argument("--json", action="store_true", help="emit the full payload")
    parser.add_argument("--direct", action="store_true",
                        help="dispatch tools in-process instead of over MCP (debugging)")
    args = parser.parse_args(argv)

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
                print(f"  {quote['status'].upper():13} PMID:{quote['record_id']}{triage}")
                print(f"                  quoted:  {quote['quote'][:88]!r}")
                if quote.get("closest_span"):
                    print(f"                  closest: {quote['closest_span'][:88]!r}")
        print(f"\nrun store: {payload['store']}")

    return 0 if payload["validation"]["passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
