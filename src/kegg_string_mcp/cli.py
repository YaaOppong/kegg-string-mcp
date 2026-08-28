"""Command line entry point for the annotation pipeline.

    gar single katG
    gar epistasis katG furA ahpC
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from kegg_string_mcp.agent import annotate_epistasis, annotate_gene


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="gar", description=__doc__)
    parser.add_argument("mode", choices=["single", "epistasis"])
    parser.add_argument("genes", nargs="+")
    parser.add_argument("--organism", default="mtu")
    parser.add_argument("--runs", type=Path, default=Path("runs"))
    parser.add_argument("--json", action="store_true", help="emit the full payload")
    args = parser.parse_args(argv)

    if args.mode == "single":
        if len(args.genes) != 1:
            parser.error("single mode takes exactly one gene")
        payload = annotate_gene(args.genes[0], args.organism, args.runs)
    else:
        if len(args.genes) < 2:
            parser.error("epistasis mode needs at least two genes")
        payload = annotate_epistasis(args.genes, args.organism, args.runs)

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
                print(f"  {quote['status'].upper():13} PMID:{quote['record_id']}  "
                      f"{quote['detail']}\n                  quoted: {quote['quote'][:90]!r}")
        print(f"\nrun store: {payload['store']}")

    return 0 if payload["validation"]["passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
