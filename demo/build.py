"""Turn raw run stores into the compact records the demo replays.

Same philosophy as tests/fixtures: real output, captured once, replayed offline.
The demo makes no API calls -- visitors have no Anthropic key and must not spend
ours -- so everything it shows is committed to the repo.

Deliberately NOT captured: the validation verdict. Only the tool results and the
model's summary are stored, and the app re-runs `validate()` at display time. The
validator has been fixed many times, several of those fixes removing false
positives, so replaying a stored verdict would show a visitor conclusions the
current code no longer reaches. Re-validating is also cheap, deterministic, and
means the demo cannot drift from the library it is demonstrating.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any

RUNS = Path(__file__).parent / "runs"


def read_store(path: Path) -> dict[str, Any]:
    """Extract the replayable parts of one run store."""
    calls: list[dict[str, Any]] = []
    turns: list[dict[str, Any]] = []
    output: dict[str, Any] | None = None
    # The store is append-only and ordered, so who asked for a call is a fact in
    # the file rather than something to infer: anything before the first decision
    # was fetched by the pipeline itself, and anything after one belongs to that
    # turn. Matching on tool name and arguments instead gets this wrong whenever
    # the model happens to request a call the pipeline already made.
    current_turn: int | None = None

    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        entry = json.loads(line)
        kind = entry.get("kind")
        if kind == "tool_result":
            calls.append({"tool": entry["tool"], "arguments": entry["arguments"],
                          "result": entry["result"],
                          "origin": "pipeline" if current_turn is None else f"turn {current_turn}"})
        elif kind == "decision":
            current_turn = entry["turn"]
            turns.append({"turn": entry["turn"], "stop_reason": entry["stop_reason"],
                          "tools_requested": entry.get("tools_requested", [])})
        elif kind == "output":
            output = entry

    if output is None:
        raise ValueError(f"{path} has no output record; the run did not finish")

    return {
        "mode": output["mode"],
        "target": output.get("gene") or ", ".join(output.get("genes", [])),
        "organism": output.get("organism", "mtu"),
        "summary": output["summary"] or "",
        "turns": turns,
        "calls": calls,
    }


def build(sources: dict[str, Path]) -> None:
    RUNS.mkdir(parents=True, exist_ok=True)
    for name, path in sources.items():
        record = read_store(path)
        record["id"] = name
        out = RUNS / f"{name}.json"
        out.write_text(json.dumps(record, indent=1, sort_keys=True), encoding="utf-8")
        print(f"  {name:24} {len(record['calls'])} tool calls, "
              f"{len(record['turns'])} turns, {out.stat().st_size // 1024} KB")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit("usage: build.py NAME=PATH [NAME=PATH ...]")
    pairs = dict(arg.split("=", 1) for arg in sys.argv[1:])
    build({name: Path(path) for name, path in pairs.items()})
