"""The agent loop. Hand-written on purpose.

The model decides which tools to call, in what order, and when it has enough --
a well-characterised gene takes one call, an obscure one takes several. Every one
of those decisions is logged to the run store, and that log is the evidence that
the orchestration is real rather than asserted. A framework would supply the loop
and take the log with it.

Division of labour, which is the whole design:

  pipeline (code)  fetches, computes set operations, writes the store, validates
  model            chooses what to look up, and interprets what came back

The model never writes to the store and never does arithmetic over record IDs.
"""

from __future__ import annotations

import json
import uuid
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable

import anthropic

from kegg_string_mcp.agent.modes import MODEL, prompt_for
from kegg_string_mcp.agent.store import RunStore

MAX_TURNS = 12
MAX_TOKENS = 16000

TOOL_SCHEMAS = [
    {
        "name": "kegg_pathways",
        "description": (
            "KEGG pathways for one gene. Accepts a KEGG gene ID (mtu:Rv1908c), a locus tag "
            "(Rv1908c), or a gene symbol (katG). Returns one record per pathway, with the KEGG "
            "pathway ID as record_id. If records is empty, read notes: it distinguishes "
            "'identifier did not resolve' from 'gene exists but maps to no pathway'."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "gene": {"type": "string"},
                "organism": {"type": "string", "default": "mtu"},
            },
            "required": ["gene"],
        },
    },
    {
        "name": "string_partners",
        "description": (
            "STRING interaction partners for one gene, with the full per-channel score "
            "breakdown. STRING's combined_score includes textmining, so check "
            "evidence_beyond_textmining before calling an interaction experimentally supported. "
            "An empty result at a high required_score is not evidence of isolation."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "gene": {"type": "string"},
                "species": {"type": "integer", "default": 83332},
                "limit": {"type": "integer", "default": 20},
                "required_score": {"type": "integer", "default": 700},
            },
            "required": ["gene"],
        },
    },
    {
        "name": "pubmed_abstracts",
        "description": (
            "PubMed abstracts for a gene, via NCBI E-utilities. Use this LAST and only when the "
            "structured tools leave a question open -- it is the weakest and noisiest channel. "
            "Unlike KEGG and STRING this is relevance-ranked text search, NOT identifier "
            "resolution: matching a query string is not the same as being about the gene, and "
            "there is no 'did not resolve' signal. Each record carries the exact retrieved text; "
            "any claim you draw from an abstract must quote a verbatim span of it."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "gene": {"type": "string"},
                "organism": {"type": "string", "default": "Mycobacterium tuberculosis"},
                "limit": {"type": "integer", "default": 10},
            },
            "required": ["gene"],
        },
    },
]


@dataclass
class LoopResult:
    text: str
    turns: int
    tool_calls: list[dict[str, Any]] = field(default_factory=list)
    stop_reason: str = ""
    usage: dict[str, int] = field(default_factory=dict)


def run_loop(
    mode: str,
    task: str,
    dispatch: Callable[[str, dict[str, Any]], dict[str, Any]],
    store: RunStore,
    client: Any | None = None,
    max_turns: int = MAX_TURNS,
) -> LoopResult:
    """Drive the model until it stops asking for tools.

    `dispatch` runs a tool and returns its envelope; injecting it keeps this loop
    testable without a network or an API key.
    """
    client = client or anthropic.Anthropic()
    messages: list[dict[str, Any]] = [{"role": "user", "content": task}]
    calls: list[dict[str, Any]] = []
    result = LoopResult(text="", turns=0)

    for turn in range(1, max_turns + 1):
        response = client.messages.create(
            model=MODEL,
            max_tokens=MAX_TOKENS,
            system=prompt_for(mode),
            tools=TOOL_SCHEMAS,
            thinking={"type": "adaptive"},
            messages=messages,
        )
        result.turns = turn
        result.stop_reason = response.stop_reason or ""

        requested = [b for b in response.content if getattr(b, "type", "") == "tool_use"]
        store.decision({
            "turn": turn,
            "stop_reason": response.stop_reason,
            "tools_requested": [{"name": b.name, "input": b.input} for b in requested],
            "usage": {"input": response.usage.input_tokens, "output": response.usage.output_tokens},
        })

        if response.stop_reason != "tool_use":
            result.text = "".join(b.text for b in response.content if getattr(b, "type", "") == "text")
            result.usage = {"input": response.usage.input_tokens,
                            "output": response.usage.output_tokens}
            result.tool_calls = calls   # was set only on the max-turns path
            return result

        messages.append({"role": "assistant", "content": response.content})

        # All tool results for one assistant turn go back in ONE user message.
        # Splitting them across messages teaches the model to stop batching calls.
        tool_results = []
        for block in requested:
            envelope = dispatch(block.name, dict(block.input))
            store.tool_result(block.name, dict(block.input), envelope)
            calls.append({"turn": turn, "tool": block.name, "input": dict(block.input),
                          "n_records": len(envelope.get("records", []))})
            tool_results.append({
                "type": "tool_result",
                "tool_use_id": block.id,
                "content": json.dumps(envelope, default=str),
            })
        messages.append({"role": "user", "content": tool_results})

    result.text = ""
    result.stop_reason = "max_turns_exhausted"
    result.tool_calls = calls
    return result


def new_store(root: Path, label: str) -> RunStore:
    run_id = f"{label}-{uuid.uuid4().hex[:8]}"
    return RunStore(path=root / f"{run_id}.jsonl", run_id=run_id)
