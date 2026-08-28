"""Tool dispatch backed by a real MCP session.

The agent used to define its own copy of the tool schemas and call the Python
clients directly. That worked, but it meant the pipeline was not a client of its
own MCP server -- and the two copies drifted, so the model driving the pipeline
saw *shorter* tool descriptions than an external MCP client would. It was missing,
among other things, the instruction that quotes are checked programmatically and
a paraphrase will not pass, which is the single most behaviour-shaping line in the
system.

Here the schemas come from `list_tools()` over the wire, so there is exactly one
definition and drift is structurally impossible rather than something a test has
to notice.

The server is spawned as `python -m kegg_string_mcp.server` rather than by the
console-script name: it works from a source checkout without the package being on
PATH, and it guarantees the child runs the same interpreter as the parent.
"""

from __future__ import annotations

import json
import os
import sys
from contextlib import asynccontextmanager
from typing import Any

from mcp import ClientSession, StdioServerParameters
from mcp.client.stdio import stdio_client

# Passed through to the child. The server needs these to identify itself to NCBI
# and STRING, and to share the caller's cache rather than starting cold.
FORWARDED_ENV = (
    "KEGG_STRING_MCP_CACHE", "KEGG_STRING_MCP_USER_AGENT",
    "STRING_CALLER_IDENTITY", "NCBI_EMAIL", "NCBI_API_KEY", "NCBI_TOOL",
)


class McpTools:
    """Dispatch callable with the same shape as `pipeline.Tools`, but over MCP."""

    def __init__(self, session: ClientSession, tools: list[Any]):
        self.session = session
        self._tools = tools

    def schemas(self) -> list[dict[str, Any]]:
        """Anthropic tool definitions built from what the server advertises."""
        return [{"name": t.name,
                 "description": t.description or "",
                 "input_schema": t.input_schema}
                for t in self._tools]

    @property
    def names(self) -> set[str]:
        return {t.name for t in self._tools}

    def _check(self, name: str, arguments: dict[str, Any]) -> list[str]:
        """Validate against the schema the server advertises.

        MCP silently drops arguments a tool does not declare, so a model asking for
        `organism="mtu"` on a tool that has no such parameter got results as though
        the constraint applied. The direct dispatch refuses that; both paths must
        behave the same, and deriving the check from the advertised schema means
        there is no second parameter table to keep in step.
        """
        schema = next((t.input_schema for t in self._tools if t.name == name), None)
        if schema is None:
            return [f"'{name}' is not a tool this server exposes"]
        properties = schema.get("properties", {})
        problems = [f"'{key}' is not a parameter of {name} "
                    f"(accepts: {', '.join(sorted(properties))})"
                    for key in arguments if key not in properties]
        for key in schema.get("required", []):
            if key not in arguments:
                problems.append(f"'{key}' is required by {name}")
        for key, value in arguments.items():
            expected = properties.get(key, {}).get("type")
            if expected == "integer" and not isinstance(value, int):
                problems.append(f"'{key}' must be an integer, got {value!r}")
            elif expected == "string" and not isinstance(value, str):
                problems.append(f"'{key}' must be a string, got {value!r}")
        return problems

    async def __call__(self, name: str, arguments: dict[str, Any]) -> dict[str, Any]:
        problems = self._check(name, arguments)
        if problems:
            return {"query": dict(arguments), "records": [], "record_ids": [],
                    "notes": [(f"Invalid argument(s), so no lookup was performed: "
                               f"{'; '.join(problems)}. An empty result here does NOT "
                               f"mean there is no data.")]}

        result = await self.session.call_tool(name, arguments)

        if getattr(result, "is_error", False):
            # A tool error is data the model can act on, not a crash. Same contract
            # as the direct dispatch, which returns an envelope for bad arguments.
            text = _first_text(result)
            return {"query": dict(arguments), "records": [], "record_ids": [],
                    "notes": [f"Tool '{name}' returned an error: {text[:400]}"]}

        structured = getattr(result, "structured_content", None)
        if isinstance(structured, dict):
            return structured

        # Fall back to the text content if a server ever omits structured output,
        # rather than silently handing the model an empty envelope.
        text = _first_text(result)
        try:
            parsed = json.loads(text)
        except ValueError:
            parsed = None
        if isinstance(parsed, dict):
            return parsed
        return {"query": dict(arguments), "records": [], "record_ids": [],
                "notes": [f"Tool '{name}' returned no structured content."]}


def _first_text(result: Any) -> str:
    for block in getattr(result, "content", []) or []:
        if getattr(block, "type", "") == "text":
            return block.text
    return ""


def child_env() -> dict[str, str]:
    env = {k: os.environ[k] for k in FORWARDED_ENV if k in os.environ}
    # The child is a Python process; without PATH and the venv it may not resolve.
    for key in ("PATH", "HOME", "VIRTUAL_ENV", "CONDA_PREFIX", "SYSTEMROOT"):
        if key in os.environ:
            env[key] = os.environ[key]
    return env


@asynccontextmanager
async def mcp_tools(command: str | None = None, args: list[str] | None = None):
    """Spawn the MCP server and yield a dispatch bound to a live session."""
    params = StdioServerParameters(
        command=command or sys.executable,
        args=args if args is not None else ["-m", "kegg_string_mcp.server"],
        env=child_env(),
    )
    async with stdio_client(params) as (read, write), ClientSession(read, write) as session:
        await session.initialize()
        listed = await session.list_tools()
        yield McpTools(session, list(listed.tools))
