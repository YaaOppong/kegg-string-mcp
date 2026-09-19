"""The two annotation modes, and the prompts that constrain them.

Both prompts share one rule, stated in the strongest terms the model will honour:
cite only identifiers that appear in a tool result, and prefer saying nothing to
saying something uncited. That rule is not trusted -- it is *checked*, by
`validate.py`, against what the tools actually returned. The prompt exists to
make the common case correct; the validator exists because prompts are not
guarantees.

The prompt text itself lives in `skills/gene-annotation/SKILL.md`, not here. That
file is also the skill an external MCP client loads, and rules that an outside
caller follows while the pipeline is scored under a second, drifting copy would
be worse than no skill at all -- the same reason the pipeline takes its tool
schemas from `list_tools()` rather than redeclaring them. This module assembles;
it does not define.
"""

from __future__ import annotations

import functools
from pathlib import Path

MODEL = "claude-opus-5"

SKILL_RELATIVE = Path("skills") / "gene-annotation" / "SKILL.md"


def _skill_file() -> Path:
    """Packaged copy first, then the source tree, so an installed wheel and a
    `pip install -e .` checkout resolve the same text."""
    packaged = Path(__file__).resolve().parent.parent / SKILL_RELATIVE
    if packaged.exists():
        return packaged
    return Path(__file__).resolve().parents[3] / SKILL_RELATIVE


@functools.lru_cache(maxsize=1)
def _skill_text() -> str:
    return _skill_file().read_text(encoding="utf-8")


def block(name: str) -> str:
    """The text between `<!-- canonical:<name> start -->` and its end marker.

    A missing or unterminated marker raises rather than returning a truncated
    prompt: a silently shortened rule list would drop the citation rule and the
    run would look normal until the validator started firing.
    """
    start = f"<!-- canonical:{name} start -->"
    end = f"<!-- canonical:{name} end -->"
    text = _skill_text()
    try:
        body = text.split(start, 1)[1].split(end, 1)[0]
    except IndexError:
        raise ValueError(
            f"{_skill_file()} has no complete canonical:{name} block"
        ) from None
    return body.strip()


def prompt_for(mode: str) -> str:
    tail = {"single": "single", "epistasis": "epistasis"}[mode]
    return block("shared") + "\n\n" + block(tail)


def __getattr__(name: str) -> str:
    # _SHARED, SINGLE_GENE and EPISTASIS were module constants before the text
    # moved into the skill. Keep them readable as attributes so tests and the
    # evaluation harness that assert on prompt content do not have to care.
    if name == "_SHARED":
        return block("shared")
    if name == "SINGLE_GENE":
        return prompt_for("single")
    if name == "EPISTASIS":
        return prompt_for("epistasis")
    raise AttributeError(name)
