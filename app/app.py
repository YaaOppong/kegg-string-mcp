"""Gradio replay viewer for the kegg-string-mcp annotation pipeline.

Runs entirely from committed fixtures: no API key, no model call, no network.
See app/replay.py for why the validation verdict is recomputed rather than stored.
"""

from __future__ import annotations

import gradio as gr

from app.replay import (
    LABELS,
    ORDERED,
    available,
    citation_rows,
    load,
    quote_rows,
    tool_call_rows,
)

REPO = "https://github.com/YaaOppong/kegg-string-mcp"

INTRO = f"""
# Does the annotation actually say what its sources say?

This is a tool that annotates *Mycobacterium tuberculosis* genes. It gives a language
model a set of lookup tools — KEGG for pathways, STRING for protein interactions,
PubMed for literature, UniProt for protein function — lets it decide which to call,
and then writes up what it found.

The interesting part is the last step. **Every identifier in the write-up is checked
against what the tools actually returned.** Not by another model judging it, but by
looking: was this exact record retrieved, for this exact gene, and if the text is
quoted, do those words really appear in the source? A model can produce a fluent,
plausible, correctly-formatted citation for something it was never shown. This catches
that.

Pick a gene below to watch one run. Nothing here is live — these are real runs captured
earlier and replayed, so the page needs no account and costs nothing. The checking,
though, runs fresh every time you load a run.

[Source code, and the reasoning behind every design decision]({REPO})
"""

VERDICT_CLEAN = """### ✅ Every citation checked out

Each identifier in the write-up above names a record a tool actually returned for this
gene, and every quoted span appears verbatim in the source it is attributed to."""

VERDICT_FAILED = """### ⚠️ The checking caught something

The write-up cites records that were **not** retrieved for this gene. See the table
below — this is the failure the project exists to catch, on a real run."""


def render(choice: str):
    replay = load(choice)
    verdict = VERDICT_CLEAN if replay.clean else VERDICT_FAILED
    quotes = quote_rows(replay)
    model_calls = sum(1 for c in replay.calls if c.get("origin", "").startswith("turn"))
    pipeline_calls = len(replay.calls) - model_calls
    note = (f"### The model made {model_calls} tool call(s) across {len(replay.turns)} turn(s)"
            + (f", after the pipeline fetched {pipeline_calls} deterministically"
               if pipeline_calls else ""))
    return (
        note,
        tool_call_rows(replay),
        replay.summary,
        verdict,
        citation_rows(replay),
        gr.update(visible=bool(quotes), value=quotes),
    )


def build() -> gr.Blocks:
    with gr.Blocks(title="Annotation with checked citations") as page:
        gr.Markdown(INTRO)

        choices = [(LABELS.get(name, name), name) for name in ORDERED if name in available()]
        picker = gr.Dropdown(choices=choices, value=choices[0][1], label="Gene")

        gr.Markdown("## 1. What the model asked for")
        gr.Markdown(
            "The model chooses which tools to call and when it has enough. A "
            "well-annotated gene takes one round; a sparsely annotated one takes "
            "several, and may fall back to literature."
        )
        turns_note = gr.Markdown()
        calls = gr.Dataframe(
            headers=["requested by", "tool", "arguments", "returned", "note from the tool"],
            wrap=True, interactive=False,
        )

        gr.Markdown("## 2. What it wrote")
        summary = gr.Markdown()

        gr.Markdown("## 3. What survived checking")
        verdict = gr.Markdown()
        citations = gr.Dataframe(
            headers=["identifier", "kind", "status", "detail"], wrap=True, interactive=False,
        )
        quotes = gr.Dataframe(
            headers=["record", "status", "quoted in the write-up", "closest text in the source"],
            wrap=True, interactive=False, visible=False, label="Quoted passages",
        )

        gr.Markdown(
            f"**CROSS-TARGET** means the record was retrieved during this run, but for a "
            f"different gene than the sentence attributes it to — a real identifier in the "
            f"wrong place, which a check of 'does this ID exist?' would miss. "
            f"**UNSUPPORTED** means no tool returned it at all.\n\n"
            f"Research use only. Not for clinical decisions. [Repository]({REPO})"
        )

        outputs = [turns_note, calls, summary, verdict, citations, quotes]
        picker.change(render, picker, outputs)
        page.load(lambda: render(choices[0][1]), None, outputs)
    return page


if __name__ == "__main__":
    build().launch(theme=gr.themes.Soft())
