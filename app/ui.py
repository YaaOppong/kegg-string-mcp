"""The page layout. One definition, used by both the server build and the static one.

app.py and page.py were near-duplicates differing only in an import line, which is
how the tool schemas drifted earlier in this project. `demo/build_pages.py` rewrites
the import for the browser build instead.
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

This annotates *Mycobacterium tuberculosis* genes. It hands a language model a set of
lookup tools — KEGG for pathways, STRING for protein interactions, PubMed for
literature — lets it decide which to call and when it has enough, then writes up what
it found.

The interesting part is the last step. **Every identifier in the write-up is checked
against what the tools actually returned.** Not by another model judging it, but by
looking: was this exact record retrieved, for this exact gene, and where text is
quoted, do those words really appear in the source? A model can produce a fluent,
plausible, correctly-formatted citation for something it was never shown. This catches
that.

Nothing here is live — these are real runs captured earlier and replayed, so the page
needs no account and costs nothing. The checking runs fresh each time you load a run.

[Source code, and the reasoning behind every design decision]({REPO})
"""

FAILED = """<div style="border-left:5px solid #b3261e;background:#fdecea;padding:1rem 1.25rem;
border-radius:6px">
<strong>⚠️ On this run, the checking caught a bad citation.</strong><br/>
The write-up below cites records that were <strong>not</strong> retrieved for this gene.
The model produced real, correctly-formatted identifiers for things it was never shown
for this gene — that is the failure this project exists to catch, and it is caught here
on a real run. Details in section 3.
</div>"""

CLEAN = """<div style="border-left:5px solid #1e8e3e;background:#e9f7ef;padding:1rem 1.25rem;
border-radius:6px">
<strong>✅ On this run, every citation checked out.</strong><br/>
Each identifier in the write-up names a record a tool actually returned for this gene,
and every quoted span appears verbatim in its source. Try one of the genes marked
⚠️ to see the checking catch something.
</div>"""


def _label(name: str) -> str:
    """Mark, in the picker itself, which runs contain a caught failure. Without it a
    visitor has to select the right gene to find the point of the demo, and most
    will not."""
    base = LABELS.get(name, name)
    return f"⚠️ {base}  — citation check FAILS" if not load(name).clean else base


def render(choice: str):
    replay = load(choice)
    model_calls = sum(1 for c in replay.calls if c.get("origin", "").startswith("turn"))
    pipeline_calls = len(replay.calls) - model_calls
    note = (f"### The model made {model_calls} tool call(s) across {len(replay.turns)} turn(s)"
            + (f", after the pipeline fetched {pipeline_calls} deterministically"
               if pipeline_calls else ""))
    quotes = quote_rows(replay)
    return (CLEAN if replay.clean else FAILED, note, tool_call_rows(replay), replay.summary,
            citation_rows(replay), gr.update(visible=bool(quotes), value=quotes))


def build() -> gr.Blocks:
    with gr.Blocks(title="Annotation with checked citations",
                   theme=gr.themes.Soft()) as page:
        gr.Markdown(INTRO)

        choices = [(_label(name), name) for name in ORDERED if name in available()]
        if not choices:
            gr.Markdown("**No runs found.** `demo/runs/` is empty — see the repository.")
            return page

        picker = gr.Dropdown(choices=choices, value=choices[0][1], label="Gene")

        # The verdict sits directly under the picker, above everything else. It is
        # the whole demonstration; putting it below the tool-call table meant a
        # visitor had to scroll to find out whether anything was caught.
        verdict = gr.HTML()

        gr.Markdown("## 1. What the model asked for")
        gr.Markdown(
            "The model chooses which tools to call and when it has enough. A well-annotated "
            "gene takes one round; a sparsely annotated one takes several, and may fall back "
            "to literature. Calls marked *pipeline* were made deterministically by the code, "
            "not chosen by the model."
        )
        turns_note = gr.Markdown()
        calls = gr.Dataframe(
            headers=["requested by", "tool", "arguments", "returned", "note from the tool"],
            wrap=True, interactive=False,
        )

        gr.Markdown("## 2. What it wrote")
        summary = gr.Markdown()

        gr.Markdown("## 3. What survived checking")
        citations = gr.Dataframe(headers=["identifier", "kind", "status", "detail"],
                                 wrap=True, interactive=False)
        quotes = gr.Dataframe(
            headers=["record", "status", "quoted in the write-up", "closest text in the source"],
            wrap=True, interactive=False, visible=False, label="Quoted passages",
        )

        gr.Markdown(
            f"**CROSS-TARGET** means the record was retrieved during this run, but for a "
            f"different gene than the sentence attributes it to — a real identifier in the "
            f"wrong place, which a check of 'does this ID exist?' would miss. **UNSUPPORTED** "
            f"means no tool returned it at all.\n\n"
            f"Research use only. Not for clinical decisions. [Repository]({REPO})"
        )

        outputs = [verdict, turns_note, calls, summary, citations, quotes]
        picker.change(render, picker, outputs)
        page.load(lambda: render(choices[0][1]), None, outputs)
    return page
