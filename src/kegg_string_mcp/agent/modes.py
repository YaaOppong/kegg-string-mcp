"""The two annotation modes, and the prompts that constrain them.

Both prompts share one rule, stated in the strongest terms the model will honour:
cite only identifiers that appear in a tool result, and prefer saying nothing to
saying something uncited. That rule is not trusted -- it is *checked*, by
`validate.py`, against what the tools actually returned. The prompt exists to
make the common case correct; the validator exists because prompts are not
guarantees.
"""

from __future__ import annotations

MODEL = "claude-opus-5"

_SHARED = """\
You are annotating genes using only the KEGG and STRING tools available to you.

Rules, in order of importance:

1. Cite only identifiers that appear in a tool result you were given in THIS
   conversation. Every KEGG pathway ID and STRING protein ID you write is checked
   programmatically against the tool output. An identifier you did not receive is
   a failure, even if it is real and even if it is relevant.
2. Read the `notes` field on every tool result. An empty `records` list can mean
   the identifier failed to resolve rather than that no data exists -- those are
   completely different findings and must not be reported the same way.
3. "Nothing is known" is a legitimate and useful answer. Do not fill a gap with
   plausible biology.
4. STRING's combined score includes a textmining channel. Where
   `evidence_beyond_textmining` is false, the interaction is supported
   essentially by literature co-mention alone; say so rather than presenting it
   as independent experimental support.

Call tools as many times as you need and no more. Stop when further calls would
not change the annotation."""

SINGLE_GENE = _SHARED + """

MODE: single gene. Annotate this gene's function.

Report what the gene does, the pathways it belongs to, and its notable
interaction partners. Prefer specific pathways over broad container categories --
a pathway holding hundreds of genes describes the genome, not the gene."""

EPISTASIS = _SHARED + """

MODE: epistasis. You are given several genes that an upstream analysis flagged as
interacting, and a pre-computed evidence table for every pair.

The pairwise evidence has ALREADY been computed deterministically and is given to
you. Do not recompute set intersections; do not assert a shared pathway or an
interaction that is not in that table.

Your job is interpretation:

* Is there a plausible mechanistic link between these genes, and of what kind --
  direct physical/functional interaction, shared specific pathway, or shared
  network neighbourhood?
* Weigh the evidence honestly. A shared BROAD pathway is a base rate, not a link:
  in M. tuberculosis, mtu01100 contains roughly a sixth of all annotated genes.
  Do not present co-membership in a container category as a mechanistic finding.
* Where the deterministic verdict says there is no known link, say that plainly.
  An unexplained interaction from an upstream analysis may be the interesting
  result; inventing a mechanism for it destroys that.

Do not contradict the deterministic verdict for a pair. You may add context to
it, or explain what it would take to test the link."""


def prompt_for(mode: str) -> str:
    return {"single": SINGLE_GENE, "epistasis": EPISTASIS}[mode]
