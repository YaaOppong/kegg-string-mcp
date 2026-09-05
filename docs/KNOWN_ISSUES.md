# Known issues

Things found, understood, and deliberately not yet fixed. Each says what breaks,
how it was found, and what the fix would be, so none of it has to be rediscovered.

## Cross-source locus fallback can pick the wrong paralogue

**Where:** `retrieval/coverage.py`, `_locus()`.

**What happens.** KEGG and UniProt disagree about gene symbols — KEGG calls
Rv3133c `devR`, UniProt has no `icl1` but holds the same protein as `icl`. When
one source fails to resolve a symbol, `assess()` retries it with the locus tag
the *other* source returned. `_locus()` takes that tag from `records[0]`, the
first record the answering source produced.

That first record is not guaranteed to be the right one. A UniProt symbol search
can match several entries — one reviewed SwissProt entry plus unreviewed TrEMBL
fragments, or genuine paralogues sharing a symbol. If a fragment sorts first, the
retry runs against the wrong locus tag and the coverage verdict describes a
different gene.

**Exposure is narrower than it looks, and was measured.** `_locus()` skips
records carrying no locus tag, and the unreviewed TrEMBL entries that pad a
symbol search almost always have none: katG returns three UniProt entries and
only the reviewed P9WIE5 carries `Rv1908c`, so the fallback lands on the right
one by construction. Across 17 genes checked — including paralogue-prone
families (esxA, esxB, mce1A, PE_PGRS56, mmpL9) — **zero** had more than one
record carrying a locus tag, so the ambiguity the code would have to resolve did
not arise once.

The fallback also only fires when a source failed to resolve, so it never
overrides a lookup that worked, and `resolved_via` records the tag used — a
wrong answer is traceable after the fact rather than silent.

This is therefore a latent risk rather than an observed one. It needs UniProt to
return two tagged entries for one symbol, which is possible in principle and did
not happen in the sample. Priority accordingly low.

**The fix, when it is worth doing.** Two parts:

1. Prefer the reviewed entry. `detail["reviewed"]` is already parsed, so pick the
   first reviewed record and fall back to `records[0]` only if none is.
2. Refuse when genuinely ambiguous. If the answering source returned records
   carrying more than one distinct locus tag, do not retry at all — mark the gene
   `unknown` and report it. That matches how the rest of the module treats a
   question it cannot answer, and never guesses between paralogues.

Together these keep the fallback useful for the `dosR` / `icl1` cases it was
built for while removing the guess.
