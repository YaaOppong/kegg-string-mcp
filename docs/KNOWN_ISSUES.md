# Known issues

Things found, understood, and deliberately not fixed. Each says what breaks, how
it was found, and what the fix would be, so none of it has to be rediscovered.

Nothing is currently open.

---

## Resolved

**Cross-source locus fallback could pick the wrong paralogue** —
`retrieval/coverage.py`, `_locus()`. When one source failed to resolve a symbol,
the retry took its locus tag from the first record the other source returned.
That record is not guaranteed to be the intended gene: a UniProt symbol search
matches several entries, and an unreviewed fragment or a paralogue sorting first
would redirect the lookup, producing a coverage verdict about a different gene.

Measured before fixing rather than assumed. `_locus()` already skipped records
with no locus tag, and the unreviewed TrEMBL entries padding a symbol search
carry none — katG returns three entries and only the reviewed P9WIE5 has
`Rv1908c`. Across 17 genes including paralogue-prone families (esxA, esxB,
mce1A, PE_PGRS56, mmpL9), **zero** returned more than one tagged record, so the
ambiguity never arose in practice.

Fixed deterministically anyway: the reviewed entry breaks a tie (SwissProt versus
TrEMBL is exactly that distinction and is already parsed), and an unbroken tie is
refused rather than guessed. The gene is flagged `needs_review` with the
competing tags in `locus_candidates`, routed nowhere, and listed in the coverage
summary for a human to settle.

Adjudicating the tie with a model was considered and rejected. The stage 1 tools
are advertised as deterministic — same input and same cache produce the same
output — which is what makes the run store auditable and the cache replayable;
a model in the middle of `_locus()` would break that for a decision with no
judgement in it, since `reviewed` and `locus_tags` already carry the answer. The
flag is the handoff instead: the agent layer can reason about an ambiguity the
tool refuses to resolve, which is the division the repo already uses.
