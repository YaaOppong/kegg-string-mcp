# Design notes

Why the tools behave the way they do. The [README](../README.md) covers what they do.

**Empty is not the same as absent.** `records: []` always comes with a note saying
which case it is: the identifier failed to resolve, or the gene resolved and genuinely
has no pathways/partners at that threshold. Silence would invite a model to fill the gap.

**STRING's headline score already includes literature.** The combined score folds in a
textmining channel, so a high score is *not* independent of literature evidence about
the same pair. Each partner carries `max_non_textmining_score` and
`evidence_beyond_textmining`, thresholded at STRING's medium-confidence band (0.4) —
not at "greater than zero", which is meaningless here. Live example: katG–embB scores
0.964 combined, of which 0.963 is textmining and the only other non-zero channel is
coexpression at 0.044.

**Symbol resolution is exact, never fuzzy.** Gene symbols resolve against the
organism's full KEGG gene list (fetched once, cached) rather than KEGG's `/find`
endpoint, which returns loose matches and would make results non-deterministic. The
match type is reported in `resolved.matched_by`.

**Cached records report their original fetch time.** A response served from disk keeps
the timestamp of the fetch that produced it. Stamping "now" onto a three-week-old
cached record would be a false provenance claim on every cache hit.

**Rate limits are honoured, not hoped for.** Per-host throttling (KEGG ~3 req/s,
STRING ~1 req/s, NCBI ~2.8 req/s), bounded retries with `Retry-After` support,
`caller_identity` on every STRING call and `tool`/`email` on every NCBI call. An NCBI
API key raises the ceiling to 10 req/s; the floor stays put, because the extra rate is
worth less than never being the client that gets your IP banned. Only 2xx responses are
cached — caching a 500 would freeze a transient failure into a permanent wrong answer.

**Credentials never reach the provenance trail.** Caller-identifying params (`email`,
`tool`, `api_key`, `caller_identity`) are stripped from the cache key, so setting
`NCBI_EMAIL` later does not invalidate a whole cache. `api_key` goes further and is
redacted from the *audit* URL too: audit URLs travel into `requests[].url` and from
there into the run store on disk, and a credential must not be written down there.

### Structured records and prose

A KEGG pathway ID and a STRING score are *structured*: the record means one thing, and
the model's only job is to relay it. An abstract is prose, so a claim drawn from one is
the model deciding what a record says — which is what the other two tools are built to
prevent. Membership in `record_ids` cannot catch that: **a real PMID attached to a
fabricated finding passes it.**

So every article record carries `quotable_text` — the exact retrieved text, with XML
markup flattened and whitespace normalised. A claim sourced from an abstract should
quote a span of that string, and the quote is then checkable by containment. Two tiers,
both deterministic, no embeddings:

| Tier | Question | Check |
|---|---|---|
| Citation | Was this record actually retrieved? | `record_id in record_ids` |
| Span | Does the cited record really say this? | `quote in record.detail["quotable_text"]` |

`has_abstract` is `false` where PubMed holds only a title; those records are returned,
flagged, and named in `notes`, because a validator needs to know which records have
nothing to quote from. The span validator itself lives in the agent layer, not here —
this tool's job is to store the exact text it can be checked against.

**PubMed search is retrieval, not resolution.** This is the one tool with no notion of a
match failing. KEGG and STRING resolve an identifier: there is a right answer, and a
miss is reportable as a miss. `esearch` is a relevance-ranked text search that always
returns something plausible, so a returned article may merely *contain* the gene string.
Results therefore report PubMed's own `querytranslation` (what it actually searched for,
after Automatic Term Mapping) and the total hit count next to the truncated set — 1382
matched, 4 returned is a materially different result than 4 matched, 4 returned. Gene and
organism are phrase-quoted, and values containing PubMed query syntax (`"`, `[`, `]`) are
rejected rather than interpolated, since they would change the question rather than the
terms.

**Literature evidence is easy to double-count.** STRING's combined score already folds
in textmining, so a PubMed abstract about the same pair may be the very source of that
score. Citing both is one line of evidence counted twice, and every `pubmed_abstracts`
result says so.


## Deterministic where you can be, model where you must be

The pipeline fetches, computes and validates; the model chooses what to look up and
interprets what came back. It never writes to the store and never does arithmetic over
record IDs -- a model doing set arithmetic over a hundred identifiers will get some of
it wrong in a way nobody can audit.

In epistasis mode every pairwise relationship is computed before the model sees it.
Shared pathways carry their size, because size decides whether sharing one means
anything: `mtu01100` ("Metabolic pathways") holds 698 of ~4,000 *M. tuberculosis*
genes, so co-membership is a base rate; `mtu00983` holds 11, and co-membership is
signal. Partner lists that hit the retrieval limit are flagged, because without a true
network degree a shared partner cannot be distinguished from what any two
well-connected proteins would share.

## Why failed quotes are ranked, not adjudicated

Failed quotes report the closest matching span and a similarity score, which separates
a quoting artefact (~0.93) from an invented claim (~0.38). That ranks failures for
human attention; it does not change the verdict.

A model deciding when a deterministic check was wrong would invert the trust model --
the validator is deterministic precisely so a model cannot argue past it, and a checker
would need its own checker. Every false positive this validator has produced had a
deterministic cause and a deterministic fix; a tolerance layer would have absorbed them
and they would still be shipping.
