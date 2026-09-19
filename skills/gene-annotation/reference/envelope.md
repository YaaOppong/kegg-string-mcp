# The result envelope

Every tool returns the same shape:

```jsonc
{
  "query":      { "gene": "katG", "organism": "mtu" },
  "resolved":   { "kegg_gene_id": "mtu:Rv1908c", "matched_by": "symbol" },
  "records":    [ { "record_id": "mtu00360", "name": "…", "url": "…",
                    "retrieved_at": "…", "cached": true, "detail": { } } ],
  "record_ids": ["mtu00360", "mtu00380", "mtu00983", "mtu01100", "mtu01110"],
  "notes":      [ "…" ],
  "requests":   [ { "url": "…", "retrieved_at": "…", "cached": true,
                    "status": 200, "content_sha256": "…" } ]
}
```

## resolved

What the gene name resolved to, and how. `matched_by` distinguishes a symbol
match from a locus-tag or alias match. Identity is resolved once, in one place,
so every tool in a run is answering about the same gene.

## record_ids

The flat citable list — the only identifiers you may write for this call. The
validator tests membership in this set rather than walking nested records,
because a set test is hard to get subtly wrong.

Membership is necessary but not sufficient for PubMed: a real PMID can carry a
fabricated finding, which is what `quotable_text` exists to catch.

## quotable_text

The exact retrieved text, held on PubMed, UniProt and corpus records. Quotes are
checked by string containment against this field. Quote from it verbatim —
including its capitalisation and hyphenation — or the check fires.

## notes

Read it on every result. An empty `records` list is a legitimate, explicit
answer, and `notes` says which kind:

- the identifier failed to resolve — a lookup failure, not absent data;
- the source was queried and holds nothing for this gene;
- the corpus does not cover the query (`corpus_search`);
- a limit was hit, so what is missing may sit below the cut.

Reporting these four the same way is the most common error this skill exists to
prevent.

## requests

Provenance for every HTTP call behind the result: URL, timestamp, cache status,
status code, and a SHA-256 of the response body. Not for citation — for
reproducing the run.
