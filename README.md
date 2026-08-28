# kegg-string-mcp

An [MCP](https://modelcontextprotocol.io) server exposing **KEGG** and **STRING** as
model-callable tools for gene annotation. Every record returned carries a stable ID
and a resolvable URL, so a downstream agent's citations can be checked
programmatically against what was actually retrieved.

Defaults target *Mycobacterium tuberculosis* H37Rv (KEGG `mtu`, NCBI taxon `83332`),
but both tools take the organism as a parameter and work for any species in KEGG and
STRING.

## Tools

| Tool | Parameters | Returns |
|---|---|---|
| `kegg_pathways` | `gene` (KEGG ID, locus tag, or symbol)<br>`organism` = `mtu` | One record per KEGG pathway: `record_id` (e.g. `mtu00360`), name, `https://www.kegg.jp/entry/…` URL |
| `string_partners` | `gene`<br>`species` = `83332`<br>`limit` = `20`<br>`required_score` = `700` | One record per interaction partner: `record_id` (e.g. `83332.Rv1909c`), preferred name, combined score, full per-channel breakdown, `https://string-db.org/network/…` URL |

Both are annotated `readOnlyHint`, expose structured output schemas, and are
deterministic: same input and same cache produce the same output.

### Result envelope

```jsonc
{
  "query":      { "gene": "katG", "organism": "mtu" },
  "resolved":   { "kegg_gene_id": "mtu:Rv1908c", "matched_by": "locus_tag_or_symbol" },
  "records":    [ { "record_id": "mtu00360", "name": "…", "url": "…",
                    "retrieved_at": "…", "cached": true, "detail": { } } ],
  "record_ids": ["mtu00360", "mtu00380", "mtu00983", "mtu01100", "mtu01110"],
  "notes":      [ "…" ],
  "requests":   [ { "url": "…", "retrieved_at": "…", "cached": true,
                    "status": 200, "content_sha256": "…" } ]
}
```

`record_ids` is the flat citable list. A citation validator checks membership in this
set — a set test, rather than walking nested records, so it is hard to get subtly wrong.

## Design notes

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
STRING ~1 req/s), bounded retries with `Retry-After` support, and `caller_identity` on
every STRING call. Only 2xx responses are cached — caching a 500 would freeze a
transient failure into a permanent wrong answer.

## Install

```bash
conda create -n kegg-string-mcp python=3.11 -y && conda activate kegg-string-mcp
pip install -e ".[dev]"
pytest
```

## Run

```bash
kegg-string-mcp          # stdio transport
```

Wire into an MCP client:

```jsonc
{
  "mcpServers": {
    "kegg-string": {
      "command": "kegg-string-mcp",
      "env": {
        "STRING_CALLER_IDENTITY": "your-name-or-project",
        "KEGG_STRING_MCP_CACHE": "~/.cache/kegg-string-mcp"
      }
    }
  }
}
```

| Variable | Default | Purpose |
|---|---|---|
| `KEGG_STRING_MCP_CACHE` | `~/.cache/kegg-string-mcp` | Response cache directory |
| `STRING_CALLER_IDENTITY` | `kegg-string-mcp` | Sent to STRING on every call, as they ask |
| `KEGG_STRING_MCP_USER_AGENT` | `kegg-string-mcp/0.1` | HTTP User-Agent |

## Tests

56 tests. **The suite never touches the network** — `tests/conftest.py` swaps in a fake
HTTP client that replays saved responses from `tests/fixtures/`, so the suite runs in
under two seconds, gives the same answer every time, works offline and in CI, and does
not hammer a free academic service. A red test means this code broke, not that an
upstream service was down.

The fixtures are *captured*, not hand-written: verbatim KEGG and STRING output pulled
from the live endpoints on 2026-08-27. That is how the four-column `/list/{organism}`
layout was caught — an invented fixture would have encoded the wrong assumption and
passed. The trade-off is that fixtures cannot detect an upstream format *change*; that
needs an occasional live check.

```bash
pytest -q
```

## Licence

This project is MIT licensed — see [LICENSE](LICENSE).

### Upstream data

This server queries public APIs on the caller's behalf and caches responses locally;
it redistributes no third-party data. Terms remain the caller's responsibility:

- **KEGG** — free for academic use; commercial use requires a licence from Pathway Solutions.
- **STRING** — CC BY 4.0, free for academic and commercial use, attribution required.

## Status

Phase 1 of a larger project. Next: an agent loop that uses these tools to annotate
genes, an append-only store written by the pipeline rather than the model, and a
validation layer that checks every citation in a generated summary against the
`record_ids` actually retrieved.
