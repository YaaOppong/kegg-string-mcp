# kegg-string-mcp

An [MCP](https://modelcontextprotocol.io) server exposing **KEGG**, **STRING** and
**PubMed** as model-callable tools for gene annotation. Every record returned carries a
stable ID and a resolvable URL, so a downstream agent's citations can be checked
programmatically against what was actually retrieved.

Defaults target *Mycobacterium tuberculosis* H37Rv (KEGG `mtu`, NCBI taxon `83332`),
but every tool takes the organism as a parameter and works for any species the
underlying service covers.

## Tools

| Tool | Parameters | Returns |
|---|---|---|
| `kegg_pathways` | `gene` (KEGG ID, locus tag, or symbol)<br>`organism` = `mtu` | One record per KEGG pathway: `record_id` (e.g. `mtu00360`), name, `https://www.kegg.jp/entry/…` URL |
| `string_partners` | `gene`<br>`species` = `83332`<br>`limit` = `20`<br>`required_score` = `700` | One record per interaction partner: `record_id` (e.g. `83332.Rv1909c`), preferred name, combined score, full per-channel breakdown, `https://string-db.org/network/…` URL |
| `pubmed_abstracts` | `gene`<br>`organism` = `Mycobacterium tuberculosis`<br>`limit` = `10` | One record per article: `record_id` (PMID, e.g. `35038342`), title, abstract, `quotable_text`, journal, year, DOI, `https://pubmed.ncbi.nlm.nih.gov/…` URL |

All three are annotated `readOnlyHint`, expose structured output schemas, and are
deterministic: same input and same cache produce the same output.

### Result envelope

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

`record_ids` is the flat citable list. A citation validator checks membership in this
set — a set test, rather than walking nested records, so it is hard to get subtly wrong.

For PubMed records that check is necessary but not sufficient, so each article also
carries `quotable_text`: the exact retrieved title-and-abstract string. See
[Structured records and prose](#structured-records-and-prose) below.

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

## Annotation pipeline

On top of the MCP tools sits an agent loop with two modes, a validator, and an
evaluation harness.

```bash
gar single katG                  # annotate one gene's function
gar epistasis katG furA ahpC     # look for mechanistic links between genes
gar eval                         # score the pipeline against a gold set
```

The organising principle is **deterministic where you can be, model where you must
be**. The pipeline fetches, computes and validates; the model chooses what to look
up and interprets what came back. It never writes to the store and never does
arithmetic over record IDs.

**Epistasis mode** computes every pairwise relationship — direct interaction, shared
pathways, shared network partners — before the model sees it. Shared pathways carry
their size, because size decides whether sharing one means anything: `mtu01100`
("Metabolic pathways") holds 698 of ~4,000 *M. tuberculosis* genes, so co-membership
is a base rate; `mtu00983` holds 11, and co-membership is signal.

### Run store

One append-only JSONL file per run, written by the pipeline and never by the model:
every tool call and its full result, every deterministic computation, and every turn
of the agent loop — which tools were chosen, and why it stopped.

```
turn 1: stop=tool_use   tools=[kegg_pathways(katG), kegg_pathways(ahpC)]
turn 2: stop=tool_use   tools=[string_partners(katG), string_partners(ahpC)]
turn 3: stop=tool_use   tools=[pubmed_abstracts("ahpC katG")]
turn 4: stop=end_turn
```

### Citation validation

Every identifier in the summary is checked against what the tools actually returned,
recorded in the store before the model saw it. Three failure classes, because they
mean different things:

- **unsupported** — a well-formed, authoritative-looking ID that no tool returned.
- **cross-target** — an ID that *was* retrieved, but for a different gene than the
  sentence attributes it to. Invisible to a global membership check.
- **quote not in source** — a claim drawn from an abstract whose quoted span does not
  appear in the retrieved text. Catches a real PMID carrying a fabricated finding.

The check is set membership and string containment, deliberately not similarity
scoring: a citation either names a record a tool returned or it does not.

```
NOT_IN_SOURCE  PMID:10609885  [likely_fabricated, similarity 0.38]
               quoted:  'KatG binds directly to AhpC in a stable complex'
               closest: 'katg is a catalase-peroxidase required for isoniazid activation'
```

Failed quotes report the closest matching span and a similarity score, which
separates a quoting artefact (~0.93) from an invented claim (~0.38). That ranks
failures for attention; it does not change the verdict. A model must not be able to
argue its way past a deterministic check.

### Corpus manifest

Runs that find papers emit `<run>.corpus.jsonl` for a downstream full-text pipeline:

```jsonc
{"pmid": "16870753", "pmcid": "PMC1538650", "in_pmc": true,
 "doi": "10.1128/AAC.00161-06", "mentions": ["katG", "ahpC"]}
```

`in_pmc` matters more than the DOI — a DOI resolves to a usually-paywalled publisher,
while a PMCID is the licit route to full text. `mentions` records the genes actually
present in the retrieved text, not the genes queried: PubMed searches `[All Fields]`,
which spans MeSH terms the model never sees, so a paper can be returned for a gene it
never names.

## Evaluation

**The reference is incomplete, and that is the point.** KEGG assigns a pathway to just
**1,171 of 4,008** *M. tuberculosis* genes — 29%. `gyrA`, one of the most studied genes
in TB, has none. So the gold set has two classes:

- **positive controls** — KEGG assigns pathways. Measures whether the pipeline
  faithfully reports what its tools returned.
- **negative controls** — KEGG assigns none, so the correct answer is "nothing found".
  Measures whether the pipeline abstains or invents. For this project, the more
  important number.

```
gene       kind expected  reported  hits  missed  cites     quotes
katG       pos  5         5         5     0       26/26     0/0
inhA       pos  3         3         3     0       27/27     3/3
rpoB       pos  1         1         1     0       25/25     3/3
gyrA       neg  0         0         0     0       25/25     7/7

Retrieval fidelity   recall 1.0   precision 1.0
Abstention           1.0
Citation integrity   citation precision 1.0   quote precision 1.0
```

Scoring reuses the validator's own `cross_target` judgement rather than re-deriving
intent from prose: a model annotating `furA` may legitimately look up its neighbour
`katG` and discuss katG's pathways as context, and counting those as claims about
furA scored a correct annotation as a fabrication.

Citation and quote precision are the numbers worth trusting — computed, not judged.

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
        "NCBI_EMAIL": "you@example.org",
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
| `NCBI_EMAIL` | *(unset)* | Sent to NCBI on every call, as they ask. Set it — they contact you before blocking |
| `NCBI_API_KEY` | *(unset)* | Optional NCBI key. Never cached, never written to the provenance trail |
| `NCBI_TOOL` | `kegg-string-mcp` | Tool name sent to NCBI |
| `KEGG_STRING_MCP_USER_AGENT` | `kegg-string-mcp/0.1` | HTTP User-Agent |

## Tests

173 tests. **The suite never touches the network** — `tests/conftest.py` swaps in a fake
HTTP client that replays saved responses from `tests/fixtures/`, so the suite runs in
under four seconds, gives the same answer every time, works offline and in CI, and does
not hammer a free academic service. A red test means this code broke, not that an
upstream service was down.

The fixtures are *captured*, not hand-written: verbatim KEGG, STRING and E-utilities
output pulled from the live endpoints on 2026-08-27/28, with each `efetch` fixture
holding the real response for exactly the PMIDs its paired `esearch` returned, so the
search → fetch path replays end to end. That is how the four-column `/list/{organism}`
layout was caught, and how two PubMed parsing traps were: titles and abstracts carry
`<i>`/`<sup>` markup, so reading `.text` instead of `itertext()` truncates an abstract
at its first italicised species name; and every entry in an article's `<ReferenceList>`
carries its own `<ArticleIdList>` — 35 of them in one fixture article — so a descendant
search for the DOI is right only by document order, and stops being right for an article
with no DOI of its own. An invented fixture would have encoded the wrong assumption and
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
- **PubMed** — records are US government works and free to use; the abstracts themselves
  are frequently under publisher copyright. This server retrieves them per query and
  caches locally for the caller; it redistributes nothing. Respect NCBI's
  [E-utilities usage policy](https://www.ncbi.nlm.nih.gov/books/NBK25497/) — set
  `NCBI_EMAIL`.

## Status

MCP server, annotation pipeline and evaluation all working, on `main`.

Deliberately out of scope: literature RAG. No embeddings, no vector store, no document
ingestion. The literature payload for a gene is ~3,000 tokens, so there is nothing to
retrieve *from* — this is structured, tool-based retrieval, and that is what it claims
to be.
