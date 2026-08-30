# kegg-string-mcp

[![CI](https://github.com/YaaOppong/kegg-string-mcp/actions/workflows/ci.yml/badge.svg)](https://github.com/YaaOppong/kegg-string-mcp/actions/workflows/ci.yml)

An [MCP](https://modelcontextprotocol.io) server exposing **KEGG**, **STRING**, **UniProt** and
**PubMed** as model-callable tools for gene annotation. Every record returned carries a
stable ID and a resolvable URL, so a downstream agent's citations can be checked
programmatically against what was actually retrieved.

Defaults target *Mycobacterium tuberculosis* H37Rv (KEGG `mtu`, NCBI taxon `83332`),
but every tool takes the organism as a parameter and works for any species the
underlying service covers.

> **Research use only.** This is a research tool. It is not a clinical decision
> support system, is not validated for diagnostic, prognostic or treatment decisions,
> and must not be used to guide patient care. It annotates genes including the drug
> resistance loci `katG`, `inhA`, `rpoB`, `pncA` and `embB`; resistance interpretation
> for clinical purposes requires validated methods and expert review. Outputs are
> generated in part by a language model and require verification against the cited
> primary sources. Provided under the MIT licence, without warranty — see
> [LICENSE](LICENSE).

## Tools

| Tool | Parameters | Returns |
|---|---|---|
| `kegg_pathways` | `gene` (KEGG ID, locus tag, or symbol)<br>`organism` = `mtu` | One record per KEGG pathway: `record_id` (e.g. `mtu00360`), name, `https://www.kegg.jp/entry/…` URL |
| `string_partners` | `gene`<br>`species` = `83332`<br>`limit` = `20`<br>`required_score` = `700` | One record per interaction partner: `record_id` (e.g. `83332.Rv1909c`), preferred name, combined score, full per-channel breakdown, `https://string-db.org/network/…` URL |
| `pubmed_abstracts` | `gene`<br>`organism` = `Mycobacterium tuberculosis`<br>`limit` = `10` | One record per article: `record_id` (PMID, e.g. `35038342`), title, abstract, `quotable_text`, journal, year, DOI, `https://pubmed.ncbi.nlm.nih.gov/…` URL |
| `uniprot_protein` | `gene`<br>`organism_id` = `83332`<br>`limit` = `3` | One record per UniProt entry: `record_id` (accession, e.g. `P9WG47`), protein name, function statements tiered by evidence code with supporting PMIDs, catalytic activity, PDB cross-refs, `quotable_text`, `https://www.uniprot.org/uniprotkb/…` URL |

All four are annotated `readOnlyHint`, expose structured output schemas, and are
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

## Demo

A zero-config replay viewer lives in [`app/`](app/): pick a gene, watch the tool calls,
read the write-up, then watch every citation in it get checked. It replays runs
committed in [`demo/runs/`](demo/runs/) — no API key, no install, no live calls — but
re-runs the validator fresh, so it shows the current verdict rather than a stored one.

It opens on a run where the checking **catches a real misattribution**: the model looked
up `katG` as context while annotating `furA`, then reported katG's pathways as furA's. A
demo where everything passes would prove nothing.

Deploying it to a Hugging Face Space: [`app/README_SPACE.md`](app/README_SPACE.md).

## Annotation pipeline

```bash
gar single katG                  # annotate one gene's function
gar epistasis katG furA ahpC     # look for mechanistic links between genes
gar eval                         # score the pipeline against a gold set
```

**The pipeline is a client of its own MCP server.** It spawns the server over stdio,
takes its tool schemas from `list_tools()`, and invokes tools with `call_tool()` — so
the schemas have exactly one definition and cannot drift from what an external MCP
client sees. (`gar --direct` dispatches in-process instead, for debugging.)

The pipeline fetches, computes and validates; the model chooses what to look up and
interprets what came back. It never writes to the store and never does arithmetic over
record IDs. In epistasis mode every pairwise relationship is computed *before* the model
sees it, with each shared pathway's size attached — `mtu01100` holds 698 of ~4,000
*M. tuberculosis* genes, so sharing it is a base rate, not a link.

Each run writes an append-only JSONL store: every tool call and its full result, every
deterministic computation, and every turn of the loop.

```
turn 1: stop=tool_use   tools=[kegg_pathways(katG), kegg_pathways(ahpC)]
turn 2: stop=tool_use   tools=[string_partners(katG), string_partners(ahpC)]
turn 3: stop=tool_use   tools=[pubmed_abstracts("ahpC katG")]
turn 4: stop=end_turn
```

### Citation validation

Every identifier in the summary is checked against what the tools actually returned,
recorded before the model saw it. Three failure classes:

- **unsupported** — an authoritative-looking ID that no tool returned.
- **cross-target** — an ID that *was* retrieved, but for a different gene than the
  sentence attributes it to. Invisible to a global membership check.
- **quote not in source** — a claim whose quoted span is absent from the retrieved
  text. Catches a real PMID carrying a fabricated finding.

```
NOT_IN_SOURCE  PMID:10609885  [likely_fabricated, similarity 0.38]
               quoted:  'KatG binds directly to AhpC in a stable complex'
               closest: 'katg is a catalase-peroxidase required for isoniazid activation'
```

Set membership and string containment, deliberately not similarity scoring: a citation
either names a record a tool returned or it does not. See
[docs/DESIGN.md](docs/DESIGN.md) for why failures are ranked but never adjudicated.

### Corpus manifest

Runs that find papers emit `<run>.corpus.jsonl` for a downstream full-text pipeline.
`in_pmc` matters more than the DOI — a DOI resolves to a usually-paywalled publisher,
a PMCID is the licit route to full text. `mentions` records the genes actually present
in the retrieved text, not the genes queried.

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

Retrieval fidelity   recall 1.0   precision 1.0        n=6 positives
Abstention           1.0                                 n=2 negatives
Citation integrity   citation precision 1.0   quote precision 1.0
```

*These figures are from a partial run: 8 of the 12 gold-set genes completed before
the run hit an API quota. The four that failed are recorded as errors and excluded
from the metrics rather than scored as misses. Run `gar eval` to reproduce over the
full set.*

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

## Running the MCP server

```bash
kegg-string-mcp                  # stdio transport
docker run -i --rm ghcr-or-local-image   # or in a container
```

The container image is built and handshake-tested on every push; see
[.github/workflows/ci.yml](.github/workflows/ci.yml).

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

246 tests. **The suite never touches the network** — `tests/conftest.py` swaps in a fake
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

MCP server, annotation pipeline and evaluation all working on `main`.
Ongoing work happens on `develop`.

Design rationale: [docs/DESIGN.md](docs/DESIGN.md).
