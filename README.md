# kegg-string-mcp

[![LLM Gene Annotation Demo](https://img.shields.io/badge/demo-LLM%20Gene%20Annotation-2ea44f)](https://yaaoppong.github.io/kegg-string-mcp/)
[![CI](https://github.com/YaaOppong/kegg-string-mcp/actions/workflows/ci.yml/badge.svg)](https://github.com/YaaOppong/kegg-string-mcp/actions/workflows/ci.yml)

An [MCP](https://modelcontextprotocol.io) server exposing **KEGG**, **STRING**, **UniProt**,
**PubMed**, the **M. tuberculosis lineage barcode** and the **WHO catalogue of resistance
mutations** as model-callable tools for gene annotation. Every record returned carries a
stable ID and a resolvable URL, so a downstream agent's citations can be checked
programmatically against what was retrieved.

Alongside the tools, a **literature retrieval arm** — BM25, dense embeddings and their
fusion — measured head to head on the same corpus and queries, with the numbers and their
limits in [docs/RETRIEVAL.md](docs/RETRIEVAL.md).

**▶ [Try it in your browser](https://yaaoppong.github.io/kegg-string-mcp/)** — pick a
gene, watch the model call tools, then watch its citations get checked. It opens on a run
where the checking catches a bad citation. No account, no install, nothing to run.

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
| `lineage_markers` | `gene`<br>`organism` = `mtu` | One record per lineage-defining SNP the gene contains: `record_id` (e.g. `tbdb:851797`), the lineage it marks, H37Rv position, allele, `https://github.com/jodyphelan/tbdb` URL |
| `resistance_variants` | `gene`<br>`drug` (optional) | Whether the gene is resistance-associated, the drugs, and per-grade counts; one record per resistance-associated variant: `record_id` (e.g. `tbdb:katG:p.Ser315Thr`), WHO grading, drug, source |

All six are annotated `readOnlyHint`, expose structured output schemas, and are
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

## LLM Gene Annotation Demo

**<https://yaaoppong.github.io/kegg-string-mcp/>**

Pick a gene, watch the tool calls in sequence, read the write-up, then watch every
citation in it get checked. It replays runs committed in [`demo/runs/`](demo/runs/) — no
API key, no install, no live calls — but re-runs the validator fresh each time, so it
shows the current verdict rather than a stored screenshot of one.

It opens on a run where the checking **catches a real misattribution**: the model looked
up `katG` as context while annotating `furA`, then reported katG's pathways as furA's.
Runs where checking fails are marked in the picker, because a demo where everything
passes proves nothing.

The page runs entirely in the browser via [Gradio-Lite](https://www.gradio.app/guides/gradio-lite)
and Pyodide — there is no server. That is only possible because the replay layer is
standard-library only, which a test enforces. `docs/index.html` is generated by
[`demo/build_pages.py`](demo/build_pages.py) on deploy and never committed: a stale copy
could show a verdict the validator no longer produces.

Running it locally, or on a Hugging Face Space instead:

```bash
pip install -e ".[demo]"
python -m app.app                       # local
```

Space deployment: [`app/README_SPACE.md`](app/README_SPACE.md).

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

Every identifier in the summary is checked against what the tools returned,
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
a PMCID is the licit route to full text. `mentions` records the genes present
in the retrieved text, not the genes queried.

## Evaluation

**The reference is incomplete, and that is the point.** KEGG assigns a pathway to just
**about 1,170 of 4,008** *M. tuberculosis* genes — 29%. `gyrA`, one of the most studied genes
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

## Literature retrieval

The tools above answer what a curated source records. For a gene they do not describe,
the answer has to come from papers — so there is a second retrieval path, with three arms
measured against each other on the same corpus and the same queries.

| Arm | Mechanism | precision@10 |
|---|---|---|
| `hybrid` | reciprocal rank fusion of the two below | **0.917** |
| `lexical` | BM25 over tokens | 0.844 |
| `dense` | MiniLM embeddings in Chroma, cosine | 0.739 |

Measured over 41 genes, 676 papers, 1,350 chunked passages, 820 gene-pair queries.
Relevance is scored deterministically — does the retrieved passage name the gene asked
about — so it is reproducible and cannot be tuned after the fact. The ranking holds when
the query set is restricted to pairs STRING has no edge for, which removes the
circularity in scoring relevance by gene names.

The dense arm is the least precise and retrieves the most different papers (Jaccard 0.158
against lexical), which is the case for keeping it and the case against replacing BM25
with it. [docs/RETRIEVAL.md](docs/RETRIEVAL.md) has the full result and, more usefully,
what it does not show.

```bash
python scripts/build_corpus.py --extended --all-genes --tag tb41   # build the corpus
python scripts/run_comparison.py data/corpus_tb41.json --tag tb41  # measure the arms
python scripts/residue.py --tag tb41                               # what stays unexplained
```

Corpora are gitignored: a few hundred PubMed abstracts under publisher copyright is bulk
redistribution. Every command above rebuilds from the tools in this repo.

**Retrieval is routed, not run on everything.** Literature is the expensive, noisy
channel, so `build_corpus.py` runs it only on genes whose structured annotation is thin —
no UniProt function, only inferred function, or no KEGG pathway — and writes the routing
decision beside the corpus so a reader can see why an annotation rested on papers. On the
41-gene set that is 22 of 41. `--all-genes` bypasses it, which the arm comparison needs
because a head-to-head on gene *pairs* wants every gene regardless of coverage.

## Install

```bash
conda create -n kegg-string-mcp python=3.11 -y && conda activate kegg-string-mcp
pip install -e ".[dev]"
pytest
```

The server and the annotation pipeline need nothing beyond the base install. The
retrieval arm is a separate extra, kept out of the default so the MCP server does not
carry a vector store it never uses:

```bash
pip install -e ".[dev,vector]"     # chromadb, rank-bm25, langgraph
pip install -e ".[dev,demo]"       # gradio, for running the demo locally
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

396 tests. **The suite never touches the network** — `tests/conftest.py` swaps in a fake
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
pytest -q                                                        # everything
pytest -q --ignore=tests/test_retrieval.py --ignore=tests/test_app.py   # base install only
```

323 of the 396 run with no optional extras installed at all, which is what CI checks on
every supported Python: importing the library must not require the vector stack. The
remaining tests belong to the `vector` and `demo` extras and run in their own CI jobs.

## Licence

This project is MIT licensed — see [LICENSE](LICENSE).

### Upstream data

This server queries public sources on the caller's behalf and caches responses locally.
It redistributes no bulk data; the exceptions are small and named below. Terms remain the
caller's responsibility:

- **KEGG** — free for academic use; commercial use requires a licence from Pathway Solutions.
- **STRING** — CC BY 4.0, free for academic and commercial use, attribution required.
- **UniProt** — CC BY 4.0, attribution required.
- **TB-Profiler / tbdb** ([jodyphelan/tbdb](https://github.com/jodyphelan/tbdb), LGPL-3.0)
  — source of the lineage barcode (`barcode.bed`, after Coll 2014 and Napier 2020) and the
  resistance catalogue (`mutations.csv`, derived from the WHO catalogue of mutations in
  *M. tuberculosis complex*). Neither file is committed; both are fetched and cached at
  run time. A small number of *derived records* are committed in `demo/runs/` so the demo
  can replay offline — 139 resistance-variant rows and 4 lineage-marker rows, each a
  gene, position or variant with its WHO grading. Anyone using the underlying catalogues
  should take them from tbdb and the WHO publication directly and check their terms.
- **PubMed** — records are US government works and free to use; the abstracts themselves
  are frequently under publisher copyright. The server retrieves them per query and caches
  locally for the caller. A small number of abstracts *are* committed, in `demo/runs/` and
  `tests/fixtures/`, because the demo cannot show quote-checking without the text it checks
  against and the tests cannot run offline without it — roughly 40 records, retained for
  research and educational use. Larger corpora built by the retrieval arm are **not**
  committed. Respect NCBI's
  [E-utilities usage policy](https://www.ncbi.nlm.nih.gov/books/NBK25497/) — set
  `NCBI_EMAIL`.

## Status

MCP server, annotation pipeline, evaluation and the retrieval arm are all on `main`.

- Design rationale: [docs/DESIGN.md](docs/DESIGN.md)
- Retrieval comparison, and what it does not show: [docs/RETRIEVAL.md](docs/RETRIEVAL.md)
- Found, understood, not fixed: [docs/KNOWN_ISSUES.md](docs/KNOWN_ISSUES.md)
