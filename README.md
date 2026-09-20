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
where the checking fires. No account, no install, nothing to run.

Defaults target *Mycobacterium tuberculosis* H37Rv (KEGG `mtu`, NCBI taxon `83332`),
but every tool takes the organism as a parameter and works for any species the
underlying service covers.

> **Research use only.** This is a research tool. It is not a clinical decision
> support system, is not validated for diagnostic, prognostic or treatment decisions,
> and must not be used to guide patient care. It annotates genes including the drug
> resistance loci `katG`, `inhA`, `rpoB`, `pncA` and `embB`; resistance interpretation
> for clinical purposes requires validated methods and expert review. Outputs are
> generated in part by a language model and require verification against the cited
> primary sources. Code provided under the MIT licence, without warranty — see
> [LICENSE](LICENSE); committed third-party data keeps its own terms — see [Licence](#licence).

## Tools

| Tool | Parameters | Returns |
|---|---|---|
| `kegg_pathways` | `gene` (KEGG ID, locus tag, or symbol)<br>`organism` = `mtu` | One record per KEGG pathway: `record_id` (e.g. `mtu00360`), name, `https://www.kegg.jp/entry/…` URL |
| `string_partners` | `gene`<br>`species` = `83332`<br>`limit` = `20`<br>`required_score` = `700` | One record per interaction partner: `record_id` (e.g. `83332.Rv1909c`), preferred name, combined score, full per-channel breakdown, `https://string-db.org/network/…` URL |
| `pubmed_abstracts` | `gene`<br>`organism` = `Mycobacterium tuberculosis`<br>`limit` = `10` | One record per article: `record_id` (PMID, e.g. `35038342`), title, abstract, `quotable_text`, journal, year, DOI, `https://pubmed.ncbi.nlm.nih.gov/…` URL |
| `uniprot_protein` | `gene`<br>`organism_id` = `83332`<br>`limit` = `3` | One record per UniProt entry: `record_id` (accession, e.g. `P9WG47`), protein name, function statements tiered by evidence code with supporting PMIDs, catalytic activity, PDB cross-refs, `quotable_text`, `https://www.uniprot.org/uniprotkb/…` URL |
| `lineage_markers` | `gene`<br>`organism` = `mtu` | One record per lineage-defining SNP the gene contains: `record_id` (e.g. `tbdb:851797`), the lineage it marks, H37Rv position, allele, `https://github.com/jodyphelan/tbdb` URL |
| `corpus_search` | `query`<br>`limit` = `5` | Abstracts from a locally built corpus, ranked by BM25 fused with dense embeddings. One record per paper: `record_id` (PMID), `quotable_text` (the whole abstract), `matched_passage` (the span that ranked), score and rank. Inert until `KEGG_STRING_MCP_CORPUS` points at a corpus from `scripts/build_corpus.py` |
| `resistance_variants` | `gene`<br>`drug` (optional) | Whether the gene is resistance-associated, the drugs, and per-grade counts; one record per resistance-associated variant: `record_id` (e.g. `tbdb:katG:p.Ser315Thr`), WHO grading, drug, source |

All seven are annotated `readOnlyHint`, expose structured output schemas, and are
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

It opens on a run where the checking **fires**: annotating `furA`, the model looked up
its neighbour `katG` for context, and the five KEGG pathways it reported for katG are
flagged because they were not retrieved for furA. Read the write-up and the flag is
conservative rather than a catch — the prose attributes those pathways to katG, which is
correct. The check compares an identifier against the gene being annotated and does not
read the sentence that cites it, so context citations surface here. That is the intended
direction of the error: a check that adjudicated intent would be the thing to distrust.
Runs where the checking fires are marked in the picker, because a demo where nothing
fires shows nothing.

**No committed run currently contains a real fabrication.** That furA run is an earlier
capture; every run re-captured after the pipeline gained identity resolution and
`corpus_search` passes, furA included, so the re-run sits beside the original as
`furA-rerun`.

The page runs entirely in the browser on **bare Pyodide, with no packages installed** —
there is no server and nothing to resolve at load. That is only possible because the
replay layer is standard-library only, which a test enforces. (An earlier attempt used
Gradio-Lite, which pulls gradio through micropip and huggingface-hub behind it, and
failed to start; `gradio` is still what runs the demo *locally*.) `docs/index.html` is generated by
[`demo/build_pages.py`](demo/build_pages.py) on deploy and never committed: a stale copy
could show a verdict the validator no longer produces.

Running it locally, or on a Hugging Face Space instead:

```bash
pip install -e ".[demo]"
python -m app.app                       # local
```

Space deployment: [`app/README_SPACE.md`](app/README_SPACE.md).

## Annotation pipeline

`gar` is the pipeline's command. `pip install -e .` puts it on your PATH alongside
`kegg-string-mcp` (the server); it runs `main()` in
[`src/kegg_string_mcp/cli.py`](src/kegg_string_mcp/cli.py), which spawns the server as a
subprocess and drives the agent loop against it.

```bash
gar single katG                  # annotate one gene's function
gar epistasis katG furA ahpC     # look for mechanistic links between genes
gar eval                         # score the pipeline against a gold set
gar table                        # turn every finished run into TSVs
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

Both modes settle the gene's names first, through `identity.py`, so every tool in a run
answers about the same gene and the aliases it accepted — and the ones it refused as
ambiguous — are on the record rather than re-derived per caller.

**Epistasis also computes why a pair might co-occur without being related**, and puts
the answer in the verdict the model is forbidden to contradict. Two genes graded
resistance-associated for the same drug are co-selected by treating with it; two genes
carrying positions that define the same lineage are inherited together by descent.
Neither is a link between the genes.

```
katG / inhA   Direct STRING interaction (combined 0.916), supported essentially only by
              literature co-mention. CONFOUND: both genes carry variants graded
              resistance-associated for isoniazid. Treating with that drug selects
              both, so they co-occur across isolates under co-selection rather than
              through any link between the genes.
```

Any scan over clinical isolates will flag katG/inhA, and the explanation is the
treatment. The lineage and catalogue lookups behind that are made by the pipeline, not
requested of the model: they are the only claims in the run that a prompt used to ask
for and nothing checked.

Each run writes an append-only JSONL store: every tool call and its full result, every
deterministic computation, and every turn of the loop.

```
turn 1: stop=tool_use   tools=[kegg_pathways(katG), kegg_pathways(ahpC)]
turn 2: stop=tool_use   tools=[string_partners(katG), string_partners(ahpC)]
turn 3: stop=tool_use   tools=[corpus_search("ahpC promoter mutations compensate katG")]
turn 4: stop=end_turn
```

### The rules, as a skill

The prompts are not in the code. [`skills/gene-annotation/SKILL.md`](skills/gene-annotation/SKILL.md)
holds them — the abstention rules, the citation rules, and the two modes — and
`agent/modes.py` reads the marked blocks out of it at import and assembles them.

It is also a loadable skill, so an agent driving these tools directly follows the same
text the scored pipeline does rather than a second copy that drifts. That is the same
argument as taking the tool schemas from `list_tools()`: one definition, or the thing
being measured and the thing being used are not the same thing.
[`reference/tools.md`](skills/gene-annotation/reference/tools.md) carries the per-tool
caveat that governs each, and
[`reference/envelope.md`](skills/gene-annotation/reference/envelope.md) the result
envelope and the four meanings of an empty `records`.

### Citation validation

Every identifier in the summary is checked against what the tools returned,
recorded before the model saw it. Three failure classes:

- **unsupported** — an authoritative-looking ID that no tool returned.
- **cross-target** — an ID that *was* retrieved, but not for the gene being annotated.
  Invisible to a global membership check. It compares the identifier against the run's
  own record of which gene returned it; it does not parse the sentence, so a neighbour
  gene cited as context is flagged too — conservative in the direction that keeps the
  check mechanical. Epistasis runs have no single target and so skip this class.
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
in the retrieved text, not the genes queried: `genes_named` for a corpus record, whose
`mentions` is the corpus-build query, and `mentions` for a live PubMed one.

### Tabular output

`gar table` reads the stores already on disk and writes four TSVs. It makes no model
call and no request, so the fan-out over many genes can be a rule per gene and one
aggregate rule here.

| File | One row per | Carries |
|---|---|---|
| `genes.tsv` | annotated gene | resolved ID, accepted and refused aliases, `matched_by` per source, per-tool counts, the resistance flag, and the validation verdict |
| `pairs.tsv` | pair in an epistasis run | the deterministic verdict, the per-channel STRING scores, shared pathways split by specificity, and the confound tokens |
| `resistance_variants.tsv` | catalogued associated variant | mutation, drug, WHO grade — the specific variants, not a count of them |
| `lineage_markers.tsv` | lineage-defining position in the locus | H37Rv position, lineage, allele |

**Built from the stores, never from the summary.** Deriving the numbers from prose
would reintroduce exactly what the citation validator exists to catch, and at a hundred
rows nobody reads the prose to notice. So the cells keep apart what prose renders
identically:

```
kegg_n = 0          KEGG was asked and assigns nothing — true of 71% of genes
kegg_n = NA         KEGG never resolved the gene
string_truncated    the partner list hit `limit`; an absent partner is a query artefact
resistance NA       absent from the WHO catalogue: never assessed, which is not negative
```

Whether an empty result is a zero or a missing value is not uniform across the tools, so
it is declared per tool rather than inferred. `lineage_markers` reports a gene with no
marker as `matched_by: "none"` — its index holds only genes that have one — and reading
that as unresolved turns an informative negative into a blank.

Counts come only from calls whose gene argument matches the target's alias set.
Annotating `furA` and looking up `katG` for context is legitimate; letting katG's
pathways land in furA's row is the cross-target error in tabular form, and unlike the
citation flag nothing downstream would show it.

## Rule annotation

An upstream rule learner classifies isolates as resistant or susceptible from
per-locus variant states, and emits rules:

```
Rv1908c=1 AND Rv2428=1  ->  R     numerosity 11   precision 0.9524
```

This layer says what the structured sources already account for, before any model
is involved. `=1` means the locus carries a qualifying variant; `=0` asserts it
matches the reference, which is a claim about every isolate the rule covers
rather than silence about the gene.

```bash
gar features LABELS.txt --annotation h37rv.bed    # do the feature labels resolve?
gar rules RULES.tsv     --annotation h37rv.bed    # annotate loci, classify rules
gar gold                --annotation h37rv.bed    # score the classifier on known rules
```

**The reference annotation is an input, never fetched.** Feature names come from
whichever annotation the variant caller used — `snpEff genes2bed
Mycobacterium_tuberculosis_h37rv`, or NCBI's H37Rv GFF3 where SnpEff is not
installed — so resolving them against a different gene list is how a locus
silently becomes the wrong locus. BED and GFF/GTF are both read, the format from
the shape of a row rather than the extension, and a GFF3's child rows (`CDS`,
`rRNA`, `exon`) are folded into their parent locus rather than counted as
separate genes.

Non-coding features are named differently by different producers. A pipeline
collapsing SnpEff `ANN` fields emits `upstream_<gene>`, `downstream_<gene>` and
`intergenic_<a>-<b>`, with a transcript version on the gene part
(`upstream_Rv1482c.1`); another names the flanking pair directly. A trailing
`.N` is tried as part of the name before it is tried as a version, because H37Rv
has ten symbols — `TB7.3`, `TB15.3`, `TB31.7` and the rest of that family —
where the number *is* the gene. Both resolve,
the prefix set is the caller's to override — a
refusal names the prefixes it tried, so a new vocabulary is a configuration
change rather than a silent loss.

Two conditions in one rule can be driven by the same variant, and a rule saying
so twice reads as epistasis. Three cases are reported, and they differ in
certainty: both conditions naming one locus, resolved spans that intersect
(H37Rv has 917 overlapping consecutive gene pairs — mostly 4 bp start/stop
junctions, but 56 of at least 50 bp), and a region feature paired with the gene
that bounds it. The last is **undetermined** rather than either answer: a caller
reporting upstream variants within a window annotates a variant inside the
flanking gene to both features, and whether any did is in the caller's distance
field. The check keys on adjacency *in the rule*, so it fires on the pairing that
can alias rather than on every gene with a near neighbour.

`upstream_X` resolves to the intergenic interval 5' of X *on X's own strand*,
and says in its note that the producer's window is usually wider: SnpEff's
default upstream window is 5,000 bp, and 3,048 of H37Rv's 3,049 intergenic
intervals are narrower than that (median 81 bp). Where X abuts or overlaps its
5' neighbour — 830 of 4,008 genes — there is no interval to resolve to, and the
answer says so rather than reading as a failed lookup. Measured against one real
2,903-label vocabulary, KEGG's curated list was missing 76 of them and disagreed
with 1,346 more on the strand suffix. KEGG stays authoritative for pathways; the
caller's annotation is authoritative for what a locus is and where it sits. Its
`sha256` is recorded, so a reference update cannot change the results invisibly.

### What a rule is classified as

Each condition gets a role from its state and what the WHO catalogue holds:

| | catalogue anchor | assessed, none associated | absent |
|---|---|---|---|
| **state 1** | `anchor` | `compensator_candidate` | `unknown` |
| **state 0** | `negated_anchor` | `negated` | `negated` |

That middle cell is what the layer is for. `rpoC` (1,401 catalogued variants, 0
graded associated), `rpoA` (507, 0) and `ahpC` (250, 0) — the canonical
compensatory loci — all carry it: in the catalogue because they recur in
resistant isolates, graded as not conferring resistance themselves. 50 of the
catalogue's 74 genes look like that.

Rules then carry every signature that applies, with precedence deciding which
leads:

| Signature | Means |
|---|---|
| `confounded:feature_overlap` | two conditions one variant could satisfy — the same observation entered twice, which reads as epistasis |
| `confounded:lineage` | every present locus marks one lineage — carried together by descent |
| `known:compensation` | anchor + a candidate assessed against the **same drug**, ideally functionally linked |
| `known:alt_route` | resistance while a canonical locus is at reference |
| `discordant` | a known resistance locus varying in isolates predicted susceptible |
| `unknown` | a locus the catalogue never assessed |
| `confounded:multidrug` | loci for **different** drugs — MDR is defined that way, so this is the diagnosis, not a relationship |
| `confounded:co_selection` | loci for the **same** drug — one drug selecting both |
| `known:resistance` | recovers known pharmacology |

Drug concordance is load-bearing. `rpoA` has the compensatory shape and STRING
supplies a co-mention edge for it, but it was only ever assessed against
rifampicin — so pairing it with `katG` is not compensation, and the verdict says
why. `gar gold` scores twelve rules of known answer, and that negative control is
the one nothing else in the suite guards.

### Output

| File | One row per | Carries |
|---|---|---|
| `loci.tsv` | distinct locus | product, pathways, partners, lineage markers, catalogue status, how many rules it appears in, and `linked_loci` |
| `rules.tsv` | rule | conditions, roles, signatures, the learner's own statistics, and the verdict |
| `questions.tsv` | gap the sources could not close | what to ask the literature, and which rules raised it |

Relationships live in a structured string (`ahpC:string_textmining_only:0.968|furA:adjacent:6bp`)
rather than a third file, because `loci.tsv` is meant to be read as a
supplementary table. A STRING edge supported only by textmining is co-mention in
papers, so the channel is carried rather than flattened into one score.

`questions.tsv` is stage A of the literature layer and costs nothing: the
questions fall out of the classification, so a run can be counted before any
retrieval happens. Each asks for a verbatim passage rather than for a verdict.

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

The tools above answer what a curated source records. KEGG assigns a pathway to 29% of
*M. tuberculosis* genes; UniProt describes many of the rest only by similarity to
something else. For those genes the answer exists, but it is in papers.

So there is a second path, and it is retrieval-augmented generation: gather the
literature first, then let the model read it, quote it, and be checked against it. Three
choices make that more than a phrase.

**What is gathered, and by whom.** `build_corpus.py` decides which genes need papers at
all — those whose structured annotation is thin — then asks PubMed for the 20 most
relevant records per gene, deduplicates by PMID, and splits each abstract into 180-word
passages. The corpus is a fixed, hashed artefact: the same file gives the same answers
tomorrow, and the run store records which file was read.

**How a passage is found.** Keyword search knows that `Rv1908c` is a rare, decisive token
and ignores a paraphrase; embeddings know that "peroxide detoxification" and "antioxidant
defence" are the same idea and have no feel for an identifier. Neither is sufficient for
literature where gene symbols and prose descriptions carry the meaning together, so the
tool runs both and merges the two ranked lists by **reciprocal rank fusion**: each paper
scores `1 / (60 + rank)` in each list, and the scores are added. Fusing on rank rather
than score is the point — a BM25 score of 27.96 and a cosine of 0.655 are not on a common
scale, and normalising them would invent a calibration neither arm has. Agreement between
two different notions of relevance outranks one arm's favourite.

**When the model reaches for it.** Never by default. Literature is the expensive, noisy
channel, so the prompt asks for it only when the structured tools leave a real question
open, and asks for `corpus_search` before live PubMed because its ranking is measured
where the corpus applies. An empty corpus result says the corpus does not cover the query,
never that no literature exists — so the agent falls back rather than concluding.

What comes back is quotable, not just readable: every claim drawn from a paper must carry
a verbatim span, and [citation validation](#citation-validation) checks that span against
the retrieved text. Retrieval supplies the evidence; the checking decides whether the
write-up actually used it.

Three arms, measured head to head on the same corpus and the same queries:

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

The candidate set bounds all of it. A corpus built from the 20 most relevant records per
gene is a gene-centric sample of the literature, so these numbers measure ranking *within*
that sample and not recall over PubMed — and a paper about a pair that ranked 21st for
both of its genes is beyond every arm's reach.

```bash
python scripts/build_corpus.py --extended --all-genes --tag tb41   # build the corpus
python scripts/run_comparison.py data/corpus_tb41.json --tag tb41  # measure the arms
python scripts/residue.py --tag tb41                               # what stays unexplained
```

Corpora are gitignored: a few hundred PubMed abstracts under publisher copyright is bulk
redistribution. Every command above rebuilds from the tools in this repo.

**One gene, one identity.** `Rv0678` and `mmpR5` are the same protein, and every stage
resolves that once, in `kegg_string_mcp/identity.py`, rather than matching whatever string
it was handed. Before it, the STRING edge between Rv0678 and Rv0676c scored 0.989 under
their symbols and read as silent under their locus tags. The same module supplies the
alias set the corpus matches on, and records which sources failed to resolve a gene — so a
failed lookup is reported as **undetermined** rather than counted as an absence of
evidence, which would make the pair look novel. See [docs/TODO.md](docs/TODO.md).

**Retrieval is routed, not run on everything.** Literature is the expensive, noisy
channel, so `build_corpus.py` runs it only on genes whose structured annotation is thin —
no UniProt function, only inferred function, or no KEGG pathway — and writes the routing
decision beside the corpus so a reader can see why an annotation rested on papers. On the
41-gene set that is 22 of 41. `--all-genes` bypasses it, which the arm comparison needs
because a head-to-head on gene *pairs* wants every gene regardless of coverage.

## Scripts and entry points

Everything runnable in the repo, and what it is for. The `gar` and `kegg-string-mcp`
commands are installed by `pip install -e .`; the rest are run with `python`.

| What | Command | Does |
|---|---|---|
| **MCP server** | `kegg-string-mcp` | Serves the seven tools over stdio. `src/kegg_string_mcp/server.py`. |
| **Annotation pipeline** | `gar single \| epistasis \| eval` | Runs the agent loop against the server and validates the result. Installed by `pip install -e .`. `src/kegg_string_mcp/cli.py`. |
| **Rule annotation** | `gar rules RULES.tsv --annotation h37rv.bed` | Resolves the rule set's feature labels against the supplied annotation, annotates each distinct locus once, and classifies every rule against the WHO catalogue. Writes `loci.tsv`, `rules.tsv`, `questions.tsv`. No model. `src/kegg_string_mcp/rules/`. |
| **Feature pre-flight** | `gar features LABELS.txt --annotation h37rv.bed` | Reports how much of a feature vocabulary resolves, and why the rest does not. Free, offline. |
| **Classifier gold set** | `gar gold --annotation h37rv.bed` | Scores the rule classifier against twelve rules of known answer, including a negative control. |
| **Run tables** | `gar table --runs runs --out tables` | Turns finished run stores into `genes.tsv`, `pairs.tsv`, `resistance_variants.tsv` and `lineage_markers.tsv`. No model, no network, re-runnable. `src/kegg_string_mcp/agent/tables.py`. |
| **Corpus build** | `python scripts/build_corpus.py --extended --all-genes --tag tb41` | Routes genes by annotation coverage, resolves every alias, fetches abstracts, dedupes by PMID, chunks to 180 words, records which genes each passage names. Writes `data/corpus_<tag>.json` and `coverage_<tag>.json`. |
| **Arm comparison** | `python scripts/run_comparison.py data/corpus_tb41.json --tag tb41` | Measures lexical, dense and hybrid retrieval over all gene pairs, then again over the pairs STRING has no edge for, which removes the circularity in scoring relevance by gene names. Writes `comparison_<tag>.json`. |
| **Unexplained residue** | `python scripts/residue.py --tag tb41` | Reads what the earlier steps wrote, adds KEGG pathway membership, and reports which pairs no source accounts for — the input to hypothesis generation. |
| **Demo capture** | `python demo/build.py NAME=path/to/run.jsonl` | Turns a raw run store into the compact record the demo replays. Deliberately drops the validation verdict, so the demo recomputes it. |
| **Demo page build** | `python demo/build_pages.py` | Generates the serverless `docs/index.html` for GitHub Pages: bare Pyodide, module bodies copied verbatim so the page cannot show a verdict the library does not produce. Run by the Pages workflow, never committed. |
| **Demo, locally** | `python -m app.app` | The same replay through gradio, for running on your own machine or a Hugging Face Space. |
| **Container handshake** | `python .github/scripts/mcp_smoke.py <image>` | Drives a real MCP stdio handshake against the built image and asserts the tool list. Runs in CI on every push. |

## Install

```bash
conda create -n kegg-string-mcp python=3.11 -y && conda activate kegg-string-mcp
pip install -e ".[dev]"
pytest
```

That puts two commands on your PATH: `kegg-string-mcp`, the MCP server, and `gar`, the
annotation pipeline.

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

465 tests. **The suite never touches the network** — `tests/conftest.py` swaps in a fake
HTTP client that replays saved responses from `tests/fixtures/`, so the suite runs in
seconds, gives the same answer every time, works offline and in CI, and does
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

384 of the 465 run with no optional extras installed at all, which is what CI checks on
every supported Python: importing the library must not require the vector stack. The
remaining tests belong to the `vector` and `demo` extras and run in their own CI jobs.

## Licence

The source code is MIT licensed — see [LICENSE](LICENSE). The MIT licence does not cover
the third-party data committed in `demo/runs/` and `tests/fixtures/`, which stays under
its owners' terms, listed below and in [NOTICE](NOTICE).

### Upstream data

This server queries public sources on the caller's behalf and caches responses locally.
It redistributes no bulk data; the exceptions are small, committed so the demo can replay
and the tests can run offline, and named below. Terms for anything the server fetches
remain the caller's responsibility:

- **KEGG** — the KEGG API is provided only for academic use by academic users. Academics
  who use KEGG to provide a service need an academic service provider licence; non-academic
  use needs a commercial licence from Pathway Solutions
  ([terms](https://www.kegg.jp/kegg/legal.html)). Committed: KEGG gene and pathway lists in
  `tests/fixtures/kegg_*.tsv` (about 550 lines) and the pathway records in `demo/runs/`.
- **STRING** — [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). Committed: API
  responses in `tests/fixtures/string_*.json` and the interaction records in `demo/runs/`.
  Cite: Szklarczyk D *et al.* The STRING database in 2023. *Nucleic Acids Res*
  51(D1):D638–D646 (2023). [doi:10.1093/nar/gkac1000](https://doi.org/10.1093/nar/gkac1000)
- **UniProt** — [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). Committed: API
  responses in `tests/fixtures/uniprot_*.json`. Cite: The UniProt Consortium. UniProt: the
  Universal Protein Knowledgebase in 2025. *Nucleic Acids Res* 53(D1):D609–D617 (2025).
  [doi:10.1093/nar/gkae1010](https://doi.org/10.1093/nar/gkae1010)
- **TB-Profiler / tbdb** ([jodyphelan/tbdb](https://github.com/jodyphelan/tbdb), LGPL-3.0)
  — source of the lineage barcode (`barcode.bed`, after Coll 2014 and Napier 2020) and the
  resistance catalogue (`mutations.csv`). Neither file is committed; both are fetched and
  cached at run time. Committed: 4 lineage-marker rows in `demo/runs/`.
- **WHO catalogue of mutations** — the resistance catalogue in tbdb is derived from WHO,
  *Catalogue of mutations in* Mycobacterium tuberculosis *complex and their association
  with drug resistance*, 2nd ed. (2023), licensed
  [CC BY-NC-SA 3.0 IGO](https://creativecommons.org/licenses/by-nc-sa/3.0/igo/):
  non-commercial use only, and adaptations must be shared under the same licence.
  Committed: 662 distinct resistance-variant rows across `demo/runs/` (for katG, rpoB,
  pncA, gyrA and gyrB), each a gene, position or variant with its WHO grading — about 1%
  of the catalogue's 49,330. These rows are not MIT licensed. Anyone using the
  catalogue should take it from the WHO publication directly.
- **PubMed** — citation metadata comes from NLM; abstract text remains under the copyright
  of the respective publishers. The server retrieves abstracts per query and caches them
  locally for the caller. Committed: 58 abstracts, in `demo/runs/`,
  `tests/fixtures/corpus_small.json` and `tests/fixtures/pubmed_efetch_*.xml`, because
  the demo cannot show quote-checking without the text it checks against. They are
  included for non-commercial research and teaching only. Larger corpora built by the
  retrieval arm are **not** committed. Respect NCBI's
  [E-utilities usage policy](https://www.ncbi.nlm.nih.gov/books/NBK25497/) — set
  `NCBI_EMAIL`.

## Status

MCP server, annotation pipeline, evaluation, the run tables and the retrieval arm are
all on `main`. Rule annotation is on `rule-annotation`: the deterministic half is
complete and scored; the literature layer is stage A only.

- The rules an external client should follow: [skills/gene-annotation/SKILL.md](skills/gene-annotation/SKILL.md)
- Design rationale: [docs/DESIGN.md](docs/DESIGN.md)
- Retrieval comparison, and what it does not show: [docs/RETRIEVAL.md](docs/RETRIEVAL.md)
- Found, understood, not fixed: [docs/KNOWN_ISSUES.md](docs/KNOWN_ISSUES.md)
