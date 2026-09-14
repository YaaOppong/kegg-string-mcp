# To do

Raised by the *M. tuberculosis* epistasis pipeline, which uses this package as
its annotation step (step 8). Findings dated 2026-09-14 against `f6c4012`,
organism H37Rv / KEGG `mtu` / taxon 83332. IDs are that report's.

The work is batched by what it changes: batch 1 changes **results**, batches 2
and 3 change the **tool**. Batch 1 is done; the rest is open.

---

## Batch 1 — done (identity and the third answer)

F1, F9, F10, F12 were four symptoms of two missing abstractions, and all four
erred the same way: making a known relationship look unexplained, which in an
epistasis run is the strongest possible claim — a novel candidate — produced by a
lookup that failed.

* **F9, F12 — no identity layer.** Gene identity was re-derived in five places
  (preferred-name keys, name-or-ID-suffix matching, a literal regex, a locus-tag
  retry, a private alias set), so the answer depended on the spelling the caller
  typed: `Rv0678`/`Rv0676c` read as silent while `mmpR5`/`mmpL5` read as 0.989.
  → `kegg_string_mcp/identity.py`, one resolution per gene, aliases filtered
  against KEGG's own ambiguity set, `Corpus.aliases` + `Passage.named_via`.
* **F1 — truncated partner lists answered as negatives.** Direct interaction was
  inferred from two ranked top-20 lists, so katG/rpoC (0.823, rank 24) and
  rpoB/rpsC (0.991, rank 23) were reported as having none.
  → `StringClient.network()` asks about the pair itself; the partner path now
  says `truncated` rather than claiming a negative.
* **Resolution tries harder before giving up.** Sources rescue each other: a
  symbol STRING refuses is retried with the locus tag UniProt or KEGG supplied,
  which is `coverage.py`'s cross-source retry generalised. The spelling that
  worked is recorded in `resolved_via`, never silently substituted. Unresolved
  now means every source refused every known spelling.
* **F10 — a resolution failure became negative evidence.**
  → `provenance.answered()` (promoted from `coverage.py`, which was the one
  module that had this right), `unresolved`/`truncated` verdicts, and a fourth
  bucket in the residue: pairs that are **undetermined** are reported in full and
  excluded from the residue and its denominator.

**Breaking, on purpose** — `retrieval.independence.gene_partner_map()` is gone.
Build a `StringEdges` with `network_edges(genes, string, identities)` (or
`partner_edges(...)` as the fallback) and pass that to `classify()`. The old
function could not be fixed in place: it discarded `resolved` before its caller
could see it, which is what made a failed lookup indistinguishable from silence.

---

## Speculative — a literature proposer for renamed genes

**No gene in the TB-41 set needs this.** Replaying the response cache, every real
gene resolves in STRING under every spelling tried, `mmpR5` and `Rv0678` both
returning `83332.Rv0678`; the only identifier that fails is `fakeGene1`, the
synthetic one from the F10 repro. F9 was never a resolution failure — it was a
lookup keyed on the preferred name. So do not build this until a real gene fails.

`identity.resolve(..., proposer=...)` is the seam, and it is unused: nothing
passes a proposer, so resolution stays deterministic and offline by default. The
proposer itself would belong in the agent layer (it needs a model and an API key,
and `identity.py` must stay importable by scripts that have neither).

The case it would answer is a gene renamed since the query set was written: no
structured source can map the old symbol, because the old symbol is exactly what
they dropped. A paper can — "Rv0678 (formerly …)" is a sentence that exists.

The contract is already enforced here: the proposer returns
`[(candidate_spelling, evidence)]`, and `identity` re-asks KEGG, UniProt and
STRING. **A name no structured source resolves is refused and recorded, and the
gene stays unresolved.** The model reads; the structured sources assert. A
hallucinated identifier therefore fails closed rather than propagating as a
resolved gene, and the evidence for an accepted one is kept in `proposals` next
to the confirming source.

If it is ever built: an agent-layer proposer using `pubmed_abstracts` (and
`corpus_search` where a corpus exists) on the unresolved symbol, with the existing
span validator requiring the quote to be verbatim in the retrieved abstract. Note
that `mmpR5`/`Rv0678` is a symbol-and-locus-tag pair, not a rename, so it is not a
test case for this; a genuine synonym pair such as `fabG1`/`mabA` is.

## Batch 2 — tool hygiene (low risk, independent of each other)

* **F2 — `gar --direct` advertises tools it cannot dispatch.** `TOOL_PARAMS` has
  4 tools; `server_tool_schemas()` advertises 6 (+`corpus_search`), so
  `lineage_markers` and `resistance_variants` return "unknown tool" under direct
  dispatch — and the epistasis prompt *requires* a lineage result.
  Fix: build the dispatch table from the server registry, and test that the
  dispatchable set equals the advertised set. `agent/pipeline.py`.
* **F3 — `resolved.accession` ignores review status.** `uniprot.py:268` takes
  `records[0]` in search order, so a TrEMBL entry can become "the" accession.
  `coverage.py` and now `identity.py` both break this tie on `reviewed`; the tool
  should too, and should expose `resolved.reviewed`.
* **F4 — no way to quiet the MCP server child.** One multi-line log entry per
  upstream request on the parent's stderr buries pipeline output over a run of
  comparisons. Fix: forward `KEGG_STRING_MCP_LOG_LEVEL` in `FORWARDED_ENV`, or
  accept an `errlog` argument. `agent/mcp_tools.py`.
* **F8 — `drug` values are undiscoverable.** The filter matches tbdb's lowercase
  full names exactly, so `INH` and catalogue-absent fluoroquinolones match
  nothing, silently. Fix: list the 18 values in the tool description, and note
  when `drug` matched nothing. `resistance.py`.
* **F11 — compensatory annotations are filtered out.** Only `Assoc w R` /
  `Assoc w R - Interim` rows are returned, but the catalogue keeps compensatory
  information as free text on `Uncertain significance` rows, so `rpoC` returns 0
  records. Step 8 currently reaches past the API into `ResistanceClient._load()`.
  Fix: a public `compensatory(gene)`; longer term, a curated citable table.

## Batch 3 — API surface and new capability

* **F5 — a stable public API.** Step 8 imports internals: `pathway_sizes()`,
  `gene_index()`, `agent.evidence.pair_evidence()`, a hand-built
  `PoliteClient(DiskCache())`, a re-implementation of `_annotated_gene_count()`,
  and `_load()`. Stage 2 is wired together only by `scripts/*.py` passing JSON
  files keyed by `--tag`. Fix: `annotate_deterministic(genes, organism)` and
  `assess_pairs(pairs)` as documented, semver-stable entry points — deliberately
  **after** batch 1, which changed the shape of exactly those calls.
* **F6 — lineage at variant level and for promoters.** `marks_position()` exists
  but is not a tool, and gene spans come from KEGG coding coordinates, so
  promoter features (inhA −15, the eis promoter) are never checked. Fix: a tool
  taking H37Rv positions plus an optional upstream window — the window is a
  scientific parameter, so it needs an explicit default recorded in provenance.
* **F7 — rule state in epistasis mode.** LCS rules say which genes are mutated
  *and* which are wild type (`phoP=1 AND phoR=0 AND gyrA=1`);
  `annotate_epistasis(genes)` loses that. Fix: optional
  `{gene: "mutated"|"wildtype"}` carried into the task text. Note this changes
  the question being asked, not just the input.

---

## Checked and not issues (from the same report)

* phoP (Rv0757) and phoR (Rv0758) genuinely have no UniProt function statement —
  only unreviewed P71814/P71815 — and genuinely have no KEGG pathway. Good
  negative controls for the gold set.
* `I6Y8F7` (Rv0678/mmpR5) really is Swiss-Prot despite the TrEMBL-shaped
  accession, so `reviewed=True` is correct.
