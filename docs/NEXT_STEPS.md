# Next steps: widening the literature search

The retrieval arm works and is measured. What it does not yet do is reach past
what STRING already encodes, and the reason is in how the corpus is built rather
than in how it is ranked.

This document is the plan for fixing that, in the order the work should happen,
with what each change must be measured against. Nothing here is started.

## The problem, stated with the numbers

On all 820 gene-pair queries the hybrid arm returns 0.31 papers naming both genes
per query. On the 379 pairs STRING holds no edge for, that falls to **0.03**.
Almost all co-mention evidence sits on pairs STRING already scores, and 49% of
those scores are textmining -- which is co-mention in papers. So the literature
arm is largely recovering the same signal STRING built its score from.

That is a property of how the corpus was selected, not of the retrievers:

* one PubMed query per gene, `"{symbol}" AND "{organism}"`, phrase-quoted;
* the 20 most relevant records by NCBI's own ranking;
* titles and abstracts only, no full text.

Each of those bounds what any ranking method can find. A paper that never writes
the configured symbol is not in the corpus, and no amount of alias-aware matching
afterwards can recover it -- labelling only works on papers already held.

## The asymmetry to fix first

`identity.py` resolves every gene to its full alias set (UniProt gene names and
locus tags, KEGG symbols, STRING preferred name) and refuses any alias that names
more than one gene. `build_corpus.py` calls it before fetching and passes the
result into `build()`.

`build()` then searches PubMed with the configured symbol alone. The alias set is
used only at the end, by `annotate_genes_named`, to label which genes each
retrieved passage mentions.

So the corpus is **retrieved by one spelling and measured across all of them**.
Every number in `docs/RETRIEVAL.md` carries that asymmetry.

## Workstream 1 -- query by every alias

One OR'd query per gene rather than one per spelling, so the 20-record cap stays
per gene rather than silently becoming 20 per alias:

```
("katG" OR "Rv1908c") AND "Mycobacterium tuberculosis"
```

Record which spelling retrieved each paper as `queried_via`, mirroring
`named_via` on the matching side. An alias set that silently changes corpus
membership is not auditable; one that records its own effect is.

Expected: a modest rise in papers per gene, concentrated on genes whose symbol
differs between sources (dosR/devR, gidB/gid, icl1/icl). Low risk -- ambiguous
aliases are already refused upstream.

## Workstream 2 -- query by protein name

The larger prize, and the one aimed at the problem above.

Gene symbols are the vocabulary of genetics: mutation surveys, resistance
catalogues, knockouts. Protein names are the vocabulary of biochemistry and
structural work: enzymology, crystallography, binding assays. These are two
literatures about the same molecule with little shared vocabulary, and the second
is where mechanistic evidence lives.

It also bears directly on the independence question. STRING's textmining channel
is built from co-mention, and symbol co-mention is what gene-name search
retrieves -- so the current corpus draws from roughly the well that produced the
score it is being compared against. Protein-name literature is a different well.

UniProt already returns `protein_name`, and `uniprot.py` already parses
`recommendedName.fullName` -- that is where "Transcriptional regulator FurA"
comes from in a run. What is missing is carrying it into the identity set and
into the query.

Three cautions:

* **Generic names must be excluded.** "Histidine kinase", "response regulator"
  and "hypothetical protein" name hundreds of proteins. Skip a name that matches
  more than one gene in the same organism -- the same instinct as refusing an
  ambiguous alias.
* **The organism clause carries the weight.** "Catalase-peroxidase" alone
  retrieves every organism's.
* **It changes what `genes_named` should mean.** A paper that says
  "catalase-peroxidase" throughout and never `katG` does not count as naming the
  gene under the current rule, so the corpus would improve while the measured
  numbers got worse. Decide the rule before scoring: either report co-mention
  twice, by symbol and by symbol-or-product, or state which was used. `named_via`
  already records which spelling matched, so both readings are available.

## Workstream 3 -- stop hand-writing the organism

The recorded query translation for katG is:

```
"katG"[All Fields] AND "Mycobacterium tuberculosis"[All Fields]
```

`[All Fields]` means PubMed did no MeSH mapping: the quotes suppressed its own
Automatic Term Mapping, which is what normally covers "M. tuberculosis", "MTB"
and the rest. `[All Fields]` still searches MeSH terms and substance lists, so an
indexed paper matches however its abstract spells the organism -- but recent
papers are not yet indexed.

Rather than hand-maintaining a synonym list, let MeSH do it while keeping the
gene term literal:

```
"katG" AND "Mycobacterium tuberculosis"[MeSH Terms]
```

Cheap to test: run katG three ways -- as now, MeSH-scoped, and with a hand-written
OR list -- and compare counts and the top-20 PMID sets. If MeSH matches the OR
list, use MeSH and let NCBI maintain the synonyms.

## Workstream 4 -- past the abstract, and past PubMed

Everything above still searches titles and abstracts. Three ways past that, in
increasing order of effort:

* **PMC full-text search.** PMC searches the body text PubMed cannot, which is
  where a mechanism is usually described rather than summarised. Restricted to
  the open-access subset, per-article licences vary, and the repo already treats
  `in_pmc` as the licit route to full text. The corpus manifest was written for
  exactly this consumer.
* **Citation-graph expansion.** Take the papers already retrieved and follow
  references and citing articles. This reaches work that shares no vocabulary
  with the query at all, which is the whole point -- a mechanistically related
  paper that names neither gene in its abstract can still be two hops away in the
  citation graph.
* **Topic and method queries.** `docs/RETRIEVAL.md` already names this as what it
  would take to test whether literature reaches past STRING: build a corpus from
  queries about mechanisms ("compensatory mutation", "efflux-mediated
  resistance", "transcriptional repression under oxidative stress") rather than
  from gene names. A corpus selected by topic cannot be accused of recovering
  gene co-mention by construction.

The last of these changes what the corpus *is*, so it belongs in its own tagged
build and its own section of the write-up.

## How each change gets measured

Every claim in this repo is tied to a specific corpus, so none of this is done in
place. Build alongside, using the tag each script already takes:

```bash
python scripts/build_corpus.py --extended --all-genes --tag tb41alias
python scripts/run_comparison.py data/corpus_tb41alias.json --tag tb41alias
python scripts/residue.py --tag tb41alias
```

Then one table, three rows, same four columns:

| Corpus | Papers | precision@10 | naming both | naming both, STRING-silent |
|---|---|---|---|---|
| `tb41` (symbol only) | 676 | 0.917 | 0.31 | 0.03 |
| `tb41alias` | | | | |
| `tb41prot` | | | | |

The last column is the one that matters. It is the weakest result in the repo and
the one a wider query should move. If it does not move, that is a finding too:
the pairs are unstudied rather than the search too narrow, which is the reading
the residue already takes.

Re-measuring is not optional. `docs/RETRIEVAL.md`, the README and the demo all
quote corpus-dependent figures, and a rebuild under the same tag would silently
invalidate them.

## Open questions, to settle before building

1. **What counts as naming a gene** once protein names are in play. Symbol only,
   symbol-or-product, or both reported separately.
2. **Whether the 20-record cap should rise** with the query widened. More
   spellings against the same cap changes which papers are dropped, and the cap
   was never measured -- a `--limit 50` build on a few genes would say whether
   pair evidence appears past rank 20.
3. **Whether the epistasis prompt should be told the co-mention count.** The
   corpus can already count papers naming both genes, deterministically, and the
   pipeline currently leaves that to one model-written query returning five
   papers. That is the same argument the repo makes for pathway intersections:
   counting belongs in code.

## What this does not change

The ranking arms, the fusion, the relevance criterion and the citation checking
all stay as they are. This is entirely about what enters the corpus. Widening the
candidate set and then measuring the same way is the only way to tell a better
corpus from a better-looking number.
