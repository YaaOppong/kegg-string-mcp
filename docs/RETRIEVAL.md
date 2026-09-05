# The retrieval comparison

Two retrieval arms over one corpus of PubMed abstracts, measured against each
other on the same queries. The comparison is the contribution; adding a vector
store to a repo is a tutorial.

Reproduce with:

```bash
NCBI_EMAIL=you@example.org python scripts/build_corpus.py --extended --all-genes --tag tb41
python scripts/run_comparison.py data/corpus_tb41.json --tag tb41
```

`--all-genes` is not optional here. In the pipeline, stage 2 runs only on genes
whose structured annotation is thin, because covering what stage 1 already
describes is wasted effort. A head-to-head on gene *pairs* needs every gene in
the corpus regardless of coverage, so the measurement bypasses the routing that
the pipeline applies. Without the flag this builds a 22-gene corpus and every
number below changes.

## The corpus

41 M. tuberculosis genes, 676 unique papers, chunked to 1,350 passages of 180
words with 40 words of overlap. Chunking is not cosmetic: the ONNX MiniLM
embedder truncates at 256 wordpieces, 80% of these abstracts run longer, and 18
of them named a gene only past the cut. Un-chunked, those genes were invisible
to the dense arm for reasons that had nothing to do with retrieval.

Corpora are gitignored. The repo redistributes no bulk PubMed content.

## The arms

| Arm | Mechanism | Fails on |
|---|---|---|
| `lexical` | BM25 over tokens | paraphrase; a query that shares no words with the passage |
| `dense` | MiniLM embeddings in Chroma, cosine | exact identifiers -- `Rv1908c` has no useful neighbourhood in embedding space |
| `hybrid` | reciprocal rank fusion of the two, smoothing 60 | inherits both, less often |

RRF fuses on rank, not score, because BM25 scores and cosine similarities are
not on a common scale and normalising them across arms invents a calibration
neither arm has.

## Relevance without hand-labelling

Each passage carries `genes_named`: the genes its own text names, computed over
the whole corpus with word boundaries. This is separate from `mentions`, which
records which *queried* terms appear -- a passage retrieved for katG names ahpC
whether or not ahpC was asked for. Conflating the two undercounted pair evidence
by a factor of three before it was caught.

So "did this return papers that discuss the gene asked about?" has an exact
answer, reproducible and impossible to tune after the fact. It is weaker than
human judgement -- a gene can be named in passing -- but it costs nothing.

## Result: all 820 pairs

| Arm | precision@10 | papers naming both genes |
|---|---|---|
| hybrid | **0.917** | 0.31 |
| lexical | 0.844 | **0.39** |
| dense | 0.739 | 0.15 |

Arm overlap (Jaccard over returned PMIDs): dense/lexical 0.158, hybrid/lexical
0.433, hybrid/dense 0.423. The dense arm is retrieving substantially different
papers, which is the case for keeping it; it is also the least precise, which is
the case against replacing BM25 with it.

## Removing the circularity

Scoring relevance by whether a passage names the queried genes is close to what
BM25 ranks on, which tilts the measurement toward the lexical arm. The fix is a
query set chosen by something neither retriever can see: STRING.

The classification below lives in `retrieval/independence.py` and is also what
the residue gate consumes (`hypothesis/residue.py`), where it does a different
job -- see "What the silent pairs turned out to be" at the end.

Every pair is classified by one STRING call per gene (820 pairs would be 820
requests against a service that asks for roughly one per second):

| Status | Pairs | What literature adds |
|---|---|---|
| `silent` | 455 (55%) | everything -- STRING returns no edge |
| `textmining_only` | 326 (40%) | nothing new; STRING's score came from this literature |
| `corroborating` | 39 (5%) | confirmation of an experimental or database channel |

The middle row is the one to sit with. For 40% of pairs STRING reports a high
combined score -- katG/pncA at 0.965 -- built almost entirely from textmining,
with every other channel below 0.11. Presenting both STRING's score and the
retrieved abstracts as evidence is one line of evidence counted twice. The repo
already refuses to do that at the tool level via `evidence_beyond_textmining`;
this applies the same rule to the retrieval comparison.

28% of the corpus postdates the STRING v12.0 release, so those papers cannot be
in any channel, textmining included.

## Result: the 455 pairs STRING is silent on

| Arm | precision@10 | papers naming both genes |
|---|---|---|
| hybrid | **0.908** | 0.04 |
| lexical | 0.854 | **0.05** |
| dense | 0.714 | 0.02 |

Two findings, one comfortable and one not.

The arm ranking is unchanged. Hybrid leads on precision, lexical on joint
evidence, dense trails on both, and overlaps move by less than 0.005. The
comparison survives removing the circularity, which is the result it needed to
survive.

Joint evidence collapses -- 0.39 papers per query to 0.05 for the lexical arm.
The pairs STRING is silent on are largely pairs this corpus is silent on too.
Almost all co-mention evidence sits on pairs STRING already scores, and most of
those scores are textmining. On this corpus, the literature arm is mostly
recovering what STRING already encodes rather than reaching past it.

That is a finding about the corpus, not about the retrievers. It was built by
querying PubMed for the genes themselves, which retrieves the well-studied
drug-resistance pairs and little else. Testing whether literature reaches past
STRING needs a corpus built to find that -- topic and method queries rather than
gene queries -- and the machinery for that measurement now exists.

## What the silent pairs turned out to be

The 455 STRING-silent pairs were constructed here as a de-biased query set: a
selection neither retriever influences, used to check that the arm ranking was
not an artefact of scoring relevance on gene names.

They are also the input to hypothesis generation, and there the polarity is
reversed. A pair no structured source connects and no paper co-mentions is not a
disappointing query -- it is a candidate. 441 of the 455 (97%) have no paper
naming both genes, against 83% across all 820 pairs, so the filter discriminates.

This makes the collapse above read differently depending on what is being asked.
As a retrieval measurement it says the corpus has little to offer on these pairs.
As a discovery filter it is the expected signature of a relationship nobody has
written down, which is the only kind worth generating a hypothesis about. Low
co-mention is necessary and nowhere near sufficient: most gene pairs are
unconnected because they are unrelated.

Neither reading changes the arm comparison. The ranking is the same on both
query sets.

## What is not measured

No reranking, no hybrid weight tuning, no query expansion. Nothing here is a
temporal holdout: every number is measured on a corpus containing everything
published to date, which is the right setup for comparing retrievers and the
wrong one for judging whether a method could have found something before it was
known. Locus tags such as
`Rv1908c` fail in every arm because the corpus text uses symbols; the exact-term
probe reports `in_corpus` alongside each arm's hit count so a zero is not
ambiguous between "retrieval missed it" and "it is not there".
