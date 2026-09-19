# The seven tools

All are read-only, expose structured output schemas, and are deterministic: same
input and same cache produce the same output. Defaults shown are for
*M. tuberculosis* H37Rv.

## kegg_pathways

`gene` (KEGG ID, locus tag, or symbol), `organism` = `mtu`

One record per pathway: `record_id` (e.g. `mtu00360`), name, entry URL.

Covers about 1,170 of 4,008 genes. `gyrA`, one of the most studied TB genes, has
none. Empty means unannotated far more often than it means uninvolved.

## uniprot_protein

`gene`, `organism_id` = `83332`, `limit` = `3`

One record per entry: `record_id` (accession, e.g. `P9WG47`), protein name,
function statements **tiered by evidence code** with supporting PMIDs, catalytic
activity, PDB cross-refs, `quotable_text`.

The evidence tier matters: an experimentally-supported function and one inferred
from homology are different claims. Check UniProt whenever KEGG is silent.

## string_partners

`gene`, `species` = `83332`, `limit` = `20`, `required_score` = `700`

One record per partner: `record_id` (e.g. `83332.Rv1909c`), preferred name,
combined score, **full per-channel breakdown**, network URL.

Read `evidence_beyond_textmining` before presenting any partner as experimental
support. If the partner list hits `limit`, true network degree is unknown and
absence of an edge is a limit of the query.

## pubmed_abstracts

`gene`, `organism` = `Mycobacterium tuberculosis`, `limit` = `10`

One record per article: `record_id` (PMID), title, abstract, `quotable_text`,
journal, year, DOI, URL.

Relevance-ranked over titles and abstracts only. Quote from `quotable_text` —
that is the exact retrieved string the validator checks against.

## corpus_search

`query` (free text), `limit` = `5`

One record per paper: `record_id` (PMID), `quotable_text` (the whole abstract),
`matched_passage` (the span that ranked), score, rank.

BM25 fused with dense embeddings; 0.917 precision@10 against 0.844 for keyword
alone. Covers only a prebuilt gene set — its `notes` says which genes it holds.
Inert until `KEGG_STRING_MCP_CORPUS` points at a corpus from
`scripts/build_corpus.py`. Empty means out of coverage, so fall back to
`pubmed_abstracts` rather than concluding anything.

## lineage_markers

`gene`, `organism` = `mtu`

One record per lineage-defining SNP the gene contains: `record_id` (e.g.
`tbdb:851797`), the lineage it marks, H37Rv position, allele, source URL.

Call for every gene. 855 of 4,008 genes contain one, so a positive does not mean
the variant in question *is* a marker — report it as a confound to test against
genotype data, never as a conclusion.

## resistance_variants

`gene`, `drug` (optional)

Whether the gene is resistance-associated, the drugs, per-grade counts, and one
record per catalogued variant: `record_id` (e.g. `tbdb:katG:p.Ser315Thr`), WHO
grading, drug, source.

A gene is resistance-associated if **any** catalogued variant is graded
associated, however many are not. Keep the three negatives apart: absent from the
catalogue, present and negative, and "Uncertain significance".
