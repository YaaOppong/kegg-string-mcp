"""Identity is resolved once, for every consumer, or the answer depends on spelling."""

from __future__ import annotations

from types import SimpleNamespace

from kegg_string_mcp.identity import MIN_ALIAS_LENGTH, GeneIdentity, IdentitySet, resolve
from kegg_string_mcp.kegg import GeneIndex


class FakeString:
    def __init__(self, hits: dict[str, dict]):
        self.hits = hits

    def resolve(self, gene, species):
        return self.hits.get(gene), SimpleNamespace(url="", retrieved_at="", cached=True,
                                                    status=200, content_sha256="")


class FakeUniProt:
    def __init__(self, records: dict[str, list]):
        self.records = records

    def protein(self, gene, organism_id=83332):
        return SimpleNamespace(records=self.records.get(gene, []), notes=[], resolved={})


def _record(accession, gene_names, locus_tags, reviewed=True):
    return SimpleNamespace(record_id=accession,
                           detail={"gene_names": gene_names, "locus_tags": locus_tags,
                                   "reviewed": reviewed})


class FakeKegg:
    def __init__(self, entries, ambiguous=(), locus_tags=()):
        self.index = GeneIndex(entries=entries, ambiguous=set(ambiguous),
                               locus_tags=set(locus_tags))

    def gene_index(self, organism):
        return self.index, SimpleNamespace(url="", retrieved_at="", cached=True,
                                           status=200, content_sha256="")


def test_one_gene_under_two_spellings_is_one_identity():
    """Rv0678 and mmpR5 are the same protein. Every consumer that keyed on the
    string it was handed disagreed about that."""
    string = FakeString({"Rv0678": {"stringId": "83332.Rv0678", "preferredName": "mmpR5"}})
    identities = resolve(["Rv0678"], string=string)
    identity = identities.get("Rv0678")
    assert identity.string_id == "83332.Rv0678"
    assert identity.locus_tag == "Rv0678"
    assert "mmpR5" in identity.aliases
    # ...and the set is addressable by either spelling.
    assert identities.get("mmpR5") is identity
    assert identities.string_id("mmpR5") == "83332.Rv0678"


def test_the_query_string_is_always_the_first_alias():
    """Downstream artefacts are keyed by what the caller asked for; an alias list
    that reorders them would rename genes between pipeline stages."""
    string = FakeString({"mmpR5": {"stringId": "83332.Rv0678", "preferredName": "mmpR5"}})
    assert resolve(["mmpR5"], string=string).get("mmpR5").aliases[0] == "mmpR5"


def test_an_unresolved_gene_is_recorded_as_unresolved_not_dropped():
    string = FakeString({})
    identities = resolve(["nosuchgene"], string=string)
    identity = identities.get("nosuchgene")
    assert identity.matched_by["string"] == "none"
    assert not identity.answered("string")
    assert identities.unresolved("string") == ["nosuchgene"]
    assert identities.unanswered() == {"nosuchgene": ["string"]}


def test_a_source_that_raises_does_not_take_the_run_down():
    """Resolving 41 genes must not be all-or-nothing: the caller's alternative is
    the ad-hoc string matching this module replaces."""

    class Broken:
        def resolve(self, gene, species):
            raise ValueError("upstream is having a day")

    identity = resolve(["katG"], string=Broken()).get("katG")
    assert identity.aliases == ["katG"]
    assert not identity.answered("string")
    assert any("retrieval failure" in n for n in identity.notes)


def test_uniprot_gene_names_and_locus_tags_become_aliases():
    uniprot = FakeUniProt({"katG": [_record("P9WIE5", ["katG"], ["Rv1908c"])]})
    identity = resolve(["katG"], uniprot=uniprot).get("katG")
    assert identity.accession == "P9WIE5"
    assert "Rv1908c" in identity.aliases


def test_the_reviewed_uniprot_entry_supplies_the_accession():
    """A TrEMBL fragment sorting first should not become 'the' accession, nor
    should its gene names be the ones matched in text."""
    uniprot = FakeUniProt({"katG": [_record("A0AFRAGMENT", ["katG"], [], reviewed=False),
                                    _record("P9WIE5", ["katG"], ["Rv1908c"])]})
    assert resolve(["katG"], uniprot=uniprot).get("katG").accession == "P9WIE5"


def test_kegg_symbols_for_the_same_gene_become_aliases():
    kegg = FakeKegg(entries={"KATG": "mtu:Rv1908c", "RV1908C": "mtu:Rv1908c"},
                    locus_tags={"RV1908C"})
    identity = resolve(["katG"], kegg=kegg).get("katG")
    assert identity.kegg_gene_id == "mtu:Rv1908c"
    assert identity.locus_tag == "Rv1908c"


def test_a_symbol_kegg_records_for_several_genes_is_refused():
    """Aliasing fixes an undercount; it must not create an overcount. A wrong
    co-mention makes a pair look explained, which is the direction of error this
    work exists to remove -- and KEGG already knows which symbols are not unique."""
    kegg = FakeKegg(entries={"KATG": "mtu:Rv1908c", "RV1908C": "mtu:Rv1908c",
                             "SHARED": "mtu:Rv1908c"},
                    ambiguous={"SHARED"}, locus_tags={"RV1908C"})
    identity = resolve(["katG"], kegg=kegg).get("katG")
    assert "SHARED" not in identity.aliases
    assert "more than one gene" in identity.rejected["SHARED"]


def test_an_alias_claimed_by_two_genes_of_the_set_is_refused_from_both():
    uniprot = FakeUniProt({"geneA": [_record("A1", ["geneA", "shared"], [])],
                           "geneB": [_record("B1", ["geneB", "shared"], [])]})
    identities = resolve(["geneA", "geneB"], uniprot=uniprot)
    assert identities.ambiguous["shared"] == ["geneA", "geneB"]
    for query in ("geneA", "geneB"):
        assert "shared" not in identities.get(query).aliases
        assert "shared" in identities.get(query).rejected


def test_short_aliases_are_refused():
    """Two-letter tokens collide with ordinary words and with other identifiers."""
    uniprot = FakeUniProt({"katG": [_record("P9WIE5", ["katG", "kg"], [])]})
    identity = resolve(["katG"], uniprot=uniprot).get("katG")
    assert "kg" not in identity.aliases
    assert identity.rejected["kg"] == f"shorter than {MIN_ALIAS_LENGTH} characters"


def test_alias_map_is_what_the_corpus_stores():
    string = FakeString({"Rv0678": {"stringId": "83332.Rv0678", "preferredName": "mmpR5"}})
    assert resolve(["Rv0678"], string=string).alias_map() == {"Rv0678": ["Rv0678", "mmpR5"]}


def test_lookup_by_identifier_as_well_as_by_name():
    identities = IdentitySet(identities={
        "katG": GeneIdentity(query="katG", string_id="83332.Rv1908c",
                             kegg_gene_id="mtu:Rv1908c", aliases=["katG"])})
    assert identities.get("83332.Rv1908c").query == "katG"
    assert identities.get("mtu:Rv1908c").query == "katG"
    assert identities.get("unknown") is None


# -- sources rescue each other ------------------------------------------------


def test_a_symbol_string_refuses_is_retried_as_the_locus_tag():
    """Nomenclature drift is not an absent gene: a source that fails on a symbol
    is retried with the locus tag another source supplied.

    The fixture is CONSTRUCTED, not observed. On the TB-41 set every real gene
    resolves in STRING under both its symbol and its locus tag -- `mmpR5` and
    `Rv0678` both return 83332.Rv0678 -- so this retry is defensive. It costs a
    resolved gene nothing (see the call-count test below) and only runs for a
    source that already failed.
    """
    string = FakeString({"Rv0678": {"stringId": "83332.Rv0678", "preferredName": "mmpR5"}})
    uniprot = FakeUniProt({"mmpR5": [_record("I6Y8F7", ["mmpR5"], ["Rv0678"])]})
    identity = resolve(["mmpR5"], string=string, uniprot=uniprot).get("mmpR5")
    assert identity.string_id == "83332.Rv0678"
    assert identity.answered("string")
    assert identity.resolved_via["string"] == "Rv0678"
    assert identity.tried["string"] == ["mmpR5", "Rv0678"]


def test_the_rescuing_spelling_is_recorded_not_silently_substituted():
    """A silent retry makes the run unreproducible by hand: the reader cannot tell
    which identifier the answer is about."""
    string = FakeString({"Rv0678": {"stringId": "83332.Rv0678", "preferredName": "mmpR5"}})
    uniprot = FakeUniProt({"mmpR5": [_record("I6Y8F7", ["mmpR5"], ["Rv0678"])]})
    identity = resolve(["mmpR5"], string=string, uniprot=uniprot).get("mmpR5")
    assert any("did not resolve 'mmpR5' but did resolve 'Rv0678'" in n for n in identity.notes)


def test_kegg_is_retried_with_what_string_supplied():
    kegg = FakeKegg(entries={"RV0678": "mtu:Rv0678"}, locus_tags={"RV0678"})
    string = FakeString({"mmpR5": {"stringId": "83332.Rv0678", "preferredName": "mmpR5"}})
    identity = resolve(["mmpR5"], string=string, kegg=kegg).get("mmpR5")
    assert identity.kegg_gene_id == "mtu:Rv0678"
    assert identity.resolved_via["kegg"] == "Rv0678"


def test_a_gene_no_source_can_resolve_is_still_unresolved():
    """The retry must not manufacture a resolution: exhausting every spelling is
    the evidence that the gene is genuinely absent, not misspelt."""
    identities = resolve(["nosuchgene"], string=FakeString({}), uniprot=FakeUniProt({}))
    identity = identities.get("nosuchgene")
    assert identity.unresolved_sources == ["string", "uniprot"]
    assert identities.unanswered() == {"nosuchgene": ["string", "uniprot"]}


def test_retrying_costs_one_extra_call_per_failing_source_not_per_gene():
    """A resolved gene must not pay for the rescue path."""
    calls = []

    class Counting(FakeString):
        def resolve(self, gene, species):
            calls.append(gene)
            return super().resolve(gene, species)

    string = Counting({"katG": {"stringId": "83332.Rv1908c", "preferredName": "katG"}})
    uniprot = FakeUniProt({"katG": [_record("P9WIE5", ["katG"], ["Rv1908c"])]})
    resolve(["katG"], string=string, uniprot=uniprot)
    assert calls == ["katG"]


# -- last resort: a proposed identifier, verified before it counts -------------


def test_a_proposed_identifier_is_accepted_only_after_a_source_confirms_it():
    """The case this exists for: a gene renamed since the query was written. The
    structured sources cannot map the old name, but a paper can -- and the paper's
    suggestion is then put back to the structured sources."""
    string = FakeString({"Rv0678": {"stringId": "83332.Rv0678", "preferredName": "mmpR5"}})
    proposals = [("Rv0678", "PMID 31932378 calls mmpR5 'Rv0678'")]
    identity = resolve(["oldName1"], string=string,
                       proposer=lambda i: proposals).get("oldName1")
    assert identity.string_id == "83332.Rv0678"
    assert identity.answered("string")
    assert identity.proposals == {"Rv0678": "PMID 31932378 calls mmpR5 'Rv0678'"}
    assert any("was confirmed by string" in n for n in identity.notes)


def test_an_unconfirmed_proposal_leaves_the_gene_unresolved():
    """The guard that makes the whole idea safe: a model may propose a name, but
    only a structured source can make it true. A hallucinated identifier fails
    closed instead of propagating as a resolved gene."""
    string = FakeString({})          # STRING confirms nothing
    identity = resolve(["oldName1"], string=string,
                       proposer=lambda i: [("Rv9999", "a confident-sounding sentence")]
                       ).get("oldName1")
    assert not identity.answered("string")
    assert identity.string_id == ""
    assert "no structured source resolved it" in identity.rejected["Rv9999"]
    assert any("remains unresolved" in n for n in identity.notes)


def test_the_proposer_is_not_called_when_everything_already_resolved():
    """It is the last resort, not a step: a resolved gene must not pay for it."""
    calls = []
    string = FakeString({"katG": {"stringId": "83332.Rv1908c", "preferredName": "katG"}})

    def proposer(identity):
        calls.append(identity.query)
        return []

    resolve(["katG"], string=string, proposer=proposer)
    assert calls == []


def test_a_proposer_that_raises_does_not_end_the_run():
    def proposer(identity):
        raise RuntimeError("the model is having a day")

    identity = resolve(["oldName1"], string=FakeString({}), proposer=proposer).get("oldName1")
    assert not identity.answered("string")
    assert any("proposer failed" in n for n in identity.notes)


def test_proposals_are_capped():
    seen = []

    def proposer(identity):
        names = [(f"Rv{i}", "guess") for i in range(10)]
        seen.extend(names)
        return names

    identity = resolve(["oldName1"], string=FakeString({}), proposer=proposer,
                       max_proposals=2).get("oldName1")
    assert len(identity.proposals) == 2
