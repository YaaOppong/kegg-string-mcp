"""The stage 2 trigger: which genes literature retrieval should run on."""

from __future__ import annotations

from types import SimpleNamespace

from kegg_string_mcp.retrieval.coverage import (
    NO_EXPERIMENTAL,
    NO_FUNCTION,
    NO_PATHWAY,
    assess,
    classify,
    route,
    summarise,
)


def _result(records, traced=True, matched_by=None):
    if matched_by is None:
        matched_by = "symbol" if records else "none"
    return SimpleNamespace(records=records, resolved={"matched_by": matched_by},
                           requests=[SimpleNamespace(url="u")] if traced else [])


def _protein(statements, experimental, traced=True, resolved=True):
    """A UniProt entry with no FUNCTION statement is still a resolved entry --
    the shape phoP (P71814) has: the record exists, its function is undescribed."""
    detail = {"function_statements": [{}] * statements,
              "has_experimental_function": experimental}
    records = [SimpleNamespace(detail=detail)] if resolved else []
    return _result(records, traced=traced,
                   matched_by="uniprot_search" if resolved else "none")


class FakeKegg:
    def __init__(self, pathways, traced=True):
        self.pathways_by_gene, self.traced = pathways, traced

    def pathways(self, gene, organism="mtu"):
        # A resolved gene with no pathway is the common case, so matched_by is
        # set even when no records come back.
        records = [SimpleNamespace(detail={"kegg_gene_id": f"mtu:{gene}"})
                   ] * self.pathways_by_gene.get(gene, 0)
        return _result(records, traced=self.traced, matched_by="symbol")


class FakeUniProt:
    def __init__(self, data, traced=True, unresolved=(), locus=None):
        self.data, self.traced, self.unresolved = data, traced, set(unresolved)
        self.locus, self.asked = locus, []

    def protein(self, gene, organism_id=83332):
        self.asked.append(gene)
        statements, experimental = self.data.get(gene, (0, False))
        result = _protein(statements, experimental, traced=self.traced,
                          resolved=gene not in self.unresolved)
        for record in result.records:
            record.detail["locus_tags"] = [self.locus] if self.locus else []
        return result


def test_well_annotated_gene_is_not_routed():
    coverage = classify("katG", kegg_pathways=3, uniprot_statements=5,
                        has_experimental_function=True)
    assert coverage.reasons == []
    assert not coverage.thin


def test_no_uniprot_function_is_the_strongest_gap():
    coverage = classify("Rv0001", 0, 0, False)
    assert coverage.reasons == [NO_FUNCTION, NO_PATHWAY]
    assert coverage.functional_gap


def test_inferred_only_annotation_counts_as_a_functional_gap():
    """A statement inferred from a rule or a homologue is not evidence about
    this protein -- the same distinction uniprot.py draws on evidence tiers."""
    coverage = classify("Rv0002", kegg_pathways=2, uniprot_statements=4,
                        has_experimental_function=False)
    assert coverage.reasons == [NO_EXPERIMENTAL]
    assert coverage.functional_gap


def test_missing_pathway_alone_is_a_weak_signal():
    """gyrA has no KEGG pathway and is among the best-characterised genes in TB.
    That is a fact about KEGG's coverage, not about the gene."""
    coverage = classify("gyrA", kegg_pathways=0, uniprot_statements=6,
                        has_experimental_function=True)
    assert coverage.reasons == [NO_PATHWAY]
    assert coverage.thin
    assert not coverage.functional_gap


def test_route_can_demand_a_functional_gap():
    coverages = [classify("gyrA", 0, 6, True),      # no pathway only
                 classify("Rv0002", 2, 4, False)]   # inferred only
    assert route(coverages) == ["gyrA", "Rv0002"]
    assert route(coverages, functional_only=True) == ["Rv0002"]


def test_assess_reads_both_sources():
    coverages = assess(["katG", "gyrA"],
                       FakeKegg({"katG": 2}),
                       FakeUniProt({"katG": (3, True), "gyrA": (2, False)}))
    assert [c.gene for c in coverages] == ["katG", "gyrA"]
    assert coverages[0].reasons == []
    assert coverages[1].reasons == [NO_EXPERIMENTAL, NO_PATHWAY]


def test_a_failed_lookup_is_not_read_as_a_missing_annotation():
    """The dangerous case: a 503 routes a well-annotated gene to literature
    retrieval, or worse, is recorded as 'nothing is known about this gene'."""
    coverages = assess(["katG"], FakeKegg({}, traced=False), FakeUniProt({}))
    coverage = coverages[0]
    assert coverage.lookup_failed
    assert not coverage.thin            # not routed
    assert not coverage.functional_gap
    assert route(coverages) == []


def test_failed_lookups_are_reported_not_absorbed():
    """A failure must not vanish into either bucket -- it is a retry list."""
    coverages = assess(["katG", "gyrA"],
                       FakeKegg({"katG": 2}),
                       FakeUniProt({"katG": (3, True)}, traced=False))
    summary = summarise(coverages)
    assert summary["lookup_failed"] == {"katG": ["uniprot"], "gyrA": ["uniprot"]}
    assert summary["well_covered"] == 0
    assert summary["thin"] == 0


def test_an_unresolved_symbol_is_not_a_missing_annotation():
    """UniProt has no entry for 'icl1'; the same protein under 'icl' or Rv0467
    is annotated. Reading that as 'nothing is known' routed a well-studied gene
    to literature retrieval on the strength of a wrong symbol."""
    coverages = assess(["icl1"], FakeKegg({"icl1": 6}),
                       FakeUniProt({}, unresolved=["icl1"]))
    coverage = coverages[0]
    assert coverage.unknown == ["uniprot"]
    assert NO_FUNCTION not in coverage.reasons
    assert not coverage.thin
    assert route(coverages) == []


def test_a_resolved_entry_with_no_function_statement_is_a_real_gap():
    """phoP resolves to P71814, which exists and describes no function. That is
    the finding the trigger is for, and must survive the unresolved guard."""
    coverages = assess(["phoP"], FakeKegg({}), FakeUniProt({"phoP": (0, False)}))
    coverage = coverages[0]
    assert coverage.unknown == []
    assert coverage.reasons == [NO_FUNCTION, NO_PATHWAY]
    assert route(coverages) == ["phoP"]


def test_one_source_failing_does_not_void_the_other():
    """KEGG resolved and reported no pathway; that stays known even though
    UniProt could not answer. The gene is still not routed, because a partial
    picture is not a basis for saying what is missing."""
    coverage = classify("x", 0, 0, False, kegg_known=True, uniprot_known=False)
    assert coverage.reasons == [NO_PATHWAY]
    assert coverage.unknown == ["uniprot"]
    assert not coverage.thin


def test_summary_counts_reasons_and_buckets():
    coverages = [classify("a", 2, 3, True),      # covered
                 classify("b", 0, 3, True),      # no pathway
                 classify("c", 0, 0, False)]     # no function + no pathway
    summary = summarise(coverages)
    assert summary["genes"] == 3
    assert summary["well_covered"] == 1
    assert summary["thin"] == 2
    assert summary["functional_gap"] == 1
    assert summary["reason_counts"] == {NO_PATHWAY: 2, NO_FUNCTION: 1}


class KeggWithLocus:
    """KEGG resolving only the locus tag, as it does for dosR (Rv3133c)."""

    def __init__(self, locus, pathways):
        self.locus, self.pathways_count, self.asked = locus, pathways, []

    def pathways(self, gene, organism="mtu"):
        self.asked.append(gene)
        if gene != self.locus:
            return _result([], matched_by="none")
        records = [SimpleNamespace(detail={"kegg_gene_id": f"mtu:{self.locus}"})
                   ] * self.pathways_count
        return _result(records, matched_by="kegg_id")


def test_a_symbol_one_source_lacks_is_retried_with_the_locus_tag():
    """KEGG calls Rv3133c devR and has no 'dosR'. UniProt resolves dosR and
    returns the locus tag, which KEGG does accept -- so a nomenclature
    difference stops looking like a missing pathway."""
    kegg = KeggWithLocus("Rv3133c", pathways=2)
    uniprot = FakeUniProt({"dosR": (3, True)}, locus="Rv3133c")
    coverage = assess(["dosR"], kegg, uniprot)[0]
    assert kegg.asked == ["dosR", "Rv3133c"]
    assert coverage.resolved_via == "Rv3133c"
    assert coverage.unknown == []
    assert coverage.reasons == []          # not routed on a naming difference


def test_the_fallback_works_in_the_other_direction_too():
    """UniProt has no 'icl1'; KEGG resolves it and supplies Rv0467."""
    kegg = KeggWithLocus("icl1", pathways=6)
    kegg.pathways = lambda gene, organism="mtu": _result(
        [SimpleNamespace(detail={"kegg_gene_id": "mtu:Rv0467"})] * 6, matched_by="symbol")
    uniprot = FakeUniProt({"Rv0467": (2, True)}, unresolved=["icl1"])
    coverage = assess(["icl1"], kegg, uniprot)[0]
    assert uniprot.asked == ["icl1", "Rv0467"]
    assert coverage.resolved_via == "Rv0467"
    assert coverage.unknown == []


def test_no_retry_when_both_sources_already_answered():
    kegg = KeggWithLocus("katG", pathways=2)
    uniprot = FakeUniProt({"katG": (3, True)}, locus="Rv1908c")
    coverage = assess(["katG"], kegg, uniprot)[0]
    assert uniprot.asked == ["katG"]        # no second call
    assert coverage.resolved_via == ""


def test_an_unresolvable_gene_still_reports_unknown():
    """The fallback must not paper over a gene neither source knows."""
    kegg = KeggWithLocus("Rv9999", pathways=0)
    uniprot = FakeUniProt({}, unresolved=["nosuchgene"])
    coverage = assess(["nosuchgene"], kegg, uniprot)[0]
    assert coverage.unknown == ["uniprot", "kegg"]
    assert not coverage.thin
