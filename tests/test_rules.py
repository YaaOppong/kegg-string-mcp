"""Feature resolution and rule parsing. No network, no model.

Every test here is a case a real vocabulary produced: 2,903 labels against a
second annotation gave 1,346 strand-suffix disagreements in both directions, 76
labels absent entirely, and a family of `Rv0063a`-style names where the trailing
letter is part of the gene rather than a strand marker.
"""

from pathlib import Path

import pytest

from kegg_string_mcp.rules import parse_annotation, parse_rules, resolve_all, vocabulary
from kegg_string_mcp.rules.annotation import normalise_locus
from kegg_string_mcp.rules.features import CODING, INTERGENIC, UNRESOLVED, Resolver

# chrom, 0-based start, end, name, score, strand, symbol
BED_ROWS = [
    ("AL123456", 99, 199, "Rv0001", ".", "+", ""),
    ("AL123456", 299, 399, "Rv0002c", ".", "-", ""),
    ("AL123456", 449, 549, "Rv0003", ".", "+", "rpoC"),
    ("AL123456", 499, 600, "Rv0004", ".", "+", ""),     # overlaps Rv0003: no gap
    ("AL123456", 699, 799, "Rv0005c", ".", "-", ""),
    ("AL123456", 899, 950, "Rv0063", ".", "+", ""),
    ("AL123456", 959, 1000, "Rv0063a", ".", "+", ""),   # a DIFFERENT gene
    ("AL123456", 1099, 1150, "Rv0070", ".", "+", ""),
    ("AL123456", 1199, 1250, "Rv0070c", ".", "-", ""),  # collides once 'c' is dropped
    ("AL123456", 1299, 1350, "Rv0080", ".", "+", "TB7.3"),  # the number IS the name
]


@pytest.fixture
def bed(tmp_path: Path) -> Path:
    path = tmp_path / "h37rv_genes.bed"
    path.write_text("\n".join("\t".join(str(c) for c in row) for row in BED_ROWS) + "\n")
    return path


@pytest.fixture
def annotation(bed: Path):
    return parse_annotation(bed)


# --- the annotation --------------------------------------------------------


def test_bed_coordinates_become_one_based_inclusive(annotation):
    """BED is 0-based half-open and GFF is 1-based inclusive. Mixing them shifts
    every start by one, which stays invisible until a variant sits on a boundary
    and lands in the wrong feature."""
    gene = annotation.gene("Rv0001")
    assert (gene.start, gene.end) == (100, 199)
    assert gene.length == 100


def test_the_format_is_detected_from_the_row_not_the_extension(tmp_path):
    """A GTF named `.bed` would otherwise be read with the wrong coordinate
    convention and every start would be one out."""
    path = tmp_path / "misnamed.bed"
    path.write_text(
        "AL123456\tena\tgene\t100\t199\t.\t+\t.\tlocus_tag=Rv0001;gene_name=dnaA\n")
    annotation = parse_annotation(path)
    assert annotation.source_format == "gff"
    # 1-based inclusive taken verbatim, not shifted as a BED start would be.
    assert (annotation.gene("Rv0001").start, annotation.gene("Rv0001").end) == (100, 199)


def test_the_annotation_records_its_own_hash(annotation, bed):
    """Every HTTP response in a run store carries a content hash. The one local
    file that decides what every feature means deserves the same, or a reference
    update changes the results with nothing in the record to say so."""
    assert len(annotation.sha256) == 64
    assert annotation.summary()["genes"] == len(BED_ROWS)
    assert annotation.summary()["path"] == str(bed)


def test_overlapping_genes_produce_no_interval(annotation):
    """Rv0003 and Rv0004 overlap, so that junction has no feature. That is a fact
    about the genome, not a parsing failure -- around 900 H37Rv pairs overlap."""
    assert annotation.intergenic("Rv0003-Rv0004") is None
    assert annotation.intergenic("Rv0002c-Rv0003") is not None


def test_intergenic_intervals_sit_between_the_flanking_genes(annotation):
    interval = annotation.intergenic("Rv0001-Rv0002c")
    assert (interval.start, interval.end) == (200, 299)
    assert interval.width == 100


def test_promoter_orientation_reads_both_flanks(annotation):
    """A gap upstream of a forward-strand gene is where its promoter would be.
    Divergent flanks qualify on both sides; convergent flanks on neither."""
    divergent = annotation.intergenic("Rv0002c-Rv0003")     # <--  -->
    assert {g.locus for g in divergent.promoter_of()} == {"Rv0002c", "Rv0003"}

    convergent = annotation.intergenic("Rv0001-Rv0002c")    # -->  <--
    assert convergent.promoter_of() == []


# --- resolving coding labels -----------------------------------------------


def test_a_label_differing_only_in_the_strand_suffix_still_resolves(annotation):
    """1,346 of 2,879 real labels disagreed with a second annotation in the
    trailing `c` alone, in both directions. A tag that differs only there is the
    same gene."""
    resolver = Resolver(annotation)

    dropped = resolver.resolve("Rv0002")        # annotation says Rv0002c
    assert dropped.kind == CODING and dropped.matched_by == "strand_suffix"
    assert dropped.gene.locus == "Rv0002c"

    added = resolver.resolve("Rv0001c")         # annotation says Rv0001
    assert added.gene.locus == "Rv0001"


def test_a_trailing_letter_that_is_not_c_names_a_different_gene(annotation):
    """Rv0063a is not Rv0063. Stripping any trailing letter collides 96 times
    across H37Rv; stripping only `c` collides zero times."""
    assert normalise_locus("Rv0063a") == "Rv0063a"
    assert normalise_locus("Rv0063c") == "Rv0063"

    resolver = Resolver(annotation)
    assert resolver.resolve("Rv0063a").gene.locus == "Rv0063a"
    # ...and the `c` form still reaches Rv0063 without being confused by Rv0063a.
    assert resolver.resolve("Rv0063c").gene.locus == "Rv0063"


def test_an_ambiguous_tag_is_refused_rather_than_guessed(annotation):
    """Rv0070 and Rv0070c both exist, so a label that normalises onto both cannot
    be resolved. A silently wrong locus is worse than a missing one because
    nothing downstream can detect it."""
    feature = Resolver(annotation).resolve("Rv0070C")
    assert feature.kind == UNRESOLVED
    assert set(feature.candidates) == {"Rv0070", "Rv0070c"}
    assert "more than one gene" in feature.note


def test_an_exact_match_wins_over_the_loose_one(annotation):
    """SnpEff falls back to a gene's name where the annotation has one, so a real
    vocabulary mixes `Rv0001` with `rpoC`. Both must resolve, exact first."""
    resolver = Resolver(annotation)
    assert resolver.resolve("Rv0070c").matched_by == "exact"
    assert resolver.resolve("rpoC").matched_by == "symbol"
    assert resolver.resolve("rpoC").gene.locus == "Rv0003"


def test_a_label_in_no_annotation_stays_unresolved(annotation):
    feature = Resolver(annotation).resolve("Rv9999")
    assert feature.kind == UNRESOLVED and feature.gene is None
    assert "no locus tag or symbol" in feature.note


# --- resolving intergenic labels -------------------------------------------


def test_an_intergenic_label_resolves_by_its_flanking_pair(annotation):
    resolver = Resolver(annotation)

    exact = resolver.resolve("Rv0001-Rv0002c")
    assert exact.kind == INTERGENIC and exact.matched_by == "exact"

    # Same gap, both flanks spelled the other way round on the suffix.
    loose = resolver.resolve("Rv0001c-Rv0002")
    assert loose.kind == INTERGENIC and loose.matched_by == "flanking"
    assert loose.interval.name == "Rv0001-Rv0002c"


def test_flanking_genes_that_exist_but_do_not_abut_say_so(annotation):
    """Distinguishable from a missing gene, and a different problem: the caller
    named a gap the annotation does not have because the genes overlap."""
    feature = Resolver(annotation).resolve("Rv0003-Rv0004")
    assert feature.kind == UNRESOLVED
    assert "not consecutive" in feature.note


def test_coverage_counts_every_outcome(annotation):
    coverage = resolve_all(
        ["Rv0001", "Rv0002", "rpoC", "Rv0001-Rv0002c", "Rv0070C", "Rv9999"], annotation)
    counts = coverage.counts()
    assert counts["total"] == 6
    assert counts[CODING] == 3          # exact, strand_suffix, symbol
    assert counts[INTERGENIC] == 1
    assert counts[UNRESOLVED] == 2      # the ambiguous one and the absent one
    assert counts["ambiguous"] == 1
    assert {f.label for f in coverage.unresolved()} == {"Rv0070C", "Rv9999"}


# --- reading the rules -----------------------------------------------------


RULES_TSV = """n_conditions\tconditions\tpredicted_class\tnumerosity\taccuracy\tcoverage\tprecision
2\tRv0001=1 AND Rv0002c=0\tR\t6\t0.8125\t8.44\t0.8125
3\tRv0001=1 AND rpoC=1 AND Rv0005c=0\tS\t11\t0.9524\t31.08\t0.9524
"""


@pytest.fixture
def rules_file(tmp_path: Path) -> Path:
    path = tmp_path / "rules.tsv"
    path.write_text(RULES_TSV)
    return path


def test_conditions_keep_their_state(rules_file):
    """`Rv0002c=0` asserts those isolates match the reference at that locus. Read
    as a bare gene set it becomes "these genes are involved", which is a
    different and wrong claim -- and it is where alternative-route patterns are."""
    rules = parse_rules(rules_file)
    assert [c.state for c in rules[0].conditions] == [1, 0]
    assert rules[0].labels == ("Rv0001", "Rv0002c")
    assert rules[0].k == 2


def test_the_learners_statistics_pass_through_prefixed(rules_file):
    """They are the scan's judgement of the rule, not anything computed here, and
    the prefix is what keeps the two from being confused in a joined table."""
    rules = parse_rules(rules_file)
    assert rules[0].extras["scan_numerosity"] == "6"
    assert rules[0].extras["scan_precision"] == "0.8125"
    assert rules[1].predicted_class == "S"


def test_a_condition_count_that_disagrees_is_reported(tmp_path):
    """The learner says how many conditions it wrote. Disagreement means the
    conjunction was split wrongly, and a rule one condition short still looks
    like a valid rule."""
    path = tmp_path / "bad.tsv"
    path.write_text("n_conditions\tconditions\tpredicted_class\n"
                    "3\tRv0001=1 AND Rv0002c=0\tR\n")
    rule = parse_rules(path)[0]
    assert any("says 3 but 2 condition(s) parsed" in p for p in rule.problems)


def test_an_unparsable_condition_is_reported_not_dropped(tmp_path):
    path = tmp_path / "odd.tsv"
    path.write_text("conditions\tpredicted_class\nRv0001=1 AND Rv0002c>0\tR\n")
    rule = parse_rules(path)[0]
    assert rule.k == 1
    assert any("unparsed condition" in p for p in rule.problems)


def test_rule_identity_ignores_condition_order_and_spelling(rules_file):
    """A conjunction is unordered, so two files listing the same rule differently
    must give the same id. And `Rv0006=1` and `Rv0006c=1` are the same locus under
    two spellings -- the same bug identity.py exists for, one level up."""
    rules = parse_rules(rules_file)
    first = rules[0]

    reordered = type(first)(row=99, conditions=tuple(reversed(first.conditions)),
                            predicted_class=first.predicted_class)
    assert reordered.rule_id() == first.rule_id()

    rename = {"Rv0001": "Rv0001", "Rv0002c": "Rv0002"}
    assert first.rule_id(rename) != first.rule_id()
    respelled = type(first)(row=100, conditions=first.conditions,
                            predicted_class=first.predicted_class)
    assert respelled.rule_id(rename) == first.rule_id(rename)


def test_the_vocabulary_is_what_the_rules_actually_use(rules_file):
    """A separate labels file may list features no rule ever mentions; this is the
    set that has to resolve for a run to mean anything."""
    assert vocabulary(parse_rules(rules_file)) == ["Rv0001", "Rv0002c", "rpoC", "Rv0005c"]


def test_a_file_with_no_conditions_column_is_refused(tmp_path):
    path = tmp_path / "wrong.tsv"
    path.write_text("gene_a\tgene_b\tp\nRv0001\tRv0002c\t1e-8\n")
    with pytest.raises(ValueError, match="no conditions column"):
        parse_rules(path)


# --- GFF3 parent/child rows ------------------------------------------------
#
# A GFF3 carries a parent `gene` row and a child `CDS` row per locus. Counting
# both put every tag in the index twice, which made the strand-suffix rule find
# two candidates that were the same gene and refuse them as ambiguous. On NCBI's
# H37Rv GFF3 that halved resolution: 7,884 "genes", 1,339 false ambiguities,
# `via strand suffix` down from 1,346 to zero.

GFF3 = """##gff-version 3
NC_000962.3\tRefSeq\tgene\t1\t1524\t.\t+\t.\tID=gene-Rv0001;locus_tag=Rv0001;gene=dnaA
NC_000962.3\tRefSeq\tCDS\t1\t1524\t.\t+\t0\tID=cds-Rv0001;locus_tag=Rv0001;gene=dnaA
NC_000962.3\tRefSeq\tgene\t2052\t3260\t.\t+\t.\tID=gene-Rv0003;locus_tag=Rv0003;gene=recF
NC_000962.3\tRefSeq\tCDS\t2052\t3260\t.\t+\t0\tID=cds-Rv0003;locus_tag=Rv0003;gene=recF
NC_000962.3\tRefSeq\tgene\t4000\t5000\t.\t-\t.\tID=gene-Rv0005;locus_tag=Rv0005;gene=gyrB
NC_000962.3\tRefSeq\tCDS\t4000\t5000\t.\t-\t0\tID=cds-Rv0005;locus_tag=Rv0005;gene=gyrB
NC_000962.3\tRefSeq\tpseudogene\t6000\t6500\t.\t+\t.\tID=gene-Rv0007;locus_tag=Rv0007
NC_000962.3\tRefSeq\tgene\t7000\t7300\t.\t+\t.\tID=gene-Rvnr01;locus_tag=Rvnr01;gene=rrs
NC_000962.3\tRefSeq\trRNA\t7000\t7300\t.\t+\t.\tID=rna-Rvnr01;locus_tag=Rvnr01;gene=rrs
"""


@pytest.fixture
def gff3(tmp_path: Path) -> Path:
    path = tmp_path / "h37rv.gff"
    path.write_text(GFF3)
    return path


def test_a_parent_and_child_row_are_one_gene(gff3):
    """Five loci, nine rows. Counting the child rows would give nine 'genes'."""
    annotation = parse_annotation(gff3)
    assert len(annotation.genes) == 5
    assert {g.locus for g in annotation.genes} == {
        "Rv0001", "Rv0003", "Rv0005", "Rv0007", "Rvnr01"}


def test_the_strand_suffix_rule_survives_a_gff3(gff3):
    """The regression itself: `Rv0003c` reaching `Rv0003` found two candidates
    that were the same gene, and was refused."""
    feature = Resolver(parse_annotation(gff3)).resolve("Rv0003c")
    assert feature.kind == CODING
    assert feature.matched_by == "strand_suffix"
    assert feature.gene.locus == "Rv0003"


def test_non_coding_and_pseudogene_loci_come_in_too(gff3):
    """`rrs` and `rrl` are named by real rule vocabularies, and a pseudogene row
    is the only row some loci have. Taking `gene` rows alone would drop the 30
    pseudogenes in H37Rv; taking CDS rows alone would drop every RNA gene."""
    annotation = parse_annotation(gff3)
    assert annotation.gene("Rv0007") is not None                 # pseudogene row
    assert Resolver(annotation).resolve("rrs").gene.locus == "Rvnr01"


def test_a_real_ambiguity_is_still_refused(tmp_path):
    """Deduplication must not weaken the refusal. Two genuinely different tags
    that collide once the suffix is stripped still have no answer."""
    path = tmp_path / "collide.gff"
    path.write_text(
        "##gff-version 3\n"
        "NC_000962.3\tRefSeq\tgene\t100\t200\t.\t+\t.\tID=g1;locus_tag=Rv0070\n"
        "NC_000962.3\tRefSeq\tCDS\t100\t200\t.\t+\t0\tID=c1;locus_tag=Rv0070\n"
        "NC_000962.3\tRefSeq\tgene\t300\t400\t.\t-\t.\tID=g2;locus_tag=Rv0070c\n"
        "NC_000962.3\tRefSeq\tCDS\t300\t400\t.\t-\t0\tID=c2;locus_tag=Rv0070c\n")
    feature = Resolver(parse_annotation(path)).resolve("Rv0070C")
    assert feature.kind == UNRESOLVED
    assert set(feature.candidates) == {"Rv0070", "Rv0070c"}


def test_a_cds_only_file_falls_back_and_says_so(tmp_path):
    """Prokka and some Ensembl bacterial dumps emit no gene rows. The file stays
    usable, and the note records that spans are translated regions."""
    path = tmp_path / "prokka.gff"
    path.write_text(
        "##gff-version 3\n"
        "NC_000962.3\tProkka\tCDS\t100\t200\t.\t+\t0\tID=c1;locus_tag=Rv0001\n"
        "NC_000962.3\tProkka\tCDS\t300\t400\t.\t+\t0\tID=c2;locus_tag=Rv0003\n")
    annotation = parse_annotation(path)
    assert len(annotation.genes) == 2
    assert any("no gene, pseudogene rows found" in n for n in annotation.notes)
    assert any("read CDS rows instead" in n for n in annotation.notes)


def test_c_coordinates_count_from_the_translation_start(tmp_path):
    """HGVS `c.-N` is relative to the ATG, not to the gene's 5' end. Three H37Rv
    gene rows begin before their CDS -- Rv0614 by 243 bp -- so counting from the
    gene start puts a promoter variant that far from where it is."""
    path = tmp_path / "leader.gff"
    path.write_text(
        "##gff-version 3\n"
        "NC_000962.3\tRefSeq\tgene\t709356\t710348\t.\t+\t.\tID=g;locus_tag=Rv0614\n"
        "NC_000962.3\tRefSeq\tCDS\t709599\t710348\t.\t+\t0\tID=c;locus_tag=Rv0614\n")
    gene = parse_annotation(path).gene("Rv0614")
    assert (gene.start, gene.end) == (709356, 710348)     # the locus's extent
    assert gene.cds_start == 709599                       # where `c.` counts from
    assert gene.coding_start == 709599


def test_a_gene_whose_cds_matches_records_no_separate_start(gff3):
    """The common case. Recording a redundant CDS start would make the note fire
    on every file."""
    assert parse_annotation(gff3).gene("Rv0001").cds_start is None
    assert parse_annotation(gff3).gene("Rv0001").coding_start == 1


# --- prefixed region names -------------------------------------------------
#
# A producer that collapses SnpEff ANN output names non-coding features
# `upstream_<gene>`, `downstream_<gene>` and `intergenic_<a>-<b>`, with a
# transcript version on the gene part. None of those is a flanking pair, so the
# whole intergenic half of such a vocabulary resolved at zero.


def test_a_transcript_version_is_stripped_from_the_gene_part(annotation):
    """`upstream_Rv0003.1` is SnpEff's naming. Locus tags carry no dot, so a
    trailing `.<digits>` is a version and nothing else."""
    from kegg_string_mcp.rules.annotation import strip_version

    assert strip_version("Rv1482c.1") == "Rv1482c"
    assert strip_version("Rv1482c") == "Rv1482c"

    feature = Resolver(annotation).resolve("upstream_Rv0003.1")
    assert feature.kind == INTERGENIC
    assert feature.interval.name == "Rv0002c-Rv0003"


def test_upstream_is_read_on_the_genes_own_strand(annotation):
    """5' is the lower coordinate on the forward strand and the higher one on the
    reverse. Reading it the same way for both puts a promoter behind the gene."""
    resolver = Resolver(annotation)

    forward = resolver.resolve("upstream_Rv0003")       # + strand
    assert forward.matched_by == "upstream"
    assert (forward.interval.start, forward.interval.end) == (400, 449)

    reverse = resolver.resolve("upstream_Rv0002c")      # - strand, diverging
    assert (reverse.interval.start, reverse.interval.end) == (400, 449)

    downstream = resolver.resolve("downstream_Rv0001")  # + strand
    assert downstream.matched_by == "downstream"
    assert downstream.interval.name == "Rv0001-Rv0002c"


def test_a_prefixed_region_says_the_window_is_wider_than_the_interval(annotation):
    """SnpEff calls an upstream variant within a window -- 5,000 bp by default --
    and 3,048 of H37Rv's 3,049 intergenic intervals are narrower than that. The
    coordinates are the intergenic part of the region, not its extent."""
    feature = Resolver(annotation).resolve("upstream_Rv0003")
    assert "wider than this interval" in feature.note


def test_a_gene_with_no_gap_upstream_gets_its_own_answer(annotation):
    """830 of H37Rv's 4,008 genes abut or overlap their 5' neighbour, so a variant
    called upstream of them lies inside that neighbour. A fact about the genome,
    and a different answer from 'not found'."""
    feature = Resolver(annotation).resolve("upstream_Rv0004")   # overlaps Rv0003
    assert feature.kind == UNRESOLVED
    assert "no intergenic interval 5' of it" in feature.note
    # The bounding gene is named -- by symbol where it has one -- rather than
    # left as "its neighbour", which is what a five-cohort run showed it saying.
    assert "rpoC (Rv0003) abuts or overlaps it" in feature.note
    assert feature.bounded_by == "Rv0003"


def test_an_intergenic_prefix_is_still_a_flanking_pair(annotation):
    resolver = Resolver(annotation)
    assert resolver.resolve("intergenic_Rv0001-Rv0002c").interval.name == "Rv0001-Rv0002c"
    # Versions on both halves.
    assert resolver.resolve("intergenic_Rv0001.1-Rv0002c.1").interval.name == "Rv0001-Rv0002c"


def test_an_unrecognised_prefix_names_what_was_tried(annotation):
    """A different caller will use different words. The refusal should make that
    a configuration change rather than someone else's bug."""
    feature = Resolver(annotation).resolve("sideways_Rv0001")
    assert feature.kind == UNRESOLVED
    assert "tried: downstream_, intergenic_, upstream_" in feature.note


def test_the_prefix_set_is_the_callers_to_override(annotation):
    """SnpEff's three are the default, not a law."""
    resolver = Resolver(annotation, prefixes={"5prime_": "upstream"})
    assert resolver.resolve("5prime_Rv0003").interval.name == "Rv0002c-Rv0003"
    # ...and the defaults are then not in play, which the refusal says.
    assert "tried: 5prime_" in resolver.resolve("upstream_Rv0003").note


def test_a_promoter_by_orientation_is_distinguishable_from_a_stated_pair(annotation):
    """One is asserted by the caller, the other inferred by us from strand. A
    table that shows them the same way hides which assumption a row rests on."""
    resolver = Resolver(annotation)
    assert resolver.resolve("Rv0001-Rv0002c").matched_by == "exact"
    assert resolver.resolve("intergenic_Rv0001c-Rv0002").matched_by == "flanking"
    assert resolver.resolve("upstream_Rv0003").matched_by == "upstream"


def test_a_symbol_that_ends_in_a_number_is_not_a_version(annotation):
    """H37Rv has ten of these -- TB7.3, TB15.3, TB31.7 and the rest of that
    family -- where the trailing number is part of the gene's name. Stripping
    first maps TB7.3 to TB7, which is nothing. Same shape as the strand-suffix
    rule: `.1` from SnpEff is a version, `.3` in TB7.3 is the gene."""
    resolver = Resolver(annotation)

    feature = resolver.resolve("TB7.3")
    assert feature.kind == CODING and feature.matched_by == "symbol"
    assert feature.gene.locus == "Rv0080"
    assert "version was stripped" not in feature.note

    # The stem alone is not a gene, and must not become one.
    assert resolver.resolve("TB7").kind == UNRESOLVED


def test_a_transcript_version_still_strips_when_nothing_else_matches(annotation):
    """SnpEff's `Rv0001.1` must keep working. The fallback is only reached once
    the name as given has failed."""
    feature = Resolver(annotation).resolve("Rv0001.1")
    assert feature.gene.locus == "Rv0001"
    assert "version was stripped" in feature.note


def test_a_label_matching_neither_reading_reports_both_attempts(annotation):
    """So the refusal says what was tried rather than leaving the reader to guess
    which reading failed."""
    feature = Resolver(annotation).resolve("Rv9999.1")
    assert feature.kind == UNRESOLVED
    assert "also tried without the trailing version, as 'Rv9999'" in feature.note


def test_a_flanking_pair_tries_the_names_as_given_first(annotation):
    """Same rule for interval names: the pair as written, then the stripped
    stems, so a versioned pair and an unversioned one both resolve."""
    resolver = Resolver(annotation)
    assert resolver.resolve("Rv0001-Rv0002c").interval.name == "Rv0001-Rv0002c"
    assert resolver.resolve("intergenic_Rv0001.1-Rv0002c.1").interval.name == "Rv0001-Rv0002c"


# --- Mycobrowser-shaped GFF ------------------------------------------------
#
# The TubercuList successor, and the fuller H37Rv annotation: 4,173 loci against
# NCBI's 4,008. Structurally different in two ways that each broke parsing.

MYCOBROWSER = """NC_000962.3\tMycobrowser_v5\tCDS\t100\t200\t.\t+\t\tLocus=Rv0001;Name=dnaA;Product=X
NC_000962.3\tMycobrowser_v5\tCDS\t300\t400\t.\t-\t\tLocus=Rv0002c;Name=dnaN;Product=Y
NC_000962.3\tMycobrowser_v5\tncRNA\t500\t560\t.\t+\t\tLocus=MTB000115;Name=ncRv10666;Product=Z
NC_000962.3\tMycobrowser_v5\trRNA\t700\t800\t.\t+\t\tLocus=MTB000019;Name=rrs;Product=16S
NC_000962.3\tMycobrowser_v5\ttRNA\t900\t980\t.\t+\t\tLocus=MTB000001;Name=alaT;Product=tRNA
NC_000962.3\tMycobrowser_v5\tpromoter\t210\t240\t.\t+\t\tLocus=MTBp0001;Name=Prv0001
"""


@pytest.fixture
def mycobrowser(tmp_path: Path) -> Path:
    path = tmp_path / "H37Rv.gff"
    path.write_text(MYCOBROWSER)
    return path


def test_the_locus_attribute_is_read_before_the_name(mycobrowser):
    """Mycobrowser writes the tag in `Locus` and the SYMBOL in `Name`. Reading
    `Name` first makes every locus its own gene symbol -- dnaA rather than
    Rv0001 -- and no locus tag in a vocabulary resolves."""
    annotation = parse_annotation(mycobrowser)
    gene = annotation.gene("Rv0001")
    assert gene is not None and gene.symbol == "dnaA"
    assert Resolver(annotation).resolve("dnaA").gene.locus == "Rv0001"


def test_rna_loci_survive_a_file_with_no_gene_rows(mycobrowser):
    """Mycobrowser emits no `gene` row anywhere, so the fallback decides what
    comes in. CDS alone would drop 141 real loci -- every RNA gene, rrs and rrl
    among them."""
    annotation = parse_annotation(mycobrowser)
    assert {g.locus for g in annotation.genes} == {
        "Rv0001", "Rv0002c", "MTB000115", "MTB000019", "MTB000001"}
    assert Resolver(annotation).resolve("rrs").gene.locus == "MTB000019"
    assert any("no gene, pseudogene rows found" in n for n in annotation.notes)


def test_a_regulatory_row_is_not_a_gene(mycobrowser):
    """`promoter`, `-35_signal` and `-10_signal` rows carry a Locus value but
    name part of a gene, not a gene. Admitting them would put spurious loci
    between real ones and change every intergenic interval around them."""
    assert parse_annotation(mycobrowser).gene("MTBp0001") is None


def test_a_tag_whose_locus_was_split_is_refused(tmp_path):
    """`Rv2306c` normalises to `Rv2306`, which exists in neither annotation --
    the locus was split into Rv2306A and Rv2306B between revisions. The trailing
    letter is part of the name, so neither half is reachable by suffix
    stripping, and guessing one would silently pick half a gene."""
    path = tmp_path / "split.gff"
    path.write_text(
        "NC_000962.3\tMycobrowser_v5\tCDS\t100\t200\t.\t+\t\tLocus=Rv2306A;Name=a\n"
        "NC_000962.3\tMycobrowser_v5\tCDS\t300\t400\t.\t+\t\tLocus=Rv2306B;Name=b\n")
    feature = Resolver(parse_annotation(path)).resolve("Rv2306c")
    assert feature.kind == UNRESOLVED
    assert feature.candidates == ()


def test_curated_fields_are_read_from_the_annotation_that_has_them(mycobrowser):
    """A Mycobrowser GFF already carries these in its attributes, so reading them
    costs nothing and needs no second file, no new tool and no fetch."""
    annotation = parse_annotation(mycobrowser)
    gene = annotation.gene("Rv0001")
    assert gene.product == "X"
    # A BED has no attributes, so the fields stay empty -- and empty means "this
    # annotation did not say", never "uncategorised".
    assert parse_annotation(mycobrowser).gene("Rv0002c").category == ""


def test_repetitive_categories_are_flagged_for_variant_calling(tmp_path):
    """PE/PPE (168 loci) and mobile elements (147) are repeated elsewhere in the
    genome, so short reads misplace and a condition on one may be an artefact."""
    from kegg_string_mcp.rules.annotation import MOBILE, PE_PPE

    path = tmp_path / "cats.gff"
    path.write_text(
        f"NC_000962.3\tMB\tCDS\t100\t200\t.\t+\t\tLocus=Rv1806;Name=PE20;"
        f"Functional_Category={PE_PPE}\n"
        f"NC_000962.3\tMB\tCDS\t300\t400\t.\t+\t\tLocus=Rv0031;Functional_Category={MOBILE}\n"
        f"NC_000962.3\tMB\tCDS\t500\t600\t.\t+\t\tLocus=Rv0002;"
        f"Functional_Category=lipid metabolism\n")
    annotation = parse_annotation(path)
    assert annotation.gene("Rv1806").repetitive == "repetitive and GC-rich"
    assert annotation.gene("Rv0031").repetitive == "a multi-copy mobile element"
    # A category that is not a mapping hazard gets no caution.
    assert annotation.gene("Rv0002").repetitive == ""
    assert annotation.gene("Rv0002").category == "lipid metabolism"
    assert any("are repetitive" in n for n in annotation.notes)


def test_the_rules_package_never_reaches_a_model():
    """The boundary this branch is for.

    Everything in `rules/` is a lookup, a count or a comparison, so a run needs
    no API key, costs nothing, and gives the same answer twice. That is what lets
    it run on a cluster with no credentials and lets a reader check any
    classification against the locus row it rests on.

    Hypothesis generation and quote extraction need a model and belong elsewhere.
    Stated as a test because a single convenient import would end the property
    silently, and nothing else would fail.
    """
    import ast

    package = Path(__file__).resolve().parents[1] / "src" / "kegg_string_mcp" / "rules"
    forbidden = ("anthropic", "kegg_string_mcp.agent.loop", "kegg_string_mcp.agent.pipeline",
                 "kegg_string_mcp.agent.modes")
    offenders = []
    for source in sorted(package.glob("*.py")):
        tree = ast.parse(source.read_text(encoding="utf-8"))
        for node in ast.walk(tree):
            names = ([a.name for a in node.names] if isinstance(node, ast.Import)
                     else [node.module or ""] if isinstance(node, ast.ImportFrom) else [])
            for name in names:
                if any(name == f or name.startswith(f + ".") for f in forbidden):
                    offenders.append(f"{source.name} imports {name}")
    assert not offenders, "; ".join(offenders)
