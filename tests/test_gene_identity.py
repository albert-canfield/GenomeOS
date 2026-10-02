# SPDX-License-Identifier: AGPL-3.0-or-later
"""The gene-identity check: its definitions, its exhaustiveness, and the resolver it reproduces."""

from __future__ import annotations

import gzip
import json
from collections import Counter
from pathlib import Path

import pytest

from genomeos.attribution import gene_identity as gi

REFERENCE = Path("data/reference")


def row(gid: str, symbol: str, gtype: str, start: int, end: int, tss: int | None = None) -> gi.GeneRow:
    return gi.GeneRow(gid, symbol, gtype, start, end, start if tss is None else tss)


def rule(target: str, start: int = 1_000_000, end: int = 1_000_200, element: str = "EH38E1") -> gi.Rule:
    return gi.Rule(
        element=element,
        chrom="chrT",
        start=start,
        end=end,
        verb="activates",
        target=target,
        cell="K562",
        strength=0.4,
        activity="activates_target",
    )


# ---- the definitions are fixed, named and carry no verdict --------------------------------------


def test_every_outcome_is_named_and_the_three_classes_are_disjoint() -> None:
    assert gi.AGREE not in gi.CAUSES and gi.DIFFER not in gi.CAUSES
    assert len(gi.OUTCOMES) == len(set(gi.OUTCOMES)) == 2 + len(gi.CAUSES)
    for name, why in gi.CAUSES.items():
        assert why and not name.startswith("other") and "residue" not in name


def test_no_outcome_name_carries_a_verdict_stem() -> None:
    """A cause is a statement about the record. It may not read as a judgement on the rule."""
    forbidden = (
        "support",
        "contradict",
        "validat",
        "confirm",
        "refut",
        "verif",
        "correct",
        "wrong",
        "true",
        "false",
        "good",
        "bad",
        "pass",
        "fail",
    )
    for name in gi.OUTCOMES:
        assert not any(s in name for s in forbidden), name


def test_the_limitations_refuse_the_reading_the_count_must_not_be_given() -> None:
    text = " ".join(gi.limitations().values())
    assert "not a measure of any feature's contribution" in text
    assert "establishes about no single rule that the rule is wrong" in text
    assert "recovered from the element's recorded window gene list" in text


def test_the_imported_reading_is_carried_word_for_word() -> None:
    assert gi.imported_readings()["lane_notopen_1678"] == (
        "1,678 rules carry a TSS further from the element than half the 1 Mb deletion window, so "
        "a repeated gene symbol resolved to another locus - a fault in the distance, not the rule"
    )


def test_the_half_window_is_imported_and_not_restated() -> None:
    from genomeos.attribution.not_open_profile import SCORER_HALF_WINDOW

    assert gi.HALF_WINDOW == SCORER_HALF_WINDOW == 500_000
    assert gi.SENSITIVITY_HALF_WINDOW == 524_288


# ---- the compiled side reproduces the project's own symbol resolution ---------------------------


def test_compiled_side_prefers_a_protein_coding_locus_then_the_lowest_tss() -> None:
    a = row("ENSG1", "DUP", "lncRNA", 100, 200, tss=100)
    b = row("ENSG2", "DUP", "protein_coding", 900, 1000, tss=900)
    c = row("ENSG3", "DUP", "protein_coding", 300, 400, tss=300)
    rows = (a, b, c)
    assert gi.compiled_locus("DUP", gi.by_symbol(rows), gi.by_id(rows)) is c


def test_compiled_side_resolves_a_bare_ensembl_id_directly() -> None:
    a = row("ENSG9", "", "lncRNA", 10, 20)
    rows = (a,)
    out = gi.compiled_locus("ENSG9", gi.by_symbol(rows), gi.by_id(rows))
    assert out is a


def test_compiled_side_returns_none_for_a_token_the_chromosome_does_not_carry() -> None:
    rows = (row("ENSG1", "A", "protein_coding", 10, 20),)
    assert gi.compiled_locus("B", gi.by_symbol(rows), gi.by_id(rows)) is None


@pytest.mark.skipif(
    not (REFERENCE / "gencode_v50_chr21.gff3.gz").exists(), reason="GENCODE v50 chr21 not on disk"
)
def test_the_compiled_resolver_agrees_with_pilot_bio_gene_tss_on_chr21() -> None:
    """The resolver lane-notopen's distance came from. This module may carry the id, never choose
    a different locus from the one the project already chooses."""
    from genomeos.attribution.pilot_bio import gene_tss

    rows = gi.annotation("chr21")
    symbols, ids = gi.by_symbol(rows), gi.by_id(rows)
    theirs = gene_tss("chr21")
    assert theirs
    for symbol, tss in theirs.items():
        mine = gi.compiled_locus(symbol, symbols, ids)
        assert mine is not None, symbol
        assert mine.tss == tss, symbol


# ---- the measured side ---------------------------------------------------------------------------


def test_window_candidates_are_the_genes_whose_body_overlaps_the_window() -> None:
    inside = row("ENSG1", "G", "protein_coding", 1_400_000, 1_400_500)
    body_reaches_in = row("ENSG2", "G", "protein_coding", 400_000, 600_000)
    outside = row("ENSG3", "G", "protein_coding", 3_000_000, 3_000_500)
    rows = (inside, body_reaches_in, outside)
    cand = gi.window_candidates("G", 1_000_100, gi.by_symbol(rows))
    assert set(c.gene_id for c in cand) == {"ENSG1", "ENSG2"}


def test_a_name_the_scorer_never_listed_is_unresolvable_and_not_a_guess() -> None:
    rows = (row("ENSG1", "G", "protein_coding", 1_000_000, 1_000_500),)
    out = gi.classify(rule("G"), gi.by_symbol(rows), gi.by_id(rows), occurrences=0)
    assert out.outcome == "target_not_in_the_recorded_window_list"
    assert out.measured_id is None and not out.resolved


def test_a_name_the_scorer_listed_twice_is_unresolvable() -> None:
    rows = (row("ENSG1", "G", "protein_coding", 1_000_000, 1_000_500),)
    out = gi.classify(rule("G"), gi.by_symbol(rows), gi.by_id(rows), occurrences=2)
    assert out.outcome == "target_named_more_than_once_in_the_recorded_window_list"


def test_two_annotated_loci_in_the_window_are_unresolvable_and_never_picked_between() -> None:
    rows = (
        row("ENSG1", "G", "protein_coding", 900_000, 900_500),
        row("ENSG2", "G", "lncRNA", 1_100_000, 1_100_500),
    )
    out = gi.classify(rule("G"), gi.by_symbol(rows), gi.by_id(rows), occurrences=1)
    assert out.outcome == "two_or_more_annotated_loci_of_that_name_in_the_scorer_window"
    assert out.measured_id is None


def test_the_same_locus_on_both_sides_agrees() -> None:
    rows = (row("ENSG1", "G", "protein_coding", 1_000_000, 1_000_500),)
    out = gi.classify(rule("G"), gi.by_symbol(rows), gi.by_id(rows), occurrences=1)
    assert out.outcome == gi.AGREE and out.compiled_id == out.measured_id == "ENSG1"


def test_a_repeated_symbol_resolved_to_a_locus_outside_the_window_differs() -> None:
    """The case the lane exists for: the project's resolver takes the first protein-coding locus of
    the symbol on the chromosome, which is not the one the scorer could have read."""
    far = row("ENSG_FAR", "G", "protein_coding", 10_000, 10_500)
    near = row("ENSG_NEAR", "G", "protein_coding", 1_020_000, 1_020_500)
    rows = (far, near)
    out = gi.classify(rule("G"), gi.by_symbol(rows), gi.by_id(rows), occurrences=1)
    assert out.outcome == gi.DIFFER
    assert out.compiled_id == "ENSG_FAR" and out.measured_id == "ENSG_NEAR"


def test_an_experimental_layer_rule_has_no_deletion_answer_to_compare() -> None:
    rows = (row("ENSG1", "G", "protein_coding", 1_000_000, 1_000_500),)
    r = rule("G", element="EH38E1_measured")
    assert r.source == gi.SOURCE_MEASURED
    out = gi.classify(r, gi.by_symbol(rows), gi.by_id(rows), occurrences=None)
    assert out.outcome == "rule_of_the_experimental_layer"


def test_an_element_with_no_window_record_is_its_own_cause() -> None:
    rows = (row("ENSG1", "G", "protein_coding", 1_000_000, 1_000_500),)
    out = gi.classify(rule("G"), gi.by_symbol(rows), gi.by_id(rows), occurrences=None)
    assert out.outcome == "element_absent_from_the_window_cache"


def test_a_token_absent_from_the_annotation_resolves_on_neither_side() -> None:
    rows = (row("ENSG1", "A", "protein_coding", 10, 20),)
    out = gi.classify(rule("B"), gi.by_symbol(rows), gi.by_id(rows), occurrences=1)
    assert out.outcome == "target_absent_from_the_chromosome_annotation"
    assert out.compiled_id is None and out.measured_id is None


def test_every_rule_lands_in_exactly_one_outcome_whatever_it_is_given() -> None:
    rows = (
        row("ENSG1", "G", "protein_coding", 1_000_000, 1_000_500),
        row("ENSG2", "G", "protein_coding", 10_000, 10_500),
        row("ENSG3", "", "lncRNA", 1_000_000, 1_000_400),
    )
    symbols, ids = gi.by_symbol(rows), gi.by_id(rows)
    tally: Counter[str] = Counter()
    for target in ("G", "ENSG3", "MISSING"):
        for occ in (None, 0, 1, 2):
            for element in ("EH38E1", "EH38E1_measured"):
                out = gi.classify(rule(target, element=element), symbols, ids, occ)
                assert out.outcome in gi.OUTCOMES
                tally[out.outcome] += 1
    assert sum(tally.values()) == 3 * 4 * 2


# ---- the compiled programs are read as they are written -----------------------------------------


def test_the_program_parser_reads_a_rule_with_its_element_block(tmp_path: Path) -> None:
    program = tmp_path / "noncoding_chrT.bio"
    program.write_text(
        "module human.noncoding.chrT\n"
        "# a comment\n"
        "element EH38E9 {\n"
        "  class: enhancer\n"
        "  locus: chrT:100-300\n"
        "  activity: represses_target\n"
        "}\n"
        "rule EH38E9 inhibits SOMEGENE { strength: 0.42; when: cell_type = K562; evidence: x }\n"
    )
    got = list(gi.rules("chrT", compiled=tmp_path))
    assert len(got) == 1
    r = got[0]
    assert (r.element, r.chrom, r.start, r.end) == ("EH38E9", "chrT", 100, 300)
    assert (r.verb, r.target, r.cell) == ("inhibits", "SOMEGENE", "K562")
    assert r.strength == 0.42 and r.effect_band == "strong"
    assert r.activity == "represses_target" and r.midpoint == 200


def test_the_parser_refuses_a_rule_it_cannot_locate(tmp_path: Path) -> None:
    program = tmp_path / "noncoding_chrT.bio"
    program.write_text("rule EH38E9 activates G { strength: 0.1; when: cell_type = K562 }\n")
    with pytest.raises(ValueError, match="no element block"):
        list(gi.rules("chrT", compiled=tmp_path))


def test_the_effect_band_is_the_rule_lines_own_strength() -> None:
    from genomeos.predict.enhancer_target import STRONG_EFFECT

    assert rule("G").effect_band == "strong"
    weak = gi.Rule("E", "chrT", 0, 2, "activates", "G", "K562", STRONG_EFFECT - 0.01, "activates_target")
    assert weak.effect_band == "weak"
    none = gi.Rule("E", "chrT", 0, 2, "activates", "G", "K562", None, "activates_target")
    assert none.effect_band is None


# ---- the window cache is streamed, and a missing record is not a zero ----------------------------


def test_window_gene_counts_streams_and_separates_absent_from_zero(tmp_path: Path) -> None:
    payload = {
        "EH38E1": {
            "id": "EH38E1",
            "genes_in_window": 3,
            "genes": [{"gene": "A"}, {"gene": "B"}, {"gene": "A"}],
        },
        "EH38E2": {"id": "EH38E2", "genes_in_window": 1, "genes": [{"gene": "C"}]},
    }
    path = tmp_path / "chrT.json.gz"
    with gzip.open(path, "wt") as fh:
        json.dump(payload, fh)
    wanted = {"EH38E1": frozenset({"A", "Z"}), "EH38E2": frozenset({"A"})}
    counts, seen, present = gi.window_gene_counts("chrT", wanted, cache=tmp_path, chunk=16)
    assert seen == 2
    assert present == {"EH38E1", "EH38E2"}
    assert counts["EH38E1"]["A"] == 2
    assert counts["EH38E1"].get("Z", 0) == 0
    assert counts.get("EH38E2", Counter()).get("A", 0) == 0


def test_window_gene_counts_on_a_chromosome_with_no_cache(tmp_path: Path) -> None:
    counts, seen, present = gi.window_gene_counts("chrNope", {}, cache=tmp_path)
    assert counts == {} and seen == 0 and present == set()


# ---- the bands are the project's own --------------------------------------------------------------


def test_the_distance_bands_are_the_projects_with_one_residual_row() -> None:
    from genomeos.attribution import executor as ex

    assert gi.distance_band(None) == "no distance"
    assert gi.distance_band(1) == ex._band(1)
    assert gi.distance_band(ex.WIDE_MAX_DISTANCE) == ex._band(ex.WIDE_MAX_DISTANCE)
    assert gi.distance_band(ex.WIDE_MAX_DISTANCE + 1) == "beyond 250 kb"


# ---- the registration is on disk before the result, and the result sums --------------------------

RESULTS = Path("data/results")


@pytest.mark.skipif(
    not (RESULTS / "gene_identity_registration.json").exists(), reason="registration not written yet"
)
def test_the_registration_fixes_the_definitions_before_any_count() -> None:
    reg = json.loads((RESULTS / "gene_identity_registration.json").read_text())
    assert reg["alphagenome_requests"] == 0 and reg["network_requests"] == 0
    assert reg["mis_resolution"] == gi.MIS_RESOLUTION
    # The registration is never edited. Amendment 1 added one cause, so the registration
    # names a subset and the amendment carries the rest.
    assert set(reg["unresolvable_symbols"]["causes"]) <= set(gi.CAUSES)
    assert reg["population"]["expected_total"] == 440589
    assert reg["sensitivity_fixed_in_advance"]["half_window"] == gi.HALF_WINDOW
    assert "before any count over the compiled programs was read" in reg["status"]
    annotation = reg["inputs_named_by_sha256"]["annotation"]
    assert annotation and all(len(e["sha256"]) == 64 for e in annotation)


@pytest.mark.skipif(not (RESULTS / "gene_identity.json").exists(), reason="result not written yet")
def test_the_result_sums_to_its_population_with_every_cause_named() -> None:
    res = json.loads((RESULTS / "gene_identity.json").read_text())
    c = res["counts"]
    assert c["sums_to_the_population"] is True
    assert c["differ"] + c["agree"] + c["unresolvable"] == res["population"]["rules"] == 440589
    assert set(c["unresolvable_by_cause"]) <= set(gi.CAUSES)
    assert sum(c["unresolvable_by_cause"].values()) == c["unresolvable"]
    assert res["population"]["reproduces_the_fixed_denominator"] is True
    assert res["alphagenome_requests"] == 0 and res["network_requests"] == 0
    assert res["the_1678"]["imported_reading"] == gi.imported_readings()["lane_notopen_1678"]
    assert "is not equated with any count here" in res["the_1678"]["the_two_populations"]
    assert "establishes about no single rule" in " ".join(res["what_this_cannot_establish"].values())


# ---- amendment 1: the compiled target is a BioLang identifier, not a gene symbol -----------------


def test_a_gene_is_found_under_its_name_and_under_its_ident() -> None:
    """GENCODE's KRTAP10-1 is written KRTAP10_1 by compile.ident, and 152 chr21 rules carry such a
    token. Both spellings must find the one gene."""
    g = row("ENSG1", "KRTAP10-1", "protein_coding", 1_000_000, 1_000_500)
    symbols = gi.by_symbol((g,))
    assert gi.compiled_locus("KRTAP10-1", symbols, gi.by_id((g,))) is g
    assert gi.compiled_locus("KRTAP10_1", symbols, gi.by_id((g,))) is g


def test_a_token_matching_two_distinct_symbols_is_refused_and_not_picked_between() -> None:
    a = row("ENSG1", "A-1", "protein_coding", 1_000_000, 1_000_500)
    b = row("ENSG2", "A_1", "protein_coding", 1_002_000, 1_002_500)
    symbols, ids = gi.by_symbol((a, b)), gi.by_id((a, b))
    assert gi.symbols_of("A_1", symbols) == frozenset({"A-1", "A_1"})
    assert gi.compiled_locus("A_1", symbols, ids) is None
    out = gi.classify(rule("A_1"), symbols, ids, occurrences=1)
    assert out.outcome == "the_token_matches_two_or_more_symbols_in_the_annotation"
    assert out.compiled_id is None and out.measured_id is None


def test_the_ident_fallback_does_not_disturb_a_token_with_one_symbol() -> None:
    a = row("ENSG1", "A-1", "protein_coding", 1_000_000, 1_000_500)
    symbols, ids = gi.by_symbol((a,)), gi.by_id((a,))
    assert gi.symbols_of("A_1", symbols) == frozenset({"A-1"})
    out = gi.classify(rule("A_1"), symbols, ids, occurrences=1)
    assert out.outcome == gi.AGREE and out.compiled_id == "ENSG1"


def test_the_recorded_names_to_look_for_are_the_real_symbols_and_the_token() -> None:
    """The cache records AlphaGenome's names, which carry the hyphen the compiled token has lost."""
    a = row("ENSG1", "A-1", "protein_coding", 1_000_000, 1_000_500)
    symbols, ids = gi.by_symbol((a,)), gi.by_id((a,))
    assert gi.recorded_names("A_1", symbols, ids) == frozenset({"A_1", "A-1"})
    assert gi.recorded_names("ENSG9", {}, {}) == frozenset({"ENSG9"})


AMENDMENT = RESULTS / "gene_identity_registration_amendment_1.json"


@pytest.mark.skipif(not AMENDMENT.exists(), reason="amendment 1 not written yet")
def test_amendment_one_accounts_for_every_cause_the_registration_does_not_name() -> None:
    reg = json.loads((RESULTS / "gene_identity_registration.json").read_text())
    am = json.loads(AMENDMENT.read_text())
    assert am["amends"] == "gene_identity_registration"
    assert set(reg["unresolvable_symbols"]["causes"]) | set(am["causes_added"]) == set(gi.CAUSES)
    assert am["counts_already_read_when_this_was_written"]
    assert am["definitions_that_did_not_move"]
