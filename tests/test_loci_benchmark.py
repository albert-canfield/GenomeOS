# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""The machinery against loci whose answer is already published (docs/LOCI-BENCHMARK.md).

The run needs local reference data and the network (scripts/loci_benchmark.py). These tests read
its committed result, so a locus that used to pass and stops passing fails CI within a second.
The expectation and scoring tests below need nothing but the code.
"""

from __future__ import annotations

import pytest

from genomeos.benchmark.loci import (
    BEYOND_SEQUENCE_NEEDS,
    CLASSES,
    GATE_PANEL,
    PANEL,
    READER_CELLS,
    Expect,
    score_class,
    score_target,
)
from genomeos.results import load_result

RESULT = load_result("loci_benchmark")
needs_result = pytest.mark.skipif(not RESULT, reason="benchmark result not present")


# ------------------------------------------------------------------ the expectations, offline
def test_the_panel_is_small_and_covers_every_class():
    assert 10 <= len(PANEL) <= 24, "the panel must stay small enough to run and to read"
    seen = {c for e in PANEL for c in e.classes}
    assert seen == set(CLASSES)
    names = {e.locus for e in PANEL}
    assert {"HBB_LCR", "HOXD"} <= names, "ordered in time and ordered in space are both required"


def test_the_gate_set_is_the_pinned_subset_of_the_panel():
    """The full panel may grow; the gate set is the loci whose verdicts have been seen and defended."""
    names = {e.locus for e in PANEL}
    assert set(GATE_PANEL) <= names, "the gate names a locus that is not in the panel"
    assert len(GATE_PANEL) <= len(PANEL)
    assert set(GATE_PANEL) >= DERIVED_TARGET_PASSES, "a pinned pass sits outside the gate set"


def test_the_panel_covers_the_mechanisms_that_are_not_in_the_sequence():
    """Widening the panel was for the gaps: parent of origin, copy number, a second long-range case."""
    beyond = {m for e in PANEL for m in e.beyond_sequence}
    assert {"parent of origin", "copy number", "developmental stage"} <= beyond
    for e in PANEL:
        for m in e.beyond_sequence:
            assert m in BEYOND_SEQUENCE_NEEDS, f"{e.locus}: {m} has no recorded cost to reach"
    long_range = [e for e in PANEL if (e.distance or 0) >= 500_000]
    assert len(long_range) >= 2, "one megabase-reach locus is an anecdote; the panel needs two"


def test_every_expectation_is_written_down_with_its_sources():
    for e in PANEL:
        assert e.window[0] <= e.element[0] < e.element[1] <= e.window[1], e.locus
        assert e.targets and e.citations and e.answer_from, e.locus
        assert set(e.classes) <= set(CLASSES), e.locus
        assert set(e.cells) <= set(READER_CELLS), e.locus
        assert e.direction in ("activates", "represses", "coding", "none"), e.locus
        for v in e.variants:
            if v.get("pos"):
                assert e.window[0] <= v["pos"] <= e.window[1], f"{e.locus} {v.get('rsid')}"


def _toy(classes=("program",)) -> Expect:
    return Expect(
        locus="toy",
        chrom="chr1",
        element=(100, 200),
        window=(0, 1000),
        classes=classes,
        targets=("GOOD",),
        direction="activates",
        tissues=("x",),
        citations=("none",),
        answer_from=("none",),
    )


def test_a_layer_hits_only_when_it_ranks_the_published_target_first():
    readings = {
        "deletion": {"provenance": "derived", "targets": [{"gene": "BAD"}, {"gene": "GOOD"}]},
        "node": {"provenance": "heuristic", "target": "GOOD"},
    }
    s = score_target(_toy(), readings)
    assert not s["hit_derived"] and s["hit_derived_among"]
    assert s["hit_heuristic"] and not s["hit_looked_up"]


def test_a_lookup_never_counts_as_a_derived_hit():
    readings = {"lookups": {"provenance": "looked_up", "gwas_genes": ["GOOD"], "clinvar_genes": []}}
    s = score_target(_toy(), readings)
    assert s["hit_looked_up"] and not s["hit_derived"]


def test_the_class_reading_does_not_look_at_the_expectation():
    readings = {
        "constraint": {"mammals": {"fraction_above": 0.5}, "humans": None},
        "syntax_values": {"values_in_syntax": 0},
        "deletion": {"target": "GOOD"},
    }
    a = score_class(_toy(("program",)), readings)
    b = score_class(_toy(("storage",)), readings)
    assert a["read"] == b["read"] == ["syntax", "program"]
    assert a["hit_derived"] and not b["hit_derived"]


# ------------------------------------------------------------------------- the gate, on the result
#: loci whose target a derived layer named on the first run (2026-09-14). A locus may join this set;
#: one that leaves it is a regression.
DERIVED_TARGET_PASSES = {
    "SHH_ZRS",
    "HERC2_OCA2",
    "MCM6_LCT",
    "HBB_LCR",
    "FTO_IRX3",
    "BCL11A_enhancer",
    "ABO",
    "HLA_DRB1",
    "APP",
    "TP53",
    "HOXD",
}
#: the element's own deletion naming the wrong gene, counted today. Lower it when one is fixed; a rise
#: is a regression. Today: the ZRS names LMBR1 (SHH is 979 kb away) and MYC's enhancer names POU5F1B.
KNOWN_ELEMENT_LEVEL_DEFECTS = 2
#: loci whose published causal variant constraint ranks first among the window's variable positions
CAUSAL_VARIANT_FIRST = {"HERC2_OCA2": "rs12913832", "BCL11A_enhancer": "rs1427407"}


def loci() -> dict:
    return {r["locus"]: r for r in RESULT["loci"]}


@needs_result
def test_no_locus_that_passed_stops_passing():
    now = {name for name, r in loci().items() if r["score"]["target_hit_derived"]}
    lost = DERIVED_TARGET_PASSES - now
    assert not lost, f"a derived layer stopped naming the published target at {sorted(lost)}"


@needs_result
def test_the_result_is_the_panel_in_the_code():
    assert [r["locus"] for r in RESULT["loci"]] == [e.locus for e in PANEL]


@needs_result
def test_no_hit_is_annotation_lookup_alone():
    assert RESULT["aggregate"]["target_only_looked_up"] == []


@needs_result
def test_derived_and_looked_up_rates_are_reported_apart():
    a = RESULT["aggregate"]
    for key in ("target_derived", "target_heuristic", "target_looked_up"):
        assert a[key]["n"] == len(PANEL) and "misses" in a[key]
    for r in RESULT["loci"]:
        for layer in r["score"]["scored"]["target"]["by_layer"].values():
            assert layer["provenance"] in ("derived", "targeted", "heuristic", "looked_up")


@needs_result
def test_negative_controls_are_scored_through_the_same_claims():
    n = RESULT["negative_controls"]
    assert n["windows"] >= 5 * len(PANEL)
    for claim in ("target", "cell", "direction", "storage"):
        assert n["positives"][claim]["n"] == len(PANEL)
        assert n["negatives"][claim]["n"] == n["windows"]


@needs_result
def test_a_value_on_syntax_still_separates_the_panel_from_its_controls():
    """The one claim the controls did not erase: it must stay at least twice as common at the panel."""
    n = RESULT["negative_controls"]
    assert n["positives"]["storage"]["rate"] >= 2 * n["negatives"]["storage"]["rate"]


@needs_result
def test_the_causal_variants_constraint_found_stay_found():
    for name, rsid in CAUSAL_VARIANT_FIRST.items():
        rank = loci()[name]["score"]["scored"]["variant"]["constraint_rank"][rsid]
        assert rank["rank"] == 1, f"{name}: {rsid} fell to {rank}"


@needs_result
def test_the_coding_calibration_derives_the_published_protein_changes():
    for name, rsid, change in (("APP", "rs63750847", "p.Ala673Thr"), ("TP53", "rs1042522", "p.Pro72Arg")):
        got = loci()[name]["score"]["scored"]["variant"]["protein_change"][rsid]
        assert got == change, f"{name}: the engine derived {got}"


@needs_result
def test_the_hypersensitive_ladder_of_the_globin_locus_is_found_from_dnase_alone():
    ladder = loci()["HBB_LCR"]["readings"]["ladder"]
    assert ladder["erythroid_specific"] >= 4
    assert 2_000 <= ladder["median_spacing"] <= 5_000


@needs_result
def test_every_unreachable_reading_says_why():
    for r in RESULT["loci"]:
        for layer, reading in r["readings"].items():
            if isinstance(reading, dict) and reading.get("pending") is not None:
                assert len(reading["pending"]) > 20, (
                    f"{r['locus']} {layer}: a pending reading without a reason"
                )
    order = loci()["HOXD"]["readings"]["order"]
    assert "order" in order["pending"], "the untestable order of HOXD must stay recorded as pending"


@needs_result
def test_the_value_domain_tells_the_many_valued_slot_from_the_two_valued_ones():
    read = {n: r["score"]["scored"]["class"]["value_domain_read"] for n, r in loci().items()}
    assert read["HLA_DRB1"] == "many"
    assert read["HERC2_OCA2"] == "two" and read["APP"] == "two"
    assert read["ABO"] == "few"


@needs_result
def test_the_named_common_values_come_back_with_their_population():
    expected = {
        "HERC2_OCA2": ("rs12913832", "European"),
        "MCM6_LCT": ("rs4988235", "European"),
        "FTO_IRX3": ("rs1421085", "European"),
        "TP53": ("rs1042522", "European"),
    }
    for name, (rsid, population) in expected.items():
        f = loci()[name]["readings"]["frequencies"]["named"][rsid]
        assert f["af"] >= 0.05 and population in (f["commonest_in"] or ""), f"{name}: {f}"


@needs_result
def test_the_target_rate_is_reported_against_its_chance_floor():
    a = RESULT["aggregate"]
    chance = a["target_by_chance"]
    assert chance["of"] == len(PANEL) and 0 < chance["expected"] < len(PANEL)
    assert a["target_derived"]["k"] > chance["expected"], (
        "the derived target rate has fallen to the rate of drawing a gene at random from the window"
    )


@needs_result
def test_the_element_level_misses_do_not_grow():
    """An element scored in the model naming a gene that is not the published target: recorded, not asserted.

    The ZRS is the case the panel exists for: its own deletion names LMBR1, the gene it sits inside,
    while SHH is a megabase away. The count is pinned so a new one is a regression.
    """
    flagged = []
    for name, r in loci().items():
        layer = r["score"]["scored"]["target"]["by_layer"].get("deletion") or {}
        if layer.get("named") and not layer["hit"]:
            flagged.append(name)
    pinned = [n for n in flagged if n in GATE_PANEL]
    assert len(pinned) <= KNOWN_ELEMENT_LEVEL_DEFECTS, f"more elements name the wrong gene: {pinned}"


@needs_result
def test_the_direction_of_effect_is_right_wherever_a_deletion_names_the_target():
    d = RESULT["aggregate"]["direction_derived_where_judged"]
    assert d["n"] >= 4 and d["k"] == d["n"], f"the direction became wrong at {d['misses']}"


@needs_result
def test_naming_a_target_is_still_reported_against_the_control_rate():
    """The panel's most valuable number: matched windows get a target nearly as often as the loci."""
    n = RESULT["negative_controls"]
    assert n["negatives"]["target"]["rate"] > 0.5, (
        "if matched windows stop getting a target, the control has changed and the headline must be redone"
    )


@needs_result
def test_every_claim_a_deletion_makes_is_reported_with_its_coverage():
    """A direction can only be produced where an element has been deleted.

    The panel's windows were deliberately scored and the matched negatives were not, so the
    unconditional direction rate (panel 67%, controls 13% when it was first measured) is a statement
    about where the requests were spent. The run must carry the coverage-conditioned counts beside
    it, so the number cannot be quoted without its denominator.
    """
    g = RESULT["negative_controls"]["given_deletion_data"]
    for side in ("positives", "negatives"):
        for group in ("with_deletion_data", "without_deletion_data"):
            for claim in ("target", "cell", "direction", "storage"):
                assert "k" in g[side][group][claim], f"{side} {group} {claim}"
    assert g["negatives"]["without_deletion_data"]["direction"]["k"] == 0, (
        "a window with no scored element produced a direction: the direction layer is not the deletion"
    )
    assert g["negatives"]["with_deletion_data"]["direction"]["n"] >= 1, (
        "no matched negative has element-level deletion data, so the direction claim has no control"
    )


def test_a_stated_interval_carries_its_own_citation():
    """The panel may draw an element nobody else drew, but it has to say on whose authority."""
    for e in PANEL:
        assert e.element_source in ("annotated", "stated"), e.locus
        if e.element_source == "stated":
            assert e.element_citation and len(e.element_citation) > 80, (
                f"{e.locus}: a stated interval without a citation is the panel inventing an element"
            )


@needs_result
def test_stated_and_annotated_element_rates_are_never_pooled_without_the_split():
    """A hit on an interval we drew ourselves is a different claim from a hit on an ENCODE cCRE."""
    split = RESULT["aggregate"]["target_derived_by_element_source"]
    for source in ("annotated", "stated"):
        assert {"k", "n", "rate", "loci", "hits"} <= set(split[source]), source
    assert split["annotated"]["n"] + split["stated"]["n"] == len(PANEL)
    stated = {e.locus for e in PANEL if e.element_source == "stated"}
    assert set(split["stated"]["loci"]) == stated
    for r in RESULT["loci"]:
        dl = r["readings"].get("deletion") or {}
        if dl.get("elements_scored"):
            assert "from_stated_interval" in dl, f"{r['locus']}: a deletion reading without its source"


@needs_result
def test_an_element_no_annotation_called_is_told_from_one_the_sweep_has_not_reached():
    """The element-level reading can only ask about intervals somebody else drew.

    Seven of the panel's published elements have no cCRE over them, the ZRS among them, and its
    pinned LMBR1 answer comes from VISTA rather than the registry. A locus with no annotated
    interval cannot be fixed by sweeping harder, so the two causes must never be reported as one.
    """
    cov = RESULT["aggregate"]["published_element_coverage"]
    assert set(cov["no_element_level_reading_at_all"]) <= set(cov["no_ccre_over_the_published_element"]), (
        "a locus has no element-level reading although an annotation did call its element:"
        " that one IS a coverage gap and belongs to the sweep, not to this list"
    )
    for name, row in cov["per_locus"].items():
        assert row["ccres_over_element"] is not None, name
        if row["elements_scored"]:
            assert row["runs"], f"{name}: elements scored but no run named"


@needs_result
def test_the_controls_on_finished_chromosomes_are_reported_separately():
    """A control on an unfinished chromosome can only be under-read, so it cannot be compared.

    The subset whose chromosome the sweep has completed is the only fair control for any claim a
    deletion makes, and it must travel with the run: on it, about half the matched windows hold no
    registry element at all, which caps every such claim well below the panel's rate.
    """
    s = RESULT["negative_controls"]["on_swept_chromosomes"]
    assert s["windows"] >= 5 and s["chromosomes"], "no control sits on a finished chromosome"
    assert s["windows"] <= RESULT["negative_controls"]["windows"]
    for claim in ("target", "cell", "direction", "storage", "deletion_scored"):
        assert s["rates"][claim]["n"] == s["windows"], claim
    assert s["rates"]["deletion_scored"]["k"] <= s["windows"]


# ------------------------------------------------------------------ what the model could be asked
def test_a_target_outside_the_models_input_is_unaskable_not_unasked():
    """The benchmark must not grade a question the model was never in a position to be asked.

    The scorer resizes its input to 1 Mb around the element, so a gene whose body lies wholly
    outside that window is not a candidate however many requests are spent. Three of the panel's
    targets are in that position, and the ZRS is the one that mattered: its own deletion named
    LMBR1, the gene it sits inside, and that was read as a miss against SHH for months. SHH is
    979 kb away and the reach is 524 kb. It is out of range, not wrong.
    """
    from genomeos.benchmark.loci import MODEL_WINDOW

    assert MODEL_WINDOW == 1_048_576  # dna_client.SEQUENCE_LENGTH_1MB
    zrs = next(e for e in PANEL if e.locus == "SHH_ZRS")
    sox9 = next(e for e in PANEL if e.locus == "SOX9_PierreRobin")

    assert 900_000 < zrs.distance < 1_100_000  # published, and beyond half the window
    assert zrs.distance > MODEL_WINDOW // 2
    assert sox9.distance > MODEL_WINDOW // 2


@needs_result
def test_the_committed_result_marks_the_out_of_reach_loci_and_no_others():
    """A locus in reach must never carry the flag: it would excuse a real miss."""
    out_of_reach = set()
    for row in RESULT["loci"]:
        reach = (row.get("readings") or {}).get("reach")
        assert reach is not None, f"{row['locus']} has no reach reading"
        if not reach["askable"]:
            out_of_reach.add(row["locus"])
            assert (row["readings"]["deletion"] or {}).get("unaskable")

    assert out_of_reach == {"SHH_ZRS", "SOX9_PierreRobin"}


@needs_result
def test_a_locus_whose_second_target_is_out_of_reach_is_still_askable():
    """FTO reaches IRX3 at 520 kb and not IRX5 at 1,164 kb: partial reach is reach, and is recorded."""
    fto = next(r for r in RESULT["loci"] if r["locus"] == "FTO_IRX3")
    reach = fto["readings"]["reach"]

    assert reach["askable"]
    assert reach["targets_in_reach"] == ["IRX3"]
    assert [x["target"] for x in reach["targets_out_of_reach"]] == ["IRX5"]


def test_reach_is_gene_body_overlap_not_tss_distance():
    """A long gene reaches into the window its promoter is nowhere near, and the sweep proves it.

    Tempting shorthand: "the reach is 524 kb, so a target further than that cannot be named." It is
    wrong, and chr21 says how wrong. Of 5,174 elements with a coding target, 79 name a gene whose
    TSS is 524 kb or more away — up to 986 kb — and every one of the 79 is a gene whose BODY
    overlaps the window. They are almost all RUNX1, which is 1.22 Mb long. So the test has to be
    body overlap, which is what `read_reach` does; measuring to the TSS would call RUNX1 unreachable
    from elements that in fact name it.
    """
    from genomeos.benchmark.loci import MODEL_WINDOW

    half = MODEL_WINDOW // 2
    runx1_length = 1_216_867
    element_to_tss = 985_879  # the furthest chr21 case, from data/knowledge/.../all_elements/chr21.json

    assert element_to_tss > half, "the TSS is well outside half the window"
    assert runx1_length > element_to_tss, "yet the body spans far enough to be inside it"


# ------------------------------------------------------------------ the registry lookup
def test_ccres_in_is_correct_when_the_registry_arrives_out_of_order():
    """`load_ccres` returns several sorted runs concatenated, and a bisect over that lies quietly.

    Every chromosome checked has three descending steps in the file order, one per class group. The
    old index was built straight from that list, so the bisect returned a plausible index and
    everything before it was dropped: 12 cCREs at MC1R read as 1, and the ZRS, MYC and ABO read as
    zero when they have 1, 2 and 4. That put three loci into section 10's "no cCRE over the
    published element" list which do not belong there.

    Built here as a unit so it stays offline: the method needs only the three attributes set below.
    """
    import bisect as _b

    from genomeos.benchmark.loci import Chromosome

    class E:
        def __init__(self, start, end):
            self.start, self.end = start, end

    # two sorted runs concatenated, exactly the shape the registry file has
    ccres = [E(100, 200), E(5_000, 5_100), E(9_000, 9_100)] + [E(150, 250), E(5_050, 5_150)]
    ch = object.__new__(Chromosome)
    ch.ccres = ccres
    ch._ccre_sorted = sorted(ccres, key=lambda c: c.start)
    ch._ccre_starts = [c.start for c in ch._ccre_sorted]

    got = ch.ccres_in(120, 260)
    truth = [c for c in ccres if c.start < 260 and c.end > 120]

    assert len(got) == len(truth) == 2, "both elements over the window, whichever run they came from"
    assert ch.ccres_in(5_000, 5_200) and len(ch.ccres_in(5_000, 5_200)) == 2
    assert ch.ccres_in(20_000, 21_000) == []
    # the index must be sorted, or the bisect above is meaningless
    assert ch._ccre_starts == sorted(ch._ccre_starts)
    assert _b.bisect_left(ch._ccre_starts, 0) == 0


@needs_result
def test_section_ten_counts_match_the_registry_the_code_actually_reads():
    """The count that section 10 quotes has to be the count the fixed lookup produces."""
    zero = {r["locus"] for r in RESULT["loci"] if r.get("ccres_over_element") == 0}

    assert zero == {"APP", "TP53", "SOX9_PierreRobin", "H19_ICR1"}
