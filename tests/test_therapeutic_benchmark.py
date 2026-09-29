# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""The therapeutic pipeline against targets whose answer is already known.

The benchmark itself needs the network and several minutes
(`scripts/therapeutic_benchmark.py`); these tests read its committed result,
so a regression that loses a known target fails in CI within a second.

Three questions are kept apart on purpose. Recovering the target and calling
the route correctly are asserted, because we believe them. Whether the
mechanism ranked first is defensible was *recorded* rather than asserted for as
long as the pipeline got it wrong, because pretending otherwise would have made
the benchmark decorative. The count is pinned: it may fall, never rise, and it
now reads 0. A pin at 0 is only worth having if the thing it counts can still
be counted, so `test_a_preferred_mechanism_is_never_merely_unrefused` checks the
property the last two defects violated rather than trusting the zero.
"""

from __future__ import annotations

import pytest

from genomeos.results import load_result

RESULT = load_result("therapeutic_benchmark")
pytestmark = pytest.mark.skipif(not RESULT, reason="benchmark result not present")

#: Known defects, counted today. Lower this when one is fixed; a rise is a regression.
#: 3 at the start: an agonist antibody offered against an activating driver, and
#: two mechanisms that headed a list on a requirement nobody had answered —
#: blocking_antibody for BRAF and adcp for PIK3CA, both at 0.25 against a
#: cytoplasmic protein whose membrane compartment is curated at 0.45, under the
#: 0.6 reachability threshold. The agonist went on 2026-09-11, the other two on
#: 2026-09-15 when a preferred mechanism was required to be established and not
#: merely unrefused. Both are still scored and still listed, now as the nearest
#: provisional mechanism with `surface_accessible` named as the open question.
KNOWN_MECHANISM_DEFECTS = 0  # was 2; the two unanswered-requirement cases closed 2026-09-15


def rows():
    return RESULT["rows"]


def test_every_known_target_is_recovered():
    lost = [r["gene"] for r in rows() if not r["recovered"]]
    assert not lost, f"the pipeline lost known targets: {lost}"


def test_the_approved_antibody_targets_are_found_as_surface_targets():
    """EGFR, ERBB2 and CD19 have approved antibodies. If these are not surface targets, nothing is."""
    surface = [r for r in rows() if r["expected"] == "surface"]
    assert len(surface) >= 3
    for r in surface:
        assert r["pass"], f"{r['gene']} has an approved antibody but was not called a surface target"
        assert (r["surface_accessibility"] or 0) > 0.5


def test_no_surface_route_is_claimed_for_an_intracellular_target():
    """The four small-molecule cases pass by *not* claiming a route that does not exist."""
    out = [r for r in rows() if r["expected"] == "out_of_scope_expected"]
    assert len(out) >= 4
    for r in out:
        assert r["pass"], f"{r['gene']}: {r['verdict']}"
        assert (r["surface_accessibility"] or 0) <= 0.5


def test_an_intracellular_driver_still_gets_the_peptide_route():
    """GenomeOS's distinctive claim: a nuclear or cytoplasmic mutation is never 'no target'.

    Scoped to point mutations, which is what it always meant and could leave
    implicit while every case was one. A peptide comes from a changed sequence,
    and for a fusion the changed sequence is the junction, which GenomeOS does
    not reconstruct: the EML4-ALK case therefore offers no peptide route and
    must not be read as this claim failing.
    """
    missense = [
        r for r in rows() if r["expected"] == "out_of_scope_expected" and r["driver_call"] == "point mutation"
    ]
    assert len(missense) >= 4
    assert all(r.get("peptide_route") for r in missense), (
        "an intracellular missense driver must still offer the peptide/HLA route"
    )


def test_a_fusion_is_judged_on_the_product_and_not_the_curated_gene():
    """The case where every database says surface and the tumour says otherwise.

    ALK is curated as a single-pass receptor with a signal peptide and a
    1,020-residue ectodomain. As the fusion's 3' partner it keeps none of it,
    and the approved drugs are small molecules acting inside the cell. The
    benchmark passes this case by *not* claiming a surface route — the same way
    the four intracellular missense cases pass — and it could not have done so
    before the orientation was read.
    """
    alk = next(r for r in rows() if r["driver_call"] == "structural variant")
    assert alk["gene"] == "ALK"
    assert alk["recovered"] and alk["pass"]
    assert alk["target_class"] == "intracellular_only"
    assert alk["surface_accessibility"] == 0.0, "scored from the product the fusion makes"
    assert not alk["best_mechanism"], "no antibody-like route exists against a cytoplasmic kinase"


def test_the_known_mechanism_defects_do_not_grow():
    flagged = [r["gene"] for r in rows() if not r.get("mechanism_sane", True)]
    assert len(flagged) <= KNOWN_MECHANISM_DEFECTS, (
        f"the top-ranked mechanism became indefensible for more targets: {flagged}"
    )
    for r in rows():
        assert "mechanism_sane" in r, "the benchmark stopped scoring the mechanism question"


def test_a_preferred_mechanism_is_never_merely_unrefused():
    """The property the two closed defects violated, checked rather than counted.

    A gate has three answers and the pipeline read two of them. Requiring the
    count of defects to be zero only helps if the pipeline can still report one,
    so this asserts the rule directly: a mechanism that heads a list has had its
    hard requirements answered. And a row with no preferred mechanism has to say
    which question is in the way, because the failure mode of the fix is a
    pipeline that goes quiet instead of one that says "not established".
    """
    for r in rows():
        if not r["recovered"]:
            continue
        assert "best_mechanism_established" in r, "the benchmark stopped recording the distinction"
        if r["best_mechanism"]:
            assert r["best_mechanism_established"], (
                f"{r['gene']}: {r['best_mechanism']} heads the list with "
                f"{r['requirements_unanswered']} unanswered"
            )
        elif r["nearest_provisional_mechanism"]:
            assert r["requirements_unanswered"], (
                f"{r['gene']}: a provisional mechanism must name the requirement that is open"
            )


def test_a_target_is_reached_by_something_other_than_a_point_mutation():
    """Six point mutations was a narrower benchmark than it read as.

    Copy number, structural variants and expression all reach the candidate
    list, and until CD19 was added no scored case turned on any of them, so the
    end-to-end evidence for three of the four routes was zero. CD19 carries no
    alteration in any tumour and is the target of four approved therapies, so
    the only measurement that can reach it is the patient's own RNA.
    """
    calls = RESULT["cases_by_driver_call"]
    assert sum(calls.values()) == RESULT["cases"]
    assert set(calls) - {"point mutation"}, "every case is still driven by a point mutation"
    cd19 = next(r for r in rows() if r["gene"] == "CD19")
    assert cd19["driver_call"] == "expression"
    assert cd19["recovered"] and cd19["pass"]
    assert cd19["best_mechanism_established"], "the approved CD19 drugs are antibody-like"


def test_the_copy_number_route_is_scored_and_not_merely_unit_tested():
    """The same target twice, so that the route is the only thing that differs.

    ERBB2 is in this benchmark under a point mutation and again under nothing
    but its amplification. Varying one thing is the whole reason the second
    case exists: CD19 changes the target and the route together, so a failure
    there cannot be attributed. Amplification is also the clinically honest
    call, since it is what trastuzumab is prescribed on.
    """
    amp = next(r for r in rows() if r["driver_call"] == "copy number")
    assert amp["gene"] == "ERBB2"
    assert amp["recovered"] and amp["pass"]
    assert amp["target_class"] == "direct_surface", "reached as itself, not as a pathway hypothesis"
    assert amp["best_mechanism_established"], "an approved antibody must not be provisional here"
    assert amp["evidence_tier"] == "observed_alteration"
    assert amp["alteration_evidence"] == 1.0, "twelve copies are an observed alteration"
    mutated = next(r for r in rows() if r["gene"] == "ERBB2" and r["driver_call"] == "point mutation")
    assert mutated["rank"] == 1
    assert not amp["outranked_by_hypotheses"], (
        "the finding this case existed to record, now closed: nothing with no measurement in this "
        "tumour outranks the target the approved therapy is prescribed on"
    )
    assert amp["rank"] <= AMPLIFIED_TARGET_RANK, (
        f"the amplified target fell below the pinned rank {AMPLIFIED_TARGET_RANK}"
    )
    # Recorded, not repaired: both routes to the same gene still produce the same
    # number, because the tier orders the candidates and is published per
    # candidate rather than averaged into the mean. What the amplification buys
    # is the place, which is what the clinic reads. The assertion allows the
    # number to rise and not to fall.
    assert amp["score"] >= mutated["score"]


#: Where the amplified ERBB2 case ranks: the case the evidence tier was built
#: for. It may not rise. It was 4, behind KDR, EGFR and PDGFRB, none of them
#: altered in that tumour and all three named only for neighbouring the mutated
#: PIK3CA. It is 1 since the evidence tier of 2026-09-27, with every score in the
#: benchmark unchanged: the tier orders the candidates and is published per
#: candidate, and it is deliberately kept out of the weighted mean. Averaging it
#: in was implemented first and measured: it put PIK3CA above ERBB2 in this very
#: tumour, because the mean is taken over the dimensions that were available and
#: so lifts the candidate with fewer of them further, and because ERBB2's raw
#: score is clipped by the poor-normal-tissue-safety cap that PIK3CA's mediocre
#: safety walks past.
AMPLIFIED_TARGET_RANK = 1  # was 4

#: Targets outranked by candidates with nothing measured about them in this
#: patient, counted today. It may fall, never rise. It read 2 while the score had
#: no dimension for the alteration — the amplified gene and the same gene as a
#: guess scored identically at 0.494 — and it reads 0 since the evidence tier,
#: with the metric also corrected: it counted "candidates with no origins
#: record", under which CD19's own route read as burial.
BURIED_SURFACE_TARGETS = 0  # was 2; the evidence tier, 2026-09-27


def test_a_recovered_target_is_not_quietly_buried_under_hypotheses():
    """Recovery was never the hard question, and three scored questions miss this one.

    A pathway-induced candidate has no origin: it is named for neighbouring
    something altered, and needs expression evidence before it means anything.
    When such a candidate outranks a target whose approved antibody exists, the
    pipeline has recovered the target without preferring it, and every verdict
    above still reads as a pass. The count is pinned rather than asserted to
    zero, because lowering it is a scoring change and this test exists to keep
    it visible until someone makes that change deliberately.
    """
    buried = [
        r for r in rows() if r["expected"] == "surface" and r["recovered"] and r["outranked_by_hypotheses"]
    ]
    assert len(buried) <= BURIED_SURFACE_TARGETS, (
        f"more approved-antibody targets are buried under evidence-free candidates than the "
        f"pinned {BURIED_SURFACE_TARGETS}: {[(r['gene'], r['outranked_by_hypotheses']) for r in buried]}"
    )
    for r in rows():
        assert r["gene"] not in r["outranked_by_hypotheses"], (
            f"{r['gene']} outranks itself: the same gene reached twice, once by its alteration and "
            "once as a hypothesis about a neighbour"
        )


def test_a_target_is_preferred_for_its_alteration_and_not_only_recovered():
    """The question the nine cases could not ask, and the trap inside asking it.

    Recovering a target and preferring it are different results. Until
    2026-09-27 nothing in the score read the alteration: `surface_accessibility`
    reads curated localisation, which describes the gene whether or not the
    tumour touched it, so twelve copies of ERBB2 and a STRING neighbour of a
    mutated gene were worth the same and three unaltered genes outranked the
    target trastuzumab is prescribed on.

    The trap is that the obvious rule — penalise a candidate with no DNA origin
    — demotes CD19, which has no origin either and is the target of four
    approved therapies. So the tier is read from the evidence about the gene in
    this patient, and CD19 sits in its middle tier rather than its bottom one.
    This test asserts both halves: no target is under an unmeasured candidate,
    and the case that carries no alteration of any kind did not pay for it.
    """
    for r in rows():
        if not r["recovered"]:
            continue
        assert r["evidence_tier"] in ("observed_alteration", "patient_measurement"), (
            f"{r['gene']}: an approved target is not a hypothesis about a neighbour"
        )
        assert not r["outranked_by_hypotheses"], (
            f"{r['gene']} is outranked by candidates measured nowhere in this tumour: "
            f"{r['outranked_by_hypotheses']}"
        )
        if r["driver_call"] != "expression":
            assert r["alteration_evidence"] == 1.0
    cd19 = next(r for r in rows() if r["gene"] == "CD19")
    assert cd19["evidence_tier"] == "patient_measurement", (
        "CD19 carries no alteration; it is reached from this patient's own RNA, which is a "
        "measurement of this gene and not an association with another one"
    )
    assert cd19["alteration_evidence"] == 0.6
    assert cd19["rank"] <= 2, "the rule aimed at hypotheses must not demote the expression route"


def test_every_route_into_the_candidate_list_is_scored():
    """The gap this test used to describe is closed, so it now asserts the closure.

    It read "the two routes no scored case turns on have to be named as
    uncovered" while copy number and structural variants were measured by unit
    tests alone. Both are scored now, and the assertion is inverted rather than
    deleted: an end-to-end case per route, so that a route cannot quietly stop
    working while its unit test keeps passing.
    """
    assert set(RESULT["cases_by_driver_call"]) == {
        "point mutation",
        "copy number",
        "structural variant",
        "expression",
    }
    assert RESULT["cases"] == len(rows()) >= 9


def test_the_benchmark_states_what_it_does_not_cover():
    """What remains out of scope is the modality, not a route: no small molecules."""
    note = RESULT["note"].lower()
    assert "small molecule" in note
    out = [r for r in rows() if r["expected"] == "out_of_scope_expected"]
    assert len(out) >= 5, "the small-molecule cases are the ones that state the boundary"
    assert all(r["approved_modality"] == "small_molecule" for r in out)


#: Targets ranked below a candidate that no modelled modality reaches with its
#: hard requirements answered, counted today. It may fall, never rise. It reads 0
#: from the mechanism gate registered 2026-09-27, and the pin is what keeps the
#: question countable: the gate closed nothing in these nine cases, because the
#: candidates that could have outranked a target were already below it on the
#: evidence tier. It exists so that a future case cannot reintroduce the defect
#: quietly.
UNREACHABLE_ABOVE_TARGET = 0

#: The ranks the two rules of 2026-09-27 registered as unmovable: every target
#: first, CD19 second behind FCRL5, which this patient's RNA measures exactly as
#: it measures CD19. A rank may improve and may not get worse.
REGISTERED_RANKS = {"CD19": 2}


def test_the_registered_ranks_did_not_move():
    """The must-not-move condition of the mechanism gate and the magnitude tiebreak.

    Both act on the order, so the order is where they have to be checked as a
    whole and not only at the case each was written for. Every target is first
    except CD19, which sits behind one gene measured in the same patient by the
    same route.
    """
    for r in rows():
        assert r["recovered"], r["gene"]
        assert r["rank"] <= REGISTERED_RANKS.get(r["gene"], 1), (
            f"{r['gene']} ({r['driver_call']}) fell to rank {r['rank']}: {r['candidates']}"
        )


def test_no_target_is_ranked_below_a_candidate_nothing_can_be_aimed_at():
    """The gate's own metric, which the three original questions could not ask.

    A gene no modelled modality reaches, ranked above a gene with an approved
    antibody, is a list nobody can act on — and until 2026-09-27 nothing in the
    ranking asked the question. The count is pinned rather than asserted to zero,
    because lowering it is an ordering change that has to be made deliberately.
    """
    above = [(r["gene"], r["outranked_by_unreachable"]) for r in rows() if r["outranked_by_unreachable"]]
    assert len(above) <= UNREACHABLE_ABOVE_TARGET, (
        f"more targets rank below candidates no modality reaches than the pinned "
        f"{UNREACHABLE_ABOVE_TARGET}: {above}"
    )
    for r in rows():
        assert "outranked_by_unreachable" in r, "the benchmark stopped asking the gate's question"


def test_the_gate_class_agrees_with_the_preferred_mechanism():
    """A falsifier of the gate, checked directly rather than trusted.

    The gate is defined as `best_mechanism` and nothing else: a mechanism whose
    hard requirement is unanswered leaves the candidate in the gated class. If the
    two ever disagree, the gate is reading something other than what was
    registered.
    """
    for r in rows():
        assert r["mechanism_reach"] in ("established_mechanism", "no_established_mechanism")
        assert r["mechanism_reach_reason"], f"{r['gene']}: a gate class with no sentence behind it"
        if r["best_mechanism"]:
            assert r["mechanism_reach"] == "established_mechanism"
            assert r["best_mechanism"] in r["mechanism_reach_reason"]
        else:
            assert r["mechanism_reach"] == "no_established_mechanism", (
                f"{r['gene']} has no preferred mechanism and is not gated"
            )


def test_the_magnitude_is_measured_and_its_absences_are_stated():
    """Twelve copies and one missense stop reading the same, without a guess.

    The amplified case publishes the count and its distance from the diploid 2.
    The fusion case publishes three absences: no copy count, no allele fraction,
    and no observed position to look a hotspot up at. An absence stated is the
    point of the rule — a default here would read like a measurement.
    """
    amp = next(r for r in rows() if r["driver_call"] == "copy number")
    assert amp["alteration_magnitude"]["copies"] == 12.0
    assert amp["alteration_magnitude"]["copies_above_diploid"] == 10.0
    assert "12 copies against the diploid 2" in amp["alteration_magnitude_reason"]

    alk = next(r for r in rows() if r["driver_call"] == "structural variant")
    mag = alk["alteration_magnitude"]
    assert mag["copies"] is None and mag["vaf"] is None and mag["hotspot"] is None
    assert mag["quantities"] == [], "a fusion record carries no measured amount"
    assert "not imputed" in alk["alteration_magnitude_reason"]

    for r in rows():
        assert "alteration_magnitude" in r and r["alteration_magnitude_reason"]
        mag = r["alteration_magnitude"]
        assert set(mag) == {"copies", "copies_above_diploid", "vaf", "hotspot", "quantities"}
        assert mag["vaf"] is None or 0.0 <= mag["vaf"] <= 1.0


#: The tenth case, pre-registered 2026-09-28 (docs/THERAPEUTICS.md) and run once
#: after the registration was committed. It is the first case chosen for a *rule*
#: rather than for a route: a HER2-positive gastroesophageal adenocarcinoma with
#: a co-amplified MYC, from TCGA stomach adenocarcinoma, where an approved
#: antibody target and a gene no modelled modality reaches sit in the same
#: evidence tier and only the mechanism gate can order them.
GATE_CASE = "HER2-positive gastroesophageal adenocarcinoma with a co-amplified MYC"

#: The four tests below wait for their artifact rather than assert against one
#: that is not there. The tenth case was registered and run on 2026-09-28, and
#: the regenerated result could not be committed with it: every re-run rewrites
#: the date, the seven counts and the manifest's write stamp, and
#: `scripts/check_staged.py` refuses a commit that removes lines a recent commit
#: added. This result is also a headline one, so landing it means a rebuild in a
#: clean checkout and the pins in `tests/test_manifest_headlines.py` and README
#: moving from 9 to 10 in the same commit. The skip is on the case's presence, so
#: these assert the moment that commit lands; nothing here is weakened, and the
#: run they were written from is in docs/THERAPEUTICS.md under 2026-09-28.
gate_case = pytest.mark.skipif(
    not any(r["case"] == GATE_CASE for r in RESULT.get("rows", [])),
    reason="the tenth case's result is registered and run but not yet committed (docs/THERAPEUTICS.md)",
)

#: What the gate did, measured: nothing, in all ten cases. `rank_without_gate` is
#: the pipeline's own ranking key with the gate removed, so a target whose two
#: ranks agree is a target the gate did not move. An equality pin rather than a
#: ceiling: a change here means the gate has finally acted on a case, which is a
#: result to read and record, not a regression to absorb quietly.
GATE_MOVED_TARGETS: list[str] = []

#: What the tiebreak did, measured: nothing, and in the tenth case it had no
#: opportunity at all. Every group the score could not separate is a pair of
#: candidates with nothing measured about either of them in that patient, which
#: is not an accident: the score ties among unmeasured hypotheses, and an
#: unmeasured candidate carries no quantity by definition, so the tiebreak's
#: opportunity set in this benchmark is empty by construction.
TIEBREAK_CHANGED_ORDER: list[str] = []


@gate_case
def test_the_gate_case_recovers_the_approved_target_over_the_undruggable_amplification():
    """The case the mechanism gate was given, and what it measured.

    Registered predictions: ERBB2 first, nothing unreachable above it, and MYC in
    the same evidence tier so that only the gate could separate them. The first two
    hold. The third holds as a fact about the tier and is beside the point as a
    test of the gate, which is the finding: MYC is in the top tier and scores
    0.246 against the target's 0.494, the lowest of the six altered candidates, so
    it never came near the target's place.
    """
    row = next(r for r in rows() if r["case"] == GATE_CASE)
    assert row["gene"] == "ERBB2" and row["driver_call"] == "copy number"
    assert row["recovered"] and row["pass"]
    assert row["rank"] == 1
    assert row["target_class"] == "direct_surface"
    assert row["best_mechanism_established"], "an approved antibody must not be provisional here"
    assert row["evidence_tier"] == "observed_alteration"
    assert not row["outranked_by_unreachable"], "the gate's own metric, in the case chosen to move it"
    assert not row["outranked_by_hypotheses"]
    assert row["alteration_magnitude"]["copies"] == 13.0
    assert row["alteration_magnitude"]["copies_above_diploid"] == 11.0
    by_gene = {q["gene"]: q for q in row["quantities_in_this_tumour"]}
    assert by_gene["MYC"]["copies"] == 8.0, "the undruggable amplification is in the list and measured"
    assert by_gene["MYC"]["gate"] == "no_established_mechanism"
    assert by_gene["MYC"]["score"] < by_gene["ERBB2"]["score"], (
        "the case's finding: the annotations that make a gene unreachable also make it score low, "
        "so an undruggable amplification does not produce the configuration the gate was written for"
    )


@gate_case
def test_the_gate_is_measured_and_not_only_defined():
    """What the gate would have to be doing to be worth its place in the key.

    `rank_without_gate` is the same key with the gate taken out, so the difference
    between it and `rank` is the gate's whole effect on that tumour. Across ten
    cases the difference is nothing.
    """
    for r in rows():
        assert r["rank_without_gate"] is not None or not r["recovered"]
    assert RESULT["gate_moved_the_target"] == GATE_MOVED_TARGETS, (
        "the gate moved a target for the first time; read the case before changing this pin"
    )


@gate_case
def test_the_tiebreak_had_an_opportunity_and_what_came_of_it():
    """A rule that did not fire is only informative if its opportunities are counted.

    Every score-tied group in every case is recorded, so "did not fire" is
    distinguishable from "was never asked". It was never asked: every tie in the
    ten cases is between candidates with nothing measured about them in that
    patient, and a candidate with no measurement carries no quantity to compare.
    """
    assert RESULT["magnitude_tiebreak_changed_order"] == TIEBREAK_CHANGED_ORDER
    for r in rows():
        assert "magnitude_tiebreaks" in r, "the benchmark stopped counting the tiebreak's opportunities"
        for group in r["magnitude_tiebreaks"]:
            assert len(group["tied_on_score"]) > 1
            if not group["decided_by_magnitude"]:
                assert not group["changed_the_order"], (
                    "a tied group was reordered without any measured quantity deciding it"
                )


@gate_case
def test_the_allele_fraction_reaches_the_ranking_even_though_the_row_has_none():
    """A registered prediction that the run falsified, kept as it was written.

    Prediction 4 of 2026-09-28 said the tenth row would publish a variant allele
    fraction of 0.49. It publishes none, and the mistake is in the registration: a
    row publishes the *target's* magnitude, the target there is reached by its
    amplification, and the 0.49 belongs to the mutated TP53 in the same tumour.
    The amendment is additive — the candidates' own quantities are published — and
    this test pins both halves so the correction cannot quietly become a claim
    that the prediction held.
    """
    row = next(r for r in rows() if r["case"] == GATE_CASE)
    assert row["alteration_magnitude"]["vaf"] is None, (
        "the target of this case is amplified, not mutated: its magnitude is a copy count"
    )
    fractions = {q["gene"]: q["vaf"] for q in row["quantities_in_this_tumour"] if q["vaf"] is not None}
    assert fractions == {"TP53": 0.49}, (
        "the first allele fraction in this benchmark, on the candidate that actually carries it"
    )
    assert RESULT["cases_where_a_candidate_carries_an_allele_fraction"] == ["ERBB2"]
