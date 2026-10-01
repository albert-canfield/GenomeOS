# SPDX-License-Identifier: AGPL-3.0-or-later
"""The gates of the next experiment's information value, on records written by the tests only.

No dataset is read, no experiment is selected and no model request is made anywhere in this file.
Every test names the behaviour it would catch: each one fails if the module answers a design with a
volume before the design has been shown to be worth running.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from genomeos.attribution import experiment_value as ev

MODULE = Path(ev.__file__)

#: The registered reading of a second-cell-type arm that did not replicate, quoted word for word
#: from docs/ATTRIBUTION.md (2026-09-28, lane-hct116). The test's point is that this wording travels
#: unchanged and is never counted as a replication.
INCONCLUSIVE_READING = (
    "in HCT116 the frozen deletion model gains +0.022 with an interval across zero, so the CRISPRi "
    "result passes its second cell type but is not replicated there"
)

REPLICATED_CONTEXT = ev.Context(
    name="the first cell type of the design's record",
    replication_status=ev.REPLICATED,
    registered_reading="replicated in one cell type",
)

INCONCLUSIVE_CONTEXT = ev.Context(
    name="HCT116",
    replication_status=ev.INCONCLUSIVE,
    registered_reading=INCONCLUSIVE_READING,
)

UNITS = "log2 fold change in the design's own readout"


# ---------------------------------------------------------------------------
# records the tests build for themselves
# ---------------------------------------------------------------------------


def an_assay(
    *,
    detection_limit: float | None = 0.10,
    available: bool = True,
    comparable: bool | None = True,
    units: str = UNITS,
) -> ev.Assay:
    return ev.Assay(
        name="the assay named by the design's record",
        units=units,
        detection_limit=detection_limit,
        available=available,
        availability_reading="the design's record says the assay is available",
        readout_comparable=comparable,
        comparability_reading="the design's record declares the readout comparability",
    )


def a_candidate(
    i: int,
    *,
    chromosome: str = "chr1",
    position: int | None = None,
    target_gene: str | None = None,
    predictions: tuple[float, ...] = (0.0, 0.5),
    prediction_units: str | None = None,
    assay: ev.Assay | None = None,
    prior_reads: tuple[str, ...] = (),
    context: ev.Context | None = None,
) -> ev.Candidate:
    assay = assay or an_assay()
    hypothesis_units = prediction_units or assay.units
    return ev.Candidate(
        candidate_id=f"c{i}",
        chromosome=chromosome,
        position=i * 5_000_000 if position is None else position,
        target_gene=f"GENE{i}" if target_gene is None else target_gene,
        assay=assay,
        hypotheses=tuple(
            ev.Hypothesis(name=f"h{n}", predicted_outcome=v, units=hypothesis_units)
            for n, v in enumerate(predictions)
        ),
        context=context or REPLICATED_CONTEXT,
        prior_reads=prior_reads,
    )


def some_units(n: int = 12, positives: int = 4, spread: int = 5_000_000) -> tuple[ev.EligibilityUnit, ...]:
    return tuple(
        ev.EligibilityUnit(
            unit_id=f"u{i}",
            chromosome="chr1",
            position=i * spread,
            target_gene=f"GENE{i}",
            label=i < positives,
        )
        for i in range(n)
    )


def some_floors(**over: object) -> ev.Floors:
    base: dict[str, object] = {
        "min_units": 10,
        "min_positive_units": 3,
        "min_label_prevalence": 0.1,
        "min_independent_loci": 3,
        "min_distinguishing_candidates": 2,
        "target_loci": 2,
        "units_per_locus": 1,
        "registered_at": "fixed by this test before any count was taken",
    }
    base.update(over)
    return ev.Floors(**base)  # type: ignore[arg-type]


def a_design(
    *,
    claim: str = ev.DEVELOPMENT_EVIDENCE,
    unit_of_analysis: str = "element-gene pairs of the design's record",
    floors: ev.Floors | None = None,
    eligibility: object = None,
    candidates: tuple[ev.Candidate, ...] | None = None,
) -> ev.Design:
    return ev.Design(
        name="a design written by the test",
        unit_of_analysis=unit_of_analysis,
        claim=claim,
        floors=floors or some_floors(),
        eligibility=some_units() if eligibility is None else eligibility,  # type: ignore[arg-type]
        candidates=(
            (
                a_candidate(1, predictions=(0.0, 0.9)),
                a_candidate(2, predictions=(0.0, 0.6)),
                a_candidate(3, predictions=(0.0, 0.4)),
                a_candidate(4, predictions=(0.0, 0.2)),
            )
            if candidates is None
            else candidates
        ),
        source_reading="records built by the test; no dataset was read",
    )


FORBIDDEN_IN_A_REFUSAL = (
    "arms",
    "volume",
    "units_requested",
    "loci_requested",
    "selections",
    "selection_probability",
)


def _keys(node: object) -> set[str]:
    """Every key appearing anywhere in a nested result, at any depth."""
    found: set[str] = set()
    if isinstance(node, dict):
        for key, value in node.items():
            found.add(str(key))
            found |= _keys(value)
    elif isinstance(node, (list, tuple)):
        for value in node:
            found |= _keys(value)
    return found


def assert_no_volume(result: dict) -> None:
    """No volume, no sample size and no selection may appear anywhere in a refused result.

    The check is on keys at any depth, not on prose: a refusal is expected to say in words that no
    volume is reported, and must not carry the field.
    """
    present = _keys(result) & set(FORBIDDEN_IN_A_REFUSAL)
    assert not present, f"a refused design reported {sorted(present)}"
    assert result["refused"], "a refused design must carry its reason"
    assert result["refused_at"] in ev.GATES
    assert json.dumps(result)  # the refusal is serialisable as it stands


# ---------------------------------------------------------------------------
# gate 1: the blinded eligibility count
# ---------------------------------------------------------------------------


def test_a_prediction_joined_to_the_eligibility_count_is_refused() -> None:
    """Catches the eligibility count being taken over records that carry model predictions."""
    design = a_design(eligibility=(a_candidate(1), a_candidate(2)))
    gate = ev.gate_eligibility(design)
    assert not gate.passed
    assert "a prediction was joined to the eligibility count" in gate.reason
    assert gate.detail["units_with_a_prediction"] == ["c1", "c2"]


def test_a_design_whose_count_is_unblinded_is_never_answered_with_a_volume() -> None:
    """Catches a volume reported for a design whose eligibility count was not blinded."""
    result = ev.plan(a_design(eligibility=(a_candidate(1), a_candidate(2))))
    assert result["refused_at"] == ev.GATE_ELIGIBILITY
    assert_no_volume(result)


def test_the_unit_of_analysis_must_be_named() -> None:
    """Catches a count reported without the population it covers."""
    gate = ev.gate_eligibility(a_design(unit_of_analysis="  "))
    assert not gate.passed
    assert "unit of analysis is not named" in gate.reason


def test_below_the_independent_locus_floor_the_design_is_refused_before_any_candidate() -> None:
    """Catches a design answered with a volume although it has too few independent loci.

    It also catches the candidates being examined before gate 1 passes: a refused design carries no
    per-candidate information value at all.
    """
    result = ev.plan(a_design(floors=some_floors(min_independent_loci=99)))
    assert result["refused_at"] == ev.GATE_ELIGIBILITY
    assert "below the registered floor for independent loci" in result["refused"]
    assert "information_value" not in result
    assert_no_volume(result)


def test_below_the_unit_and_prevalence_floors_the_design_is_refused() -> None:
    """Catches a design proceeding below its registered coverage and prevalence floors."""
    thin = ev.plan(a_design(eligibility=some_units(n=4, positives=2), floors=some_floors()))
    assert thin["refused_at"] == ev.GATE_ELIGIBILITY
    assert "below the registered floor for units" in thin["refused"]

    rare = ev.plan(
        a_design(eligibility=some_units(n=12, positives=3), floors=some_floors(min_label_prevalence=0.5))
    )
    assert rare["refused_at"] == ev.GATE_ELIGIBILITY
    assert "below the registered floor for label prevalence" in rare["refused"]


def test_every_eligibility_count_names_its_population() -> None:
    """Catches a bare number reported with no population beside it."""
    gate = ev.gate_eligibility(a_design())
    assert gate.passed
    assert len(gate.counts) == 4
    for count in gate.counts:
        assert count.population.strip()


def test_a_count_cannot_be_built_without_a_population() -> None:
    with pytest.raises(ValueError, match="must name the population"):
        ev.Count(12, "")


def test_a_gate_result_cannot_be_built_without_a_reason() -> None:
    """Catches a refusal returned with no reason."""
    with pytest.raises(ValueError, match="must state its reason"):
        ev.GateResult(ev.GATE_ELIGIBILITY, False, "   ")


# ---------------------------------------------------------------------------
# the registered locus convention
# ---------------------------------------------------------------------------


def test_independent_loci_follow_the_convention_not_the_observation_count() -> None:
    """Catches independent loci being inferred from how many observations there are.

    Five records inside one 1 Mb window on one chromosome are one locus, not five.
    """
    records = tuple(
        ev.EligibilityUnit(f"u{i}", "chr1", 1_000_000 + i * 100_000, f"GENE{i}", True) for i in range(5)
    )
    loci = ev.independent_loci(records)
    assert len(loci) == 1
    assert loci[0] == ("u0", "u1", "u2", "u3", "u4")


def test_the_same_target_gene_is_one_locus_wherever_the_records_lie() -> None:
    records = (
        ev.EligibilityUnit("u1", "chr1", 1_000, "SAME", True),
        ev.EligibilityUnit("u2", "chr7", 90_000_000, "SAME", True),
    )
    assert ev.independent_loci(records) == [("u1", "u2")]


def test_records_farther_apart_than_the_window_are_separate_loci() -> None:
    records = (
        ev.EligibilityUnit("u1", "chr1", 1_000_000, "A", True),
        ev.EligibilityUnit("u2", "chr1", 1_000_000 + ev.LOCUS_WINDOW + 1, "B", True),
    )
    assert ev.independent_loci(records) == [("u1",), ("u2",)]


def test_the_window_boundary_is_inclusive_and_the_grouping_is_transitive() -> None:
    records = (
        ev.EligibilityUnit("u1", "chr1", 0, "A", True),
        ev.EligibilityUnit("u2", "chr1", ev.LOCUS_WINDOW, "B", True),
        ev.EligibilityUnit("u3", "chr1", 2 * ev.LOCUS_WINDOW, "C", True),
    )
    assert ev.independent_loci(records) == [("u1", "u2", "u3")]


def test_a_different_chromosome_is_a_different_locus() -> None:
    records = (
        ev.EligibilityUnit("u1", "chr1", 1_000, "A", True),
        ev.EligibilityUnit("u2", "chr2", 1_000, "B", True),
    )
    assert ev.independent_loci(records) == [("u1",), ("u2",)]


def test_the_locus_window_is_the_registered_one_megabase() -> None:
    """Pins the span, so it cannot be widened or narrowed after a count has been seen.

    The same span is registered for the second-cell-type count in
    genomeos/attribution/cell2.py (INDEPENDENT_LOCUS_SPAN, 2026-10-01); the project states one
    convention and the two must not diverge.
    """
    assert ev.LOCUS_WINDOW == 1_000_000


def test_the_locus_convention_states_what_it_is_not() -> None:
    """Catches the grouping being presented as biological independence."""
    text = ev.LOCUS_CONVENTION
    assert "operational grouping" in text
    assert "NOT established biological independence" in text
    assert f"{ev.LOCUS_WINDOW:,} bp" in text
    assert "upper bound" in text
    assert ev.plan(a_design())["locus_convention"] == text


# ---------------------------------------------------------------------------
# gate 2: assay feasibility
# ---------------------------------------------------------------------------


def test_an_unstated_detection_limit_refuses_the_design() -> None:
    """Catches a design proceeding without saying what its assay can see."""
    design = a_design(
        candidates=(
            a_candidate(1, assay=an_assay(detection_limit=None)),
            a_candidate(2),
        )
    )
    result = ev.plan(design)
    assert result["refused_at"] == ev.GATE_FEASIBILITY
    assert "state no positive detection limit" in result["refused"]
    assert_no_volume(result)


def test_an_undeclared_readout_comparability_refuses_the_design() -> None:
    """Catches the failure that stopped the paired-enhancer survey: comparability never declared."""
    design = a_design(candidates=(a_candidate(1, assay=an_assay(comparable=None)), a_candidate(2)))
    result = ev.plan(design)
    assert result["refused_at"] == ev.GATE_FEASIBILITY
    assert "do not declare" in result["refused"]
    assert "comparable" in result["refused"]
    assert_no_volume(result)


def test_an_unavailable_assay_drops_the_candidate_and_keeps_the_design() -> None:
    design = a_design(
        candidates=(
            a_candidate(1, predictions=(0.0, 0.9)),
            a_candidate(2, predictions=(0.0, 0.6)),
            a_candidate(3, predictions=(0.0, 0.4), assay=an_assay(available=False)),
        )
    )
    gate = ev.gate_feasibility(design)
    assert gate.passed
    assert gate.detail["feasible"] == ["c1", "c2"]
    assert "c3" in gate.detail["unavailable"]


def test_a_readout_declared_not_comparable_drops_the_candidate() -> None:
    design = a_design(
        candidates=(
            a_candidate(1, predictions=(0.0, 0.9)),
            a_candidate(2, predictions=(0.0, 0.6)),
            a_candidate(3, predictions=(0.0, 0.4), assay=an_assay(comparable=False)),
        )
    )
    gate = ev.gate_feasibility(design)
    assert gate.passed
    assert gate.detail["readout_not_comparable"] == ["c3"]
    assert "c3" not in gate.detail["feasible"]


def test_a_design_with_no_feasible_candidate_is_refused() -> None:
    design = a_design(candidates=(a_candidate(1, assay=an_assay(available=False)),))
    result = ev.plan(design)
    assert result["refused_at"] == ev.GATE_FEASIBILITY
    assert_no_volume(result)


# ---------------------------------------------------------------------------
# gate 3: distinguishable hypotheses
# ---------------------------------------------------------------------------


def test_two_hypotheses_with_identical_predictions_are_indistinguishable() -> None:
    """Catches two hypotheses predicting the same outcome passing the distinguishability gate."""
    value = ev.information_value(a_candidate(1, predictions=(0.42, 0.42)))
    assert value.separation == 0
    assert not value.separates
    assert "predict the same outcome" in value.refusal
    assert "not a limit of the assay" in value.refusal


def test_a_design_of_identical_predictions_is_refused_without_a_volume() -> None:
    design = a_design(
        candidates=(
            a_candidate(1, predictions=(0.42, 0.42)),
            a_candidate(2, predictions=(0.1, 0.1)),
        )
    )
    result = ev.plan(design)
    assert result["refused_at"] == ev.GATE_DISTINGUISHABILITY
    assert result["gates"][-1]["detail"]["identical_predictions"] == ["c1", "c2"]
    assert_no_volume(result)


def test_a_separation_below_the_detection_limit_is_refused() -> None:
    """Catches a candidate kept although the assay could not report its separation."""
    value = ev.information_value(
        a_candidate(1, predictions=(0.0, 0.05), assay=an_assay(detection_limit=0.10))
    )
    assert value.separation == pytest.approx(0.05)
    assert not value.separates
    assert "below the assay's detection limit" in value.refusal


def test_a_separation_at_the_detection_limit_is_kept() -> None:
    """Fixes the boundary: at the limit is reportable, below it is not."""
    value = ev.information_value(
        a_candidate(1, predictions=(0.0, 0.10), assay=an_assay(detection_limit=0.10))
    )
    assert value.separates
    assert value.refusal == ""


def test_undetectable_candidates_do_not_reach_an_arm() -> None:
    design = a_design(
        candidates=(
            a_candidate(1, predictions=(0.0, 0.05)),
            a_candidate(2, predictions=(0.0, 0.02)),
        )
    )
    result = ev.plan(design)
    assert result["refused_at"] == ev.GATE_DISTINGUISHABILITY
    assert sorted(result["gates"][-1]["detail"]["below_the_detection_limit"]) == ["c1", "c2"]
    assert_no_volume(result)


def test_a_single_hypothesis_cannot_be_separated() -> None:
    value = ev.information_value(a_candidate(1, predictions=(0.9,)))
    assert not value.separates
    assert "has nothing to be separated from" in value.refusal


def test_predictions_must_be_in_the_assay_units() -> None:
    """Catches a separation compared with a detection limit on a different scale."""
    value = ev.information_value(
        a_candidate(1, predictions=(0.0, 0.9), prediction_units="percent of control")
    )
    assert not value.separates
    assert "not on one scale" in value.refusal


def test_every_information_value_records_the_predicted_outcomes_and_the_limit() -> None:
    result = ev.plan(a_design())
    assert result["refused"] is None
    for record in result["information_value"]:
        assert record["detection_limit"] == pytest.approx(0.10)
        assert len(record["predicted"]) == 2
        assert all("predicted_outcome" in p for p in record["predicted"])
        assert isinstance(record["separates"], bool)


# ---------------------------------------------------------------------------
# gate 4: prior exposure and replication
# ---------------------------------------------------------------------------


def test_a_candidate_already_read_is_development_evidence() -> None:
    """Catches a re-read candidate being labelled fresh validation."""
    read = a_candidate(1, prior_reads=("read on 2026-09-20 by the design's own development work",))
    assert read.evidence_role == ev.DEVELOPMENT_EVIDENCE
    value = ev.information_value(read)
    assert value.evidence_role == ev.DEVELOPMENT_EVIDENCE
    assert value.may_be_reported_as_fresh_validation is False
    assert value.as_dict()["may_be_reported_as_fresh_validation"] is False


def test_a_fresh_validation_claim_is_refused_when_a_candidate_was_read() -> None:
    """Catches held-out data being re-read and reported as fresh validation."""
    design = a_design(
        claim=ev.FRESH_VALIDATION,
        candidates=(
            a_candidate(1, predictions=(0.0, 0.9), prior_reads=("read while the model was built",)),
            a_candidate(2, predictions=(0.0, 0.6)),
            a_candidate(3, predictions=(0.0, 0.4)),
        ),
    )
    result = ev.plan(design)
    assert result["refused_at"] == ev.GATE_FRESHNESS
    assert "is not fresh validation" in result["refused"]
    assert_no_volume(result)


def test_a_fresh_validation_claim_passes_when_nothing_was_read() -> None:
    result = ev.plan(a_design(claim=ev.FRESH_VALIDATION))
    assert result["refused"] is None
    assert "arms" in result


def test_an_inconclusive_context_is_not_counted_as_a_replication() -> None:
    """Catches an inconclusive second-cell-type arm counted as a replication."""
    assert INCONCLUSIVE_CONTEXT.counts_as_replication is False
    design = a_design(
        candidates=(
            a_candidate(1, predictions=(0.0, 0.9), context=INCONCLUSIVE_CONTEXT),
            a_candidate(2, predictions=(0.0, 0.6), context=INCONCLUSIVE_CONTEXT),
            a_candidate(3, predictions=(0.0, 0.4), context=INCONCLUSIVE_CONTEXT),
        )
    )
    gate = ev.gate_freshness(design, [ev.information_value(c) for c in design.candidates])
    assert gate.passed
    assert gate.detail["counted_as_replication"] == []
    assert gate.detail["contexts"]["HCT116"]["counts_as_replication"] is False
    assert gate.detail["contexts"]["HCT116"]["registered_reading"] == INCONCLUSIVE_READING


def test_a_replication_claim_is_refused_when_every_context_is_inconclusive() -> None:
    design = a_design(
        claim=ev.REPLICATION,
        candidates=(
            a_candidate(1, predictions=(0.0, 0.9), context=INCONCLUSIVE_CONTEXT),
            a_candidate(2, predictions=(0.0, 0.6), context=INCONCLUSIVE_CONTEXT),
        ),
    )
    result = ev.plan(design)
    assert result["refused_at"] == ev.GATE_FRESHNESS
    assert "stays inconclusive and is not counted as a replication" in result["refused"]
    assert INCONCLUSIVE_READING in result["refused"]
    assert_no_volume(result)


def test_a_context_cannot_assert_a_status_outside_the_registered_set() -> None:
    """Catches an inconclusive reading promoted by the caller rather than derived."""
    with pytest.raises(ValueError, match="replication status must be one of"):
        ev.Context("somewhere", "replicated-ish", "a reading")
    with pytest.raises(ValueError, match="must carry its registered reading"):
        ev.Context("somewhere", ev.REPLICATED, "   ")


# ---------------------------------------------------------------------------
# the two arms
# ---------------------------------------------------------------------------


def test_every_selection_records_its_selection_probability() -> None:
    """Catches a selected candidate reported without the probability that selected it."""
    result = ev.plan(a_design())
    assert result["refused"] is None
    for arm in ev.ARMS:
        selections = result["arms"][arm]["selections"]
        assert selections
        for s in selections:
            assert isinstance(s["selection_probability"], float)
            assert 0.0 < s["selection_probability"] <= 1.0
            assert s["scheme"].strip()


def test_a_selection_cannot_be_recorded_without_its_probability() -> None:
    for bad in (None, 0.0, -0.1, 1.5, True):
        with pytest.raises(ValueError):
            ev.Selection("c1", ev.REPRESENTATIVE, bad, "a scheme", True)  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="must name the scheme"):
        ev.Selection("c1", ev.REPRESENTATIVE, 0.5, "  ", True)


def test_the_two_arms_cannot_be_merged() -> None:
    """Catches the representative and disagreement arms counted as one sample."""
    design = a_design()
    values = [ev.information_value(c) for c in design.candidates]
    representative = ev.select_representative(design, values)
    disagreement = ev.select_disagreement(design, values)
    assert ev.pool(representative) == len(representative)
    assert ev.pool(disagreement) == len(disagreement)
    with pytest.raises(ValueError, match="cannot be merged"):
        ev.pool(list(representative) + list(disagreement))


def test_the_arms_are_reported_separately_and_answer_different_questions() -> None:
    result = ev.plan(a_design())
    arms = result["arms"]
    assert set(arms) == set(ev.ARMS)
    assert arms[ev.REPRESENTATIVE]["question"] == "how common a function is"
    assert arms[ev.DISAGREEMENT]["question"] == "which model is right"
    assert all(s["representative_of_population"] for s in arms[ev.REPRESENTATIVE]["selections"])
    assert not any(s["representative_of_population"] for s in arms[ev.DISAGREEMENT]["selections"])
    schemes = {arm: {s["scheme"] for s in arms[arm]["selections"]} for arm in ev.ARMS}
    assert schemes[ev.REPRESENTATIVE].isdisjoint(schemes[ev.DISAGREEMENT])
    assert result["merge_refused"] == ev.MERGE_REFUSED


def test_the_disagreement_arm_takes_the_largest_separations() -> None:
    design = a_design()
    values = [ev.information_value(c) for c in design.candidates]
    selected = [s.candidate_id for s in ev.select_disagreement(design, values)]
    assert selected == ["c1", "c2"]
    assert all(s.selection_probability == 1.0 for s in ev.select_disagreement(design, values))


def test_the_representative_arm_records_the_inclusion_probability_of_the_draw() -> None:
    """Four eligible loci, two drawn, one candidate per locus: 2/4 times 1/1."""
    design = a_design()
    values = [ev.information_value(c) for c in design.candidates]
    selections = ev.select_representative(design, values)
    assert len(selections) == design.floors.target_loci
    for s in selections:
        assert s.selection_probability == pytest.approx(2 / 4)
    assert ev.select_representative(design, values) == selections


def test_the_representative_draw_splits_the_probability_inside_a_locus() -> None:
    """Two candidates in one 1 Mb window are one locus, so each carries half the locus's draw."""
    design = a_design(
        floors=some_floors(target_loci=1),
        candidates=(
            a_candidate(1, chromosome="chr1", position=1_000_000, target_gene="A", predictions=(0.0, 0.9)),
            a_candidate(2, chromosome="chr1", position=1_500_000, target_gene="B", predictions=(0.0, 0.8)),
            a_candidate(3, chromosome="chr5", position=9_000_000, target_gene="C", predictions=(0.0, 0.7)),
        ),
    )
    values = [ev.information_value(c) for c in design.candidates]
    selections = ev.select_representative(design, values)
    assert len(selections) == 1
    chosen = selections[0]
    expected = (1 / 2) * (1 / 2) if chosen.candidate_id in {"c1", "c2"} else (1 / 2) * 1.0
    assert chosen.selection_probability == pytest.approx(expected)


# ---------------------------------------------------------------------------
# volume, last
# ---------------------------------------------------------------------------


def test_a_volume_appears_only_when_every_gate_has_passed() -> None:
    """Catches a sample size returned while a gate is unmet."""
    passing = ev.plan(a_design())
    assert passing["refused"] is None
    assert [g["gate"] for g in passing["gates"]] == list(ev.GATES)
    assert all(g["passed"] for g in passing["gates"])
    for arm in ev.ARMS:
        vol = passing["arms"][arm]["volume"]
        assert vol["units_requested"]["value"] == 2
        assert vol["units_requested"]["population"].strip()
        assert vol["loci_available"]["value"] == 4
        assert "not a power calculation" in vol["basis"]

    for refused in (
        a_design(floors=some_floors(min_units=99)),
        a_design(candidates=(a_candidate(1, assay=an_assay(comparable=None)), a_candidate(2))),
        a_design(candidates=(a_candidate(1, predictions=(0.2, 0.2)), a_candidate(2, predictions=(0.3, 0.3)))),
        a_design(
            claim=ev.FRESH_VALIDATION,
            candidates=(
                a_candidate(1, predictions=(0.0, 0.9), prior_reads=("read before",)),
                a_candidate(2, predictions=(0.0, 0.6)),
            ),
        ),
    ):
        assert_no_volume(ev.plan(refused))


def test_the_gates_run_in_the_registered_order() -> None:
    assert ev.GATES == (
        ev.GATE_ELIGIBILITY,
        ev.GATE_FEASIBILITY,
        ev.GATE_DISTINGUISHABILITY,
        ev.GATE_FRESHNESS,
    )


def test_the_module_sets_no_floor_value() -> None:
    """Catches a threshold living in the module, where it could be moved after an outcome."""
    with pytest.raises(TypeError):
        ev.Floors()  # type: ignore[call-arg]


# ---------------------------------------------------------------------------
# nothing is hard-coded
# ---------------------------------------------------------------------------


def test_no_dataset_or_cell_line_name_is_hard_coded_in_the_module() -> None:
    """Catches a dataset name or number written into the module instead of read from the record."""
    source = MODULE.read_text(encoding="utf-8")
    for name in (
        "HCT116",
        "K562",
        "GM12878",
        "Jurkat",
        "WTC11",
        "HepG2",
        "IMR-90",
        "ENCODE",
        "Gasperini",
        "Schraivogel",
        "Nasser",
        "Fulco",
        "Xie",
        "Morris",
        "GTEx",
        "IGVF",
        "lentiMPRA",
        "AlphaGenome",
        "VISTA",
        "ClinVar",
        "Perturb-seq",
        "FlowFISH",
        "CRISPRi",
    ):
        assert name not in source, f"{name!r} is hard-coded in {MODULE.name}"


def test_every_reading_in_the_result_comes_from_the_records() -> None:
    """The context reading in the result is the record's, word for word."""
    design = a_design(
        candidates=(
            a_candidate(1, predictions=(0.0, 0.9), context=INCONCLUSIVE_CONTEXT),
            a_candidate(2, predictions=(0.0, 0.6), context=REPLICATED_CONTEXT),
        )
    )
    result = ev.plan(design)
    readings = {r["context"]: r["context_registered_reading"] for r in result["information_value"]}
    assert readings["HCT116"] == INCONCLUSIVE_READING
    assert readings[REPLICATED_CONTEXT.name] == REPLICATED_CONTEXT.registered_reading
