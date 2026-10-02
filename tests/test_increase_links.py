# SPDX-License-Identifier: AGPL-3.0-or-later
"""The versioned measured-layer extractor: v1 unchanged, v2 additive, and the R2 prohibition enforced.

Almost everything here is synthetic. The assertions are about the registration's own content, about the
two versions' behaviour on rows built in the test, and about the arithmetic of additivity. The point of
them is that the default, the conflict rule, the eligibility predicate or the R2 wording cannot be moved
without a test failing.

The two tests that read the cached benchmark are marked and SKIP WITH THE NAME OF THE MISSING INPUT when
it is not on the machine, rather than passing quietly: a test that cannot run says so, which is the rule
the project adopted after a test passed because an unrelated extra happened to install a package.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from genomeos.attribution import cell2, crispri, fresh
from genomeos.attribution import increase_links as il
from genomeos.attribution import increases as inc
from genomeos.attribution import measured as ms
from genomeos.attribution import repress2 as rp

# ---- rows built here, never read from disk --------------------------------------------------------


def _pair(**over: object) -> dict[str, object]:
    base: dict[str, object] = {
        "gene": "GENE1",
        "cell": "K562",
        "dataset": "D1",
        "regulated": False,
        "effect_size": 0.4,
        "outcome": ms.INCREASE,
        "split": ms.TRAINING,
        "p_adjusted": 0.01,
    }
    base.update(over)
    return base


def _decrease(**over: object) -> dict[str, object]:
    base: dict[str, object] = {"regulated": True, "effect_size": -0.4, "outcome": ms.DECREASE}
    base.update(over)
    return _pair(**base)


def _row(pairs: list[dict[str, object]]) -> dict[str, object]:
    return {
        "id": "e1",
        "chrom": "chr21",
        "start": 1_000,
        "end": 1_500,
        "measured": {"crispri": {"pairs": pairs}},
    }


# ---- the versions, and which one a call that names none gets --------------------------------------


def test_the_two_version_names_are_the_ones_the_registration_fixed() -> None:
    assert (ms.EXTRACTOR_V1, ms.EXTRACTOR_V2) == (il.EXTRACTOR_V1, il.EXTRACTOR_V2)
    assert ms.EXTRACTORS == il.EXTRACTORS == ("v1", "v2")


def test_v1_is_the_default_and_this_is_what_keeps_every_pinned_count_still() -> None:
    assert ms.DEFAULT_EXTRACTOR == il.DEFAULT_EXTRACTOR == ms.EXTRACTOR_V1
    row = _row([_decrease(), _pair(gene="GENE2")])
    assert ms.rule_links(row) == ms.rule_links(row, ms.EXTRACTOR_V1)


def test_an_unknown_version_is_refused_rather_than_silently_treated_as_one() -> None:
    with pytest.raises(ValueError, match="unknown extractor"):
        ms.rule_links(_row([_decrease()]), "v3")


# ---- v1 is unchanged -----------------------------------------------------------------------------


def test_v1_still_raises_nothing_from_a_significant_increase() -> None:
    """The structural cause lane-increase established, asserted against the default once more."""
    assert ms.rule_links(_row([_pair()])) == []
    assert ms.rule_links(_row([_pair()]), ms.EXTRACTOR_V1) == []


def test_v1_raises_an_activates_link_from_a_decrease_and_its_action_is_unchanged() -> None:
    links = ms.rule_links(_row([_decrease()]))
    assert links == [("GENE1", "activates", 0.4, "K562", ms.TRAINING)]


def test_under_v1_no_link_anywhere_can_carry_the_inhibits_action() -> None:
    """`regulated` is awarded to a significant decrease only, so v1 cannot reach the branch."""
    row = _row([_decrease(), _pair(gene="GENE2"), _decrease(gene="GENE3", effect_size=-0.9)])
    assert {a for _g, a, _s, _c, _sp in ms.rule_links(row)} == {"activates"}


# ---- v2 reaches the branch that was already written ----------------------------------------------


def test_v2_raises_an_inhibits_link_from_a_significant_increase() -> None:
    links = ms.rule_links(_row([_pair()]), ms.EXTRACTOR_V2)
    assert links == [("GENE1", "inhibits", 0.4, "K562", ms.TRAINING)]


def test_v2_takes_its_action_from_the_same_line_v1_does() -> None:
    """Not a new direction rule: the increase arm passes through `_link` exactly as the decrease arm."""
    assert ms._link([_pair(effect_size=0.7)], "GENE1")[1] == "inhibits"
    assert ms._link([_decrease(effect_size=-0.7)], "GENE1")[1] == "activates"


def test_a_links_strength_is_a_magnitude_with_no_sign_in_it() -> None:
    up = ms.rule_links(_row([_pair(effect_size=0.6)]), ms.EXTRACTOR_V2)
    down = ms.rule_links(_row([_decrease(effect_size=-0.6)]), ms.EXTRACTOR_V2)
    assert up[0][2] == down[0][2] == 0.6
    assert up[0][1] == "inhibits" and down[0][1] == "activates"


def test_a_strength_above_one_is_capped_as_v1_caps_it() -> None:
    assert ms.rule_links(_row([_pair(effect_size=3.5)]), ms.EXTRACTOR_V2)[0][2] == 1.0


def test_an_increase_link_found_only_in_heldout_pairs_keeps_the_heldout_split() -> None:
    row = _row([_pair(split=ms.HELDOUT)])
    assert ms.rule_links(row, ms.EXTRACTOR_V2)[0][4] == ms.HELDOUT


def test_a_training_increase_pair_sets_the_number_and_a_heldout_one_never_does() -> None:
    row = _row([_pair(effect_size=0.2), _pair(effect_size=0.9, split=ms.HELDOUT)])
    link = ms.rule_links(row, ms.EXTRACTOR_V2)[0]
    assert link[2] == 0.2 and link[4] == ms.TRAINING


def test_two_cells_that_measured_one_gene_are_two_links_as_r1_requires() -> None:
    row = _row([_pair(cell="K562"), _pair(cell="WTC11", effect_size=0.8)])
    links = sorted(ms.rule_links(row, ms.EXTRACTOR_V2), key=lambda x: x[3])
    assert [(x[0], x[3], x[2]) for x in links] == [("GENE1", "K562", 0.4), ("GENE1", "WTC11", 0.8)]


def test_a_non_significant_pair_raises_nothing_under_either_version() -> None:
    for outcome in (ms.NULL_INFORMATIVE, ms.NULL_INCONCLUSIVE, ms.MISSING):
        row = _row([_pair(outcome=outcome)])
        assert ms.rule_links(row) == []
        assert ms.rule_links(row, ms.EXTRACTOR_V2) == []


def test_the_predicate_reads_the_cached_outcome_field_and_not_the_effect_size() -> None:
    """A positive effect size with a non-significant outcome is not an increase-derived candidate."""
    row = _row([_pair(outcome=ms.NULL_INFORMATIVE, effect_size=0.9)])
    assert ms.rule_links(row, ms.EXTRACTOR_V2) == []
    assert il.INCREASE is inc.INCREASE is ms.INCREASE


# ---- additive: v2 is v1 plus the new links, by count and by content ------------------------------


def test_v2_is_v1_plus_the_new_links_by_count_and_by_content() -> None:
    row = _row(
        [
            _decrease(gene="A"),
            _pair(gene="B", effect_size=0.5),
            _decrease(gene="C", effect_size=-0.7),
            _pair(gene="D", effect_size=0.3, cell="WTC11"),
        ]
    )
    v1 = ms.rule_links(row)
    v2 = ms.rule_links(row, ms.EXTRACTOR_V2)
    assert v2[: len(v1)] == v1, "v1's links must survive unchanged and in order inside v2's output"
    assert len(v2) == len(v1) + 2
    assert all(link in v2 for link in v1), "no link v1 emitted may be dropped"
    assert {x[1] for x in v2[len(v1) :]} == {"inhibits"}


def test_an_element_with_no_crispri_block_has_no_links_under_either_version() -> None:
    row = {"id": "e1", "measured": {}}
    assert ms.rule_links(row) == [] and ms.rule_links(row, ms.EXTRACTOR_V2) == []
    assert ms.rule_link_conflicts(row) == []


# ---- the conflict rule ---------------------------------------------------------------------------


def test_a_gene_cell_with_both_signs_emits_no_rule_at_all() -> None:
    """Not an activates rule, not an inhibits rule, nothing. The sign is never chosen."""
    row = _row([_decrease(effect_size=-0.9), _pair(effect_size=0.2)])
    assert ms.rule_links(row, ms.EXTRACTOR_V2) == []
    # and the larger magnitude did not win, which is exactly what must not happen
    assert ms.rule_links(row, ms.EXTRACTOR_V2) != [("GENE1", "activates", 0.9, "K562", ms.TRAINING)]


def test_the_conflict_is_listed_with_both_observations_and_their_effect_sizes() -> None:
    row = _row([_decrease(effect_size=-0.9), _pair(effect_size=0.2)])
    conflicts = ms.rule_link_conflicts(row)
    assert len(conflicts) == 1
    c = conflicts[0]
    assert (c["gene"], c["cell"], c["emits"]) == ("GENE1", "K562", "no rule")
    assert c["element"] == {"chrom": "chr21", "start": 1_000, "end": 1_500}
    assert {o["outcome"] for o in c["observations"]} == {ms.DECREASE, ms.INCREASE}
    assert sorted(o["effect_size"] for o in c["observations"]) == [-0.9, 0.2]


def test_a_conflict_in_one_cell_does_not_suppress_a_clean_link_in_another() -> None:
    row = _row(
        [
            _decrease(cell="K562", effect_size=-0.9),
            _pair(cell="K562", effect_size=0.2),
            _pair(cell="WTC11", effect_size=0.5),
        ]
    )
    assert ms.rule_links(row, ms.EXTRACTOR_V2) == [("GENE1", "inhibits", 0.5, "WTC11", ms.TRAINING)]
    assert [c["cell"] for c in ms.rule_link_conflicts(row)] == ["K562"]


def test_v1_is_not_subject_to_the_conflict_rule_because_v1_cannot_see_the_increase() -> None:
    """The one place v2 can remove a link v1 emitted. v1 itself is untouched, which is the point."""
    row = _row([_decrease(effect_size=-0.9), _pair(effect_size=0.2)])
    assert ms.rule_links(row) == [("GENE1", "activates", 0.9, "K562", ms.TRAINING)]


def test_the_conflict_rule_is_registered_as_never_resolved_by_choosing_a_sign() -> None:
    assert "NEVER resolved by choosing a sign" in il.CONFLICT_RULE
    assert "larger magnitude" in il.CONFLICT_RULE
    assert "0 of 48" not in il.CONFLICT_IS_NOT_A_FINDING  # it is not quoted as a finding
    assert "not a finding" in il.CONFLICT_IS_NOT_A_FINDING


# ---- R2: binding in code, not only in prose ------------------------------------------------------


def test_an_increase_derived_link_records_the_measured_net_direction() -> None:
    [d] = ms.rule_links_detail(_row([_pair()]), ms.EXTRACTOR_V2)
    assert d["action"] == il.ACTION == "inhibits"
    assert d["evidence_status"] == il.EVIDENCE_STATUS == "observed"
    assert d["outcome"] == il.OUTCOME_TEXT == "increase on knockdown"
    assert d["derived_from"] == ms.INCREASE


def test_an_increase_derived_links_molecular_role_is_unresolved_and_not_an_empty_string() -> None:
    [d] = ms.rule_links_detail(_row([_pair()]), ms.EXTRACTOR_V2)
    assert d["molecular_role"] is None
    assert il.MOLECULAR_ROLE is None and il.MOLECULAR_ROLE_TEXT == "unresolved"


def test_a_planted_silencer_label_makes_the_check_fail() -> None:
    """The R2 test the approval requires. If this stops failing, the prohibition has been weakened."""
    [d] = ms.rule_links_detail(_row([_pair()]), ms.EXTRACTOR_V2)
    d["molecular_role"] = "silencer"
    with pytest.raises(ValueError, match="silencer"):
        il.check_no_mechanism_claim([d])


def test_a_planted_repressor_label_makes_the_check_fail() -> None:
    [d] = ms.rule_links_detail(_row([_pair()]), ms.EXTRACTOR_V2)
    d["molecular_role"] = "repressor"
    with pytest.raises(ValueError, match="repressor"):
        il.check_no_mechanism_claim([d])


def test_text_claiming_a_mechanism_makes_the_check_fail() -> None:
    [d] = ms.rule_links_detail(_row([_pair()]), ms.EXTRACTOR_V2)
    d["outcome"] = "the element represses the gene"
    with pytest.raises(ValueError, match="represses"):
        il.check_no_mechanism_claim([d])


def test_the_permitted_phrase_is_the_observation_and_passes() -> None:
    records = ms.rule_links_detail(_row([_pair()]), ms.EXTRACTOR_V2)
    il.check_no_mechanism_claim(records)  # must not raise
    assert records[0]["outcome"] == "increase on knockdown"


def test_the_detail_call_refuses_a_mechanism_claim_on_its_way_out() -> None:
    """The check is not optional: `rule_links_detail` runs it before returning."""
    row = _row([_pair()])
    assert ms.rule_links_detail(row, ms.EXTRACTOR_V2)  # the clean path returns records
    assert "check_no_mechanism_claim" in ms.rule_links_detail.__doc__


def test_the_forbidden_words_are_named_and_include_silencer_and_repressor() -> None:
    assert "silencer" in il.FORBIDDEN_OF_AN_INCREASE_LINK
    assert "repressor" in il.FORBIDDEN_OF_AN_INCREASE_LINK
    assert "never a silencer, never a repressor" in il.R2_WORDING


def test_a_decrease_derived_link_is_not_policed_by_the_increase_prohibition() -> None:
    """The prohibition is about what an increase may be called, and it applies to increases only."""
    records = ms.rule_links_detail(_row([_decrease()]), ms.EXTRACTOR_V2)
    assert records[0]["derived_from"] == ms.DECREASE
    assert records[0]["outcome"] == "decrease on knockdown"


# ---- the by-construction property, stated by the links themselves --------------------------------


def test_every_link_carries_the_by_construction_statement_of_itself() -> None:
    row = _row([_decrease(gene="A"), _pair(gene="B")])
    records = ms.rule_links_detail(row, ms.EXTRACTOR_V2)
    assert len(records) == 2
    for d in records:
        assert d["by_construction"] is il.BY_CONSTRUCTION
        assert "NOT evidence that the model covered" in d["by_construction"]
        assert "212 of" in d["by_construction"]


def test_the_by_construction_wording_is_the_one_already_on_the_record() -> None:
    assert il.DECREASE_ARM_IS_NOT_MODEL_COVERAGE.startswith("the 212-of-212")
    assert "EXTRACTOR'S OWN CONSTRUCTION" in il.DECREASE_ARM_IS_NOT_MODEL_COVERAGE


# ---- the floors are imports, and this lane applies neither ---------------------------------------


def test_both_floors_are_the_imported_ones_and_nobody_chose_them_here() -> None:
    assert il.POSITIVE_FLOOR is fresh.POSITIVE_FLOOR
    assert il.LOCUS_FLOOR is fresh.LOCUS_FLOOR
    assert fresh.LOCUS_FLOOR == cell2.POOLED_LOCUS_FLOOR
    assert (il.POSITIVE_FLOOR, il.LOCUS_FLOOR) == (30, 20)
    assert "no gate is taken by this lane" in il.FLOORS["not_applied_here"]


def test_the_registration_carries_the_inherited_readings_word_for_word() -> None:
    assert il.INHERITED_INCREASE_GO is inc.GO
    assert il.INHERITED_REPRESS2_NO_GO is rp.GATE_NO_GO
    assert il.STRUCTURAL_CAUSE is rp.MEASURED_AXIS_IS_NOT_THE_LINK_DIRECTION


def test_the_registration_says_v1_stays_the_default_and_why() -> None:
    reg = il.registration()["extractor_is_versioned_not_replaced"]
    assert reg["default"] == "v1"
    assert "5176" in reg["why_v1_stays_the_default"]
    assert "byte-identical" in reg["why_v1_stays_the_default"]


def test_the_registration_fixes_the_downstream_names_before_anything_claims_them() -> None:
    assert "_v2" in il.V2_CARRIES_NEW_NAMES
    assert "extractor: v2" in il.V2_CARRIES_NEW_NAMES


def test_the_registration_states_that_no_direction_is_read_here() -> None:
    reg = il.registration()
    assert any("no direction is read" in s for s in reg["what_this_lane_does_not_do"])
    assert any("21 ANSWERABLE links of 48 MEASURED" in s for s in reg["what_this_lane_does_not_do"])


def test_the_registration_result_on_disk_was_written_with_its_code_committed() -> None:
    p = Path("data/results/increase_links_registration.json")
    if not p.exists():
        pytest.skip(f"the registration result is not on this machine: {p}")
    d = json.loads(p.read_text())
    cc = d["result_manifest"]["code_cleanliness"]
    assert cc["own_code_is_committed"] is True
    assert cc["foreign_uncommitted_code_on_the_counting_path"] == []
    assert d["extractor_is_versioned_not_replaced"]["default"] == "v1"


# ---- the cached benchmark, read only if it is here -----------------------------------------------


def _benchmark_present() -> Path | None:
    for name in ms.CRISPRI_FILES:
        if not (crispri.KNOWLEDGE / name).exists():
            return None
    return crispri.KNOWLEDGE


def test_regulated_is_awarded_to_a_significant_decrease_only() -> None:
    """The fact the whole lane rests on, checked against the tables rather than taken on trust.

    If a `Regulated TRUE` row ever carried a non-negative effect size, v1's `inhibits` branch would
    already be reachable and the premise of this lane would be wrong.
    """
    if _benchmark_present() is None:
        pytest.skip(
            f"the CRISPRi benchmark tables are not on this machine: {crispri.KNOWLEDGE}/"
            f"{{{', '.join(ms.CRISPRI_FILES)}}}"
        )
    pairs, _invalid = ms.load_crispri(split=ms.ALL)
    regulated = [p for p in pairs if p.regulated]
    assert regulated, "the tables hold no regulated pair at all, so nothing here was tested"
    assert all(p.effect_size < 0 for p in regulated)
    assert all(p.outcome == ms.DECREASE for p in regulated)
    increases = [p for p in pairs if p.outcome == ms.INCREASE]
    assert increases and all(not p.regulated for p in increases)
