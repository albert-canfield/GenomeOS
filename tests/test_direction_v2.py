# SPDX-License-Identifier: AGPL-3.0-or-later
"""The registration of the versioned direction rule: the floor is imported, the default is v1."""

from __future__ import annotations

import pytest

from genomeos.attribution import direction_v2 as dv
from genomeos.attribution import increase_links as il
from genomeos.predict import enhancer_target as et


def test_default_version_is_v1() -> None:
    assert dv.DIRECTION_RULE_DEFAULT == dv.V1
    assert (dv.V1, dv.V2) == (1, 2)


def test_magnitude_floor_is_the_imported_object() -> None:
    """The threshold is `predict.enhancer_target.MIN_EFFECT` itself, not a copied literal."""
    assert dv.MAGNITUDE_FLOOR is et.MIN_EFFECT


def test_cells_are_the_imported_tuple() -> None:
    """The cell names v2 can read a second value for are the project's, not a restated list."""
    assert dv.CELLS is et.CELLS


def _word_lists(ns: dict[str, object]) -> list[str]:
    """Every module-level collection of strings that overlaps the increase-link prohibition.

    The prohibition lives in `increase_links.check_no_mechanism_claim` and its own tuple. A second
    tuple here would shadow it: a reader fixing one would leave the other, so none may exist.
    """
    out = []
    forbidden = set(il.FORBIDDEN_OF_AN_INCREASE_LINK)
    for name, value in ns.items():
        if not isinstance(value, (tuple, list, set, frozenset)):
            continue
        words = {v.lower() for v in value if isinstance(v, str)}
        if words & forbidden:
            out.append(name)
    return out


def test_module_defines_no_word_list_that_could_shadow_the_checker() -> None:
    assert _word_lists(vars(dv)) == []


def test_the_shadow_scan_would_catch_one_planted() -> None:
    """The counterfactual: the scan above is not vacuous, it flags a planted list."""
    planted = {"FORBIDDEN": ("silencer", "repressor")}
    assert _word_lists(planted) == ["FORBIDDEN"]


def test_the_imported_mechanism_check_is_live() -> None:
    """A planted increase-derived record is rejected by the imported checker, not by a local list."""
    with pytest.raises(ValueError):
        il.check_no_mechanism_claim(
            [{"derived_from": il.INCREASE, "gene": "G", "cell": "K562", "note": "a silencer"}]
        )
    il.check_no_mechanism_claim([{"derived_from": il.INCREASE, "gene": "G", "cell": "K562"}])


def test_the_registration_states_what_may_not_be_claimed() -> None:
    assert "never independent evidence" in dv.NOT_VALIDATION
    assert "not an improvement in accuracy" in dv.NOT_AN_IMPROVEMENT
    assert "the sign rule did not resolve" in dv.UNRESOLVED_MEANS
    assert "never means the element has no action" in dv.UNRESOLVED_MEANS
    assert "0 by construction" in dv.FLIPS_ARE_STRUCTURALLY_IMPOSSIBLE
    assert "not detected open at the reader's registered call, never closed" in dv.WORDING


# ---- v1: the committed rule, unchanged -----------------------------------------------------------


def test_v1_is_the_committed_expression() -> None:
    """Exactly what compile.py held: `activates` iff the prediction says so, `inhibits` otherwise."""
    for action, expected in (("activates", dv.ACTIVATES), ("represses", dv.INHIBITS), (None, dv.INHIBITS)):
        pc = {"action": action, "tissue": "K562", "gene": "G"}
        assert dv.compiled_action(pc) == expected
        assert dv.compiled_action(pc, dv.V1) == expected
        assert dv.direction_call(pc, dv.V1).action == expected


def test_v1_opens_no_file() -> None:
    """The default path reads `pc` and nothing else: no chrom, no element id, no cache."""
    assert dv.direction_call({"action": "activates"}, dv.V1).values == ()


def test_an_unknown_version_is_refused() -> None:
    with pytest.raises(ValueError):
        dv.direction_call({"action": "activates"}, 3)


# ---- the unresolved value is first class ----------------------------------------------------------


def _unresolved_call() -> dv.DirectionCall:
    return dv.direction_call({"action": "represses", "tissue": ""}, dv.V2)


def test_unresolved_is_neither_direction() -> None:
    call = _unresolved_call()
    assert call.resolved is False
    assert call.token is dv.UNRESOLVED
    assert call.token != dv.ACTIVATES
    assert call.token != dv.INHIBITS
    assert call.token not in (dv.ACTIVATES, dv.INHIBITS)
    assert call.reason == "context_unknown"


def test_reading_a_direction_off_an_unresolved_call_raises() -> None:
    with pytest.raises(dv.UnresolvedDirectionError):
        _ = _unresolved_call().action
    with pytest.raises(dv.UnresolvedDirectionError):
        dv.compiled_action({"action": "represses", "tissue": ""}, dv.V2)


def test_an_unresolved_call_says_molecular_action_unknown() -> None:
    assert _unresolved_call().molecular_action == dv.MOLECULAR_ACTION_UNKNOWN
    assert "did not resolve" in dv.MOLECULAR_ACTION_UNKNOWN


def test_the_unresolved_value_cannot_be_read_as_a_boolean() -> None:
    """`if token:` would read it as something; it raises instead."""
    with pytest.raises(dv.UnresolvedDirectionError):
        bool(dv.UNRESOLVED)


def test_the_planted_fall_through_a_naive_reader_would_take() -> None:
    """The counterfactual for UNRESOLVED_IS_FIRST_CLASS, planted rather than argued.

    A reader written as `activates if token == ACTIVATES else inhibits` - the shape compile.py used -
    would quietly call an unresolved rule `inhibits`. The test asserts that this is what such a
    reader does, and that the supported accessor raises instead, so the two cannot be confused.
    """
    call = _unresolved_call()
    naive = dv.ACTIVATES if call.token == dv.ACTIVATES else dv.INHIBITS
    assert naive == dv.INHIBITS  # the defect, if anyone writes it this way again
    with pytest.raises(dv.UnresolvedDirectionError):
        _ = call.action


# ---- v2: the sign in the rule's own assigned cell --------------------------------------------------


def _row(**kw: object) -> dict:
    row = {
        "gene": "G",
        "max_drop_log2fc": 0.0,
        "max_drop_tissue": "",
        "max_rise_log2fc": 0.0,
        "max_rise_tissue": "",
        "by_cell": {},
    }
    row.update(kw)
    return row


def _call(row: dict, tissue: str = "K562", action: str = "represses") -> dv.DirectionCall:
    return dv.direction_call({"action": action, "tissue": tissue, "gene": "G"}, dv.V2, row=row)


def test_v2_resolves_a_rise_in_the_cell_to_inhibits() -> None:
    call = _call(_row(max_rise_log2fc=0.5, max_rise_tissue="K562", by_cell={"K562": 0.3}))
    assert call.resolved and call.action == dv.INHIBITS
    assert call.values == (("max_rise", 0.5), ("by_cell", 0.3))


def test_v2_resolves_a_fall_in_the_cell_to_activates() -> None:
    call = _call(_row(max_drop_log2fc=-0.5, max_drop_tissue="K562", by_cell={"K562": -0.3}))
    assert call.resolved and call.action == dv.ACTIVATES


def test_v2_is_unresolved_when_the_cell_points_both_ways() -> None:
    call = _call(_row(max_rise_log2fc=0.5, max_rise_tissue="K562", by_cell={"K562": -0.3}))
    assert not call.resolved and call.reason == "signs_disagree_in_cell"


def test_v2_is_unresolved_when_one_value_only_is_retained() -> None:
    """The defect itself: one track, which is what v1 resolves on."""
    call = _call(_row(max_rise_log2fc=0.5, max_rise_tissue="occipital lobe"), tissue="occipital lobe")
    assert not call.resolved and call.reason == "one_value_only"


def test_v2_is_unresolved_when_two_fields_carry_one_track() -> None:
    """AMENDMENT_1: by_cell can be the very track the maximum selected."""
    call = _call(_row(max_rise_log2fc=0.1453, max_rise_tissue="K562", by_cell={"K562": 0.1453}))
    assert not call.resolved and call.reason == "one_track_seen_twice"


def test_v2_is_unresolved_below_the_imported_floor() -> None:
    call = _call(_row(max_rise_log2fc=0.09, max_rise_tissue="K562", by_cell={"K562": 0.05}))
    assert not call.resolved and call.reason == "below_magnitude_floor"


def test_v2_resolves_at_the_floor_exactly() -> None:
    call = _call(_row(max_rise_log2fc=dv.MAGNITUDE_FLOOR, max_rise_tissue="K562", by_cell={"K562": 0.05}))
    assert call.resolved and call.action == dv.INHIBITS


def test_v2_is_unresolved_on_a_zero_value() -> None:
    call = _call(_row(max_rise_log2fc=0.5, max_rise_tissue="K562", by_cell={"K562": 0.0}))
    assert not call.resolved and call.reason == "zero_value_in_cell"


def test_v2_is_unresolved_when_the_cell_has_no_retained_value() -> None:
    call = _call(_row(max_rise_log2fc=0.5, max_rise_tissue="A375"), tissue="HepG2")
    assert not call.resolved and call.reason == "no_value_for_cell"


def test_v2_is_unresolved_without_a_row_for_the_target() -> None:
    call = dv.direction_call({"action": "represses", "tissue": "K562", "gene": "G"}, dv.V2)
    assert not call.resolved and call.reason == "no_row_for_target"


def test_a_placeholder_zero_belongs_to_no_cell() -> None:
    """A gene row carries 0.0 with an empty tissue for a direction no track took; it is not read."""
    assert dv.retained_values(_row(by_cell={"K562": 0.2}), "") == ()
    assert dv.retained_values(_row(), "K562") == ()


def test_v2_reads_no_other_cell() -> None:
    """The tracks of other cells are not evidence about this rule's action in its own cell."""
    row = _row(max_rise_log2fc=0.5, max_rise_tissue="K562", by_cell={"K562": 0.3, "HepG2": -0.9})
    assert dv.retained_values(row, "K562") == (("max_rise", 0.5), ("by_cell", 0.3))


def test_every_reason_the_code_records_is_registered() -> None:
    rows = [
        ({"action": "represses", "tissue": ""}, None),
        ({"action": "represses", "tissue": "K562", "gene": "G"}, None),
        ({"action": "represses", "tissue": "K562", "gene": "G"}, _row(max_rise_log2fc=0.5)),
        (
            {"action": "represses", "tissue": "K562", "gene": "G"},
            _row(max_rise_log2fc=0.5, max_rise_tissue="K562"),
        ),
        (
            {"action": "represses", "tissue": "K562", "gene": "G"},
            _row(max_rise_log2fc=0.4, max_rise_tissue="K562", by_cell={"K562": 0.4}),
        ),
        (
            {"action": "represses", "tissue": "K562", "gene": "G"},
            _row(max_rise_log2fc=0.4, max_rise_tissue="K562", by_cell={"K562": 0.0}),
        ),
        (
            {"action": "represses", "tissue": "K562", "gene": "G"},
            _row(max_rise_log2fc=0.4, max_rise_tissue="K562", by_cell={"K562": -0.2}),
        ),
        (
            {"action": "represses", "tissue": "K562", "gene": "G"},
            _row(max_rise_log2fc=0.09, max_rise_tissue="K562", by_cell={"K562": 0.02}),
        ),
    ]
    seen = set()
    for pc, row in rows:
        call = dv.direction_call(pc, dv.V2, row=row) if row else dv.direction_call(pc, dv.V2)
        assert call.reason in dv.UNRESOLVED_REASONS
        seen.add(call.reason)
    assert seen == set(dv.UNRESOLVED_REASONS)


# ---- the seam in the compiler ----------------------------------------------------------------------


def _chry_elements() -> list[dict]:
    from genomeos.attribution import compile as cp

    try:
        return cp._attributed("chrY")
    except FileNotFoundError:  # pragma: no cover - a machine without the attribution results
        return []


def test_the_compiler_default_is_v1_and_writes_the_same_bytes() -> None:
    """v1 is the default and naming it changes nothing: the same program, byte for byte."""
    from genomeos.attribution import compile as cp

    if not _chry_elements():
        pytest.skip("no attributed elements for chrY on this machine")
    assert cp.compile_chromosome("chrY") == cp.compile_chromosome("chrY", direction_version=dv.V1)


def test_the_compiler_refuses_to_write_a_direction_that_did_not_resolve() -> None:
    """V2_NOT_WRITABLE in force: the grammar has no token for it, so the compiler stops.

    Writing `inhibits` or `activates` for an unresolved call would be the defect v2 exists to stop,
    so there is deliberately no way to compile one.
    """
    from genomeos.attribution import compile as cp

    elements = _chry_elements()
    if not elements:
        pytest.skip("no attributed elements for chrY on this machine")
    unresolved = [
        e
        for e in elements
        if not dv.direction_call(e["predicted_coding"], dv.V2, chrom="chrY", element_id=e["id"]).resolved
    ]
    if not unresolved:
        pytest.skip("every chrY rule resolves under v2 on this machine")
    with pytest.raises(dv.UnresolvedDirectionError):
        cp.compile_chromosome("chrY", direction_version=dv.V2)
