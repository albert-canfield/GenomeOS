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
