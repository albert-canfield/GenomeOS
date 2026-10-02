"""The mesoderm diagnosis registration and its variant builder.

The registration is checked as committed code: a closed list that cannot grow, every entry carrying a
citation, every exclusion named. The variant builder is checked hermetically, on a planted
base text, so these tests say nothing about the state of the working checkout.
"""

from __future__ import annotations

import itertools
from math import comb

import pytest

from genomeos.lang import parse_file
from scripts.mesoderm_diagnosis import BOUND_LAST_FRACTION, BOUND_SHARE, REGISTRATION, variant_text

PLANTED = """\
module development.planted
gene A { max: 10; basal: 0.05; produces: Ap }
protein Ap { half_life: 2 }
rule Ap inhibits A { strength: 1.0; threshold: 2.0; hill: 3 }
rule Ap activates A { strength: 1.0; threshold: 2.0; hill: 3 }
"""


@pytest.fixture(scope="module")
def reg() -> dict:
    return REGISTRATION


def test_closed_list_size_and_pair_count(reg: dict) -> None:
    ids = [e["id"] for e in reg["closed_list"]]
    assert len(ids) == len(set(ids)) == reg["closed_list_n"] == 10
    assert reg["pairs"]["n"] == 10
    assert reg["pairs"]["n_choose_2"] == comb(10, 2) == 45
    assert reg["pairs"]["runnable"] + reg["pairs"]["not_runnable"] == 45
    assert len(list(itertools.combinations(ids, 2))) == 45


def test_every_entry_cites_a_source(reg: dict) -> None:
    for e in reg["closed_list"]:
        assert e["citation"].strip(), e["id"]
        assert e["acts_in_mesoderm_specification"].strip(), e["id"]
        assert any(c.isdigit() for c in e["citation"]), e["id"]


def test_every_entry_is_a_rule_a_species_or_a_declared_non_run(reg: dict) -> None:
    for e in reg["closed_list"]:
        if not e.get("runnable", True):
            assert e["no_run_reason"].strip(), e["id"]
            continue
        assert e.get("rules") or e.get("remove"), e["id"]


def test_exclusions_and_eliminations_are_named(reg: dict) -> None:
    assert [x["id"] for x in reg["excluded"]] == ["X1", "X2", "X3", "X4", "X5"]
    for x in reg["excluded"]:
        assert "EXCLUDED" in x["reason"], x["id"]
    (elim,) = reg["eliminated_in_advance"]
    assert elim["item"] == "the activator combination rule"
    assert elim["shas"] == ["93caf61", "267cc88", "a6fc5c4"]
    assert elim["withdrawn_attribution"] == "3f7b0f5"
    assert elim["runs_spent"] == 0


def test_bound_is_the_weak_one_and_not_the_census_range(reg: dict) -> None:
    assert reg["bound"]["share"] == BOUND_SHARE == 0.01
    assert reg["bound"]["last_fraction"] == BOUND_LAST_FRACTION == 0.1
    # 120 cells, so the bound needs two cells, not floating-point dust
    assert reg["bound"]["cells_needed"] == 2
    assert reg["bound"]["share"] < 0.694  # the CS7 mesoderm range is not the bound


def test_variant_text_comments_out_exactly_the_named_rule() -> None:
    out = variant_text(PLANTED, [{"remove": ["rule Ap inhibits A"]}])
    assert "# removed for this run: rule Ap inhibits A" in out
    assert "\nrule Ap activates A" in out  # the other rule, whose prefix differs, is untouched


def test_variant_text_refuses_a_removal_that_does_not_match_exactly_once() -> None:
    with pytest.raises(ValueError, match="matched 0 lines"):
        variant_text(PLANTED, [{"remove": ["rule Ap represses A"]}])
    doubled = PLANTED + "rule Ap inhibits A { strength: 1.0; threshold: 2.0; hill: 3 }\n"
    with pytest.raises(ValueError, match="matched 2 lines"):
        variant_text(doubled, [{"remove": ["rule Ap inhibits A"]}])


def test_a_species_two_entries_both_add_is_declared_once() -> None:
    sp = {"gene": "WNT3A", "protein": "Wnt3a", "note": "planted"}
    one = {"species": [sp], "rules": [{"src": "Wnt3a", "dst": "A", "action": "activates", "cite": "x"}]}
    two = {"species": [sp], "rules": [{"src": "Ap", "dst": "WNT3A", "action": "activates", "cite": "y"}]}
    out = variant_text(PLANTED, [one, two])
    assert out.count("gene WNT3A ") == 1
    assert out.count("protein Wnt3a ") == 1
    assert "rule Wnt3a activates A" in out and "rule Ap activates WNT3A" in out


def test_no_variant_reverts_the_sign_correction(reg: dict) -> None:
    """The correction at 5cbce26 is not negotiable: no entry removes or flips it."""
    for e in reg["closed_list"]:
        for target in e.get("remove", []):
            assert "SOX17" not in target or "Sox17" in target  # only Sox17->TBXT, never Tbxt->SOX17
            assert target != "rule Tbxt activates SOX17"
        for rule in e.get("rules", []):
            assert not (rule["src"] == "Tbxt" and rule["dst"] == "SOX17")


def test_the_module_still_states_the_corrected_sign() -> None:
    module = parse_file("data/demo/gastrulation.bio")
    signs = {(r.source, r.target): r.action.name for r in module.rules}
    assert signs[("Tbxt", "SOX17")] == "ACTIVATE"
    assert signs[("Sox17", "TBXT")] == "INHIBIT"
