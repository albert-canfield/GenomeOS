# SPDX-License-Identifier: Apache-2.0
"""The MaBoSS adapter (runtime/maboss.py): chosen by name, a translation that refuses what it cannot
carry, and an honest refusal when the `maboss` extra is absent. ROADMAP §5 item 13, area C."""

import sys
from itertools import product

import pytest

from genomeos.runtime import maboss
from genomeos.runtime.boolean import ENGINES, BooleanNetwork, runtime

# An engine test (it imports only the engine), so it travels with the biolang package, whose own suite
# accepts a skip only when the reason says the extra is not installed; the compose tests do the same.
# tests/test_maboss_score.py pins, under requires_extra("maboss"), that this flag is True wherever the
# extra is installed, so these cannot skip quietly in CI.
needs_maboss = pytest.mark.skipif(
    not maboss.available(), reason="cmaboss not installed (uv sync --extra maboss)"
)

FAURE = "data/models/mammalian_cell_cycle.bnet"
RING = "A, !C\nB, A & !B\nC, (B | A) & !A\nD, 1\n"


def _same_rules(a: BooleanNetwork, b: BooleanNetwork) -> bool:
    """Same nodes and the same truth table for every rule."""
    if set(a.rules) != set(b.rules):
        return False
    for bits in product((False, True), repeat=len(a.rules)):
        state = dict(zip(a.rules, bits, strict=True))
        if a.step(state) != b.step(state):
            return False
    return True


def test_the_stand_in_stays_the_default_and_engines_are_named():
    net = BooleanNetwork.from_bnet(RING)
    assert ENGINES[0] == "standin" and "maboss" in ENGINES
    assert runtime(net) is net
    with pytest.raises(ValueError, match="unknown Boolean engine"):
        runtime(net, "boolnet")
    with pytest.raises(TypeError, match="no simulation settings"):
        runtime(net, sample_count=10)


def test_a_missing_extra_is_reported_never_substituted(monkeypatch):
    """Without cmaboss the adapter says which extra to install; the stand-in is not run in its place."""
    monkeypatch.setitem(sys.modules, "cmaboss", None)  # importing it now raises ImportError
    assert maboss.available() is False
    net = BooleanNetwork.from_bnet(RING)
    with pytest.raises(maboss.EngineUnavailableError, match=r"uv sync --extra maboss"):
        runtime(net, "maboss")
    with pytest.raises(ImportError, match=r"--extra maboss"):
        maboss.MaBoSSRuntime(maboss.MaBoSSModel.from_network(net))


def test_the_subset_translates_and_the_rest_is_refused():
    assert maboss.logic_of("A & !B | C") == "((A & !B) | C)"
    assert maboss.logic_of("!(A | B)") == "!(A | B)"
    assert maboss.logic_of("1") == "1" and maboss.logic_of("0") == "0"
    for bad in ("2", "A & ", "node & A"):
        with pytest.raises(maboss.TranslationError):
            maboss.logic_of(bad)
    net = BooleanNetwork.from_bnet(RING)
    with pytest.raises(maboss.TranslationError, match="not nodes"):
        maboss.MaBoSSModel.from_network(net, initial={"Z": True})
    with pytest.raises(maboss.TranslationError, match="both fixed and random"):
        maboss.MaBoSSModel.from_network(net, initial={"A": True}, random_nodes=["A"])
    odd = BooleanNetwork()
    odd.add("logic", "logic")
    with pytest.raises(maboss.TranslationError, match="not one MaBoSS can read"):
        maboss.MaBoSSModel.from_network(odd)


def test_a_network_goes_to_maboss_and_back_unchanged():
    net = BooleanNetwork.from_file(FAURE)
    model = maboss.MaBoSSModel.from_network(net, initial={"CycD": True, "Rb": True}, random_nodes=["p27"])
    back, initial, random_nodes, rate = maboss.to_network(model)
    assert _same_rules(net, back)
    assert initial["CycD"] and initial["Rb"] and not initial["E2F"] and "p27" not in initial
    assert random_nodes == ["p27"] and rate == 1.0
    small = BooleanNetwork.from_bnet(RING)
    back, _, _, _ = maboss.to_network(maboss.MaBoSSModel.from_network(small))
    assert _same_rules(small, back)


def test_to_network_refuses_where_maboss_and_the_stand_in_differ():
    cfg = "A.istate = 0; B.istate = 1;\n"
    ok = "Node A { logic = B; }\nNode B { logic = !A; }\n"
    assert maboss.to_network(maboss.MaBoSSModel(ok, cfg))[3] == 1.0
    cases = {
        "unequal rates": ("Node A { logic = B; rate_up = @logic ? 2 : 0; }\nNode B { logic = !A; }\n", cfg),
        "is not @logic": ("Node A { logic = B; rate_up = B ? 1 : 0; }\nNode B { logic = !A; }\n", cfg),
        "fields not translated": ("Node A { logic = B; max_rate = 3; }\nNode B { logic = !A; }\n", cfg),
        "outside the subset": ("Node A { logic = B ^ A; }\nNode B { logic = !A; }\n", cfg),
        "joint initial": (ok, "[A, B].istate = 0.5 [0, 1], 0.5 [1, 0];\n"),
        "neither 0, 1 nor random": (ok, "A.istate = 0.3; B.istate = 1;\n"),
        "move one way only": (
            "Node A { logic = B; rate_up = 0; }\nNode B { logic = !A; }\n",
            "B.istate = 1;\n",
        ),
    }
    for reason, (bnd, c) in cases.items():
        with pytest.raises(maboss.TranslationError, match=reason):
            maboss.to_network(maboss.MaBoSSModel(bnd, c))


def test_a_held_node_has_no_rates_and_starts_at_its_value():
    model = maboss.MaBoSSModel.from_network(BooleanNetwork.from_bnet(RING), random_nodes=["A"])
    held = model.held({"A": True, "C": False})
    fields = maboss.parse_bnd(held.bnd)
    assert fields["A"]["rate_up"] == "0" and fields["A"]["rate_down"] == "0"
    assert fields["A"]["logic"] == maboss.parse_bnd(model.bnd)["A"]["logic"]
    assert fields["B"] == maboss.parse_bnd(model.bnd)["B"]
    cfg = maboss.parse_cfg(held.cfg)
    assert cfg["istate"]["A"] == "1" and cfg["istate"]["C"] == "0" and not cfg["joint_istate"]
    net, initial, random_nodes, _ = maboss.to_network(held)
    assert net.rules["A"] == "1" and net.rules["C"] == "0" and initial["A"] and random_nodes == []
    with pytest.raises(maboss.TranslationError, match="not nodes"):
        model.held({"Z": True})


def test_the_istate_form_that_crashes_cmaboss_is_refused_before_it_reaches_it():
    with pytest.raises(maboss.TranslationError, match="crashes"):
        maboss.MaBoSSModel("Node A { logic = A; }\n", "A.istate = 0.5 [0], 0.5 [1];\n")
    maboss.MaBoSSModel("Node A { logic = A; }\n", "[A].istate = 0.5 [0], 0.5 [1];\n")


@needs_maboss
def test_the_faure_cell_cycle_on_maboss_matches_the_stand_in_fixed_point():
    """No growth factor: one quiescent fixed point (Rb, p27, Cdh1), reached with probability 1 on
    MaBoSS, as the stand-in's synchronous attractor finds it. With CycD on, the cyclins keep
    cycling: no single output state holds the population at t = 50."""
    net = BooleanNetwork.from_file(FAURE)
    rest = runtime(net, "maboss", initial={"CycD": False}, sample_count=2000, max_time=50).run()
    fixed = net.attractor({"CycD": False})
    assert len(fixed) == 1 and rest.states == {
        tuple(sorted(n for n, v in fixed[0].items() if v)): pytest.approx(1.0)
    }
    cycling = runtime(
        net, "maboss", initial={"CycD": True, "Rb": True, "p27": True, "Cdh1": True}, sample_count=2000
    ).run()
    assert max(cycling.states.values()) < 0.5 and cycling.nodes["CycD"] == pytest.approx(1.0)
    assert 0.1 < cycling.nodes["CycB"] < 0.9 and cycling.settings["cmaboss"] == maboss.version()


@needs_maboss
def test_a_held_node_never_moves_on_maboss():
    model = maboss.MaBoSSModel.from_network(
        BooleanNetwork.from_bnet("A, !B\nB, A\n"), random_nodes=["A", "B"], sample_count=500, max_time=5
    )
    assert maboss.MaBoSSRuntime(model.held({"A": True})).run().nodes["A"] == pytest.approx(1.0)
    assert maboss.MaBoSSRuntime(model.held({"A": False})).run().nodes["A"] == pytest.approx(0.0)
