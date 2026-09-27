"""Area H: a network knockout's consequence, ranked against matched random knockouts."""

import pytest

from genomeos.forge.network_experiment import (
    CONTROL_DRAWS,
    NETWORK_KNOCKOUT_KNOWN_CASE,
    _matched_sets,
    parse_knockouts,
    run,
    run_experiment,
)
from genomeos.ir import Experiment
from genomeos.lang import parse
from genomeos.runtime.boolean import BooleanNetwork

FAURE = "data/models/mammalian_cell_cycle.bnet"
CYCLING = {"CycD": True}  # the fixture's start: growth factor on, every other node off

RING = """
module net.knockout.ring
gene a { max: 200; basal: 0.2; produces: A }
gene b { max: 200; basal: 0.2; produces: B }
gene c { max: 200; basal: 0.2; produces: C }
protein A {}
protein B {}
protein C {}
rule A inhibits b { strength: 1.0; threshold: 1.0; hill: 2 }
rule B inhibits c { strength: 1.0; threshold: 1.0; hill: 2 }
rule C inhibits a { strength: 1.0; threshold: 1.0; hill: 2 }
param mrna_half_life = 0.693 h {}
param protein_half_life = 0.139 h {}
param translation_rate = 5 {}
"""


def _faure() -> BooleanNetwork:
    return BooleanNetwork.from_file(FAURE)


def test_the_known_case_cycd_off_turns_the_cycle_into_the_g1_fixed_point():
    assert NETWORK_KNOCKOUT_KNOWN_CASE["knockout"] == "CycD"
    assert "Faure" in NETWORK_KNOCKOUT_KNOWN_CASE["source"]
    r = run(_faure(), ["CycD"], start=CYCLING)
    assert r.unperturbed.startswith("cycle") and r.perturbed == "fixed point" and r.changed_kind
    on = {n for n, v in r.levels.items() if v == 1.0}
    assert {"Rb", "p27", "Cdh1"} <= on
    assert r.levels["CycA"] == 0.0 and r.levels["CycB"] == 0.0 and r.levels["CycE"] == 0.0
    assert r.consequence > 0.5  # over half the scale, averaged over the nine nodes read in both runs
    assert "CycD" not in r.per_node  # the knocked-out node is not part of its own consequence
    assert r.p_value == pytest.approx(0.1) and r.exhaustive and r.draws == 9  # the floor at nine draws


def test_the_falsifier_an_unread_node_scores_zero_at_p_one():
    net = BooleanNetwork.from_bnet("A, !C\nB, !A\nC, !B\nSpectator, A\n")  # no rule reads Spectator
    r = run(net, ["Spectator"], start={"A": True})
    assert r.consequence == 0.0 and r.p_value == 1.0
    assert r.per_node == {"A": 0.0, "B": 0.0, "C": 0.0} and r.unperturbed == r.perturbed
    assert not r.changed_kind


def test_the_control_never_draws_a_chosen_node_or_edge_and_matches_the_size():
    net = _faure()
    r = run(net, ["CycD", "Cdh1>CycB"], start=CYCLING, draws=40)
    free_nodes = [n for n in net.nodes if n != "CycD"]
    free_edges = [(s, t) for t in net.nodes for s in sorted(net.inputs_of(t)) if (s, t) != ("Cdh1", "CycB")]
    matched, exhaustive = _matched_sets(free_nodes, free_edges, 1, 1, 40, 0)
    assert len(matched) == 40 and not exhaustive  # far more matched knockouts exist, so this is a sample
    for nodes, edges in matched:
        assert len(nodes) == 1 and len(edges) == 1
        assert "CycD" not in nodes and ("Cdh1", "CycB") not in edges
    ge = sum(1 for c in r.control if c >= r.consequence - 1e-12)
    assert r.p_value == (1 + ge) / (1 + r.draws) and r.draws == 40


def test_every_matched_knockout_is_run_when_few_exist_and_the_registered_default_is_200():
    assert CONTROL_DRAWS == 200
    net = _faure()
    single = run(net, ["Rb"], start=CYCLING)
    assert single.exhaustive and single.draws == 9  # nine other nodes: every matched knockout was run
    pair = run(net, ["Rb", "CycA"], start=CYCLING, draws=12)
    assert not pair.exhaustive and pair.draws == 12


def test_an_edge_knockout_only_blinds_its_target():
    """Cdh1>CycB stops CycB reading Cdh1; CycA and UbcH10 still read it, so it is not a node knockout."""
    net = _faure()
    edge = run(net, ["Cdh1>CycB"], start=CYCLING)
    node = run(net, ["Cdh1"], start=CYCLING)
    assert edge.consequence != node.consequence
    assert "Cdh1" in edge.per_node and "Cdh1" not in node.per_node
    assert edge.levels["Cdh1"] > 0.0  # the edge's source still runs its own rule
    assert node.levels["Cdh1"] == 0.0


def test_a_knockout_that_names_nothing_is_refused():
    net = _faure()
    with pytest.raises(ValueError, match="names nothing"):
        run(net, ["NoSuchNode"], start=CYCLING)
    with pytest.raises(ValueError, match="names nothing"):
        run(net, ["CycD>CycB"], start=CYCLING)  # CycB's rule does not read CycD


def test_knockouts_are_read_as_nodes_and_edges():
    assert parse_knockouts(["A", "B>C", " A ", "B > C", ""]) == (["A"], [("B", "C")])


def test_an_experiment_block_knockout_list_runs_the_same_network_knockout():
    ex = Experiment(name="cycd_off", knockouts=["CycD"])
    a = run_experiment(_faure(), ex, start=CYCLING)
    b = run(_faure(), ["CycD"], start=CYCLING)
    assert a.to_dict() == b.to_dict()


def test_a_continuous_node_is_held_at_zero_all_run_and_a_cut_edge_stops_the_oscillation():
    m = parse(RING)
    off = run(m, ["a"], start={"A": 10.0}, draws=5, hours=60.0)
    assert off.levels["a"] == 0.0  # exactly zero: inside the integrator's sub-steps, not only between
    assert off.baseline["a"] > 1.0 and off.unperturbed == "oscillating" and off.perturbed == "steady"
    assert off.changed_kind and off.consequence > 1.0
    cut = run(m, ["A>b"], start={"A": 10.0}, draws=5, hours=60.0)
    assert cut.changed_kind and cut.perturbed == "steady"
    assert cut.levels["a"] > 1.0  # gene a still runs; only b stopped reading A
    assert cut.levels["b"] > cut.baseline["b"]  # b lost its repressor
