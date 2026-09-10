"""Task 2.4: Boolean network engine reproduces a published logical model."""

from genomeos.ir import Action
from genomeos.runtime.boolean import BooleanNetwork


def test_boolean_repressilator_cycles():
    net = BooleanNetwork.from_bnet("A, !C\nB, !A\nC, !B\n")
    att = net.attractor({"A": True, "B": False, "C": False})
    assert len(att) == 6  # the 3-node ring has a period-6 synchronous cycle
    assert len(net.attractors(samples=50)) >= 1


def test_faure_mammalian_cell_cycle():
    net = BooleanNetwork.from_file("data/models/mammalian_cell_cycle.bnet")
    assert len(net.nodes) == 10
    # no growth factor: quiescent G1 fixed point with Rb, p27 and Cdh1 on
    rest = net.attractor({"CycD": False})
    assert len(rest) == 1
    s = rest[0]
    assert s["Rb"] and s["p27"] and s["Cdh1"] and not s["CycB"] and not s["CycA"]
    # growth factor present: a cyclic attractor (the cell cycle), never a fixed point
    cyc = net.attractor({"CycD": True})
    assert len(cyc) > 1
    assert any(x["CycB"] for x in cyc) and any(x["CycE"] for x in cyc)  # both cyclins take their turn
    atts = net.attractors(samples=300)
    assert any(len(a) == 1 for a in atts) and any(len(a) > 1 for a in atts)


def test_boolean_to_bioir():
    net = BooleanNetwork.from_bnet("A, !C\nB, A & !B\n")
    m = net.to_module("bool.test")
    acts = {(r.source, r.target): r.action for r in m.rules}
    assert acts[("C", "A")] is Action.INHIBIT
    assert acts[("A", "B")] is Action.ACTIVATE
    assert acts[("B", "B")] is Action.INHIBIT
    assert m.confidence_report()["node"] == 0.5
