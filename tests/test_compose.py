"""Task 2.5: process-bigraph composition spike."""

import pytest

from genomeos.runtime import compose

SRC = open("data/demo/repressilator.bio").read()  # noqa: SIM115


@pytest.mark.skipif(not compose.AVAILABLE, reason="process-bigraph not installed (uv sync --extra compose)")
def test_network_and_ageing_run_on_one_clock():
    sim = compose.build_composite(SRC, cells=50)
    sim.run(24.0)  # 24 hours
    state = sim.state
    levels = state["levels"]
    assert any(v > 0 for v in levels.values())  # the repressilator has moved off zero
    assert abs(state["years"] - 24 / (24 * 365.25)) < 1e-9
    assert abs(state["global_time"] - 24.0) < 1e-9


@pytest.mark.skipif(not compose.AVAILABLE, reason="process-bigraph not installed (uv sync --extra compose)")
def test_body_process_grows_the_worm_embryo_on_the_composite_clock():
    sim = compose.build_body_composite(
        "data/organisms/celegans/embryo.bio", means=True
    )  # a path: imports resolve
    sim.run(14.0)  # 14 hours = 840 min of organism time, past hatching on the source axis
    st = sim.state
    assert st["minutes"] == 840.0 and st["global_time"] == 14.0
    assert st["cells"] == 610 and st["deaths"] == 110 and st["born"] == 1439
    assert st["fates"]["Neuron"] == 226 and st["fates"]["Intestine"] == 34 and st["fates"]["UNKNOWN"] == 0


@pytest.mark.skipif(not compose.AVAILABLE, reason="process-bigraph not installed (uv sync --extra compose)")
def test_network_gates_an_organism_through_a_factor():
    toy = """
module toy.gated
import bio.std.development
organism G { root: R; cell_type: Zygote; resolution: populations }
cell_type Ready { parent: Blastomere }
decision grow { action: divide; when: cell_type = Zygote, count = <64; fraction: 1.0; after: 30 min }
decision ready { action: differentiate; when: cell_type = Zygote, TetR_high = present; to: Ready }
"""
    on = compose.build_body_composite(toy, SRC, gates={"TetR_high": {"species": "TetR", "threshold": 1.0}})
    on.run(4.0)
    assert on.state["levels"]["TetR"] >= 1.0 and on.state["fates"]["Ready"] > 0
    off = compose.build_body_composite(toy, SRC, gates={"TetR_high": {"species": "TetR", "threshold": 1e6}})
    off.run(4.0)
    assert off.state["fates"]["Ready"] == 0 and off.state["fates"]["Zygote"] == 64
