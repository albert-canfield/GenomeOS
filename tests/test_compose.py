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
