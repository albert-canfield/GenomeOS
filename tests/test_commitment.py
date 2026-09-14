"""Commitment, competence and lateral inhibition: the statistics, their nulls and the reference model."""

import json
import random
from collections import Counter
from pathlib import Path

import pytest

from genomeos.organism import commitment as cm
from genomeos.organism import lateral

RESULT = Path("data/results/celegans_commitment.json")


def test_curveball_keeps_both_margins():
    rng = random.Random(1)
    rows = [0b1011, 0b0110, 0b1100, 0b0001, 0b1111]
    out = cm.curveball(rows, rng, trades=200)
    assert [r.bit_count() for r in out] == [r.bit_count() for r in rows]
    cols = lambda rs: [sum(r >> k & 1 for r in rs) for k in range(4)]  # noqa: E731
    assert cols(out) == cols(rows) and out != rows


def test_programmes_mutual_information_and_bins():
    assert cm.programmes_on({"ELT-2", "HLH-1", "X"}) == {"intestine", "muscle"}
    assert cm.time_bin(0) == 0 and cm.time_bin(45) == 1 and cm.time_bin(199) == 4 and cm.time_bin(500) == 4
    # time independent of fate within each stratum: no information; time determines fate: one bit
    same = [("s", "0", {"gut": 0.5, "skin": 0.5}), ("s", "1", {"gut": 0.5, "skin": 0.5})]
    assert cm._soft_cmi(same) == 0.0
    split = [("s", "0", {"gut": 1.0}), ("s", "1", {"skin": 1.0})]
    assert abs(cm._soft_cmi(split) - 1.0) < 1e-9
    assert cm._axis("ABala", "ABalp") == "a/p" and cm._axis("ABal", "ABar") == "l/r"
    assert cm._tv(Counter(neuron=2), Counter(neuron=1, glia=1)) == 0.5


def test_lateral_inhibition_needs_noise_to_break_symmetry():
    equal = lateral.two_cells(noise=0.0)
    assert not equal.diverged and equal.winner is None and equal.delta[0] == equal.delta[1]
    tilted = lateral.two_cells(start=(0.5, 0.51))
    assert tilted.diverged and tilted.winner == 1  # a small head start is amplified, not erased
    stats = lateral.selection_statistics(noise=0.02, runs=60)
    assert not stats["deterministic_equal_start_diverges"] and stats["diverged"] >= 55
    assert 0.25 < stats["first_cell_win_share"] < 0.75  # the winner is chosen by the noise


@pytest.mark.skipif(not RESULT.exists(), reason="run scripts/celegans_commitment.py")
def test_distilled_commitment_result():
    r = json.loads(RESULT.read_text())
    assert r["cells_with_complete_lifetimes"] > 700
    st = r["states"]
    assert st["series"]["factors_per_cell"][0] < st["series"]["factors_per_cell"][-1]
    assert set(st["p_shuffled_time"]) >= {"discreteness", "fate_silhouette", "discreteness_excess"}
    pop = r["sisters"]["pop1_anterior_over_posterior"]
    assert pop["anterior_higher"] > 0.8 * pop["pairs"]  # Lin et al. 1998: POP-1 high in the anterior sister
    assert r["lateral_inhibition_model"]["deterministic_equal_start_diverges"] is False
