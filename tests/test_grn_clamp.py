# SPDX-License-Identifier: Apache-2.0
"""A clamped species in the gene-network runtime is held at its clamp value, not only between steps.

Found by lane-h on 2026-09-28 (8e705f7): `NetworkRuntime.run(clamp=...)` re-applied the clamp only after
a whole step, so inside the Runge-Kutta stages a clamped mRNA still took its basal and regulated rate
and its protein was translated off those stages. On the three-gene ring below, clamping `a.mRNA` to 0
left protein A at 1.87 against an unperturbed 19.2: a knockdown at a tenth, not a knockout. A signal
clamped to a non-zero level had the smaller version: a clamped protein decayed inside each stage, and
an external signal no gene of the module makes decayed inside the stages it was read in.

These tests watch the states the integrator actually evaluates, every stage of every step, and not
only the recorded trajectory, because the recorded points were already right while the stages were not.
"""

from __future__ import annotations

import math

import pytest

from genomeos.lang import parse
from genomeos.runtime.grn import NetworkRuntime

RING = """
module net.clamp.ring
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

# a signal the module reads and no gene of it makes, the way gastrulation reads its NODAL gradient
SIGNAL = """
module net.clamp.signal
gene g { max: 10; basal: 0.0; produces: G }
protein G {}
protein S {}
rule S activates g { strength: 1.0; threshold: 1.0; hill: 2 }
param mrna_half_life = 0.693 h {}
param protein_half_life = 0.139 h {}
param translation_rate = 5 {}
"""


def _watched(vm: NetworkRuntime) -> list[dict[str, float]]:
    """Record every state the integrator evaluates a derivative at (all four RK stages per step)."""
    seen: list[dict[str, float]] = []
    inner = vm._derivatives

    def spy(state: dict[str, float]) -> dict[str, float]:
        seen.append(dict(state))
        return inner(state)

    vm._derivatives = spy  # type: ignore[method-assign]
    return seen


CASES = [
    (RING, {"a.mRNA": 0.0}, {"A": 10.0}),  # a knockout: an mRNA held at zero
    (RING, {"B": 3.0}, {"A": 10.0}),  # a declared protein held at a non-zero level
    (SIGNAL, {"S": 2.0}, {}),  # an external signal no gene makes, held at a non-zero level
]


@pytest.mark.xfail(strict=True, reason="the clamp is re-applied only between whole steps (fixed next)")
@pytest.mark.parametrize(("source", "clamp", "initial"), CASES, ids=["mRNA-zero", "protein-3", "signal-2"])
def test_a_clamped_species_equals_its_clamp_at_every_stage_and_every_recorded_point(
    source: str, clamp: dict[str, float], initial: dict[str, float]
) -> None:
    vm = NetworkRuntime(parse(source), seed=0)
    seen = _watched(vm)
    traj = vm.run(hours=5.0, dt=0.02, initial=initial, record_every=5, clamp=clamp)

    assert len(seen) == 4 * 250  # RK4: four stages per step, and every one is checked
    for name, value in clamp.items():
        worst = max(seen, key=lambda st: abs(st.get(name, 0.0) - value)).get(name, 0.0)
        assert all(st.get(name, 0.0) == value for st in seen), (name, value, worst)
        if name in traj.levels:
            assert all(x == value for x in traj.levels[name])


@pytest.mark.xfail(strict=True, reason="the clamped mRNA is transcribed inside the stages (fixed next)")
def test_a_protein_whose_mrna_is_clamped_at_zero_only_decays() -> None:
    """Zero production: protein A follows its own decay path from its start and nothing else.

    RK4 on dA/dt = -delta * A multiplies A by the same factor every step, so the expected value is
    exact to rounding rather than to the integrator's truncation error.
    """
    dt, steps, a0 = 0.02, 300, 10.0
    vm = NetworkRuntime(parse(RING), seed=0)
    traj = vm.run(hours=steps * dt, dt=dt, initial={"A": a0}, record_every=1, clamp={"a.mRNA": 0.0})

    z = math.log(2) / 0.139 * dt
    factor = 1 - z + z**2 / 2 - z**3 / 6 + z**4 / 24
    for n, got in enumerate(traj.levels["A"]):
        assert got == pytest.approx(a0 * factor**n, rel=1e-9, abs=1e-300), (n, got)
    assert traj.levels["A"][-1] < 1e-8  # against 1.87 before the fix


def test_an_unclamped_run_is_unchanged_by_the_clamp_machinery() -> None:
    """No clamp, no change: the ring still oscillates around the 19.2 lane-h measured."""
    traj = NetworkRuntime(parse(RING), seed=0).run(hours=60.0, dt=0.02, initial={"A": 10.0}, record_every=10)
    xs = traj.levels["A"][len(traj.levels["A"]) // 2 :]
    assert sum(xs) / len(xs) == pytest.approx(19.2136, abs=1e-3)
    assert traj.peaks("A") >= 3


# The fix landed (2026-09-28, the commit after the census): the two strict expected-failure markers above
# are kept as the record of what was stated before the runtime changed, and lifted here so both tests
# now have to pass. A strict marker left in place would fail them as an unexpected pass.
for _test in (
    test_a_clamped_species_equals_its_clamp_at_every_stage_and_every_recorded_point,
    test_a_protein_whose_mrna_is_clamped_at_zero_only_decays,
):
    _test.pytestmark = [m for m in _test.pytestmark if m.name != "xfail"]  # type: ignore[attr-defined]


# ---- the declared dynamics against the independently calculated decay curve (2026-09-29) ----------
#
# What an mRNA clamp at zero promises is that the mRNA is zero; it does not remove the protein already
# made. With the only production path gone, dP/dt = -k P with k = ln 2 / the declared protein half-life,
# so P(t) = P0 exp(-k t) from the moment the clamp starts. The expected values below are that formula,
# computed here without the runtime; the tolerance is RK4's own leading truncation error on this
# equation, a relative n z^5 / 120 after n steps of z = k dt, doubled. Before 2026-09-28 (51c10b8) the
# runtime failed this by eight orders of magnitude at 5 h: the clamped mRNA was transcribed inside the
# stages, and the protein settled at k_tl r dt / (2k) instead of decaying (the ring's 1.87 at dt 0.02,
# 0.96 at dt 0.01), a steady leak proportional to the step, not a property of the model.

DECAY = """
module net.clamp.decay
gene g { max: 20; basal: 0.0; produces: P }
protein P {}
param mrna_half_life = 0.693 h {}
param protein_half_life = 0.139 h {}
param translation_rate = 5 {}
"""


def _decay_cases() -> list[tuple[str, str, str, float]]:
    return [(DECAY, "g.mRNA", "P", 5.0), (RING, "a.mRNA", "A", 60.0)]


@pytest.mark.parametrize("dt", [0.005, 0.01, 0.02])
@pytest.mark.parametrize(("source", "mrna", "protein", "hours"), _decay_cases(), ids=["decay", "ring"])
def test_a_zero_mrna_clamp_leaves_its_protein_on_the_declared_exponential_decay(
    source: str, mrna: str, protein: str, hours: float, dt: float
) -> None:
    k = math.log(2) / 0.139  # the declared protein_half_life, the only decay either module states
    # t0: the module's own state after 20 h unclamped (144 protein half-lives), so P0 is what the gene
    # made and the mRNA is non-zero when the clamp starts
    before = NetworkRuntime(parse(source), seed=0).run(hours=20.0, dt=dt, record_every=10**9).final()
    p0 = before[protein]
    assert p0 > 1.0 and before[mrna] > 1.0

    vm = NetworkRuntime(parse(source), seed=0)
    seen = _watched(vm)
    traj = vm.run(hours=hours, dt=dt, initial=dict(before), record_every=1, clamp={mrna: 0.0})

    # production is k_tl * mRNA at each stage the integrator evaluates, and the mRNA is zero at all of them
    assert all(st[mrna] == 0.0 for st in seen)
    z = k * dt
    for n, (t, got) in enumerate(zip(traj.times, traj.levels[protein], strict=True)):
        expected = p0 * math.exp(-k * t)
        assert abs(got / expected - 1) <= 2 * n * z**5 / 120, (t, got, expected)
    # not removed at once: one step after t0 the protein is still P0 exp(-k dt), and after one declared
    # half-life it is half of P0
    assert traj.levels[protein][1] > 0.9 * p0
    half = round(0.139 / dt)
    assert traj.levels[protein][half] == pytest.approx(p0 * math.exp(-k * half * dt), rel=1e-4)
