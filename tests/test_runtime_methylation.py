"""The CpG methylation state machine: what must hold whatever the rates are.

Conservation, determinism under a seed, the exact answers the mechanism has when a rate is 0 or 1,
the neighbour dependence the whole prediction rests on, and the refusal to run on UNKNOWN rates.
"""

from pathlib import Path

import pytest

from genomeos.ir import UNKNOWN
from genomeos.lang import parse
from genomeos.runtime.methylation import (
    DEFAULT_CONTEXTS,
    PARAM_PREFIX,
    PARAMETERS,
    CpGContext,
    MethylationParams,
    UnknownParametersError,
    divisions_for_level,
    expected_trajectory,
    fidelity_for_level,
    initial_distribution,
    level,
    run_module,
    simulate,
    steady_state,
)

SOLO = CpGContext("solo_wcgw", 0, "W")
DENSE = CpGContext("dense", 8, "S")


def params(**over: float) -> MethylationParams:
    base = dict(
        maintenance_fidelity_dense=0.97,
        maintenance_fidelity_solo_wcgw=0.90,
        maintenance_fidelity_solo_scgs=0.95,
        de_novo_dense=0.05,
        de_novo_solo=0.01,
        tet_erasure=0.0,
        neighbour_scale=2.0,
    )
    base.update(over)
    return MethylationParams(**base)  # type: ignore[arg-type]


# ---- conservation -------------------------------------------------------------------


def test_the_three_states_always_hold_every_cpg():
    p = params(tet_erasure=0.02)
    for ctx in DEFAULT_CONTEXTS:
        dist = initial_distribution(0.8, hemi=0.1)
        for _ in range(200):
            from genomeos.runtime.methylation import division

            dist = division(dist, ctx, p)
            assert abs(sum(dist) - 1.0) < 1e-12
            assert all(x >= -1e-15 for x in dist)
            assert 0.0 <= level(dist) <= 1.0


def test_no_cpg_is_created_or_lost_in_the_stochastic_run():
    run = simulate(params(tet_erasure=0.01), divisions=40, cpgs=200, seed=3)
    for counts in run.final_states.values():
        assert sum(counts.values()) == 200


# ---- determinism --------------------------------------------------------------------


def test_the_same_seed_gives_the_same_run_and_a_different_seed_does_not():
    a = simulate(params(), divisions=30, cpgs=100, seed=7)
    b = simulate(params(), divisions=30, cpgs=100, seed=7)
    c = simulate(params(), divisions=30, cpgs=100, seed=8)
    assert a.levels == b.levels and a.final_states == b.final_states
    assert a.levels != c.levels


# ---- the exact answers the mechanism has -------------------------------------------


def test_perfect_maintenance_with_no_de_novo_and_no_erasure_keeps_the_state():
    p = params(
        maintenance_fidelity_dense=1.0,
        maintenance_fidelity_solo_wcgw=1.0,
        maintenance_fidelity_solo_scgs=1.0,
        de_novo_dense=0.0,
        de_novo_solo=0.0,
        tet_erasure=0.0,
    )
    for ctx in DEFAULT_CONTEXTS:
        assert expected_trajectory(ctx, p, 500, initial_level=1.0) == pytest.approx([1.0] * 501)
        assert expected_trajectory(ctx, p, 500, initial_level=0.0) == pytest.approx([0.0] * 501)
        assert expected_trajectory(ctx, p, 500, initial_level=0.4) == pytest.approx([0.4] * 501)
    run = simulate(p, divisions=50, cpgs=100, seed=1, initial_level=0.6)
    for xs in run.levels.values():
        assert xs == pytest.approx([xs[0]] * 51)


def test_no_maintenance_halves_the_level_at_every_division():
    p = params(
        maintenance_fidelity_dense=0.0,
        maintenance_fidelity_solo_wcgw=0.0,
        maintenance_fidelity_solo_scgs=0.0,
        de_novo_dense=0.0,
        de_novo_solo=0.0,
    )
    xs = expected_trajectory(DENSE, p, 6, initial_level=1.0)
    assert xs == pytest.approx([1.0, 0.5, 0.25, 0.125, 0.0625, 0.03125, 0.015625])


def test_erasure_alone_removes_methylation_even_with_perfect_maintenance():
    p = params(
        maintenance_fidelity_dense=1.0,
        maintenance_fidelity_solo_wcgw=1.0,
        maintenance_fidelity_solo_scgs=1.0,
        de_novo_dense=0.0,
        de_novo_solo=0.0,
        tet_erasure=0.05,
    )
    xs = expected_trajectory(DENSE, p, 200, initial_level=1.0)
    assert xs[-1] < 0.01 and all(b <= a + 1e-12 for a, b in zip(xs, xs[1:], strict=False))


# ---- the neighbour dependence, which is what the prediction rests on ----------------


def test_an_isolated_cpg_erodes_faster_than_a_clustered_one_and_the_gap_widens():
    p = params()
    solo = expected_trajectory(SOLO, p, 100, initial_level=0.95)
    dense = expected_trajectory(DENSE, p, 100, initial_level=0.95)
    gaps = [d - s for d, s in zip(dense, solo, strict=True)]
    assert gaps[0] == pytest.approx(0.0)
    assert gaps[10] < gaps[30] < gaps[60] and gaps[-1] > 0.5


def test_maintenance_and_de_novo_rise_with_the_neighbour_count():
    p = params()
    fid = [p.maintenance(CpGContext(f"n{n}", n, "W")) for n in range(0, 12)]
    den = [p.de_novo(CpGContext(f"n{n}", n, "W")) for n in range(0, 12)]
    assert fid == sorted(fid) and den == sorted(den)
    assert fid[0] == pytest.approx(p.maintenance_fidelity_solo_wcgw)
    assert fid[-1] == pytest.approx(p.maintenance_fidelity_dense, abs=1e-3)
    # only the solo end depends on the flanking bases
    assert p.maintenance(CpGContext("s", 0, "S")) > p.maintenance(CpGContext("w", 0, "W"))
    assert p.maintenance(CpGContext("s", 20, "S")) == pytest.approx(
        p.maintenance(CpGContext("w", 20, "W")), abs=1e-5
    )


def test_the_sampled_mean_tracks_the_exact_mean():
    run = simulate(params(), divisions=60, cpgs=3000, seed=11, initial_level=0.95)
    for key, xs in run.levels.items():
        if key.endswith("gap"):
            continue
        assert xs[-1] == pytest.approx(run.expected[key][-1], abs=0.03)


# ---- the steady state ---------------------------------------------------------------


def test_the_steady_state_is_a_fixed_point():
    from genomeos.runtime.methylation import division

    p = params(tet_erasure=0.01)
    for ctx in (SOLO, DENSE):
        st, n = steady_state(ctx, p)
        assert 0 < n < 200_000
        dist = initial_distribution(1.0)
        for _ in range(n + 50):
            dist = division(dist, ctx, p)
        assert level(dist) == pytest.approx(st, abs=1e-6)
        again = division(dist, ctx, p)
        assert level(again) == pytest.approx(st, abs=1e-6)


def test_the_solo_steady_state_is_below_the_dense_one():
    p = params()
    assert steady_state(SOLO, p)[0] < steady_state(DENSE, p)[0]


# ---- the inverses the application script fits with ---------------------------------


def test_fidelity_for_level_round_trips():
    p = params()
    for target_f in (0.90, 0.95, 0.99):
        got = expected_trajectory(
            SOLO, params(maintenance_fidelity_solo_wcgw=target_f), 80, initial_level=0.9
        )[-1]
        back = fidelity_for_level(got, 80, p, SOLO, knob="solo", initial_level=0.9)
        assert back == pytest.approx(target_f, abs=1e-4)


def test_an_unreachable_level_is_refused_rather_than_fitted():
    p = params(de_novo_dense=0.05, de_novo_solo=0.05)
    # no fidelity in 0..1 can hold a dense CpG below the level de novo methylation alone sustains
    assert fidelity_for_level(0.001, 200, p, DENSE, knob="dense", initial_level=0.9) is None
    assert fidelity_for_level(0.999, 200, p, DENSE, knob="dense", initial_level=0.9) is None


def test_divisions_for_level_counts_the_divisions_a_loss_needs():
    p = params()
    n = divisions_for_level(0.5, p, SOLO, initial_level=0.95)
    assert n is not None and expected_trajectory(SOLO, p, n, initial_level=0.95)[-1] <= 0.5
    assert expected_trajectory(SOLO, p, n - 1, initial_level=0.95)[-1] > 0.5
    # a level below the steady state is never reached, and the engine says so
    assert divisions_for_level(0.01, p, DENSE, initial_level=0.95) is None


# ---- refusing to invent a rate ------------------------------------------------------


def test_the_standard_library_says_unknown_where_the_literature_does_not_establish_a_rate():
    m = parse("module t\nimport bio.std.methylation\n")
    unknown = sorted(k for k, p in m.parameters.items() if p.value is UNKNOWN)
    assert unknown == [
        PARAM_PREFIX + n
        for n in (
            "de_novo_dense",
            "de_novo_solo",
            "maintenance_fidelity_solo_scgs",
            "maintenance_fidelity_solo_wcgw",
            "neighbour_scale",
            "tet_erasure",
        )
    ]
    assert m.parameters[PARAM_PREFIX + "maintenance_fidelity_dense"].value == 0.97
    assert all(p.evidence.source for p in m.parameters.values())
    with pytest.raises(UnknownParametersError) as e:
        MethylationParams.from_module(m)
    assert "tet_erasure" in str(e.value) and "neighbour_scale" in str(e.value)


def test_a_caller_may_bind_the_unknown_rates_and_the_run_records_them():
    m = parse("module t\nimport bio.std.methylation\n")
    p = MethylationParams.from_module(
        m,
        assume={
            "maintenance_fidelity_solo_wcgw": 0.9,
            "maintenance_fidelity_solo_scgs": 0.95,
            "de_novo_dense": 0.05,
            "de_novo_solo": 0.01,
            "tet_erasure": 0.0,
            "neighbour_scale": 2.0,
        },
    )
    assert p.maintenance_fidelity_dense == 0.97  # the one number the literature does establish
    assert len(p.assumptions) == 6 and any("tet_erasure" in a for a in p.assumptions)
    assert {q.name for q in p.parameters()} == set(PARAMETERS)


def test_the_demo_program_binds_them_and_runs():
    run = run_module(parse(Path("data/demo/methylation_erosion.bio").read_text()))
    assert run.divisions == 100 and run.seed == 0 and run.cpgs == 500
    solo = run.levels[PARAM_PREFIX + "solo_wcgw"]
    dense = run.levels[PARAM_PREFIX + "dense"]
    assert solo[0] == pytest.approx(0.95, abs=0.01) and solo[-1] < 0.3 < dense[-1]
    assert run.summary()["assumptions"]


def test_an_unknown_parameter_survives_the_ir_round_trip():
    from genomeos.ir import Module

    m = parse(
        'module t\nparam k = unknown /division { evidence: experimental "not measured"; confidence: 0 }\n'
        'param j = 0.5 /division { evidence: experimental "measured"; confidence: 0.7 }\n'
    )
    assert m.parameters["k"].value is UNKNOWN and m.parameters["j"].value == 0.5
    back = Module.from_dict(m.to_dict())
    assert back.parameters["k"].value is UNKNOWN and back.parameters["j"].value == 0.5


def test_a_rate_outside_zero_to_one_is_a_modelling_error():
    with pytest.raises(ValueError):
        params(maintenance_fidelity_dense=1.4)
    with pytest.raises(ValueError):
        params(neighbour_scale=0.0)
    with pytest.raises(ValueError):
        initial_distribution(1.2)
    with pytest.raises(ValueError):
        CpGContext("x", 0, "N")


# ---- the claim, as it stands in the committed result --------------------------------


def test_the_registered_claim_and_its_verdict_are_on_the_record():
    from genomeos.results import load_result

    r = load_result("methylation_erosion")
    if r is None:  # the run needs the WGBS caches, which are not in the repository
        pytest.skip("methylation_erosion has not been run on this machine")
    assert r["claim"]["name"] == "solo_wcgw_erosion_ordering"
    assert r["claim"]["registered"] <= r["date"]
    v = r["verdict"]
    assert set(v) == {"gap", "gap_over_dense", "gap_within_50kb"}
    # the ordering that was predicted: the long-cultured lines separate the two classes more
    assert v["gap"]["holds"] and v["gap_within_50kb"]["holds"]
    assert v["gap"]["min_long_cultured"] > v["gap"]["max_intact"]
    # and the direction the engine insists on: a solo CpG is never the better-kept one in a line
    for cell in r["claim"]["long_cultured"]:
        assert r["measured"][cell]["gap"] > 0
