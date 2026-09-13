"""Haematopoiesis: a mechanism module at population resolution, decided by named factors and fitted to
the measured outputs; knockouts remove the lineages the genetics says they should."""

from pathlib import Path

from genomeos.ir import to_minutes
from genomeos.lang import parse_file
from genomeos.organism.experiment import run_experiment
from genomeos.organism.haematopoiesis import COMPARTMENTS, MATURE, mutants_bio, plan, to_bio
from genomeos.runtime.body import Body

THREE_YEARS = to_minutes(3, "yr")


def test_plan_solves_the_flows_from_the_outputs_upward():
    p = plan()
    for c in COMPARTMENTS:
        if c.get("feeds"):
            assert p[c["id"]]["output"] == MATURE[c["feeds"]]["output"]
    assert p["CMP"]["output"] == p["MEP"]["inflow"] + p["GMP"]["inflow"]
    assert abs(p["MEP"]["share"] + p["GMP"]["share"] - 1.0) < 1e-9
    assert 100 < p["HSC"]["implied_division_days"] < 600  # Catlin 2011: ~280 d; Lee-Six 2018: 2-20 months
    for cid, c in p.items():
        if cid != "HSC":
            assert (
                0 < c["growth_fraction"] < c["step_out"] * 1.1
            )  # growth about balances the outflow per step
    bio = to_bio()
    assert "organism Haematopoiesis" in bio and "GATA1 = present" in bio and "decision HSC_to_MPP" in bio
    assert bio == Path("data/organisms/human/haematopoiesis.bio").read_text()
    assert mutants_bio() == Path("data/organisms/human/haematopoiesis_mutants.bio").read_text()


def test_steady_state_reproduces_pools_and_turnover():
    m = parse_file("data/organisms/human/haematopoiesis.bio")
    body = Body(m).run(until=THREE_YEARS)
    alive = {c.cell_type: c.count for c in body.alive_at(body.time)}
    p = plan()
    for cid, c in p.items():
        assert 0.8 <= alive[cid] / c["pool"] <= 1.2, (cid, alive[cid], c["pool"])
    for mid, mm in MATURE.items():
        assert 0.8 <= alive[mid] / (mm["output"] * mm["lifespan"]) <= 1.2, (mid, alive[mid])
    expected = sum(mm["output"] for mm in MATURE.values())
    assert 0.9 <= body.turnover_per_day() / expected <= 1.1
    assert all(a["ok"] for a in body.check_asserts()), body.check_asserts()
    assert body.unknown == {} and body.summary()["populations"] == 16


def test_knockouts_remove_the_dependent_lineages():
    m = parse_file("data/organisms/human/haematopoiesis_mutants.bio")
    by = {e.name: e for e in m.experiments}
    assert set(by) == {"gata1", "klf1", "spi1", "cebpa", "irf8", "pax5", "notch1", "ikzf1", "tal1"}
    for name in ("gata1", "spi1", "tal1"):
        r = run_experiment(m, by[name])
        assert all(a["ok"] for a in r.asserts()), (name, r.asserts())
        assert r.experiment.expect and r.experiment.evidence.source
    gata1 = run_experiment(m, by["gata1"])
    alive = {c.cell_type: c.count for c in gata1.mutant.alive_at(gata1.mutant.time)}
    assert alive.get("Erythrocyte", 0) == 0 and alive["Neutrophil"] > 1e9
    assert "MEP" not in alive and "CMP" in alive  # GATA1 is required at the CMP -> MEP step
    tal1 = run_experiment(m, by["tal1"])
    assert tal1.mutant.count_at(tal1.mutant.time) < 10  # no stem-cell pool, no blood
