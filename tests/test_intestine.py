"""The small-intestinal module: gated by a closed sourced list, fitted in one place, mutants registered.

Three things are tested, and they are the three the registration (scripts/intestine_registration.py,
committed before any number was fitted) can be held to:

1. the module gates only on the registered factors, at the registered steps, and the generated
   programs on disk are exactly what the generator writes;
2. the wild-type run lands on its designed pools and on the published daily output — the closure of
   the one fitted quantity, which is not a prediction and is labelled so;
3. every registered knockout removes exactly the lineages the published phenotype removes and keeps
   the ones it keeps, including the accumulation of secretory progenitors registered for SPDEF.

The last test also pins the two consequences of the formalism that the registration did not predict:
a blocked branch inflates the open ones, because the growth term was solved for the wild-type
outflows. They are facts about a module made of flows, not claims about the intestine.
"""

from __future__ import annotations

import sys
from pathlib import Path

from genomeos.ir import to_minutes
from genomeos.lang import parse_file
from genomeos.organism.experiment import run_experiment
from genomeos.organism.intestine import (
    COMPARTMENTS,
    DAILY_OUTPUT,
    MATURE,
    TOTAL_EPITHELIUM,
    mutants_bio,
    plan,
    to_bio,
)
from genomeos.runtime.body import Body

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

import intestine_registration as reg  # noqa: E402

ONE_YEAR = to_minutes(1, "yr")
MODULE = "data/organisms/human/intestine.bio"
MUTANTS = "data/organisms/human/intestine_mutants.bio"


def test_the_gates_are_the_registered_closed_list_and_nothing_else():
    used = {f for c in COMPARTMENTS for f in c["requires"]} | {
        f for d in MATURE.values() for f in d["requires"]
    }
    reg.check_factors(used)  # raises on a factor the registration does not name
    assert used == reg.factors()
    for c in COMPARTMENTS:
        for f in c["requires"]:
            reg.check_gate(f, c["id"])
    for m, d in MATURE.items():
        for f in d["requires"]:
            reg.check_gate(f, m)
    assert set(MATURE) == set(reg.MATURE)
    assert to_bio() == Path(MODULE).read_text()
    assert mutants_bio() == Path(MUTANTS).read_text()


def test_the_plan_solves_the_flows_and_fits_exactly_one_quantity():
    p = plan()
    # the output of the amplifying compartment is the whole measured output, and what it needs from
    # the stem compartment is that divided by the amplification
    assert abs(p["TA"]["output"] - DAILY_OUTPUT) / DAILY_OUTPUT < 1e-12
    assert abs(p["TA"]["inflow"] * p["TA"]["amplification"] - DAILY_OUTPUT) / DAILY_OUTPUT < 1e-12
    assert p["ISC"]["output"] == p["TA"]["inflow"]
    assert abs(sum(p[m]["output"] for m in MATURE) - DAILY_OUTPUT) / DAILY_OUTPUT < 1e-12
    # the one fitted quantity closes the published count, and it is the only residence without a source
    assert abs(p["_fit"]["total_designed"] - TOTAL_EPITHELIUM) / TOTAL_EPITHELIUM < 1e-9
    fitted = [m for m, d in MATURE.items() if d["residence"] is None]
    assert fitted == ["Enterocyte", "EnteroendocrineCell"]
    villus = p["_fit"]["villus_residence_days"]
    assert 3.0 < villus < p["_fit"]["implied_mean_residence_days"]  # between the mouse 3 d and the mean
    for m in fitted:
        assert p[m]["residence"] == villus
    # growth balances outflow per step everywhere: a pass-through compartment needs none
    assert p["SecretoryProgenitor"]["growth_fraction"] == 0.0
    assert 0.0 < p["TA"]["growth_fraction"] < 1.0
    assert p["ISC"]["implied_crypts"] > 0  # reported, checked against nothing


def test_the_wild_type_lands_on_its_designed_pools_and_the_published_output():
    m = parse_file(MODULE)
    body = Body(m).run(until=ONE_YEAR)
    alive = {c.cell_type: c.count for c in body.alive_at(body.time)}
    p = plan()
    for cid in [c["id"] for c in COMPARTMENTS] + list(MATURE):
        assert 0.9 <= alive[cid] / p[cid]["pool"] <= 1.1, (cid, alive[cid], p[cid]["pool"])
    assert 0.95 <= body.turnover_per_day() / DAILY_OUTPUT <= 1.05
    assert all(a["ok"] for a in body.check_asserts()), body.check_asserts()
    assert body.unknown == {}


def test_every_registered_knockout_removes_its_lineages_and_keeps_the_rest():
    m = parse_file(MUTANTS)
    by = {e.name: e for e in m.experiments}
    assert set(by) == {x["name"] for x in reg.MUTANTS}
    wild = plan()
    for entry in reg.MUTANTS:
        r = run_experiment(m, by[entry["name"]])
        assert all(a["ok"] for a in r.asserts()), (entry["name"], r.asserts())
        alive = {c.cell_type: c.count for c in r.mutant.alive_at(r.mutant.time)}
        for lost in entry["lost"]:
            assert alive.get(lost, 0.0) == 0.0, (entry["name"], lost)
        for kept in entry["kept"]:
            assert alive.get(kept, 0.0) >= reg.KEPT_FLOOR, (entry["name"], kept)
        assert r.experiment.expect and r.experiment.evidence.source
    # SPDEF: the registered extra, an accumulation of secretory progenitors above the wild-type pool
    spdef = run_experiment(m, by["spdef"])
    alive = {c.cell_type: c.count for c in spdef.mutant.alive_at(spdef.mutant.time)}
    assert alive["SecretoryProgenitor"] > 2.0 * wild["SecretoryProgenitor"]["pool"]
    # the two consequences of the formalism, pinned so they cannot be mistaken for claims: a blocked
    # branch sends its share to the open ones, because the growth term was solved for the wild type
    atoh1 = run_experiment(m, by["atoh1"])
    kept = {c.cell_type: c.count for c in atoh1.mutant.alive_at(atoh1.mutant.time)}
    assert kept["Enterocyte"] > 1.5 * wild["Enterocyte"]["pool"]  # the paper says "unaffected"
    gfi1 = run_experiment(m, by["gfi1"])
    kept = {c.cell_type: c.count for c in gfi1.mutant.alive_at(gfi1.mutant.time)}
    assert kept["EnteroendocrineCell"] > 5.0 * wild["EnteroendocrineCell"]["pool"]  # supernumerary
