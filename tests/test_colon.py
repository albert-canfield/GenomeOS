"""The colonic module: gated by a closed list sourced FOR THE COLON, fitted in one place, 7 mutants.

The tests are the ones scripts/colon_registration.py (committed before any number was fitted) can be
held to:

1. the module gates only on the registered factors, at the registered steps, and the programs on
   disk are exactly what the generator writes;
2. the wild-type run lands on its designed pools and on the published daily output - the closure of
   the one fitted quantity, which is not a prediction;
3. every registered knockout removes exactly the lineages its published COLON phenotype removes and
   keeps the ones it keeps, with the two CONTESTED gates counted apart: holding their pattern is
   agreement with one named paper, not a reproduction of the colon;
4. the registered readings about the formalism, including R2's number (51.0), which was computed
   from the registered shares before the run.

The measured 20% goblet share is pinned here as a CONSTRUCTION and not as a check, because the share
was inferred as the output share and all three mature residences are the one fitted number.
"""

from __future__ import annotations

import sys
from pathlib import Path

from genomeos.ir import to_minutes
from genomeos.lang import parse_file
from genomeos.organism.colon import (
    COMPARTMENTS,
    DAILY_OUTPUT,
    MATURE,
    TOTAL_EPITHELIUM,
    mutants_bio,
    plan,
    to_bio,
)
from genomeos.organism.experiment import run_experiment
from genomeos.runtime.body import Body

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

import colon_registration as reg  # noqa: E402

ONE_YEAR = to_minutes(1, "yr")
MODULE = "data/organisms/human/colon.bio"
MUTANTS = "data/organisms/human/colon_mutants.bio"


def test_the_gates_are_the_registered_closed_list_and_nothing_else():
    used = {f for c in COMPARTMENTS for f in c["requires"]} | {
        f for d in MATURE.values() for f in d["requires"]
    }
    reg.check_factors(used)  # raises on a factor the registration does not name
    assert used == reg.factors()
    assert len(reg.factors()) == 7
    for c in COMPARTMENTS:
        for f in c["requires"]:
            reg.check_gate(f, c["id"])
    for m, d in MATURE.items():
        for f in d["requires"]:
            reg.check_gate(f, m)
    assert set(MATURE) == set(reg.MATURE)
    # ASCL2 and SPDEF gate the small-intestinal module and are excluded here, for want of a colon
    # phenotype and not for want of genetics: the point of a per-organ closed list.
    excluded = {e["factor"] for e in reg.EXCLUDED}
    assert {"ASCL2", "SPDEF", "REG4", "POU2F3", "HES1", "NOX1"} <= excluded
    assert not (excluded & reg.factors())
    assert all(g["source"] and g["quote"] for g in reg.GATES)
    assert to_bio() == Path(MODULE).read_text()
    assert mutants_bio() == Path(MUTANTS).read_text()


def test_the_plan_solves_the_flows_and_fits_exactly_one_quantity():
    p = plan()
    assert abs(p["ColonTA"]["output"] - DAILY_OUTPUT) / DAILY_OUTPUT < 1e-12
    assert abs(p["ColonTA"]["inflow"] * p["ColonTA"]["amplification"] - DAILY_OUTPUT) / DAILY_OUTPUT < 1e-12
    assert p["ColonStem"]["output"] == p["ColonTA"]["inflow"]
    assert abs(sum(p[m]["output"] for m in MATURE) - DAILY_OUTPUT) / DAILY_OUTPUT < 1e-12
    # the one fitted quantity closes the published count, and every mature lineage shares it: the colon
    # has no resident type, so there is no second residence to fit
    assert abs(p["_fit"]["total_designed"] - TOTAL_EPITHELIUM) / TOTAL_EPITHELIUM < 1e-9
    assert [m for m, d in MATURE.items() if d["residence"] is None] == list(MATURE)
    surface = p["_fit"]["surface_residence_days"]
    assert 0.0 < surface < p["_fit"]["implied_mean_residence_days"]
    for m in MATURE:
        assert p[m]["residence"] == surface
    assert p["ColonSecretoryProgenitor"]["growth_fraction"] == 0.0
    assert 0.0 < p["ColonTA"]["growth_fraction"] < 1.0
    assert p["ColonStem"]["implied_crypts"] > 0  # reported, checked against nothing


def test_the_measured_goblet_share_returns_by_construction_and_is_not_a_check():
    """Katz et al. 2002 measured 20% of colonic epithelial cells as mature goblet cells, and that
    share was used as the share of daily output. With all three mature residences equal to the one
    fitted number, it therefore comes back exactly - as a share of the MATURE cells. As a share of
    every epithelial cell the crypt compartments make it smaller, and Katz counted total epithelial
    cells in sections, so neither figure may be reported as agreement."""
    p = plan()
    assert abs(p["_fit"]["goblet_share_of_mature"] - reg.MEASURED[3]["value"]) < 1e-12
    assert p["_fit"]["goblet_share_of_all"] < 0.16  # the same pool over the whole epithelium
    assert any("by construction" in s.lower() for s in reg.NOT_CLAIMED)


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
    assert len(reg.MUTANTS) == 7
    for entry in reg.MUTANTS:
        r = run_experiment(m, by[entry["name"]])
        assert all(a["ok"] for a in r.asserts()), (entry["name"], r.asserts())
        alive = {c.cell_type: c.count for c in r.mutant.alive_at(r.mutant.time)}
        for lost in entry["lost"]:
            assert alive.get(lost, 0.0) == 0.0, (entry["name"], lost)
        for kept in entry["kept"]:
            assert alive.get(kept, 0.0) >= reg.KEPT_FLOOR, (entry["name"], kept)
        assert r.experiment.expect and r.experiment.evidence.source


def test_the_two_contested_gates_are_counted_as_agreement_and_not_as_reproduction():
    """Registered reading R1. TCF7L2 and KLF4 have two fetched colon papers that disagree, so the
    module's gate can agree with one of them; the experiment text must name the other."""
    contested = [x for x in reg.MUTANTS if x["match"] == "contested"]
    assert {x["name"] for x in contested} == {"tcf7l2", "klf4"}
    assert len([x for x in reg.MUTANTS if x["match"] == "partial"]) == 3
    assert len([x for x in reg.MUTANTS if x["match"] == "full"]) == 2
    m = parse_file(MUTANTS)
    by = {e.name: e for e in m.experiments}
    for x in contested:
        assert x["agrees_with"] and x["disagrees_with"] != x["agrees_with"]
        text = by[x["name"]].expect
        assert "CONTESTED" in text and x["disagrees_with"].split(" et al")[0] in text
    for g in reg.GATES:
        if g["match"] == "contested":
            assert g["contradicted_by"] and g["contradicting_quote"]


def test_the_registered_readings_about_the_formalism():
    """R2 was a number written before the run: with the goblet route blocked the secretory progenitor
    backs up until its one open outflow carries the whole inflow, so the enteroendocrine pool stands
    at 1 / (0.004 / 0.204) = 51.0 times the wild type. R3 registered only a direction, because the
    transit-amplifying compartment carries a growth term and a cap. Both are properties of a module
    made of flows, not claims about the colon."""
    p = plan()
    predicted = next(r["number"] for r in reg.READINGS if r["id"] == "R2")
    assert abs(p["_fit"]["eec_fold_if_goblet_blocked"] - predicted) < 1e-9
    m = parse_file(MUTANTS)
    by = {e.name: e for e in m.experiments}
    for name in ("gfi1", "sox9", "klf4"):
        r = run_experiment(m, by[name])
        alive = {c.cell_type: c.count for c in r.mutant.alive_at(r.mutant.time)}
        fold = alive["ColonEnteroendocrineCell"] / p["ColonEnteroendocrineCell"]["pool"]
        assert abs(fold - predicted) / predicted < 0.01, (name, fold)
        assert alive["ColonSecretoryProgenitor"] / p["ColonSecretoryProgenitor"]["pool"] > 0.99 * predicted
    # R3: the direction only. Lo et al. 2016 report the absorptive lineage remaining, not expanding.
    atoh1 = run_experiment(m, by["atoh1"])
    kept = {c.cell_type: c.count for c in atoh1.mutant.alive_at(atoh1.mutant.time)}
    assert kept["Colonocyte"] > 1.5 * p["Colonocyte"]["pool"]
