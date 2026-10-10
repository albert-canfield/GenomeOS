"""The gastric module: a closed list sourced FOR THE STOMACH, one fitted number, 5 mutants.

The tests are the ones scripts/stomach_registration.py (committed before any number was fitted, f345ffe)
can be held to:

1. the module gates only on the registered factors, at the registered steps, and the programs on disk
   are exactly what the generator writes;
2. the wild-type run lands on its designed pools, and the ONE fitted quantity - the isthmus residence
   - falls inside an upper bound FETCHED BEFORE THE FIT (registered reading R2);
3. every registered knockout removes exactly the lineages its published STOMACH phenotype removes and
   keeps the ones it keeps, with the CONTESTED gate counted apart: holding its pattern is agreement
   with one experiment of one paper, not a reproduction of the stomach;
4. the registered readings about the formalism, including R3's numbers (the standing composition and
   the output composition diverge by an order of magnitude) and R5's (nothing rests on the one
   residence that was not fetched).

Two things are pinned here as properties of the formalism rather than as checks: the measured standing
composition, which comes back BY CONSTRUCTION because it is what the module was built from, and the
published daily output, which was used as given. Neither may be reported as agreement. And one thing
the registration did NOT foresee is pinned as such - the root compartment backs up whenever any of its
outflows is gated off, so a knockout of a small lineage raises every other lineage - in the way
tests/test_intestine.py pins what that lane's registration got wrong.
"""

from __future__ import annotations

import sys
from pathlib import Path

from genomeos.ir import to_minutes
from genomeos.lang import parse_file
from genomeos.organism.experiment import run_experiment
from genomeos.organism.stomach import (
    COMPARTMENTS,
    DAILY_OUTPUT,
    MATURE,
    ROOT_CAP,
    TOTAL_EPITHELIUM,
    mutants_bio,
    plan,
    to_bio,
)
from genomeos.organism.stomach import MUTANTS as GEN_MUTANTS
from genomeos.runtime.body import Body

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

import stomach_registration as reg  # noqa: E402

ONE_YEAR = to_minutes(1, "yr")
MODULE = "data/organisms/human/stomach.bio"
MUTANTS = "data/organisms/human/stomach_mutants.bio"


def test_the_gates_are_the_registered_closed_list_and_nothing_else():
    used = {f for c in COMPARTMENTS for f in c["requires"]} | {
        f for d in MATURE.values() for f in d["requires"]
    }
    reg.check_factors(used)  # raises on a factor the registration does not name
    assert used == reg.factors()
    assert len(reg.factors()) == 5
    for c in COMPARTMENTS:
        for f in c["requires"]:
            reg.check_gate(f, c["id"])
    for m, d in MATURE.items():
        for f in d["requires"]:
            reg.check_gate(f, m)
    assert set(MATURE) == set(reg.MATURE)
    # ATOH1 gates the secretory branch of BOTH gut modules and is excluded here for want of any stomach
    # sentence; SOX2 and GATA4 are candidates the task named and each is excluded on what was fetched
    # about it; XBP1, HNF4A, SPDEF, FOXQ1 and KLF4 are the stomach's own pattern, a factor whose loss
    # spoils a cell that is still there.
    excluded = {e["factor"] for e in reg.EXCLUDED}
    assert {"ATOH1", "SOX2", "GATA4", "XBP1", "HNF4A", "SPDEF", "FOXQ1", "KLF4"} <= excluded
    assert not (excluded & reg.factors())
    assert all(g["source"] and g["quote"] for g in reg.GATES)
    assert all(e["source"] and e["reason"] for e in reg.EXCLUDED)
    assert to_bio() == Path(MODULE).read_text()
    assert mutants_bio() == Path(MUTANTS).read_text()


def test_the_pit_lineage_has_no_gate_and_carries_most_of_the_turnover():
    """Registered reading R7. The lineage with the largest share of the stomach's daily output is the
    one no knockout but the root's touches, because every fetched pit-cell factor leaves the cell
    present and empty or hyperplastic. That is a statement about the record, not about the stomach."""
    assert MATURE["GastricPitCell"]["requires"] == []
    p = plan()
    assert p["_fit"]["output_share"]["GastricPitCell"] > 0.9
    takes_the_pit = [name for name, _f, _s, _e, lost, _k in GEN_MUTANTS if "GastricPitCell" in lost]
    assert takes_the_pit == ["rbpj"]  # the root knockout, which takes every lineage with it


def test_the_plan_solves_from_the_measured_composition_and_fits_exactly_one_quantity():
    p = plan()
    # the output shares are DERIVED from the measured standing counts and residences, not assumed
    for m in MATURE:
        assert abs(p[m]["pool"] / p[m]["residence"] - p[m]["output"]) / p[m]["output"] < 1e-12
    assert abs(sum(p[m]["output"] for m in MATURE) - DAILY_OUTPUT) / DAILY_OUTPUT < 1e-12
    # the neck compartment carries the chief route's whole throughput and has exactly one outflow
    assert abs(p["GastricNeckCell"]["output"] - p["GastricChiefCell"]["output"]) < 1e-6
    assert p["GastricNeckCell"]["growth_fraction"] == 0.0
    # the root's outflow is the whole daily output, and its residence is the one fitted number
    assert abs(p["GastricIsthmusStem"]["output"] - DAILY_OUTPUT) / DAILY_OUTPUT < 1e-12
    assert abs(p["_fit"]["total_designed"] - TOTAL_EPITHELIUM) / TOTAL_EPITHELIUM < 1e-9
    assert [c["id"] for c in COMPARTMENTS if c["residence"] is None] == ["GastricIsthmusStem"]
    assert len(reg.FITTED) == 1
    assert p["_fit"]["implied_units"] > 0  # reported, checked against nothing


def test_the_one_fitted_number_falls_inside_a_bound_fetched_before_the_fit():
    """Registered reading R2. The knob had one degree of freedom and a bound it could have missed:
    Karam & Leblond part I put the isthmal granule-free cell's turnover at less than 2.6 days, and
    that sentence was fetched and committed before the fit was run. The 2.50 d the registration states
    is arithmetic on the registered inputs, not insight; the comparison with the bound is the content."""
    p = plan()
    fitted = p["_fit"]["isthmus_residence_days"]
    r2 = next(r for r in reg.READINGS if r["id"] == "R2")
    assert abs(fitted - r2["number"]) < r2["tolerance"]
    assert fitted < r2["bound"] == p["_fit"]["granule_free_turnover_bound_days"]


def test_the_standing_composition_returns_by_construction_and_is_not_a_check():
    """The measured per-unit counts ARE the module's composition, so they come back exactly. Registered
    reading R3 is what is not by construction: with four residences a factor of 60 apart, the standing
    composition and the output composition diverge by an order of magnitude, which the colonic module
    could not show because all three of its mature residences were one fitted number."""
    p = plan()
    fit = p["_fit"]
    pit, chief = MATURE["GastricPitCell"], MATURE["GastricChiefCell"]
    ratio = (p["GastricPitCell"]["pool"] / pit["unit_count"]) / (
        p["GastricChiefCell"]["pool"] / chief["unit_count"]
    )
    assert abs(ratio - 1.0) < 1e-9  # by construction, not a check
    r3 = next(r for r in reg.READINGS if r["id"] == "R3")
    assert abs(fit["standing_share"]["GastricPitCell"] - r3["pit_standing"]) < r3["tolerance"]
    assert abs(fit["output_share"]["GastricPitCell"] - r3["pit_output"]) < r3["tolerance"]
    assert abs(fit["standing_share"]["GastricChiefCell"] - r3["chief_standing"]) < r3["tolerance"]
    assert abs(fit["output_share"]["GastricChiefCell"] - r3["chief_output"]) < r3["tolerance"]
    # the divergence itself: the standing share says almost nothing about the share of the turnover
    assert fit["output_share"]["GastricPitCell"] / fit["standing_share"]["GastricPitCell"] > 4.0
    assert fit["standing_share"]["GastricChiefCell"] / fit["output_share"]["GastricChiefCell"] > 10.0
    assert any("by construction" in s.lower() for s in reg.NOT_CLAIMED)


def test_nothing_in_the_module_rests_on_the_one_residence_that_was_not_fetched():
    """Registered reading R5. No number was fetched in this lane for the chief cell's residence, so the
    registration stated in advance that the fit barely moves between 50 and 1000 days: a longer chief
    residence raises the chief pool and lowers the chief output and the neck pool in one proportion."""
    r5 = next(r for r in reg.READINGS if r["id"] == "R5")
    low = plan(50.0)["_fit"]["isthmus_residence_days"]
    high = plan(1000.0)["_fit"]["isthmus_residence_days"]
    assert abs(low - r5["fitted_at_50_days"]) < r5["tolerance"]
    assert abs(high - r5["fitted_at_1000_days"]) < r5["tolerance"]
    assert abs(low - high) / low < 0.01
    for r in (50.0, 1000.0):
        assert abs(plan(r)["_fit"]["total_designed"] - TOTAL_EPITHELIUM) / TOTAL_EPITHELIUM < 1e-9


def test_the_wild_type_lands_on_its_designed_pools_and_returns_the_output_it_was_given():
    m = parse_file(MODULE)
    body = Body(m).run(until=ONE_YEAR)
    alive = {c.cell_type: c.count for c in body.alive_at(body.time)}
    p = plan()
    for cid in [c["id"] for c in COMPARTMENTS] + list(MATURE):
        assert 0.9 <= alive[cid] / p[cid]["pool"] <= 1.1, (cid, alive[cid], p[cid]["pool"])
    # the total is the fit's own target, so this is the fit closing and not a prediction; the daily
    # output was used as given, so its return is a construction and not a check either
    assert 0.95 <= sum(alive.values()) / TOTAL_EPITHELIUM <= 1.05
    assert 0.95 <= body.turnover_per_day() / DAILY_OUTPUT <= 1.05
    assert all(a["ok"] for a in body.check_asserts()), body.check_asserts()
    assert body.unknown == {}


def test_every_registered_knockout_removes_its_lineages_and_keeps_the_rest():
    m = parse_file(MUTANTS)
    by = {e.name: e for e in m.experiments}
    assert set(by) == {x["name"] for x in reg.MUTANTS}
    assert len(reg.MUTANTS) == 5
    for entry in reg.MUTANTS:
        r = run_experiment(m, by[entry["name"]])
        assert all(a["ok"] for a in r.asserts()), (entry["name"], r.asserts())
        alive = {c.cell_type: c.count for c in r.mutant.alive_at(r.mutant.time)}
        for lost in entry["lost"]:
            assert alive.get(lost, 0.0) == 0.0, (entry["name"], lost)
        for kept in entry["kept"]:
            assert alive.get(kept, 0.0) >= reg.KEPT_FLOOR, (entry["name"], kept)
        assert r.experiment.expect and r.experiment.evidence.source
    # neurog3 is the one single-lineage knockout with NO kept lineage registered, because Lee et al.
    # 2002 report intestinal metaplasia and so do not report the other lineages as normal
    assert [x["name"] for x in reg.MUTANTS if not x["kept"]] == ["rbpj", "neurog3"]


def test_the_contested_gate_is_counted_as_agreement_and_not_as_reproduction():
    """Registered reading R1. This registration widens "contested" from two papers to two EXPERIMENTS:
    Kim & Shivdasani 2011's pharmacological arm leaves no Ki67+ cells in most gastric glands while its
    own genetic arm leaves the units alive and 47% fewer, and one gate cannot satisfy both."""
    contested = [x for x in reg.MUTANTS if x["match"] == "contested"]
    assert {x["name"] for x in contested} == {"rbpj"}
    assert len([x for x in reg.MUTANTS if x["match"] == "partial"]) == 3
    assert len([x for x in reg.MUTANTS if x["match"] == "full"]) == 1
    m = parse_file(MUTANTS)
    by = {e.name: e for e in m.experiments}
    for x in contested:
        assert x["agrees_with"] and x["disagrees_with"] != x["agrees_with"]
        assert "CONTESTED" in by[x["name"]].expect
    for g in reg.GATES:
        if g["match"] == "contested":
            assert g["contradicted_by"] and g["contradicting_quote"]
        if g["match"] == "partial":
            assert g["not_expressible"]


def test_blocking_the_chief_route_has_nowhere_to_send_the_flow():
    """Registered reading R4, and the difference from the colon. The colonic secretory progenitor had a
    second outflow and settled at a predictable 51 times its pool; the gastric neck compartment has
    exactly one, so nothing leaves it and it accumulates for as long as the run lasts. Only the
    DIRECTION was registered, and the magnitude is a property of the run length, not of the stomach."""
    p = plan()
    m = parse_file(MUTANTS)
    by = {e.name: e for e in m.experiments}
    r = run_experiment(m, by["mist1"])
    alive = {c.cell_type: c.count for c in r.mutant.alive_at(r.mutant.time)}
    assert alive.get("GastricChiefCell", 0.0) == 0.0
    assert alive["GastricNeckCell"] > 10.0 * p["GastricNeckCell"]["pool"]  # direction only
    # the accumulation is at the wild-type chief inflow rate, which is what makes it a run-length
    # property: the pool is of the order of that inflow times the time spent in the Steady stage
    assert alive["GastricNeckCell"] < 365.0 * p["GastricChiefCell"]["output"]
    # and because the neck route is upstream of nothing else, the other three lineages do not move
    for other in ("GastricPitCell", "GastricParietalCell", "GastricEnteroendocrineCell"):
        assert abs(alive[other] / p[other]["pool"] - 1.0) < 0.01


def test_what_the_registration_did_not_foresee_the_root_backs_up_behind_any_blocked_route():
    """NOT registered in advance, and pinned here for the next module in the way
    tests/test_intestine.py pins what that lane's registration got wrong.

    The root has no homeostatic feedback, only the cap it carries as the colonic transit-amplifying
    compartment did. Its self-renewal term was solved for the WILD-TYPE outflow, so gating off any one
    route - however small - leaves that term too large and the root grows until the cap stops it,
    carrying every surviving lineage up with it. Deleting the parietal lineage, 3.8% of the daily
    output, puts the other three lineages at about 2.5 times their wild-type pools; deleting the
    endocrine lineage, 0.09% of the output, puts them at about 1.1 times. Neither multiple is a claim
    about the stomach: the direction is the formalism's and the magnitude is the cap's."""
    p = plan()
    m = parse_file(MUTANTS)
    by = {e.name: e for e in m.experiments}
    esrrg = run_experiment(m, by["esrrg"])
    alive = {c.cell_type: c.count for c in esrrg.mutant.alive_at(esrrg.mutant.time)}
    assert alive.get("GastricParietalCell", 0.0) == 0.0
    assert 2.0 < alive["GastricIsthmusStem"] / p["GastricIsthmusStem"]["pool"] <= ROOT_CAP
    for other in ("GastricPitCell", "GastricChiefCell", "GastricEnteroendocrineCell"):
        assert 1.5 < alive[other] / p[other]["pool"] <= ROOT_CAP
    # the smallest route of all still compounds over a year, by much less
    ascl1 = run_experiment(m, by["ascl1"])
    small = {c.cell_type: c.count for c in ascl1.mutant.alive_at(ascl1.mutant.time)}
    assert small.get("GastricEnteroendocrineCell", 0.0) == 0.0
    assert 1.0 < small["GastricPitCell"] / p["GastricPitCell"]["pool"] < 1.3
    # the cap is what keeps the root finite, and it is the one structural borrowing from the colonic
    # module: in the stomach the root IS the amplifying compartment
    assert ROOT_CAP == 3.0
