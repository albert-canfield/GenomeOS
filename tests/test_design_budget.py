"""Area H: a `design` run under a budget never answers over it, and says what the limit excluded."""

import pytest

import genomeos.organism.forge as forge
from genomeos.lang import parse, parse_file
from genomeos.organism.forge import (
    BUDGET_KNOWN_CASE,
    Budget,
    BudgetedDesignResult,
    run_design,
    run_designs,
)

CELEGANS = "data/organisms/celegans/designs.bio"
BLOOD = "data/organisms/human/haematopoiesis_designs.bio"


def _check_budget_holds(r, limit):
    rep = r.budget_report()
    over = [c for c in r.candidates if c.spend > limit + 1e-12]
    assert r.best is None or r.best.spend <= limit + 1e-12
    assert all(not c.within_budget for c in over)
    assert sorted(e["perturbation"] for e in rep["excluded"]) == sorted(c.label() for c in over)
    b = r.best
    better = {
        c.label() for c in over if c.feasible and (b is None or not b.feasible or c.loss < b.loss - 1e-12)
    }
    assert set(rep["excluded_better"]) == better and rep["limit_cost_the_target"] == bool(better)


def test_known_case_pop1_at_limit_one_and_wild_type_at_limit_zero():
    assert "Lin et al. 1995" in BUDGET_KNOWN_CASE["source"]
    m = parse_file(CELEGANS)
    one = run_designs(m, ["two_intestinal_founders"], budget=Budget(1))[0]
    assert one.solved and one.best.label() == "-POP-1" and one.best.spend == 1
    assert one.budget_report()["excluded"] == []
    zero = run_designs(m, ["two_intestinal_founders"], budget=Budget(0))[0]
    assert not zero.solved and zero.best.label() == "wild type" and zero.best.spend == 0
    rep = zero.budget_report()
    pop1 = next(e for e in rep["excluded"] if e["perturbation"] == "-POP-1")
    assert pop1["loss"] == 0 and pop1["spend"] == 1
    assert "-POP-1" in rep["excluded_better"] and rep["limit_cost_the_target"]
    assert "knockout: POP-1" not in zero.to_bio()
    assert "budget  0 perturbations: bought wild type" in zero.format()


@pytest.mark.parametrize("limit", [0, 1, 2, 3])
def test_no_celegans_design_answers_over_any_limit(limit):
    m = parse_file(CELEGANS)
    for r in run_designs(m, budget=Budget(limit)):
        _check_budget_holds(r, limit)


def test_no_blood_design_answers_over_any_limit(monkeypatch):
    """The same check on the haematopoiesis designs, whose organism runs take about a second each.

    Their designs move no knobs, so every run is of the parsed program itself and depends only on the
    knockouts and additions; one run per perturbation is shared across the four limits. Runs of any
    other module (a knob search's copy) are never cached."""
    m = parse_file(BLOOD)
    assert all(not d.vary for d in m.designs)
    real = forge.Body
    runs: dict[tuple, object] = {}

    class SharedBody:
        def __init__(self, module, seed=None, means=True, knockouts=frozenset(), adds=frozenset()):
            self.args = (module, seed, means, frozenset(knockouts), frozenset(adds))

        def run(self, until):
            module, *key = self.args
            if module is not m:
                return real(module, seed=key[0], means=key[1], knockouts=key[2], adds=key[3]).run(until=until)
            k = (*key, until)
            if k not in runs:
                runs[k] = real(m, seed=key[0], means=key[1], knockouts=key[2], adds=key[3]).run(until=until)
            return runs[k]

    monkeypatch.setattr(forge, "Body", SharedBody)
    for limit in (0, 1, 2, 3):
        for d in m.designs:
            _check_budget_holds(run_design(m, d, seed=0, budget=Budget(limit)), limit)


def test_no_budget_changes_nothing():
    m = parse_file(CELEGANS)
    plain = run_designs(m)
    assert all("budget" not in r.to_dict() and not isinstance(r, BudgetedDesignResult) for r in plain)
    wide = run_designs(m, budget=Budget(10))
    for a, b in zip(plain, wide, strict=True):
        assert [c.label() for c in a.candidates] == [c.label() for c in b.candidates]
        assert a.to_bio() == b.to_bio()


def test_priced_perturbation_moves_the_answer():
    m = parse_file(CELEGANS)
    r = run_designs(m, ["two_intestinal_founders"], budget=Budget(1, costs={"POP-1": 2}))[0]
    assert r.best.label() != "-POP-1" and not r.solved
    assert r.budget_report()["excluded"] == [
        {"perturbation": "-POP-1", "spend": 2.0, "loss": 0.0, "feasible": True}
    ]


TWO_KNOBS = """
module toy.budget
import bio.std.development
organism T { root: R; cell_type: Zygote; resolution: populations }
timer cycle { duration: 60 min; when: cell_type = Zygote }
decision grow { action: divide; when: cell_type = Zygote; fraction: 1.0 }
design pace {
  vary: timer cycle duration 10..120
  vary: decision grow fraction 0.1..1.0
  target: count at 300 min in 60..70
  until: 300 min
}
"""


@pytest.mark.parametrize("limit", [0, 1])
def test_knob_moves_are_counted_against_the_budget(limit):
    m = parse(TWO_KNOBS)
    r = run_design(m, m.designs[0], seed=1, iterations=25, restarts=2, budget=Budget(limit))
    assert r.best is not None and len(r.best.moved) <= limit
    for c in r.candidates:
        assert c.spend == len(c.knockouts) + len(c.adds) + len(c.moved)
    if limit == 0:
        assert r.best.values == {"timer:cycle.duration": 60.0, "decision:grow.fraction": 1.0}


def test_negative_limit_is_refused():
    with pytest.raises(ValueError):
        Budget(-1)
