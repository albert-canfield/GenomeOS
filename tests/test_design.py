"""BioForge over organisms: design blocks find the perturbation that reaches a goal under constraints."""

from genomeos.lang import parse, parse_file
from genomeos.organism.forge import Knob, distance, run_design, run_designs


def test_distance_is_zero_when_the_assert_holds_and_grows_with_the_gap():
    assert distance({"ok": True, "value": 5}, "count at 1 min = 5") == 0.0
    assert distance({"ok": False, "value": 8}, "count at 1 min = 10") == 0.2
    assert distance({"ok": False, "value": 8}, "count at 1 min >= 10") == 0.2
    assert distance({"ok": False, "value": 12}, "count at 1 min <= 10") == 0.2
    assert distance({"ok": False, "value": 30}, "count at 1 min in 10..20") == 1.0
    assert distance({"ok": False, "error": "cannot parse"}, "nonsense") == 10.0


def test_celegans_designs_find_pop1_and_pie1():
    m = parse_file("data/organisms/celegans/designs.bio")
    r1, r2 = run_designs(m)
    assert r1.design.name == "two_intestinal_founders" and r1.solved and r1.best.knockouts == ["POP-1"]
    assert r1.evaluations == 7 and r1.candidates[0].label() == "-POP-1"
    assert r2.solved and r2.best.knockouts == ["PIE-1"]
    bio = r1.to_bio()
    assert "experiment design_two_intestinal_founders {" in bio and "knockout: POP-1" in bio
    assert 'evidence: predicted "BioForge search' in bio and "assert: type EPrecursor at 100 min = 2" in bio
    d = r1.to_dict()
    assert d["solved"] and d["candidates"][0]["knockouts"] == ["POP-1"] and "FAIL" in r1.format()


def test_continuous_knob_tunes_a_timer_to_a_count():
    src = """
module toy.design
import bio.std.development
organism T { root: R; cell_type: Zygote; resolution: populations }
timer cycle { duration: 60 min; when: cell_type = Zygote }
decision grow { action: divide; when: cell_type = Zygote; fraction: 1.0 }
design pace { vary: timer cycle duration 10..120; target: count at 300 min in 30..34; until: 300 min }
"""
    m = parse(src)
    r = run_design(m, m.designs[0], seed=1, iterations=25, restarts=2)
    assert r.best is not None and r.best.feasible and r.solved
    duration = r.best.values["timer:cycle.duration"]
    assert 55 <= duration <= 75  # five doublings in 300 minutes: one every 60 minutes, give or take
    assert Knob.parse("timer cycle duration 10..120").path == "timer:cycle.duration"
    assert "knobs: timer:cycle.duration" in r.to_bio()


def test_infeasible_design_reports_no_candidate():
    src = """
module toy.none
import bio.std.development
organism T { root: R; cell_type: Zygote; factors: F }
decision stay { action: quiesce; when: cell = R }
design impossible { knockout_any_of: F; target: count at 10 min = 5; keep: count at 10 min = 1 }
"""
    m = parse(src)
    r = run_design(m, m.designs[0])
    assert not r.solved and r.best is not None and r.best.feasible and r.best.loss == 4 / 5
    assert r.to_experiment() is not None  # feasible but not solved: the least-bad experiment is still emitted
    m2 = parse(src.replace("keep: count at 10 min = 1", "keep: count at 10 min = 3"))
    r2 = run_design(m2, m2.designs[0])
    assert r2.to_experiment() is None and "no feasible candidate" in r2.to_bio()
