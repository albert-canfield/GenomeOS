"""The human body as counted populations: the same Body runtime at population resolution."""

from genomeos.ir import to_minutes
from genomeos.lang import parse_file
from genomeos.runtime.body import Body

PROGRAM = "data/organisms/human/body.bio"


def test_human_body_grows_to_adult_counts_and_turns_over():
    m = parse_file(PROGRAM)
    assert m.organism.resolution == "populations" and len(m.stages) == 7
    body = Body(m).run(until=to_minutes(20, "yr"))
    day = lambda d: to_minutes(d, "d")  # noqa: E731
    assert body.count_at(day(1)) == 2 and body.count_at(day(3)) == 8  # cleavage: 2-cell day 1, 8-cell day 3
    assert 30 <= body.count_at(day(5)) <= 130  # blastocyst
    assert 1e2 <= body.count_at(day(21)) <= 1e4  # gastrulation done
    assert 1e12 <= body.count_at(day(266)) <= 1e13  # birth
    total = body.count_at(body.time)
    assert 2.0e13 <= total <= 3.5e13, total  # Sender & Milo 3.0e13; Hatton 36e12 for a 70 kg male
    s = body.summary()
    assert 2.0e11 <= s["turnover_per_day"] <= 4.5e11, s["turnover_per_day"]  # 3.3e11 per day published
    alive = {c.cell_type: c.count for c in body.alive_at(body.time)}
    assert alive["Erythrocyte"] >= 2.0e13 and alive["NeuronHuman"] <= 1.5e11 and alive["Myocyte"] <= 4e8
    assert alive["Trophectoderm"] < 100 and alive["Hypoblast"] < 100  # extra-embryonic, not counted in growth
    assert s["populations"] == 18 and s["cells_born"] == 18 and body.unknown == {}
    assert all(c["ok"] for c in body.check_asserts()), body.check_asserts()
    rep = body.uncertainty().to_dict()
    assert rep["organism"]["label"] in ("low", "medium") and rep["cellular"]["label"] == "low"
    assert rep["organism"]["items"] > 10_000  # every daily decision counted
    tree = body.tree(depth=2)
    assert any("[Erythrocyte]" in line for line in tree)


def test_population_semantics_in_isolation():
    src = """
module toy.pop
import bio.std.development
organism P { root: R; resolution: populations; cell_type: Zygote; assert: count at 30 min = 8 }
stage Grow { from: 0 min; to: 31 min }
stage Hold { from: 31 min }
cell_type A { parent: Blastomere }
cell_type B { parent: Blastomere }
decision blast { action: differentiate; when: cell_type = Zygote; to: Blastomere }
decision double { action: divide; when: cell_type = Blastomere, stage = Grow; fraction: 1.0; after: 10 min }
decision split { action: differentiate; when: cell_type = Blastomere, stage = Hold; to: A; fraction: 0.25 }
decision rest { action: differentiate; when: cell_type = Blastomere, stage = Hold; to: B; fraction: 1.0 }
decision cap { action: quiesce; when: cell_type = B, count = >=12 }
decision b_grow { action: divide; when: cell_type = B, stage = Hold; fraction: 0.5; after: 5 min }
decision b_loss { action: die; when: cell_type = B, stage = Hold; fraction: 0.2; after: 5 min }
"""
    body = Body(parse_file.__globals__["parse"](src)).run(until=200)
    assert body.count_at(30) == 8 and body.check_asserts()[0]["ok"]
    a = next(c for c in body.cells.values() if c.cell_type == "A")
    b = next(c for c in body.cells.values() if c.cell_type == "B")
    assert a.count == 2 and a.population and b.population and b.name == "R"
    # B: +50% and -20% every 5 min, capped at 12 by quiescence, regrowing after each loss
    assert 6 <= b.count <= 12.5 and body.culled > 0 and body.history[-1][0] <= 200
    assert body.count_at(31) == 8  # history answers past totals for populations
