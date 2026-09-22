"""Space in the Body: positions, fields stepped between events, gradient signals, placement of daughters,
contact inhibition, migration, and the French flag grown from one cell."""

from genomeos.ir import to_minutes
from genomeos.lang import parse, parse_file
from genomeos.organism.experiment import run_experiment
from genomeos.runtime.body import Body

HOURS = 60.0


def test_flag_grows_from_one_cell_into_three_bands():
    m = parse_file("data/demo/flag_organism.bio")
    body = Body(m).run(until=to_minutes(200, "h"))
    assert body.spatial and body.count_at(body.time) == 30 and body.blocked_divisions > 0
    rows = body.type_map()
    assert len(rows) == 1 and rows[0].startswith("BBB") and rows[0].endswith("RRR") and "W" in rows[0]
    bands = body.bands_along_x()
    assert [b[0] for b in bands] == ["Blue", "White", "Red"]  # positional information in order
    profile = body.fields["Morphogen"].profile_x()
    assert profile[0] > profile[10] > profile[29]  # a decaying gradient from the source edge
    assert all(a["ok"] for a in body.check_asserts()), body.check_asserts()
    # the founder's site stays occupied by its first-daughter line; every blue cell sits nearest the source
    xs = sorted(x for _, x, _, t in body.positions() if t == "Blue")
    assert xs == list(range(7)) and (0, 0) in body.occupied


def test_daughters_take_free_sites_and_deaths_free_them():
    src = """
module toy.space
import bio.std.development
organism S { root: A; cell_type: Zygote; space: 4 x 1; origin: 0,0 }
decision first { action: divide; when: cell = A; daughters: A1, A2; direction: +x; after: 10 min }
decision second { action: divide; when: cell = A2; daughters: B1, B2; direction: +x; after: 10 min }
decision third { action: divide; when: cell = B2; daughters: C1, C2; direction: +x; after: 10 min }
decision fourth { action: divide; when: cell = C2; daughters: D1, D2; direction: +x; after: 10 min }
decision gone { action: die; when: cell = B1; after: 25 min }
"""
    body = Body(parse(src)).run(until=100)
    # A1 keeps the founder's site; each second daughter takes the next free site to the right; C2 finds
    # no room at 40 min, waits a cycle, and by 50 min B1 has died (at 45) so D2 takes its freed site
    sites = {name: (x, y) for name, x, y, _ in body.positions()}
    assert sites == {"A1": (0, 0), "D2": (1, 0), "C1": (2, 0), "D1": (3, 0)}
    assert body.cells["B1"].dies_at == 45 and body.cells["C2"].divides_at == 50
    assert body.count_at(100) == 4 and body.blocked_divisions == 1


def test_migration_toward_a_source_and_along_a_direction():
    src = """
module toy.move
import bio.std.development
organism M { root: W; cell_type: Zygote; space: 12 x 1; origin: 11,0; sense: 10 min }
field Attractant { diffusion: 0.4; decay: 0.02; source: 0,0 = 1.0 }
cell_type Near { parent: Blastomere }
signal Close { mode: gradient; field: Attractant; threshold: 1.0; sets: Close = present }
decision chase {
  action: migrate; when: cell = W, Close = absent; toward: Attractant; steps: 1; after: 60 min
}
decision arrive { action: differentiate; when: cell = W, Close = present; to: Near }
"""
    body = Body(parse(src)).run(until=to_minutes(100, "h"))
    w = body.cells["W"]
    assert w.x is not None and w.x < 11  # it walked up the gradient
    assert (
        w.cell_type == "Near" and w.factors.get("Close") == "present"
    )  # and stopped where the signal is high
    steps = """
module toy.steps
import bio.std.development
organism T { root: R; cell_type: Zygote; space: 10 x 1; origin: 9,0 }
decision left { action: migrate; when: cell = R; direction: -x; steps: 3 }
"""
    b2 = Body(parse(steps)).run(until=10)
    assert b2.cells["R"].x == 6 and b2.occupied == {(6, 0): "R"}


def test_gradient_signal_knockout_removes_the_high_fate():
    m = parse(
        parse_file.__globals__["Path"]("data/demo/flag_organism.bio").read_text()
        + '\nexperiment noHigh { knockout: High; until: 200 h; expect: "no blue band" }\n',
        base_dir=parse_file.__globals__["Path"]("data/demo"),
    )
    r = run_experiment(m, m.experiments[0])
    fates = r.mutant.summary()["fates"]
    assert "Blue" not in fates and fates.get("White", 0) >= 10 and r.wild_type.summary()["fates"]["Blue"] == 7
