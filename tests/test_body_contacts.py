# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""BioLang v0.4 §7.4: neighbours, contact amounts, per-cell networks, seeded noise and replicate asserts.

The falsifier is area E's reproduction of Collier et al. 1996: two equivalent cells must not diverge
without noise, must diverge in nearly every run with it, and must split the winner about evenly across
seeds. A runtime that picks the same winner every run has smuggled in an order, and the tests below
check that swapping the order the cells are created in changes nothing.
"""

from __future__ import annotations

from collections import Counter
from pathlib import Path

import pytest

from genomeos.lang import parse, parse_file
from genomeos.runtime.body import Body, evaluate_replicate_asserts, read_contacts, replicate

ROOT = Path(__file__).resolve().parent.parent
COLLIER = ROOT / "data" / "demo" / "lateral_inhibition.bio"
PAIR = ("Z1.ppp", "Z4.aaa")


def _anchor(body: Body) -> list[str]:
    return [n for n in PAIR if body.cells[n].cell_type == "AnchorCell"]


def test_two_equal_cells_never_diverge_without_noise():
    m = parse_file(COLLIER)
    m.parameters["noise"].value = 0.0
    for seed in range(5):
        b = Body(m, seed=seed).run(until=40 * 60)
        a, z = (b.cells[n].levels["Dp"] for n in PAIR)
        assert a == z and _anchor(b) == []  # the symmetric state is a steady state


def test_noise_decides_and_the_winner_is_not_an_order():
    m = parse_file(COLLIER)
    swapped = parse(COLLIER.read_text().replace("daughters: Z1.ppp, Z4.aaa", "daughters: Z4.aaa, Z1.ppp"))
    winners, same = Counter(), 0
    for seed in range(30):
        first = _anchor(Body(m, seed=seed).run(until=40 * 60))
        again = _anchor(Body(swapped, seed=seed).run(until=40 * 60))
        assert len(first) == 1  # exactly one anchor cell in every run
        winners[first[0]] += 1
        same += first == again
    assert same == 30  # creation order changes nothing: the seed and the cell decide
    assert 8 <= winners["Z1.ppp"] <= 22  # both cells win, about evenly


def test_noise_without_a_seed_is_refused():
    m = parse_file(COLLIER)
    m.organism.seed = None
    with pytest.raises(ValueError, match="needs a seed"):
        Body(m)


TABLE_PROGRAM = """
module toy.contacts
import bio.std.development
organism O { root: Z; cell_type: Progenitor; contacts: contacts.tsv }
stage S { from: 0 min }
cell_type Inner { parent: PostMitotic }
decision split { action: divide; when: cell = Z; daughters: A, B; after: 1 min }
decision split_b { action: divide; when: cell = B; daughters: Ba, Bp; after: 1 min }
signal touching { mode: contact; reads: amount; from: cell_type = Progenitor; sets: senders }
decision inner { action: differentiate; when: cell = A, senders = >=2; to: Inner }
"""


def test_neighbours_come_from_a_time_resolved_table_and_amounts_count_senders(tmp_path):
    (tmp_path / "contacts.tsv").write_text(
        "time\tcell\tcell\tarea\n0\tA\tB\t1\n2\tA\tBa\t1\n2\tA\tBp\t0.5\n2\tBa\tBp\t1\n"
    )
    snaps = read_contacts(str(tmp_path / "contacts.tsv"))
    assert [t for t, _ in snaps] == [0.0, 2.0] and snaps[1][1]["A"] == {"Ba": 1.0, "Bp": 0.5}
    prog = tmp_path / "toy.bio"
    prog.write_text(TABLE_PROGRAM)
    b = Body(parse_file(prog)).run(until=1.5)
    assert b.neighbours(b.cells["A"]) == {"B": 1.0} and b.summary()["neighbours_from"] == "table"
    b.run(until=5)
    assert b.neighbours(b.cells["A"]) == {"Ba": 1.0, "Bp": 0.5}
    assert b.context(b.cells["A"])["senders"] == "1.5"  # area-weighted count of touching senders
    assert b.cells["A"].cell_type == "Progenitor"  # 1.5 is below the rule's 2: nothing is stated


def test_a_missing_contact_table_is_a_compile_error(tmp_path):
    prog = tmp_path / "toy.bio"
    prog.write_text(TABLE_PROGRAM)
    with pytest.raises(SyntaxError, match="contact table"):
        parse_file(prog)


def test_daughters_go_along_the_axis_their_names_imply():
    src = """
module toy.axis
import bio.std.development
organism O { root: Z; cell_type: Progenitor; space: 5 x 5; origin: 2,2; placement: names }
stage S { from: 0 min }
decision split { action: divide; when: cell = Z; daughters: Zl, Zr; after: 1 min }
"""
    b = Body(parse(src)).run(until=2)
    assert (b.cells["Zl"].x, b.cells["Zl"].y) == (2, 2) and (b.cells["Zr"].x, b.cells["Zr"].y) == (2, 3)
    b2 = Body(parse(src.replace("; placement: names", ""))).run(until=2)
    assert (b2.cells["Zr"].x, b2.cells["Zr"].y) == (3, 2)  # nearest free site, as before


def test_replicate_asserts_catch_a_runtime_that_always_picks_the_same_winner():
    smuggled = (
        COLLIER.read_text()
        .replace("when: Dp = >=0.6;", "when: cell = Z1.ppp;")
        .replace("when: Np = >=0.6;", "when: cell = Z4.aaa;")
    )
    m = parse(smuggled)
    bodies = replicate(m, 40 * 60, range(10))
    exactly, share = evaluate_replicate_asserts(bodies, m.organism.asserts)
    assert exactly["ok"] and exactly["winners"] == {"Z1.ppp": 10}  # exactly one anchor cell in every run...
    assert not share["ok"] and share["value"] == 100.0  # ...but always the same one: an order, not a choice


def test_a_revised_fate_is_counted_and_drops_its_old_name():
    """Naming the anchor cell by rule while the circuit decides otherwise: the fate is revised later,
    which is counted (until `commitment` can refuse it) and the cell stops answering to its old name."""
    m = parse(COLLIER.read_text().replace("when: Dp = >=0.6;", "when: cell = Z1.ppp;"))
    revised = 0
    for seed in range(10):
        b = Body(m, seed=seed).run(until=40 * 60)
        z1 = b.cells["Z1.ppp"]
        if z1.cell_type == "VentralUterine":
            revised += 1
            assert z1.terminal_name == "" and b.summary()["revised_fates"] >= 1
    assert revised > 0
