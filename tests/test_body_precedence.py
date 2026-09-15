# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""Four ways the Body decided quietly where it should refuse, report or keep still.

Found by area E using the contact runtime hard, and all of the same family as the fate-precedence
bug of §7.3: a runtime that chooses silently is a runtime whose results cannot be read.

1. Two decisions sharing an id: the one without a `cell` clause was dropped from every candidate
   list and never fired. Shared ids stay legal (mechanism is stated before the generated lookup,
   as the worm's founders are), so they are now reachable and reported instead.
2. A settled fate could be taken back at a later decision point by a decision that had already lost
   the precedence contest — a `cell_network` step gave every losing rule a second bite.
3. A birth makes its neighbours read their contacts again; a death did not.
4. `asymmetric: X -> D` created X in the keeper when the mother had none.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from genomeos.lang import parse, parse_file
from genomeos.runtime.body import Body

ROOT = Path(__file__).resolve().parent.parent
WORM = ROOT / "data" / "organisms" / "celegans" / "embryo_factors.bio"
PLAIN_WORM = ROOT / "data" / "organisms" / "celegans" / "embryo.bio"  # the same worm, no factor rules

SHARED_ID = """
module toy
import bio.std.development
organism T { root: P0; cell_type: Blastomere; factors: M; observe: fates }
stage S { from: 0 min }
timer cycle { duration: 10 min }
decision split { action: divide; when: cell = P0; daughters: Pa, Pp; asymmetric: M -> a, N -> p }
decision mark { action: express; when: cell = P0; sets: K }
decision mark { action: express; when: cell_type = Blastomere; sets: G }
"""

CONTACT = """
module toy
import bio.std.development
cell_type Touched { parent: PostMitotic }
cell_type Orphan { parent: PostMitotic }
organism T { root: P0; cell_type: Zygote; space: 4 x 1; placement: names; observe: fates }
stage S { from: 0 min }
timer cycle { duration: 10 min; when: generation = 0 }
signal Touch { mode: contact; ligand: L; reads: amount; from: cell_type = Blastomere; sets: Lext }
decision start { action: differentiate; when: cell_type = Zygote; to: Blastomere }
decision split { action: divide; when: cell = P0; daughters: Pa, Pp }
decision touched { action: differentiate; when: cell = Pa, Lext = >=1; to: Touched }
decision orphaned { action: differentiate; when: cell = Pa, cell_type = Touched, Lext = <1; to: Orphan }
decision kill { action: die; when: cell = Pp; after: 30 min }
"""


def test_a_decision_sharing_an_id_is_reachable_and_the_sharing_is_reported():
    body = Body(parse(SHARED_ID), seed=None).run(until=40)
    root = body.cells["P0"]
    assert "K" in root.factors  # the cell-specific one
    assert "G" in root.factors  # the general one, which used to be dropped from every candidate list
    assert body.summary()["duplicate_decision_ids"] == {"mark": 2}


def test_asymmetric_segregates_a_factor_and_never_creates_one():
    body = Body(parse(SHARED_ID), seed=None).run(until=40)
    assert body.cells["Pa"].factors.get("M") == "present"  # the mother had it; the keeper keeps it
    assert "M" not in body.cells["Pp"].factors
    assert "N" not in body.cells["Pp"].factors  # the mother had none, so no daughter may have it
    assert "N" not in body.cells["Pa"].factors


def test_a_death_makes_its_neighbours_read_their_contacts_again():
    body = Body(parse(CONTACT), seed=None).run(until=60)
    kinds = {c.name: c.cell_type for c in body.alive_at(60)}
    assert kinds == {"Pa": "Orphan"}  # Touched while Pp lived, Orphan once it died


@pytest.fixture(scope="module")
def worm():
    return parse_file(WORM)


def _fate_score(module, cadence: int):
    from genomeos.organism.diff import compare
    from genomeos.organism.reference import ReferenceLineage

    module.organism.cell_network = cadence
    body = Body(module, seed=None).run(until=800)
    diff = compare(body, ReferenceLineage.load(), until=800)
    return diff.fates_correct, diff.fates_checked, body.summary()


def test_a_network_cadence_does_not_change_which_fates_are_taken(worm):
    """Area E measured the worm's Sulston fate score falling when a network was stepped between
    events: every step re-decided every cell, and a rule that had lost on precedence took the fate
    back. `fates: first` means one fate per decision point, and a later decision point is not an
    appeal, so the score must not depend on the cadence."""
    still = _fate_score(worm, 0)
    stepped = _fate_score(worm, 6)
    assert stepped[0] == still[0] and stepped[1] == still[1]
    assert stepped[2]["network_steps"] > 0  # the network really ran
    assert stepped[2]["overruled_fates"] > 0  # and losing rules really tried


READS = """
module toy
import bio.std.development
cell_type Deep { parent: PostMitotic }
cell_type Shallow { parent: PostMitotic }
organism T { root: P0; cell_type: Blastomere; factors: M; observe: fates }
stage S { from: 0 min }
timer cycle { duration: 10 min }
decision d0 { action: divide; when: cell = P0; daughters: Pa, Pp; asymmetric: M -> a }
decision d1 { action: divide; when: generation = 1 }
decision deep { action: differentiate; when: generation = 2, M.exposure(lineage) >= 15; to: Deep }
decision shallow { action: differentiate; when: generation = 2, M.exposure(lineage) = <15, M.mean(lineage) = <1, M.mean(cell) = <1; to: Shallow }
"""  # noqa: E501


def test_exposure_and_mean_are_read_over_the_window_the_program_names():
    """v0.4 §7.2a: the three readings of a factor give different answers, so each names its window.
    M segregates to Pa at the first division, so the two halves of the lineage differ by construction:
    the Pa side carries it for the whole path, the Pp side only for the 10 minutes P0 held it."""
    body = Body(parse(READS), seed=None).run(until=40)
    kinds = {c.name: c.cell_type for c in body.alive_at(40)}
    assert kinds == {"Pal": "Deep", "Par": "Deep", "Ppl": "Shallow", "Ppr": "Shallow"}
    read = {n: body.context(body.cells[n], 40) for n in ("P0", "Pa", "Pal", "Ppl")}
    assert read["Pal"]["M.exposure(lineage)"] == "40"  # 10 in P0, 10 in Pa, 20 of its own
    assert read["Ppl"]["M.exposure(lineage)"] == "10"  # only what P0 carried before the division
    assert read["P0"]["M.exposure(lineage)"] == "10"  # a cell stops accumulating when it divides
    assert float(read["Pal"]["M.mean(lineage)"]) == 1.0
    assert float(read["Ppl"]["M.mean(lineage)"]) == 0.25


def test_a_program_that_asks_for_no_integrated_read_computes_none():
    """A program pays for an integrated read only by naming one. `embryo_factors.bio` asks for 33 of
    them since area E rewrote its fate rules on this runtime, so the program that pins this is the one
    without factor rules at all."""
    assert Body(parse_file(PLAIN_WORM), seed=None)._reads == []
    assert Body(parse_file(WORM), seed=None)._reads  # and the one that asks does get them


def test_an_equal_or_higher_precedence_reading_may_still_change_a_fate():
    """The refusal is precedence, not a freeze: a new reading by a rule that outranks the one that
    settled the fate still applies, which is how a contact signal changes a cell's mind."""
    source = CONTACT.replace("to: Orphan }", "to: Orphan; priority: 9 }").replace(
        "to: Touched }", "to: Touched; priority: 1 }"
    )
    body = Body(parse(source), seed=None).run(until=60)
    assert {c.name: c.cell_type for c in body.alive_at(60)} == {"Pa": "Orphan"}
    assert body.summary()["overruled_fates"] == 0
