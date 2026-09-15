# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""BioLang v0.4 §7.2a: `competence`, `commitment` and a perturbation that arrives at a stated time.

The falsifier is the published series (Fukushige & Krause 2005; Yuzyuk et al. 2009) written as
data/organisms/celegans/plasticity.bio: forcing HLH-1 inside the window converts the embryo, forcing it
after the window does nothing, removing MES-2 reopens the window, and cells that have already
differentiated never respond. The tests below check the runtime pieces separately and then that each
construct is load-bearing: strip it and an arm of the series stops reproducing.
"""

from __future__ import annotations

import copy
from pathlib import Path

import pytest

from genomeos.lang import parse, parse_file
from genomeos.lang.parser import BioLangError
from genomeos.runtime.body import Body, evaluate_assert

ROOT = Path(__file__).resolve().parent.parent
PROGRAM = ROOT / "data" / "organisms" / "celegans" / "plasticity.bio"

TOY = """
module toy
import bio.std.development
cell_type Muscle { parent: PostMitotic }
organism Toy { root: P0; cell_type: Zygote; observe: fates }
stage Early { from: 0 min; to: 100 min }
stage Late { from: 100 min }
timer cycle { duration: 20 min }
decision grow { action: divide; when: cell_type = any }
competence window { allows: Muscle; closes: at 50 min; closed_by: MES-2 }
decision make_muscle { action: differentiate; when: X = present; to: Muscle; competence: window }
"""


def toy(extra: str = "") -> object:
    return parse(TOY + extra)


# ---- the constructs compile, and refuse what they cannot do -------------------------


def test_a_window_must_say_when_it_closes():
    with pytest.raises(BioLangError, match="not a window"):
        parse(TOY.replace("closes: at 50 min; ", ""))


def test_a_window_that_gates_nothing_is_refused():
    with pytest.raises(BioLangError, match="governs no decision"):
        parse(TOY.replace("competence: window", "priority: 1"))


def test_a_decision_may_not_need_an_undeclared_window():
    with pytest.raises(BioLangError, match="undeclared competence"):
        parse(TOY.replace("competence: window", "competence: other"))


def test_a_window_must_allow_the_fate_it_gates():
    with pytest.raises(BioLangError, match="does not allow"):
        parse(TOY.replace("allows: Muscle", "allows: PostMitotic"))


def test_clauses_specified_but_not_implemented_are_compile_errors():
    for clause in ("maintain: X >= 0.2", "excludes: Y", "hysteresis: enter 0.8, leave 0.3"):
        with pytest.raises(BioLangError, match="not implemented"):
            parse(TOY + f"\ncommitment c {{ establish: cell_type = Muscle; {clause} }}\n")


def test_an_integrated_read_must_name_its_window():
    parse(TOY + "\ncommitment c { establish: ELT-2.exposure(lineage) >= 0.8 }\n")  # compiles
    for bad in ("ELT-2.exposure >= 0.8", "ELT-2.mean(path) >= 0.8", "ELT-2.exposure(both) >= 0.8"):
        with pytest.raises(BioLangError, match="must name its window"):
            parse(TOY + f"\ncommitment c {{ establish: {bad} }}\n")


# ---- what the runtime does with them ------------------------------------------------


def test_a_factor_can_be_forced_at_a_time_and_is_then_inherited():
    module = toy()
    body = Body(module, seed=None, add_at={"X": 30.0}).run(until=90)
    assert body.summary()["forced"] == {"X": len(body.alive_at(30.0))}
    assert all(c.factors.get("X") == "present" for c in body.alive_at(90))


def test_the_window_closes_and_the_fate_is_refused_not_ignored():
    module = toy()
    inside = Body(module, seed=None, adds={"MES-2"}, add_at={"X": 30.0}).run(until=90)
    outside = Body(module, seed=None, adds={"MES-2"}, add_at={"X": 60.0}).run(until=90)
    assert any(c.cell_type == "Muscle" for c in inside.alive_at(90))
    assert not any(c.cell_type == "Muscle" for c in outside.alive_at(90))
    assert outside.summary()["outside_competence"] > 0  # refused and counted, which is not UNKNOWN


def test_without_the_closing_machinery_the_window_stays_open():
    module = toy()
    late = Body(module, seed=None, add_at={"X": 60.0}).run(until=90)  # no MES-2 at all
    assert any(c.cell_type == "Muscle" for c in late.alive_at(90))


def test_a_committed_cell_refuses_a_later_fate_and_its_daughters_inherit_the_lock():
    extra = (
        "\ncell_type Other { parent: PostMitotic }"
        "\ncommitment locked { establish: cell_type = Muscle; inherit: daughters; release: never }"
        "\ndecision make_other { action: differentiate; when: Y = present; to: Other; priority: 5 }\n"
    )
    module = parse(TOY + extra)
    body = Body(module, seed=None, adds={"MES-2"}, add_at={"X": 30.0, "Y": 70.0}).run(until=90)
    committed = body.summary()["committed"]
    assert list(committed) == ["Muscle"] and committed["Muscle"] > 0
    assert body.summary()["refused_committed"] > 0
    assert not any(c.cell_type == "Other" for c in body.alive_at(90))
    assert all(c.committed == "Muscle" for c in body.alive_at(90))  # inherited, not re-established


# ---- the published series, and what each construct is worth --------------------------


@pytest.fixture(scope="module")
def series():
    return parse_file(PROGRAM)


def arms(module) -> dict[str, bool]:
    out = {}
    for ex in module.experiments:
        body = Body(
            module,
            seed=None,
            knockouts=set(ex.knockouts),
            adds=set(ex.adds),
            add_at=dict(ex.add_at),
            environment=ex.environment,
        ).run(until=ex.until)
        out[ex.name] = all(evaluate_assert(body, a)["ok"] for a in ex.asserts)
    return out


def test_the_published_series_is_reproduced(series):
    assert arms(series) == {ex.name: True for ex in series.experiments}


def test_removing_competence_loses_the_closed_window(series):
    stripped = copy.deepcopy(series)
    stripped.competences.clear()
    for d in stripped.decisions:
        d.competence = ""
    result = arms(stripped)
    assert result["hlh1_late"] is False  # the late factor converts the embryo again
    assert result["hlh1_late_pha4"] is False


def test_removing_commitment_loses_the_cells_that_already_differentiated(series):
    stripped = copy.deepcopy(series)
    stripped.commitments.clear()
    result = arms(stripped)
    assert result["hlh1_terminal_mes2"] is False
    assert result["hlh1_early"] is False  # and the converted cells do not stay converted


# ---- the hole in the last arm, and what closed it ------------------------------------


def _born_after_the_induction(module) -> dict[str, str]:
    ex = next(e for e in module.experiments if e.name == "hlh1_terminal_mes2")
    body = Body(
        module,
        seed=None,
        knockouts=set(ex.knockouts),
        adds=set(ex.adds),
        add_at=dict(ex.add_at),
    ).run(until=ex.until)
    return {n: c.cell_type for n, c in body.cells.items() if c.born > 700}


def test_a_cell_born_after_the_induction_is_protected_by_what_its_parent_committed_to(series):
    """The documented hole in this arm: `commitment` establishes when a cell differentiates, so the
    twelve cells born after the factor arrives had nothing to inherit and took the forced fate. They
    are protected now, because their parents reach a terminal fate before dividing and the lock is
    inherited — and the ablation is per block, because a coarse one named the wrong cause once:
    either block alone is enough, since both establish on a terminal type."""
    late = _born_after_the_induction(series)
    assert len(late) == 12
    assert [t for t in late.values() if t == "Muscle"] == []
    names = [c.name for c in series.commitments]
    assert len(names) > 1, "the per-block ablation below needs more than one block to be meaningful"
    for drop in names:
        one = copy.deepcopy(series)
        one.commitments = [c for c in one.commitments if c.name != drop]
        assert [t for t in _born_after_the_induction(one).values() if t == "Muscle"] == [], drop
    both = copy.deepcopy(series)
    both.commitments.clear()
    assert set(_born_after_the_induction(both).values()) == {"Muscle"}  # all twelve, with none left


def test_a_precursor_that_never_differentiates_can_commit_on_a_sustained_read():
    """The hole is closed on this program but the shape of it is real: a cell that runs a programme
    for a long time without differentiating has nothing to pass on. A program can say so itself, with
    no new construct — `establish` on an integrated read is §7.2a's own example, and it compiles since
    the reads landed. Here the precursor divides for ever and never takes a fate."""
    source = """
module toy
import bio.std.development
cell_type Gut { parent: PostMitotic }
cell_type Muscle { parent: PostMitotic }
organism T { root: P0; cell_type: Blastomere; factors: ELT-2; observe: fates }
stage S { from: 0 min }
timer cycle { duration: 20 min }
decision grow { action: divide; when: cell_type = Blastomere }
decision forced { action: differentiate; when: HLH-1 = present; to: Muscle; priority: 9 }
"""
    precursor = (
        "commitment gut_programme { programme: Gut; establish: ELT-2.exposure(lineage) >= 30; "
        "inherit: daughters; release: never }\n"
    )
    for source_text, expected in ((source, "Muscle"), (source + precursor, "Blastomere")):
        body = Body(parse(source_text), seed=None, add_at={"HLH-1": 70.0}).run(until=120)
        kinds = {c.cell_type for c in body.alive_at(120)}
        assert kinds == {expected}, source_text[-80:]
