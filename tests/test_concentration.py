# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
# ruff: noqa: E501  (BioLang source kept one block per line)
"""A compartment's absolute volume, and the concentration thresholds it makes expressible.

BioLang v0.4 §4.1 and §10 decision 2, resolved by Albert on 2026-09-19: a `compartment` gains an
absolute volume *beside* the fraction it already carries. The fraction is what the containment tree
needs and what the committed programs state; the absolute volume is what a threshold can be divided
by, and a fraction of a cell cannot express cells of different volume, which is why §9's registered
test ("run the gate both ways on cells of different volume") could not be run at all.

What is demonstrated here is that one claim: the same rule, written once in molecules and once in
nanomolar, is the *same* statement in the cell the two were calibrated on and a *different* statement
in a cell of another volume. The second half is the part the language could not say before.
"""

from __future__ import annotations

import pytest

from genomeos.ir import UNKNOWN, Module, molecules_in, to_femtolitres, to_molar
from genomeos.lang import BioLangError, parse
from genomeos.runtime import NetworkRuntime
from genomeos.runtime.located import LocatedRuntime

# One cell, one fixed repressor, and one target regulated twice: T's threshold is a concentration and
# U's is the amount that concentration comes to in 94 fL (100 nM * 94 fL * Avogadro = 5661 molecules).
# Only the volume changes between the two instances of this template.
CELL = """
module demo.two_units

regime counted {{ units: copies }}

compartment Cytosol {{ translation: yes; volume: 1.0; absolute_volume: {volume}
  evidence: experimental "a cytosol"; confidence: 0.7 }}

protein Rp {{ location: Cytosol; initial: 20000; half_life: 1e6 h }}
gene T {{ location: Cytosol; basal: 1; max: 1000; produces: Tp }}
protein Tp {{ location: Cytosol }}
gene U {{ location: Cytosol; basal: 1; max: 1000; produces: Up }}
protein Up {{ location: Cytosol }}

rule Rp inhibits T {{ threshold: 100 nM; hill: 8; strength: 1.0 }}
rule Rp inhibits U {{ threshold: 5661; hill: 8; strength: 1.0 }}
"""

REPRESSED, FREE = 5.0, 900.0  # basal alone is 1; basal + max is 1001


def _run(volume: str) -> dict[str, float]:
    module = parse(CELL.format(volume=volume))
    res = LocatedRuntime(module).run(hours=48.0, dt=0.05)
    return {"concentration": res.final("Tp@Cytosol"), "amount": res.final("Up@Cytosol")}


# ---- the field ---------------------------------------------------------------------------


def test_a_compartment_states_an_absolute_volume_beside_its_fraction():
    m = parse(
        "module t\n"
        "compartment Cytosol { volume: 0.54; absolute_volume: 2.7 pL\n"
        '  evidence: curated "x"; confidence: 0.5 }\n'
        "compartment Mito { parent: Cytosol; absolute_volume: 300 um3 }\n"
    )
    assert m.entities["Cytosol"].volume == 0.54  # the fraction keeps its meaning and its name
    assert m.entities["Cytosol"].absolute_volume_fl == 2700.0  # 2.7 pL
    assert m.entities["Mito"].absolute_volume_fl == 300.0  # 1 um3 is exactly 1 fL
    assert m.entities["Mito"].volume is UNKNOWN


def test_an_absolute_volume_survives_the_json_round_trip():
    m = parse("module t\ncompartment C { volume: 1.0; absolute_volume: 94 fL }\n")
    back = Module.from_dict(m.to_dict())
    assert back.entities["C"].absolute_volume_fl == 94.0
    blank = Module.from_dict(parse("module t\ncompartment C { }\n").to_dict())
    assert blank.entities["C"].absolute_volume_fl is UNKNOWN


def test_a_volume_without_a_unit_is_refused_because_the_fraction_is_already_a_bare_number():
    with pytest.raises(BioLangError, match="needs a number and a unit"):
        parse("module t\ncompartment C { absolute_volume: 94 }\n")
    with pytest.raises(BioLangError, match="unknown volume unit"):
        parse("module t\ncompartment C { absolute_volume: 94 litres }\n")
    with pytest.raises(BioLangError, match="greater than zero"):
        parse("module t\ncompartment C { absolute_volume: -1 fL }\n")


def test_unknown_is_a_statement_a_program_may_make():
    m = parse("module t\ncompartment C { volume: 0.54; absolute_volume: unknown }\n")
    assert m.entities["C"].absolute_volume_fl is UNKNOWN


def test_the_unit_tables_agree_with_the_conversion_they_are_for():
    assert to_femtolitres(1.0, "um3") == to_femtolitres(1.0, "fL") == 1.0
    assert to_femtolitres(1.0, "L") == 1e15
    assert to_molar(100.0, "nM") == pytest.approx(1e-7)
    # 100 nM in 94 fL, the number the committed erythrocyte's volume gives
    assert molecules_in(to_molar(100.0, "nM"), 94.0) == pytest.approx(5661.0, rel=1e-3)


# ---- what it makes possible -------------------------------------------------------------


def test_a_concentration_threshold_runs_and_is_converted_at_the_compartment():
    module = parse(CELL.format(volume="94 fL"))
    vm = LocatedRuntime(module)
    assert vm.thresholds[("Rp inhibits T", "Cytosol")] == pytest.approx(5661.0, rel=1e-3)
    assert vm.thresholds[("Rp inhibits U", "Cytosol")] == 5661.0  # an amount passes through untouched


def test_the_two_units_agree_in_the_cell_they_were_calibrated_on():
    small = _run("94 fL")
    # 20,000 molecules is 353 nM here, above 100 nM however the threshold is written, so both are off
    assert small["concentration"] < REPRESSED and small["amount"] < REPRESSED
    # they agree to the 3e-5 by which the program's readable 5661 differs from the exact 5660.81
    assert small["concentration"] == pytest.approx(small["amount"], rel=1e-4)


def test_the_two_units_disagree_across_cells_of_different_volume():
    """§9's gate, which could not be run before: the same 20,000 molecules are 353 nM in a red-cell
    volume and 35.3 nM in a cell ten times larger, so a 100 nM threshold represses in one and not in
    the other while a threshold of 5,661 molecules represses in both."""
    small, large = _run("94 fL"), _run("940 fL")
    assert large["amount"] == pytest.approx(small["amount"], rel=1e-9)  # an amount does not notice
    assert large["concentration"] > FREE  # the concentration does
    assert large["concentration"] / large["amount"] > 100.0


# ---- what is refused rather than guessed ------------------------------------------------


def test_a_concentration_needs_a_volume_to_divide_by_and_says_so_at_compile_time():
    with pytest.raises(BioLangError, match="declares no absolute_volume"):
        parse(
            "module t\n"
            "compartment Cytosol { translation: yes; volume: 1.0 }\n"
            "protein Rp { location: Cytosol; initial: 10 }\n"
            "gene T { location: Cytosol; basal: 1; max: 10; produces: Tp }\n"
            "protein Tp { location: Cytosol }\n"
            "rule Rp inhibits T { threshold: 100 nM }\n"
        )


def test_a_threshold_unit_that_is_not_a_concentration_is_refused():
    with pytest.raises(BioLangError, match="amount .no unit. or a concentration"):
        parse(
            "module t\ncompartment C { absolute_volume: 1 pL }\n"
            "gene T { location: C; basal: 1 }\nprotein P { location: C }\n"
            "rule P inhibits T { threshold: 5 mg }\n"
        )


def test_the_unlocated_runtime_refuses_a_concentration_because_it_has_no_places():
    module = parse(
        "module t\ngene T { basal: 1; max: 10; produces: Tp }\nprotein Tp { }\nprotein Rp { }\n"
        "rule Rp inhibits T { threshold: 100 nM }\n"
    )
    with pytest.raises(ValueError, match="needs a compartment with an absolute_volume"):
        NetworkRuntime(module)


def test_a_concentration_needs_counted_units_not_arbitrary_ones():
    module = parse(CELL.format(volume="94 fL").replace("units: copies", "units: au"))
    with pytest.raises(ValueError, match="units: copies"):
        LocatedRuntime(module)
