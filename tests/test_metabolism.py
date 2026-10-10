# SPDX-License-Identifier: Apache-2.0
"""Stage 3: energy and mitochondria (v0.4 §6, §9.5). The language surface first.

`heteroplasmy:` on a compartment is the one key stage 3 adds. Each refusal here is a modelling error
that would otherwise run: a variant on a compartment with no mitochondrial genome, a fraction that is
not a fraction of copies, two variants claiming more copies than there are.
"""

from __future__ import annotations

import pytest

from genomeos.lang import BioLangError, parse

MITO = "compartment Mito { genome: chrM; copies: 1000; heteroplasmy: %s }"


def compile_one(block: str):
    return parse("module t\n" + block)


def test_heteroplasmy_is_a_fraction_of_copies_per_variant():
    m = compile_one(MITO % "m.3243A>G = 0.7, m.8344A>G = 0.1")
    mito = m.entities["Mito"]
    assert mito.heteroplasmy == {"m.3243A>G": 0.7, "m.8344A>G": 0.1}
    assert mito.mitochondrial and mito.copies == 1000


def test_a_compartment_without_heteroplasmy_states_none():
    m = compile_one(
        "compartment Cytosol { volume: 0.7 }\ncompartment Mito { parent: Cytosol; genome: chrM }\n"
        "compartment Nucleus { parent: Cytosol; genome: nuclear }"
    )
    assert m.entities["Mito"].heteroplasmy == {} and m.entities["Mito"].mitochondrial
    assert not m.entities["Nucleus"].mitochondrial  # `nuclear` excludes chrM


@pytest.mark.parametrize(
    "block, why",
    [
        ("compartment Cytosol { heteroplasmy: m.3243A>G = 0.7 }", "does not read chrM"),
        ("compartment N { genome: nuclear; heteroplasmy: m.3243A>G = 0.7 }", "does not read chrM"),
        (MITO % "m.3243A>G = 1.2", "within 0..1"),
        (MITO % "m.3243A>G = -0.1", "within 0..1"),
        (MITO % "m.3243A>G = 0.7, m.8344A>G = 0.4", "sum to"),
        (MITO % "m.3243A>G = 0.2, m.3243A>G = 0.3", "twice"),
        (MITO % "3243A>G = 0.7", "expected `m.<position>"),
        (MITO % "m.3243A>A = 0.7", "expected `m.<position>"),
        (MITO % "m.3243A>G", "expected `m.<position>"),
        (MITO % "m.0A>G = 0.1", "expected `m.<position>"),
    ],
)
def test_heteroplasmy_refuses_what_is_not_a_fraction_of_mitochondrial_copies(block, why):
    with pytest.raises(BioLangError, match=why):
        compile_one(block)


def test_a_pool_rate_keeps_its_unit_and_is_held_per_hour():
    """Before stage 3 the unit was dropped, so `60 /min` ran as 60 per hour; no runtime read it then."""
    m = compile_one(
        "pool A { size: 1; regenerates: 2e9 /h }\npool B { size: 1; regenerates: 60 /min }\n"
        "pool C { size: 1; regenerates: 24 /d }\npool D { size: 1; regenerates: 5 }"
    )
    assert [m.entities[p].regeneration for p in "ABCD"] == [2e9, 3600.0, 1.0, 5.0]
    with pytest.raises(BioLangError, match="not a rate"):
        compile_one("pool P { size: 1; regenerates: 5 /fortnight }")


def test_processes_are_still_named_not_checked_at_compile_time():
    """Stage 2's and stage 4's demos write `regenerates: from Glycolysis, OXPHOS` with no mitochondrion;
    they must still compile. The stage 3 runtime refuses them when asked to run them."""
    m = compile_one(
        "compartment Cytosol { volume: 1.0 }\n"
        "pool ATP { location: Cytosol; size: 3e9; regenerates: from Glycolysis, OXPHOS }"
    )
    assert m.entities["ATP"].regenerates_from == ["Glycolysis", "OXPHOS"]
