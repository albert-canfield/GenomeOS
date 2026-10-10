# SPDX-License-Identifier: Apache-2.0
"""Stage 3: energy and mitochondria (v0.4 §6, §9.5). The language surface first.

`heteroplasmy:` on a compartment is the one key stage 3 adds. Each refusal here is a modelling error
that would otherwise run: a variant on a compartment with no mitochondrial genome, a fraction that is
not a fraction of copies, two variants claiming more copies than there are.
"""

from __future__ import annotations

import pytest

from genomeos.ir import Allocation, Cost, Pool
from genomeos.lang import BioLangError, parse, parse_file
from genomeos.runtime.economy import Economy
from genomeos.runtime.metabolism import Metabolism, MetabolismError

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


# ---- the runtime: refusals rather than defaults ---------------------------------------------------

DEMO = "data/demo/stage3_energy.bio"
RBC = "data/demo/stage3_erythrocyte.bio"

CELL = """module t
import bio.std.core_metabolism
compartment Cytosol { volume: 0.7 }
compartment Mito { parent: Cytosol; genome: chrM; copies: 10 }
pool ATP { location: Cytosol; size: 1e6; regenerates: from %s }
pool Glucose { location: Cytosol; size: 1; regenerates: 100 /h }
pool O2 { location: Cytosol; size: 1; regenerates: 600 /h }
param oxphos_capacity_per_copy = 100 molecules/h { confidence: 0.1 }
event work { rate: 320 /h; cost: ATP 1 per event }
"""


def cell(processes: str = "Glycolysis, TCA, OXPHOS", text: str = CELL) -> Metabolism:
    return Metabolism.from_module(parse("" + text % processes))


@pytest.mark.parametrize(
    "processes, why",
    [
        ("Glycolysis, Photosynthesis", "runs Glycolysis, TCA, OXPHOS only"),
        ("OXPHOS", "not regenerated from Glycolysis"),
        ("Glycolysis, TCA", "TCA without OXPHOS"),
    ],
)
def test_the_processes_a_pool_names_are_refused_unless_the_module_can_run_them(processes, why):
    with pytest.raises(MetabolismError, match=why):
        cell(processes)


@pytest.mark.parametrize(
    "drop, why",
    [
        ("compartment Mito { parent: Cytosol; genome: chrM; copies: 10 }\n", "needs a mitochondrion"),
        ("pool O2 { location: Cytosol; size: 1; regenerates: 600 /h }\n", "no pool O2"),
        ("pool Glucose { location: Cytosol; size: 1; regenerates: 100 /h }\n", "no pool Glucose"),
        (
            "param oxphos_capacity_per_copy = 100 molecules/h { confidence: 0.1 }\n",
            "does not invent a capacity",
        ),
        ("import bio.std.core_metabolism\n", "holds no stoichiometry"),
    ],
)
def test_what_the_program_does_not_state_is_refused_not_defaulted(drop, why):
    assert drop in CELL
    with pytest.raises(MetabolismError, match=why):
        cell(text=CELL.replace(drop, ""))


def test_stage_2_and_stage_4_demos_compile_but_their_atp_pool_cannot_be_run():
    """They write §5.1's example, `regenerates: from Glycolysis, OXPHOS`, in a cell with no
    mitochondrion and no stoichiometry: compiled as before, refused when asked to run."""
    for path in ("data/demo/stage2_pools.bio", "data/demo/stage4_division.bio"):
        with pytest.raises(MetabolismError):
            Metabolism.from_module(parse_file(path))


def test_a_small_cell_by_hand():
    """10 copies x 100 pairs/h = 1000 pairs/h; 320 ATP/h at 32 per glucose is 10 glucose, 120 pairs,
    60 O2. With no oxygen it needs 160 glucose, over the ceiling of 100: 200 ATP of 320."""
    m = cell()
    hi = m.solve("high")
    assert (hi.supply, hi.glucose, hi.o2, hi.lactate, hi.limit) == (320.0, 10.0, 60.0, 0.0, "demand")
    assert hi.atp_by_process == {"Glycolysis": 20.0, "TCA": 20.0, "OXPHOS": 280.0}
    anoxic = Metabolism.from_module(m.module, supply={"O2": 0.0}).solve("high")
    assert (anoxic.supply, anoxic.glucose, anoxic.lactate, anoxic.factors["work"]) == (
        200.0,
        100.0,
        200.0,
        0.625,
    )
    assert m.available(10.0, "high") == 1e6 + (min(100, 1000 / 12) * 32 + (100 - 1000 / 12) * 2 - 320) * 10


def test_glycolysis_alone_reads_no_oxygen():
    m = cell("Glycolysis")
    assert m.capacity == 0.0
    for o2 in (600.0, 60.0, 0.0):
        f = Metabolism.from_module(m.module, supply={"O2": o2, "Glucose": 1000.0}).solve("low")
        assert (f.supply, f.glucose, f.o2, f.atp_by_process["OXPHOS"]) == (320.0, 160.0, 0.0, 0.0)
        capped = Metabolism.from_module(m.module, supply={"O2": o2}).solve("low")  # 100 glucose/h
        assert (capped.supply, capped.glucose, capped.factors["work"]) == (200.0, 100.0, 0.625)


def test_heteroplasmy_cuts_capacity_not_flux_so_a_reserve_hides_it():
    """Rossignol's biochemical threshold: the flux holds until the copies left cannot carry the
    demand. Here 120 pairs/h are needed of 1000, so the threshold is at h = 0.88."""
    base = cell().solve("high").atp_by_process["OXPHOS"]
    for h, same in ((0.5, True), (0.85, True), (0.9, False), (1.0, False)):
        text = CELL.replace("copies: 10 }", f"copies: 10; heteroplasmy: m.3243A>G = {h} }}")
        got = Metabolism.from_module(parse(text % "Glycolysis, TCA, OXPHOS")).solve("high")
        assert (got.atp_by_process["OXPHOS"] == base) is same, h


def test_the_economy_is_unchanged_without_a_capacity_and_reads_one_when_given():
    pools = {"ATP": Pool(id="ATP", kind="pool", size=100.0)}
    costs = {"a": [Cost(pool="ATP", amount=1.0)], "b": [Cost(pool="ATP", amount=1.0)]}
    allocations = {"ATP": Allocation(id="x", kind="allocation", pool="ATP", policy="proportional")}
    work = {"a": 100.0, "b": 100.0}
    stage2 = Economy(pools=pools, costs=costs, allocations=allocations)
    assert (
        stage2.scale(work) == {"a": 0.5, "b": 0.5}
        and stage2.report(work)["pools"]["ATP"]["capacity"] == 100.0
    )
    stage3 = Economy(pools=pools, costs=costs, allocations=allocations, capacity={"ATP": 400.0})
    assert (
        stage3.scale(work) == {"a": 1.0, "b": 1.0}
        and stage3.report(work)["pools"]["ATP"]["capacity"] == 400.0
    )
    with pytest.raises(ValueError, match="no pool declares"):
        Economy(pools=pools, costs=costs, capacity={"GTP": 1.0}).check()


def test_the_demo_cell_reads_its_declared_numbers():
    m = Metabolism.from_module(parse_file(DEMO))
    assert (m.copies, m.mutant_copies, m.capacity, m.demand()) == (1000.0, 700.0, 6e9, 1e10)
    assert m.yields("high")["o2"] == m.yields("low")["o2"] == 6.0
    rbc = Metabolism.from_module(parse_file(RBC))
    assert rbc.processes == ("Glycolysis",) and rbc.shunt == 0.19


def test_a_supply_names_a_declared_pool_or_is_refused():
    m = cell()
    with pytest.raises(MetabolismError, match="no pool declares"):
        Metabolism.from_module(m.module, supply={"Oxygen": 0.0})
    with pytest.raises(MetabolismError, match="cannot be negative"):
        Metabolism.from_module(m.module, supply={"O2": -1.0})
