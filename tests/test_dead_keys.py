# SPDX-License-Identifier: Apache-2.0
"""The engine's two allegedly dead keys: `policy: optimise` and `regenerates:` (ROADMAP §5 item 13 A).

WHY THIS EXISTS (2026-10-10, lane-deadkeys). The open debt reads: "`policy: optimise` compiles and
the runtime refuses it, and `regenerates:` is stored and read by nothing". This file is the lane's
decision, written down before the implementation, and the first-hand read behind it. Only the
engine's own keys are involved: no biosample, nothing in area J, no node model, no AlphaGenome
request, no BioForge calibration.

THE DECISION, PER KEY.

**`policy: optimise` -> outcome 2: refused at CHECK time, not at run time.** First-hand: the grammar
lists it among the allocation policies and documents `objective:`; the parser requires an objective
for it; the IR carries it in `REGIME_ALLOCATIONS`; and `runtime/economy.py` `_share` raises
`EconomyError` on it as "a declared modelling device". So the roadmap's sentence is right about this
key: a program could be written, pass `bio check`, and then die on its first allocation. No committed
`.bio` program writes it (grep for `optimise` over `data/organisms`, `data/demo` and `genomeos/std`:
zero hits on 2026-10-10), so refusing it at check time changes no program's behaviour. The engine
does not solve an objective, so a program that names one is refused where the mistake is -- at the
check -- with the policy named. Removing the key from the language is NOT done here: that is Albert's
decision, and the refusal is what makes the gap honest in the meantime. The run-time refusal stays,
because the IR can still be built in Python without the parser (`tests/test_economy.py` pins it).

**`regenerates:` -> THE ROADMAP'S SENTENCE IS WRONG ABOUT THIS KEY.** It is read in two of its three
positions, and the one dead position is a language decision, not an engine bug. First-hand, the key
has three positions:

  1. `regenerates: from <Process>, <Process>` -> READ. `lang/parser.py` stores `regenerates_from` and
     `runtime/metabolism.py` reads it: it is how the module finds the one regenerated pool and the
     processes that refill it, and it refuses zero such pools and more than one.
  2. `regenerates: <rate> /h` on the `Glucose` or `O2` pool -> READ. `metabolism._rate` takes the
     pool's `regeneration` as the glucose uptake ceiling and the oxygen supply, the two numbers that
     bound oxidation; `scripts/metabolism_gate.py` reads the same two off the demo programs.
  3. `regenerates: <rate> /h` on any other pool -> READ BY NOTHING. `runtime/economy.py` holds a
     pool's demand against `capacity_of()`, which is a capacity its CALLER passed or else the pool's
     `size`; no code path anywhere reads `Pool.regeneration` for a pool that is not Glucose or O2.
     Exactly one committed program writes one: `data/demo/stage2_pools.bio`'s
     `pool Ribosomes { size: 5e6; regenerates: 1e5 /h }`.

Outcome for position 3: **outcome 3, a language decision for Albert**, recorded and pinned rather
than changed. Both ways of closing it change what is already committed. Making the rate the capacity
would hold that demo's ribosome demand against 1e5/h instead of 5e6, a fiftyfold tighter ceiling, so
a committed program's burden would move; refusing the rate at check time would make a committed demo
fail `bio check`. And the choice between them IS the open question -- what a per-hour refill rate
means for a pool whose demand is held against a standing stock, i.e. what interval a stage-2 draw
covers. The engine should not answer that quietly. The test below is the falsifier that keeps the
position honest: it fails the day anything starts reading the rate.
"""

from __future__ import annotations

import pytest

from genomeos.ir import UNKNOWN
from genomeos.lang import BioLangError, parse
from genomeos.runtime.economy import Economy, EconomyError
from genomeos.runtime.metabolism import Metabolism

CELL = """module t
import bio.std.core_metabolism
compartment Cytosol { volume: 0.7 }
compartment Mito { parent: Cytosol; genome: chrM; copies: 10 }
pool ATP { location: Cytosol; size: 1e6; regenerates: from %s }
pool Glucose { location: Cytosol; size: 1; regenerates: %s /h }
pool O2 { location: Cytosol; size: 1; regenerates: %s /h }
param oxphos_capacity_per_copy = 100 molecules/h { confidence: 0.1 }
event work { rate: 320 /h; cost: ATP 1 per event }
"""

#: A pool divided by `optimise`, with a demand above its capacity so the policy is actually applied.
OPTIMISE = """module t
pool R { size: 1000 }
allocation a { pool: R; policy: optimise%s }
event work { rate: 1500 /h; cost: R 1 per event }
"""

#: A pool that is neither Glucose nor O2, with both a stock and a refill rate: position 3 above.
STOCK = """module t
compartment Cytosol { volume: 0.7 }
pool Ribosomes { location: Cytosol; size: 100; %s returns_as: Ribosomes }
event work { rate: 1000 /h; cost: Ribosomes 1 per event }
"""


def cell(processes: str = "Glycolysis, TCA, OXPHOS", glucose: str = "100", o2: str = "600") -> Metabolism:
    return Metabolism.from_module(parse(CELL % (processes, glucose, o2)))


# ---- `policy: optimise` ---------------------------------------------------------------------


def test_the_debt_optimise_compiles_today_and_only_the_runtime_refuses_it() -> None:
    """THIS TEST IS THE DEBT ITSELF, committed before the fix so the gap is on the record as a fact
    and not as a sentence. A program naming a policy the engine cannot run passes the check and dies
    later. The next commit on this file replaces it with the check-time refusal the decision chose."""
    module = parse(OPTIMISE % "; objective: growth")

    assert module.entities["a"].policy == "optimise", "the check accepts it today"
    with pytest.raises(EconomyError, match="does not solve an objective"):
        Economy.from_module(module).scale({"work": 1500.0})


def test_optimise_still_needs_an_objective_whatever_else_changes() -> None:
    """Pinned here as well as in `tests/test_pools.py`: whichever refusal fires first, the missing
    parameter keeps its own message, because "you named no objective" and "the engine solves no
    objective" are different mistakes."""
    with pytest.raises(BioLangError, match="declare its objective"):
        parse(OPTIMISE % "")


# ---- `regenerates:` -------------------------------------------------------------------------


def test_regenerates_from_is_read_and_the_processes_it_names_reach_the_runtime() -> None:
    """Position 1. The roadmap says nothing reads `regenerates:`; the module is built out of it."""
    assert cell().processes == ("Glycolysis", "TCA", "OXPHOS")
    assert cell(processes="Glycolysis").processes == ("Glycolysis",)
    with pytest.raises(ValueError, match="no pool is regenerated"):
        Metabolism.from_module(parse("module t\npool ATP { size: 1e6 }"))


def test_a_regenerates_rate_on_glucose_or_o2_is_read_as_the_supply_that_bounds_oxidation() -> None:
    """Position 2. Both numbers move the run, so both are read: the glucose rate is the uptake
    ceiling and the O2 rate is the oxygen supply, and the arm's limit names whichever binds."""
    assert (cell().glucose_uptake, cell().o2_supply) == (100.0, 600.0)
    assert cell(glucose="7").glucose_uptake == 7.0
    assert cell(o2="6").solve("high").limit == "oxygen"
    assert cell(glucose="1", o2="6000").solve("high").limit == "glucose"


@pytest.mark.parametrize("rate", ["", "regenerates: 1e5 /h;", "regenerates: 1e9 /h;", "regenerates: 1 /h;"])
def test_a_regenerates_rate_on_any_other_pool_is_read_by_nothing(rate: str) -> None:
    """Position 3, the one dead position, kept as a falsifier rather than fixed (see the module
    docstring: closing it is Albert's decision). The pool is neither Glucose nor O2, so no metabolism
    reads it, and the economy holds its demand against its stock: 1000 draws on a stock of 100 are
    cut to a tenth whatever the declared refill, and the whole report is identical. The day anything
    reads the rate -- an interval capacity, a refusal, anything -- this test fails and that is the
    point, because it would also move `data/demo/stage2_pools.bio`."""
    module = parse(STOCK % rate)
    pool = module.entities["Ribosomes"]
    assert (pool.regeneration is UNKNOWN) == (rate == ""), "the parser stores the rate either way"

    economy = Economy.from_module(module)

    assert economy.capacity_of("Ribosomes") == 100.0
    assert economy.scale({"work": 1000.0})["work"] == pytest.approx(0.1)
    assert economy.report({"work": 1000.0}) == Economy.from_module(parse(STOCK % "")).report({"work": 1000.0})
