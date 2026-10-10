# SPDX-License-Identifier: AGPL-3.0-or-later
"""Run stage 3's gate against its registration (runtime/metabolism_gate.py). No network, no quota.

    uv run python scripts/metabolism_gate.py [--no-save]

The design is `PRE_REGISTRATION`, committed in `genomeos/runtime/metabolism_gate.py` (e09edf7) before
the runtime, the standard-library module, the two demo programs and this file existed. Nothing here may
change it, and the script refuses to run when a demo program declares a number the registration does
not. Every prediction it is held to was worked by hand into that file; this script only measures.

Registered there, before any run, and therefore not findings of this script: the hierarchy of
Buttgereit & Brand is reproduced only as a priority order written into the program (proportional and
competitive hand every consumer one scalar); the Pasteur clause holds by construction of the routing;
the heteroplasmy threshold is the declared reserve; every cell here collapses without glucose.
"""

from __future__ import annotations

import argparse
import copy
import sys
import time
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from genomeos import manifest as mf  # noqa: E402
from genomeos.ir import Allocation, Compartment, Event, Pool  # noqa: E402
from genomeos.lang import BioLangError, parse, parse_file  # noqa: E402
from genomeos.runtime.division import Division  # noqa: E402
from genomeos.runtime.economy import Economy  # noqa: E402
from genomeos.runtime.metabolism import ARMS, Metabolism, MetabolismError  # noqa: E402
from genomeos.runtime.metabolism_gate import (  # noqa: E402
    ATP_STOCK,
    BPG_SHUNT,
    CAPACITY_PER_COPY,
    CLOSED_FORM_RTOL,
    CONSUMERS,
    DEMAND_PER_H,
    DEMO_HETEROPLASMY,
    DIVISION_COST_ATP,
    GLUCOSE_UPTAKE_PER_H,
    HETEROPLASMY_SWEEP,
    INTERVAL_H,
    MACROMOLECULE,
    MITO_COPIES,
    O2_SUPPLY_PER_H,
    PRE_REGISTRATION,
    PREDICTED_ANOXIC_FACTORS,
    PREDICTED_AVAILABLE,
    PREDICTED_GLUCOSE,
    PREDICTED_MAX_SUPPLY,
    PREDICTED_O2_PER_GLUCOSE,
    PREDICTED_OXPHOS_SHARE,
    PREDICTED_PAIRS_PER_GLUCOSE,
    PREDICTED_PASTEUR_RATIO,
    PREDICTED_RBC_GLUCOSE,
    PREDICTED_RBC_LACTATE,
    PREDICTED_RBC_SHUNT_FLUX,
    PREDICTED_RBC_YIELD,
    PREDICTED_SUPPLY,
    PREDICTED_THRESHOLD,
    PREDICTED_YIELD_AEROBIC,
    PREDICTED_YIELD_ANAEROBIC,
    PREDICTED_YIELD_SPEC_LOW_END,
    PRIORITY_ORDER,
    PROGRAM,
    RBC_CONSUMERS,
    RBC_DEMAND_PER_H,
    RBC_GLUCOSE_UPTAKE_PER_H,
    RBC_PROGRAM,
    SUPPLY_SWEEP,
    UNCAPPED_GLUCOSE_PER_H,
)

RESULT = "metabolism_gate"
ENTRY = "scripts/metabolism_gate.py"
STAGE4_PROGRAM = "data/demo/stage4_division.bio"
STD = "genomeos/std/core_metabolism.bio"
OWN_CODE = (
    "scripts/metabolism_gate.py",
    "genomeos/runtime/metabolism.py",
    "genomeos/runtime/metabolism_gate.py",
    "genomeos/runtime/economy.py",
    STD,
    PROGRAM,
    RBC_PROGRAM,
    "tests/test_metabolism.py",
    "tests/test_metabolism_gate.py",
)


def close(got: float, want: float, rtol: float = CLOSED_FORM_RTOL) -> bool:
    """Within a relative `rtol`; a zero must be exactly zero."""
    return abs(got - want) <= rtol * max(abs(want), 1e-300)


# ---- the programs against the registration ---------------------------------------------------------


def disagreements(cell: Any, rbc: Any) -> list[str]:
    """Every declared number of the two programs that the registration does not hold."""
    out: list[str] = []

    def want(label: str, got: Any, expected: Any) -> None:
        if got != expected:
            out.append(f"{label}: program {got!r}, registered {expected!r}")

    for module, consumers, label in ((cell, CONSUMERS, PROGRAM), (rbc, RBC_CONSUMERS, RBC_PROGRAM)):
        events = {e.id: e for e in module.events if isinstance(e, Event) and e.costs}
        want(f"{label} consumers", sorted(events), sorted(consumers))
        for name, spec in consumers.items():
            if name in events:
                e = events[name]
                want(f"{label} {name} rate", (e.rate, e.rate_unit), (spec["rate_per_h"], "/h"))
                want(
                    f"{label} {name} cost",
                    [(c.pool, c.amount, c.per) for c in e.costs],
                    [("ATP", spec["atp_per_event"], "event")],
                )
    ent = cell.entities
    want("Glucose uptake", ent["Glucose"].regeneration, GLUCOSE_UPTAKE_PER_H)
    want("O2 supply", ent["O2"].regeneration, O2_SUPPLY_PER_H)
    want("ATP stock", ent["ATP"].size, ATP_STOCK)
    want("ATP processes", ent["ATP"].regenerates_from, ["Glycolysis", "TCA", "OXPHOS"])
    want("copies", ent["Mitochondrion"].copies, MITO_COPIES)
    want("heteroplasmy", ent["Mitochondrion"].heteroplasmy, {"m.3243A>G": DEMO_HETEROPLASMY})
    want("capacity per copy", cell.parameters["oxphos_capacity_per_copy"].value, CAPACITY_PER_COPY)
    alloc = [e for e in ent.values() if isinstance(e, Allocation)]
    want(
        "allocation",
        [(a.pool, a.policy, tuple(a.order)) for a in alloc],
        [("ATP", "priority", PRIORITY_ORDER)],
    )
    want("red cell glucose uptake", rbc.entities["Glucose"].regeneration, RBC_GLUCOSE_UPTAKE_PER_H)
    want("red cell processes", rbc.entities["ATP"].regenerates_from, ["Glycolysis"])
    want("red cell shunt", rbc.parameters["bpg_shunt"].value, BPG_SHUNT)
    want(
        "red cell mitochondria",
        [c.id for c in rbc.entities.values() if isinstance(c, Compartment) and c.mitochondrial],
        [],
    )
    return out


# ---- the arms ----------------------------------------------------------------------------------------


def at(
    module: Any, o2: float = 1.0, glucose: float = 1.0, glucose_ceiling: float = GLUCOSE_UPTAKE_PER_H
) -> Metabolism:
    return Metabolism.from_module(
        module, supply={"O2": O2_SUPPLY_PER_H * o2, "Glucose": glucose_ceiling * glucose}
    )


def with_heteroplasmy(module: Any, h: float) -> Any:
    m = copy.deepcopy(module)
    m.entities["Mitochondrion"].heteroplasmy = {"m.3243A>G": h} if h else {}
    return m


def with_policy(module: Any, policy: str) -> Any:
    m = copy.deepcopy(module)
    for e in m.entities.values():
        if isinstance(e, Allocation) and e.pool == "ATP":
            e.policy, e.order = policy, []
    return m


def with_param(module: Any, name: str, value: float) -> Any:
    m = copy.deepcopy(module)
    m.parameters[name].value = value
    return m


def clause_yields(cell: Any) -> dict[str, Any]:
    m = Metabolism.from_module(cell)
    k = m.constants
    got: dict[str, Any] = {}
    ok = True
    for arm in ARMS:
        y = m.yields(arm)
        full = y["glycolytic"] + y["tca"] + y["oxphos"]
        got[arm] = {"aerobic": full, "anaerobic": y["glycolytic"], "o2": y["o2"], "pairs": y["pairs"]}
        ok &= full == PREDICTED_YIELD_AEROBIC[arm] and y["glycolytic"] == PREDICTED_YIELD_ANAEROBIC
        ok &= y["o2"] == PREDICTED_O2_PER_GLUCOSE and y["pairs"] == PREDICTED_PAIRS_PER_GLUCOSE
    # §6's "30" by hand from the module's constants: Hinkle's 2.5 and 1.5, glycerol 3-phosphate shuttle
    spec_low = (
        k["glycolysis_atp"]
        + k["glycolysis_pyruvate"] * k["tca_atp"]
        + k["glycolysis_pyruvate"] * (k["tca_nadh"] * k["po_nadh"] + k["tca_fadh2"] * k["po_fadh2"])
        + k["glycolysis_nadh"] * k["po_fadh2"]
    )
    got["p_o_2.5_1.5_glycerol_phosphate"] = spec_low
    ok &= spec_low == PREDICTED_YIELD_SPEC_LOW_END
    return {"ok": bool(ok), "measured": got}


def clause_budget(cell: Any) -> dict[str, Any]:
    out, ok = {}, True
    for arm in ARMS:
        f = at(cell).solve(arm)
        y = f.yield_aerobic
        out[arm] = {
            "supply": f.supply,
            "glucose": f.glucose,
            "o2": f.o2,
            "lactate": f.lactate,
            "oxphos_share": f.oxphos_share,
            "limit": f.limit,
        }
        ok &= f.supply == DEMAND_PER_H and close(f.glucose, DEMAND_PER_H / y) and close(f.o2, 6.0 * f.glucose)
        ok &= f.lactate == 0.0 and close(f.oxphos_share, (y - 4.0) / y)
    return {"ok": bool(ok), "measured": out}


def sweep(module: Any, arm: str) -> list[Any]:
    return [at(module, o2, g).solve(arm) for o2, g in SUPPLY_SWEEP]


def clause_pasteur(cell: Any) -> dict[str, Any]:
    out, ok = {}, True
    oxygen_steps = [i for i, (_, g) in enumerate(SUPPLY_SWEEP) if g == 1.0]
    for arm in ARMS:
        runs = sweep(cell, arm)
        share = [r.oxphos_share for r in runs]
        glucose = [r.glucose for r in runs]
        out[arm] = {
            "supply": [r.supply for r in runs],
            "glucose": glucose,
            "oxphos_share": share,
            "lactate": [r.lactate for r in runs],
        }
        ox_share = [share[i] for i in oxygen_steps]
        ox_glucose = [glucose[i] for i in oxygen_steps]
        ok &= all(b <= a for a, b in zip(ox_share, ox_share[1:], strict=False)) and ox_share[-1] == 0.0
        ok &= all(b >= a for a, b in zip(ox_glucose, ox_glucose[1:], strict=False))
        for i, r in enumerate(runs):
            ok &= close(r.supply, PREDICTED_SUPPLY[arm][i]) and close(r.glucose, PREDICTED_GLUCOSE[arm][i])
            ok &= close(r.oxphos_share, PREDICTED_OXPHOS_SHARE[arm][i])
    return {"ok": bool(ok), "measured": out}


def clause_pasteur_uncapped(cell: Any) -> dict[str, Any]:
    out, ok = {}, True
    for arm in ARMS:
        aerobic = at(cell, 1.0, 1.0, UNCAPPED_GLUCOSE_PER_H).solve(arm)
        anoxic = at(cell, 0.0, 1.0, UNCAPPED_GLUCOSE_PER_H).solve(arm)
        ratio = anoxic.glucose / aerobic.glucose
        out[arm] = {
            "aerobic_glucose": aerobic.glucose,
            "anoxic_glucose": anoxic.glucose,
            "ratio": ratio,
            "anoxic_supply": anoxic.supply,
        }
        ok &= aerobic.supply == anoxic.supply == DEMAND_PER_H and close(ratio, PREDICTED_PASTEUR_RATIO[arm])
    return {"ok": bool(ok), "measured": out}


def clause_glycolysis_only(rbc: Any) -> dict[str, Any]:
    out, ok = {}, True
    for arm in ARMS:
        base = at(rbc, 1.0, 1.0, RBC_GLUCOSE_UPTAKE_PER_H).solve(arm)
        per_o2 = [
            at(rbc, o2, 1.0, RBC_GLUCOSE_UPTAKE_PER_H).solve(arm).readings()
            for o2, g in SUPPLY_SWEEP
            if g == 1.0
        ]
        noshunt = at(with_param(rbc, "bpg_shunt", 0.0), 1.0, 1.0, RBC_GLUCOSE_UPTAKE_PER_H).solve(arm)
        out[arm] = {
            "supply": base.supply,
            "yield": base.yield_anaerobic,
            "glucose": base.glucose,
            "lactate": base.lactate,
            "shunt_flux": base.shunt_flux,
            "atp_by_process": base.atp_by_process,
            "identical_at_every_oxygen_level": all(r == per_o2[0] for r in per_o2),
            "glucose_without_shunt": noshunt.glucose,
            "ratio_with_to_without_shunt": noshunt.glucose / base.glucose,
        }
        ok &= base.atp_by_process["TCA"] == 0.0 and base.atp_by_process["OXPHOS"] == 0.0
        ok &= close(base.yield_anaerobic, PREDICTED_RBC_YIELD) and close(base.glucose, PREDICTED_RBC_GLUCOSE)
        ok &= close(base.lactate, PREDICTED_RBC_LACTATE) and close(base.shunt_flux, PREDICTED_RBC_SHUNT_FLUX)
        ok &= base.supply == RBC_DEMAND_PER_H and out[arm]["identical_at_every_oxygen_level"]
        ok &= close(out[arm]["ratio_with_to_without_shunt"], PREDICTED_RBC_YIELD / 2.0)
    return {"ok": bool(ok), "measured": out}


def clause_glucose_removed(cell: Any, rbc: Any) -> dict[str, Any]:
    out, ok = {}, True
    for arm in ARMS:
        r = at(rbc, 1.0, 0.0, RBC_GLUCOSE_UPTAKE_PER_H).solve(arm)
        c = at(cell, 1.0, 0.0).solve(arm)
        out[arm] = {"red_cell_supply": r.supply, "red_cell_factors": r.factors, "nucleated_supply": c.supply}
        ok &= (
            r.supply == 0.0
            and all(v == 0.0 for v in r.factors.values())
            and len(r.factors) == len(RBC_CONSUMERS)
        )
    return {
        "ok": bool(ok),
        "measured": out,
        "reported_not_scored": "the nucleated cell supplies 0 too: no other fuel",
    }


def first_falls(runs: list[Any]) -> dict[str, int | None]:
    return {c: next((i for i, r in enumerate(runs) if r.factors[c] < 1.0), None) for c in CONSUMERS}


def clause_hierarchy(cell: Any) -> dict[str, Any]:
    out: dict[str, Any] = {}
    ok = True
    anoxic = SUPPLY_SWEEP.index((0.0, 1.0))
    for arm in ARMS:
        runs = sweep(cell, arm)
        falls = first_falls(runs)
        out[arm] = {"priority": {"factors": [r.factors for r in runs], "first_fall_step": falls}}
        macro = [falls[c] for c in MACROMOLECULE]
        ordered = (
            None not in (*macro, falls["na_k_pumping"], falls["ca_pumping"])
            and max(macro) < falls["na_k_pumping"] < falls["ca_pumping"]
        )
        ok &= bool(ordered) and runs[anoxic].factors == PREDICTED_ANOXIC_FACTORS
        for policy in ("proportional", "competitive"):
            scalar = sweep(with_policy(cell, policy), arm)
            equal = all(len(set(r.factors.values())) == 1 for r in scalar)
            out[arm][policy] = {
                "factors": [r.factors for r in scalar],
                "identical_at_every_step": equal,
                "first_fall_step": first_falls(scalar),
            }
            ok &= equal
    return {
        "ok": bool(ok),
        "measured": out,
        "reading": (
            "reproduced as an input, not predicted: the order appears under the priority the program "
            "writes, and under proportional and competitive every consumer falls together"
        ),
    }


def clause_threshold(cell: Any) -> dict[str, Any]:
    out: dict[str, Any] = {}
    ok = True
    for arm in ARMS:
        runs = {h: at(with_heteroplasmy(cell, h)).solve(arm) for h in HETEROPLASMY_SWEEP}
        base = runs[0.0]
        star = PREDICTED_THRESHOLD[arm]
        rows = {
            h: {
                "oxphos_atp": r.atp_by_process["OXPHOS"],
                "glucose": r.glucose,
                "supply": r.supply,
                "limit": r.limit,
            }
            for h, r in runs.items()
        }
        below = all(
            r.atp_by_process["OXPHOS"] == base.atp_by_process["OXPHOS"] for h, r in runs.items() if h < star
        )
        above = all(
            r.atp_by_process["OXPHOS"] < base.atp_by_process["OXPHOS"] and r.glucose > base.glucose
            for h, r in runs.items()
            if h > star
        )
        last_flat = max(h for h in HETEROPLASMY_SWEEP if h < star)
        first_drop = min(h for h in HETEROPLASMY_SWEEP if h > star)
        full_load = close(runs[1.0].supply, PREDICTED_SUPPLY[arm][SUPPLY_SWEEP.index((0.0, 1.0))])
        # falsifier: no reserve, capacity exactly what the demand needs
        y = base.yield_aerobic
        tight = with_param(cell, "oxphos_capacity_per_copy", 12.0 * DEMAND_PER_H / (y * MITO_COPIES))
        t_runs = {h: at(with_heteroplasmy(tight, h)).solve(arm) for h in HETEROPLASMY_SWEEP}
        t_base = t_runs[0.0].atp_by_process["OXPHOS"]
        no_plateau = all(r.atp_by_process["OXPHOS"] < t_base for h, r in t_runs.items() if h > 0)
        out[arm] = {
            "sweep": rows,
            "registered_threshold": star,
            "last_flat_point": last_flat,
            "first_lower_point": first_drop,
            "flat_below": below,
            "lower_above": above,
            "full_load_equals_anoxia": full_load,
            "no_reserve": {
                "oxphos_atp": {h: r.atp_by_process["OXPHOS"] for h, r in t_runs.items()},
                "every_h_above_0_lower": no_plateau,
            },
        }
        ok &= below and above and full_load and no_plateau
    return {"ok": bool(ok), "measured": out}


def clause_checkpoint(cell: Any) -> dict[str, Any]:
    division = Division.from_module(parse_file(ROOT / STAGE4_PROGRAM))
    cost = [(c.pool, c.amount, c.per) for c in division.event.costs]
    pool: Pool = cell.entities["ATP"]

    def economy(capacity: float | None) -> Economy:
        cap = {} if capacity is None else {"ATP": capacity}
        return Economy(
            pools={"ATP": pool}, costs={division.event.id: list(division.event.costs)}, capacity=cap
        )

    met = {"size": "2"}  # stage 4's own checkpoint, `when: size >= 2`, met
    stock_ready, stock_why = division.ready(met, economy(None))
    out: dict[str, Any] = {"division_cost": cost, "standing_stock": {"ready": stock_ready, "why": stock_why}}
    ok = cost == [("ATP", DIVISION_COST_ATP, "division")] and not stock_ready
    for arm in ARMS:
        m = at(cell)
        available = m.available(INTERVAL_H, arm)
        ready, why = division.ready(met, economy(available))
        starved = at(cell, 1.0, 0.0).available(INTERVAL_H, arm)
        s_ready, s_why = division.ready(met, economy(starved))
        out[arm] = {
            "max_supply": m.solve(arm).max_supply,
            "available_over_interval": available,
            "ready": ready,
            "available_without_glucose": starved,
            "ready_without_glucose": s_ready,
            "why_without_glucose": s_why,
        }
        ok &= close(m.solve(arm).max_supply, PREDICTED_MAX_SUPPLY[arm]) and close(
            available, PREDICTED_AVAILABLE[arm]
        )
        ok &= ready and not s_ready and starved == ATP_STOCK
    return {"ok": bool(ok), "measured": out}


REFUSED_AT_COMPILE = {
    "heteroplasmy without chrM": "compartment C { heteroplasmy: m.3243A>G = 0.7 }",
    "fraction above 1": "compartment M { genome: chrM; heteroplasmy: m.3243A>G = 1.2 }",
    "fractions summing past 1": (
        "compartment M { genome: chrM; heteroplasmy: m.3243A>G = 0.7, m.8344A>G = 0.4 }"
    ),
    "not a substitution": "compartment M { genome: chrM; heteroplasmy: m.3243delA = 0.2 }",
    "no m. prefix": "compartment M { genome: chrM; heteroplasmy: 3243A>G = 0.2 }",
}
CELL = """module t
import bio.std.core_metabolism
compartment Cytosol { volume: 0.7 }
compartment Mito { parent: Cytosol; genome: chrM; copies: 10 }
pool ATP { location: Cytosol; size: 1e6; regenerates: from Glycolysis, TCA, OXPHOS }
pool Glucose { location: Cytosol; size: 1; regenerates: 100 /h }
pool O2 { location: Cytosol; size: 1; regenerates: 600 /h }
param oxphos_capacity_per_copy = 100 molecules/h { confidence: 0.1 }
event work { rate: 320 /h; cost: ATP 1 per event }
"""
PROCESSES = "from Glycolysis, TCA, OXPHOS"
REFUSED_AT_RUN = {
    "unknown process": CELL.replace(PROCESSES, "from Glycolysis, Photosynthesis"),
    "no Glycolysis": CELL.replace(PROCESSES, "from OXPHOS"),
    "TCA without OXPHOS": CELL.replace(PROCESSES, "from Glycolysis, TCA"),
    "OXPHOS without a mitochondrion": CELL.replace(
        "compartment Mito { parent: Cytosol; genome: chrM; copies: 10 }\n", ""
    ),
    "no Glucose pool": CELL.replace("pool Glucose { location: Cytosol; size: 1; regenerates: 100 /h }\n", ""),
    "no O2 pool": CELL.replace("pool O2 { location: Cytosol; size: 1; regenerates: 600 /h }\n", ""),
    "no capacity": CELL.replace("param oxphos_capacity_per_copy = 100 molecules/h { confidence: 0.1 }\n", ""),
    "no stoichiometry": CELL.replace("import bio.std.core_metabolism\n", ""),
}


def clause_refusals() -> dict[str, Any]:
    out: dict[str, Any] = {"compile": {}, "run": {}}
    for name, block in REFUSED_AT_COMPILE.items():
        try:
            parse("module t\n" + block)
            out["compile"][name] = "ACCEPTED"
        except BioLangError as e:
            out["compile"][name] = f"refused: {str(e)[:100]}"
    for name, text in REFUSED_AT_RUN.items():
        assert text != CELL, name
        try:
            Metabolism.from_module(parse(text))
            out["run"][name] = "ACCEPTED"
        except MetabolismError as e:
            out["run"][name] = f"refused: {str(e)[:100]}"
    control = Metabolism.from_module(parse(CELL)).solve("high").supply  # the unbroken program runs
    out["control_runs"] = control == 320.0
    ok = out["control_runs"] and all(
        v.startswith("refused") for part in ("compile", "run") for v in out[part].values()
    )
    return {"ok": bool(ok), "measured": out}


def manifest() -> dict[str, Any]:
    return {
        "sources": [
            {
                "accession": (
                    "PubMed abstracts PMID 7492307 (Buttgereit & Brand 1995), 15620362 (Hinkle "
                    "2005), 12467494 (Rossignol et al. 2003), 10477268 (Mulquiney, Bubb & Kuchel 1999)"
                ),
                "version": (
                    "read 2026-10-10 and quoted in genomeos/runtime/metabolism_gate.py; no dataset is read"
                ),
            }
        ],
        "inputs": [
            mf.input_entry(PROGRAM, partition=None, role="the nucleated cell, every number declared"),
            mf.input_entry(RBC_PROGRAM, partition=None, role="the glycolysis-only cell"),
            mf.input_entry(STD, partition=None, role="the stoichiometry"),
            mf.input_entry(STAGE4_PROGRAM, partition=None, role="stage 4's division cost, read only"),
            mf.input_entry("genomeos/runtime/metabolism_gate.py", partition=None, role="the registration"),
            mf.input_entry("genomeos/runtime/metabolism.py", partition=None, role="the instrument"),
        ],
        "assembly": "n/a: no genomic coordinate is read; m.3243A>G is a declared label, not a lookup",
        "coordinates": "n/a: no genomic coordinate is read",
        "parameters": {
            "demand_per_h": DEMAND_PER_H,
            "glucose_uptake_per_h": GLUCOSE_UPTAKE_PER_H,
            "o2_supply_per_h": O2_SUPPLY_PER_H,
            "mito_copies": MITO_COPIES,
            "capacity_per_copy": CAPACITY_PER_COPY,
            "supply_sweep": [list(s) for s in SUPPLY_SWEEP],
            "heteroplasmy_sweep": list(HETEROPLASMY_SWEEP),
            "bpg_shunt": BPG_SHUNT,
            "interval_h": INTERVAL_H,
            "closed_form_rtol": CLOSED_FORM_RTOL,
        },
        "exclusions": [f"{k}: {v}" for k, v in PRE_REGISTRATION["not_askable"].items()],
        "partitions": "n/a: nothing is fitted, so there is nothing to hold out",
        "code_cleanliness": mf.code_cleanliness(ENTRY, OWN_CODE, ROOT),
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--no-save", action="store_true")
    args = ap.parse_args(argv)
    t0 = time.time()
    cell, rbc = parse_file(ROOT / PROGRAM), parse_file(ROOT / RBC_PROGRAM)
    problems = disagreements(cell, rbc)
    if problems:
        print("REFUSED: a demo program disagrees with its registration:\n  " + "\n  ".join(problems))
        return 2
    clauses = {
        "1_yields": clause_yields(cell),
        "2_budget": clause_budget(cell),
        "3_pasteur": clause_pasteur(cell),
        "4_pasteur_uncapped": clause_pasteur_uncapped(cell),
        "5_glycolysis_only": clause_glycolysis_only(rbc),
        "6_glucose_removed": clause_glucose_removed(cell, rbc),
        "7_hierarchy": clause_hierarchy(cell),
        "8_threshold": clause_threshold(cell),
        "9_resource_checkpoint": clause_checkpoint(cell),
        "10_refusals": clause_refusals(),
    }
    assert set(clauses) == set(PRE_REGISTRATION["clauses"]), "the clauses run are the clauses registered"
    met = [k for k, v in clauses.items() if v["ok"]]
    passed = len(met) == len(clauses)
    verdict = {
        "clauses_met": f"{len(met)} of {len(clauses)}",
        "passed": passed,
        "gate_as_worded": (
            "(i) the budget's magnitude and its split are NOT ASKABLE (no measured turnover or shares held; "
            "§10 decision 4 stays open), and Buttgereit & Brand's hierarchy is REPRODUCED AS AN INPUT, NOT "
            "PREDICTED; (ii) lowering oxygen shifts production to glycolysis and raises glucose uptake, BY "
            "CONSTRUCTION of the routing, and removing glucose collapses the glycolysis-only cell (and every "
            "cell, there being no other fuel); (iii) the red cell runs on glycolysis alone, and its ATP and "
            "2,3-BPG ranges are NOT ASKABLE (not held; a concentration is not a stoichiometric quantity)"
        ),
    }
    for k, v in clauses.items():
        print(f"{'ok  ' if v['ok'] else 'FAIL'} {k}")
    print(f"\n{verdict['clauses_met']} registered clauses met\n{verdict['gate_as_worded']}")
    out = {
        "result": RESULT,
        "result_manifest": manifest(),
        "pre_registration": PRE_REGISTRATION,
        "programs": [PROGRAM, RBC_PROGRAM],
        "clauses": clauses,
        **verdict,
        "seconds": round(time.time() - t0, 2),
    }
    if args.no_save:
        print("\nnot saved (--no-save)")
    else:
        from genomeos.results import save_result

        print(f"\nsaved {save_result(RESULT, out)}")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
