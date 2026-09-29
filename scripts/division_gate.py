# SPDX-License-Identifier: AGPL-3.0-or-later
"""Run stage 4's gate against its registration (runtime/division_gate.py). No network, no quota.

    uv run python scripts/division_gate.py

The design is `PRE_REGISTRATION`, committed in `genomeos/runtime/division_gate.py` before this file
existed. Nothing here may change it. In particular the registration says in advance that the gate as
§8 words it — *"dilution must separate a stable protein from a short-lived one"* — cannot be met,
because partitioning multiplies every species by the same one-half and the division term cancels out
of the ratio. This script's job is to measure that rather than to argue it: the duplicating arm is
the falsifier, and a run in which the two arms' separation differed would mean the partitioning is
indexed by protein, which is a bug of the shape that looks like a finding.
"""

from __future__ import annotations

import argparse
import math
import random
import statistics
import time
from dataclasses import replace
from pathlib import Path
from typing import Any

from genomeos.ir import UNKNOWN
from genomeos.lang import parse_file
from genomeos.results import save_result
from genomeos.runtime.division import Division
from genomeos.runtime.division_gate import (
    BINOMIAL_MEAN_TOL,
    BINOMIAL_VAR_RTOL,
    CLOSED_FORM_RTOL,
    DIVISIONS,
    INTERVAL_H,
    LOW_COPIES,
    PRE_REGISTRATION,
    PREDICTED_COMPRESSION,
    PREDICTED_DILUTION_ONLY,
    PREDICTED_SEPARATION_DECADES,
    PREDICTED_STEADY_RATIO_DIVIDING,
    PREDICTED_STEADY_RATIO_STATIC,
    SEEDS,
    SHORT_HALF_LIFE_H,
    STABLE_HALF_LIFE_H,
    START_AMOUNT,
)
from genomeos.runtime.economy import Economy

PROGRAM = Path("data/demo/stage4_division.bio")
STABLE, SHORT, UNDEGRADED = "STABLEp", "SHORTp", "UNDEGRADEDp"
MET = {"size": "2"}  # the size checkpoint as the program writes it: `when: size >= 2`
UNMET = {"size": "1"}
START = dict.fromkeys((STABLE, SHORT, UNDEGRADED), START_AMOUNT)


def _close(got: float, want: float, rtol: float = CLOSED_FORM_RTOL) -> bool:
    return abs(got - want) <= rtol * max(abs(want), 1e-300)


def arm(division: Division, context: dict[str, str], half_lives: dict[str, Any] | None = None) -> dict:
    """One chase: `DIVISIONS` cycles of decay then division, from `START`, following one daughter."""
    d = division if half_lives is None else replace(division, half_lives=half_lives)
    got = d.run(START, DIVISIONS, INTERVAL_H, context)
    got["closed_form"] = {
        s: d.expected(START_AMOUNT, d.half_lives[s], got["divisions_fired"], INTERVAL_H, DIVISIONS)
        for s in START
    }
    got["matches_closed_form"] = all(_close(got["final"][s], got["closed_form"][s]) for s in START)
    got["separation_decades"] = math.log10(got["final"][STABLE]) - math.log10(got["final"][SHORT])
    return got


def binomial_arm(division: Division, copies: float) -> dict:
    """`SEEDS` independent single divisions of one species, and what they did to it."""
    drawn, conserved, whole = [], 0, 0
    for seed in range(SEEDS):
        got = division.split({"X": copies}, random.Random(seed))
        drawn.append(got.a["X"])
        conserved += got.conserves({"X": copies})
        whole += float(got.a["X"]).is_integer()
    return {
        "copies": copies,
        "seeds": SEEDS,
        "mean": statistics.fmean(drawn),
        "variance": statistics.pvariance(drawn),
        "analytic_mean": copies / 2.0,
        "analytic_variance": copies / 4.0,
        "conserved": conserved,
        "whole_daughters": whole,
        "distinct_values": len(set(drawn)),
    }


def resource_arm(division: Division, module: Any) -> dict:
    """The other checkpoint §8 names: the division's own cost, against stage 2's pools.

    **Not part of the verdict, and not because it is inconvenient.** The registration's clause 8 is
    the `when` arm alone; this arm was written while building and is reported as what it is, an
    observation the registration did not anticipate. It is the more interesting half of the run.
    """
    economy = Economy.from_module(module)
    declared = division.ready(MET, economy)
    rich_pools = {p: replace(pool, size=float(pool.size) * 1e3) for p, pool in economy.pools.items()}
    rich = division.ready(MET, replace(economy, pools=rich_pools))
    poor_pools = {p: replace(pool, size=1.0) for p, pool in economy.pools.items()}
    poor = division.ready(MET, replace(economy, pools=poor_pools))
    return {
        "contributes_to_verdict": False,
        "as_declared": {"ready": declared[0], "why": declared[1]},
        "pools_x1000": {"ready": rich[0], "why": rich[1]},
        "unaffordable": {"ready": poor[0], "why": poor[1]},
        "discriminates": rich[0] and not poor[0],
        "finding": (
            "the program declares §5.2's own example, `event divide { cost: ATP 1e10 }`, against "
            "§5.1's own ATP pool of 3e9, so the resource checkpoint refuses every division: one "
            "division costs more ATP than the cell holds. That is not a bug in either number. A "
            "pool is a standing stock and a division cost is a draw over a whole cycle, and nothing "
            "in the engine converts between them: `regenerates` takes a rate, this pool declares "
            "none (`from Glycolysis, OXPHOS`), and `Economy` never integrates a rate over an "
            "interval. So the resource form of §8's checkpoint is only askable against a standing "
            "stock today, and is reported rather than scored. What would make it askable is stage "
            "3, which is what supplies ATP over time"
        ),
    }


def verdict(runs: dict[str, Any]) -> dict[str, Any]:
    """The registered bars, applied without interpretation."""
    part, dup = runs["partitioned"], runs["duplicating"]
    equal, unmet = runs["equal_half_lives"], runs["checkpoint_unmet"]
    low, high = runs["binomial_low"], runs["binomial_high"]
    steady = runs["steady_state"]

    c1 = part["matches_closed_form"] and dup["matches_closed_form"]
    c2 = (
        part["final"][UNDEGRADED] == START_AMOUNT * PREDICTED_DILUTION_ONLY
        and dup["final"][UNDEGRADED] == START_AMOUNT
    )
    c3 = _close(part["separation_decades"], dup["separation_decades"]) and _close(
        part["separation_decades"], PREDICTED_SEPARATION_DECADES
    )
    c4 = _close(steady["ratio_dividing"], PREDICTED_STEADY_RATIO_DIVIDING) and _close(
        steady["ratio_static"], PREDICTED_STEADY_RATIO_STATIC
    )
    c5 = part["conserved"] and low["conserved"] == SEEDS and high["conserved"] == SEEDS
    c6 = (
        abs(low["mean"] - low["analytic_mean"]) <= BINOMIAL_MEAN_TOL
        and abs(low["variance"] / low["analytic_variance"] - 1.0) <= BINOMIAL_VAR_RTOL
        and low["whole_daughters"] == SEEDS
        and high["variance"] == 0.0
    )
    c7 = equal["final"][STABLE] == equal["final"][SHORT]
    # clause 8 as registered is the `when` arm alone. The resource arm below it is reported and
    # deliberately excluded, the way gate (a) excluded its unasked magnitude clause (§9.2).
    c8 = unmet["divisions_fired"] == 0 and all(
        _close(unmet["final"][s], unmet["closed_form"][s]) for s in START
    )

    clauses = {
        "1_closed_form": c1,
        "2_dilution_only": c2,
        "3_dilution_separates_nothing": c3,
        "4_compression": c4,
        "5_conservation": c5,
        "6_binomial": c6,
        "7_equal_half_lives": c7,
        "8_checkpoint": c8,
    }
    passed = all(clauses.values())
    return {
        "clauses": clauses,
        "passed": passed,
        "separation_partitioned_decades": part["separation_decades"],
        "separation_duplicating_decades": dup["separation_decades"],
        "dilution_contributes_decades": part["separation_decades"] - dup["separation_decades"],
        "resource_checkpoint_excluded": runs["resource_checkpoint"]["finding"],
        "measured_turnover": (
            "not askable: this project holds no protein turnover table, so no clause compares these "
            "declared half-lives with a measurement (registered)"
        ),
        "gate_as_worded": (
            "REFUTED, as registered before the run: dilution contributes "
            f"{part['separation_decades'] - dup['separation_decades']:.1e} decades of the "
            f"{part['separation_decades']:.6f} that separate a stable protein from a short-lived "
            f"one, and in a steady state division compresses the separation "
            f"{PREDICTED_COMPRESSION:.2f}-fold. What separates them is turnover, which separates "
            "them just as well in a cell that never divides"
        ),
        "verdict": (
            "stage 4's mechanism passes all eight registered clauses. §8's gate as worded is not "
            "passed by them and cannot be: it is answered with a negative, derived before the run"
            if passed
            else "stage 4 FAILS: " + ", ".join(n for n, ok in clauses.items() if not ok)
        ),
    }


def manifest() -> dict[str, Any]:
    """What this result was made from: the program, the registration and the two engine modules.

    No public dataset appears here and that is the finding rather than an omission — the half-lives
    are declared by the program, not read from a turnover table, because this project holds none.
    """
    from genomeos import manifest as mf

    return {
        "sources": [
            {
                "accession": "no external dataset",
                "version": "n/a: the half-lives are declared by "
                "data/demo/stage4_division.bio, not measured; §10 decision 6 is open",
            },
        ],
        "inputs": [
            mf.input_entry(str(PROGRAM)),
            mf.input_entry("genomeos/runtime/division.py"),
            mf.input_entry("genomeos/runtime/division_gate.py"),
        ],
        "assembly": "n/a: no genomic coordinate is read",
        "coordinates": "n/a: no genomic coordinate is read",
        "parameters": {
            "divisions": DIVISIONS,
            "interval_h": INTERVAL_H,
            "stable_half_life_h": STABLE_HALF_LIFE_H,
            "short_half_life_h": SHORT_HALF_LIFE_H,
            "start_amount": START_AMOUNT,
            "low_copies": LOW_COPIES,
            "seeds": SEEDS,
            "closed_form_rtol": CLOSED_FORM_RTOL,
            "binomial_mean_tol": BINOMIAL_MEAN_TOL,
            "binomial_var_rtol": BINOMIAL_VAR_RTOL,
        },
        "exclusions": [
            "the resource form of the checkpoint: measured and reported, excluded from the verdict "
            "because the registration's clause 8 is the `when` form alone",
            "any comparison with a measured protein turnover set: registered as not askable",
        ],
        "partitions": "n/a: nothing is fitted, so there is nothing to hold out",
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--no-save", action="store_true")
    args = ap.parse_args(argv)

    t0 = time.time()
    module = parse_file(str(PROGRAM))
    partitioning = Division.from_module(module)
    duplicating = replace(partitioning, event=replace(partitioning.event, partition="duplicate"))
    equal = {STABLE: STABLE_HALF_LIFE_H, SHORT: STABLE_HALF_LIFE_H, UNDEGRADED: UNKNOWN}

    runs = {
        "partitioned": arm(partitioning, MET),
        "duplicating": arm(duplicating, MET),
        "equal_half_lives": arm(partitioning, MET, equal),
        "checkpoint_unmet": arm(partitioning, UNMET),
        "binomial_low": binomial_arm(partitioning, float(LOW_COPIES)),
        "binomial_high": binomial_arm(partitioning, START_AMOUNT),
        "resource_checkpoint": resource_arm(partitioning, module),
        "steady_state": {
            "ratio_static": partitioning.steady_state(1.0, STABLE_HALF_LIFE_H, None)
            / partitioning.steady_state(1.0, SHORT_HALF_LIFE_H, None),
            "ratio_dividing": partitioning.steady_state(1.0, STABLE_HALF_LIFE_H, INTERVAL_H)
            / partitioning.steady_state(1.0, SHORT_HALF_LIFE_H, INTERVAL_H),
        },
    }
    runs["steady_state"]["compression"] = (
        runs["steady_state"]["ratio_static"] / runs["steady_state"]["ratio_dividing"]
    )
    got = verdict(runs)

    print(f"after {DIVISIONS} divisions {INTERVAL_H:g} h apart, from {START_AMOUNT:.0e} of each:")
    print(f"{'protein':>14} {'partitioned':>14} {'duplicating':>14}")
    for s in (STABLE, SHORT, UNDEGRADED):
        print(f"{s:>14} {runs['partitioned']['final'][s]:14.6g} {runs['duplicating']['final'][s]:14.6g}")
    print(
        f"\nseparation stable:short  partitioned {runs['partitioned']['separation_decades']:.6f} "
        f"decades, duplicating {runs['duplicating']['separation_decades']:.6f} decades"
    )
    print(
        f"steady state             {runs['steady_state']['ratio_static']:.4f} without division, "
        f"{runs['steady_state']['ratio_dividing']:.4f} with it "
        f"(compression {runs['steady_state']['compression']:.4f})"
    )
    print(
        f"binomial, {LOW_COPIES} copies    mean {runs['binomial_low']['mean']:.4f} "
        f"(analytic {runs['binomial_low']['analytic_mean']}), variance "
        f"{runs['binomial_low']['variance']:.4f} (analytic {runs['binomial_low']['analytic_variance']}), "
        f"{runs['binomial_low']['conserved']} of {SEEDS} conserved"
    )
    print(f"\nresource checkpoint (reported, not scored): {runs['resource_checkpoint']['finding']}")
    print(f"\n{got['gate_as_worded']}\n{got['verdict']}")

    out = {
        "result": "division_gate",
        "result_manifest": manifest(),
        "pre_registration": PRE_REGISTRATION,
        "program": str(PROGRAM),
        "runs": runs,
        **got,
        "seconds": round(time.time() - t0, 2),
    }
    if args.no_save:
        print("\nnot saved (--no-save)")
    else:
        print(f"\nsaved {save_result(out['result'], out)}")
    return 0 if got["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
