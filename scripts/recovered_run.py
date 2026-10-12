# SPDX-License-Identifier: AGPL-3.0-or-later
"""lane-recovered step 3: THE ONE RUN. Read the registered endpoint and report it, whatever it says.

    uv run python scripts/recovered_run.py [--no-save]

The registration is `data/results/recovered_registration.json` (`266f031`), written and committed before
any P, D or Delta MSE of this population was read. This script re-checks every BOUND field against
`genomeos/attribution/recovered.py`, runs the HARNESS, and then reads ONE number: the held-out Delta MSE
of the multiplicative composition at M2h's fixed constant +0.2146 over the genes lane-scalefit's scaling
makes fittable. `magcal`'s `bootstrap_delta`, `permuted_control`, `decide`, `reading` and
`combined_reading` are called UNCHANGED and nothing in this lane is fitted: magcal's scale-fitting
and scale-bootstrapping functions are never called, and a test holds their names out of this file.

ONE RUN. The reading is whatever the registered words say for the outcome, including a WORSE or an
INCONCLUSIVE one, which DECISION step 6 requires to be reported first and unsoftened.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import magcal as mc  # noqa: E402
from genomeos.attribution import recovered as rc  # noqa: E402
from genomeos.results import save_result  # noqa: E402

RESULT = "recovered_result"
REGISTRATION = Path("data/results/recovered_registration.json")
POPULATION = Path("data/results/recovered_population.json")
MULTRNA = Path("data/results/multrna_result.json")
MAGCAL = Path("data/results/magcal_result.json")
FRESHRNA = Path("data/results/freshrna_result.json")
SCALEFIT_CENSUS = Path("data/results/scalefit_census.json")
OWN_CODE = (
    "genomeos/attribution/recovered.py",
    "scripts/recovered_population.py",
    "scripts/recovered_register.py",
    "scripts/recovered_run.py",
    "tests/test_recovered.py",
)


# ---- the rows -------------------------------------------------------------------------------------------


def rows_of_7h() -> list[dict[str, Any]]:
    """7h's committed multiplicative rows, in the committed order, as magcal's functions expect them."""
    table = json.loads(MULTRNA.read_text())["genes_table"]
    i_s = table["columns"].index("symbol")
    i_p = table["columns"].index("P")
    i_d = table["columns"].index("D")
    return [{"symbol": r[i_s], "P": r[i_p], "D": r[i_d]} for r in table["rows"]]


def pairs(rows: list[dict]) -> int:
    return sum(len(r["D"]) for r in rows)


def sums(rows: list[dict]) -> dict[str, float]:
    """mean(D^2), mean(P^2), mean(SCALE^2 P^2) and mean(D*P): DESCRIPTIVES, deciding nothing."""
    dp = pp = dd = 0.0
    n = 0
    for r in rows:
        for p, d in zip(r["P"], r["D"], strict=True):
            dp += p * d
            pp += p * p
            dd += d * d
            n += 1
    return {
        "mean_D2": dd / n,
        "mean_P2": pp / n,
        "mean_scaled_P2": rc.SCALE * rc.SCALE * pp / n,
        "mean_D_times_P": dp / n,
        "pairs": n,
    }


# ---- the harness, run before any value of this population is read --------------------------------------


def harness(rows: list[dict], pop: dict) -> dict[str, Any]:
    """recovered.HARNESS, the seven checks, in order."""
    m2h = json.loads(MAGCAL.read_text())
    mult = m2h["per_form"]["multiplicative"]
    committed = json.loads(MULTRNA.read_text())
    by_symbol = {r["symbol"]: r for r in rows}
    held_pop = pop["symbols"]["held_out"]
    m2h_train = {r["symbol"] for r in rows if not mc.is_held_out(r["symbol"])}
    m2h_held = [r for r in rows if mc.is_held_out(r["symbol"])]

    out: dict[str, Any] = {
        "1_population_record_established": (
            pop["outcome"] == "POPULATION_ESTABLISHED" and pop["population_harness"]["agrees"]
        ),
        "2_population_sha256_recomputes": (
            rc.sha256_of_symbols(pop["symbols"]["all"]) == pop["population_sha256"]
            and rc.sha256_of_symbols(held_pop) == pop["held_out_sha256"]
        ),
        "3_7h_table_shape": {
            "genes": len(rows),
            "pairs": pairs(rows),
            "genes_committed": 16509,
            "pairs_committed": 214617,
            "thirteen_each": all(len(r["P"]) == 13 and len(r["D"]) == 13 for r in rows),
            "no_symbol_twice": len(by_symbol) == len(rows),
        },
        "4_delta_mse_at_scale_1_reproduces_7h": {
            "recomputed": mc.delta_mse(rows, 1.0),
            "committed": committed["multiplicative"]["delta_mse"],
            "tolerance": rc.TOLERANCE,
        },
        "5_scale_and_interval_are_m2h_s": (
            m2h["scales"]["multiplicative"] == rc.SCALE
            and m2h["scale_cis"]["multiplicative"] == list(rc.SCALE_CI)
        ),
        "6_m2h_split_reproduces": {
            "held_out_genes": len(m2h_held),
            "held_out_genes_committed": mult["split"]["held_out_genes"],
            "held_out_pairs": pairs(m2h_held),
            "held_out_pairs_committed": 108264,
            "training_genes": len(rows) - len(m2h_held),
            "training_genes_committed": mult["split"]["training_genes"],
        },
        "7_this_population_s_held_out_half_is_held_out_under_m2h_s_split": (
            all(mc.is_held_out(s) for s in held_pop) and not any(s in m2h_train for s in held_pop)
        ),
    }
    s3 = out["3_7h_table_shape"]
    s4 = out["4_delta_mse_at_scale_1_reproduces_7h"]
    s6 = out["6_m2h_split_reproduces"]
    out["agrees"] = bool(
        out["1_population_record_established"]
        and out["2_population_sha256_recomputes"]
        and s3["genes"] == s3["genes_committed"]
        and s3["pairs"] == s3["pairs_committed"]
        and s3["thirteen_each"]
        and s3["no_symbol_twice"]
        and abs(s4["recomputed"] - s4["committed"]) <= rc.TOLERANCE
        and out["5_scale_and_interval_are_m2h_s"]
        and s6["held_out_genes"] == s6["held_out_genes_committed"]
        and s6["held_out_pairs"] == s6["held_out_pairs_committed"]
        and s6["training_genes"] == s6["training_genes_committed"]
        and out["7_this_population_s_held_out_half_is_held_out_under_m2h_s_split"]
    )
    return out


# ---- the one reading -----------------------------------------------------------------------------------


def side(rows: list[dict], label: str) -> dict[str, Any]:
    """One side's Delta MSE at the fixed scale with its interval, plus the descriptives."""
    d = mc.delta_mse(rows, rc.SCALE)
    lo, hi = mc.bootstrap_delta(rows, rc.SCALE, rc.BOOTSTRAP["resamples"], rc.BOOTSTRAP["seed"])
    s = sums(rows)
    return {
        "side": label,
        "genes": len(rows),
        "pairs": s["pairs"],
        "delta_mse": d,
        "ci": [lo, hi],
        "variance_share": d / s["mean_D2"],
        "delta_mse_at_scale_1": mc.delta_mse(rows, 1.0),
        **{k: v for k, v in s.items() if k != "pairs"},
    }


def per_cell(rows: list[dict]) -> dict[str, float]:
    """DESCRIPTIVES: the held-out Delta MSE per reader cell at the fixed scale. Decides nothing."""
    out: dict[str, float] = {}
    for j, cell in enumerate(rc.CELLS):
        dp = pp = 0.0
        for r in rows:
            p, d = r["P"][j], r["D"][j]
            dp += p * d
            pp += p * p
        n = len(rows)
        out[cell] = (2.0 * rc.SCALE * dp - rc.SCALE * rc.SCALE * pp) / n if n else float("nan")
    return out


def run() -> dict:
    reg = json.loads(REGISTRATION.read_text())
    for field, value in rc.registered_fields().items():
        if reg[field] != value:
            raise SystemExit(f"the run differs from the registration at {field}; nothing is read")

    pop = json.loads(POPULATION.read_text())
    rows = rows_of_7h()
    hc = harness(rows, pop)
    base: dict[str, Any] = {
        "lane": rc.LANE,
        "date": time.strftime("%Y-%m-%d"),
        "registration": str(REGISTRATION),
        "step": "3 of 3: THE ONE RUN",
        "question": rc.QUESTION,
        "form": rc.FORM,
        "form_words": rc.FORM_WORDS,
        "scale": rc.SCALE,
        "scale_rounded": rc.SCALE_ROUNDED,
        "scale_ci_carried_from_m2h": list(rc.SCALE_CI),
        "scale_is_refitted": False,
        "harness_check": hc,
    }
    if not hc["agrees"]:
        return {
            **base,
            "outcome": "HARNESS_DEFECT",
            "reading": mc.reading("HARNESS_DEFECT", what="see harness_check"),
            "rna_read": "none: the harness missed and nothing was read",
            "alphagenome_requests": 0,
            "money": "none: no request is sent and no key is read",
        }

    wanted = set(pop["symbols"]["all"])
    held_symbols = set(pop["symbols"]["held_out"])
    mine = [r for r in rows if r["symbol"] in wanted]
    held = [r for r in mine if r["symbol"] in held_symbols]
    train = [r for r in mine if r["symbol"] not in held_symbols]

    if len(held) < rc.MIN_GENES:
        return {
            **base,
            "outcome": "NOT_TESTABLE",
            "reading": mc.reading(
                "NOT_TESTABLE",
                form=rc.FORM_WORDS,
                g=len(held),
                min=rc.MIN_GENES,
                scale="the constant +0.2146",
            ),
            "rna_read": "the population was too small and no comparison was read",
            "alphagenome_requests": 0,
            "money": "none: no request is sent and no key is read",
        }

    # ---- THE REGISTERED ENDPOINT -----------------------------------------------------------------------
    held_side = side(held, "held out (the registered endpoint)")
    train_side = side(train, "training (descriptive, decides nothing)")
    all_side = side(mine, "the whole population (descriptive, decides nothing)")
    lo, hi = held_side["ci"]
    scale_lo, scale_hi = rc.SCALE_CI
    outcome = mc.decide(lo, hi, scale_lo, scale_hi)
    control = mc.permuted_control(train, held, rc.PERMUTATIONS["count"], rc.PERMUTATIONS["seed"])
    control_p = (1 + sum(1 for c in control if c >= held_side["delta_mse"])) / (1 + len(control))

    scale_clause = (
        f"the scale is the constant {rc.SCALE_ROUNDED:+.4f} taken from ac5580b, M2h's fitted scale for"
        f" this form [98.33% CI {scale_lo:+.4f}, {scale_hi:+.4f}] over its own 8,181 training genes, and"
        f" NOTHING IS FITTED IN THIS LANE"
    )
    reading = rc.full_reading(
        outcome,
        delta=held_side["delta_mse"],
        mean_d2=held_side["mean_D2"],
        held=len(held),
        total=pop["re_derived_counts"]["runnable_after"],
        form=rc.FORM_WORDS,
        d=f"{held_side['delta_mse']:+.5g}",
        ci=rc.INTERVAL_WORDS,
        lo=f"{lo:+.4g}",
        hi=f"{hi:+.4g}",
        n=held_side["pairs"],
        g=len(held),
        p=round(control_p, 6),
        scale=scale_clause,
    )

    m2h = json.loads(MAGCAL.read_text())
    m2i = json.loads(FRESHRNA.read_text())
    return {
        **base,
        "outcome": outcome,
        "reading": reading,
        "primary": held_side,
        "decision_inputs": {
            "held_out_ci": [lo, hi],
            "scale_ci": [scale_lo, scale_hi],
            "which_governs": rc.WHICH_GOVERNS,
            "which_scale_interval": rc.WHICH_SCALE_INTERVAL,
            "calibration_removes_was_unreachable": not (scale_lo <= 0.0 <= scale_hi),
        },
        "control_p": control_p,
        "control_delta_mse_median": sorted(control)[len(control) // 2],
        "control_delta_mse_max": max(control),
        "control_difference": rc.CONTROL_DIFFERENCE,
        "variance_share": held_side["variance_share"],
        "descriptives_not_verdicts": {
            "words": rc.DESCRIPTIVES,
            "training_side": train_side,
            "whole_population": all_side,
            "held_out_per_cell_delta_mse": per_cell(held),
        },
        "context_verbatim": {
            "m2h_multiplicative_reading": m2h["per_form"]["multiplicative"]["reading"],
            "m2h_held_out_delta_mse": m2h["per_form"]["multiplicative"]["held_out_delta_mse"],
            "m2h_held_out_ci": m2h["per_form"]["multiplicative"]["held_out_ci"],
            "m2i_verdict_sentence": m2i["verdict_sentence"],
            "m2i_variance_fraction": m2i["variance_fraction"],
            "scalefit_verdict": json.loads(SCALEFIT_CENSUS.read_text())["verdict"],
            "note": "quoted, not re-decided, and never strengthened",
        },
        "exposure": rc.EXPOSURE,
        "falsifier": rc.FALSIFIER,
        "stop_rule": rc.STOP_RULE,
        "limitations": list(rc.LIMITATIONS),
        "rna_read": (
            "7h's committed multiplicative P and D rows for this population's genes, the same rows items"
            " 7h, M2h and M2i read; no new RNA table and no new biosample"
        ),
        "alphagenome_requests": 0,
        "money": "none: no request is sent and no key is read",
    }


def manifest() -> dict:
    inputs = [
        mf.input_entry(REGISTRATION, partition=None),
        mf.input_entry(POPULATION, partition=None),
        mf.input_entry(MULTRNA, partition=None),
        mf.input_entry(MAGCAL, partition=None),
        mf.input_entry(FRESHRNA, partition=None),
        mf.input_entry(SCALEFIT_CENSUS, partition=None),
    ]
    return {
        "sources": [
            {"accession": "data/results/recovered_registration.json (step 2)", "version": "266f031"},
            {"accession": "data/results/recovered_population.json (step 1)", "version": "99e2ea5"},
            {
                "accession": "data/results/multrna_result.json (7h, the P and D rows read here)",
                "version": "pinned by sha256",
            },
            {
                "accession": "data/results/magcal_result.json (M2h, the fixed scale and the split)",
                "version": "ac5580b",
            },
            {"accession": "data/results/freshrna_result.json (M2i, context)", "version": "pinned"},
            {"accession": "data/results/scalefit_census.json (lane-scalefit)", "version": "7f65544"},
        ],
        "inputs": inputs,
        "input_count": len(inputs),
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "form": rc.FORM,
            "population_form": rc.POPULATION_FORM,
            "scale": rc.SCALE,
            "scale_is_refitted": False,
            "interval": list(rc.INTERVAL),
            "bootstrap_resamples": rc.BOOTSTRAP["resamples"],
            "bootstrap_seed": rc.BOOTSTRAP["seed"],
            "permutations": rc.PERMUTATIONS,
            "min_genes": rc.MIN_GENES,
        },
        "exclusions": [
            "the 536 recovered-and-runnable genes with no row in 7h's committed multiplicative table:"
            " they have no P and no D and cannot be read",
            "the training half of this population is read as a descriptive only; the endpoint is the"
            " held-out half, the genes whose RNA did not enter M2h's fit of the scale",
        ],
        "partitions": {
            "held_out": "the registered endpoint: genes with magcal.is_held_out(symbol) true",
            "training": "descriptive only: genes M2h fitted its scale on",
        },
        "code_cleanliness": mf.code_cleanliness(__file__, OWN_CODE),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-save", action="store_true")
    args = ap.parse_args()
    t0 = time.time()
    out = run()
    out["seconds"] = round(time.time() - t0, 1)
    print("OUTCOME", out["outcome"], flush=True)
    print("harness agrees", out["harness_check"]["agrees"], flush=True)
    if "primary" in out:
        p = out["primary"]
        print(
            f"held-out Delta MSE {p['delta_mse']:+.6f} [{p['ci'][0]:+.6f}, {p['ci'][1]:+.6f}]"
            f" over {p['pairs']} pairs, {p['genes']} genes; variance share {p['variance_share']:.4%}",
            flush=True,
        )
        print("control p", out["control_p"], flush=True)
    print("READING:", out["reading"], flush=True)
    if not args.no_save:
        print("saved", save_result(RESULT, out, manifest=manifest()))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
