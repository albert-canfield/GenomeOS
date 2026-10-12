# SPDX-License-Identifier: AGPL-3.0-or-later
"""lane-scalehet, the run: is M2h's magnitude one global scale, or does it vary by region?

    uv run python scripts/scalehet_run.py

It refuses to start unless data/results/scalehet_registration.json exists, every BOUND field in it
equals the committed module's, and every pinned input's sha256 still matches the file on disk; and it
refuses to overwrite its own result. Then, exactly once: the per-chromosome fitted scales on the
admitted chromosomes, the global fitted scale, Q and R, the registered permutation null at its
registered count and seed, the two one-sided p values, the registered outcome through
`scalehet.decide_heterogeneity`, and the registered sentence for that outcome with its numbers filled.

No held-out Delta MSE of this lane is read, no per-region model is fitted whatever the outcome, no
file of magcal's, scalesplit's or multrna's is edited, no request is sent and no key is read.
"""

from __future__ import annotations

import json
import os
import sys
import time
from datetime import date
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
os.chdir(ROOT)

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import scalehet as sh  # noqa: E402
from genomeos.results import save_result  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from scalehet_harness import (  # noqa: E402
    MAGCAL_RESULT,
    MR_RESULT,
    OWN_CODE,
    READERS,
    SCALESPLIT_RESULT,
    SOURCE_FILES,
    load_rows,
)

RESULT = "scalehet_result"
OUT = "data/results/scalehet_result.json"
REGISTRATION = "data/results/scalehet_registration.json"
HARNESS_RESULT = "data/results/scalehet_harness.json"


def refuse_unless_registered() -> dict[str, Any]:
    """The registration must exist, bind this module word for word, and pin these inputs."""
    if Path(OUT).exists():
        raise SystemExit(f"refused: {OUT} already exists; the run is not re-run (STOP_RULE)")
    if not Path(REGISTRATION).exists():
        raise SystemExit(f"refused: {REGISTRATION} does not exist; register and commit it first")
    if not Path(HARNESS_RESULT).exists():
        raise SystemExit(f"refused: {HARNESS_RESULT} does not exist; run and commit the harness first")
    reg = json.loads(Path(REGISTRATION).read_text())
    if json.loads(Path(HARNESS_RESULT).read_text())["outcome"] != "HARNESS_REPRODUCES":
        raise SystemExit("refused: the committed harness did not reproduce M2h's numbers")
    drift = [f for f, v in sh.registered_fields().items() if reg.get(f) != v]
    if drift:
        raise SystemExit(f"refused: the registration differs from the module in {drift}")
    pinned = reg["pinned_sha256"]
    moved = [p for p, want in pinned.items() if mf.input_entry(p)["sha256"] != want]
    if moved:
        raise SystemExit(f"refused: a pinned input moved since the registration: {moved}")
    return reg


def observed(rows: list[dict[str, Any]]) -> tuple[dict[str, Any], list[str], list[str], list[int]]:
    """The admitted chromosomes' scales, the global scale, Q and R, through the registered fit."""
    counts = sh.gene_counts(rows)
    adm, exc = sh.admitted_chromosomes(counts)
    dps, pps = sh.gene_sums(rows)
    obs = sh.observed_statistics(dps, pps, [r["chrom"] for r in rows], adm)
    return obs, adm, exc, [counts[c] for c in adm]


def per_chromosome_table(
    rows: list[dict[str, Any]],
    obs: dict[str, Any],
    adm: list[str],
    exc: list[str],
    null: dict[str, Any],
) -> list[dict[str, Any]]:
    """DESCRIPTIVES: every chromosome's own fitted scale with its counts and its null band.

    The bands carry NO correction over the admitted chromosomes and are NOT a test of any one of them.
    An excluded chromosome's scale is reported and FLAGGED, and enters no statistic.
    """
    counts = sh.gene_counts(rows)
    dps, pps = sh.gene_sums(rows)
    out: list[dict[str, Any]] = []
    for i, c in enumerate(adm):
        lo, hi = sh.percentile_range(null["per_group"][i])
        a = obs["scales"][i]
        out.append(
            {
                "chrom": c,
                "admitted": True,
                "genes": counts[c],
                "pairs": counts[c] * len(sh.CELLS),
                "dp": obs["dp"][i],
                "pp": obs["pp"][i],
                "scale": a,
                "deviation_from_global": a - obs["global_scale"],
                "weighted_squared_deviation": obs["pp"][i] * (a - obs["global_scale"]) ** 2,
                "null_scale_98_33_range": [lo, hi],
                "inside_its_own_null_band": bool(lo <= a <= hi),
                "inside_m2h_scale_interval": bool(
                    sh.M2H_COMMITTED["scale_ci"][0] <= a <= sh.M2H_COMMITTED["scale_ci"][1]
                ),
            }
        )
    for c in exc:
        idx = [j for j, r in enumerate(rows) if r["chrom"] == c]
        out.append(
            {
                "chrom": c,
                "admitted": False,
                "flag": (
                    f"BELOW MIN_GENES = {sh.MIN_GENES} at {counts[c]} genes: reported only, and it"
                    " enters no statistic, no null and no verdict"
                ),
                "genes": counts[c],
                "pairs": counts[c] * len(sh.CELLS),
                "dp": sum(dps[j] for j in idx),
                "pp": sum(pps[j] for j in idx),
                "scale": sh.scale_from_sums(sum(dps[j] for j in idx), sum(pps[j] for j in idx)),
                "null_scale_98_33_range": None,
            }
        )
    return out


def context() -> dict[str, Any]:
    """DESCRIPTIVES, second paragraph: the committed readings, copied verbatim, not re-decided."""
    mres = json.loads(Path(MR_RESULT).read_text())
    gres = json.loads(Path(MAGCAL_RESULT).read_text())
    sres = json.loads(Path(SCALESPLIT_RESULT).read_text())
    return {
        "seven_h": {"source": MR_RESULT, "outcome": mres["outcome"], "reading": mres["reading"]},
        "m2h": {
            "source": MAGCAL_RESULT,
            "outcome": gres["per_form"]["multiplicative"]["outcome"],
            "reading": gres["per_form"]["multiplicative"]["reading"],
            "committed": dict(sh.M2H_COMMITTED),
        },
        "scalesplit": {
            "source": SCALESPLIT_RESULT,
            "outcome": sres["outcome"],
            "transport_verdict": sres["transport_verdict"],
            "reading": sres["reading"],
        },
        "note": (
            "quoted and not re-decided. Nothing in this lane revises M2h's about 2% of the variance or"
            " scalesplit's chromosome-held-out reading, and a scale remains a MAGNITUDE CORRECTION and"
            " never an explanation of expression"
        ),
    }


def manifest(entries: list[dict[str, Any]], adm: list[str], exc: list[str]) -> dict[str, Any]:
    return {
        "sources": [
            {"accession": "data/results/multrna_result.json (7h, lane-multrna)", "version": "pinned"},
            {"accession": "data/results/magcal_result.json (M2h, lane-magcal)", "version": "pinned"},
            {
                "accession": "data/results/scalesplit_result.json (lane-scalesplit)",
                "version": "pinned; its post hoc per-chromosome spread is the pattern tested here",
            },
            {
                "accession": "ENCODE polyA plus RNA-seq gene quantifications, V29, 13 reader biosamples",
                "version": "read by 7d, 7h and M2h; this lane reads their committed rows, not the files",
            },
        ],
        "inputs": entries,
        "input_count": len(entries),
        "readers": READERS,
        "reader_count": len(READERS),
        "assembly": "GRCh38; genes as 7h matched them, by gene_id",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "form": sh.FORM,
            "grouping": "by chromosome, as 7h's committed table names it",
            "statistic": ("Q = sum over admitted chromosomes of pp_c * (a_c - a_g)^2; R = max a_c - min a_c"),
            "null": ("the chromosome labels permuted among the genes, each group's gene count held exactly"),
            "min_genes": sh.MIN_GENES,
            "min_chromosomes": sh.MIN_CHROMOSOMES,
            "admitted_chromosomes": adm,
            "permutations": sh.PERMUTATIONS,
            "control_alpha": sh.CONTROL_ALPHA,
            "form_interval": list(sh.FORM_INTERVAL),
            "tolerance": sh.TOLERANCE,
            "free_parameters": (
                "one closed-form scale per admitted chromosome and one global scale; the null refits"
                " every one of them on each draw, so nothing fitted here is unpriced"
            ),
        },
        "exclusions": [
            f"{c}: fewer than MIN_GENES = {sh.MIN_GENES} genes in 7h's committed table, by the floor"
            " fixed in the registration before any value of this lane was read"
            for c in exc
        ]
        or ["none: every chromosome in 7h's table cleared the registered gene floor"],
        "partitions": {
            "grouping": "by chromosome; a gene is wholly in one group with all thirteen of its pairs",
            "held_out": (
                "none, by registration (scalehet.NO_SPLIT): the statistic is a dispersion of FITTED"
                " scales and the permutation null prices the in-sample fitting exactly"
            ),
        },
        "code_cleanliness": mf.code_cleanliness(__file__, OWN_CODE),
    }


def main() -> int:
    t0 = time.time()
    reg = refuse_unless_registered()
    rows = load_rows()
    obs, adm, exc, sizes = observed(rows)
    dps, pps = sh.gene_sums(rows)
    out: dict[str, Any]
    if obs.get("degenerate"):
        numbers: dict[str, Any] = {"what": obs["degenerate"]}
        outcome = sh.decide_heterogeneity(1.0, 1.0, len(adm), obs["degenerate"])
        out = {
            "outcome": outcome,
            "degenerate": obs["degenerate"],
            "p_q": 1.0,
            "p_r": 1.0,
            "admitted_count": len(adm),
            "reading_numbers": numbers,
            "reading": sh.full_reading(outcome, **numbers),
        }
    else:
        keep = [i for i, r in enumerate(rows) if r["chrom"] in set(adm)]
        null = sh.permutation_null(
            [dps[i] for i in keep],
            [pps[i] for i in keep],
            sizes,
            sh.PERMUTATIONS["count"],
            sh.PERMUTATIONS["seed"],
        )
        degenerate = null["degenerate"]
        q_lo = q_hi = float("nan")
        if not degenerate:
            q_lo, q_hi = sh.percentile_range(null["q"])
            if q_hi - q_lo == 0.0:
                degenerate = "the null Q distribution has zero width: the permutation gave no spread"
        if degenerate:
            numbers = {"what": degenerate}
            outcome = sh.decide_heterogeneity(1.0, 1.0, len(adm), degenerate)
            out = {
                "outcome": outcome,
                "degenerate": degenerate,
                "p_q": 1.0,
                "p_r": 1.0,
                "admitted_count": len(adm),
                "reading_numbers": numbers,
                "reading": sh.full_reading(outcome, **numbers),
            }
        else:
            r_lo, r_hi = sh.percentile_range(null["r"])
            p_q = sh.permutation_p(obs["q"], null["q"])
            p_r = sh.permutation_p(obs["r"], null["r"])
            outcome = sh.decide_heterogeneity(p_q, p_r, len(adm), None)
            i_min = min(range(len(adm)), key=lambda i: obs["scales"][i])
            i_max = max(range(len(adm)), key=lambda i: obs["scales"][i])
            q_sorted = sorted(null["q"])
            r_sorted = sorted(null["r"])
            numbers = {
                "n_admitted": len(adm),
                "a_global": obs["global_scale"],
                "q_obs": obs["q"],
                "q_med": q_sorted[len(q_sorted) // 2],
                "q_lo": q_lo,
                "q_hi": q_hi,
                "p_q": p_q,
                "r_obs": obs["r"],
                "a_min": obs["scales"][i_min],
                "a_max": obs["scales"][i_max],
                "chrom_min": adm[i_min],
                "chrom_max": adm[i_max],
                "r_med": r_sorted[len(r_sorted) // 2],
                "r_lo": r_lo,
                "r_hi": r_hi,
                "p_r": p_r,
                "m2h_scale": sh.M2H_COMMITTED["scale"],
                "count": sh.PERMUTATIONS["count"],
                "seed": sh.PERMUTATIONS["seed"],
            }
            table = per_chromosome_table(rows, obs, adm, exc, null)
            scales = list(obs["scales"])
            out = {
                "outcome": outcome,
                "degenerate": None,
                "q_observed": obs["q"],
                "r_observed": obs["r"],
                "global_scale": obs["global_scale"],
                "p_q": p_q,
                "p_r": p_r,
                "admitted": adm,
                "admitted_count": len(adm),
                "excluded_below_min_genes": exc,
                "genes_admitted": sum(sizes),
                "pairs_admitted": sum(sizes) * len(sh.CELLS),
                "null_q": null["q"],
                "null_r": null["r"],
                "null_q_summary": {
                    "median": numbers["q_med"],
                    "percentile_98_33": [q_lo, q_hi],
                    "min": q_sorted[0],
                    "max": q_sorted[-1],
                },
                "null_r_summary": {
                    "median": numbers["r_med"],
                    "percentile_98_33": [r_lo, r_hi],
                    "min": r_sorted[0],
                    "max": r_sorted[-1],
                },
                "null_global_scale_is_invariant": (
                    "by construction: permuting the labels redistributes the same genes, so the total"
                    " sums and therefore the global fitted scale are the same in every draw"
                ),
                "per_chromosome": table,
                "reading_numbers": numbers,
                "reading": sh.full_reading(outcome, **numbers),
                "post_hoc_deciding_nothing": {
                    "rule": (
                        "every number in this block decides nothing and can change no verdict; the"
                        " per-chromosome null bands carry NO correction over the admitted chromosomes"
                        " and are not a test of any chromosome"
                    ),
                    "scale_quantiles": sh.quantiles(scales),
                    "scale_min": min(scales),
                    "scale_max": max(scales),
                    "negative_scales": sum(1 for a in scales if a < 0.0),
                    "inside_m2h_scale_interval": sum(
                        1 for c in table if c["admitted"] and c["inside_m2h_scale_interval"]
                    ),
                    "inside_own_null_band": sum(
                        1 for c in table if c["admitted"] and c["inside_its_own_null_band"]
                    ),
                    "p_q_against_0_0167": (
                        "0.0167 is the 0.05/3 Bonferroni divisor M2h's 98.33% interval came from; the"
                        " verdict uses 0.05 and only 0.05 (scalehet.ALPHA)"
                    ),
                    "q_rejects_at_0_0167": bool(p_q <= 0.05 / 3),
                    "r_rejects_at_0_0167": bool(p_r <= 0.05 / 3),
                    "global_scale_all_24": sh.scale_from_sums(sum(dps), sum(pps)),
                    "m2h_committed_scale": sh.M2H_COMMITTED["scale"],
                    "m2h_committed_scale_ci": sh.M2H_COMMITTED["scale_ci"],
                },
            }
    out.update(
        {
            "lane": sh.LANE,
            "date": date.today().isoformat(),
            "step": "4 of 4: the run, after the registration was committed",
            "registration": REGISTRATION,
            "registration_date": reg["date"],
            "harness": HARNESS_RESULT,
            "question": sh.QUESTION,
            "what_is_not_asked": sh.WHAT_IS_NOT_ASKED,
            "grouping": sh.GROUPING,
            "no_split": sh.NO_SPLIT,
            "statistic": list(sh.STATISTIC),
            "null_description": list(sh.NULL),
            "permutations": sh.PERMUTATIONS,
            "alpha": sh.ALPHA,
            "interval": sh.INTERVAL,
            "which_governs": sh.WHICH_GOVERNS,
            "not_vacuous": sh.NOT_VACUOUS,
            "falsifier": sh.FALSIFIER,
            "prior_exposure": list(sh.PRIOR_EXPOSURE),
            "prior_exposure_clause": sh.PRIOR_EXPOSURE_CLAUSE,
            "stop_rule": sh.STOP_RULE,
            "limitations": list(sh.LIMITATIONS),
            "descriptives": list(sh.DESCRIPTIVES),
            "context_verbatim": context(),
            "no_regional_model_was_fitted": (
                "whatever the outcome above: no per-region scale, no per-chromosome model and no cause"
                " of a difference was fitted, proposed or named by this run. That is DECISION step (8),"
                " fixed before any statistic was computed"
            ),
            "cells": list(sh.CELLS),
            "readers": READERS,
            "alphagenome_requests": 0,
            "money": "none: no request is sent and no key is read",
            "seconds": round(time.time() - t0, 1),
        }
    )
    entries = [
        mf.input_entry(p, partition="pinned by sha256 in the registration")
        for p in (*SOURCE_FILES, HARNESS_RESULT, REGISTRATION)
    ]
    print("saved", save_result(RESULT, out, manifest=manifest(entries, adm, exc)))
    print(out["outcome"])
    print(out["reading"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
