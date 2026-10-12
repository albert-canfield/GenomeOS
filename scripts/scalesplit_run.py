# SPDX-License-Identifier: AGPL-3.0-or-later
"""lane-scalesplit, the registered run: the scale refitted on training CHROMOSOMES, read on held-out
chromosomes, on both complementary folds.

    uv run python scripts/scalesplit_run.py

Runs only from the committed registration (data/results/scalesplit_registration.json); refuses if any
registered field differs from genomeos/attribution/scalesplit.py or any pinned file's sha256 differs
from the file on disk.

Order, as registered (scalesplit.DECISION):
  1. 7h's committed genes table is read and the five HARNESS points are checked, M2h's committed
     gene-level numbers among them. Any miss stops the run here, before anything of this lane is split
     by chromosome, fitted or read (HARNESS_DEFECT);
  2. the chromosome partition is applied (scalesplit.SPLIT, the parity of the chromosome name);
  3. per fold, one scale is REFITTED on the TRAINING CHROMOSOMES' pairs alone and bootstrapped on
     those chromosomes' genes;
  4. only then is the held-out Delta MSE computed at that scale, with its interval and the within-gene
     cell permutation control, and M2h's four-way table applied as written by magcal.decide;
  5. the transport verdict is scalesplit.TRANSPORT_RULE on the two folds' outcomes, and the registered
     reading puts it first;
  6. the post hoc of scalesplit.POST_HOC is computed last and decides nothing.

Every function that fits, bootstraps, permutes, decides or words an outcome is magcal's, imported
read-only. No file of magcal's, multrna's or combine2's is edited. No model request, no network, no
money, 0 AlphaGenome requests.
"""

from __future__ import annotations

import hashlib
import json
import os
import sys
import time
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
os.chdir(ROOT)

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import magcal as mc  # noqa: E402
from genomeos.attribution import scalesplit as ss  # noqa: E402
from genomeos.results import save_result  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from scalesplit_harness import (  # noqa: E402
    MAGCAL_RESULT,
    MR_RESULT,
    OWN_CODE,
    SOURCE_FILES,
    chromosome_checks,
    load_rows,
    reproduce_m2h,
    table_checks,
)

RESULT = "scalesplit_result"
REGISTRATION = "data/results/scalesplit_registration.json"
#: step 1's committed result; read because check_registration re-hashes every path the registration
#: pins, and the registration pins it. Declared in inputs() for that reason (manifest contract).
HARNESS_RESULT = "data/results/scalesplit_harness.json"
#: M2i's replication, quoted as context and not re-decided
M2I = {
    "lane": "lane-freshrna",
    "commit": "f6e71c6",
    "held_out_delta_mse": "+0.03727 [+0.03202, +0.04299]",
    "note": "nothing fitted, on RNA fresh to this project; quoted, not re-decided here",
}


def check_registration() -> dict[str, Any]:
    """Refuse to run unless every BOUND field and every pinned sha256 is what was registered."""
    reg_path = Path(REGISTRATION)
    if not reg_path.exists():
        raise SystemExit(f"refused: {REGISTRATION} does not exist; register before running")
    reg = json.loads(reg_path.read_text())
    bad = []
    for field, name in ss.BOUND.items():
        want = getattr(ss, name)
        want = list(want) if isinstance(want, tuple) else want
        if json.dumps(reg.get(field), sort_keys=True) != json.dumps(want, sort_keys=True):
            bad.append(f"registered field {field} differs from scalesplit.{name}")
    for path, sha in reg["pinned_sha256"].items():
        p = Path(path)
        got = hashlib.sha256(p.read_bytes()).hexdigest() if p.exists() else "missing"
        if got != sha:
            bad.append(f"{path}: sha256 {got} != registered {sha}")
    if bad:
        raise SystemExit("refused:\n  " + "\n  ".join(bad))
    return reg


def harness(rows: list[dict[str, Any]]) -> dict[str, Any]:
    """HARNESS, all five points, through the committed harness step's own functions."""
    tbl, p1 = table_checks(rows)
    rep, p2 = reproduce_m2h(rows)
    chr_, p3 = chromosome_checks(rows)
    return {
        "table_checks": tbl,
        "m2h_reproduction": rep,
        "chromosome_checks": chr_,
        "problems": p1 + p2 + p3,
    }


def _sides(rs: list[dict[str, Any]], scale: float) -> dict[str, Any]:
    """DESCRIPTIVES: the size the calibration keeps, as numbers. magcal's sums, nothing recomputed."""
    dp, pp, dd, n = mc._sums(rs)
    return {
        "genes": len(rs),
        "pairs": n,
        "mean_D2": dd / n if n else None,
        "mean_P2": pp / n if n else None,
        "mean_scaled_P2": scale * scale * pp / n if n else None,
        "mean_D_times_P": dp / n if n else None,
    }


def score_fold(fold: str, rows: list[dict[str, Any]]) -> dict[str, Any]:
    """One registered fold: the refitted scale, the held-out comparison, the control, the verdict.

    The structure of scripts/magcal_run.py:score_form with magcal's chromosome-blind split replaced by
    scalesplit.split_rows_by_chromosome. Every fit, bootstrap, permutation, decision and word is
    magcal's.
    """
    words = ss.FOLD_WORDS[fold]
    train, held = ss.split_rows_by_chromosome(rows, fold)
    scale = mc.fit_scale(train)
    base: dict[str, Any] = {
        "fold_words": words,
        "split": {
            "training_chromosomes": list(ss.FOLD_TRAINING[fold]),
            "held_out_chromosomes": list(ss.FOLD_HELD_OUT[fold]),
            "training_genes": len(train),
            "held_out_genes": len(held),
            "rule": ss.SPLIT,
            "refitted": ss.REFITTED,
        },
    }
    if scale is None:
        clause = mc.scale_clause("not_fittable")
        return {
            **base,
            "outcome": "NOT_FITTABLE",
            "scale": None,
            "scale_clause": clause,
            "reading": mc.reading("NOT_FITTABLE", form=words, scale=clause),
        }
    alo, ahi = mc.bootstrap_scale(train, ss.BOOTSTRAP_RESAMPLES, ss.BOOTSTRAP_SEED)
    tn = sum(len(r["D"]) for r in train)
    clause = mc.scale_clause(
        "fitted",
        a=f"{scale:+.4g}",
        ci=ss.FORM_INTERVAL_WORDS,
        alo=f"{alo:+.4g}",
        ahi=f"{ahi:+.4g}",
        tg=f"{len(train):,}",
        tn=f"{tn:,}",
    )
    negative = mc.scale_is_negative(alo, ahi)
    if negative:
        clause += "; " + mc.scale_clause(
            "negative",
            a=f"{scale:+.4g}",
            ci=ss.FORM_INTERVAL_WORDS,
            alo=f"{alo:+.4g}",
            ahi=f"{ahi:+.4g}",
            tg=f"{len(train):,}",
        )
    out: dict[str, Any] = {
        **base,
        "scale": scale,
        "scale_ci": [alo, ahi],
        "scale_interval_includes_zero": bool(alo <= 0.0 <= ahi),
        "scale_is_negative": bool(negative),
        "scale_clause": clause,
        "training_delta_mse_at_refitted_scale": mc.delta_mse(train, scale),
        "training_delta_mse_at_scale_1": mc.delta_mse(train, 1.0),
    }
    if len(held) < ss.MIN_GENES:
        out["outcome"] = "NOT_TESTABLE"
        out["reading"] = mc.reading(
            "NOT_TESTABLE", form=words, g=f"{len(held):,}", min=ss.MIN_GENES, scale=clause
        )
        return out
    dm = mc.delta_mse(held, scale)
    lo, hi = mc.bootstrap_delta(held, scale, ss.BOOTSTRAP_RESAMPLES, ss.BOOTSTRAP_SEED)
    perm = mc.permuted_control(train, held, ss.PERMUTATIONS["count"], ss.PERMUTATIONS["seed"])
    p_ctrl = (1 + sum(1 for x in perm if x >= dm)) / (1 + len(perm))
    kind = mc.decide(lo, hi, alo, ahi)
    n = sum(len(r["D"]) for r in held)
    m2h = ss.M2H_COMMITTED
    out.update(
        {
            "outcome": kind,
            "held_out_delta_mse": dm,
            "held_out_ci": [lo, hi],
            "control_p": p_ctrl,
            "reading": mc.reading(
                kind,
                form=words,
                d=f"{dm:+.4g}",
                ci=ss.FORM_INTERVAL_WORDS,
                lo=f"{lo:+.4g}",
                hi=f"{hi:+.4g}",
                n=f"{n:,}",
                g=f"{len(held):,}",
                p=f"{p_ctrl:.3g}",
                scale=clause,
            ),
            "magnitude_against_m2h": [
                ss.magnitude_finding(
                    lo, hi, m2h["held_out_ci"][0], m2h["held_out_ci"][1], f"{fold} held-out Delta MSE"
                ),
                ss.magnitude_finding(
                    alo, ahi, m2h["scale_ci"][0], m2h["scale_ci"][1], f"{fold} fitted scale"
                ),
            ],
            "descriptives_not_verdicts": {
                "held_out_delta_mse_at_scale_1": mc.delta_mse(held, 1.0),
                "held_out_ci_at_scale_1": list(
                    mc.bootstrap_delta(held, 1.0, ss.BOOTSTRAP_RESAMPLES, ss.BOOTSTRAP_SEED)
                ),
                "training_side": _sides(train, scale),
                "held_out_side": _sides(held, scale),
                "held_out_delta_mse_by_cell": {
                    c: sum(2 * scale * r["P"][j] * r["D"][j] - scale * scale * r["P"][j] ** 2 for r in held)
                    / len(held)
                    for j, c in enumerate(ss.CELLS)
                },
                "predicted_range_log2_quantiles_held_out": mc.quantiles(
                    [max(r["P"]) - min(r["P"]) for r in held]
                ),
                "scaled_predicted_range_log2_quantiles_held_out": mc.quantiles(
                    [abs(scale) * (max(r["P"]) - min(r["P"])) for r in held]
                ),
                "control_delta_mse_median": sorted(perm)[len(perm) // 2],
                "control_delta_mse_95th_percentile": sorted(perm)[int(0.95 * len(perm)) - 1],
            },
        }
    )
    return out


def post_hoc(rows: list[dict[str, Any]]) -> dict[str, Any]:
    """POST_HOC: the per-chromosome scales, deciding nothing. No interval, no comparison, no verdict."""
    chroms = sorted({r["chrom"] for r in rows})
    per: dict[str, Any] = {}
    for c in chroms:
        own = [r for r in rows if r["chrom"] == c]
        others = [r for r in rows if r["chrom"] != c]
        per[c] = {
            "genes": len(own),
            "pairs": sum(len(r["D"]) for r in own),
            "scale_on_this_chromosome_alone": mc.fit_scale(own),
            "scale_on_the_other_23": mc.fit_scale(others),
            "below_min_genes": len(own) < ss.MIN_GENES,
            "chromosome_set": ss.chromosome_set(c),
        }
    vals = sorted(
        v["scale_on_this_chromosome_alone"]
        for v in per.values()
        if v["scale_on_this_chromosome_alone"] is not None
    )
    m2h_lo, m2h_hi = ss.M2H_COMMITTED["scale_ci"]
    inside = [
        c
        for c, v in per.items()
        if v["scale_on_this_chromosome_alone"] is not None
        and m2h_lo <= v["scale_on_this_chromosome_alone"] <= m2h_hi
    ]
    return {
        "what_this_is": ss.POST_HOC,
        "per_chromosome": per,
        "spread_of_the_within_chromosome_scales": {
            "n": len(vals),
            "min": vals[0] if vals else None,
            "q1": vals[len(vals) // 4] if vals else None,
            "median": vals[len(vals) // 2] if vals else None,
            "q3": vals[(3 * len(vals)) // 4] if vals else None,
            "max": vals[-1] if vals else None,
            "inside_m2h_scale_interval": sorted(inside),
            "inside_m2h_scale_interval_count": len(inside),
            "m2h_scale_interval": [m2h_lo, m2h_hi],
        },
        "deciding_nothing": (
            "no interval is placed on any per-chromosome scale, no comparison between chromosomes is"
            " registered, and none of this can change a fold's verdict or the transport verdict"
        ),
    }


def manifest(entries: list[dict[str, Any]]) -> dict[str, Any]:
    return {
        "sources": [
            {"accession": "data/results/multrna_result.json (7h, lane-multrna)", "version": "pinned"},
            {"accession": "data/results/magcal_result.json (M2h, lane-magcal)", "version": "pinned"},
            {
                "accession": "ENCODE polyA plus RNA-seq gene quantifications, V29, 13 reader biosamples",
                "version": "read by 7d, 7h and M2h; this lane reads their committed rows, not the files",
            },
        ],
        "inputs": entries,
        "input_count": len(entries),
        "assembly": "GRCh38; genes as 7h matched them, by gene_id",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "form": ss.FORM,
            "folds": list(ss.FOLDS),
            "split": ss.SPLIT,
            "refitted": ss.REFITTED,
            "min_genes": ss.MIN_GENES,
            "form_interval": list(ss.FORM_INTERVAL),
            "interval_decision": ss.INTERVAL,
            "bootstrap": ss.BOOTSTRAP,
            "permutations": ss.PERMUTATIONS,
            "control_alpha": ss.CONTROL_ALPHA,
            "tolerance": ss.TOLERANCE,
            "m2h_tolerance": ss.M2H_TOLERANCE,
            "free_parameters": "one real scale per fold; no per-gene and no per-cell parameter",
        },
        "exclusions": ["none of this lane's own: the population is 7h's, unchanged"],
        "partitions": {
            "rule": ss.SPLIT,
            "odd_held": {
                "training": "chr2, chr4, chr6, chr8, chr10, chr12, chr14, chr16, chr18, chr20, chr22,"
                " chrX, chrY: the scale is refitted here only",
                "held_out": "chr1, chr3, chr5, chr7, chr9, chr11, chr13, chr15, chr17, chr19, chr21:"
                " the registered endpoint is read here only",
            },
            "even_and_sex_held": {
                "training": "chr1, chr3, chr5, chr7, chr9, chr11, chr13, chr15, chr17, chr19, chr21:"
                " the scale is refitted here only",
                "held_out": "chr2, chr4, chr6, chr8, chr10, chr12, chr14, chr16, chr18, chr20, chr22,"
                " chrX, chrY: the registered endpoint is read here only",
            },
        },
        "code_cleanliness": mf.code_cleanliness(__file__, OWN_CODE),
    }


def main() -> int:
    t0 = time.time()
    reg = check_registration()
    entries = [
        mf.input_entry(REGISTRATION, partition="this lane's registration, re-checked field by field"),
        *(
            mf.input_entry(p, partition="pinned by sha256; the source of the committed P and D")
            for p in SOURCE_FILES
        ),
        mf.input_entry(
            HARNESS_RESULT,
            partition="step 1's committed harness, pinned by sha256 and re-hashed by check_registration",
        ),
    ]
    rows = load_rows()
    hc = harness(rows)
    base: dict[str, Any] = {
        "lane": ss.LANE,
        "registration": REGISTRATION,
        "registration_date": reg.get("date"),
        "question": ss.QUESTION,
        "why_stronger": ss.WHY_STRONGER,
        "prior_exposure": ss.PRIOR_EXPOSURE,
        "which_split_is_whose": ss.WHICH_SPLIT_IS_WHOSE,
        "refitted": ss.REFITTED,
        "interval_decision": ss.INTERVAL,
        "harness_check": hc,
        "alphagenome_requests": 0,
        "money": "none: no request is sent and no key is read",
    }
    if hc["problems"]:
        out = {
            **base,
            "outcome": "HARNESS_DEFECT",
            "reading": mc.reading("HARNESS_DEFECT", what="; ".join(hc["problems"])),
            "seconds": round(time.time() - t0, 1),
        }
        print("saved", save_result(RESULT, out, manifest=manifest(entries)))
        print(out["reading"])
        return 1
    per_fold = {f: score_fold(f, rows) for f in ss.FOLDS}
    outcomes = {f: per_fold[f]["outcome"] for f in ss.FOLDS}
    transport = ss.transport_verdict(outcomes)
    magnitude = [m for f in ss.FOLDS for m in per_fold[f].get("magnitude_against_m2h", [])]
    reading = ss.registered_reading(transport, [per_fold[f]["reading"] for f in ss.FOLDS], magnitude)
    mres = json.loads(Path(MR_RESULT).read_text())
    gres = json.loads(Path(MAGCAL_RESULT).read_text())
    out = {
        **base,
        "outcome": transport,
        "transport_verdict": transport,
        "transport_sentence": ss.TRANSPORT[transport],
        "transport_rule": ss.TRANSPORT_RULE,
        "reading": reading,
        "outcomes": outcomes,
        "scales": {f: per_fold[f].get("scale") for f in ss.FOLDS},
        "scale_cis": {f: per_fold[f].get("scale_ci") for f in ss.FOLDS},
        "held_out_delta_mse": {f: per_fold[f].get("held_out_delta_mse") for f in ss.FOLDS},
        "held_out_cis": {f: per_fold[f].get("held_out_ci") for f in ss.FOLDS},
        "per_fold": per_fold,
        "magnitude_comparison": ss.MAGNITUDE_COMPARISON,
        "magnitude_findings": magnitude,
        "post_hoc_deciding_nothing": post_hoc(rows),
        "which_governs": ss.WHICH_GOVERNS,
        "not_vacuous": ss.NOT_VACUOUS,
        "falsifier": ss.FALSIFIER,
        "prior_exposure_clause": ss.PRIOR_EXPOSURE_CLAUSE,
        "stop_rule": ss.STOP_RULE,
        "limitations": list(ss.LIMITATIONS),
        "context_verbatim": {
            "seven_h_reading": mres["reading"],
            "m2h_reading": gres["per_form"]["multiplicative"]["reading"],
            "m2h_committed": dict(ss.M2H_COMMITTED),
            "m2i_replication": dict(M2I),
            "note": (
                "7h's and M2h's readings are their own, quoted and not re-decided here; M2i's numbers are"
                " quoted from its commit message and are not recomputed"
            ),
        },
        "cells": list(ss.CELLS),
        "seconds": round(time.time() - t0, 1),
    }
    print("saved", save_result(RESULT, out, manifest=manifest(entries)))
    print(reading)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
