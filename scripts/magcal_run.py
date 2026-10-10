# SPDX-License-Identifier: AGPL-3.0-or-later
"""lane-magcal, the registered run: one scale per form, fitted on training genes, read on held-out ones.

    uv run python scripts/magcal_run.py

Runs only from the committed registration (data/results/magcal_registration.json); refuses if any
registered field differs from genomeos/attribution/magcal.py or any pinned file's sha256 differs from
the file on disk.

Order, as registered (magcal.DECISION):
  1. the three committed genes tables are read and the HARNESS checked on all of them -- each form's
     gene count, pair count and uncalibrated Delta MSE against its committed result. Any miss stops the
     run here, before anything is split, fitted or read (HARNESS_DEFECT);
  2. the split is applied (magcal.SPLIT, lane-saturable's sha256 of the symbol);
  3. per form, one scale is fitted on the TRAINING pairs alone and bootstrapped on training genes;
  4. only then is the held-out Delta MSE computed at that scale, with its interval and the within-gene
     cell permutation control, and the four-way table of magcal.WHICH_GOVERNS applied as written.

No model request, no network, no money, 0 AlphaGenome requests.
"""

from __future__ import annotations

import hashlib
import json
import sys
import time
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import magcal as mc  # noqa: E402
from genomeos.results import RESULTS_DIR, save_result  # noqa: E402

RESULT = "magcal_result"
REGISTRATION = RESULTS_DIR / "magcal_registration.json"
C2_RESULT = RESULTS_DIR / "combine2_result.json"
MR_RESULT = RESULTS_DIR / "multrna_result.json"
C2_REGISTRATION = RESULTS_DIR / "combine2_registration.json"
MR_REGISTRATION = RESULTS_DIR / "multrna_registration.json"
#: everything the registration pins, in the order it pins it
SOURCE_FILES = (C2_RESULT, MR_RESULT, C2_REGISTRATION, MR_REGISTRATION)
OWN_CODE = (
    "genomeos/attribution/magcal.py",
    "scripts/magcal_register.py",
    "scripts/magcal_run.py",
    "tests/test_magcal.py",
)
#: HARNESS: the committed numbers each form's table must reproduce, read from the committed results
COMMITTED = {
    "mean": {"genes": 4725, "pairs": 61425, "delta_mse": -0.44348590854013453},
    "competition": {"genes": 2501, "pairs": 32513, "delta_mse": -0.8101520857173696},
    "multiplicative": {"genes": 16509, "pairs": 214617, "delta_mse": -0.5361286982098615},
}


def check_registration() -> dict[str, Any]:
    """Refuse to run unless every BOUND field and every pinned sha256 is what was registered."""
    if not REGISTRATION.exists():
        raise SystemExit(f"refused: {REGISTRATION} does not exist; register before running")
    reg = json.loads(REGISTRATION.read_text())
    bad = []
    for field, name in mc.BOUND.items():
        want = getattr(mc, name)
        want = list(want) if isinstance(want, tuple) else want
        if json.dumps(reg.get(field), sort_keys=True) != json.dumps(want, sort_keys=True):
            bad.append(f"registered field {field} differs from magcal.{name}")
    for path, sha in reg["pinned_sha256"].items():
        p = Path(path)
        got = hashlib.sha256(p.read_bytes()).hexdigest() if p.exists() else "missing"
        if got != sha:
            bad.append(f"{path}: sha256 {got} != registered {sha}")
    if bad:
        raise SystemExit("refused:\n  " + "\n  ".join(bad))
    return reg


def load_rows() -> dict[str, list[dict[str, Any]]]:
    """Each form's committed genes table as rows of {gene, symbol, chrom, n_rules, n_inhibiting, P, D}."""
    c2res = json.loads(C2_RESULT.read_text())
    mres = json.loads(MR_RESULT.read_text())
    tables = {
        "mean": c2res["genes_tables"]["mean"],
        "competition": c2res["genes_tables"]["competition"],
        "multiplicative": mres["genes_table"],
    }
    out: dict[str, list[dict[str, Any]]] = {}
    for form, tab in tables.items():
        ix = {c: i for i, c in enumerate(tab["columns"])}
        out[form] = [
            {
                "gene": r[ix["gene"]],
                "symbol": r[ix["symbol"]],
                "chrom": r[ix["chrom"]],
                "n_rules": r[ix["n_rules"]],
                "n_inhibiting": r[ix["n_inhibiting"]],
                "P": list(r[ix["P"]]),
                "D": list(r[ix["D"]]),
            }
            for r in tab["rows"]
        ]
    return out


def harness(rows: dict[str, list[dict[str, Any]]]) -> dict[str, Any]:
    """HARNESS: counts, the uncalibrated Delta MSE against the committed one, row shapes, uniqueness."""
    problems: list[str] = []
    per: dict[str, Any] = {}
    for form in mc.FORMS:
        rs = rows[form]
        want = COMMITTED[form]
        pairs = sum(len(r["D"]) for r in rs)
        got = mc.delta_mse(rs, 1.0)
        per[form] = {
            "genes": len(rs),
            "pairs": pairs,
            "committed_delta_mse": want["delta_mse"],
            "recomputed_delta_mse_at_scale_1": got,
            "abs_difference": abs(got - want["delta_mse"]),
        }
        if len(rs) != want["genes"]:
            problems.append(f"{form}: {len(rs)} genes, committed {want['genes']}")
        if pairs != want["pairs"]:
            problems.append(f"{form}: {pairs} pairs, committed {want['pairs']}")
        if abs(got - want["delta_mse"]) > mc.TOLERANCE:
            problems.append(
                f"{form}: recomputed Delta MSE {got:.6f} differs from committed"
                f" {want['delta_mse']:.6f} by more than {mc.TOLERANCE}"
            )
        if any(len(r["P"]) != len(mc.CELLS) or len(r["D"]) != len(mc.CELLS) for r in rs):
            problems.append(f"{form}: a row does not have {len(mc.CELLS)} P and D values")
        if any(not r["symbol"] for r in rs):
            problems.append(f"{form}: a row has an empty symbol")
        syms = [r["symbol"] for r in rs]
        if len(set(syms)) != len(syms):
            problems.append(f"{form}: {len(syms) - len(set(syms))} symbols occur more than once")
    return {"per_form": per, "problems": problems}


def _sides(rs: list[dict[str, Any]], scale: float) -> dict[str, Any]:
    dp, pp, dd, n = mc._sums(rs)
    return {
        "genes": len(rs),
        "pairs": n,
        "mean_D2": dd / n if n else None,
        "mean_P2": pp / n if n else None,
        "mean_scaled_P2": scale * scale * pp / n if n else None,
        "mean_D_times_P": dp / n if n else None,
    }


def score_form(form: str, rs: list[dict[str, Any]]) -> dict[str, Any]:
    """The registered scale, the held-out comparison, the control and the verdict for one form."""
    words = mc.FORM_WORDS[form]
    train, held = mc.split_rows(rs)
    scale = mc.fit_scale(train)
    base: dict[str, Any] = {
        "form_words": words,
        "split": {"training_genes": len(train), "held_out_genes": len(held), "rule": mc.SPLIT},
    }
    if scale is None:
        clause = mc.scale_clause("not_fittable")
        return {
            **base,
            "outcome": "NOT_FITTABLE",
            "scale": None,
            "reading": mc.reading("NOT_FITTABLE", form=words, scale=clause),
        }
    alo, ahi = mc.bootstrap_scale(train, mc.BOOTSTRAP["resamples"], mc.BOOTSTRAP["seed"])
    tn = sum(len(r["D"]) for r in train)
    clause = mc.scale_clause(
        "fitted",
        a=f"{scale:+.4g}",
        ci=mc.FORM_INTERVAL_WORDS,
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
            ci=mc.FORM_INTERVAL_WORDS,
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
        "training_delta_mse_at_fitted_scale": mc.delta_mse(train, scale),
        "training_delta_mse_at_scale_1": mc.delta_mse(train, 1.0),
    }
    if len(held) < mc.MIN_GENES:
        out["outcome"] = "NOT_TESTABLE"
        out["reading"] = mc.reading(
            "NOT_TESTABLE", form=words, g=f"{len(held):,}", min=mc.MIN_GENES, scale=clause
        )
        return out
    dm = mc.delta_mse(held, scale)
    lo, hi = mc.bootstrap_delta(held, scale, mc.BOOTSTRAP["resamples"], mc.BOOTSTRAP["seed"])
    perm = mc.permuted_control(train, held, mc.PERMUTATIONS["count"], mc.PERMUTATIONS["seed"])
    p_ctrl = (1 + sum(1 for x in perm if x >= dm)) / (1 + len(perm))
    kind = mc.decide(lo, hi, alo, ahi)
    n = sum(len(r["D"]) for r in held)
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
                ci=mc.FORM_INTERVAL_WORDS,
                lo=f"{lo:+.4g}",
                hi=f"{hi:+.4g}",
                n=f"{n:,}",
                g=f"{len(held):,}",
                p=f"{p_ctrl:.3g}",
                scale=clause,
            ),
            "descriptives_not_verdicts": {
                "held_out_delta_mse_at_scale_1": mc.delta_mse(held, 1.0),
                "held_out_ci_at_scale_1": list(
                    mc.bootstrap_delta(held, 1.0, mc.BOOTSTRAP["resamples"], mc.BOOTSTRAP["seed"])
                ),
                "training_side": _sides(train, scale),
                "held_out_side": _sides(held, scale),
                "held_out_delta_mse_by_cell": {
                    c: sum(2 * scale * r["P"][j] * r["D"][j] - scale * scale * r["P"][j] ** 2 for r in held)
                    / len(held)
                    for j, c in enumerate(mc.CELLS)
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


def manifest(entries: list[dict[str, Any]]) -> dict[str, Any]:
    return {
        "sources": [
            {"accession": "data/results/combine2_result.json (7d, lane-combine2)", "version": "pinned"},
            {"accession": "data/results/multrna_result.json (7h, lane-multrna)", "version": "pinned"},
            {
                "accession": "ENCODE polyA plus RNA-seq gene quantifications, V29, 13 reader biosamples",
                "version": "read by 7d and 7h; this lane reads their committed rows, not the files",
            },
        ],
        "inputs": entries,
        "input_count": len(entries),
        "assembly": "GRCh38; genes as 7d and 7h matched them, by gene_id",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "forms": list(mc.FORMS),
            "split": mc.SPLIT,
            "min_genes": mc.MIN_GENES,
            "form_interval": list(mc.FORM_INTERVAL),
            "bootstrap": mc.BOOTSTRAP,
            "permutations": mc.PERMUTATIONS,
            "control_alpha": mc.CONTROL_ALPHA,
            "tolerance": mc.TOLERANCE,
            "free_parameters": "one real scale per form; no per-gene and no per-cell parameter",
        },
        "exclusions": ["none of this lane's own: each form's population is 7d's or 7h's, unchanged"],
        "partitions": {
            "rule": mc.SPLIT,
            "training": "even first 32 bits of sha256(symbol): the scale is fitted here only",
            "held_out": "odd first 32 bits of sha256(symbol): the registered endpoint is read here only",
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
    ]
    rows = load_rows()
    hc = harness(rows)
    base = {
        "lane": mc.LANE,
        "registration": "data/results/magcal_registration.json",
        "registration_date": reg.get("date"),
        "question": mc.QUESTION,
        "prior_exposure": mc.PRIOR_EXPOSURE,
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
    per_form = {f: score_form(f, rows[f]) for f in mc.FORMS}
    reading = mc.combined_reading(per_form)
    c2res = json.loads(C2_RESULT.read_text())
    mres = json.loads(MR_RESULT.read_text())
    out = {
        **base,
        "reading": reading,
        "outcomes": {f: per_form[f]["outcome"] for f in mc.FORMS},
        "scales": {f: per_form[f]["scale"] for f in mc.FORMS},
        "scale_cis": {f: per_form[f].get("scale_ci") for f in mc.FORMS},
        "per_form": per_form,
        "which_governs": mc.WHICH_GOVERNS,
        "not_vacuous": mc.NOT_VACUOUS,
        "falsifier": mc.FALSIFIER,
        "prior_exposure_clause": mc.PRIOR_EXPOSURE_CLAUSE,
        "limitations": list(mc.LIMITATIONS),
        "context_verbatim": {
            "seven_d_reading": c2res["reading"],
            "seven_d_control_p": {f: c2res["per_form"][f]["control_p"] for f in ("mean", "competition")},
            "seven_h_reading": mres["reading"],
            "seven_h_control_p": mres["multiplicative"]["control_p"],
            "note": (
                "7d's and 7h's readings are their own, on the uncalibrated forms, quoted and not"
                " re-decided here; their shuffled-reader control is not this lane's control"
            ),
        },
        "cells": list(mc.CELLS),
        "seconds": round(time.time() - t0, 1),
    }
    print("saved", save_result(RESULT, out, manifest=manifest(entries)))
    print(reading)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
