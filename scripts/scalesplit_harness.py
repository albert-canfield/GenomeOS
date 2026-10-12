# SPDX-License-Identifier: AGPL-3.0-or-later
"""lane-scalesplit, step 1: the harness, run and committed BEFORE the registration.

    uv run python scripts/scalesplit_harness.py

It answers one question and no other: reading 7h's committed genes table through M2h's own functions,
at M2h's own resamples, seeds and permutation count, do M2h's COMMITTED GENE-LEVEL NUMBERS COME BACK?
That is scalesplit.HARNESS points (1) to (5). It applies M2h's sha256-of-symbol split, which is M2h's
and not this lane's (scalesplit.WHICH_SPLIT_IS_WHOSE), and it fits NO scale on any chromosome set and
computes NO Delta MSE on any chromosome subset: the first chromosome-fold number of this design is
computed by scripts/scalesplit_run.py, after the registration is committed.

Nothing of magcal's, multrna's or combine2's is edited or re-run. No model request, no network, no
money, 0 AlphaGenome requests.
"""

from __future__ import annotations

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

RESULT = "scalesplit_harness"
#: everything read, repository-relative (d72f959)
MR_RESULT = "data/results/multrna_result.json"
MAGCAL_RESULT = "data/results/magcal_result.json"
MAGCAL_REGISTRATION = "data/results/magcal_registration.json"
MR_REGISTRATION = "data/results/multrna_registration.json"
SOURCE_FILES = (MR_RESULT, MAGCAL_RESULT, MAGCAL_REGISTRATION, MR_REGISTRATION)
OWN_CODE = (
    "genomeos/attribution/scalesplit.py",
    "scripts/scalesplit_harness.py",
    "scripts/scalesplit_register.py",
    "scripts/scalesplit_run.py",
    "tests/test_scalesplit.py",
)


def load_rows() -> list[dict[str, Any]]:
    """7h's committed genes table as magcal's rows: {gene, symbol, chrom, n_rules, n_inhibiting, P, D}.

    The same columns scripts/magcal_run.py:load_rows reads, for the multiplicative form alone.
    """
    tab = json.loads(Path(MR_RESULT).read_text())["genes_table"]
    ix = {c: i for i, c in enumerate(tab["columns"])}
    return [
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


def table_checks(rows: list[dict[str, Any]]) -> tuple[dict[str, Any], list[str]]:
    """HARNESS (1), (2) and (3): the table is 7h's, at 7h's committed uncalibrated Delta MSE."""
    problems: list[str] = []
    want = ss.SEVEN_H_COMMITTED
    pairs = sum(len(r["D"]) for r in rows)
    got = mc.delta_mse(rows, 1.0)
    out = {
        "genes": len(rows),
        "committed_genes": want["genes"],
        "pairs": pairs,
        "committed_pairs": want["pairs"],
        "committed_delta_mse": want["delta_mse"],
        "recomputed_delta_mse_at_scale_1": got,
        "abs_difference": abs(got - want["delta_mse"]),
        "tolerance": ss.TOLERANCE,
    }
    if len(rows) != want["genes"]:
        problems.append(f"{len(rows)} genes, 7h committed {want['genes']}")
    if pairs != want["pairs"]:
        problems.append(f"{pairs} pairs, 7h committed {want['pairs']}")
    if abs(got - want["delta_mse"]) > ss.TOLERANCE:
        problems.append(
            f"recomputed Delta MSE at scale 1 {got:.6f} differs from 7h's committed"
            f" {want['delta_mse']:.6f} by more than {ss.TOLERANCE}"
        )
    if any(len(r["P"]) != len(ss.CELLS) or len(r["D"]) != len(ss.CELLS) for r in rows):
        problems.append(f"a row does not have {len(ss.CELLS)} P and D values")
    if any(not r["symbol"] for r in rows):
        problems.append("a row has an empty symbol")
    if any(not r["chrom"] for r in rows):
        problems.append("a row has an empty chromosome")
    syms = [r["symbol"] for r in rows]
    if len(set(syms)) != len(syms):
        problems.append(f"{len(syms) - len(set(syms))} symbols occur more than once")
    return out, problems


def reproduce_m2h(rows: list[dict[str, Any]]) -> tuple[dict[str, Any], list[str]]:
    """HARNESS (4): M2h's committed gene-level numbers, recomputed through M2h's own functions.

    M2h's split, M2h's resamples, M2h's seeds, M2h's permutation count. This is the check that makes
    every later difference attributable to the SPLIT and to nothing else.
    """
    problems: list[str] = []
    want = ss.M2H_COMMITTED
    train, held = mc.split_rows(rows)
    scale = mc.fit_scale(train)
    alo, ahi = mc.bootstrap_scale(train, mc.BOOTSTRAP["resamples"], mc.BOOTSTRAP["seed"])
    dm = mc.delta_mse(held, scale)
    lo, hi = mc.bootstrap_delta(held, scale, mc.BOOTSTRAP["resamples"], mc.BOOTSTRAP["seed"])
    perm = mc.permuted_control(train, held, mc.PERMUTATIONS["count"], mc.PERMUTATIONS["seed"])
    p_ctrl = (1 + sum(1 for x in perm if x >= dm)) / (1 + len(perm))
    outcome = mc.decide(lo, hi, alo, ahi)
    got = {
        "training_genes": len(train),
        "held_out_genes": len(held),
        "scale": scale,
        "scale_ci": [alo, ahi],
        "held_out_delta_mse": dm,
        "held_out_ci": [lo, hi],
        "control_p": p_ctrl,
        "outcome": outcome,
    }
    for k in ("training_genes", "held_out_genes", "outcome"):
        if got[k] != want[k]:
            problems.append(f"M2h {k}: recomputed {got[k]!r}, committed {want[k]!r}")
    for k in ("scale", "held_out_delta_mse", "control_p"):
        if abs(float(got[k]) - float(want[k])) > ss.M2H_TOLERANCE:
            problems.append(f"M2h {k}: recomputed {got[k]!r}, committed {want[k]!r}")
    for k in ("scale_ci", "held_out_ci"):
        for g, w in zip(got[k], want[k], strict=True):
            if abs(g - w) > ss.M2H_TOLERANCE:
                problems.append(f"M2h {k}: recomputed {got[k]!r}, committed {want[k]!r}")
    return {
        "split_rule": "M2h's, sha256 of the gene symbol (magcal.SPLIT); NOT this lane's split",
        "m2h_tolerance": ss.M2H_TOLERANCE,
        "recomputed": got,
        "committed": {k: want[k] for k in got},
        "reproduces": not problems,
    }, problems


def chromosome_checks(rows: list[dict[str, Any]]) -> tuple[dict[str, Any], list[str]]:
    """HARNESS (5): the two registered chromosome sets cover the table exactly once each.

    Gene COUNTS only. No scale and no Delta MSE is computed on any chromosome subset here.
    """
    problems: list[str] = []
    odd = set(ss.ODD_CHROMOSOMES)
    even = set(ss.EVEN_AND_SEX_CHROMOSOMES)
    table = sorted({r["chrom"] for r in rows})
    per_chrom = {c: sum(1 for r in rows if r["chrom"] == c) for c in table}
    n_odd = sum(v for c, v in per_chrom.items() if c in odd)
    n_even = sum(v for c, v in per_chrom.items() if c in even)
    if odd & even:
        problems.append(f"the two chromosome sets are not disjoint: {sorted(odd & even)}")
    unnamed = [c for c in table if ss.chromosome_set(c) == "unnamed"]
    if unnamed:
        problems.append(f"chromosomes in the table named by neither set: {unnamed}")
    extra = sorted((odd | even) - set(table))
    if extra:
        problems.append(f"chromosomes named by a set but absent from the table: {extra}")
    if n_odd + n_even != len(rows):
        problems.append(f"the two sets hold {n_odd} + {n_even} genes, the table has {len(rows)}")
    return {
        "chromosomes_in_table": table,
        "genes_per_chromosome": per_chrom,
        "odd_chromosomes": list(ss.ODD_CHROMOSOMES),
        "even_and_sex_chromosomes": list(ss.EVEN_AND_SEX_CHROMOSOMES),
        "genes_on_odd": n_odd,
        "genes_on_even_and_sex": n_even,
        "sets_are_disjoint": not (odd & even),
        "union_is_the_table": set(table) == (odd | even),
        "below_min_genes": {c: v for c, v in per_chrom.items() if v < ss.MIN_GENES},
    }, problems


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
            "split_reproduced": "M2h's, sha256 of the gene symbol; this lane's chromosome split is not run",
            "min_genes": ss.MIN_GENES,
            "form_interval": list(ss.FORM_INTERVAL),
            "bootstrap": {"resamples": ss.BOOTSTRAP_RESAMPLES, "seed": ss.BOOTSTRAP_SEED},
            "permutations": ss.PERMUTATIONS,
            "tolerance": ss.TOLERANCE,
            "m2h_tolerance": ss.M2H_TOLERANCE,
        },
        "exclusions": ["none: the population is 7h's, unchanged"],
        "partitions": {
            "reproduced": "M2h's sha256-of-symbol training and held-out genes, to recover its numbers",
            "this_lanes_own": "not computed here: the chromosome folds are run only after registration",
        },
        "code_cleanliness": mf.code_cleanliness(__file__, OWN_CODE),
    }


def main() -> int:
    t0 = time.time()
    entries = [
        mf.input_entry(
            p,
            partition=(
                "pinned by sha256; the source of the committed P and D"
                if p == MR_RESULT
                else "pinned by sha256; the committed numbers that must come back, quoted not re-decided"
            ),
        )
        for p in SOURCE_FILES
    ]
    rows = load_rows()
    tbl, p1 = table_checks(rows)
    rep, p2 = reproduce_m2h(rows)
    chr_, p3 = chromosome_checks(rows)
    problems = p1 + p2 + p3
    out = {
        "lane": ss.LANE,
        "step": "1 of 4: the harness, before the registration",
        "question": ss.QUESTION,
        "why_stronger": ss.WHY_STRONGER,
        "harness": ss.HARNESS,
        "which_split_is_whose": ss.WHICH_SPLIT_IS_WHOSE,
        "prior_exposure": ss.PRIOR_EXPOSURE,
        "table_checks": tbl,
        "m2h_reproduction": rep,
        "chromosome_checks": chr_,
        "problems": problems,
        "outcome": "HARNESS_REPRODUCES" if not problems else "HARNESS_DEFECT",
        "reading": (
            "M2h's committed gene-level numbers REPRODUCE from 7h's committed genes table through M2h's"
            " own functions at M2h's own seeds: training 8,181 genes, held out 8,328, fitted scale"
            " +0.2146, held-out Delta MSE +0.04339, control p = 0.000999, outcome BETTER. Nothing of this"
            " lane's chromosome split is fitted or read here"
            if not problems
            else mc.reading("HARNESS_DEFECT", what="; ".join(problems))
        ),
        "nothing_of_this_lane_was_read": (
            "no scale was fitted on any chromosome set and no Delta MSE was computed on any chromosome"
            " subset by this script; only gene COUNTS per chromosome, which depend on no measured value"
        ),
        "alphagenome_requests": 0,
        "money": "none: no request is sent and no key is read",
        "seconds": round(time.time() - t0, 1),
    }
    print("saved", save_result(RESULT, out, manifest=manifest(entries)))
    print(out["outcome"], "|", out["reading"][:400])
    return 0 if not problems else 1


if __name__ == "__main__":
    raise SystemExit(main())
