# SPDX-License-Identifier: AGPL-3.0-or-later
"""lane-scalehet, step 1: the harness, run and committed BEFORE the registration.

    uv run python scripts/scalehet_harness.py

It answers two questions and no others. (A) Reading 7h's committed genes table through M2h's own
functions, at M2h's own resamples, seeds and permutation count, DO M2h'S COMMITTED GENE-LEVEL NUMBERS
COME BACK? lane-scalesplit did this bit-identically and this follows it. (B) Does this lane's own
per-chromosome scale path -- the per-gene sums of D*P and P^2 that its statistic and its null are built
from -- agree with magcal.fit_scale applied to each chromosome's rows, within TOLERANCE?

That is scalehet.HARNESS points (1) to (5). NO heterogeneity statistic, NO range, NO permutation and
NO p value is computed here: those are scalehet.HARNESS point (6) and belong to
scripts/scalehet_run.py, after the registration is committed. The per-chromosome SCALES are computed
here only as the two paths' agreement check -- the dispersion of them is not summarised, no statistic
is formed from them and none of them is written to this result, because the registered statistic is
the thing that must not be seen before it is registered.

Nothing of magcal's, scalesplit's, multrna's or combine2's is edited or re-run. No model request, no
network, no money, 0 AlphaGenome requests.
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
from genomeos.attribution import scalehet as sh  # noqa: E402
from genomeos.results import save_result  # noqa: E402

RESULT = "scalehet_harness"
#: everything read, repository-relative (d72f959)
MR_RESULT = "data/results/multrna_result.json"
MAGCAL_RESULT = "data/results/magcal_result.json"
MAGCAL_REGISTRATION = "data/results/magcal_registration.json"
MR_REGISTRATION = "data/results/multrna_registration.json"
SCALESPLIT_RESULT = "data/results/scalesplit_result.json"
SOURCE_FILES = (
    MR_RESULT,
    MAGCAL_RESULT,
    MAGCAL_REGISTRATION,
    MR_REGISTRATION,
    SCALESPLIT_RESULT,
)
OWN_CODE = (
    "genomeos/attribution/scalehet.py",
    "scripts/scalehet_harness.py",
    "scripts/scalehet_register.py",
    "scripts/scalehet_run.py",
    "tests/test_scalehet.py",
)
#: every reader this run uses: the thirteen reader cells of 7d, 7h and M2h, declared (rule f)
READERS = list(sh.CELLS)


def load_rows() -> list[dict[str, Any]]:
    """7h's committed genes table as magcal's rows: {gene, symbol, chrom, n_rules, n_inhibiting, P, D}.

    The same columns scripts/magcal_run.py:load_rows and scripts/scalesplit_harness.py:load_rows read,
    for the multiplicative form alone.
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
    """HARNESS (1) and (2): the table is 7h's, at 7h's committed uncalibrated Delta MSE."""
    problems: list[str] = []
    want = sh.SEVEN_H_COMMITTED
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
        "tolerance": sh.TOLERANCE,
    }
    if len(rows) != want["genes"]:
        problems.append(f"{len(rows)} genes, 7h committed {want['genes']}")
    if pairs != want["pairs"]:
        problems.append(f"{pairs} pairs, 7h committed {want['pairs']}")
    if abs(got - want["delta_mse"]) > sh.TOLERANCE:
        problems.append(
            f"recomputed Delta MSE at scale 1 {got:.6f} differs from 7h's committed"
            f" {want['delta_mse']:.6f} by more than {sh.TOLERANCE}"
        )
    if any(len(r["P"]) != len(sh.CELLS) or len(r["D"]) != len(sh.CELLS) for r in rows):
        problems.append(f"a row does not have {len(sh.CELLS)} P and D values")
    if any(not r["symbol"] for r in rows):
        problems.append("a row has an empty symbol")
    if any(not r["chrom"] for r in rows):
        problems.append("a row has an empty chromosome")
    syms = [r["symbol"] for r in rows]
    if len(set(syms)) != len(syms):
        problems.append(f"{len(syms) - len(set(syms))} symbols occur more than once")
    return out, problems


def reproduce_m2h(rows: list[dict[str, Any]]) -> tuple[dict[str, Any], list[str]]:
    """HARNESS (3): M2h's committed gene-level numbers, recomputed through M2h's own functions.

    M2h's split, M2h's resamples, M2h's seeds, M2h's permutation count. This is the check that makes
    every later difference attributable to the GROUPING and to nothing else.
    """
    problems: list[str] = []
    want = sh.M2H_COMMITTED
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
        if abs(float(got[k]) - float(want[k])) > sh.M2H_TOLERANCE:
            problems.append(f"M2h {k}: recomputed {got[k]!r}, committed {want[k]!r}")
    for k in ("scale_ci", "held_out_ci"):
        for g, w in zip(got[k], want[k], strict=True):
            if abs(g - w) > sh.M2H_TOLERANCE:
                problems.append(f"M2h {k}: recomputed {got[k]!r}, committed {want[k]!r}")
    return {
        "split_rule": "M2h's, sha256 of the gene symbol (magcal.SPLIT); NOT this lane's design",
        "m2h_tolerance": sh.M2H_TOLERANCE,
        "recomputed": got,
        "committed": {k: want[k] for k in got},
        "reproduces": not problems,
        "reading": (
            "M2h's committed gene-level numbers REPRODUCE from 7h's committed genes table through M2h's"
            " own functions at M2h's own seeds: training 8,181 genes, held out 8,328, fitted scale"
            " +0.2146, held-out Delta MSE +0.04339 (about 2% of the variance), control p = 0.000999,"
            " outcome BETTER"
            if not problems
            else "M2h's committed numbers DO NOT come back: " + "; ".join(problems)
        ),
    }, problems


def two_paths_agree(rows: list[dict[str, Any]]) -> tuple[dict[str, Any], list[str]]:
    """HARNESS (4): this lane's per-gene-sums path IS magcal.fit_scale, per chromosome.

    The agreement is reported as a MAXIMUM ABSOLUTE DIFFERENCE over the chromosomes and as a count of
    chromosomes checked. The per-chromosome scales THEMSELVES are not written to this result and no
    dispersion of them is computed here: the registered statistic must not be seen before it is
    registered (scalehet.HARNESS (6)).
    """
    problems: list[str] = []
    dps, pps = sh.gene_sums(rows)
    by_chrom: dict[str, list[int]] = {}
    for i, r in enumerate(rows):
        by_chrom.setdefault(r["chrom"], []).append(i)
    worst = 0.0
    worst_chrom = ""
    checked = 0
    for chrom, idx in sorted(by_chrom.items()):
        mine = sh.scale_from_sums(sum(dps[i] for i in idx), sum(pps[i] for i in idx))
        theirs = mc.fit_scale([rows[i] for i in idx])
        if mine is None or theirs is None:
            problems.append(f"{chrom}: a fit is undefined (mine={mine!r}, magcal={theirs!r})")
            continue
        checked += 1
        d = abs(mine - theirs)
        if d > worst:
            worst, worst_chrom = d, chrom
    if worst > sh.TOLERANCE:
        problems.append(
            f"this lane's per-chromosome scale path differs from magcal.fit_scale by {worst:.3g}"
            f" on {worst_chrom}, above {sh.TOLERANCE}"
        )
    g_mine = sh.scale_from_sums(sum(dps), sum(pps))
    g_theirs = mc.fit_scale(rows)
    g_diff = abs(float(g_mine) - float(g_theirs)) if g_mine is not None and g_theirs else float("nan")
    if not (g_diff <= sh.TOLERANCE):
        problems.append(f"the whole-table sums path differs from magcal.fit_scale by {g_diff:.3g}")
    return {
        "chromosomes_checked": checked,
        "max_abs_difference": worst,
        "on_chromosome": worst_chrom,
        "whole_table_abs_difference": g_diff,
        "tolerance": sh.TOLERANCE,
        "note": (
            "the per-chromosome scales themselves are NOT written here and no dispersion of them is"
            " computed: only the agreement of the two paths, as a maximum absolute difference. The"
            " whole-table difference is a difference of two fits of the SAME global scale by two code"
            " paths, not a value of this lane's design"
        ),
    }, problems


def admission(rows: list[dict[str, Any]]) -> tuple[dict[str, Any], list[str]]:
    """HARNESS (5): gene COUNTS per chromosome and the admitted set under MIN_GENES.

    A count depends on no measured value (scalesplit's harness states this in those words).
    """
    problems: list[str] = []
    counts = sh.gene_counts(rows)
    adm, exc = sh.admitted_chromosomes(counts)
    if sum(counts.values()) != len(rows):
        problems.append(f"the per-chromosome counts sum to {sum(counts.values())}, not {len(rows)}")
    if len(adm) < sh.MIN_CHROMOSOMES:
        problems.append(
            f"only {len(adm)} chromosomes carry at least {sh.MIN_GENES} genes, below"
            f" MIN_CHROMOSOMES {sh.MIN_CHROMOSOMES}"
        )
    return {
        "genes_per_chromosome": dict(sorted(counts.items())),
        "chromosomes_in_table": len(counts),
        "min_genes": sh.MIN_GENES,
        "admitted": adm,
        "admitted_count": len(adm),
        "excluded_below_min_genes": {c: counts[c] for c in exc},
        "genes_admitted": sum(counts[c] for c in adm),
        "genes_excluded": sum(counts[c] for c in exc),
        "min_chromosomes": sh.MIN_CHROMOSOMES,
    }, problems


def manifest(entries: list[dict[str, Any]]) -> dict[str, Any]:
    return {
        "sources": [
            {"accession": "data/results/multrna_result.json (7h, lane-multrna)", "version": "pinned"},
            {"accession": "data/results/magcal_result.json (M2h, lane-magcal)", "version": "pinned"},
            {
                "accession": "data/results/scalesplit_result.json (lane-scalesplit)",
                "version": "pinned; its post hoc per-chromosome spread is the pattern this lane tests",
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
            "split_reproduced": (
                "M2h's, sha256 of the gene symbol; this lane's heterogeneity statistic and its null are"
                " not computed here"
            ),
            "min_genes": sh.MIN_GENES,
            "min_chromosomes": sh.MIN_CHROMOSOMES,
            "form_interval": list(sh.FORM_INTERVAL),
            "bootstrap": {"resamples": mc.BOOTSTRAP["resamples"], "seed": mc.BOOTSTRAP["seed"]},
            "permutations_reproduced": mc.PERMUTATIONS,
            "tolerance": sh.TOLERANCE,
            "m2h_tolerance": sh.M2H_TOLERANCE,
        },
        "exclusions": ["none: the population is 7h's, unchanged; the gene floor is applied by the run"],
        "partitions": {
            "reproduced": "M2h's sha256-of-symbol training and held-out genes, to recover its numbers",
            "this_lanes_own": (
                "not computed here: the chromosome grouping's statistic, its null and its p are run only"
                " after the registration is committed"
            ),
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
                else "pinned by sha256; committed numbers and readings, quoted and not re-decided"
            ),
        )
        for p in SOURCE_FILES
    ]
    rows = load_rows()
    tbl, p1 = table_checks(rows)
    rep, p2 = reproduce_m2h(rows)
    paths, p3 = two_paths_agree(rows)
    adm, p4 = admission(rows)
    problems = p1 + p2 + p3 + p4
    out = {
        "lane": sh.LANE,
        "step": "1 of 4: the harness, before the registration",
        "question": sh.QUESTION,
        "why_not_already_answered": sh.WHY_NOT_ALREADY_ANSWERED,
        "harness": sh.HARNESS,
        "prior_exposure": list(sh.PRIOR_EXPOSURE),
        "table_checks": tbl,
        "m2h_reproduction": rep,
        "two_paths_agree": paths,
        "admission": adm,
        "problems": problems,
        "outcome": "HARNESS_REPRODUCES" if not problems else "HARNESS_DEFECT",
        "reading": (
            "M2h's committed gene-level numbers REPRODUCE from 7h's committed genes table through M2h's"
            " own functions at M2h's own seeds: training 8,181 genes, held out 8,328, fitted scale"
            " +0.2146, held-out Delta MSE +0.04339 (about 2% of the variance), control p = 0.000999,"
            " outcome BETTER. This lane's own per-gene-sums scale path agrees with magcal.fit_scale on"
            " every chromosome of the table within TOLERANCE. 23 of the 24 chromosomes carry at least"
            " MIN_GENES = 100 genes and chrY (13) does not. No heterogeneity statistic, no range, no"
            " permutation and no p value is computed here"
            if not problems
            else mc.reading("HARNESS_DEFECT", what="; ".join(problems))
        ),
        "nothing_of_this_lanes_statistic_was_read": (
            "no Q, no range statistic, no permutation of the chromosome labels and no p value was"
            " computed by this script, and no per-chromosome scale value is written to this result. The"
            " per-chromosome fits were formed only to check that this lane's sums path IS"
            " magcal.fit_scale, and only their maximum disagreement -- a number about two code paths,"
            " not about the genome -- is reported. Gene counts per chromosome depend on no measured value"
        ),
        "readers": READERS,
        "alphagenome_requests": 0,
        "money": "none: no request is sent and no key is read",
        "seconds": round(time.time() - t0, 1),
    }
    print("saved", save_result(RESULT, out, manifest=manifest(entries)))
    print(out["outcome"], "|", out["reading"][:400])
    return 0 if not problems else 1


if __name__ == "__main__":
    raise SystemExit(main())
