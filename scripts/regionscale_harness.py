# SPDX-License-Identifier: AGPL-3.0-or-later
"""lane-regionscale, step 1: the harness, run and committed BEFORE the registration.

    uv run python scripts/regionscale_harness.py

It answers the questions regionscale.HARNESS (1) to (6) names and no others.

(A) FIRST, AND REPORTED BEFORE ANYTHING ELSE: reading 7h's committed genes table through M2h's OWN
functions, at M2h's own resamples, seeds and permutation count, DO M2h'S COMMITTED NUMBERS COME BACK?
lane-scalesplit and lane-scalehet each did this bit-identically and this follows them.

(B) Is this lane's one-scale-per-gene metric path the SAME path as magcal's? With every gene given the
same scale, `delta_mse_rowwise` must equal `magcal.delta_mse` within TOLERANCE, and the paired
bootstrap must return a DELTA interval of EXACTLY [0, 0] while its per-arm interval equals
`magcal.bootstrap_delta`'s within TOLERANCE. The per-arm agreement is not bit-for-bit and the result
says so and says why: the two accumulate a gene's term in a different ORDER, which is a floating-point
difference of order 1e-17 and not an algebraic one.

(C) Is there an HONEST HELD-OUT DESIGN at these gene counts at all? The gene counts per chromosome per
half, the admitted set under ADMISSION, and whether it clears MIN_REGIONS. A COUNT DEPENDS ON NO
MEASURED VALUE (regionscale.GENE_COUNTS_ARE_NOT_MEASUREMENTS), which is why this is computed here.

(D) The variance denominator of "about X% of the variance": the no-rule baseline's own held-out MSE,
which is the denominator M2h's committed "about 2%" is already stated against.

NO VALUE OF THIS LANE'S DESIGN IS COMPUTED HERE: no per-chromosome scale is fitted, no per-region
Delta MSE, no DELTA, no bootstrap of DELTA and no permutation. The constant scale used for the path
reduction is M2h's OWN committed scale, not a value of this design. That is regionscale.HARNESS (7),
and the registered statistic is the thing that must not be seen before it is registered.

Nothing of magcal's, scalehet's, hetblock's, scalesplit's or multrna's is edited or re-run. No model
request, no network, no money, 0 AlphaGenome requests.
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
from genomeos.attribution import regionscale as rs  # noqa: E402
from genomeos.results import save_result  # noqa: E402

RESULT = "regionscale_harness"
#: everything read, repository-relative (d72f959)
MR_RESULT = "data/results/multrna_result.json"
MAGCAL_RESULT = "data/results/magcal_result.json"
MAGCAL_REGISTRATION = "data/results/magcal_registration.json"
SCALEHET_RESULT = "data/results/scalehet_result.json"
HETBLOCK_RESULT = "data/results/hetblock_result.json"
SOURCE_FILES = (
    MR_RESULT,
    MAGCAL_RESULT,
    MAGCAL_REGISTRATION,
    SCALEHET_RESULT,
    HETBLOCK_RESULT,
)
OWN_CODE = (
    "genomeos/attribution/regionscale.py",
    "scripts/regionscale_harness.py",
    "scripts/regionscale_register.py",
    "scripts/regionscale_run.py",
    "tests/test_regionscale.py",
)
#: every reader this run uses: the thirteen reader cells of 7d, 7h and M2h, declared (rule f)
READERS = list(rs.CELLS)


def load_rows() -> list[dict[str, Any]]:
    """7h's committed genes table as magcal's rows: {gene, symbol, chrom, n_rules, n_inhibiting, P, D}.

    The same columns scripts/magcal_run.py:load_rows and scripts/scalehet_harness.py:load_rows read,
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
    want = rs.SEVEN_H_COMMITTED
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
        "tolerance": rs.TOLERANCE,
    }
    if len(rows) != want["genes"]:
        problems.append(f"{len(rows)} genes, 7h committed {want['genes']}")
    if pairs != want["pairs"]:
        problems.append(f"{pairs} pairs, 7h committed {want['pairs']}")
    if abs(got - want["delta_mse"]) > rs.TOLERANCE:
        problems.append(
            f"recomputed Delta MSE at scale 1 {got:.6f} differs from 7h's committed"
            f" {want['delta_mse']:.6f} by more than {rs.TOLERANCE}"
        )
    if any(len(r["P"]) != len(rs.CELLS) or len(r["D"]) != len(rs.CELLS) for r in rows):
        problems.append(f"a row does not have {len(rs.CELLS)} P and D values")
    if any(not r["symbol"] for r in rows):
        problems.append("a row has an empty symbol")
    if any(not r["chrom"] for r in rows):
        problems.append("a row has an empty chromosome")
    syms = [r["symbol"] for r in rows]
    if len(set(syms)) != len(syms):
        problems.append(f"{len(syms) - len(set(syms))} symbols occur more than once")
    return out, problems


def reproduce_m2h(rows: list[dict[str, Any]]) -> tuple[dict[str, Any], list[str]]:
    """HARNESS (3): M2h's committed numbers, recomputed through M2h's OWN functions. FIRST.

    M2h's split, M2h's resamples, M2h's seeds, M2h's permutation count. This is the check that makes
    every later difference attributable to the REGIONAL UNIT and to nothing else, and the task of this
    lane requires it to be reported before anything else is.
    """
    problems: list[str] = []
    want = rs.M2H_COMMITTED
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
        if abs(float(got[k]) - float(want[k])) > rs.M2H_TOLERANCE:
            problems.append(f"M2h {k}: recomputed {got[k]!r}, committed {want[k]!r}")
    for k in ("scale_ci", "held_out_ci"):
        for g, w in zip(got[k], want[k], strict=True):
            if abs(g - w) > rs.M2H_TOLERANCE:
                problems.append(f"M2h {k}: recomputed {got[k]!r}, committed {want[k]!r}")
    return {
        "split_rule": "M2h's, sha256 of the gene symbol (magcal.SPLIT); carried UNCHANGED by this lane",
        "m2h_tolerance": rs.M2H_TOLERANCE,
        "recomputed": got,
        "committed": {k: want[k] for k in got},
        "reproduces": not problems,
        "reading": (
            "M2h'S COMMITTED NUMBERS REPRODUCE from 7h's committed genes table through M2h's own"
            " functions at M2h's own seeds: training 8,181 genes, held out 8,328, fitted scale"
            " +0.2146 [98.33% CI +0.2007, +0.2283], held-out Delta MSE +0.04339 [+0.03790, +0.04894]"
            " (about 2% of the variance), control p = 0.000999, outcome BETTER"
            if not problems
            else "M2h's committed numbers DO NOT come back: " + "; ".join(problems)
        ),
    }, problems


def one_path(rows: list[dict[str, Any]]) -> tuple[dict[str, Any], list[str]]:
    """HARNESS (4): this lane's one-scale-per-gene path IS magcal's, in metric and in interval.

    The constant scale is M2h's OWN committed scale, so no value of this lane's design is formed. No
    per-chromosome scale is fitted here.
    """
    problems: list[str] = []
    a = rs.GLOBAL_SCALE_COMMITTED
    _train, held = mc.split_rows(rows)
    dps, pps, ns = rs.gene_sums(held)
    const = [a] * len(held)
    mine = rs.delta_mse_rowwise(dps, pps, ns, const)
    theirs = mc.delta_mse(held, a)
    d_metric = abs(mine - theirs)
    if d_metric > rs.TOLERANCE:
        problems.append(
            f"delta_mse_rowwise at a constant scale differs from magcal.delta_mse by {d_metric:.3g},"
            f" above {rs.TOLERANCE}"
        )
    boot = rs.bootstrap_paired(dps, pps, ns, const, const, rs.BOOTSTRAP["resamples"], rs.BOOTSTRAP["seed"])
    mc_lo, mc_hi = mc.bootstrap_delta(held, a, mc.BOOTSTRAP["resamples"], mc.BOOTSTRAP["seed"])
    arm_lo, arm_hi = boot["region_arm_ci"]
    d_arm = max(abs(arm_lo - mc_lo), abs(arm_hi - mc_hi))
    if d_arm > rs.TOLERANCE:
        problems.append(
            f"the paired bootstrap's per-arm interval differs from magcal.bootstrap_delta's by"
            f" {d_arm:.3g}, above {rs.TOLERANCE}: mine [{arm_lo!r}, {arm_hi!r}],"
            f" magcal [{mc_lo!r}, {mc_hi!r}]"
        )
    if boot["delta_ci"] != [0.0, 0.0] or boot["delta_ci_nominal"] != [0.0, 0.0]:
        problems.append(
            f"two identical arms did not give a DELTA interval of exactly [0, 0]:"
            f" {boot['delta_ci']} and {boot['delta_ci_nominal']}"
        )
    # the two arms' intervals must also agree with each other, since both arms are the same vector
    if boot["region_arm_ci"] != boot["global_arm_ci"]:
        problems.append("two identical arms gave two different per-arm intervals")
    return {
        "constant_scale_used": a,
        "whose_scale_it_is": (
            "M2h's own committed fitted scale, quoted; NOT a per-region value and NOT a value of this"
            " lane's design"
        ),
        "metric_mine": mine,
        "metric_magcal": theirs,
        "metric_abs_difference": d_metric,
        "tolerance": rs.TOLERANCE,
        "arm_interval_mine": [arm_lo, arm_hi],
        "arm_interval_magcal": [mc_lo, mc_hi],
        "arm_interval_max_abs_difference": d_arm,
        "delta_interval_of_two_identical_arms": boot["delta_ci"],
        "delta_interval_nominal_of_two_identical_arms": boot["delta_ci_nominal"],
        "note": (
            "with every gene given one and the same scale, this lane's metric IS magcal.delta_mse and"
            " this lane's per-arm bootstrap IS magcal.bootstrap_delta WITHIN TOLERANCE, because the"
            " draw sequence is magcal's own. It is not BIT-for-bit, and the reason is recorded rather"
            " than smoothed over: magcal accumulates a gene's term cell by cell while this lane"
            " accumulates sum(P D) and sum(P^2) first -- its permutation null must re-weight those two"
            " group sums at a thousand scales -- so the two differ in floating-point summation ORDER"
            " and not algebraically, by the number in"
            " arm_interval_max_abs_difference, which is of order 1e-17 against TOLERANCE = 1e-4. Two"
            " identical arms nonetheless give a DELTA interval of EXACTLY [0, 0], because the paired"
            " draw cancels the two arms' terms index for index; that is the property which makes a"
            " non-zero DELTA interval attributable to the scales and not to the resampling"
        ),
    }, problems


def design_feasibility(rows: list[dict[str, Any]]) -> tuple[dict[str, Any], list[str]]:
    """HARNESS (5): the gene COUNTS, the admitted set, and whether an honest design exists at all.

    A COUNT DEPENDS ON NO MEASURED VALUE (regionscale.GENE_COUNTS_ARE_NOT_MEASUREMENTS); this is the
    one thing the harness is permitted to establish about the design, and it is the thing that decides
    whether the lane can proceed or must report "a per-region scale cannot be tested out of sample at
    this gene count" and stop.
    """
    problems: list[str] = []
    counts = rs.split_counts(rows)
    adm, fall = rs.admitted_regions(counts)
    testable = rs.is_testable(adm)
    tr_tot = sum(d["train"] for d in counts.values())
    ho_tot = sum(d["held"] for d in counts.values())
    if tr_tot + ho_tot != len(rows):
        problems.append(f"the per-chromosome half-counts sum to {tr_tot + ho_tot}, not {len(rows)}")
    if tr_tot != rs.M2H_COMMITTED["training_genes"]:
        problems.append(f"{tr_tot} training genes, M2h committed {rs.M2H_COMMITTED['training_genes']}")
    if ho_tot != rs.M2H_COMMITTED["held_out_genes"]:
        problems.append(f"{ho_tot} held-out genes, M2h committed {rs.M2H_COMMITTED['held_out_genes']}")
    adm_tr = sorted(counts[c]["train"] for c in adm)
    return {
        "genes_per_chromosome_per_half": {c: counts[c] for c in sorted(counts)},
        "chromosomes_in_table": len(counts),
        "min_genes_each_half": rs.MIN_GENES,
        "min_regions": rs.MIN_REGIONS,
        "admitted": adm,
        "admitted_count": len(adm),
        "fallback_to_global_both_arms": {c: counts[c] for c in fall},
        "fallback_count": len(fall),
        "held_out_genes_in_admitted": sum(counts[c]["held"] for c in adm),
        "held_out_genes_in_fallback": sum(counts[c]["held"] for c in fall),
        "training_genes_in_admitted": sum(counts[c]["train"] for c in adm),
        "admitted_training_genes_min": adm_tr[0] if adm_tr else None,
        "admitted_training_genes_median": adm_tr[len(adm_tr) // 2] if adm_tr else None,
        "admitted_training_genes_max": adm_tr[-1] if adm_tr else None,
        "the_held_out_design_is_testable": testable,
        "what_that_means": (
            "the design CAN be made honest at these gene counts: a region's scale is fitted on its"
            " training genes alone and scored on its held-out genes, and enough regions clear both"
            " floors for the comparison to be made"
            if testable
            else "the design CANNOT be made honest at these gene counts, and the registered finding is"
            " 'a per-region scale cannot be tested out of sample at this gene count'"
        ),
        "a_count_depends_on_no_measured_value": rs.GENE_COUNTS_ARE_NOT_MEASUREMENTS,
    }, problems


def variance_denominator(rows: list[dict[str, Any]]) -> tuple[dict[str, Any], list[str]]:
    """HARNESS (6): the denominator of "about X% of the variance", and M2h's own "about 2%" confirmed.

    The no-rule baseline's held-out MSE is the mean of D^2 over the held-out pairs, because D is the
    measured deviation FROM that baseline. M2h's committed reading already states its Delta MSE as
    about 2% of the variance, so this is a re-derivation of a quantity already published and not a
    newly exposed one.
    """
    problems: list[str] = []
    _train, held = mc.split_rows(rows)
    dd = 0.0
    n = 0
    for r in held:
        for d in r["D"]:
            dd += d * d
            n += 1
    var = dd / n if n else float("nan")
    pct = 100.0 * rs.M2H_COMMITTED["held_out_delta_mse"] / var if var else float("nan")
    if not (1.0 <= pct <= 3.0):
        problems.append(
            f"M2h's committed held-out Delta MSE is {pct:.3f}% of the held-out variance, which is not"
            " the 'about 2%' its committed reading states"
        )
    return {
        "held_out_pairs": n,
        "held_out_variance": var,
        "what_it_is": (
            "the no-rule baseline's own held-out MSE, the mean of D^2 over the held-out gene-cell"
            " pairs; D is the measured deviation from that baseline, so this is the denominator every"
            " 'about X% of the variance' in this lane is taken against"
        ),
        "m2h_held_out_delta_mse": rs.M2H_COMMITTED["held_out_delta_mse"],
        "m2h_percent_of_variance": pct,
        "m2h_committed_wording": "about 2% of the variance",
        "confirms": not problems,
    }, problems


def manifest(entries: list[dict[str, Any]]) -> dict[str, Any]:
    return {
        "sources": [
            {"accession": "data/results/multrna_result.json (7h, lane-multrna)", "version": "pinned"},
            {"accession": "data/results/magcal_result.json (M2h, lane-magcal)", "version": "pinned"},
            {
                "accession": "data/results/scalehet_result.json (lane-scalehet)",
                "version": "pinned; its rejection of a single global scale is why this lane exists",
            },
            {
                "accession": "data/results/hetblock_result.json (lane-hetblock2)",
                "version": "pinned; its locus-block re-test of that rejection, carried word for word",
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
            "form": rs.FORM,
            "split_rule": "magcal.SPLIT, sha256 of the gene symbol, carried UNCHANGED",
            "regional_unit": "chromosome, inherited from scalehet and hetblock2",
            "min_genes_each_half": rs.MIN_GENES,
            "min_regions": rs.MIN_REGIONS,
            "form_interval": list(rs.FORM_INTERVAL),
            "nominal_interval": list(rs.NOMINAL_INTERVAL),
            "bootstrap": {
                "resamples": rs.BOOTSTRAP["resamples"],
                "seed": rs.BOOTSTRAP["seed"],
            },
            "permutations_not_run_here": rs.PERMUTATIONS,
            "tolerance": rs.TOLERANCE,
            "m2h_tolerance": rs.M2H_TOLERANCE,
        },
        "exclusions": [
            "none: the population is 7h's, unchanged. The ADMISSION floors are applied by the run, and"
            " a chromosome that does not clear them is NOT excluded -- its held-out genes stay in the"
            " scored set and are scored by the global scale under both arms"
        ],
        "partitions": {
            "reproduced": "M2h's sha256-of-symbol training and held-out genes, to recover its numbers",
            "this_lanes_own": (
                "the same partition, taken WITHIN each chromosome: a region's scale is fitted on its"
                " training genes and scored on its held-out genes. No scale of this design is fitted"
                " here and no held-out value under it is read"
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
    rep, p_rep = reproduce_m2h(rows)
    tbl, p_tbl = table_checks(rows)
    path, p_path = one_path(rows)
    feas, p_feas = design_feasibility(rows)
    var, p_var = variance_denominator(rows)
    problems = p_rep + p_tbl + p_path + p_feas + p_var
    ok = not problems
    reading = (
        rep["reading"]
        + ". This lane's one-scale-per-gene metric IS magcal.delta_mse at a constant scale, differing by"
        f" {path['metric_abs_difference']:.3g} against TOLERANCE {rs.TOLERANCE}, and its per-arm"
        " bootstrap IS magcal.bootstrap_delta, differing by"
        f" {path['arm_interval_max_abs_difference']:.3g} -- the REGISTERED claim is agreement within"
        " TOLERANCE, because the two accumulate a gene's term in a different ORDER and that costs about"
        " 1e-17 on a small table; on THIS table the observed difference is exactly 0. Two identical arms"
        " give a DELTA interval of exactly [0, 0]. "
        + (
            f"{feas['admitted_count']} of the {feas['chromosomes_in_table']} chromosomes carry at least"
            f" MIN_GENES = {rs.MIN_GENES} TRAINING genes AND at least {rs.MIN_GENES} HELD-OUT genes,"
            f" above MIN_REGIONS = {rs.MIN_REGIONS}, SO AN HONEST HELD-OUT DESIGN FOR A PER-REGION"
            " SCALE EXISTS at these gene counts and the lane may proceed to register it; the"
            f" remaining {feas['fallback_count']}"
            f" ({', '.join(sorted(feas['fallback_to_global_both_arms'])) or 'none'}) keep their"
            f" {feas['held_out_genes_in_fallback']} held-out genes in the scored set and are scored by"
            " the global scale under BOTH arms"
            if feas["the_held_out_design_is_testable"]
            else f"only {feas['admitted_count']} chromosomes carry MIN_GENES in BOTH halves, below"
            f" MIN_REGIONS = {rs.MIN_REGIONS}, SO A PER-REGION SCALE CANNOT BE TESTED OUT OF SAMPLE AT"
            " THIS GENE COUNT and that is the finding"
        )
        + f". M2h's committed held-out Delta MSE is {var['m2h_percent_of_variance']:.2f}% of the"
        " held-out variance, confirming its committed 'about 2% of the variance'. No per-chromosome"
        " scale, no per-region Delta MSE, no DELTA, no bootstrap of DELTA and no permutation was"
        " computed here"
    )
    out = {
        "lane": rs.LANE,
        "step": "1 of 4: the harness, before the registration",
        "question": rs.QUESTION,
        "why_not_already_answered": rs.WHY_NOT_ALREADY_ANSWERED,
        "harness": rs.HARNESS,
        "prior_exposure": list(rs.PRIOR_EXPOSURE),
        "m2h_reproduction_FIRST": rep,
        "table_checks": tbl,
        "one_path": path,
        "design_feasibility": feas,
        "variance_denominator": var,
        "problems": problems,
        "outcome": "HARNESS_REPRODUCES" if ok else "HARNESS_DEFECT",
        "reading": reading if ok else rs.reading("HARNESS_DEFECT", what="; ".join(problems)),
        "nothing_of_this_lanes_design_was_read": (
            "no per-chromosome scale was fitted, no per-region Delta MSE was formed, no DELTA, no"
            " paired bootstrap of DELTA and no permutation of the chromosome labels was computed, and"
            " no value of this lane's design is written to this result. The constant scale used to"
            " check that this lane's metric IS magcal's is M2h's OWN committed scale. The gene counts"
            " depend on no measured value, and the variance denominator is the quantity M2h's"
            " committed 'about 2%' is already stated against"
        ),
        "readers": READERS,
        "alphagenome_requests": 0,
        "money": "none: no request is sent and no key is read",
        "seconds": round(time.time() - t0, 1),
    }
    print("saved", save_result(RESULT, out, manifest=manifest(entries)))
    print(out["outcome"], "|", out["reading"][:700])
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
