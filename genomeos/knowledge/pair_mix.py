# SPDX-License-Identifier: AGPL-3.0-or-later
"""The mixed calibration of area J's tip-noise pair null: one side a real library, the other random.

GRAMMAR-BY-COMPARISON §23 read 154 of 860 library pairs as more alike than age and shared ancestry
predict under Pagel's lambda on binomially scaled deviations (`pair_null`, model `lambda_binomial`),
after the null passed its calibration on two sets of random gene-set pairs. Its note names what that
calibration did not check: random gene sets carry small deviations from their age expectation (their
sides fit lambda 0.65 to 0.77, the library sides 0.84 to 0.90), so the null was never tried on a pair
where a side carries a deviation as large and as tree-shaped as a real library's. A library whose
deviation happens to line up with another stratum make-up's age expectation could read as alike to
almost anything, and the random-pair calibration would not see it.

THE MIXED PAIRS. For every testable library pair (a, b), with a's and b's exclusive members as §14
tested them, two mixed pairs:

1. first library real: a's exclusive members against a random gene set drawn to the origin-stratum
   make-up of b's exclusive members;
2. second library real: a random gene set drawn to the make-up of a's exclusive members against b's
   exclusive members.

So one side carries a real library's deviation and the other a random one, and each library keeps the
side (and so the bank of draws) it had when §23 read the pair. The random set is §14's draw, per
stratum without replacement, from the stratum's genes LESS the library side's members, so the two sides
share no gene, as the library pair's exclusive sides share none. Where a stratum's pool less those
members is smaller than the make-up asks (the Homo, Hominoidea and Haplorrhini strata hold one gene
each), the draw takes what is left and the shortfall is recorded.

Each mixed pair is tested exactly as §23 tested a library pair (`pair_null.pair_test`, the same banks,
seed and fits). The null is calibrated on mixed pairs when, in each direction separately (860 mixed
pairs each, as each of §23's calibration sets was 860), at most 10% have p at or below 0.05 and at most
1% have a Benjamini-Hochberg q over that direction at or below 0.05: §23's two ceilings. Before any
mixed pair is computed the run reproduces §14's statistic and §23's read, pair by pair, so that what is
calibrated here is the null the 154 were read under. `scripts/pair_mix_register.py` registers this and
the consequence for the 154 under each outcome; `scripts/pair_mix_run.py` reads it once.
"""

from __future__ import annotations

from collections import defaultdict
from typing import Any

import numpy as np

from genomeos.knowledge import pair_null as pn
from genomeos.knowledge import species_null as sn

MODEL = "lambda_binomial"  # the null §23 selected and read the library pairs under
MIX_SEED = 20261017  # the mixed draws; §22 used 20261010-14, §23 20261015-16
DIRECTIONS = ("first_library_real", "second_library_real")
DIRECTION_NAMES = {
    "first_library_real": "the first library real, a random set of the second's make-up beside it",
    "second_library_real": "a random set of the first library's make-up, the second library real",
}


# ------------------------------------------------------------------------------------- the draw


def matched_rows_excluding(
    d: sn.PairInputs, makeup_rows: np.ndarray, exclude: np.ndarray, rng: np.random.Generator
) -> tuple[np.ndarray, int]:
    """A random gene set of `makeup_rows`'s origin-stratum make-up drawn from each stratum's genes less
    `exclude` (§14's draw, without replacement per stratum, strata in sorted order). Returns the rows
    and how many genes the make-up asked for that the reduced pools could not give."""
    strata = d.origin_of_row[makeup_rows]
    banned = np.asarray(exclude, dtype=np.int64)
    out = []
    short = 0
    for s in sorted(set(strata.tolist()) - {""}):
        k = int((strata == s).sum())
        pool = d.pools[s]
        if len(banned):
            pool = pool[~np.isin(pool, banned)]
        take = min(k, len(pool))
        short += k - take
        if take:
            out.append(rng.choice(pool, size=take, replace=False))
    rows = np.sort(np.concatenate(out)) if out else np.array([], dtype=np.int64)
    return rows, short


def mixed_pairs(d: sn.PairInputs, pairs: list[tuple], seed: int = MIX_SEED) -> list[dict[str, Any]]:
    """The two mixed pairs of every testable library pair, in testable-pair order. One generator; per
    pair, first the random set beside the first library (the second's make-up), then the random set
    beside the second library (the first's make-up)."""
    rng = np.random.default_rng(seed)
    out = []
    for a, b, oa, ob in pairs:
        ra, rb = sn.rows_of(d, oa), sn.rows_of(d, ob)
        beside_a, short_a = matched_rows_excluding(d, rb, ra, rng)
        beside_b, short_b = matched_rows_excluding(d, ra, rb, rng)
        out.append(
            {
                "libraries": [a, b],
                "first_library_real": (ra, beside_a, short_a),
                "second_library_real": (beside_b, rb, short_b),
            }
        )
    return out


# ------------------------------------------------------------------------------------- the test


def lambda_summary(values: list[float]) -> dict[str, Any] | None:
    """n, min, quartiles and max of a list of fitted lambdas."""
    if not values:
        return None
    q = np.quantile(values, [0.0, 0.25, 0.5, 0.75, 1.0])
    return {
        "n": len(values),
        "min": float(q[0]),
        "q1": float(q[1]),
        "median": float(q[2]),
        "q3": float(q[3]),
        "max": float(q[4]),
    }


def calibrate_mixed(
    d: sn.PairInputs, tn: pn.TipNoise, mixed: list[dict[str, Any]], progress: pn.Progress = None
) -> dict[str, Any]:
    """Every mixed pair under the null §23 read; per direction the calibration figures (with BH over
    that direction), pooled figures, and the fitted lambda of library and random sides."""
    rows: dict[str, list[dict[str, Any]]] = {k: [] for k in DIRECTIONS}
    for k, m in enumerate(mixed):
        for direction in DIRECTIONS:
            s1, s2, short = m[direction]
            t = pn.pair_test(d, tn, MODEL, s1, s2)
            rows[direction].append(
                {
                    "libraries": m["libraries"],
                    "p_pair": t.get("p_pair"),
                    "r": t.get("r_exclusive"),
                    "r_age_expected": t.get("r_age_expected"),
                    "lambda": t.get("lambda"),
                    "random_genes": int(len(s2) if direction == "first_library_real" else len(s1)),
                    "random_genes_short": short,
                }
            )
        if progress and k % 100 == 0:
            progress("mixed pairs", k, len(mixed))
    out: dict[str, Any] = {}
    for direction in DIRECTIONS:
        tested = [r for r in rows[direction] if r["p_pair"] is not None]
        ps = [r["p_pair"] for r in tested]
        out[direction] = pn.summarise(ps)
        for r, q in zip(tested, sn.bh(ps), strict=True):
            r["q_pair"] = q
    out["pooled"] = pn.summarise(
        [r["p_pair"] for k in DIRECTIONS for r in rows[k] if r["p_pair"] is not None]
    )
    lib, rand = [], []
    for direction in DIRECTIONS:
        real = 0 if direction == "first_library_real" else 1
        for r in rows[direction]:
            if r["lambda"]:
                lib.append(r["lambda"][real])
                rand.append(r["lambda"][1 - real])
    out["lambda_library_sides"] = lambda_summary(lib)
    out["lambda_random_sides"] = lambda_summary(rand)
    out["random_genes_short"] = {
        "mixed_pairs_short": sum(1 for k in DIRECTIONS for r in rows[k] if r["random_genes_short"]),
        "genes_short": sum(r["random_genes_short"] for k in DIRECTIONS for r in rows[k]),
    }
    out["rows"] = rows
    return out


def passes(mix: dict[str, Any]) -> bool:
    """The registered verdict: both directions within both ceilings."""
    return all(bool((mix.get(k) or {}).get("passes")) for k in DIRECTIONS)


# ------------------------------------------------------------------------------- reproduction


def reproduce_read(
    d: sn.PairInputs,
    tn: pn.TipNoise,
    pairs: list[tuple],
    stored_read: dict[str, Any],
    progress: pn.Progress = None,
) -> dict[str, Any]:
    """Every library pair's p under the null §23 read, recomputed and set against §23's stored count of
    draws at or above, and the pairs at q at or below 0.05 recomputed against §23's."""
    stored = {tuple(r["libraries"]): r for r in stored_read["pairs"] if r.get("p_pair") is not None}
    got: dict[tuple, float] = {}
    mism = 0
    for k, (a, b, oa, ob) in enumerate(pairs):
        t = pn.pair_test(d, tn, MODEL, sn.rows_of(d, oa), sn.rows_of(d, ob))
        if progress and k % 100 == 0:
            progress("library pairs, reproduced", k, len(pairs))
        if t.get("p_pair") is None:
            continue
        got[(a, b)] = t["p_pair"]
        s = stored.get((a, b))
        if s is None or s["draws_at_or_above"] != t["draws_at_or_above"]:
            mism += 1
    keys = list(got)
    qs = sn.bh([got[k] for k in keys])
    surv = {k for k, q in zip(keys, qs, strict=True) if q <= pn.Q_ALPHA}
    surv_stored = {k for k, r in stored.items() if r.get("q_pair") is not None and r["q_pair"] <= pn.Q_ALPHA}
    return {
        "n_tested": len(got),
        "n_stored": len(stored),
        "n_mismatch_read": mism,
        "same_pairs": set(got) == set(stored),
        "n_survivors_recomputed": len(surv),
        "n_survivors_stored": len(surv_stored),
        "same_survivors": surv == surv_stored,
    }


# ------------------------------------------------------------------------------- descriptive


def per_library(mix: dict[str, Any]) -> list[dict[str, Any]]:
    """For each library, its mixed pairs as the real side and how many have p at or below 0.05."""
    seen: dict[str, list[float]] = defaultdict(list)
    for direction in DIRECTIONS:
        real = 0 if direction == "first_library_real" else 1
        for r in mix["rows"][direction]:
            if r["p_pair"] is not None:
                seen[r["libraries"][real]].append(r["p_pair"])
    return [
        {"library": lib, "mixed_pairs": len(ps), "at_p_005": sum(1 for p in ps if p <= pn.P_ALPHA)}
        for lib, ps in sorted(seen.items())
    ]


def among_read(mix: dict[str, Any], survivors: set[tuple]) -> dict[str, Any]:
    """Of §23's pairs at q at or below 0.05, how many have one or both mixed pairs at p at or below 0.05
    (one random draw per mixed pair: too noisy to judge any one pair, and it decides nothing)."""
    hits: dict[tuple, int] = defaultdict(int)
    for direction in DIRECTIONS:
        for r in mix["rows"][direction]:
            key = tuple(r["libraries"])
            if key in survivors and r["p_pair"] is not None and r["p_pair"] <= pn.P_ALPHA:
                hits[key] += 1
    return {
        "pairs": len(survivors),
        "one_or_both_mixed_at_p_005": sum(1 for k in survivors if hits[k] >= 1),
        "both_mixed_at_p_005": sum(1 for k in survivors if hits[k] == 2),
    }


# -------------------------------------------------------------------------------------- outcomes
# The exact wording of every outcome, fixed before the run. The run fills the numbers into the
# sentence its decision rule selects and changes no other word.

MIX_OUTCOMES = {
    "not_reproduced": "No reading: the null §23 read was not reproduced, so a calibration computed here "
    "would not be a calibration of that null. {n_mismatch_14} of {n_tested} testable library pairs "
    "differ from §14's stored r_exclusive by more than 0.0001 (the testable pairs being §14's: "
    "{same_pairs_14}), {n_mismatch_read} have a count of null draws at or above their r other than the "
    "count §23 stored (the pairs being §23's: {same_pairs_23}), and {n_23_recomputed} pairs, against "
    "§23's {n_23}, have a Benjamini-Hochberg q at or below 0.05 (the same pairs: {same_survivors}). No "
    "mixed pair was computed. §23's reading stands as registered, with its caveat that the null was not "
    "checked when a side carries a library-sized deviation.",
    "calibrated": "The mixed calibration passes, and §23's reading stands as read: its {n_23} of "
    "{n_tested} library pairs ({k_23_14} of §14's {n_14} among them) need no re-reading. Each testable "
    "library pair's two libraries were set, one at a time, against a random gene set drawn to the "
    "partner's origin-stratum make-up and sharing no gene with the library, so one side of each mixed "
    "pair carries a real library's deviation and the other a random one. Under Pagel's lambda on "
    "binomially scaled deviations, the null §23 read, {a_hits} of {n_a} ({a_share}) mixed pairs with the "
    "first library real and {b_hits} of {n_b} ({b_share}) with the second library real have p at or "
    "below 0.05, and {a_q} and {b_q} have a Benjamini-Hochberg q at or below 0.05, within the registered "
    "ceilings of 10% and 1% in each direction. The library sides fit lambda {lam_lib} (quartiles) and "
    "the random sides of the mixed pairs {lam_rand}. This checks the null and reads no pair anew: "
    "co-occurrence across species is not a dependency, and no library pair is called a dependency.",
    "not_calibrated": "No reading: the mixed calibration fails, and §23's {n_23} of {n_tested} library "
    "pairs are withdrawn to no reading. Each testable library pair's two libraries were set, one at a "
    "time, against a random gene set drawn to the partner's origin-stratum make-up and sharing no gene "
    "with the library, so one side of each mixed pair carries a real library's deviation and the other a "
    "random one. Under Pagel's lambda on binomially scaled deviations, the null §23 read, {a_hits} of "
    "{n_a} ({a_share}) mixed pairs with the first library real and {b_hits} of {n_b} ({b_share}) with "
    "the second library real have p at or below 0.05, and {a_q} and {b_q} have a Benjamini-Hochberg q at "
    "or below 0.05, against registered ceilings of 10% and 1% in each direction: the null finds a real "
    "library alike to random genes more often than it should when one side carries a library-sized "
    "deviation, so a library pair's p under it is not evidence. Withdrawn with the {n_23} is §23's "
    "finding that {lost_23_14} of §14's {n_14} do not survive this null; §14's {n_14} of {n_tested} "
    "pairs stand as written, with their caveat. The library sides fit lambda {lam_lib} (quartiles) and "
    "the random sides of the mixed pairs {lam_rand}. No library pair is called a dependency.",
}
POOLED_SENTENCE = (
    " Pooled over the {n_pool} mixed pairs, {pool_hits} ({pool_share}) have p at or below 0.05 and "
    "{pool_q} have a Benjamini-Hochberg q, over all of them, at or below 0.05; the pooled figure decides "
    "nothing."
)
MIX_RULE = {
    "1_reproduce": "recompute §14's r_exclusive for every testable pair (pair_null.reproduce) and every "
    "library pair's p under the null §23 read (Pagel's lambda on binomially scaled deviations, §23's "
    "banks and fits); if the testable pairs are not §14's, any r_exclusive differs from the stored one "
    "by more than 0.0001, the pairs with a p are not §23's, any pair's count of null draws at or above "
    "its r differs from the count §23 stored, or the set of pairs with a Benjamini-Hochberg q at or "
    "below 0.05 is not §23's, select not_reproduced and compute no mixed pair",
    "2_mixed": "else draw the mixed pairs (MIX_SEED, in testable-pair order; per pair, first the random "
    "set beside the first library, drawn to the make-up of the second library's exclusive members, then "
    "the random set beside the second library, drawn to the make-up of the first's; each from the "
    "stratum pools less the library side's members) and test each under the same null with each library "
    "on the side it had in §23's read; the first-library-real and second-library-real mixed pairs are "
    "two sets, each as large as the number of testable pairs",
    "3_verdict": "select calibrated when, in EACH direction, at most 10% of the mixed pairs have p at or "
    "below 0.05 AND at most 1% have a Benjamini-Hochberg q, computed over that direction's mixed pairs, "
    "at or below 0.05; otherwise select not_calibrated",
    "consequence": "calibrated: §23's 154 pairs stand as read and no re-reading is needed; "
    "not_calibrated: §23's 154 pairs are withdrawn to no reading, together with §23's finding about "
    "§14's 35, and §14's 35 stand as written with their caveat; not_reproduced: nothing is decided "
    "about §23, whose reading stands as registered with its caveat",
    "ties": "'at or below' and 'at most' include equality",
    "always": "the pooled sentence is appended to calibrated and not_calibrated; it, every lambda "
    "figure, the shortfall of the random draws and every per-library or per-pair figure are "
    "descriptive and decide nothing",
}


def quartile_range(s: dict[str, Any] | None) -> str | None:
    return f"{s['q1']:.3f} to {s['q3']:.3f}" if s else None


def select_mix(
    repro14: dict[str, Any],
    repro23: dict[str, Any],
    stored_14: dict[str, Any],
    stored_23_figures: dict[str, Any],
    mix: dict[str, Any] | None = None,
) -> tuple[str, str, dict[str, Any]]:
    """Apply MIX_RULE to what the run computed; return (outcome, filled wording, the figures used).

    `stored_23_figures` is §23's own `figures` block: n_surv, k_14 and lost_14 as §23 recorded them.
    `mix` is None exactly when rule 1 decided.
    """
    f: dict[str, Any] = {
        "n_tested": repro14["n_tested"],
        "n_mismatch_14": repro14["n_mismatch"],
        "same_pairs_14": "yes" if repro14["same_pairs"] else "no",
        "n_mismatch_read": repro23["n_mismatch_read"],
        "same_pairs_23": "yes" if repro23["same_pairs"] else "no",
        "n_23_recomputed": repro23["n_survivors_recomputed"],
        "same_survivors": "yes" if repro23["same_survivors"] else "no",
        "n_23": stored_23_figures["n_surv"],
        "k_23_14": stored_23_figures["k_14"],
        "lost_23_14": stored_23_figures["lost_14"],
        "n_14": len(stored_14["significant_q05"]),
    }
    reproduced = (
        not repro14["n_mismatch"]
        and repro14["same_pairs"]
        and not repro23["n_mismatch_read"]
        and repro23["same_pairs"]
        and repro23["same_survivors"]
    )
    if not reproduced:
        if mix is not None:
            raise ValueError("rule 1 failed, so no mixed pair may have been computed")
        key = "not_reproduced"
        return key, MIX_OUTCOMES[key].format(**f), f
    if mix is None:
        raise ValueError("the read reproduced, so the mixed pairs must have been computed")
    a, b, pool = mix["first_library_real"], mix["second_library_real"], mix["pooled"]
    f.update(
        {
            "n_a": a["pairs"],
            "a_hits": a["at_or_below_005"],
            "a_share": pn.share(a["share_at_or_below_005"]),
            "a_q": a["q_at_or_below_005"],
            "n_b": b["pairs"],
            "b_hits": b["at_or_below_005"],
            "b_share": pn.share(b["share_at_or_below_005"]),
            "b_q": b["q_at_or_below_005"],
            "lam_lib": quartile_range(mix["lambda_library_sides"]),
            "lam_rand": quartile_range(mix["lambda_random_sides"]),
            "n_pool": pool["pairs"],
            "pool_hits": pool["at_or_below_005"],
            "pool_share": pn.share(pool["share_at_or_below_005"]),
            "pool_q": pool["q_at_or_below_005"],
        }
    )
    key = "calibrated" if passes(mix) else "not_calibrated"
    return key, MIX_OUTCOMES[key].format(**f) + POOLED_SENTENCE.format(**f), f
