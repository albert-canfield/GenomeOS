# SPDX-License-Identifier: AGPL-3.0-or-later
"""A pair null with tip noise for area J, calibrated on random gene-set pairs before any pair is read.

GRAMMAR-BY-COMPARISON §22 re-tested §14's library pairs under a Brownian null along Ensembl Compara's
species tree and got no reading: 205 of 860 pairs of random age-matched gene sets had p at or below
0.05, against a registered ceiling of 10%. Its note, which decides nothing, names the likely cause:
the simulated deviations swamp the age expectation (null means -0.14 to +0.04 where the age-expected
r is -0.96 to +0.99). A Brownian rate fitted by REML on Grafen branch lengths, where the tips of a
cherry are 1/354 of the root height apart, reads every difference between two sister species as fast
evolution, so the fitted rate is large and the simulated motion lives in the deep clades. A random
gene set's deviation from its age expectation is mostly species-by-species sampling noise; under that
null it looks like nothing the null draws, and its p is small whenever its age expectation alone makes
the pair alike. This module holds the null §22 named next: the same statistic, the same age
expectation, the deviation drawn with a tip-noise term.

THE NULL FAMILY. Each side's profile P = E + D, E the exact mean of §14's matched draw
(`species_null.age_expectation`), D the side's deviation. The null holds E and replaces D by an
independent draw from a model fitted to D itself:

    D = h * y,   y ~ N(mu 1, sigma^2 V(lambda)),   V(lambda) = lambda R + (1 - lambda) I

R is the tree's correlation matrix (the Brownian covariance on Grafen lengths, which is already a
correlation because the tree is ultrametric with the root at height 1), so V is Pagel's lambda and,
on this tree, exactly Brownian motion plus independent white noise at the tips. lambda in [0, 1],
sigma^2 and mu are fitted by REML for each side on its own deviation, lambda on a grid of 1,001
points. Two members, in a registered order:

1. `lambda_binomial`: h is the binomial sd of the side's origin-stratum make-up at each species,
   sqrt(sum_k n_k p_k (1 - p_k) / sum_k n_k), p_k the share of stratum k's genes present there. A
   share cannot move at a species where every stratum it is made of is all present or all absent
   (104 species for a Eukaryota-born set, 341 for a Catarrhini-born one), and D is exactly zero there
   for random sets and libraries alike; the model draws nothing there and scales the rest to the room
   each share has.
2. `lambda`: h = 1, Pagel's lambda on the deviations as they stand, the null §22 named.

CALIBRATION, ON RANDOM PAIRS ONLY. Every testable pair contributes one random pair of gene sets, each
drawn to its side's origin-stratum make-up (§14's draw); the first set uses §22's own seed and order,
so it is the very set of random pairs §22's Brownian null failed on. A null passes when at most 10% of
the random pairs have p at or below 0.05 and at most 1% have a Benjamini-Hochberg q (over the random
pairs) at or below 0.05. The selected null is the first in the order that passes; it must pass again
on a second, independent set before any library pair is read under it. `scripts/pair_null_register.py`
registers this, `scripts/pair_null_run.py` reads it once, in that order, and stops where it fails.
"""

from __future__ import annotations

import math
from collections.abc import Callable
from dataclasses import dataclass
from typing import Any

import numpy as np

from genomeos.knowledge import species_null as sn

P_SIMS = 50_000  # draws per bank, as §22; the smallest attainable p is 1 / (P_SIMS + 1)
P_SEED = 20261015  # the banks
CAL1_SEED = sn.A_CAL_SEED  # §22's own random pairs
CAL2_SEED = 20261016  # the confirmation set
HALVES_SEED = sn.A_HALVES_SEED  # §22's own halves
LAMBDA_GRID = 1001  # lambda in 0, 0.001, ..., 1
P_ALPHA = 0.05
Q_ALPHA = 0.05
CAL_MAX_SHARE = 0.10  # random pairs at p <= 0.05 above this share: not calibrated (§22's ceiling)
CAL_MAX_Q_SHARE = 0.01  # random pairs at q <= 0.05 above this share: not calibrated in its tail
HALVES_MIN_SHARE = 0.50  # libraries seen against their own halves below this share: low power
CHUNK = 5_000  # draws per block in the null's Pearson
SCALE_EPS = 1e-12  # a species where h is at or below this holds no deviation

MODELS = ("lambda_binomial", "lambda")  # the registered order
MODEL_NAMES = {
    "lambda_binomial": "Pagel's lambda on binomially scaled deviations",
    "lambda": "Pagel's lambda on the deviations as they stand",
}


# ------------------------------------------------------------------------------------------ model


@dataclass(frozen=True)
class Fit:
    """One side's REML fit: lambda, the rate sigma^2, the root mu (on the scaled axis), species used."""

    lam: float
    sigma2: float
    mu: float
    columns: int


def reml_profile(
    y: np.ndarray, evals: np.ndarray, evecs: np.ndarray, lams: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(restricted log likelihood, sigma^2, mu) at each lambda in `lams`, sigma^2 profiled out.

    y ~ N(mu 1, sigma^2 V), V = lambda R + (1 - lambda) I, R = evecs diag(evals) evecs'. The restricted
    log likelihood, up to a constant, is -1/2 [(n - 1) log sigma^2 + log|V| + log(1' V^-1 1)], with
    sigma^2 = r' V^-1 r / (n - 1) and r the residual from the GLS mean. At lambda = 1 sigma^2 is the
    Brownian REML rate (`species_null.Brownian.rate`).
    """
    n = len(y)
    yt = evecs.T @ y
    xt = evecs.T @ np.ones(n)
    lam = np.asarray(lams, dtype=float)[:, None]
    v = lam * evals[None, :] + (1.0 - lam)
    w = 1.0 / v
    xwx = (w * xt * xt).sum(axis=1)
    mu = (w * xt * yt).sum(axis=1) / xwx
    s2 = (w * (yt[None, :] - mu[:, None] * xt[None, :]) ** 2).sum(axis=1) / (n - 1)
    with np.errstate(divide="ignore"):
        ll = -0.5 * ((n - 1) * np.log(s2) + np.log(v).sum(axis=1) + np.log(xwx))
    return ll, s2, mu


def reml_lambda(y: np.ndarray, evals: np.ndarray, evecs: np.ndarray, grid: int = LAMBDA_GRID) -> Fit:
    """The REML fit of lambda (on a grid over [0, 1], the first maximum when two tie), sigma^2 and mu."""
    n = len(y)
    if n < 3:
        return Fit(0.0, 0.0, float(y.mean()) if n else 0.0, n)
    lams = np.linspace(0.0, 1.0, grid)
    ll, s2, mu = reml_profile(y, evals, evecs, lams)
    if not np.any(s2 > 0):
        return Fit(0.0, 0.0, float(mu[0]), n)
    k = int(np.argmax(ll))
    return Fit(float(lams[k]), float(s2[k]), float(mu[k]), n)


def binomial_scale(d: sn.PairInputs, rows: np.ndarray) -> np.ndarray:
    """The binomial sd of a set's origin-stratum make-up at each species (members without an origin
    are left out, as the age expectation leaves them out)."""
    strata = d.origin_of_row[rows]
    strata = strata[strata != ""]
    if not len(strata):
        return np.zeros(len(d.species))
    names, counts = np.unique(strata.astype(str), return_counts=True)
    total = np.zeros(len(d.species))
    for s, c in zip(names, counts, strict=True):
        m = d.stratum_mean[str(s)]
        total += c * m * (1.0 - m)
    return np.sqrt(total / counts.sum())


def scale_of(d: sn.PairInputs, rows: np.ndarray, model: str) -> np.ndarray:
    if model == "lambda_binomial":
        return binomial_scale(d, rows)
    if model == "lambda":
        return np.ones(len(d.species))
    raise ValueError(f"unknown model {model!r}")


@dataclass
class Side:
    """One side under the null: Y = fixed + p * b + q * w, with b a unit-diagonal Brownian draw along
    the tree and w independent unit white noise, both from the side's own bank."""

    fixed: np.ndarray
    p: np.ndarray
    q: np.ndarray
    fit: Fit


class TipNoise:
    """Two independent banks of (Brownian, white) draws on one tree, and the per-side REML fits."""

    def __init__(self, cov: np.ndarray, sims: int = P_SIMS, seed: int = P_SEED):
        sd = np.sqrt(np.diag(cov))
        self.corr = cov / np.outer(sd, sd)
        self.n = cov.shape[0]
        self.sims = sims
        chol = np.linalg.cholesky(self.corr)
        rng = np.random.default_rng(seed)
        banks = []
        for _ in range(2):
            b = (rng.standard_normal((sims, self.n)) @ chol.T).astype(np.float32)
            w = rng.standard_normal((sims, self.n)).astype(np.float32)
            banks.append((b, w))
        (self.b1, self.w1), (self.b2, self.w2) = banks
        self._eig: dict[bytes, tuple[np.ndarray, np.ndarray]] = {}

    def eig(self, mask: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """The eigendecomposition of R restricted to the species in `mask` (cached by mask)."""
        key = np.packbits(mask).tobytes()
        if key not in self._eig:
            idx = np.flatnonzero(mask)
            evals, evecs = np.linalg.eigh(self.corr[np.ix_(idx, idx)])
            self._eig[key] = (np.clip(evals, 1e-12, None), evecs)
        return self._eig[key]

    def side(self, f: np.ndarray, dev: np.ndarray, h: np.ndarray) -> Side:
        """Fit the side's deviation on the species where h is positive; hold f plus the fitted root."""
        mask = h > SCALE_EPS
        if not mask.any():
            z = np.zeros(self.n)
            return Side(f.astype(float), z, z, Fit(0.0, 0.0, 0.0, 0))
        evals, evecs = self.eig(mask)
        fit = reml_lambda(dev[mask] / h[mask], evals, evecs)
        hm = np.where(mask, h, 0.0)
        s = math.sqrt(fit.sigma2)
        return Side(f + hm * fit.mu, s * math.sqrt(fit.lam) * hm, s * math.sqrt(1.0 - fit.lam) * hm, fit)

    def null_r(self, a: Side, b: Side, chunk: int = CHUNK) -> np.ndarray:
        """Pearson r of the two sides' simulated profiles, one value per draw (side a on bank one)."""
        out = np.empty(self.sims)
        for lo in range(0, self.sims, chunk):
            hi = min(lo + chunk, self.sims)
            ya = a.fixed + self.b1[lo:hi] * a.p + self.w1[lo:hi] * a.q
            yb = b.fixed + self.b2[lo:hi] * b.p + self.w2[lo:hi] * b.q
            ya -= ya.mean(axis=1, keepdims=True)
            yb -= yb.mean(axis=1, keepdims=True)
            num = np.einsum("ij,ij->i", ya, yb)
            den = np.sqrt(np.einsum("ij,ij->i", ya, ya) * np.einsum("ij,ij->i", yb, yb))
            with np.errstate(divide="ignore", invalid="ignore"):
                out[lo:hi] = num / den
        return out


# ------------------------------------------------------------------------------------ the test


def pair_test(
    d: sn.PairInputs, tn: TipNoise, model: str, rows_a: np.ndarray, rows_b: np.ndarray
) -> dict[str, Any]:
    """§14's statistic for one pair, and where it falls under the tip-noise null `model`."""
    pa, pb = sn.profile(d, rows_a), sn.profile(d, rows_b)
    ea, eb = sn.age_expectation(d, rows_a), sn.age_expectation(d, rows_b)
    r_obs = sn.pearson(pa - d.genome, pb - d.genome)
    if r_obs is None:
        return {"r_exclusive": None}
    sa = tn.side(ea - d.genome, pa - ea, scale_of(d, rows_a, model))
    sb = tn.side(eb - d.genome, pb - eb, scale_of(d, rows_b, model))
    null = tn.null_r(sa, sb)
    finite = null[np.isfinite(null)]
    exceed = int((finite >= r_obs).sum())
    return {
        "r_exclusive": r_obs,
        "r_age_expected": sn.pearson(ea - d.genome, eb - d.genome),
        "lambda": [sa.fit.lam, sb.fit.lam],
        "rate": [sa.fit.sigma2, sb.fit.sigma2],
        "species_free_to_move": [sa.fit.columns, sb.fit.columns],
        "null_mean": float(finite.mean()) if len(finite) else None,
        "null_sd": float(finite.std(ddof=1)) if len(finite) > 1 else None,
        "draws_at_or_above": exceed,
        "p_pair": (exceed + 1) / (len(null) + 1),
    }


def random_pairs(
    d: sn.PairInputs, pairs: list[tuple], seed: int
) -> list[tuple[str, str, np.ndarray, np.ndarray]]:
    """One random pair of gene sets per testable pair, each drawn to its side's origin-stratum make-up,
    in the order and with the generator `species_null.null_a` uses: §22's seed gives §22's pairs."""
    rng = np.random.default_rng(seed)
    out = []
    for a, b, oa, ob in pairs:
        ra = sn.matched_rows(d, sn.rows_of(d, oa), rng)
        rb = sn.matched_rows(d, sn.rows_of(d, ob), rng)
        out.append((a, b, ra, rb))
    return out


def summarise(ps: list[float]) -> dict[str, Any]:
    """The calibration figures of a list of random-pair p values, and whether they pass."""
    n = len(ps)
    hits = sum(1 for p in ps if p <= P_ALPHA)
    q_hits = sum(1 for q in sn.bh(ps) if q <= Q_ALPHA) if ps else 0
    share = hits / n if n else None
    q_share = q_hits / n if n else None
    return {
        "pairs": n,
        "at_or_below_005": hits,
        "share_at_or_below_005": share,
        "q_at_or_below_005": q_hits,
        "q_share_at_or_below_005": q_share,
        "passes": bool(n and share <= CAL_MAX_SHARE and q_share <= CAL_MAX_Q_SHARE),
    }


Progress = Callable[[str, int, int], None] | None


def calibrate(
    d: sn.PairInputs, tn: TipNoise, rand: list[tuple], model: str, progress: Progress = None
) -> dict[str, Any]:
    """Every random pair under the null `model`; random pairs only, no library pair is read."""
    rows = []
    for k, (a, b, ra, rb) in enumerate(rand):
        t = pair_test(d, tn, model, ra, rb)
        rows.append(
            {
                "for": [a, b],
                "p_pair": t.get("p_pair"),
                "r": t.get("r_exclusive"),
                "r_age_expected": t.get("r_age_expected"),
                "lambda": t.get("lambda"),
            }
        )
        if progress and k % 100 == 0:
            progress(f"calibration {model}", k, len(rand))
    out = summarise([r["p_pair"] for r in rows if r["p_pair"] is not None])
    lams = [x for r in rows for x in (r["lambda"] or [])]
    out["lambda_quartiles"] = [float(v) for v in np.quantile(lams, [0.25, 0.5, 0.75])] if lams else None
    out["rows"] = rows
    return out


def brownian_hits(d: sn.PairInputs, bm: sn.Brownian, rand: list[tuple]) -> int:
    """§22's Brownian null recomputed on a set of random pairs: how many have p at or below 0.05."""
    hits = 0
    for _, _, ra, rb in rand:
        p = sn.pair_test(d, bm, ra, rb).get("p_tree")
        hits += int(p is not None and p <= P_ALPHA)
    return hits


def choose_null(cal1: dict[str, dict[str, Any]]) -> str | None:
    """The first null in the registered order that passes the first calibration set."""
    for model in MODELS:
        if (cal1.get(model) or {}).get("passes"):
            return model
    return None


def reproduce(d: sn.PairInputs, pairs: list[tuple], stored_14: dict[str, Any]) -> dict[str, Any]:
    """§14's r_exclusive recomputed for every testable pair and set against the stored value."""
    stored = {tuple(r["libraries"]): r for r in stored_14["pairs"] if r.get("r_exclusive") is not None}
    mism = 0
    seen = set()
    for a, b, oa, ob in pairs:
        r = sn.pearson(
            sn.profile(d, sn.rows_of(d, oa)) - d.genome, sn.profile(d, sn.rows_of(d, ob)) - d.genome
        )
        if r is None:
            continue
        seen.add((a, b))
        if (a, b) not in stored or abs(r - stored[(a, b)]["r_exclusive"]) > 1e-4:
            mism += 1
    return {
        "n_tested": len(seen),
        "n_tested_14": len(stored),
        "n_mismatch": mism,
        "same_pairs": seen == set(stored),
    }


def read_pairs(
    d: sn.PairInputs, tn: TipNoise, model: str, pairs: list[tuple], progress: Progress = None
) -> dict[str, Any]:
    """Every testable library pair under the selected null, once, with Benjamini-Hochberg q."""
    rows: list[dict[str, Any]] = []
    for k, (a, b, oa, ob) in enumerate(pairs):
        row: dict[str, Any] = {"libraries": [a, b], "exclusive_members": [len(oa), len(ob)]}
        row.update(pair_test(d, tn, model, sn.rows_of(d, oa), sn.rows_of(d, ob)))
        rows.append(row)
        if progress and k % 100 == 0:
            progress("library pairs", k, len(pairs))
    tested = [r for r in rows if r.get("p_pair") is not None]
    for r, q in zip(tested, sn.bh([r["p_pair"] for r in tested]), strict=True):
        r["q_pair"] = q
    return {"pairs": rows, "tested": len(tested)}


def halves(
    d: sn.PairInputs, tn: TipNoise, model: str, sets: dict[str, set[str]], seed: int = HALVES_SEED
) -> dict[str, Any]:
    """Each library against its own random half (§22's split, by its seed): the power check."""
    rng = np.random.default_rng(seed)
    hv = []
    for lib in sorted(sets):
        r = sn.rows_of(d, sets[lib])
        perm = rng.permutation(r)
        h1, h2 = np.sort(perm[: len(r) // 2]), np.sort(perm[len(r) // 2 :])
        t = pair_test(d, tn, model, h1, h2)
        hv.append({"library": lib, "p_pair": t.get("p_pair"), "r": t.get("r_exclusive")})
    ps = [h["p_pair"] for h in hv if h["p_pair"] is not None]
    hits = sum(1 for p in ps if p <= P_ALPHA)
    return {
        "libraries": len(ps),
        "at_or_below_005": hits,
        "share_at_or_below_005": hits / len(ps) if ps else None,
        "rows": hv,
    }


# -------------------------------------------------------------------------------------- outcomes
# The exact wording of every outcome, fixed before the run. The run fills the numbers into the
# sentence its decision rule selects and changes no other word.

P_OUTCOMES = {
    "not_reproduced": "No reading: §14's observed statistic did not reproduce. {n_mismatch} of "
    "{n_tested} testable library pairs differ from the stored r_exclusive by more than 0.0001, or the "
    "testable pairs are not §14's {n_tested_14}, so the pairs this null would judge are not the pairs "
    "§14 judged. No calibration was run and no library pair was read.",
    "null_not_calibrated": "No reading: neither tip-noise null passes the registered calibration (at "
    "most 10% of random age-matched pairs at p at or below 0.05 and at most 1% at q at or below 0.05), "
    "so a library pair's p under either would not be evidence, and no library pair was read under "
    "either. §14's {n_14} of {n_tested} pairs stand as written, with their caveat.",
    "not_confirmed": "No reading: {null_name} passed the registered calibration on the first {n_cal} "
    "random pairs but not on {n_cal2} independent random pairs drawn the same way, where {cal2_hits} "
    "({cal2_share}) have p at or below 0.05 and {cal2_q} have q at or below 0.05, against ceilings of "
    "10% and 1%. No library pair was read. §14's {n_14} of {n_tested} pairs stand as written, with "
    "their caveat.",
    "none_survive": "None of the {n_tested} testable library pairs is more alike than age and shared "
    "ancestry predict under {null_name} along the Ensembl Compara species tree (Grafen branch lengths; "
    "lambda, rate and root fitted by REML for each side): no pair has a Benjamini-Hochberg q at or "
    "below 0.05, and that includes all {n_14} pairs §14 called more alike than age at the same false "
    "discovery rate. The null passed the registered calibration on two independent sets of random "
    "age-matched pairs ({cal1_hits} of {n_cal} and {cal2_hits} of {n_cal2} at p at or below 0.05) and "
    "sees {halves_hits} of {n_halves} libraries against their own random halves at p at or below 0.05. "
    "On this model §14's {n_14} pairs are not distinguished from what age and shared ancestry give, "
    "and no library pair is called a dependency.",
    "none_survive_low_power": "None of the {n_tested} testable library pairs is more alike than age "
    "and shared ancestry predict under {null_name} along the Ensembl Compara species tree (Grafen "
    "branch lengths; lambda, rate and root fitted by REML for each side): no pair has a "
    "Benjamini-Hochberg q at or below 0.05, including the {n_14} pairs §14 called more alike than age. "
    "This is weak evidence against them: the null passed the registered calibration ({cal1_hits} of "
    "{n_cal} and {cal2_hits} of {n_cal2} random age-matched pairs at p at or below 0.05) but sees only "
    "{halves_hits} of {n_halves} libraries against their own random halves at p at or below 0.05, "
    "below the registered half, so the test has little power and a pair's failure does not show that "
    "it is age and shared ancestry. No library pair is called a dependency.",
    "some_survive": "{n_surv} of {n_tested} testable library pairs are more alike than age and shared "
    "ancestry predict under {null_name} along the Ensembl Compara species tree (Grafen branch lengths; "
    "lambda, rate and root fitted by REML for each side; Benjamini-Hochberg q at or below 0.05): "
    "{k_14} of the {n_14} pairs §14 called more alike than age, and {k_new} that §14 did not call. The "
    "other {lost_14} of §14's {n_14} do not survive this null. The null passed the registered "
    "calibration on two independent sets of random age-matched pairs ({cal1_hits} of {n_cal} and "
    "{cal2_hits} of {n_cal2} at p at or below 0.05) and sees {halves_hits} of {n_halves} libraries "
    "against their own random halves at p at or below 0.05. This is agreement of orthologue presence "
    "across species beyond age and this model of shared ancestry: co-occurrence across species is not "
    "a dependency, it does not show that one library acts on the other, and an orthologue's presence "
    "is not function. No library pair is called a dependency.",
}
CAL_SENTENCE = (
    " The first calibration set is §22's own {n_cal} random pairs, on which §22's Brownian null "
    "recomputed here puts {bm_c1_hits} at p at or below 0.05 ({bm_c1_hits_22} in §22). On it, Pagel's "
    "lambda on binomially scaled deviations puts {cal_bin_hits} ({cal_bin_share}) at p at or below "
    "0.05 and {cal_bin_q} at q at or below 0.05, and Pagel's lambda on the deviations as they stand "
    "puts {cal_lam_hits} ({cal_lam_share}) and {cal_lam_q}; the registered order tries the scaled null "
    "first."
)
P_RULE = {
    "1_reproduce": "if the testable pairs are not §14's (same pairs, same count) or any pair's "
    "recomputed r_exclusive differs from the stored one by more than 0.0001, select not_reproduced "
    "and run no calibration",
    "2_calibration": "else compute both nulls on the first calibration set (one random pair of gene "
    "sets per testable pair, each set drawn to its side's origin-stratum make-up with §22's seed and "
    "order, so the random pairs §22 drew); a null passes when at most 10% of these pairs have p at or "
    "below 0.05 and at most 1% have a Benjamini-Hochberg q, over the calibration pairs, at or below "
    "0.05; the selected null is the first that passes in the registered order (lambda_binomial, then "
    "lambda); if neither passes, select null_not_calibrated and read no library pair",
    "3_confirmation": "else compute the selected null on a second calibration set drawn the same way "
    "with its own seed; if it does not pass the same two ceilings there, select not_confirmed and read "
    "no library pair",
    "4_none": "else read every testable library pair once under the selected null, with "
    "Benjamini-Hochberg q over them; if no pair has q at or below 0.05, select none_survive_low_power "
    "when fewer than half of the libraries' random halves have p at or below 0.05, otherwise "
    "none_survive",
    "5_some": "else select some_survive",
    "ties": "'at or below' and 'at most' include equality; 'fewer than half' is strict",
    "always": "the calibration sentence is appended to every outcome except not_reproduced and "
    "decides nothing beyond rules 2 and 3",
}


def share(x: float | None) -> str | None:
    return f"{x:.1%}" if x is not None else None


def select(
    repro: dict[str, Any],
    stored_14: dict[str, Any],
    cal1: dict[str, dict[str, Any]] | None = None,
    bm_c1_hits: int | None = None,
    bm_c1_hits_22: int | None = None,
    chosen: str | None = None,
    cal2: dict[str, Any] | None = None,
    read: dict[str, Any] | None = None,
    hv: dict[str, Any] | None = None,
) -> tuple[str, str, dict[str, Any]]:
    """Apply P_RULE to what the run computed; return (outcome, filled wording, the figures used).

    The run computes in the rule's order and stops where a step fails, so each later argument is
    None exactly when an earlier step decided.
    """
    called = {tuple(r["libraries"]) for r in stored_14["significant_q05"]}
    f: dict[str, Any] = {
        "n_tested": repro["n_tested"],
        "n_tested_14": repro["n_tested_14"],
        "n_mismatch": repro["n_mismatch"],
        "n_14": len(called),
    }
    if repro["n_mismatch"] or not repro["same_pairs"]:
        key = "not_reproduced"
        return key, P_OUTCOMES[key].format(**f), f
    if cal1 is None:
        raise ValueError("the pairs reproduced, so the first calibration set must have been computed")
    cb, cl = cal1["lambda_binomial"], cal1["lambda"]
    f.update(
        {
            "n_cal": cb["pairs"],
            "bm_c1_hits": bm_c1_hits,
            "bm_c1_hits_22": bm_c1_hits_22,
            "cal_bin_hits": cb["at_or_below_005"],
            "cal_bin_share": share(cb["share_at_or_below_005"]),
            "cal_bin_q": cb["q_at_or_below_005"],
            "cal_lam_hits": cl["at_or_below_005"],
            "cal_lam_share": share(cl["share_at_or_below_005"]),
            "cal_lam_q": cl["q_at_or_below_005"],
            "selected": choose_null(cal1),
        }
    )
    if chosen != f["selected"]:
        raise ValueError(f"the run selected {chosen!r}, the rule selects {f['selected']!r}")
    tail = CAL_SENTENCE.format(**f)
    if chosen is None:
        key = "null_not_calibrated"
        return key, P_OUTCOMES[key].format(**f) + tail, f
    c1 = cal1[chosen]
    if cal2 is None:
        raise ValueError("a null passed the first set, so the confirmation set must have been computed")
    f.update(
        {
            "null_name": MODEL_NAMES[chosen],
            "cal1_hits": c1["at_or_below_005"],
            "cal1_q": c1["q_at_or_below_005"],
            "n_cal2": cal2["pairs"],
            "cal2_hits": cal2["at_or_below_005"],
            "cal2_share": share(cal2["share_at_or_below_005"]),
            "cal2_q": cal2["q_at_or_below_005"],
        }
    )
    if not cal2["passes"]:
        key = "not_confirmed"
        return key, P_OUTCOMES[key].format(**f) + tail, f
    if read is None or hv is None:
        raise ValueError("the null was confirmed, so the pairs and the halves must have been read")
    rows = [r for r in read["pairs"] if r.get("p_pair") is not None]
    surv = {tuple(r["libraries"]) for r in rows if r["q_pair"] <= Q_ALPHA}
    f.update(
        {
            "n_surv": len(surv),
            "k_14": len(surv & called),
            "k_new": len(surv - called),
            "lost_14": len(called - surv),
            "halves_hits": hv["at_or_below_005"],
            "n_halves": hv["libraries"],
        }
    )
    if not surv:
        low = (hv.get("share_at_or_below_005") or 0) < HALVES_MIN_SHARE
        key = "none_survive_low_power" if low else "none_survive"
    else:
        key = "some_survive"
    text = P_OUTCOMES[key].format(**f) + tail
    f["survivors"] = sorted(list(p) for p in surv)
    f["survivors_among_14"] = sorted(list(p) for p in surv & called)
    return key, text, f
