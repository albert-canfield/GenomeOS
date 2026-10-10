# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""Gate (b) asked again, of the per-demander `saturable` policy, registered before any run.

Stage 2's gate (b) failed algebraically (BIOLANG-v0.4-ECONOMY.md §9.2): `proportional` and
`competitive` hand every demander of a pool one factor, so a baseline that fits a constant absorbs
the pool layer exactly (log10 MAE 0.952005 against 0.952005 on 15,159 genes). The record named the
one change that could make the gate answerable, a share that saturates in each demander's OWN demand,
and sized it: measured protein rises 0.515 [0.499, 0.531] decades per decade of transcript, so such a
policy would have to supply 0.485 decades per decade of compression. 91e2a85 built that policy
(`saturable`, K = `half_saturation`). This script asks gate (b)'s question of it, on gate (b)'s join,
with its one parameter K fitted on a training half of the genes chosen by a hash rule fixed here, and
its error read only on the other half.

    uv run python scripts/saturable_gate.py register   # writes the registration; reads no protein value
    uv run python scripts/saturable_gate.py dry        # the real transcript arm with a synthetic protein
                                                       # arm: exercises the run's code, writes nothing
    uv run python scripts/saturable_gate.py run        # once, after the registration is committed

`register` writes data/results/saturable_gate_registration.json and `run` writes
data/results/saturable_gate_result.json. `run` refuses unless the committed registration's binding
words equal REGISTRATION below and every file it pins has the sha256 it pinned.
"""

from __future__ import annotations

import hashlib
import json
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from genomeos.ir import Allocation, Cost, Pool  # noqa: E402
from genomeos.runtime.abundance_gate import (  # noqa: E402
    ELONGATION_AA_PER_S,
    MEASURED_TOTAL_PROTEINS,
    MILO_BAND,
    MIN_COPIES,
    MIN_IMPROVEMENT_LOG10,
    PROTEIN_HALF_LIFE_H,
    RIBOSOME_CAPACITY,
    TRANSCRIPT_SOURCE,
)
from genomeos.runtime.economy import Economy  # noqa: E402

REGISTRATION_NAME = "saturable_gate_registration"
RESULT_NAME = "saturable_gate_result"
REGISTRATION_FILE = ROOT / "data" / "results" / f"{REGISTRATION_NAME}.json"
PAXDB_FILE = ROOT / "data" / "knowledge" / "paxdb" / "9606-WHOLE_ORGANISM-integrated.txt"
GTEX_FILE = ROOT / "data" / "knowledge" / "gtex" / "median_tpm_distilled.tsv"
PROTEINS = ROOT / "data" / "knowledge" / "proteins"

#: the code this lane answers for (the cleanliness block's own code)
OWN_CODE = [
    "scripts/saturable_gate.py",
    "tests/test_saturable_gate.py",
    "genomeos/runtime/economy.py",
    "genomeos/lang/parser.py",
    "genomeos/lang/grammar.py",
    "genomeos/ir/model.py",
]
#: what the registration pins by sha256 and `run` refuses to run past a change in
PINNED_CODE = ["scripts/saturable_gate.py", "genomeos/runtime/economy.py"]

SEED = 20261010
BOOTSTRAP = 2000
#: log10 K from -1.00 to 7.00 in steps of 0.01 (TPM units): from saturating nearly every gene to
#: saturating none, since 1e7 TPM is ten times the whole transcriptome
LOG10_K_GRID = tuple(round(-1.0 + 0.01 * i, 2) for i in range(801))
#: the comparator's slope, 0.000 to 1.500 in steps of 0.001
SLOPE_GRID = tuple(round(0.001 * j, 3) for j in range(1501))
#: gate (b)'s registered slope of log10 copies on log10 TPM, as the record states it (§9.2)
REGISTERED_SLOPE = 0.515
REGISTERED_SLOPE_95 = (0.499, 0.531)
#: gate (b) drove the pools 1e4 times harder before anything was short; the gene-blind arms do the same
GENE_BLIND_DEMAND_MULTIPLIER = 1.0e4
CLOSED_FORM_TOLERANCE_LOG10 = 1e-9
GENE_BLIND_TOLERANCE_LOG10 = 1e-9
MIN_HELD_OUT = 1000

PLANTED_N = 4000
PLANTED_K = 30.0
PLANTED_SD = 0.15
PLANTED_SHIFT = 3.0

REGISTRATION: dict[str, Any] = {
    "question": (
        "does the per-demander `saturable` allocation policy (91e2a85), with its one parameter K fitted"
        " on training genes only, predict absolute protein copies from transcript abundance on held-out"
        " genes better than gate (b)'s one-to-one baseline fitted on the same training genes?"
    ),
    "form": (
        "prediction_i = log10 TPM_i + log10 f_i + c, where f_i is the factor `Economy.scale` returns"
        " for gene i under a program with ONE pool, Ribosomes, size 5e6 (§5.1, the value gate (b)"
        " used), each gene charging 1 per chain with its work rate equal to its TPM (gate (b)'s"
        " convention), and `allocation { pool: Ribosomes; policy: saturable; half_saturation: K }`,"
        " parsed from .bio text. So f_i = (K / (K + TPM_i)) * cut, cut = min(1, 5e6 / sum_j TPM_j K /"
        " (K + TPM_j)) over every analysed gene. The ATP pool gate (b) also declared is left out: a"
        " gene-blind pool in the per-gene minimum could only clip, and gate (b) found neither pool"
        " binding at real demand"
    ),
    "parameters": {
        "K": (
            "FITTED on the training genes only: the value on LOG10_K_GRID (log10 K from -1.00 to 7.00"
            " in 0.01 steps, TPM units) minimising training log10 MAE, searched from the top so a tie"
            " keeps the larger K (closer to no saturation). No independent measurement supplies a"
            " per-gene half-saturation demand in TPM units, so none is cited"
        ),
        "c": "FITTED on the training genes only: the median training residual, the L1-optimal constant",
        "capacity": "5e6 ribosomes, §5.1 (HeLa, BioNumbers), as gate (b); not fitted",
        "baseline_c": "the median training residual of log10 copies minus log10 TPM; gate (b)'s baseline",
    },
    "data": (
        "gate (b)'s join, rebuilt by gate (b)'s own readers (scripts/abundance_gate.py: paxdb, gtex,"
        " lengths): PaxDb integrated human whole-organism ppm times 5e9 proteins per cell, against GTEx"
        " v8 median TPM in Cells - Cultured fibroblasts, keeping a gene with at least 10 copies, TPM"
        " above 0 and a residue count (15,159 genes in gate (b)). The copy floor reads every gene's"
        " protein value, held-out genes included, as gate (b)'s denominator did; it selects the"
        " population and fits nothing"
    ),
    "split": (
        "held out when int(sha256(symbol as UTF-8).hexdigest()[:8], 16) is odd, training when even;"
        " fixed here, before any protein value of either half was read, and blind to every value"
    ),
    "circularity_guard": (
        "no protein value of a held-out gene enters any fit: K, c, the baseline's constant, the"
        " comparator's slope and intercept, the declared-capacity arm's constant and the gene-blind"
        " arms' constants are all fitted on training genes alone. Held-out TPM values do enter the"
        " pool's served total (the cut), because the program serves every gene from one pool; that is"
        " a transcript value, not an outcome. tests/test_saturable_gate.py replaces every held-out"
        " protein value and asserts every fitted number is unchanged"
    ),
    "metric": (
        "log10 MAE in protein copies per cell over the held-out genes (gate (b)'s primary), for the"
        " saturable arm and the baseline, each with its training-fitted constant"
    ),
    "bar": (
        f"improvement = baseline MAE minus saturable MAE on held-out genes, at least"
        f" {MIN_IMPROVEMENT_LOG10} (gate (b)'s MIN_IMPROVEMENT_LOG10), with the 95% percentile interval"
        f" of {BOOTSTRAP} bootstrap resamples of the held-out genes (seed {SEED}, every fitted number"
        " held fixed) above 0, AND gate (b)'s secondary in band (total protein per cell from the"
        " declared ribosome capacity within one order of magnitude of Milo 2013; computed from TPM and"
        " residue counts, no allocation in it, so it is the same number under any policy)"
    ),
    "slope": (
        "reported, not a verdict: the OLS slope of the saturable arm's held-out predictions on held-out"
        " log10 TPM (no protein value in it), against the registered 0.515 [0.499, 0.531] of measured"
        " protein on transcript; inside or outside that interval is stated. The baseline's is 1 by"
        " construction"
    ),
    "comparator": (
        "reported, not a verdict: a two-parameter straight line, prediction = b log10 TPM + a, b on"
        " SLOPE_GRID (0 to 1.5, steps 0.001, a tie keeping the b nearest 1) and a the median training"
        " residual, fitted on training genes. It has the saturable arm's parameter count and no"
        " mechanism, so it says whether an improvement is specific to the saturable form. Its"
        " held-out MAE minus the saturable arm's gets the same bootstrap"
    ),
    "gene_blind_arms": (
        "`proportional` and `competitive` through the same runtime, pool driven short by multiplying"
        " every work rate by 1e4 as gate (b) did, constant fitted on training: their held-out"
        " improvement over the baseline must be 0 within 1e-9, as gate (b) found"
    ),
    "declared_capacity_arm": (
        "reported, not a verdict: K set to the capacity itself (5e6), the form §9.2 named,"
        " capacity / (capacity + demand_i), with only c fitted on training"
    ),
    "instrument_checks": (
        "INSTRUMENT_FAILED, and no reading made, if: the runtime's factors differ from the closed form"
        " above by more than 1e-9 in log10 on any gene; either gene-blind arm's held-out improvement"
        " differs from 0 by more than 1e-9; or fewer than 1,000 genes are held out"
    ),
    "verdict_rule": {
        "PASS": "improvement >= 0.05 AND its 95% interval's lower bound > 0 AND the secondary in band",
        "NO_IMPROVEMENT_DETECTED": "anything else that is not INSTRUMENT_FAILED",
    },
    "readings": {
        "PASS": (
            "On held-out genes, the saturable policy's log10 MAE is {sat} against the fitted"
            " baseline's {base}: an improvement of {imp} [{lo}, {hi}], which clears gate (b)'s"
            " registered bar of 0.05."
        ),
        "NO_IMPROVEMENT_DETECTED": (
            "On held-out genes, the saturable policy's log10 MAE is {sat} against the fitted"
            " baseline's {base}: an improvement of {imp} [{lo}, {hi}]; no improvement detected at"
            " gate (b)'s registered bar of 0.05."
        ),
        "INSTRUMENT_FAILED": "INSTRUMENT_FAILED: {why}. No reading is made.",
        "comparator_line_better": (
            "A two-parameter straight line does better than the saturable policy on the same held-out"
            " genes ({line} against {sat}, difference {d} [{lo}, {hi}]), so the improvement is not"
            " specific to the saturable form."
        ),
        "comparator_saturable_better": (
            "The saturable policy's held-out log10 MAE is lower than a two-parameter straight line's"
            " ({sat} against {line}, difference {d} [{lo}, {hi}])."
        ),
        "comparator_not_separated": (
            "The saturable policy and a two-parameter straight line are not separated on held-out"
            " genes ({sat} against {line}, difference {d} [{lo}, {hi}])."
        ),
        "slope": (
            "The saturable policy's held-out predictions rise {p} decades per decade of transcript,"
            " {where} the registered 0.515 [0.499, 0.531] of measured protein."
        ),
        "boundary_top": (
            "The fitted K is at the top of the grid: no saturation within it was preferred, and the"
            " policy reduces to the baseline."
        ),
        "boundary_bottom": (
            "The fitted K is at the bottom of the grid: every gene is saturated and the policy predicts"
            " nearly one constant."
        ),
        "binding_no": (
            "At real demand under the fitted K, the ribosome pool serves {served} against a capacity"
            " of 5e6, so the capacity does not bind; every factor below 1 comes from K, not from the"
            " pool."
        ),
        "binding_yes": (
            "At real demand under the fitted K, the ribosome pool serves {served} against a capacity"
            " of 5e6, so the capacity binds and every factor carries a pool cut of {cut} beside K."
        ),
    },
    "falsifier": (
        "the policy fails gate (b) on this join if its held-out improvement is below 0.05 or its"
        " interval reaches 0. A pass is falsified as a statement about the saturable FORM, though not"
        " as the registered reading, if the comparator line does better on the same genes"
    ),
    "planted_contrasts": {
        "what": (
            f"synthetic transcriptomes where the answer is known: {PLANTED_N} genes p0000..p3999,"
            f" log10 TPM uniform on [-1, 4] (numpy default_rng seed {SEED}), Gaussian noise of sd"
            f" {PLANTED_SD} on log10 copies, shift {PLANTED_SHIFT}; run through the same split, fits,"
            " runtime and verdict as the real data"
        ),
        "saturable_truth": (
            f"log10 copies = log10(TPM * {PLANTED_K} / ({PLANTED_K} + TPM)) + shift + noise. Expected:"
            " PASS; the fitted log10 K within 0.3 of log10 30; both gene-blind arms' improvement 0"
            " within 1e-9, i.e. they cannot pass where the answer is per-gene saturation"
        ),
        "null_truth": (
            "log10 copies = log10 TPM + shift + noise. Expected: NO_IMPROVEMENT_DETECTED (the fitted"
            " form must not manufacture a pass where nothing saturates)"
        ),
        "straight_truth": (
            "log10 copies = 0.5 log10 TPM + shift + noise. Expected: the comparator reading is"
            " comparator_line_better, which shows the comparator can say an improvement is not"
            " specific to the saturable form"
        ),
    },
    "expectation_stated_before_the_run": (
        "the bar is gate (b)'s, against a baseline whose slope is 1, and gate (b) measured 0.515. Any"
        " fitted compression is expected to beat a slope-1 baseline on such data: the planted straight"
        " truth (slope 0.5, no saturation in it) PASSES the bar with the line comparator better, and"
        " the dry run's synthetic protein arm, unrelated to TPM, also passes. So a PASS here is weak"
        " evidence for the saturable form; the comparator reading is the one that can tell the form"
        " from compression in general, and it is reported beside the verdict whatever it says"
    ),
    "cannot_do": [
        "a held-out improvement from a fitted per-gene compression is also what noise in the"
        " transcript arm would produce: GTEx fibroblast TPM against PaxDb's whole-organism average is"
        " two populations, and measurement noise in x attenuates a slope toward 0 with no biology in"
        " it. This gate cannot tell a saturable allocation from that attenuation; the comparator line"
        " only says whether the saturable shape adds anything beyond any compression",
        "a pass licenses the wording 'on held-out genes, the policy's log10 MAE is X against Y' and no"
        " statement about how a cell allocates ribosomes",
        "K is fitted, in TPM units; it is not a measured ribosome affinity and is not comparable to one",
        "the covered-set trap of gate (b) stands: PaxDb holds what could be quantified",
    ],
    "exploration_preceded_this": (
        "before writing this the lane had read gate (b)'s code, registration, amendment and tests, and"
        " the published aggregate numbers in BIOLANG-v0.4-ECONOMY.md §9.2 and ROADMAP-HISTORY.md (the"
        " two MAEs, the slope and its interval, the factor 1.754e-4, the coverage counts, 1.36e10). It"
        " had read the header and two rows of the distilled GTEx file (DDX11L1 and WASH7P TPM) and no"
        " row of the PaxDb file. No per-gene protein value, no split, no fit and no error of this gate"
        " had been computed. Before registering, `dry` ran the code on the real transcript arm with a"
        " synthetic protein arm (PaxDb not opened) to check the code runs; nothing from it is kept"
    ),
    "money": "none; no request is sent, no key is read; both arms are the cached files gate (b) fetched",
}


#: outside the binding words, which no amendment changes; each entry says what changed and why
AMENDMENTS: list[dict[str, str]] = [
    {
        "number": "1",
        "date": "2026-10-10",
        "what_happened": (
            "run 1, against the registration committed in 9d2ab1f, computed every reading and was then"
            " refused at the write by the result contract (genomeos/results.py): the manifest's"
            " `partitions` was a string, and the contract takes a dict or 'n/a: why'. save_result"
            " quarantined the payload at data/quarantine/results/saturable_gate_result.json. The run"
            " printed nothing before the write, and the quarantined file was not opened, so no number"
            " of run 1 was seen before this amendment"
        ),
        "change": (
            "the result manifest's `partitions` becomes a dict naming the split rule; `run` also"
            " compares its payload with run 1's quarantined payload key for key and records whether"
            " they are identical. Nothing in REGISTRATION, the data, the fits, the split, the bar or"
            " the readings changes; the pins are re-taken because the script's bytes changed"
        ),
        "rule": (
            "the reading is run 2's only if run 2's payload is identical to run 1's; if they differ,"
            " that is recorded and run 1's quarantined payload is the registered run"
        ),
    },
]
QUARANTINED_RUN_1 = ROOT / "data" / "quarantine" / "results" / f"{RESULT_NAME}.json"


# ---- the rules fixed above, as code -----------------------------------------------------------


def is_held_out(symbol: str) -> bool:
    """The split, fixed at registration: odd first 32 bits of sha256(symbol) are held out."""
    return int(hashlib.sha256(symbol.encode("utf-8")).hexdigest()[:8], 16) % 2 == 1


def closed_form_log_factor(tpm: np.ndarray, log10_k: float, capacity: float) -> np.ndarray:
    """log10 of `saturable`'s factor for every gene, as 91e2a85's runtime defines it."""
    k = 10.0**log10_k
    f = k / (k + tpm)
    served = float(np.sum(tpm * f))
    cut = capacity / served if served > capacity else 1.0
    return np.log10(f * cut)


def _mae(err: np.ndarray) -> float:
    return float(np.mean(np.abs(err)))


def fit_baseline(y: np.ndarray, x: np.ndarray) -> float:
    return float(np.median(y - x))


def fit_saturable(
    y_tr: np.ndarray, x_tr: np.ndarray, tpm_all: np.ndarray, train: np.ndarray, capacity: float
) -> dict[str, float]:
    """K on LOG10_K_GRID and c by median, minimising training MAE. Receives training y only."""
    best: dict[str, float] | None = None
    for lk in reversed(LOG10_K_GRID):  # from the top: a tie keeps the larger K
        p = x_tr + closed_form_log_factor(tpm_all, lk, capacity)[train]
        c = float(np.median(y_tr - p))
        m = _mae(y_tr - p - c)
        if best is None or m < best["train_mae"]:
            best = {"log10_k": lk, "c": c, "train_mae": m}
    assert best is not None
    return best


def fit_line(y_tr: np.ndarray, x_tr: np.ndarray) -> dict[str, float]:
    """The comparator: slope on SLOPE_GRID (ties keep the slope nearest 1), intercept by median."""
    best: dict[str, float] | None = None
    for b in sorted(SLOPE_GRID, key=lambda v: (abs(v - 1.0), v)):
        p = b * x_tr
        a = float(np.median(y_tr - p))
        m = _mae(y_tr - p - a)
        if best is None or m < best["train_mae"]:
            best = {"b": b, "a": a, "train_mae": m}
    assert best is not None
    return best


def declared_program(k: float) -> tuple[Pool, Allocation]:
    """The pool and allocation as a `.bio` program states them, parsed, so K enters by the language."""
    from genomeos.lang import parse

    text = (
        "module saturable_gate\n\n"
        "compartment Cytosol { volume: 1.0; translation: yes }\n"
        f"pool Ribosomes {{ location: Cytosol; size: {RIBOSOME_CAPACITY!r} }}\n"
        f"allocation ribosome_share {{ pool: Ribosomes; policy: saturable; half_saturation: {k!r} }}\n"
    )
    m = parse(text)
    return m.entities["Ribosomes"], m.entities["ribosome_share"]


def runtime_log_factors(
    symbols: list[str], tpm: np.ndarray, policy: str, k: float | None = None, multiplier: float = 1.0
) -> np.ndarray:
    """The REAL runtime's factor per gene, one pool, 1 per chain, work = TPM times `multiplier`."""
    if policy == "saturable":
        pool, al = declared_program(float(k))
    else:
        pool = Pool(id="Ribosomes", kind="pool", size=RIBOSOME_CAPACITY)
        al = Allocation(id="a", kind="allocation", pool="Ribosomes", policy=policy)
    econ = Economy(
        pools={"Ribosomes": pool},
        costs={s: [Cost(pool="Ribosomes", amount=1.0, per="chain")] for s in symbols},
        allocations={"Ribosomes": al},
    )
    econ.check()
    got = econ.scale({s: float(t) * multiplier for s, t in zip(symbols, tpm, strict=True)})
    return np.log10(np.array([got[s] for s in symbols]))


def _boot(err_a: np.ndarray, err_b: np.ndarray) -> tuple[float, float]:
    """95% percentile interval of mean|a| - mean|b| over bootstrap resamples of the held-out genes."""
    rng = np.random.default_rng(SEED)
    n = len(err_a)
    aa, bb = np.abs(err_a), np.abs(err_b)
    got = np.empty(BOOTSTRAP)
    for i in range(BOOTSTRAP):
        idx = rng.integers(0, n, size=n)
        got[i] = aa[idx].mean() - bb[idx].mean()
    got.sort()
    return float(got[int(0.025 * BOOTSTRAP)]), float(got[int(0.975 * BOOTSTRAP) - 1])


def _ols_slope(y: np.ndarray, x: np.ndarray) -> float:
    xm = x - x.mean()
    den = float(np.sum(xm * xm))
    return float(np.sum(xm * (y - y.mean())) / den) if den else float("nan")


def _r(v: float) -> float:
    return round(float(v), 6)


def run_gate(symbols: list[str], tpm: np.ndarray, y: np.ndarray) -> dict[str, Any]:
    """Every registered reading for one join. `y` is log10 copies; nothing about held-out y is fitted."""
    tpm = np.asarray(tpm, dtype=float)
    y = np.asarray(y, dtype=float)
    x = np.log10(tpm)
    held = np.array([is_held_out(s) for s in symbols])
    train = ~held
    n_tr, n_ho = int(train.sum()), int(held.sum())
    y_tr, x_tr = y[train], x[train]

    # ---- every fit, on training genes alone
    c_base = fit_baseline(y_tr, x_tr)
    sat = fit_saturable(y_tr, x_tr, tpm, train, RIBOSOME_CAPACITY)
    line = fit_line(y_tr, x_tr)
    k_hat = 10.0 ** sat["log10_k"]
    lf_decl = closed_form_log_factor(tpm, math.log10(RIBOSOME_CAPACITY), RIBOSOME_CAPACITY)
    c_decl = float(np.median(y_tr - x_tr - lf_decl[train]))
    blind_fit: dict[str, tuple[np.ndarray, float]] = {}
    for policy in ("proportional", "competitive"):
        lf = runtime_log_factors(symbols, tpm, policy, multiplier=GENE_BLIND_DEMAND_MULTIPLIER)
        blind_fit[policy] = (lf, float(np.median(y_tr - x_tr - lf[train])))

    # ---- the runtime, not the closed form, makes the saturable arm's predictions
    lf_rt = runtime_log_factors(symbols, tpm, "saturable", k=k_hat)
    lf_cf = closed_form_log_factor(tpm, sat["log10_k"], RIBOSOME_CAPACITY)
    closed_form_gap = float(np.max(np.abs(lf_rt - lf_cf)))

    # ---- held-out errors, every parameter frozen
    e_base = y[held] - (x[held] + c_base)
    e_sat = y[held] - (x[held] + lf_rt[held] + sat["c"])
    e_line = y[held] - (line["b"] * x[held] + line["a"])
    e_decl = y[held] - (x[held] + lf_decl[held] + c_decl)
    mae_base, mae_sat, mae_line = _mae(e_base), _mae(e_sat), _mae(e_line)
    improvement = mae_base - mae_sat
    lo, hi = _boot(e_base, e_sat)
    d_line = mae_line - mae_sat  # positive: the saturable arm is better
    l_lo, l_hi = _boot(e_line, e_sat)

    blind: dict[str, Any] = {}
    for policy, (lf, cb) in blind_fit.items():
        mae_b = _mae(y[held] - (x[held] + lf[held] + cb))
        blind[policy] = {
            "distinct_factors": len({round(float(v), 12) for v in lf}),
            "log10_factor": _r(lf[0]),
            "heldout_mae": _r(mae_b),
            "improvement_over_baseline": float(mae_base - mae_b),
        }

    pred_sat_ho = x[held] + lf_rt[held] + sat["c"]
    slope_sat = _ols_slope(pred_sat_ho, x[held])
    inside = REGISTERED_SLOPE_95[0] <= slope_sat <= REGISTERED_SLOPE_95[1]

    served = float(np.sum(tpm * (k_hat / (k_hat + tpm))))
    binds = served > RIBOSOME_CAPACITY
    readings = REGISTRATION["readings"]
    if binds:
        bind_reading = readings["binding_yes"].format(
            served=f"{served:.4g}", cut=f"{RIBOSOME_CAPACITY / served:.4g}"
        )
    else:
        bind_reading = readings["binding_no"].format(served=f"{served:.4g}")

    # ---- instrument checks, then the readings
    why = []
    if closed_form_gap > CLOSED_FORM_TOLERANCE_LOG10:
        why.append(f"runtime and closed form differ by {closed_form_gap:.3e} in log10")
    for policy, row in blind.items():
        if abs(row["improvement_over_baseline"]) > GENE_BLIND_TOLERANCE_LOG10:
            why.append(f"{policy} moved the held-out error by {row['improvement_over_baseline']:.3e}")
    if n_ho < MIN_HELD_OUT:
        why.append(f"only {n_ho} genes held out")

    passed_primary = improvement >= MIN_IMPROVEMENT_LOG10 and lo > 0.0
    fmt = {"sat": f"{mae_sat:.6f}", "base": f"{mae_base:.6f}", "imp": f"{improvement:.6f}"}
    fmt |= {"lo": f"{lo:.6f}", "hi": f"{hi:.6f}"}
    if l_hi < 0.0:
        comp_key = "comparator_line_better"
    elif l_lo > 0.0:
        comp_key = "comparator_saturable_better"
    else:
        comp_key = "comparator_not_separated"
    comp = readings[comp_key].format(
        sat=f"{mae_sat:.6f}", line=f"{mae_line:.6f}", d=f"{d_line:.6f}", lo=f"{l_lo:.6f}", hi=f"{l_hi:.6f}"
    )
    boundary = None
    if sat["log10_k"] == LOG10_K_GRID[-1]:
        boundary = readings["boundary_top"]
    elif sat["log10_k"] == LOG10_K_GRID[0]:
        boundary = readings["boundary_bottom"]

    return {
        "split": {"training": n_tr, "held_out": n_ho, "rule": REGISTRATION["split"]},
        "fitted_on_training": {
            "baseline_c": _r(c_base),
            "saturable_log10_k": sat["log10_k"],
            "saturable_k_tpm": _r(k_hat),
            "saturable_c": _r(sat["c"]),
            "saturable_training_mae": _r(sat["train_mae"]),
            "line_b": line["b"],
            "line_a": _r(line["a"]),
            "line_training_mae": _r(line["train_mae"]),
            "declared_capacity_c": _r(c_decl),
            "gene_blind_c": {p: _r(cb) for p, (_, cb) in blind_fit.items()},
        },
        "held_out": {
            "baseline_mae": _r(mae_base),
            "saturable_mae": _r(mae_sat),
            "improvement": _r(improvement),
            "improvement_95": [_r(lo), _r(hi)],
            "bar": MIN_IMPROVEMENT_LOG10,
            "primary_clears_bar": bool(passed_primary),
            "line_mae": _r(mae_line),
            "line_minus_saturable": _r(d_line),
            "line_minus_saturable_95": [_r(l_lo), _r(l_hi)],
            "declared_capacity_arm_mae": _r(_mae(e_decl)),
        },
        "slope": {
            "saturable_predictions_on_log10_tpm_heldout": _r(slope_sat),
            "registered_measured": REGISTERED_SLOPE,
            "registered_measured_95": list(REGISTERED_SLOPE_95),
            "inside_registered_interval": bool(inside),
            "line_b": line["b"],
            "baseline": 1.0,
            "reading": readings["slope"].format(
                p=f"{slope_sat:.3f}", where="inside" if inside else "outside"
            ),
        },
        "binding": {
            "capacity": RIBOSOME_CAPACITY,
            "served_total_at_fitted_k": served,
            "demand_total": float(np.sum(tpm)),
            "capacity_binds": bool(binds),
            "genes_above_k": int(np.sum(tpm > k_hat)),
            "reading": bind_reading,
        },
        "gene_blind_arms": blind,
        "instrument": {
            "runtime_vs_closed_form_max_log10_gap": closed_form_gap,
            "failed": bool(why),
            "why": why,
        },
        "comparator_reading": comp,
        "comparator_key": comp_key,
        "boundary": boundary,
        "_fmt": fmt,
        "_primary": bool(passed_primary),
    }


def verdict(gate: dict[str, Any], secondary_in_band: bool) -> tuple[str, str]:
    """The registered verdict and its reading, word for word."""
    if gate["instrument"]["failed"]:
        key = "INSTRUMENT_FAILED"
        return key, REGISTRATION["readings"][key].format(why="; ".join(gate["instrument"]["why"]))
    key = "PASS" if gate["_primary"] and secondary_in_band else "NO_IMPROVEMENT_DETECTED"
    return key, REGISTRATION["readings"][key].format(**gate["_fmt"])


def planted(case: str) -> tuple[list[str], np.ndarray, np.ndarray]:
    """A synthetic transcriptome whose answer is known (REGISTRATION['planted_contrasts'])."""
    rng = np.random.default_rng(SEED)
    symbols = [f"p{i:04d}" for i in range(PLANTED_N)]
    x = rng.uniform(-1.0, 4.0, size=PLANTED_N)
    tpm = 10.0**x
    noise = rng.normal(0.0, PLANTED_SD, size=PLANTED_N)
    if case == "saturable_truth":
        y = np.log10(tpm * PLANTED_K / (PLANTED_K + tpm)) + PLANTED_SHIFT + noise
    elif case == "null_truth":
        y = x + PLANTED_SHIFT + noise
    elif case == "straight_truth":
        y = 0.5 * x + PLANTED_SHIFT + noise
    else:
        raise ValueError(case)
    return symbols, tpm, y


def planted_outcomes() -> dict[str, Any]:
    out: dict[str, Any] = {}
    for case in ("saturable_truth", "null_truth", "straight_truth"):
        g = run_gate(*planted(case))
        key, _ = verdict(g, secondary_in_band=True)
        out[case] = {
            "verdict": key,
            "fitted_log10_k": g["fitted_on_training"]["saturable_log10_k"],
            "improvement": g["held_out"]["improvement"],
            "improvement_95": g["held_out"]["improvement_95"],
            "comparator_key": g["comparator_key"],
            "gene_blind_improvement": {
                p: row["improvement_over_baseline"] for p, row in g["gene_blind_arms"].items()
            },
        }
    return out


# ---- the real join ---------------------------------------------------------------------------


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def real_join(synthetic_protein: bool = False) -> dict[str, Any]:
    """Gate (b)'s join by gate (b)'s readers. `synthetic_protein` (the dry run) never opens PaxDb."""
    from abundance_gate import gtex, lengths, paxdb  # noqa: PLC0415 - gate (b)'s script

    rna = gtex()
    length = lengths(set(rna))
    if synthetic_protein:
        rng = np.random.default_rng(SEED)
        protein_ppm = {s: float(10.0 ** rng.uniform(-2, 3)) for s in sorted(rna)}
    else:
        protein_ppm = paxdb()
    kept: list[str] = []
    for sym in sorted(set(protein_ppm) & set(rna)):
        copies = protein_ppm[sym] * 1e-6 * MEASURED_TOTAL_PROTEINS
        if copies < MIN_COPIES or rna[sym]["declared_tissue"] <= 0.0 or sym not in length:
            continue
        kept.append(sym)
    tpm = np.array([rna[s]["declared_tissue"] for s in kept])
    y = np.log10(np.array([protein_ppm[s] * 1e-6 * MEASURED_TOTAL_PROTEINS for s in kept]))
    # gate (b)'s secondary, from TPM and residue counts only
    weighted_len = float(sum(length[s] * rna[s]["declared_tissue"] for s in kept) / tpm.sum())
    chains_per_s = RIBOSOME_CAPACITY * ELONGATION_AA_PER_S / weighted_len
    total_model = chains_per_s * PROTEIN_HALF_LIFE_H * 3600.0 / math.log(2.0)
    in_band = MILO_BAND[0] / 10.0 <= total_model <= MILO_BAND[1] * 10.0
    return {
        "symbols": kept,
        "tpm": tpm,
        "y": y,
        "rna": rna,
        "secondary": {
            "expression_weighted_mean_residues": round(weighted_len, 1),
            "model_total_proteins_per_cell": total_model,
            "milo_2013_band": list(MILO_BAND),
            "within_one_order_of_magnitude": bool(in_band),
            "policy_independent": True,
        },
    }


def _binding_words() -> dict[str, Any]:
    return json.loads(json.dumps(REGISTRATION))


def _manifest(inputs: list[dict[str, Any]], partitions: str) -> dict[str, Any]:
    from genomeos import manifest as mf

    # amendment 1: the contract takes a dict or 'n/a: why', and run 1 passed the split rule's text
    if not partitions.startswith("n/a:"):
        partitions = {"rule": partitions, "training": "even", "held_out": "odd"}  # type: ignore[assignment]

    return {
        "sources": [
            {"accession": "PaxDb 9606-WHOLE_ORGANISM-integrated", "version": "as fetched by gate (b)"},
            {"accession": "GTEx v8 gene median TPM, distilled by gate (b)", "version": "v8"},
            {"accession": "the packaged proteome, data/knowledge/proteins", "version": "on disk"},
        ],
        "inputs": inputs,
        "assembly": "n/a: gene-level abundances joined by symbol",
        "coordinates": "n/a: no genomic intervals",
        "parameters": {
            "seed": SEED,
            "bootstrap": BOOTSTRAP,
            "log10_k_grid": [LOG10_K_GRID[0], LOG10_K_GRID[-1], 0.01],
            "slope_grid": [SLOPE_GRID[0], SLOPE_GRID[-1], 0.001],
            "ribosome_capacity": RIBOSOME_CAPACITY,
            "min_improvement_log10": MIN_IMPROVEMENT_LOG10,
            "declared_tissue": TRANSCRIPT_SOURCE["tissue"],
            "min_copies": MIN_COPIES,
        },
        "exclusions": [
            "gate (b)'s: below 10 copies per cell, zero TPM in the declared tissue, no residue count",
        ],
        "partitions": partitions,
        "code_cleanliness": mf.code_cleanliness(__file__, OWN_CODE),
    }


def register() -> int:
    from genomeos import manifest as mf
    from genomeos.results import save_result

    pins = {p: _sha(ROOT / p) for p in PINNED_CODE}
    for f in (PAXDB_FILE, GTEX_FILE):
        pins[str(f.relative_to(ROOT))] = _sha(f)
    out = {
        "result": REGISTRATION_NAME,
        "date": "2026-10-10",
        "lane": "lane-saturable",
        "roadmap_item": (
            "area A, Next 1: a per-demander saturable allocation policy, registered against gate (b)"
        ),
        "policy_commit": "91e2a85",
        "registration": REGISTRATION,
        "planted_outcomes_at_registration": planted_outcomes(),
        "pinned_sha256": pins,
        "what_the_run_will_write": f"data/results/{RESULT_NAME}.json",
        "protein_values_read": 0,
        "amendments": AMENDMENTS,
    }
    inputs = [
        mf.input_entry(PAXDB_FILE.relative_to(ROOT), partition="pinned by sha256; no value read"),
        mf.input_entry(GTEX_FILE.relative_to(ROOT), partition="pinned by sha256; no value read"),
    ]
    print("saved", save_result(REGISTRATION_NAME, out, manifest=_manifest(inputs, "n/a: a registration")))
    for case, row in out["planted_outcomes_at_registration"].items():
        print(case, row["verdict"], row["fitted_log10_k"], row["improvement"], row["comparator_key"])
    return 0


def _check_registration() -> dict[str, Any]:
    reg = json.loads(REGISTRATION_FILE.read_text())
    if reg["registration"] != _binding_words():
        raise SystemExit("REFUSED: the registration's binding words differ from REGISTRATION in this script")
    for path, sha in reg["pinned_sha256"].items():
        if _sha(ROOT / path) != sha:
            raise SystemExit(f"REFUSED: {path} has changed since the registration pinned it")
    return reg


def dry() -> int:
    """The real transcript arm, a synthetic protein arm: the run's code end to end, nothing written."""
    j = real_join(synthetic_protein=True)
    g = run_gate(j["symbols"], j["tpm"], j["y"])
    key, _ = verdict(g, j["secondary"]["within_one_order_of_magnitude"])
    print("dry run on synthetic protein values (PaxDb not opened):", len(j["symbols"]), "genes")
    print(
        " split",
        g["split"]["training"],
        g["split"]["held_out"],
        "instrument failed:",
        g["instrument"]["failed"],
    )
    print(" verdict key on synthetic values, meaningless:", key)
    return 0


def run() -> int:
    from genomeos import manifest as mf
    from genomeos.results import save_result

    _check_registration()
    j = real_join()
    g = run_gate(j["symbols"], j["tpm"], j["y"])
    key, reading = verdict(g, j["secondary"]["within_one_order_of_magnitude"])
    fmt = g.pop("_fmt")
    g.pop("_primary")
    out = {
        "result": RESULT_NAME,
        "date": "2026-10-10",
        "lane": "lane-saturable",
        "registration": f"data/results/{REGISTRATION_NAME}.json",
        "analysed_genes": len(j["symbols"]),
        **g,
        "secondary": j["secondary"],
        "verdict": key,
        "reading": reading,
        "reading_numbers": fmt,
    }
    residues = [PROTEINS / f"{s}.json" for s in sorted(j["rna"]) if (PROTEINS / f"{s}.json").exists()]
    inputs = [
        mf.input_entry(REGISTRATION_FILE.relative_to(ROOT), partition=None),
        mf.input_entry(PAXDB_FILE.relative_to(ROOT), partition="training and held-out, split by symbol hash"),
        mf.input_entry(GTEX_FILE.relative_to(ROOT), partition="training and held-out, split by symbol hash"),
        mf.files_entry("data/knowledge/proteins/<symbol>.json for every GTEx symbol that has one", residues),
    ]
    # amendment 1: run 1 was quarantined unread; its payload must equal this run's, key for key
    if QUARANTINED_RUN_1.exists():
        q = json.loads(QUARANTINED_RUN_1.read_text())
        differ = sorted(k for k in out if q.get(k) != json.loads(json.dumps(out[k])))
        out["run_1_quarantined_payload"] = {
            "path": str(QUARANTINED_RUN_1.relative_to(ROOT)),
            "identical_to_this_run": not differ,
            "keys_that_differ": differ,
        }
        inputs.append(mf.input_entry(QUARANTINED_RUN_1.relative_to(ROOT), partition=None))
    print("saved", save_result(RESULT_NAME, out, manifest=_manifest(inputs, REGISTRATION["split"])))
    if "run_1_quarantined_payload" in out:
        print("run 1 payload identical:", out["run_1_quarantined_payload"]["identical_to_this_run"])
    for line in (reading, g["comparator_reading"], g["slope"]["reading"], g["binding"]["reading"]):
        print(line)
    if g["boundary"]:
        print(g["boundary"])
    return 0


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else ""
    raise SystemExit({"register": register, "dry": dry, "run": run}.get(cmd, lambda: 2)())
