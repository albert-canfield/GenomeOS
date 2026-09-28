# SPDX-License-Identifier: Apache-2.0
# Part of the BioLang engine (language, IR, VM, standard library); see LICENSING.md.
"""Pre-registration for stage 4's gate, written before the instrument and before any run.

§8 states the stage and the gate: *"A `divide` event partitions contents (binomially for low copy
numbers) instead of duplicating them, guarded by a size or resource checkpoint, so growth costs
something. Gate: over several divisions dilution must separate a stable protein from a short-lived
one."*

**The gate as worded cannot be met, and the reason is arithmetic that can be written down before the
instrument exists.** Partitioning multiplies every species in the cell by the same one-half. It
carries no index over proteins, exactly as stage 2's `Economy._share` carried no index over genes
(§9.2). In a chase, after `n` divisions of interval `T`, a protein of half-life `t` stands at

    A_n = A_0 · 2^(-n) · 2^(-nT/t)

and the ratio between a stable protein and a short-lived one is

    log10(A_stable / A_short) = n · T · log10(2) · (1/t_short − 1/t_stable)

in which **the division term has cancelled**. Dilution therefore separates nothing: it removes both
proteins by the same factor, and what separates them is turnover, which separates them just as well
in a cell that never divides. Worse for the gate's wording: in a steady state with synthesis on, a
level is `s / (k_degradation + k_division)`, so dividing *adds a constant to both denominators* and
**compresses** the separation. With the constants below the stable:short ratio falls from 40.0 in a
cell that does not divide to 4.545… in one that divides daily — a compression of 8.80×.

This is registered here, before the run, as a prediction the instrument must reproduce, and as the
reason the verdict will not be able to say "stage 4's gate passed" whatever the numbers are. A run
that *did* show dilution separating the two would mean the partitioning is indexed by protein, which
would be a bug of exactly the shape that looks like the finding.

**What dilution does do, and what this gate therefore asks instead.** A protein the program gives no
half-life is never degraded, so in a non-dividing cell it never leaves. Division is the only removal
it has. That is the one separation dilution makes on its own — between what turnover can remove and
what only division can — and it is askable here, so it is the gate's positive clause. The clauses
below are the ones this project can ask with what it holds; each can fail, and two are falsifier arms
whose whole purpose is to catch an implementation that passes the first clause for the wrong reason.

**What is NOT askable, registered as unaskable rather than attempted.** Whether the half-lives used
here are the half-lives of any real protein. `bio.std.human_turnover` holds **cell** lifespans
(Sender & Milo 2021), not protein half-lives, and nothing under `data/knowledge` carries a protein
turnover table; §10 decision 6 (which turnover set, and whether a mouse one is acceptable) is open
and this registration does not close it. The constants below are therefore **declared, not
measured**, and no clause compares them with a measurement. The same applies to the partitioning
model itself: `binomial` is the null for independent molecules, and Huh & Paulsson 2011 show that
clustering makes real partitioning noisier than binomial. This project holds no partitioning-error
dataset, so the binomial arm is checked against **its own analytic mean and variance** and against
conservation, never against a measured spread.

**And what this stage must not be said to have done for stage 2.** §9.2 named what would let gate (b)
move: a term that differs per gene, worth 0.485 decades per decade of compression. A per-gene
half-life is such a term, and stage 4 makes one expressible. It does not follow that gate (b) would
now pass: that is a different run, on a proteome, and it needs the measured turnover set this
registration has just said the project does not hold. Nothing here is evidence about gate (b).
"""

from __future__ import annotations

import math
from typing import Any

#: how many divisions "over several divisions" means. Fixed here so the horizon is not chosen after
#: seeing which horizon makes the separation look best.
DIVISIONS = 8

#: hours between divisions. A day, which is the order of a cultured mammalian cell cycle.
INTERVAL_H = 24.0

#: the stable protein's half-life, in hours. Declared, not measured (see the docstring).
STABLE_HALF_LIFE_H = 240.0

#: the short-lived protein's half-life, in hours. Declared, not measured.
SHORT_HALF_LIFE_H = 6.0

#: the amount each protein starts the chase with. Far above the regime threshold, so the
#: deterministic arm is deterministic.
START_AMOUNT = 1e6

#: copies for the stochastic arm, below the regime's default `threshold: 50`, so partitioning is
#: binomial there and deterministic above.
LOW_COPIES = 20

#: seeds the binomial arm runs, fixed before the run.
SEEDS = 2000

#: the separation the closed form predicts after DIVISIONS divisions, in decades of log10. It is the
#: same number in the partitioned arm and the duplicating arm, which is the whole point.
PREDICTED_SEPARATION_DECADES = (
    DIVISIONS * INTERVAL_H * math.log10(2.0) * (1.0 / SHORT_HALF_LIFE_H - 1.0 / STABLE_HALF_LIFE_H)
)

#: steady-state stable:short ratio in a cell that never divides, and in one that divides every
#: INTERVAL_H hours. The second is smaller: division compresses the separation.
PREDICTED_STEADY_RATIO_STATIC = STABLE_HALF_LIFE_H / SHORT_HALF_LIFE_H
PREDICTED_STEADY_RATIO_DIVIDING = (1.0 / SHORT_HALF_LIFE_H + 1.0 / INTERVAL_H) / (
    1.0 / STABLE_HALF_LIFE_H + 1.0 / INTERVAL_H
)
PREDICTED_COMPRESSION = PREDICTED_STEADY_RATIO_STATIC / PREDICTED_STEADY_RATIO_DIVIDING

#: the factor a protein with no declared half-life falls by over DIVISIONS divisions: dilution alone.
PREDICTED_DILUTION_ONLY = 2.0**-DIVISIONS

#: relative tolerance for a clause that compares a run against a closed form. Tight, because both
#: sides are the same arithmetic in floating point and anything looser would hide a real difference.
CLOSED_FORM_RTOL = 1e-9

#: bars for the binomial arm, fixed before the run. The mean bar is five standard errors of the mean
#: at SEEDS seeds (sqrt((N/4)/SEEDS) = 0.05 for N = LOW_COPIES); the variance bar is 15% of N/4.
BINOMIAL_MEAN_TOL = 0.25
BINOMIAL_VAR_RTOL = 0.15

PRE_REGISTRATION: dict[str, Any] = {
    "written": (
        "2026-09-28, before genomeos/runtime/division.py and scripts/division_gate.py existed and "
        "before any number was produced; committed in this file before either of them"
    ),
    "claim": (
        "a divide event that partitions contents removes a protein the program never gives a "
        "half-life to, which no amount of turnover can remove, and removes it by exactly 2^-n over "
        "n divisions; and partitioning is conservative, so no molecule is created by dividing"
    ),
    "counter_claim_registered_in_advance": (
        "dilution does NOT separate a stable protein from a short-lived one. Partitioning "
        "multiplies every species by the same one-half, so the division term cancels out of the "
        f"stable:short ratio, and the instrument must return the same {PREDICTED_SEPARATION_DECADES:.6f} "
        "decades whether contents are partitioned or duplicated. In a steady state division "
        f"compresses the separation from {PREDICTED_STEADY_RATIO_STATIC:.4f} to "
        f"{PREDICTED_STEADY_RATIO_DIVIDING:.4f}, a factor of {PREDICTED_COMPRESSION:.4f}. §8's gate "
        "as worded is therefore refuted rather than passed, and the verdict says so"
    ),
    "design": (
        "one program (data/demo/stage4_division.bio) with three proteins that differ only in the "
        f"half-life they declare: a stable one at {STABLE_HALF_LIFE_H} h, a short-lived one at "
        f"{SHORT_HALF_LIFE_H} h, and one that declares none at all. A chase (no synthesis) runs "
        f"{DIVISIONS} divisions {INTERVAL_H} h apart from {START_AMOUNT:.0e} of each, in two arms: "
        "contents partitioned, and contents duplicated (the pre-stage-4 behaviour, which is the "
        "falsifier arm for every dilution claim). A separate stochastic arm partitions "
        f"{LOW_COPIES} copies over {SEEDS} seeds"
    ),
    "clauses": {
        "1_closed_form": (
            "every amount in the partitioned arm equals A0 * 2^-n * 2^(-nT/t) and every amount in "
            f"the duplicating arm equals A0 * 2^(-nT/t), to a relative {CLOSED_FORM_RTOL:g}"
        ),
        "2_dilution_only": (
            "the protein with no declared half-life stands at exactly "
            f"{PREDICTED_DILUTION_ONLY} of its start in the partitioned arm and at exactly 1.0 of "
            "it in the duplicating arm. This is the positive clause: the one removal that is "
            "dilution's and turnover's alone cannot be"
        ),
        "3_dilution_separates_nothing": (
            "the stable:short separation is identical in the two arms to a relative "
            f"{CLOSED_FORM_RTOL:g}. If it is not, partitioning is indexed by protein, which is a bug"
        ),
        "4_compression": (
            "the steady-state stable:short ratio with division equals "
            f"{PREDICTED_STEADY_RATIO_DIVIDING:.6f} and without it {PREDICTED_STEADY_RATIO_STATIC:.6f}"
        ),
        "5_conservation": (
            "at every division, in both the deterministic and the binomial arm, the two daughters' "
            "amounts sum to the parent's exactly; no species and no seed is exempt"
        ),
        "6_binomial": (
            f"partitioning {LOW_COPIES} copies over {SEEDS} seeds has mean within "
            f"{BINOMIAL_MEAN_TOL} of {LOW_COPIES / 2} and variance within {BINOMIAL_VAR_RTOL:.0%} of "
            f"{LOW_COPIES / 4}, every daughter is a whole number, and the same species at "
            f"{START_AMOUNT:.0e} copies (above the regime threshold) has variance exactly 0"
        ),
        "7_equal_half_lives": (
            "falsifier: with both half-lives set equal, the two proteins are identical at every "
            "division in both arms, difference exactly 0. A difference would mean the separation "
            "came from bookkeeping rather than from turnover"
        ),
        "8_checkpoint": (
            "falsifier: with the checkpoint's condition unmet, 0 divisions fire, and every amount "
            "equals the decay-only value. A model that divides anyway is not checkpointed, and "
            "growth would then cost nothing, which is the thing §8 asks the checkpoint to prevent"
        ),
    },
    "not_askable": {
        "measured_turnover": (
            "whether these half-lives are any protein's. bio.std.human_turnover holds cell "
            "lifespans, not protein half-lives, and no protein turnover table is under "
            "data/knowledge. §10 decision 6 stays open; what would close it is a digitised set "
            "(Mathieson et al. 2018, Nat Commun 9:689, or Schwanhäusser et al. 2011, Nature "
            "473:337) with its cell type and method stated"
        ),
        "partitioning_noise": (
            "whether real partitioning is binomial. Huh & Paulsson 2011 show clustering makes it "
            "noisier; this project holds no partitioning-error dataset, so the binomial arm is "
            "checked against its own analytic moments and never against a measurement"
        ),
        "gate_b": (
            "whether a per-gene half-life would move stage 2's gate (b). That is a proteome-scale "
            "run needing the turnover set above; nothing here is evidence about it"
        ),
    },
    "must_not_move": {
        "body_runtime": (
            "genomeos/runtime/body.py is not touched and `partition` defaults to `duplicate`, so "
            "area E's numbers stand: 1,439 of 1,439 fates, terminal 496 of 555, 45 ambiguous "
            "decision points to 0, 522 of 555 under §7.5"
        ),
        "stage_2": "burden factors 1.0,1.0,1.0,1.0,0.50,0.25,0.05; abundance MAE 0.952005, improvement 0.0",
        "stage_1": "13 of 13 chrM, 1,123 of 1,123 imported, 0 with the routes closed",
        "engine_package": "12 of 12 checks",
    },
}
