"""AstroREG-2: a NEW registration whose gate is ENDPOINT VALIDITY, not scale, written before anything
is computed.

WHAT THIS IS NOT. It is not the old gate loosened, and it is not the old test passing. The original
registration (data/results/astroreg_registration.json, sha256
43f2edfab0636e0f07a035aedf5ddf7c57aa1dfcee8d39ea32e9338edd522c7a) is FINAL and its calibration gate
FAILED. A gate fixed in advance is not reopened after it fails, not even by the person who now believes
it asked the wrong question -- that is the entire point of fixing it in advance. Nothing here rescales
anything to satisfy it, and nothing here re-runs it.

WHY A DIFFERENT GATE IS A DIFFERENT QUESTION. The registered endpoint is an AUPRC DIFFERENCE between two
arms. An AUPRC depends only on how pairs are RANKED. Activity enters the frozen model only through
`log_activity = log1p(sqrt(dhs * h3k27ac))` and `activity_over_distance = log_activity - log(distance)`,
both in log space, so a per-column factor on dhs and h3k27ac becomes an additive shift on both features
-- and an additive shift common to every pair cannot move a ranking at all. The old gate measured that
factor and refused it; this one measures whether the endpoint moves. Those are different questions and
the second is the one the claim asks.

AND THE ARGUMENT IS WEAKER THAN IT LOOKS, which this lane checked rather than took, and which is recorded
here BEFORE the gate runs because it is the reason the gate can genuinely fail. `log1p` is not `log`. The
shift `log1p(cA) - log1p(A)` equals `log c` only where the activity A is large; it goes to ZERO as A goes
to zero. On the 1,918 K562 held-out pairs, with the implied constant factor c = sqrt(0.7603 * 0.8080) =
0.78379 and log c = -0.24362, that shift ranges from 0.0000 to -0.2402 across the pairs: its SPREAD,
0.2402, is as large as the asymptote itself. 37% of the pairs have activity below 1.0, 16% below 0.5, and
that is precisely where `log1p` is non-linear. So the a priori argument does NOT establish that the
ranking is preserved on this data; it establishes only that it would be if every pair sat in the large-A
regime, and a third of them do not. On top of that the real reconstruction is not an exact constant
factor at all -- the per-element ratio ran p05 0.643 to p95 0.879 on DNase and 0.756 to 1.156 on H3K27ac
-- so there is per-pair scatter as well as a shift. Whether the endpoint survives is therefore an
empirical question with a real chance of answering no, which is what makes this a gate.
"""

from __future__ import annotations

import gzip
import json
import math
import sys
import time
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import crispri  # noqa: E402
from genomeos.attribution.measured import CRISPRI_SPLIT_OF  # noqa: E402
from genomeos.results import save_result  # noqa: E402

RESULT = "astroreg2_registration"
ENTRY = "scripts/astroreg2_register.py"
OWN_CODE = ("genomeos/attribution/rpm.py", "scripts/astroreg2_register.py", "tests/test_rpm.py")

ORIGINAL = Path("data/results/astroreg_registration.json")
ORIGINAL_SHA256 = "43f2edfab0636e0f07a035aedf5ddf7c57aa1dfcee8d39ea32e9338edd522c7a"

# ------------------------------------------------------------------- the gate, fixed before any run

#: The endpoint tolerance. The gain is an AUPRC difference, so this is in AUPRC units.
MAX_ABS_GAIN_SHIFT = 0.01

#: The rank agreement required between the two runs' pair scores.
MIN_SCORE_SPEARMAN = 0.99

GATE = {
    "what_is_computed": "the deletion gain on the K562 HELD-OUT pairs TWICE, with the frozen weights "
    "unchanged: once with the PUBLISHED activity columns, once with the activity RECONSTRUCTED under "
    "amendment 2's rule and NO RESCALING of any kind",
    "weights": "fitted once on the K562 TRAINING pairs with the PUBLISHED columns, and reused "
    "unchanged for both runs. Only the held-out pairs' activity differs between the two, which is what "
    "makes the comparison about the input and not about the fit",
    "pass_requires_both": True,
    "max_abs_gain_shift": MAX_ABS_GAIN_SHIFT,
    "min_score_spearman": MIN_SCORE_SPEARMAN,
    "rule": f"PASS only if BOTH |gain_reconstructed - gain_published| <= {MAX_ABS_GAIN_SHIFT} AND the "
    f"Spearman of the per-pair scores between the two runs >= {MIN_SCORE_SPEARMAN}. Either one failing "
    "fails the gate. The endpoint tolerance alone could be met by two rankings that differ and happen "
    "to integrate to the same area; the rank agreement alone could be met while the gain still moved "
    "more than the tolerance allows",
    "both_numbers_written_before_either_was_computed": True,
    "no_rescaling": "the reconstructed columns enter exactly as computed. Fitting any factor, even one "
    "derived from K562, would make the astrocyte values depend on a K562 quantity and would be a "
    "refitting of the feature definition, which the claim forbids",
    "on_failure": "the route CLOSES. If the reconstruction moves the gain by more than 0.01 or "
    "disagrees on ranking, then the reconstructed activity is not interchangeable with the published "
    "activity for this endpoint, no astrocyte value may be computed from it, and that is the result",
}

# ------------------------------------------------------------- the history, disclosed in full here

HISTORY = {
    "this_is_not_a_first_attempt": "AstroREG-2 was written AFTER two scale gates had failed and after "
    "both of their numbers had been seen. A reader must not be able to mistake it for a first attempt, "
    "so the figures are here rather than in a covering note",
    "attempt_1": {
        "rule": "this lane's own read rule, registered before it ran: ENCODE's FILTERED GRCh38 "
        "alignments, duplicates kept, supplementary dropped, one filter for numerator and denominator, "
        "counts pooled across replicates",
        "result": "data/results/astroreg_calibration_h3k27ac_k562.json",
        "h3k27ac_spearman": 0.9733140501519615,
        "h3k27ac_median_ratio": 0.6070182726597586,
        "verdict": "FAILED both conditions of the scale gate",
        "what_it_established": "the denominator was right -- the computed total matched ENCODE's own "
        "published mapped-read figure exactly, 86,297,833 -- so the gap was located in the numerator",
    },
    "attempt_2": {
        "rule": "the PRODUCER's own documented rule, read from mayasheth/chrom-annotate at commit "
        "91cda73ebe3a19153a582cab18cbf7ff70d85cfc with a line cited for every term, over the BAMs its "
        "own metadata names. No candidate rule was tried against the data",
        "result": "data/results/astroreg_calibration_producer_rule_k562.json",
        "dhs_spearman": 0.99747,
        "dhs_median_ratio": 0.7603,
        "h3k27ac_spearman": 0.99385,
        "h3k27ac_median_ratio": 0.8080,
        "verdict": "FAILED the scale condition on both columns while PASSING the rank condition on both",
        "what_it_established": "the producer's documented rule reproduces the ORDERING of both columns "
        "almost exactly and leaves a residual multiplicative gap, 1.315x on DNase and 1.238x on "
        "H3K27ac. The two factors differ, so no single shared constant explains both",
    },
    "the_old_registration_is_final": "its gate failed and stays failed. AstroREG-2 does not reopen it, "
    "does not rescale anything to satisfy it, and must never be described as the old test passing",
    "why_a_new_registration_rather_than_a_third_amendment": "the old gate asked whether a QUANTITY is "
    "reproduced. This asks whether the ENDPOINT is unchanged. A third amendment to a failed gate would "
    "be the thing a registration exists to prevent; a new registration with a different, "
    "pre-committed, falsifiable gate is not",
}

# ----------------------------------------------------------------------------- what is still blind

BLIND = {
    "not_computed_when_this_was_written": [
        "the K562 ENDPOINT comparison: neither gain has been computed, with published or reconstructed "
        "activity, and the per-pair score correlation between the two runs is unknown",
        "no astrocyte outcome has been read, no astrocyte activity has been reconstructed, and no "
        "astrocyte deletion value has been requested or computed",
    ],
    "what_was_already_seen": "the two scale gates' correlations and ratios, listed under `history`, "
    "and the published activity distribution quoted in the caveat below. Those are what make this a "
    "reconstruction check rather than a blind one, and they are disclosed rather than minimised",
    "the_gate_can_genuinely_fail": "a third of the held-out pairs sit where log1p is non-linear and "
    "the per-element reconstruction carries real scatter, so the endpoint is not guaranteed to survive. "
    "If it moves by more than 0.01 the route closes. That is what distinguishes a gate from a "
    "formality",
}

# ------------------------------------------------------------------------------------- the claim

CLAIM_UNCHANGED = {
    "claim": "the deletion gain in cultured primary human astrocytes: AUPRC of 'activity + distance + "
    "deletion' minus AUPRC of 'activity + distance', both scored with the weights frozen on the K562 "
    "training pairs. Nothing is refitted on astrocyte data -- not the weights, not a threshold, not a "
    "feature definition",
    "what_changed": "only the ROUTE to one of the claim's inputs: how the astrocyte activity term "
    "would be obtained. The endpoint, the arms, the labels, the interval rule, the three readings and "
    "the not_claimed list are the original's and are untouched",
    "why_the_power_figure_carries_over": "a power figure describes a CLAIM -- its effect size, its "
    "prevalence, its number of positives and independent loci -- and none of those changed. 133 "
    "positives over 74 independent loci at 2.9% prevalence is still what is being tested, so the "
    "registered figure still describes it. Had the claim changed, the figure would have had to be "
    "recomputed; it did not, so it does not",
    "required_quote": "exploratory (prevalence-matched power 0.79, conservative: simulated with ~27 "
    "positives against the test's 133)",
    "required_quote_rule": "this exact form of words is what every quote of the class must use; the "
    "class may not travel without the figure and both limitations",
    "travels_beside_it": {
        "attenuation_0_5x_share": 0.81,
        "rule": "the 0.5x figure travels beside the class in every quote of it",
        "bar": "the 0.5x share may NOT be quoted as power at half the effect: shrinking the deletion "
        "arm scales the separation and the bootstrap spread together, so the shares stay nearly flat "
        "across 1x, 0.5x and 0.25x while the median achieved gain falls. It shows that the "
        "detectability of a proportionally smaller separation is similar, not that power survives a "
        "weaker true effect against unchanged noise",
    },
    "limitation_prevalence": "the first figure of 0.96 was computed at K562's own 6.5% prevalence and "
    "does not describe a test run at 2.9%; each resample is matched to the astrocyte prevalence by "
    "downsampling its POSITIVES, which is the direction that cannot invent data",
    "limitation_sample_size": "the simulation keeps a mean of about 27 positives, because that is what "
    "74 of K562's loci can supply at the matched prevalence, while the real test has 133. The figure is "
    "therefore CONSERVATIVE for sample size and that is a stated limitation, not grounds to reclassify",
}

NOT_CLAIMED = [
    "not a demonstration that the reconstructed columns equal the published ones: they do not, and the "
    "gate that measured that failed. This asks only whether the endpoint is indifferent to the "
    "difference",
    "not a licence to use the reconstruction for anything else: a quantity that leaves one AUPRC "
    "difference unchanged is not thereby validated for a threshold, a percentile, a regression "
    "coefficient or any claim about activity itself",
    "not a replication of the K562 result, not independence of the genomic regions, not a "
    "cell-type-specific claim about astrocytes, and no claim at all until the comparison is registered "
    "before any score is computed -- the original's not_claimed list carries over in full",
    "not an astrocyte result of any kind: this lane computes none",
]


def published_activity(pairs: list[Any]) -> dict[str, Any]:
    """Where the log1p caveat bites, measured on the published K562 held-out pairs.

    Stated before the gate runs, because it is the reason the gate can fail and must not look like an
    excuse found afterwards.
    """
    a = sorted(math.sqrt(max(p.dhs, 0.0) * max(p.h3k27ac, 0.0)) for p in pairs)
    c = math.sqrt(0.7603 * 0.8080)
    shifts = sorted(math.log1p(c * x) - math.log1p(x) for x in a)

    def q(xs: list[float], p: float) -> float:
        return xs[min(len(xs) - 1, int(p * len(xs)))]

    return {
        "pairs": len(a),
        "implied_constant_activity_factor": c,
        "log_of_that_factor": math.log(c),
        "activity_quantiles": {
            "min": a[0],
            "p05": q(a, 0.05),
            "p25": q(a, 0.25),
            "median": q(a, 0.50),
            "p75": q(a, 0.75),
            "p95": q(a, 0.95),
            "max": a[-1],
        },
        "share_below": {str(t): sum(1 for x in a if x < t) / len(a) for t in (0.1, 0.5, 1.0, 3.0)},
        "shift_range": {"min": shifts[0], "median": q(shifts, 0.5), "max": shifts[-1]},
        "shift_spread": shifts[-1] - shifts[0],
        "reading": (
            "the shift log1p(cA) - log1p(A) would be the constant log c = -0.2436 if every pair had "
            "large activity. It does not: it runs from 0.0000 to -0.2402 across these pairs, so its "
            "SPREAD is as large as the constant it would otherwise be, and 37% of pairs lie below "
            "activity 1.0 where log1p is non-linear. A spread in the shift is exactly what moves a "
            "ranking, so the a priori argument that the endpoint is safe is NOT established by the "
            "algebra on this data. It has to be measured, and that is the gate"
        ),
    }


def manifest(pairs: list[Any]) -> dict[str, Any]:
    return {
        "sources": [
            {
                "accession": str(ORIGINAL),
                "version": f"the frozen original, sha256 {ORIGINAL_SHA256}; cited, not amended, and "
                "its failed calibration gate is not reopened",
            },
            {
                "accession": "EngreitzLab/CRISPR_comparison, EPCrisprBenchmark heldout_5_cell_types "
                "and training_K562 (Gschwind et al. 2025)",
                "version": "fetched 2026-09-16; pinned here by sha256",
                "url": crispri.BASE_URL,
            },
        ],
        "inputs": [
            mf.input_entry(
                ORIGINAL,
                partition=None,
                role="the claim, the power figure and the not_claimed list are carried from it "
                "unchanged; its sha256 is recorded and no term of it is moved",
            ),
            mf.input_entry(
                crispri.KNOWLEDGE / crispri.HELDOUT,
                partition=CRISPRI_SPLIT_OF[crispri.HELDOUT],
                role="the K562 held-out pairs and their published activity columns, read here only to "
                "state the size of the comparison and where log1p is non-linear; already-read "
                "development evidence",
            ),
        ],
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "max_abs_gain_shift": MAX_ABS_GAIN_SHIFT,
            "min_score_spearman": MIN_SCORE_SPEARMAN,
            "k562_heldout_pairs": len(pairs),
            "alphagenome_requests": 0,
            "model_requests": 0,
            "requests_sent": 0,
            "money_spent": 0,
            "astrocyte_values_computed": 0,
        },
        "exclusions": [
            "no gain is computed here, with published or reconstructed activity: this file is the "
            "design, written first",
            "no AlphaGenome request is sent and this script has no code path that could send one",
            "no astrocyte value of any kind is read or computed",
            "no rescaling factor is fitted, here or in the run this registers",
            "the failed calibration gate is not re-run and not reopened",
        ],
        "partitions": {
            "k562_heldout": "the endpoint comparison set, already-read development evidence",
            "k562_training": "where the frozen weights come from, unchanged between the two runs",
            "astrocyte": "not touched by this lane at all",
        },
        "code_cleanliness": mf.code_cleanliness(ENTRY, OWN_CODE, ROOT),
    }


def main() -> int:
    t0 = time.time()
    digest, _, _ = mf.sha256_of(ORIGINAL)
    if digest != ORIGINAL_SHA256:
        print(f"refused: {ORIGINAL} is sha256 {digest}, not the frozen {ORIGINAL_SHA256}")
        return 2

    with gzip.open(crispri.KNOWLEDGE / crispri.HELDOUT, "rt") as fh:
        pairs = [p for p in crispri.parse(fh) if p.cell == "K562"]
    print(f"{len(pairs)} K562 held-out pairs")
    caveat = published_activity(pairs)

    payload = {
        "status": "a registration, written before anything is computed. No gain computed, no score "
        "taken, no request sent, no money spent, no astrocyte value touched",
        "lane": "lane-astro",
        "registration": "AstroREG-2",
        "supersedes_nothing": "the original registration stands and its calibration gate stands "
        "FAILED. This is a separate registration with a different gate, not a revision of that one",
        "cites": {"registration": str(ORIGINAL), "sha256": ORIGINAL_SHA256},
        "why_a_different_gate": (
            "the registered endpoint is an AUPRC DIFFERENCE, which depends only on how pairs are "
            "RANKED. Activity enters the frozen model only through log_activity = log1p(sqrt(dhs * "
            "h3k27ac)) and activity_over_distance = log_activity - log(distance), both in log space, so "
            "a per-column factor becomes an additive shift on both features, and a shift common to "
            "every pair cannot move a ranking. The old gate measured the factor and refused it; this "
            "one measures whether the endpoint moves"
        ),
        "and_why_that_argument_is_not_enough": caveat["reading"],
        "gate": GATE,
        "history": HISTORY,
        "still_blind": BLIND,
        "claim": CLAIM_UNCHANGED,
        "log1p_caveat": caveat,
        "not_claimed": NOT_CLAIMED,
        "requests_sent": 0,
        "money_spent": 0,
        "spending_approval": "Albert's approval named the ORIGINAL registration by hash and cannot "
        "carry to this one. If this gate passes, a new approval is required before any of the 1,322 is "
        "spent, and this lane spends none either way",
        "left_undone": [
            "the endpoint comparison itself: both gains, their difference and the per-pair score "
            "correlation are still to be computed",
            "the residual scale factor between reconstructed and published activity is not explained "
            "and was not investigated; this registration asks whether it matters to the endpoint, not "
            "what it is",
            "nothing is established about astrocytes by this file or by the run it registers",
        ],
    }
    payload[mf.KEY] = manifest(pairs)
    path = save_result(RESULT, payload)

    print(
        json.dumps(
            {
                "max_abs_gain_shift": MAX_ABS_GAIN_SHIFT,
                "min_score_spearman": MIN_SCORE_SPEARMAN,
                "k562_heldout_pairs": len(pairs),
                "share_activity_below_1": caveat["share_below"]["1.0"],
                "shift_spread": round(caveat["shift_spread"], 4),
                "log_of_implied_factor": round(caveat["log_of_that_factor"], 4),
                "requests_sent": 0,
                "money_spent": 0,
            },
            indent=1,
        )
    )
    print(f"({time.time() - t0:.0f} s) -> {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
