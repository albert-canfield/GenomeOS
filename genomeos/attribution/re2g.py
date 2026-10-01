# SPDX-License-Identifier: AGPL-3.0-or-later
"""The paired comparison against ENCODE-rE2G on identical held-out CRISPR pairs, registered before
any comparator score was read.

Everything in this module that fixes the comparison -- the claim, the endpoint, the comparator and the
order it was preferred in, the join rule, the gate, the populations, the statistic, the reading, the
four independence statements and the falsifier -- was written and committed before any analysis code
ran, and is not to be changed afterwards. `scripts/re2g_comparator.py` records what could and could
not be obtained under the lane's download budget; it reads none of the quantities below.

The project's earlier figure against ENCODE-rE2G was a band ("where does our number sit against their
published interval"), not a paired comparison: the two models were never scored on the same resamples,
so no difference between them was ever measured. That band was withdrawn on 2026-09-27 for a separate
reason (docs/ROADMAP.md: the held-out comparison flatters this project, and the published held-out
model reads DNase only while the activity term here reads H3K27ac). This module is what would replace
it: a difference, with an interval, on identical pairs and identical resamples.

docs/ATTRIBUTION.md, "The paired comparison against ENCODE-rE2G, registered before any comparator
score was read".
"""

from __future__ import annotations

from collections.abc import Callable
from typing import Any

from genomeos.attribution import crispri

#: The one thing this comparison asserts. Nothing wider is registered here: no statement about the
#: deletion feature's value in general, and none about any pair set other than the one named.
CLAIM = (
    "On the Gschwind et al. 2026 held-out CRISPR benchmark (doi:10.1038/s41586-026-10781-4, "
    "Supplementary Table 3; the 4,378 pairs, 190 positives, 157.39 weighted, already described in "
    "data/results/crispri_published.json), the frozen 'activity + distance + deletion' model ranks "
    "regulated pairs differently from ENCODE-rE2G when both are scored on the identical pairs"
)

#: Weighted AUPRC by the benchmark's own estimator (`crispri.benchmark_auprc`), which reproduced the
#: published distance-to-TSS figures exactly (0.4359 unweighted on training, 0.3631 weighted on
#: held-out) before any model of this project was scored with it. The positive is the benchmark's own
#: `Regulated` column as published, and the weight its own `direct_vs_indirect_negative`. No
#: relabelling, no filtering, no re-derivation of either.
ENDPOINT = (
    "weighted AUPRC using the benchmark's own pair weights (direct_vs_indirect_negative), by the "
    "benchmark's own estimator; positive = the benchmark's Regulated label exactly as published, with "
    "no relabelling"
)

#: The order the comparator was preferred in, fixed before looking for any of them. The order exists so
#: that the cheapest and most faithful source wins on its merits and not on the number it produces.
COMPARATOR_PREFERENCE = (
    "1. a per-pair ENCODE-rE2G score carried in the benchmark table itself, if present",
    "2. otherwise the ENCODE portal's rE2G prediction file for each benchmark biosample -- model "
    "version, accession and sha256 recorded -- joined to the pairs by the benchmark's own overlap and "
    "gene rule (JOIN_RULE), with unpredicted pairs scored 0 as the benchmark itself does",
)

#: Preference 1 is settled by reading the column names of the held-out table and nothing else: the 42
#: columns of EPCrisprBenchmark_combined_data.heldout_5_cell_types.GRCh38.tsv.gz carry the measurement,
#: the power columns, the chromatin annotations and the weight, and no predictor's score. So preference
#: 2 applies, and the join rule below is the one that governs.
COMPARATOR_IN_TABLE = (
    "absent: the held-out table's 42 columns carry the measurement, the five power columns, the "
    "chromatin annotations and the weight, and no column holding any predictor's score, ENCODE-rE2G's "
    "included. Preference 1 is therefore unavailable and preference 2 governs"
)

#: The commit the join rule is read at: the tag v1.0.0, which the paper's data availability names, and
#: not the moving branch. `crisprComparisonMergeFunctions.R` is byte-identical at v1.0.0 and at main
#: (sha256 fee907d41ff9a72740dd9cf57ace795536fa9adb67f8ff462c94777a6fea916c, 506 lines), and
#: `createPredConfig.R` carries the same aggregate_function and fill_value at both, so the rule below is
#: the one the published comparison used.
BENCHMARK_COMMIT = "50587422e6b11259ead6fbc6f867681c788f39b7"
BENCHMARK_TAG = "v1.0.0"

#: The join, taken from the benchmark's own merge code and not invented here:
#: EngreitzLab/CRISPR_comparison, workflow/scripts/crisprComparisonMergeFunctions.R,
#: `combineSingleExptPred` (steps 1 to 3) with the defaults of `createPredConfig.R`. Stated in full
#: before any join was performed, because a join rule chosen after seeing a number is not a rule.
JOIN_RULE = {
    "key": (
        "a pair joins a prediction only within the same cell type and the same target gene: the "
        "experiment's (CellType, chrom, measuredGeneSymbol) against the prediction's "
        "(ExperimentCellType, chr, TargetGene), gene matched by symbol. This is the pipeline's own "
        "trick of folding cell type and gene into the GRanges seqname so that findOverlaps cannot "
        "cross either"
    ),
    "overlap": (
        "any overlap of at least one base between the perturbed element [chromStart, chromEnd] and the "
        "predicted element [start, end]. Both sides are read from BED-derived tables into the same "
        "interval convention, as the pipeline does, so the convention cancels"
    ),
    "several_overlaps": (
        "where one perturbed element overlaps several predicted elements, the score is aggregated by "
        "the pipeline's documented default aggregate_function, sum "
        "(workflow/scripts/createPredConfig.R). 'max' is registered here as a named sensitivity and is "
        "NOT what the gate is judged on: the gate is judged on sum alone, so that the aggregation is "
        "not chosen by which value passes"
    ),
    "unpredicted": (
        "a pair with no overlapping prediction takes the pipeline's fill_value, 0, and is kept in the "
        "pair set -- the benchmark's own rule (fillMissingPredictions, Prediction = 0). A sensitivity "
        "analysis drops the unpredicted pairs instead of scoring them 0, and is reported beside the "
        "primary figure, never in place of it"
    ),
    "biosample": (
        "one rE2G prediction set per benchmark cell type, named by the accession Supplementary Table "
        "12 gives for that sample; where that table lists several annotations for one cell type the "
        "ambiguity is recorded and the comparison is not run until it is resolved"
    ),
}

#: The gate. The published figure is the pooled one and the only one the source gives: Supplementary
#: Table 3, sheet "Held-out benchmarks", has six predictors times three metrics for a single "Held-out"
#: dataset and no per-cell-type held-out AUPRC at all. So the gate can only be run on all 4,378 pairs.
GATE_TARGET = 0.556151
GATE_INTERVAL = (0.467852, 0.631224)
GATE_TOLERANCE = 0.02

#: The comparator files preference 2 names, one per benchmark cell type, as Supplementary Table 12 gives
#: them: the annotation accession for that sample, the smallest full "element gene links" file in it, and
#: that file's size. The thresholded files in the same annotations are a tenth of a percent of the size
#: but hold only the pairs above the paper's 70%-recall threshold, so they cannot carry this endpoint:
#: scoring every other pair 0 censors the score distribution the AUPRC is taken over.
COMPARATOR_FILES = {
    "K562": {
        "annotation": "ENCSR627ANP",
        "sample": "K562_ENCSR000EOT",
        "file": "ENCFF970QAX",
        "bytes": 333743457,
    },
    "GM12878": {
        "annotation": "ENCSR604WPK",
        "sample": "GM12878_ENCSR000EMT",
        "file": "ENCFF260CIO",
        "bytes": 301527890,
    },
    "HCT116": {
        "annotation": "ENCSR575WDI",
        "sample": "HCT116_ENCSR000ENM",
        "file": "ENCFF407YWT",
        "bytes": 297785157,
    },
    "Jurkat": {
        "annotation": "ENCSR805BXG",
        "sample": "Jurkat__Clone_E6-1_ENCSR000EOS",
        "file": "ENCFF681NXF",
        "bytes": 304665637,
    },
    "WTC11": {
        "annotation": "ENCSR464GSY",
        "sample": "WTC11_ENCSR785ZUI",
        "file": "ENCFF071ZZU",
        "bytes": 331445001,
    },
}

#: The lane's download budget, and what the gate would cost against it. Recorded here because the gap is
#: the reason the comparison is registered and not run, and because a later lane must not quietly spend
#: more than was authorised.
DOWNLOAD_BUDGET_BYTES = 1_000_000_000
COMPARATOR_BYTES_REQUIRED = sum(v["bytes"] for v in COMPARATOR_FILES.values())

#: What a K562-only prediction file can attain as a pooled weighted AUPRC over all 4,378 held-out pairs,
#: with the other 2,460 pairs forced to 0 by the benchmark's own fill rule: no skill inside K562 at one
#: end, perfect separation inside K562 at the other. Computed from the benchmark's real labels and real
#: weights with synthetic scores, so no comparator score is involved; `scripts/re2g_k562_only_range.py`
#: reproduces it. The published target 0.556151 falls INSIDE this range, which is why the one file that
#: fits the download budget cannot serve the gate: passing it would not mean the scores are the
#: published ones. A test holds the bracket, so this cannot be forgotten.
K562_ONLY_RANGE = (0.064, 0.6341)

#: Why the gate is registered and not run. Availability is not the obstacle: every file is public on the
#: ENCODE portal, needs no login and costs nothing. The obstacle is that the only gate value the source
#: publishes is the pooled one, and the pooled figure cannot be reproduced from one cell type's
#: predictions.
BLOCKER = {
    "the_gate_is_pooled_only": (
        "Supplementary Table 3's held-out sheet gives one 'Held-out' dataset, six predictors times three "
        "metrics, and no per-cell-type held-out AUPRC. So there is no published figure to gate the K562 "
        "stratum against, and the gate can only be run on all 4,378 pairs"
    ),
    "the_pooled_gate_needs_all_five_cell_types": (
        "held-out K562 carries 100.79 of the 157.39 weighted positives, 64.0%. A K562-only prediction "
        "file forces the other 2,460 pairs to 0, so 36.0% of the positive weight goes into one tie at "
        "the bottom and the quantity computed is not the published one, which used predictions in all "
        "five cell types. The reason that disqualifies it is not that it must fail the gate but that it "
        f"could pass it for the wrong reason: on the pooled population a K562-only file can attain "
        f"anything from {K562_ONLY_RANGE[0]} (no skill inside K562) to {K562_ONLY_RANGE[1]} (perfect "
        f"separation inside K562), and that range brackets the published {GATE_TARGET:.4f}. A figure "
        "within 0.02 of the target would therefore be no evidence that the scores in hand are the "
        "published ones, which is the only thing the gate exists to establish"
    ),
    "bytes": (
        f"the five full prediction files total {COMPARATOR_BYTES_REQUIRED:,} bytes "
        f"({COMPARATOR_BYTES_REQUIRED / 1e9:.2f} GB), against an authorised budget of one public download "
        f"under {DOWNLOAD_BUDGET_BYTES / 1e9:.0f} GB. The single K562 file alone fits the budget but "
        "cannot serve the gate, for the reason above; the set that could serve the gate does not fit "
        "the budget"
    ),
    "hct116_is_ambiguous": (
        "Supplementary Table 12 lists 16 HCT116 rE2G annotations. HCT116_ENCSR000ENM is recorded above as "
        "the untreated ENCODE HCT116 DNase sample, but which one the benchmark used is not established "
        "from the table alone, and JOIN_RULE['biosample'] forbids running the comparison until it is"
    ),
    "not_a_paywall": (
        "no file here is behind a login or a payment. The paper's baseline-predictor bundle on Synapse "
        "(syn58896208) would need an account, and is not used for that reason"
    ),
}

GATE = (
    "before any comparison, reproduce ENCODE-rE2G's published held-out weighted AUPRC "
    f"({GATE_TARGET:.4f}, interval {GATE_INTERVAL[0]:.4f} to {GATE_INTERVAL[1]:.4f}, Supplementary "
    "Table 3, sheet 'Held-out benchmarks', pooled over all 4,378 pairs) from the per-pair comparator "
    f"scores. A miss of more than {GATE_TOLERANCE} stops the lane: the scores are not the published "
    "ones. The miss is recorded and nothing is compared"
)

#: Populations, with the count beside each name. Measured from the held-out table on 2026-10-01, before
#: any comparator score existed: pairs, positives, weighted positives.
POPULATIONS = {
    "primary": {
        "name": "held-out K562 pairs",
        "pairs": 1918,
        "positives": 118,
        "weighted_positives": 100.79,
        "deletion_feature": "available",
        "why": "the stratum where the deletion feature exists",
    },
    "secondary": {
        "name": "all held-out pairs, pooled",
        "pairs": 4378,
        "positives": 190,
        "weighted_positives": 157.39,
        "deletion_feature": "available on the K562 and GM12878 pairs only",
        "why": "the population the published comparator figure is stated on, and the only one it is",
    },
    "gm12878": {
        "name": "held-out GM12878 pairs",
        "pairs": 68,
        "positives": 16,
        "weighted_positives": 14.30,
        "deletion_feature": "available",
        "why": "its own stratum, reported apart and never pooled into the primary",
    },
}

#: The strata that may never carry a number for this comparison, because the deletion feature does not
#: exist in them. `crispri.gain_where_available` is what enforces it; these are the names it refuses.
FEATURE_UNAVAILABLE = ("HCT116", "Jurkat", "WTC11")

#: The statistic. One seed, shared by both models, so that every draw resamples the same clusters for
#: ours and for the comparator and the difference is paired within the draw.
DRAWS = 2000
SEED = 0
CLUSTER = "chromosome"
STATISTIC = (
    f"paired delta weighted AUPRC = ours minus ENCODE-rE2G, on identical resamples: {CLUSTER}-cluster "
    f"bootstrap, {DRAWS} draws requested, one seed ({SEED}) shared by both models so each draw is a "
    "paired difference. The result states clusters, draws_requested and draws_dropped "
    "(crispri._interval_provenance), because a dropped draw is a draw the interval does not rest on"
)

#: The reading, fixed in advance. These three strings are the only wording this comparison may be
#: reported in, and `reading()` is the only thing that chooses between them.
READS_BETTER = "ranks better than ENCODE-rE2G on these pairs"
NO_DIFFERENCE = "no difference detected"
READS_WORSE = "ranks worse"

#: What is true about exposure whatever the number comes out as. All four are stated in the result,
#: not left to a reader to infer, because each of them weakens the comparison in a different way.
INDEPENDENCE = {
    "a_comparator_trained_on_this_compendium": (
        "ENCODE-rE2G was trained on K562 CRISPR pairs from the same compendium, and these held-out "
        "pairs were held out by its own authors"
    ),
    "b_our_weights_frozen": (
        "our weights were fitted on the training split and frozen before the held-out pairs were scored"
    ),
    "c_reused_benchmark": (
        "this held-out set has already been scored by this project, so it is a reused benchmark: a "
        "comparison of two frozen models, not fresh validation"
    ),
    "d_deletion_feature_exposure": (
        "the deletion feature comes from a model trained on ENCODE K562 tracks (features, not labels). "
        "Whether its training saw these CRISPR outcomes is not established, and is not claimed either "
        "way"
    ),
}

#: What the comparison is not allowed to be used for. The second clause holds whatever the number is.
FALSIFIER = (
    "if the primary interval covers 0 or lies below it, README may not say the deletion model performs "
    "above ENCODE-rE2G. The earlier 'above the published interval' wording is withdrawn either way, "
    "and was in fact already withdrawn on 2026-09-27 for a separate reason"
)


def gate(value: float | None) -> dict[str, Any]:
    """Whether the comparator scores reproduce the published figure closely enough to be compared.

    A miss larger than `GATE_TOLERANCE` means the scores in hand are not the published ones, so the
    lane stops: `passed` is False and `may_compare` is False. `value` of None is a gate that could not
    be run at all, which is also not a pass.
    """
    if value is None:
        return {
            "target": GATE_TARGET,
            "published_interval": list(GATE_INTERVAL),
            "tolerance": GATE_TOLERANCE,
            "reproduced": None,
            "miss": None,
            "passed": False,
            "may_compare": False,
            "reason": "the gate was not run: no per-pair comparator score was obtained",
        }
    # The miss is judged on the same figure the result reports, so the verdict rests on the number a
    # reader can see rather than on a float representation a reader cannot. This decides the boundary
    # case only: a miss of exactly the tolerance passes, as the tolerance is written inclusively.
    miss = round(abs(value - GATE_TARGET), 4)
    passed = miss <= GATE_TOLERANCE
    return {
        "target": GATE_TARGET,
        "published_interval": list(GATE_INTERVAL),
        "tolerance": GATE_TOLERANCE,
        "reproduced": round(value, 4),
        "miss": miss,
        "passed": passed,
        "may_compare": passed,
        "reason": (
            "the comparator scores reproduce the published held-out figure within the tolerance"
            if passed
            else f"missed the published figure by {miss}, above {GATE_TOLERANCE}: the scores are "
            "not the published ones, so nothing is compared"
        ),
    }


def reading(ci95: list[float] | None) -> str | None:
    """The registered wording for an interval, and nothing else.

    None when there is no interval to read: a comparison without an interval has no reading, and is
    not allowed to borrow one from its point estimate.
    """
    if not ci95 or len(ci95) != 2:
        return None
    lo, hi = ci95
    if lo > 0:
        return READS_BETTER
    if hi < 0:
        return READS_WORSE
    return NO_DIFFERENCE


def paired_delta(
    ours: list[float],
    theirs: list[float],
    pairs: list[crispri.Pair],
    seed: int = SEED,
    draws: int = DRAWS,
) -> dict[str, Any]:
    """Ours minus ENCODE-rE2G on identical resamples, with the registered reading attached.

    Both models are resampled by the same `random.Random(seed)` walk over chromosomes inside
    `crispri.weighted_gain`, so each draw holds one paired difference rather than two independent ones.
    """
    out = crispri.weighted_gain(ours, theirs, pairs, weighted=True, seed=seed, n=draws)
    out["delta_auprc"] = out.pop("gain")
    out["reading"] = reading(out.get("ci95"))
    out["seed"] = seed
    out["paired"] = "identical resamples; one seed shared by both models"
    return out


def delta_where_available(cell: str, delta: Callable[[], dict[str, Any]]) -> dict[str, Any]:
    """The paired delta for a stratum, or a refusal naming why no number may stand there.

    A stratum in `FEATURE_UNAVAILABLE` has no deletion value, so the difference between the two models
    there would be an artefact of the model form rather than a measurement. It reports null and the
    reason, never a number -- the same rule `crispri.gain_where_available` applies to the deletion gain.
    """
    return crispri.gain_where_available(cell not in FEATURE_UNAVAILABLE, delta)
