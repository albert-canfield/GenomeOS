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


# --- applying the registered rule -------------------------------------------------------------------
# Everything above was committed before any comparator score was read. What follows implements the rule
# already registered above; it does not restate or alter it.

#: The columns the join needs from an ENCODE rE2G `element gene links` BED, by header name. The file
#: carries a `#chr`-prefixed header line, nineteen columns, with the model's own features between the
#: identity and the score; only the identity and the score are read.
PREDICTION_COLUMNS = ("#chr", "start", "end", "TargetGene", "CellType", "Score")


def _overlaps(a_start: int, a_end: int, b_start: int, b_end: int) -> bool:
    """The pipeline's own overlap test, which is what JOIN_RULE['overlap'] names.

    `combineSingleExptPred` puts both sides' raw BED numbers into `IRanges(start, end)` and calls
    `findOverlaps`, so the comparison is inclusive on both ends: the intervals meet when neither starts
    after the other finishes. Reproduced here rather than replaced by a half-open test, because the
    difference is one base at each boundary and the point is to join as the benchmark joined.
    """
    return a_start <= b_end and b_start <= a_end


def join_predictions(
    pairs: list[crispri.Pair],
    path: Any,
    cell: str,
    aggregate: str = "sum",
) -> tuple[dict[int, float], int]:
    """Scores for the pairs of one cell type, from one rE2G prediction file, by the registered rule.

    Returns the mapping of index in `pairs` to aggregated score and the number of overlapping prediction
    rows it was built from. The mapping holds only the pairs that joined: a pair absent from it took no
    prediction and is the caller's to fill with 0, the benchmark's own fill rule. `aggregate` is "sum"
    (the pipeline's documented default, and what the gate is judged on) or "max" (the registered
    sensitivity).
    """
    import gzip

    if aggregate not in ("sum", "max"):
        raise ValueError(f"aggregate must be sum or max, not {aggregate!r}")
    # (chrom, gene) -> the pairs of this cell type there, so each prediction row costs one dict lookup
    index: dict[tuple[str, str], list[int]] = {}
    for i, p in enumerate(pairs):
        if p.cell == cell:
            index.setdefault((p.chrom, p.gene), []).append(i)
    out: dict[int, float] = {}
    hits = 0
    with gzip.open(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        col = {name: header.index(name) for name in PREDICTION_COLUMNS}
        c_chr, c_start, c_end = col["#chr"], col["start"], col["end"]
        c_gene, c_score = col["TargetGene"], col["Score"]
        width = len(header)
        for line in fh:
            row = line.rstrip("\n").split("\t")
            if len(row) != width:
                continue
            candidates = index.get((row[c_chr], row[c_gene]))
            if not candidates:
                continue
            start, end = int(row[c_start]), int(row[c_end])
            score = float(row[c_score])
            for i in candidates:
                p = pairs[i]
                if not _overlaps(p.start, p.end, start, end):
                    continue
                hits += 1
                if i not in out:
                    out[i] = score
                elif aggregate == "sum":
                    out[i] += score
                else:
                    out[i] = max(out[i], score)
    return out, hits


# --- the second registration: like-for-like, written before anything was computed -------------------

#: A separate, like-for-like paired comparison, registered in full before any figure of it was computed.
#: The first comparison's activity term reads H3K27ac, which the benchmark's own positive selection is
#: total on, while the comparator's published held-out model reads DNase only. This one matches OUR
#: feature set to the comparator's input instead -- `dnase + distance + deletion` against ENCODE-rE2G --
#: so that the shared input is DNase in both and the selection-correlated input is in neither. Nothing
#: of the comparator is refitted, reweighted or adjusted: the like-for-like is achieved by changing our
#: features, never theirs.
SECOND_REGISTRATION = {
    "claim": (
        "On the Gschwind et al. 2026 held-out CRISPR benchmark, the frozen 'dnase + distance + "
        "deletion' model ranks regulated pairs differently from ENCODE-rE2G when both are scored on "
        "the identical pairs and both read DNase rather than H3K27ac"
    ),
    "why_it_exists": (
        "all 190 held-out positives carry H3K27ac at the tested element and 1,438 non-regulated pairs "
        "do not, so an H3K27ac-reading model is handed a separation by the benchmark's positive "
        "selection that ENCODE-rE2G's DNase-only held-out model cannot use. This comparison removes "
        "H3K27ac from our side. It is the comparison that may be cited about ENCODE-rE2G; the H3K27ac "
        "one is reported with its qualifiers and never cited alone"
    ),
    "models": {
        "ours": "crispri.DNASE_FEATURES['dnase + distance + deletion'] (log_distance, log_dnase, "
        "dnase_over_distance, top_target, deletion_drop), weights fitted on the covered training pairs "
        "and frozen before the held-out pairs are scored",
        "baseline_for_decomposition": "crispri.DNASE_FEATURES['dnase + distance']",
        "comparator": "the same reconstructed ENCODE-rE2G per-pair scores as the first comparison, "
        "unchanged: same five portal files, same join rule, same fill convention. Nothing of the "
        "comparator is refitted, reweighted or adjusted",
    },
    "endpoint": ENDPOINT,
    "population": (
        "primary and only registered population: the 1,918 held-out K562 pairs, where the deletion "
        "feature is that cell line's own. The pooled population is reported as descriptive context, "
        "with the same unavailable-feature refusal for HCT116, Jurkat and WTC11"
    ),
    "gate": "the same gate, already passed: the comparator scores reproduce the published pooled "
    f"weighted AUPRC {GATE_TARGET:.4f} within {GATE_TOLERANCE}. The comparator is not re-derived, so "
    "the gate is not re-run; its verdict carries over because the scores are identical",
    "statistic": STATISTIC,
    "reading": {
        "lower_bound_above_zero": READS_BETTER,
        "interval_covers_zero": NO_DIFFERENCE,
        "upper_bound_below_zero": READS_WORSE,
        "no_other_wording": "the same three strings, chosen by `reading()` and nothing else",
    },
    "prior_exposure": {
        "stated_because_it_is_real": (
            "unpaired figures for this variant already exist in the committed "
            "crispri_published.json post_hoc_positive_filter.dnase_only_diagnostic and have been seen, "
            "so this registration is NOT blind. It is a pre-specification of a paired statistic that "
            "does not exist yet, not a first look at the variant"
        ),
        "heldout_pooled_weighted_already_seen": {
            "dnase + distance": 0.4757,
            "dnase + distance + deletion": 0.6393,
            "deletion_gain": {"gain": 0.1636, "ci95": [0.1015, 0.2368], "resamples": 200},
            "against_encode_re2g": "the file already carries the phrase 'above the published "
            "interval', unpaired, against the published 0.5562",
        },
        "heldout_k562_already_seen": {
            "note": "the coordinator's brief said no K562 figure exists for this variant; one does, and "
            "understating the exposure would be wrong. It is not the registered quantity: it covers "
            "1,744 covered pairs rather than all 1,918, is unweighted, and uses the average-precision "
            "estimator rather than the benchmark's",
            "pairs": 1744,
            "positives": 114,
            "dnase + distance": 0.5021,
            "dnase + distance + deletion": 0.6813,
            "deletion_gain_gain": 0.1792,
        },
        "what_does_not_exist": (
            "no paired delta against ENCODE-rE2G for this variant, in any population; and no weighted "
            "benchmark-estimator figure on all 1,918 K562 pairs. Those are what this registration fixes"
        ),
        "the_old_figures_are_exposure_not_findings": (
            "the 200-resample intervals above are the old coarse ones. Everything reported under this "
            f"registration is recomputed at {DRAWS} draws; the figures above are never quoted as "
            "results of it"
        ),
    },
    "falsifier": (
        "if the K562 interval covers 0 or lies below it, README may not say the deletion model ranks "
        "better than ENCODE-rE2G at all, because this is the like-for-like comparison and the H3K27ac "
        "one may not be cited alone"
    ),
    "independence_and_exposure": INDEPENDENCE,
    "also_true_of_this_comparison": (
        "DNase remains a shared input, read for both models from the benchmark's own table for ours and "
        "from the ENCODE portal predictions for theirs; the reconstruction still falls "
        f"{GATE_TOLERANCE} or less short of the published figure, and the actual shortfall is reported "
        "with the result; and this is still a reused benchmark and a comparison of two frozen models, "
        "not fresh validation"
    ),
}
