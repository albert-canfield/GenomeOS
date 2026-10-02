# SPDX-License-Identifier: AGPL-3.0-or-later
"""Registration of the astrocyte CRISPRi test: every term fixed before any score exists.

The eligibility lane (`genomeos/attribution/fresh.py`, result `fresh_crispri_eligibility.json`)
established that Green et al. 2025's CRISPRi screen in cultured primary human astrocytes clears the
floors set before it was opened: 7,759 element-gene pairs, 158 measured positives over 86 independent
loci on 22 chromosomes. This module registers what would be done with it. Nothing here reads an
effect size, and no model request is made by anything in this module.

Why a registration before a run. The project's CRISPRi result is a K562 result. Every attempt to
find a second cell type has either reused a set already scored here or been too small to carry an
interval. This is the first fresh, non-K562, adequately sized set, and it is therefore the only
chance to be wrong in public: if the design is written after the numbers, the test cannot fail. So
the labels, the claim, the interval, the power rule and the readings are fixed here, in code, and the
readings travel word for word into every later quote.

The one substantive thing this module refuses. The sweep's per-element cache
(`crispri.ElementCache`) keeps four cell lines' values per gene -- GM12878, HepG2, IMR-90, K562 --
and astrocyte is not among them, so there is no astrocyte deletion value to reuse. `astrocyte` does
survive in the cache as a *winning track name*, on 224 of chr21's 348,033 gene rows. Those 224 are
not reusable, and the reason is registered here rather than discovered later: a value survives there
only because it was the most extreme of the 371 tracks in its window, so reusing them would import
selection on a correlate of the outcome.
"""

from __future__ import annotations

import random
from collections import defaultdict
from typing import Any

from genomeos.attribution import cell2, crispri

#: The eligibility lane counted 158 HITS over 86 independent loci. The label rule registered below is
#: narrower than "hit": a positive is a significant *decrease*, the benchmark's own construction, so
#: 25 of those 158 are increases and are held apart. The positive class is therefore 133 pairs over
#: 74 loci, and **74 is what the interval and the power calculation are sized against** -- not 86.
#: The two figures are both kept here so that nothing quietly inherits the wider one: a number
#: computed over all hits may not be quoted for a test run on the decreases.
ELIGIBILITY_HITS = 158
ELIGIBILITY_LOCI = 86
REGISTERED_POSITIVES = 133
REGISTERED_INCREASES_HELD_APART = 25
REGISTERED_LOCI = 74
REGISTERED_PAIRS = 7_759
REGISTERED_CHROMOSOMES = 22

#: The reconciliation, stated rather than left to be noticed. Both floors still clear: 133 positives
#: against a floor of 30, and 74 loci against a floor of 20.
LOCUS_RECONCILIATION = (
    f"the eligibility lane reported {ELIGIBILITY_HITS} positives over {ELIGIBILITY_LOCI} independent "
    f"loci, counting every pair the authors call a hit. The label rule registered here is narrower: a "
    f"positive is a significant DECREASE, so {REGISTERED_INCREASES_HELD_APART} increases are held "
    f"apart and the positive class is {REGISTERED_POSITIVES} pairs over {REGISTERED_LOCI} loci on "
    f"{REGISTERED_CHROMOSOMES} chromosomes. Both floors still clear ({REGISTERED_POSITIVES} against "
    f"30, {REGISTERED_LOCI} against 20). Every interval and the power figure are sized on "
    f"{REGISTERED_LOCI}; {ELIGIBILITY_LOCI} is the eligibility figure and may not be quoted for this "
    "test"
)

#: AlphaGenome's astrocyte output exists. Established at 0 requests from the sweep's own cache, where
#: `astrocyte` appears as an exact winning-track name; it is therefore a track the model returned,
#: not a track this lane hopes for. What does NOT exist is a retained astrocyte value per gene.
ASTROCYTE_TRACK_EXPOSED = (
    "AlphaGenome returned an astrocyte output: `astrocyte` is an exact winning-track name in the "
    "sweep's own per-element cache, on 224 of chr21's 348,033 gene rows (0.0644%), over 218 of its "
    "12,158 elements (1.79%), among 316 distinct winning names and 371 tracks scored per window. "
    "crispri.MODEL_CELLS governs which four cell lines' values are RETAINED per gene, not which "
    "tracks the model exposes: the astrocyte values were computed and discarded at write time, the "
    "same defect as the per-cell collapse fixed at 6c44443. So the earlier no-go was about the wrong "
    "thing -- the input is not missing from the model, it is missing from the cache"
)

#: Why the 224 surviving astrocyte rows may not be reused as the astrocyte deletion value. Registered
#: in these terms before any of them is looked at.
CACHE_REUSE_REFUSED = (
    "a value survives there only because it was the most extreme of the 371 tracks in its window, so "
    "reusing them would import selection on a correlate of the outcome. The 224 rows are counted as "
    "evidence that the track exists and are never read as values; every astrocyte deletion value "
    "used by the test must come from a fresh request on the same code path"
)

#: The label rules, fixed before any effect, p-value or fold change from the screen is read. They
#: follow the benchmark's own construction (`Regulated` = significant AND effect < 0) so that the
#: positive class means in astrocytes what it means in the K562 pairs the weights were frozen on.
LABELS = {
    "positive": (
        "the authors' own call of a significant DECREASE: Supplementary Table 3F's "
        "`Hit&Downregulated`. This is the benchmark's construction (significant AND effect below "
        "zero), so the positive class means the same thing in both cell contexts"
    ),
    "primary_negative": (
        "a pair the authors do NOT call a hit AND which they flag `WellPowered_at_FC_0.25`. The 0.25 "
        "power convention is the one the frozen model was scored under on the K562 benchmark, so an "
        "underpowered non-hit is not evidence of absence and is excluded rather than counted as a "
        "negative"
    ),
    "sensitivity_negative": (
        "the same rule at `WellPowered_at_FC_0.15`, reported as a sensitivity analysis beside the "
        "primary and never in place of it"
    ),
    "increases_held_apart": (
        "a pair the authors call a hit whose change is an INCREASE is neither a positive nor a "
        "negative. It is counted and reported separately and is never folded into the negatives, "
        "because an enhancer whose silencing raises a gene is a different phenomenon from one the "
        "model is asked to rank"
    ),
}

#: The primary claim. Frozen weights, nothing refitted, and the arms are the two the K562 result
#: already compares.
PRIMARY_CLAIM = (
    "the deletion gain in cultured primary human astrocytes: AUPRC of 'activity + distance + "
    "deletion' minus AUPRC of 'activity + distance', both scored with the weights frozen on the K562 "
    "training pairs. Nothing is refitted on astrocyte data -- not the weights, not a threshold, not "
    "a feature definition. The gain is the only primary quantity; an AUPRC on its own is not a claim"
)

#: The interval. Two clusterings are computed and the WIDER is the one quoted, so that whichever
#: grouping is the weaker guarantee governs the published tail.
INTERVAL_RULE = (
    f"a cluster bootstrap at {crispri.BOOTSTRAPS} draws over the {REGISTERED_LOCI} registered "
    f"independent loci, AND the chromosome bootstrap over {REGISTERED_CHROMOSOMES} clusters. Both "
    "are reported and the WIDER of the two is the interval that is quoted, because the weaker "
    "guarantee must govern the published tail. crispri.MIN_CLUSTERS_FOR_AN_INTERVAL applies to each: "
    f"below {crispri.MIN_CLUSTERS_FOR_AN_INTERVAL} clusters there is no ci95 at all, only "
    "'interval unreliable: N clusters'. The locus grouping is cell2's operational convention and is "
    "not established biological independence, so the chromosome interval is the more conservative of "
    "the two and is expected to be the wider"
)

#: The readings, fixed now so that no outcome can choose its own wording later.
READINGS = {
    "lower_bound_above_zero": "the deletion gain is detected in cultured fetal astrocytes",
    "interval_covers_zero": "no gain detected",
    "upper_bound_below_zero": "the deletion feature lowers ranking here",
    "never": (
        "an interval covering zero is never reported as 'no effect': it is 'no gain detected'. The "
        "distinction is the whole difference between an absent effect and an underpowered test"
    ),
}

#: The power rule, by a registered method, computed before any astrocyte score exists. The S8 classes
#: decide what the result may be called, and the class travels with every quote of it.
POWER_RULE = (
    f"resample the K562 held-out loci with replacement up to {REGISTERED_LOCI} loci, at K562's own "
    "observed gain, and report the share of resamples whose interval excludes zero. Classes from S8, "
    "fixed here: at or above 0.8 the test is confirmatory; 0.5 to 0.8 it is exploratory and is "
    "labelled so in every quote of it; below 0.5 the test is REFUSED and not run. The figure is "
    "computed before any astrocyte pair is scored, so it cannot be chosen after an outcome"
)

#: Amendment 1, registered 2026-10-02 before the recomputation it describes, after the reviewer found
#: the first power figure optimistic in two independent ways. The first figure, 0.96, kept K562's own
#: prevalence and assumed the whole K562 gain transfers; neither holds for the test actually proposed.
POWER_AMENDMENT_1 = (
    "the first power figure (0.96) was optimistic in two independent ways and is superseded, not "
    "revised. (a) PREVALENCE. It resampled K562 at K562's prevalence, 114 positives in 1,744 pairs "
    "= 6.5%. The astrocyte test runs at 133 in 4,580 = 2.9%, less than half. An AUPRC gain is "
    "prevalence-dependent -- both the attainable gain and its spread move with the positive rate -- "
    "so a figure computed at 6.5% does not describe a test run at 2.9%. Each resample is therefore "
    "matched to the astrocyte prevalence by downsampling its POSITIVES, keeping every negative, "
    "which is the direction that cannot invent data. (b) EFFECT SIZE. It assumed the full K562 gain "
    "transfers, the most favourable assumption available. Power is therefore also reported at 0.5x "
    "and 0.25x, by shrinking the deletion arm's score toward the no-deletion arm, "
    "score(k) = without + k * (with - without), and the ACHIEVED gain is reported at each k so that "
    "the attenuation is visible rather than asserted. "
    "The S8 class is set by the prevalence-matched figure at the FULL effect, and the 0.5x figure "
    "travels beside it in every quote of the class. Below 0.5 the test is refused and nothing is "
    "spent. This method was registered before it was run"
)

#: The two prevalences the amendment turns on, fixed from the label counts before the recomputation.
ASTROCYTE_PREVALENCE = 133 / (133 + 4447)
K562_PREVALENCE_NOTE = (
    "K562's own prevalence is computed from the covered held-out K562 pairs at run time and reported "
    "beside the astrocyte figure, so the size of the mismatch is on the record and not implied"
)

#: The attenuation factors reported. 1.0 sets the class; 0.5 travels beside it in every quote.
ATTENUATIONS = (1.0, 0.5, 0.25)

#: What the attenuation arm does NOT establish. Found on running it, recorded rather than presented as
#: a result: the scheme shrinks the signal and its own sampling spread in the same proportion, so the
#: share of intervals excluding zero barely moves while the achieved gain falls as intended. A reader
#: must not take the 0.5x share as evidence that power survives a halved true effect.
ATTENUATION_LIMITATION = (
    "the attenuation reduces the achieved gain as intended, but it does NOT answer the question it "
    "was asked. Shrinking the deletion arm toward the no-deletion arm scales the score separation, "
    "and the cluster bootstrap's spread scales with it, so the signal-to-noise ratio is close to "
    "preserved and the share of intervals excluding zero stays nearly flat across 1x, 0.5x and 0.25x "
    "even though the median achieved gain falls. The arm therefore shows that the DETECTABILITY of a "
    "proportionally smaller separation is similar; it does NOT show that power survives a weaker true "
    "effect against unchanged measurement noise, which is what a conventional effect-size sensitivity "
    "means. The 0.5x share may not be quoted as power at half the effect. A scheme that held the "
    "noise fixed while lowering the signal would be needed, and is not registered or run here"
)

#: Term 3: what the deletion value does if `astrocyte` turns out to be several AlphaGenome tracks.
#: Registered now so that no rule is invented at scoring time, when an outcome is visible.
TRACK_ROSTER_RULE = (
    "the astrocyte deletion value uses the frozen K562 rule UNCHANGED -- whatever the legacy "
    "per-cell-name collapse did for K562, it does for astrocyte -- because the frozen K562 values "
    "came through that same rule and the two definitions must match. Whether `astrocyte` is one "
    "AlphaGenome track or several is not known at registration time: the cache records a winning "
    "track name, not the roster. So the FULL astrocyte track roster from the first response is "
    "written to the result BEFORE any pair is scored. If it is several tracks, that fact is "
    "recorded; it is not resolved by a new rule invented after the data are in view"
)

#: Term 2: the activity term's inputs, provenanced. Read from the data files' own headers and
#: confirmed against the ENCODE portal; both experiments are Homo sapiens astrocyte.
ACTIVITY_INPUTS = {
    "h3k27ac": {
        "file_accession": "ENCFF970DKF",
        "experiment": "ENCSR000AOQ",
        "assay": "ChIP-seq",
        "target": "H3K27ac-human",
        "biosample": "Homo sapiens astrocyte",
        "output_type": "replicated peaks",
        "file_type": "bed narrowPeak",
        "assembly": "GRCh38",
        "status": "released",
        "date_created": "2020-09-30T04:40:52.048936+00:00",
        "portal_md5sum": "67285e9f88b12069df30875dbfe016d3",
        "portal_file_size": 1_639_983,
        "biological_replicates": [1, 2],
        "local_derivation": "data/knowledge/epigenome/peaks/astrocyte_H3K27ac_chr*.bed.gz, 24 "
        "chromosomes, each carrying the accession in its own first line",
    },
    "dnase": {
        "file_accession": "ENCFF874OPW",
        "experiment": "ENCSR000EPM",
        "assay": "DNase-seq",
        "target": None,
        "biosample": "Homo sapiens astrocyte",
        "output_type": "peaks",
        "file_type": "bed narrowPeak",
        "assembly": "GRCh38",
        "status": "released",
        "date_created": "2020-11-18T18:31:31.509548+00:00",
        "portal_md5sum": "12869af50a461f22f4101a7bde0fe12d",
        "portal_file_size": 3_221_892,
        "biological_replicates": [1],
        "local_derivation": "data/results/dnase_astrocyte_chr*.bed.gz, 24 chromosomes, each "
        "carrying the accession in its own first line",
    },
    "correction": (
        "an earlier left-undone item of this lane said these files carried no accession anywhere in "
        "the repository. That was wrong: each per-chromosome file records its ENCODE accession in its "
        "own first line. What was missing was a manifest entry, not the provenance"
    ),
    "why_it_matters": (
        "the activity term is half the model: 'activity + distance' is one of the two arms the gain "
        "is a difference between, so an unprovenanced activity input would leave the comparison "
        "unauditable on the side that is not even being tested"
    ),
}

#: Everything already read that touches this comparison. Stated so that a reader can discount the
#: result correctly rather than being told it is clean.
EXPOSURE = (
    "the screen's labels have never been read here: established by grep over docs/, genomeos/, "
    "scripts/, tests/ and the knowledge READMEs for Voineagu, AstroREG, PsychENCODE, EGrf, CROP-seq, "
    "NHA, the DOI stem and PMID 41413662, with no hit, and the study is absent from the project's "
    "CRISPRi exposure ledger",
    "astrocyte CHROMATIN features have already been read here: the node reader lanes hold astrocyte "
    "DNase for 24 chromosomes and the epigenome store holds astrocyte H3K27ac for 24 chromosomes. "
    "These are the activity term's own inputs, so the activity arm is not naive to this cell type",
    "AlphaGenome was trained on ENCODE tracks genome-wide, which include astrocyte chromatin. The "
    "labels are unseen; the genomic regions are not",
    "the elements were selected on PsychENCODE astrocyte activity, which flatters the activity term "
    "in both arms equally and so leaves the gain comparison less affected",
    "the frozen K562 values came through the same legacy per-cell collapse as this arm would (one "
    "track per cell name, pre-6c44443), so the two definitions match and the identity check below is "
    "the precedent for saying so",
)

#: Carried word for word from the eligibility lane's own text. Not to be shortened or softened.
NOT_CLAIMED = (
    "not a replication of the K562 result: a different assay readout (single-cell CROP-seq against "
    "the benchmark's mixed FlowFISH and scRNA-seq screens), a different element selection "
    "(PsychENCODE astrocyte enhancers), a different hit threshold and a different cell lineage, so a "
    "difference cannot be attributed to the cell type",
    "not independence of the genomic regions: the locus grouping is cell2's operational convention, "
    "and AlphaGenome was trained on ENCODE tracks genome-wide, which include astrocyte chromatin; "
    "the labels are unseen, the regions are not",
    "not a cell-type-specific claim about astrocytes: two donor lines of cultured fetal-derived "
    "astrocytes are not primary brain tissue",
    "no claim at all until the comparison is registered before any score is computed",
)

#: The same code path, term (a). The claim is checkable by reading one line, and is also tested.
SAME_CODE_PATH = (
    "the astrocyte deletion value is computed by exactly the function that produced the frozen K562 "
    "values: crispri.annotate with `cells` extended to include astrocyte, nothing else changed, "
    "weights frozen. crispri._annotate_one reaches the deletion features through "
    "`if p.cell in cells`, so `cells` is a membership gate and nothing more: extending the tuple "
    "cannot alter the value of any pair whose cell was already in it. That is why adding astrocyte "
    "leaves every K562 pair untouched, and the identity check recomputes 20 K562 pairs through the "
    "extended path to show it rather than assert it"
)

#: The build, settled in the eligibility lane and carried here so no join can start without it.
BUILD = (
    "GRCh38/hg38, confirmed against the project's own GENCODE table on five of the screen's measured "
    "genes before any join: HSPB1 within 8 bp of the annotated TSS, FTH1 6 bp, ATP2B4 102 bp, with "
    "RHOC and CRABP2 inside the gene body, against roughly 500 kb of displacement had it been hg19. "
    "No liftover is needed and no pair is lost to one"
)

#: What the money buys, in the terms the decision is made in.
COST_FRAMING = (
    "the figure Albert sees is the cost of the information, not the volume: the total requests "
    "needed, the covered share of positives and of negatives, and the expected requests per "
    "positive, each quoted against the HCT116 precedent of 705 requests for 363 covered pairs. "
    f"{REGISTERED_POSITIVES} positives over {REGISTERED_LOCI} loci is what the money buys"
)

#: The HCT116 arm, for the comparison the cost is quoted against.
HCT116_PRECEDENT = {"requests": 705, "covered_pairs": 363, "regulated": 34}


def cluster_bootstrap(
    a: list[float],
    b: list[float],
    labels: list[bool],
    clusters: list[Any],
    seed: int = 0,
    n: int = crispri.BOOTSTRAPS,
) -> dict[str, Any]:
    """AUPRC(`a`) - AUPRC(`b`), with a 95% interval from resampling whole `clusters`.

    `crispri.gain_interval` resamples chromosomes and is left untouched, because other lanes are in
    that file; this is the same procedure over an arbitrary clustering, so the locus interval and the
    chromosome interval are produced by one function and differ only in what a cluster is.

    Below `crispri.MIN_CLUSTERS_FOR_AN_INTERVAL` clusters no `ci95` is returned at all: the key is
    absent and `unreliable` says how many clusters there were. A caller cannot mistake a wide
    interval for a refused one.
    """
    members: dict[Any, list[int]] = defaultdict(list)
    for i, c in enumerate(clusters):
        members[c].append(i)
    keys = sorted(members, key=str)
    point = (crispri.average_precision(a, labels) or 0) - (crispri.average_precision(b, labels) or 0)
    out: dict[str, Any] = {
        "point": point,
        "clusters": len(keys),
        "requested_draws": n,
        "enough_clusters": len(keys) >= crispri.MIN_CLUSTERS_FOR_AN_INTERVAL,
    }
    if len(keys) < crispri.MIN_CLUSTERS_FOR_AN_INTERVAL:
        out["unreliable"] = f"interval unreliable: {len(keys)} clusters"
        return out
    rng = random.Random(seed)
    diffs = []
    for _ in range(n):
        idx = [i for c in (rng.choice(keys) for _ in keys) for i in members[c]]
        lab = [labels[i] for i in idx]
        if not any(lab):
            continue
        diffs.append(
            (crispri.average_precision([a[i] for i in idx], lab) or 0)
            - (crispri.average_precision([b[i] for i in idx], lab) or 0)
        )
    diffs.sort()
    out["kept_draws"] = len(diffs)
    out["dropped_draws_with_no_positive"] = n - len(diffs)
    if len(diffs) >= crispri.MIN_CLUSTERS_FOR_AN_INTERVAL:
        lo = diffs[int(0.025 * (len(diffs) - 1))]
        hi = diffs[int(0.975 * (len(diffs) - 1))]
        out["ci95"] = [lo, hi]
    else:
        out["unreliable"] = f"interval unreliable: {len(diffs)} usable draws"
    return out


def wider(*intervals: dict[str, Any]) -> dict[str, Any]:
    """Of the intervals that have a `ci95`, the one with the widest span; the registered rule.

    An interval with no `ci95` cannot be compared on width, so it cannot win on width either; if none
    of them has one, the first is returned so the refusal travels rather than being dropped.
    """
    with_ci = [i for i in intervals if "ci95" in i]
    if not with_ci:
        return dict(intervals[0])
    return dict(max(with_ci, key=lambda i: i["ci95"][1] - i["ci95"][0]))


def reading(interval: dict[str, Any]) -> str:
    """The registered reading of an interval, word for word. No other wording is permitted."""
    if "ci95" not in interval:
        return interval.get("unreliable", "interval unreliable")
    lo, hi = interval["ci95"]
    if lo > 0:
        return READINGS["lower_bound_above_zero"]
    if hi < 0:
        return READINGS["upper_bound_below_zero"]
    return READINGS["interval_covers_zero"]


def power_class(share: float) -> dict[str, Any]:
    """The S8 class of a power share, and whether the test may be run at all."""
    if share >= 0.8:
        return {"share": share, "class": "confirmatory", "run": True}
    if share >= 0.5:
        return {
            "share": share,
            "class": "exploratory",
            "run": True,
            "label_required": "exploratory, and labelled so in every quote of it",
        }
    return {
        "share": share,
        "class": "refused",
        "run": False,
        "why": "below 0.5 the test is refused and not run, by the rule registered before the figure",
    }


def attenuate(with_d: list[float], without: list[float], k: float) -> list[float]:
    """The deletion arm shrunk toward the no-deletion arm by `k`: `without + k * (with - without)`.

    At k=1 it is the deletion arm unchanged; at k=0 the two arms are identical and the gain is nil.
    This attenuates the *score separation*, which is why the achieved AUPRC gain is reported at each
    k rather than assumed to scale with it.
    """
    return [b + k * (a - b) for a, b in zip(with_d, without, strict=True)]


def match_prevalence(indices: list[int], labels: list[bool], target: float, rng: random.Random) -> list[int]:
    """`indices` thinned to `target` prevalence by dropping POSITIVES, keeping every negative.

    Downsampling positives is the only direction that cannot invent data: raising the negative count
    would mean resampling negatives that are not there. If the subset is already at or below the
    target prevalence it is returned unchanged, because thinning further would overshoot.
    """
    pos = [i for i in indices if labels[i]]
    neg = [i for i in indices if not labels[i]]
    if not neg or not pos:
        return list(indices)
    want = int(round(target * len(neg) / (1.0 - target)))
    if want >= len(pos):
        return list(indices)
    kept = rng.sample(pos, max(want, 0))
    return kept + neg


def label_of(hit: bool, downregulated: bool, well_powered: bool) -> str:
    """One pair's registered class: positive, increase_held_apart, negative or excluded."""
    if hit and downregulated:
        return "positive"
    if hit and not downregulated:
        return "increase_held_apart"
    if well_powered:
        return "negative"
    return "excluded_underpowered"


def terms() -> dict[str, Any]:
    """Every registered term, for the result file. This is the document Albert rules on."""
    return {
        "what_is_registered": (
            "the design of a first test of the frozen model's deletion feature in a cell context "
            "this project has never scored, on a set whose outcomes it has never read. Written "
            "before any astrocyte score exists and before any request is bought"
        ),
        "eligibility": {
            "source": "data/results/fresh_crispri_eligibility.json",
            "pairs": REGISTERED_PAIRS,
            "eligibility_hits": ELIGIBILITY_HITS,
            "eligibility_loci": ELIGIBILITY_LOCI,
            "registered_positives": REGISTERED_POSITIVES,
            "registered_increases_held_apart": REGISTERED_INCREASES_HELD_APART,
            "registered_loci": REGISTERED_LOCI,
            "chromosomes": REGISTERED_CHROMOSOMES,
            "locus_reconciliation": LOCUS_RECONCILIATION,
            "locus_rule": cell2.INDEPENDENT_LOCUS_RULE,
            "locus_rule_is_not": "established biological independence",
        },
        "feasibility_gate": {
            "astrocyte_track_exposed": ASTROCYTE_TRACK_EXPOSED,
            "cache_reuse_refused": CACHE_REUSE_REFUSED,
            "rule": (
                "if an input the frozen model needs is missing, the model is not applied and no "
                "substitute is put in its place: a no-go is recorded instead"
            ),
        },
        "build": BUILD,
        "labels": LABELS,
        "primary_claim": PRIMARY_CLAIM,
        "interval_rule": INTERVAL_RULE,
        "readings": READINGS,
        "power_rule": POWER_RULE,
        "power_amendment_1": POWER_AMENDMENT_1,
        "astrocyte_prevalence": ASTROCYTE_PREVALENCE,
        "k562_prevalence_note": K562_PREVALENCE_NOTE,
        "attenuations_reported": list(ATTENUATIONS),
        "attenuation_limitation": ATTENUATION_LIMITATION,
        "track_roster_rule": TRACK_ROSTER_RULE,
        "activity_inputs": ACTIVITY_INPUTS,
        "same_code_path": SAME_CODE_PATH,
        "exposure": list(EXPOSURE),
        "not_claimed": list(NOT_CLAIMED),
        "cost_framing": COST_FRAMING,
        "hct116_precedent": HCT116_PRECEDENT,
        "requests_sent_by_this_lane": 0,
        "money_spent_by_this_lane": 0,
    }
