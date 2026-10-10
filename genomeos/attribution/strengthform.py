# SPDX-License-Identifier: AGPL-3.0-or-later
"""lane-strengthform: the strength-dependent activator-combination form, on held-out enhancer pairs.

WHAT R3 LEFT OPEN. On 2026-09-28 lane-combine registered a choice between four fixed activator
rules and ran it once (`7f3e597` registered, `93caf61` result): no rule survived the measured
enhancer pairs of Bothma et al. 2015, the runtime kept the mean, and the choice was recorded as
undetermined. The reason given was that the measured pairs change with enhancer strength, which no
fixed function of the terms reproduces, and a STRENGTH-DEPENDENT form was named there and never
tested. `runtime/grn.py` already carries it as the opt-in rule `competition`, A = W / (K + W)
(Bothma et al. 2015, eLife 4:e07956, Methods eqs. 7-9): weak engagement adds, strong engagement
saturates. This module tests that form, and nothing in the engine is changed to do it.

WHY NOT ON LANE-COMBINE'S OWN PAIRS. The five Bothma pairs' outcomes were read by lane-combine
before it wrote its rule, and it said so. A form chosen in knowledge of those outcomes is not
independent of them, so this lane tests on pairs lane-combine never read: the three-way
combinatorial reporter assay of Martinez-Ara M, Comoglio F, van Steensel B (2024), "Large-scale
analysis of the integration of enhancer-enhancer signals by promoters", eLife 12:e91994
(doi 10.7554/eLife.91994), GEO GSE240586, in mouse E14tg2a embryonic stem cells. Its processed table
holds every enhancer-enhancer-promoter combination of 59 enhancers and 20 control fragments over 8
promoters. Neither the publication nor the accession appears anywhere in docs/ATTRIBUTION.md or
docs/ROADMAP.md (checked by search before the registration), so nothing of it has been scored here.

EXPOSURE, STATED RATHER THAN CLAIMED AWAY (`EXPLORATION_PRECEDED_THIS`). Before the form was fixed
this lane read the paper's own summary of its finding, read the Methods, and read from the table its
DESIGN columns, its promoter baselines, its control-control rows and its single-enhancer boost
indices. It read no activity of any enhancer-enhancer combination. The form's shape and its one
declared constant come from Bothma et al. 2015 and from 7d's committed THETA, not from this dataset.

WHAT IS DECLARED AND NOT FITTED. The form has one free number per promoter, the transcriptional
efficiency r with the promoter fully engaged. It is fixed by the same declaration 7d made
(`attribution/combine2.THETA` = 0.5): the promoter's STRONGEST SINGLE enhancer engages it half the
time, so with K = 1 that enhancer's term is 1 and r = 2 * (S_max - B). Every other single's term
then follows from its own measured boost index, and the pair prediction follows from the two terms.
Nothing is fitted on any pair, and the form reproduces each single it is given exactly.

The comparators are the rules lane-combine compared, each as the parameter-free prediction its own
algebra gives in these units (`predicted_boost`), plus the no-rule baseline of 7d and the
multiplicative form of 7f/7h. B is the promoter's measured baseline and b_i the measured boost index
of single i, so S_i / B = 2 ** b_i:

  no_rule             P / B = 1                       (the promoter alone, 7d's no-rule baseline)
  mean                P / B = (2**b1 + 2**b2) / 2     (the runtime's retained rule)
  max                 P / B = 2**max(b1, b2)
  additive_envelope   max(b1, b2) <= boost <= log2(2**b1 + 2**b2 - 1)   (sum_capped and or)
  multiplicative      P / B = 2**(b1 + b2)            (Zhou et al. 2024, the log link of 7f)
  competition         P / B = 1 + 2*(2**bmax - 1) * W/(1 + W)           (this lane's form)

`additive_envelope` is an INTERVAL because sum_capped's cap and or's product are in the engine's
[0, 1] drive units, for which this assay gives no measured mapping; the interval is the whole of
what those two rules say without a unit assumption, and it is the envelope lane-combine itself used.
A rule's error inside its interval is zero.
"""

from __future__ import annotations

import csv
import gzip
import math
import random
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from genomeos.runtime import grn

LANE = "lane-strengthform"
ROADMAP_ITEM = (
    "area I Next 1 (ROADMAP §5 item 11, review R3): the strength-dependent combination form R3 left"
    " open, registered before it is run"
)

#: the machine-local store this lane reads, and the one file in it (scripts/strengthform_fetch.py)
STORE = Path("data/cache/strengthform")
TABLE = STORE / "GSE240586_EEP_normalized_data_MMA2023030806.tsv.gz"
TABLE_SHA256 = "8a0160c0022546187831710230d7a521ca32f4d7db602ceb9ddad403665d191c"
TABLE_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE240nnn/GSE240586/suppl/"
    "GSE240586_EEP_normalized_data_MMA2023030806.tsv.gz"
)
SOURCE = {
    "accession": "GEO GSE240586",
    "publication": (
        "Martinez-Ara M, Comoglio F, van Steensel B (2024) Large-scale analysis of the integration of"
        " enhancer-enhancer signals by promoters. eLife 12:e91994, doi 10.7554/eLife.91994"
    ),
    "assay": (
        "three-way combinatorial reporter assay (enhancer position 1, enhancer position 2, promoter),"
        " barcoded, in mouse E14tg2a embryonic stem cells"
    ),
    "file": str(TABLE),
    "url": TABLE_URL,
    "sha256": TABLE_SHA256,
    "prior_use_in_this_project": "none: neither the publication nor the accession appears in docs/",
}

#: the form under test, as `runtime/grn.py` already names it (an opt-in rule; the engine is unchanged)
FORM = "competition"
#: the one declared constant: the promoter's strongest single engages it this share of the time. The
#: same number 7d declared (`attribution/combine2.THETA`), not fitted here or there.
THETA = 0.5
#: K in A = W / (K + W). With THETA = 0.5 the strongest single's term is 1, so K = 1 fixes the units.
HALF_SATURATION = 1.0
#: every comparator, in the order they are reported
COMPARATORS = ("no_rule", "mean", "max", "additive_envelope", "multiplicative")
INTERVAL_RULES = ("additive_envelope",)

#: the columns read before the run: the design, the promoter baseline, and the two single boosts,
#: which are the authors' own position- and strand-matched definition of a single enhancer's effect
DESIGN_COLUMNS = (
    "id",
    "class",
    "Promoter",
    "frag1",
    "strand1",
    "frag2",
    "strand2",
    "Baseline",
    "Boost1",
    "Boost2",
)
#: the one column the run opens: the measured activity of the combination
OUTCOME_COLUMN = "activity_all"
#: columns this lane never reads, because they are the authors' own expectations and comparisons and
#: reading them would make this test a comparison with their analysis rather than with their data
COLUMNS_NEVER_READ = (
    "AdditiveCoop",
    "Additivebaseline",
    "MiltipCoop",
    "Coop",
    "Coop2",
    "R1",
    "R2",
    "boost",
    "activity_br1",
    "activity_br2",
    "activity_br3",
)

#: the observed quantity: log2 of the combination's activity over its promoter's baseline
ENDPOINT = (
    "the squared error of the observed boost index log2(activity_all / Baseline) of an"
    " enhancer-enhancer-promoter combination against each form's predicted boost index, in log2 units"
)
#: a combination is informative only under all of these, each counted in the result
INFORMATIVE = (
    "class EE (both fragments are enhancers, by the table's own class column)",
    "the promoter's Baseline > 0",
    "both single boost indices present (not NA) and finite",
    "both single boost indices > 0: both singles activate above the promoter's baseline",
    "the promoter's largest single boost index > 0",
    "the combination's activity_all > 0",
)
#: resampling: the unordered enhancer pair, ignoring promoter, position and orientation, because the
#: same pair of enhancers recurs across all eight promoters and both orientations of each position
CLUSTER = "the unordered pair of enhancer fragment names"
BOOTSTRAP_DRAWS = 2000
BOOTSTRAP_SEED = 20261010
CI_LEVEL = 0.95
#: a comparator is read as better or worse only when the interval excludes zero
READING_WORDS = ("better", "worse", "no difference detected")
FALSIFIER = (
    "The claim that the strength-dependent form is better than the rules lane-combine compared is"
    " FALSIFIED unless the form's mean squared error is below EVERY comparator's (no_rule, mean, max,"
    " additive_envelope, multiplicative) with the 95% cluster-bootstrap interval of the difference"
    " excluding zero. A single comparator it does not beat is enough to falsify it, and a comparator"
    " that beats it with its interval excluding zero is reported as 'worse'."
)
#: reported, deciding nothing: what lane-combine's own parameter-free envelope says about these pairs
SECONDARY = (
    "the share of informative combinations above the additive expectation (supra-additive, which no"
    " monotone rule of the terms can produce) and below the stronger single (which contradicts every"
    " monotone rule, the Bothma snail case lane-combine used), reported as counts and deciding nothing"
)

DECISIONS = (
    "the population is Martinez-Ara et al. 2024 (GSE240586), not lane-combine's Bothma pairs, because"
    " lane-combine read those outcomes before writing its rule and said so",
    "the form is the opt-in `competition` rule already in runtime/grn.py; the engine is not changed,"
    " and a test pins its header and every existing form byte for byte",
    f"the one free number per promoter is fixed by declaration, THETA = {THETA}, the constant 7d"
    " declared (attribution/combine2.THETA); nothing is fitted on any enhancer pair",
    "the singles are the authors' own Boost1 and Boost2 columns, position- and strand-matched, so no"
    " reconstruction of a single's effect enters the comparison",
    "the no-rule baseline, the mean, max, the sum_capped/or envelope and the multiplicative form are"
    " all scored on the same combinations, so no form is read on a population another cannot reach",
)


def _f(value: str) -> float | None:
    """A table cell as a float, or None when it is NA or unparseable."""
    try:
        x = float(value)
    except (TypeError, ValueError):
        return None
    return x if math.isfinite(x) else None


@dataclass(frozen=True)
class Combination:
    """One enhancer-enhancer-promoter row: its design, its two single boosts, its outcome."""

    promoter: str
    frag1: str
    strand1: str
    frag2: str
    strand2: str
    baseline: float
    b1: float
    b2: float
    activity: float | None = None

    @property
    def cluster(self) -> tuple[str, str]:
        return tuple(sorted((self.frag1, self.frag2)))  # type: ignore[return-value]

    @property
    def observed_boost(self) -> float:
        if self.activity is None or self.activity <= 0 or self.baseline <= 0:
            raise ValueError("no observed boost: the activity was not read, or is not positive")
        return math.log2(self.activity / self.baseline)


def read_table(path: Path = TABLE, *, read_outcome: bool) -> dict[str, Any]:
    """The design, the baselines, the singles and (only if `read_outcome`) the pair activities.

    `read_outcome` False never touches `OUTCOME_COLUMN` of an EE row, so the census and the
    registration can be computed and committed with no pair outcome read.
    """
    baselines: dict[str, float] = {}
    rows: list[Combination] = []
    counts = {
        "rows": 0,
        "class_EE": 0,
        "class_not_EE": 0,
        "no_baseline": 0,
        "single_missing": 0,
        "single_not_above_baseline": 0,
        "no_promoter_maximum": 0,
        "activity_not_positive": 0,
        "informative": 0,
    }
    singles: dict[tuple[str, int, str, str], float] = {}
    with gzip.open(path, "rt", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            counts["rows"] += 1
            promoter = row["Promoter"]
            baseline = _f(row["Baseline"])
            if baseline is not None:
                baselines[promoter] = baseline
            if row["class"] != "EE":
                counts["class_not_EE"] += 1
                continue
            counts["class_EE"] += 1
            b1, b2 = _f(row["Boost1"]), _f(row["Boost2"])
            if b1 is not None:
                singles[(promoter, 1, row["frag1"], row["strand1"])] = b1
            if b2 is not None:
                singles[(promoter, 2, row["frag2"], row["strand2"])] = b2
            if baseline is None or baseline <= 0:
                counts["no_baseline"] += 1
                continue
            if b1 is None or b2 is None:
                counts["single_missing"] += 1
                continue
            if not (b1 > 0 and b2 > 0):
                counts["single_not_above_baseline"] += 1
                continue
            activity = _f(row[OUTCOME_COLUMN]) if read_outcome else None
            if read_outcome and (activity is None or activity <= 0):
                counts["activity_not_positive"] += 1
                continue
            rows.append(
                Combination(
                    promoter=promoter,
                    frag1=row["frag1"],
                    strand1=row["strand1"],
                    frag2=row["frag2"],
                    strand2=row["strand2"],
                    baseline=baseline,
                    b1=b1,
                    b2=b2,
                    activity=activity,
                )
            )
    maxima = promoter_maxima(singles)
    kept = [c for c in rows if maxima.get(c.promoter, 0.0) > 0]
    counts["no_promoter_maximum"] = len(rows) - len(kept)
    counts["informative"] = len(kept)
    return {
        "baselines": baselines,
        "singles": singles,
        "promoter_maxima": maxima,
        "combinations": kept,
        "counts": counts,
    }


def promoter_maxima(singles: dict[tuple[str, int, str, str], float]) -> dict[str, float]:
    """Each promoter's largest single boost index, over both positions and both orientations.

    This is what THETA is declared about: that single, and no other, half-engages its promoter.
    """
    out: dict[str, float] = {}
    for (promoter, _position, _frag, _strand), boost in singles.items():
        out[promoter] = max(out.get(promoter, float("-inf")), boost)
    return out


def terms(b1: float, b2: float, bmax: float) -> tuple[float, float]:
    """The two activator terms of the strength-dependent form, from the measured singles.

    r / B = 2 * (2**bmax - 1) by the declaration (the strongest single's term is 1 at THETA = 0.5),
    y_i = (S_i - B) / r is that single's own occupancy, and t_i = y_i / (1 - y_i) its term. A single
    at the promoter's maximum gives t = 1 exactly, and the form reproduces every single it is given.
    """
    r_rel = 2.0 * (2.0**bmax - 1.0)
    if not r_rel > 0:
        raise ValueError(f"the promoter's largest single boost must be > 0, got {bmax!r}")
    out = []
    for b in (b1, b2):
        y = (2.0**b - 1.0) / r_rel
        if not 0 < y < 1:
            raise ValueError(f"a single's occupancy must lie in (0, 1), got {y!r} from boost {b!r}")
        out.append(y / (1.0 - y))
    return out[0], out[1]


def predicted_boost(form: str, b1: float, b2: float, bmax: float) -> float | tuple[float, float]:
    """A form's prediction of a pair's boost index: one number, or an interval for an interval rule.

    Every form is computed through the engine's own `grn.combine_activators`, on terms in the units
    the engine takes, so no form is reimplemented here.
    """
    linear = [2.0**b1 - 1.0, 2.0**b2 - 1.0]  # t_i in units of B, for the fixed rules
    if form == "no_rule":
        return 0.0
    if form == "mean":
        return math.log2(1.0 + grn.combine_activators(linear, "mean"))
    if form == "max":
        return math.log2(1.0 + grn.combine_activators(linear, "max"))
    if form == "additive_envelope":
        lower = math.log2(1.0 + grn.combine_activators(linear, "max"))
        upper = math.log2(1.0 + sum(linear))
        return (lower, upper)
    if form == "multiplicative":
        gained = grn.combine_activators([b1 * math.log(2), b2 * math.log(2)], "multiplicative", gain=1.0)
        return math.log2(gained)
    if form == FORM:
        t1, t2 = terms(b1, b2, bmax)
        drive = grn.combine_activators([t1, t2], FORM, HALF_SATURATION)
        return math.log2(1.0 + 2.0 * (2.0**bmax - 1.0) * drive)
    raise ValueError(f"unknown form {form!r}")


def squared_error(observed: float, prediction: float | tuple[float, float]) -> float:
    """The squared error in log2 units; zero inside an interval rule's interval."""
    if isinstance(prediction, tuple):
        lower, upper = prediction
        if lower <= observed <= upper:
            return 0.0
        edge = lower if observed < lower else upper
        return (observed - edge) ** 2
    return (observed - prediction) ** 2


def mean_squared_error(errors: list[float]) -> float:
    return sum(errors) / len(errors) if errors else float("nan")


def errors_per_form(combinations: list[Combination], maxima: dict[str, float]) -> dict[str, list[float]]:
    """Every form's per-combination squared error, on exactly the same combinations."""
    out: dict[str, list[float]] = {name: [] for name in (FORM, *COMPARATORS)}
    for c in combinations:
        observed = c.observed_boost
        bmax = maxima[c.promoter]
        for name in out:
            out[name].append(squared_error(observed, predicted_boost(name, c.b1, c.b2, bmax)))
    return out


def cluster_bootstrap(
    clusters: list[tuple[str, str]],
    errors: dict[str, list[float]],
    *,
    draws: int = BOOTSTRAP_DRAWS,
    seed: int = BOOTSTRAP_SEED,
    level: float = CI_LEVEL,
) -> dict[str, dict[str, float]]:
    """Per comparator, the interval of MSE(comparator) - MSE(the form), resampling whole clusters.

    Positive means the strength-dependent form is closer to the measured pairs than that comparator.
    """
    index: dict[tuple[str, str], list[int]] = {}
    for i, key in enumerate(clusters):
        index.setdefault(key, []).append(i)
    keys = sorted(index)
    rng = random.Random(seed)
    deltas: dict[str, list[float]] = {name: [] for name in COMPARATORS}
    for _ in range(draws):
        picked: list[int] = []
        for _ in keys:
            picked.extend(index[keys[rng.randrange(len(keys))]])
        own = sum(errors[FORM][i] for i in picked) / len(picked)
        for name in COMPARATORS:
            deltas[name].append(sum(errors[name][i] for i in picked) / len(picked) - own)
    lo_q, hi_q = (1 - level) / 2, 1 - (1 - level) / 2
    out: dict[str, dict[str, float]] = {}
    for name, values in deltas.items():
        values.sort()
        lo = values[min(len(values) - 1, int(lo_q * len(values)))]
        hi = values[min(len(values) - 1, int(hi_q * len(values)))]
        point = mean_squared_error(errors[name]) - mean_squared_error(errors[FORM])
        out[name] = {"delta_mse": point, "ci_low": lo, "ci_high": hi, "draws": len(values)}
    return out


def word(interval: dict[str, float]) -> str:
    """The registered reading word for one comparator, from its interval alone."""
    if interval["ci_low"] > 0 and interval["ci_high"] > 0:
        return "better"
    if interval["ci_low"] < 0 and interval["ci_high"] < 0:
        return "worse"
    return "no difference detected"


#: the reading's wording, fixed before the run and rendered by `reading` whatever the outcome is
READING_TEMPLATE = (
    "the strength-dependent combination form (promoter competition, Bothma et al. 2015 eLife 4:e07956"
    " Methods eqs. 7-9, THETA = 0.5 declared and nothing fitted on any pair), run once on the held-out"
    " enhancer-enhancer-promoter combinations of Martinez-Ara et al. 2024 (GSE240586, mouse E14tg2a"
    " embryonic stem cells): {verdict}; over {n} combinations ({clusters} unordered enhancer pairs,"
    " {promoters} promoters) its mean squared error of the measured boost index is {mse} log2 squared,"
    " and against each comparator: {per_comparator}. {falsified}"
)
VERDICT_BETTER = "better than every rule lane-combine compared"
VERDICT_NOT_BETTER = "not better than every rule lane-combine compared"
FALSIFIED_YES = (
    "The claim that the strength-dependent form is better than those rules is FALSIFIED on this population."
)
FALSIFIED_NO = (
    "The claim that the strength-dependent form is better than those rules is NOT FALSIFIED on this"
    " population."
)


def reading(
    intervals: dict[str, dict[str, float]],
    *,
    own_mse: float,
    n: int,
    clusters: int,
    promoters: int,
) -> dict[str, Any]:
    """The registered reading, rendered from the intervals alone, whatever they are."""
    words = {name: word(intervals[name]) for name in COMPARATORS}
    better = all(words[name] == "better" for name in COMPARATORS)
    parts = []
    for name in COMPARATORS:
        row = intervals[name]
        parts.append(
            f"{name}: {words[name]} (Delta MSE {row['delta_mse']:+.4f}"
            f" [95% CI {row['ci_low']:+.4f}, {row['ci_high']:+.4f}])"
        )
    not_beaten = [name for name in COMPARATORS if words[name] != "better"]
    verdict = VERDICT_BETTER if better else f"{VERDICT_NOT_BETTER} ({', '.join(not_beaten)})"
    sentence = READING_TEMPLATE.format(
        verdict=verdict,
        n=n,
        clusters=clusters,
        promoters=promoters,
        mse=f"{own_mse:.4f}",
        per_comparator="; ".join(parts),
        falsified=FALSIFIED_NO if better else FALSIFIED_YES,
    )
    return {
        "sentence": sentence,
        "words": words,
        "better_than_every_comparator": better,
        "falsified": not better,
    }


def registered_fields() -> dict[str, Any]:
    """Every binding word of this lane, as the registration writes it and the run re-checks it."""
    return {
        "form": FORM,
        "form_source": (
            "Bothma et al. 2015, eLife 4:e07956, Methods eqs. 7-9, as runtime/grn.py's opt-in"
            " `competition` rule: A = W / (K + W), W the sum of the activator terms"
        ),
        "theta": THETA,
        "half_saturation": HALF_SATURATION,
        "population": SOURCE,
        "endpoint": ENDPOINT,
        "informative": list(INFORMATIVE),
        "comparators": list(COMPARATORS),
        "interval_rules": list(INTERVAL_RULES),
        "statistic": (
            "MSE(comparator) - MSE(the form) in log2 units, with a cluster bootstrap over"
            f" {CLUSTER} ({BOOTSTRAP_DRAWS} draws, seed {BOOTSTRAP_SEED}) and a"
            f" {int(CI_LEVEL * 100)}% percentile interval"
        ),
        "cluster": CLUSTER,
        "bootstrap_draws": BOOTSTRAP_DRAWS,
        "bootstrap_seed": BOOTSTRAP_SEED,
        "ci_level": CI_LEVEL,
        "reading_words": list(READING_WORDS),
        "reading_template": READING_TEMPLATE,
        "verdict_words": [VERDICT_BETTER, VERDICT_NOT_BETTER],
        "falsifier": FALSIFIER,
        "secondary_deciding_nothing": SECONDARY,
        "decisions": list(DECISIONS),
        "columns_read_before_the_run": list(DESIGN_COLUMNS),
        "column_the_run_opens": OUTCOME_COLUMN,
        "columns_never_read": list(COLUMNS_NEVER_READ),
    }
