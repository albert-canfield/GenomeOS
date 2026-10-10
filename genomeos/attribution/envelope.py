# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""lane-envelope (ROADMAP area I, 2026-10-10): the additive envelope as its own combination form.

WHAT WAS ASKED. lane-strengthform (`7d23bc3` code, `10ca354` registration, `fc328c5` result) ran the
strength-dependent form (promoter competition) once on held-out enhancer pairs and found exactly one
comparator it did not beat: `additive_envelope`, worse by Delta MSE -0.0463 [-0.0488, -0.0440]. The
envelope has never been held against measured RNA as its own form, so this lane was asked to register
it on 7d's population (combine2.py) with 7h's wording, after establishing from the FITS alone --
before any RNA value is opened -- whether the design can distinguish it from the forms already read.

THE DEFINITION, READ FIRST-HAND FROM THE LANE THAT USED IT. `strengthform.predicted_boost`
("additive_envelope") returns an INTERVAL, not a number:

    lower = log2(1 + grn.combine_activators(linear, "max"))        # the engine's `max` rule
    upper = log2(1 + sum(linear))                                   # the uncapped additive sum
    linear = [2**b1 - 1, 2**b2 - 1]                                 # t_i in units of the baseline B

so in boost units max(b1, b2) <= boost <= log2(2**b1 + 2**b2 - 1), and
`strengthform.squared_error` gives an interval rule ZERO error anywhere between the two ends.
The upper end is additivity above the baseline, P = (S1 - B) + (S2 - B) + B: in the engine's [0, 1]
drive units that is the uncapped sum of the terms, and with the cap `sum_capped` carries it is
A = min(1, W), W = sum of s_i*H(x_i). That IS `grn.combine_activators`'s `sum_capped`, the rule
`grn.ACTIVATOR_RULES` has held since 2026-09-28. The lower end is `grn`'s `max`, one of the four
rules lane-combine compared on 2026-09-28. No free number is added here, and the engine is not
touched; `identity_checks` verifies both identities numerically rather than asserting them
(`IDENTITY`).

CONSEQUENCE FOR THE RNA QUESTION (`RNA_SIDE`). The envelope's upper end as a form, with its cap, is
not merely close to the form lane-capsumrna ran against RNA -- it is the SAME engine rule,
`sum_capped`, at the same conventions. So lane-capsumrna's census (`43dd4ac`, `8994fca`, `f74dca5`;
data/results/capsumrna_census.json) IS this form's census on 7d's population, and this lane reads its
committed counts rather than re-deriving them: of 17,983 multi-element genes combine2.classify admits
the form fits 2,288 on its uncapped root, on every one of which its prediction in the thirteen reader
cells is the runtime's mean's (largest difference 6.7e-16 in log2, because the mean is linear under
the reader's hold and so is a sum below its cap), 2,165 genes have a unique fit and 3 of those a
prediction other than the mean's. 7d's DECISION 2 (`combine2.MIN_GENES` = 100) makes 3 genes
NOT_TESTABLE, so no registration is written and NO RNA VALUE IS OPENED. The 126 genes with a second,
capped root are the ambiguity lane-rootpick is registered to decide; this lane picks no root and
reports the count as that lane's.

WHAT THIS LANE ADDS BESIDE THE IDENTITY (`PAIR_SIDE`, descriptive and deciding nothing). Because the
envelope is an interval, "the rule that beat the strength-dependent form" is not one combination form,
and the win is not evidence for either end of it. On lane-strengthform's own already-read,
already-registered combinations, run through its own frozen code, this lane reports each end as a
POINT form beside the interval, and the share of combinations that fall inside the interval and so
score zero. These are descriptive statistics of a comparator's definition on a population whose
outcomes that lane read and committed; they register nothing, they decide nothing here, and no
registered reading of any lane is restated or revised by them.

No model request, no key, no money. No new biosample: the thirteen reader cells unchanged.
"""

from __future__ import annotations

import json
import math
import random
from pathlib import Path
from typing import Any

from genomeos.attribution import combine2 as c2
from genomeos.attribution import strengthform as sf
from genomeos.runtime import grn

DATE = "2026-10-10"
LANE = "lane-envelope"
ROADMAP_ITEM = (
    "area I: the additive envelope (lane-strengthform's one unbeaten comparator) as its own combination"
    " form on 7d's population, established from the fits before any RNA value is opened"
)

#: the engine rule the envelope's UPPER end is, with the cap it names: A = min(1, W)
FORM = "sum_capped"
#: the engine rule the envelope's LOWER end is: A = max(t_i)
LOWER_FORM = "max"
#: the form lane-strengthform tested, against which the pair-side point forms are differenced
REFERENCE_FORM = sf.FORM

MIN_GENES = c2.MIN_GENES
CELLS = c2.CELLS

#: lane-capsumrna's committed census: the same engine rule, the same conventions, on 7d's population
CAPSUM_CENSUS = Path("data/results/capsumrna_census.json")
CAPSUM_COMMITS = ("43dd4ac", "8994fca", "f74dca5")
STRENGTHFORM_RESULT = Path("data/results/strengthform_result.json")
STRENGTHFORM_COMMITS = ("7d23bc3", "10ca354", "fc328c5")

DEFINITION = (
    "lane-strengthform's `additive_envelope` (strengthform.COMPARATORS, strengthform.INTERVAL_RULES) is"
    " the INTERVAL max(b1, b2) <= boost <= log2(2**b1 + 2**b2 - 1), computed in that lane's"
    " predicted_boost as lower = log2(1 + grn.combine_activators(linear, 'max')) and"
    " upper = log2(1 + sum(linear)) on linear = [2**b1 - 1, 2**b2 - 1], with zero error anywhere"
    " between the ends (strengthform.squared_error). Its UPPER end is additivity above the baseline,"
    " P = (S1 - B) + (S2 - B) + B, which in the engine's [0, 1] drive units is the uncapped sum of the"
    " terms; with the cap sum_capped carries it is A = min(1, W), W = sum of s_i*H(x_i), which IS"
    " grn.combine_activators's 'sum_capped' (grn.ACTIVATOR_RULES). Its LOWER end is grn's 'max'. The"
    " form registered on RNA here is the upper end with its cap, as the brief fixes it; no free number"
    " is added and runtime/grn.py is not changed"
)

IDENTITY = (
    "two identities, verified numerically by identity_checks rather than asserted: (1) in the engine's"
    " drive units, min(1, sum(t)) equals grn.combine_activators(t, 'sum_capped') on random term"
    " vectors; (2) in the assay's boost units, log2(1 + sum(2**b_i - 1)) equals the upper end"
    " strengthform.predicted_boost('additive_envelope', ...) returns, and log2(1 + max(2**b_i - 1))"
    " equals its lower end, which is also what grn's 'max' gives. So the envelope's two ends are two"
    " rules the engine already carries, and the upper end is the form lane-capsumrna ran against RNA"
)

RNA_SIDE = (
    "because the envelope's upper end with its cap is the SAME engine rule at the SAME conventions as"
    " lane-capsumrna's form (grn 'sum_capped' at combine2's CONVENTIONS), that lane's committed census"
    " is this form's census on 7d's population and is read here, not re-derived: the counts come from"
    " data/results/capsumrna_census.json (commits 43dd4ac, 8994fca, f74dca5). The structural reason is"
    " combine2.STRUCTURE's: under the reader's hold the mean is linear in the element terms, and so is"
    " a sum below its cap, so on its uncapped root the capped sum gives the mean's level in every"
    " reader cell. No RNA value is opened by this lane"
)

PAIR_SIDE = (
    "descriptive, registering nothing and deciding nothing: on lane-strengthform's own combinations,"
    " read and committed by that lane (data/results/strengthform_result.json, fc328c5) and recomputed"
    " through its own frozen code, each end of the envelope as a POINT form beside the interval, with"
    " the share of combinations inside the interval (which score zero under an interval rule) and the"
    " same cluster bootstrap that lane registered (seed and draws from strengthform). paired_bootstrap"
    " is checked against that lane's committed intervals for its five registered comparators before"
    " any point form is reported, so the procedure is the registered one. No registered reading is"
    " restated or revised"
)

DECISION = (
    "fixed before any count is read, and unchanged by any of them. (1) The form registered on RNA is"
    " the envelope's upper end with its cap, which identity_checks must show to be grn's 'sum_capped';"
    " if the identity fails, nothing is registered and the failure is the result. (2) Given the"
    " identity, the form's census on 7d's population is lane-capsumrna's committed census, and the"
    " genes on which the form could be told apart from what 7d has already read are the genes whose"
    " fit is unique and whose prediction is not the mean's. (3) 7d's DECISION 2 applies unchanged:"
    " fewer than combine2.MIN_GENES = 100 such genes is NOT_TESTABLE -- no registration is written, no"
    " comparison is read, and no RNA value is opened. (4) The genes with a second, capped root are"
    " lane-rootpick's registered decision; this lane picks no root and reports their count as that"
    " lane's question"
)
DECISIONS = ("not_testable_below_7d_minimum", "register_on_the_genes_that_differ", "identity_failed")

CENSUS_SENTENCE = (
    "the additive envelope, lane-strengthform's one unbeaten comparator, is an INTERVAL between two"
    " rules runtime/grn.py already carries: its lower end is 'max' and its upper end, with the cap it"
    " names, is 'sum_capped' (both identities exact to {identity_worst} in the checks run here). As its"
    " own combination form on 7d's population the upper end is therefore the form lane-capsumrna ran:"
    " of the {multi_element} multi-element genes combine2.classify admits it fits {fitted} on its"
    " uncapped root, on every one of which its prediction in the {cells} reader cells is the runtime's"
    " mean's (largest difference {largest} in log2), {unique} genes have a unique fit and {differ} of"
    " those a prediction other than the mean's. {differ} is below 7d's DECISION 2 minimum of"
    " {min_genes} genes, so the envelope as its own form is NOT_TESTABLE on this design, nothing is"
    " registered and no RNA value was opened. The {two_roots} genes with a second, capped root are"
    " lane-rootpick's registered decision and no root is picked here"
)

PAIR_SENTENCE = (
    "descriptive and deciding nothing, on lane-strengthform's committed combinations: of {n}"
    " combinations {inside} ({inside_pct}%) fall inside the envelope and score zero under an interval"
    " rule. As POINT forms the ends are not the interval: the upper end's MSE is {upper_mse} against"
    " the interval's {interval_mse} and the strength-dependent form's {reference_mse}, a difference"
    " from that form of {upper_delta} [95% CI {upper_lo}, {upper_hi}] ({upper_word}), while the"
    " interval's difference is {interval_delta} [95% CI {interval_lo}, {interval_hi}] ({interval_word})"
    " and the lower end's is {lower_delta} [95% CI {lower_lo}, {lower_hi}] ({lower_word}). So the"
    " margin lane-strengthform's registered reading reports against additive_envelope is the"
    " interval's, not either end's as a form; that reading stands as registered and is not revised"
)

#: the term vectors identity_checks draws, fixed here so the check is the same run to run
IDENTITY_DRAWS = 20_000
IDENTITY_SEED = 20261010
IDENTITY_TOLERANCE = 1e-12


def identity_checks(*, draws: int = IDENTITY_DRAWS, seed: int = IDENTITY_SEED) -> dict[str, Any]:
    """The envelope's two ends against the two engine rules, numerically. No data is read."""
    rng = random.Random(seed)
    worst_drive = 0.0
    worst_upper = 0.0
    worst_lower = 0.0
    worst_lower_rule = 0.0
    for _ in range(draws):
        n = rng.randint(2, 6)
        terms = [rng.random() for _ in range(n)]
        closed = min(1.0, sum(terms))
        worst_drive = max(worst_drive, abs(closed - grn.combine_activators(terms, FORM)))
        b1, b2 = rng.random() * 3.0, rng.random() * 3.0
        linear = [2.0**b1 - 1.0, 2.0**b2 - 1.0]
        lower, upper = sf.predicted_boost("additive_envelope", b1, b2, max(b1, b2))
        worst_upper = max(worst_upper, abs(math.log2(1.0 + sum(linear)) - upper))
        worst_lower = max(worst_lower, abs(math.log2(1.0 + max(linear)) - lower))
        by_rule = math.log2(1.0 + grn.combine_activators(linear, LOWER_FORM))
        worst_lower_rule = max(worst_lower_rule, abs(by_rule - lower))
    worst = max(worst_drive, worst_upper, worst_lower, worst_lower_rule)
    return {
        "draws": draws,
        "seed": seed,
        "upper_end_is_sum_capped_in_drive_units": worst_drive,
        "upper_end_is_the_uncapped_additive_sum_in_boost_units": worst_upper,
        "lower_end_is_the_stronger_single_in_boost_units": worst_lower,
        "lower_end_is_grn_max": worst_lower_rule,
        "worst": worst,
        "tolerance": IDENTITY_TOLERANCE,
        "holds": worst <= IDENTITY_TOLERANCE,
    }


def paired_bootstrap(
    clusters: list[tuple[str, str]],
    errors: dict[str, list[float]],
    *,
    reference: str = REFERENCE_FORM,
    draws: int = sf.BOOTSTRAP_DRAWS,
    seed: int = sf.BOOTSTRAP_SEED,
    level: float = sf.CI_LEVEL,
) -> dict[str, dict[str, float]]:
    """MSE(form) - MSE(reference) per form, resampling whole clusters.

    lane-strengthform's own `cluster_bootstrap` with the reference form a parameter, so the five
    comparators it registered come back identical under the same seed (`check_against_committed`).
    """
    index: dict[tuple[str, str], list[int]] = {}
    for i, key in enumerate(clusters):
        index.setdefault(key, []).append(i)
    keys = sorted(index)
    rng = random.Random(seed)
    names = [name for name in errors if name != reference]
    deltas: dict[str, list[float]] = {name: [] for name in names}
    for _ in range(draws):
        picked: list[int] = []
        for _ in keys:
            picked.extend(index[keys[rng.randrange(len(keys))]])
        own = sum(errors[reference][i] for i in picked) / len(picked)
        for name in names:
            deltas[name].append(sum(errors[name][i] for i in picked) / len(picked) - own)
    lo_q, hi_q = (1 - level) / 2, 1 - (1 - level) / 2
    out: dict[str, dict[str, float]] = {}
    for name, values in deltas.items():
        values.sort()
        lo = values[min(len(values) - 1, int(lo_q * len(values)))]
        hi = values[min(len(values) - 1, int(hi_q * len(values)))]
        point = sf.mean_squared_error(errors[name]) - sf.mean_squared_error(errors[reference])
        out[name] = {"delta_mse": point, "ci_low": lo, "ci_high": hi, "draws": len(values)}
    return out


def point_errors(combinations: list[sf.Combination], maxima: dict[str, float]) -> dict[str, list[float]]:
    """Each end of the envelope as a point form, and which combinations fall inside the interval.

    The two ends come from `strengthform.predicted_boost` itself, so neither is reimplemented.
    """
    upper: list[float] = []
    lower: list[float] = []
    inside: list[float] = []
    for c in combinations:
        lo, hi = sf.predicted_boost("additive_envelope", c.b1, c.b2, maxima[c.promoter])
        upper.append((c.observed_boost - hi) ** 2)
        lower.append((c.observed_boost - lo) ** 2)
        inside.append(1.0 if lo <= c.observed_boost <= hi else 0.0)
    return {"envelope_upper_point": upper, "envelope_lower_point": lower, "inside": inside}


def check_against_committed(
    intervals: dict[str, dict[str, float]], result: Path = STRENGTHFORM_RESULT
) -> dict[str, Any]:
    """`paired_bootstrap` must reproduce lane-strengthform's committed intervals, to 1e-12."""
    committed = json.loads(result.read_text())["intervals"]
    problems = []
    worst = 0.0
    for name, row in committed.items():
        if name not in intervals:
            problems.append(f"{name}: not recomputed")
            continue
        for key in ("delta_mse", "ci_low", "ci_high"):
            diff = abs(intervals[name][key] - row[key])
            worst = max(worst, diff)
            if diff > 1e-12:
                problems.append(f"{name}.{key}: {intervals[name][key]!r} vs committed {row[key]!r}")
    return {
        "source": str(result),
        "comparators": sorted(committed),
        "worst": worst,
        "problems": problems,
    }


def capsum_counts(path: Path = CAPSUM_CENSUS) -> dict[str, Any]:
    """The committed census of the same engine rule on 7d's population, read field by field."""
    census = json.loads(path.read_text())
    sel = census["reader_selection"]
    return {
        "source": str(path),
        "commits": list(CAPSUM_COMMITS),
        "form": census["form"],
        "eligible_for_a_fit": sel["eligible_for_a_fit"],
        "capped_sum_status": sel["capped_sum_status"],
        "unique_fit": sel["unique_fit"],
        "two_roots": sel["two_roots"],
        "size_check": census["size_check"],
        "sentence": census["sentence"],
        "decision": census["decision"],
        "rna_opened": census["size_check"]["rna_opened"],
    }


def decide(identity: dict[str, Any], differ: int) -> str:
    """DECISION, from the identity and the count of genes that could tell the form apart."""
    if not identity["holds"]:
        return "identity_failed"
    if differ < MIN_GENES:
        return "not_testable_below_7d_minimum"
    return "register_on_the_genes_that_differ"
