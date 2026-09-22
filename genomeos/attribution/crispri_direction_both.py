# SPDX-License-Identifier: AGPL-3.0-or-later
"""The upward arm: the half of `action` that the passing direction result could not test.

`crispri_direction.py` asked whether the deletion layer gets the **sign** of an enhancer-gene effect
right, and it passed: 41 of 44 held-out pairs, 0.9318, Wilson 95% [0.8177, 0.9765], against a
constant-sign caller worth 0.5361 and a post-hoc magnitude-matched 0.6537. That result named its own
limit in its last paragraph. **All 44 answerable pairs carried a measured decrease.** Zero upward. So
nothing in it tested the `represses` half of `action`, and the only three upward calls the model made
were all wrong, 0 of 3. A claim that has only been tested where it succeeds is not yet a claim.

This module removes that limit, and it does so for **zero requests**. The costed extension was 30
AlphaGenome deletion requests for K562 and 35 with GM12878 — the distinct elements behind the 33 and
38 in-reach upward held-out pairs. That costing is reproduced exactly here (`costing`). It is also
unnecessary, and the reason is worth stating precisely, because it was a fact about a *derived table*
being mistaken for a fact about the measurement.

**The sweep did not discard the other genes; one of its two outputs did.** The compact table
`data/knowledge/alphagenome/all_elements/<chrom>.json` keeps one `predicted` and one
`predicted_coding` gene per element, and that is the table `crispri_direction.py` joins to. But the
run that built it wrote a second output, the per-element response cache
`data/knowledge/alphagenome/elements/<chrom>.json.gz`, and `scripts/enhancer_targets_all.py` scored
it with `threshold=0.0` (`worker_scorer`), so that cache carries **every gene in the scorer's 1 Mb
window with its signed log2 fold change on each of K562, HepG2, GM12878 and IMR-90's own track**,
uncensored. The quantity 30 requests would have purchased — the predicted log2 fold change for a
*named* gene on deleting a *named* element in a *named* cell line — is already on disk for every
element the sweep ever touched. Reading it needs a new reader over `genes`, not a new request.

That turns "the sweep stores one target per element, so the upward arm is unanswerable" from a
limitation of the measurement into a limitation of one projection of it, and it grows the answerable
held-out set from 44 one-signed pairs to 152 two-signed ones.

**Why both signs change the test, and the trap that replaces the old one.** With one measured sign,
agreement was the model's own rate of emitting that sign, so 0.5 was the wrong chance level and the
registered baseline was the sweep's marginal down-rate. With both signs present that particular trap
is gone — but a second one takes its place immediately, and reporting raw agreement against 0.5 would
walk straight into it. The two-signed set is **not balanced**: it is roughly 116 down against 36 up,
so a caller that says "down" and nothing else scores about 0.76, not 0.5. Raw agreement against 0.5
would therefore manufacture a success exactly as before, one step further along.

So the registered primary statistic for the combined set is **balanced accuracy**, the unweighted
mean of the two per-sign sensitivities. Its chance level is 0.5 for any class mix, and a constant-sign
caller — of either sign — scores exactly 0.5 on it by construction. That is the honest form of "both
signs present, so 0.5 is available at last". Raw agreement is reported too, and never without the
majority-class rate printed beside it.

No AlphaGenome request is made anywhere in this module. `requests` is 0 and is asserted, not claimed.
"""

from __future__ import annotations

import gzip
import json
import random
from collections import defaultdict
from pathlib import Path
from typing import Any

from genomeos.attribution.crispri import HELDOUT, TRAINING
from genomeos.attribution.crispri_direction import (
    ELEMENTS,
    MAGNITUDE_EDGES,
    in_reach,
    magnitude_bin,
    marginal_down_rate,
    overlapping,
    rows,
    signed,
    wilson,
)

CACHE = Path("data/knowledge/alphagenome/elements")
CELLS = ("K562", "GM12878")
CHANCE = 0.5
CHANCE_MARGIN = 0.15  # the same margin the first half registered: a pass needs 0.65
MARGINAL_MARGIN = 0.10  # and the same 0.10 over the model's own habit
STRONG = 0.1  # the magnitude the first half found the model did not miss above (28 of 28)
WEAK = 0.05  # and below which every one of its three errors sat
BOOTSTRAPS = 2000
SEED = 20260922
# the 5c842ca headline, fixed here so the re-derivation check cannot be relaxed after the fact
PRIOR = {"k": 41, "n": 44, "rate": 0.9318}

EVIDENCE = (
    "experimental: signed EffectSize of CRISPRi enhancer-gene screens, ENCODE benchmark "
    "(EngreitzLab/CRISPR_comparison, Gschwind et al. 2025); predicted: signed AlphaGenome deletion "
    "log2 fold change for the named gene on the named cell line's own track, read from the "
    "per-element response cache of the finished genome-wide sweep (threshold=0.0, uncensored). "
    "No request"
)

PREREGISTERED = (
    "On the held-out arm, K562 and GM12878, reading the per-element response cache rather than the "
    "one-target-per-element compact table, so that upward measured pairs become answerable: "
    "(1) THE UPWARD ARM ALONE, measured increases only, sensitivity_up clears 0.5 by 0.15 and clears "
    "the model's marginal UP-rate over the whole sweep by 0.10 with the Wilson lower bound above it; "
    "(2) THE COMBINED TWO-SIGNED SET, primary statistic balanced accuracy = mean of the two per-sign "
    "sensitivities, whose chance level is 0.5 for any class mix and on which a constant-sign caller "
    "scores exactly 0.5: PASS at balanced accuracy >= 0.65 with the stratified bootstrap 95% lower "
    "bound above 0.5, FAIL at a lower bound at or below 0.5, UNDECIDABLE in between; raw agreement is "
    "reported only beside the majority-class rate, never against 0.5; "
    "(3) THE MAGNITUDE-MATCHED SUB-TEST, both arms restricted to |predicted log2fc| > 0.1, the band in "
    "which the first half was 28 of 28, registered so that a failing upward arm cannot be excused as "
    "small magnitudes; "
    "(4) THE RE-DERIVATION CHECK, the 44 pairs of 5c842ca must reappear with identical signs or the "
    "run is refused; and the downward arm re-read over all 116 answerable downward pairs, which "
    "removes the top-target selection the 41/44 was drawn under. "
    "Expected: the upward arm at or below chance. Requests: 0."
)

FALSIFIES = (
    "The pass of 5c842ca is restated over both signs and does NOT survive if either: balanced "
    "accuracy on the combined set has a bootstrap lower bound at or below 0.5; or the upward arm "
    "comes in at or below chance (sensitivity_up <= 0.5) -- including when balanced accuracy is high, "
    "because a high balanced accuracy carried entirely by the downward half means the layer reads "
    "'down', not direction. In either case the 0.9318 headline is restated as a statement about "
    "DOWNWARD effects only, in the same place and with the same prominence, and no subset is searched "
    "for in which the upward arm holds. A third way to lose it: the downward arm collapsing on the "
    "116 unselected pairs, which would show the 41/44 was the top-target selection rather than the "
    "layer."
)

EXPECTED = (
    "Stated plainly before scoring, and expected to be borne out: the upward arm comes back AT OR "
    "BELOW CHANCE, sensitivity_up somewhere between 0.2 and 0.45. The model is not shy of emitting "
    "'up' -- its marginal up-rate over the whole sweep is about 0.46, so it says 'up' on nearly half "
    "of all elements -- it simply says it in the wrong places, and the only three up-calls it made on "
    "the one-signed set were all wrong, 0 of 3. Balanced accuracy is therefore expected between 0.60 "
    "and 0.70, above 0.5 but carried by the downward half, which lands the combined result in the "
    "UNDECIDABLE band and makes the honest headline a downward-only one. If that is what happens it "
    "is reported as the result, not as a caveat."
)

LADDER = {
    "<= 0.55": (
        "the layer does not read direction at all. The 0.9318 was the measured side's constancy plus "
        "the magnitude selection the top-target rule imposed, and 5c842ca becomes a statement about "
        "downward effects only"
    ),
    "0.55 to 0.65": (
        "direction is faintly readable and not usable. UNDECIDABLE: the represses half of action "
        "stays unsupported, and every action word GenomeOS writes on a predicted rise is unbacked"
    ),
    "0.65 to 0.80": (
        "the layer reads direction on both signs with the upward half materially weaker. 5c842ca "
        "survives as a two-signed claim, but its headline number must be quoted as the balanced "
        "accuracy and never again as 0.9318"
    ),
    "> 0.80": (
        "the layer reads direction on both signs. 5c842ca survives in full, and the one-sidedness "
        "caveat on measurement 2 is discharged rather than narrowed"
    ),
}


def costing(name: str = HELDOUT) -> dict[str, Any]:
    """What buying the upward arm would have cost, by the coordinator's own rule, reproduced.

    One deletion request per distinct element behind an in-reach upward significant held-out pair.
    Reported whatever it comes to, and reported beside the fact that it need not be spent.
    """
    up = [r for r in signed(rows(name)) if float(r["EffectSize"]) > 0 and in_reach(r)]
    out: dict[str, Any] = {}
    for label, cells in (("K562", ("K562",)), ("K562+GM12878", CELLS)):
        s = [r for r in up if r["CellType"] in cells]
        els = {(r["chrom"], r["chromStart"], r["chromEnd"]) for r in s}
        out[label] = {"pairs": len(s), "requests": len(els)}
    out["spent"] = 0
    out["why_zero"] = (
        "the per-element response cache of the same sweep was written with threshold=0.0, so every "
        "gene in the scorer's 1 Mb window already carries a signed per-cell log2 fold change. The "
        "one-target-per-element limit belongs to the compact all_elements table, not to what was "
        "measured. These requests would have re-bought an answer that is on disk"
    )
    return out


def cached_genes(chrom: str, element_id: str, archive: dict[str, Any]) -> list[dict[str, Any]]:
    """Every gene the sweep scored for this element: the archive first, then a loose per-element file."""
    hit = archive.get(element_id)
    if hit is None:
        p = CACHE / chrom / f"{element_id}.json"
        if p.exists():
            hit = json.loads(p.read_text())
    return (hit or {}).get("genes") or []


def value_for_gene(
    els: list[dict[str, Any]], chrom: str, gene: str, cell: str, archive: dict[str, Any]
) -> dict[str, Any] | None:
    """The sweep's signed predicted change for THIS gene in THIS cell, whether or not it is the top target.

    The same resolution rule the first half registered, moved from the compact table to the cache:
    among every overlapping scored element that carries a value for the gene on the cell's own track,
    take the largest absolute change, ties broken by element id. Every match is kept so that pairs
    whose overlapping elements disagree in sign are counted rather than hidden by the choice.
    """
    matches = []
    for e in els:
        for g in cached_genes(chrom, e["id"], archive):
            if g["gene"] != gene:
                continue
            v = (g.get("by_cell") or {}).get(cell)
            if v is not None:
                matches.append({"value": float(v), "element": e["id"]})
    if not matches:
        return None
    best = min(matches, key=lambda m: (-abs(m["value"]), m["element"]))
    signs = {1 if m["value"] > 0 else (-1 if m["value"] < 0 else 0) for m in matches}
    return {**best, "matches": len(matches), "matches_disagree": len(signs) > 1}


def top_target_of(els: list[dict[str, Any]], gene: str, cell: str) -> bool:
    """Whether the compact table could have answered this pair: the old route's filter, kept as a stratum."""
    return any(
        e.get(k) and e[k].get("gene") == gene and (e.get(bc) or {}).get(cell) is not None
        for e in els
        for k, bc in (("predicted_coding", "predicted_coding_by_cell"), ("predicted", "predicted_by_cell"))
    )


def collect(table: list[dict[str, str]], cells: tuple[str, ...] = CELLS) -> dict[str, Any]:
    """Join one arm's signed pairs to the per-element cache, a chromosome at a time."""
    pairs = signed(table)
    by_chrom: dict[str, list[dict[str, str]]] = defaultdict(list)
    for r in pairs:
        by_chrom[r["chrom"]].append(r)
    answered: list[dict[str, Any]] = []
    unanswered: list[dict[str, Any]] = []
    strata: dict[str, int] = defaultdict(int)
    for chrom in sorted(by_chrom):
        relevant = [r for r in by_chrom[chrom] if r["CellType"] in cells and in_reach(r)]
        for r in by_chrom[chrom]:
            cell = r["CellType"]
            sign = "up" if float(r["EffectSize"]) > 0 else "down"
            strata[f"signed::{cell}::{sign}"] += 1
            if cell not in cells:
                strata[f"refused_no_model_line::{cell}::{sign}"] += 1
            elif not in_reach(r):
                strata[f"out_of_reach::{cell}::{sign}"] += 1
        p = ELEMENTS / f"{chrom}.json"
        if not relevant or not p.exists():
            continue
        els_all = sorted(json.loads(p.read_text()), key=lambda e: e["start"])
        starts = [e["start"] for e in els_all]
        arc_p = CACHE / f"{chrom}.json.gz"
        archive: dict[str, Any] = {}
        if arc_p.exists():
            with gzip.open(arc_p, "rt") as fh:
                archive = json.load(fh)
        for r in relevant:
            cell, gene = r["CellType"], r["measuredGeneSymbol"]
            effect = float(r["EffectSize"])
            sign = "up" if effect > 0 else "down"
            ov = overlapping(els_all, starts, int(r["chromStart"]), int(r["chromEnd"]))
            if not ov:
                strata[f"element_not_in_the_sweep::{cell}::{sign}"] += 1
                continue
            strata[f"covered::{cell}::{sign}"] += 1
            row = {
                "cell": cell,
                "chrom": chrom,
                "element": [int(r["chromStart"]), int(r["chromEnd"])],
                "gene": gene,
                "measured": round(effect, 4),
                "measured_sign": sign,
                "distance_to_tss": int(float(r["distanceToTSS"])),
                "dataset": r["Dataset"],
                "was_top_target": top_target_of(ov, gene, cell),
            }
            m = value_for_gene(ov, chrom, gene, cell, archive)
            if m is None:
                strata[f"gene_absent_from_the_scored_window::{cell}::{sign}"] += 1
                unanswered.append(row)
                continue
            strata[f"answerable::{cell}::{sign}"] += 1
            strata[f"newly_answerable::{cell}::{sign}"] += not row["was_top_target"]
            answered.append(
                {
                    **row,
                    "predicted": round(m["value"], 4),
                    "predicted_sign": "down" if m["value"] < 0 else ("up" if m["value"] > 0 else "zero"),
                    "matches": m["matches"],
                    "matches_disagree": m["matches_disagree"],
                }
            )
        del els_all, starts, archive
    return {
        "answered": answered,
        "unanswered": unanswered,
        "strata": dict(sorted(strata.items())),
        "signed": len(pairs),
    }


def usable(answered: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """A predicted value of exactly 0.0 is no direction; such pairs leave the denominator and are counted."""
    return [r for r in answered if r["predicted_sign"] != "zero"]


def sensitivity(answered: list[dict[str, Any]], sign: str) -> dict[str, Any]:
    """Agreement within one measured sign: the share of measured-`sign` pairs the model also calls `sign`."""
    arm = [r for r in usable(answered) if r["measured_sign"] == sign]
    return wilson(sum(r["predicted_sign"] == sign for r in arm), len(arm))


def balanced(answered: list[dict[str, Any]], seed: int = SEED, reps: int = BOOTSTRAPS) -> dict[str, Any]:
    """Balanced accuracy and its stratified bootstrap interval.

    The mean of the two per-sign sensitivities. Chance is 0.5 for any class mix, and a caller that
    emits one sign and nothing else scores exactly 0.5 whichever sign it picks, which is the whole
    reason this and not raw agreement is the registered primary. The bootstrap resamples within each
    sign arm separately, so the class mix is held fixed and only the within-arm uncertainty is drawn.
    """
    u = usable(answered)
    down = [r for r in u if r["measured_sign"] == "down"]
    up = [r for r in u if r["measured_sign"] == "up"]
    if not down or not up:
        return {"rate": None, "ci95": None, "n_down": len(down), "n_up": len(up), "both_signs": False}
    hit = {
        "down": [r["predicted_sign"] == "down" for r in down],
        "up": [r["predicted_sign"] == "up" for r in up],
    }
    point = (sum(hit["down"]) / len(down) + sum(hit["up"]) / len(up)) / 2
    rng = random.Random(seed)
    draws = []
    for _ in range(reps):
        a = sum(rng.choice(hit["down"]) for _ in down) / len(down)
        b = sum(rng.choice(hit["up"]) for _ in up) / len(up)
        draws.append((a + b) / 2)
    draws.sort()
    lo = draws[int(0.025 * reps)]
    hi = draws[min(reps - 1, int(0.975 * reps))]
    return {
        "rate": round(point, 4),
        "ci95": [round(lo, 4), round(hi, 4)],
        "n_down": len(down),
        "n_up": len(up),
        "both_signs": True,
        "method": f"stratified bootstrap, {reps} draws, seed {seed}",
    }


def raw_agreement(answered: list[dict[str, Any]]) -> dict[str, Any]:
    """Plain sign agreement, and beside it the majority-class rate that makes 0.5 the wrong comparison."""
    u = usable(answered)
    k = sum(r["predicted_sign"] == r["measured_sign"] for r in u)
    down = sum(r["measured_sign"] == "down" for r in u)
    up = len(u) - down
    return {
        **wilson(k, len(u)),
        "measured_down": down,
        "measured_up": up,
        "predicted_down": sum(r["predicted_sign"] == "down" for r in u),
        "predicted_up": sum(r["predicted_sign"] == "up" for r in u),
        "majority_class_rate": round(max(down, up) / len(u), 4) if u else None,
        "majority_class_note": (
            "the constant-sign caller's score on this set. Raw agreement is compared to THIS, not to "
            "0.5; balanced accuracy is the statistic that compares to 0.5 honestly"
        ),
        "predicted_zero_excluded": len(answered) - len(u),
    }


def marginal_up_rate(
    marginal: dict[str, Any], answered: list[dict[str, Any]], cells: tuple[str, ...]
) -> float | None:
    """The constant-sign caller's score on the UPWARD arm: one minus the sweep's down-rate, cell-weighted."""
    arm = [r for r in usable(answered) if r["measured_sign"] == "up"]
    if not arm:
        return None
    w: dict[str, int] = defaultdict(int)
    for r in arm:
        w[r["cell"]] += 1
    num = sum(w[c] * (1 - (marginal[c]["rate"] or 0)) for c in w if c in cells)
    return num / sum(w.values())


def judge_up(head: dict[str, Any], marginal: float | None) -> dict[str, Any]:
    """The upward arm, by the same two margins the first half registered for the downward one."""
    rate, ci = head.get("rate"), head.get("ci95")
    if rate is None or not ci:
        return {"verdict": "refused", "why": "no upward pair with a non-zero predicted sign"}
    if rate <= CHANCE:
        return {
            "verdict": "failed",
            "why": (
                f"upward sensitivity {rate} is at or below chance {CHANCE}. The layer cannot call an "
                f"increase, and the passing direction headline is a statement about downward effects only"
            ),
        }
    if rate <= CHANCE + CHANCE_MARGIN:
        return {
            "verdict": "failed",
            "why": f"upward sensitivity {rate} does not clear {CHANCE + CHANCE_MARGIN}",
        }
    if marginal is None:
        return {"verdict": "undecidable", "why": "no marginal up-rate to compare against"}
    if rate - marginal >= MARGINAL_MARGIN and ci[0] > marginal:
        return {
            "verdict": "passed",
            "why": (
                f"upward sensitivity {rate} clears both {CHANCE + CHANCE_MARGIN} and the "
                f"marginal up-rate {round(marginal, 4)}"
            ),
        }
    return {
        "verdict": "undecidable",
        "why": (
            f"upward sensitivity {rate} clears {CHANCE + CHANCE_MARGIN} but not the model's own "
            f"up-habit {round(marginal, 4)} by {MARGINAL_MARGIN} with the interval clear of it"
        ),
    }


def judge_combined(ba: dict[str, Any]) -> dict[str, Any]:
    """The combined set, on balanced accuracy, against a chance level that is 0.5 for real this time."""
    rate, ci = ba.get("rate"), ba.get("ci95")
    if rate is None or not ci:
        return {"verdict": "refused", "why": "the answerable set does not carry both measured signs"}
    band = next(
        (
            k
            for k, lo, hi in (
                ("<= 0.55", 0, 0.55),
                ("0.55 to 0.65", 0.55, 0.65),
                ("0.65 to 0.80", 0.65, 0.80),
                ("> 0.80", 0.80, 2),
            )
            if lo < rate <= hi
        ),
        "> 0.80",
    )
    reading = LADDER[band]
    if ci[0] <= CHANCE:
        return {
            "verdict": "failed",
            "why": f"balanced accuracy {rate}, 95% lower bound {ci[0]} at or below {CHANCE}; {reading}",
            "band": band,
        }
    if rate >= CHANCE + CHANCE_MARGIN:
        return {
            "verdict": "passed",
            "why": (
                f"balanced accuracy {rate} >= {CHANCE + CHANCE_MARGIN} with lower bound "
                f"{ci[0]} above {CHANCE}; {reading}"
            ),
            "band": band,
        }
    return {
        "verdict": "undecidable",
        "why": f"balanced accuracy {rate} is above chance but under {CHANCE + CHANCE_MARGIN}; {reading}",
        "band": band,
    }


def by_magnitude(answered: list[dict[str, Any]]) -> dict[str, Any]:
    """The first half's sharpest finding, re-read over the new pairs: agreement against |predicted log2fc|.

    It found every error below 0.055 and 28 of 28 above 0.1. That is restated here as a prediction
    which can fail, per measured sign, so that the ladder is not only re-read on the arm that made it.
    """
    out: dict[str, Any] = {}
    for sign in ("down", "up", "both"):
        arm = [r for r in usable(answered) if sign == "both" or r["measured_sign"] == sign]
        bins: dict[str, dict[str, int]] = {}
        for r in arm:
            lo = str(MAGNITUDE_EDGES[magnitude_bin(r["predicted"])])
            slot = bins.setdefault(lo, {"agree": 0, "pairs": 0})
            slot["pairs"] += 1
            slot["agree"] += r["predicted_sign"] == r["measured_sign"]
        weak = [r for r in arm if abs(r["predicted"]) <= WEAK]
        strong = [r for r in arm if abs(r["predicted"]) > STRONG]
        out[sign] = {
            "bins": dict(sorted(bins.items(), key=lambda t: float(t[0]))),
            f"at_or_below_{WEAK}": wilson(
                sum(r["predicted_sign"] == r["measured_sign"] for r in weak), len(weak)
            ),
            f"above_{STRONG}": wilson(
                sum(r["predicted_sign"] == r["measured_sign"] for r in strong), len(strong)
            ),
            "largest_error_magnitude": round(
                max(
                    (abs(r["predicted"]) for r in arm if r["predicted_sign"] != r["measured_sign"]),
                    default=0.0,
                ),
                4,
            ),
        }
    return out


def magnitude_matched(answered: list[dict[str, Any]]) -> dict[str, Any]:
    """The registered sub-test: both arms above |0.1|, the band the first half did not miss in.

    The upward pairs are not selected for being the element's top target, so their predicted
    magnitudes run smaller than the 44's. If the upward arm fails, small magnitudes are a candidate
    explanation. This removes it: if the upward arm fails even where the downward arm was 28 of 28,
    the failure is about sign and not about magnitude, and it is registered before either is known.
    """
    strong = [r for r in usable(answered) if abs(r["predicted"]) > STRONG]
    return {
        "threshold": STRONG,
        "down": sensitivity(strong, "down"),
        "up": sensitivity(strong, "up"),
        "balanced": balanced(strong),
        "note": (
            "registered in advance so that a failing upward arm cannot be excused as small magnitudes, "
            "and a passing one cannot be credited to them"
        ),
    }


def rederivation_check(answered: list[dict[str, Any]], prior_path: Path) -> dict[str, Any]:
    """The 44 of 5c842ca must reappear here with the same signs, or the run is refused.

    The cache route and the compact-table route read the same numbers by different paths. Where the
    measured gene IS the element's top target the two must agree exactly. A mismatch means one of the
    two readers is wrong, and there is no version of that in which this run's numbers mean anything.
    """
    if not prior_path.exists():
        return {"checked": False, "why": "the first half's result file is not on disk"}
    prior = json.loads(prior_path.read_text())["answered"]
    here = {(r["chrom"], tuple(r["element"]), r["gene"], r["cell"]): r for r in answered}
    missing, changed = [], []
    for p in prior:
        key = (p["chrom"], tuple(p["element"]), p["gene"], p["cell"])
        q = here.get(key)
        if q is None:
            missing.append(key)
        elif q["predicted_sign"] != p["predicted_sign"]:
            changed.append({"pair": key, "then": p["predicted"], "now": q["predicted"]})
    return {
        "checked": True,
        "prior_pairs": len(prior),
        "prior_headline": PRIOR,
        "missing_here": len(missing),
        "sign_changed": len(changed),
        "changed": changed[:10],
        "passes": not missing and not changed,
    }


def direction_both(
    cells: tuple[str, ...] = CELLS, prior: Path = Path("data/results/crispri_direction.json")
) -> dict[str, Any]:
    """The whole registered measurement, for zero requests."""
    held = collect(rows(HELDOUT), cells)
    train = collect(rows(TRAINING), cells)
    answered = held["answered"]
    marginal = marginal_down_rate(cells)

    up_head = sensitivity(answered, "up")
    down_head = sensitivity(answered, "down")
    up_marginal = marginal_up_rate(marginal, answered, cells)
    ba = balanced(answered)

    # the downward arm without the top-target selection the 41/44 was drawn under
    selected = [r for r in answered if r["was_top_target"]]
    unselected = [r for r in answered if not r["was_top_target"]]

    verdict_up = judge_up(up_head, up_marginal)
    verdict_combined = judge_combined(ba)
    survives = (
        verdict_combined["verdict"] == "passed"
        and (up_head.get("rate") or 0) > CHANCE
        and verdict_up["verdict"] != "failed"
    )
    return {
        "result": "crispri_direction_both",
        "evidence": EVIDENCE,
        "preregistered": PREREGISTERED,
        "falsifies": FALSIFIES,
        "expected_before_scoring": EXPECTED,
        "ladder": LADDER,
        "requests": 0,
        "costing_of_the_route_not_taken": costing(),
        "cells": list(cells),
        "upward_arm": {
            **up_head,
            "marginal_up_rate": None if up_marginal is None else round(up_marginal, 4),
            **verdict_up,
        },
        "downward_arm": down_head,
        "combined": {
            "balanced_accuracy": ba,
            "raw_agreement": raw_agreement(answered),
            **verdict_combined,
        },
        "pass_of_5c842ca_survives_restatement_over_both_signs": survives,
        "why_survives": (
            "both the combined set and the upward arm clear their registered bars"
            if survives
            else (
                "restated over both signs it does not survive as a direction claim: "
                + verdict_combined["why"]
                + " | upward arm: "
                + verdict_up["why"]
            )
        ),
        "magnitude_matched_sub_test": magnitude_matched(answered),
        "magnitude_ladder": by_magnitude(answered),
        "rederivation_of_the_44": rederivation_check(answered, prior),
        "top_target_selection": {
            "note": (
                "the 41/44 was drawn from pairs whose measured gene is the element's top predicted "
                "target, which selects for a large predicted change. These are the same pairs split "
                "by that selection, both signs pooled"
            ),
            "was_top_target": raw_agreement(selected),
            "newly_answerable": raw_agreement(unselected),
            "balanced_newly_answerable": balanced(unselected),
        },
        "per_cell": {
            c: {
                "up": sensitivity([r for r in answered if r["cell"] == c], "up"),
                "down": sensitivity([r for r in answered if r["cell"] == c], "down"),
                "balanced": balanced([r for r in answered if r["cell"] == c]),
                "raw": raw_agreement([r for r in answered if r["cell"] == c]),
            }
            for c in cells
            if any(r["cell"] == c for r in answered)
        },
        "marginal_down_rate": {
            c: {k: v for k, v in marginal[c].items() if k != "by_magnitude"} for c in cells
        },
        "training_arm_fitted_on_not_a_test": {
            "up": sensitivity(train["answered"], "up"),
            "down": sensitivity(train["answered"], "down"),
            "balanced": balanced(train["answered"]),
            "strata": train["strata"],
            "note": "the 2026-09-16 lane fitted logistic weights on this arm; reported for shape only",
        },
        "coverage": {
            "heldout": held["strata"],
            "heldout_signed_pairs": held["signed"],
            "unanswered_pairs": len(held["unanswered"]),
            "unanswered": held["unanswered"],
        },
        "answered": answered,
    }
