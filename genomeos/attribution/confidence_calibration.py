# SPDX-License-Identifier: AGPL-3.0-or-later
"""Do the project's own confidences mean anything? The measured layer as a reliability table.

GenomeOS puts an evidence kind and a confidence on every stated fact, and the Evidence explorer
counts them: 960,096 compiled facts, 880,754 of them `predicted`, mean confidence around 0.6. The
promise a reader takes from that is a probability - a fact stated at 0.4 is right about 40% of the
time - and nobody has ever checked it. This module checks it, on the only elements where it can be
checked at all.

**What the number on a compiled fact actually is.** `attribution/compile.py` writes

    confidence = round(min(PREDICTED_CAP, max(0.05, |log2 fold change|)), 2)

on the element block and on the rule beside it. It is a clipped effect size, exactly the quantity
`target_calibration` found unusable as a probability for a different question (the probability that a
CRISPRi screen calls a pair regulated). Nothing has ever said what it means for *this* question, which
is the compiled claim itself: deleting this element moves this gene.

**What a measurement says about that claim, and it is not the same thing in four assays.** Since
2026-09-17 `attribution/measured.py` attaches real outcomes to compiled elements. Each assay answers
a different question and its "agreement" is a different event, so the table is reported per assay and
`WHY_POOLING_IS_REFUSED` says why one number over all four would be arithmetic rather than a reading.
Only CRISPRi on the named gene asks the compiled claim's own question, so only there is the *level* of
the confidence interpretable at all; elsewhere the table can say whether the number ORDERS the facts
and nothing about what it means (`WHAT_AGREEMENT_MEANS`, `level_interpretable`).

**Three verdicts, named in code before any rate was computed** (`UNINFORMATIVE`,
`ORDERED_BUT_MISCALIBRATED`, `CALIBRATED`, and `LEVEL_WITHHELD` for the assays whose agreement is not
the compiled claim). The bands are not chosen here either: `CONFIDENCE_BANDS` is imported unchanged
from `target_calibration`, where it was fixed on 2026-09-17 for a different quantity and before this
lane existed. A band the rates were drawn to fit is not a band.

**Two denominators, never one.** The layer separates the agreement rate over pairs where the predicted
gene was tested (128 elements, the number a careless reader quotes) from the rate over all matched
elements (1,505), and both travel through every table here. The wide one is a mixture of "the
prediction was wrong" and "the screen never asked", so it is bounded above by the screen's own choice
of genes and can never be calibrated against anything; it exists so the narrow one cannot stand alone.

**Scope is the result, not a caveat on it.** 4.33% of compiled elements are measured at all, and the
measured slice is longer, more GC-rich and closer to a promoter than the rest. A reliability table
over it is a statement about measured elements. `SCOPE` and `FALSIFIES_TRANSFER` say so and say what
would show the transfer failing.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field
from statistics import median
from typing import Any

from genomeos.attribution import measured
from genomeos.attribution.compile import PREDICTED_CAP
from genomeos.attribution.measured import AGREES, ASSAYS, DISAGREES
from genomeos.attribution.panel_background import TSS_EDGES
from genomeos.attribution.target_calibration import CONFIDENCE_BANDS, wilson
from genomeos.compare import Strata, imbalance, input_presence, standardised
from genomeos.evidence import WEAK

# --- what is being tested -------------------------------------------------------------------------
STATED_CONFIDENCE_FORMULA = (
    "round(min(PREDICTED_CAP, max(0.05, |log2 fold change|)), 2), as attribution/compile.py writes it "
    f"on the element block and on the rule beside it; PREDICTED_CAP = {PREDICTED_CAP}, so no compiled "
    "element can state more than that however large the predicted effect"
)
CLAIM_UNDER_TEST = (
    "the compiled element block and its rule say that deleting this element moves this gene, and the "
    "number beside them is read by every reader of the Evidence explorer as the probability that it "
    "does. The question here is whether that reading survives contact with a measurement"
)

# --- the bands, and they were not chosen here -------------------------------------------------------
# imported unchanged from target_calibration, where they were fixed for a different quantity on
# 2026-09-17, before this lane existed. Half-open, [lo, hi), exactly as that module reads them.
BANDS = CONFIDENCE_BANDS
MIN_FOR_A_BAND = 30  # below this a band is printed and described, never judged
MIN_BANDS_FOR_A_VERDICT = 3  # fewer judged bands than this and the assay gets no verdict at all
MIN_SHARE_CONSISTENT = 0.7  # target_calibration's 7-of-10, as a share because the band count varies
# the satmut lane's bar, imported from the lane that set it rather than restated: two copies of a
# number that cite each other are folklore, and this one cited a home that did not exist
MIN_FOR_A_COMPARISON = measured.MIN_FOR_A_COMPARISON
HIGH_CONFIDENCE = WEAK  # 0.5, the line the Program tab and the Evidence explorer already draw

# --- the failure patterns, stated before the rates were computed -------------------------------------
UNINFORMATIVE = (
    "the observed agreement rate does not rise with the stated confidence: the highest judged band's "
    "rate is at or below the lowest judged band's, or their 95% Wilson intervals overlap. If this is "
    "the pattern, the number on a compiled fact orders nothing, and a reader who preferred a 0.6 fact "
    "to a 0.2 one gained nothing by it"
)
ORDERED_BUT_MISCALIBRATED = (
    "the rate rises - the highest judged band is above the lowest and their intervals are disjoint - "
    "but the band's mean stated confidence lies outside the 95% Wilson interval of its observed rate "
    f"in more than {1 - MIN_SHARE_CONSISTENT:.0%} of the judged bands. The number ranks and its level "
    "is wrong. This is the ordinary outcome for a scoring system read as a probability, and it is a "
    "useful result: a rank is worth having. The offset is reported as the median signed gap, and a "
    "level is a property of the population and not of the scorer - target_calibration's curve "
    "transferred while its level did not, because the base rate of the new population differed"
)
CALIBRATED = (
    "the rate rises AND the band's mean stated confidence lies inside the 95% Wilson interval of its "
    f"observed rate in at least {MIN_SHARE_CONSISTENT:.0%} of the judged bands. Then the number may be "
    "read as a probability, within the scope below and nowhere else"
)
LEVEL_WITHHELD = (
    "the assay's agreement is not the compiled claim's own event, so the ordering may be read and the "
    "level may not: a rate of 0.19 against a stated 0.25 is not a miscalibration when the rate being "
    "measured is 'a 200 bp sequence drives a reporter' and the claim is 'deleting the element in its "
    "own chromosome moves this gene'. The level verdict is withheld rather than computed and ignored"
)
PATTERNS = {
    "uninformative": UNINFORMATIVE,
    "ordered_but_miscalibrated": ORDERED_BUT_MISCALIBRATED,
    "calibrated": CALIBRATED,
    "level_withheld_ordering_only": LEVEL_WITHHELD,
}

# --- what "agreement" is, assay by assay ------------------------------------------------------------
WHAT_AGREEMENT_MEANS: dict[str, dict[str, Any]] = {
    "crispri": {
        "question": "silencing this element in its own chromosome, does the predicted gene fall?",
        "is_the_compiled_claim": True,
        "level_interpretable": True,
        "why": (
            "the screen perturbs the element where it lives and measures the very gene the deletion "
            "named, which is the compiled claim asked directly. A negative here is a strong "
            "contradiction and the rate may be read as a probability"
        ),
    },
    "lentimpra": {
        "question": "out of its chromosome, does a 200 bp copy of this sequence drive a reporter?",
        "is_the_compiled_claim": False,
        "level_interpretable": False,
        "why": (
            "episomal: it measures the sequence and not the locus, and it never sees the predicted "
            "gene. A silence is a weak contradiction, and the base rate of 'active in a reporter' has "
            "no reason to equal the base rate of 'deleting this moves that gene'. Ordering only"
        ),
    },
    "vista": {
        "question": "in a transgenic mouse embryo at e11.5, is this sequence an enhancer?",
        "is_the_compiled_claim": False,
        "level_interpretable": False,
        "why": (
            "a different organism, one developmental stage and a reporter construct. A negative is a "
            "negative in a mouse embryo, which is not a measurement of a human cell's transcription, "
            "and the assay's own positives are enriched by how its sequences were chosen. Ordering only"
        ),
    },
    "satmut": {
        "question": "which single bases inside this element change a reporter's activity?",
        "is_the_compiled_claim": False,
        "level_interpretable": False,
        "why": (
            "it cannot disagree at all: substituting one base at a time never deletes the element and "
            "never measures the predicted gene, so its agreement rate is 1.0 by construction wherever "
            "any base is functional. A reliability table over it would measure nothing twice over, and "
            "at four matched elements it would not measure it. See measured.SATMUT_CANNOT_DISAGREE"
        ),
    },
}
WHY_POOLING_IS_REFUSED = (
    "one agreement rate over all four assays would be a rate whose event changes with the denominator. "
    "A CRISPRi negative on the named gene contradicts the compiled claim; a lentiMPRA silence says a "
    "200 bp copy does not drive a reporter out of its chromosome; a VISTA negative says a mouse embryo "
    "did not stain at e11.5; saturation mutagenesis cannot disagree at all and so contributes "
    "agreements and never disagreements. Pooling them would also pool their sizes - lentiMPRA's 17,869 "
    "matched elements against CRISPRi's 128 tested pairs is 140 to 1 - so the pooled number would be "
    "the reporter's number wearing the perturbation's name. Every table here is per assay and the "
    "result holds no pooled rate"
)

# --- scope -------------------------------------------------------------------------------------------
SCOPE = (
    "This is a reliability table over MEASURED compiled elements, and that is not a random sample of "
    "the compiled genome. 4.33% of compiled elements are measured at all (19,072 of 440,377) and only "
    "5.98% were ever eligible, because an assay's footprint had to reach them. The measured slice is "
    "longer, more GC-rich and closer to a coding TSS than the rest - the assays chose candidate "
    "regulatory sequence, which is the same property the deletion scores highly - so the bands' rates "
    "are conditional on an element having been chosen for an assay. The table says how often the "
    "project's confidence is borne out ON MEASURED ELEMENTS. It says nothing directly about the "
    "414,053 elements no assay ever covered, and the imbalance between the two arms is printed in the "
    "result rather than described here, so the size of the extrapolation is a number"
)
FALSIFIES_TRANSFER = (
    "The transfer from measured elements to the compiled genome fails if any of these is observed. "
    "(a) A screen that tiles a region without choosing candidate sequence - an unbiased tiling rather "
    "than a cCRE list - finds the top band's rate far below what it is here, which is what selection "
    "on testability looks like from outside. (b) The agreement rate moves with the measured slice's "
    "own covariates once they are held fixed: if the standardised difference between the high and low "
    "confidence arms collapses toward zero when length, GC and distance to a coding TSS are matched, "
    "then the band was reading the covariate and not the confidence, and on unmeasured elements the "
    "covariate distribution is different. (c) A new assay extends the footprint into elements that "
    "look unlike the measured ones - far from a promoter, AT-rich, short - and the bands do not hold "
    "there. (d) The prevalence shift target_calibration already found: a population whose base rate of "
    "agreement differs from this one's needs its own intercept, so a band quoted without its "
    "population's rate beside it has stopped being a measurement of anything"
)
TWO_DENOMINATORS = (
    "over_all_matched_elements is a mixture of 'the prediction was wrong' and 'the screen never asked "
    "about that gene', so it is bounded above by the screen's own choice of genes and cannot be "
    "calibrated against a confidence; it is carried beside the narrow rate so that the narrow rate - "
    "computed over a subset the prediction itself selected - can never stand alone. For lentiMPRA and "
    "VISTA the two coincide by construction, since every matched element is tested, and the result "
    "says so rather than printing one number and letting it look like a choice"
)


# --- the elements ------------------------------------------------------------------------------------
@dataclass(slots=True)
class Element:
    """One compiled element: the confidence it states, its covariates, and what each assay said.

    `eligible` is the second denominator - the assays whose footprint touched it at one shared base -
    and it is carried per element so that every band can print compiled, eligible and measured in that
    order rather than measured alone.
    """

    id: str
    chrom: str
    stated_confidence: float
    length: float
    gc: float | None
    nearest_coding_tss: float | None
    eligible: frozenset[str] = frozenset()
    verdicts: dict[str, str] = field(default_factory=dict)

    @property
    def band(self) -> str:
        return band_of(self.stated_confidence)

    @property
    def covariates_present(self) -> bool:
        return self.gc is not None and self.nearest_coding_tss is not None


def stated_confidence(element: dict[str, Any]) -> float | None:
    """The number the compiler writes on this element's block and rule, recomputed identically.

    Recomputed rather than re-read, because the committed programs are generated text that a rebuild
    replaces; the formula in `compile.py` is the fact. A element with no predicted target states
    nothing and is None, not zero.
    """
    pc = element.get("predicted_coding") or {}
    if not pc.get("gene"):
        return None
    raw = pc.get("confidence")
    if raw is None:
        raw = abs(float(pc.get("log2_fold_change") or 0.0))
    return round(min(PREDICTED_CAP, max(0.05, float(raw))), 2)


def band_of(c: float) -> str:
    """The band a confidence falls in, half-open [lo, hi), as target_calibration reads its own."""
    for lo, hi in BANDS:
        if lo <= c < hi:
            return f"{lo:g}-{hi if hi <= 1 else 1:g}"
    return "1"


BAND_NAMES = tuple(f"{lo:g}-{hi if hi <= 1 else 1:g}" for lo, hi in BANDS)
BAND_EDGES = {f"{lo:g}-{hi if hi <= 1 else 1:g}": (lo, hi) for lo, hi in BANDS}


def elements_of(
    chrom: str,
    attributed: list[dict[str, Any]],
    layer: measured.Layer,
    measured_rows: list[dict[str, Any]],
    covariates: dict[str, dict[str, Any]],
) -> list[Element]:
    """Every compiled element of one chromosome, with its confidence, covariates, footprint and verdicts."""
    by_id = {r["id"]: r for r in measured_rows}
    out: list[Element] = []
    for e in attributed:
        conf = stated_confidence(e)
        if conf is None:
            continue
        cov = covariates.get(e["id"]) or {}
        row = by_id.get(e["id"])
        verdicts = {}
        if row is not None:
            for assay in ASSAYS:
                v = row["agreement"].get(assay)
                if v is not None and assay in row["measured"]:
                    verdicts[assay] = v
        gc = cov.get("gc")
        tss = cov.get("nearest_coding_tss")
        out.append(
            Element(
                id=e["id"],
                chrom=chrom,
                stated_confidence=conf,
                length=float(e["end"] - e["start"]),
                gc=None if gc is None else float(gc),
                nearest_coding_tss=None if tss is None else float(tss),
                eligible=frozenset(layer.touched_by(e["start"], e["end"])),
                verdicts=verdicts,
            )
        )
    return out


# --- the reliability table ---------------------------------------------------------------------------
def _medians(rows: list[Element]) -> dict[str, float | None]:
    """The three covariates every compared group in this project prints, as medians."""
    gc = [r.gc for r in rows if r.gc is not None]
    tss = [r.nearest_coding_tss for r in rows if r.nearest_coding_tss is not None]
    return {
        "median_length": round(median([r.length for r in rows]), 1) if rows else None,
        "median_gc": round(median(gc), 4) if gc else None,
        "median_nearest_coding_tss": round(median(tss), 1) if tss else None,
        "elements_without_covariates": sum(1 for r in rows if not r.covariates_present),
    }


def _labels(rows: list[Element], assay: str, narrow: bool) -> list[Element]:
    """The rows one denominator is computed over.

    narrow: the assay actually asked the compiled claim's question of this element. For CRISPRi that
    is the pairs where the predicted gene was tested; `predicted_gene_not_tested` is excluded, because
    an unasked question is not a wrong answer. For lentiMPRA and VISTA every matched element is tested,
    so the narrow and the wide set are the same rows and the result says so by name.
    """
    matched = [r for r in rows if assay in r.verdicts]
    if not narrow:
        return matched
    return [r for r in matched if r.verdicts[assay] in (AGREES, DISAGREES)]


def table(
    population: list[Element], assay: str, narrow: bool, bands: tuple[str, ...] = BAND_NAMES
) -> dict[str, Any]:
    """One reliability table: per band the stated confidence against the observed agreement rate.

    Three denominators per band and in that order - the compiled elements whose confidence falls in
    the band, those an assay's footprint ever touched, and those this assay measured - so a band's rate
    is read after, and never instead of, the coverage that produced it. Each band also prints length,
    GC and median distance to a coding TSS, because a band is a group being compared with the others.
    """
    scored = _labels(population, assay, narrow)
    by_band: dict[str, list[Element]] = defaultdict(list)
    for r in scored:
        by_band[r.band].append(r)
    compiled_by_band: dict[str, list[Element]] = defaultdict(list)
    for r in population:
        compiled_by_band[r.band].append(r)

    rows: list[dict[str, Any]] = []
    for name in bands:
        lo, hi = BAND_EDGES[name]
        here = by_band.get(name, [])
        all_here = compiled_by_band.get(name, [])
        eligible_here = sum(1 for r in all_here if assay in r.eligible)
        n = len(here)
        k = sum(1 for r in here if r.verdicts[assay] == AGREES)
        mean_stated = sum(r.stated_confidence for r in here) / n if n else None
        w_lo, w_hi = wilson(k, n) if n else (None, None)
        rows.append(
            {
                "band": name,
                "stated_range": [lo, min(hi, 1.0)],
                # 1. compiled, 2. eligible, 3. measured - always in this order
                "compiled_elements_in_the_band": len(all_here),
                "in_this_assays_footprint": eligible_here,
                "measured_by_this_assay": n,
                "coverage_of_the_band": round(n / len(all_here), 5) if all_here else None,
                "share_of_the_eligible_raised": (round(n / eligible_here, 4) if eligible_here else None),
                "mean_stated_confidence": round(mean_stated, 4) if mean_stated is not None else None,
                "agrees": k,
                "observed": round(k / n, 4) if n else None,
                "ci95": [round(w_lo, 4), round(w_hi, 4)] if n else None,
                "stated_inside_the_interval": (
                    bool(w_lo <= mean_stated <= w_hi) if n and mean_stated is not None else None
                ),
                "gap_observed_minus_stated": (
                    round(k / n - mean_stated, 4) if n and mean_stated is not None else None
                ),
                "judged": n >= MIN_FOR_A_BAND,
                **_medians(here),
            }
        )
    judged = [r for r in rows if r["judged"]]
    n_total = sum(r["measured_by_this_assay"] for r in rows)
    k_total = sum(r["agrees"] for r in rows)
    ece = (
        round(
            sum(
                abs(r["gap_observed_minus_stated"]) * r["measured_by_this_assay"] / n_total
                for r in rows
                if r["measured_by_this_assay"]
            ),
            4,
        )
        if n_total
        else None
    )
    return {
        "denominator": (
            "the predicted gene was tested by this assay"
            if narrow
            else "every element this assay measured, whether or not it asked about the predicted gene"
        ),
        "denominator_name": (
            "where_the_predicted_gene_was_tested" if narrow else "over_all_matched_elements"
        ),
        "elements": n_total,
        "agrees": k_total,
        "observed_overall": round(k_total / n_total, 4) if n_total else None,
        "mean_stated_confidence_overall": (
            round(sum(r.stated_confidence for r in scored) / n_total, 4) if n_total else None
        ),
        "bands_populated": sum(1 for r in rows if r["measured_by_this_assay"]),
        "bands_judged": len(judged),
        "min_for_a_band": MIN_FOR_A_BAND,
        "expected_calibration_error": ece,
        "rows": rows,
    }


def verdict(tbl: dict[str, Any], level_interpretable: bool) -> dict[str, Any]:
    """Which of the pre-stated patterns this table matches, and the numbers the pattern rests on.

    Ordering is read first and separately from level, because they fail independently: a scorer that
    ranks correctly and sits five points high is a different object from one that ranks nothing. The
    level is only read where the assay's agreement IS the compiled claim's event; elsewhere it is
    withheld by name rather than computed and quietly ignored.
    """
    judged = [r for r in tbl["rows"] if r["judged"]]
    populated = [r for r in tbl["rows"] if r["measured_by_this_assay"]]
    out: dict[str, Any] = {
        "bands_judged": len(judged),
        "min_bands_for_a_verdict": MIN_BANDS_FOR_A_VERDICT,
        "level_interpretable": level_interpretable,
        # description, never a verdict: the gaps over every populated band, judged or not, so a refusal
        # does not also hide the offset. A refused pattern means no claim is made FROM these numbers
        "described_not_judged": {
            "bands_populated": len(populated),
            "gap_per_populated_band": {r["band"]: r["gap_observed_minus_stated"] for r in populated},
            "median_gap_over_populated_bands": (
                round(median([r["gap_observed_minus_stated"] for r in populated]), 4) if populated else None
            ),
        },
    }
    if len(judged) < MIN_BANDS_FOR_A_VERDICT:
        short = {
            r["band"]: {
                "measured": r["measured_by_this_assay"],
                "short_of_the_bar_by": MIN_FOR_A_BAND - r["measured_by_this_assay"],
            }
            for r in populated
            if not r["judged"]
        }
        out["pattern"] = "refused"
        out["bands_populated_but_below_the_bar"] = short
        out["reading"] = (
            f"{len(judged)} band{' holds' if len(judged) == 1 else 's hold'} at least "
            f"{MIN_FOR_A_BAND} measured elements, below the {MIN_BANDS_FOR_A_VERDICT} this lane fixed "
            "before computing anything. The table is printed and no verdict is drawn from it: a curve "
            "through fewer than three points is a line by construction. The bands that are populated "
            "and below the bar are "
            "named with the size of the shortfall, because a bar missed by one is still a bar and "
            "moving it after seeing the counts is how a pre-registration stops being one"
        )
        return out
    bottom, top = judged[0], judged[-1]
    disjoint = bottom["ci95"][1] < top["ci95"][0]
    monotone = top["observed"] > bottom["observed"]
    rises = monotone and disjoint
    consistent = sum(1 for r in judged if r["stated_inside_the_interval"])
    share = consistent / len(judged)
    gaps = sorted(r["gap_observed_minus_stated"] for r in judged)
    out.update(
        {
            "lowest_judged_band": bottom["band"],
            "lowest_observed": bottom["observed"],
            "highest_judged_band": top["band"],
            "highest_observed": top["observed"],
            "point_estimates_rise": monotone,
            "intervals_disjoint": disjoint,
            "ordering": (
                "rises with the stated confidence"
                if rises
                else "the point estimates rise but the intervals overlap"
                if monotone
                else "does not rise"
            ),
            "bands_where_the_stated_confidence_is_inside_the_interval": consistent,
            "share_consistent": round(share, 4),
            "min_share_consistent": MIN_SHARE_CONSISTENT,
            "median_signed_offset_observed_minus_stated": round(median(gaps), 4),
            "offset_range": [gaps[0], gaps[-1]],
        }
    )
    if not rises:
        out["pattern"] = "uninformative"
    elif not level_interpretable:
        out["pattern"] = "level_withheld_ordering_only"
    elif share >= MIN_SHARE_CONSISTENT:
        out["pattern"] = "calibrated"
    else:
        out["pattern"] = "ordered_but_miscalibrated"
    out["pattern_text"] = PATTERNS[out["pattern"]]
    return out


# --- is the band reading the confidence or the covariate? ----------------------------------------------
def _row(e: Element, assay: str) -> dict[str, Any]:
    return {
        "id": e.id,
        "length": e.length,
        "gc": e.gc,
        "nearest_coding_tss": e.nearest_coding_tss,
        "stated_confidence": e.stated_confidence,
        "agrees": e.verdicts.get(assay) == AGREES,
    }


# the edges come from the lane that declared them rather than a fifth copy of the same tuple
STRATA = Strata(length=(200.0, 400.0), gc=(0.35, 0.45, 0.55), nearest_coding_tss=TSS_EDGES)


def carries_information(population: list[Element], assay: str, narrow: bool) -> dict[str, Any]:
    """Does a high stated confidence still agree more often once the covariates are held fixed?

    `input_presence` runs first and over the whole population, because the question it answers comes
    before any stratum: the stated confidence is on every compiled element with a target, so it is
    FREE, while a verdict from this assay is on 4.33% of them, so it is BOUGHT. A coverage artefact can
    only live in the bought one, and saying which is which costs nothing and stops 'conditioning on
    coverage changed nothing' from reading as a test that was passed when no test was administered.

    The comparison itself is refused below MIN_FOR_A_COMPARISON targets, as the satmut lane refused at
    four: a standardised difference over a handful of elements is noise carrying a p-value.
    """
    scored = [e for e in _labels(population, assay, narrow) if e.covariates_present]
    unmeasured = [e for e in population if assay not in e.verdicts and e.covariates_present]
    presence = input_presence(
        [{**_row(e, assay), "assay_verdict": e.verdicts.get(assay)} for e in scored],
        [{**_row(e, assay), "assay_verdict": None} for e in unmeasured],
        {
            "a stated confidence on the compiled fact": "stated_confidence",
            f"a verdict from {assay} over this element": "assay_verdict",
        },
    )
    high = [_row(e, assay) for e in scored if e.stated_confidence >= HIGH_CONFIDENCE]
    low = [_row(e, assay) for e in scored if e.stated_confidence < HIGH_CONFIDENCE]
    out: dict[str, Any] = {
        "split": (
            f"stated confidence at or above {HIGH_CONFIDENCE} against below it, which is the line "
            "genomeos/evidence.py already draws as WEAK and the Program tab already prints; it was not "
            "chosen here"
        ),
        "denominator_name": (
            "where_the_predicted_gene_was_tested" if narrow else "over_all_matched_elements"
        ),
        "input_presence": presence,
        "targets_high_confidence": len(high),
        "controls_low_confidence": len(low),
        "imbalance": imbalance(high, low, STRATA),
    }
    if len(high) < MIN_FOR_A_COMPARISON or len(low) < MIN_FOR_A_COMPARISON:
        out["comparison_refused"] = (
            f"{len(high)} high-confidence and {len(low)} low-confidence elements, at least one of them "
            f"below {MIN_FOR_A_COMPARISON}: the groups are described above and not compared, because a "
            "standardised difference over this many elements is noise carrying a p-value. The "
            "footprint is the finding, exactly as it was for the base-level assay at four targets"
        )
    else:
        out["standardised"] = standardised(high, low, STRATA, hit="agrees")
    return out


# --- is the measured slice the compiled genome? ---------------------------------------------------------
SCOPE_STRATA = Strata(
    length=(200.0, 400.0),
    gc=(0.35, 0.45, 0.55),
    nearest_coding_tss=TSS_EDGES,
    stated_confidence=(0.15, 0.25, 0.5),
)


def measured_slice(population: list[Element]) -> dict[str, Any]:
    """What the measured elements are, beside what the unmeasured ones are: the size of the transfer.

    No standardised difference is computed here and that is deliberate. The question is not whether
    one arm out-rates the other on some hit; it is how far apart the two arms are on the covariates
    and on the stated confidence itself, which is `imbalance` and nothing more. A difference would
    need a hit both arms have, and the unmeasured arm has no measurement by definition - that is the
    whole problem, and computing a number over it would be the 'conditioned on was-this-measured'
    mistake compare.py's docstring names.
    """
    with_cov = [e for e in population if e.covariates_present]
    seen = [e for e in with_cov if e.verdicts]
    unseen = [e for e in with_cov if not e.verdicts]
    rows_seen = [
        {
            "length": e.length,
            "gc": e.gc,
            "nearest_coding_tss": e.nearest_coding_tss,
            "stated_confidence": e.stated_confidence,
        }
        for e in seen
    ]
    rows_unseen = [
        {
            "length": e.length,
            "gc": e.gc,
            "nearest_coding_tss": e.nearest_coding_tss,
            "stated_confidence": e.stated_confidence,
        }
        for e in unseen
    ]
    band_seen: dict[str, int] = defaultdict(int)
    band_unseen: dict[str, int] = defaultdict(int)
    for e in seen:
        band_seen[e.band] += 1
    for e in unseen:
        band_unseen[e.band] += 1
    return {
        "compiled_elements_with_a_stated_confidence": len(population),
        "elements_without_covariates": len(population) - len(with_cov),
        "measured_by_any_assay": len(seen),
        "never_measured": len(unseen),
        "coverage": round(len(seen) / len(with_cov), 5) if with_cov else None,
        "imbalance": imbalance(rows_seen, rows_unseen, SCOPE_STRATA),
        "band_counts_measured": {b: band_seen.get(b, 0) for b in BAND_NAMES},
        "band_counts_never_measured": {b: band_unseen.get(b, 0) for b in BAND_NAMES},
        "coverage_by_band": {
            b: (
                round(band_seen.get(b, 0) / (band_seen.get(b, 0) + band_unseen.get(b, 0)), 5)
                if (band_seen.get(b, 0) + band_unseen.get(b, 0))
                else None
            )
            for b in BAND_NAMES
        },
        "reading": (
            "the two arms' medians are the size of the extrapolation. A coverage that rises with the "
            "stated confidence is the selection this lane cannot undo: the bands would then be "
            "comparing differently-selected populations, and coverage_by_band is where that shows"
        ),
    }


# --- the report -----------------------------------------------------------------------------------------
def assay_report(population: list[Element], assay: str) -> dict[str, Any]:
    """One assay: the two denominators, both tables, both verdicts, and the covariate check."""
    meaning = WHAT_AGREEMENT_MEANS[assay]
    matched = [e for e in population if assay in e.verdicts]
    narrow_rows = _labels(population, assay, True)
    coincide = len(narrow_rows) == len(matched)
    out: dict[str, Any] = {
        "assay": assay,
        "source": measured.SOURCES[assay],
        "what_agreement_means": meaning,
        "elements_in_the_footprint": sum(1 for e in population if assay in e.eligible),
        "elements_measured": len(matched),
        "elements_where_the_predicted_gene_was_tested": len(narrow_rows),
        "denominators_coincide": coincide,
        "two_denominators": TWO_DENOMINATORS,
    }
    if assay in measured.BASE_LEVEL:
        out["refused"] = (
            f"{len(matched)} matched elements, and the assay cannot disagree: "
            f"{measured.SATMUT_CANNOT_DISAGREE}. A reliability table needs an event that can come out "
            "false, and this one cannot, so no table is computed. The group is described in "
            "`the_measured_slice` with the rest"
        )
        out["agrees"] = sum(1 for e in matched if e.verdicts[assay] == AGREES)
        out["stated_confidences"] = sorted(e.stated_confidence for e in matched)
        out["medians"] = _medians(matched)
        return out
    for narrow in (True, False):
        key = "where_the_predicted_gene_was_tested" if narrow else "over_all_matched_elements"
        if narrow is False and coincide:
            out[key] = {
                "identical_to": "where_the_predicted_gene_was_tested",
                "why": (
                    "every element this assay measured was tested by it, so the wide and the narrow "
                    "denominator are the same rows. The key is kept so the narrow rate is never the "
                    "only rate on the page"
                ),
            }
            continue
        tbl = table(population, assay, narrow)
        level = meaning["level_interpretable"] and narrow
        out[key] = {
            "table": tbl,
            "verdict": verdict(tbl, level),
            "carries_information": carries_information(population, assay, narrow),
        }
        if meaning["level_interpretable"] and not narrow:
            out[key]["level_withheld"] = (
                "the wide denominator counts 'the screen never asked' as a non-agreement, so its rate "
                "is bounded above by the share of elements whose predicted gene the screen chose to "
                "test. A confidence cannot be calibrated against a bound; the ordering can still be "
                "read, and is"
            )
    return out


LEVEL_IS_A_PROPERTY_OF_THE_POPULATION = (
    "the same stated number, measured against two populations, is off in OPPOSITE DIRECTIONS. That is "
    "not an inconsistency in the tables; it is what a level is. target_calibration found the same thing "
    "for a different quantity - the curve transferred and the level did not, because the base rate of "
    "the new population differed - and it is the reason no single offset can be added to the compiler's "
    "formula to fix it. An offset quoted without the population it was measured on is not a number"
)


def the_answer(per_assay: dict[str, Any]) -> dict[str, Any]:
    """The answer assembled from the computed tables, so prose never has to retype a figure.

    Every field here is lifted out of `per_assay`; nothing is recomputed and nothing is typed. The
    offsets are collected across assays because their SIGNS are the finding: a level measured on one
    population does not transfer to another, and that is visible only when both are printed together.
    """
    rows: dict[str, Any] = {}
    for assay, rep in per_assay.items():
        if "refused" in rep:
            rows[assay] = {"table": "refused", "why": rep["refused"]}
            continue
        for key in ("where_the_predicted_gene_was_tested", "over_all_matched_elements"):
            block = rep.get(key) or {}
            if "identical_to" in block:
                continue
            t, v = block["table"], block["verdict"]
            rows[f"{assay}:{key}"] = {
                "elements": t["elements"],
                "observed_overall": t["observed_overall"],
                "mean_stated_confidence_overall": t["mean_stated_confidence_overall"],
                "expected_calibration_error": t["expected_calibration_error"],
                "pattern": v["pattern"],
                "ordering": v.get("ordering"),
                "level_interpretable": v["level_interpretable"],
                "median_gap_over_populated_bands": v["described_not_judged"][
                    "median_gap_over_populated_bands"
                ],
                "gap_per_populated_band": v["described_not_judged"]["gap_per_populated_band"],
                "standardised_difference": (
                    (block["carries_information"].get("standardised") or {}).get("matched", {}) or {}
                ).get("difference"),
                "standardised_p": (
                    (block["carries_information"].get("standardised") or {}).get("matched", {}) or {}
                ).get("p_one_sided"),
            }
    ordered = [k for k, r in rows.items() if r.get("ordering") == "rises with the stated confidence"]
    not_ordered = [k for k, r in rows.items() if r.get("ordering") and k not in ordered]
    populated = inside = 0
    for rep in per_assay.values():
        for key in ("where_the_predicted_gene_was_tested", "over_all_matched_elements"):
            block = rep.get(key) or {}
            if not isinstance(block, dict) or "table" not in block:
                continue
            for r in block["table"]["rows"]:
                if not r["measured_by_this_assay"]:
                    continue
                populated += 1
                inside += bool(r["stated_inside_the_interval"])
    return {
        "per_denominator": rows,
        # the tables that exist, which is not the number of entries: an assay that cannot disagree has
        # an entry saying so and no table, and it must not inflate the denominator of a headline
        "tables": sum(1 for r in rows.values() if "elements" in r),
        "entries_with_no_table_at_all": [k for k, r in rows.items() if "elements" not in r],
        "bands_populated_across_every_table": populated,
        "bands_where_the_stated_level_is_inside_its_own_interval": inside,
        "bands_where_it_is_not": populated - inside,
        "tables_where_the_rate_rises_with_the_confidence": ordered,
        "tables_where_it_does_not": not_ordered,
        "tables_refused": [k for k, r in rows.items() if r.get("pattern") in ("refused", None)],
        "level_is_a_property_of_the_population": LEVEL_IS_A_PROPERTY_OF_THE_POPULATION,
        "offsets_by_table": {
            k: r["median_gap_over_populated_bands"]
            for k, r in rows.items()
            if r.get("median_gap_over_populated_bands") is not None
        },
    }


def report(population: list[Element], per_chromosome: dict[str, Any]) -> dict[str, Any]:
    """The whole answer: what is being tested, what would count as a failure, and what happened."""
    per_assay = {a: assay_report(population, a) for a in ASSAYS}
    return {
        "question": (
            "GenomeOS states a confidence on every fact. Is it a probability? For the compiled facts a "
            "measurement touches, the observed agreement rate is binned by the stated confidence and "
            "reported per assay, with an interval"
        ),
        "claim_under_test": CLAIM_UNDER_TEST,
        "stated_confidence_formula": STATED_CONFIDENCE_FORMULA,
        "bands": {
            "edges": [list(b) for b in BANDS],
            "half_open": "[lo, hi)",
            "where_they_come_from": (
                "imported unchanged from genomeos/attribution/target_calibration.CONFIDENCE_BANDS, "
                "fixed there on 2026-09-17 for a different quantity and before this lane existed. They "
                "were not chosen after seeing a rate, and they are not adjusted by anything here"
            ),
            "min_for_a_band": MIN_FOR_A_BAND,
            "min_bands_for_a_verdict": MIN_BANDS_FOR_A_VERDICT,
            "min_share_consistent": MIN_SHARE_CONSISTENT,
        },
        "pre_stated_patterns": PATTERNS,
        "why_pooling_is_refused": WHY_POOLING_IS_REFUSED,
        "two_denominators": TWO_DENOMINATORS,
        "per_chromosome": per_chromosome,
        "the_measured_slice": measured_slice(population),
        "per_assay": per_assay,
        "the_answer": the_answer(per_assay),
        "no_pooled_rate": (
            "there is deliberately no pooled agreement rate in this result. See why_pooling_is_refused"
        ),
        "scope": SCOPE,
        "falsifies_transfer": FALSIFIES_TRANSFER,
    }


def summary_line(got: dict[str, Any], assay: str) -> str:
    """One line per assay for the console, eligible before measured and both denominators named."""
    a = got["per_assay"][assay]
    if "refused" in a:
        return f"{assay}: {a['elements_measured']} measured, no table ({a['refused'][:60]}...)"
    bits = []
    for key in ("where_the_predicted_gene_was_tested", "over_all_matched_elements"):
        block = a.get(key) or {}
        if "identical_to" in block:
            bits.append(f"{key}: identical to the narrow denominator")
            continue
        t, v = block["table"], block["verdict"]
        bits.append(
            f"{key}: n={t['elements']}, observed {t['observed_overall']}, stated "
            f"{t['mean_stated_confidence_overall']}, {t['bands_judged']} bands judged, "
            f"{v['pattern']}"
        )
    return (
        f"{assay}: {a['elements_in_the_footprint']} in the footprint, {a['elements_measured']} "
        f"measured; " + "; ".join(bits)
    )
