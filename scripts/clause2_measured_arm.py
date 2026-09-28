# SPDX-License-Identifier: AGPL-3.0-or-later
"""Milestone 1.3's clause 2 asked of measurement rather than of the model: what an assay says about the
elements inside a constrained-unknown block.

    uv run python scripts/clause2_measured_arm.py [--chroms chr21] [--no-save]

Four lanes have now put clause 2's second part against four controls -- a length-matched draw
(`d717b28`), a length- and coding-TSS-density-matched draw (`e1dbcf3`/`0c8b82d`), a length- and
element-count-matched draw (`da5764e`/`faeb0da`) and a control-free per-element estimator -- and the
real unknown stays between 18 and 29 points below every one of them. Each lane ended with the same
sentence: beyond the descriptive controls, only measurement separates "these blocks hold no element
that regulates a coding gene" from "the deletion model cannot name a target for sequence like this".

This is that arm, and it is not the one `clause2_matched_control.PRE_REGISTRATION` left registered.
That entry (`second_route_not_run`) registered a *designed* experiment -- 284,001 unsynthesised oligos
over 878 blocks -- and correctly said nothing measured exists for it. The arm run here is the one the
project can actually run today: the measurements it already holds, asked of the same elements, over
the same blocks, against the same windows, with the same statistic and the same interval method.

What is fixed by `clause2_matched_control` and reused unchanged: the target sets, the window draw,
the coding-TSS decile matching, the seed, the draw count, the per-block matched difference, the two
bootstraps and the reproduction gate. What this script registers, because that one leaves it open:
which measured outcome counts as "regulates a coding gene", which denominator the measured arm uses,
the power floor below which no interval is read, and -- fixed before the first count -- how each of
the three possible outcomes is to be read, including the likely one, that measurement is absent.

Model requests: none. It reads the stored all-element archive one chromosome at a time and the four
assay tables through `genomeos.attribution.measured`; the per-element response cache is never opened.
"""

from __future__ import annotations

import argparse
import importlib.util
import math
import time
from pathlib import Path
from typing import Any

from genomeos import manifest as mf
from genomeos.attribution import measured as ms
from genomeos.attribution import mpra, organise, vista
from genomeos.results import save_result

_spec = importlib.util.spec_from_file_location(
    "clause2_matched_control", Path(__file__).resolve().parent / "clause2_matched_control.py"
)
mc = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(mc)
cut = mc.cut  # scripts/constrained_unknown_targets.py, through the matched-control module

RESULT = "clause2_measured_arm"

# --- everything below is `clause2_matched_control`'s, imported and not restated ---------------------
DRAWS = mc.DRAWS
SEED = mc.SEED
UNMATCHED_TRIES = mc.UNMATCHED_TRIES
BOOTSTRAP = mc.BOOTSTRAP

# --- what this lane fixes ---------------------------------------------------------------------------
#: the model question the gate reproduces, asked on the same windows at no extra draw
MODEL_QUESTION = "names_a_coding_gene"
#: the measured question the lane turns on: an assay measured this element to move a protein-coding gene
PRIMARY_QUESTION = "measured_to_regulate_a_coding_gene"
#: the denominator of the primary: an assay tested this element against at least one coding gene and
#: could therefore have said yes. Without it a never-measured element would count as a measured "no".
DENOMINATOR_QUESTION = "crispri_tested_against_a_coding_gene"
#: `measured.MIN_FOR_A_COMPARISON`, adopted unchanged: below this many blocks no interval is read
MIN_BLOCKS = ms.MIN_FOR_A_COMPARISON
#: the effect the measured arm would have to detect to overturn the model arm, in proportion, from
#: `0c8b82d`'s committed primary (-27.25 points)
EFFECT_TO_DETECT = 0.2725
#: (z_{0.975} + z_{0.80})^2, the constant of the one-sample power formula used for "coverage needed"
POWER_CONSTANT = (1.959964 + 0.841621) ** 2
#: `0c8b82d`'s committed matched primary, which this run must reproduce on its own windows to the digit
REPRODUCE = {
    "matched_difference_points": -27.25,
    "ci95_over_blocks": [-30.91, -23.58],
    "n_blocks_compared": 531,
    "blocks_carrying_an_element": 531,
    "windows_drawn": 44081,
    "pooled": {"blocks_yes": 158, "windows_yes": 12876, "windows_carrying": 22228},
}

PRE_REGISTRATION = {
    "registered": "2026-09-28",
    "registration_reused": (
        "`scripts/clause2_matched_control.PRE_REGISTRATION` (e1dbcf3, and its dated section in "
        "docs/ATTRIBUTION.md) settles the target sets, the control draw, the covariate and its bins, "
        "the per-block estimator, the interval method, the reproduction gate and the pass rule for the "
        "model arm. None of it is redesigned here: this script imports those functions rather than "
        "restating them, and the model question is carried through the same draw so the gate is free. "
        "Its `second_route_not_run` entry registered a DESIGNED measured arm (284,001 unsynthesised "
        "oligos) and said nothing measured exists for it; that is still true and that experiment is "
        "still unrun. This lane runs the different arm it left open: the measurements the project "
        "already holds, over the same blocks and the same windows"
    ),
    "targets": (
        "unchanged from e1dbcf3: the real unknown, constrained_unknown blocks with the copies out "
        "(882 blocks, 531 carrying a scored element) from organise.blocks on all 24 chromosomes; the "
        "neutral tier as the secondary target set"
    ),
    "control": (
        "unchanged from e1dbcf3 and drawn by the same function with the same seed: windows of each "
        "block's exact length, placed uniformly inside the span of the chromosome's scored-element "
        "midpoints, rejected on overlap with any organiser block and rejected unless the count of "
        "GENCODE protein-coding TSSs within 524,288 bp of the window's midpoint falls in the block's "
        "own decile; 50 accepted windows per block, 20,000 tries. The unmatched draw (d717b28's, 4,000 "
        "tries, no density test) is reported beside it"
    ),
    "measurement": (
        "`genomeos.attribution.measured` and nothing else, read but never edited by this lane. A "
        "measurement is OF an element only at RECIPROCAL_OVERLAP = 0.5 both ways, the module's own "
        "rule. Four assays, each under its own name: CRISPRi (the ENCODE enhancer-gene benchmark, both "
        "files, R2's five outcomes kept apart and R5's split field carried); lentiMPRA (ENCODE4 "
        "ENCSR106SZM, R6's per-tile reading -- the per-cell label is the median of the matched tiles "
        "under REPORTER_LABEL_RULE, never the strongest tile, and a tile-conflicting cell is counted "
        "under its own name); VISTA (transgenic mouse e11.5, positive and negative alike); saturation "
        "mutagenesis (GSE126550), which is base-level and is counted and reported apart, never pooled "
        "into an element-level verdict"
    ),
    "questions": {
        MODEL_QUESTION: (
            "the model arm, carried on the same windows so the gate costs no extra draw: does the "
            "block (or window) hold an element whose deletion names a protein-coding gene "
            "(`predicted_coding`)"
        ),
        "eligible_for_an_assay": (
            "`Layer.touched_by`: some assay's tested interval shares at least one base with the "
            "element. The weakest bar, and the one that says where the screens were pointed"
        ),
        "measured_by_an_element_level_assay": (
            "CRISPRi, lentiMPRA or VISTA measures this element at reciprocal overlap 0.5"
        ),
        "measured_by_any_assay": "the three above or saturation mutagenesis",
        DENOMINATOR_QUESTION: (
            "CRISPRi tested this element against at least one gene whose GENCODE symbol is "
            "protein-coding, in any cell, in either split. This and only this is the denominator of "
            "the primary: an element no screen ever pointed at is not a measured negative"
        ),
        PRIMARY_QUESTION: (
            "of the elements CRISPRi tested against a coding gene, at least one coding gene moved "
            "significantly on silencing -- a significant DECREASE (the benchmark's own `Regulated`) or "
            "a significant INCREASE. Both signs count, because R2 established that the benchmark's "
            "`Regulated` is a decrease only and that a significant increase is a measured regulatory "
            "effect and not a null. The decrease-only form is reported beside it"
        ),
        "measured_regulated_decrease_only": "the benchmark's own `Regulated`, restricted to coding genes",
        "measured_null_well_powered_on_every_coding_gene": (
            "CRISPRi tested the element against at least one coding gene, no coding gene moved "
            "significantly either way, and at least one coding gene's null carries power >= 0.8 at a "
            "20% effect. This is the only measured outcome that can say 'this element regulates no "
            "coding gene'; an underpowered null cannot and is counted apart"
        ),
        "measured_active_but_no_gene_named": (
            "lentiMPRA calls the element active in at least one cell under REPORTER_LABEL_RULE, or "
            "VISTA calls it a positive enhancer. Registered as a SEPARATE question and never folded "
            "into the primary: both assays measure that the sequence drives transcription and neither "
            "names the gene it drives, so neither can answer a clause that asks for a gene"
        ),
    },
    "primary": (
        "the real unknown against its density- and length-matched windows on "
        f"`{PRIMARY_QUESTION}`, by e1dbcf3's estimator with e1dbcf3's denominator replaced by the "
        "measured one: the mean over blocks that hold at least one element CRISPRi tested against a "
        "coding gene AND have at least one accepted window holding one, of (block yes, 0 or 1) minus "
        "(share of that block's TESTED windows saying yes), in points"
    ),
    "interval": (
        "unchanged from e1dbcf3: 95% percentile bootstrap over blocks, 10,000 resamples, numpy "
        "default_rng(20260913), with a bootstrap over the chromosomes reported beside it and the "
        "pooled rates printed as d717b28 printed them"
    ),
    "power_floor": (
        f"no interval is read below {MIN_BLOCKS} compared blocks. That is "
        "`measured.MIN_FOR_A_COMPARISON`, the bar the measured layer already sets for a comparison of "
        "its own, adopted here unchanged rather than chosen"
    ),
    "gate": (
        "before any measured figure is read, the model question carried through the same matched draw "
        "must reproduce 0c8b82d's committed primary to the digit: -27.25 points, 95% -30.91 to -23.58 "
        "over blocks, 531 compared blocks, 44,081 windows drawn, 158 blocks yes, 12,876 window-yes of "
        "22,228 carrying windows. If it does not, these are not the same windows and no measured "
        "figure is reported"
    ),
    "secondary": [
        "the same statistic on `measured_regulated_decrease_only`",
        "the neutral tier under its own deciles, same statistic and interval",
        "the same-instrument estimator: e1dbcf3's denominator left alone, every block carrying a scored "
        "element in, an unmeasured element counting as a no. Reported as DESCRIPTIVE ONLY and never as "
        "the primary, because a rate that counts absence of measurement as a negative is the error this "
        "project has had to withdraw three times",
        "the coverage table at reciprocal overlap 0.25 and 0.75 beside 0.5 "
        "(`measured.OVERLAP_SENSITIVITY`), on the counts only, so the rule can be seen rather than "
        "trusted",
    ],
    "coverage_is_the_first_result": (
        "reported before any difference and regardless of every outcome: of the 882 real-unknown blocks "
        "and the 531 that carry a scored element, how many carry an element that any assay is eligible "
        "to have measured, how many carry one it did measure, how many carry one CRISPRi tested against "
        "a coding gene, and how many carry one measured to move it. The same four counts for the matched "
        "windows and for the neutral tier"
    ),
    "readings": {
        "model_failed": (
            "the blocks DO hold elements measured to regulate a coding gene, at a rate not below their "
            "matched windows: the primary's 95% interval over blocks reaches or lies above 0 with at "
            f"least {MIN_BLOCKS} compared blocks. Then the failure is the MODEL'S. The sequence "
            "regulates coding genes and the deletion model could not name the target, so clause 2's "
            "second part failed on the instrument and not on the sequence; the clause stays not met, "
            "and its stated reason changes from 'these blocks hold nothing that regulates a coding "
            "gene' to 'the deletion model cannot name a target for sequence like this'. The four "
            "descriptive controls keep their figures and lose their reading"
        ),
        "wording_wrong": (
            "the blocks hold NO element measured to regulate a coding gene while their matched windows "
            f"do: the primary's 95% interval lies wholly below 0 with at least {MIN_BLOCKS} compared "
            "blocks. Then the model was right about these blocks. What is wrong is clause 2's WORDING: "
            "a milestone cannot ask for a gene and a tissue from sequence that regulates no gene, and "
            "the clause is to be rewritten rather than re-tested. The blocks' constraint then wants an "
            "account that is not regulation of a coding gene"
        ),
        "cannot_decide": (
            f"fewer than {MIN_BLOCKS} real-unknown blocks hold an element CRISPRi tested against a "
            "coding gene, or fewer than that many also have a matched window holding one. This is the "
            "likely case and it is registered here as a FINISHED LANE and a real result, not a "
            "failure: no interval is reported, both readings above stay open, and the lane's output is "
            "the coverage that exists and the coverage that would be needed -- the real-unknown blocks "
            "with no measured element at all, the blocks short of the floor, and the number of blocks "
            "an assay would have to reach for 80% power against an effect the size of the model arm's "
            "own -27.25 points"
        ),
    },
    "coverage_needed_formula": (
        "n_for_80_percent_power = ceil((z_0.975 + z_0.80)^2 * s^2 / d^2) with d = 0.2725, the model "
        "arm's committed matched difference, and s the standard deviation of THIS RUN'S per-block "
        "model-arm differences on the same target set, which the same draw produces. Fixed here before "
        "any count so the shortfall is not chosen after seeing it"
    ),
    "what_this_arm_cannot_do": (
        "CRISPRi is the only one of the four assays that names a gene, so it is the only one that can "
        "answer this clause, and it was pointed at candidate regulatory sequence near expressed genes "
        "in a handful of cell lines -- not at constrained unknown sequence. A block with no CRISPRi "
        "element is not evidence either way. lentiMPRA and VISTA can say the sequence does something "
        "and cannot say what to; saturation mutagenesis can say which bases matter inside an element "
        "and, by SATMUT_CANNOT_DISAGREE, can never say an element does nothing. None of that is a "
        "defect of this run; it is the size of the hole the lane is measuring"
    ),
    "cost": (
        "0 AlphaGenome requests. The per-element response cache is never opened; the all-element "
        "archive and one chromosome's measured layer are held at a time"
    ),
}


# ---- the measured verdict of one element -----------------------------------------------------------


def coding_symbols(chroms: list[str]) -> set[str]:
    """Every GENCODE protein-coding gene symbol, from the same annotation `coding_tss` reads."""
    from genomeos.genome import Annotation
    from genomeos.genome.annotation import default_gencode

    out: set[str] = set()
    for c in chroms:
        ann = Annotation.from_gff3(default_gencode({c}), {c})
        out |= {g.symbol for g in ann.protein_coding() if g.symbol}
    return out


def verdict(
    layer: ms.Layer, start: int, end: int, coding: set[str], fraction: float = ms.RECIPROCAL_OVERLAP
) -> dict[str, Any]:
    """What measurement says about one element, each assay under its own name and nothing pooled."""
    m = layer.for_element(start, end, fraction=fraction)
    cr, lm, vs, sm = (m.get(k) for k in ("crispri", "lentimpra", "vista", "satmut"))
    tested = set(cr["genes_tested"]) & coding if cr else set()
    decrease = set(cr["genes_regulated"]) & coding if cr else set()
    increase = set(cr["genes_increased"]) & coding if cr else set()
    powered = set(cr["genes_no_effect_well_powered"]) & coding if cr else set()
    weak = set(cr["genes_no_effect_underpowered"]) & coding if cr else set()
    active = bool(lm and lm["cells_active"])
    positive = bool(vs and vs["positive"])
    return {
        "eligible_for_an_assay": bool(layer.touched_by(start, end)),
        "measured_by_an_element_level_assay": bool(cr or lm or vs),
        "measured_by_any_assay": bool(cr or lm or vs or sm),
        "crispri": bool(cr),
        "lentimpra": bool(lm),
        "vista": bool(vs),
        "satmut": bool(sm),
        DENOMINATOR_QUESTION: bool(tested),
        PRIMARY_QUESTION: bool(decrease or increase),
        "measured_regulated_decrease_only": bool(decrease),
        "measured_increase_only": bool(increase and not decrease),
        "measured_null_well_powered_on_every_coding_gene": bool(
            tested and not (decrease or increase) and powered
        ),
        "measured_null_underpowered_only": bool(tested and not (decrease or increase or powered) and weak),
        "measured_active_but_no_gene_named": bool(active or positive),
        "lentimpra_active": active,
        ms.TILES_CONFLICT: bool(lm and lm["cells_conflicting"]),
        "vista_positive": positive,
        "vista_negative": bool(vs and vs["negative"] and not positive),
        "crispri_tested_against_any_gene": bool(cr),
        "crispri_coding_genes_tested": len(tested),
        "crispri_coding_genes_moved": len(decrease | increase),
    }


#: the questions carried through the window draw, in the order the result prints them
QUESTIONS = (
    MODEL_QUESTION,
    "eligible_for_an_assay",
    "measured_by_an_element_level_assay",
    "measured_by_any_assay",
    DENOMINATOR_QUESTION,
    PRIMARY_QUESTION,
    "measured_regulated_decrease_only",
    "measured_null_well_powered_on_every_coding_gene",
    "measured_active_but_no_gene_named",
)
#: the flags counted per element in the coverage table but not carried through the draw
COVERAGE_FLAGS = (
    "crispri",
    "lentimpra",
    "vista",
    "satmut",
    "measured_increase_only",
    "measured_null_underpowered_only",
    "lentimpra_active",
    ms.TILES_CONFLICT,
    "vista_positive",
    "vista_negative",
)


def predicates() -> dict[str, Any]:
    """The questions as predicates over an element dict whose verdict is already attached."""
    out: dict[str, Any] = {MODEL_QUESTION: cut.moves_coding}
    for q in QUESTIONS:
        if q != MODEL_QUESTION:
            out[q] = lambda e, q=q: bool(e["_measured"][q])
    return out


def attach(
    els: list[dict[str, Any]], layer: ms.Layer, coding: set[str], fraction: float = ms.RECIPROCAL_OVERLAP
) -> None:
    """Attach each element's measured verdict once, so a window costs a dict lookup and not a scan."""
    for e in els:
        e["_measured"] = verdict(layer, e["start"], e["end"], coding, fraction=fraction)


# ---- coverage ---------------------------------------------------------------------------------------

COVERAGE_TOTALS = (
    "blocks",
    "blocks_carrying_a_scored_element",
    "scored_elements",
    "blocks_with_no_measured_element",
    "crispri_coding_pairs_tested",
)


def coverage(rows: list[dict[str, Any]], els_by_block: dict[str, list[dict[str, Any]]]) -> dict[str, Any]:
    """The lane's first result: how much of this block set measurement reaches at all."""
    flags = list(QUESTIONS[1:]) + list(COVERAGE_FLAGS)
    blocks = dict.fromkeys(flags, 0)
    elements = dict.fromkeys(flags, 0)
    n_elements = untouched = pairs = 0
    for r in rows:
        got = els_by_block.get(r["block"], [])
        n_elements += len(got)
        untouched += not any(e["_measured"]["measured_by_any_assay"] for e in got)
        pairs += sum(e["_measured"]["crispri_coding_genes_tested"] for e in got)
        for f in flags:
            hits = sum(1 for e in got if e["_measured"][f])
            elements[f] += hits
            blocks[f] += hits > 0
    return {
        "blocks": len(rows),
        "blocks_carrying_a_scored_element": sum(1 for r in rows if r["elements"]),
        "scored_elements": n_elements,
        "blocks_with_no_measured_element": untouched,
        "crispri_coding_pairs_tested": pairs,
        "blocks_with": blocks,
        "elements_with": elements,
    }


def pool_coverage(per_chrom: list[dict[str, Any]]) -> dict[str, Any]:
    if not per_chrom:
        return {}
    out: dict[str, Any] = {k: sum(c[k] for c in per_chrom) for k in COVERAGE_TOTALS}
    for key in ("blocks_with", "elements_with"):
        out[key] = {f: sum(c[key][f] for c in per_chrom) for f in per_chrom[0][key]}
    return out


# ---- the statistic, e1dbcf3's with the measured denominator ------------------------------------------


def restricted(records: list[dict[str, Any]], denom: str = DENOMINATOR_QUESTION) -> list[dict[str, Any]]:
    """Blocks the denominator question reaches, that also have a window it reaches."""
    return [r for r in records if r["yes"][denom] and r["windows_yes"][denom]]


def summarise_measured(
    by_chrom: dict[str, list[dict[str, Any]]],
    q: str,
    denom: str = DENOMINATOR_QUESTION,
    seed: int = SEED,
    n_boot: int = BOOTSTRAP,
) -> dict[str, Any]:
    """`clause2_matched_control.summarise` step for step, with the measured denominator in place of
    "carries a scored element": a block's window rate is over the windows the assay could have
    measured, not over every window carrying an element."""
    import numpy as np

    def diffs(recs: list[dict[str, Any]]) -> list[float]:
        return [
            float(r["yes"][q]) - r["windows_yes"][q] / r["windows_yes"][denom]
            for r in restricted(recs, denom)
        ]

    chroms = [c for c in by_chrom if restricted(by_chrom[c], denom)]
    per = {c: np.array(diffs(by_chrom[c])) for c in chroms}
    allrec = [r for c in by_chrom for r in by_chrom[c]]
    comp = restricted(allrec, denom)
    d = np.concatenate([per[c] for c in chroms]) if chroms else np.array([])
    out: dict[str, Any] = {
        "question": q,
        "denominator": denom,
        "target_blocks": len(allrec),
        "blocks_the_denominator_reaches": sum(1 for r in allrec if r["yes"][denom]),
        "blocks_reached_without_a_reached_window": sum(
            1 for r in allrec if r["yes"][denom] and not r["windows_yes"][denom]
        ),
        "n_blocks_compared": len(comp),
        "windows_the_denominator_reaches": sum(r["windows_yes"][denom] for r in comp),
        "power_floor": MIN_BLOCKS,
        "below_the_power_floor": len(comp) < MIN_BLOCKS,
    }
    if not len(d) or len(comp) < MIN_BLOCKS:
        out["matched_difference_points"] = None
        out["why"] = (
            f"{len(comp)} compared blocks is below the registered floor of {MIN_BLOCKS}: no interval is read"
            if len(d)
            else "no block of this set holds an element the denominator question reaches"
        )
        return out
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(d), size=(n_boot, len(d)))
    boots = d[idx].mean(axis=1)
    sums = np.array([per[c].sum() for c in chroms])
    ns = np.array([len(per[c]) for c in chroms])
    cidx = rng.integers(0, len(chroms), size=(n_boot, len(chroms)))
    cboots = sums[cidx].sum(axis=1) / ns[cidx].sum(axis=1)
    by = sum(1 for r in comp if r["yes"][q])
    wy = sum(r["windows_yes"][q] for r in comp)
    wd = sum(r["windows_yes"][denom] for r in comp)
    out.update(
        matched_difference_points=round(100 * float(d.mean()), 2),
        ci95_over_blocks=[round(100 * float(x), 2) for x in np.percentile(boots, [2.5, 97.5])],
        ci95_over_chromosomes=[round(100 * float(x), 2) for x in np.percentile(cboots, [2.5, 97.5])],
        block_rate=round(by / len(comp), 4),
        pooled={
            "blocks_yes": by,
            "windows_yes": wy,
            "windows_reached": wd,
            "block_rate": round(by / len(comp), 4),
            "window_rate": round(wy / wd, 4),
            "difference_in_points": round(100 * (by / len(comp) - wy / wd), 2),
        },
    )
    return out


def coverage_needed(
    records: list[dict[str, Any]], compared_now: int, model_q: str = MODEL_QUESTION
) -> dict[str, Any]:
    """How many blocks an assay would have to reach, by the formula registered above."""
    import numpy as np

    d = np.array(mc.differences(records, model_q))
    s = float(d.std(ddof=1)) if len(d) > 1 else None
    need = math.ceil(POWER_CONSTANT * s**2 / EFFECT_TO_DETECT**2) if s else None
    return {
        "blocks_compared_now": compared_now,
        "power_floor": MIN_BLOCKS,
        "blocks_short_of_the_floor": max(0, MIN_BLOCKS - compared_now),
        "model_arm_sd_of_per_block_differences": round(s, 4) if s is not None else None,
        "effect_to_detect": EFFECT_TO_DETECT,
        "n_for_80_percent_power": need,
        "blocks_short_of_80_percent_power": max(0, need - compared_now) if need else None,
        "share_of_the_882_real_unknown_blocks_needed": round(need / 882, 4) if need else None,
    }


#: Where the power question is actually answered, after the statistical review of 2026-09-28 found
#: that `coverage_needed` above could not answer it. `coverage_needed` is kept exactly as it was,
#: because it is the record of what this run computed and what a7f207f reported from it.
POWER_LANE = "scripts/clause2_design_power.py"


def reporting_floor_and_power_assumptions(
    records: list[dict[str, Any]], compared_now: int, model_q: str = MODEL_QUESTION
) -> dict[str, Any]:
    """The same three numbers `coverage_needed` returns, separated into what each one is.

    Added 2026-09-28 (lane-design) after a statistical review of milestone 1.3's notes. Nothing above
    is changed or removed: `coverage_needed` stays as the record of what was computed and published.
    This function exists because that record's field names let three different things be read as one:

    * a **reporting floor**, `measured.MIN_FOR_A_COMPARISON`, which is a rule about when a number may
      be printed and carries no statement about power at all. The "19 blocks short" that was quoted
      is this and nothing else: 20 minus the 1 block compared;
    * two **assumptions borrowed from the model arm** -- the effect to detect and the dispersion --
      used to size an experiment whose endpoint is a measurement. Both are model output. Nothing
      licenses carrying them across, and they are labelled here as assumptions wherever they appear;
    * a **provisional calculation** from those assumptions, which is not a power result for the
      measured experiment and must not be quoted as the size of it.
    """
    base = coverage_needed(records, compared_now, model_q)
    return {
        "reporting_floor": {
            "floor_blocks": MIN_BLOCKS,
            "blocks_compared_now": compared_now,
            "blocks_short_of_the_floor": base["blocks_short_of_the_floor"],
            "what_it_is": (
                "`measured.MIN_FOR_A_COMPARISON`: below this many compared blocks no interval is "
                "printed. A rule about reporting"
            ),
            "what_it_is_not": (
                "a power calculation, an effect size, a sample size, or a statement that this many "
                "blocks would decide anything"
            ),
        },
        "model_derived_assumptions": {
            "effect_to_detect": {
                "value": EFFECT_TO_DETECT,
                "is": "AN ASSUMPTION, not an estimate for this endpoint",
                "source": (
                    "the model arm's committed matched difference (0c8b82d): a difference in how "
                    "often a deletion model names a target. It is not an estimate of a difference in "
                    "measured regulation"
                ),
            },
            "dispersion": {
                "value": base["model_arm_sd_of_per_block_differences"],
                "is": "AN ASSUMPTION, not an estimate for this endpoint",
                "source": "the standard deviation of this run's per-block MODEL-arm differences",
            },
            "endpoint_they_were_used_for": PRIMARY_QUESTION,
            "why_that_is_a_mismatch": (
                "both inputs describe the model's behaviour; the experiment they size has a measured "
                "outcome. A model's effect size and a model's dispersion are not transferable to it"
            ),
        },
        "provisional_calculation": {
            "n_if_those_assumptions_held": base["n_for_80_percent_power"],
            "is_a_power_result_for_the_measured_experiment": False,
            "why_not": (
                "its two inputs are model output, and 80% power is in any case the probability of "
                "detecting an effect of an assumed size, never a guarantee that an experiment decides "
                "the clause"
            ),
            "coincidence_to_note": (
                "this returned 20, which is also the reporting floor, so the floor's '19 short' and "
                "the calculation's '19 short' looked like one number confirming another. They are "
                "unrelated quantities that happened to agree"
            ),
        },
        "where_the_power_question_is_answered": POWER_LANE,
        "record_of_what_was_computed": base,
    }


# ---- the run ------------------------------------------------------------------------------------


def window_coverage(by_chrom: dict[str, list[dict[str, Any]]]) -> dict[str, Any]:
    recs = [r for c in by_chrom for r in by_chrom[c]]
    return {
        "windows_drawn": sum(r["drawn"] for r in recs),
        "windows_carrying_a_scored_element": sum(r["carrying"] for r in recs),
        "windows_with": {q: sum(r["windows_yes"][q] for r in recs) for q in QUESTIONS},
        "blocks_with_at_least_one_such_window": {
            q: sum(1 for r in recs if r["windows_yes"][q]) for q in QUESTIONS
        },
    }


def measured_reading(primary: dict[str, Any]) -> dict[str, Any]:
    """Which of the three registered readings this run supports, by the registered rules alone."""
    if primary.get("matched_difference_points") is None:
        key = "cannot_decide"
    else:
        key = "wording_wrong" if primary["ci95_over_blocks"][1] < 0 else "model_failed"
    return {
        "outcome": key,
        "decides_clause_2": key != "cannot_decide",
        "text": PRE_REGISTRATION["readings"][key],
    }


def drop_boots(s: dict[str, Any]) -> dict[str, Any]:
    s.pop("_boots", None)
    return s


def collect(chroms: list[str]) -> dict[str, Any]:
    t0 = time.time()
    blocks: dict[str, list[dict[str, Any]]] = {}
    tss: dict[str, list[int]] = {}
    for c in chroms:
        blocks[c] = organise.blocks(c)
        tss[c] = mc.coding_tss(c)
    coding = coding_symbols(chroms)
    print(f"{len(coding)} protein-coding symbols over {len(chroms)} chromosomes", flush=True)
    edges = {
        name: mc.decile_edges(
            [
                mc.tss_count(tss[c], (b["start"] + b["end"]) // 2)
                for c in chroms
                for b in mc.target_sets(blocks[c])[name]
            ]
        )
        for name in ("real_unknown", "neutral")
    }
    preds = predicates()
    runs: dict[str, dict[str, list[dict[str, Any]]]] = {
        k: {}
        for k in ("unmatched_real_unknown", "matched_real_unknown", "unmatched_neutral", "matched_neutral")
    }
    cov: dict[str, list[dict[str, Any]]] = {"real_unknown": [], "neutral": []}
    other = [f for f in ms.OVERLAP_SENSITIVITY if f != ms.RECIPROCAL_OVERLAP]
    cov_sens: dict[float, dict[str, list[dict[str, Any]]]] = {
        f: {"real_unknown": [], "neutral": []} for f in other
    }
    layers: dict[str, dict[str, int]] = {}
    for c in chroms:
        els = cut.elements_of(c)
        if not els:
            continue
        els.sort(key=lambda e: e["start"])
        got = cut.read_chromosome(c, els)
        if not got:
            continue
        layer = ms.Layer.load(c)
        layers[c] = layer.counts()
        attach(els, layer, coding)
        by_block = {r["block"]: cut.inside(r, els) for r in got}
        sets = mc.target_sets(got)
        t = tss[c]
        for name in ("real_unknown", "neutral"):
            cov[name].append(coverage(sets[name], by_block))
            runs[f"unmatched_{name}"][c] = mc.matched_windows(
                sets[name], got, els, preds, None, max_tries=UNMATCHED_TRIES
            )
            e = edges[name]
            runs[f"matched_{name}"][c] = mc.matched_windows(
                sets[name],
                got,
                els,
                preds,
                lambda m, e=e, t=t: mc.bin_of(mc.tss_count(t, m), e),
                raw_fn=lambda m, t=t: mc.tss_count(t, m),
            )
        for f in other:  # the overlap rule reported at three values, on the counts only
            attach(els, layer, coding, fraction=f)
            for name in ("real_unknown", "neutral"):
                cov_sens[f][name].append(coverage(sets[name], by_block))
        del els, got, by_block, layer
        print(f"{c}: done ({time.time() - t0:.0f} s)", flush=True)

    model = drop_boots(mc.summarise(runs["matched_real_unknown"], MODEL_QUESTION))
    gate = {
        "want": REPRODUCE,
        "got": {
            "matched_difference_points": model["matched_difference_points"],
            "ci95_over_blocks": model["ci95_over_blocks"],
            "n_blocks_compared": model["n_blocks_compared"],
            "blocks_carrying_an_element": model["blocks_carrying_an_element"],
            "windows_drawn": model["windows_drawn"],
            "pooled": {k: model["pooled"][k] for k in ("blocks_yes", "windows_yes", "windows_carrying")},
        },
    }
    gate["passed"] = gate["got"] == gate["want"]
    out: dict[str, Any] = {
        "result": RESULT,
        "chromosomes": chroms,
        "registration": PRE_REGISTRATION,
        "measured_layer_counts": layers,
        "coverage": {
            "reciprocal_overlap": ms.RECIPROCAL_OVERLAP,
            "real_unknown": pool_coverage(cov["real_unknown"]),
            "neutral": pool_coverage(cov["neutral"]),
            "windows_matched_real_unknown": window_coverage(runs["matched_real_unknown"]),
            "windows_unmatched_real_unknown": window_coverage(runs["unmatched_real_unknown"]),
            "sensitivity_by_overlap": {
                str(f): {n: pool_coverage(v[n]) for n in v} for f, v in cov_sens.items()
            },
        },
        "gate_reproduces_0c8b82d": gate,
        "model_arm_on_these_windows": model,
        "seconds": round(time.time() - t0, 1),
    }
    if not gate["passed"]:
        out["reading"] = (
            "the reproduction gate failed: these are not 0c8b82d's windows, no measured figure is reported"
        )
        return out

    allrec = [r for c in runs["matched_real_unknown"] for r in runs["matched_real_unknown"][c]]
    primary = summarise_measured(runs["matched_real_unknown"], PRIMARY_QUESTION)
    out["primary"] = primary
    out["coverage_needed"] = coverage_needed(allrec, primary["n_blocks_compared"])
    # beside it, never in place of it: the same three numbers with what each one is (lane-design)
    out["reporting_floor_and_power_assumptions"] = reporting_floor_and_power_assumptions(
        allrec, primary["n_blocks_compared"]
    )
    out["primary_reading"] = measured_reading(primary)
    out["secondary"] = {
        "decrease_only": summarise_measured(runs["matched_real_unknown"], "measured_regulated_decrease_only"),
        "neutral_primary": summarise_measured(runs["matched_neutral"], PRIMARY_QUESTION),
        "unmatched_primary": summarise_measured(runs["unmatched_real_unknown"], PRIMARY_QUESTION),
        "same_instrument_descriptive_only": {
            q: drop_boots(mc.summarise(runs["matched_real_unknown"], q))
            for q in (
                PRIMARY_QUESTION,
                "measured_by_an_element_level_assay",
                "measured_active_but_no_gene_named",
            )
        },
    }
    out["seconds"] = round(time.time() - t0, 1)
    return out


@mf.depends_on_models("alphagenome")  # the model arm is carried through the same draw (R9)
def manifest(chroms: list[str]) -> dict[str, Any]:
    base = mc.manifest(chroms)
    inputs = list(base["inputs"])
    for name in ms.CRISPRI_FILES:
        p = ms.CRISPRI_KNOWLEDGE / name
        if p.exists():
            inputs.append(mf.input_entry(p, partition=ms.CRISPRI_SPLIT_OF[name]))
    for acc in mpra.FILES.values():
        p = mpra.KNOWLEDGE / f"{acc}.bed.gz"
        if p.exists():
            inputs.append(mf.input_entry(p, partition=None))
    p = vista.locus_path(vista.KNOWLEDGE)
    if p.exists():
        inputs.append(mf.input_entry(p, partition=None))
    base["inputs"] = inputs
    base["sources"] = list(base["sources"]) + [
        {
            "accession": "EngreitzLab/CRISPR_comparison EPCrisprBenchmark (Gschwind et al. 2025)",
            "version": "main, as fetched; pinned by sha256",
        },
        {
            "accession": f"ENCODE {mpra.LIBRARY} ({', '.join(mpra.FILES.values())})",
            "version": "as fetched 2026-09-12; pinned by sha256",
        },
        {"accession": "VISTA Enhancer Browser locus table (vista-data)", "version": "main; pinned by sha256"},
        {
            "accession": "this repository, data/results/clause2_matched_control.json at 0c8b82d",
            "version": "git 0c8b82d, the matched primary the gate reproduces",
        },
    ]
    base["parameters"] = {
        **base["parameters"],
        "reciprocal_overlap": ms.RECIPROCAL_OVERLAP,
        "overlap_sensitivity": list(ms.OVERLAP_SENSITIVITY),
        "mpra_active_log2": ms.MPRA_ACTIVE,
        "reporter_label_rule": ms.REPORTER_LABEL_RULE,
        "power_for_negatives": ms.POWER_FOR_NEGATIVES,
        "well_powered": ms.WELL_POWERED,
        "min_blocks_for_a_comparison": MIN_BLOCKS,
        "effect_to_detect": EFFECT_TO_DETECT,
        "primary_question": PRIMARY_QUESTION,
        "denominator_question": DENOMINATOR_QUESTION,
    }
    base["exclusions"] = list(base["exclusions"]) + [
        "a measurement counts for an element only at reciprocal overlap >= 0.5, both ways",
        "a CRISPRi gene counts only if its GENCODE symbol is protein-coding",
        "saturation mutagenesis is base-level: counted, reported, never pooled into an element verdict",
        "lentiMPRA and VISTA name no gene and never enter the primary",
        "blocks no CRISPRi coding pair reaches leave the primary; they are the coverage result",
    ]
    return base


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--chroms", default="", help="comma-separated; default every chromosome with a table")
    ap.add_argument("--no-save", action="store_true", help="print without writing the result")
    args = ap.parse_args(argv)
    chroms = (
        args.chroms.split(",")
        if args.chroms
        else sorted((p.stem for p in cut.ELEMENTS.glob("chr*.json")), key=lambda c: (len(c), c))
    )
    out = collect(chroms)
    if args.no_save or args.chroms:
        print("not saved (a partial run never overwrites the genome's)")
    else:
        print(f"saved {save_result(RESULT, out, manifest=manifest(chroms))}")
    ru = out["coverage"]["real_unknown"]
    print(
        f"coverage: {ru['blocks']} real-unknown blocks, "
        f"{ru['blocks_carrying_a_scored_element']} carry a scored element, "
        f"{ru['blocks_with']['measured_by_any_assay']} carry a measured one, "
        f"{ru['blocks_with'][DENOMINATOR_QUESTION]} carry one CRISPRi tested against a coding gene, "
        f"{ru['blocks_with'][PRIMARY_QUESTION]} carry one measured to move it"
    )
    print("gate passed:", out["gate_reproduces_0c8b82d"]["passed"])
    if "primary" in out:
        p = out["primary"]
        print(f"primary: n={p['n_blocks_compared']} compared blocks, {p.get('matched_difference_points')}")
        print("reading:", out["primary_reading"]["outcome"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
