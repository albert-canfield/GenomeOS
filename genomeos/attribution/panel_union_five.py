# SPDX-License-Identifier: AGPL-3.0-or-later
"""The fifth covariate, the one that runs the other way: segmental duplication put into the union.

`panel_union.py` took four covariates out of the human panel's GC- and replication-timing-matched
background at once -- coding-gene introns, exons of any gene, conserved elements, promoter proximity
-- and reported where that union sits: fossil 98.7%, regulatory 141.0%, neutral 55.5% of the offset,
against largest singles of 56.5 / 65.6 / 30.1 and sums of 101.8 / 129.6 / 61.3. Its own "left undone"
named what it had omitted:

> **Segmental duplication, the fifth covariate, which runs the other way.** The first section
> measured it and it *widens* the offset when removed -- fossil 29%, regulatory 34%, neutral 33%. It
> is not in this union, so the union above is a union of the four that close the offset and says
> nothing about what a fifth that opens it would do to the total.

This module is that five-way arm. Nothing is re-derived: the background is rebuilt through the
panel's own `build_background` by way of `panel_background.tier_ratios`, the decomposition through
`panel_background.explained`, the degeneracy flag at `panel_leftover.DEGENERATE_MIN` through
`panel_union.block_share_after`, the gene model through `panel_union.covariate_model`, the window
comparison through `panel_union.union_windows` and `panel_union.comparisons`, and the four-way
reading through `panel_union.where_it_sits` itself, so that the four-way row in this result is that
lane's own code reading that lane's own arm and cannot drift from it.

**A widening covariate must not be averaged in as though it narrowed the offset.** The four have a
positive `explains`: the sequence they hold is quieter than the rest of the background, so taking it
out raises the background's rate and the tier's ratio falls towards 1. Segmental duplication has a
negative one. Every share it appears in is carried with its sign, the narrowing singles and the
widening singles are summed separately as well as together, and there is no mean anywhere that mixes
them. The consequence is stated rather than smoothed: *a union that contains a widening term can be
SMALLER than the same union without it, because the two kinds of covariate move the background's rate
in opposite directions and part of what the four bought is spent putting the fifth back.* If that is
what the arms do, it is the result.

Which way a covariate runs is not asserted here either. `which_way_each_covariate_runs` reads it off
the background's own kilobases -- the recurring rate where the covariate is against the rate where it
is not, which `panel_background.feature_splits` already computes -- and says whether that rate and
the sign of the arm agree. A covariate noisier than the rest of the background must widen the offset
when it is removed; if the two disagree, the result says so instead of picking one.
"""

from __future__ import annotations

from array import array
from pathlib import Path
from typing import Any

from genomeos.attribution.human_panel import BACKGROUND_BIN, CACHE, TIERS, merge_intervals
from genomeos.attribution.panel_background import (
    CELL,
    CLAIM_WITH_COVERAGE,
    CLAIM_WITHOUT_COVERAGE,
    CONSTRAINT_CACHE,
    COVERAGE_MEASURE,
    SPLIT_AT,
    TRACKS,
    Rebuilt,
    background_bins,
    block_composition,
    composition,
    cover_cells,
    explained,
    feature_splits,
    local_tracks,
    rebuild,
    reproduces_the_panel,
    span_bp,
    tier_ratios,
)
from genomeos.attribution.panel_background import SOURCES as BACKGROUND_SOURCES
from genomeos.attribution.panel_leftover import DEGENERATE_MIN, bases
from genomeos.attribution.panel_union import (
    ABOVE_SUM,
    BELOW_MAX,
    BETWEEN,
    BLOCKS,
    CLAIMS,
    DEGENERATE_WHY,
    MEASURED_EMPTY,
    NEAR,
    NEAR_MAX,
    NEAR_SUM,
    NO_GENCODE,
    NO_OVERSHOOT,
    NO_TRACK,
    OVERSHOOT,
    STRATIFICATIONS,
    UNASSESSED,
    block_share_after,
    comparisons,
    covariate_model,
    intersect,
    union_windows,
    where_it_sits,
)
from genomeos.attribution.panel_union import COVARIATES as FOUR
from genomeos.attribution.panel_union import FROM_LANE as FOUR_FROM_LANE
from genomeos.attribution.panel_union import SOURCES as FOUR_SOURCES
from genomeos.attribution.panel_union import UNION as UNION_FOUR
from genomeos.attribution.panel_union import covariate_intervals as four_intervals
from genomeos.attribution.panel_union import pairwise_overlap as four_pairwise
from genomeos.attribution.panel_union import why_absent as four_why_absent

# The fifth covariate, by the definition of the lane that measured it: `panel_background` reads it
# from the same per-chromosome BED, under the same name, and found it widens the offset.
FIFTH = "segmental_duplication"
COVARIATES = (*FOUR, FIFTH)
# Which of them was measured running the other way. This is a declaration of what the module expects,
# not of what it found: `which_way_each_covariate_runs` reads the direction off the background and
# reports a disagreement rather than taking this tuple's word for it.
WIDENING_WHEN_MEASURED = (FIFTH,)
UNION = "all_five"
ARMS = (*(f"without_{c}" for c in COVARIATES), f"without_{UNION_FOUR}", f"without_{UNION}")
UNIONS: dict[str, tuple[str, ...]] = {UNION_FOUR: FOUR, UNION: COVARIATES}

FROM_LANE = {
    **FOUR_FROM_LANE,
    FIFTH: "panel_background.local_tracks: UCSC genomicSuperDups, per-chromosome BED",
}
SOURCES = {**FOUR_SOURCES, FIFTH: BACKGROUND_SOURCES[FIFTH]}

# The floor `panel_background.reads_explained` draws between "accounts for none" and a share, used
# here for the same job on a difference between two unions. `tests/test_panel_union_five.py` asserts
# it is still that function's own default, so the two cannot drift apart.
MOVED_AT_ALL = 0.02
# Two per-kilobase rates within this relative distance of each other are not called one quieter than
# the other; the rates themselves are reported beside the verdict either way.
SAME_WITHIN = 0.02

# -- the sign, as data ---------------------------------------------------------------------------
SIGN_RULE = (
    "a covariate that widens the offset when it is removed is carried with a negative sign and is "
    "never averaged into a share with the covariates that narrow it: the narrowing singles and the "
    "widening singles are summed separately as well as together, and no mean mixes them"
)
WITH_A_WIDENING_TERM = (
    "a union that contains a widening term can explain LESS than the same union without it. The four "
    "narrowing covariates hold sequence quieter than the rest of the background, so removing them "
    "raises the background's rate and every tier's ratio falls towards 1; the widening one holds "
    "sequence noisier than the rest, so removing it lowers that rate again and the ratio rises away "
    "from 1. Part of what the four bought is spent putting the fifth back, and a five-way union "
    "smaller than the four-way one is that arithmetic and not an anomaly to smooth"
)
FIFTH_LOWERS = (
    "the fifth covariate LOWERS the union: adding a term that widens the offset on its own takes "
    "back part of what the four narrowing covariates bought, which is what a signed sum predicts"
)
FIFTH_RAISES = (
    "the fifth covariate RAISES the union, which its own arm does not predict: the exclusions "
    "interact, and the sign of a single arm is not the sign of its contribution to a union"
)
FIFTH_MOVES_NOTHING = (
    "the fifth covariate moves the union by less than the floor this project draws under a share, so "
    "the five-way union and the four-way union are the same answer"
)
NOT_COMPARABLE = "one of the two unions has no share, so the difference between them is not a number"

# -- which way a covariate runs, read off the background's own kilobases --------------------------
QUIETER_THAN_THE_REST = (
    "measured: the kilobases of the background holding it are quieter than the kilobases that do "
    "not, so taking it out raises the background's rate and a tier's ratio falls towards 1"
)
NOISIER_THAN_THE_REST = (
    "measured: the kilobases of the background holding it are noisier than the kilobases that do "
    "not, so taking it out lowers the background's rate and a tier's ratio rises away from 1 -- this "
    "is a covariate that runs the other way, and removing it widens the offset rather than closing it"
)
SAME_AS_THE_REST = "measured: the kilobases holding it read the same as the kilobases that do not"
RATE_NOT_ASSESSED = "the covariate has no kilobase on either side on this chromosome"
AGREE = "the rate and the arm agree"
DISAGREE = (
    "the rate and the arm disagree: the kilobases holding this covariate read one way and the "
    "rebuilt background reads the other, so the arm is not a restatement of the rate and neither is "
    "reported as though it settled the other"
)


# ================================================================================================
# the five interval sets
# ================================================================================================
def covariate_intervals(model: dict[str, Any], tracks: dict[str, Any]) -> dict[str, list[tuple[int, int]]]:
    """The five covariates as merged interval sets: `panel_union`'s four, and the fifth beside them."""
    ivs: dict[str, list[tuple[int, int]]] = dict(four_intervals(model, tracks))
    ivs[FIFTH] = tracks.get(FIFTH, [])
    return ivs


def why_absent(name: str, model: dict[str, Any], tracks: dict[str, Any]) -> str:
    """Measured-and-absent or never-looked, as two different sentences, never as one zero."""
    if name == FIFTH:
        return NO_TRACK if FIFTH in tracks.get("missing", []) else MEASURED_EMPTY
    return four_why_absent(name, model, tracks)


def pairwise_overlap(ivs: dict[str, list[tuple[int, int]]], length: int) -> dict[str, Any]:
    """Bases shared by each pair of the five, with the four-way lane's own totals beside them.

    The six pairs among the four, and the four-way union and sum, are `panel_union.pairwise_overlap`'s
    output unchanged -- the same function on the same intervals -- so the 168.4 Mb that lane reported
    double-counted is this lane's number for the same thing rather than a second computation of it.
    The four pairs the fifth adds, and the five-way union and sum, are computed here the same way.
    """
    four = four_pairwise({n: ivs.get(n, []) for n in FOUR}, length)
    pairs: dict[str, Any] = dict(four["pairs"])
    for a in FOUR:
        key = f"{a}_and_{FIFTH}"
        if not ivs.get(a) or not ivs.get(FIFTH):
            pairs[key] = {
                "shared_bases": None,
                "bases_a": None,
                "bases_b": None,
                "share_of_the_smaller": None,
            }
            continue
        ba, bb = bases(ivs[a]), bases(ivs[FIFTH])
        shared = ba + bb - bases(ivs[a] + ivs[FIFTH])
        smaller = min(ba, bb)
        pairs[key] = {
            "shared_bases": shared,
            "bases_a": ba,
            "bases_b": bb,
            "share_of_the_smaller": round(shared / smaller, 4) if smaller else None,
        }
    names = [n for n in COVARIATES if ivs.get(n)]
    union = bases(merge_intervals([iv for n in names for iv in ivs[n]]))
    summed = sum(bases(ivs[n]) for n in names)
    return {
        "pairs": pairs,
        "covariates_with_bases": names,
        "union_bases": union,
        "sum_of_the_singles_bases": summed,
        "bases_counted_more_than_once_by_the_sum": summed - union,
        "the_sum_over_the_union": round(summed / union, 4) if union else None,
        "chromosome_bases": length,
        # the four-way lane's own numbers, from its own function, for the row this one is read against
        "the_four_way": {
            "covariates_with_bases": four["covariates_with_bases"],
            "union_bases": four["union_bases"],
            "sum_of_the_singles_bases": four["sum_of_the_singles_bases"],
            "bases_counted_more_than_once_by_the_sum": four["bases_counted_more_than_once_by_the_sum"],
        },
    }


# ================================================================================================
# the arms and the arrays every table reads
# ================================================================================================
def arm_intervals(ivs: dict[str, list[tuple[int, int]]]) -> dict[str, list[tuple[int, int]]]:
    """The five singles and the two unions: the four-way as `panel_union` builds it, and the five-way."""
    arms: dict[str, list[tuple[int, int]]] = {n: ivs[n] for n in COVARIATES}
    for name, members in UNIONS.items():
        arms[name] = merge_intervals([iv for n in members for iv in ivs[n]])
    return arms


def arm_arrays(rb: Rebuilt, arms: dict[str, list[tuple[int, int]]]) -> dict[str, array]:
    """A coverage array per arm, per arm-inside-a-block, plus the tiers and the canonical-CDS cut.

    This is `panel_union.union_arrays` over an arms dict given to it rather than over that module's
    own four, because two unions are wanted here and it builds one. The test asserts the arrays it
    produces for the four and their union are identical to that function's, so the generalisation
    cannot become a second definition of an arm.
    """
    cells = rb.length // CELL + 2
    blocks = merge_intervals([iv for t in TIERS for iv in rb.tier_intervals(t)])
    sets: dict[str, list[tuple[int, int]]] = {"cds_canonical": rb.cds_all, BLOCKS: blocks}
    sets.update(arms)
    for name, a in arms.items():
        sets[f"{BLOCKS}_and_{name}"] = intersect(blocks, a, rb.length) if a else []
    for t in TIERS:
        sets[f"in_{t}"] = rb.tier_intervals(t)
    return {k: cover_cells(v, cells) for k, v in sets.items()}


def decompose(
    rb: Rebuilt,
    arms: dict[str, list[tuple[int, int]]],
    rows: list[dict[str, Any]],
    model: dict[str, Any],
    tracks: dict[str, Any],
) -> dict[str, Any]:
    """Every tier's ratio with each covariate out, then with the four out, then with all five out.

    The four-way arm is rebuilt here rather than read from the committed four-way result, so that the
    two unions are measured on the same chromosomes by the same call on the same day and their
    difference is a difference between arms and not between runs.
    """
    out: dict[str, Any] = {"as_the_panel_builds_it": tier_ratios(rb, rb.cds_all)}
    for name, excluded in arms.items():
        key = f"without_{name}"
        members = UNIONS.get(name, (name,))
        if not excluded:
            out[key] = {
                "assessed": False,
                "degenerate": False,
                "why": "; ".join(sorted({why_absent(n, model, tracks) for n in members if not arms[n]})),
            }
            continue
        share = block_share_after(rows, name)
        degenerate = share is not None and share >= DEGENERATE_MIN
        entry: dict[str, Any] = {
            "assessed": not degenerate,
            "block_share_of_the_remaining_background": share,
            "degenerate": degenerate,
            "degenerate_threshold": DEGENERATE_MIN,
            "background_bases_left": BACKGROUND_BIN * len(rows) - sum(r[name] for r in rows),
            **tier_ratios(rb, merge_intervals(rb.cds_all + excluded)),
        }
        if name in UNIONS:
            entry["covariates_included"] = [n for n in members if arms[n]]
            entry["covariates_left_out"] = {n: why_absent(n, model, tracks) for n in members if not arms[n]}
            entry["widening_covariates_included"] = [
                n for n in members if n in WIDENING_WHEN_MEASURED and arms[n]
            ]
        if degenerate:
            entry["why"] = DEGENERATE_WHY
        out[key] = entry
    return out


# ================================================================================================
# which way each covariate runs, read off the background rather than asserted
# ================================================================================================
def which_way_it_runs(with_per_kb: float | None, without_per_kb: float | None) -> str:
    """The direction the background's own kilobases say, before any arm is rebuilt."""
    if with_per_kb is None or without_per_kb is None or not without_per_kb:
        return RATE_NOT_ASSESSED
    rel = (with_per_kb - without_per_kb) / without_per_kb
    if abs(rel) < SAME_WITHIN:
        return SAME_AS_THE_REST
    return NOISIER_THAN_THE_REST if rel > 0 else QUIETER_THAN_THE_REST


def which_way_each_covariate_runs(splits: dict[str, Any], exp: dict[str, Any]) -> dict[str, Any]:
    """Per covariate: the rate where it is, the rate where it is not, and whether the arm agrees.

    The rate is descriptive and the arm is the measurement; they are two different readings of the
    same covariate and this says whether they point the same way. `exp` is the per-tier decomposition,
    so the arm's sign is taken on every tier the chromosome carries rather than on a chosen one.
    """
    out: dict[str, Any] = {}
    for name in COVARIATES:
        row = splits.get(name) or {}
        reads = which_way_it_runs(row.get("with_per_kb"), row.get("without_per_kb"))
        signs = {
            t: exp[t][f"without_{name}"]["explains"]
            for t in exp
            if isinstance(exp[t].get(f"without_{name}"), dict)
            and exp[t][f"without_{name}"].get("explains") is not None
        }
        widens = [t for t, v in signs.items() if v < -MOVED_AT_ALL]
        narrows = [t for t, v in signs.items() if v > MOVED_AT_ALL]
        entry: dict[str, Any] = {
            "recurring_per_kb_where_it_is": row.get("with_per_kb"),
            "recurring_per_kb_where_it_is_not": row.get("without_per_kb"),
            "share_of_background_bases": row.get("share_of_background_bases"),
            "reads": reads,
            "declared_as_widening": name in WIDENING_WHEN_MEASURED,
            "tiers_whose_arm_widens_the_offset": sorted(widens),
            "tiers_whose_arm_narrows_the_offset": sorted(narrows),
            "explains_by_tier_with_its_sign": {t: signs[t] for t in sorted(signs)},
        }
        if reads in (RATE_NOT_ASSESSED, SAME_AS_THE_REST) or not signs:
            entry["the_rate_and_the_arm"] = UNASSESSED
        elif (reads == NOISIER_THAN_THE_REST) == (len(widens) >= len(narrows)):
            entry["the_rate_and_the_arm"] = AGREE
        else:
            entry["the_rate_and_the_arm"] = DISAGREE
        out[name] = entry
    return out


# ================================================================================================
# where the five-way union sits, with the widening term kept apart from the narrowing ones
# ================================================================================================
def signed_singles(exp: dict[str, Any]) -> dict[str, float]:
    """Each single arm's share with its sign, for the arms that were assessed and not degenerate."""
    return {
        n: exp[f"without_{n}"]["explains"]
        for n in COVARIATES
        if isinstance(exp.get(f"without_{n}"), dict) and exp[f"without_{n}"].get("explains") is not None
    }


def _position(union: float, lo: float, hi: float) -> tuple[float | None, str]:
    """Where a union sits between one end of an interval and the other, and the name for it."""
    pos = (union - lo) / (hi - lo) if hi > lo else None
    if union < lo:
        return pos, BELOW_MAX
    if union > hi:
        return pos, ABOVE_SUM
    if pos is None or pos <= NEAR:
        return pos, NEAR_MAX
    if pos >= 1 - NEAR:
        return pos, NEAR_SUM
    return pos, BETWEEN


def where_the_five_sit(exp: dict[str, Any]) -> dict[str, Any]:
    """The five-way union beside the four-way one, with every share carried with its sign.

    The bracket a four-way union had to sit in -- between the largest single and the sum of the
    singles -- is not a bracket for five, because a negative term makes the sum of the singles smaller
    than the largest of them. So two intervals are reported, not one: the narrowing singles' interval,
    which is the four-way lane's construction on the arms that behave as it assumed, and the signed
    sum of all five, which is where a reader who adds every arm up would land. Neither is called the
    union; the union is measured.
    """
    singles = signed_singles(exp)
    narrowing = {n: v for n, v in singles.items() if v > 0}
    widening = {n: v for n, v in singles.items() if v < 0}
    row5 = exp.get(f"without_{UNION}") or {}
    row4 = exp.get(f"without_{UNION_FOUR}") or {}
    five = row5.get("explains") if isinstance(row5, dict) else None
    four = row4.get("explains") if isinstance(row4, dict) else None
    out: dict[str, Any] = {
        "union_explains": five,
        "union_assessed": bool(isinstance(row5, dict) and row5.get("assessed")),
        "union_degenerate": bool(isinstance(row5, dict) and row5.get("degenerate")),
        "sign_rule": SIGN_RULE,
        "what_a_union_with_a_widening_term_means": WITH_A_WIDENING_TERM,
        "per_single_with_its_sign": {n: singles[n] for n in COVARIATES if n in singles},
        "narrowing_singles": [n for n in COVARIATES if n in narrowing],
        "widening_singles": [n for n in COVARIATES if n in widening],
        "singles_that_move_nothing": [n for n in COVARIATES if n in singles and singles[n] == 0],
        "singles_not_in_any_sum": [n for n in COVARIATES if n not in singles],
        "largest_narrowing_single": round(max(narrowing.values()), 4) if narrowing else None,
        "sum_of_the_narrowing_singles": round(sum(narrowing.values()), 4) if narrowing else None,
        "sum_of_the_widening_singles": round(sum(widening.values()), 4) if widening else None,
        "signed_sum_of_all_five": round(sum(singles.values()), 4) if singles else None,
        "the_widening_singles_are_not_averaged_in": True,
        "four_way_union_explains": four,
    }
    out["union_overshoots_the_offset"] = bool(five is not None and five > 1)
    if out["union_overshoots_the_offset"]:
        out["overshoot_reads"] = OVERSHOOT
    else:
        out["overshoot_reads"] = NO_OVERSHOOT if five is not None else UNASSESSED
    if five is None or four is None:
        out["the_fifth_covariate_moves_the_union_by"] = None
        out["the_fifth_covariate_reads"] = NOT_COMPARABLE
    else:
        delta = round(five - four, 4)
        out["the_fifth_covariate_moves_the_union_by"] = delta
        moved = FIFTH_LOWERS if delta < 0 else FIFTH_RAISES
        out["the_fifth_covariate_reads"] = FIFTH_MOVES_NOTHING if abs(delta) < MOVED_AT_ALL else moved
    lo, hi = out["largest_narrowing_single"], out["sum_of_the_narrowing_singles"]
    out["near_threshold"] = NEAR
    if five is None or lo is None or hi is None:
        out["position_among_the_narrowing_singles"] = None
        out["reads"] = UNASSESSED
        return out
    pos, reads = _position(five, lo, hi)
    out["position_among_the_narrowing_singles"] = round(pos, 4) if pos is not None else None
    out["reads"] = reads
    signed = out["signed_sum_of_all_five"]
    out["union_against_the_signed_sum_of_all_five"] = round(five - signed, 4) if signed is not None else None
    return out


# ================================================================================================
# the windows, which are `panel_union`'s own with the fifth covariate counted into them
# ================================================================================================
def add_the_fifth(rows: list[dict[str, Any]], feats: dict[str, array]) -> list[dict[str, Any]]:
    """Count the fifth covariate into windows `panel_union.union_windows` built counting four.

    The window machinery -- the subsampling step, the coverage, the GC, the timing, the tier and the
    hit -- is that function's, unchanged and uncopied. Only the two annotation fields it derives from
    its own four covariates are extended, and a window is covariate-free here when not one base of any
    of the FIVE lies in it, which is the rule the five-way exclusion arm applies.
    """
    arr = feats.get(FIFTH)
    for r in rows:
        held = span_bp(arr, r["start"], r["start"] + CELL) if arr is not None else 0
        r["covariates_held"] += 1 if held else 0
        r["covariate_bases"] += held
        r["holds_a_widening_covariate"] = bool(held)
    return rows


def window_cost(rows: list[dict[str, Any]], cost: dict[str, Any]) -> dict[str, Any]:
    """`union_windows`' own cost line, with the counts that named four re-counted over five."""
    out = {k: v for k, v in cost.items() if k != "windows_free_of_all_four"}
    out["windows_free_of_all_five"] = sum(1 for r in rows if not r["covariates_held"])
    out["windows_holding_at_least_one"] = sum(1 for r in rows if r["covariates_held"])
    out["windows_holding_a_widening_covariate"] = sum(1 for r in rows if r["holds_a_widening_covariate"])
    return out


# ================================================================================================
# one chromosome
# ================================================================================================
def analyse(
    chrom: str,
    cache: Path = CACHE,
    results_dir: Path | None = None,
    gff: Path | None = None,
    tracks_dir: Path = TRACKS,
    constraint_dir: Path = CONSTRAINT_CACHE,
) -> dict[str, Any]:
    """Take all five covariates out of one chromosome's background, and the four, in one run."""
    rb = rebuild(chrom, cache, results_dir)
    model = covariate_model(chrom, gff)
    tracks = local_tracks(chrom, tracks_dir, constraint_dir)
    ivs = covariate_intervals(model, tracks)
    present = [n for n in COVARIATES if ivs[n]]
    if not present:
        return {"chrom": chrom, "assessed": False, "why": NO_GENCODE if not model else NO_TRACK}
    arms = arm_intervals(ivs)
    feats = arm_arrays(rb, arms)
    tss = model.get("tss", [])
    rows, dropped = background_bins(rb, feats, tss)
    names = [*COVARIATES, UNION_FOUR, UNION]
    absent = [n for n in COVARIATES if not ivs[n]]
    groups: dict[str, Any] = {
        "background": composition(rows, names, absent),
        "background_free_of_all_five": composition([r for r in rows if r[UNION] < SPLIT_AT], names, absent),
    }
    for n in names:
        groups[f"kilobases_holding_{n}"] = composition([r for r in rows if r[n] >= SPLIT_AT], names, absent)
    for t in TIERS:
        groups[f"tier_{t}"] = block_composition(rb, t, feats, tss, names, absent)
    offset = decompose(rb, arms, rows, model, tracks)
    exp = {t: explained(offset, t) for t in TIERS}
    splits = feature_splits(rows, names, absent)
    win, cost = union_windows(rb, feats, tss)
    add_the_fifth(win, feats)
    cost = window_cost(win, cost)
    comp = comparisons(win)
    comp["control_windows_free_of_all_five"] = comp.pop("control_windows_free_of_all_four")
    matched = [
        v[STRATIFICATIONS[-1]]["targets_matched"]
        for v in comp.values()
        if isinstance(v, dict) and v.get("assessed")
    ]
    return {
        "chrom": chrom,
        "assessed": True,
        "assemblies": rb.panel.n,
        "evidence": SOURCES,
        "covariate_definitions": FROM_LANE,
        "sign_rule": SIGN_RULE,
        "what_a_union_with_a_widening_term_means": WITH_A_WIDENING_TERM,
        # one of two sentences and never a third, carried as data so no reader can upgrade it
        "claim_available": CLAIM_WITH_COVERAGE if matched and all(matched) else CLAIM_WITHOUT_COVERAGE,
        "reproduces_the_panel": reproduces_the_panel(rb, offset["as_the_panel_builds_it"]),
        "chromosome_bases": rb.length,
        "eligible": {
            "chromosome_bases": rb.length,
            "kilobases_in_the_panel_alignment": len(rows) + sum(dropped.values()),
            "covariates_declared": list(COVARIATES),
            "covariates_with_bases_on_this_chromosome": present,
            "covariates_without_bases": {n: why_absent(n, model, tracks) for n in absent},
            "covariates_declared_widening": list(WIDENING_WHEN_MEASURED),
            "tiers_declared": list(TIERS),
            "tiers_with_blocks_on_this_chromosome": [t for t in TIERS if rb.tier_intervals(t)],
            "windows_in_the_panel_alignment": cost["cells_in_panel_blocks"],
            "genes_in_the_annotation": model.get("genes", 0),
            "coding_tss_read": len(tss),
        },
        "overlap": pairwise_overlap(ivs, rb.length),
        "groups": groups,
        "feature_splits": splits,
        "which_way_each_covariate_runs": which_way_each_covariate_runs(splits, exp),
        "offset": offset,
        "explained": exp,
        "where_the_union_sits": {
            t: {
                # the four-way row is `panel_union`'s own reading of its own arm, not a copy of it
                "the_four_way_as_it_stands": where_it_sits(exp[t]),
                "the_five_way": where_the_five_sit(exp[t]),
            }
            for t in TIERS
        },
        "comparisons": comp,
        "coverage": {
            "kilobases_in_panel_blocks": len(rows) + sum(dropped.values()),
            "kilobases_measured": len(rows),
            "kilobases_not_measured": dropped,
            "kilobases_free_of_all_five": sum(1 for r in rows if r[UNION] < SPLIT_AT),
            "tracks_absent": tracks.get("missing", []),
            "gencode_read": bool(model),
            "coverage_measure": COVERAGE_MEASURE,
            "windows": cost,
        },
    }


__all__ = [
    "ARMS",
    "CLAIMS",
    "COVARIATES",
    "FIFTH",
    "STRATIFICATIONS",
    "UNION",
    "UNION_FOUR",
    "analyse",
    "arm_arrays",
    "arm_intervals",
    "covariate_intervals",
    "decompose",
    "pairwise_overlap",
    "where_the_five_sit",
    "which_way_each_covariate_runs",
    "why_absent",
]
