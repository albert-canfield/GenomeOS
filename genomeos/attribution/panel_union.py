# SPDX-License-Identifier: AGPL-3.0-or-later
"""All four covariates out of the panel's background at once: where the union sits, and whether it can.

Two lanes have measured what explains the offset by which every non-coding tier reads above the human
panel's GC- and replication-timing-matched background, and neither measured the other's covariate:

- `panel_background.py` took out exons of any gene, phastCons conserved elements and sequence within
  2 kb of a coding TSS, one at a time and then together: 41% (fossil), 56% (regulatory), 25%
  (neutral). Segmental duplication runs the other way and masks about a third of the offset.
- `panel_leftover.py` took out the classes of the 2,006 Mb the budget never cut into a block, and
  found that coding-gene introns alone account for 56% / 66% / 30% -- more than the first lane's three
  covariates together.

**Nobody has removed all four at once.** The union lies somewhere between the larger single share and
the sum of the singles, because the four overlap heavily -- a conserved element inside a coding intron
is one base in both arms and two bases in the sum -- and quoting either end as "explained" would be
wrong. This module computes the union arm, beside each single arm and beside the sum, by the same
method on the same chromosomes.

Nothing here invents a third way of asking. The background is rebuilt through the panel's own
`build_background`, the arms go through `panel_background.tier_ratios`, the decomposition through
`panel_background.explained`, and the two-group comparison through `compare.standardised` after
`compare.input_presence` has said whether each input was bought or free.

**Degeneracy is the live risk, and both earlier lanes hit it.** Taking every gene body out of this
background leaves a remainder that is 99.8% block sequence, which is a tier average standing in for a
control. The union arm removes less than that but a great deal more than any single arm, so every arm
here carries `block_share_of_the_remaining_background` and the `degenerate` flag at
`panel_leftover.DEGENERATE_MIN`, a degenerate arm is excluded from every reported share, and if the
union arm is degenerate then that is the result: the question cannot be answered this way.
"""

from __future__ import annotations

from array import array
from collections import Counter
from pathlib import Path
from typing import Any

from genomeos.attribution.human_panel import (
    BACKGROUND_BIN,
    CACHE,
    MIN_RECURRING,
    TIERS,
    merge_intervals,
    rt_over,
    rt_stratum,
)
from genomeos.attribution.panel_background import (
    CELL,
    CLAIM_WITH_COVERAGE,
    CLAIM_WITHOUT_COVERAGE,
    CONSTRAINT_CACHE,
    COVERAGE_MEASURE,
    PROMOTER,
    SPLIT_AT,
    TRACKS,
    WINDOW_CAP,
    Rebuilt,
    background_bins,
    block_composition,
    composition,
    cover_cells,
    coverage_over,
    explained,
    feature_splits,
    informing,
    local_tracks,
    nearest,
    rebuild,
    reproduces_the_panel,
    span_bp,
    tier_ratios,
)
from genomeos.attribution.panel_leftover import (
    COVERAGE_FULL,
    DEGENERATE_MIN,
    LADDER,
    bases,
    subtract,
)
from genomeos.compare import input_presence, standardised

# The four covariates, each one exactly as the lane that measured it defined it. The names are that
# lane's names, so that a row here and a row there can be set side by side without a translation.
COVARIATES = ("coding_gene_intron", "exon_any", "conserved_elements", "promoter")
FROM_LANE = {
    "coding_gene_intron": "panel_leftover.classify: a protein-coding gene's span, minus every exon of it",
    "exon_any": "panel_background.gene_features: every exon of every transcript of every gene",
    "conserved_elements": "panel_background.local_tracks: UCSC phastConsElements100way",
    "promoter": f"panel_background._ivs: within {PROMOTER} bases of a protein-coding gene's TSS",
}
UNION = "all_four"
ARMS = (*(f"without_{c}" for c in COVARIATES), f"without_{UNION}")
BLOCKS = "in_a_block"  # the merged tier blocks, as one interval set: the thing an arm must not leave alone
STRATIFICATIONS = tuple(LADDER)

SOURCES = {
    "gene_model": "GENCODE 50 (GRCh38.p14): gene, ncRNA_gene and pseudogene spans, exons, CDS, coding TSS",
    "conserved_elements": "UCSC phastConsElements100way, the cache attribution/constraint.py writes",
    "background": "human_panel.build_background, rebuilt from the local store, nothing re-derived",
    "blocks": "the committed human_panel_<chrom> blocks and their tiers",
    "coverage": COVERAGE_MEASURE,
}

NO_TRACK = "no track for this chromosome"
NO_GENCODE = "GENCODE not read for this chromosome"
MEASURED_EMPTY = "read, and empty on this chromosome"
DEGENERATE_WHY = (
    "the remaining background is almost entirely block sequence, so this arm is a tier average against "
    "a tier and is reported without being counted as evidence"
)

# Where the union sits between the largest single arm and the sum of the singles, as a name rather
# than as a reader's arithmetic. NEAR_MAX means the covariates are so nested that adding the other
# three to the largest buys almost nothing; NEAR_SUM means they are nearly disjoint in their effect.
NEAR_MAX = "the union is near the largest single arm: the covariates are nested, and the sum double-counts"
NEAR_SUM = "the union is near the sum of the singles: the covariates act on nearly separate sequence"
BETWEEN = "the union sits between the largest single arm and the sum, as overlapping covariates do"
BELOW_MAX = "the union explains LESS than its largest single arm: the arms interfere and do not add"
ABOVE_SUM = "the union explains MORE than the sum of its singles: the arms reinforce rather than overlap"
UNASSESSED = "not assessed"
NEAR = 0.15  # within this share of one end of the interval, the union is called near that end

# An explained share above 1 does not mean the offset is more than explained. It means the exclusion
# carried the ratio PAST 1 and out the other side: the tier now reads BELOW the background the arm
# leaves. That is a different fact from "explained", it is the fact a reader is most likely to
# mis-quote, and it is a flag in the result rather than a sentence in a document.
OVERSHOOT = (
    "the union does not close this tier's offset, it crosses it: with all four covariates gone the tier "
    "reads BELOW the background the exclusion leaves, so the share above 100% is over-correction and not "
    "explanation, and the covariates removed were quieter than the tier rather than merely as quiet"
)
NO_OVERSHOOT = "the union moves this tier's ratio towards the background without crossing it"

# The one sentence a coverage stratum licenses and the one it does not, carried from the lanes that
# wrote them so that no reader can upgrade the weaker by paraphrase. "Nothing explains it" is not
# available at any sample size and is not a constant in this module.
CLAIMS = (CLAIM_WITHOUT_COVERAGE, CLAIM_WITH_COVERAGE)


# ================================================================================================
# the four interval sets, each read once
# ================================================================================================
def covariate_model(chrom: str, gff: Path | None = None) -> dict[str, Any]:
    """The gene-model pieces the four covariates are cut from, parsing the annotation once.

    `panel_background.gene_features` and `panel_leftover.gene_model` each parse GENCODE for their own
    lane, and calling both would parse it twice per chromosome for pieces that come from one walk of
    the same genes. This walks it once and produces both lanes' pieces; `tests/test_panel_union.py`
    asserts that what comes out equals what those two functions return on the same annotation, so the
    saving cannot become a drift.
    """
    from genomeos.genome.annotation import Annotation, default_gencode

    gff = gff or default_gencode({chrom})
    if gff is None:
        return {}
    ann = Annotation.from_gff3(gff, {chrom})
    coding: list[tuple[int, int]] = []
    coding_exons: list[tuple[int, int]] = []
    exon_any: list[tuple[int, int]] = []
    tss: list[int] = []
    genes = 0
    for g in ann.genes.values():
        if g.locus.chrom != chrom:
            continue
        genes += 1
        exons = [(x.start, x.end) for t in g.transcripts.values() for x in t.exons]
        exon_any += exons
        if g.type == "protein_coding":
            coding.append((g.locus.start, g.locus.end))
            coding_exons += exons
            tss.append(g.locus.start if "MINUS" not in str(g.locus.strand) else g.locus.end)
    return {
        "genes": genes,
        "coding_spans": merge_intervals(coding),
        "coding_exons": merge_intervals(coding_exons),
        "exon_any": merge_intervals(exon_any),
        "tss": sorted(tss),
    }


def covariate_intervals(model: dict[str, Any], tracks: dict[str, Any]) -> dict[str, list[tuple[int, int]]]:
    """The four covariates as merged interval sets, by the definitions `FROM_LANE` names."""
    if not model:
        return {n: [] for n in COVARIATES}
    return {
        "coding_gene_intron": subtract(model["coding_spans"], model["coding_exons"]),
        "exon_any": model["exon_any"],
        "conserved_elements": tracks.get("conserved_elements", []),
        "promoter": merge_intervals([(max(0, p - PROMOTER), p + PROMOTER) for p in model["tss"]]),
    }


def why_absent(name: str, model: dict[str, Any], tracks: dict[str, Any]) -> str:
    """Measured-and-absent or never-looked, as two different sentences, never as one zero."""
    if name == "conserved_elements" and "conserved_elements" in tracks.get("missing", []):
        return NO_TRACK
    if name != "conserved_elements" and not model:
        return NO_GENCODE
    return MEASURED_EMPTY


def intersect(a: list[tuple[int, int]], b: list[tuple[int, int]], length: int) -> list[tuple[int, int]]:
    """Every base in both, through `panel_leftover.subtract` so there is one interval arithmetic here."""
    return subtract(a, subtract([(0, length)], b))


def pairwise_overlap(ivs: dict[str, list[tuple[int, int]]], length: int) -> dict[str, Any]:
    """How many bases each pair of covariates shares -- the number that decides where the union sits.

    Two covariates that never touch make the union the sum of their effects; two that are nested make
    it the larger of them. This is the cheap half of the answer and it is computed before any
    background is rebuilt, so the expensive half can be read against it.
    """
    out: dict[str, Any] = {}
    names = [n for n in COVARIATES if ivs.get(n)]
    for i, a in enumerate(names):
        for b in names[i + 1 :]:
            shared = bases(ivs[a]) + bases(ivs[b]) - bases(ivs[a] + ivs[b])
            smaller = min(bases(ivs[a]), bases(ivs[b]))
            out[f"{a}_and_{b}"] = {
                "shared_bases": shared,
                "bases_a": bases(ivs[a]),
                "bases_b": bases(ivs[b]),
                "share_of_the_smaller": round(shared / smaller, 4) if smaller else None,
            }
    for i, a in enumerate(COVARIATES):
        for b in COVARIATES[i + 1 :]:
            out.setdefault(
                f"{a}_and_{b}",
                {"shared_bases": None, "bases_a": None, "bases_b": None, "share_of_the_smaller": None},
            )
    union = merge_intervals([iv for n in names for iv in ivs[n]])
    summed = sum(bases(ivs[n]) for n in names)
    return {
        "pairs": out,
        "covariates_with_bases": names,
        "union_bases": bases(union),
        "sum_of_the_singles_bases": summed,
        "bases_counted_more_than_once_by_the_sum": summed - bases(union),
        "chromosome_bases": length,
    }


# ================================================================================================
# the arrays every table reads
# ================================================================================================
def union_arrays(
    rb: Rebuilt, ivs: dict[str, list[tuple[int, int]]]
) -> tuple[dict[str, array], dict[str, list[tuple[int, int]]]]:
    """A coverage array per covariate, per arm, and per arm-inside-a-block, plus the tiers and the CDS cut.

    The arm-inside-a-block arrays are what the degeneracy check needs: an arm takes bases out of the
    blocks as well as out of the background around them, and the share of the REMAINING background
    that is block sequence cannot be computed without knowing how much of each the arm took.
    """
    cells = rb.length // CELL + 2
    blocks = merge_intervals([iv for t in TIERS for iv in rb.tier_intervals(t)])
    arms: dict[str, list[tuple[int, int]]] = {n: ivs[n] for n in COVARIATES}
    arms[UNION] = merge_intervals([iv for n in COVARIATES for iv in ivs[n]])
    sets: dict[str, list[tuple[int, int]]] = {"cds_canonical": rb.cds_all, BLOCKS: blocks}
    sets.update(arms)
    for name, a in arms.items():
        sets[f"{BLOCKS}_and_{name}"] = intersect(blocks, a, rb.length) if a else []
    for t in TIERS:
        sets[f"in_{t}"] = rb.tier_intervals(t)
    return {k: cover_cells(v, cells) for k, v in sets.items()}, arms


def block_share_after(rows: list[dict[str, Any]], arm: str) -> float | None:
    """Of the background bases an arm leaves, what share lies inside an unknown block.

    `panel_leftover.block_share_left` computes this for its own arms from the fact that its classes
    and the blocks partition the sequence, which is true of the never-cut classes and false of these
    four: a conserved element or an exon can sit inside a block, so what an arm leaves has to be read
    as (block bases minus the block bases the arm took) over (all bases minus the bases the arm took).
    Same quantity, same line at `DEGENERATE_MIN`, computed the one way that is right here. The
    denominator is the whole kilobase the background bins by, as it is there.
    """
    span = BACKGROUND_BIN * len(rows)
    if not span:
        return None
    gone = sum(r[arm] for r in rows)
    block_gone = sum(r[f"{BLOCKS}_and_{arm}"] for r in rows)
    left = span - gone
    if left <= 0:
        return None
    return round(max(0, sum(r[BLOCKS] for r in rows) - block_gone) / left, 4)


def decompose(
    rb: Rebuilt,
    arms: dict[str, list[tuple[int, int]]],
    rows: list[dict[str, Any]],
    model: dict[str, Any],
    tracks: dict[str, Any],
) -> dict[str, Any]:
    """Every tier's ratio with each covariate, and then all four, taken out of the background.

    Each arm is rebuilt through the panel's own `build_background` by way of
    `panel_background.tier_ratios`, which is the same call the two earlier lanes made. A degenerate
    arm keeps its numbers, carries the flag and the block share beside them, and is marked
    `assessed: false` so that `panel_background.explained` gives it no share.
    """
    out: dict[str, Any] = {"as_the_panel_builds_it": tier_ratios(rb, rb.cds_all)}
    for name, excluded in arms.items():
        key = f"without_{name}"
        if not excluded:
            missing = [n for n in COVARIATES if not arms[n]] if name == UNION else [name]
            out[key] = {
                "assessed": False,
                "degenerate": False,
                "why": "; ".join(sorted({why_absent(n, model, tracks) for n in missing})),
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
        if name == UNION:
            entry["covariates_included"] = [n for n in COVARIATES if arms[n]]
            entry["covariates_left_out"] = {
                n: why_absent(n, model, tracks) for n in COVARIATES if not arms[n]
            }
        if degenerate:
            entry["why"] = DEGENERATE_WHY
        out[key] = entry
    return out


# ================================================================================================
# where the union sits, per tier
# ================================================================================================
def where_it_sits(exp: dict[str, Any]) -> dict[str, Any]:
    """The union's share beside the largest single's and the sum of the singles', with a name for it.

    Only assessed, non-degenerate singles enter the sum, and which they were is in the result: a sum
    over three covariates and a union over four are not the same comparison and must not read as one.
    """
    singles = {
        n: exp[f"without_{n}"]["explains"]
        for n in COVARIATES
        if isinstance(exp.get(f"without_{n}"), dict) and exp[f"without_{n}"].get("explains") is not None
    }
    union_row = exp.get(f"without_{UNION}") or {}
    union = union_row.get("explains") if isinstance(union_row, dict) else None
    out: dict[str, Any] = {
        "union_explains": union,
        "union_assessed": bool(isinstance(union_row, dict) and union_row.get("assessed")),
        "union_degenerate": bool(isinstance(union_row, dict) and union_row.get("degenerate")),
        "singles_in_the_sum": sorted(singles),
        "singles_not_in_the_sum": [n for n in COVARIATES if n not in singles],
        "per_single": {n: singles[n] for n in sorted(singles)},
        "largest_single": round(max(singles.values()), 4) if singles else None,
        "sum_of_the_singles": round(sum(singles.values()), 4) if singles else None,
    }
    # a share above 1 is a ratio carried past 1, not an offset more than explained; it is named here
    # rather than left for a reader to infer from a percentage that looks like good news
    out["union_overshoots_the_offset"] = bool(union is not None and union > 1)
    if out["union_overshoots_the_offset"]:
        out["overshoot_reads"] = OVERSHOOT
    else:
        out["overshoot_reads"] = NO_OVERSHOOT if union is not None else UNASSESSED
    lo, hi = out["largest_single"], out["sum_of_the_singles"]
    if union is None or lo is None or hi is None:
        out["position_between_the_largest_and_the_sum"] = None
        out["reads"] = UNASSESSED
        return out
    pos = (union - lo) / (hi - lo) if hi > lo else None
    out["position_between_the_largest_and_the_sum"] = round(pos, 4) if pos is not None else None
    out["near_threshold"] = NEAR
    if union < lo:
        out["reads"] = BELOW_MAX
    elif union > hi:
        out["reads"] = ABOVE_SUM
    elif pos is None or pos <= NEAR:
        out["reads"] = NEAR_MAX
    elif pos >= 1 - NEAR:
        out["reads"] = NEAR_SUM
    else:
        out["reads"] = BETWEEN
    return out


# ================================================================================================
# the two-group comparison, where both arms are real groups of sequence
# ================================================================================================
def union_windows(
    rb: Rebuilt, feats: dict[str, array], tss: list[int], cap: int = WINDOW_CAP
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """100-base windows carrying their covariates, their tier, their coverage and the hit.

    The hit is "does this window carry a recurring event", the unit `panel_background.window_rows`
    chose so that a rate can be standardised as a share. A window is covariate-free when NOT ONE base
    of any of the four lies in it -- the same rule the exclusion arm applies, which removes bases and
    not windows -- and the count of windows that are partly covered is reported rather than folded in.
    Windows with no coverage value are kept here so `input_presence` can count what is present before
    a stratum is applied; they are dropped when the stratified comparison runs, and counted there.
    """
    panel = rb.panel
    per_block = informing(panel)
    ev_cell: Counter = Counter()
    for ev in rb.evs:
        if ev.minor >= MIN_RECURRING:
            ev_cell[ev.start // CELL] += 1
    cells = sorted(
        {
            c
            for i in range(len(panel.block_start))
            for c in range(panel.block_start[i] // CELL, (panel.block_end[i] - 1) // CELL + 1)
        }
    )
    step = max(1, -(-len(cells) // cap))
    rows: list[dict[str, Any]] = []
    thin = 0
    for c in cells[::step]:
        s, e = c * CELL, (c + 1) * CELL
        cov, aligned = coverage_over(panel, per_block, s, e)
        gc = panel.gc_fraction([(s, e)])
        if aligned < CELL / 2 or gc is None:
            thin += 1
            continue
        held = {n: span_bp(feats[n], s, e) for n in COVARIATES if n in feats}
        tier = next((t for t in TIERS if span_bp(feats[f"in_{t}"], s, e) >= CELL / 2), None)
        rows.append(
            {
                "start": s,
                "length": aligned,
                "gc": round(gc, 4),
                "rt": rt_stratum(rt_over(rb.rt, [(s, e)]), rb.edges) or 0,
                "tss": nearest(tss, s + CELL // 2),
                "coverage": round(cov, 4) if cov is not None else None,
                "recurring": bool(ev_cell.get(c, 0)),
                "tier": tier,
                # the annotation input: how many of the four the window carries, which is present on
                # every window by construction and is what `input_presence` is asked about
                "covariates_held": sum(1 for v in held.values() if v),
                "covariate_bases": sum(held.values()),
                "coding": span_bp(feats["cds_canonical"], s, e) > 0,
            }
        )
    return rows, {
        "step": step,
        "cells_in_panel_blocks": len(cells),
        "windows": len(rows),
        "windows_too_thin_or_without_gc": thin,
        "windows_without_a_coverage_value": sum(1 for r in rows if r["coverage"] is None),
        "windows_free_of_all_four": sum(1 for r in rows if not r["covariates_held"]),
        "windows_holding_at_least_one": sum(1 for r in rows if r["covariates_held"]),
        "windows_in_a_block": sum(1 for r in rows if r["tier"]),
        "coverage_measure": COVERAGE_MEASURE,
    }


INPUTS = {
    "the panel's alignment over the window": "coverage",
    "the local annotation of the window": "covariates_held",
}


def compare_groups(targets: list[dict[str, Any]], controls: list[dict[str, Any]]) -> dict[str, Any]:
    """One two-group comparison: `input_presence` first, then the three stratifications.

    Presence before the stratum, always, and in that order: a coverage stratum only bites a claim
    whose input had to be bought, so saying which BEFORE running it stops "conditioning changed
    nothing" from reading as a test passed when no test was administered.
    """
    out: dict[str, Any] = {
        "targets": len(targets),
        "controls": len(controls),
        "input_presence": input_presence(targets, controls, INPUTS),
    }
    t = [r for r in targets if r["coverage"] is not None]
    c = [r for r in controls if r["coverage"] is not None]
    out["targets_without_a_coverage_value"] = len(targets) - len(t)
    out["controls_without_a_coverage_value"] = len(controls) - len(c)
    out["targets_below_full_coverage"] = sum(1 for r in t if r["coverage"] < COVERAGE_FULL)
    out["controls_below_full_coverage"] = sum(1 for r in c if r["coverage"] < COVERAGE_FULL)
    if not t or not c:
        out["assessed"] = False
        out["why"] = "one of the two arms has no window with a coverage value"
        return out
    out["assessed"] = True
    for name, strata in LADDER.items():
        out[name] = standardised(t, c, strata, hit="recurring")
    return out


def comparisons(rows: list[dict[str, Any]]) -> dict[str, Any]:
    """The same question asked from both ends, with the covariate-free sequence as the control.

    The exclusion arm above rebuilds a background with the four covariates gone; the window comparison
    is that background made of real sequence, and it is run both ways round the same control:

    - *each tier against the covariate-free background*, which is what the tier ratio would read if
      the background were clean -- the union arm's own question, with neither arm a tier average
      because the control is sequence the tiers do not contain;
    - *the covariate-carrying background against the covariate-free background*, which is whether the
      four are quieter than what removing them leaves. If they are not, no union of them can explain
      an offset, however much sequence they cover.
    """
    pool = [r for r in rows if not r["coding"]]
    free = [r for r in pool if not r["covariates_held"] and not r["tier"]]
    held = [r for r in pool if r["covariates_held"] and not r["tier"]]
    out: dict[str, Any] = {
        "coding_windows_excluded": sum(1 for r in rows if r["coding"]),
        "control_windows_free_of_all_four": len(free),
        "background_windows_holding_at_least_one": len(held),
        "the_covariates_against_the_covariate_free_background": compare_groups(held, free),
    }
    for t in TIERS:
        targets = [r for r in pool if r["tier"] == t]
        out[f"{t}_against_the_covariate_free_background"] = (
            compare_groups(targets, free)
            if targets and free
            else {"assessed": False, "targets": len(targets), "why": "one of the two arms has no window"}
        )
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
    """Take all four covariates out of one chromosome's background at once, and say where that sits."""
    rb = rebuild(chrom, cache, results_dir)
    model = covariate_model(chrom, gff)
    tracks = local_tracks(chrom, tracks_dir, constraint_dir)
    ivs = covariate_intervals(model, tracks)
    present = [n for n in COVARIATES if ivs[n]]
    if not present:
        return {"chrom": chrom, "assessed": False, "why": NO_GENCODE if not model else NO_TRACK}
    feats, arms = union_arrays(rb, ivs)
    tss = model.get("tss", [])
    rows, dropped = background_bins(rb, feats, tss)
    names = [*COVARIATES, UNION]
    absent = [n for n in COVARIATES if not ivs[n]]
    groups: dict[str, Any] = {
        "background": composition(rows, names, absent),
        "background_free_of_all_four": composition([r for r in rows if r[UNION] < SPLIT_AT], names, absent),
    }
    for n in names:
        groups[f"kilobases_holding_{n}"] = composition([r for r in rows if r[n] >= SPLIT_AT], names, absent)
    for t in TIERS:
        groups[f"tier_{t}"] = block_composition(rb, t, feats, tss, names, absent)
    offset = decompose(rb, arms, rows, model, tracks)
    exp = {t: explained(offset, t) for t in TIERS}
    win, win_cost = union_windows(rb, feats, tss)
    comp = comparisons(win)
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
        "claim_available": CLAIM_WITH_COVERAGE if matched and all(matched) else CLAIM_WITHOUT_COVERAGE,
        # nothing below means anything if the rebuild is not the panel's own: this is the same check
        # `panel_background` makes, against the same committed ratios, before anything is read
        "reproduces_the_panel": reproduces_the_panel(rb, offset["as_the_panel_builds_it"]),
        "chromosome_bases": rb.length,
        # what could have been measured, before what was
        "eligible": {
            "chromosome_bases": rb.length,
            "kilobases_in_the_panel_alignment": len(rows) + sum(dropped.values()),
            "covariates_declared": list(COVARIATES),
            "covariates_with_bases_on_this_chromosome": present,
            "covariates_without_bases": {n: why_absent(n, model, tracks) for n in absent},
            "tiers_declared": list(TIERS),
            "tiers_with_blocks_on_this_chromosome": [t for t in TIERS if rb.tier_intervals(t)],
            "windows_in_the_panel_alignment": win_cost["cells_in_panel_blocks"],
            "genes_in_the_annotation": model.get("genes", 0),
            "coding_tss_read": len(tss),
        },
        "overlap": pairwise_overlap(ivs, rb.length),
        "groups": groups,
        "feature_splits": feature_splits(rows, names, absent),
        "offset": offset,
        "explained": exp,
        "where_the_union_sits": {t: where_it_sits(exp[t]) for t in TIERS},
        "comparisons": comp,
        "coverage": {
            "kilobases_in_panel_blocks": len(rows) + sum(dropped.values()),
            "kilobases_measured": len(rows),
            "kilobases_not_measured": dropped,
            "kilobases_free_of_all_four": sum(1 for r in rows if r[UNION] < SPLIT_AT),
            "tracks_absent": tracks.get("missing", []),
            "gencode_read": bool(model),
            "windows": win_cost,
        },
    }


__all__ = [
    "ARMS",
    "COVARIATES",
    "STRATIFICATIONS",
    "UNION",
    "analyse",
    "block_share_after",
    "covariate_intervals",
    "covariate_model",
    "pairwise_overlap",
    "where_it_sits",
]
