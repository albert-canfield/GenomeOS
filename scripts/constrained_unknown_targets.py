# SPDX-License-Identifier: AGPL-3.0-or-later
"""What the completed sweep says about the real unknown: the constrained-unknown blocks that are not copies.

    uv run python scripts/constrained_unknown_targets.py [--chroms chr21,chr22]

Milestone 1.3 asks for the constrained-unknown blocks attributed to a gene and a tissue. The sweep
finished on 2026-09-16, so every ENCODE element inside a node has been deleted in AlphaGenome on every
chromosome, and the answer is now arithmetic over tables already on disk: no model request.

For every UNKNOWN block (the organiser's join of tier, human axis and copy flag) this counts the
elements that lie inside it, how many of them move a gene, and what they name. It reports the
constrained-unknown tier with the copies taken out -- the organiser's "real unknown" -- split by the
case the two axes make, and it reports it **against the other tiers**, because a rate with no control
is what this project has had to withdraw three times:

- per element, the share that moves a gene, which does not depend on how long a block is;
- per block, the share carrying at least one element that moves a gene, inside length deciles, since a
  longer block holds more elements for no biological reason.

The result is `constrained_unknown_targets`. It names genes but claims nothing about them. Since
2026-09-27 it carries its own control (`matched_random_control`): the same block-level question asked
of 50 random windows of each block's length outside the organiser's blocks, lifted from
`attribution.unknown_scoring` and checked against its chr21 figures to the digit before any genome-wide
figure is drawn. The 87% it used to quote was the locus benchmark's, a different instrument.
"""

from __future__ import annotations

import argparse
import json
import time
from collections import defaultdict
from pathlib import Path
from statistics import median
from typing import Any

from genomeos.attribution import organise
from genomeos.results import load_result, save_result

ELEMENTS = Path("data/knowledge/alphagenome/all_elements")
MIN_LOG2 = 0.1  # the sweep's own threshold for "this deletion moved a gene"
DECILES = 10
RANDOM_DRAWS = 50  # windows per block, as `unknown_scoring.matched_random_windows`
SEED = 20260913  # the same seed, so the lift can be checked against chr21 to the digit
CHANCE_BAND = 5.0  # points either side of the random-window rate that read as "at chance"

# Registered 2026-09-27, before the control below was run on any chromosome but chr21 (where the
# unknown-scoring lane measured it). The 87% this script used to quote was the locus benchmark's,
# a different instrument; this is the script's own.
CONTROL_REGISTRATION = {
    "registered": "2026-09-27",
    "lift_check": (
        "the lifted function, handed chr21's 446 UNKNOWN blocks and unknown_scoring's predicate (names a "
        "coding gene), must reproduce 7f7c8c1 to the digit: 22,200 windows drawn, 2 undrawable, 16,706 "
        "carrying an element, 10,908 naming; 285 blocks carrying an element, 140 naming. Any difference "
        "means the lift is wrong and no genome-wide figure is reported"
    ),
    "reproduce_exactly": (
        "every field this script already wrote on 2026-09-16 (real_unknown 331 of 531 blocks with an "
        "element move, 0.6234; moves per element 0.2476; the other tiers; by case; the length deciles; "
        "331 named blocks with their targets): the tables are unchanged and the window reading cannot "
        "change a yes or no (see scripts/onetarget_consumers.py)"
    ),
    "primary": (
        "real unknown (882 blocks, copies out), the script's own question: does a block carrying a "
        "scored element carry one that moves a gene. Against 50 windows per block of the block's length, "
        "drawn inside the chromosome's scored-element span and rejected on overlap with any organiser "
        "block, pooled over the 24 chromosomes"
    ),
    "readings": {
        "below_chance": (
            f"block rate more than {CHANCE_BAND} points below the random rate: a lead at a real-unknown "
            "block is rarer than at random sequence of the same length outside unknown space. Registered as "
            "a live possibility and as the expected one (chr21, coding targets: 49.1% against 65.3%), not "
            "as a failure: it says what the 331 leads are worth"
        ),
        "at_chance": f"within {CHANCE_BAND} points: a lead carries no block-level information",
        "above_chance": f"more than {CHANCE_BAND} points above: the first above-chance reading of the leads",
    },
    "secondary": [
        "the same control with 'names a coding gene' as the question, the chr21 instrument's",
        "the neutral tier as a second target set, since it is the script's own comparison",
        "the window reading's genes at the bar per block, beside the head gene per element",
    ],
    "caveat": (
        "the random windows avoid all UNKNOWN space, so 'chance' is sequence outside it of the same "
        "length, which is gene-richer than the unknown; the tier comparison inside length deciles stays "
        "the within-unknown control"
    ),
}


def elements_of(chrom: str) -> list[dict[str, Any]]:
    """The chromosome's scored elements, or nothing if the sweep never wrote it."""
    p = ELEMENTS / f"{chrom}.json"
    if not p.exists():
        return []
    return json.loads(p.read_text())


def moves(element: dict[str, Any]) -> bool:
    """Did deleting this element move a gene, by the sweep's own threshold?"""
    pred = element.get("predicted") or {}
    log2 = pred.get("log2_fold_change")
    return bool(pred.get("gene")) and log2 is not None and abs(log2) >= MIN_LOG2


def inside(block: dict[str, Any], els: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Elements whose midpoint falls in the block: an element is not split between two blocks."""
    lo, hi = block["start"], block["end"]
    return [e for e in els if lo <= (e["start"] + e["end"]) // 2 < hi]


def read_chromosome(chrom: str, els: list[dict[str, Any]] | None = None) -> list[dict[str, Any]]:
    """Every UNKNOWN block of one chromosome with the sweep's answers attached."""
    els = elements_of(chrom) if els is None else els
    if not els:
        return []
    els.sort(key=lambda e: e["start"])
    rows = []
    for block in organise.blocks(chrom):
        got = inside(block, els)
        moving = [e for e in got if moves(e)]
        rows.append(
            {
                "chrom": chrom,
                "block": f"{chrom}:{block['start']}-{block['end']}",
                "start": block["start"],
                "end": block["end"],
                "length": block["length"],
                "tier": block.get("tier"),
                "case": block.get("case"),
                "copy": bool(block.get("copy")),
                "elements": len(got),
                "moving": len(moving),
                "targets": sorted({(e["predicted"] or {}).get("gene") for e in moving} - {None}),
                "cells": sorted({(e["predicted"] or {}).get("tissue") for e in moving} - {None}),
                "strongest": max(
                    (abs((e["predicted"] or {}).get("log2_fold_change", 0)) for e in moving), default=0
                ),
            }
        )
    return rows


def rate(rows: list[dict[str, Any]]) -> dict[str, Any]:
    """Per-element and per-block rates for a set of blocks, with the counts they came from."""
    elements = sum(r["elements"] for r in rows)
    moving = sum(r["moving"] for r in rows)
    with_any = [r for r in rows if r["elements"]]
    return {
        "blocks": len(rows),
        "blocks_with_an_element": len(with_any),
        "elements": elements,
        "moving": moving,
        "moves_per_element": round(moving / elements, 4) if elements else None,
        "blocks_with_a_moving_element": sum(1 for r in rows if r["moving"]),
        "share_of_blocks_with_an_element_that_move": (
            round(sum(1 for r in with_any if r["moving"]) / len(with_any), 4) if with_any else None
        ),
        "median_length_kb": round(median([r["length"] for r in rows]) / 1000, 1) if rows else None,
    }


def by_length_decile(target: list[dict[str, Any]], control: list[dict[str, Any]]) -> dict[str, Any]:
    """The per-block comparison inside length deciles, so block size cannot produce the difference."""
    pool = sorted(r["length"] for r in control + target)
    if not pool:
        return {}
    cuts = [pool[int(len(pool) * i / DECILES)] for i in range(1, DECILES)]

    def decile(length: int) -> int:
        return sum(length >= c for c in cuts)

    buckets: dict[int, dict[str, list]] = defaultdict(lambda: {"target": [], "control": []})
    for r in target:
        buckets[decile(r["length"])]["target"].append(r)
    for r in control:
        buckets[decile(r["length"])]["control"].append(r)
    rows, weighted, weight = {}, 0.0, 0
    for d in sorted(buckets):
        t, c = buckets[d]["target"], buckets[d]["control"]
        if not t or not c:
            continue
        t_rate = sum(1 for r in t if r["moving"]) / len(t)
        c_rate = sum(1 for r in c if r["moving"]) / len(c)
        rows[str(d)] = {
            "target_blocks": len(t),
            "control_blocks": len(c),
            "target": round(t_rate, 4),
            "control": round(c_rate, 4),
            "difference": round(t_rate - c_rate, 4),
        }
        weighted += (t_rate - c_rate) * len(t)
        weight += len(t)
    return {
        "per_decile": rows,
        "deciles_compared": len(rows),
        "difference_weighted_by_target_blocks": round(weighted / weight, 4) if weight else None,
    }


# ------------------------------------------------ the control: length-matched random windows (2026-09-27)


def matched_random_windows(
    targets: list[dict[str, Any]],
    exclude: list[dict[str, Any]],
    rows: list[dict[str, Any]],
    predicates: dict[str, Any],
    seed: int = SEED,
    draws: int = RANDOM_DRAWS,
) -> dict[str, Any]:
    """The block-level question asked of random windows of each target block's length, on one chromosome.

    Lifted from `attribution.unknown_scoring.matched_random_windows` (7f7c8c1) with its draw unchanged,
    so the same seed, blocks and rows give the same windows to the digit: for every target block, in
    start order, `draws` windows of its length placed uniformly inside the span of the rows' midpoints
    and rejected on overlap with any block of `exclude`; 4,000 tries at most, and a block longer than
    the span is undrawable and counted. An element belongs to the window its midpoint falls in. What it
    generalises is the question: `predicates` maps a name to a test on one element, and a block or a
    window says yes to it when any element inside does, so several questions share one set of windows.

    Rates are over blocks and windows carrying at least one element; a window with none was never asked.
    Counts are returned rather than rates so chromosomes pool by adding.
    """
    import bisect
    import random

    out: dict[str, Any] = {"blocks": len(targets), "windows_drawn": 0}
    if not targets or not rows:
        out.update(
            undrawable_blocks=0,
            undrawable_lengths=[],
            windows_carrying_an_element=0,
            blocks_carrying_an_element=0,
            **{name: {"blocks": 0, "windows": 0} for name in predicates},
        )
        return out
    mids = sorted(((r["start"] + r["end"]) // 2, i) for i, r in enumerate(rows))
    keys = [m for m, _ in mids]
    lo, hi = keys[0], keys[-1]
    unknown = sorted((b["start"], b["end"]) for b in exclude)
    u_starts = [s for s, _ in unknown]

    def overlaps(s: int, t: int) -> bool:
        i = max(0, bisect.bisect_right(u_starts, s) - 1)
        while i < len(unknown) and unknown[i][0] < t:
            if unknown[i][1] > s:
                return True
            i += 1
        return False

    def inside_window(s: int, t: int) -> list[dict[str, Any]]:
        return [rows[i] for _, i in mids[bisect.bisect_left(keys, s) : bisect.bisect_left(keys, t)]]

    rng = random.Random(seed)
    drawn = carry = 0
    undrawable: list[int] = []
    yes = dict.fromkeys(predicates, 0)
    for b in sorted(targets, key=lambda b: b["start"]):
        length, got, tries = b["length"], 0, 0
        while got < draws and tries < 4000:
            tries += 1
            top = hi - length
            if top <= lo:
                break
            s = rng.randrange(lo, top)
            if overlaps(s, s + length):
                continue
            got += 1
            drawn += 1
            sel = inside_window(s, s + length)
            if not sel:
                continue
            carry += 1
            for name, test in predicates.items():
                yes[name] += any(test(r) for r in sel)
        if not got:
            undrawable.append(length)
    with_element = [sel for sel in (inside_window(b["start"], b["end"]) for b in targets) if sel]
    out.update(
        span_bp=hi - lo,
        windows_drawn=drawn,
        undrawable_blocks=len(undrawable),
        undrawable_lengths=sorted(undrawable, reverse=True)[:5],
        windows_carrying_an_element=carry,
        blocks_carrying_an_element=len(with_element),
    )
    for name, test in predicates.items():
        out[name] = {
            "blocks": sum(1 for sel in with_element if any(test(r) for r in sel)),
            "windows": yes[name],
        }
    return out


def pooled(per_chrom: dict[str, dict[str, Any]], names: list[str]) -> dict[str, Any]:
    """The per-chromosome counts added up, then the rates and the registered reading."""
    tot = {
        k: sum(v.get(k, 0) for v in per_chrom.values())
        for k in ("blocks", "windows_drawn", "undrawable_blocks", "windows_carrying_an_element")
    }
    tot["blocks_carrying_an_element"] = sum(
        v.get("blocks_carrying_an_element", 0) for v in per_chrom.values()
    )
    b, w = tot["blocks_carrying_an_element"], tot["windows_carrying_an_element"]
    for name in names:
        yb = sum(v[name]["blocks"] for v in per_chrom.values())
        yw = sum(v[name]["windows"] for v in per_chrom.values())
        diff = round(100 * (yb / b - yw / w), 2) if b and w else None
        tot[name] = {
            "blocks_yes": yb,
            "block_rate": round(yb / b, 4) if b else None,
            "windows_yes": yw,
            "random_window_rate": round(yw / w, 4) if w else None,
            "difference_in_points": diff,
            "reading": (
                None
                if diff is None
                else "below_chance"
                if diff < -CHANCE_BAND
                else "above_chance"
                if diff > CHANCE_BAND
                else "at_chance"
            ),
        }
    return tot


def moves_coding(element: dict[str, Any]) -> bool:
    """Does the deletion name a coding gene -- the chr21 instrument's question (`names_coding`)?"""
    return bool(element.get("predicted_coding"))


def lift_check() -> dict[str, Any]:
    """The lifted function against 7f7c8c1 on chr21, with that module's own blocks and rows."""
    from genomeos.attribution import unknown_scoring as us

    blocks = us.unknown_blocks("chr21")
    rows = us.scored_elements("chr21")
    got = matched_random_windows(blocks, blocks, rows, {"names_coding": moves_coding}, seed=us.SEED)
    want = (load_result("unknown_scoring_chr21") or {})["matched_random_windows"]
    pairs = {
        "blocks": (got["blocks"], want["blocks"]),
        "windows_drawn": (got["windows_drawn"], want["windows_drawn"]),
        "undrawable_blocks": (got["undrawable_blocks"], want["undrawable_blocks"]),
        "undrawable_lengths": (got["undrawable_lengths"], want["undrawable_lengths"]),
        "span_bp": (got["span_bp"], want["span_bp"]),
        "windows_carrying_an_element": (
            got["windows_carrying_an_element"],
            want["windows_carrying_an_element"],
        ),
        "blocks_carrying_an_element": (
            got["blocks_carrying_an_element"],
            want["unknown_blocks_carrying_an_element"],
        ),
        "blocks_naming": (
            got["names_coding"]["blocks"],
            want["from_the_named_gene"]["unknown_blocks_naming_a_coding_gene"],
        ),
        "windows_naming": (
            got["names_coding"]["windows"],
            want["from_the_named_gene"]["random_windows_naming_a_coding_gene"],
        ),
    }
    return {
        "against": "unknown_scoring_chr21.matched_random_windows (7f7c8c1)",
        "fields": {k: {"lifted": a, "stored": b} for k, (a, b) in pairs.items()},
        "passed": all(a == b for a, b in pairs.values()),
    }


def genes_at_bar_per_block(
    chrom: str, rows: list[dict[str, Any]], els: list[dict[str, Any]]
) -> dict[str, int]:
    """Secondary: over the real-unknown blocks with a mover, the head genes against every gene at the bar
    in the cache's window. One chromosome's archive, dropped before the next."""
    from genomeos.attribution.targets import ElementResponses

    r = ElementResponses()
    c = {
        "blocks": 0,
        "moving_elements": 0,
        "head_genes": 0,
        "genes_at_bar": 0,
        "not_cached": 0,
        "yes_no_disagree": 0,
    }
    for row in rows:
        if not row["moving"]:
            continue
        s, t = (int(x) for x in row["block"].split(":")[1].split("-"))
        got = [e for e in els if s <= (e["start"] + e["end"]) // 2 < t]
        heads = {(e.get("predicted") or {}).get("gene") for e in got if moves(e)} - {None}
        genes = set(heads)
        for e in got:
            w = r.at_bar(chrom, e["id"], MIN_LOG2)
            if w is None:
                c["not_cached"] += 1
                continue
            c["yes_no_disagree"] += bool(w) != moves(e)
            genes.update(g for g, _ in w)
        c["blocks"] += 1
        c["moving_elements"] += row["moving"]
        c["head_genes"] += len(heads)
        c["genes_at_bar"] += len(genes)
    return c


def collect(chroms: list[str], control: bool = True, window: bool = True) -> dict[str, Any]:
    t0 = time.time()
    rows: list[dict[str, Any]] = []
    # the lift is checked before any genome-wide figure is drawn; if it fails, none is reported
    lift = lift_check() if control else None
    control = bool(lift and lift["passed"])
    questions = {"moves_a_gene": moves, "names_a_coding_gene": moves_coding}
    per_chrom: dict[str, dict[str, dict[str, Any]]] = {"real_unknown": {}, "neutral": {}}
    at_bar: dict[str, int] = {}
    for chrom in chroms:
        els = elements_of(chrom)
        els.sort(key=lambda e: e["start"])
        got = read_chromosome(chrom, els)
        rows.extend(got)
        print(f"{chrom}: {len(got)} blocks, {sum(r['elements'] for r in got)} elements inside", flush=True)
        if control and got:
            sets = {
                "real_unknown": [r for r in got if r["tier"] == "constrained_unknown" and not r["copy"]],
                "neutral": [r for r in got if r["tier"] == "neutral"],
            }
            for name, target in sets.items():
                per_chrom[name][chrom] = matched_random_windows(target, got, els, questions)
        if window and got:
            real = [r for r in got if r["tier"] == "constrained_unknown" and not r["copy"]]
            for k, v in genes_at_bar_per_block(chrom, real, els).items():
                at_bar[k] = at_bar.get(k, 0) + v

    real_unknown = [r for r in rows if r["tier"] == "constrained_unknown" and not r["copy"]]
    copies = [r for r in rows if r["tier"] == "constrained_unknown" and r["copy"]]
    others = {t: [r for r in rows if r["tier"] == t] for t in ("neutral", "fossil", "regulatory")}
    by_case = {
        case: rate([r for r in real_unknown if r["case"] == case])
        for case in sorted({r["case"] for r in real_unknown if r["case"]})
    }
    named = [r for r in real_unknown if r["targets"]]
    return {
        "result": "constrained_unknown_targets",
        "chromosomes": chroms,
        "min_log2": MIN_LOG2,
        "real_unknown": rate(real_unknown),
        "constrained_unknown_copies": rate(copies),
        "other_tiers": {t: rate(v) for t, v in others.items()},
        "by_case": by_case,
        "against_neutral_in_length_deciles": by_length_decile(real_unknown, others["neutral"]),
        "blocks_with_a_named_target": len(named),
        "named": sorted(
            (
                {
                    "block": r["block"],
                    "case": r["case"],
                    "elements": r["elements"],
                    "moving": r["moving"],
                    "targets": r["targets"][:6],
                    "cells": r["cells"][:4],
                    "strongest_abs_log2": round(r["strongest"], 3),
                }
                for r in named
            ),
            key=lambda r: -r["strongest_abs_log2"],
        ),
        **(
            {
                "matched_random_control": {
                    "registration": CONTROL_REGISTRATION,
                    "lift_check": lift,
                    "draws_per_block": RANDOM_DRAWS,
                    "seed_per_chromosome": SEED,
                    "real_unknown": pooled(per_chrom["real_unknown"], list(questions)),
                    "neutral": pooled(per_chrom["neutral"], list(questions)),
                    "per_chromosome_real_unknown": {
                        c: {
                            "blocks_carrying_an_element": v["blocks_carrying_an_element"],
                            "windows_carrying_an_element": v["windows_carrying_an_element"],
                            "undrawable_blocks": v["undrawable_blocks"],
                            **{q: v[q] for q in questions},
                        }
                        for c, v in per_chrom["real_unknown"].items()
                    },
                }
            }
            if lift is not None
            else {}
        ),
        **({"genes_at_bar_real_unknown": at_bar} if window else {}),
        "reading": reading(per_chrom, list(questions))
        if control
        else (
            "no random-window control on this run"
            if lift is None
            else "the lift check failed: no control reported"
        ),
        "seconds": round(time.time() - t0, 1),
    }


def reading(per_chrom: dict[str, dict[str, dict[str, Any]]], questions: list[str]) -> str:
    """The registered reading of the primary, in words, replacing the borrowed 87%."""
    p = pooled(per_chrom["real_unknown"], questions)["moves_a_gene"]
    words = {
        "below_chance": "below chance: a lead at a real-unknown block is rarer than at random sequence of "
        "the same length outside unknown space",
        "at_chance": "at chance: a lead carries no block-level information",
        "above_chance": "above chance: the first above-chance reading of the leads",
    }[p["reading"]]
    return (
        f"{p['blocks_yes']} real-unknown blocks carrying a scored element carry one that moves a gene "
        f"({100 * p['block_rate']:.1f}%), against {100 * p['random_window_rate']:.1f}% of "
        "length-matched random windows outside the organiser's blocks "
        f"({p['difference_in_points']:+.2f} points), which reads "
        f"{words}. A named target here is a lead, not a finding; the random windows sample gene-richer "
        "sequence than the unknown, so the tier comparison inside length deciles stays the "
        "within-unknown control"
    )


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--chroms", default="", help="comma-separated; default every chromosome with a table")
    ap.add_argument(
        "--no-save",
        action="store_true",
        help="print without writing the result: a run on a few chromosomes must not overwrite the genome's",
    )
    ap.add_argument("--no-control", action="store_true", help="skip the matched random windows")
    ap.add_argument("--no-window", action="store_true", help="skip the cache's genes at the bar")
    args = ap.parse_args(argv)
    chroms = (
        args.chroms.split(",")
        if args.chroms
        else sorted((p.stem for p in ELEMENTS.glob("chr*.json")), key=lambda c: (len(c), c))
    )
    out = collect(chroms, control=not args.no_control, window=not args.no_window)
    before = load_result(out["result"]) or {}
    if before and not args.chroms:
        # the registered reproduction: every field the 2026-09-16 run wrote, bar its reading and timing
        same = {k: before[k] == out.get(k) for k in before if k not in ("reading", "seconds", "date")}
        out["reproduced_from_2026_09_16"] = {"fields": same, "all": all(same.values())}
        print(f"reproduced 2026-09-16 fields: {out['reproduced_from_2026_09_16']['all']}", flush=True)
    if args.no_save:
        print("\nnot saved (--no-save)")
    else:
        print(f"\nsaved {save_result(out['result'], out)}")
    ru, ne = out["real_unknown"], out["other_tiers"]["neutral"]
    print(f"real unknown: {ru['blocks']} blocks, {ru['elements']} elements, moves {ru['moves_per_element']}")
    print(f"neutral tier: {ne['blocks']} blocks, {ne['elements']} elements, moves {ne['moves_per_element']}")
    deciles = out["against_neutral_in_length_deciles"]
    print(f"per block inside length deciles: {deciles.get('difference_weighted_by_target_blocks')}")
    print(f"blocks with a named target: {out['blocks_with_a_named_target']} of {ru['blocks']}")
    ctl = out.get("matched_random_control")
    if ctl:
        print(f"lift check passed: {ctl['lift_check']['passed']}")
        for tier in ("real_unknown", "neutral"):
            print(tier, json.dumps({k: v for k, v in ctl[tier].items()}))
    print(out["reading"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
