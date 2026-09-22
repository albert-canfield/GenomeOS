# SPDX-License-Identifier: AGPL-3.0-or-later
"""The union axis: registered before it is fitted, and separated from the pairs that suggested it.

    uv run python scripts/union_axis.py --counts   # population sizes only; no label is read
    uv run python scripts/union_axis.py            # the measurement the registration asks for

On 2026-09-22 a lane crossed the two axes it had already measured -- the predicted K562 deletion
drop and the track the sweep named the element's target on -- found that their union bands 128 of
128 on-gate pairs at 1.0000 [0.9709, 1.0], and refused to adopt it, because an axis chosen after
seeing it win is priced by the search that found it and not by the interval it reports.

This lane adopts it the right way round. `target_calibration.PREREGISTERED_UNION` is fixed in the
module before a single number below is computed, and it separates three questions that the 128 of
128 runs together:

- the LEVEL (is the union stratum's rate inside 0.9-1) can only be read on the same 285 on-gate
  pairs that suggested the axis, so it is reported as a description and never as a test;
- the SEARCH (could a crossing this clean arise from two uninformative axes) is priced by a
  permutation of the labels over the whole family of partitions the crossing offers;
- the SEPARATION (does either arm carry information where it was never looked at) is tested on
  three populations that took no part in choosing the axis: the off-gate pairs read through the
  sweep's own window, the held-out pairs in the five cell lines that are not K562, and the union's
  own `neither` stratum read on the drop axis that was registered first.

Reads the cached benchmark tables (data/knowledge/crispri), the all-enhancer deletion table
(data/knowledge/alphagenome/all_elements), the per-element response cache
(data/knowledge/alphagenome/elements, 775 MB, one chromosome held at a time) and GENCODE v50. No
model request, no network. Saves data/results/union_axis.json.
"""

from __future__ import annotations

import argparse
import json
import random
import time
from collections import Counter
from typing import Any

from genomeos.attribution import crispri
from genomeos.attribution import target_calibration as tc
from genomeos.attribution import targets as tg
from genomeos.results import save_result

PERMUTATIONS = 20_000
SEED = 20260922


def on_gate(rows: list[crispri.Pair]) -> list[crispri.Pair]:
    return [p for p in rows if p.features["top_target"]]


def off_gate(rows: list[crispri.Pair]) -> list[crispri.Pair]:
    return [p for p in rows if not p.features["top_target"]]


def prepare(cache: crispri.ElementCache | None) -> tuple[list[crispri.Pair], list[crispri.Pair], Any]:
    """Both benchmark arms annotated the same way, with or without the sweep's own window.

    Both arms get `add_features`, which the 2026-09-22 driver called on the training arm only in
    its windowed pass; without it `scored` drops every held-out pair and the windowed table is
    training-only. That is corrected here and reported.
    """
    training, heldout = crispri.load(crispri.TRAINING), crispri.load(crispri.HELDOUT)
    table = crispri.DeletionTable()
    crispri.annotate(training, table, cache)
    crispri.annotate(heldout, table, cache)
    responses = tg.ElementResponses() if cache is not None else None
    tc.add_features(training, table, responses=responses)
    tc.add_features(heldout, table, responses=responses)
    return training, heldout, table


def k562(rows: list[crispri.Pair]) -> list[crispri.Pair]:
    return [p for p in rows if p.cell == tc.CELL]


def other_cells(rows: list[crispri.Pair]) -> list[crispri.Pair]:
    return [p for p in rows if p.cell != tc.CELL]


# ------------------------------------------------------------------------------------------
# Counts: the populations and their sizes, with no label read
# ------------------------------------------------------------------------------------------


def sizes(pairs: list[crispri.Pair], table: crispri.DeletionTable) -> dict[str, Any]:
    cells = tc.crossing_counts(pairs, table)
    return {
        "pairs": len(pairs),
        "by_crossing_cell": cells,
        "union": sum(v for k, v in cells.items() if k != "neither"),
        "neither": cells["neither"],
    }


def population_sizes(compact: tuple[Any, ...], windowed: tuple[Any, ...]) -> dict[str, Any]:
    training, heldout, table = compact
    wtrain, wheld, wtable = windowed
    tr, he = tc.scored(training), tc.scored(k562(heldout))
    wtr, whe = tc.scored(wtrain), tc.scored(k562(wheld))
    others = tc.scored(other_cells(heldout))
    return {
        "on the gate, K562, the 285 that suggested the axis": sizes(on_gate(tr) + on_gate(he), table),
        "off the gate, through the sweep's own window": sizes(off_gate(wtr) + off_gate(whe), wtable),
        "held out in the five cell lines that are not K562, on the gate": sizes(on_gate(others), table),
        "held-out cells": dict(sorted(Counter(p.cell for p in others).items())),
        "datasets of the 285": dict(sorted(Counter(p.dataset for p in on_gate(tr) + on_gate(he)).items())),
    }


# ------------------------------------------------------------------------------------------
# The three questions the registration separates
# ------------------------------------------------------------------------------------------


def described(pairs: list[crispri.Pair], table: crispri.DeletionTable, key: Any) -> dict[str, Any]:
    return tc.population_bands(pairs, lambda p: key(p, table))


def price_the_search(
    pairs: list[crispri.Pair], table: crispri.DeletionTable, permutations: int = PERMUTATIONS
) -> dict[str, Any]:
    """How often two uninformative axes hand the same search a stratum banded 0.9-1.

    The cells are held fixed and the labels are permuted, so the null is that neither axis carries
    information about regulation while both keep exactly the sizes they have. The search is the
    whole family it really was: every way of grouping the four cells of the crossing into a chosen
    stratum and a remainder. The statistic is whether the best such stratum clears the adoption
    rule -- at least `MIN_POOLED_FOR_A_BAND` pairs and a Wilson interval inside 0.9-1.
    """
    cells = [tc.pair_union_cell(p, table) for p in pairs]
    labels = [bool(p.regulated) for p in pairs]
    families = [tuple(sorted(g)) for n in range(1, len(tc.UNION_CELLS)) for g in _subsets(tc.UNION_CELLS, n)]
    # the membership of each grouping never changes under a permutation of the labels, so it is
    # resolved once and only the counting runs 20,000 times
    searched = [
        (fam, [i for i, c in enumerate(cells) if c in fam])
        for fam in families
        if sum(1 for c in cells if c in fam) >= tc.MIN_POOLED_FOR_A_BAND
    ]

    def best(ls: list[bool]) -> tuple[float, str]:
        top, name = -1.0, ""
        for fam, idx in searched:
            lo, _ = tc.wilson(sum(ls[i] for i in idx), len(idx))
            if lo > top:
                top, name = lo, " + ".join(fam)
        return top, name

    observed, chosen = best(labels)
    rng = random.Random(SEED)
    shuffled = list(labels)
    hits = 0
    for _ in range(permutations):
        rng.shuffle(shuffled)
        if best(shuffled)[0] >= 0.9:
            hits += 1
    return {
        "what": (
            "the labels permuted with the crossing held fixed; the statistic is the best Wilson "
            "lower bound any grouping of the four cells reaches under the adoption rule"
        ),
        "families_searched": len(families),
        "permutations": permutations,
        "best_grouping_observed": chosen,
        "best_lower_bound_observed": round(observed, 4),
        "permutations_reaching_0.9": hits,
        "p_value": round((hits + 1) / (permutations + 1), 5),
    }


def _subsets(items: tuple[str, ...], n: int) -> list[tuple[str, ...]]:
    if n == 0:
        return [()]
    if not items:
        return []
    head, rest = items[0], items[1:]
    return [(head, *s) for s in _subsets(rest, n - 1)] + _subsets(rest, n)


def heterogeneity(pairs: list[crispri.Pair], table: crispri.DeletionTable) -> dict[str, Any]:
    """The union's `neither` stratum read on the drop axis that was registered first.

    A stratum given one number has to be one population. `neither` is every pair below the drop
    threshold and off the K562 track, which the drop axis already reports as three strata with
    different rates, so this is the union axis's own falsifier applied to itself.
    """
    rest = [p for p in pairs if tc.pair_union_stratum(p, table) == "neither"]
    return {
        "pairs": len(rest),
        "by drop stratum": tc.population_bands(rest, tc.drop_stratum, min_pooled=1),
        "pooled": tc.observed_band(sum(p.regulated for p in rest), len(rest)),
    }


def eqtl_arm() -> dict[str, Any]:
    """The independent assay: GTEx cis-eQTLs, cached, no request."""
    from genomeos.attribution import eqtl
    from genomeos.genome import Annotation, default_gencode

    hits = eqtl.load_hits()
    symbols: dict[str, str] = {}
    for c in tc.CHROMS:
        path = default_gencode({c})
        if path is None:
            continue
        ann = Annotation.from_gff3(path, {c})
        symbols.update({gid.split(".")[0]: g.symbol for gid, g in ann.genes.items()})
    out = tc.eqtl_separation(hits, symbols, progress=lambda s: print(f"  {s}", flush=True))
    out["evidence"] = eqtl.EVIDENCE
    out["caveats"] = (
        "the distilled elements are a pre-selected index of 9,286, not a draw from the sweep; the "
        "label is association in bulk GTEx tissue, not a CRISPRi call; GTEx has no K562, so the "
        "track arm is asked a cross-context question. This prices direction, never a band's level."
    )
    return out


def verdict(out: dict[str, Any]) -> dict[str, Any]:
    """The registered adoption rule, applied without reinterpretation."""
    p = out["the search that found it, priced by permutation"]["p_value"]
    off = out["separation, off the gate through the sweep's own window"]["by crossing cell"]
    off_union = out["separation, off the gate through the sweep's own window"]["by union stratum"]
    het = out["the union's own neither stratum, on the drop axis registered first"]
    eq = out["separation, GTEx cis-eQTLs, an independent assay"]["any gene"]

    track_arm = tc.separated(off.get("track only", {}), off.get("neither", {}))
    union_arm = tc.separated(off_union.get("union", {}), off_union.get("neither", {}))
    homogeneous = all(
        v["band"] == het["pooled"]["band"]
        for v in het["by drop stratum"].values()
        if v.get("pairs", 0) >= tc.MIN_POOLED_FOR_A_BAND
    )
    a = p < 0.01
    b = bool(track_arm or union_arm or eq.get("separated"))
    return {
        "(a) the search is priced below 0.01": a,
        "(b) a population that took no part in the choice separates": b,
        "   off the gate, the track arm": track_arm,
        "   off the gate, the union stratum": union_arm,
        "   GTEx eQTLs, union against neither": eq.get("separated"),
        "(c) the neither stratum is homogeneous enough for one number": homogeneous,
        "adopted": a and b,
        "bands the neither stratum": a and b and homogeneous,
        "what that means": (
            "the union stratum is banded and `neither` is banded 'not calibrated here'"
            if (a and b and not homogeneous)
            else "nothing is banded from this axis"
            if not (a and b)
            else "both strata are banded"
        ),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--counts", action="store_true", help="population sizes only; no label read")
    ap.add_argument("--permutations", type=int, default=PERMUTATIONS)
    ap.add_argument("--no-sweep", action="store_true", help="skip the genome walk")
    args = ap.parse_args()
    t0 = time.time()

    compact = prepare(None)
    print(f"compact arm annotated ({time.time() - t0:.0f} s); reading the response cache...")
    windowed = prepare(crispri.ElementCache())
    print(f"windowed arm annotated ({time.time() - t0:.0f} s)")

    pops = population_sizes(compact, windowed)
    print("\nPOPULATION SIZES (no label read)")
    print(json.dumps(pops, indent=1))

    training, heldout, table = compact
    wtrain, wheld, wtable = windowed
    tr, he = tc.scored(training), tc.scored(k562(heldout))
    gate = on_gate(tr) + on_gate(he)
    weights = {
        "target": tc.fit(on_gate(tr), tc.TARGET_FEATURES),
        "sweep": tc.fit(tr, tc.SWEEP_FEATURES),
    }

    if args.counts:
        if not args.no_sweep:
            print("\nGENOME CENSUS (no weight fitted, no label read)")
            census = tc.crossing_census(progress=lambda s: print(f"  {s}", flush=True))
            print(json.dumps(census["genome_wide"], indent=1))
        print(f"\ncounts only, nothing measured ({time.time() - t0:.0f} s)")
        return 0

    out: dict[str, Any] = {
        "preregistered": tc.PREREGISTERED_UNION,
        "evidence": tc.EVIDENCE,
        "scope": tc.SCOPE,
        "what": (
            "the union of the two measured axes, registered before it was fitted; the level is a "
            "description of the same 285 pairs that suggested it, the search that found it is "
            "priced by permutation, and the separation is tested only where the axis was never "
            "looked at"
        ),
        "population_sizes": pops,
        "min_pooled_for_a_band": tc.MIN_POOLED_FOR_A_BAND,
        "union_drop_threshold": tc.UNION_DROP,
        "weights": {k: [round(x, 6) for x in v] for k, v in weights.items()},
    }

    # (1) the level: a description of the pairs that suggested the axis, never a test
    out["the level, described on the same 285 pairs"] = {
        "by crossing cell": described(gate, table, tc.pair_union_cell),
        "by union stratum": described(gate, table, tc.pair_union_stratum),
        "caveat": (
            "these are the pairs the axis was chosen on. The interval is the interval of the "
            "chosen subset and is biased upward by the choice; it prices no genome-wide target."
        ),
    }

    # (2) the search that found it
    out["the search that found it, priced by permutation"] = price_the_search(gate, table, args.permutations)

    # (3) separation, only where the axis was never looked at
    off = off_gate(tc.scored(wtrain)) + off_gate(tc.scored(k562(wheld)))
    others = on_gate(tc.scored(other_cells(heldout)))
    out["separation, off the gate through the sweep's own window"] = {
        "by crossing cell": described(off, wtable, tc.pair_union_cell),
        "by union stratum": described(off, wtable, tc.pair_union_stratum),
    }
    out["separation, held out in the cell lines that are not K562"] = {
        "by crossing cell": described(others, table, tc.pair_union_cell),
        "by union stratum": described(others, table, tc.pair_union_stratum),
        "transfer": (
            "the drop arm is read in each pair's own screen cell line, so this arm is a transfer "
            "claim as well as an independence claim; the track arm is the element's, not the cell's"
        ),
    }
    out["the union's own neither stratum, on the drop axis registered first"] = heterogeneity(gate, table)
    out["separation, GTEx cis-eQTLs, an independent assay"] = eqtl_arm()

    # the adoption rule, applied exactly as registered
    out["the adoption rule, as registered"] = verdict(out)

    for k, v in out.items():
        if isinstance(v, dict) and k.startswith(("the ", "separation")):
            print(f"\n{k}")
            print(json.dumps(v, indent=1)[:2500])

    rule = out["the adoption rule, as registered"]
    print("\nTHE ADOPTION RULE, AS REGISTERED")
    print(json.dumps(rule, indent=1))

    # the census is counted whatever the verdict: how much of the genome the axis can reach at all
    if not args.no_sweep:
        out["coverage"] = tc.crossing_census(progress=lambda s: print(f"  {s}", flush=True))

    if not rule["adopted"]:
        out["banding"] = (
            "nothing. Registered clause (b) was not met: no population that took no part in "
            "choosing this axis separates on it, so the axis is described and not adopted, and no "
            "swept target is re-banded from it."
        )
        print(f"\n{out['banding']}")
    elif not args.no_sweep:
        stratum_band = dict(described(gate, table, tc.pair_union_stratum))
        if not rule["bands the neither stratum"]:
            stratum_band["neither"] = {
                **stratum_band["neither"],
                "band": tc.NOT_CALIBRATED,
                "why": "registered clause (c): this stratum is not one population",
            }
        out["stratum_band_applied"] = stratum_band
        print("\nwalking the genome on the union axis...")
        out["rebanding, union axis"] = tc.reband(
            weights, stratum_band, progress=lambda s: print(f"  {s}", flush=True), axis="union"
        )
        for label, v in out["rebanding, union axis"]["genome_wide"].items():
            print(f"\n[union] {label}: {v['targets']} targets")
            print(f"  strata:    {json.dumps(v['by_stratum'])}")
            print(f"  re-banded: {json.dumps(v['rebanded'])}")
            print(f"  moves:     {json.dumps(v['moves'])}")

    path = save_result("union_axis", out)
    print(f"\nsaved {path} ({time.time() - t0:.0f} s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
