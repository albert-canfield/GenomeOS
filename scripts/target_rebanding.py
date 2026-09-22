# SPDX-License-Identifier: AGPL-3.0-or-later
"""Re-band the genome from a curve fitted on the population each band is quoted for.

    uv run python scripts/target_rebanding.py --strata   # the population table only, no sweep
    uv run python scripts/target_rebanding.py            # the population table and the re-banding

The published band table prices 593,765 swept targets in eight bands from a curve fitted on 245
CRISPRi pairs, 86.8% of them priced by 121 training and 14 held-out pairs. A band is a claim about
a rate, so this gives each population the band its own measured rate can carry -- and where the
rate's 95% interval spans more than one published band, or the population holds fewer than 30
pooled pairs, the targets are banded `not calibrated here` and carry the count and the interval
instead of a number. Registered in full as `target_calibration.PREREGISTERED_BANDS`.

Both axes named in the registration are computed: the predicted K562 deletion drop (the registered
one) and the track the sweep named the target on (the alternative, which replaces it if it is the
cleaner separator). The off-gate population -- the genes `targets.ElementResponses` has made
askable and the sweep does not band -- gets its own table from its own pairs.

Reads the cached benchmark tables (data/knowledge/crispri), the all-enhancer deletion table
(data/knowledge/alphagenome/all_elements), the per-element response cache
(data/knowledge/alphagenome/elements, 775 MB, one chromosome held at a time) and GENCODE v50. No
model request, no network. Saves data/results/target_rebanding.json.
"""

from __future__ import annotations

import argparse
import json
import time
from typing import Any

from genomeos.attribution import crispri
from genomeos.attribution import target_calibration as tc
from genomeos.attribution import targets as tg
from genomeos.results import save_result


def on_gate(rows: list[crispri.Pair]) -> list[crispri.Pair]:
    return [p for p in rows if p.features["top_target"]]


def off_gate(rows: list[crispri.Pair]) -> list[crispri.Pair]:
    return [p for p in rows if not p.features["top_target"]]


def tissue_axis(table: crispri.DeletionTable) -> Any:
    """The alternative axis for a benchmark pair: the track its gene was named on."""

    def key(p: crispri.Pair) -> str:
        el = tc.matched_element(table.overlapping(p.chrom, p.start, p.end), p.gene)
        if el is None:
            return "another track"
        for k in ("predicted", "predicted_coding"):
            q = el.get(k)
            if q and q.get("gene") == p.gene:
                return tc.CELL if (q.get("tissue") or "?") == tc.CELL else "another track"
        return "another track"

    return key


def populations(
    training: list[crispri.Pair], heldout: list[crispri.Pair], table: crispri.DeletionTable
) -> dict[str, Any]:
    """The band each population's own measured rate can carry, on both registered axes."""
    train = tc.scored(training)
    held = tc.scored([p for p in heldout if p.cell == tc.CELL])
    pooled_on = on_gate(train) + on_gate(held)
    pooled_off = off_gate(train) + off_gate(held)
    return {
        "on the gate, by predicted K562 drop (the registered axis)": tc.population_bands(
            pooled_on, tc.drop_stratum
        ),
        "on the gate, by the track the target was named on (the alternative axis)": tc.population_bands(
            pooled_on, tissue_axis(table)
        ),
        "off the gate, by predicted K562 drop (the population the sweep does not band)": tc.population_bands(
            pooled_off, tc.drop_stratum
        ),
        "pooled_on_gate_pairs": len(pooled_on),
        "pooled_off_gate_pairs": len(pooled_off),
    }


def show(name: str, t: dict[str, Any]) -> None:
    print(f"\n{name}")
    for s, v in t.items():
        if not isinstance(v, dict) or "band" not in v:
            continue
        if not v.get("pairs"):
            print(f"  {s:<24} no pairs -> {v['band']}")
            continue
        print(
            f"  {s:<24} n={v['pairs']:<5} observed={v['observed_rate']:.4f} "
            f"{v['ci95']} -> {v['band']}  ({v['why']})"
        )


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--strata", action="store_true", help="the population table only; no sweep")
    args = ap.parse_args()
    t0 = time.time()

    training, heldout = crispri.load(crispri.TRAINING), crispri.load(crispri.HELDOUT)
    table = crispri.DeletionTable()

    # the compact arm: exactly what the published curve was fitted through
    crispri.annotate(training, table)
    crispri.annotate(heldout, table)
    tc.add_features(training, table)
    tc.add_features(heldout, table)
    pops = populations(training, heldout, table)
    train = tc.scored(training)
    weights = {
        "target": tc.fit(on_gate(train), tc.TARGET_FEATURES),
        "sweep": tc.fit(train, tc.SWEEP_FEATURES),
    }

    for k, v in pops.items():
        if isinstance(v, dict):
            show(k, v)

    # the off-gate population read through the sweep's own window, where its magnitude is real
    print("\nreading the response cache for the off-gate magnitude (one chromosome at a time)...")
    tr2, he2 = crispri.load(crispri.TRAINING), crispri.load(crispri.HELDOUT)
    cache = crispri.ElementCache()
    crispri.annotate(tr2, table, cache)
    crispri.annotate(he2, table, cache)
    tc.add_features(tr2, table, responses=tg.ElementResponses())
    windowed = populations(tr2, he2, table)
    off_key = "off the gate, by predicted K562 drop (the population the sweep does not band)"
    pops["off the gate, through the sweep's own window (its magnitude is real here)"] = windowed[off_key]
    show("off the gate, through the sweep's own window", windowed[off_key])

    out: dict[str, Any] = {
        "preregistered": tc.PREREGISTERED_BANDS,
        "evidence": tc.EVIDENCE,
        "scope": tc.SCOPE,
        "what": (
            "the published band table kept as written, beside a band drawn from the population it "
            "is quoted for; a population whose measured rate has an interval spanning more than one "
            "published band, or fewer than 30 pooled pairs, is banded 'not calibrated here'"
        ),
        "min_pooled_for_a_band": tc.MIN_POOLED_FOR_A_BAND,
        "populations": pops,
        "weights": {k: [round(x, 6) for x in v] for k, v in weights.items()},
    }
    if args.strata:
        print(f"\npopulations only, nothing re-banded ({time.time() - t0:.0f} s)")
        print(json.dumps(out["populations"], indent=1)[:2000])
        return 0

    # both axes are re-banded, because the registration says the cleaner one replaces the other and
    # which is cleaner is not known until the intervals are read
    for axis, key in (
        ("drop", "on the gate, by predicted K562 drop (the registered axis)"),
        ("track", "on the gate, by the track the target was named on (the alternative axis)"),
    ):
        print(f"\nre-banding the sweep on the {axis} axis, chromosome by chromosome...")
        out[f"rebanding, {axis} axis"] = tc.reband(
            weights, pops[key], progress=lambda s: print(f"  {s}", flush=True), axis=axis
        )
        for label, v in out[f"rebanding, {axis} axis"]["genome_wide"].items():
            print(f"\n[{axis}] {label}: {v['targets']} targets")
            print(f"  published: {json.dumps(v['published_bands'])}")
            print(f"  re-banded: {json.dumps(v['rebanded'])}")
            print(f"  moves:     {json.dumps(v['moves'])}")
            print(f"  transitions: {json.dumps(v['transitions'])}")

    path = save_result("target_rebanding", out)
    print(f"\nsaved {path} ({time.time() - t0:.0f} s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
