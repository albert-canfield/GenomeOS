# SPDX-License-Identifier: AGPL-3.0-or-later
"""The link-level direction test data/results/direction_link_registration.json fixed, taken after that
registration was committed (data/results/direction_link.json).

    uv run --frozen python scripts/direction_link.py
    uv run --frozen python scripts/direction_link.py --chroms chr21 --result direction_link_chr21

In the registered order:

1. **The eligibility join**, blind to every predicted value: one link per (attributed element, gene,
   cell) whose reciprocally overlapping cached CRISPRi pairs hold a significant effect, split into the
   decreases arm and the increases arm by the committed outcome labels, with the links whose pairs hold
   both counted on their own as contested.
2. **The answerability pass**: for each link the cached signed deletion value for its OWN gene on its
   OWN cell track. A cell the cache does not carry is unanswerable and nothing is substituted for it; an
   absent value and an exactly-zero value are counted and excluded from every denominator.
3. **The gate** on each arm separately against both imported floors.
4. If and only if the gate passes: both arms' rates apart, balanced accuracy, the difference with its
   baseline and excess, the cluster bootstrap intervals over independent loci, and the sign shuffle.
5. **The reading**, in the words `direction_link.ESTABLISHED` and `direction_link.NOT_DETECTED`
   registered before this run, and the section saying what the result cannot establish.

Nothing is changed. `measured.Layer.near`, `measured.reciprocal_overlap` and `targets.attributed` are
called exactly as committed; no compiled rule is read anywhere, because a rule-level test of these
increases is impossible (1 of 48 has a rule naming its gene, 0 have one gated on its cell). No
threshold, floor, predicate, seed or registered reading is moved, nothing is refitted, nothing is
downloaded and no model request is made.

**How the per-element cache is read.** The sweep's archives are streamed, not loaded: each
`<chrom>.json.gz` is decompressed in chunks and only the records of the elements this run needs are
materialised, so peak memory is one record rather than one chromosome. `json.load` on these archives is
the heavy read the project gates; this reader does not take it.
"""

from __future__ import annotations

import argparse
import gzip
import json
import re
import sys
import time
from collections import Counter, defaultdict
from datetime import date
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import crispri  # noqa: E402
from genomeos.attribution import direction_link as dl  # noqa: E402
from genomeos.attribution import measured as ms  # noqa: E402
from genomeos.attribution.targets import attributed  # noqa: E402
from genomeos.results import RESULTS_DIR, save_result  # noqa: E402

RESULT = "direction_link"
REGISTRATION = "direction_link_registration"
CHROMS = tuple(f"chr{c}" for c in [*range(1, 23), "X", "Y"])
CACHE_DIR = Path("data/knowledge/alphagenome/elements")

#: This lane's own files. Everything else uncommitted in this shared checkout belongs to another lane.
OWN_CODE = (
    "genomeos/attribution/direction_link.py",
    "scripts/direction_link_register.py",
    "scripts/direction_link.py",
    "tests/test_direction_link.py",
)

#: The same per-chromosome results `targets.attributed` reads, each digested by its own path.
INPUT_GLOBS = (
    "enhancer_targets_chr*.json",
    "enhancer_targets_all_chr*.json",
    "constrained_targets_chr*.json",
)

#: A top-level record of the per-element cache starts with its own id, which is what anchors the scan.
RECORD = re.compile(rb'"([^"]{1,64})":\s*\{"id":')


def _object_end(buf: bytes, start: int) -> int | None:
    """The index just past the balanced JSON object that begins at `buf[start] == '{'`, or None.

    Strings and their escapes are respected, so a brace inside a gene or tissue name cannot close a
    record early. None means the object is not yet complete in `buf`.
    """
    depth = 0
    in_string = False
    escaped = False
    for i in range(start, len(buf)):
        c = buf[i]
        if in_string:
            if escaped:
                escaped = False
            elif c == 0x5C:  # backslash
                escaped = True
            elif c == 0x22:  # quote
                in_string = False
            continue
        if c == 0x22:
            in_string = True
        elif c == 0x7B:  # {
            depth += 1
        elif c == 0x7D:  # }
            depth -= 1
            if depth == 0:
                return i + 1
    return None


def cached_records(path: Path, wanted: set[str], chunk: int = 1 << 23) -> dict[str, dict[str, Any]]:
    """The records of `wanted` out of one per-element cache archive, streamed.

    The archive is a single JSON object keyed by element id. It is decompressed in chunks and only a
    wanted record is ever accumulated, so a 700 MB chromosome costs one record of memory and not one
    chromosome. Nothing is written, nothing is cached and no value is read here: this returns the
    records, and the caller takes the one gene and one cell the registration allows it.
    """
    out: dict[str, dict[str, Any]] = {}
    if not wanted or not path.exists():
        return out
    remaining = set(wanted)
    buf = b""
    holding: str | None = None
    held = b""
    with gzip.open(path, "rb") as fh:
        while remaining:
            block = fh.read(chunk)
            if not block:
                break
            if holding is not None:
                held += block
                end = _object_end(held, 0)
                if end is None:
                    continue
                out[holding] = json.loads(held[:end])
                remaining.discard(holding)
                buf = held[end:]
                holding, held = None, b""
            else:
                buf += block
            while remaining:
                hit = next((m for m in RECORD.finditer(buf) if m.group(1).decode() in remaining), None)
                if hit is None:
                    break
                start = buf.index(b"{", hit.end() - len(b'{"id":'))
                end = _object_end(buf, start)
                name = hit.group(1).decode()
                if end is None:
                    holding, held = name, buf[start:]
                    buf = b""
                    break
                out[name] = json.loads(buf[start:end])
                remaining.discard(name)
                buf = buf[end:]
            if len(buf) > (1 << 20):  # keep only enough tail to span a record boundary
                buf = buf[-(1 << 20) :]
    return out


def overlap_index(
    elements: list[dict[str, Any]], layer: ms.Layer
) -> tuple[dict[int, list[str]], dict[str, tuple[int, int]]]:
    """Which attributed elements each cached pair is a measurement *of*, by the committed rule.

    The same two calls `measured.rows` makes, from the element side: `Layer.near` over
    `measured.REACH` and `measured.reciprocal_overlap` at `measured.RECIPROCAL_OVERLAP`. The answer is
    indexed by the pair's position in `layer.crispri`, so a pair covered by several elements is one
    pair. This is lane-increase's own `overlap_index`, unchanged.
    """
    where: dict[int, list[str]] = defaultdict(list)
    at = {id(p): i for i, p in enumerate(layer.crispri)}
    spans: dict[str, tuple[int, int]] = {}
    for e in elements:
        start, end = e["start"], e["end"]
        spans[e["id"]] = (start, end)
        for p in layer.near("crispri", start, end):
            if ms.reciprocal_overlap(start, end, p.start, p.end) >= ms.RECIPROCAL_OVERLAP:
                where[at[id(p)]].append(e["id"])
    return where, spans


def links_of(chrom: str) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Every eligible link of one chromosome, and the contested ones. Blind to every predicted value."""
    layer = ms.Layer.load(chrom)
    pairs = layer.crispri
    if not pairs:
        return [], []
    elements = attributed(chrom, RESULTS_DIR)
    at_element, spans = overlap_index(elements, layer)
    grouped: dict[tuple[str, str, str], list[ms.CrispriPair]] = defaultdict(list)
    for i, p in enumerate(pairs):
        if p.outcome not in (ms.DECREASE, ms.INCREASE):
            continue
        for e in at_element.get(i) or []:
            grouped[(e, p.gene, p.cell)].append(p)
    links: list[dict[str, Any]] = []
    contested: list[dict[str, Any]] = []
    for (e, gene, cell), gp in sorted(grouped.items()):
        down = [p for p in gp if p.outcome == ms.DECREASE]
        up = [p for p in gp if p.outcome == ms.INCREASE]
        start, end = spans[e]
        row = {
            "element": e,
            "chrom": chrom,
            "start": start,
            "end": end,
            "gene": gene,
            "cell": cell,
            "decreases": len(down),
            "increases": len(up),
            "splits": sorted({p.split for p in gp}),
            "datasets": sorted({p.dataset for p in gp}),
            "largest_abs_effect_size": round(max(abs(p.effect_size) for p in gp), 4),
        }
        arm = dl.arm_of(bool(down), bool(up))
        if arm is None:
            contested.append(row)
            continue
        links.append({**row, "arm": arm, "measured_sign": dl.MEASURED_SIGN[arm]})
    return links, contested


def answer(links: list[dict[str, Any]], cell: str) -> dict[str, Any]:
    """The cached signed value for each link's own gene on its own cell track. No substitution."""
    by_chrom: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for r in links:
        by_chrom[r["chrom"]].append(r)
    answered: list[dict[str, Any]] = []
    zero: list[dict[str, Any]] = []
    absent: list[dict[str, Any]] = []
    archives: list[str] = []
    for chrom in sorted(by_chrom):
        path = CACHE_DIR / f"{chrom}.json.gz"
        rows = by_chrom[chrom]
        records = cached_records(path, {r["element"] for r in rows})
        if path.exists():
            archives.append(path.as_posix())
        for r in rows:
            rec = records.get(r["element"])
            value = None
            for g in (rec or {}).get("genes") or []:
                if g.get("gene") == r["gene"]:
                    value = (g.get("by_cell") or {}).get(cell)
                    break
            if value is None:
                absent.append({**r, "why": "no cached value for this gene on this cell's track"})
                continue
            sign = dl.sign_of(value)
            row = {**r, "predicted_value": round(float(value), 6), "predicted_sign": sign}
            if sign is None:
                zero.append(row)
                continue
            answered.append({**row, "agrees": sign == r["measured_sign"]})
    return {"answered": answered, "zero": zero, "absent": absent, "archives": sorted(archives)}


def _counts(rows: list[dict[str, Any]], field: str) -> dict[str, int]:
    c: Counter[str] = Counter()
    for r in rows:
        v = r[field]
        for x in v if isinstance(v, list) else [v]:
            c[x] += 1
    return dict(sorted(c.items()))


def arm_payload(rows: list[dict[str, Any]], arm: str) -> dict[str, Any]:
    """One arm, on its own: the denominator, the rate, the breakdown, never pooled with the other."""
    return {
        "arm": arm,
        "measured_sign": dl.MEASURED_SIGN[arm],
        "answered_links": len(rows),
        "agreements": sum(1 for r in rows if r["agrees"]),
        "agreement_rate": None if not rows else round(dl.rate(rows), 4),
        "independent_loci": len(set(dl.loci_of(rows))) if rows else 0,
        "by_split": _counts(rows, "splits"),
        "by_chromosome": _counts(rows, "chrom"),
        "by_dataset": _counts(rows, "datasets"),
        "predicted_sign_breakdown": {
            "down": sum(1 for r in rows if r["predicted_sign"] == -1),
            "up": sum(1 for r in rows if r["predicted_sign"] == 1),
        },
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--chroms", default=",".join(CHROMS))
    ap.add_argument("--result", default=RESULT, help="a run over fewer chromosomes writes its own name")
    args = ap.parse_args()
    chroms = [c for c in args.chroms.split(",") if c]
    if args.result == RESULT and tuple(chroms) != CHROMS:
        raise SystemExit(f"{RESULT} is the genome-wide reading; a run over {len(chroms)} needs --result")

    started = time.time()
    links: list[dict[str, Any]] = []
    contested: list[dict[str, Any]] = []
    for chrom in chroms:
        got, bad = links_of(chrom)
        links.extend(got)
        contested.extend(bad)

    eligible = {
        "call": dl.ELIGIBILITY,
        "links": len(links),
        "by_arm": {a: sum(1 for r in links if r["arm"] == a) for a in dl.ARMS},
        "by_cell": _counts(links, "cell"),
        "by_arm_and_cell": {a: _counts([r for r in links if r["arm"] == a], "cell") for a in dl.ARMS},
        "contested": len(contested),
        "contest_rule": dl.CONTEST_RULE,
        "contested_rows": contested[:20],
    }

    cell = dl.PRIMARY_CELL
    in_cell = [r for r in links if r["cell"] == cell]
    unanswerable = [r for r in links if r["cell"] not in dl.CACHED_CELLS]
    other_cached = [r for r in links if r["cell"] != cell and r["cell"] in dl.CACHED_CELLS]
    got = answer(in_cell, cell)
    answered = got["answered"]

    answerability = {
        "cell_read": cell,
        "call": dl.CELL_AVAILABILITY,
        "cached_cell_tracks": list(dl.CACHED_CELLS),
        "cells_not_cached": list(dl.CELLS_NOT_CACHED),
        "links_in_the_cell_read": len(in_cell),
        "answered": len(answered),
        "predicted_zero_excluded": len(got["zero"]),
        "absent": len(got["absent"]),
        "breakdown_reconciles": len(answered) + len(got["zero"]) + len(got["absent"]) == len(in_cell),
        "sign_rule": dl.SIGN_CALL,
        "unanswerable_because_the_cache_carries_no_such_cell": {
            "links": len(unanswerable),
            "by_arm_and_cell": {
                a: _counts([r for r in unanswerable if r["arm"] == a], "cell") for a in dl.ARMS
            },
            "nothing_substituted": (
                "no other cell, gene, cache or pooled value stands in for any of these links; they are "
                "absent from the test rather than negative in it"
            ),
        },
        "links_in_another_cached_cell_not_read": {
            "links": len(other_cached),
            "by_arm_and_cell": {
                a: _counts([r for r in other_cached if r["arm"] == a], "cell") for a in dl.ARMS
            },
            "why": dl.PRIMARY_CELL_CALL,
        },
        "by_arm": {
            a: {
                "in_the_cell_read": sum(1 for r in in_cell if r["arm"] == a),
                "answered": sum(1 for r in answered if r["arm"] == a),
                "predicted_zero_excluded": sum(1 for r in got["zero"] if r["arm"] == a),
                "absent": sum(1 for r in got["absent"] if r["arm"] == a),
            }
            for a in dl.ARMS
        },
    }

    by_arm = {a: [r for r in answered if r["arm"] == a] for a in dl.ARMS}
    gates = {a: dl.gate(by_arm[a], a) for a in dl.ARMS}
    gate_met = all(g["meets_both_floors"] for g in gates.values())
    gate = {
        "call": dl.GATE_CALL,
        "per_arm": gates,
        "both_arms_meet_both_floors": gate_met,
        "reading": dl.GATE_PASS if gate_met else dl.GATE_NO_GO,
    }

    payload: dict[str, Any] = {
        "result": args.result,
        "date": date.today().isoformat(),
        "lane": "lane-direction",
        "registered_before_this_run": f"data/results/{REGISTRATION}.json",
        "level": "link level, not rule level",
        "why_not_rule_level": dl.registration()["why_not_rule_level"],
        "binding_wording_r2": dl.registration()["binding_wording_r2"],
        "follows_from": dl.registration()["follows_from"],
        "prior_direction_lanes": dl.PRIOR_DIRECTION,
        "chromosomes": chroms,
        "eligibility": eligible,
        "answerability": answerability,
        "gate": gate,
        "arms": {a: arm_payload(by_arm[a], a) for a in dl.ARMS},
        "arms_never_pooled": dl.registration()["arms_never_pooled"],
        "decrease_arm_is_not_model_coverage": dl.DECREASE_ARM_IS_NOT_MODEL_COVERAGE,
        "split_rule": dl.SPLIT_RULE,
        "not_balanced": dl.NOT_BALANCED,
        "floors_are_on_counts": dl.FLOORS_ARE_ON_COUNTS,
        "locus_convention": dl.registration()["locus_convention"],
        "cannot_establish": dl.CANNOT_ESTABLISH,
        "requests": 0,
        "money": "none: every input was already on disk",
        "no_requests": dl.NO_REQUESTS,
        "nothing_changed": dl.registration()["nothing_changed"],
    }

    if not gate_met:
        payload["comparison"] = {
            "taken": False,
            "why": dl.GATE_NO_GO,
            "short_by": {a: gates[a]["short_by"] for a in dl.ARMS},
        }
    else:
        d = dl.difference(by_arm)
        ba = dl.balanced_accuracy(answered)
        base = dl.baseline(answered)
        boot = dl.bootstrap(answered)
        control = dl.shuffle_control(answered, d)
        read = dl.reading(boot["ci95"]["balanced_accuracy"], control["p_one_sided"])
        payload["comparison"] = {
            "taken": True,
            "statistic_the_reading_is_taken_on": dl.TEST_STATISTIC,
            "balanced_accuracy": {
                "rate": round(ba, 4),
                "chance": dl.CHANCE,
                "ci95_over_independent_loci": boot["ci95"]["balanced_accuracy"],
                "why_this_statistic": dl.BALANCED_CALL,
            },
            "difference_between_the_arms": {
                "statistic": dl.PRIMARY_STATISTIC,
                "value": round(d, 4),
                "ci95_over_independent_loci": boot["ci95"]["difference"],
                "marginal_predicted_down_rate": round(dl.marginal_down(answered), 4),
                "baseline": round(base, 4),
                "baseline_call": dl.BASELINE,
                "excess_over_the_baseline": round(dl.excess(answered), 4),
                "excess_ci95_over_independent_loci": boot["ci95"]["excess"],
                "excess_call": dl.EXCESS_CALL,
                "why_this_is_not_the_reading": dl.COLLAPSE,
            },
            "per_arm_ci95_over_independent_loci": {a: boot["ci95"][a] for a in dl.ARMS},
            "bootstrap": {k: v for k, v in boot.items() if k != "ci95"},
            "control": control,
            "reading": read,
        }

    paths: set[Path] = {Path(f"data/results/{REGISTRATION}.json")}
    for name in ms.CRISPRI_FILES:
        p = crispri.KNOWLEDGE / name
        if p.exists():
            paths.add(p)
    for glob in INPUT_GLOBS:
        for chrom in chroms:
            paths.update(RESULTS_DIR.glob(glob.replace("chr*", chrom)))
    inputs = [mf.input_entry(p, partition=None) for p in sorted(paths, key=lambda q: q.as_posix())]
    if got["archives"]:
        inputs.append(
            mf.files_entry(
                "per_element_response_cache",
                got["archives"],
                partition="the finished deletion sweep's own per-element cache, read for a sign only",
            )
        )

    payload["seconds"] = round(time.time() - started, 1)
    payload["result_manifest"] = {
        "sources": [
            {
                "accession": "ENCODE CRISPRi enhancer-gene benchmark "
                "(EngreitzLab/CRISPR_comparison, Gschwind et al.)",
                "version": "the two benchmark tables cached under data/knowledge, digested in inputs",
            },
            {
                "accession": "the finished genome-wide AlphaGenome deletion sweep's per-element "
                "response cache, written with threshold=0.0",
                "version": "the per-chromosome archives this run streamed, digested as a group",
            },
            {
                "accession": "the attributed elements of the chromosomes this run covers",
                "version": "re-enumerated in this run from the results on disk, not read from a file",
            },
        ],
        "inputs": inputs,
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "eligibility_predicate": dl.ELIGIBILITY,
            "measured_sign": dl.MEASURED_SIGN_CALL,
            "predicted_direction": dl.PREDICTED_CALL,
            "sign_rule_zero_and_absent": dl.SIGN_CALL,
            "primary_statistic": dl.TEST_STATISTIC,
            "chance": dl.CHANCE,
            "difference_between_the_arms": dl.PRIMARY_STATISTIC,
            "why_the_difference_is_not_the_reading": dl.COLLAPSE,
            "baseline": dl.BASELINE,
            "independent_locus_rule": dl.LOCUS_RULE,
            "independent_locus_span": dl.LOCUS_SPAN,
            "locus_floor": dl.LOCUS_FLOOR,
            "link_floor": dl.POSITIVE_FLOOR,
            "reciprocal_overlap": ms.RECIPROCAL_OVERLAP,
            "reach": ms.REACH,
            "bootstrap_unit": "independent locus",
            "bootstrap_draws": dl.DRAWS,
            "bootstrap_seed": dl.SEED,
            "shuffles": dl.SHUFFLES,
            "shuffle_seed": dl.SHUFFLE_SEED,
            "alpha_one_sided": dl.CONTROL_ALPHA,
            "cells_read": [cell],
            "cached_cell_tracks": list(dl.CACHED_CELLS),
            "cells_not_cached": list(dl.CELLS_NOT_CACHED),
            "chromosomes": chroms,
            "aggregate_function": "mean of a 0/1 agreement per answered link, per arm",
            "fill_value": "none: an absent or exactly-zero predicted value is excluded and counted",
        },
        "exclusions": [
            dl.CONTEST_RULE,
            dl.SIGN_CALL,
            dl.CELL_AVAILABILITY,
            "no pooled raw agreement rate over the two arms is reported as a result",
            "no pair is excluded for its magnitude: no effect-size bar is applied anywhere",
        ],
        "partitions": {
            "decreases": "links holding a significant decrease and no significant increase; sign -1",
            "increases": "links holding a significant increase and no significant decrease; sign +1",
            "contested": dl.CONTEST_RULE,
            "training_and_heldout": dl.SPLIT_RULE,
        },
        "code_cleanliness": mf.code_cleanliness(Path(__file__).resolve(), OWN_CODE),
    }

    path = save_result(args.result, payload)
    print(f"{args.result}: {path}")
    print(f"  eligible links: {eligible['links']} ({eligible['by_arm']}), contested {eligible['contested']}")
    print(
        f"  cell read: {cell}; answered {answerability['answered']} of {len(in_cell)}, "
        f"zero {answerability['predicted_zero_excluded']}, absent {answerability['absent']}"
    )
    no_cell = answerability["unanswerable_because_the_cache_carries_no_such_cell"]["links"]
    print(f"  unanswerable because the cache carries no such cell: {no_cell}")
    for a in dl.ARMS:
        p = payload["arms"][a]
        print(
            f"  {a}: {p['agreements']}/{p['answered_links']} = {p['agreement_rate']}, "
            f"{p['independent_loci']} independent loci -> "
            f"{'floors met' if gates[a]['meets_both_floors'] else '; '.join(gates[a]['short_by'])}"
        )
    if payload["comparison"]["taken"]:
        c = payload["comparison"]
        print(
            f"  balanced accuracy {c['balanced_accuracy']['rate']} "
            f"{c['balanced_accuracy']['ci95_over_independent_loci']} against {dl.CHANCE}"
        )
        print(
            f"  difference {c['difference_between_the_arms']['value']} "
            f"{c['difference_between_the_arms']['ci95_over_independent_loci']}, "
            f"baseline {c['difference_between_the_arms']['baseline']}, "
            f"excess {c['difference_between_the_arms']['excess_over_the_baseline']}"
        )
        print(
            f"  control p {c['control']['p_one_sided']} (balanced accuracy), "
            f"{c['control']['p_one_sided_difference']} (difference)"
        )
        print(f"  established: {c['reading']['established']}")
    else:
        print("  comparison not taken: gate no-go")


if __name__ == "__main__":
    main()
