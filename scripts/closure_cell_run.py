#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""The M3 run: the cell-level readout, measured once, against the committed registration.

Four gates before a count is read, each one a refusal and not a warning:

1. `data/results/closure_cell_registration.json` must exist and must be committed -- the file's
   own `code_cleanliness` is not enough, so the gate asks git whether the registration's bytes are
   in HEAD and unmodified;
2. the frozen population on disk must still hash to the digest the registration pinned;
3. every streamed file's md5 must equal the md5 ENCODE publishes for it;
4. no file may have more than `MAX_MALFORMED_SHARE` of its rows outside the published layout.

Then the eight guide-quantification files are streamed once each, nothing of them reaches disk, and
the arms' readings are computed by the decision rule the registration fixed. The reading printed
and written is the registered wording with the numbers filled in, and no other word of it changes.
"""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import closure_cell as cc  # noqa: E402
from genomeos.results import load_result, save_result  # noqa: E402
from genomeos.therapeutics import depmap  # noqa: E402

NAME = "closure_cell_result"
ENTRY = "scripts/closure_cell_run.py"
OWN_CODE = (
    "genomeos/attribution/closure_cell.py",
    "scripts/closure_cell_inventory.py",
    "scripts/closure_cell_register.py",
    "scripts/closure_cell_run.py",
    "tests/test_closure_cell.py",
)
REGISTRATION = "closure_cell_registration"
REGISTRATION_FILE = "data/results/closure_cell_registration.json"
MAX_MALFORMED_SHARE = 0.001
ARM_A = "attributed_essential_neighbourhood_not"
ARM_B = "attributed_essential_neighbourhood_essential"
ARM_C = "attributed_not_neighbourhood_essential"
ARM_D = "attributed_not_neighbourhood_not"


def committed(path: str) -> dict[str, Any]:
    """Whether `path` is tracked, in HEAD and unmodified, as git answers and not as a file says."""

    def git(*args: str) -> str:
        return subprocess.run(
            ["git", *args], cwd=ROOT, capture_output=True, text=True, check=False
        ).stdout.strip()

    tracked = git("ls-files", "--error-unmatch", path) == path
    dirty = bool(git("status", "--porcelain", "--", path))
    head_blob = git("rev-parse", f"HEAD:{path}")
    return {
        "path": path,
        "tracked": tracked,
        "modified_or_untracked": dirty,
        "in_head": bool(head_blob),
        "head_blob": head_blob or None,
        "committed": tracked and not dirty and bool(head_blob),
    }


def gate(log) -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    """The registration, the frozen population and the gate's own record, or a refusal."""
    registration = load_result(REGISTRATION)
    if registration is None:
        raise SystemExit(f"no {REGISTRATION} result: the registration must be written and committed first")
    state = committed(REGISTRATION_FILE)
    if not state["committed"]:
        raise SystemExit(
            f"the registration is not committed ({json.dumps(state)}): a test registered after its "
            "outcome was read is not registered"
        )
    population = cc.load_population()
    digest = hashlib.sha256(cc.POPULATION.read_bytes()).hexdigest()
    pinned = registration["population"]["sha256"]
    if digest != pinned:
        raise SystemExit(f"the population has moved: {digest} on disk against {pinned} registered")
    print(f"gate: registration at {state['head_blob']}, population {digest}", file=log, flush=True)
    return (
        registration,
        population,
        {
            "registration": state,
            "registration_sha256": hashlib.sha256((ROOT / REGISTRATION_FILE).read_bytes()).hexdigest(),
            "population_sha256": digest,
            "population_matches_the_registration": True,
        },
    )


def stream_all(index: cc.GuideIndex, log) -> dict[str, dict[str, Any]]:
    """Every file of both experiments, streamed once, with the md5 and layout gates applied."""
    held: dict[str, dict[str, Any]] = {}
    for side, experiment in (("screen", cc.SCREEN), ("control", cc.CONTROL)):
        for replicate, accession in sorted(experiment["files"].items()):
            print(f"streaming {side} replicate {replicate}: {accession}", file=log, flush=True)
            out = cc.stream_guide_counts(accession, index, log=log)
            if not out["md5_matches_encode"]:
                raise SystemExit(
                    f"{accession}: streamed md5 {out['md5']} is not the md5 ENCODE publishes "
                    f"({cc.FILE_MD5.get(accession)})"
                )
            census = out["census"]
            share = census["malformed"] / census["rows"] if census["rows"] else 1.0
            if share > MAX_MALFORMED_SHARE:
                raise SystemExit(
                    f"{accession}: {census['malformed']} of {census['rows']} rows are outside the "
                    f"published layout ({share:.4%} > {MAX_MALFORMED_SHARE:.1%})"
                )
            out["side"] = side
            out["replicate"] = replicate
            held[f"{side}_{replicate}"] = out
    return held


def readouts(population: list[dict[str, Any]], held: dict[str, dict[str, Any]]) -> dict[str, Any]:
    """Each element's pooled readout, and the per-replicate-pair ones beside it."""
    screens = [v for v in held.values() if v["side"] == "screen"]
    controls = [v for v in held.values() if v["side"] == "control"]
    screen_total = sum(v["library_total"] for v in screens)
    control_total = sum(v["library_total"] for v in controls)
    rows: list[dict[str, Any]] = []
    dropped = {"no_guide_seen": 0, "under_min_guides": 0, "under_min_control_count": 0}
    for element in population:
        eid = element["id"]
        screen_count = sum(v["per_element"].get(eid, 0.0) for v in screens)
        control_count = sum(v["per_element"].get(eid, 0.0) for v in controls)
        guides_seen = max([v["guides_per_element"].get(eid, 0) for v in (*screens, *controls)], default=0)
        if guides_seen == 0:
            dropped["no_guide_seen"] += 1
            continue
        if guides_seen < cc.MIN_GUIDES:
            dropped["under_min_guides"] += 1
            continue
        if control_count < cc.MIN_CONTROL_COUNT:
            dropped["under_min_control_count"] += 1
            continue
        pooled = cc.log2_ratio(screen_count, screen_total, control_count, control_total)
        if pooled is None:
            dropped["no_guide_seen"] += 1
            continue
        pairs = []
        for replicate in sorted(cc.SCREEN["files"]):
            s = held.get(f"screen_{replicate}")
            c = held.get(f"control_{replicate}")
            if not s or not c:
                continue
            value = cc.log2_ratio(
                s["per_element"].get(eid, 0.0),
                s["library_total"],
                c["per_element"].get(eid, 0.0),
                c["library_total"],
            )
            if value is not None:
                pairs.append(value)
        rows.append(
            {
                **element,
                "screen_count": screen_count,
                "control_count": control_count,
                "guides_seen": guides_seen,
                "log2_ratio": pooled,
                "per_pair": pairs,
            }
        )
    return {
        "rows": rows,
        "dropped": dropped,
        "screen_library_total": screen_total,
        "control_library_total": control_total,
    }


def arm_table(rows: list[dict[str, Any]], cut: float | None) -> dict[str, Any]:
    """Every arm's size, readout and depletion rate, each with its denominator beside it."""
    by_arm: dict[str, list[dict[str, Any]]] = {}
    for row in rows:
        if row.get("arm"):
            by_arm.setdefault(row["arm"], []).append(row)
    table: dict[str, Any] = {}
    for arm, members in sorted(by_arm.items()):
        values = [m["log2_ratio"] for m in members]
        depleted = sum(1 for v in values if cut is not None and v <= cut)
        table[arm] = {
            "measured_elements": len(members),
            "mean_log2_ratio": cc.mean(values),
            "median_log2_ratio": cc.median(values),
            "elements_at_or_below_the_floor_s_5th_percentile": depleted,
            "rate_with_denominator": f"{depleted} of {len(members)}",
            "mean_guides_seen": cc.mean([float(m["guides_seen"]) for m in members]),
            "mean_pooled_control_count": cc.mean([float(m["control_count"]) for m in members]),
            "mean_prediction_magnitude": cc.mean([abs(float(m["effect_k562"])) for m in members]),
            "mean_distance_to_nearest_coding_tss": cc.mean(
                [float(m["nearest_coding_distance"]) for m in members]
            ),
        }
    return table


def comparisons(by_arm: dict[str, list[float]]) -> dict[str, Any]:
    """The registered comparisons, each one-sided for the first arm below the floor."""
    floor = by_arm.get(ARM_D, [])
    out: dict[str, Any] = {}
    for arm, values in sorted(by_arm.items()):
        if arm == ARM_D:
            continue
        out[f"{arm}_below_{ARM_D}"] = cc.mann_whitney_u(values, floor)
    return out


def select(table: dict[str, Any], tests: dict[str, Any], min_per_arm: int) -> tuple[str, dict[str, Any]]:
    """The outcome the registration's decision rule selects, applied in the order it fixes."""

    def n(arm: str) -> int:
        return int((table.get(arm) or {}).get("measured_elements") or 0)

    def med(arm: str) -> float | None:
        return (table.get(arm) or {}).get("median_log2_ratio")

    def p(arm: str) -> float | None:
        return (tests.get(f"{arm}_below_{ARM_D}") or {}).get("p_one_sided_a_below_b")

    def separates(arm: str) -> bool:
        m, f, pv = med(arm), med(ARM_D), p(arm)
        return m is not None and f is not None and m < f and pv is not None and pv < 0.05

    why: dict[str, Any] = {
        "sizes": {a: n(a) for a in (ARM_A, ARM_B, ARM_C, ARM_D, cc.POSITIVE_CONTROL)},
        "separates": {a: separates(a) for a in (ARM_A, ARM_B, ARM_C, cc.POSITIVE_CONTROL)},
        "min_per_arm": min_per_arm,
    }
    if min(n(ARM_A), n(ARM_D), n(cc.POSITIVE_CONTROL)) < min_per_arm:
        return "stopped_for_power", {**why, "rule": "1_stop"}
    if not separates(cc.POSITIVE_CONTROL):
        return "readout_not_shown_to_work", {**why, "rule": "2_instrument"}
    if separates(ARM_A):
        return "chain_supported", {**why, "rule": "3_chain"}
    if separates(ARM_C):
        return "proximity_explains_it", {**why, "rule": "4_proximity"}
    return "no_difference_detected", {**why, "rule": "5_negative"}


def fill(wording: str, table: dict[str, Any], tests: dict[str, Any], min_per_arm: int) -> str:
    """The registered sentence with the numbers in it and not one other word changed."""

    def n(arm: str) -> int:
        return int((table.get(arm) or {}).get("measured_elements") or 0)

    def med(arm: str) -> str:
        v = (table.get(arm) or {}).get("median_log2_ratio")
        return "not computed" if v is None else f"{v:+.3f}"

    def pv(arm: str) -> str:
        v = (tests.get(f"{arm}_below_{ARM_D}") or {}).get("p_one_sided_a_below_b")
        return "not computed" if v is None else (f"{v:.1e}" if v < 1e-3 else f"{v:.4f}")

    rate = (table.get(ARM_A) or {}).get("rate_with_denominator", "not computed")
    return wording.format(
        n_a=n(ARM_A),
        n_b=n(ARM_B),
        n_c=n(ARM_C),
        n_d=n(ARM_D),
        n_p=n(cc.POSITIVE_CONTROL),
        median_a=med(ARM_A),
        median_b=med(ARM_B),
        median_c=med(ARM_C),
        median_d=med(ARM_D),
        median_p=med(cc.POSITIVE_CONTROL),
        p_ad=pv(ARM_A),
        p_cd=pv(ARM_C),
        p_pd=pv(cc.POSITIVE_CONTROL),
        rate_a=rate.split(" of ")[0] if " of " in rate else rate,
        min_per_arm=min_per_arm,
    )


def result_manifest(read_now: list[str], streamed: dict[str, dict[str, Any]]) -> dict[str, Any]:
    groups = {
        "the committed registration and census of this test": [
            p for p in read_now if p.startswith("data/results/closure_cell")
        ],
        "the frozen population and the distilled K562 dependency column": [
            p for p in read_now if p.startswith("data/knowledge/closure_cell/")
        ],
    }
    inputs = [
        mf.files_entry(label, sorted(paths), partition="the registered inputs of this run")
        for label, paths in groups.items()
        if paths
    ]
    grouped = {p for paths in groups.values() for p in paths}
    inputs += [
        mf.input_entry(p, partition="read by the run")
        for p in sorted(p for p in set(read_now) - grouped if (ROOT / p).exists())
    ]
    inputs += [
        {
            "path": out["url"],
            "sha256": out["sha256"],
            "md5": out["md5"],
            "md5_matches_encode": out["md5_matches_encode"],
            "bytes": out["streamed_bytes"],
            "partition": f"{out['side']} replicate {out['replicate']}: streamed over HTTP and "
            "discarded row by row, never written to disk; the digest is over the bytes read",
        }
        for out in streamed.values()
    ]
    return {
        "sources": [
            {"accession": cc.SCREEN_SOURCE, "version": "ENCODE portal, released; " + cc.ENCODE_LICENCE},
            {
                "accession": "DepMap CRISPRGeneEffect (Chronos), K562 = " + cc.K562_MODEL,
                "version": f"{depmap.RELEASE}, doi:{depmap.DOI}; {depmap.LICENCE}",
            },
            {
                "accession": "the project's own frozen population and committed registration",
                "version": f"{REGISTRATION_FILE} as committed, digested in inputs",
            },
        ],
        "inputs": inputs,
        "input_count": len(inputs),
        "assembly": "GRCh38: the elements, the guide quantifications and the library design are all "
        "on one build, and no coordinate is lifted",
        "coordinates": {
            "base": 0,
            "interval": "half-open",
            "note": "a guide counts for an element when its 3 bp target interval overlaps the "
            "element's interval; a guide inside two population elements counts for both",
        },
        "parameters": {
            "min_guides": cc.MIN_GUIDES,
            "min_pooled_control_count": cc.MIN_CONTROL_COUNT,
            "pseudocount": cc.PSEUDOCOUNT,
            "depletion_percentile": cc.DEPLETION_PERCENTILE,
            "max_malformed_share": MAX_MALFORMED_SHARE,
            "alpha": 0.05,
        },
        "exclusions": [
            "every element outside the frozen population",
            "every element whose guides were not seen, or seen fewer than MIN_GUIDES times, or "
            "whose pooled control count is under MIN_CONTROL_COUNT",
            "the non-targeting guides from every arm: they are the assay's zero and not an arm",
        ],
        "partitions": {
            "arms": "fixed by the census before any count was read and pinned by the "
            "registration's population digest",
            "floor": ARM_D,
            "positive_control": cc.POSITIVE_CONTROL,
        },
        "code_cleanliness": mf.code_cleanliness(ENTRY, OWN_CODE, ROOT),
    }


def main() -> int:
    mf.trace_begin()
    started = time.time()
    log = sys.stdout
    registration, population, gate_record = gate(log)
    rows = population["elements"]
    index = cc.index_from_elements(rows)
    streamed = stream_all(index, log)
    measured = readouts(rows, streamed)
    floor_values = [
        r["log2_ratio"] for r in measured["rows"] if r.get("arm") == ARM_D and r["log2_ratio"] is not None
    ]
    cut = cc.percentile(floor_values, cc.DEPLETION_PERCENTILE)
    table = arm_table(measured["rows"], cut)
    by_arm: dict[str, list[float]] = {}
    for row in measured["rows"]:
        if row.get("arm"):
            by_arm.setdefault(row["arm"], []).append(row["log2_ratio"])
    tests = comparisons(by_arm)
    min_per_arm = int(registration["thresholds"]["min_per_arm"])
    outcome, why = select(table, tests, min_per_arm)
    reading = fill(registration["outcome_wordings"][outcome], table, tests, min_per_arm)
    non_targeting = {
        side: {
            "guides": sum(v["non_targeting"]["guides"] for v in streamed.values() if v["side"] == side),
            "count": sum(v["non_targeting"]["count"] for v in streamed.values() if v["side"] == side),
        }
        for side in ("screen", "control")
    }
    assay_zero = cc.log2_ratio(
        non_targeting["screen"]["count"],
        measured["screen_library_total"],
        non_targeting["control"]["count"],
        measured["control_library_total"],
    )
    per_pair: dict[str, Any] = {}
    for arm in sorted({r["arm"] for r in measured["rows"] if r.get("arm")}):
        members = [r for r in measured["rows"] if r.get("arm") == arm]
        means = [cc.mean(m["per_pair"]) for m in members if m["per_pair"]]
        per_pair[arm] = {
            "elements_with_at_least_one_pair": len(means),
            "median_of_the_per_pair_means": cc.median([m for m in means if m is not None]),
        }
    record: dict[str, Any] = {
        "question": registration["question"],
        "cell": cc.CELL,
        "registered": {
            "result": REGISTRATION,
            "sha256": gate_record["registration_sha256"],
            "committed_blob": gate_record["registration"]["head_blob"],
            "population_sha256": gate_record["population_sha256"],
        },
        "gate": gate_record,
        "streamed": {
            key: {
                "accession": v["accession"],
                "side": v["side"],
                "replicate": v["replicate"],
                "streamed_bytes": v["streamed_bytes"],
                "sha256": v["sha256"],
                "md5": v["md5"],
                "md5_matches_encode": v["md5_matches_encode"],
                "census": v["census"],
                "library_total": v["library_total"],
                "non_targeting": v["non_targeting"],
            }
            for key, v in sorted(streamed.items())
        },
        "population_elements": len(rows),
        "measured_elements": len(measured["rows"]),
        "not_measured": measured["dropped"],
        "depletion_cut": cut,
        "arms": table,
        "comparisons": tests,
        "per_replicate_pair": per_pair,
        "assay_zero_from_the_non_targeting_guides": {
            "log2_ratio": assay_zero,
            "guides_and_counts": non_targeting,
            "what_it_is": "the same ratio over the pooled non-targeting guides of the same library: "
            "the assay's own zero, not an arm and in no comparison",
        },
        "outcome": outcome,
        "decision": why,
        "reading": reading,
        "reading_is_the_registered_wording": "the sentence above is the registration's "
        f"outcome_wordings[{outcome!r}] with its numbers filled in and no other word changed",
        "cannot_establish": registration["cannot_establish"],
        "alphagenome_requests": 0,
        "money": "none",
        "seconds": round(time.time() - started, 1),
    }
    record["result_manifest"] = result_manifest(mf.traced_reads(), streamed)
    path = save_result(NAME, record)
    print(f"\nwrote {path}", file=log)
    print(f"outcome: {outcome}", file=log)
    print(reading, file=log)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
