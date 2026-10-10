#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""The tip-noise pair null's run: calibration on random pairs, then, only if it passes, §14's pairs.

    uv run python scripts/pair_null_run.py

Refuses to start unless `data/results/pair_null_registration.json` is committed and unmodified, the
code it pins is the code on disk, the tree topology hashes to what it records, and every input result
hashes to its recorded digest. Then it follows the registered rule in order and stops where a step
fails: reproduce §14's statistic; calibrate both tip-noise nulls on §22's random pairs (and recompute
§22's Brownian null on them, for the record); confirm the selected null on a second random set; only
then read every library pair once and the libraries' halves. It fills the registered wording and
writes `pair_null_result`. No request; about half an hour.
"""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from genomeos import manifest as mf  # noqa: E402
from genomeos.knowledge import pair_null as pn  # noqa: E402
from genomeos.knowledge import species_null as sn  # noqa: E402
from genomeos.knowledge.homology import library_members  # noqa: E402
from genomeos.results import load_result, save_result  # noqa: E402

NAME = "pair_null_result"
ENTRY = "scripts/pair_null_run.py"
REGISTRATION = "pair_null_registration"
REGISTRATION_FILE = f"data/results/{REGISTRATION}.json"
OWN_CODE = (
    "genomeos/knowledge/pair_null.py",
    "scripts/pair_null_register.py",
    "scripts/pair_null_run.py",
    "tests/test_pair_null.py",
)


def git(*args: str) -> str:
    return subprocess.run(
        ["git", *args], cwd=ROOT, capture_output=True, text=True, check=False
    ).stdout.strip()


def committed(path: str) -> dict[str, Any]:
    tracked = git("ls-files", "--error-unmatch", path) == path
    dirty = bool(git("status", "--porcelain", "--", path))
    head_blob = git("rev-parse", f"HEAD:{path}")
    return {
        "path": path,
        "tracked": tracked,
        "modified_or_untracked": dirty,
        "head_blob": head_blob or None,
        "committed": tracked and not dirty and bool(head_blob),
    }


def gate(registration: dict[str, Any], species: list[str]) -> dict[str, Any]:
    """Every registered pin, checked; a refusal names the first that moved."""
    state = committed(REGISTRATION_FILE)
    if not state["committed"]:
        raise SystemExit(f"the registration is not committed ({json.dumps(state)}): nothing is read")
    for path, blob in registration["frozen_code"]["blobs"].items():
        now = committed(path)
        if not now["committed"] or now["head_blob"] != blob:
            raise SystemExit(f"{path} is not the registered code ({now['head_blob']} against {blob})")
    for entry in registration["result_manifest"]["inputs"]:
        digest = mf.sha256_of(Path(entry["path"]))[0]
        if digest != entry["sha256"]:
            raise SystemExit(f"{entry['path']} moved since registration ({digest} against {entry['sha256']})")
    root, _ = sn.species_tree(species, with_human=False)
    if sn.topology_digest(root) != registration["tree_as_built"]["topology_sha256"]:
        raise SystemExit("the tree topology moved since registration")
    return {
        "registration": state,
        "registration_sha256": hashlib.sha256((ROOT / REGISTRATION_FILE).read_bytes()).hexdigest(),
        "code_blobs_match": True,
        "inputs_match": True,
        "topology_matches": True,
    }


def log(stage: str, k: int, n: int) -> None:
    print(f"  {stage} {k + 1}/{n}  {time.strftime('%H:%M:%S')}", flush=True)


def compact(cal: dict[str, Any]) -> dict[str, Any]:
    """A calibration block with its rows cut to the figures a reader checks."""
    rows = [{"for": r["for"], "p_pair": r["p_pair"], "r": r["r"], "lambda": r["lambda"]} for r in cal["rows"]]
    return {k: v for k, v in cal.items() if k != "rows"} | {"rows": rows}


def named_rows(read: dict[str, Any], registration: dict[str, Any]) -> list[dict[str, Any]]:
    want = [tuple(p) for p in registration["readings_retested"]["named_rows_reported_individually"]]
    by = {tuple(r["libraries"]): r for r in read["pairs"]}
    keep = ("r_exclusive", "r_age_expected", "lambda", "rate", "null_mean", "null_sd", "p_pair", "q_pair")
    return [{"libraries": list(p), **{k: by.get(p, {}).get(k) for k in keep}} for p in want]


def quartiles(values: list[float]) -> list[float] | None:
    return [float(v) for v in np.quantile(values, [0.25, 0.5, 0.75])] if values else None


def result_manifest(registration: dict[str, Any]) -> dict[str, Any]:
    reg = registration["result_manifest"]
    return {
        "sources": reg["sources"],
        "inputs": [
            *[mf.input_entry(Path(e["path"]), partition=e["partition"]) for e in reg["inputs"]],
            mf.input_entry(
                Path(REGISTRATION_FILE), partition="the committed registration this run reads against"
            ),
        ],
        "assembly": reg["assembly"],
        "coordinates": reg["coordinates"],
        "parameters": reg["parameters"],
        "exclusions": [
            "libraries with fewer than 8 members with a profile, and pairs with fewer than 8 members "
            "exclusive to either side, exactly as §14 and §22 excluded them",
            "genes with no origin stratum, from the age expectation, the binomial scale and the matched "
            "draws only",
            "every library pair and every library half, unless the selected null passed both "
            "calibration sets (the registered rule stops before them otherwise)",
        ],
        "partitions": "n/a: no training or held-out split; the calibration sets are random gene sets "
        "drawn for that purpose and nothing is fitted across pairs",
        "code_cleanliness": mf.code_cleanliness(ENTRY, OWN_CODE, ROOT),
    }


def main() -> int:
    mf.trace_begin()
    registration = load_result(REGISTRATION)
    if registration is None:
        raise SystemExit(f"no {REGISTRATION}: register first")
    presence = load_result("origin_presence_genome_wide")
    gate_record = gate(registration, presence["species"])
    print(f"gate: registration blob {gate_record['registration']['head_blob']}", flush=True)
    origin = load_result("origin_genome_wide")
    members = library_members() or {}
    stored_14 = load_result("profiling_genome_wide")["null"]
    s22 = load_result("species_null_result")
    bm_c1_hits_22 = s22["A_library_pairs"]["grafen"]["calibration"]["at_or_below_005"]

    d = sn.pair_inputs(presence, origin)
    sets = sn.library_sets(members, d)
    pairs = sn.testable_pairs(sets)
    out: dict[str, Any] = {"result": NAME, "date": time.strftime("%Y-%m-%d"), "registration": gate_record}
    timing: dict[str, float] = {}

    # rule 1: reproduce §14's statistic
    repro = pn.reproduce(d, pairs, stored_14)
    out["reproduction"] = repro
    cal1 = cal2 = read = hv = None
    chosen = None
    bm_hits = None
    if not repro["n_mismatch"] and repro["same_pairs"]:
        root, _ = sn.species_tree(d.species, with_human=False)
        t = sn.arrays(root, d.species)
        cov = sn.covariance(t, sn.grafen_lengths(t))
        rand1 = pn.random_pairs(d, pairs, pn.CAL1_SEED)

        # §22's Brownian null on its own random pairs, for the record
        t0 = time.time()
        bm = sn.Brownian(cov)
        bm_hits = pn.brownian_hits(d, bm, rand1)
        del bm
        timing["brownian_first_set"] = round(time.time() - t0, 1)
        print(f"Brownian on the first set: {bm_hits} of {len(rand1)} (§22: {bm_c1_hits_22})", flush=True)

        # rule 2: both tip-noise nulls on the first calibration set, random pairs only
        tn = pn.TipNoise(cov)
        cal1 = {}
        for model in pn.MODELS:
            t0 = time.time()
            cal1[model] = pn.calibrate(d, tn, rand1, model, progress=log)
            timing[f"first_set_{model}"] = round(time.time() - t0, 1)
            c = cal1[model]
            print(
                f"first set, {model}: {c['at_or_below_005']} of {c['pairs']} at p <= 0.05, "
                f"{c['q_at_or_below_005']} at q <= 0.05, passes {c['passes']}",
                flush=True,
            )
        chosen = pn.choose_null(cal1)

        # rule 3: the selected null on an independent random set
        if chosen is not None:
            t0 = time.time()
            rand2 = pn.random_pairs(d, pairs, pn.CAL2_SEED)
            cal2 = pn.calibrate(d, tn, rand2, chosen, progress=log)
            timing["second_set"] = round(time.time() - t0, 1)
            print(
                f"second set, {chosen}: {cal2['at_or_below_005']} of {cal2['pairs']} at p <= 0.05, "
                f"{cal2['q_at_or_below_005']} at q <= 0.05, passes {cal2['passes']}",
                flush=True,
            )

            # rules 4 and 5: only now is any library pair computed under the null
            if cal2["passes"]:
                t0 = time.time()
                read = pn.read_pairs(d, tn, chosen, pairs, progress=log)
                hv = pn.halves(d, tn, chosen, sets)
                timing["library_pairs_and_halves"] = round(time.time() - t0, 1)

    key, text, figures = pn.select(repro, stored_14, cal1, bm_hits, bm_c1_hits_22, chosen, cal2, read, hv)
    out["outcome"] = key
    out["reading"] = text
    out["figures"] = figures
    if cal1 is not None:
        out["first_calibration_set"] = {m: compact(c) for m, c in cal1.items()}
        out["first_calibration_set"]["brownian_section_22_recomputed"] = {
            "at_or_below_005": bm_hits,
            "section_22_recorded": bm_c1_hits_22,
        }
    if cal2 is not None:
        out["second_calibration_set"] = compact(cal2)
    if read is not None and hv is not None:
        out["library_pairs"] = read
        out["named_rows"] = named_rows(read, registration)
        out["halves"] = hv
        lib_lams = [x for r in read["pairs"] for x in (r.get("lambda") or [])]
        rand_lams = [x for r in cal1[chosen]["rows"] for x in (r["lambda"] or [])]
        out["descriptive"] = {
            "lambda_quartiles_library_sides": quartiles(lib_lams),
            "lambda_quartiles_random_sides_first_set": quartiles(rand_lams),
            "null_mean_range_over_pairs": [
                min(r["null_mean"] for r in read["pairs"] if r.get("null_mean") is not None),
                max(r["null_mean"] for r in read["pairs"] if r.get("null_mean") is not None),
            ],
        }
    out["seconds"] = timing
    out["result_manifest"] = result_manifest(registration)
    path = save_result(NAME, out)
    print(f"\nsaved {path}")
    print(f"\n{key}: {text}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
