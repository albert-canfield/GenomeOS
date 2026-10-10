#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""The mixed calibration's run: reproduce §23's read, then test every mixed pair once under its null.

    uv run python scripts/pair_mix_run.py

Refuses to start unless `data/results/pair_mix_registration.json` is committed and unmodified, the code
it pins is the code on disk, the tree topology hashes to what it records, and every input result hashes
to its recorded digest. Then it follows the registered rule in order: reproduce §14's statistic and
every library pair's count of null draws under §23's null; only if both reproduce, draw the mixed pairs
(each library against a random set of its partner's make-up) and test each once; apply the ceilings in
each direction. It fills the registered wording and writes `pair_mix_result`. No request; a few minutes.
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
from genomeos.knowledge import pair_mix as pm  # noqa: E402
from genomeos.knowledge import pair_null as pn  # noqa: E402
from genomeos.knowledge import species_null as sn  # noqa: E402
from genomeos.knowledge.homology import library_members  # noqa: E402
from genomeos.results import load_result, save_result  # noqa: E402

NAME = "pair_mix_result"
ENTRY = "scripts/pair_mix_run.py"
REGISTRATION = "pair_mix_registration"
REGISTRATION_FILE = f"data/results/{REGISTRATION}.json"
OWN_CODE = (
    "genomeos/knowledge/pair_mix.py",
    "scripts/pair_mix_register.py",
    "scripts/pair_mix_run.py",
    "tests/test_pair_mix.py",
)
CONSEQUENCE = {
    "calibrated": "stand_as_read",
    "not_calibrated": "withdrawn_to_no_reading",
    "not_reproduced": "nothing_decided",
}


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
    for k, v in registration["null_parameters"]["checked"].items():
        if getattr(pn, k) != v:
            raise SystemExit(f"pair_null.{k} is {getattr(pn, k)!r}, registered {v!r}")
    return {
        "registration": state,
        "registration_sha256": hashlib.sha256((ROOT / REGISTRATION_FILE).read_bytes()).hexdigest(),
        "code_blobs_match": True,
        "inputs_match": True,
        "topology_matches": True,
        "null_parameters_match": True,
    }


def log(stage: str, k: int, n: int) -> None:
    print(f"  {stage} {k + 1}/{n}  {time.strftime('%H:%M:%S')}", flush=True)


def compact(mix: dict[str, Any]) -> dict[str, Any]:
    """The mixed calibration with its rows cut to the figures a reader checks."""
    keep = (
        "libraries",
        "p_pair",
        "q_pair",
        "r",
        "r_age_expected",
        "lambda",
        "random_genes",
        "random_genes_short",
    )
    rows = {k: [{f: r.get(f) for f in keep} for r in mix["rows"][k]] for k in pm.DIRECTIONS}
    return {k: v for k, v in mix.items() if k != "rows"} | {"rows": rows}


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
            "exclusive to either side, exactly as §14, §22 and §23 excluded them",
            "genes with no origin stratum, from the age expectation, the binomial scale and the matched "
            "draws only",
            "the library side's own members, from the pools the random set beside it is drawn from",
            "every mixed pair, unless §14's statistic and §23's read reproduced first (the registered "
            "rule stops before them otherwise)",
        ],
        "partitions": "n/a: no training or held-out split; the random sets are drawn for the "
        "calibration and nothing is fitted across pairs",
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
    s23 = load_result("pair_null_result")
    stored_read = s23["library_pairs"]
    survivors_23 = {
        tuple(r["libraries"])
        for r in stored_read["pairs"]
        if r.get("q_pair") is not None and r["q_pair"] <= pn.Q_ALPHA
    }

    d = sn.pair_inputs(presence, origin)
    sets = sn.library_sets(members, d)
    pairs = sn.testable_pairs(sets)
    out: dict[str, Any] = {"result": NAME, "date": time.strftime("%Y-%m-%d"), "registration": gate_record}
    timing: dict[str, float] = {}

    # rule 1: reproduce §14's statistic and §23's read
    repro14 = pn.reproduce(d, pairs, stored_14)
    root, _ = sn.species_tree(d.species, with_human=False)
    t = sn.arrays(root, d.species)
    tn = pn.TipNoise(sn.covariance(t, sn.grafen_lengths(t)))
    t0 = time.time()
    repro23 = pm.reproduce_read(d, tn, pairs, stored_read, progress=log)
    timing["reproduce_read"] = round(time.time() - t0, 1)
    out["reproduction"] = {"section_14": repro14, "section_23": repro23}
    print(f"reproduction: §14 {json.dumps(repro14)}; §23 {json.dumps(repro23)}", flush=True)
    reproduced = (
        not repro14["n_mismatch"]
        and repro14["same_pairs"]
        and not repro23["n_mismatch_read"]
        and repro23["same_pairs"]
        and repro23["same_survivors"]
    )

    # rules 2 and 3: only now is any mixed pair drawn
    mix = None
    if reproduced:
        t0 = time.time()
        mixed = pm.mixed_pairs(d, pairs, pm.MIX_SEED)
        mix = pm.calibrate_mixed(d, tn, mixed, progress=log)
        timing["mixed_pairs"] = round(time.time() - t0, 1)
        for k in (*pm.DIRECTIONS, "pooled"):
            c = mix[k]
            print(
                f"{k}: {c['at_or_below_005']} of {c['pairs']} at p <= 0.05, "
                f"{c['q_at_or_below_005']} at q <= 0.05, within ceilings {c['passes']}",
                flush=True,
            )

    key, text, figures = pm.select_mix(repro14, repro23, stored_14, s23["figures"], mix)
    out["outcome"] = key
    out["reading"] = text
    out["section_23_consequence"] = CONSEQUENCE[key]
    out["section_23_consequence_registered"] = registration["mixed_calibration"]["what_happens_to_the_154"][
        key
    ]
    out["figures"] = figures
    if mix is not None:
        out["mixed_calibration"] = compact(mix)
        out["descriptive"] = {
            "per_library_as_the_real_side": pm.per_library(mix),
            "section_23_pairs_with_mixed_pairs_at_p_005": pm.among_read(mix, survivors_23),
            "lambda_library_sides_section_23": s23["descriptive"]["lambda_quartiles_library_sides"],
            "lambda_random_sides_section_23_first_set": s23["descriptive"][
                "lambda_quartiles_random_sides_first_set"
            ],
            "note": "descriptive only; none of these decides",
        }
    out["seconds"] = timing
    out["result_manifest"] = result_manifest(registration)
    path = save_result(NAME, out)
    print(f"\nsaved {path}")
    print(f"\n{key}: {text}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
