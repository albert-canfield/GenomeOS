#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""The species-tree null's run: area J's two readings, read once against the committed registration.

    uv run python scripts/species_null_run.py

Refuses to start unless `data/results/species_null_registration.json` is committed and unmodified,
the module and this script are the blobs it pins, the tree files and pruned topologies hash to what
it records, and every input result hashes to its recorded digest. Then it computes null A (library
pairs, Brownian deviations along the tree, Grafen and equal branch lengths) and null B (library
coherence, stochastic Dollo along the tree against birth-matched real genes), applies each decision
rule, fills the registered wording, and writes `species_null_result`. No request; about a minute.
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
from genomeos.knowledge import species_null as sn  # noqa: E402
from genomeos.knowledge.homology import library_members  # noqa: E402
from genomeos.results import load_result, save_result  # noqa: E402

NAME = "species_null_result"
ENTRY = "scripts/species_null_run.py"
REGISTRATION = "species_null_registration"
REGISTRATION_FILE = f"data/results/{REGISTRATION}.json"
OWN_CODE = (
    "genomeos/knowledge/species_null.py",
    "scripts/species_null_register.py",
    "scripts/species_null_run.py",
    "tests/test_species_null.py",
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
    pins = registration["frozen_code"]["blobs"]
    for path, blob in pins.items():
        now = committed(path)
        if not now["committed"] or now["head_blob"] != blob:
            raise SystemExit(f"{path} is not the registered code ({now['head_blob']} against {blob})")
    for entry in registration["result_manifest"]["inputs"]:
        digest = mf.sha256_of(Path(entry["path"]))[0]
        if digest != entry["sha256"]:
            raise SystemExit(f"{entry['path']} moved since registration ({digest} against {entry['sha256']})")
    for label, with_human in (("without_human_for_A", False), ("with_human_for_B", True)):
        root, _ = sn.species_tree(species, with_human=with_human)
        want = registration["tree_as_built"][label]["topology_sha256"]
        if sn.topology_digest(root) != want:
            raise SystemExit(f"the {label} topology moved since registration")
    return {
        "registration": state,
        "registration_sha256": hashlib.sha256((ROOT / REGISTRATION_FILE).read_bytes()).hexdigest(),
        "code_blobs_match": True,
        "inputs_match": True,
        "topologies_match": True,
    }


def log(stage: str, k: int, n: int) -> None:
    print(f"  {stage} {k + 1}/{n}", flush=True)


def run_a(presence: dict, origin: dict, members: dict) -> tuple[dict, dict, dict]:
    d = sn.pair_inputs(presence, origin)
    sets = sn.library_sets(members, d)
    root, _ = sn.species_tree(d.species, with_human=False)
    t = sn.arrays(root, d.species)
    out = {}
    for label, lengths in (("grafen", sn.grafen_lengths(t)), ("equal", sn.equal_lengths(t))):
        t0 = time.time()
        bm = sn.Brownian(sn.covariance(t, lengths))
        out[label] = sn.null_a(d, bm, sets, progress=log)
        out[label]["seconds"] = round(time.time() - t0, 1)
        print(f"null A ({label}): {out[label]['tested']} pairs in {out[label]['seconds']} s", flush=True)
    return out["grafen"], out["equal"], {"libraries": len(sets)}


def run_b(presence: dict, members: dict) -> tuple[list[dict], dict]:
    species, genes, bits = sn.presence_bits(presence)
    root, _ = sn.species_tree(species, with_human=True)
    t = sn.arrays(root, [*species, sn.HUMAN])
    t0 = time.time()
    dl = sn.Dollo(t, bits, human_column=len(species))
    sims = dl.simulate(sn.B_SIMS, sn.B_SIM_SEED)
    row = {g: i for i, g in enumerate(genes)}
    library_rows = {}
    for lib, syms in members.items():
        rows = np.array(sorted(row[s] for s in syms if s in row), dtype=np.int64)
        if len(rows) >= sn.MIN_MEMBERS:
            library_rows[lib] = rows
    rows_out = sn.null_b(dl, bits, sims, library_rows, progress=log)
    model = {
        "births_by_path_node": {t.names[int(u)]: int((dl.birth == j).sum()) for j, u in enumerate(dl.path)},
        "nodes_off_the_human_path": int((~dl.on_path).sum()),
        "loss_probability_range": [float(dl.p[~dl.on_path].min()), float(dl.p[~dl.on_path].max())],
        "theta_median": float(np.median(dl.theta)),
        "mean_simulated_species_per_gene": float(sims.sum(axis=2).mean()),
        "mean_observed_species_per_gene": float(bits.sum(axis=1).mean()),
        "seconds": round(time.time() - t0, 1),
    }
    return rows_out, model


def named_rows(grafen: dict, registration: dict) -> list[dict]:
    want = [
        tuple(p)
        for p in registration["readings_retested"]["A_library_pairs"]["named_rows_reported_individually"]
    ]
    by = {tuple(r["libraries"]): r for r in grafen["pairs"]}
    keep = ("r_exclusive", "r_age_expected", "tree_null_mean", "tree_null_sd", "p_tree", "q_tree")
    return [{"libraries": list(p), **{k: by.get(p, {}).get(k) for k in keep}} for p in want]


def result_manifest(registration: dict) -> dict[str, Any]:
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
            "exclusive to either side, exactly as §14 and §19 excluded them",
            "genes with no origin stratum, from the age expectation and the matched draws only",
        ],
        "partitions": "n/a: no training or held-out split; every figure is a null comparison on the "
        "full species and library sets",
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
    stored_19 = load_result("profiling_genes")

    grafen, equal, a_info = run_a(presence, origin, members)
    key_a, text_a, figures_a = sn.select_a(grafen, equal, stored_14)
    b_rows, model = run_b(presence, members)
    key_b, text_b, figures_b = sn.select_b(b_rows, stored_19)

    compact_equal = [
        {"libraries": r["libraries"], "p_tree": r.get("p_tree"), "q_tree": r.get("q_tree")}
        for r in equal["pairs"]
    ]
    out = {
        "result": NAME,
        "date": time.strftime("%Y-%m-%d"),
        "registration": gate_record,
        "A_library_pairs": {
            "outcome": key_a,
            "reading": text_a,
            "figures": figures_a,
            "named_rows": named_rows(grafen, registration),
            "grafen": {k: v for k, v in grafen.items() if k != "pairs"} | {"pairs": grafen["pairs"]},
            "equal_lengths": {k: v for k, v in equal.items() if k not in ("pairs", "halves")}
            | {"pairs": compact_equal, "halves": {k: v for k, v in equal["halves"].items() if k != "rows"}},
            "libraries": a_info["libraries"],
        },
        "B_library_coherence": {
            "outcome": key_b,
            "reading": text_b,
            "figures": figures_b,
            "libraries": b_rows,
            "model": model,
        },
    }
    out["result_manifest"] = result_manifest(registration)
    path = save_result(NAME, out)
    print(f"\nsaved {path}")
    print(f"\nA ({key_a}): {text_a}\n\nB ({key_b}): {text_b}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
