# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""Census of every caller of `NetworkRuntime.run(clamp=...)` and the figures that depend on it.

lane-h found on 2026-09-28 (8e705f7) that the clamp was re-applied only between whole steps, so
inside the Runge-Kutta stages a clamped species moved: a clamped mRNA was transcribed and its protein
translated off it, a clamped signal decayed. This script lists the callers (found by grep over
genomeos/, scripts/, tests/ and the web client, stated below by hand and checked against the source
at run time) and measures each dependent figure on the runtime as it stands, under a label:

    uv run python scripts/grn_clamp_census.py before   # the runtime before the fix
    uv run python scripts/grn_clamp_census.py after    # the runtime after it

Each run keeps the other label's measurements, so the result holds old and new side by side.
Writes `data/results/grn_clamp_census.json`. No network, no AlphaGenome.
"""

from __future__ import annotations

import json
import re
import subprocess
import sys
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tests"))

from genomeos import manifest as mf  # noqa: E402
from genomeos.lang import parse  # noqa: E402
from genomeos.results import RESULTS_DIR, save_result  # noqa: E402
from genomeos.runtime.gastrulation import run_gastrulation  # noqa: E402
from genomeos.runtime.grn import NetworkRuntime  # noqa: E402
from genomeos.web.server import Api  # noqa: E402

NAME = "grn_clamp_census"
# one file per label after the first: the before record stays byte for byte as committed, and
# `after` writes data/results/grn_clamp_census_after.json beside it
if __name__ == "__main__" and len(sys.argv) > 1 and sys.argv[1] != "before":
    NAME = f"{NAME}_{sys.argv[1]}"
ROOT = Path(__file__).resolve().parents[1]

CALLERS: list[dict[str, Any]] = [
    {
        "caller": "genomeos/runtime/gastrulation.py run_gastrulation",
        "clamp": "{'Nodal': nodal_max * exp(-x / decay_length)}, one value per cell, non-zero",
        "kind": "external signal held at a non-zero level; Nodal is a declared protein (gene NODAL max 0)",
        "reached_from": ["genomeos develop gastrulation (cli)", "web /api/develop model=gastrulation"],
        "dependents": [
            "tests/test_gastrulation.py::test_three_layers_in_order_and_plausible_proportions "
            "(cells=60, hours=30, dt=0.05): order endoderm/mesoderm/ectoderm, every proportion > 0.05, "
            "|proportion - expected| < 0.25",
            "docs/PROGRESS.md row 4.3 and docs/LESSONS.md line 58: qualitative (order, within expectation)",
        ],
    },
    {
        "caller": "genomeos/web/server.py Api._integrate <- Api.cell_run(silence=True)",
        "clamp": "{'<gene>.mRNA': 0.0} for every gene the cell type silences",
        "kind": "mRNA knockout at zero, but only of genes the runtime already silences (derivative "
        "-delta_m * x, which is 0 at 0), so the stages never left zero",
        "reached_from": ["web Cell tab 'Run and watch', the 'hold at zero' toggle (/api/cell/run silence=1)"],
        "dependents": [
            "tests/test_cell_view.py::"
            "test_a_gene_the_cell_does_not_express_stays_at_zero_without_being_held_down "
            "(data/demo/cell_context.bio, Neuron2, 20 h): silence is a no-op on the levels",
        ],
    },
]

NOT_CALLERS: list[dict[str, str]] = [
    {
        "where": "genomeos/runtime/body.py Body._network_step",
        "why": "its own Euler loop over vm._derivatives that skips clamped species "
        "(`if s in clamp: continue`); "
        "Euler has one stage, so the signal is held exactly; unaffected by the defect or the fix",
    },
    {
        "where": "genomeos/forge/network_experiment.py",
        "why": "deliberately avoids clamp=: a node knockout zeroes basal, max and produces instead (8e705f7)",
    },
    {
        "where": "BioLang experiment blocks (genomeos/bio.py, LocatedRuntime)",
        "why": "knockouts go to LocatedRuntime(knockouts=...) or the network experiment; "
        "the v0.3 `clamps` clause "
        "named in docs/ARCHITECTURE.md is not parsed anywhere, so no program reaches clamp=",
    },
    {"where": "scripts/", "why": "no script passes clamp="},
    {"where": "genomeos/web/static/index.html", "why": "reaches clamp only through /api/cell/run silence"},
]


def _grep_callers() -> list[str]:
    """Every `clamp=` in Python source outside this script and the runtime's own signature."""
    out = subprocess.run(
        ["git", "grep", "-n", "clamp=", "--", "genomeos", "scripts", "tests"],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=False,
    ).stdout
    keep = []
    for line in out.splitlines():
        path = line.split(":", 1)[0]
        if path in ("scripts/grn_clamp_census.py", "tests/test_grn_clamp.py") or not path.endswith(".py"):
            continue
        if re.search(r"clamp=\w|clamp=\{", line):
            keep.append(line.strip())
    return keep


def _stage_deviation(source: str, clamp: dict[str, float], initial: dict[str, float]) -> float:
    """Largest distance of a clamped species from its clamp over every RK stage of a 5 h run."""
    vm = NetworkRuntime(parse(source), seed=0)
    inner, worst = vm._derivatives, [0.0]

    def spy(state: dict[str, float]) -> dict[str, float]:
        for k, v in clamp.items():
            worst[0] = max(worst[0], abs(state.get(k, 0.0) - v))
        return inner(state)

    vm._derivatives = spy  # type: ignore[method-assign]
    vm.run(hours=5.0, dt=0.02, initial=initial, record_every=5, clamp=clamp)
    return worst[0]


def measure() -> dict[str, Any]:
    from test_grn_clamp import RING, SIGNAL

    def mean2(xs: list[float]) -> float:
        h = len(xs) // 2
        return sum(xs[h:]) / len(xs[h:])

    ring = parse(RING)
    free = NetworkRuntime(ring, seed=0).run(hours=60, dt=0.02, initial={"A": 10.0}, record_every=10)
    held = NetworkRuntime(ring, seed=0).run(
        hours=60, dt=0.02, initial={"A": 10.0}, record_every=10, clamp={"a.mRNA": 0.0}
    )
    v: dict[str, Any] = {
        "ring_A_second_half_mean_unperturbed": mean2(free.levels["A"]),
        "ring_A_second_half_mean_a_mRNA_clamped_0": mean2(held.levels["A"]),
        "stage_max_deviation": {
            "ring a.mRNA clamped 0": _stage_deviation(RING, {"a.mRNA": 0.0}, {"A": 10.0}),
            "ring B clamped 3": _stage_deviation(RING, {"B": 3.0}, {"A": 10.0}),
            "signal S clamped 2": _stage_deviation(SIGNAL, {"S": 2.0}, {}),
        },
    }
    for label, kw in {
        "test (cells 60, 30 h, dt 0.05)": {"cells": 60, "hours": 30, "dt": 0.05},
        "default = cli and web (cells 120, 40 h, dt 0.05)": {},
    }.items():
        r = run_gastrulation(**kw)
        v[f"gastrulation {label}"] = {
            "proportions": r.proportions(),
            "boundaries": [i for i in range(1, len(r.fates)) if r.fates[i] != r.fates[i - 1]],
            "order": [f for i, f in enumerate(r.fates) if i == 0 or f != r.fates[i - 1]],
        }
    api = Api(ROOT)
    for sil in (False, True):
        o = api.cell_run("data/demo/cell_context.bio", "Neuron2", hours=20, silence=sil)
        v[f"cell_run Neuron2 20 h silence={sil} final"] = {s: xs[-1] for s, xs in o["levels"].items()}
    return v


def main(label: str) -> None:
    p = RESULTS_DIR / f"{NAME}.json"
    old = json.loads(p.read_text()) if p.exists() else {}
    measured = dict(old.get("measured", {}))
    measured[label] = measure()
    payload = {
        "question": "which callers of NetworkRuntime.run(clamp=...) and which recorded figures or "
        "pinned tests "
        "depend on it, and what each reads on the runtime before and after the clamp is held inside every "
        "Runge-Kutta stage",
        "callers": CALLERS,
        "not_callers": NOT_CALLERS,
        "grep_clamp_equals": _grep_callers(),
        "measured": measured,
    }
    manifest = {
        "sources": [
            {"accession": "data/demo/gastrulation.bio", "version": "in-repo demo module"},
            {"accession": "data/demo/cell_context.bio", "version": "in-repo demo module"},
            {"accession": "tests/test_grn_clamp.py RING and SIGNAL", "version": "in-repo fixtures"},
        ],
        "inputs": [
            mf.input_entry("data/demo/gastrulation.bio"),
            mf.input_entry("data/demo/cell_context.bio"),
            mf.input_entry("genomeos/runtime/grn.py"),
            mf.input_entry("tests/test_grn_clamp.py"),
        ],
        "assembly": "n/a: simulated networks, no genome coordinates",
        "coordinates": "n/a: no genomic intervals",
        "parameters": {"label": label, "ring": "60 h, dt 0.02, A0 10", "stage_runs": "5 h, dt 0.02"},
        "exclusions": [],
        "partitions": "n/a: no evaluation split",
    }
    save_result(NAME, payload, manifest=manifest)
    print(json.dumps(measured[label], indent=1))


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "current")
