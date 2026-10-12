# SPDX-License-Identifier: AGPL-3.0-or-later
"""Score the Boolean engines against a published MaBoSS figure: ROADMAP §5 item 13, area C.

    uv run --extra maboss python scripts/maboss_score.py                  # all three engines
    uv run python scripts/maboss_score.py --engines standin standin_async  # without the extra

The model is Cohen et al. 2015's metastasis model in its published MaBoSS files
(data/models/cohen2015_metastasis.bnd and .cfg, checked against the sha256 recorded in
genomeos/molecules/maboss_reference.py before anything runs), and the reference is the four pie charts
of the paper's Fig 3. The references, the phenotype rule, the mutant rule, the tolerance and what each
engine is expected to do were fixed in that module before this script's first scored run.

Three engines, each in every panel:

    maboss          runtime/maboss.py, the published cfg as published (seed 100, 50,000 trajectories)
    standin_async   runtime/boolean.py's step_async under an exponential clock, 50,000 trajectories,
                    seed 100, on the network maboss.to_network translates (it refuses where the two
                    processes differ)
    standin         runtime/boolean.py's synchronous attractor from each of the 16 initial states

A panel passes on an engine only when every check passes; the model passes when all four do. MaBoSS
and standin_async are also compared with each other (total variation over the output states). No
money, no key: the model files are the paper's public supplement.
"""

from __future__ import annotations

import argparse
import hashlib
import random
import sys
import time
from collections import Counter
from itertools import product
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.molecules import maboss_reference as ref  # noqa: E402
from genomeos.results import save_result  # noqa: E402
from genomeos.runtime import maboss  # noqa: E402
from genomeos.runtime.boolean import BooleanNetwork  # noqa: E402

RESULT = "maboss_scores"
OWN_CODE = [
    "scripts/maboss_score.py",
    "genomeos/molecules/maboss_reference.py",
    "genomeos/runtime/maboss.py",
    "genomeos/runtime/boolean.py",
    "data/models/cohen2015_metastasis.bnd",
    "data/models/cohen2015_metastasis.cfg",
]

State = tuple[str, ...]


def published() -> maboss.MaBoSSModel:
    """The published pair, refused unless both files are byte for byte the paper's."""
    for kind in ("bnd", "cfg"):
        digest = hashlib.sha256(Path(ref.ROOT_FILES[kind]).read_bytes()).hexdigest()
        if digest != ref.SOURCE[f"{kind}_sha256"]:
            raise SystemExit(f"{ref.ROOT_FILES[kind]} is not the published file: sha256 {digest}")
    return maboss.MaBoSSModel.from_files(ref.ROOT_FILES["bnd"], ref.ROOT_FILES["cfg"])


def _key(state: dict[str, bool]) -> State:
    return tuple(sorted(o for o in ref.OUTPUTS if state[o]))


def run_maboss(model: maboss.MaBoSSModel) -> tuple[dict[State, float], dict, float]:
    t0 = time.time()
    rt = maboss.MaBoSSRuntime(model)
    res = rt.run()
    return res.states, {"window": list(res.window), **res.settings}, round(time.time() - t0, 2)


def run_standin_async(
    net: BooleanNetwork, initial: dict[str, bool], random_nodes: list[str], rate: float, spec: dict
) -> tuple[dict[State, float], dict, float]:
    """The stand-in's step_async under an exponential clock of rate (nodes x rate): MaBoSS's process
    when every rate is equal. A trajectory that reaches a fixed point stops there (it is absorbing)."""
    t0 = time.time()
    rng = random.Random(spec["seed"])
    nodes = net.nodes
    clock = len(nodes) * rate
    counts: Counter[State] = Counter()
    stopped_at_fixed_point = 0
    for _ in range(spec["trajectories"]):
        state = {n: bool(initial.get(n, False)) for n in nodes}
        for n in random_nodes:
            state[n] = rng.random() < 0.5
        t, steps = 0.0, 0
        while True:
            t += rng.expovariate(clock)
            if t > spec["t_end"]:
                break
            state = net.step_async(state, rng)
            steps += 1
            if steps % len(nodes) == 0 and net.step(state) == state:
                stopped_at_fixed_point += 1
                break
        counts[_key(state)] += 1
    n = spec["trajectories"]
    probs = {k: v / n for k, v in counts.items()}
    detail = {**spec, "clock_rate": clock, "stopped_at_a_fixed_point_before_t_end": stopped_at_fixed_point}
    return probs, detail, round(time.time() - t0, 2)


def run_standin_sync(
    net: BooleanNetwork, initial: dict[str, bool], random_nodes: list[str]
) -> tuple[dict[State, float], dict, float]:
    """The synchronous attractor from every initial state the random nodes allow, equally weighted."""
    t0 = time.time()
    if len(random_nodes) > 16:
        raise SystemExit(f"{len(random_nodes)} random nodes: too many initial states to enumerate")
    probs: dict[State, float] = {}
    lengths: Counter[int] = Counter()
    combos = list(product((False, True), repeat=len(random_nodes)))
    for bits in combos:
        att = net.attractor({**initial, **dict(zip(random_nodes, bits, strict=True))})
        lengths[len(att)] += 1
        for s in att:
            probs[_key(s)] = probs.get(_key(s), 0.0) + 1.0 / len(combos) / len(att)
    detail = {
        "initial_states": len(combos),
        "attractor_lengths": {str(k): v for k, v in sorted(lengths.items())},
    }
    return probs, detail, round(time.time() - t0, 2)


def _states_out(states: dict[State, float]) -> dict[str, float]:
    return {
        (" -- ".join(k) or maboss.NIL): round(v, 6) for k, v in sorted(states.items(), key=lambda kv: -kv[1])
    }


def run_panel(base: maboss.MaBoSSModel, panel: str, engines: list[str]) -> dict:
    held = ref.CONDITIONS[panel]
    model = base.held(held) if held else base
    entry: dict = {"held": {k: int(v) for k, v in held.items()}, "figure": ref.FIGURE[panel], "engines": {}}
    translated = None
    if {"standin", "standin_async"} & set(engines):
        translated = maboss.to_network(model)
        entry["standin_translation"] = {
            "nodes": len(translated[0].nodes),
            "random_nodes": translated[2],
            "common_rate": translated[3],
            "constant_rules": sorted(n for n, r in translated[0].rules.items() if r in ("0", "1")),
        }
    finals: dict[str, dict[State, float]] = {}
    for engine in engines:
        if engine == "maboss":
            states, detail, secs = run_maboss(model)
        elif engine == "standin_async":
            net, initial, random_nodes, rate = translated  # type: ignore[misc]
            states, detail, secs = run_standin_async(
                net, initial, random_nodes, rate, ref.SETTINGS["standin_async"]
            )
        elif engine == "standin":
            net, initial, random_nodes, _ = translated  # type: ignore[misc]
            states, detail, secs = run_standin_sync(net, initial, random_nodes)
        else:
            raise SystemExit(f"unknown engine {engine!r}; choose from {ref.ENGINES}")
        finals[engine] = states
        probs = ref.phenotypes(states)
        checks = ref.score(panel, probs)
        entry["engines"][engine] = {
            "phenotypes": {k: round(v, 6) for k, v in sorted(probs.items())},
            "checks": checks,
            "passed": sum(c["pass"] for c in checks),
            "of": len(checks),
            "states": _states_out(states),
            "run": detail,
            "seconds": secs,
        }
    if "maboss" in finals and "standin_async" in finals:
        tv = ref.total_variation(finals["maboss"], finals["standin_async"])
        entry["agreement"] = {
            "total_variation": round(tv, 6),
            "bound": ref.AGREEMENT,
            "agree": tv <= ref.AGREEMENT,
        }
    return entry


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=(__doc__ or "").split("\n\n")[0])
    ap.add_argument("--engines", nargs="+", default=list(ref.ENGINES))
    args = ap.parse_args(argv)

    t_start = time.time()
    base = published()
    panels = {p: run_panel(base, p, args.engines) for p in ref.CONDITIONS}
    table = {}
    for engine in args.engines:
        passed = sum(panels[p]["engines"][engine]["passed"] for p in panels)
        total = sum(panels[p]["engines"][engine]["of"] for p in panels)
        table[engine] = {
            "checks_passed": f"{passed}/{total}",
            "passes": passed == total,
            "per_panel": {
                p: f"{panels[p]['engines'][engine]['passed']}/{panels[p]['engines'][engine]['of']}"
                for p in panels
            },
            "expected": ref.EXPECTED[engine],
        }
    if "agreement" in next(iter(panels.values())):
        table["agreement"] = {
            "all_panels_agree": all(panels[p]["agreement"]["agree"] for p in panels),
            "total_variation": {p: panels[p]["agreement"]["total_variation"] for p in panels},
            "expected": ref.EXPECTED["agreement"],
        }
    out = {
        "result": RESULT,
        "question": (
            "Does MaBoSS, behind runtime/maboss.py, reproduce the phenotype probabilities Cohen et al. "
            "2015 print in Fig 3 from their published MaBoSS files, and what do the stand-in's two update "
            "schemes give, at a tolerance fixed before the runs?"
        ),
        "source": ref.SOURCE,
        "tolerance": ref.TOLERANCE,
        "agreement_threshold": ref.AGREEMENT,
        "phenotype_rule": {
            "outputs": list(ref.OUTPUTS),
            "precedence": ["HS if none", *ref.PRECEDENCE, "CCA"],
        },
        "settings": ref.SETTINGS,
        "table": table,
        "panels": panels,
        "seconds": round(time.time() - t_start, 1),
        "money": "none: the paper's public supplement, no key",
    }
    paths = [ref.ROOT_FILES["bnd"], ref.ROOT_FILES["cfg"]]
    man = {
        "sources": [
            {
                "accession": ref.SOURCE["model_files"],
                "version": f"S1 File of doi:{ref.SOURCE['doi']}, "
                f"zip sha256 {ref.SOURCE['model_zip_sha256']}; references read from {ref.SOURCE['figure']}",
            }
        ],
        "inputs": [mf.input_entry(p, partition=None) for p in paths],
        "input_count": len(paths),
        "assembly": "n/a: a logical model, no genome coordinates",
        "coordinates": "n/a: a logical model, no genome coordinates",
        "parameters": {
            "tolerance": ref.TOLERANCE,
            "agreement": ref.AGREEMENT,
            "engines": args.engines,
            "settings": ref.SETTINGS,
        },
        "exclusions": [],
        "partitions": "n/a: nothing was fitted; the references are printed in the paper's figure",
        "code_cleanliness": mf.code_cleanliness(__file__, OWN_CODE),
    }
    print("saved", save_result(RESULT, out, manifest=man))
    for engine, row in table.items():
        print(engine, row)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
