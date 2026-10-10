# SPDX-License-Identifier: AGPL-3.0-or-later
"""Score the SBML engines against BioModels' curation figures: ROADMAP §5 item 13, area C.

    uv run python scripts/roadrunner_score.py            # both engines, every model, needs the extra
    uv run python scripts/roadrunner_score.py --engines standin

Each model in genomeos/molecules/sbml_reference.py is run on the in-house stand-in (runtime/sbml.py,
RK4 at the dt its committed run used) and on libRoadRunner (runtime/roadrunner.py, CVODE at its
defaults), on one time grid, and every feature is checked against the tolerance fixed in that module
before these runs. A model passes on an engine only when every check passes. The two engines'
trajectories are also compared with each other (`AGREEMENT`), because a model the stand-in "passed"
could still disagree with the standard simulator. Schoeberl2002, the stand-in's recorded failure, is
also run on the stand-in at a tenth of its dt, as a diagnostic that is not scored.

The SBML files are BioModels' public downloads, cached git-ignored in data/knowledge/biomodels (the
repressilator is committed in data/models). No money, no key.
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.molecules import biomodels  # noqa: E402
from genomeos.molecules import sbml_reference as ref  # noqa: E402
from genomeos.results import save_result  # noqa: E402
from genomeos.runtime.sbml import SbmlModel, runtime  # noqa: E402

RESULT = "roadrunner_scores"
COMMITTED_MODEL = {"BIOMD0000000012": Path("data/models/BIOMD0000000012.xml")}
SCHOEBERL = "BIOMD0000000019"
SCHOEBERL_DEFAULT = "50 ng/ml (x1 = 4962, the SBML default)"
OWN_CODE = [
    "scripts/roadrunner_score.py",
    "genomeos/molecules/sbml_reference.py",
    "genomeos/molecules/biomodels.py",
    "genomeos/runtime/sbml.py",
    "genomeos/runtime/roadrunner.py",
]


def model_path(model_id: str) -> Path:
    return COMMITTED_MODEL.get(model_id) or biomodels.fetch(model_id)


def simulate(path: Path, engine: str, spec: dict, dose: dict, dt: float) -> tuple[Any, Any, float]:
    model = SbmlModel.from_file(path)
    rt = runtime(model, engine)
    every = max(1, round(spec["record_dt"] / dt))
    t0 = time.time()
    traj = rt.run(duration=spec["duration"], dt=dt, record_every=every, initial=dose["initial"] or None)
    return model, traj, round(time.time() - t0, 2)


def score_run(model, traj, spec: dict, dose: dict) -> list[dict]:
    checks: list[dict] = []
    for feature in dose["features"]:
        xs = ref.series(model, traj, feature[1], dose["initial"])
        checks += ref.score_feature(feature, traj.times, xs, spec["axes"])
    return checks


def run_model(model_id: str, engines: list[str]) -> dict[str, Any]:
    spec = ref.REFERENCES[model_id]
    path = model_path(model_id)
    entry: dict[str, Any] = {"name": spec["name"], "figure": spec["figure"], "file": str(path), "doses": {}}
    for dose_name, dose in spec["doses"].items():
        runs: dict[str, Any] = {}
        trajs: dict[str, Any] = {}
        for engine in engines:
            dt = spec["standin_dt"]  # the stand-in's step; for libRoadRunner only the grid it reports on
            try:
                model, traj, seconds = simulate(path, engine, spec, dose, dt)
            except Exception as e:  # an engine that cannot run the model fails it, with the reason
                runs[engine] = {"pass": False, "error": f"{type(e).__name__}: {e}"[:400]}
                print(f"{model_id} {dose_name} {engine}: ERROR {e}", flush=True)
                continue
            checks = score_run(model, traj, spec, dose)
            trajs[engine] = traj
            passed = sum(c["pass"] for c in checks)
            runs[engine] = {
                "pass": bool(checks) and passed == len(checks),
                "checks_passed": passed,
                "checks": len(checks),
                "dt": dt,
                "seconds": seconds,
                "detail": checks,
            }
            print(f"{model_id} {dose_name} {engine}: {passed}/{len(checks)} in {seconds}s", flush=True)
        if len(trajs) == 2:
            runs["agreement"] = ref.agreement(trajs["standin"], trajs["roadrunner"])
        entry["doses"][dose_name] = runs
    if model_id == SCHOEBERL and "standin" in engines:
        dose = spec["doses"][SCHOEBERL_DEFAULT]
        dt = spec["standin_dt"] / 10
        model, traj, seconds = simulate(path, "standin", spec, dose, dt)
        checks = score_run(model, traj, spec, dose)
        passed = sum(c["pass"] for c in checks)
        entry["diagnostic_standin_dt_over_10"] = {
            "dt": dt,
            "seconds": seconds,
            "checks_passed": passed,
            "checks": len(checks),
            "scored": False,
            "detail": checks,
        }
        print(f"{model_id} diagnostic stand-in dt={dt}: {passed}/{len(checks)} in {seconds}s", flush=True)
    has_features = any(d["features"] for d in spec["doses"].values())
    for engine in engines:
        scored = [d[engine] for d in entry["doses"].values() if engine in d]
        entry[f"{engine}_passes"] = (
            (bool(scored) and all(r["pass"] for r in scored)) if has_features else None
        )
    if len(engines) == 2:
        entry["engines_agree"] = all(
            d.get("agreement", {}).get("agree", False) for d in entry["doses"].values()
        )
    return entry


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=(__doc__ or "").split("\n\n")[0])
    ap.add_argument("--engines", nargs="+", default=["standin", "roadrunner"])
    ap.add_argument("--models", nargs="+", default=list(ref.REFERENCES))
    args = ap.parse_args(argv)

    t_start = time.time()
    settings = None
    if "roadrunner" in args.engines:
        from genomeos.runtime.roadrunner import RoadRunnerRuntime

        settings = RoadRunnerRuntime(SbmlModel.from_file(model_path("BIOMD0000000084"))).settings()
    models = {m: run_model(m, args.engines) for m in args.models}
    table = {
        m: {e: entry.get(f"{e}_passes") for e in args.engines} | {"engines_agree": entry.get("engines_agree")}
        for m, entry in models.items()
    }
    out = {
        "result": RESULT,
        "question": (
            "Does libRoadRunner, behind the stand-in's interface, reproduce Schoeberl2002 and the models the "
            "stand-in already ran, against BioModels' curation figures, at a tolerance fixed before the runs?"
        ),
        "tolerance": ref.TOLERANCE,
        "agreement_threshold": ref.AGREEMENT,
        "roadrunner": settings,
        "table": table,
        "models": models,
        "seconds": round(time.time() - t_start, 1),
        "money": "none: public BioModels files, no key",
    }
    paths = [model_path(m) for m in args.models]
    man = {
        "sources": [
            {
                "accession": "https://www.ebi.ac.uk/biomodels",
                "version": "the SBML main files as downloaded, sha256 per input; references read from the "
                "entries' curation_image files (genomeos/molecules/sbml_reference.py)",
            }
        ],
        "inputs": [mf.input_entry(p, partition=None) for p in paths],
        "input_count": len(paths),
        "assembly": "n/a: kinetic models, no genome coordinates",
        "coordinates": "n/a: kinetic models, no genome coordinates",
        "parameters": {"tolerance": ref.TOLERANCE, "agreement": ref.AGREEMENT, "engines": args.engines},
        "exclusions": [
            {
                "reason": "McClean2007: its curation figure is a parameter scan, not a time course; "
                "compared on both engines, not scored",
                "count": 1,
            }
        ],
        "partitions": "n/a: nothing was fitted; the references are read from figures, not trained on",
        "code_cleanliness": mf.code_cleanliness(__file__, OWN_CODE),
    }
    print("saved", save_result(RESULT, out, manifest=man))
    for m, row in table.items():
        print(m, row)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
