# SPDX-License-Identifier: AGPL-3.0-or-later
"""M1b run: the registered network-level test of the compiled regulator rules, once.

    uv run python scripts/regulators_run.py

Refuses unless data/results/regulators_registration.json exists, the analysis code
(genomeos/attribution/regulators.py and this file) has the sha256 the registration pinned, the test
inputs have theirs, and every knockdown file was fetched with ENCODE's md5. Refuses a second run.
Applies the registered conditions in order, then `regulators.analyse`, and writes
data/results/regulators_result.json with the registered reading verbatim.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import math
import sys
import time
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import regulators as rg  # noqa: E402
from genomeos.results import RESULTS_DIR, save_result  # noqa: E402

RESULT = "regulators_result"
REG = RESULTS_DIR / "regulators_registration.json"
LEDGER = rg.KD_DIR / "fetched.json"
OWN_CODE = (
    "genomeos/attribution/regulators.py",
    "scripts/regulators_inventory.py",
    "scripts/regulators_census.py",
    "scripts/regulators_register.py",
    "scripts/regulators_fetch.py",
    "scripts/regulators_run.py",
    "tests/test_regulators.py",
)


def sha256(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def tpm(acc: str) -> dict[str, float]:
    with gzip.open(rg.KD_DIR / f"{acc}.tsv.gz", "rt") as fh:
        return rg.parse_quantification(fh)


def checks(reg: dict) -> list[str]:
    problems = []
    for p, h in reg["frozen_code"]["sha256"].items():
        if sha256(Path(p)) != h:
            problems.append(f"{p} differs from the registered code")
    ti = reg["test_inputs"]
    if sha256(Path(ti["path"])) != ti["sha256"]:
        problems.append("the test inputs differ from the registered ones")
    ledger = json.loads(LEDGER.read_text()) if LEDGER.exists() else {}
    for u in reg["units"]:
        for e in u["experiments"]:
            for f in e["knockdown_files"] + e["control_files"]:
                got = ledger.get(f["accession"])
                if not got or got["md5"] != f["md5sum"]:
                    problems.append(f"{f['accession']} not fetched with ENCODE's md5")
                elif sha256(rg.KD_DIR / f"{f['accession']}.tsv.gz") != got["gz_sha256"]:
                    problems.append(f"{f['accession']} changed since it was fetched")
    return problems


def direction_signs(key, accession, universe, arm, kd, ctrl, directed) -> list[dict]:
    """For directed rows: the share of genes whose signed change has the predicted sign (activates ->
    lower on knockdown, inhibits -> higher), beside the share among the universe's genes outside the arm."""
    signed = {
        g: math.log2((kd.get(g, 0.0) + rg.PSEUDO_TPM) / (ctrl.get(g, 0.0) + rg.PSEUDO_TPM)) for g in universe
    }
    others = [g for g in universe if g not in arm]
    out = []
    for direction, want in (("activates", -1), ("inhibits", 1)):
        gs = [g for g in directed.get(direction, []) if g in signed]
        if gs:
            out.append(
                {
                    "unit": key,
                    "accession": accession,
                    "direction": direction,
                    "genes": len(gs),
                    "share_with_predicted_sign": sum(1 for g in gs if signed[g] * want > 0) / len(gs),
                    "share_among_other_universe_genes": sum(1 for g in others if signed[g] * want > 0)
                    / max(1, len(others)),
                }
            )
    return out


def main() -> int:
    t0 = time.time()
    if (RESULTS_DIR / f"{RESULT}.json").exists():
        raise SystemExit("refused: the registered test has already run")
    reg = json.loads(REG.read_text())
    problems = checks(reg)
    if problems:
        raise SystemExit("refused before any value was read:\n" + "\n".join(problems))
    inputs = json.loads(gzip.decompress(Path(reg["test_inputs"]["path"]).read_bytes()))
    tss = {g: (v[0], v[1]) for g, v in inputs["tss"].items()}
    cache: dict[str, dict[str, float]] = {}

    def mean_of(files: list[dict]) -> dict[str, float]:
        reps = []
        for f in files:
            if f["accession"] not in cache:
                cache[f["accession"]] = tpm(f["accession"])
            reps.append(cache[f["accession"]])
        return rg.mean_tpm(reps)

    units, rows, signs, motif_rows = [], [], [], []
    excluded: Counter = Counter()
    for u in reg["units"]:
        key = f"{u['factor']} {u['cell']}"
        fgene = inputs["factor_gene"][u["factor"]]
        arm_all = set(inputs["arms"][key])
        arm_motif_all = set(inputs["arms_with_motif"][key])
        views = []
        for e in u["experiments"]:
            kd, ctrl = mean_of(e["knockdown_files"]), mean_of(e["control_files"])
            ratio = rg.knockdown_ratio(fgene, kd, ctrl)
            row = {
                "factor": u["factor"],
                "cell": u["cell"],
                "accession": e["accession"],
                "assay": e["assay"],
                "knockdown_ratio": ratio,
            }
            if ratio is None:
                row["excluded"] = "knockdown_not_assessable"
            elif ratio > rg.KD_MAX_RATIO:
                row["excluded"] = "knockdown_below_threshold"
            else:
                universe = rg.universe_of(inputs["rule_genes"][u["cell"]], ctrl, fgene, tss)
                arm = arm_all & set(universe)
                row.update({"universe": len(universe), "rule_arm": len(arm)})
                if len(arm) < rg.MIN_ARM or len(universe) - len(arm) < rg.MIN_ARM:
                    row["excluded"] = "arm_too_small"
                else:
                    near = rg.near_peaks(universe, tss, inputs["peaks"][key])
                    base = rg.baseline_arm(universe, len(arm), near)
                    view = rg.experiment_view(universe, ctrl, rg.response(kd, ctrl, universe), arm, base)
                    views.append(view)
                    row["baseline_arm"] = len(base)
                    row["baseline_overlap_with_rule_arm"] = len(base & arm)
                    marm = arm_motif_all & set(universe)
                    if len(marm) >= rg.MIN_ARM and len(universe) - len(marm) >= rg.MIN_ARM:
                        lab = [g in marm for g in universe]
                        motif_rows.append(
                            {
                                "unit": key,
                                "accession": e["accession"],
                                "genes": len(marm),
                                "auroc": rg.stratified_auroc(view["y"], lab, view["strata"]),
                            }
                        )
                    signs += direction_signs(
                        key, e["accession"], universe, arm, kd, ctrl, inputs["directed"].get(key, {})
                    )
            if "excluded" in row:
                excluded[row["excluded"]] += 1
            rows.append(row)
        if views:
            units.append({"factor": u["factor"], "cell": u["cell"], "experiments": views})
    for unit in units:
        same = [x for x in units if x["cell"] == unit["cell"]]
        unit["arms"] = {x["factor"]: set(inputs["arms"][f"{x['factor']} {x['cell']}"]) for x in same}
    analysis = rg.analyse(units)
    k562 = [x for x in analysis.get("per_unit", []) if x["cell"] == "K562"]
    k562_only = None
    if k562:
        diffs = [x["auroc_rules"] - x["auroc_baseline"] for x in k562]
        lo, hi = rg.bootstrap_ci(diffs)
        k562_only = {"units": len(k562), "mean_paired_difference": sum(diffs) / len(diffs), "ci95": [lo, hi]}
    out = {
        "lane": rg.LANE,
        "registration": "regulators_registration",
        "question": reg["question"],
        "outcome": analysis["outcome"],
        "reading": analysis["reading"],
        "analysis": {k: v for k, v in analysis.items() if k not in ("outcome", "reading")},
        "experiments": rows,
        "experiments_excluded": dict(excluded),
        "descriptives_not_verdicts": {
            "rule_arm_with_motif": motif_rows,
            "k562_units_alone": k562_only,
            "direction_signs": signs,
        },
        "limitations": reg["limitations"],
        "files_read": len(cache),
        "alphagenome_requests": 0,
        "money": "none: no request is sent and no key is read",
        "seconds": round(time.time() - t0, 1),
    }
    man = {
        "sources": [
            {
                "accession": "https://www.encodeproject.org",
                "version": "files pinned by md5 in the registration",
            },
            {
                "accession": "data/results/regulators_registration.json",
                "version": "the committed registration",
            },
        ],
        "inputs": [
            mf.input_entry(REG, partition=None),
            mf.input_entry(Path(reg["test_inputs"]["path"]), partition=None),
            mf.input_entry(LEDGER, partition=None),
        ],
        "input_count": 3,
        "assembly": "GRCh38; the knockdown files' own GENCODE annotation per unit",
        "coordinates": "n/a: per-gene expression",
        "parameters": rg.TEST_CONSTANTS,
        "exclusions": [{"reason": k, "count": v} for k, v in sorted(excluded.items())],
        "partitions": "n/a: the endpoint is outside every training set this project fitted",
        "code_cleanliness": mf.code_cleanliness(__file__, OWN_CODE),
    }
    print("saved", save_result(RESULT, out, manifest=man))
    print(out["outcome"])
    print(out["reading"])
    print(json.dumps(out["experiments_excluded"]), k562_only)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
