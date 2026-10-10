# SPDX-License-Identifier: AGPL-3.0-or-later
"""The registered run: the multi-target compiled choice on M1c's 1,918 held-out K562 CRISPRi pairs.

    uv run python scripts/targetchoice_run.py

Runs only from the committed registration (data/results/targetchoice_registration.json) and refuses if
any registered field differs from genomeos/attribution/targetchoice.py, if its own code or the code it
calls differs from what the registration pinned, or if a pinned input's sha256 moved.

Order, as registered (targetchoice.STOP_RULE):
  1. the population is read and its ids and sha256 checked against the registration (M1c's population);
  2. the anchor's columns are computed (crispri.annotate) and every choice's columns are computed for the
     fit rows and the population, one chromosome at a time; no label is read;
  3. the predictor-side harness (targetchoice.harness) is checked; a defect stops the run here;
  4. the anchor is recomputed and must equal the committed arm-1 figures;
  5. the one-target columns must equal M1c's committed s2_target columns pair for pair, and the
     one-target reading must equal M1c's committed s2_target reading;
  6. only then is the multi-target choice (and the top-k choice) read and decided.

The run reads its own registration and declares it among its inputs (M1c's first write was quarantined
for not doing so). No model request, no network, no money: every input is on this machine.
"""

from __future__ import annotations

import hashlib
import json
import sys
import time
from collections import Counter
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import crispri  # noqa: E402
from genomeos.attribution import signaltrace as st  # noqa: E402
from genomeos.attribution import targetchoice as tc  # noqa: E402
from genomeos.attribution.measured import CRISPRI_SPLIT_OF  # noqa: E402
from genomeos.results import RESULTS_DIR, save_result  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import targetchoice_census as census  # noqa: E402

RESULT = "targetchoice_result"
REGISTRATION = RESULTS_DIR / "targetchoice_registration.json"
M1C_REGISTRATION = RESULTS_DIR / "signaltrace_registration.json"
M1C_RESULT = RESULTS_DIR / "signaltrace_result.json"
CENSUS = Path(tc.CENSUS)
OWN_CODE = (
    "genomeos/attribution/targetchoice.py",
    "scripts/targetchoice_census.py",
    "scripts/targetchoice_register.py",
    "scripts/targetchoice_run.py",
    "tests/test_targetchoice.py",
)
#: the code the run must equal, pinned by sha256 in the registration; the tests are not on the run's path
FROZEN = OWN_CODE[:4]
#: the committed code the choices run, not this lane's; pinned by sha256 in the registration
STAGE_CODE = (
    "genomeos/attribution/crispri.py",
    "genomeos/attribution/compile.py",
    "genomeos/attribution/signaltrace.py",
    "genomeos/attribution/readergate.py",
    "genomeos/attribution/targets.py",
    "genomeos/predict/enhancer_target.py",
    "genomeos/genome/annotation.py",
    "genomeos/ir/model.py",
    "genomeos/lang/parser.py",
)


def load_pairs() -> tuple[list[crispri.Pair], list[crispri.Pair]]:
    """The training pairs and the population: the held-out file's K562 pairs, in file order."""
    training, heldout = crispri.load(crispri.TRAINING), crispri.load(crispri.HELDOUT)
    return training, [p for p in heldout if p.cell == tc.CELL]


def chromosomes_of(training: list[crispri.Pair], population_ids: list[str]) -> list[str]:
    seen = {p.chrom for p in training} | {i.split(":", 1)[0] for i in population_ids}
    return [c for c in census.CHROMS if c in seen]


def code_sha256(paths: tuple[str, ...]) -> dict[str, str]:
    return {p: hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in paths}


def inputs(chroms: list[str], n_training: int, n_population: int) -> list[dict[str, Any]]:
    """Every file the run reads but its own registration, pinned by sha256 as the registration pins them."""
    return [
        mf.input_entry(
            crispri.KNOWLEDGE / crispri.TRAINING,
            partition=CRISPRI_SPLIT_OF[crispri.TRAINING],
            pairs=n_training,
            role="fit rows: the covered K562 training pairs",
        ),
        mf.input_entry(
            crispri.KNOWLEDGE / crispri.HELDOUT,
            partition=CRISPRI_SPLIT_OF[crispri.HELDOUT],
            pairs=n_population,
            role="the population: its K562 pairs, scored",
        ),
        *census.shared_inputs(chroms),
        mf.input_entry(
            M1C_REGISTRATION, partition=None, role="M1c's registration: the population's ids and sha256"
        ),
        mf.input_entry(
            M1C_RESULT,
            partition=None,
            role="M1c's result: the s2_target columns and reading the one-target choice must reproduce",
        ),
        mf.input_entry(CENSUS, partition="training", role="the training census that fixed TOP_K"),
    ]


def registration_entry() -> dict[str, Any]:
    return mf.input_entry(
        REGISTRATION, partition=None, role="this run's registration, read before any reading"
    )


def check_registration() -> dict[str, Any]:
    """The committed registration; refused if a bound field or a pinned code file differs."""
    reg = json.loads(REGISTRATION.read_text())
    bad = [f for f, v in tc.registered_fields().items() if reg.get(f) != v]
    if bad:
        raise SystemExit(f"refused: the registration and targetchoice.py disagree on {', '.join(bad)}")
    moved = [p for p, h in code_sha256(FROZEN).items() if reg["frozen_code"]["sha256"].get(p) != h]
    moved += [p for p, h in code_sha256(STAGE_CODE).items() if reg["stage_code"]["sha256"].get(p) != h]
    if moved:
        raise SystemExit(f"refused: code differs from what the registration pinned: {moved}")
    return reg


def manifest(entries: list[dict[str, Any]]) -> dict[str, Any]:
    base = census.manifest(entries)
    return {
        **base,
        "sources": [
            {
                "accession": "EngreitzLab/CRISPR_comparison EPCrisprBenchmark training_K562 and"
                " heldout_5_cell_types (Gschwind et al.)",
                "version": "fetched 2026-09-16; pinned here by sha256",
                "url": crispri.BASE_URL,
            },
            *base["sources"][1:],
        ],
        "parameters": {
            "cell": tc.CELL,
            "choices": {c: {"columns": list(tc.columns(c))} for c in tc.read_choices(tc.TOP_K)},
            "min_effect": tc.MIN_EFFECT,
            "top_k": tc.TOP_K,
            "fit": "crispri.logistic_fit, lam 1e-3, 25 rounds, on the covered K562 training pairs",
            "bootstrap": tc.BOOTSTRAP,
        },
        "exclusions": [
            st.PREDICTED_ONLY,
            "held-out pairs in cell types other than K562: the anchor is a K562 result",
            "a pair a choice does not keep is kept in the population with 0 in that choice's two columns",
        ],
        "partitions": {
            "training": "the fit rows only: the covered K562 training pairs; TOP_K was fixed on them",
            "heldout": "the scored population: the held-out file's K562 pairs; no held-out pair enters a fit",
        },
        "code_cleanliness": mf.code_cleanliness(__file__, OWN_CODE),
    }


def pairs_table(pop: list[crispri.Pair], choices: tuple[str, ...]) -> dict[str, Any]:
    cols = ["pair", "regulated", "rank", *[f"{c}_{x}" for c in choices for x in ("kept", "drop")]]
    rows = [
        [
            st.pair_id(p),
            p.regulated,
            p.features["rank"],
            *[round(p.features[f"{c}_{x}"], 6) for c in choices for x in ("kept", "drop")],
        ]
        for p in pop
    ]
    return {"columns": cols, "rows": rows}


def descriptives(lookups: dict[str, tc.ChoiceLookup], pop: list[crispri.Pair]) -> dict[str, Any]:
    reached: set[tuple[str, str]] = set()
    why: Counter = Counter()
    for p in pop:
        cl = lookups.get(p.chrom)
        els = cl.overlapping(p.start, p.end) if cl else []
        reached.update((p.chrom, el) for el in els)
        if not p.features[f"{tc.MULTI}_kept"]:
            why[
                "no compiled element overlaps"
                if not els
                else "the overlapping elements list other genes only"
            ] += 1
    sizes = [len(lookups[c].candidates[el]) for c, el in sorted(reached)]
    ranks = Counter(
        "not listed" if p.features["rank"] is None else tc.size_band(p.features["rank"]) for p in pop
    )
    counts: Counter = Counter()
    for cl in lookups.values():
        counts.update(cl.counts)
    return {
        "compiled_elements_reached_by_a_population_pair": tc.element_census(sizes),
        "population_by_the_measured_genes_position_among_candidates": {
            b: ranks.get(b, 0) for b in ("1", "2", "3", "4-5", "6-10", "11+", "not listed")
        },
        "population_not_kept_by_the_multi_target_choice_because": dict(why.most_common()),
        "compiled_elements_on_the_chromosomes_read": dict(counts),
    }


def main() -> int:
    t0 = time.time()
    if (RESULTS_DIR / f"{RESULT}.json").exists():
        raise SystemExit("refused: the registered run has already been written")
    reg = check_registration()
    top_k = reg["top_k"]
    choices = tc.read_choices(top_k)
    training, pop = load_pairs()
    ids = st.population_ids(pop)
    chroms = chromosomes_of(training, ids)
    entries = inputs(chroms, len(training), len(pop))
    base: dict[str, Any] = {
        "lane": tc.LANE,
        "registration": "data/results/targetchoice_registration.json",
        "registration_date": reg.get("date"),
        "question": tc.QUESTION,
        "alphagenome_requests": 0,
        "money": "none: no request is sent and no key is read",
    }

    def defect(what: str) -> int:
        out = {
            **base,
            "outcome": "HARNESS_DEFECT",
            "reading": tc.reading("HARNESS_DEFECT", what=what),
            "seconds": round(time.time() - t0, 1),
        }
        print("saved", save_result(RESULT, out, manifest=manifest([*entries, registration_entry()])))
        print(out["reading"])
        return 1

    pinned = {e["path"]: e["sha256"] for e in reg["result_manifest"]["inputs"]}
    moved = [e["path"] for e in entries if pinned.get(e["path"]) != e["sha256"]]
    if moved:
        return defect(f"pinned inputs moved since the registration: {moved}")
    if st.ids_sha256(ids) != reg["population_sha256"] or ids != reg["population_ids"]:
        return defect("the population's ids or sha256 differ from the registration's")
    if json.loads(M1C_REGISTRATION.read_text())["population_sha256"] != reg["population_sha256"]:
        return defect("the registration's population is not M1c's")

    census.annotate(training + pop)
    fit_rows = [p for p in training if p.covered]
    lookups = census.lookups_for(chroms)
    for p in training + pop:
        cl = lookups.get(p.chrom)
        p.features.update(tc.trace_pair(p, cl, top_k=top_k if tc.TOPK in choices else None))
        p.features["rank"] = tc.rank_of(p, cl)
    problems = tc.harness(training + pop, choices)
    if problems:
        return defect(f"{len(problems)} harness problems: " + "; ".join(problems[:20]))
    anchor = st.anchor_reading(fit_rows, pop)
    if not anchor["reproduced"]:
        return defect(
            f"the anchor did not reproduce: recomputed {anchor['recomputed']}, committed {tc.ANCHOR}"
        )
    stray = tc.comparator_mismatches(pop, json.loads(M1C_RESULT.read_text())["pairs_table"])
    if stray:
        return defect(f"{len(stray)} one-target columns differ from M1c's s2_target: {stray[:10]}")
    one = tc.score_one(fit_rows, pop, anchor)
    if not tc.comparator_reproduced(one):
        shown = {k: one[k] for k in ("population", "gain_over_activity_distance", "against_anchor")}
        return defect(f"the one-target choice did not reproduce M1c's s2_target: {shown}")

    scored = tc.score(fit_rows, pop, anchor, one, choices, top_k)
    out = {
        **base,
        "population": {
            "pairs": len(pop),
            "positives": sum(p.regulated for p in pop),
            "sha256": reg["population_sha256"],
            "chromosomes": len({p.chrom for p in pop}),
        },
        "fit_rows": {"pairs": len(fit_rows), "positives": sum(p.regulated for p in fit_rows)},
        "anchor": {k: anchor[k] for k in ("recomputed", "committed", "reproduced")},
        "one_target": {
            **{k: v for k, v in one.items() if k != "scores"},
            "reproduces_m1c_s2_target": True,
        },
        "outcome": scored["outcome"],
        "reading": scored["reading"],
        "falsified": scored["falsified"],
        "readings_verbatim": {c: r["reading"] for c, r in scored["choices"].items()},
        "choices": scored["choices"],
        "top_k": top_k,
        "descriptives_not_verdicts": descriptives(lookups, pop),
        "limitations": list(tc.LIMITATIONS),
        "pairs_table": pairs_table(pop, choices),
        "seconds": round(time.time() - t0, 1),
    }
    print("saved", save_result(RESULT, out, manifest=manifest([*entries, registration_entry()])))
    print(json.dumps({"outcome": out["outcome"], "reading": out["reading"]}, indent=1))
    for c, r in out["readings_verbatim"].items():
        print(f"{c}: {r}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
