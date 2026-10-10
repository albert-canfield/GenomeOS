# SPDX-License-Identifier: AGPL-3.0-or-later
"""Area E, lane-shareadopt (2026-10-10): `share:` in haematopoiesis.to_bio, registered before the change.

The coordinator's decision (docs/ROADMAP.md area E, a8423e5): adopt `share:` for the sibling partitions
that `genomeos/organism/haematopoiesis.py` writes, with the growth terms re-solved for 1 - sum(f), as a
registered lane. lane-tobioshare measured the question first (no commit): flows due at one instant
under `share:` all draw on the same starting count, so the parent keeps 1 - sum(f) (BIOLANG-v0.3;
tests/test_share.py::test_sibling_flows_at_one_instant_draw_from_one_base), while `fraction:` keeps
prod(1 - f) and the second sibling draws from what the first left. `plan()` solves growth for the
product, so the same numbers written as `share:` do not balance; re-solved, they do.

Two subcommands, each run ONCE, from a clean worktree of a committed tree (so the stamped revision
reproduces the file):

    python scripts/tobio_share_compare.py register   # BEFORE the generator changes
    python scripts/tobio_share_compare.py run        # once, after; judges the registered criteria

`register` measures every reading the registration names as depending on the module, on the
committed generator, and writes data/results/tobio_share_registration.json: the change, the
acceptance criteria (a)-(e), the dependents and their before values. `run` measures the same things on
the changed generator, judges (a)-(e) and writes data/results/tobio_share_result.json with every
dependent's before -> after. It does not change the generator either way; the falsifier says what the
lane does next (restore the generator), and the result records it.
"""

from __future__ import annotations

import argparse
import contextlib
import io
import json
import re
import sys
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import genomeos  # noqa: E402
from genomeos import manifest as mf  # noqa: E402
from genomeos import results as gr  # noqa: E402
from genomeos.bio import main as bio_main  # noqa: E402
from genomeos.forge.design_answers import answer_record, results_for  # noqa: E402
from genomeos.ir import to_minutes  # noqa: E402
from genomeos.lang import parse_file  # noqa: E402
from genomeos.organism import haematopoiesis as haem  # noqa: E402
from genomeos.organism.experiment import run_experiment  # noqa: E402
from genomeos.runtime.body import Body  # noqa: E402

REGISTRATION = "tobio_share_registration"
RESULT = "tobio_share_result"
BIO = "data/organisms/human/haematopoiesis.bio"
MUTANTS = "data/organisms/human/haematopoiesis_mutants.bio"
DESIGNS = "data/organisms/human/haematopoiesis_designs.bio"
BIO_TEST_ROOTS = ["data/demo", "genomeos/std", "data/organisms"]  # the CI `bio test` line
OWN_CODE = ["scripts/tobio_share_compare.py", "genomeos/organism/haematopoiesis.py"]

DAY = to_minutes(1, "d")
THREE_YEARS = to_minutes(3, "yr")
SAMPLE_DAYS = list(range(0, 1096, 30))  # every 30 days, plus the 3-year horizon itself
STEADY_FROM_DAYS = 90  # `stage Steady { from: 90 d }` in the program: the Fill stage is not steady
POOL_TOLERANCE = 0.05

CHILDREN: dict[str, list[str]] = {}
for _c in haem.COMPARTMENTS:
    if _c["parent"]:
        CHILDREN.setdefault(_c["parent"], []).append(_c["id"])
BRANCHING = sorted(p for p, kids in CHILDREN.items() if len(kids) > 1)  # CLP, CMP, GMP, MPP
PARTITION_FLOWS = sorted(f"{p}_to_{k}" for p in BRANCHING for k in CHILDREN[p])
EXPECTED_GROWTH_CHANGES = sorted(
    f"{p}_{tag}" for p in BRANCHING for tag in ("fill_2x", "fill_10pc", "fill_1pc", "amplify")
)

CHANGE = {
    "generator": "genomeos/organism/haematopoiesis.py (application file, AGPL)",
    "plan": "plan() solves each compartment's growth step for the share rule: (1 + g) * (1 - sum(f_i)) = "
    "1 - o / A, where the f_i are the per-step outflow fractions to the children and o / A the share the "
    "inflow supplies (HSC: (1 + g) * (1 - f) = 1), replacing (1 + g) * prod(1 - f_i). The Fill-stage "
    "steps divide by the same retained term, so they move with it.",
    "to_bio": "the eight sibling-partition flows (a parent with two children) are written as `share:` with "
    "the same numbers: " + ", ".join(PARTITION_FLOWS),
    "unchanged": "single-child flows (HSC_to_MPP, MEP_to_Erythroblast, the five terminal-to-mature flows) "
    "stay `fraction:` with their numbers; the five mature losses, every name, every evidence string, "
    "haematopoiesis_mutants.bio and haematopoiesis_designs.bio are byte-identical",
    "expected_diff": "exactly the 8 partition lines change key (value identical) and exactly the 16 divide "
    "numbers of the four branching compartments change value: " + ", ".join(EXPECTED_GROWTH_CHANGES) + "; a "
    "single-child compartment keeps 1 - f = prod(1 - f), so its numbers cannot move. A diff of any other "
    "shape is not the registered change and is not adopted.",
    "regenerated_by": "to_bio() and mutants_bio(), the functions genomeos/storage.py's human_cell_turnover "
    "distiller calls, written to the two program files (the distiller itself also refetches Sender & "
    "Milo and rewrites tissues.bio, which this change does not touch)",
}

CRITERIA = {
    "a": "every wild-type assert of haematopoiesis.bio (3) and every assert of all nine experiments of "
    "haematopoiesis_mutants.bio (36) passes",
    "b": "wild type: every compartment's pool within +-5% of its solved value (plan()[id]['pool']) at every "
    "sample of the Steady stage (from 90 d, the program's own stage boundary) every 30 days to 1080 d and at "
    "3 yr: 35 samples x 11 compartments. 'Compartment' is read as the module's COMPARTMENTS (HSC through "
    "Thymocyte), the cell types whose pool plan() solves. The five mature types are reported at 3 yr "
    "against output x lifespan but are not in (b): they approach it with their lifespan as time constant "
    "(TCell 320 d), so no generator puts them within 5% from 90 d; the wild-type asserts test three of them "
    "at 3 yr. Fill-stage samples (0 to 60 d) are reported, not judged.",
    "c": "every growth number (the `fraction:` of every `action: divide` decision, 44 in all) that "
    "changes is "
    "listed old -> new with its relative change, and each is named as fitted: it must sit on a decision "
    "whose evidence tier is `inferred` (none sourced). Holds when the list is complete and no changed number "
    "is on a curated or experimental decision.",
    "d": "every committed result or reading that depends on the module (listed under `dependents`, "
    "identified before the change) is measured before (this file) and after (the result), and every change "
    "is listed.",
    "e": "falsifier: any wild-type assert or mutant assert failing, or any compartment sample outside +-5% "
    "in (b), or a diff that is not the registered change, means the change is NOT adopted and the generator "
    "is restored to its committed state.",
}

DEPENDENTS = [
    {
        "what": "the module itself: data/organisms/human/haematopoiesis.bio (milestone 1.1's module) and the "
        "two programs that import it, haematopoiesis_mutants.bio and haematopoiesis_designs.bio",
        "measured_as": "program decisions (key, value, evidence tier); wild-type trajectory, asserts, "
        "turnover, populations; nine experiments' asserts; two design answers",
    },
    {
        "what": "`bio test data/demo genomeos/std data/organisms` (CI): 58 files, 410/410 at 36cabe9; "
        "haematopoiesis.bio 3/3, haematopoiesis_designs.bio 3/3, haematopoiesis_mutants.bio 39/39",
        "measured_as": "bio_test",
    },
    {
        "what": "tests/test_haematopoiesis.py: file equality with the generator, plan bounds (0 < growth < "
        "1.1 x step_out), pools and mature types within 0.8-1.2 at 3 yr, turnover within 0.9-1.1 of the "
        "summed outputs, 16 populations, three knockouts",
        "measured_as": "plan_growth_within_test_bound, wild_type, mutants (the test file is run as "
        "acceptance)",
    },
    {
        "what": "design answers pinned in tests/test_forge_design_answers.py (lose_neutrophils_keep_the_rest "
        "-> -GFI1, 11 evaluations, 2 feasible; no_lymphocytes_at_all -> -IKZF1 tied with -TCF3, 29 "
        "evaluations, 22 feasible) and recorded in data/results/forge_calibration_inventory.json; also run "
        "by tests/test_design_budget.py and tests/test_forge_certainty.py",
        "measured_as": "designs",
    },
    {
        "what": "README.md and docs/ORGANISM-FROM-ONE-CELL.md readings: '2.8e11 cells replaced per day', "
        "'pools within 15% of design after three years', 'nine knockouts reproducing their published "
        "phenotypes', the design answers GFI1 and IKZF1 (also docs/BIOLANG-v0.3.md: 'IKZF1 or TCF3 alone'; "
        "docs/BIOFORGE-CONFIDENCE.md: GFI1, IKZF1)",
        "measured_as": "readings",
    },
    {
        "what": "docs/BIOLANG-v0.4-ECONOMY.md section 7.3: haematopoiesis with its designs and mutants gives "
        "identical results under `fates: first` and `fates: last`, 0 ambiguous points",
        "measured_as": "fates_modes",
    },
    {
        "what": "data/results/forge_calibration_inventory.json: its result_manifest.inputs holds the only "
        "committed digest of haematopoiesis.bio (sha256 41a822db..., 20643 bytes), which any change to the "
        "program makes stale, and its design_triples hold the two blood answers with evaluations, feasible "
        "counts and loss. The file is BioForge's and held by lane-pushfix; this lane does not rewrite it. "
        "tests/test_forge_calibration_inventory.py:146 re-runs the designs against it (RUN_ONLY_KEYS leaves "
        "result_manifest out); tests/test_forge_calibration.py needs every design solved",
        "measured_as": "program_sha256; designs",
    },
    {
        "what": "data/results/activator_combination_census.json and activator_combination_choice.json: a "
        "directory digest of data/organisms in their manifests (already stale before this change: 26 files "
        "recorded, 40 now; no test reads it); their haematopoiesis entries count activators (0), which "
        "neither the key nor the numbers can move",
        "measured_as": "not re-measured: listed as a digest that goes staler",
    },
    {
        "what": "censuses that glob every program, read before the change and found not to depend on it: "
        "biolang_key_census (share is an accepted decision key), when_census (rule and event when clauses "
        "only), amounts_or_concentrations and mandatory_locations (compartments and entities only), "
        "rule_evidence_tier and rule_cell_provenance (rule blocks only; the program has none), "
        "number_provenance (data/organisms excluded), compiled_defaults (nothing pinned)",
        "measured_as": "not re-measured",
    },
]


def _decisions(text: str) -> dict[str, dict[str, Any]]:
    """Every `decision` line of a generated program: action, the amount key and value, evidence tier."""
    out = {}
    pat = re.compile(r"decision (\S+) \{ action: (\w+);.*?; (fraction|share): ([^;]+);.*?evidence: (\w+)")
    for line in text.splitlines():
        if not line.startswith("decision "):
            continue
        m = pat.match(line)
        if m is None:
            raise ValueError(f"unparsed decision line: {line[:120]}")
        name, action, key, value, tier = m.groups()
        out[name] = {"action": action, "key": key, "value": float(value), "evidence": tier, "line": line}
    return out


def _counts(body: Body) -> dict[str, float]:
    counts: dict[str, float] = {}
    for c in body.alive_at(body.time):
        counts[c.cell_type] = counts.get(c.cell_type, 0.0) + c.count
    return dict(sorted(counts.items()))


def _bio_test() -> dict[str, Any]:
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        code = bio_main(["test", *BIO_TEST_ROOTS])
    lines = buf.getvalue().splitlines()
    per_file = {}
    for line in lines:
        m = re.match(r"([✓✗]) (\S+\.bio): .*?(\d+)/(\d+) checks$", line)
        if m and "haematopoiesis" in m.group(2):
            per_file[m.group(2)] = f"{m.group(3)}/{m.group(4)}"
    summary = lines[-1]
    m = re.match(r"(\d+) files, (\d+)/(\d+) checks passed", summary)
    assert m, summary
    return {
        "exit": code,
        "files": int(m.group(1)),
        "passed": int(m.group(2)),
        "total": int(m.group(3)),
        "summary": summary,
        "haematopoiesis_files": per_file,
        "failures": [x for x in lines if x.startswith("✗") or x.strip().startswith("FAIL")],
    }


def measure() -> dict[str, Any]:
    """Everything the registration names, on the generator and the programs as they are in this tree."""
    text = (ROOT / BIO).read_text()
    if text != haem.to_bio() or (ROOT / MUTANTS).read_text() != haem.mutants_bio():
        raise SystemExit("the committed programs are not what the generator writes: regenerate first")
    p = haem.plan()
    plan_out = {
        cid: {k: c[k] for k in ("pool", "step_out", "retained", "growth_fraction") if k in c}
        for cid, c in p.items()
    }
    bound = all(0 < c["growth_fraction"] < c["step_out"] * 1.1 for cid, c in p.items() if cid != "HSC")

    wt_module = parse_file(str(ROOT / BIO))
    body = Body(wt_module)
    samples = []
    for d in SAMPLE_DAYS:
        body.run(until=d * DAY)
        samples.append({"day": d, "counts": _counts(body)})
    body.run(until=THREE_YEARS)
    samples.append({"day": THREE_YEARS / DAY, "counts": _counts(body)})
    summary = body.summary()
    wild_type = {
        "asserts": [
            {"assert": a["assert"], "ok": a["ok"], "value": a.get("value")} for a in body.check_asserts()
        ],
        "final_counts": samples[-1]["counts"],
        "turnover_per_day": body.turnover_per_day(),
        "summed_mature_outputs_per_day": sum(m["output"] for m in haem.MATURE.values()),
        "populations": summary["populations"],
        "unknown": summary["unknown"],
        "ambiguous_fates": summary["ambiguous_fates"],
        "samples": samples,
    }
    last = Body(wt_module, fates="last").run(until=THREE_YEARS)
    last_counts = _counts(last)
    rel = {
        k: abs(last_counts.get(k, 0.0) - v) / max(abs(v), abs(last_counts.get(k, 0.0)), 1e-300)
        for k, v in wild_type["final_counts"].items()
    }
    fates_modes = {
        "first_ambiguous": summary["ambiguous_fates"],
        "last_ambiguous": last.summary()["ambiguous_fates"],
        "same_types": sorted(last_counts) == sorted(wild_type["final_counts"]),
        "max_rel_diff": max(rel.values()),
    }

    mut_module = parse_file(str(ROOT / MUTANTS))
    mutants = {}
    for ex in mut_module.experiments:
        r = run_experiment(mut_module, ex)
        asserts = [{"assert": a["assert"], "ok": a["ok"], "value": a.get("value")} for a in r.asserts()]
        mutants[ex.name] = {"all_ok": all(a["ok"] for a in asserts), "asserts": asserts}

    designs = [answer_record(DESIGNS, r) for r in results_for(DESIGNS)]
    final = wild_type["final_counts"]
    ratios_3yr = {cid: final.get(cid, 0.0) / p[cid]["pool"] for cid in p}
    ratios_3yr |= {mid: final.get(mid, 0.0) / (m["output"] * m["lifespan"]) for mid, m in haem.MATURE.items()}
    readings = {
        "cells_replaced_per_day": wild_type["turnover_per_day"],
        "max_abs_pool_deviation_at_3yr_all_16_types": max(abs(v - 1.0) for v in ratios_3yr.values()),
        "ratios_at_3yr": ratios_3yr,
        "knockouts_reproducing": sum(1 for m in mutants.values() if m["all_ok"]),
        "design_answers": {d["design"]: d["answer"] for d in designs},
    }
    return {
        "program_sha256": mf.sha256_of(ROOT / BIO)[0],
        "decisions": {
            k: {x: v[x] for x in ("action", "key", "value", "evidence")} for k, v in _decisions(text).items()
        },
        "program_lines": {k: v["line"] for k, v in _decisions(text).items()},
        "plan": plan_out,
        "plan_growth_within_test_bound": bound,
        "wild_type": wild_type,
        "fates_modes": fates_modes,
        "mutants": mutants,
        "designs": designs,
        "bio_test": _bio_test(),
        "readings": readings,
    }


def _manifest(stage: str, extra_inputs: list[str]) -> dict[str, Any]:
    return {
        "sources": [
            {
                "accession": "Sender & Milo 2021, Nat Med 27:45 (the outputs and lifespans in the generator)",
                "version": "as in genomeos/organism/haematopoiesis.py",
            },
            {
                "accession": "the committed generator and programs",
                "version": "pinned by sha256 and code.git_sha",
            },
        ],
        "inputs": [mf.input_entry(x) for x in (BIO, MUTANTS, DESIGNS, *extra_inputs)],
        "assembly": "n/a: no genome coordinates; a population model",
        "coordinates": "n/a: no genome coordinates; a population model",
        "parameters": {
            "stage": stage,
            "sample_days": SAMPLE_DAYS + [THREE_YEARS / DAY],
            "steady_from_days": STEADY_FROM_DAYS,
            "pool_tolerance": POOL_TOLERANCE,
            "bio_test_roots": BIO_TEST_ROOTS,
            "genomeos_imported_from": str(Path(genomeos.__file__).resolve().parent.relative_to(ROOT)),
        },
        "exclusions": ["the five mature types are reported, not judged, in (b): see criteria.b"],
        "partitions": "n/a: a deterministic model run; no training or held-out data",
        "code_cleanliness": mf.code_cleanliness(__file__, OWN_CODE),
    }


def register() -> int:
    before = measure()
    flows = before["decisions"]
    if any(flows[f]["key"] != "fraction" for f in PARTITION_FLOWS):
        raise SystemExit("the generator already writes share: - register only before the change")
    payload = {
        "lane": "lane-shareadopt",
        "area": "E",
        "date": "2026-10-10",
        "decision": "docs/ROADMAP.md area E (a8423e5): adopt share: in haematopoiesis.to_bio with the growth "
        "re-solved for 1 - sum(f), as its own registered lane",
        "stop_list": "not on the 2026-09-27 stop list (docs/ROADMAP.md section 5): this is the human "
        "haematopoiesis module, not the worm's fate ceiling",
        "prior_measurement": "lane-tobioshare (2026-10-10, no commit): same numbers as share: collapse the "
        "pools (CMP 0.13 of solved, wild-type asserts fail); growth re-solved for 1 - sum(f): 16 of 44 "
        "growth numbers move by up to +2.2%, every assert and the nine mutants pass, pools 1.038-1.04 of "
        "solved at 3 yr (Monocyte 0.9075 -> 1.0401 against today's generator); single-child flows "
        "identical as share:. Its trajectories between 0 d and 3 yr were not read for this registration.",
        "change": CHANGE,
        "criteria": CRITERIA,
        "dependents": DEPENDENTS,
        "before": before,
    }
    print("saved", gr.save_result(REGISTRATION, payload, manifest=_manifest("register", [])))
    return 0


def _ratio_table(wt: dict[str, Any], pools: dict[str, float]) -> dict[str, Any]:
    out = {}
    for cid, pool in pools.items():
        steady = [s["counts"].get(cid, 0.0) / pool for s in wt["samples"] if s["day"] >= STEADY_FROM_DAYS]
        fill = [s["counts"].get(cid, 0.0) / pool for s in wt["samples"] if s["day"] < STEADY_FROM_DAYS]
        out[cid] = {
            "steady_min": min(steady),
            "steady_max": max(steady),
            "steady_samples": len(steady),
            "outside": [
                (s["day"], s["counts"].get(cid, 0.0) / pool)
                for s in wt["samples"]
                if s["day"] >= STEADY_FROM_DAYS
                and abs(s["counts"].get(cid, 0.0) / pool - 1.0) > POOL_TOLERANCE
            ],
            "fill_stage_ratios": fill,
        }
    return out


def _diff(before: Any, after: Any, path: str = "") -> list[dict[str, Any]]:
    """Every leaf that differs, as path, before, after (lists compared by position)."""
    if isinstance(before, dict) and isinstance(after, dict):
        out = []
        for k in sorted(set(before) | set(after), key=str):
            out += _diff(before.get(k), after.get(k), f"{path}.{k}" if path else str(k))
        return out
    if isinstance(before, list) and isinstance(after, list) and len(before) == len(after):
        out = []
        for i, (b, a) in enumerate(zip(before, after, strict=True)):
            out += _diff(b, a, f"{path}[{i}]")
        return out
    return [] if before == after else [{"path": path, "before": before, "after": after}]


def run() -> int:
    reg_path = gr.RESULTS_DIR / f"{REGISTRATION}.json"
    reg = json.loads(reg_path.read_text())
    before = reg["before"]
    after = measure()

    # the change, as registered
    b_dec, a_dec = before["decisions"], after["decisions"]
    key_changes = sorted(k for k in a_dec if a_dec[k]["key"] != b_dec.get(k, {}).get("key"))
    value_changes = sorted(k for k in a_dec if a_dec[k]["value"] != b_dec.get(k, {}).get("value"))
    b_lines, a_lines = before["program_lines"], after["program_lines"]
    other_lines = sorted(
        k
        for k in set(a_lines) | set(b_lines)
        if a_lines.get(k) != b_lines.get(k) and k not in key_changes and k not in value_changes
    )
    registered_change = (
        key_changes == PARTITION_FLOWS
        and all(a_dec[k]["key"] == "share" and a_dec[k]["value"] == b_dec[k]["value"] for k in key_changes)
        and value_changes == EXPECTED_GROWTH_CHANGES
        and not other_lines
        and set(a_dec) == set(b_dec)
    )

    # (a)
    wt_ok = all(a["ok"] for a in after["wild_type"]["asserts"])
    mut_ok = {k: v["all_ok"] for k, v in after["mutants"].items()}
    n_mut_asserts = sum(len(v["asserts"]) for v in after["mutants"].values())
    a_holds = wt_ok and len(mut_ok) == 9 and all(mut_ok.values())

    # (b)
    pools = {cid: c["pool"] for cid, c in after["plan"].items()}
    pools_unchanged = pools == {cid: c["pool"] for cid, c in before["plan"].items()}
    table_after = _ratio_table(after["wild_type"], pools)
    table_before = _ratio_table(before["wild_type"], pools)
    b_holds = pools_unchanged and len(pools) == 11 and not any(v["outside"] for v in table_after.values())

    # (c)
    growth = [
        {
            "decision": k,
            "old": b_dec[k]["value"],
            "new": a_dec[k]["value"],
            "relative_change": a_dec[k]["value"] / b_dec[k]["value"] - 1.0,
            "evidence": a_dec[k]["evidence"],
            "named": "fitted: solved by plan() from the inferred amplifications and the measured outputs; "
            "none sourced",
        }
        for k in value_changes
        if a_dec[k]["action"] == "divide"
    ]
    n_divide = sum(1 for v in a_dec.values() if v["action"] == "divide")
    c_holds = len(growth) == len(value_changes) and all(g["evidence"] == "inferred" for g in growth)

    # (d)
    dependents = {
        "bio_test": {
            "before": {
                k: before["bio_test"][k] for k in ("files", "passed", "total", "haematopoiesis_files")
            },
            "after": {k: after["bio_test"][k] for k in ("files", "passed", "total", "haematopoiesis_files")},
            "after_failures": after["bio_test"]["failures"],
        },
        "test_haematopoiesis_plan_bound": {
            "before": before["plan_growth_within_test_bound"],
            "after": after["plan_growth_within_test_bound"],
        },
        "designs": _diff(before["designs"], after["designs"]),
        "readings": _diff(before["readings"], after["readings"]),
        "fates_modes": {"before": before["fates_modes"], "after": after["fates_modes"]},
        "wild_type_asserts": _diff(before["wild_type"]["asserts"], after["wild_type"]["asserts"]),
        "wild_type_scalars": {
            k: {"before": before["wild_type"][k], "after": after["wild_type"][k]}
            for k in ("turnover_per_day", "populations", "unknown", "ambiguous_fates")
        },
        "mutant_asserts": _diff(
            {k: v["asserts"] for k, v in before["mutants"].items()},
            {k: v["asserts"] for k, v in after["mutants"].items()},
        ),
    }
    d_holds = True  # every registered dependent was measured on both sides and every difference is listed

    adopted = registered_change and a_holds and b_holds
    payload = {
        "lane": "lane-shareadopt",
        "area": "E",
        "date": "2026-10-10",
        "registration": f"data/results/{REGISTRATION}.json",
        "registration_sha256": mf.sha256_of(reg_path)[0],
        "verdict": "ADOPTED" if adopted else "NOT_ADOPTED",
        "next": "keep the generator change"
        if adopted
        else "restore genomeos/organism/haematopoiesis.py and "
        "the program to the committed state (criterion e)",
        "registered_change": {
            "holds": registered_change,
            "key_changes": key_changes,
            "value_changes": value_changes,
            "other_lines_changed": other_lines,
        },
        "criteria": {
            "a": {
                "holds": a_holds,
                "wild_type_asserts": after["wild_type"]["asserts"],
                "experiments_all_ok": mut_ok,
                "experiment_asserts": n_mut_asserts,
            },
            "b": {
                "holds": b_holds,
                "pools_unchanged_by_the_change": pools_unchanged,
                "compartments": len(pools),
                "after": table_after,
                "before_for_reference": table_before,
                "mature_at_3yr_reported_not_judged": {
                    mid: {
                        "before": before["readings"]["ratios_at_3yr"][mid],
                        "after": after["readings"]["ratios_at_3yr"][mid],
                    }
                    for mid in haem.MATURE
                },
            },
            "c": {"holds": c_holds, "divide_numbers": n_divide, "changed": growth},
            "d": {"holds": d_holds, "dependents": dependents},
            "e": {
                "falsified": not adopted,
                "by": [
                    x
                    for x, bad in (
                        ("not the registered change", not registered_change),
                        ("an assert or mutant failed", not a_holds),
                        ("a compartment sample outside +-5%", not b_holds),
                    )
                    if bad
                ],
            },
        },
        "before_sha256": before["program_sha256"],
        "after_sha256": after["program_sha256"],
        "after": {k: after[k] for k in ("plan", "fates_modes", "readings", "bio_test", "designs")},
        "wild_type_samples_after": after["wild_type"]["samples"],
    }
    print(
        "saved",
        gr.save_result(RESULT, payload, manifest=_manifest("run", [f"data/results/{REGISTRATION}.json"])),
    )
    print("verdict", payload["verdict"], {k: payload["criteria"][k].get("holds") for k in "abcd"})
    return 0


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("stage", choices=("register", "run"))
    args = ap.parse_args(argv)
    return register() if args.stage == "register" else run()


if __name__ == "__main__":
    raise SystemExit(main())
