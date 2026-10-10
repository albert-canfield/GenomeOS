# SPDX-License-Identifier: AGPL-3.0-or-later
"""M2d, the registered run: signaltrace's bridge stage re-run with combine2's joint fit.

    uv run python scripts/s4joint_run.py

Runs only from the committed registration (data/results/s4joint_registration.json) and refuses if any
registered field differs from genomeos/attribution/s4joint.py, if its own code or the stage code it calls
differs from what the registration pinned, or if a pinned input's sha256 moved.

Order, as registered (s4joint.STOP_RULE):
  1. the population is read and checked against signaltrace's pinned list; combine2's registration is
     checked against combine2.py;
  2. the anchor's columns are computed (crispri.annotate), signaltrace's stages and the joint stages are
     computed for the fit rows and the population, one chromosome at a time; no label is read;
  3. signaltrace's harness, the joint harness and the reproduction of signaltrace's pairs table are
     checked; a defect stops the run here;
  4. the anchor is recomputed and must equal the committed arm-1 figures, and s2_target and
     s3_compiled_value must equal signaltrace's committed readings; a mismatch stops the run;
  5. only then is every joint stage read and decided.

No model request, no network, no money: every input is on this machine and pinned by sha256.
"""

from __future__ import annotations

import json
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import context_evidence as ce  # noqa: E402
from genomeos.attribution import crispri  # noqa: E402
from genomeos.attribution import regulators as regs  # noqa: E402
from genomeos.attribution import s4joint as sj  # noqa: E402
from genomeos.attribution import signaltrace as st  # noqa: E402
from genomeos.evidence import COMPILED_DIR  # noqa: E402
from genomeos.results import RESULTS_DIR, save_result  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import signaltrace_run as sr  # noqa: E402
from runrules_census import load_rates  # noqa: E402

RESULT = "s4joint_result"
REGISTRATION = RESULTS_DIR / "s4joint_registration.json"
ST_REGISTRATION = RESULTS_DIR / "signaltrace_registration.json"
ST_RESULT = RESULTS_DIR / "signaltrace_result.json"
C2_REGISTRATION = RESULTS_DIR / "combine2_registration.json"
OWN_CODE = (
    "genomeos/attribution/s4joint.py",
    "scripts/s4joint_register.py",
    "scripts/s4joint_run.py",
    "tests/test_s4joint.py",
)
#: the code the run must equal, pinned by sha256 in the registration; the tests are not on the run's path
FROZEN = OWN_CODE[:3]
#: the committed code the stages run, not this lane's; pinned by sha256 in the registration
STAGE_CODE = (
    *sr.STAGE_CODE,
    "genomeos/attribution/signaltrace.py",
    "genomeos/attribution/combine2.py",
    "scripts/signaltrace_run.py",
)


def inputs(chroms: list[str], training: list[crispri.Pair], pop: list[crispri.Pair]) -> list[dict[str, Any]]:
    """Every file the run reads but its own registration, pinned by sha256 (the registration pins these)."""
    return [
        *sr.inputs(chroms, training, pop),
        mf.input_entry(
            ST_REGISTRATION, partition=None, role="signaltrace's registration: the population ids"
        ),
        mf.input_entry(
            ST_RESULT, partition=None, role="signaltrace's committed readings and pairs table: the comparator"
        ),
        mf.input_entry(C2_REGISTRATION, partition=None, role="combine2's registration: the two forms"),
    ]


def check_registration() -> dict[str, Any]:
    """The committed registration; refused if a bound field or a pinned code file differs."""
    reg = json.loads(REGISTRATION.read_text())
    bad = [f for f, v in sj.registered_fields().items() if reg.get(f) != v]
    if bad:
        raise SystemExit(f"refused: the registration and s4joint.py disagree on {', '.join(bad)}")
    moved = [p for p, h in sr.code_sha256(FROZEN).items() if reg["frozen_code"]["sha256"].get(p) != h]
    moved += [p for p, h in sr.code_sha256(STAGE_CODE).items() if reg["stage_code"]["sha256"].get(p) != h]
    if moved:
        raise SystemExit(f"refused: code differs from what the registration pinned: {moved}")
    return reg


def manifest(entries: list[dict[str, Any]]) -> dict[str, Any]:
    return {
        "sources": [
            {
                "accession": "EngreitzLab/CRISPR_comparison EPCrisprBenchmark training_K562 and"
                " heldout_5_cell_types (Gschwind et al.)",
                "version": "fetched 2026-09-16; pinned here by sha256",
                "url": crispri.BASE_URL,
            },
            {
                "accession": "ENCODE SCREEN cCREs scored by AlphaGenome deletion (all-element table and"
                " per-element response cache)",
                "version": "the 2026-09 all-element sweep; pinned here by sha256",
            },
            {
                "accession": "the compiled non-coding programs (data/knowledge/compiled)",
                "version": "pinned by sha256",
            },
            {
                "accession": "data/results/gene_rates.json (mouse NIH 3T3, e921fae)",
                "version": "pinned by sha256",
            },
            {
                "accession": regs.CLUSTER_SOURCE["name"],
                "version": regs.CLUSTER_SOURCE["last_modified"],
                "url": regs.CLUSTER_SOURCE["url"],
            },
            {
                "accession": "ENCODE DNase-seq narrowPeak sets read by reader v1",
                "version": "pinned by sha256",
            },
            {
                "accession": "data/results/signaltrace_registration.json and signaltrace_result.json (M1c,"
                " 9637584 and 1103719, promoted d1df762)",
                "version": "pinned by sha256",
            },
            {
                "accession": "data/results/combine2_registration.json (lane-combine2, the two forms)",
                "version": "pinned by sha256",
            },
        ],
        "inputs": entries,
        "input_count": len(entries),
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "cell": sj.CELL,
            "forms": list(sj.FORMS),
            "theta": sj.THETA,
            "stages": {
                k: {"follows": v["follows"], "columns": list(st.columns(k))} for k, v in sj.STAGES.items()
            },
            "chains": {f: list(c) for f, c in sj.CHAINS.items()},
            "branches": {f: list(b) for f, b in sj.BRANCHES.items()},
            "fit": "crispri.logistic_fit, lam 1e-3, 25 rounds, on the covered K562 training pairs",
            "bootstrap": sj.BOOTSTRAP,
            "closed_form_tolerance": sj.CLOSED_FORM_TOLERANCE,
            "rate_tier": "gene_rates.json, else its genome median (borrowed)",
            "selection": "compile.SELECTION_READER_CELLS on the predicted-only program, as signaltrace's s4",
        },
        "exclusions": [
            sj.PREDICTED_ONLY,
            "held-out pairs in cell types other than K562: the anchor is a K562 result",
            "a pair a stage drops is kept in the population with 0 in that stage's two columns",
            "a gene combine2.classify excludes or a form's fit leaves not runnable: its pairs are dropped at"
            " s4_joint for that form",
        ],
        "partitions": {
            "training": "the fit rows only: the covered K562 training pairs",
            "heldout": "the scored population: the held-out file's K562 pairs; no held-out pair enters a fit",
        },
        "code_cleanliness": mf.code_cleanliness(__file__, OWN_CODE),
    }


def reproduction_mismatch(pop: list[crispri.Pair], committed: dict[str, Any]) -> list[str]:
    """Population pairs whose recomputed signaltrace columns differ from signaltrace's pairs table."""
    cols = committed["columns"]
    rows = {r[0]: r for r in committed["rows"]}
    stage_cols = [(i, c) for i, c in enumerate(cols) if c not in ("pair", "regulated")]
    bad = []
    for p in pop:
        row = rows.get(st.pair_id(p))
        if row is None:
            bad.append(f"{st.pair_id(p)}: not in signaltrace's pairs table")
            continue
        for i, c in stage_cols:
            if round(p.features[c], 6) != row[i]:
                bad.append(f"{st.pair_id(p)} {c}: {round(p.features[c], 6)} against {row[i]}")
    if len(rows) != len(pop):
        bad.append(f"signaltrace's pairs table has {len(rows)} pairs, the population {len(pop)}")
    return bad


def _quantiles(xs: list[float]) -> dict[str, float]:
    s = sorted(xs)
    return {str(q): s[min(len(s) - 1, int(q * len(s)))] for q in (0.1, 0.25, 0.5, 0.75, 0.9)} if s else {}


def descriptives(joints: dict[str, sj.JointLookup], st_result: dict[str, Any]) -> dict[str, Any]:
    genes = [g for j in joints.values() for g in j.genes.values()]
    per_form = {}
    for form in sj.FORMS:
        ok = [g for g in genes if g.forms[form]["ok"]]
        per_form[form] = {
            "reached_genes": len(genes),
            "runnable": len(ok),
            "not_runnable_by_reason": dict(
                Counter(g.forms[form]["reason"] for g in genes if not g.forms[form]["ok"]).most_common()
            ),
            "runnable_by_tier": dict(Counter(g.forms[form]["tier"] for g in ok)),
            "q_quantiles_runnable": _quantiles([g.forms[form]["q"] for g in ok]),
            "rule_count_quantiles_runnable": _quantiles([float(len(g.elements)) for g in ok]),
            "k562_hold_runs": sum(g.forms[form]["hold_runs"] for g in ok),
            "max_abs_run_minus_closed_form_log2": max(
                (g.forms[form]["max_abs_run_minus_closed_form_log2"] for g in ok), default=None
            ),
        }
    return {
        "per_form": per_form,
        "reached_genes_rule_count_quantiles": _quantiles([float(len(g.elements)) for g in genes]),
        "per_chromosome": {c: dict(j.counts) for c, j in joints.items()},
        "signaltrace_context_verbatim": {
            "source": "data/results/signaltrace_result.json",
            "s4_bridge": st_result["readings_verbatim"]["s4_bridge"],
            "chain_outcome": st_result["outcome"],
            "chain_reading": st_result["reading"],
        },
    }


def pairs_table(pop: list[crispri.Pair]) -> dict[str, Any]:
    stages = list(sj.JOINT_STAGES)
    cols = ["pair", "regulated", *[f"{s}_{c}" for s in stages for c in ("kept", "drop")]]
    rows = [
        [
            st.pair_id(p),
            p.regulated,
            *[round(p.features[f"{s}_{c}"], 6) for s in stages for c in ("kept", "drop")],
        ]
        for p in pop
    ]
    return {"columns": cols, "rows": rows}


def main() -> int:
    t0 = time.time()
    if (RESULTS_DIR / f"{RESULT}.json").exists():
        raise SystemExit("refused: the registered run has already been written")
    reg = check_registration()
    training, pop = sr.load_pairs()
    chroms = sr.chromosomes(training + pop)
    entries = inputs(chroms, training, pop)
    own = mf.input_entry(REGISTRATION, partition=None, role="this run's registration")
    base: dict[str, Any] = {
        "lane": sj.LANE,
        "registration": "data/results/s4joint_registration.json",
        "registration_date": reg.get("date"),
        "question": sj.QUESTION,
        "alphagenome_requests": 0,
        "money": "none: no request is sent and no key is read",
    }

    def defect(what: str) -> int:
        out = {
            **base,
            "outcome": "HARNESS_DEFECT",
            "reading": st.reading("HARNESS_DEFECT", what=what),
            "seconds": round(time.time() - t0, 1),
        }
        print("saved", save_result(RESULT, out, manifest=manifest([*entries, own])))
        print(out["reading"])
        return 1

    pinned = {e["path"]: e["sha256"] for e in reg["result_manifest"]["inputs"]}
    moved = [e["path"] for e in entries if pinned.get(e["path"]) != e["sha256"]]
    if moved:
        return defect(f"pinned inputs moved since the registration: {moved}")
    st_reg = json.loads(ST_REGISTRATION.read_text())
    ids = st.population_ids(pop)
    if (
        st.ids_sha256(ids) != reg["population_record"]["sha256"]
        or reg["population_record"]["sha256"] != st_reg["population_sha256"]
        or ids != st_reg["population_ids"]
    ):
        return defect("the population differs from signaltrace's pinned list")
    if [f for f, v in st.registered_fields().items() if st_reg.get(f) != v]:
        return defect("signaltrace's registration and signaltrace.py disagree")
    c2_bad = sj.combine2_mismatch(json.loads(C2_REGISTRATION.read_text()))
    if c2_bad:
        return defect(f"combine2's registration and combine2.py disagree on {c2_bad}")

    crispri.annotate(training + pop, crispri.DeletionTable(), crispri.ElementCache())
    fit_rows = [p for p in training if p.covered]
    rates, borrowed, prov, lent = load_rates()
    clusters = regs.load_clusters(regs.CLUSTERS, chroms)
    readers = ce.Readers(RESULTS_DIR)
    by_chrom: dict[str, list[crispri.Pair]] = defaultdict(list)
    for p in training + pop:
        by_chrom[p.chrom].append(p)
    lookups: dict[str, st.ChromLookup] = {}
    joints: dict[str, sj.JointLookup] = {}
    for chrom in chroms:
        text = (Path(COMPILED_DIR) / f"noncoding_{chrom}.bio").read_text()
        lookups[chrom] = st.chromosome_lookup(
            chrom, text, by_chrom[chrom], rates, borrowed, prov, lent, clusters, readers, RESULTS_DIR
        )
        try:
            joints[chrom] = sj.joint_lookup(
                chrom, text, lookups[chrom], rates, borrowed, prov, lent, readers, RESULTS_DIR
            )
        except ValueError as e:
            return defect(f"{chrom}: the joint fit could not read a reached gene: {e}")
        del text
        print(chrom, dict(joints[chrom].counts), f"{time.time() - t0:.0f} s", flush=True)
    whys: dict[str, dict[str, str]] = {}
    for p in training + pop:
        f, why = st.trace_pair(p, lookups.get(p.chrom))
        g, why2 = sj.trace_joint(p, lookups.get(p.chrom), joints.get(p.chrom))
        p.features.update(f)
        p.features.update(g)
        whys[st.pair_id(p)] = {**why, **why2}
    st_result = json.loads(ST_RESULT.read_text())
    problems = st.harness(lookups, training + pop) + sj.harness(lookups, joints, training + pop)
    problems += reproduction_mismatch(pop, st_result["pairs_table"])
    if problems:
        return defect(f"{len(problems)} harness problems: " + "; ".join(problems[:20]))
    anchor = st.anchor_reading(fit_rows, pop)
    if not anchor["reproduced"]:
        return defect(
            f"the anchor did not reproduce: recomputed {anchor['recomputed']}, committed {st.ANCHOR}"
        )
    comp_rows, comp_scores, comp_weights = sj.comparators(fit_rows, pop, anchor)
    bad = sj.comparator_mismatch(comp_rows, st_result["stages"])
    if bad:
        return defect("the comparators did not reproduce signaltrace's readings: " + "; ".join(bad[:10]))

    scored = sj.score_forms(fit_rows, pop, anchor, comp_rows, comp_scores, comp_weights)
    pop_whys = [whys[st.pair_id(p)] for p in pop]
    for stage, row in scored["stages"].items():
        row["population_dropped_because"] = st.dropped_because(pop_whys, stage)
    out = {
        **base,
        "population": {
            "pairs": len(pop),
            "positives": sum(p.regulated for p in pop),
            "sha256": reg["population_record"]["sha256"],
            "chromosomes": len({p.chrom for p in pop}),
        },
        "fit_rows": {"pairs": len(fit_rows), "positives": sum(p.regulated for p in fit_rows)},
        "anchor": {k: anchor[k] for k in ("recomputed", "committed", "reproduced")},
        "comparators_reproduced": True,
        "outcome": {f: {"bridge": v["bridge"], "chain": v["chain"]} for f, v in scored["per_form"].items()},
        "reading": scored["reading"],
        "per_form": scored["per_form"],
        "readings_verbatim": {s: r["reading"] for s, r in scored["stages"].items()},
        "stages": scored["stages"],
        "chains": {f: list(c) for f, c in sj.CHAINS.items()},
        "branches": {f: list(b) for f, b in sj.BRANCHES.items()},
        "descriptives_not_verdicts": descriptives(joints, st_result),
        "limitations": list(sj.LIMITATIONS),
        "pairs_table": pairs_table(pop),
        "seconds": round(time.time() - t0, 1),
    }
    print("saved", save_result(RESULT, out, manifest=manifest([*entries, own])))
    print(json.dumps({"outcome": out["outcome"], "reading": out["reading"]}, indent=1))
    for s, r in out["readings_verbatim"].items():
        print(f"{s}: {r}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
