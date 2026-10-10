# SPDX-License-Identifier: AGPL-3.0-or-later
"""M2e, the registered run: signaltrace's bridge stage under the multiplicative composition.

    uv run python scripts/saturate_run.py

Runs only from the committed registration (data/results/saturate_registration.json) and refuses if any
registered field differs from genomeos/attribution/saturate.py, if its own code or the stage code it calls
differs from what the registration pinned, or if a pinned input's sha256 moved.

Order, as registered (saturate.STOP_RULE):
  1. the population is read and checked against signaltrace's pinned list;
  2. the anchor's columns are computed (crispri.annotate), signaltrace's stages and the multiplicative
     stages are computed for the fit rows and the population, one chromosome at a time; no label is read;
  3. signaltrace's harness, the multiplicative harness and the reproduction of signaltrace's pairs table
     are checked; a defect stops the run here; the census (a count, no label) is taken;
  4. the anchor is recomputed and must equal the committed arm-1 figures, and s2_target and
     s3_compiled_value must equal signaltrace's committed readings; a mismatch stops the run;
  5. only then is every multiplicative stage read and decided.

No model request, no network, no money: every input is on this machine and pinned by sha256.
"""

from __future__ import annotations

import json
import math
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
from genomeos.attribution import saturate as sa  # noqa: E402
from genomeos.attribution import signaltrace as st  # noqa: E402
from genomeos.evidence import COMPILED_DIR  # noqa: E402
from genomeos.results import RESULTS_DIR, save_result  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import s4joint_run as sjr  # noqa: E402
import signaltrace_run as sr  # noqa: E402
from runrules_census import load_rates  # noqa: E402

RESULT = "saturate_result"
REGISTRATION = RESULTS_DIR / "saturate_registration.json"
ST_REGISTRATION = RESULTS_DIR / "signaltrace_registration.json"
ST_RESULT = RESULTS_DIR / "signaltrace_result.json"
M2D_RESULT = RESULTS_DIR / "s4joint_result.json"
OWN_CODE = (
    "genomeos/attribution/saturate.py",
    "scripts/saturate_register.py",
    "scripts/saturate_run.py",
    "tests/test_saturate.py",
)
#: the code the run must equal, pinned by sha256 in the registration; the tests are not on the run's path
FROZEN = OWN_CODE[:3]
#: the committed code the stages run, not this lane's; pinned by sha256 in the registration
STAGE_CODE = (
    *sr.STAGE_CODE,
    "genomeos/attribution/signaltrace.py",
    "genomeos/attribution/combine2.py",
    "genomeos/attribution/s4joint.py",
    "scripts/signaltrace_run.py",
    "scripts/s4joint_run.py",
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
        mf.input_entry(M2D_RESULT, partition=None, role="M2d's committed census and readings, quoted"),
    ]


def check_registration() -> dict[str, Any]:
    """The committed registration; refused if a bound field or a pinned code file differs."""
    reg = json.loads(REGISTRATION.read_text())
    bad = [f for f, v in sa.registered_fields().items() if reg.get(f) != v]
    if bad:
        raise SystemExit(f"refused: the registration and saturate.py disagree on {', '.join(bad)}")
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
                " 9637584 and 1103719)",
                "version": "pinned by sha256",
            },
            {
                "accession": "data/results/s4joint_result.json (M2d, 86ddf73): its census and readings",
                "version": "pinned by sha256",
            },
        ],
        "inputs": entries,
        "input_count": len(entries),
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "cell": sa.CELL,
            "form": sa.FORM,
            "form_source": "Zhou et al. 2024, Cell Genomics 100672, doi 10.1016/j.xgen.2024.100672",
            "free_parameters": "none (saturate.PARAMETER_SOURCE)",
            "exponent_limit": sa.EXPONENT_LIMIT,
            "stages": {
                k: {"follows": v["follows"], "columns": list(st.columns(k))} for k, v in sa.STAGES.items()
            },
            "chain": list(sa.CHAIN),
            "branches": list(sa.BRANCHES),
            "fit": "crispri.logistic_fit, lam 1e-3, 25 rounds, on the covered K562 training pairs",
            "bootstrap": sa.BOOTSTRAP,
            "closed_form_tolerance": sa.CLOSED_FORM_TOLERANCE,
            "rate_tier": "gene_rates.json, else its genome median (borrowed)",
            "selection": "compile.SELECTION_READER_CELLS on the predicted-only program, as signaltrace's s4",
        },
        "exclusions": [
            sa.PREDICTED_ONLY,
            "held-out pairs in cell types other than K562: the anchor is a K562 result",
            "a pair a stage drops is kept in the population with 0 in that stage's two columns",
            "a gene combine2.classify excludes or the multiplicative fit leaves not runnable: its pairs are"
            f" dropped at {sa.S4}",
        ],
        "partitions": {
            "training": "the fit rows only: the covered K562 training pairs",
            "heldout": "the scored population: the held-out file's K562 pairs; no held-out pair enters a fit",
        },
        "code_cleanliness": mf.code_cleanliness(__file__, OWN_CODE),
    }


def _summed_beta(rho: dict[str, float]) -> float:
    return sum(-math.log(r) for r in rho.values() if r < 1.0)


def descriptives(joints: dict[str, sj.JointLookup], st_result: dict[str, Any], m2d: dict[str, Any]) -> dict:
    genes = [g for j in joints.values() for g in j.genes.values()]
    admitted = [g for g in genes if not g.excluded]
    ok = [g for g in genes if g.forms[sa.FORM]["ok"]]
    q = sjr._quantiles
    return {
        "summed_drop_quantiles_admitted": q([sa.summed_drop(g.rho) for g in admitted]),
        "summed_beta_quantiles_admitted": q([_summed_beta(g.rho) for g in admitted]),
        "rule_count_quantiles_reached": q([float(len(g.elements)) for g in genes]),
        "rule_count_quantiles_identifiable": q([float(len(g.elements)) for g in ok]),
        "identifiable_by_tier": dict(sorted(Counter(g.forms[sa.FORM]["tier"] for g in ok).items())),
        "gene_runs": sum(g.forms[sa.FORM]["gene_runs"] for g in ok),
        "k562_hold_runs": sum(g.forms[sa.FORM]["hold_runs"] for g in ok),
        "max_abs_run_minus_closed_form_log2": max(
            (g.forms[sa.FORM]["max_abs_run_minus_closed_form_log2"] for g in ok), default=None
        ),
        "per_chromosome": {c: dict(j.counts) for c, j in joints.items()},
        "context_verbatim": {
            "signaltrace_source": "data/results/signaltrace_result.json",
            "signaltrace_s4_bridge": st_result["readings_verbatim"]["s4_bridge"],
            "signaltrace_chain_outcome": st_result["outcome"],
            "signaltrace_chain_reading": st_result["reading"],
            "m2d_source": "data/results/s4joint_result.json",
            "m2d_per_form": {
                f: {k: m2d["per_form"][f][k] for k in ("bridge", "bridge_reading", "chain", "chain_reading")}
                for f in ("mean", "competition")
            },
        },
    }


def pairs_table(pop: list[crispri.Pair]) -> dict[str, Any]:
    stages = list(sa.MULT_STAGES)
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
        "lane": sa.LANE,
        "registration": "data/results/saturate_registration.json",
        "registration_date": reg.get("date"),
        "question": sa.QUESTION,
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
            joints[chrom] = sa.joint_lookup(
                chrom, text, lookups[chrom], rates, borrowed, prov, lent, readers, RESULTS_DIR
            )
        except ValueError as e:
            return defect(f"{chrom}: the multiplicative fit could not read a reached gene: {e}")
        del text
        print(chrom, dict(joints[chrom].counts), f"{time.time() - t0:.0f} s", flush=True)
    whys: dict[str, dict[str, str]] = {}
    for p in training + pop:
        f, why = st.trace_pair(p, lookups.get(p.chrom))
        g, why2 = sa.trace(p, lookups.get(p.chrom), joints.get(p.chrom))
        p.features.update(f)
        p.features.update(g)
        whys[st.pair_id(p)] = {**why, **why2}
    st_result = json.loads(ST_RESULT.read_text())
    m2d = json.loads(M2D_RESULT.read_text())
    problems = st.harness(lookups, training + pop) + sa.harness(lookups, joints, training + pop)
    problems += sjr.reproduction_mismatch(pop, st_result["pairs_table"])
    if problems:
        return defect(f"{len(problems)} harness problems: " + "; ".join(problems[:20]))
    census = sa.census(joints, m2d)
    print("census:", census["sentence"], flush=True)
    anchor = st.anchor_reading(fit_rows, pop)
    if not anchor["reproduced"]:
        return defect(
            f"the anchor did not reproduce: recomputed {anchor['recomputed']}, committed {st.ANCHOR}"
        )
    comp_rows, comp_scores, comp_weights = sj.comparators(fit_rows, pop, anchor)
    bad = sj.comparator_mismatch(comp_rows, st_result["stages"])
    if bad:
        return defect("the comparators did not reproduce signaltrace's readings: " + "; ".join(bad[:10]))

    scored = sa.score(fit_rows, pop, anchor, comp_rows, comp_scores, comp_weights)
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
        "census": census,
        "outcome": {"bridge": scored["bridge"], "chain": scored["chain"]},
        "reading": scored["reading"],
        "bridge_reading": scored["bridge_reading"],
        "chain_reading": scored["chain_reading"],
        "stage_where_it_stops": scored["stage_where_it_stops"],
        "readings_verbatim": {s: r["reading"] for s, r in scored["stages"].items()},
        "stages": scored["stages"],
        "chain": list(sa.CHAIN),
        "branches": list(sa.BRANCHES),
        "by_construction": list(sa.BY_CONSTRUCTION),
        "descriptives_not_verdicts": descriptives(joints, st_result, m2d),
        "limitations": list(sa.LIMITATIONS),
        "pairs_table": pairs_table(pop),
        "seconds": round(time.time() - t0, 1),
    }
    print("saved", save_result(RESULT, out, manifest=manifest([*entries, own])))
    summary = {"census": census["sentence"], "outcome": out["outcome"], "reading": out["reading"]}
    print(json.dumps(summary, indent=1))
    for s, r in out["readings_verbatim"].items():
        print(f"{s}: {r}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
