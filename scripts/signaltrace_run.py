# SPDX-License-Identifier: AGPL-3.0-or-later
"""M1c, the registered run: the held-out K562 CRISPRi pairs traced through each compiled stage.

    uv run python scripts/signaltrace_run.py

Runs only from the committed registration (data/results/signaltrace_registration.json) and refuses if
any registered field differs from genomeos/attribution/signaltrace.py, if its own code or the stage code
it calls differs from what the registration pinned, or if a pinned input's sha256 moved.

Order, as registered (signaltrace.STOP_RULE):
  1. the population is read and its sha256 checked against the registration;
  2. the anchor's columns are computed (crispri.annotate) and every compiled stage is computed for the
     fit rows and the population, one chromosome at a time; no label is read;
  3. the predictor-side harness (signaltrace.harness) is checked; a defect stops the run here;
  4. the anchor is recomputed and must equal the committed arm-1 figures; a mismatch stops the run;
  5. only then is every stage beyond the anchor read and decided.

No model request, no network, no money: every input is on this machine and pinned by sha256.
"""

from __future__ import annotations

import hashlib
import json
import math
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import compile as cp  # noqa: E402
from genomeos.attribution import context_evidence as ce  # noqa: E402
from genomeos.attribution import crispri  # noqa: E402
from genomeos.attribution import regulators as regs  # noqa: E402
from genomeos.attribution import signaltrace as st  # noqa: E402
from genomeos.attribution.measured import CRISPRI_SPLIT_OF  # noqa: E402
from genomeos.evidence import COMPILED_DIR  # noqa: E402
from genomeos.genome import reader as rd  # noqa: E402
from genomeos.results import RESULTS_DIR, save_result  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from runrules_census import RATE_TABLE, load_rates  # noqa: E402

RESULT = "signaltrace_result"
REGISTRATION = RESULTS_DIR / "signaltrace_registration.json"
ANCHOR_RESULT = RESULTS_DIR / "crispri_published_v2.json"
CHROMS = [f"chr{c}" for c in [*range(1, 23), "X", "Y"]]
OWN_CODE = (
    "genomeos/attribution/signaltrace.py",
    "scripts/signaltrace_register.py",
    "scripts/signaltrace_run.py",
    "tests/test_signaltrace.py",
)
#: the code the run must equal, pinned by sha256 in the registration; the tests are not on the run's path
FROZEN = OWN_CODE[:3]
#: the committed code the stages run, not this lane's; pinned by sha256 in the registration
STAGE_CODE = (
    "genomeos/attribution/crispri.py",
    "genomeos/attribution/compile.py",
    "genomeos/attribution/bridge.py",
    "genomeos/attribution/readergate.py",
    "genomeos/attribution/regulators.py",
    "genomeos/attribution/runrules.py",
    "genomeos/attribution/context_evidence.py",
    "genomeos/runtime/grn.py",
    "genomeos/ir/model.py",
    "genomeos/lang/parser.py",
    "scripts/runrules_census.py",
)


def load_pairs() -> tuple[list[crispri.Pair], list[crispri.Pair]]:
    training, heldout = crispri.load(crispri.TRAINING), crispri.load(crispri.HELDOUT)
    return training, [p for p in heldout if p.cell == st.CELL]


def chromosomes(pairs: list[crispri.Pair]) -> list[str]:
    seen = {p.chrom for p in pairs}
    return [c for c in CHROMS if c in seen]


def code_sha256(paths: tuple[str, ...]) -> dict[str, str]:
    return {p: hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in paths}


def peak_files(chroms: list[str]) -> list[str]:
    """The reader's DNase-seq files gate_program reads for these chromosomes, where present."""
    out = []
    for b, _ in cp.reader_cells():
        for c in chroms:
            path = RESULTS_DIR / rd.peaks_path(b, c).name
            if path.exists():
                out.append(str(path))
    return out


def inputs(chroms: list[str], training: list[crispri.Pair], pop: list[crispri.Pair]) -> list[dict[str, Any]]:
    """Every file the run reads, pinned by sha256 (the registration pins the same list)."""
    out = [
        mf.input_entry(
            crispri.KNOWLEDGE / crispri.TRAINING,
            partition=CRISPRI_SPLIT_OF[crispri.TRAINING],
            pairs=len(training),
            role="fit rows: the covered K562 training pairs",
        ),
        mf.input_entry(
            crispri.KNOWLEDGE / crispri.HELDOUT,
            partition=CRISPRI_SPLIT_OF[crispri.HELDOUT],
            pairs=len(pop),
            role="the population: its K562 pairs, scored",
        ),
        mf.input_entry(crispri.ELEMENTS, partition=None, role="DeletionTable: elements by overlap (anchor)"),
        mf.input_entry(
            crispri.ELEMENT_CACHE, partition=None, role="ElementCache: every gene's deletion value"
        ),
        mf.input_entry(ANCHOR_RESULT, partition=None, role="the anchor's committed arm-1 figures"),
        mf.input_entry(RATE_TABLE, partition=None, role="the bridge's rates (borrowed tier)"),
        mf.input_entry(regs.CLUSTERS, partition=None, role="ENCODE TF ChIP-seq peak clusters (s6_regulator)"),
    ]
    out += [
        mf.input_entry(Path(COMPILED_DIR) / f"noncoding_{c}.bio", partition=None, role="compiled program")
        for c in chroms
    ]
    out.append(
        mf.files_entry(
            "reader v1 DNase-seq narrowPeak sets, thirteen biosamples (compile.gate_program's reader)",
            peak_files(chroms),
            partition=None,
        )
    )
    return out


def check_registration() -> dict[str, Any]:
    """The committed registration; refused if a bound field or a pinned code file differs."""
    reg = json.loads(REGISTRATION.read_text())
    bad = [f for f, v in st.registered_fields().items() if reg.get(f) != v]
    if bad:
        raise SystemExit(f"refused: the registration and signaltrace.py disagree on {', '.join(bad)}")
    moved = [p for p, h in code_sha256(FROZEN).items() if reg["frozen_code"]["sha256"].get(p) != h]
    moved += [p for p, h in code_sha256(STAGE_CODE).items() if reg["stage_code"]["sha256"].get(p) != h]
    if moved:
        raise SystemExit(f"refused: code differs from what the registration pinned: {moved}")
    return reg


def pairs_table(pop: list[crispri.Pair]) -> dict[str, Any]:
    stages = list(st.STAGES)
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


def closed_form(run: dict[str, Any]) -> float:
    return -math.log2(run["rho"]) if run["state"] == "open" else 0.0


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
        ],
        "inputs": entries,
        "input_count": len(entries),
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "cell": st.CELL,
            "stages": {
                k: {"follows": v["follows"], "columns": list(st.columns(k))} for k, v in st.STAGES.items()
            },
            "chain": list(st.CHAIN),
            "branches": list(st.BRANCHES),
            "fit": "crispri.logistic_fit, lam 1e-3, 25 rounds, on the covered K562 training pairs",
            "bootstrap": st.BOOTSTRAP,
            "closed_form_tolerance": st.CLOSED_FORM_TOLERANCE,
            "bridge_tier": "borrowed: gene_rates.json and its genome median",
            "selection": "compile.SELECTION_READER_CELLS for s4_bridge, s5_reader, s6_regulator and s7_run;"
            " the predicted program as compiled for s4a_bridge_as_compiled",
        },
        "exclusions": [
            st.PREDICTED_ONLY,
            "held-out pairs in cell types other than K562: the anchor is a K562 result",
            "a pair a stage drops is kept in the population with 0 in that stage's two columns",
        ],
        "partitions": "training: the fit rows only (covered K562 training pairs); heldout: the scored"
        " population (its K562 pairs); no held-out pair enters a fit",
        "code_cleanliness": mf.code_cleanliness(__file__, OWN_CODE),
    }


def descriptives(lookups: dict[str, st.ChromLookup], pop: list[crispri.Pair]) -> dict[str, Any]:
    reasons_reader: Counter = Counter()
    reasons_compiled: Counter = Counter()
    runs: list[dict[str, Any]] = []
    for look in lookups.values():
        reasons_reader.update(look.why_reader.values())
        reasons_compiled.update(look.why_compiled.values())
        runs += list(look.runs.values())
    measured = [
        p for p in pop if p.chrom in lookups and lookups[p.chrom].measured_names(p.start, p.end, p.gene)
    ]
    return {
        "bridge_unresolved_reasons_reader_cells": dict(reasons_reader.most_common()),
        "bridge_unresolved_reasons_as_compiled": dict(reasons_compiled.most_common()),
        "rules_reached": sum(look.counts["rules_reached"] for look in lookups.values()),
        "runs": len(runs),
        "run_states_k562": dict(Counter(r["state"] for r in runs)),
        "run_max_abs_difference_from_closed_form_log2": max(
            (abs(r["change"] - closed_form(r)) for r in runs), default=None
        ),
        "experimental_rule_lines_removed": sum(
            look.counts["experimental_rule_lines_removed"] for look in lookups.values()
        ),
        "population_pairs_the_experimental_layer_names_in_k562": len(measured),
        "of_them_regulated": sum(p.regulated for p in measured),
        "per_chromosome": {c: dict(look.counts) for c, look in lookups.items()},
    }


def main() -> int:
    t0 = time.time()
    if (RESULTS_DIR / f"{RESULT}.json").exists():
        raise SystemExit("refused: the registered run has already been written")
    reg = check_registration()
    training, pop = load_pairs()
    chroms = chromosomes(training + pop)
    entries = inputs(chroms, training, pop)
    base: dict[str, Any] = {
        "lane": st.LANE,
        "registration": "data/results/signaltrace_registration.json",
        "registration_date": reg.get("date"),
        "question": st.QUESTION,
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
        print("saved", save_result(RESULT, out, manifest=manifest(entries)))
        print(out["reading"])
        return 1

    pinned = {e["path"]: e["sha256"] for e in reg["result_manifest"]["inputs"]}
    moved = [e["path"] for e in entries if pinned.get(e["path"]) != e["sha256"]]
    if moved:
        return defect(f"pinned inputs moved since the registration: {moved}")
    ids = st.population_ids(pop)
    if st.ids_sha256(ids) != reg["population_sha256"] or ids != reg["population_ids"]:
        return defect("the population's sha256 differs from the registration's")

    crispri.annotate(training + pop, crispri.DeletionTable(), crispri.ElementCache())
    fit_rows = [p for p in training if p.covered]
    rates, borrowed, prov, lent = load_rates()
    clusters = regs.load_clusters(regs.CLUSTERS, chroms)
    readers = ce.Readers(RESULTS_DIR)
    by_chrom: dict[str, list[crispri.Pair]] = defaultdict(list)
    for p in training + pop:
        by_chrom[p.chrom].append(p)
    lookups: dict[str, st.ChromLookup] = {}
    for chrom in chroms:
        text = (Path(COMPILED_DIR) / f"noncoding_{chrom}.bio").read_text()
        lookups[chrom] = st.chromosome_lookup(
            chrom, text, by_chrom[chrom], rates, borrowed, prov, lent, clusters, readers, RESULTS_DIR
        )
        print(chrom, dict(lookups[chrom].counts), f"{time.time() - t0:.0f} s", flush=True)
    whys: dict[str, dict[str, str]] = {}
    for p in training + pop:
        f, why = st.trace_pair(p, lookups.get(p.chrom))
        p.features.update(f)
        whys[st.pair_id(p)] = why
    problems = st.harness(lookups, training + pop)
    if problems:
        return defect(f"{len(problems)} harness problems: " + "; ".join(problems[:20]))
    anchor = st.anchor_reading(fit_rows, pop)
    if not anchor["reproduced"]:
        return defect(
            f"the anchor did not reproduce: recomputed {anchor['recomputed']}, committed {st.ANCHOR}"
        )

    scored = st.score(fit_rows, pop, anchor)
    pop_whys = [whys[st.pair_id(p)] for p in pop]
    for stage, row in scored["stages"].items():
        row["population_dropped_because"] = st.dropped_because(pop_whys, stage)
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
        "outcome": scored["outcome"],
        "reading": scored["reading"],
        "stage_where_it_stops": scored["stage_where_it_stops"],
        "readings_verbatim": {s: r["reading"] for s, r in scored["stages"].items()},
        "stages": scored["stages"],
        "chain": scored["chain"],
        "branches": scored["branches"],
        "descriptives_not_verdicts": descriptives(lookups, pop),
        "limitations": list(st.LIMITATIONS),
        "pairs_table": pairs_table(pop),
        "seconds": round(time.time() - t0, 1),
    }
    print("saved", save_result(RESULT, out, manifest=manifest(entries)))
    print(json.dumps({"outcome": out["outcome"], "reading": out["reading"]}, indent=1))
    for s, r in out["readings_verbatim"].items():
        print(f"{s}: {r}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
