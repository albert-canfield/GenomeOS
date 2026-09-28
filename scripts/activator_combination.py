# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""The activator-combination census and choice test (lane-combine, 2026-09-28).

`genomeos/runtime/grn.py` (and `runtime/located.py`, which repeats the formula) combine a gene's
activators by their mean. This script answers, by computing and not by reading code:

    uv run python scripts/activator_combination.py static     # every program: genes with >= 2 activators
    uv run python scripts/activator_combination.py gastrulation   # the one program it moves, per rule
    uv run python scripts/activator_combination.py collate --scratch DIR   # per-rule test and bio-test runs
    uv run python scripts/activator_combination.py choose     # the registered choice test, run once

`static` parses every .bio program under data/demo, data/organisms, genomeos/std and the 24 compiled
chromosomes in data/knowledge/compiled, and for every cell context a program's rules name counts the
genes a run there would give two or more activating rules, and whether any of them can be simulated
(declares `max_rate`). `collate` reads the scratch outputs of the per-rule runs made with
`scripts/activator_combination_patch.py` (pytest plugin and `bio test` wrapper) and writes which tests
and program self-tests pass or fail under each candidate rule. Writes
`data/results/activator_combination_census.json`. No network, no AlphaGenome.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from collections import Counter
from pathlib import Path

from genomeos.ir import Action
from genomeos.lang.parser import parse_file
from genomeos.results import RESULTS_DIR, save_result

ROOT = Path(__file__).resolve().parents[1]
NAME = "activator_combination_census"
PROGRAM_DIRS = ["data/demo", "data/organisms", "genomeos/std"]
COMPILED_DIR = ROOT / "data/knowledge/compiled"
RULES = ("mean", "sum_capped", "max", "or")


def contexts_of(m) -> list[dict[str, str]]:
    """The empty context and every distinct `when` a rule of the module states."""
    seen = {(): {}}
    for r in m.rules:
        if r.when:
            seen.setdefault(tuple(sorted(r.when.items())), dict(r.when))
    return list(seen.values())


def census_module(m) -> dict:
    """Genes with two or more active activators, per context the module's rules name."""
    genes = {g.id: g for g in m.genes()}
    pairs = []
    ctxs = contexts_of(m)
    for ctx in ctxs:
        try:
            active = m.active_rules(ctx)
        except Exception:  # noqa: BLE001 - a context the module refuses is reported, not guessed
            continue
        acts = Counter(r.target for r in active if r.action is Action.ACTIVATE and r.target in genes)
        for g, n in acts.items():
            if n >= 2:
                pairs.append(
                    {
                        "context": ctx,
                        "gene": g,
                        "activators": sorted(
                            r.source for r in active if r.action is Action.ACTIVATE and r.target == g
                        ),
                        "simulable": "max_rate" in genes[g].attrs,
                    }
                )
    return {"contexts": len(ctxs), "pairs": pairs}


def static() -> dict:
    programs = {}
    for d in PROGRAM_DIRS:
        for p in sorted((ROOT / d).rglob("*.bio")):
            rel = str(p.relative_to(ROOT))
            try:
                m = parse_file(p)
            except Exception as e:  # noqa: BLE001
                programs[rel] = {"error": f"{type(e).__name__}: {e}"[:200]}
                continue
            c = census_module(m)
            # a gene counted once however many contexts give it several activators
            genes = sorted({(x["gene"], x["simulable"]) for x in c["pairs"]})
            programs[rel] = {
                "contexts": c["contexts"],
                "genes_with_2plus_activators": [g for g, _ in genes],
                "simulable": [g for g, s in genes if s],
                "pairs": c["pairs"] if len(c["pairs"]) <= 40 else c["pairs"][:40],
                "pair_count": len(c["pairs"]),
            }
    compiled = {}
    for p in sorted(COMPILED_DIR.glob("*.bio")):
        t0 = time.time()
        m = parse_file(p)
        genes = {g.id: g for g in m.genes()}
        cells = sorted({r.when.get("cell_type") for r in m.rules if r.when.get("cell_type")})
        n_pairs = n_multi = n_multi_sim = rules_on_multi = 0
        for cell in cells:
            active = m.active_rules({"cell_type": cell})
            acts = Counter(r.target for r in active if r.action is Action.ACTIVATE and r.target in genes)
            n_pairs += len(acts)
            multi = [g for g, n in acts.items() if n >= 2]
            n_multi += len(multi)
            rules_on_multi += sum(acts[g] for g in multi)
            n_multi_sim += sum(1 for g in multi if "max_rate" in genes[g].attrs)
        compiled[p.name] = {
            "cell_contexts": len(cells),
            "gene_cell_pairs_with_an_activator": n_pairs,
            "gene_cell_pairs_with_2plus_activators": n_multi,
            "activating_rules_on_those": rules_on_multi,
            "of_which_simulable": n_multi_sim,
            "seconds": round(time.time() - t0, 1),
        }
        print(p.name, compiled[p.name], file=sys.stderr, flush=True)
    return {"programs": programs, "compiled": compiled}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("step", choices=["static", "gastrulation", "collate", "choose"])
    ap.add_argument("--scratch", help="directory holding the per-rule runs (collate)")
    ap.add_argument("--out", help="write here instead of data/results")
    args = ap.parse_args(argv)
    path = Path(args.out) if args.out else RESULTS_DIR / f"{NAME}.json"
    doc = json.loads(path.read_text()) if path.exists() else {}
    if args.step == "static":
        doc["static"] = static()
    elif args.step == "choose":
        doc["choice_test"] = choose()
    elif args.step == "gastrulation":
        doc["gastrulation"] = gastrulation()
    else:
        doc["dynamic"] = collate(Path(args.scratch))
    doc.pop("result_manifest", None)
    doc["question"] = (
        "which programs, tests and results depend on how a gene's activators combine, what each "
        "candidate rule does to them, and whether measured evidence picks a rule without tuning"
    )
    if args.out:
        path.write_text(json.dumps(doc, indent=1, sort_keys=True))
    else:
        save_result(NAME, doc, manifest=manifest("choice_test" in doc))


def manifest(with_crispri: bool) -> dict:
    from genomeos import manifest as mf

    inputs = [mf.input_entry(d) for d in PROGRAM_DIRS]
    inputs += [
        mf.input_entry(str(COMPILED_DIR.relative_to(ROOT))),
        mf.input_entry("genomeos/runtime/grn.py"),
        mf.input_entry("genomeos/runtime/located.py"),
        mf.input_entry(RESULTS_DIR / "gastrulation_census_comparison_sign_corrected.json"),
    ]
    sources = [
        {"accession": "in-repo .bio programs and compiled chromosomes", "version": "this revision"},
        {"accession": "Bothma et al. 2015 eLife 4:e07956, PMC4532966", "version": "published"},
        {"accession": "Zhou et al. 2024 Cell Genomics (bioRxiv 2023.04.26.538501)", "version": "published"},
    ]
    partitions = "n/a: no evaluation split read"
    if with_crispri:
        from genomeos.attribution import measured

        name = next(n for n, sp in measured.CRISPRI_SPLIT_OF.items() if sp == measured.TRAINING)
        inputs.append(mf.input_entry(measured.CRISPRI_KNOWLEDGE / name, partition=measured.TRAINING))
        sources.append({"accession": "EPCrisprBenchmark training_K562 GRCh38", "version": name})
        partitions = {"read": [measured.TRAINING], "heldout": "not read; development_only enforced"}
    return {
        "sources": sources,
        "inputs": inputs,
        "assembly": "GRCh38 for the CRISPRi table; otherwise n/a: simulated networks",
        "coordinates": "n/a: element intervals only key elements within one table, never compared",
        "parameters": {"rules": list(RULES), "gastrulation_fixture": "60 cells, 30 h, dt 0.05"},
        "exclusions": ["CRISPRi element-gene-cells with conflicting outcomes across rows"],
        "partitions": partitions,
    }
    return 0


def gastrulation() -> dict:
    """The gastrulation model under each rule: the census run, the test fixture, and the census verdict
    judged exactly as `scripts/gastrulation_census.py compare` judges it (intervals from the committed
    comparison, tolerance CENSUS_TOLERANCE). Descriptive: the census is not a choice criterion."""
    from genomeos.runtime import gastrulation as g
    from genomeos.runtime import grn

    committed = json.loads((RESULTS_DIR / "gastrulation_census_comparison_sign_corrected.json").read_text())
    intervals = {k: v["interval"] for k, v in committed["human_cs7_judgement"].items()}
    fixture = {"cells": 60, "hours": 30.0, "dt": 0.05}
    out = {"intervals_human_cs7": intervals, "tolerance": g.CENSUS_TOLERANCE, "per_rule": {}}
    saved = grn.ACTIVATOR_COMBINATION
    try:
        for rule in RULES:
            grn.ACTIVATOR_COMBINATION = rule
            census = {
                k: round(v, 4) for k, v in g.run_gastrulation(**g.CENSUS_MODEL_RUN).proportions().items()
            }
            fix = {k: round(v, 4) for k, v in g.run_gastrulation(**fixture).proportions().items()}
            judge = {}
            for layer, (lo, hi) in intervals.items():
                m = census.get(layer, 0.0)
                gap = round(max(lo - m, m - hi, 0.0), 4)
                judge[layer] = {"model": m, "outside_by": gap, "pass": gap <= g.CENSUS_TOLERANCE}
            verdict = "not contradicted" if all(v["pass"] for v in judge.values()) else "falsified"
            out["per_rule"][rule] = {
                "census_run": census,
                "test_fixture_60_cells_30h": fix,
                "census_judgement": judge,
                "census_verdict": verdict,
            }
            print(rule, census, fix, verdict, file=sys.stderr, flush=True)
    finally:
        grn.ACTIVATOR_COMBINATION = saved
    return out


def crispri_training() -> dict:
    """Single-removal outcomes per (gene, cell) in the CRISPRi training split, development only.

    An element-gene-cell measured in several rows with different outcomes is counted as conflicting
    and left out of the per-gene counts. Decreases and increases are the benchmark's significant
    calls split by the sign of EffectSize (measured.OUTCOMES)."""
    from genomeos.attribution import measured

    pairs, invalid = measured.load_crispri(split=measured.TRAINING)
    pairs = measured.development_only(pairs)
    by_elem: dict[tuple, set] = {}
    effect: dict[tuple, list[float]] = {}
    for p in pairs:
        k = (p.gene, p.cell, p.chrom, p.start, p.end)
        by_elem.setdefault(k, set()).add(p.outcome)
        effect.setdefault(k, []).append(p.effect_size)
    genes: dict[tuple, Counter] = {}
    mag: dict[tuple, Counter] = {}
    conflicting = 0
    for k, outs in by_elem.items():
        if len(outs) > 1:
            conflicting += 1
            continue
        (o,) = outs
        gc = k[:2]
        genes.setdefault(gc, Counter())[o] += 1
        if o in (measured.DECREASE, measured.INCREASE):
            mag.setdefault(gc, Counter())[o] += abs(sum(effect[k]) / len(effect[k]))
    dec = measured.DECREASE
    inc = measured.INCREASE
    two = [gc for gc, c in genes.items() if c[dec] >= 2]
    dist = Counter(min(c[dec], 5) for c in genes.values())
    d_sum = sum(mag[gc][dec] for gc in two)
    i_sum = sum(mag.get(gc, Counter())[inc] for gc in two)
    return {
        "rows": len(pairs),
        "rows_invalid": invalid,
        "element_gene_cells": len(by_elem),
        "element_gene_cells_conflicting": conflicting,
        "gene_cells": len(genes),
        "gene_cells_by_significant_decreases_capped_at_5": {str(k): v for k, v in sorted(dist.items())},
        "gene_cells_with_2plus_decreases": len(two),
        "of_those_with_an_increase": sum(1 for gc in two if genes[gc][inc] >= 1),
        "increase_over_decrease_magnitude_on_those": round(i_sum / d_sum, 4) if d_sum else None,
        "gene_cells_with_an_increase": sum(1 for c in genes.values() if c[inc] >= 1),
        "examples_2plus_decreases": sorted(
            (
                {"gene": g, "cell": c, "decreases": genes[(g, c)][dec], "increases": genes[(g, c)][inc]}
                for g, c in two
            ),
            key=lambda x: (-x["decreases"], x["gene"]),
        )[:15],
    }


def choose() -> dict:
    """The registered choice test (grn.COMBINATION_TEST), run once."""
    from genomeos.runtime import grn

    t = grn.COMBINATION_TEST
    crispri = crispri_training()
    contradicted: dict[str, list[str]] = {r: [] for r in RULES}
    excluded = []
    for name, case in t["cases"].items():
        c = case["contradicts"]
        if name == "crispri_training_two_decreases":
            c = (
                ("max",)
                if crispri["gene_cells_with_2plus_decreases"] >= t["crispri_min_gene_cells_for_max"]
                else ()
            )
        if t["exclude_cases_contradicting_every_rule"] and set(c) >= set(RULES):
            excluded.append(name)
            continue
        for r in c:
            contradicted[r].append(name)
    survivors = [r for r in RULES if r != "mean" and not contradicted[r]]
    if contradicted["mean"] and len(survivors) == 1:
        decision = f"switch to {survivors[0]}"
    else:
        decision = "mean retained, undetermined"
    return {
        "registration": "docs/DESIGN-MINIMAL-CELL.md, 'The activator combination rule: a choice test "
        "(registered 2026-09-28)'; grn.COMBINATION_TEST",
        "crispri_training": crispri,
        "excluded_as_contradicting_every_rule": excluded,
        "contradicted_by": contradicted,
        "surviving_alternatives": survivors,
        "decision": decision,
        "not_judged": t["not_judged"],
    }


def collate(scratch: Path) -> dict:
    """Per rule: the full pytest outcome, the tests that evaluated a gene with >= 2 activators, and the
    `bio test` totals, from the scratch runs made with scripts/activator_combination_patch.py."""
    import re

    out = {}
    for rule in RULES:
        log = (scratch / f"pytest_{rule}.log").read_text()
        summary = [ln for ln in log.splitlines() if re.search(r"\d+ passed", ln)][-1].strip("= ")
        failed = sorted({m.group(1) for m in re.finditer(r"^FAILED (\S+)", log, re.M)})
        xpassed = sorted({m.group(1) for m in re.finditer(r"^XPASS (\S+)", log, re.M)})
        sens = json.loads((scratch / f"sens_{rule}.json").read_text())
        sens = sens.get("sensitive", sens)
        bio = (scratch / f"bio_{rule}.log").read_text().strip().splitlines()[-1]
        out[rule] = {
            "pytest_summary": summary,
            "failed": failed,
            "xpassed": xpassed,
            "tests_reaching_a_gene_with_2plus_activators": sorted(k for k in sens if "::" in k),
            "bio_test": bio,
        }
    return out


if __name__ == "__main__":
    raise SystemExit(main())
