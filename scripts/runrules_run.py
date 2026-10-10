# SPDX-License-Identifier: AGPL-3.0-or-later
"""M1, the registered run: the compiled rules, run, against measured RNA in 54 GTEx tissues.

    uv run python scripts/runrules_run.py

Runs only from the committed registration (data/results/runrules_registration.json) and refuses if
any registered field differs from the constants in genomeos/attribution/runrules.py that apply it.

Order, as registered (runrules.DECISION):
  1. every (gene, panel tissue) pair the bridge resolves is bridged on its own and RUN by
     `runtime.grn.NetworkRuntime` with its element clamped present and absent;
  2. the run is checked against its closed form at every pair; a pair off by more than the tolerance
     stops the run here, before the measured RNA is opened;
  3. only then is the GTEx table opened, the population assembled, and the comparison scored.

No model request, no network, no money: every input is on this machine and pinned by sha256.
"""

from __future__ import annotations

import json
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import bridge  # noqa: E402
from genomeos.attribution import runrules as rr  # noqa: E402
from genomeos.evidence import COMPILED_DIR  # noqa: E402
from genomeos.lang.parser import parse  # noqa: E402
from genomeos.results import RESULTS_DIR, save_result  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from runrules_census import OWN_CODE, RATE_TABLE, load_rates  # noqa: E402

RESULT = "runrules_result"
REGISTRATION = RESULTS_DIR / "runrules_registration.json"
CHROMS = [f"chr{c}" for c in [*range(1, 23), "X", "Y"]]
#: the registration fields that must equal the module's constants, or the run refuses
BOUND = {
    "question": "QUESTION",
    "reading": "READING",
    "population": "POPULATION",
    "measured": "MEASURED",
    "baseline": "BASELINE",
    "metric": "METRIC",
    "control": "CONTROL",
    "decision": "DECISION",
    "falsifier": "FALSIFIER",
    "stop_rule": "STOP_RULE",
    "outcomes": "OUTCOMES",
    "limitations": "LIMITATIONS",
    "min_reference_tissues": "MIN_REFERENCE_TISSUES",
    "min_population": "MIN_POPULATION",
    "pseudocount": "PSEUDOCOUNT",
    "bootstrap": "BOOTSTRAP",
    "permutations": "PERMUTATIONS",
    "control_alpha": "CONTROL_ALPHA",
    "runtime": "RUNTIME",
    "closed_form_tolerance": "CLOSED_FORM_TOLERANCE",
    "gtex_tissues": "GTEX_V8_TISSUES",
}


def check_registration() -> dict:
    """The committed registration, refused if any bound field differs from the code that applies it."""
    reg = json.loads(REGISTRATION.read_text())
    bad = []
    for field, const in BOUND.items():
        want = json.loads(json.dumps(getattr(rr, const)))
        if reg.get(field) != want:
            bad.append(field)
    if bad:
        raise SystemExit(f"refused: the registration and runrules.py disagree on {', '.join(bad)}")
    return reg


def run_all(chroms: list[str]) -> tuple[list[dict], dict[str, set[str]], Counter, list[Path]]:
    """Every resolved panel pair, run; each gene's rule tissues; the counts on the way."""
    rates, borrowed, prov, lent = load_rates()
    gene_tissues: dict[str, set[str]] = defaultdict(set)
    mechanisms: Counter = Counter()
    found: dict[tuple[str, str], list[dict]] = defaultdict(list)
    counts: Counter = Counter()
    paths = []
    for chrom in chroms:
        path = Path(COMPILED_DIR) / f"noncoding_{chrom}.bio"
        paths.append(path)
        m = parse(path.read_text(), path.stem)
        for g, cells in rr.rule_tissues(m).items():
            gene_tissues[g] |= cells
        present = sorted({r.when.get("cell_type") for r in m.rules} & set(rr.PANEL))
        for cell in present:
            ctx = {"cell_type": cell}
            for g, k in rr.mechanisms_in(m, ctx).items():
                mechanisms[(g, cell)] += k
            b = bridge.parameterize(
                m,
                ctx,
                build=False,
                rates=rates,
                borrowed_rate=borrowed,
                rate_provenance=prov,
                borrowed_provenance=lent,
            )
            for i in b.unresolved:
                counts[f"panel_unresolved_{i.reason}"] += 1
            for key, s in b.strengths.items():
                gene = key[1]
                pair = rr.run_pair(rr.one_gene_module(m, gene), gene, ctx, rates, borrowed, prov, lent)
                counts["pairs_run"] += 1
                if pair is None or abs(pair["strength"] - s) > bridge.STRENGTH_TOLERANCE:
                    counts["small_module_disagrees_with_whole_program"] += 1
                    continue
                pair["chrom"] = chrom
                found[(gene, cell)].append(pair)
        print(chrom, counts["pairs_run"], flush=True)
    pairs = []
    for (gene, cell), rows in sorted(found.items()):
        if mechanisms[(gene, cell)] != 1 or len(rows) != 1:
            counts["excluded_several_mechanisms_over_programs"] += len(rows)
            continue
        pairs.append(rows[0])
    counts["resolved_pairs_single_over_programs"] = len(pairs)
    return pairs, gene_tissues, counts, paths


def manifest(paths: list[Path], gtex_read: bool) -> dict:
    inputs = (
        [mf.input_entry(p, partition=None) for p in paths]
        + [mf.input_entry(RATE_TABLE, partition=None), mf.input_entry(REGISTRATION, partition=None)]
        + (
            [mf.input_entry(rr.GTEX, partition="measured outcome, opened after the closed-form check")]
            if gtex_read
            else []
        )
    )
    return {
        "sources": [
            {"accession": "the 24 compiled non-coding programs", "version": "pinned by sha256 in inputs"},
            {
                "accession": "data/results/gene_rates.json (mouse NIH 3T3, e921fae)",
                "version": "pinned by sha256",
            },
            {
                "accession": "GTEx v8 gene median TPM, GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9",
                "version": "the cached file, pinned by sha256 in inputs and in the registration",
            },
        ],
        "inputs": inputs,
        "input_count": len(inputs),
        "assembly": "GRCh38 for the compiled programs; GTEx v8 rows are matched by gene symbol only",
        "coordinates": "n/a: genes are matched by symbol, tissues by name (runrules.tissue_label)",
        "parameters": {
            "runtime": rr.RUNTIME,
            "closed_form_tolerance": rr.CLOSED_FORM_TOLERANCE,
            "min_reference_tissues": rr.MIN_REFERENCE_TISSUES,
            "min_population": rr.MIN_POPULATION,
            "pseudocount": rr.PSEUDOCOUNT,
            "bootstrap": rr.BOOTSTRAP,
            "permutations": rr.PERMUTATIONS,
            "control_alpha": rr.CONTROL_ALPHA,
        },
        "exclusions": [
            "pairs with several mechanisms in their tissue, within a program or pooled over programs",
            "pairs the bridge does not resolve on the borrowed tier",
            "genes whose symbol is not on exactly one GTEx row",
            f"genes with fewer than {rr.MIN_REFERENCE_TISSUES} panel tissues free of any compiled rule",
        ],
        "partitions": "n/a: no training partition; every compiled panel pair that qualifies is scored",
        "code_cleanliness": mf.code_cleanliness(__file__, OWN_CODE),
    }


def main() -> int:
    t0 = time.time()
    reg = check_registration()
    pairs, gene_tissues, counts, paths = run_all(CHROMS)
    bad = [p for p in pairs if abs(p["m"] - p["closed_form"]) > rr.CLOSED_FORM_TOLERANCE]
    base = {
        "lane": rr.LANE,
        "registration": "data/results/runrules_registration.json",
        "registration_date": reg.get("date"),
        "question": rr.QUESTION,
        "structure": rr.STRUCTURE,
        "run_counts": dict(counts),
        "closed_form_check": {
            "pairs": len(pairs),
            "beyond_tolerance": len(bad),
            "max_abs_difference_log2": max((abs(p["m"] - p["closed_form"]) for p in pairs), default=None),
        },
        "alphagenome_requests": 0,
        "money": "none: no request is sent and no key is read",
    }
    if bad or counts["small_module_disagrees_with_whole_program"]:
        out = {
            **base,
            "outcome": "HARNESS_DEFECT",
            "reading": rr.reading(
                "HARNESS_DEFECT",
                tol=rr.CLOSED_FORM_TOLERANCE,
                bad=len(bad) + counts["small_module_disagrees_with_whole_program"],
                n=len(pairs),
            ),
            "seconds": round(time.time() - t0, 1),
        }
        print("saved", save_result(RESULT, out, manifest=manifest(paths, gtex_read=False)))
        print(out["reading"])
        return 1
    tissues, tpm, gtex_counts = rr.load_gtex(rr.GTEX)
    pop, why = rr.assemble(pairs, gene_tissues, tpm, tissues)
    scored = rr.score(pop)
    cols = (
        "gene",
        "symbol",
        "tissue",
        "element",
        "chrom",
        "action",
        "rho",
        "m",
        "D",
        "tier",
        "on_in_rule_tissue",
    )
    out = {
        **base,
        "gtex": gtex_counts,
        "population_exclusions": dict(why),
        **scored,
        "pairs_table": {"columns": list(cols), "rows": [[p[c] for c in cols] for p in pop]},
        "seconds": round(time.time() - t0, 1),
    }
    print("saved", save_result(RESULT, out, manifest=manifest(paths, gtex_read=True)))
    print(json.dumps({k: out[k] for k in ("outcome", "reading", "pairs", "genes", "tissues")}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
