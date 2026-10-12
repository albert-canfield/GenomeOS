#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
"""execscale demo: one compiled chr21 programme RUN under the calibrated magnitude.

    uv run python scripts/execscale_demo.py        # write data/results/execscale_demo.json

WHAT THIS IS. The census (`execscale_census.py`, record `81d3949`) wrote down the committed state
before the calibrated magnitude existed; the parameter is now in the runtime, opt-in, default 1.0.
This is the one demonstration that the decoded genome RUNS at that magnitude: chr21 is compiled twice
by `attribution.compile.compile_chromosome`, once with `magnitude="compiled"` (what the compiler has
always emitted) and once with `magnitude="calibrated"` (the same programme plus one
`param magnitude_correction` block), and the calibrated programme is parsed and integrated by
`runtime.grn.NetworkRuntime` under the multiplicative composition, which reads the magnitude off the
programme itself.

WHAT IT CLAIMS: that it runs, and what it emits. NOTHING ABOUT RNA. No measured expression is read
here, nothing is fitted, and no comparison to a measurement is made or implied. The RNA evidence for
the magnitude is M2h's (`ac5580b`) and M2i's (`f6e71c6`), already registered, and is not re-read.

WHAT THE RUN HAD TO DECLARE, because the compiled programme does not carry it. A compiled chromosome
is RULES WITHOUT KINETICS: its gene stubs declare no `max_rate` and no `basal_rate`, and its rule
sources are elements, which are not species the runtime integrates. So this demo declares, as its own
inputs and as measurements of nothing: `max_rate` 1.0 and `basal_rate` 0.0 a.u. per hour on every
gene, every rule source held at 1.0 a.u. by `clamp`, `multiplicative_gain` left at the runtime's
default 1.0, and the kinetic parameters left at `grn.DEFAULTS`, which are arbitrary units. Every level
below is therefore in arbitrary units and is a statement about the COMPOSITION, not about expression.
The two runs differ in the magnitude and in nothing else: same programme rules, same declarations,
same clamp, same step size.

Input paths go to `input_entry` REPO-RELATIVE, so the record does not pin this machine's checkout.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import statistics
import sys
import time
import warnings
from pathlib import Path
from typing import Any

from genomeos import manifest as mf
from genomeos.attribution import compile as compile_mod
from genomeos.lang import parse
from genomeos.results import save_result
from genomeos.runtime import grn

RESULT = "execscale_demo"
OWN_CODE = (
    "scripts/execscale_demo.py",
    "scripts/execscale_census.py",
    "tests/test_execscale.py",
    "genomeos/runtime/grn.py",
    "genomeos/attribution/compile.py",
)
CHROM = "chr21"
#: what the demo itself declares, because the compiled programme carries no kinetics. Arbitrary units,
#: measurements of nothing, and identical between the two runs so that the magnitude is the only
#: difference between them.
DECLARED = {
    "max_rate": 1.0,
    "basal_rate": 0.0,
    "source_level": 1.0,
    "hours": 1.0,
    "dt": 0.05,
    "multiplicative_gain": grn.DEFAULT_MULTIPLICATIVE_GAIN,
    "units": "arbitrary concentration units and hours (grn.DEFAULTS); nothing here is measured",
}
#: The cell this run is in. A compiled rule names a cell in a `when: cell_type = ...` clause, so a run
#: with no context has NO active rule at all -- the first attempt at this demo composed nothing for
#: exactly that reason and is why the active-rule count is reported below. K562 is the cell the most
#: chr21 rules name (372 of 5,176), and it is a DECLARED CHOICE of this demo and not a finding: what a
#: compiled cell label is worth is `compile.py`'s own caveat and `data/results/argmaxcell.json`'s
#: measurement (the compiled cell is the argmax over the scorer's RNA track axis, the measured cell on
#: 612 of 3,228 K562 elements), and nothing here rests on it. Any other cell would run that cell's
#: rules the same way.
CONTEXT = {"cell_type": "K562"}
CLAIM = (
    "one compiled chr21 programme, emitted with the calibrated magnitude by"
    " compile_chromosome(magnitude='calibrated'), parses and RUNS under the multiplicative composition"
    " with the runtime reading the magnitude off the programme. What is reported is that it ran and"
    " what it emitted: the composed drive and the mRNA levels it produced, in arbitrary units, beside"
    " the same programme run at the compiled magnitude of 1.0"
)
DOES_NOT_CLAIM = (
    "this demo says NOTHING about RNA. No measured expression is read, nothing is fitted, no"
    " prediction is scored and no comparison to a measurement is made. Whether the calibrated"
    " composition is closer to measured RNA is M2h's result (`ac5580b`) and M2i's (`f6e71c6`) and is"
    " not re-read here. The magnitude is also a global AVERAGE over a magnitude that varies by"
    " chromosome (hetblock2, `d9e3ddd`), so running it on chr21 applies an average and not a"
    " chr21 number. The levels are arbitrary units from declared rates and a declared clamp, and are"
    " not expression"
)


def compile_both(chrom: str = CHROM) -> dict[str, Any]:
    """Compile the chromosome at both magnitudes and prove the option changed only the param block."""
    plain = compile_mod.compile_chromosome(chrom, magnitude=compile_mod.MAGNITUDE_COMPILED)
    calibrated = compile_mod.compile_chromosome(chrom, magnitude=compile_mod.MAGNITUDE_CALIBRATED)
    added = compile_mod.magnitude_lines(compile_mod.MAGNITUDE_CALIBRATED)
    a, b = plain.splitlines(), calibrated.splitlines()
    # The calibrated programme must be the compiled one with the block INSERTED and nothing else
    # touched: b == a[:i] + added + a[i:] for exactly one i. Checked line by line rather than by
    # string surgery, so a rule line that moved cannot hide inside a substring replacement.
    where = [i for i in range(len(a) + 1) if a[:i] + added + a[i:] == b]
    if len(where) != 1:
        raise SystemExit(
            "the calibrated emission differs from the compiled one by more than the magnitude block; "
            "the demo refuses to run a programme whose rules may have moved"
        )
    only_in_calibrated = b[where[0] : where[0] + len(added)]
    return {
        "chrom": chrom,
        "compiled_sha256": hashlib.sha256(plain.encode()).hexdigest(),
        "calibrated_sha256": hashlib.sha256(calibrated.encode()).hexdigest(),
        "compiled_lines": len(a),
        "calibrated_lines": len(b),
        "lines_added": len(b) - len(a),
        "magnitude_block_lines": len(added),
        "inserted_at_line": where[0] + 1,
        "added_lines_are_the_magnitude_block": only_in_calibrated == added,
        "declared_tests_identical": [ln for ln in a if ln.startswith("# test:")]
        == [ln for ln in b if ln.startswith("# test:")],
        "declared_tests": [ln for ln in b if ln.startswith("# test:")],
        "text": calibrated,
        "plain_text": plain,
    }


def _prepare(text: str) -> Any:
    """Parse the programme and put the demo's own declared kinetics on it."""
    module = parse(text)
    for gene in module.genes():
        gene.attrs["max_rate"] = DECLARED["max_rate"]
        gene.attrs["basal_rate"] = DECLARED["basal_rate"]
    return module


def one_run(text: str, correction: float | None) -> dict[str, Any]:
    """Integrate the programme once. `correction` None reads the programme's own declaration."""
    module = _prepare(text)
    t0 = time.time()
    rt = grn.NetworkRuntime(
        module,
        context=CONTEXT,
        activator_combination="multiplicative",
        magnitude_correction=correction,
    )
    if not rt.active_rules:
        # The first attempt at this demo ran with no context and composed nothing: every compiled rule
        # carries `when: cell_type = ...`, so a contextless run has no activator and A = 1 for every
        # gene whatever the magnitude is. A run that composes nothing must not be filed as a run.
        raise SystemExit(f"no rule of the programme is active in {CONTEXT}; nothing was composed")
    sources = sorted({r.source for r in rt.active_rules})
    clamp = dict.fromkeys(sources, DECLARED["source_level"])
    with warnings.catch_warnings():
        # the runtime WARNS about every regulator it holds no state for; they are all clamped above,
        # and a warning is not an error. The count of unresolved entries is reported instead.
        warnings.simplefilter("ignore", grn.UnresolvedModelWarning)
        traj = rt.run(DECLARED["hours"], dt=DECLARED["dt"], clamp=clamp, record_every=10_000)
        unresolved = rt.diagnose(clamp)
    acts = {g.id: len(rt.activators[g.id]) for g in rt.genes if rt.activators[g.id]}
    final = traj.final()
    mrna = {k: v for k, v in final.items() if k.endswith(".mRNA")}
    driven = {k: v for k, v in mrna.items() if v > 0.0}
    return {
        "magnitude_correction": rt.magnitude_correction,
        "read_from": "the programme's own param" if correction is None else "the caller's argument",
        "seconds": round(time.time() - t0, 2),
        "context": dict(CONTEXT),
        "genes": len(rt.genes),
        "active_rules": len(rt.active_rules),
        "genes_with_at_least_one_activator": len(acts),
        "most_activators_on_one_gene": max(acts.values()) if acts else 0,
        "clamped_sources": len(clamp),
        "unresolved_entries": len(unresolved),
        "species": len(traj.species),
        "mrna_species": len(mrna),
        "genes_with_a_level_above_zero": len(driven),
        "mrna_levels": {
            "min": min(mrna.values()) if mrna else None,
            "median": statistics.median(mrna.values()) if mrna else None,
            "max": max(mrna.values()) if mrna else None,
            "sum": sum(mrna.values()),
        },
        "highest_genes": sorted(mrna.items(), key=lambda kv: -kv[1])[:10],
        "_mrna": mrna,
    }


def compare(plain: dict[str, Any], calibrated: dict[str, Any]) -> dict[str, Any]:
    """How the two runs differ, as the composition's own numbers. No measurement enters this."""
    a, b = plain["_mrna"], calibrated["_mrna"]
    if set(a) != set(b):
        raise SystemExit("the two runs do not hold the same species; the comparison is refused")
    moved = [k for k in a if a[k] != b[k]]
    ratios = [b[k] / a[k] for k in a if a[k] > 0]
    moved_ratios = [b[k] / a[k] for k in moved if a[k] > 0]
    return {
        "the_direction_is_a_property_of_the_form_and_not_a_finding": (
            "under the multiplicative rule A = exp(c * gain * W) with every term t_i >= 0, so W >= 0"
            " and A >= 1: a correction below 1 can only move the drive TOWARD 1. That every moved"
            " species moved DOWN is therefore arithmetic, not evidence, and says nothing about"
            " whether the lower value is closer to anything measured"
        ),
        "moved_level_ratio_calibrated_over_compiled": {
            "n": len(moved_ratios),
            "min": min(moved_ratios) if moved_ratios else None,
            "median": statistics.median(moved_ratios) if moved_ratios else None,
            "max": max(moved_ratios) if moved_ratios else None,
        },
        "species_unmoved_have_no_activator_in_this_cell": (
            "a gene with no active activator rule composes A = exp(0) = 1 under every magnitude, so"
            " its level is identical in the two runs; that is what the unmoved count is"
        ),
        "species_compared": len(a),
        "species_whose_level_moved": len(moved),
        "species_unmoved": len(a) - len(moved),
        "level_ratio_calibrated_over_compiled": {
            "n": len(ratios),
            "min": min(ratios) if ratios else None,
            "median": statistics.median(ratios) if ratios else None,
            "max": max(ratios) if ratios else None,
        },
        "log_level_ratio_median": (statistics.median([math.log(r) for r in ratios]) if ratios else None),
        "calibrated_is_lower_for": sum(1 for k in a if b[k] < a[k]),
        "calibrated_is_higher_for": sum(1 for k in a if b[k] > a[k]),
    }


def manifest(programmes: dict[str, Any], read_now: list[str]) -> dict:
    """The record's manifest, with the inputs taken from what the audit hook actually recorded.

    The compiler reads some forty files for one chromosome -- the budget, the domains, the unknown
    blocks, the reader tracks, the cCRE and repeat tables, GENCODE -- and a hand-written list of four
    was wrong on the first attempt, which the registry refused. So the declaration is built from
    `mf.traced_reads()`: every file this run opened under data/ is declared and hashed, grouped by
    tree, and nothing is declared that was not read.
    """
    groups = {
        "the compiler's committed per-chromosome results for this chromosome": [
            p for p in read_now if p.startswith("data/results/")
        ],
        "the local element tables and knowledge files the compiler's readers open": [
            p for p in read_now if p.startswith("data/knowledge/")
        ],
        "the reference annotation read while compiling": [
            p for p in read_now if p.startswith("data/reference/")
        ],
        "cached tables read while compiling": [p for p in read_now if p.startswith("data/cache/")],
    }
    inputs = [
        mf.files_entry(label, sorted(paths), partition=CHROM) for label, paths in groups.items() if paths
    ]
    grouped = {p for paths in groups.values() for p in paths}
    inputs += [
        mf.input_entry(p, partition=CHROM)
        for p in sorted(set(read_now) - grouped)
        if p.startswith("data/") and (Path(p).exists())
    ]
    return {
        "sources": [
            {"accession": "the compiler's own chr21 inputs in data/results", "version": "as committed"},
            {"accession": "M2h's fitted scale, ac5580b", "version": "0.2146, read and not refitted"},
        ],
        "inputs": inputs,
        "input_count": len(inputs),
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "chrom": CHROM,
            "emission": compile_mod.EMISSION_DEFAULT,
            "magnitudes": [compile_mod.MAGNITUDE_COMPILED, compile_mod.MAGNITUDE_CALIBRATED],
            "activator_combination": "multiplicative",
            "context": dict(CONTEXT),
            "declared_by_the_demo": DECLARED,
            "calibrated_programme_sha256": programmes["calibrated_sha256"],
        },
        "exclusions": [
            "no model request, no key, no network, no money; no RNA is read, nothing is fitted,"
            " and no committed programme is written"
        ],
        "code_cleanliness": mf.code_cleanliness(__file__, OWN_CODE),
        "partitions": {"chromosome": CHROM},
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--print", action="store_true", help="print the record instead of saving it")
    args = ap.parse_args(argv)
    mf.trace_begin()  # the manifest's inputs are what this run opened, not what anyone remembered
    programmes = compile_both()
    plain = one_run(programmes["plain_text"], grn.DEFAULT_MAGNITUDE_CORRECTION)
    calibrated = one_run(programmes["text"], None)
    if calibrated["magnitude_correction"] != grn.MEASURED_MAGNITUDE_CORRECTION:
        raise SystemExit("the calibrated run did not read the programme's magnitude; nothing is saved")
    diff = compare(plain, calibrated)
    for run in (plain, calibrated):
        run.pop("_mrna")
    payload = {
        "question": (
            "does one compiled chr21 programme RUN under the calibrated magnitude, and what does it emit?"
        ),
        "claim": CLAIM,
        "does_not_claim": DOES_NOT_CLAIM,
        "magnitude": {
            "value": grn.MEASURED_MAGNITUDE_CORRECTION,
            "source_commit": "ac5580b",
            "it_is_a_global_average": (
                "+0.2146 is a global AVERAGE over a magnitude that is not one global number: hetblock2"
                " (`d9e3ddd`) measured it varying by chromosome, surviving a locus-block null"
            ),
            "default_in_the_runtime": grn.DEFAULT_MAGNITUDE_CORRECTION,
        },
        "programmes": {k: v for k, v in programmes.items() if k not in ("text", "plain_text")},
        "run_at_the_compiled_magnitude": plain,
        "run_at_the_calibrated_magnitude": calibrated,
        "difference": diff,
        "declared_by_the_demo": DECLARED,
    }
    if args.print:
        print(json.dumps(payload, indent=1, default=str))
        return 0
    print("saved", save_result(RESULT, payload, manifest=manifest(programmes, mf.traced_reads())))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
