# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""Retrospective benchmark: does the therapeutic pipeline recover known targets?

Six tumours whose driver and whose approved therapy are public knowledge are
run through `genomeos therapeutic`, and the result is compared with what is
actually approved for that alteration. The point is falsifiability: until a
pipeline is asked to reproduce something already known, its rankings are
assertions.

Two different questions are scored separately, because conflating them would
flatter the pipeline:

  1. **Is the target recovered?** Does the driver's gene appear among the
     ranked candidates at all? A pipeline that loses ERBB2 in an
     ERBB2-amplified tumour is broken regardless of what it says about
     mechanism.

  2. **Is the route called correctly?** GenomeOS models antibody-like
     modalities: blocking antibodies, ADCC, engagers, conjugates,
     radionuclides, TCR and TCR-mimic. It models **no small molecules**.
     For BRAF, KRAS, PIK3CA and IDH1 the approved drug is a small-molecule
     inhibitor of an intracellular protein, so the correct answer is "no
     surface route, and this is outside what GenomeOS models". Counting those
     as failures would be dishonest; counting them as successes without
     saying why would be worse. They are scored as `out_of_scope_expected`,
     and the benchmark fails if the pipeline claims a surface route for them.

Run: `uv run python scripts/therapeutic_benchmark.py`
Writes: data/results/therapeutic_benchmark.json
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

from genomeos.results import save_result
from genomeos.therapeutics import analyse_vcf

BENCH = Path("data/demo/benchmark")

#: Each case: the alteration, the gene an oncologist would name, and what is
#: actually approved for it. `modality` is the class of the approved drug.
CASES: tuple[dict, ...] = (
    {
        "case": "EGFR L858R, lung adenocarcinoma",
        "vcf": "egfr_l858r_lung.vcf",
        "gene": "EGFR",
        "approved": (
            "osimertinib (small molecule); cetuximab and panitumumab are approved "
            "antibodies against EGFR in colorectal carcinoma"
        ),
        "modality": "antibody",
        "expect": "surface",
        "why": (
            "a single-pass receptor with a large extracellular domain and an approved "
            "antibody, so a surface route must be found"
        ),
    },
    {
        "case": "ERBB2-amplified breast carcinoma",
        "vcf": "erbb2_amplified_breast.vcf",
        "cnv": "erbb2_amplified_breast.cnv",
        "gene": "ERBB2",
        "approved": "trastuzumab, pertuzumab (antibodies); trastuzumab deruxtecan (antibody-drug conjugate)",
        "modality": "antibody",
        "expect": "surface",
        "why": "the canonical surface target in oncology; if any case must work, it is this one",
    },
    {
        "case": "BRAF V600E, melanoma",
        "vcf": "braf_v600e_melanoma.vcf",
        "gene": "BRAF",
        "approved": "vemurafenib, dabrafenib (small molecules)",
        "modality": "small_molecule",
        "expect": "out_of_scope_expected",
        "why": "a cytoplasmic kinase; the approved drug is a small molecule, a class GenomeOS does not model",
    },
    {
        "case": "KRAS G12C, lung adenocarcinoma",
        "vcf": "kras_g12c_lung.vcf",
        "gene": "KRAS",
        "approved": "sotorasib, adagrasib (covalent small molecules)",
        "modality": "small_molecule",
        "expect": "out_of_scope_expected",
        "why": (
            "membrane-anchored on the inner leaflet, so no outward-facing epitope; the "
            "mutant peptide route is the one GenomeOS can reason about"
        ),
    },
    {
        "case": "PIK3CA H1047R, breast carcinoma",
        "vcf": "pik3ca_h1047r_breast.vcf",
        "gene": "PIK3CA",
        "approved": "alpelisib (small molecule)",
        "modality": "small_molecule",
        "expect": "out_of_scope_expected",
        "why": "intracellular lipid kinase",
    },
    {
        "case": "IDH1 R132H, glioma and AML",
        "vcf": "idh1_r132h_glioma.vcf",
        "gene": "IDH1",
        "approved": "ivosidenib (small molecule)",
        "modality": "small_molecule",
        "expect": "out_of_scope_expected",
        "why": (
            "cytosolic enzyme; a neomorphic mutation whose product is a metabolite, "
            "which GenomeOS does not model at all"
        ),
    },
)


def run_case(case: dict, net: bool, log) -> dict:
    vcf = BENCH / case["vcf"]
    cnv = str(BENCH / case["cnv"]) if case.get("cnv") else None
    t0 = time.time()
    a = analyse_vcf(str(vcf), cnv=cnv, top_genes=8, net=net, indirect=True, log=None)
    gene = case["gene"]
    ranked = [c.gene for c in a["candidates"]]
    hit = next((c for c in a["candidates"] if c.gene == gene), None)
    out = {
        "case": case["case"],
        "gene": gene,
        "approved": case["approved"],
        "approved_modality": case["modality"],
        "expected": case["expect"],
        "why": case["why"],
        "candidates": ranked,
        "recovered": hit is not None,
        "rank": (ranked.index(gene) + 1) if hit is not None else None,
        "seconds": round(time.time() - t0, 1),
    }
    if hit is not None:
        best = hit.best_mechanism
        out.update(
            {
                "target_class": hit.target_class,
                "score": None if hit.scores.overall is None else round(hit.scores.overall, 3),
                "surface_accessibility": hit.scores.value("surface_accessibility"),
                "best_mechanism": best.mechanism if best else None,
                "best_compatibility": round(best.compatibility, 3) if best else None,
                "peptide_route": bool(hit.neoantigen and hit.neoantigen.peptides),
            }
        )
        # A third question, scored apart from the other two: is the mechanism the
        # pipeline puts first defensible? Two ways it is not. An agonist antibody
        # against an activating oncogenic driver would push the pathway the tumour
        # already over-drives. And a mechanism that needs an extracellular epitope
        # should not head the list for a protein with no established outward-facing
        # part, however low its compatibility.
        mech = out.get("best_mechanism") or ""
        acc = hit.scores.value("surface_accessibility") or 0
        problems = []
        if "agonist" in mech:
            problems.append(
                f"top mechanism is {mech} against an activating driver, which would drive the pathway"
            )
        if acc <= 0.5 and mech in ("targeted_radionuclide", "adc", "blocking_antibody", "adcc", "adcp"):
            problems.append(f"top mechanism {mech} needs an extracellular epitope, accessibility is {acc}")
        out["mechanism_concerns"] = problems
        out["mechanism_sane"] = not problems

        surface = acc > 0.5
        if case["expect"] == "surface":
            out["verdict"] = "recovered as a surface target" if surface else "RECOVERED BUT NO SURFACE ROUTE"
            out["pass"] = surface
        else:
            out["verdict"] = (
                "correctly not a surface target; approved drug is out of scope"
                if not surface
                else "CLAIMS A SURFACE ROUTE THAT DOES NOT EXIST"
            )
            out["pass"] = not surface
    else:
        out["verdict"] = "TARGET LOST"
        out["pass"] = False
    print(f"  {gene}: {out['verdict']} (rank {out['rank']}, {out['seconds']}s)", file=log, flush=True)
    return out


def main() -> int:
    net = "--offline" not in sys.argv
    log = sys.stdout
    print(f"therapeutic benchmark: {len(CASES)} cases, network={'on' if net else 'off'}", file=log)
    rows = [run_case(c, net, log) for c in CASES]
    recovered = sum(1 for r in rows if r["recovered"])
    passed = sum(1 for r in rows if r["pass"])
    sane = sum(1 for r in rows if r.get("mechanism_sane"))
    surface_cases = [r for r in rows if r["expected"] == "surface"]
    result = {
        "result": "therapeutic_benchmark",
        "cases": len(rows),
        "targets_recovered": recovered,
        "verdicts_correct": passed,
        "top_mechanism_defensible": sane,
        "surface_cases": len(surface_cases),
        "surface_cases_recovered": sum(1 for r in surface_cases if r["pass"]),
        "evidence": "derived: GenomeOS therapeutic pipeline against publicly approved targets",
        "note": (
            "Two questions scored apart: whether the driver's gene is recovered at all, and whether the "
            "route is called correctly. GenomeOS models antibody-like modalities and no small molecules, "
            "so for a target whose approved drug is a small-molecule inhibitor the correct answer is that "
            "no surface route exists; those cases pass by not claiming one. A third question, whether "
            "the mechanism ranked first is defensible, is scored separately and is where the pipeline "
            "currently falls short: see mechanism_concerns per row."
        ),
        "rows": rows,
    }
    save_result("therapeutic_benchmark", result)
    print(
        f"\n{recovered}/{len(rows)} targets recovered; {passed}/{len(rows)} verdicts correct; "
        f"{result['surface_cases_recovered']}/{len(surface_cases)} approved antibody targets "
        f"found as surface targets; {sane}/{len(rows)} top mechanisms defensible",
        file=log,
    )
    return 0 if passed == len(rows) else 1


if __name__ == "__main__":
    raise SystemExit(main())
