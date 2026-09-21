# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""Retrospective benchmark: does the therapeutic pipeline recover known targets?

Seven tumours whose target and whose approved therapy are public knowledge are
run through `genomeos therapeutic`, and the result is compared with what is
actually approved for that alteration. The point is falsifiability: until a
pipeline is asked to reproduce something already known, its rankings are
assertions.

Six of the seven targets are reached through a point mutation, which for a long
time was the whole benchmark and was a narrower test than it read as: the
copy-number, structural and expression routes into the candidate list were
covered by unit tests and by nothing that scored the pipeline end to end. The
seventh, CD19, carries no alteration of any kind and is reached only from the
patient's own RNA. Each case records the `call` that puts its target on the
list, so what the benchmark does not cover is visible in the result.

Two different questions are scored separately, because conflating them would
flatter the pipeline:

  1. **Is the target recovered?** Does the target's gene appear among the
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

DEMO = Path("data/demo")
BENCH = DEMO / "benchmark"

#: Each case: the alteration, the gene an oncologist would name, and what is
#: actually approved for it. `modality` is the class of the approved drug, and
#: `call` is the kind of measurement that puts the target on the list at all —
#: for six of the seven cases a point mutation, which is why the seventh exists.
CASES: tuple[dict, ...] = (
    {
        "case": "EGFR L858R, lung adenocarcinoma",
        "vcf": "egfr_l858r_lung.vcf",
        "call": "point mutation",
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
        "call": "point mutation",
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
        "call": "point mutation",
        "gene": "BRAF",
        "approved": "vemurafenib, dabrafenib (small molecules)",
        "modality": "small_molecule",
        "expect": "out_of_scope_expected",
        "why": "a cytoplasmic kinase; the approved drug is a small molecule, a class GenomeOS does not model",
    },
    {
        "case": "KRAS G12C, lung adenocarcinoma",
        "vcf": "kras_g12c_lung.vcf",
        "call": "point mutation",
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
        "call": "point mutation",
        "gene": "PIK3CA",
        "approved": "alpelisib (small molecule)",
        "modality": "small_molecule",
        "expect": "out_of_scope_expected",
        "why": "intracellular lipid kinase",
    },
    {
        "case": "IDH1 R132H, glioma and AML",
        "vcf": "idh1_r132h_glioma.vcf",
        "call": "point mutation",
        "gene": "IDH1",
        "approved": "ivosidenib (small molecule)",
        "modality": "small_molecule",
        "expect": "out_of_scope_expected",
        "why": (
            "cytosolic enzyme; a neomorphic mutation whose product is a metabolite, "
            "which GenomeOS does not model at all"
        ),
    },
    {
        # The case that is not a point mutation. Six drivers above are single
        # base changes, so until this one existed the benchmark never asked the
        # pipeline to reach a target through anything but a VCF line, and the
        # copy-number, structural and expression routes were covered by unit
        # tests alone. CD19 is the sharpest version of that question: it is not
        # mutated, not amplified and not rearranged in any tumour, so no DNA
        # event can reach it, and it is nonetheless the target of four approved
        # therapies. The measurement that puts it on the list is the patient's
        # own RNA against the healthy-tissue atlas.
        "case": "CD19 in a B-cell lymphoma, reached from the tumour's own RNA",
        "dir": "expression",
        "vcf": "bcell_lymphoma.vcf",
        "rna": "bcell_lymphoma_rna.tsv",
        "scan": 6,
        "call": "expression",
        "gene": "CD19",
        "approved": (
            "tafasitamab (Fc-engineered antibody); blinatumomab (CD19xCD3 engager); "
            "tisagenlecleucel and axicabtagene ciloleucel (CAR-T); loncastuximab "
            "tesirine (antibody-drug conjugate)"
        ),
        "modality": "antibody",
        "expect": "surface",
        "why": (
            "the tumour's only DNA driver is TP53; CD19 carries no alteration at all, so a "
            "pipeline that starts from alterations cannot propose the target of the two "
            "best-known cell therapies. Every approved CD19 drug is antibody-like, which is "
            "the one modality class GenomeOS models, so this case tests the route and not "
            "only the recovery"
        ),
    },
)


def run_case(case: dict, net: bool, log) -> dict:
    folder = DEMO / case.get("dir", "benchmark")
    vcf = folder / case["vcf"]
    cnv = str(folder / case["cnv"]) if case.get("cnv") else None
    sv = str(folder / case["sv"]) if case.get("sv") else None
    rna = str(folder / case["rna"]) if case.get("rna") else None
    t0 = time.time()
    a = analyse_vcf(
        str(vcf),
        cnv=cnv,
        sv=sv,
        rna=rna,
        scan_expression=case.get("scan", 0),
        top_genes=8,
        net=net,
        indirect=True,
        log=None,
    )
    gene = case["gene"]
    ranked = [c.gene for c in a["candidates"]]
    hit = next((c for c in a["candidates"] if c.gene == gene), None)
    out = {
        "case": case["case"],
        "gene": gene,
        "driver_call": case["call"],
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
        # A preferred mechanism is one whose hard requirements were answered.
        # The nearest one whose requirement is merely open is recorded beside
        # it, with the requirement named, so that a row saying "no preferred
        # mechanism" says which question is in the way rather than going quiet.
        nearest = hit.best_provisional_mechanism
        out.update(
            {
                "target_class": hit.target_class,
                "score": None if hit.scores.overall is None else round(hit.scores.overall, 3),
                "surface_accessibility": hit.scores.value("surface_accessibility"),
                "best_mechanism": best.mechanism if best else None,
                "best_compatibility": round(best.compatibility, 3) if best else None,
                "best_mechanism_established": bool(best and best.established),
                "nearest_provisional_mechanism": nearest.mechanism if nearest else None,
                "requirements_unanswered": list(nearest.provisional_requirements) if nearest else [],
                "peptide_route": bool(hit.neoantigen and hit.neoantigen.peptides),
            }
        )
        # A third question, scored apart from the other two: is the mechanism the
        # pipeline puts first defensible? Three ways it is not. An agonist antibody
        # against an activating oncogenic driver would push the pathway the tumour
        # already over-drives. A mechanism that needs an extracellular epitope
        # should not head the list for a protein with no established outward-facing
        # part, however low its compatibility. And a mechanism whose hard
        # requirement was never answered should not head it either: not refused is
        # not established, and the third check is what keeps a row that reports no
        # preferred mechanism from passing this question by going quiet.
        mech = out.get("best_mechanism") or ""
        acc = hit.scores.value("surface_accessibility") or 0
        problems = []
        if "agonist" in mech:
            problems.append(
                f"top mechanism is {mech} against an activating driver, which would drive the pathway"
            )
        if acc <= 0.5 and mech in ("targeted_radionuclide", "adc", "blocking_antibody", "adcc", "adcp"):
            problems.append(f"top mechanism {mech} needs an extracellular epitope, accessibility is {acc}")
        if best is not None and not best.established:
            open_ = ", ".join(best.provisional_requirements) or "an unnamed requirement"
            problems.append(f"top mechanism {mech} heads the list with {open_} unanswered")
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
    by_call: dict[str, int] = {}
    for r in rows:
        by_call[r["driver_call"]] = by_call.get(r["driver_call"], 0) + 1
    result = {
        "result": "therapeutic_benchmark",
        "cases": len(rows),
        "targets_recovered": recovered,
        "verdicts_correct": passed,
        "top_mechanism_defensible": sane,
        "cases_by_driver_call": by_call,
        "surface_cases": len(surface_cases),
        "surface_cases_recovered": sum(1 for r in surface_cases if r["pass"]),
        "evidence": "derived: GenomeOS therapeutic pipeline against publicly approved targets",
        "note": (
            "Two questions scored apart: whether the target's gene is recovered at all, and whether the "
            "route is called correctly. GenomeOS models antibody-like modalities and no small molecules, "
            "so for a target whose approved drug is a small-molecule inhibitor the correct answer is that "
            "no surface route exists; those cases pass by not claiming one. A third question, whether "
            "the mechanism ranked first is defensible, is scored separately: a mechanism is a concern if "
            "it is an agonist against an activating driver, if it needs an extracellular epitope the "
            "protein has no established outward-facing part for, or if it heads the list with a hard "
            "requirement unanswered rather than answered. A row with no preferred mechanism names the "
            "open requirement in requirements_unanswered rather than falling silent. "
            "cases_by_driver_call says which measurement reaches each target: six point mutations and "
            "one expression call, CD19, which carries no alteration of any kind. Copy number and "
            "structural variants reach the ranking and are covered by unit tests, but no scored case "
            "here turns on one, so the benchmark does not yet measure those two routes."
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
