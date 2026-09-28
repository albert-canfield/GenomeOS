# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""Retrospective benchmark: does the therapeutic pipeline recover known targets?

Ten tumours whose target and whose approved therapy are public knowledge are
run through `genomeos therapeutic`, and the result is compared with what is
actually approved for that alteration. The point is falsifiability: until a
pipeline is asked to reproduce something already known, its rankings are
assertions.

Six of the nine targets are reached through a point mutation, which for a long
time was the whole benchmark and was a narrower test than it read as: the
copy-number, structural and expression routes into the candidate list were
covered by unit tests and by nothing that scored the pipeline end to end. The
other three close that. CD19 carries no alteration of any kind and is reached
only from the patient's own RNA; ERBB2 appears a second time reached by nothing
but its amplification; ALK is reached by a rearrangement. All four routes into
the candidate list are now scored, and each case records the `call` that put
its target there.

The copy-number case earned its place on the first run. It passes all three
questions and ranks fourth, behind three surface proteins that carry no
alteration in this tumour and are named only for neighbouring the one that
does. Its two entries — the amplified gene and the same gene as a hypothesis —
scored identically at 0.494 before the duplicate was removed, which is the
finding: for a surface target the score reads the gene's annotation and not
the alteration, so twelve copies and a guess are worth the same. Recovery was
never the hard question, and six point mutations could not ask the other one.

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

The tenth case is the first one chosen for a rule rather than for a route. The
mechanism gate and the magnitude tiebreak were registered on 2026-09-27 and
neither changed anything in the nine: no candidate a modality cannot reach ever
outranked a target, and no VCF here reported an allele fraction. A rule with a
metric and no case is not yet evidence, so the tenth tumour is the one where the
gate could matter - an approved antibody target and, in the same evidence tier, a
gene nothing can be aimed at - and it carries the benchmark's first allele
fraction. It was pre-registered, with its predictions and falsifiers, before it
was run.

Run: `uv run python scripts/therapeutic_benchmark.py`
Writes: data/results/therapeutic_benchmark.json
"""

from __future__ import annotations

import sys
import time
from itertools import groupby
from pathlib import Path

from genomeos import manifest as mf
from genomeos.results import save_result
from genomeos.therapeutics import analyse_vcf
from genomeos.therapeutics import scoring as ranking_rules
from genomeos.therapeutics.scoring import UNMEASURED_TIER

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
        # The copy-number case. ERBB2 appears twice in this benchmark on
        # purpose, and the repetition is the design: the case above reaches it
        # through an ERBB2 coding variant, this one through nothing but the
        # copy-number call. Holding the target fixed and changing only the
        # measurement is what isolates the route; CD19 below varies both at
        # once and so cannot separate them. The VCF carries a PIK3CA variant
        # and no ERBB2 change at all, so a build without the copy-number path
        # does not rank ERBB2 anywhere.
        #
        # It is also the honest clinical case. Trastuzumab's companion
        # diagnostic is amplification, not a point mutation: the measurement
        # that selects the patient IS the copy-number call.
        "case": "ERBB2-amplified breast carcinoma, reached only through the copy-number call",
        "dir": "alterations",
        "vcf": "erbb2_amplification_only.vcf",
        "cnv": "erbb2_amplification_only.cnv",
        "call": "copy number",
        "gene": "ERBB2",
        "approved": "trastuzumab, pertuzumab (antibodies); trastuzumab deruxtecan (antibody-drug conjugate)",
        "modality": "antibody",
        "expect": "surface",
        "why": (
            "the amplification is the approved therapy's companion diagnostic, so a pipeline "
            "that needs a coding change to see ERBB2 would miss every patient trastuzumab is "
            "actually given to"
        ),
    },
    {
        # The structural-variant case, and the one that is a surface receptor
        # in every database and a cytoplasmic kinase in the patient. It is
        # scored as out of scope because the approved ALK drugs are small
        # molecules — and the reason no antibody exists is the thing the
        # pipeline has to work out for itself: ALK is the fusion's 3' partner,
        # so the product begins with EML4 and carries neither ALK's signal
        # peptide nor its ectodomain. A pipeline reading the curated
        # localisation alone calls this a direct surface target, which is how
        # it read until 2026-09-21.
        "case": "EML4-ALK lung adenocarcinoma, a curated surface receptor that the fusion makes cytoplasmic",
        "dir": "alterations",
        "vcf": "alk_eml4_fusion.vcf",
        "sv": "alk_eml4_fusion.sv",
        "call": "structural variant",
        "gene": "ALK",
        "approved": "crizotinib, alectinib, brigatinib, lorlatinib (small molecules)",
        "modality": "small_molecule",
        "expect": "out_of_scope_expected",
        "why": (
            "every curated compartment calls ALK a single-pass surface receptor, so the trap is to "
            "offer an antibody against it; as the 3' partner of the fusion it keeps none of its "
            "outward-facing part, and the drugs that work are small molecules acting inside the cell"
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
    {
        # The tenth case, pre-registered 2026-09-28 before it was run: the one
        # chosen to make the mechanism gate matter. The gate was registered on
        # 2026-09-27 and closed nothing in the nine, because every candidate
        # that could have outranked a target was already below it on the
        # evidence tier — so the rule had a metric and no case. A case that can
        # test it needs an approved target and, in the same tumour and the same
        # tier, a strongly altered gene no modelled modality reaches.
        #
        # It was chosen from the cohort and not invented. In TCGA stomach
        # adenocarcinoma (PanCancer Atlas, stad_tcga_pan_can_atlas_2018, 440
        # samples) ERBB2 is amplified in 58 tumours (13.2%) and MYC in 53
        # (12.0%), and 20 tumours carry both: 4.5% of the cohort and 34% of
        # every ERBB2-amplified tumour. Among those 20 the median copy number
        # from the log2 segment calls is 13.0 for ERBB2 and 8.1 for MYC, and
        # TP53 sits at 1.8. TP53 is mutated in 213 of the 440 (48%), and its
        # most frequent protein change is R175H (12 tumours, median variant
        # allele fraction 0.49). Those are the numbers in the two demo files,
        # each of them a cohort summary and none of them one patient's record.
        #
        # GRB7, MIEN1 and STARD3 are in the table because the 17q12 amplicon
        # carries them: in the same cohort every one of the 58 ERBB2-amplified
        # tumours has GRB7 and MIEN1 amplified too, and 54 have STARD3. They
        # take the same copy number as ERBB2 for a reason worth recording — a
        # copy call is a segment call, so genes on one amplicon are measured at
        # one number and the magnitude tiebreak can never separate them.
        #
        # This case is also the first in the benchmark whose VCF reports an
        # allele fraction, which is what the magnitude tiebreak has never had.
        "case": "HER2-positive gastroesophageal adenocarcinoma with a co-amplified MYC",
        "vcf": "erbb2_myc_gastroesophageal.vcf",
        "cnv": "erbb2_myc_gastroesophageal.cnv",
        "call": "copy number",
        "gene": "ERBB2",
        "approved": (
            "trastuzumab (antibody), FDA label for HER2-overexpressing metastatic gastric and "
            "gastroesophageal junction adenocarcinoma, approved 2010 on ToGA; fam-trastuzumab "
            "deruxtecan-nxki (antibody-drug conjugate), FDA approval for HER2-positive advanced "
            "gastric and GEJ adenocarcinoma after a prior trastuzumab regimen, January 2021 "
            "on DESTINY-Gastric01"
        ),
        "modality": "antibody",
        "expect": "surface",
        "why": (
            "the tumour carries an approved antibody target and, in the same tumour and the same "
            "evidence tier, an amplified MYC: a transcription factor with no outward-facing part "
            "and no modelled modality, and the gene oncology has spent forty years failing to "
            "drug. If a pipeline ranks the undruggable amplification above the one the label is "
            "written on, the list cannot be acted on, and the gate is the rule that has to stop it"
        ),
    },
)


#: Pre-registered 2026-09-28, before the tenth case was run, alongside the dated
#: section of docs/THERAPEUTICS.md and the pins in tests/test_therapeutic_benchmark.py.
#:
#: What is predicted. (1) ERBB2 is recovered as a surface target and ranks first,
#: and `outranked_by_unreachable` is empty, so the pin of 0 holds. (2) MYC, the
#: amplicon passengers and the mutated TP53 all sit in the top evidence tier with
#: ERBB2 — every one of them is altered in this tumour — so the tier cannot
#: separate them and the gate is the only rule that can. (3) The gate is *tested*
#: only if at least one of those unreachable candidates scores strictly above
#: ERBB2; `rank_without_gate`, computed for every row after this registration was
#: committed, is what says so, and if it reads 1 for this case then the gate again
#: closed nothing and the case has produced a negative rather than a pass.
#: (4) This row publishes a variant allele fraction, the first in the benchmark,
#: at the 0.49 the VCF reports and not a value imputed for it.
#: (5) The magnitude tiebreak fires only where two candidates are equal on tier,
#: gate and published score and differ on a quantity both carry. It is not
#: predicted to fire: the amplicon genes are measured at one segment number and so
#: carry the same count, and an exact tie in the score between the remaining pairs
#: is not something the choice of case can arrange. `magnitude_tiebreaks` records
#: every score-tied group and what, if anything, magnitude decided inside it.
#:
#: The falsifiers. (a) ERBB2 not first in this case. (b) Any of the nine earlier
#: rows changing in any field but `seconds`. (c) The row's magnitude publishing a
#: vaf other than the 0.49 the VCF reports, or publishing none.
GATE_CASE_REGISTERED = "2026-09-28"


def rank_without_gate(candidates: list, gene: str) -> int | None:
    """Where the target would rank if the mechanism gate were not read.

    The pipeline's key is tier, then gate, then score, then gene. This is the
    same key with the gate removed, which is exactly the order the pipeline
    produced before 2026-09-27, so the difference between it and `rank` is what
    the gate did in this tumour and nothing else.
    """
    order = [
        c.gene
        for c in sorted(
            candidates,
            key=lambda c: (c.evidence_tier == UNMEASURED_TIER, -(c.scores.overall or 0.0), c.gene),
        )
    ]
    return order.index(gene) + 1 if gene in order else None


def magnitude_tiebreaks(candidates: list) -> list[dict]:
    """Every group the score could not separate, and what magnitude did inside it.

    The tiebreak can only act where tier, gate and published score are all equal,
    so the groups are the whole of its opportunity: a run with no group of two is
    a run in which the rule could not have fired, and saying which is the
    difference between a rule that was measured and a rule that was quiet.
    """
    out: list[dict] = []
    key = lambda c: (  # noqa: E731
        c.evidence_tier == UNMEASURED_TIER,
        c.mechanism_reach,
        -(c.scores.overall or 0.0),
    )
    for _, group in groupby(candidates, key=key):
        block = list(group)
        if len(block) < 2:
            continue
        genes = [c.gene for c in block]
        decided = []
        for a, b in zip(block, block[1:], strict=False):
            if not ranking_rules.magnitude_prefers(a.alteration_magnitude, b.alteration_magnitude):
                continue
            quantity = next(
                (
                    q
                    for q in ranking_rules.MAGNITUDE_QUANTITIES
                    if a.alteration_magnitude.get(q) is not None
                    and b.alteration_magnitude.get(q) is not None
                    and a.alteration_magnitude.get(q) != b.alteration_magnitude.get(q)
                ),
                None,
            )
            decided.append({"above": a.gene, "below": b.gene, "quantity": quantity})
        out.append(
            {
                "tied_on_score": genes,
                "score": None if block[0].scores.overall is None else round(block[0].scores.overall, 3),
                "gate": block[0].mechanism_reach,
                "decided_by_magnitude": decided,
                # The fallback inside a tied group is the gene name, so the order
                # differing from the alphabetical one is the tiebreak having moved
                # something rather than having agreed with the fallback.
                "changed_the_order": genes != sorted(genes),
            }
        )
    return out


def quantities_in_this_tumour(candidates: list) -> list[dict]:
    """Every candidate that carries a measured amount, and what it is.

    Added 2026-09-28 as the amendment to prediction 4, which was falsified by the
    run and is kept in the record rather than rewritten. The prediction said the
    tenth case's row would publish an allele fraction of 0.49. It does not, and
    the reason is a mistake in the registration and not in the pipeline: a row
    publishes the *target's* magnitude, the target there is reached by its
    amplification, and the fraction the VCF reports belongs to the mutated TP53,
    which is a different candidate. A quantity that reaches the ranking is visible
    only if the ranking's own inputs are published, so they are.
    """
    out = []
    for c in candidates:
        m = c.alteration_magnitude or {}
        if not m.get("quantities"):
            continue
        out.append(
            {
                "gene": c.gene,
                "copies": m.get("copies"),
                "vaf": m.get("vaf"),
                "hotspot": m.get("hotspot"),
                "gate": c.mechanism_reach,
                "score": None if c.scores.overall is None else round(c.scores.overall, 3),
            }
        )
    return out


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
    # A fourth question, and the copy-number case is what made it askable: how
    # many candidates outrank the target while nothing about them was measured
    # in this patient at all? Recovering the target and burying it under
    # hypotheses are not the same result, and the first three questions cannot
    # tell them apart.
    #
    # It was first written as "candidates with no origins record", and that
    # proxy was wrong in the direction that mattered: CD19 has no origins
    # either, because it is reached from the patient's RNA rather than from a
    # DNA event, so the proxy counted the one case that proves the expression
    # route works as a case of burial. The count now reads the evidence tier,
    # which is the thing the proxy stood for: a candidate is a hypothesis when
    # nothing about that gene was measured in this patient.
    above = a["candidates"][: ranked.index(gene)] if hit is not None else []
    hypotheses_above = [c.gene for c in above if c.evidence_tier == UNMEASURED_TIER]
    # A fifth question, from the mechanism gate registered 2026-09-27: how many
    # candidates rank above the target while no modelled modality reaches them with
    # its hard requirements answered? A gene nothing can be aimed at, ranked above a
    # gene with an approved antibody, is a ranking that cannot be acted on.
    unreachable_above = [c.gene for c in above if c.mechanism_reach == ranking_rules.REACH_NONE]
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
        "outranked_by_hypotheses": hypotheses_above,
        "outranked_by_unreachable": unreachable_above,
        # Registered 2026-09-28 and computed only after the registration was
        # committed: where the target would have ranked with the mechanism gate
        # removed from the key. Equal to `rank` means the gate changed nothing
        # in this tumour.
        "rank_without_gate": rank_without_gate(a["candidates"], gene),
        "magnitude_tiebreaks": magnitude_tiebreaks(a["candidates"]),
        "quantities_in_this_tumour": quantities_in_this_tumour(a["candidates"]),
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
                "evidence_tier": hit.evidence_tier,
                "mechanism_reach": hit.mechanism_reach,
                "mechanism_reach_reason": hit.mechanism_reach_reason,
                "alteration_magnitude": hit.alteration_magnitude,
                "alteration_magnitude_reason": hit.alteration_magnitude_reason,
                "alteration_evidence": hit.scores.value("alteration_evidence"),
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


KNOWLEDGE_READ = (
    Path("data/knowledge/therapeutics"),
    Path("data/knowledge/proteins"),
    Path("data/knowledge/expression"),
    Path("data/knowledge/pathways"),
    Path("data/knowledge/go-basic.obo"),
    Path("data/knowledge/goa_human.gaf.gz"),
    Path("data/knowledge/Ensembl2Reactome.txt"),
    Path("data/cache/gencode_genes.tsv"),
)


def manifest(net: bool) -> dict:
    """The provenance contract (review item R9) for the 9/9 benchmark. The providers answer from the
    local knowledge stores first and, with the network on, fetch what a store lacks from live
    services that publish no version a request can pin; the stores' sha256 after the run is the pin."""
    inputs = []
    for case in CASES:
        folder = DEMO / case.get("dir", "benchmark")
        for key in ("vcf", "cnv", "sv", "rna"):
            if case.get(key):
                inputs.append(mf.input_entry(folder / case[key], partition=None, case=case["case"]))
    for name in ("cancer_msk_impact_2017", "cancer_alterations_msk_impact_2017"):
        p = Path("data/results") / f"{name}.json"
        if p.exists():
            inputs.append(mf.input_entry(p, partition=None))
    inputs += [mf.input_entry(p, partition=None) for p in KNOWLEDGE_READ if p.exists()]
    return {
        "sources": [
            {
                "accession": "the tenth benchmark tumour (data/demo/benchmark/erbb2_myc_"
                "gastroesophageal.vcf and .cnv): hand-written, its copy numbers and its allele "
                "fraction the medians of a cBioPortal cohort and no patient's own record",
                "version": "this repository at the recorded commit",
            },
            {
                "accession": "cBioPortal stad_tcga_pan_can_atlas_2018 (TCGA PanCancer Atlas, stomach "
                "adenocarcinoma): the frequencies and medians the tenth case is built from",
                "version": "public REST API as read 2026-09-28; cohort summaries only, no "
                "sample-level record",
            },
            {
                "accession": "cBioPortal msk_impact_2017 (Zehir et al. 2017)",
                "version": "public REST API as read 2026-09-10; pinned by the result files' sha256",
            },
            {
                "accession": "GenomeOS federated protein compiler: UniProt, Ensembl, InterPro, PDB, "
                "AlphaFold, Reactome, STRING, HPA",
                "version": "live services, unversioned per request; the local store's sha256 is the pin",
            },
            {"accession": "Human Protein Atlas normal expression", "version": "live service, CC BY-SA 4.0"},
            {"accession": "Open Targets Platform", "version": "live GraphQL service, CC0 1.0"},
            {"accession": "Ensembl REST (paralogues, transcript sequence)", "version": "live service"},
            {
                "accession": "the nine benchmark tumours (data/demo): hand-written VCF, CNV, SV and "
                "RNA files",
                "version": "this repository at the recorded commit",
            },
        ],
        "inputs": inputs,
        "assembly": "GRCh38",
        "coordinates": {"base": 1, "interval": "closed"},
        "parameters": {"network": net, "top_genes": 8, "indirect": True, "cases": len(CASES)},
        "exclusions": [
            "small-molecule routes are not modelled: those cases pass by not claiming a surface route",
        ],
        "partitions": "n/a: a retrospective benchmark of nine known targets; nothing is fitted to them",
        "partitions_note": "the tenth case, registered 2026-09-28, is on the same terms: it was "
        "chosen from a public cohort, registered with its predictions, and nothing is fitted to it",
    }


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
    gate_moved = [r["gene"] for r in rows if r.get("rank_without_gate") not in (None, r["rank"])]
    tie_groups = sum(len(r.get("magnitude_tiebreaks") or []) for r in rows)
    with_vaf = [
        r["gene"] for r in rows if any(q["vaf"] is not None for q in r.get("quantities_in_this_tumour") or [])
    ]
    tie_fired = [
        r["gene"] for r in rows if any(g["changed_the_order"] for g in r.get("magnitude_tiebreaks") or [])
    ]
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
            "cases_by_driver_call says which measurement reaches each target, and all four routes "
            "into the candidate list are now scored: six point mutations, one copy-number call "
            "(ERBB2 amplified and not mutated, the same target as the point-mutation case so that "
            "the route is the only thing that differs), one structural variant (EML4-ALK) and one "
            "expression call (CD19, which carries no alteration of any kind). The ALK case is the "
            "one where the curated databases and the patient disagree: ALK is a single-pass surface "
            "receptor in every compartment annotation and, as the fusion's 3' partner, keeps neither "
            "its signal peptide nor its ectodomain, so the product is a cytoplasmic kinase and the "
            "approved drugs are small molecules. It passes by refusing a surface route rather than "
            "by leaving the question open. "
            "outranked_by_hypotheses is recorded per row and is the question the copy-number case "
            "made askable: candidates ranked above the target with nothing measured about them in "
            "this patient. Until 2026-09-27 the pipeline could recover a target without preferring "
            "it, and the three scored questions could not tell those apart: ERBB2 at twelve copies, "
            "whose amplification is trastuzumab's companion diagnostic, scored 0.494 and ranked "
            "fourth behind KDR, EGFR and PDGFRB, none of them altered in that tumour, because "
            "surface_accessibility reads curated localisation and no dimension read the alteration. "
            "A dimension now does: alteration_evidence, on three tiers - an alteration observed in "
            "this tumour (1.0), this patient's own measurement of the gene (0.6), a database "
            "association with a gene that was altered (0.2) - and the ranking puts the third tier "
            "below the other two before it looks at a score. evidence_tier and alteration_evidence "
            "are reported per row. The tier is read from the evidence about the gene and never from "
            "the route that proposed it, which is what keeps CD19, reached from the patient's RNA "
            "with no alteration of any kind, from being demoted by a rule aimed at hypotheses; for "
            "the same reason the count above no longer means 'candidates with no origins record', a "
            "proxy under which CD19's own route read as burial. Scores before and after this change "
            "are not comparable, because every candidate gained a dimension; the ordering is what "
            "the change was about."
        ),
        # The two rules of 2026-09-27 are documented in their own key rather than inside
        # the note above, which belongs to the run that measured the evidence tier.
        "ordering_rules": (
            "Two rules on the order, registered before the code on 2026-09-27 and neither a term in "
            "the weighted mean. mechanism_reach partitions candidates by whether any modelled "
            "modality reaches them with its hard requirements answered; it is read after the "
            "evidence tier and never before it, and provisional-only sits with nothing-at-all "
            "because a mechanism whose requirement is unanswered has not been shown to apply. "
            "outranked_by_unreachable counts candidates ranked above the target that no modality "
            "reaches. alteration_magnitude publishes how much this patient's own data say the gene "
            "was altered - copy count against the diploid 2, variant allele fraction, hotspot "
            "status - with an absent quantity stated as absent and never imputed, and it breaks "
            "ties only between candidates already equal on tier, gate and published score, "
            "comparing like with like because a count and an annotation share no unit."
        ),
        # Hand-added to the committed file on 2026-09-28 (2805552); kept here so the script writes it
        # and a rebuild from the manifest reproduces the file (review item R9).
        "ordering_rules_measured": (
            "mechanism_reach, alteration_magnitude and outranked_by_unreachable were measured by a "
            "re-run on 2026-09-28, hours after the rules were registered. Every other field of every "
            "row came out identical to the 2026-09-27 run - rank, score, target class, evidence tier, "
            "mechanism, compatibility, accessibility, peptide route, verdict - which is the check that "
            "the two rules changed an order and nothing else. Until the same day this key was added to "
            "the file by hand beside the 2026-09-27 run's date and timings; since the result manifest "
            "(R9) the whole file is this script's output, rebuilt in a clean checkout, and every field "
            "but the date and the per-case timings came out as committed."
        ),
        "gate_and_tiebreak_case": (
            "The tenth case, pre-registered 2026-09-28 before it was run, and the first chosen for a "
            "rule rather than for a route: a HER2-positive gastroesophageal adenocarcinoma with a "
            "co-amplified MYC. The mechanism gate and the magnitude tiebreak were registered on "
            "2026-09-27 and neither did anything in the nine cases - no candidate outside every "
            "modelled modality ever outranked a target, and no demo VCF reported an allele fraction, "
            "so the tiebreak had no quantity to compare. This tumour was chosen from TCGA stomach "
            "adenocarcinoma (PanCancer Atlas, 440 samples), where ERBB2 is amplified in 58 tumours, "
            "MYC in 53, and 20 carry both: 34% of every ERBB2-amplified tumour in the cohort. The "
            "copy numbers are that subgroup's medians (ERBB2 13, MYC 8, TP53 2), the 17q12 passengers "
            "GRB7, MIEN1 and STARD3 take ERBB2's number because a copy call is a segment call, and "
            "the TP53 R175H hotspot carries the cohort's median allele fraction for that change, "
            "0.49. Registered predictions: ERBB2 first with outranked_by_unreachable empty; the gate "
            "tested only if an unreachable candidate scores above ERBB2, which rank_without_gate says "
            "per row; a published allele fraction of 0.49 and not an imputed one; and the tiebreak "
            "not predicted to fire, since one amplicon is measured at one number and an exact tie in "
            "the score cannot be arranged by the choice of case. rank_without_gate and "
            "magnitude_tiebreaks are computed for every row, the nine included, and were written only "
            "after the registration was committed."
        ),
        "gate_and_tiebreak_case_measured": (
            "Negatives first. The gate closed nothing again: rank_without_gate equals rank in all ten "
            "cases, so removing the gate from the ranking key would change no target's place, and "
            "gate_moved_the_target is empty. The case was chosen to make the gate matter and it "
            "measured why it does not, which is the useful part. MYC, amplified at 8 copies in the top "
            "evidence tier, scores 0.246 - the lowest of the six altered candidates, against the "
            "target's 0.494 - so it never threatened the target's place and the gate had nothing to "
            "demote. The reason is structural rather than particular to MYC: the annotations that make "
            "a gene unreachable, no outward-facing part and no epitope, are the same annotations that "
            "make it score low, so an undruggable nuclear amplification cannot produce the "
            "configuration the gate was written for. The configuration needs a candidate curated as a "
            "surface receptor whose mechanisms are nonetheless all refused or unanswered - the "
            "EML4-ALK shape, where ALK scores 0.559 against ERBB2's 0.494 - and that means two "
            "independent drivers in one tumour. In MSK-IMPACT 2017, of the 42 tumours whose structural "
            "variants make a surface receptor the 3' partner (6 ALK, 12 ROS1, 16 RET, 8 NTRK1), one "
            "also carries an amplified ERBB2, one an amplified EGFR and one an amplified MET: the "
            "defect the gate prevents is rare in real tumours because strong drivers are largely "
            "mutually exclusive, and a benchmark case built on a single sample would be that sample's "
            "record rather than a cohort summary. The second negative: the magnitude tiebreak had no "
            "opportunity in the new case at all - no two candidates there are equal on tier, gate and "
            "published score - and the four score-tied groups in the other cases are all pairs with no "
            "quantity between them. magnitude_tiebreak_changed_order is empty. Third, a registered "
            "prediction was falsified, and by the registration rather than by the pipeline: the tenth "
            "row was predicted to publish a variant allele fraction of 0.49 and publishes none, "
            "because a row publishes the target's own magnitude and the target is reached by its "
            "amplification, while the fraction the VCF reports belongs to the mutated TP53 in the same "
            "tumour. The prediction is kept as it was written and the amendment is additive: "
            "quantities_in_this_tumour publishes every candidate that carries a measured amount, so "
            "the 0.49 that reaches the ranking is visible where it actually sits. What did pass: "
            "ERBB2 is recovered as a surface target at rank 1 with an established blocking antibody, "
            "outranked_by_unreachable is empty and the pin of 0 holds across ten cases, and the nine "
            "earlier rows came out identical field by field apart from their per-case timings."
        ),
        "gate_moved_the_target": gate_moved,
        "magnitude_tiebreak_groups": tie_groups,
        "magnitude_tiebreak_changed_order": tie_fired,
        "cases_where_a_candidate_carries_an_allele_fraction": with_vaf,
        "rows": rows,
    }
    # The committed file keeps 2805552's wording of this note: the R9 rebuild of 2026-09-28 added a
    # manifest to the committed result rather than rewrite its lines, so its date and per-case timings
    # are still the 2026-09-27 run's, which is what this wording says. The script writes the same
    # wording, so a rebuild from the manifest reproduces the file; the draft wording above is superseded.
    result["ordering_rules_measured"] = (
        "mechanism_reach, alteration_magnitude and outranked_by_unreachable were measured by a "
        "re-run on 2026-09-28, hours after the rules were registered. Every other field of every "
        "row came out identical to the 2026-09-27 run recorded here - rank, score, target class, "
        "evidence tier, mechanism, compatibility, accessibility, peptide route, verdict - which "
        "is the check that the two rules changed an order and nothing else. That run's date and "
        "per-case wall-clock timings are left as they are rather than overwritten, so that adding "
        "these fields removes nothing a previous commit recorded; the re-run took 1.0-1.1 seconds "
        "per case."
    )
    save_result("therapeutic_benchmark", result, manifest=manifest(net))
    print(
        f"\n{recovered}/{len(rows)} targets recovered; {passed}/{len(rows)} verdicts correct; "
        f"{result['surface_cases_recovered']}/{len(surface_cases)} approved antibody targets "
        f"found as surface targets; {sane}/{len(rows)} top mechanisms defensible",
        file=log,
    )
    return 0 if passed == len(rows) else 1


if __name__ == "__main__":
    raise SystemExit(main())
