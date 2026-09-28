# SPDX-License-Identifier: AGPL-3.0-or-later
"""Item 12 S3 / item 13, the coherence pilot's second gate: the debugger run on the validation
chromosomes with every source but one, its revised labels scored on the held-out source against the
unchanged labels, the independent-block variant and distance to TSS; validated corrections per
compute-hour beside, with errors and abstentions; S4's target accuracy, role accuracy and coverage
beside that.

    uv run python scripts/pilot_biological_gate.py

Registered 2026-09-29 before it was run (genomeos/attribution/pilot_bio.py; docs/ATTRIBUTION.md "Item
13 pilot, gate 2 registered"). An internal development result on withheld sources, never a
validation. No model request; the per-element response cache is never opened. Writes
data/results/pilot_biological_gate.json.
"""

from __future__ import annotations

import resource
import sys
import time
from collections import defaultdict
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import correctness as cx  # noqa: E402
from genomeos.attribution import holdout as ho  # noqa: E402
from genomeos.attribution import measured as ms  # noqa: E402
from genomeos.attribution import pilot as pl  # noqa: E402
from genomeos.attribution import pilot_bio as pb  # noqa: E402
from genomeos.results import save_result  # noqa: E402

NAME = "pilot_biological_gate"
EXAMPLES_PER_STATUS = 3


@mf.depends_on_models("alphagenome")  # the unchanged labels and the pilot's target prior: the compiled layer
def manifest(parameters: dict[str, Any]) -> dict[str, Any]:
    inputs = [mf.input_entry(p, partition=None) for p in ho.compiled_programs()]
    for name in ms.CRISPRI_FILES:
        p = ms.CRISPRI_KNOWLEDGE / name
        if p.exists():
            inputs.append(mf.input_entry(p, partition=ms.CRISPRI_SPLIT_OF[name]))
    for acc in ho.mpra.FILES.values():
        p = ho.mpra.KNOWLEDGE / f"{acc}.bed.gz"
        if p.exists():
            inputs.append(mf.input_entry(p, partition=None))
    for p in (ho.vista.locus_path(ho.vista.KNOWLEDGE), ms.satmut_knowledge.DATA_PATH):
        if p.exists():
            inputs.append(mf.input_entry(p, partition=None))
    inputs += [mf.input_entry(p, partition=None) for p in sorted(ho.GTEX_KNOWLEDGE.glob("hits_*.tsv"))]
    inputs += [mf.input_entry(p, partition=None) for p in ho.gencode_paths()]
    for c in pb.VALIDATION:
        inputs.append(mf.input_entry(ROOT / pb.CCRE_TEMPLATE.format(chrom=c), partition="validation"))
        for s in pl.BREADTH_BIOSAMPLES:
            p = ROOT / pb.PEAK_TEMPLATE.format(sample=s, chrom=c)
            if p.exists():
                inputs.append(mf.input_entry(p, partition="validation"))
    return {
        "sources": [
            {
                "accession": "EngreitzLab/CRISPR_comparison EPCrisprBenchmark (Gschwind et al. 2025)",
                "version": "pinned by sha256",
            },
            {"accession": f"ENCODE {ho.mpra.LIBRARY} lentiMPRA", "version": "pinned by sha256"},
            {"accession": "VISTA Enhancer Browser locus table", "version": "pinned by sha256"},
            {"accession": "GEO GSE126550 (Kircher et al. 2019)", "version": "pinned by sha256"},
            {
                "accession": "GTEx v8 single-tissue cis-eQTL, distilled by genomeos/attribution/eqtl.py",
                "version": "pinned by sha256",
            },
            {"accession": "GENCODE v50 (per chromosome)", "version": "pinned by sha256"},
            {
                "accession": "ENCODE cCREs v3 (downloads.wenglab.org/V3/GRCh38-cCREs.bed), per chromosome",
                "version": "pinned by sha256",
            },
            {
                "accession": "ENCODE H3K27ac replicated peaks, 13 biosamples, per chromosome",
                "version": "pinned by sha256",
            },
            {
                "accession": "this repository: data/knowledge/compiled/noncoding_<chrom>.bio",
                "version": "pinned by sha256",
            },
        ],
        "inputs": inputs,
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": parameters,
        "exclusions": [
            "only the validation chromosomes are run and scored; the development ones never are",
            "CRISPRi increases, underpowered nulls and missing effects are not the pilot's evidence; counted",
            "GTEx non-associations are not evidence for the pilot and are counted",
            "the CRISPRi benchmark's held-out file is never evidence (R5)",
            "WTC11_DC_TAP is not held out: its only scorable endpoint (increase) is below the floors on the "
            "validation chromosomes",
            "the AlphaGenome per-element response cache is not opened; the model is not called",
        ],
        "partitions": {
            "validation": "chromosomes "
            + ", ".join(pb.VALIDATION)
            + ": the pilot is run and scored here only",
            "development": "the other chromosomes: timing and debugging only, never scored",
            "held_out_source": "each source in turn; the rest, masked by locus, is the only evidence a "
            "labelling reads",
            ms.TRAINING: "the CRISPRi benchmark's training file (K562): evidence when another source is out",
            ms.HELDOUT: "the CRISPRi benchmark's held-out file: evaluation only, never evidence",
        },
    }


def main() -> None:
    t0, c0 = time.time(), time.process_time()
    fixed = pb.Fixed.load(pb.VALIDATION)
    units = ho.all_units()
    refs = ho.crispri_references()
    distance = ho.distance_labels()
    unchanged = ho.unchanged_labels()
    load_cpu = time.process_time() - c0
    cache: dict[str, pl.Outcome] = {}
    # the unchanged labels' S4 claims on the validation chromosomes, read once
    unchanged_claims = [
        c for c in cx.compiled_claims(chroms=pb.VALIDATION) if c.axis in (cx.TARGET, cx.CONTEXT, cx.ACTIVITY)
    ]
    per_source: dict[str, Any] = {}
    scores: list[dict[str, Any]] = []
    comparisons: list[dict[str, Any]] = []
    pilot_cpu = 0.0
    scoring_cpu = 0.0
    totals: dict[str, int] = defaultdict(int)
    examples: list[dict[str, Any]] = []
    for src in pb.HELD_OUT:
        ts = time.process_time()
        view = ho.evidence(src, units, refs)
        resp = pb.responsiveness(view)
        held = [u for u in units[src] if u.chrom in pb.VALIDATION]
        solved: dict[str, dict[str, pb.Solved]] = {m: {} for m in ("joint", "independent", "prior")}
        skipped_total: dict[str, int] = defaultdict(int)
        obs_count = 0
        for chrom in pb.VALIDATION:
            obs, skipped = pb.view_obs(view, fixed, chrom)
            obs_count += len(obs)
            for k, v in skipped.items():
                skipped_total[k] += v
            rel = pb.relevant_blocks([u for u in held if u.chrom == chrom], fixed.blocks[chrom])
            for mode in solved:
                solved[mode][chrom] = pb.solve(mode, chrom, rel, obs, fixed, resp, cache)
        cpu_modes = {m: round(sum(s.cpu for s in v.values()), 2) for m, v in solved.items()}
        src_pilot_cpu = time.process_time() - ts
        pilot_cpu += src_pilot_cpu
        reads = pb.READS_ALWAYS | frozenset(view.sources)
        pilot = pb.labels_from(pb.PILOT, solved["joint"], fixed, reads, src, note="the joint debugger")
        independent = pb.labels_from(
            pb.INDEPENDENT, solved["independent"], fixed, reads, src, note="one block at a time, no coupling"
        )
        prior = pb.labels_from(pb.PRIOR, solved["prior"], fixed, pb.READS_ALWAYS, src, note="priors only")
        pilot_cov = pb.labels_from(
            pb.PILOT_COVERAGE,
            solved["joint"],
            fixed,
            reads,
            src,
            restrict=unchanged,
            note="at unchanged coverage",
        )
        tsc = time.process_time()
        for ep in [ep for s, ep in (*pb.PRIMARY, *pb.SECONDARY) if s == src]:
            res = ho.evaluate(
                [pilot, independent, prior, pilot_cov, unchanged, distance],
                src,
                ep,
                units=units,
                chroms=pb.VALIDATION,
                references=refs,
                keep_draws=True,
            )
            by = {r["labels"]: r for r in res}
            for a, b in pb.COMPARISONS:
                d = ho.difference(by[a], by[b])
                d["primary"] = (src, ep) in pb.PRIMARY
                comparisons.append(d)
            for r in res:
                r.pop("_draws", None)
                r["primary"] = (src, ep) in pb.PRIMARY
                scores.append(r)
            print(src, ep, [(r["labels"], r.get("value"), r["coverage"]["share"]) for r in res], flush=True)
        scoring_cpu += time.process_time() - tsc
        vep = pb.VALIDATION_ENDPOINT[ho.assay_of(src)]
        val_joint = pb.validate(solved["joint"], held, vep)
        val_ind = pb.validate(solved["independent"], held, vep)
        for k in pb.METRIC_KEYS:
            totals[k] += val_joint[k]
        # S4 beside: the pilot's committed labels and the unchanged labels on the same blocks
        s4_pilot = cx.judge(pb.claims_of(solved["joint"]), pilot, units=units, sources=[src], references=refs)
        s4_unch = cx.judge(
            pb.unchanged_claims_on(solved["joint"], unchanged_claims),
            cx.unchanged_labels(),
            units=units,
            sources=[src],
            references=refs,
        )
        # returned neighbourhoods: the first disputed ones of each status, in genome order
        picked: dict[str, int] = defaultdict(int)
        for chrom in pb.VALIDATION:
            for hood, res_ in solved["joint"][chrom].hoods:
                st = res_.status
                if res_.disputed and st in ("resolved", "abstained") and picked[st] < EXAMPLES_PER_STATUS:
                    picked[st] += 1
                    examples.append({"held_out": src, **pb.neighbourhood_report(hood, res_)})
        s4 = {"pilot": s4_pilot.to_dict(), "unchanged_on_same_blocks": s4_unch.to_dict()}
        for rep in s4.values():
            rep["judged"] = rep["judged"][:30]
        per_source[src] = {
            "family": ho.FAMILY[ho.assay_of(src)],
            "units_on_validation_chromosomes": len(held),
            "sources_in_view": sorted(view.sources),
            "masked_in_view": view.masked,
            "observations_on_validation_chromosomes": obs_count,
            "not_used": dict(skipped_total),
            "neighbourhoods": {m: sum(len(s.hoods) for s in v.values()) for m, v in solved.items()},
            "solved_fresh": {m: sum(s.solved for s in v.values()) for m, v in solved.items()},
            "reused": {m: sum(s.reused for s in v.values()) for m, v in solved.items()},
            "capped": {m: sum(s.capped for s in v.values()) for m, v in solved.items()},
            "cpu_seconds_by_mode": cpu_modes,
            "pilot_cpu_seconds": round(src_pilot_cpu, 2),
            "corrections": {"pilot": val_joint, "independent_block": val_ind},
            "s4": s4,
        }
        brief = {k: val_joint[k] for k in ("committed", "validated", "errors", "neutral", "untested")}
        print(src, "corrections", brief, round(src_pilot_cpu, 1), "cpu s", flush=True)
    verdict = pb.assess(comparisons, scores)
    hours = pb.cpu_hours(pilot_cpu)
    metric = {
        **{k: totals[k] for k in pb.METRIC_KEYS},
        "pilot_cpu_hours": round(hours, 4),
        "validated_per_cpu_hour": round(totals["validated"] / hours, 2) if hours else None,
        "errors_per_cpu_hour": round(totals["errors"] / hours, 2) if hours else None,
        "scoring_cpu_seconds": round(scoring_cpu, 1),
        "loading_cpu_seconds": round(load_cpu, 1),
        "within_budget": hours <= pb.BUDGET_CPU_HOURS,
    }
    peak_mb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1024 * 1024)
    payload = {
        "question": "held out one source at a time, do the pilot's jointly revised labels predict its own "
        "endpoint better than the unchanged labels, the same revision one block at a time, and distance "
        "to TSS, on the validation chromosomes",
        "status": pb.STATUS,
        "may_be_called": "an internal development result on withheld sources",
        "may_not_be_called": ["a validation", *ho.MAY_NOT_BE_CALLED],
        "registration": pb.registration(),
        "verdict": verdict,
        "stop_rule": pb.STOP_RULE,
        "metric": metric,
        "sources": per_source,
        "scores": scores,
        "comparisons": comparisons,
        "neighbourhood_examples": examples,
        "alphagenome_requests": 0,
        "per_element_response_cache_opened": False,
        "compute": {
            "total_cpu_seconds": round(time.process_time() - c0, 1),
            "wall_seconds": round(time.time() - t0, 1),
            "peak_rss_mb": round(peak_mb, 1),
            "measured_with": "time.process_time and getrusage, this process",
        },
    }
    p = save_result(NAME, payload, manifest=manifest(pb.registration()))
    print("verdict", verdict["reading"], "metric", metric, flush=True)
    print("saved", p, flush=True)


if __name__ == "__main__":
    main()
