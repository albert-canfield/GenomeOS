# SPDX-License-Identifier: AGPL-3.0-or-later
"""Item 13 C4, the diagnostic run first: remove every AlphaGenome-derived input and ask which of the
project's conclusions survive on experimental evidence alone.

    uv run python scripts/c4_ablation.py

Registered 2026-09-28 before it was run (genomeos/attribution/ablation.py, docs/ATTRIBUTION.md "Item 13
C4 registered"). No model request; the per-element response cache is never opened. The model's own
labels are read from the compiled programs (data/knowledge/compiled), which is what the project states.
Writes data/results/c4_alphagenome_ablation.json.
"""

from __future__ import annotations

import time
from collections import Counter
from pathlib import Path
from typing import Any

from genomeos import manifest as mf
from genomeos.attribution import ablation as ab
from genomeos.attribution import holdout as ho
from genomeos.attribution import measured as ms
from genomeos.results import RESULTS_DIR, save_result

NAME = "c4_alphagenome_ablation"


def _counts(c: Counter) -> list[dict[str, Any]]:
    return [
        {"class": k[0], "reason": k[1], "n": n}
        for k, n in sorted(c.items(), key=lambda kv: (-kv[1], str(kv[0])))
    ]


def _axes(side: dict[str, Counter]) -> dict[str, list[dict[str, Any]]]:
    out = {}
    for axis, c in side.items():
        out[axis] = [
            {"value": v, "class": cls, "reason": reason, "n": n}
            for (v, cls, reason), n in sorted(c.items(), key=lambda kv: (-kv[1], str(kv[0])))
        ]
    return out


def _by_class(side: dict[str, Counter]) -> dict[str, dict[str, int]]:
    """Axis values per class, over every axis of one block kind."""
    out: dict[str, Counter] = {}
    for axis, c in side.items():
        agg: Counter = Counter()
        for (_, cls, reason), n in c.items():
            agg[cls if reason is None else f"{cls}/{reason}"] += n
        out[axis] = agg
    return {a: dict(sorted(c.items())) for a, c in out.items()}


@mf.depends_on_models("alphagenome")  # it reads the compiled predicted layer to classify it
def manifest(programs: list[Path], headlines: list[dict[str, Any]]) -> dict[str, Any]:
    inputs = [mf.input_entry(p, partition=None) for p in programs]
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
    inputs += [
        mf.input_entry(p, partition=None)
        for c in ho.CHROMS
        if (p := RESULTS_DIR / f"ccres_{c}.bed.gz").exists()
    ]
    named = (
        {"manifest_headlines"}
        | {h["result"] for h in headlines}
        | {r["arm"]["result"] for r in headlines if r.get("arm")}
    )
    inputs += [mf.input_entry(RESULTS_DIR / f"{n}.json", partition=None) for n in sorted(named)]
    return {
        "sources": [
            {
                "accession": "this repository: the 24 compiled programs "
                "data/knowledge/compiled/noncoding_<chrom>.bio (generated, untracked)",
                "version": "pinned by sha256",
            },
            {
                "accession": "EngreitzLab/CRISPR_comparison EPCrisprBenchmark (Gschwind et al. 2025)",
                "version": "main, as fetched; pinned by sha256",
            },
            {
                "accession": f"ENCODE {ho.mpra.LIBRARY} lentiMPRA",
                "version": "as fetched 2026-09-12; pinned by sha256",
            },
            {"accession": "VISTA Enhancer Browser locus table", "version": "pinned by sha256"},
            {"accession": "GEO GSE126550 (Kircher et al. 2019)", "version": "pinned by sha256"},
            {
                "accession": "GTEx v8 single-tissue cis-eQTL, distilled by genomeos/attribution/eqtl.py",
                "version": "pinned by sha256",
            },
            {"accession": "GENCODE v50 and ENCODE cCREs v3 (per chromosome)", "version": "pinned by sha256"},
            {
                "accession": "this repository: the headline results and their experimental arms",
                "version": "pinned by sha256",
            },
        ],
        "inputs": inputs,
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "reciprocal_overlap": ms.RECIPROCAL_OVERLAP,
            "power_for_negatives": ms.POWER_FOR_NEGATIVES,
            "well_powered": ms.WELL_POWERED,
            "mpra_active_log2": ms.MPRA_ACTIVE,
            "classification_rule": ab.CLASSIFICATION_RULE,
            "registered": ab.REGISTERED,
        },
        "exclusions": [
            "CRISPRi pairs the benchmark marks as not a valid connection are not read",
            "the AlphaGenome per-element response cache is not opened; the model is not called",
        ],
        "partitions": {
            ms.TRAINING: "the CRISPRi benchmark's training file (K562); read here as measurement, not fitted",
            ms.HELDOUT: "the CRISPRi benchmark's held-out file; read here as measurement, reported apart",
        },
    }


def main() -> None:
    t0 = time.time()
    programs = ho.compiled_programs()
    pairs, invalid = ms.load_crispri()
    crispri = ab.CrispriIndex(pairs)
    other = ab.OtherEndpoints(ho.all_units())
    headlines = ab.headline_census()
    print("headlines", [(h["result"], h["cls"], h["reason"]) for h in headlines], flush=True)
    census = ab.program_census(programs, crispri, other)
    place = ab.placement(pairs, census["predicted_intervals"], ab.ccre_intervals())
    n_pred = sum(census["elements"].values())
    rules_total = sum(census["rules"].values())
    head_counts = Counter(h["cls"] for h in headlines)
    payload = {
        "question": "with every AlphaGenome-derived input removed, which of the project's conclusions "
        "survive on experimental evidence alone",
        "registration": ab.registration(),
        "alphagenome_requests": 0,
        "per_element_response_cache_opened": False,
        "headlines": {"results": headlines, "by_class": dict(sorted(head_counts.items()))},
        "elements": {
            "predicted_elements": n_pred,
            "link_by_class": _counts(census["elements"]),
            "link_by_detail": dict(census["element_detail"].most_common()),
            "other_endpoints_beside_links_that_do_not_survive": {
                k: dict(v.most_common())
                for k, v in sorted(census["beside"].items(), key=lambda kv: str(kv[0]))
            },
            "axis_values_by_class": {
                "predicted_blocks": _by_class(census["axes"]["predicted"]),
                "measured_blocks": _by_class(census["axes"]["measured"]),
            },
            "axis_values": {
                "predicted_blocks": _axes(census["axes"]["predicted"]),
                "measured_blocks": _axes(census["axes"]["measured"]),
            },
            "region_blocks": dict(census["regions"]),
            "region_blocks_class": ab.NEVER_MODEL_DEPENDENT,
        },
        "rules": {
            "total": rules_total,
            "by_class": _counts(census["rules"]),
            "predicted_by_detail": dict(census["rule_detail"].most_common()),
            "experimental": dict(census["experimental_rules"]),
        },
        "placement": {
            "what": "measured regulatory links (CRISPRi significant decreases as element, gene, cell, split) "
            "and whether a model-compiled element met the overlap rule to carry them; for the rest, whether "
            "an ENCODE cCRE would, with no model involved",
            **place,
        },
        "crispri_invalid_pairs_not_read": invalid,
        "programs": len(programs),
        "status": "internal development reading; every input has been read by this project before",
        "seconds": round(time.time() - t0, 1),
    }
    p = save_result(NAME, payload, manifest=manifest(programs, headlines))
    print("saved", p, flush=True)


if __name__ == "__main__":
    main()
