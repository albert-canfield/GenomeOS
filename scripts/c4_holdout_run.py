# SPDX-License-Identifier: AGPL-3.0-or-later
"""Item 13 C4, the first run of the leave-one-source-out harness: every holdable source held out in
turn, its endpoint predicted by the remaining sources (`rest`) and by the two baselines (distance to
TSS; the unchanged labels), each scored with a 95% locus-bootstrap interval, and the three paired
differences.

    uv run python scripts/c4_holdout_run.py

Registered 2026-09-28 before any score was read (genomeos/attribution/holdout.py, docs/ATTRIBUTION.md
"Item 13 C4 registered"). An internal development benchmark only. No model request; the per-element
response cache is never opened. Writes data/results/c4_holdout_scores.json.
"""

from __future__ import annotations

import time
from typing import Any

from genomeos import manifest as mf
from genomeos.attribution import holdout as ho
from genomeos.attribution import measured as ms
from genomeos.results import save_result

NAME = "c4_holdout_scores"
#: the paired differences the first run reports for every scored endpoint, registered with the rest
COMPARISONS = ((ho.REST, ho.DISTANCE), (ho.REST, ho.UNCHANGED), (ho.UNCHANGED, ho.DISTANCE))


@mf.depends_on_models("alphagenome")  # the unchanged labels are the compiled predicted layer
def manifest(programs: list, parameters: dict[str, Any]) -> dict[str, Any]:
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
    return {
        "sources": [
            {
                "accession": "EngreitzLab/CRISPR_comparison EPCrisprBenchmark (Gschwind et al. 2025)",
                "version": "main, as fetched; pinned by sha256",
            },
            {
                "accession": f"ENCODE {ho.mpra.LIBRARY} lentiMPRA ({', '.join(ho.mpra.FILES.values())})",
                "version": "as fetched 2026-09-12; pinned by sha256",
            },
            {"accession": "VISTA Enhancer Browser locus table", "version": "pinned by sha256"},
            {"accession": "GEO GSE126550 (Kircher et al. 2019)", "version": "pinned by sha256"},
            {
                "accession": "GTEx v8 single-tissue cis-eQTL, distilled by genomeos/attribution/eqtl.py",
                "version": "pinned by sha256",
            },
            {"accession": "GENCODE v50 (per chromosome)", "version": "pinned by sha256"},
            {
                "accession": "this repository: the 24 compiled programs "
                "data/knowledge/compiled/noncoding_<chrom>.bio (generated, untracked)",
                "version": "pinned by sha256",
            },
        ],
        "inputs": inputs,
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": parameters,
        "exclusions": [
            "CRISPRi pairs the benchmark marks as not a valid connection are not read",
            "CRISPRi underpowered nulls and missing effects are counted per endpoint, never scored (R2)",
            "GTEx associations with a gene that is not protein-coding in GENCODE v50, or not in it, are "
            "counted, never scored",
            "the CRISPRi benchmark's held-out file is never evidence (R5); it is held out and scored like "
            "any source",
            "the AlphaGenome per-element response cache is not opened; the model is not called",
        ],
        "partitions": {
            "held_out_source": "each source in turn is the evaluation partition; the rest, masked by "
            "locus, is the only evidence a labelling may read",
            ms.TRAINING: "the CRISPRi benchmark's training file (K562): evidence when another source is "
            "held out",
            ms.HELDOUT: "the CRISPRi benchmark's held-out file: evaluation only, never evidence",
        },
    }


def main() -> None:
    t0 = time.time()
    units = ho.all_units()
    refs = ho.crispri_references()
    distance = ho.distance_labels()
    unchanged = ho.unchanged_labels()
    per_source: dict[str, Any] = {}
    scores: list[dict[str, Any]] = []
    comparisons: list[dict[str, Any]] = []
    for src in ho.sources(units):
        t = time.time()
        view = ho.evidence(src, units, refs)
        rest = ho.rest_labels(view)
        per_source[src] = {
            "family": ho.FAMILY[ho.assay_of(src)],
            "units": len(units[src]),
            "siblings_removed": [s for s in view.removed if s != src],
            "masked_in_view": view.masked,
            "heldout_file_pairs_never_evidence": view.heldout_file_pairs_dropped,
            "sources_in_view": sorted(view.sources),
            "relation_to_its_family": ho.relation_counts(src, units),
        }
        for endpoint in ho.ENDPOINTS[ho.assay_of(src)]:
            res = ho.evaluate(
                [rest, distance, unchanged], src, endpoint, units=units, references=refs, keep_draws=True
            )
            by = {r["labels"]: r for r in res}
            for a, b in COMPARISONS:
                comparisons.append(ho.difference(by[a], by[b]))
            for r in res:
                r.pop("_draws", None)
                scores.append(r)
            print(
                src,
                endpoint,
                [(r["labels"], r["status"], r.get("value"), r.get("interval")) for r in res],
                round(time.time() - t, 1),
                "s",
                flush=True,
            )
    params = {**ho.registration(), "comparisons": [f"{a} - {b}" for a, b in COMPARISONS]}
    payload = {
        "question": "held out one source at a time, how well do the remaining sources, distance to TSS "
        "and the unchanged labels predict its own endpoint",
        "status": ho.STATUS,
        "reuse": ho.REUSE,
        "may_be_called": ho.MAY_BE_CALLED,
        "may_not_be_called": list(ho.MAY_NOT_BE_CALLED),
        "registration": ho.registration(),
        "comparisons_registered": [f"{a} - {b}" for a, b in COMPARISONS],
        "sources": per_source,
        "gtex_excluded": ho.gtex_excluded(),
        "scores": scores,
        "comparisons": comparisons,
        "alphagenome_requests": 0,
        "per_element_response_cache_opened": False,
        "seconds": round(time.time() - t0, 1),
    }
    p = save_result(NAME, payload, manifest=manifest(ho.compiled_programs(), params))
    print("saved", p, flush=True)


if __name__ == "__main__":
    main()
