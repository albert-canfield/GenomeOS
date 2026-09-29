# SPDX-License-Identifier: AGPL-3.0-or-later
"""Item 12 S4, applied once: the unchanged compiled labels judged claim by claim against every held-out
source, with target accuracy, role accuracy and coverage reported apart, per axis, with denominators.

    uv run python scripts/s4_correctness_run.py

Registered 2026-09-29 before it was run (genomeos/attribution/correctness.py, docs/ATTRIBUTION.md "S4
registered"). An internal development benchmark only. No model request; the per-element response cache
is never opened. Writes data/results/attribution_correctness.json.

`--compiled DIR` and `--chroms` exist so the pipeline can be exercised on a synthetic program; the
registered run uses neither.
"""

from __future__ import annotations

import argparse
import time
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

from genomeos import manifest as mf
from genomeos.attribution import correctness as co
from genomeos.attribution import holdout as ho
from genomeos.attribution import measured as ms
from genomeos.results import save_result

NAME = "attribution_correctness"


@mf.depends_on_models("alphagenome")  # the unchanged labels are the compiled predicted layer
def manifest(programs: list[Path], parameters: dict[str, Any]) -> dict[str, Any]:
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
            "the measured twins (`*_measured` blocks and experimental rules) are the evidence and are not "
            "read as labels",
            "`unknown`, `unassigned` and unchosen `|` alternatives are counted as not claims, never judged",
            "CRISPRi pairs the benchmark marks as not a valid connection are not read",
            "observation kinds holdout does not load (contact, conservation, allelic readouts, the model, "
            "annotation and registry) are in the table and judge nothing in this run",
            "the AlphaGenome per-element response cache is not opened; the model is not called",
        ],
        "partitions": {
            "held_out_source": "every holdout source judges the claims; the unchanged labels read none of "
            "them, and holdout.check_provenance is called for each",
            ms.TRAINING: "the CRISPRi benchmark's training file: judges like any source",
            ms.HELDOUT: "the CRISPRi benchmark's held-out file: judges (evaluation), never evidence for a "
            "labelling; verdicts decided with one of its pairs are counted per axis",
        },
    }


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--compiled", type=Path, default=ho.COMPILED_DIR)
    ap.add_argument("--chroms", nargs="*", default=None)
    ap.add_argument("--dry", action="store_true", help="print the quantities, write nothing")
    a = ap.parse_args(argv)
    t0 = time.time()
    units = ho.all_units()
    refs = ho.crispri_references()
    labels = co.unchanged_labels()
    nc: dict[str, Counter] = defaultdict(Counter)
    claims = co.compiled_claims(a.compiled, a.chroms, nc)
    report = co.judge(claims, labels, units=units, references=refs, not_claims=nc)
    body = report.to_dict()
    for q, block in body["quantities"].items():
        for axis, s in block.items():
            print(q, axis, f"{s['numerator']:,} / {s['denominator']:,}", s["beside"], flush=True)
    if a.dry:
        return
    params = {
        "registration": co.registration(),
        "holdout": ho.registration(),
        "compiled": str(a.compiled),
        "chroms": a.chroms or "all",
    }
    payload = {
        "question": "of the claims the unchanged compiled labels state, on each of five axes, how many can "
        "any held-out observation judge, and of those how many does it establish or refute",
        "status": co.STATUS,
        "reuse": co.REUSE,
        "may_be_called": co.MAY_BE_CALLED,
        "may_not_be_called": list(co.MAY_NOT_BE_CALLED),
        "expected_before_the_run": co.EXPECTED,
        "registration": co.registration(),
        "units_per_source": {s: len(units[s]) for s in ho.sources(units)},
        **body,
        "alphagenome_requests": 0,
        "per_element_response_cache_opened": False,
        "seconds": round(time.time() - t0, 1),
    }
    programs = ho.compiled_programs(a.compiled)
    if a.chroms:
        programs = [p for p in programs if p.stem.removeprefix("noncoding_") in set(a.chroms)]
    p = save_result(NAME, payload, manifest=manifest(programs, params))
    print("saved", p, round(time.time() - t0, 1), "s", flush=True)


if __name__ == "__main__":
    main()
