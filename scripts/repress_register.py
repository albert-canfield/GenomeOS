# SPDX-License-Identifier: AGPL-3.0-or-later
"""Register the repression lane before any count, table or direction is read
(data/results/repress2_registration.json).

    uv run --frozen python scripts/repress_register.py

It reads nothing but its own code. The mechanism, the prediction it makes, the refutation condition
that would overturn it, the blinded eligibility predicate, both floors and the three sets of reading
words all come from `genomeos.attribution.repress2`, so the registration and the code that applies it
cannot drift apart. No count is taken here, no measured pair is opened, no direction is read and no
model request is made.

The inputs the run will read are listed by their own paths with their digests, so a later run cannot
quietly read a different file. Each path is entered with `manifest.input_entry` one at a time;
`manifest.files_entry` is not used, because that path has a known defect reporting a group absent
without checking a file.
"""

from __future__ import annotations

import sys
from datetime import date
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import context_evidence as ce  # noqa: E402
from genomeos.attribution import crispri  # noqa: E402
from genomeos.attribution import measured as ms  # noqa: E402
from genomeos.attribution import repress2 as rp  # noqa: E402
from genomeos.results import save_result  # noqa: E402

RESULT = "repress2_registration"
#: The readings this lane's run will read, each named by its own result file and path.
DESCRIBES = (
    "data/results/not_open_profile.json",
    "data/results/context_evidence.json",
    "data/results/context_evidence_baserate.json",
)


def inputs() -> list[dict[str, Any]]:
    """The named files this registration can digest before any run, one entry per path.

    The compiled rules are re-enumerated from the per-chromosome results on disk, so the run digests
    those itself; here the three readings this lane describes, the track metadata behind the reader
    mapping and the two benchmark tables behind the measured layer are frozen by their own paths.
    """
    out: list[dict[str, Any]] = []
    for name in DESCRIBES:
        p = Path(name)
        if p.exists():
            out.append(mf.input_entry(p, partition="the reading this lane describes"))
    if ce.TRACK_METADATA.exists():
        out.append(mf.input_entry(ce.TRACK_METADATA, partition="the reader mapping table"))
    for name in ms.CRISPRI_FILES:
        p = crispri.KNOWLEDGE / name
        if p.exists():
            out.append(mf.input_entry(p, partition=ms.CRISPRI_SPLIT_OF[name]))
    return out


def main() -> None:
    payload = rp.registration()
    payload["date"] = date.today().isoformat()
    payload["describes"] = list(DESCRIBES)
    payload["what_the_run_will_write"] = "data/results/repress2_population.json"
    payload["order"] = [
        "1. gate 1: the blinded eligibility count, before any direction is read",
        "2. if below a floor, record the no-go with the count and continue only as the descriptive part",
        "3. the registered prediction, tested as the repression share by effect band and reader state",
        "4. the 0-request internal check: sign concordance across the four retained cells",
    ]
    entries = inputs()
    payload["frozen_inputs"] = [{k: e[k] for k in ("path", "sha256", "bytes") if k in e} for e in entries]
    payload["result_manifest"] = {
        "sources": [
            {
                "accession": "the compiled non-coding programs of the 24 chromosomes, as "
                "genomeos.attribution.compile emits them",
                "version": "the attributed element results on this machine, digested in frozen_inputs",
            },
            {
                "accession": "ENCODE CRISPRi enhancer-gene benchmark (EngreitzLab/CRISPR_comparison, "
                "Gschwind et al.), as genomeos.attribution.measured already caches it",
                "version": "the assay cache on this machine; the run digests the layer it reads",
            },
        ],
        "inputs": entries,
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "independent_locus_span": rp.LOCUS_SPAN,
            "locus_floor": rp.LOCUS_FLOOR,
            "link_floor": rp.POSITIVE_FLOOR,
            "strong_effect": rp.nop.STRONG_EFFECT,
            "retained_cells": list(rp.RETAINED_CELLS),
            "aggregate_function": "count",
            "fill_value": 0,
        },
        "exclusions": [
            "no rule is excluded: every rule the compiler emits is enumerated and the ones outside a "
            "reported population are counted in their own row",
            "a rule whose activity axis holds more than one term is in neither the repression nor the "
            "activation population and is counted separately",
        ],
        "partitions": {
            "gate_1_cell_matched": "the measured pair's cell equals the rule's own `when: cell_type`",
            "gate_1_cell_ignored": "element and gene only, reported beside it and never pooled with it",
            "descriptive": "the assessable rules of the predicted layer, by effect band and reader state",
        },
    }
    path = save_result(RESULT, payload)
    print(f"{RESULT}: {path}")
    print(f"  floors: {rp.LOCUS_FLOOR} independent loci, {rp.POSITIVE_FLOOR} eligible measured links")
    print(f"  inputs digested: {len(entries)}")
    print("  nothing counted, no direction read, 0 model requests")


if __name__ == "__main__":
    main()
