# SPDX-License-Identifier: AGPL-3.0-or-later
"""lane-hetblock, step 3: the registration, committed BEFORE any statistic of this lane is read.

    uv run python scripts/hetblock_register.py

It writes data/results/hetblock_registration.json: every BOUND field of
`genomeos/attribution/hetblock.py` imported from the module itself -- the question, the exchangeable
unit and why it is not the element-sharing component, the selection rule, the null that holds every
chromosome's gene count exactly, the statistic and the secondary statistic, the draw count and seed,
the gene floor, ALPHA with the stricter number reported, the degeneracy rule, the double exposure, the
decision table with BOTH OUTCOMES WRITTEN OUT, the falsifier, the stop rule and the limitations -- plus
the sha256 of every input the run will read.

It refuses to run unless data/results/hetblock_harness.json exists and reads HARNESS_REPRODUCES, so the
registration cannot be written before lane-scalehet's committed figures have been shown to come back;
and it refuses to overwrite a registration that already exists.

NO STATISTIC OF THIS LANE IS COMPUTED HERE. The only numbers that enter are the committed figures of
other lanes, the address-level counts the harness already committed, and sha256 digests.
"""

from __future__ import annotations

import json
import os
import sys
import time
from datetime import date
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
os.chdir(ROOT)

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import hetblock as hb  # noqa: E402
from genomeos.results import save_result  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from hetblock_harness import OWN_CODE, READERS, SOURCE_FILES  # noqa: E402

RESULT = "hetblock_registration"
OUT = "data/results/hetblock_registration.json"
HARNESS_RESULT = "data/results/hetblock_harness.json"


def refuse_unless_the_harness_reproduced() -> dict[str, Any]:
    """The harness must exist and read HARNESS_REPRODUCES: scalehet's figures come back first."""
    if Path(OUT).exists():
        raise SystemExit(f"refused: {OUT} already exists; a registration is written once")
    if not Path(HARNESS_RESULT).exists():
        raise SystemExit(f"refused: {HARNESS_RESULT} does not exist; run and commit the harness first")
    h = json.loads(Path(HARNESS_RESULT).read_text())
    if h["outcome"] != "HARNESS_REPRODUCES":
        raise SystemExit(f"refused: the committed harness reads {h['outcome']!r}: {h['problems']}")
    if not h["scalehet_reproduction"]["all_bit_identical"]:
        raise SystemExit("refused: the harness did not reproduce scalehet's figures bit-identically")
    if not h["singleton_block_equivalence"]["bit_identical"]:
        raise SystemExit("refused: the singleton-block null is not scalehet's null bit-identically")
    if h["degeneracy_check"]["degenerate"]:
        raise SystemExit(f"refused: the governing unit is degenerate: {h['degeneracy_check']['degenerate']}")
    return h


def manifest(entries: list[dict[str, Any]]) -> dict[str, Any]:
    return {
        "sources": [
            {
                "accession": "data/results/hetblock_harness.json (this lane, step 2)",
                "version": "pinned; it reproduced lane-scalehet's committed figures bit-identically",
            },
            {
                "accession": "data/results/scalehet_result.json (lane-scalehet)",
                "version": "pinned; the reading this lane re-tests",
            },
            {"accession": "data/results/multrna_result.json (7h, lane-multrna)", "version": "pinned"},
            {
                "accession": "the 24 compiled non-coding programs and the GENCODE gene table",
                "version": "pinned by sha256; addresses only",
            },
        ],
        "inputs": entries,
        "input_count": len(entries),
        "readers": READERS,
        "reader_count": len(READERS),
        "assembly": "GRCh38; genes as 7h matched them, by gene_id",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "registration_only": (
                "no statistic of this lane is computed here: the file fixes the words and the numbers"
                " the run is not allowed to differ from"
            ),
            "bound_fields": len(hb.BOUND),
            "block_span": hb.BLOCK_SPAN,
            "robust_span": hb.ROBUST_SPAN,
            "permutations": hb.PERMUTATIONS,
            "control_alpha": hb.CONTROL_ALPHA,
            "stricter_alpha": hb.STRICTER_ALPHA,
            "min_genes": hb.MIN_GENES,
            "min_blocks": hb.MIN_BLOCKS,
            "max_block_share": hb.MAX_BLOCK_SHARE,
        },
        "exclusions": ["none: the population is scalehet's, which is 7h's"],
        "partitions": {"none": "a registration partitions nothing"},
        "code_cleanliness": mf.code_cleanliness(__file__, OWN_CODE),
    }


def main() -> int:
    t0 = time.time()
    h = refuse_unless_the_harness_reproduced()
    pinned = [*SOURCE_FILES, HARNESS_RESULT]
    entries = [mf.input_entry(p, partition="pinned by sha256") for p in pinned]
    out: dict[str, Any] = {
        "result": RESULT,
        "date": date.today().isoformat(),
        "lane": hb.LANE,
        "step": "3 of 4: the registration, before any statistic of this lane is read",
        "harness": HARNESS_RESULT,
        "harness_outcome": h["outcome"],
        "scalehet_figures_reproduce": h["scalehet_reproduction"]["reading"],
        "governing_unit": h["governing_unit"]["governing_unit"],
        "governing_unit_reading": h["governing_unit"]["reading"],
        "element_sharing_graph": h["element_sharing_graph"],
        "locus_block_structure": h["locus_block_structure"],
        "singleton_block_equivalence": h["singleton_block_equivalence"],
        "both_outcomes_are_written_here": {
            "if_it_survives": hb.OUTCOMES["SURVIVES_BLOCK_NULL"],
            "if_it_does_not": hb.OUTCOMES["WITHDRAWN_NOT_ESTABLISHED"],
            "note": (
                "BOTH ARE WRITTEN BEFORE EITHER IS SEEN, in the words they will be reported in, and the"
                " second one WITHDRAWS lane-scalehet's rejection to NOT ESTABLISHED and says so plainly."
                " The withdrawal is the more valuable answer if that is what the data says, and it is the"
                " direction scalehet itself named as the one its null was biased against"
            ),
        },
        "pinned_sha256": {p: e["sha256"] for p, e in zip(pinned, entries, strict=True)},
        "readers": READERS,
        "alphagenome_requests": 0,
        "money": "none: no request is sent and no key is read",
        "nothing_of_this_lanes_statistic_is_in_this_file": (
            "no Q, no R, no p and no block null of this lane appears here. The only numbers are other"
            " lanes' committed figures, the address-level counts the harness committed, and sha256"
            " digests. The run refuses to start if any field below differs from the committed module or"
            " if any pinned sha256 has moved"
        ),
        **hb.registered_fields(),
        "seconds": round(time.time() - t0, 1),
    }
    print("saved", save_result(RESULT, out, manifest=manifest(entries)))
    print("registered", len(hb.BOUND), "bound fields; governing unit", out["governing_unit"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
