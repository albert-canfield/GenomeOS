# SPDX-License-Identifier: AGPL-3.0-or-later
"""The genome-wide census of what `compile.py` emits per element now and at the bar.

    uv run python scripts/emit7g_census.py

ROADMAP §5 item 7g's next step: "the compiler emits every coding target at its bar (a `compile.py`
change, registered against the compiled programs' census tests)". This script reads the 24 compiled
programs as they stand (data/knowledge/compiled, the files M1 ran and every census of this series
reads) and the per-element response cache, and writes data/results/emit7g_census.json: the rules
emitted now, the rules the emission would write, the targets per element before and after, and how many
elements gain a target. It is run and committed BEFORE anything in `compile.py` moves.

Nothing is fitted, nothing is tuned, no label is read: the bar is the compiler's own object
(`predict.enhancer_target.MIN_EFFECT`) and the candidate list is `targetchoice.candidate_list`, the list
registered at b746a23. No request is sent, no key is read. 0 model requests, no network, no money.
"""

from __future__ import annotations

import hashlib
import json
import sys
import time
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import crispri  # noqa: E402
from genomeos.attribution import emit7g as e7  # noqa: E402
from genomeos.attribution.targets import ElementResponses  # noqa: E402
from genomeos.results import save_result  # noqa: E402

RESULT = "emit7g_census"
OWN_CODE = (
    "genomeos/attribution/emit7g.py",
    "scripts/emit7g_census.py",
    "tests/test_emit7g.py",
)


def code_sha256(paths: tuple[str, ...]) -> dict[str, str]:
    return {p: hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in paths}


def gencode_path(chrom: str) -> Path:
    from genomeos.genome.annotation import default_gencode

    path = default_gencode({chrom})
    if path is None:
        raise SystemExit(f"refused: no GENCODE annotation for {chrom} on this machine")
    return Path(path)


def coding_symbols(chrom: str) -> set[str]:
    """The chromosome's protein-coding symbols in GENCODE v50 (targetchoice.CANDIDATES' own set)."""
    from genomeos.genome import Annotation

    ann = Annotation.from_gff3(gencode_path(chrom), {chrom})
    return {g.symbol for g in ann.protein_coding() if g.locus.chrom == chrom}


def program_path(chrom: str) -> Path:
    return e7.COMPILED_DIR / f"noncoding_{chrom}.bio"


def chromosomes() -> list[str]:
    """Every chromosome of `emit7g.CHROMS` whose compiled program is on this machine, or a refusal.

    The whole store or none: a census of a partial store would report a smaller emission as if it were
    the whole genome's, which is exactly the defect the glob census of 2026-10-11 was opened on.
    """
    missing = [c for c in e7.CHROMS if not program_path(c).exists()]
    if missing:
        raise SystemExit(f"refused: {len(missing)} compiled programs are not on this machine: {missing}")
    return list(e7.CHROMS)


def inputs(chroms: list[str]) -> list[dict[str, Any]]:
    return [
        mf.input_entry(
            crispri.ELEMENT_CACHE,
            partition=None,
            role="ElementResponses: every gene's predicted deletion response at an element",
        ),
        *[mf.input_entry(program_path(c), partition=None, role="compiled program") for c in chroms],
        mf.files_entry(
            "GENCODE v50 per-chromosome GFF3 (protein-coding symbols, targetchoice.CANDIDATES)",
            [str(gencode_path(c)) for c in chroms],
            partition=None,
        ),
    ]


def manifest(entries: list[dict[str, Any]]) -> dict[str, Any]:
    return {
        "sources": [
            {
                "accession": "ENCODE SCREEN cCREs scored by AlphaGenome deletion (per-element response"
                " cache)",
                "version": "the 2026-09 all-element sweep; pinned here by sha256",
            },
            {
                "accession": "the compiled non-coding programs (data/knowledge/compiled)",
                "version": "untracked, as built in this checkout; pinned by sha256",
            },
            {"accession": "GENCODE v50 basic annotation, per chromosome", "version": "pinned by sha256"},
        ],
        "inputs": entries,
        "input_count": len(entries),
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "min_effect": e7.BAR,
            "bar": e7.BAR_NOTE,
            "candidate_list": "targetchoice.candidate_list, the list registered at b746a23",
        },
        "exclusions": [
            "experimental rule lines are counted and never replaced: what an assay measured is not a"
            " model's attribution (signaltrace.PREDICTED_ONLY names the same lines)",
            "no label of any screen is read: this census has no outcome column at all",
        ],
        "partitions": {"none": "a census of the compiled store; there is no split and no fit"},
        "code_cleanliness": mf.code_cleanliness(__file__, OWN_CODE),
    }


def main() -> int:
    t0 = time.time()
    chroms = chromosomes()
    entries = inputs(chroms)
    responses = ElementResponses(crispri.ELEMENT_CACHE)
    per_chrom: dict[str, Any] = {}
    for chrom in chroms:
        per_chrom[chrom] = e7.chromosome_census(
            chrom, program_path(chrom).read_text(), coding_symbols(chrom), responses
        )
        row = per_chrom[chrom]
        print(
            chrom,
            row["rules_predicted_now"],
            "->",
            row["rules_predicted_after"],
            f"(+{row['rules_added']}), {row['elements_gaining_a_target']} elements gain",
            f"{time.time() - t0:.0f} s",
            flush=True,
        )
    total = e7.combine(per_chrom)
    out = {
        "lane": e7.LANE,
        "item": e7.ITEM,
        "population": e7.POPULATION,
        "bar": e7.BAR_NOTE,
        "min_effect": e7.BAR,
        "emission_now": e7.EMISSION_NOW,
        "emission_at_the_bar": e7.EMISSION_AFTER,
        "rule_fields": list(e7.RULE_FIELDS),
        "genome_wide": total,
        "per_chromosome": per_chrom,
        "frozen_code": {"sha256": code_sha256(OWN_CODE)},
        "alphagenome_requests": 0,
        "money": "none: no request is sent and no key is read",
        "network": "none",
        "seconds": round(time.time() - t0, 1),
    }
    print("saved", save_result(RESULT, out, manifest=manifest(entries)))
    print(json.dumps(total, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
