# SPDX-License-Identifier: AGPL-3.0-or-later
"""The training census of the multi-target compiled choice, run before the registration.

    uv run python scripts/targetchoice_census.py

Reads the TRAINING split only (crispri.TRAINING; the covered pairs are the fit rows) and nothing of the
held-out file: how many compiled elements hold 2 or more candidate targets at the compiler's bar
(targetchoice.CANDIDATES), where the measured training gene sits among them, and the k of the secondary
top-k choice by targetchoice.CENSUS_RULE. Writes data/results/targetchoice_census.json. The k it names
is then written into targetchoice.TOP_K, committed, and only after that is the registration made.

No request is sent, no key is read. 0 model requests, no network, no money.
"""

from __future__ import annotations

import hashlib
import json
import sys
import time
from collections import defaultdict
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import crispri  # noqa: E402
from genomeos.attribution import targetchoice as tc  # noqa: E402
from genomeos.attribution.measured import CRISPRI_SPLIT_OF  # noqa: E402
from genomeos.attribution.targets import ElementResponses  # noqa: E402
from genomeos.evidence import COMPILED_DIR  # noqa: E402
from genomeos.results import save_result  # noqa: E402

RESULT = "targetchoice_census"
CHROMS = [f"chr{c}" for c in [*range(1, 23), "X", "Y"]]
OWN_CODE = (
    "genomeos/attribution/targetchoice.py",
    "scripts/targetchoice_census.py",
    "tests/test_targetchoice.py",
)


def chromosomes(pairs: list[crispri.Pair]) -> list[str]:
    seen = {p.chrom for p in pairs}
    return [c for c in CHROMS if c in seen]


def code_sha256(paths: tuple[str, ...]) -> dict[str, str]:
    return {p: hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in paths}


def gencode_path(chrom: str) -> Path:
    from genomeos.genome.annotation import default_gencode

    path = default_gencode({chrom})
    if path is None:
        raise SystemExit(f"refused: no GENCODE annotation for {chrom} on this machine")
    return Path(path)


def coding_symbols(chrom: str) -> set[str]:
    """The chromosome's protein-coding symbols in GENCODE v50 (targetchoice.CANDIDATES)."""
    from genomeos.genome import Annotation

    ann = Annotation.from_gff3(gencode_path(chrom), {chrom})
    return {g.symbol for g in ann.protein_coding() if g.locus.chrom == chrom}


def program_path(chrom: str) -> Path:
    return Path(COMPILED_DIR) / f"noncoding_{chrom}.bio"


def lookups_for(chroms: list[str], log: Any = print) -> dict[str, tc.ChoiceLookup]:
    """Every chromosome's compiled elements with their candidate targets, one chromosome at a time."""
    t0 = time.time()
    responses = ElementResponses(crispri.ELEMENT_CACHE)
    out: dict[str, tc.ChoiceLookup] = {}
    for chrom in chroms:
        out[chrom] = tc.chromosome_choice(
            chrom, program_path(chrom).read_text(), coding_symbols(chrom), responses
        )
        log(chrom, dict(out[chrom].counts), f"{time.time() - t0:.0f} s", flush=True)
    return out


def annotate(pairs: list[crispri.Pair]) -> None:
    """The anchor's columns (top_target, deletion_drop) and coverage, by the anchor's own code."""
    crispri.annotate(pairs, crispri.DeletionTable(), crispri.ElementCache())


def shared_inputs(chroms: list[str]) -> list[dict[str, Any]]:
    """The inputs the census and the run both read: the deletion tables, the compiled programs, GENCODE."""
    return [
        mf.input_entry(crispri.ELEMENTS, partition=None, role="DeletionTable: elements by overlap (anchor)"),
        mf.input_entry(
            crispri.ELEMENT_CACHE,
            partition=None,
            role="ElementCache and ElementResponses: every gene's deletion response at an element",
        ),
        *[mf.input_entry(program_path(c), partition=None, role="compiled program") for c in chroms],
        mf.files_entry(
            "GENCODE v50 per-chromosome GFF3 (protein-coding symbols, targetchoice.CANDIDATES)",
            [str(gencode_path(c)) for c in chroms],
            partition=None,
        ),
    ]


def inputs(chroms: list[str], training: list[crispri.Pair]) -> list[dict[str, Any]]:
    return [
        mf.input_entry(
            crispri.KNOWLEDGE / crispri.TRAINING,
            partition=CRISPRI_SPLIT_OF[crispri.TRAINING],
            pairs=len(training),
            role="the training split: its covered pairs are the fit rows; the census reads nothing else",
        ),
        *shared_inputs(chroms),
    ]


def manifest(entries: list[dict[str, Any]]) -> dict[str, Any]:
    return {
        "sources": [
            {
                "accession": "EngreitzLab/CRISPR_comparison EPCrisprBenchmark training_K562"
                " (Gschwind et al.)",
                "version": "fetched 2026-09-16; pinned here by sha256",
                "url": crispri.BASE_URL,
            },
            {
                "accession": "ENCODE SCREEN cCREs scored by AlphaGenome deletion (all-element table and"
                " per-element response cache)",
                "version": "the 2026-09 all-element sweep; pinned here by sha256",
            },
            {
                "accession": "the compiled non-coding programs (data/knowledge/compiled)",
                "version": "pinned by sha256",
            },
            {"accession": "GENCODE v50 basic annotation, per chromosome", "version": "pinned by sha256"},
        ],
        "inputs": entries,
        "input_count": len(entries),
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "min_effect": tc.MIN_EFFECT,
            "k_tried": list(tc.K_TRIED),
            "census_rule": tc.CENSUS_RULE,
            "fit": "crispri.logistic_fit, lam 1e-3, 25 rounds, on the covered K562 training pairs",
        },
        "exclusions": [
            "the held-out file is not read: the census and the choice of k use the training split only",
            "experimental rule lines of the compiled programs (signaltrace.PREDICTED_ONLY)",
        ],
        "partitions": {
            "training": "the fit rows: the covered K562 training pairs; every reading here is in-sample",
        },
        "code_cleanliness": mf.code_cleanliness(__file__, OWN_CODE),
    }


def main() -> int:
    t0 = time.time()
    training = crispri.load(crispri.TRAINING)
    chroms = chromosomes(training)
    entries = inputs(chroms, training)
    annotate(training)
    fit_rows = [p for p in training if p.covered]
    lookups = lookups_for(chroms)
    body = tc.census(lookups, fit_rows)
    per_chrom: dict[str, Any] = defaultdict(dict)
    for c, cl in lookups.items():
        per_chrom[c] = dict(cl.counts)
    out = {
        "lane": tc.LANE,
        "question": "the training census before the registration: " + tc.QUESTION,
        "training_split": tc.TRAINING_SPLIT,
        "candidates": tc.CANDIDATES,
        "min_effect": tc.MIN_EFFECT,
        **body,
        "per_chromosome": dict(per_chrom),
        "frozen_code": {"sha256": code_sha256(OWN_CODE[:2])},
        "heldout_read": False,
        "alphagenome_requests": 0,
        "money": "none: no request is sent and no key is read",
        "seconds": round(time.time() - t0, 1),
    }
    print("saved", save_result(RESULT, out, manifest=manifest(entries)))
    summary = {k: out[k] for k in ("top_k", "secondary")}
    summary["reached"] = out["compiled_elements_reached_by_a_fit_row"]
    summary["one"] = {k: out["one_target_in_sample"][k] for k in ("kept", "kept_regulated")}
    summary["by_k"] = {
        k: {x: v[x] for x in ("kept", "kept_regulated", "average_precision_in_sample")}
        for k, v in out["by_k_in_sample"].items()
    }
    print(json.dumps(summary, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
