# SPDX-License-Identifier: AGPL-3.0-or-later
"""Build increment 2 of the response map and write it to the result registry.

    uv run --frozen python scripts/response_map_increment2.py
    uv run --frozen python scripts/response_map_increment2.py --chroms chr21 \
        --result response_map_increment2_chr21

The population was counted before anything was built, by `scripts/response_map_coverage.py`
(`data/results/response_map_coverage.json`, lane-rmap2): 110 of the 473 independent loci of the
measured layer's perturbation assay carry a measured CRISPRi perturbation, a compiled rule for the
element it attaches to and a reader state. This script arranges those loci - it selects nothing, fits
nothing, scores nothing, computes no verdict, calls no model and downloads nothing. Every input was
already on disk.

`genomeos.response_map2.build` refuses rather than write a payload whose built loci plus its named
exclusions do not come back to that count, and `response_map.check` - the same gate increment 1
passes - refuses an assertion that cannot name its file, its record key and that file's sha256, or a
model output labelled `observed`.
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos import response_map2 as r2  # noqa: E402
from genomeos.attribution import cell2  # noqa: E402
from genomeos.attribution import context_evidence as ce  # noqa: E402
from genomeos.attribution import measured as ms  # noqa: E402
from genomeos.genome import reader  # noqa: E402
from genomeos.results import RESULTS_DIR, save_result  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]

#: This lane's own files. Everything else uncommitted in this shared checkout belongs to another lane.
OWN_CODE = (
    "genomeos/response_map2.py",
    "scripts/response_map_increment2.py",
    "tests/test_response_map2.py",
)

#: The per-chromosome inputs the compiler and the measured layer read, as globs, so the manifest names
#: them without a hand-written list of several hundred paths going stale. One entry per file, never a
#: group: `scripts/manifest_rebuild.py` resolves an input by its `path` and a grouped entry carries a
#: label there, so a grouped input would read as absent without any file being checked.
INPUT_GLOBS = (
    "enhancer_targets_all_chr*.json",
    "constrained_targets_chr*.json",
    "enhancer_targets_chr*.json",
)


def inputs_of(payload: dict, chroms: list[str]) -> list[dict]:
    """Every file this run read, by its own path and sha256.

    The files the assertions name are already hashed on the payload, so they are recorded from there
    rather than hashed twice; the rest are the per-chromosome tables the compiler read to enumerate
    the elements.
    """
    out: list[dict] = []
    seen: set[str] = set()
    for s in payload["sources"]:
        if not s["available"]:
            continue
        p = ROOT / s["path"]
        if not p.is_file() or s["path"] in seen:
            continue
        seen.add(s["path"])
        out.append(mf.input_entry(s["path"], partition=None, sha256_on_the_map=s["sha256"]))
    paths: list[Path] = [ce.TRACK_METADATA]
    for f in ms.CRISPRI_FILES:
        paths.append(ms.CRISPRI_KNOWLEDGE / f)
    for cell in ce.READER_TERMS:
        paths += [RESULTS_DIR / reader.peaks_path(cell, c).name for c in chroms]
    for glob in INPUT_GLOBS:
        for c in chroms:
            paths += sorted(RESULTS_DIR.glob(glob.replace("chr*", c)))
    for p in paths:
        rel = Path(p).as_posix()
        if rel in seen or not (ROOT / rel).is_file():
            continue
        seen.add(rel)
        out.append(mf.input_entry(rel, partition=ms.CRISPRI_SPLIT_OF.get(Path(rel).name)))
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--chroms", default=",".join(r2.CHROMS))
    ap.add_argument(
        "--result",
        default=r2.RESULT,
        help=(
            "the result name to write. A run over fewer than all 24 chromosomes writes under a name "
            "of its own, so a partial build can never overwrite the genome-wide one"
        ),
    )
    ap.add_argument("--summary", action="store_true", help="print the reading and write nothing")
    args = ap.parse_args()
    chroms = [c for c in args.chroms.split(",") if c]
    if args.result == r2.RESULT and tuple(chroms) != r2.CHROMS:
        raise SystemExit(
            f"{r2.RESULT} is the genome-wide build; a run over {len(chroms)} chromosomes needs --result"
        )

    started = time.time()
    payload = r2.build(ROOT, chroms)
    print(r2.summary(payload))
    if args.summary:
        return

    payload["seconds"] = round(time.time() - started, 1)
    payload["alphagenome_requests"] = 0
    payload["money"] = "none: every input was already on disk"
    payload["code_cleanliness"] = ce.code_cleanliness(Path(__file__).resolve(), OWN_CODE)
    payload["result_manifest"] = {
        "sources": [
            {
                "accession": "EngreitzLab/CRISPR_comparison EPCrisprBenchmark (Gschwind et al. 2025)",
                "version": "the two local benchmark tables under data/knowledge/crispri",
            },
            {
                "accession": "ENCODE DNase-seq narrowPeak, GRCh38, released, 13 biosamples",
                "version": "as reader v1 cached them under data/results/dnase_*_chr*.bed.gz",
            },
            {"accession": "AlphaGenome output track metadata", "version": str(ce.TRACK_METADATA)},
            {
                "accession": "the attributed elements of the 24 chromosomes",
                "version": (
                    "the three deletion runs' own tables, read through genomeos.attribution.targets "
                    "in this process, each element's assertion naming the file it came from"
                ),
            },
            {
                "accession": "response_map_coverage",
                "version": (
                    f"{r2.COVERAGE}, sha256 {r2.COVERAGE_SHA256}: the count of the population, taken "
                    "before anything was built"
                ),
            },
        ],
        "inputs": inputs_of(payload, chroms),
        "inputs_are_recorded_one_file_per_entry": (
            "scripts/manifest_rebuild.py resolves an input by its `path`, and a grouped entry carries "
            "a label there, so a grouped input would read as absent without any file being checked"
        ),
        "assembly": payload["conventions"]["assembly"],
        "coordinates": payload["conventions"]["coordinates"],
        "parameters": {
            "independent_locus_span": cell2.INDEPENDENT_LOCUS_SPAN,
            "reciprocal_overlap": ms.RECIPROCAL_OVERLAP,
            "openness_call": ce.OPENNESS_CALL,
            "states": list(ce.STATES),
            "statuses": list(r2.STATUSES),
            "page_size": r2.PAGE,
            "chromosomes": chroms,
        },
        "exclusions": [
            "the population is the loci the count already fixed; this run selects none of its own, "
            "and every addable locus is either built or named under an exclusion cause in "
            "`reconciliation.loci_excluded`",
            "the measured layer's other three assays are not perturbations of the native locus and no "
            "locus of theirs is here, which is stated in `coverage.what_it_does_not_cover` rather "
            "than left to be inferred",
            "no new cut-off is introduced: the overlap rule, the locus span and the openness call are "
            "all imported from where they were registered",
        ],
        "partitions": {
            "loci": "one entry per independent locus under the imported 1 Mb grouping",
            "by_status": "observed, predicted, inferred, unknown over the assertions",
            "reconciliation": "the built loci against the count, with every exclusion by cause",
        },
        "code_cleanliness": payload["code_cleanliness"],
    }
    path = save_result(args.result, payload, compact=True)
    print(f"{args.result}: {path}")


if __name__ == "__main__":
    main()
