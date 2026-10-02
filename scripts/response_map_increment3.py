# SPDX-License-Identifier: AGPL-3.0-or-later
"""Build increment 3 of the response map and save it, with its manifest.

    uv run --frozen python scripts/response_map_increment3.py
    uv run --frozen python scripts/response_map_increment3.py --summary

The population was counted before this was built, at `data/results/response_map_increment3_count.json`
(committed first, on purpose), and `genomeos.response_map3.build` refuses to run over a set of a
different size or to return a payload whose loci do not reconcile against it. 0 model requests, no
downloads, no money: every input was already on disk. Every open of the per-element response cache is
recorded by patching `builtins.open` for the whole run, because `targets.run_elements` reads its path
out of a result file's field and no reading of the source can see it.
"""

from __future__ import annotations

import argparse
import builtins
import sys
import time
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos import response_map3 as rm3  # noqa: E402
from genomeos.attribution import cell2, eqtl, mpra  # noqa: E402
from genomeos.attribution import context_evidence as ce  # noqa: E402
from genomeos.attribution import measured as ms  # noqa: E402
from genomeos.attribution import targets as tg  # noqa: E402
from genomeos.genome import reader  # noqa: E402
from genomeos.results import RESULTS_DIR, save_result  # noqa: E402

RESULT = rm3.RESULT
ROOT = Path(__file__).resolve().parents[1]
ENTRY = "scripts/response_map_increment3.py"
OWN_CODE = (
    "genomeos/response_map3.py",
    "scripts/response_map_increment3.py",
    "tests/test_response_map3.py",
)
INPUT_GLOBS = (
    "budget_chr*.json",
    "budget_axes_chr*.json",
    "variation_chr*.json",
    "duplication_chr*.json",
    "domains_chr*.json",
    "unknown_chr*.json",
    "reader_*_chr*.json",
    "ccres_chr*.bed.gz",
    "enhancer_targets_chr*.json",
    "enhancer_targets_all_chr*.json",
    "constrained_targets_chr*.json",
)

#: Files this build opens that no glob of `data/results` reaches, each with why. An input is what the
#: writer opened, not what it said it read: the audit hook in `genomeos.manifest` records every open
#: under `data/` and `save_result` refuses a result whose reads exceed its declared inputs. That is
#: also what keeps this list from going stale, so it is written out here rather than shared with the
#: count's own list: if either drifts, the result is refused rather than published.
OPENED_BESIDE_THE_RESULTS = {
    # `attribution.targets` reads a run's elements out of the table its result's `elements_where` field
    # points at, so the path comes out of a result file at run time and no reading of the source names
    # it. About 1.4 GB over 26 chromosomes.
    "the per-chromosome element tables reached through a result file's `elements_where` pointer": [
        f"data/knowledge/alphagenome/all_elements/{c}.json" for c in rm3.CHROMS
    ],
    # the measured layer loads all four assays per chromosome. The lentiMPRA files are declared below
    # as sources of this map's own reporter assertions; VISTA and saturation mutagenesis contribute no
    # assertion here and are opened all the same.
    "the measured layer's other assay sources, opened by `compile._measured_rows`": [
        "data/knowledge/vista/locus.tsv.gz",
        "data/knowledge/satmut/elements.tsv.gz",
    ],
    # increment 1's map is read whole by `response_map2.map1_keys`, which is how this build knows which
    # locus increment 1 already covered; everything that view reads is opened with it
    "increment 1's own view, read by `response_map2.map1_keys`": [
        "data/results/attribution_correctness_v3.json",
        "data/results/crispri_direction.json",
        "data/results/discovery_review.json",
        "data/results/loci_benchmark.json",
        "data/knowledge/ReactomePathways.txt",
        "data/knowledge/compiled/noncoding_chr11.bio",
        "data/knowledge/hic_contact/K562_4DNFITUOMFUQ_5000.json",
        "data/organisms/human/erythrocyte.bio",
        "data/cache/rates/schofield2018_TableS2_halflives.xlsx",
    ],
}


class CacheAudit:
    """Every open of a file under the per-element response cache, recorded rather than asserted."""

    def __init__(self, root: Path = tg.ELEMENT_CACHE) -> None:
        self.root = Path(root).resolve()
        self.opens: list[str] = []
        self._real = builtins.open

    def __enter__(self) -> CacheAudit:
        audit = self

        def patched(file, *a, **kw):  # type: ignore[no-untyped-def]
            try:
                p = Path(file).resolve()
                if p == audit.root or audit.root in p.parents:
                    audit.opens.append(p.as_posix())
            except (TypeError, ValueError, OSError):
                pass
            return audit._real(file, *a, **kw)

        builtins.open = patched  # type: ignore[assignment]
        return self

    def __exit__(self, *exc: object) -> None:
        builtins.open = self._real  # type: ignore[assignment]

    def block(self) -> dict[str, Any]:
        root = self.root
        return {
            "per_element_response_cache_root": root.relative_to(ROOT).as_posix()
            if root.is_relative_to(ROOT)
            else root.as_posix(),
            "opens": len(self.opens),
            "files": sorted(set(self.opens))[:50],
            "how_this_is_known": (
                "recorded, not inferred: `builtins.open` was patched for the whole run and every open "
                "of a path under the cache root was appended. No reading of the source could see it, "
                "because `targets.run_elements` reads its path out of a result file's field"
            ),
        }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--chroms", default=",".join(rm3.CHROMS))
    ap.add_argument("--result", default=RESULT)
    ap.add_argument("--summary", action="store_true", help="print the summary and save nothing")
    args = ap.parse_args()
    chroms = [c for c in args.chroms.split(",") if c]
    name = args.result
    if name == RESULT and tuple(chroms) != tuple(rm3.CHROMS):
        raise SystemExit(
            f"{RESULT} is the genome-wide map; a run over {len(chroms)} chromosomes needs --result"
        )

    started = time.time()
    with CacheAudit() as audit:
        payload = rm3.build(ROOT, chroms)
    payload["seconds"] = round(time.time() - started, 1)
    payload["per_element_response_cache"] = audit.block()
    payload["alphagenome_requests"] = 0
    payload["money"] = "none: every input was already on disk"
    payload["code_cleanliness"] = mf.code_cleanliness(ENTRY, OWN_CODE, ROOT)

    if args.summary:
        print(rm3.summary(payload))
        return

    inputs = [mf.input_entry(ce.TRACK_METADATA, partition=None)]
    for f in ms.CRISPRI_FILES:
        p = ms.CRISPRI_KNOWLEDGE / f
        if p.exists():
            inputs.append(mf.input_entry(p, partition=ms.CRISPRI_SPLIT_OF[f]))
    seen = {str(i["path"]) for i in inputs}
    paths = [Path(rm3.COUNT_RESULT)]
    paths += [
        p
        for cell in ce.READER_TERMS
        for chrom in chroms
        if (p := RESULTS_DIR / reader.peaks_path(cell, chrom).name).exists()
    ]
    for g in INPUT_GLOBS:
        for chrom in chroms:
            paths += sorted(RESULTS_DIR.glob(g.replace("chr*", chrom)))
    for acc in mpra.FILES.values():
        p = mpra.KNOWLEDGE / f"{acc}.bed.gz"
        if p.exists():
            paths.append(p)
    paths += sorted(eqtl.KNOWLEDGE.glob("hits_*.tsv"))
    p = eqtl.KNOWLEDGE / "distil_summary.json"
    if p.exists():
        paths.append(p)
    for group in OPENED_BESIDE_THE_RESULTS.values():
        # repository-relative, so `manifest_rebuild.py` can check the bytes its own
        # worktree reads rather than an absolute path outside the linked stores
        paths += [Path(x) for x in group]
    for q in paths:
        if q.as_posix() not in seen:
            seen.add(q.as_posix())
            inputs.append(mf.input_entry(q, partition=None))

    payload["result_manifest"] = {
        "sources": [
            {
                "accession": "EngreitzLab/CRISPR_comparison (Gschwind et al. 2025)",
                "version": "the two local benchmark tables under data/knowledge/crispri",
            },
            {
                "accession": "ENCODE DNase-seq narrowPeak, GRCh38, released, 13 biosamples",
                "version": "as reader v1 cached them under data/results/dnase_*_chr*.bed.gz",
            },
            {
                "accession": "ENCODE4 lentiMPRA, joint library ENCSR106SZM",
                "version": ", ".join(f"{c} {a}" for c, a in sorted(mpra.FILES.items())),
            },
            {
                "accession": "GTEx v8 single-tissue cis-eQTLs, significant variant-gene pairs",
                "version": (
                    "the local distillation under data/knowledge/gtex, kept only inside the three "
                    "sampled element sets; see data/knowledge/gtex/distil_summary.json"
                ),
            },
            {
                "accession": "the compiled non-coding programs of the chromosomes",
                "version": (
                    "the compiler's own element and measured-layer enumeration, run in this process "
                    "from the results on disk, not read back from a .bio file"
                ),
            },
        ],
        "inputs": inputs,
        "inputs_opened_beside_the_declared_results": {
            "why_they_are_listed": (
                "an input is what the writer opened, not what it said it read. The audit hook in "
                "`genomeos.manifest` records every open under `data/` and `save_result` refuses a "
                "result whose reads exceed its declared inputs, which is how these came to be named"
            ),
            "groups": {k: len(v) for k, v in OPENED_BESIDE_THE_RESULTS.items()},
            "the_pointer_case": (
                "`attribution.targets` reads a run's elements out of the table its result's "
                "`elements_where` field points at, so the path comes out of a result file at run time "
                "and no reading of the source could name it"
            ),
        },
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
            "mpra_active_log2": ms.MPRA_ACTIVE,
            "eqtl_window_bp": payload["conventions"]["eqtl_window_bp"],
            "page_size": rm3.PAGE,
            "chromosomes": chroms,
        },
        "exclusions": [
            "every locus that left the map did so under one of the named causes in "
            "`reconciliation.loci_excluded_by_cause`, and the build refuses to return a payload whose "
            "built and excluded loci do not sum to the population counted before it",
            "VISTA and saturation mutagenesis contribute nothing: no element of this map carries one, "
            "and both are named with what they measure rather than left out silently",
            "an assertion that cannot name its source path, record key and sha256 is excluded under a "
            "cause and never shown",
        ],
        "partitions": {
            "kinds_of_evidence": "what each locus carries beyond the perturbation and the reader",
            "observed_by_assay": "`observed` split by what measured it, never pooled",
        },
        "code_cleanliness": mf.code_cleanliness(ENTRY, OWN_CODE, ROOT),
    }

    path = save_result(name, payload)
    print(f"{name}: {path}")
    print(rm3.summary(payload))
    print(f"  per-element cache opens: {payload['per_element_response_cache']['opens']}")
    print(f"  pages of {rm3.PAGE}: {max(1, -(-len(payload['loci']) // rm3.PAGE))}")


if __name__ == "__main__":
    main()
