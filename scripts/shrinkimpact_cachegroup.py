#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""Why four results of 2026-10-02 declare 22 per-element cache archives, and what it changed.

    uv run python scripts/shrinkimpact_cachegroup.py

Run it in a clean worktree with the git-ignored stores linked read-only
(`scripts/link_stores_read_only.py`). Nothing is rewritten.

THE FOUR: `direction_link`, `hct116_count`, `reptest_answerable` and `rescore_cost`, each declaring
a `per_element_response_cache` group of 22 members against a 24-chromosome complement, chr13 and
chrY the two it does not name.

HOW THAT GROUP IS BUILT, in all four: the run groups its own rows by chromosome, walks
`sorted(by_chrom)`, and appends `data/knowledge/alphagenome/elements/<chrom>.json.gz` to the group
IF `path.exists()`. That `.exists()` is the defect this lane investigates: had an archive been
absent, every row on that chromosome would have been counted unanswerable and the declaration would
simply not have mentioned the file. So there are exactly two things to establish, and this result
establishes both by measurement rather than by argument:

1. WAS AN ARCHIVE THE RUN NEEDED ABSENT? Every one of the 24 archives is listed with its sha256 and
   byte count as it is on disk now, and every declared member's digest is checked against the
   manifest's. A declared member whose bytes differ, or an archive missing from the store, is named.

2. DID THE RUN NEED chr13 OR chrY AT ALL? The upstream population is the cached ENCODE CRISPRi
   benchmark: `measured.Layer.load(chrom).crispri`, which `links_of` reads before anything else and
   returns empty on. The number of cached CRISPRi pairs and of eligible links is recomputed for all
   24 chromosomes, so "the group is 22 because the benchmark has rows on 22 chromosomes" is either
   shown or refuted.

AND FOR `direction_link` SPECIFICALLY, the counts are recomputed, not just the scope: its whole
`eligibility` block and the four `answerability` counts (`links_in_the_cell_read`, `answered`,
`predicted_zero_excluded`, `absent`) are re-derived with the complete stores and compared leaf by
leaf. `absent` is the count that would have moved had an archive been missing, so it is the one
worth the 22 archive streams.

THE LIMIT, NAMED AND NOT BURIED: the other three results' own counts are NOT recomputed here. Each
runs a different downstream pipeline over the same upstream population, and this result speaks about
their DECLARED GROUP and the store behind it, not about their arm figures. What is established for
them is the scope of the group and the integrity of the files it names; what is not is left to be
said in those words.

NO NETWORK, NO MODEL, NO MONEY.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
import time
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import measured as ms  # noqa: E402
from genomeos.attribution import shrinkimpact as si  # noqa: E402
from genomeos.results import save_result  # noqa: E402

NAME = "shrinkimpact_per_element_cache_group"
ENTRY = "scripts/shrinkimpact_cachegroup.py"
OWN_CODE = (
    "genomeos/attribution/shrinkimpact.py",
    "scripts/shrinkimpact_cachegroup.py",
    "tests/test_shrinkimpact.py",
)
CHROMS = (*(f"chr{i}" for i in range(1, 23)), "chrX", "chrY")
CACHE_DIR = Path("data/knowledge/alphagenome/elements")
GROUP_LABEL = "per_element_response_cache"
THE_FOUR = ("direction_link", "hct116_count", "reptest_answerable", "rescore_cost")


def _direction_link_module() -> Any:
    """`scripts/direction_link.py` as a module, so its own `links_of` and `answer` are reused."""
    path = ROOT / "scripts/direction_link.py"
    spec = importlib.util.spec_from_file_location("direction_link", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def archive_census() -> dict[str, Any]:
    """Every one of the 24 per-chromosome archives as it is on disk now, with its digest."""
    out: dict[str, Any] = {}
    for chrom in CHROMS:
        path = ROOT / CACHE_DIR / f"{chrom}.json.gz"
        if not path.exists():
            out[chrom] = {"present": False}
            continue
        digest = hashlib.sha256()
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1 << 22), b""):
                digest.update(chunk)
        out[chrom] = {
            "present": True,
            "sha256": digest.hexdigest(),
            "bytes": path.stat().st_size,
        }
    return out


def declared_group(name: str) -> dict[str, Any]:
    """The `per_element_response_cache` group of one committed result, by its declared members."""
    manifest = json.loads((ROOT / f"data/results/{name}.json").read_text())["result_manifest"]
    entry = next(
        (e for e in manifest["inputs"] if isinstance(e, dict) and e.get("path") == GROUP_LABEL),
        None,
    )
    if entry is None:
        raise RuntimeError(f"{name} declares no {GROUP_LABEL} group")
    members = {Path(m["path"]).name.split(".")[0]: m for m in entry["members"]}
    return {"files": entry["files"], "members": members}


def check_declared_against_disk(name: str, census: dict[str, Any]) -> dict[str, Any]:
    """Does every file the result declared still hold the bytes it declared, and is any missing?"""
    group = declared_group(name)
    same: list[str] = []
    different: list[dict[str, Any]] = []
    gone: list[str] = []
    for chrom, member in sorted(group["members"].items()):
        on_disk = census.get(chrom) or {"present": False}
        if not on_disk["present"]:
            gone.append(chrom)
        elif on_disk["sha256"] == member["sha256"]:
            same.append(chrom)
        else:
            different.append(
                {
                    "chrom": chrom,
                    "declared_sha256": member["sha256"],
                    "on_disk_sha256": on_disk["sha256"],
                }
            )
    undeclared_but_present = sorted(
        c for c in CHROMS if c not in group["members"] and (census.get(c) or {}).get("present")
    )
    return {
        "declared_members": group["files"],
        "declared_chromosomes": sorted(group["members"]),
        "declared_and_byte_identical_on_disk": len(same),
        "declared_but_different_bytes": different,
        "declared_but_missing_from_the_store": gone,
        "not_declared_though_present_in_the_store": undeclared_but_present,
    }


def benchmark_population() -> dict[str, Any]:
    """The upstream population per chromosome: cached CRISPRi pairs, eligible links, contested."""
    direction = _direction_link_module()
    pairs: dict[str, int] = {}
    links: dict[str, int] = {}
    contested: dict[str, int] = {}
    all_links: list[dict[str, Any]] = []
    for chrom in CHROMS:
        layer = ms.Layer.load(chrom)
        pairs[chrom] = len(layer.crispri or ())
        got, bad = direction.links_of(chrom)
        links[chrom] = len(got)
        contested[chrom] = len(bad)
        all_links.extend(got)
    return {
        "cached_crispri_pairs_per_chromosome": pairs,
        "eligible_links_per_chromosome": links,
        "contested_per_chromosome": contested,
        "chromosomes_with_an_eligible_link": sorted(c for c, n in links.items() if n),
        "links": all_links,
        "module": direction,
    }


def direction_link_counts(direction: Any, links: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """`direction_link`'s own eligibility and answerability counts, recomputed and compared."""
    import genomeos.attribution.direction_link as dl

    committed = json.loads((ROOT / "data/results/direction_link.json").read_text())
    counts = direction._counts
    eligible = {
        "links": len(links),
        "by_arm": {a: sum(1 for r in links if r["arm"] == a) for a in dl.ARMS},
        "by_cell": counts(links, "cell"),
        "by_arm_and_cell": {a: counts([r for r in links if r["arm"] == a], "cell") for a in dl.ARMS},
    }
    was = committed["eligibility"]
    reports = [
        si.compare(
            {k: was[k] for k in eligible},
            eligible,
            label="direction_link eligibility: links, by_arm, by_cell, by_arm_and_cell",
        )
    ]
    cell = dl.PRIMARY_CELL
    in_cell = [r for r in links if r["cell"] == cell]
    got = direction.answer(in_cell, cell)
    answerability = {
        "links_in_the_cell_read": len(in_cell),
        "answered": len(got["answered"]),
        "predicted_zero_excluded": len(got["zero"]),
        "absent": len(got["absent"]),
    }
    reports.append(
        si.compare(
            {k: committed["answerability"][k] for k in answerability},
            answerability,
            label="direction_link answerability: the four counts, absent among them",
        )
    )
    reports.append(
        si.compare(
            {"archives": sorted(declared_group("direction_link")["members"])},
            {"archives": sorted(Path(p).name.split(".")[0] for p in got["archives"])},
            label="direction_link: the archives the rerun opened against the ones it declared",
        )
    )
    return reports


def recount() -> dict[str, Any]:
    """The whole thing, as the payload of its own result."""
    census = archive_census()
    population = benchmark_population()
    direction = population.pop("module")
    links = population.pop("links")
    per_result = {name: check_declared_against_disk(name, census) for name in THE_FOUR}
    with_links = population["chromosomes_with_an_eligible_link"]
    reports = direction_link_counts(direction, links)
    for name, found in sorted(per_result.items()):
        reports.append(
            si.compare(
                {"chromosomes": found["declared_chromosomes"]},
                {"chromosomes": with_links},
                label=f"{name}: the declared group's chromosomes against the benchmark's own",
            )
        )
    return {
        "question": "four results of 2026-10-02 declare a 22-member per-element cache group "
        "against a 24-chromosome complement, and the group is built with an `.exists()` filter, "
        "so the declaration would have shrunk silently had an archive been absent. Was one, and "
        "did the runs need chr13 or chrY at all?",
        "corrects": "nothing: this is an additive check on four committed results, none of which "
        "is touched here",
        "the_four": list(THE_FOUR),
        "comparisons": reports,
        "no_number_moves": si.every_reading_is_no_movement(reports),
        "every_archive_as_it_is_on_disk_now": census,
        "the_upstream_population": population,
        "per_result_declaration_against_the_store": per_result,
        "limit_named": "the own counts of hct116_count, reptest_answerable and rescore_cost are "
        "NOT recomputed here: each runs a different downstream pipeline over the same upstream "
        "population, and this result speaks about their declared group and the files behind it. "
        "`direction_link`'s own counts ARE recomputed, eligibility and all four answerability "
        "counts, because `absent` is the figure that would have moved had an archive been missing.",
        "wording": "a recomputed count is what the complete stores give; it does not by itself "
        "overturn a registered reading, and no registered test is re-run here",
        "stores": "a clean worktree of the committed tree with data/cache, data/reference and "
        "data/knowledge cloned in at 0444 (scripts/link_stores_read_only.py)",
        "alphagenome_requests": 0,
        "model_requests": 0,
        "money": "none: every byte read here is already on this machine",
    }


def result_manifest(read_now: list[str]) -> dict[str, Any]:
    """The contract's block, inputs repository-relative and grouped by what the code means."""
    groups = {
        "the finished sweep's per-element response archives, every one of the 24": [
            p for p in read_now if p.startswith(f"{CACHE_DIR.as_posix()}/") and p.endswith(".json.gz")
        ],
        "the four committed results whose declarations are under examination": [
            p for p in read_now if p.startswith("data/results/") and Path(p).stem in THE_FOUR
        ],
        "the cached ENCODE CRISPRi benchmark and the sweep's committed per-chromosome results": [
            p
            for p in read_now
            if (p.startswith("data/results/") and Path(p).stem not in THE_FOUR)
            or p.startswith("data/knowledge/crispri/")
        ],
    }
    inputs = [
        mf.files_entry(label, sorted(paths), partition="not an evaluation partition: the recount's inputs")
        for label, paths in groups.items()
        if paths
    ]
    grouped = {p for paths in groups.values() for p in paths}
    other = sorted(p for p in set(read_now) - grouped if (ROOT / p).exists())
    inputs += [mf.input_entry(p, partition="read by the recount") for p in other]
    return {
        "sources": [
            {
                "accession": ", ".join(f"data/results/{n}.json" for n in THE_FOUR),
                "version": "as committed 2026-10-02; opened read-only and not rewritten",
            },
            {
                "accession": "ENCODE CRISPRi enhancer-gene benchmark "
                "(EngreitzLab/CRISPR_comparison, Gschwind et al.)",
                "version": "the benchmark tables cached under data/knowledge, digested in inputs",
            },
            {
                "accession": "the finished genome-wide AlphaGenome deletion sweep's per-element "
                "response cache, written with threshold=0.0",
                "version": "all 24 per-chromosome archives on this machine, digested in the result",
            },
            {
                "accession": "the project's own all-enhancer deletion sweep (enhancer_targets_all_<chrom>)",
                "version": "as committed on this machine",
            },
        ],
        "inputs": inputs,
        "input_count": len(inputs),
        "assembly": "GRCh38: the sweep's elements and the benchmark's tested intervals are on one "
        "build, and no coordinate is lifted here",
        "coordinates": {
            "base": 0,
            "interval": "half-open",
            "note": "the committed join's own convention, unchanged: reciprocal overlap at "
            f"{ms.RECIPROCAL_OVERLAP} over the candidates within {ms.REACH} bp",
        },
        "parameters": {
            "chromosomes": list(CHROMS),
            "reciprocal_overlap": ms.RECIPROCAL_OVERLAP,
            "reach": ms.REACH,
            "results_examined": list(THE_FOUR),
            "group_label": GROUP_LABEL,
        },
        "exclusions": [
            "the own counts of hct116_count, reptest_answerable and rescore_cost: three different "
            "downstream pipelines, not recomputed here, and said so in `limit_named`",
            "hct116_count's separate `elements_hct116_per_element_opens` group: a partial cache of "
            "a different shape over five chromosomes, not the 24-chromosome archive set",
            "every registered test: re-running one on a changed population needs a new "
            "registration and is not this recount's to do",
        ],
        "partitions": "n/a: a check of four declarations against the store behind them",
        "code_cleanliness": mf.code_cleanliness(ENTRY, OWN_CODE, ROOT),
    }


def main() -> int:
    mf.trace_begin()
    started = time.time()
    print("checking the per-element cache group of the four results of 2026-10-02", flush=True)
    record = recount()
    record["seconds"] = round(time.time() - started, 1)
    record["result_manifest"] = result_manifest(mf.traced_reads())
    path = save_result(NAME, record)
    print(f"wrote {path}", flush=True)
    for report in record["comparisons"]:
        print("  " + report["reading"], flush=True)
    print(f"no_number_moves = {record['no_number_moves']}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
