# SPDX-License-Identifier: AGPL-3.0-or-later
"""The blinded HCT116 count `genomeos/attribution/hct116_count.py` registered, and the pooled arms
against the two imported floors (`data/results/hct116_count.json`).

    uv run --frozen python scripts/hct116_count.py

In the registered order, `hct116_count.ORDER`:

1. **The eligibility join**, lane-direction's own committed `links_of`, blind to every predicted value.
   Nothing of it is changed, and the script reconciles its own eligible counts against the figures
   `data/results/direction_link.json` published before this lane existed.
2. **The blinded presence pass**: for each HCT116 link, whether the partial per-element cache
   `data/knowledge/alphagenome/elements_hct116` holds an entry for that link's own gene on the HCT116
   track; and for each K562 link, the same question of the finished sweep's own archives. Presence, never
   a value: the reader returns a `bool` and the gene names, and no value has a path out of it.
3. **The pooled arms against both imported floors**, per arm, by lane-direction's own `gate`.
4. **The verdict** in the words `hct116_count.GATE_PASS` and `hct116_count.GATE_NO_GO`, and nothing
   beyond it.

**No direction is read anywhere.** Not a predicted value, not a sign, not a magnitude, not an aggregate
of any of them. That is structural and not careful: see `hct116_count.BLINDING` and the three tests in
`tests/test_hct116_count.py`.

**How the caches are read.** The main sweep's `<chrom>.json.gz` archives are STREAMED, not loaded:
`cached_records` decompresses in chunks and materialises one record at a time, so peak memory is one
record and not one chromosome. That reader is lane-direction's own, imported from
`scripts/direction_link.py` by `importlib.util.spec_from_file_location` rather than copied, the way
`scripts/clause2_measured_arm.py` imports its own predecessor. The partial HCT116 cache is one small
JSON file per element, and **every per-element open this run makes is recorded in the result by path**.

**What the manifest declares.** `targets.run_elements` takes a path out of a result file's
`elements_where` FIELD, so no grep over this script can see which element tables it opens. This script
therefore resolves those fields itself and declares both the pointer results and the tables behind them,
alongside the benchmark tables, the two cache groups and lane-direction's own result.
"""

from __future__ import annotations

import importlib.util
import sys
import time
from collections import Counter, defaultdict
from datetime import date
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import crispri  # noqa: E402
from genomeos.attribution import direction_link as dl  # noqa: E402
from genomeos.attribution import hct116_count as hc  # noqa: E402
from genomeos.attribution import measured as ms  # noqa: E402
from genomeos.attribution import targets as tg  # noqa: E402
from genomeos.results import RESULTS_DIR, load_result, save_result  # noqa: E402

RESULT = "hct116_count"
CHROMS = tuple(f"chr{c}" for c in [*range(1, 23), "X", "Y"])
HCT116_CACHE = Path(hc.CACHE_ROOT)
MAIN_CACHE = Path("data/knowledge/alphagenome/elements")

#: This lane's own files. Everything else uncommitted in this shared checkout belongs to another lane.
OWN_CODE = (
    "genomeos/attribution/hct116_count.py",
    "scripts/hct116_count.py",
    "tests/test_hct116_count.py",
)

INPUT_GLOBS = (
    "enhancer_targets_chr*.json",
    "enhancer_targets_all_chr*.json",
    "constrained_targets_chr*.json",
)


def _predecessor() -> Any:
    """lane-direction's computing script, imported rather than copied: its `links_of` is the committed
    eligibility join and its `cached_records` is the streaming reader this run reuses unchanged."""
    path = Path(__file__).resolve().parent / "direction_link.py"
    spec = importlib.util.spec_from_file_location("lane_direction_script", path)
    if spec is None or spec.loader is None:  # pragma: no cover - the file is committed beside this one
        raise SystemExit(f"{path} is not importable; lane-direction's reader is required")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _counts(rows: list[dict[str, Any]], field: str) -> dict[str, int]:
    c: Counter[str] = Counter()
    for r in rows:
        v = r[field]
        for x in v if isinstance(v, list) else [v]:
            c[x] += 1
    return dict(sorted(c.items()))


def hct116_presence(rows: list[dict[str, Any]]) -> dict[str, Any]:
    """The blinded presence pass over the partial cache. Presence and the path opened, never a value."""
    present: list[dict[str, Any]] = []
    no_entry: list[dict[str, Any]] = []
    no_file: list[dict[str, Any]] = []
    opens: list[str] = []
    for r in sorted(rows, key=lambda q: (q["chrom"], q["element"], q["gene"])):
        got = hc.answer_presence(HCT116_CACHE, r["chrom"], r["element"], r["gene"], hc.CELL)
        if got["opened"]:
            opens.append(got["opened"])
        row = {
            "arm": r["arm"],
            "cell": r["cell"],
            "chrom": r["chrom"],
            "element": r["element"],
            "start": r["start"],
            "end": r["end"],
            "gene": r["gene"],
            "has_answer": got["has_answer"],
            "why": got["why"],
            "opened": got["opened"],
        }
        if got["has_answer"]:
            present.append({**r, **row})
        elif got["why"] == hc.NO_FILE:
            no_file.append(row)
        else:
            no_entry.append(row)
    return {
        "present": present,
        "no_entry": no_entry,
        "no_file": no_file,
        "per_element_opens": sorted(set(opens)),
    }


def k562_presence(rows: list[dict[str, Any]], records_of: Any) -> dict[str, Any]:
    """The same presence question of the finished sweep's own archives, for the primary cell.

    The archives are streamed by lane-direction's `cached_records`; presence here is `the record holds an
    entry for this gene on this cell's track`, and the entry is never read. lane-direction reported 0
    exactly-zero values over all 230 of its in-cell links, so this presence count is expected to
    reproduce its `answered` counts exactly, and the result records whether it does.
    """
    by_chrom: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for r in rows:
        by_chrom[r["chrom"]].append(r)
    present: list[dict[str, Any]] = []
    absent: list[dict[str, Any]] = []
    archives: list[str] = []
    for chrom in sorted(by_chrom):
        path = MAIN_CACHE / f"{chrom}.json.gz"
        want = by_chrom[chrom]
        records = records_of(path, {r["element"] for r in want})
        if path.exists():
            archives.append(path.as_posix())
        for r in want:
            names = hc.genes_with_an_answer(records.get(r["element"]), hc.PRIMARY_CELL)
            (present if r["gene"] in names else absent).append(r)
    return {"present": present, "absent": absent, "archives": sorted(archives)}


def declared_element_tables(chroms: list[str]) -> tuple[list[Path], list[str]]:
    """Every element table this run's `targets.attributed` actually opens, and every pointer that is
    dangling, resolved out of the `elements_where` FIELD a grep cannot see."""
    tables: list[Path] = []
    missing: list[str] = []
    for name in tg.RUNS:
        for chrom in chroms:
            r = load_result(f"{name}_{chrom}", RESULTS_DIR) or {}
            if "elements" in r:
                continue
            where = r.get("elements_where")
            if not where:
                continue
            p = Path(where)
            if not p.is_absolute() and not p.exists():
                p = RESULTS_DIR.parent.parent / where
            if p.exists():
                tables.append(p)
            else:
                missing.append(str(where))
    return sorted(set(tables), key=lambda q: q.as_posix()), sorted(set(missing))


def main() -> None:
    started = time.time()
    if not HCT116_CACHE.exists():
        raise SystemExit(
            f"{HCT116_CACHE} is not on this machine: that is this lane's result, and nothing is "
            "substituted for it"
        )
    pre = _predecessor()
    chroms = list(CHROMS)

    links: list[dict[str, Any]] = []
    contested: list[dict[str, Any]] = []
    for chrom in chroms:
        got, bad = pre.links_of(chrom)
        links.extend(got)
        contested.extend(bad)

    eligible = {
        "call": dl.ELIGIBILITY,
        "join": "lane-direction's own committed `links_of`, imported and unchanged",
        "links": len(links),
        "by_arm": {a: sum(1 for r in links if r["arm"] == a) for a in dl.ARMS},
        "by_cell": _counts(links, "cell"),
        "by_arm_and_cell": {a: _counts([r for r in links if r["arm"] == a], "cell") for a in dl.ARMS},
        "contested": len(contested),
        "contest_rule": dl.CONTEST_RULE,
    }
    published = hc.INHERITED_FIGURES
    eligible["reconciles_with_lane_direction"] = {
        "file": published["file"],
        "links": eligible["links"] == published["eligible_links"],
        "by_arm": eligible["by_arm"] == published["eligible_by_arm"],
        "by_cell": eligible["by_cell"] == published["eligible_by_cell"],
        "by_arm_and_cell": eligible["by_arm_and_cell"] == published["eligible_by_arm_and_cell"],
        "published": {
            "links": published["eligible_links"],
            "by_arm": published["eligible_by_arm"],
            "by_cell": published["eligible_by_cell"],
            "by_arm_and_cell": published["eligible_by_arm_and_cell"],
        },
    }

    mine = [r for r in links if r["cell"] == hc.CELL]
    primary = [r for r in links if r["cell"] == hc.PRIMARY_CELL]
    got = hct116_presence(mine)
    k = k562_presence(primary, pre.cached_records)

    cache_chroms = sorted(p.name for p in HCT116_CACHE.iterdir() if p.is_dir())
    cache_files = sorted(p.as_posix() for p in HCT116_CACHE.rglob("*.json"))
    count = {
        "cell": hc.CELL,
        "cache": hc.CACHE,
        "cache_shape": hc.CACHE_SHAPE,
        "cache_chromosomes": cache_chroms,
        "cache_element_files": len(cache_files),
        "presence_call": hc.PRESENCE_CALL,
        "blinding": hc.BLINDING,
        "presence_is_an_upper_bound": hc.PRESENCE_IS_AN_UPPER_BOUND,
        "eligible_links_in_this_cell": len(mine),
        "eligible_by_arm": {a: sum(1 for r in mine if r["arm"] == a) for a in dl.ARMS},
        "links_with_an_answer": len(got["present"]),
        "links_with_an_answer_by_arm": {a: sum(1 for r in got["present"] if r["arm"] == a) for a in dl.ARMS},
        "answer_presence_gap": {
            "call": hc.ANSWER_PRESENCE_GAP,
            "links": len(mine) - len(got["present"]),
            "no_cached_element_file": {
                "links": len(got["no_file"]),
                "by_arm": {a: sum(1 for r in got["no_file"] if r["arm"] == a) for a in dl.ARMS},
                "why": hc.NO_FILE,
            },
            "element_cached_but_no_entry_for_the_pairs_own_gene": {
                "links": len(got["no_entry"]),
                "by_arm": {a: sum(1 for r in got["no_entry"] if r["arm"] == a) for a in dl.ARMS},
                "why": hc.NO_ENTRY,
            },
            "breakdown_reconciles": len(got["present"]) + len(got["no_entry"]) + len(got["no_file"])
            == len(mine),
        },
        "by_chromosome": _counts(mine, "chrom"),
        "links_with_an_answer_by_chromosome": _counts(got["present"], "chrom"),
        "per_element_opens": got["per_element_opens"],
        "per_element_opens_count": len(got["per_element_opens"]),
        "links": [
            {kk: r[kk] for kk in ("arm", "chrom", "element", "start", "end", "gene", "has_answer", "why")}
            for r in sorted(
                [*got["present"], *got["no_entry"], *got["no_file"]],
                key=lambda q: (q["arm"], q["chrom"], q["element"], q["gene"]),
            )
        ],
    }

    primary_side = {
        "cell": hc.PRIMARY_CELL,
        "cache": dl.CACHE,
        "streamed_not_loaded": (
            "the archives are decompressed in chunks and one record is materialised at a time, by "
            "lane-direction's own `cached_records` imported from scripts/direction_link.py; json.load "
            "on these archives is the heavy read the project gates and this run does not take it"
        ),
        "eligible_links": len(primary),
        "links_with_an_answer": len(k["present"]),
        "links_with_an_answer_by_arm": {a: sum(1 for r in k["present"] if r["arm"] == a) for a in dl.ARMS},
        "absent": len(k["absent"]),
        "archives_read": k["archives"],
        "reconciles_with_lane_directions_answered": {
            "published": published["answered"],
            "published_predicted_zero_excluded": published["predicted_zero_excluded"],
            "matches": {
                a: sum(1 for r in k["present"] if r["arm"] == a) == published["answered"][a] for a in dl.ARMS
            },
            "why_it_should": (
                "lane-direction's answered links are the present-and-not-exactly-zero ones and it "
                "reported 0 exactly-zero over all 230 of its in-cell links, so a presence count must "
                "reproduce its answered counts exactly; where it does not, this lane's pooled counts "
                "are not comparable with that lane's and the mismatch is the finding"
            ),
        },
    }

    k_by_arm = {a: [r for r in k["present"] if r["arm"] == a] for a in dl.ARMS}
    mine_by_arm = {a: [r for r in got["present"] if r["arm"] == a] for a in dl.ARMS}
    gates = {a: hc.pooled_gate(k_by_arm[a], mine_by_arm[a], a) for a in dl.ARMS}
    read = hc.verdict(gates)

    payload: dict[str, Any] = {
        "result": RESULT,
        "date": date.today().isoformat(),
        "lane": "lane-hct",
        "registration": hc.registration(),
        "eligibility": eligible,
        "this_lanes_count": count,
        "primary_cell_side": primary_side,
        "pooled_gate": {
            "call": dl.GATE_CALL,
            "pooling": hc.POOLING,
            "floors": hc.FLOORS,
            "per_arm": gates,
            "both_arms_meet_both_floors": read["both_pooled_arms_meet_both_floors"],
        },
        "verdict": read,
        "inherited_gate_no_go_unchanged": hc.INHERITED_GATE_NO_GO,
        "named_not_read_by_lane_direction": hc.NAMED_NOT_READ,
        "nothing_registered_or_run_here": (
            "no direction test is registered and none is run in this lane; no direction, value or sign "
            "is read; no floor, predicate or registered reading of any lane is moved; lane-direction's "
            "no-go stands as that lane's reading of that lane's population"
        ),
        "requests": 0,
        "money": "none: every input was already on disk",
        "no_requests": hc.NO_REQUESTS,
    }

    paths: set[Path] = {Path("data/results/direction_link.json")}
    for name in ms.CRISPRI_FILES:
        p = crispri.KNOWLEDGE / name
        if p.exists():
            paths.add(p)
    # `measured.Layer.load` opens the whole layer, so this run OPENS the lentiMPRA, VISTA and
    # saturation-mutagenesis tables as well, although its counting path reads only `layer.crispri`. An
    # input is what a writer opened and not what it meant to read (the audit hook in genomeos/manifest.py
    # records every open under data/), so they are declared, enumerated exactly as
    # `measured.result_manifest` enumerates them rather than named by hand. Imported here and not at the
    # top, the way that function imports its own helper, so the declaration cannot disturb a line of the
    # import block another lane's diff is reading.
    from genomeos.attribution import mpra, vista
    from genomeos.knowledge import satmut as satmut_knowledge

    for acc in mpra.FILES.values():
        p = mpra.KNOWLEDGE / f"{acc}.bed.gz"
        if p.exists():
            paths.add(p)
    for p in (vista.locus_path(vista.KNOWLEDGE), satmut_knowledge.DATA_PATH):
        if p.exists():
            paths.add(p)
    for glob in INPUT_GLOBS:
        for chrom in chroms:
            paths.update(RESULTS_DIR.glob(glob.replace("chr*", chrom)))
    tables, dangling = declared_element_tables(chroms)
    paths.update(tables)
    inputs = [mf.input_entry(p, partition=None) for p in sorted(paths, key=lambda q: q.as_posix())]
    if got["per_element_opens"]:
        inputs.append(
            mf.files_entry(
                "elements_hct116_per_element_opens",
                got["per_element_opens"],
                partition="the partial HCT116 per-element cache, read for answer PRESENCE only",
            )
        )
    if k["archives"]:
        inputs.append(
            mf.files_entry(
                "per_element_response_cache",
                k["archives"],
                partition="the finished deletion sweep's own per-element cache, read for presence only",
            )
        )

    payload["seconds"] = round(time.time() - started, 1)
    payload["result_manifest"] = {
        "sources": [
            {
                "accession": "ENCODE CRISPRi enhancer-gene benchmark "
                "(EngreitzLab/CRISPR_comparison, Gschwind et al.)",
                "version": "the two benchmark tables cached under data/knowledge, digested in inputs",
            },
            {
                "accession": "the finished genome-wide AlphaGenome deletion sweep's per-element "
                "response cache",
                "version": "the per-chromosome archives this run streamed, digested as a group",
            },
            {
                "accession": "the partial HCT116 per-element response cache "
                "data/knowledge/alphagenome/elements_hct116, written by an earlier run",
                "version": "every per-element file this run opened, digested as a group",
            },
            {
                "accession": "the attributed elements of the chromosomes this run covers",
                "version": "re-enumerated in this run from the results on disk, with the element "
                "tables behind the `elements_where` pointers declared by path",
            },
        ],
        "inputs": inputs,
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "order": hc.ORDER,
            "eligibility_predicate": dl.ELIGIBILITY,
            "presence_call": hc.PRESENCE_CALL,
            "blinding": hc.BLINDING,
            "presence_is_an_upper_bound": hc.PRESENCE_IS_AN_UPPER_BOUND,
            "answer_presence_gap": hc.ANSWER_PRESENCE_GAP,
            "independent_locus_rule": hc.LOCUS_RULE,
            "independent_locus_span": hc.LOCUS_SPAN,
            "not_biological_independence": hc.NOT_BIOLOGICAL_INDEPENDENCE,
            "pooling": hc.POOLING,
            "locus_floor": hc.LOCUS_FLOOR,
            "link_floor": hc.POSITIVE_FLOOR,
            "floors_imported_from": hc.FLOORS,
            "reciprocal_overlap": ms.RECIPROCAL_OVERLAP,
            "reach": ms.REACH,
            "cells_counted": [hc.CELL, hc.PRIMARY_CELL],
            "chromosomes": chroms,
            "aggregate_function": "a count of links, and a count of independent loci under cell2.group",
            "fill_value": "none: a link with no entry is counted as having no answer, never defaulted",
            "element_tables_opened": [p.as_posix() for p in tables],
            "element_pointers_dangling": dangling,
            "why_the_tables_are_declared": (
                "targets.run_elements takes a path out of a result file's `elements_where` FIELD, so no "
                "grep over this script can see which element tables it opens; this run resolves those "
                "fields itself and declares the tables by path alongside the pointer results"
            ),
        },
        "exclusions": [
            dl.CONTEST_RULE,
            hc.PRESENCE_IS_AN_UPPER_BOUND,
            "no predicted value, direction or sign is read, recorded or reported anywhere",
            "no cell, gene, value or cache is substituted for a link this cache cannot answer",
            "the two arms are never pooled together to reach a floor",
        ],
        "partitions": {
            "decreases": "links holding a significant decrease and no significant increase",
            "increases": "links holding a significant increase and no significant decrease",
            "contested": dl.CONTEST_RULE,
        },
        "code_cleanliness": mf.code_cleanliness(Path(__file__).resolve(), OWN_CODE),
    }

    path = save_result(RESULT, payload)
    print(f"{RESULT}: {path}")
    print(
        f"  eligible links: {eligible['links']} ({eligible['by_arm']}), "
        f"{hc.CELL} {eligible['by_cell'].get(hc.CELL, 0)}"
    )
    print(
        f"  {hc.CELL}: {count['links_with_an_answer']} of {count['eligible_links_in_this_cell']} links "
        f"have an answer {count['links_with_an_answer_by_arm']}; gap "
        f"{count['answer_presence_gap']['links']} "
        f"(no file {len(got['no_file'])}, no entry for the pair's own gene {len(got['no_entry'])}); "
        f"{count['per_element_opens_count']} per-element opens"
    )
    print(
        f"  {hc.PRIMARY_CELL}: {primary_side['links_with_an_answer']} of {len(primary)} "
        f"{primary_side['links_with_an_answer_by_arm']}, reconciles "
        f"{primary_side['reconciles_with_lane_directions_answered']['matches']}"
    )
    for a in dl.ARMS:
        g = gates[a]
        print(
            f"  pooled {a}: {g['links_with_an_answer']} links "
            f"({g['from_the_primary_cell']} + {g['from_this_lanes_cell']}), "
            f"{g['independent_loci']} independent loci -> "
            f"{'floors met' if g['meets_both_floors'] else '; '.join(g['short_by'])}"
        )
    print(f"  both pooled arms meet both floors: {read['both_pooled_arms_meet_both_floors']}")


if __name__ == "__main__":
    main()
