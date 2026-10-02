# SPDX-License-Identifier: AGPL-3.0-or-later
"""How many of extractor v2's increase-derived links are answerable in their own cell, and nothing else.

    uv run --frozen python scripts/reptest_answerable.py
    uv run --frozen python scripts/reptest_answerable.py --chroms chr21 --result reptest_answerable_chr21

The predicate, the causes, the floors and the verdict rule are `genomeos.attribution.reptest`, written
and committed before this count. In the registered order:

1. **The population**, read out of the extractor and not rebuilt beside it: every record
   `measured.rule_links_detail(row, measured.EXTRACTOR_V2)` returns whose `derived_from` is a
   significant increase, over the rows `measured.rows` builds from `targets.attributed`. The decrease
   arm is enumerated in the same pass, because the two arms' denominators are what the collapse
   identity's premise is checked against. Blind to every cached value.
2. **The answerability pass**: one cause per link out of `reptest.CAUSES`, exhaustive and disjoint,
   from the cached signed deletion value for that link's OWN gene on that link's OWN cell track.
3. **The verdict** against both imported floors on the ANSWERABLE links and their independent loci.
4. **Stop** if either floor is missed. No direction is read, no agreement rate is computed and no
   registration of a test is written. If both clear, the result says a test may be registered and this
   script still reads no direction: the registration is a separate commit.

**No agreement is ever computed here, even where it would be free.** The cached value's SIGN is read,
because that is what distinguishes an answerable link from one whose value is exactly zero, and it is
never compared with the link's measured sign. That comparison is the direction test and it is gated.

**How the per-element cache is read.** Two readers, both committed behaviour and neither widened.
The per-chromosome archive is STREAMED, decompressed in chunks with only a wanted record ever
materialised, so a 77 MB chromosome costs one record of memory rather than one chromosome: that reader
is `scripts/direction_link.py`'s, carried unchanged. Elements the archive does not hold are then looked
up as loose per-element files under `data/knowledge/alphagenome/elements/<chrom>/<id>.json`, which is
what `targets.ElementResponses.element` does and what the streaming reader on its own does NOT do. Both
are the same element's own record for the same gene on the same cell, so reading the loose file is
reading what is already on disk and is not a substitution; every loose file opened is named in the
manifest. `json.load` over a whole archive, the heavy read the project gates, is not taken.

Nothing is changed. `measured.rows`, `measured.rule_links`, `measured.rule_links_detail`,
`measured.rule_link_conflicts` and `targets.attributed` are called exactly as committed. No threshold,
floor, predicate or registered reading is moved, nothing is refitted, nothing is downloaded and no
model request is made.
"""

from __future__ import annotations

import argparse
import gzip
import json
import re
import sys
import time
from collections import Counter, defaultdict
from datetime import date
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import (  # noqa: E402
    crispri,  # noqa: E402
    mpra,
    vista,
)
from genomeos.attribution import increase_links as il  # noqa: E402
from genomeos.attribution import measured as ms  # noqa: E402
from genomeos.attribution import reptest as rt  # noqa: E402
from genomeos.attribution.targets import attributed  # noqa: E402
from genomeos.knowledge import satmut as satmut_knowledge  # noqa: E402
from genomeos.results import RESULTS_DIR, save_result  # noqa: E402

RESULT = "reptest_answerable"
CHROMS = tuple(f"chr{c}" for c in [*range(1, 23), "X", "Y"])
CACHE_DIR = Path("data/knowledge/alphagenome/elements")

#: This lane's own files. Everything else uncommitted in this shared checkout belongs to another lane.
OWN_CODE = (
    "genomeos/attribution/reptest.py",
    "scripts/reptest_answerable.py",
    "tests/test_reptest.py",
)

#: The same per-chromosome results `targets.attributed` reads, each digested by its own path. The
#: `enhancer_targets_all` summaries are POINTERS: the table they name under `elements_where` is what the
#: run really opens, so `pointer_tables` collects those paths and declares them too.
INPUT_GLOBS = (
    "enhancer_targets_chr*.json",
    "enhancer_targets_all_chr*.json",
    "constrained_targets_chr*.json",
)

#: A top-level record of the per-element cache starts with its own id, which is what anchors the scan.
RECORD = re.compile(rb'"([^"]{1,64})":\s*\{"id":')


def _object_end(buf: bytes, start: int) -> int | None:
    """The index just past the balanced JSON object that begins at `buf[start] == '{'`, or None.

    Strings and their escapes are respected, so a brace inside a gene or tissue name cannot close a
    record early. None means the object is not yet complete in `buf`. Carried unchanged from
    `scripts/direction_link.py`, so the two lanes read the archive identically.
    """
    depth = 0
    in_string = False
    escaped = False
    for i in range(start, len(buf)):
        c = buf[i]
        if in_string:
            if escaped:
                escaped = False
            elif c == 0x5C:  # backslash
                escaped = True
            elif c == 0x22:  # quote
                in_string = False
            continue
        if c == 0x22:
            in_string = True
        elif c == 0x7B:  # {
            depth += 1
        elif c == 0x7D:  # }
            depth -= 1
            if depth == 0:
                return i + 1
    return None


def cached_records(path: Path, wanted: set[str], chunk: int = 1 << 23) -> dict[str, dict[str, Any]]:
    """The records of `wanted` out of one per-element cache archive, streamed.

    The archive is a single JSON object keyed by element id. It is decompressed in chunks and only a
    wanted record is ever accumulated, so a 77 MB chromosome costs one record of memory and not one
    chromosome. Carried unchanged from `scripts/direction_link.py`. Nothing is written and no value is
    read here: this returns the records, and the caller takes the one gene and one cell the
    registration allows it.
    """
    out: dict[str, dict[str, Any]] = {}
    if not wanted or not path.exists():
        return out
    remaining = set(wanted)
    buf = b""
    holding: str | None = None
    held = b""
    with gzip.open(path, "rb") as fh:
        while remaining:
            block = fh.read(chunk)
            if not block:
                break
            if holding is not None:
                held += block
                end = _object_end(held, 0)
                if end is None:
                    continue
                out[holding] = json.loads(held[:end])
                remaining.discard(holding)
                buf = held[end:]
                holding, held = None, b""
            else:
                buf += block
            while remaining:
                hit = next((m for m in RECORD.finditer(buf) if m.group(1).decode() in remaining), None)
                if hit is None:
                    break
                start = buf.index(b"{", hit.end() - len(b'{"id":'))
                end = _object_end(buf, start)
                name = hit.group(1).decode()
                if end is None:
                    holding, held = name, buf[start:]
                    buf = b""
                    break
                out[name] = json.loads(buf[start:end])
                remaining.discard(name)
                buf = buf[end:]
            if len(buf) > (1 << 20):  # keep only enough tail to span a record boundary
                buf = buf[-(1 << 20) :]
    return out


def links_of(chrom: str) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, int]]:
    """The extractor's own links on one chromosome, under v2, with the conflicts and the v1 check.

    `measured.rule_links_detail` is called as committed and nothing is re-derived from the pairs: the
    gene, the cell, the action, the split, the strength, the outcome and the by-construction property
    are the extractor's. The row supplies the element's interval, which the extractor's record does not
    carry and the locus grouping needs. Blind to every cached value.
    """
    layer = ms.Layer.load(chrom)
    elements = attributed(chrom, RESULTS_DIR)
    links: list[dict[str, Any]] = []
    conflicts: list[dict[str, Any]] = []
    counts = {"rows": 0, "v1_links": 0, "v2_links": 0}
    for row in ms.rows(chrom, elements, layer):
        if not row["measured"].get("crispri"):
            continue
        detail = ms.rule_links_detail(row, rt.EXTRACTOR)
        counts["rows"] += 1
        counts["v1_links"] += len(ms.rule_links(row, ms.EXTRACTOR_V1))
        counts["v2_links"] += len(detail)
        conflicts.extend(ms.rule_link_conflicts(row))
        for r in detail:
            links.append(
                {
                    **r,
                    "chrom": chrom,
                    "start": int(row["start"]),
                    "end": int(row["end"]),
                    "element": r["element_id"],
                }
            )
    il.check_no_mechanism_claim(links)
    return links, conflicts, counts


def answerability(links: list[dict[str, Any]]) -> dict[str, Any]:
    """One cause per link, out of the cache, for that link's own gene on its own cell track.

    The archive is streamed for the elements of links whose cell the cache carries at all; an element
    the archive does not hold is then looked up as a loose per-element file, which is what the
    committed reader does. No other element, gene, cell or cache is consulted.
    """
    by_chrom: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for r in links:
        by_chrom[r["chrom"]].append(r)
    out: list[dict[str, Any]] = []
    archives: list[str] = []
    loose: list[str] = []
    for chrom in sorted(by_chrom):
        rows = by_chrom[chrom]
        wanted = {r["element"] for r in rows if r["cell"] in rt.CACHED_CELLS}
        path = CACHE_DIR / f"{chrom}.json.gz"
        records = cached_records(path, wanted)
        if wanted and path.exists():
            archives.append(path.as_posix())
        for element in sorted(wanted - set(records)):
            p = CACHE_DIR / chrom / f"{element}.json"
            if p.exists():
                records[element] = json.loads(p.read_text())
                loose.append(p.as_posix())
        for r in rows:
            cause, value = rt.cause_of(records.get(r["element"]), r["gene"], r["cell"])
            row = {**r, "answerability": cause}
            if value is not None:
                row["cached_value"] = round(value, 6)
            out.append(row)
    return {"links": out, "archives": sorted(archives), "loose_files": sorted(loose)}


def _counts(rows: list[dict[str, Any]], field: str) -> dict[str, int]:
    c: Counter[str] = Counter()
    for r in rows:
        v = r[field]
        for x in v if isinstance(v, list) else [v]:
            c[str(x)] += 1
    return dict(sorted(c.items()))


def arm_payload(rows: list[dict[str, Any]], arm: str) -> dict[str, Any]:
    """One arm's population and its answerability, apart from the other. No agreement is computed."""
    answerable = [r for r in rows if r["answerability"] == rt.ANSWERABLE]
    causes = _counts(rows, "answerability")
    return {
        "arm": arm,
        "measured_links": len(rows),
        "answerable_in_its_own_cell": len(answerable),
        "independent_loci_of_the_answerable_links": rt.independent_loci(answerable),
        "causes": causes,
        "causes_sum_to_the_population": rt.causes_reconcile(causes, len(rows)),
        "by_cell": _counts(rows, "cell"),
        "by_split": _counts(rows, "split"),
        "answerable_by_cell": _counts(answerable, "cell"),
        "answerable_by_split": _counts(answerable, "split"),
        "answerable_by_chromosome": _counts(answerable, "chrom"),
        "no_agreement_computed": (
            "the cached value's sign is read only to tell an answerable link from one whose value is "
            "exactly zero. It is never compared with the link's measured sign here: that comparison is "
            "the direction test and it is gated on the floors"
        ),
    }


def pointer_tables(chroms: list[str]) -> list[Path]:
    """The tables the `enhancer_targets_all` summaries point at through `elements_where`.

    A declared input must be the bytes the run really opened. These summaries are pointers: each names
    a table under `data/knowledge/alphagenome/all_elements`, and `targets.run_elements` opens that
    table. `save_result` refuses a result whose traced reads exceed its declared inputs, and it exists
    because the 24 summaries declare about 39 KB of pointers and read about 1.4 GB of tables.
    """
    out: list[Path] = []
    for chrom in chroms:
        p = RESULTS_DIR / f"enhancer_targets_all_{chrom}.json"
        if not p.exists():
            continue
        where = json.loads(p.read_text()).get("elements_where")
        if not where:
            continue
        q = Path(where)
        if not q.is_absolute() and not q.exists():
            q = RESULTS_DIR.parent.parent / where
        if q.exists():
            out.append(q)
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--chroms", default=",".join(CHROMS))
    ap.add_argument("--result", default=RESULT, help="a run over fewer chromosomes writes its own name")
    args = ap.parse_args()
    chroms = [c for c in args.chroms.split(",") if c]
    if args.result == RESULT and tuple(chroms) != CHROMS:
        raise SystemExit(f"{RESULT} is the genome-wide count; a run over {len(chroms)} needs --result")

    started = time.time()
    links: list[dict[str, Any]] = []
    conflicts: list[dict[str, Any]] = []
    extracted = {"rows": 0, "v1_links": 0, "v2_links": 0}
    for chrom in chroms:
        got, bad, counts = links_of(chrom)
        links.extend(got)
        conflicts.extend(bad)
        for k, v in counts.items():
            extracted[k] += v

    got = answerability(links)
    scored = got["links"]
    increases = [r for r in scored if r["derived_from"] == ms.INCREASE]
    decreases = [r for r in scored if r["derived_from"] == ms.DECREASE]
    answerable_up = [r for r in increases if r["answerability"] == rt.ANSWERABLE]
    answerable_down = [r for r in decreases if r["answerability"] == rt.ANSWERABLE]
    verdict = rt.verdict(answerable_up, len(increases))

    payload: dict[str, Any] = {
        "result": args.result,
        "date": date.today().isoformat(),
        "lane": rt.LANE,
        "question": (
            "how many of extractor v2's increase-derived links are ANSWERABLE IN THEIR OWN CELL, "
            "counted before any direction test and reported beside the measured count"
        ),
        "registered_in_code_before_this_count": "genomeos/attribution/reptest.py",
        "registration": rt.registration(),
        "binding_wording_r2": il.R2_WORDING,
        "count_is_not_coverage": rt.COUNT_IS_NOT_COVERAGE,
        "chromosomes": chroms,
        "population": {
            "call": rt.POPULATION_IS_THE_EXTRACTORS_OUTPUT,
            "extractor": rt.EXTRACTOR,
            "rows_with_a_crispri_measurement": extracted["rows"],
            "v1_links": extracted["v1_links"],
            "v2_links": extracted["v2_links"],
            "increase_derived_links": len(increases),
            "decrease_derived_links": len(decreases),
            "v1_links_removed_by_v2": extracted["v1_links"] - len(decreases),
            "additive": extracted["v2_links"] == extracted["v1_links"] + len(increases),
            "conflicts": len(conflicts),
            "conflict_rows": conflicts,
            "conflict_is_not_a_finding": il.CONFLICT_IS_NOT_A_FINDING,
            "by_construction": rt.POPULATION_IS_NOT_COVERAGE,
        },
        "agrees_with_lane_direction": {
            "what_is_compared": (
                "this lane counted the extractor's own output; lane-direction enumerated an equivalent "
                "population by hand from the cached pairs. Both counts are named so that the agreement "
                "is a finding rather than an assumption, and neither is restated as the other's"
            ),
            "lane_direction_figures": dict(rt.INHERITED_FIGURES),
            "increase_links_here": len(increases),
            "increase_links_there": rt.INHERITED_FIGURES["lane_direction_increase_arm_measured"],
            "increase_links_agree": len(increases)
            == rt.INHERITED_FIGURES["lane_direction_increase_arm_measured"],
            "answerable_increase_links_here": len(answerable_up),
            "answerable_increase_links_there": rt.INHERITED_FIGURES["lane_direction_increase_arm_answered"],
            "answerable_increase_links_agree": len(answerable_up)
            == rt.INHERITED_FIGURES["lane_direction_increase_arm_answered"],
        },
        "answerability": {
            "call": rt.ANSWERABLE_CALL,
            "predicted_direction": rt.PREDICTED_CALL,
            "sign_rule_zero_and_absent": rt.SIGN_CALL,
            "causes": {c: rt.CAUSE_TEXT[c] for c in rt.CAUSES},
            "causes_are_exhaustive": rt.CAUSES_ARE_EXHAUSTIVE,
            "cached_cell_tracks": list(rt.CACHED_CELLS),
            "cells_not_cached": list(rt.CELLS_NOT_CACHED),
            "archives_streamed": got["archives"],
            "loose_per_element_files_opened": got["loose_files"],
        },
        "arms": {
            "increases": arm_payload(increases, "increases"),
            "decreases": arm_payload(decreases, "decreases"),
            "reported_apart": (
                "the two arms are reported apart and are never pooled to reach a floor. The decreases "
                "arm is here because the collapse identity's premise is checked against both answered "
                "denominators, and its own count is the extractor's construction and NOT model "
                "coverage: " + rt.DECREASE_ARM_IS_NOT_MODEL_COVERAGE
            ),
        },
        "verdict": verdict,
        "statistic": rt.collapse_applies(len(answerable_down), len(answerable_up)),
        "increase_links": increases,
    }

    paths: set[Path] = set()
    for name in ("increase_links_registration", "increase_population", "direction_link"):
        p = RESULTS_DIR / f"{name}.json"
        if p.exists():
            paths.add(p)
    for name in ms.CRISPRI_FILES:
        p = crispri.KNOWLEDGE / name
        if p.exists():
            paths.add(p)
    for acc in mpra.FILES.values():
        p = mpra.KNOWLEDGE / f"{acc}.bed.gz"
        if p.exists():
            paths.add(p)
    for p in (vista.locus_path(vista.KNOWLEDGE), satmut_knowledge.DATA_PATH):
        if Path(p).exists():
            paths.add(Path(p))
    for glob in INPUT_GLOBS:
        for chrom in chroms:
            paths.update(RESULTS_DIR.glob(glob.replace("chr*", chrom)))
    paths.update(pointer_tables(chroms))
    inputs = [mf.input_entry(p, partition=None) for p in sorted(paths, key=lambda q: q.as_posix())]
    cache_files = got["archives"] + got["loose_files"]
    if cache_files:
        inputs.append(
            mf.files_entry(
                "per_element_response_cache",
                cache_files,
                partition="the finished deletion sweep's own per-element cache, read for one gene on "
                "one cell track per link and for nothing else",
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
                "response cache, written with threshold=0.0",
                "version": "the per-chromosome archives and loose per-element files this run opened, "
                "digested as a group",
            },
            {
                "accession": "the attributed elements of the chromosomes this run covers",
                "version": "re-enumerated in this run from the results on disk and the tables their "
                "`elements_where` fields point at, not read from one file",
            },
        ],
        "inputs": inputs,
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "extractor": rt.EXTRACTOR,
            "population": rt.POPULATION_IS_THE_EXTRACTORS_OUTPUT,
            "answerable_predicate": rt.ANSWERABLE_CALL,
            "predicted_direction": rt.PREDICTED_CALL,
            "sign_rule_zero_and_absent": rt.SIGN_CALL,
            "causes": list(rt.CAUSES),
            "answerable_links_floor": rt.POSITIVE_FLOOR,
            "independent_loci_floor": rt.LOCUS_FLOOR,
            "floors": dict(rt.FLOORS),
            "gate_call": rt.GATE_CALL,
            "independent_locus_rule": rt.LOCUS_RULE,
            "independent_locus_span": rt.LOCUS_SPAN,
            "not_biological_independence": rt.NOT_BIOLOGICAL_INDEPENDENCE,
            "reciprocal_overlap": ms.RECIPROCAL_OVERLAP,
            "reach": ms.REACH,
            "statistic_a_test_would_be_registered_on": rt.TEST_STATISTIC,
            "chance": rt.CHANCE,
            "why_a_between_arm_difference_is_not_it": rt.COLLAPSE,
            "cached_cell_tracks": list(rt.CACHED_CELLS),
            "cells_not_cached": list(rt.CELLS_NOT_CACHED),
            "chromosomes": chroms,
            "aggregate_function": "a count per cause; no rate, no mean and no interval is taken",
            "fill_value": "none: an unanswerable link is counted under its cause and nothing is "
            "substituted for it",
            "model_requests": 0,
        },
        "exclusions": [
            rt.CAUSE_TEXT[rt.CELL_NOT_CACHED],
            rt.CAUSE_TEXT[rt.VALUE_EXACTLY_ZERO],
            rt.SIGN_CALL,
            il.CONFLICT_RULE,
            "no agreement rate, interval, control or direction is computed by this run under any "
            "outcome of the gate",
            "no pair is excluded for its magnitude: no effect-size bar is applied anywhere",
        ],
        "partitions": {
            "increases": "links the extractor derived from a significant increase under v2; action "
            "`inhibits`, outcome `increase on knockdown`, molecular role unresolved",
            "decreases": "links the extractor derived from a benchmark `Regulated` pair, which is a "
            "significant decrease only; v1's own links, unchanged under v2",
            "conflicts": il.CONFLICT_RULE,
            "training_and_heldout": "the benchmark's own split, one per table, as the extractor "
            "assigns it to a link",
        },
        "code_cleanliness": mf.code_cleanliness(Path(__file__).resolve(), OWN_CODE),
    }

    path = save_result(args.result, payload)
    print(f"{args.result}: {path}")
    print(
        f"  population under {rt.EXTRACTOR}: {extracted['v2_links']} links = "
        f"{len(decreases)} decrease-derived + {len(increases)} increase-derived, "
        f"v1 {extracted['v1_links']}, conflicts {len(conflicts)}"
    )
    for arm in ("increases", "decreases"):
        p = payload["arms"][arm]
        print(
            f"  {arm}: {p['answerable_in_its_own_cell']} answerable of {p['measured_links']} measured, "
            f"{p['independent_loci_of_the_answerable_links']} independent loci, "
            f"causes reconcile {p['causes_sum_to_the_population']}"
        )
        for cause, n in p["causes"].items():
            print(f"      {n:>4}  {cause}")
    v = payload["verdict"]
    print(
        f"  verdict: answerable {v['answerable_links']} of {v['population_measured_links']} measured, "
        f"{v['independent_loci_of_the_answerable_links']} loci; floors "
        f"{rt.POSITIVE_FLOOR} links / {rt.LOCUS_FLOOR} loci -> "
        f"{'both floors met' if v['meets_both_floors'] else 'NO-GO'}"
    )
    for s in v["short_by"]:
        print(f"      {s}")
    print(f"  statistic: {payload['statistic']['implication']}")


if __name__ == "__main__":
    main()
