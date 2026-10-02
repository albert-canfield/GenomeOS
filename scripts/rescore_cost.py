# SPDX-License-Identifier: AGPL-3.0-or-later
"""What a re-score would cost to make the increase-derived direction test answerable. 0 requests.

    uv run --frozen python scripts/rescore_cost.py
    uv run --frozen python scripts/rescore_cost.py --chroms chr21 --result rescore_cost_chr21

The method, the population rules, the request unit, the window facts, the purchasability rule and the
ceiling rule are `genomeos.attribution.rescore_cost`, written and committed before this costing. In the
registered order:

1. **The extractor's own links**, both directions, re-enumerated genome-wide from
   `measured.rule_links_detail(row, measured.EXTRACTOR_V2)` over the rows `measured.rows` builds, blind
   to every cached value. The population is re-enumerated rather than read out of lane-reptest's result
   because that result dumps only the increase arm, and this costing must cost the decrease arm too.
   lane-reptest's figures are read as well and asserted to agree, so neither count can drift.
2. **The two populations, each by its own rule**, neither stated in terms of the no-go:
   the CELL population is every link whose own cell is not one of the tracks the sweep kept, computed
   against `CACHED_CELLS` and not read from `CELLS_NOT_CACHED`; the WINDOW population is every link
   whose tested gene's transcription start lies outside the window the sweep scored at its element.
3. **Each link's tested gene's TSS**, from the benchmark pairs behind that link, found by the committed
   overlap rule. A link with several tested intervals keeps them all and is placed inside the window if
   ANY of them puts the TSS there, so no pair is chosen.
4. **The window each link would need**, against the four sequence lengths the client accepts.
5. **The purchasability of each**, exhaustive and disjoint, and the request count by DISTINCT element.
6. **The ceiling** of the purchase against both imported floors, the two arms reported apart.
7. **The precedents**, read from the results that record them rather than typed in.

**Nothing is bought and nothing is sent.** The per-element cache is read for the population's elements,
one gene per link, exactly as lane-reptest read it: the per-chromosome archive is streamed with only a
wanted record materialised, and an element the archive does not hold is looked up as a loose per-element
file. The second, partial cache under `data/knowledge/alphagenome/elements_hct116` is OPENED HERE for one
purpose only: to count which purchasable elements already carry a delivered answer, so the costing can
say what the project has already paid for. No value from it enters any finding, no sign is read from it,
and the registered predicate still refuses it - a re-score under that predicate must write into the
sweep's own cache, which is why those elements are costed at one request each all the same.

**No direction is read.** Where a delivered answer is reported it is reported as "carries a value that
is not exactly zero" - the answerability predicate's own question - and never as a value, a sign or a
comparison. R2 is binding: `increase_links.check_no_mechanism_claim` is called on every link record.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
import time
from collections import Counter
from datetime import date
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import crispri, mpra, vista  # noqa: E402
from genomeos.attribution import direction_link as dl  # noqa: E402
from genomeos.attribution import increase_links as il  # noqa: E402
from genomeos.attribution import measured as ms  # noqa: E402
from genomeos.attribution import reptest as rt  # noqa: E402
from genomeos.attribution import rescore_cost as rc  # noqa: E402
from genomeos.attribution.targets import attributed  # noqa: E402
from genomeos.knowledge import satmut as satmut_knowledge  # noqa: E402
from genomeos.results import RESULTS_DIR, save_result  # noqa: E402

RESULT = "rescore_cost"
CHROMS = tuple(f"chr{c}" for c in [*range(1, 23), "X", "Y"])
CACHE_DIR = Path("data/knowledge/alphagenome/elements")
SECOND_CACHE = Path("data/knowledge/alphagenome/elements_hct116")
INHERITED = Path("data/results/reptest_answerable.json")
PRECEDENT_CRISPRI = Path("data/results/crispri_published_v2.json")
PRECEDENT_ASTROREG = Path("data/results/astroreg_request_plan.json")

#: The same per-chromosome results `targets.attributed` reads, each digested by its own path. The
#: `enhancer_targets_all` summaries are POINTERS: the table they name under `elements_where` is what the
#: run really opens, so `pointer_tables` collects those paths and declares them too.
INPUT_GLOBS = (
    "enhancer_targets_chr*.json",
    "enhancer_targets_all_chr*.json",
    "constrained_targets_chr*.json",
)

#: This lane's own files. Everything else uncommitted in this shared checkout belongs to another lane.
OWN_CODE = (
    "genomeos/attribution/rescore_cost.py",
    "scripts/rescore_cost.py",
    "tests/test_rescore_cost.py",
)


def streaming_reader() -> Any:
    """lane-reptest's own archive reader, loaded from its committed runner rather than copied.

    The reader is `scripts/reptest_answerable.py`'s `cached_records`, which is
    `scripts/direction_link.py`'s carried unchanged. Loading it keeps one implementation of how the
    archive is read, so this lane and that one cannot come to read the cache differently.
    """
    path = Path(__file__).resolve().parents[1] / "scripts" / "reptest_answerable.py"
    spec = importlib.util.spec_from_file_location("reptest_answerable_reader", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.cached_records


def links_of(chrom: str) -> list[dict[str, Any]]:
    """The extractor's own links on one chromosome under v2, both directions. Blind to every value.

    `measured.rule_links_detail` is called exactly as committed and nothing is re-derived from the
    pairs. The row supplies the element's interval, which the extractor's record does not carry and both
    the window test and the locus grouping need. Carried from `scripts/reptest_answerable.py`, so the
    two lanes enumerate one population.
    """
    layer = ms.Layer.load(chrom)
    elements = attributed(chrom, RESULTS_DIR)
    links: list[dict[str, Any]] = []
    for row in ms.rows(chrom, elements, layer):
        if not row["measured"].get("crispri"):
            continue
        for r in ms.rule_links_detail(row, rt.EXTRACTOR):
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
    return links


def benchmark_pairs() -> dict[tuple[str, str, str], list[crispri.Pair]]:
    """Every valid benchmark pair, keyed by (chromosome, gene, cell). The tables are read as committed."""
    out: dict[tuple[str, str, str], list[crispri.Pair]] = {}
    for name in ms.CRISPRI_FILES:
        for p in crispri.load(name, crispri.KNOWLEDGE):
            out.setdefault((p.chrom, p.gene, p.cell), []).append(p)
    return out


def pairs_behind(
    link: dict[str, Any], pairs: dict[tuple[str, str, str], list[crispri.Pair]]
) -> list[crispri.Pair]:
    """The benchmark pairs behind a link, by the committed overlap rule. All of them, none chosen.

    The link's gene and cell are its own and the element interval is the row's, so the candidates are
    the pairs of that (chromosome, gene, cell) whose tested interval reciprocally overlaps the element at
    `measured.RECIPROCAL_OVERLAP`, which is the rule the extractor's own layer used to raise the link.
    Several pairs are kept rather than resolved: a link is inside the window if ANY of its pairs puts the
    TSS inside, and the width it needs is the smallest any of them needs, so no pair is picked.
    """
    return [
        p
        for p in pairs.get((link["chrom"], link["gene"], link["cell"]), [])
        if ms.measures(link["start"], link["end"], p.start, p.end) and p.tss is not None
    ]


def records_for(
    links: list[dict[str, Any]], read_archive: Any
) -> tuple[dict[str, dict[str, Any]], list[str], list[str]]:
    """The sweep's own record for each link's element, with the files opened named for the manifest."""
    by_chrom: dict[str, set[str]] = {}
    for r in links:
        by_chrom.setdefault(r["chrom"], set()).add(r["element"])
    records: dict[str, dict[str, Any]] = {}
    archives: list[str] = []
    loose: list[str] = []
    for chrom in sorted(by_chrom):
        wanted = by_chrom[chrom]
        path = CACHE_DIR / f"{chrom}.json.gz"
        got = read_archive(path, wanted)
        if wanted and path.exists():
            archives.append(path.as_posix())
        records.update(got)
        for element in sorted(wanted - set(got)):
            p = CACHE_DIR / chrom / f"{element}.json"
            if p.exists():
                records[element] = json.loads(p.read_text())
                loose.append(p.as_posix())
    return records, sorted(archives), sorted(loose)


def delivered_elsewhere(chrom: str, element: str, gene: str, cell: str) -> dict[str, Any]:
    """Whether the partial second cache already holds a delivered answer for this link's own cell.

    One purpose: to say what the project has already paid for. The answer is a presence and a NOT-ZERO
    fact, which is the answerability predicate's own question; no value, no sign and no comparison is
    taken from it, and the registered predicate still refuses this cache, so the link is costed at one
    request all the same.
    """
    p = SECOND_CACHE / chrom / f"{element}.json"
    if not p.exists():
        return {"second_cache_holds_this_element": False, "path": p.as_posix(), "opened": False}
    row = json.loads(p.read_text())
    entry = next((g for g in row.get("genes") or [] if g.get("gene") == gene), None)
    value = (entry.get("by_cell") or {}).get(cell) if entry else None
    return {
        "second_cache_holds_this_element": True,
        "path": p.as_posix(),
        "opened": True,
        "names_this_gene": entry is not None,
        "carries_this_cells_track": value is not None,
        "value_is_not_exactly_zero": None if value is None else float(value) != 0.0,
        "no_value_is_taken_from_it": (
            "the registered predicate names the sweep's own per-element cache and refuses this one, so "
            "this link is costed at one request regardless. What is recorded here is that the request "
            "has been delivered once already, which is a fact about money spent and not a value"
        ),
    }


def geometry(
    links: list[dict[str, Any]],
    records: dict[str, dict[str, Any]],
    pairs: dict[tuple[str, str, str], list[crispri.Pair]],
) -> list[dict[str, Any]]:
    """One row per link: its scored window, its tested gene's TSS, and what it would need. No value read."""
    out: list[dict[str, Any]] = []
    for r in links:
        rec = records.get(r["element"])
        behind = pairs_behind(r, pairs)
        row: dict[str, Any] = {
            "chrom": r["chrom"],
            "element": r["element"],
            "gene": r["gene"],
            "cell": r["cell"],
            "split": r["split"],
            "derived_from": r["derived_from"],
            "pairs_behind_the_link_with_a_tss": len(behind),
            "cell_track_is_one_the_sweep_kept": r["cell"] in rt.CACHED_CELLS,
            "element_record_in_the_sweeps_cache": rec is not None,
        }
        if rec is None or not behind:
            row["window_width_needed"] = None
            row["smallest_supported_window_that_would_contain_the_tss"] = None
            row["the_scored_window_reaches_the_tss"] = None
            row["gene_named_in_the_scored_record"] = (
                None if rec is None else r["gene"] in {g.get("gene") for g in rec.get("genes") or []}
            )
            row["why_no_geometry"] = (
                "the sweep's cache holds no record of this element, so the window it was scored at is "
                "not on disk"
                if rec is None
                else "no benchmark pair behind this link carries a startTSS, so the tested gene's "
                "position is not in the table"
            )
            out.append(row)
            continue
        lo, hi = rc.window_bounds(rec["start"], rec["end"], rc.WINDOW_USED_BY_THE_SWEEP)
        tss = sorted({int(p.tss) for p in behind})  # type: ignore[arg-type]
        centre = rc.window_centre(rec["start"], rec["end"])
        needed = min(rc.window_needed(rec["start"], rec["end"], t) for t in tss)
        names = {g.get("gene") for g in rec.get("genes") or []}
        row.update(
            {
                "scored_interval": [rec["start"], rec["end"]],
                "window_scored": [lo, hi],
                "genes_in_window_of_the_scored_record": rec.get("genes_in_window"),
                "tss": tss,
                "offsets_of_tss_from_the_window_centre": [t - centre for t in tss],
                "benchmark_distance_to_tss": sorted({p.distance for p in behind}),
                "window_width_needed": needed,
                "smallest_supported_window_that_would_contain_the_tss": rc.smallest_supported_window(needed),
                "the_scored_window_reaches_the_tss": any(lo <= t < hi for t in tss),
                "gene_named_in_the_scored_record": r["gene"] in names,
            }
        )
        out.append(row)
    return out


def cost_cell_population(rows: list[dict[str, Any]]) -> dict[str, Any]:
    """The CELL population costed: every link in a cell the sweep kept no track of, both directions."""
    mine = [r for r in rows if not r["cell_track_is_one_the_sweep_kept"]]
    for r in mine:
        r["purchasability"] = rc.purchasability(
            bool(r["gene_named_in_the_scored_record"]),
            r["window_width_needed"] if r["window_width_needed"] is not None else 1 << 62,
            False,
        )
    counts = Counter(r["purchasability"] for r in mine)
    elements = {r["element"] for r in mine}
    purchasable_elements = {r["element"] for r in mine if r["purchasability"] == rc.PURCHASABLE_CELL_TRACK}
    by_direction: dict[str, dict[str, Any]] = {}
    for arm, label in ((ms.INCREASE, "increase_derived"), (ms.DECREASE, "decrease_derived")):
        side = [r for r in mine if r["derived_from"] == arm]
        by_direction[label] = {
            "links": len(side),
            "by_cell": dict(sorted(Counter(r["cell"] for r in side).items())),
            "purchasability": dict(sorted(Counter(r["purchasability"] for r in side).items())),
            "purchasable_links": sum(1 for r in side if r["purchasability"] == rc.PURCHASABLE_CELL_TRACK),
            "elements": len({r["element"] for r in side}),
            "requests_if_this_direction_were_bought_alone": rc.requests_for([r["element"] for r in side]),
        }
    return {
        "population_rule": rc.POPULATION_RULE_CELL_ARM,
        "cells_not_covered_is_computed_not_inherited": rc.CELLS_NOT_COVERED_IS_COMPUTED_NOT_INHERITED,
        "cells_not_covered": rc.cells_not_covered([r["cell"] for r in rows]),
        "cells_the_sweep_kept": list(rt.CACHED_CELLS),
        "cells_not_cached_constant_names_only": list(dl.CELLS_NOT_CACHED),
        "links": len(mine),
        "elements": len(elements),
        "requests_under_the_unit": rc.requests_for([r["element"] for r in mine]),
        "requests_serving_no_link_the_test_could_score": len(elements - purchasable_elements),
        "purchasability": {k: counts[k] for k in rc.PURCHASABILITY if counts[k]},
        "purchasability_sums_to_the_population": sum(counts.values()) == len(mine),
        "by_direction": by_direction,
        "directions_apart": (
            "the two directions are reported apart because they are different propositions: the floors "
            "that refused the test are the increase arm's, so a request that makes a decrease link "
            "answerable costs the same and moves neither of them. The population rule admits both and "
            "the purchase would buy both; the costing refuses to pool them into one yield"
        ),
        "what_this_is_not": (
            "not a coverage figure and not an answerable count. " + rt.COUNT_IS_NOT_COVERAGE
        ),
        "links_detail": mine,
    }


def cost_window_population(rows: list[dict[str, Any]]) -> dict[str, Any]:
    """The WINDOW population counted: every link whose tested gene the scored window does not reach."""
    mine = [r for r in rows if r["the_scored_window_reaches_the_tss"] is False]
    widths = sorted(r["window_width_needed"] for r in mine if r["window_width_needed"] is not None)
    over = [w for w in widths if w > rc.MAX_SUPPORTED_WINDOW]
    return {
        "population_rule": rc.POPULATION_RULE_WINDOW_ARM,
        "closed_on_a_service_limit": rc.WINDOW_ARM_IS_CLOSED_ON_A_SERVICE_LIMIT,
        "links": len(mine),
        "elements": len({r["element"] for r in mine}),
        "requests_available": 0,
        "cost": "NO AVAILABLE REQUEST: the sweep already asks at the largest width the client accepts, "
        f"{rc.MAX_SUPPORTED_WINDOW:,} bases, so no wider request exists to send at any price",
        "by_direction": {
            label: {
                "links": sum(1 for r in mine if r["derived_from"] == arm),
                "widths_above_the_service_maximum": sum(
                    1
                    for r in mine
                    if r["derived_from"] == arm and (r["window_width_needed"] or 0) > rc.MAX_SUPPORTED_WINDOW
                ),
                "widths_at_or_below_the_service_maximum": sum(
                    1
                    for r in mine
                    if r["derived_from"] == arm
                    and r["window_width_needed"] is not None
                    and r["window_width_needed"] <= rc.MAX_SUPPORTED_WINDOW
                ),
            }
            for arm, label in ((ms.INCREASE, "increase_derived"), (ms.DECREASE, "decrease_derived"))
        },
        "window_width_needed": {
            "n": len(widths),
            "min": widths[0] if widths else None,
            "median": widths[len(widths) // 2] if widths else None,
            "max": widths[-1] if widths else None,
            "above_the_service_maximum": len(over),
            "at_or_below_the_service_maximum": len(widths) - len(over),
            "all_widths": widths,
            "one_width_does_not_fit_all": (
                "the width is computed per link from that link's own element-to-TSS offset and the "
                "widths are listed in full, because a single figure would hide that most are beyond "
                "the service maximum and some are not"
            ),
        },
        "a_width_at_or_below_the_maximum_is_not_a_purchase": (
            "a link the scored window did not reach whose TSS nonetheless lies within the largest "
            "accepted width is not a link a request would fix: the window already covered the place, so "
            "the gap is not the window. Those links are counted apart and are not costed as purchases"
        ),
        "links_detail": mine,
    }


def naming_gap_population(rows: list[dict[str, Any]]) -> dict[str, Any]:
    """The links the scored window DID reach and whose record names no gene by the benchmark's symbol.

    A third population, counted because it costs nothing and because it is where a purchase would be
    wasted. The window already covered the place, so no request changes what the record names; the gap is
    between the benchmark's gene symbol and the symbol the scorer emitted. Whether the registered
    predicate may match a gene under a second symbol is a question about the predicate and not a cost,
    and this lane does not amend it, does not match an alias and counts none of these as answerable.
    """
    mine = [
        r
        for r in rows
        if r["the_scored_window_reaches_the_tss"] is True and not r["gene_named_in_the_scored_record"]
    ]
    return {
        "population_rule": (
            "every measured significant link the extractor emits under v2 whose tested gene's "
            "transcription start lies INSIDE the window the sweep scored at that element and whose "
            "cached record names no gene by the benchmark's symbol. The rule names no count and no "
            "outcome"
        ),
        "links": len(mine),
        "elements": len({r["element"] for r in mine}),
        "requests_available": 0,
        "cost": "NO REQUEST WOULD HELP: the window already reached the place, so the gap is not the "
        "window and not the cell. It is a symbol that does not match, and a request returns the same "
        "symbols again",
        "by_direction": {
            label: sum(1 for r in mine if r["derived_from"] == arm)
            for arm, label in ((ms.INCREASE, "increase_derived"), (ms.DECREASE, "decrease_derived"))
        },
        "by_cell": dict(sorted(Counter(r["cell"] for r in mine).items())),
        "it_is_not_an_answerable_link": (
            "none of these is counted as answerable here. The registered predicate asks for a gene entry "
            "whose name is exactly the link's gene, and this lane neither amends that predicate nor "
            "matches an alias; it reports the count and the symbols the record does name so that the "
            "question can be put to whoever owns the predicate"
        ),
        "links_detail": [
            {
                "chrom": r["chrom"],
                "element": r["element"],
                "gene": r["gene"],
                "cell": r["cell"],
                "derived_from": r["derived_from"],
                "offsets_of_tss_from_the_window_centre": r["offsets_of_tss_from_the_window_centre"],
                "genes_in_window_of_the_scored_record": r["genes_in_window_of_the_scored_record"],
            }
            for r in mine
        ],
    }


def pointer_tables(chroms: list[str]) -> list[Path]:
    """The tables the `enhancer_targets_all` summaries point at through `elements_where`.

    A declared input must be the bytes the run really opened. These summaries are pointers: each names a
    table under `data/knowledge/alphagenome/all_elements`, and `targets.run_elements` opens that table.
    `save_result` refuses a result whose traced reads exceed its declared inputs, and it exists because
    the 24 summaries declare about 39 KB of pointers and read about 1.4 GB of tables.
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
        raise SystemExit(f"{RESULT} is the genome-wide costing; a run over {len(chroms)} needs --result")

    started = time.time()
    read_archive = streaming_reader()
    links: list[dict[str, Any]] = []
    for chrom in chroms:
        links.extend(links_of(chrom))
    increases = [r for r in links if r["derived_from"] == ms.INCREASE]
    decreases = [r for r in links if r["derived_from"] == ms.DECREASE]

    inherited = json.loads(INHERITED.read_text())
    arm = inherited["arms"]["increases"]
    inherited_figures = {
        "increase_derived_links": inherited["population"]["increase_derived_links"],
        "answerable_in_its_own_cell": arm["answerable_in_its_own_cell"],
        "independent_loci_of_the_answerable_links": arm["independent_loci_of_the_answerable_links"],
        **{c: arm["causes"].get(c, 0) for c in rt.CAUSES if c != rt.ANSWERABLE},
    }
    disagree = {
        k: (v, rc.INHERITED_FIGURES[k])
        for k, v in inherited_figures.items()
        if k in rc.INHERITED_FIGURES and v != rc.INHERITED_FIGURES[k]
    }
    if args.result == RESULT and (disagree or len(increases) != inherited_figures["increase_derived_links"]):
        raise SystemExit(
            "the result this costing is taken against does not hold the figures registered for it: "
            f"{disagree or (len(increases), inherited_figures['increase_derived_links'])}. Nothing is "
            "costed against a denominator that moved"
        )
    answerable_increases = [r for r in inherited["increase_links"] if r["answerability"] == rt.ANSWERABLE]

    records, archives, loose = records_for(links, read_archive)
    pairs = benchmark_pairs()
    rows = geometry(links, records, pairs)
    cell_arm = cost_cell_population(rows)
    window_arm = cost_window_population(rows)
    naming_arm = naming_gap_population(rows)

    purchasable = [r for r in cell_arm["links_detail"] if r["purchasability"] == rc.PURCHASABLE_CELL_TRACK]
    for r in purchasable:
        r["already_delivered_once"] = delivered_elsewhere(r["chrom"], r["element"], r["gene"], r["cell"])
    purchasable_keys = {(r["element"], r["gene"], r["cell"]) for r in purchasable}
    up = [r for r in purchasable if r["derived_from"] == ms.INCREASE]
    down = [r for r in purchasable if r["derived_from"] == ms.DECREASE]
    hypothetical = [r for r in increases if (r["element"], r["gene"], r["cell"]) in purchasable_keys]
    requests = cell_arm["requests_under_the_unit"]
    loci_now = len(set(dl.loci_of(answerable_increases)))
    loci_ceiling = len(set(dl.loci_of(answerable_increases + hypothetical)))
    verdict = rc.ceiling(len(answerable_increases), len(up), loci_ceiling)
    naming_increases = [
        r
        for r in increases
        if (r["element"], r["gene"], r["cell"])
        in {(n["element"], n["gene"], n["cell"]) for n in naming_arm["links_detail"]}
    ]
    generous_links = hypothetical + naming_increases
    generous = rc.ceiling(
        len(answerable_increases),
        len(up) + len(naming_increases),
        len(set(dl.loci_of(answerable_increases + generous_links))),
    )
    generous["rule"] = rc.MOST_GENEROUS_CEILING_RULE
    generous["naming_gap_increase_links_assumed_closed_at_no_cost"] = len(naming_increases)
    all_increase_loci = len(set(dl.loci_of(increases)))

    crispri_precedent = json.loads(PRECEDENT_CRISPRI.read_text())
    astroreg_precedent = json.loads(PRECEDENT_ASTROREG.read_text())
    hct = crispri_precedent["second_cell_type_hct116"]
    hct_requests = hct["requests"]["delivered_in_the_carried_run"]["sent"]
    hct_pairs = hct["as_carried"]["covered_pairs"]
    wtc11_cost = crispri_precedent["second_cell_type"]["cost"]["WTC11"]["requests_frozen_feature"]
    astro = astroreg_precedent["summary"]
    precedents = {
        "HCT116 second cell type, delivered": {
            "requests": hct_requests,
            "covered_pairs": hct_pairs,
            "requests_per_covered_pair": round(hct_requests / hct_pairs, 4),
            "read_from": rc.PRECEDENT_SOURCES["HCT116 second cell type, delivered"],
        },
        "WTC11 second cell type, costed and not bought": {
            "requests": wtc11_cost,
            "read_from": rc.PRECEDENT_SOURCES["WTC11 second cell type, costed and not bought"],
        },
        "AstroREG request plan, authorised and not sent": {
            "requests": astro["requests"],
            "requests_serving_at_least_one_scored_pair": astro["requests_serving_at_least_one_scored_pair"],
            "requests_serving_no_scored_pair": astro["requests_serving_no_scored_pair"],
            "read_from": rc.PRECEDENT_SOURCES["AstroREG request plan, authorised and not sent"],
        },
        "this costing, beside them": {
            "requests": requests,
            "answerable_positives_at_the_ceiling_increase_arm": len(up),
            "answerable_positives_at_the_ceiling_decrease_arm": len(down),
            "requests_serving_no_link_the_test_could_score": cell_arm[
                "requests_serving_no_link_the_test_could_score"
            ],
            "share_of_the_WTC11_arm_costed_but_not_bought": round(requests / wtc11_cost, 4),
            "what_it_would_buy": (
                "a larger answerable count on the increase arm that is still below an imported floor, "
                "so no testable number: see the ceiling"
                if not verdict["meets_both_floors_at_the_ceiling"]
                else "an increase arm whose ceiling clears both imported floors"
            ),
        },
        "precedents_are_read_not_quoted": rc.PRECEDENTS_ARE_READ_NOT_QUOTED,
    }

    already = [r for r in purchasable if r["already_delivered_once"].get("carries_this_cells_track")]
    payload: dict[str, Any] = {
        "result": args.result,
        "date": date.today().isoformat(),
        "lane": rc.LANE,
        "question": (
            "what it would cost, in model requests, to make the direction test on extractor v2's "
            "increase-derived links answerable, counted at zero requests over a population stated "
            "without reference to the no-go, with the two directions apart and the ceiling of the "
            "purchase against both imported floors"
        ),
        "registered_in_code_before_this_costing": "genomeos/attribution/rescore_cost.py",
        "registration": rc.registration(),
        "headline": verdict["reading"],
        "follows_a_no_go_whose_outcome_was_never_read": rc.FOLLOWS_A_NO_GO_WHOSE_OUTCOME_WAS_NEVER_READ,
        "the_shortfall_may_not_be_bought": rc.THE_SHORTFALL_MAY_NOT_BE_BOUGHT,
        "binding_wording_r2": rc.R2_WORDING,
        "no_direction_is_read": rc.NO_DIRECTION_IS_READ,
        "chromosomes": chroms,
        "population": {
            "call": rt.POPULATION_IS_THE_EXTRACTORS_OUTPUT,
            "extractor": rt.EXTRACTOR,
            "v2_links": len(links),
            "increase_derived_links": len(increases),
            "decrease_derived_links": len(decreases),
            "by_construction": rt.POPULATION_IS_NOT_COVERAGE,
            "count_is_not_coverage": rt.COUNT_IS_NOT_COVERAGE,
        },
        "inherited": {
            "lane": "lane-reptest",
            "result": rc.INHERITED_RESULT,
            "figures_as_read": inherited_figures,
            "figures_as_registered": dict(rc.INHERITED_FIGURES),
            "figures_agree": not disagree,
            "increase_links_re_enumerated_here": len(increases),
            "no_go_carried_verbatim": rc.INHERITED_NO_GO,
            "answerable_increase_links_now": len(answerable_increases),
            "independent_loci_now": loci_now,
        },
        "request_unit": {
            "unit": rc.REQUEST_UNIT,
            "read_from": dict(rc.REQUEST_UNIT_READ_FROM),
            "checked_not_carried_over": rc.REQUEST_UNIT_IS_CHECKED_NOT_CARRIED_OVER,
            "counted_over": "distinct element ids, never links",
        },
        "window": {
            "used_by_the_sweep": rc.WINDOW_USED_BY_THE_SWEEP,
            "supported_sequence_lengths": list(rc.SUPPORTED_SEQUENCE_LENGTHS),
            "max_supported_window": rc.MAX_SUPPORTED_WINDOW,
            "the_sweep_already_asks_at_the_maximum": rc.WINDOW_USED_BY_THE_SWEEP == rc.MAX_SUPPORTED_WINDOW,
            "is_at_the_maximum": rc.WINDOW_IS_AT_THE_MAXIMUM,
            "centring": rc.WINDOW_CENTRING,
            "tss_source": rc.TSS_IS_THE_BENCHMARK_S,
            "necessary_not_sufficient": rc.WINDOW_TEST_IS_NECESSARY_NOT_SUFFICIENT,
        },
        "cell_population": cell_arm,
        "window_population": window_arm,
        "naming_gap_population": naming_arm,
        "populations_overlap": {
            "links_in_both": sum(
                1
                for r in rows
                if not r["cell_track_is_one_the_sweep_kept"]
                and r["the_scored_window_reaches_the_tss"] is False
            ),
            "what_it_means": (
                "a link in an uncovered cell whose tested gene the scored window also fails to reach is "
                "in both populations and is purchasable in neither: the cell request would be sent and "
                "could not place that gene in the window. It is counted in both and costed as "
                "unpurchasable"
            ),
        },
        "cost": {
            "requests_total": requests,
            "what_the_requests_are_for": (
                "every element carrying a link of the CELL population, whichever direction that link's "
                "measurement had and whether or not the link turns out purchasable. The population rule "
                "is applied as stated and the requests it implies are counted in full, because buying "
                "only the purchasable subset would be buying the gap"
            ),
            "per_request_and_per_answerable_positive": rc.per_request_figures(
                requests,
                len(up),
                len(down),
                loci_ceiling - loci_now,
                cell_arm["requests_serving_no_link_the_test_could_score"],
            ),
            "already_delivered_once_in_the_second_cache": {
                "links": len(already),
                "elements": sorted({r["element"] for r in already}),
                "what_it_means": (
                    "the request has been delivered once already, into a cache the registered "
                    "predicate refuses. The link is costed at one request all the same, because a "
                    "re-score under that predicate must write into the sweep's own cache. This is "
                    "money the project has already spent, not an answer it may use"
                ),
            },
            "window_population_requests": 0,
            "window_population_cost": window_arm["cost"],
        },
        "ceiling": verdict,
        "most_generous_ceiling": generous,
        "if_every_currently_unanswerable_increase_link_became_answerable": {
            "links": len(increases),
            "independent_loci": all_increase_loci,
            "clears_the_links_floor": len(increases) >= rc.POSITIVE_FLOOR,
            "clears_the_loci_floor": all_increase_loci >= rc.LOCUS_FLOOR,
            "margin": [
                f"links {len(increases)} against the imported floor {rc.POSITIVE_FLOOR}: "
                f"{len(increases) - rc.POSITIVE_FLOOR:+d}",
                f"independent loci {all_increase_loci} against the imported floor {rc.LOCUS_FLOOR}: "
                f"{all_increase_loci - rc.LOCUS_FLOOR:+d}",
            ],
            "is_it_reachable": (
                "no. The hypothetical clears both floors and is recorded for completeness, and it "
                "cannot be reached: the links the widest accepted window cannot reach are unpurchasable "
                "at any price and a link the window already reached is not a purchase. The CEILING is "
                "the figure that bears on a decision"
            ),
        },
        "what_cannot_be_established_without_a_request": {
            "figure": "whether a purchasable element's request would return a value that is not "
            "exactly zero on that link's own cell track",
            "why": (
                "the sweep kept no value for that cell at that element, so there is nothing on disk to "
                "read. The ceiling therefore assumes every purchasable request returns a non-zero "
                "value, which is why it is a ceiling and not a count. It does not change the headline: "
                "the ceiling is an upper bound and the gate is taken on it"
            ),
            "it_is_named_and_not_estimated": True,
        },
        "precedents": precedents,
        "the_registration_that_would_spend_this_is_a_new_one": (
            "nothing here registers a purchase or a test. A purchase would need its own registration, "
            "with its population rule, its floors and its readings fixed BEFORE any request is sent, "
            "and this lane does not write it. This lane costs it, and the decision is the owner's"
        ),
    }

    paths: set[Path] = set()
    for p in (INHERITED, PRECEDENT_CRISPRI, PRECEDENT_ASTROREG):
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
    cache_files = archives + loose
    if cache_files:
        inputs.append(
            mf.files_entry(
                "per_element_response_cache",
                cache_files,
                partition="the finished deletion sweep's own per-element cache, read for the scored "
                "interval and which genes each record names, and for nothing else",
            )
        )
    second = sorted(
        p
        for p in {
            r["already_delivered_once"]["path"]
            for r in purchasable
            if r["already_delivered_once"].get("opened")
        }
        if p
    )
    if second:
        inputs.append(
            mf.files_entry(
                "second_partial_cache_opened_to_count_delivered_requests",
                second,
                partition="the partial second cache, opened ONLY to count which purchasable elements "
                "already carry a delivered answer. No value, sign or comparison is taken from it and "
                "the registered predicate still refuses it",
            )
        )

    payload["seconds"] = round(time.time() - started, 1)
    payload["result_manifest"] = {
        "sources": [
            {
                "accession": "ENCODE CRISPRi enhancer-gene benchmark "
                "(EngreitzLab/CRISPR_comparison, Gschwind et al.)",
                "version": "the two benchmark tables cached under data/knowledge, digested in inputs; "
                "the startTSS and distanceToTSS columns are read and no annotation is joined in",
            },
            {
                "accession": "the finished genome-wide AlphaGenome deletion sweep's per-element "
                "response cache, written with threshold=0.0",
                "version": "the per-chromosome archives and loose per-element files this run opened, "
                "digested as a group; read for the scored interval and the gene names only",
            },
            {
                "accession": "the partial second per-element cache under "
                "data/knowledge/alphagenome/elements_hct116",
                "version": "the files this run opened, digested as a group; opened to count delivered "
                "requests and for nothing else",
            },
            {
                "accession": "the attributed elements of the chromosomes this run covers",
                "version": "re-enumerated in this run from the results on disk and the tables their "
                "`elements_where` fields point at, not read from one file",
            },
            {
                "accession": "alphagenome 0.9.0, the version uv.lock pins (wheel sha256 "
                "a4f35884341ae85b5d2cf088dfe0304961de7ae4a590d4538653069673de32f4)",
                "version": "`alphagenome/models/dna_client.py` SUPPORTED_SEQUENCE_LENGTHS, recorded as "
                "a committed literal and cited; the package is an optional extra and is not imported "
                "by this run",
            },
        ],
        "inputs": inputs,
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "what_this_is": rc.WHAT_THIS_IS,
            "follows_a_no_go_whose_outcome_was_never_read": (rc.FOLLOWS_A_NO_GO_WHOSE_OUTCOME_WAS_NEVER_READ),
            "the_shortfall_may_not_be_bought": rc.THE_SHORTFALL_MAY_NOT_BE_BOUGHT,
            "population_rule_cell_arm": rc.POPULATION_RULE_CELL_ARM,
            "population_rule_window_arm": rc.POPULATION_RULE_WINDOW_ARM,
            "cells_not_covered_is_computed_not_inherited": (rc.CELLS_NOT_COVERED_IS_COMPUTED_NOT_INHERITED),
            "window_arm_is_closed_on_a_service_limit": rc.WINDOW_ARM_IS_CLOSED_ON_A_SERVICE_LIMIT,
            "extractor": rt.EXTRACTOR,
            "request_unit": rc.REQUEST_UNIT,
            "request_unit_read_from": dict(rc.REQUEST_UNIT_READ_FROM),
            "window_used_by_the_sweep": rc.WINDOW_USED_BY_THE_SWEEP,
            "supported_sequence_lengths": list(rc.SUPPORTED_SEQUENCE_LENGTHS),
            "max_supported_window": rc.MAX_SUPPORTED_WINDOW,
            "window_is_at_the_maximum": rc.WINDOW_IS_AT_THE_MAXIMUM,
            "window_centring": rc.WINDOW_CENTRING,
            "purchasability": list(rc.PURCHASABILITY),
            "purchasability_is_exhaustive": rc.PURCHASABILITY_IS_EXHAUSTIVE,
            "ceiling_rule": rc.CEILING_RULE,
            "ceiling_is_not_a_count": rc.CEILING_IS_NOT_A_COUNT,
            "links_floor": rc.POSITIVE_FLOOR,
            "independent_loci_floor": rc.LOCUS_FLOOR,
            "floors_are_imported": rc.FLOORS_ARE_IMPORTED,
            "independent_locus_rule": rt.LOCUS_RULE,
            "independent_locus_span": rt.LOCUS_SPAN,
            "not_biological_independence": rt.NOT_BIOLOGICAL_INDEPENDENCE,
            "reciprocal_overlap": ms.RECIPROCAL_OVERLAP,
            "reach": ms.REACH,
            "chromosomes": chroms,
            "aggregate_function": "a count per purchasability cause and a request count per distinct "
            "element; no rate, no mean, no interval and no direction is taken",
            "fill_value": "none: a link whose cost cannot be established is named and not estimated",
            "model_requests": 0,
            "alphagenome_requests": 0,
        },
        "exclusions": [
            rc.NO_DIRECTION_IS_READ,
            rc.THE_SHORTFALL_MAY_NOT_BE_BOUGHT,
            rc.WINDOW_TEST_IS_NECESSARY_NOT_SUFFICIENT,
            rc.CEILING_IS_NOT_A_COUNT,
            rt.CAUSE_TEXT[rt.CELL_NOT_CACHED],
            "no floor is moved, recomputed, reinterpreted or pooled around",
            "no value, sign or comparison is taken from the partial second cache",
            "no pair is excluded for its magnitude: no effect-size bar is applied anywhere",
        ],
        "partitions": {
            "increase_derived": "links the extractor derived from a significant increase under v2; "
            "action `inhibits`, outcome `increase on knockdown`, molecular role unresolved",
            "decrease_derived": "links the extractor derived from a benchmark `Regulated` pair, which "
            "is a significant decrease only; v1's own links, unchanged under v2",
            "cell_population": rc.POPULATION_RULE_CELL_ARM,
            "window_population": rc.POPULATION_RULE_WINDOW_ARM,
            "naming_gap_population": "every link whose tested gene's TSS the scored window reaches and "
            "whose cached record names no gene by the benchmark's symbol; costed at no request, because "
            "the window already reached the place",
            "training_and_heldout": "the benchmark's own split, one per table, as the extractor "
            "assigns it to a link",
        },
        "code_cleanliness": mf.code_cleanliness(Path(__file__).resolve(), OWN_CODE),
    }

    path = save_result(args.result, payload)
    print(f"{args.result}: {path}")
    print(
        f"  window: the sweep asks at {rc.WINDOW_USED_BY_THE_SWEEP:,}; the service's maximum is "
        f"{rc.MAX_SUPPORTED_WINDOW:,}"
    )
    print(
        f"  population under {rt.EXTRACTOR}: {len(links)} links = {len(decreases)} decrease-derived + "
        f"{len(increases)} increase-derived"
    )
    print(
        f"  CELL population ({', '.join(cell_arm['cells_not_covered'])}): {cell_arm['links']} links on "
        f"{cell_arm['elements']} elements -> {requests} requests"
    )
    for label, block in cell_arm["by_direction"].items():
        print(
            f"      {label}: {block['links']} links, {block['purchasable_links']} purchasable, cells "
            f"{block['by_cell']}"
        )
    print(f"  WINDOW population: {window_arm['links']} links -> {window_arm['cost']}")
    print(
        f"  NAMING-GAP population: {naming_arm['links']} links "
        f"{naming_arm['by_direction']} -> {naming_arm['cost']}"
    )
    c = payload["ceiling"]
    print(
        f"  ceiling on the increase arm: {c['links_at_the_ceiling']} links / "
        f"{c['independent_loci_at_the_ceiling']} loci against floors {rc.POSITIVE_FLOOR}/"
        f"{rc.LOCUS_FLOOR} -> "
        f"{'both met at the ceiling' if c['meets_both_floors_at_the_ceiling'] else 'NO-GO AT THE CEILING'}"
    )
    for s in c["short_by"]:
        print(f"      {s}")
    g = payload["most_generous_ceiling"]
    print(
        f"  most generous ceiling (purchase AND the naming gap closed free): "
        f"{g['links_at_the_ceiling']} links / {g['independent_loci_at_the_ceiling']} loci -> "
        f"{'both met' if g['meets_both_floors_at_the_ceiling'] else 'STILL A NO-GO'}"
    )
    for s in g["short_by"]:
        print(f"      {s}")
    r = payload["cost"]["per_request_and_per_answerable_positive"]
    print(
        f"  requests per answerable positive at the ceiling: increase arm "
        f"{r['requests_per_answerable_positive_at_the_ceiling']['increase_arm_only']}, decrease arm "
        f"{r['requests_per_answerable_positive_at_the_ceiling']['decrease_arm_only']}, pooled "
        f"{r['requests_per_answerable_positive_at_the_ceiling']['both_arms_pooled']}"
    )
    print(
        f"  precedents: HCT116 {hct_requests} requests / {hct_pairs} covered pairs; AstroREG "
        f"{astro['requests']} of which {astro['requests_serving_no_scored_pair']} served no pair"
    )
    print(f"  headline: {payload['headline']}")


if __name__ == "__main__":
    main()
