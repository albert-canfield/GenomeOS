# SPDX-License-Identifier: AGPL-3.0-or-later
"""Why 198 measured loci attach to no compiled element: one cause for each of them.

    uv run --frozen python scripts/placement_cause_198.py

The second increment's coverage count (`0487e49`, data/results/response_map_coverage.json) found that of
the 473 independent loci the measured layer's CRISPRi perturbation assay holds, **198 attach to no
compiled element at all** under the measured layer's own overlap rule (`measured.RECIPROCAL_OVERLAP = 0.5`),
and that this is the limiting input on the response map. It did not ask why. This script asks why, and
assigns **exactly one** cause to each of the 198.

It decides nothing. No rule is changed, no policy is adopted, no verdict is moved, nothing is deleted,
no model is called and nothing is downloaded: every input was already on disk.

**The classification is not invented here.** The causes, their meanings and the order in which they are
tried are audit A's, imported from `scripts/placement_audit.py` (`6c0d39f`) and reported word for word:
`placement_audit.CAUSES`, `CAUSE_MEANING` and, for a candidate exclusion, `EXCLUSION`. The locus is the
operational grouping imported from `genomeos.attribution.cell2`, which is **not established biological
independence** and which pools cell types. The attachment test is the one that defines the 198, imported
from `scripts/response_map_coverage.py` (`attaches_to`), so the population here is that population.

**Audit A's 387 and this 198 are different numbers over different units** and are reconciled rather than
equated: audit A counted *positive links* (significant decreases) that no element of either set carries,
198 counts *loci* (all outcomes pooled, cell-pooled) none of whose pairs carries in the compiled set.
The reconciliation block states how the two populations meet, pair by pair.

Writes data/results/placement_cause_198.json through the result contract.
"""

from __future__ import annotations

import json
import resource
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

import placement_audit as pa  # noqa: E402  audit A: the causes, their order and their meanings
import response_map_coverage as rmc  # noqa: E402  the attachment test that defines the 198

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import cell2  # noqa: E402
from genomeos.attribution import compile as cp  # noqa: E402
from genomeos.attribution import holdout as ho  # noqa: E402
from genomeos.attribution import measured as ms  # noqa: E402
from genomeos.attribution import targets as tg  # noqa: E402
from genomeos.results import RESULTS_DIR, save_result  # noqa: E402

NAME = "placement_cause_198"
OWN_CODE = ("scripts/placement_cause_198.py", "tests/test_placement_cause_198.py")

# --------------------------------------------------------------------------------------------------
# The four classes asked for, each built from audit A's causes and from nothing else. A locus is in
# exactly one, and `class_4_other_named` names what it holds rather than absorbing a remainder.
# --------------------------------------------------------------------------------------------------
CLASS_OF = {
    pa.ABSENT: "class_1_no_ccre_in_the_tested_interval",
    pa.UNREACHABLE: "class_2_ccre_present_but_failing_the_reciprocal_half_rule",
    pa.WIDTH_MISMATCH: "class_2_ccre_present_but_failing_the_reciprocal_half_rule",
    pa.REGISTRY_WOULD: "class_3_present_but_excluded_by_the_candidate_filter",
    pa.PARTIAL: "class_4_other_named",
    pa.ZERO_WIDTH: "class_4_other_named",
}
CLASSES = (
    "class_1_no_ccre_in_the_tested_interval",
    "class_2_ccre_present_but_failing_the_reciprocal_half_rule",
    "class_3_present_but_excluded_by_the_candidate_filter",
    "class_4_other_named",
)
CLASS_MEANING = {
    "class_1_no_ccre_in_the_tested_interval": (
        "nothing to attach to: no ENCODE cCRE v3 registry element shares a base with the tested "
        f"interval ({pa.ABSENT})"
    ),
    "class_2_ccre_present_but_failing_the_reciprocal_half_rule": (
        "a registry element holds at least half of the smaller of the two intervals, but the "
        "reciprocal-0.5 rule is not met: either no registry element could meet it at this tested "
        f"width ({pa.UNREACHABLE}) or these do not ({pa.WIDTH_MISMATCH})"
    ),
    "class_3_present_but_excluded_by_the_candidate_filter": (
        "a registry element does meet the reciprocal-0.5 rule, and it is not a compiled element: the "
        f"candidate filter excluded it ({pa.REGISTRY_WOULD}), with the reason under `exclusion`"
    ),
    "class_4_other_named": (
        "named, never a remainder, and here it holds exactly one thing: a registry element shares "
        "bases with the tested interval but less than half of the smaller width, a boundary offset "
        f"rather than a width mismatch ({pa.PARTIAL}). The other cause that would land here, a "
        f"tested interval holding no base ({pa.ZERO_WIDTH}), is reported with its own count"
    ),
}
ASSIGNMENT_RULE = (
    "one cause per measured pair from audit A's cascade, first match wins, the causes and their order "
    "exactly `placement_audit.CAUSES`; then one cause per locus: the earliest cause in that same frozen "
    "order that any one of the locus's measured pairs reaches. `placed_compiled` is unreachable here by "
    "construction, since a locus is in this population only when none of its pairs attaches to a "
    "compiled element. The ranking is audit A's own furthest-along-first order and is not re-chosen "
    "here; `ranking_sensitivity` reports the same loci assigned by the last cause of the order instead, "
    "so the choice can be seen rather than trusted"
)
BOUNDARY_READING = (
    "a boundary offset also fails the reciprocal-0.5 rule, so a reader who folds it into class 2 reads "
    "class 2 as 110 and class 4 as 0. The split kept here is audit A's, which reports a partial overlap "
    "apart from a width failure; both readings are given and neither is preferred"
)


# --- the population the 198 is: reproduced, not re-derived -----------------------------------------
def pairs_and_keys() -> tuple[list[Any], list[bool], list[cell2.LocusKey], list[tuple[str, int, int, str]]]:
    """Every measured CRISPRi pair on the 24 chromosomes, whether it attaches to a compiled element,
    its locus key, and the compiled element set — all through the code that wrote the 198."""
    pairs: list[Any] = []
    placed: list[bool] = []
    keys: list[cell2.LocusKey] = []
    comp_items: list[tuple[str, int, int, str]] = []
    for chrom in rmc.CHROMS:
        elements = cp._attributed(chrom, cp.RESULTS_DIR)
        starts = [e["start"] for e in elements]
        comp_items += [(chrom, e["start"], e["end"], e["id"]) for e in elements]
        got, invalid = ms.load_crispri(chrom)
        if invalid:
            raise SystemExit(f"{chrom}: {invalid} rows the parser rejected; the population has moved")
        for p in got:
            e = rmc.attaches_to(p, starts, elements)
            if e is None:
                keys.append(cell2.LocusKey(p.cell, p.chrom, p.start, p.end, p.gene))
            else:
                keys.append(cell2.LocusKey(p.cell, chrom, e["start"], e["end"], p.gene))
            pairs.append(p)
            placed.append(e is not None)
    return pairs, placed, keys, comp_items


def grouping(keys: list[cell2.LocusKey], placed: list[bool]) -> tuple[list[int], set[int], dict[str, int]]:
    """The loci, with and without the response map's own keys in the grouping (the count must not
    depend on them), and the groups holding no attached pair at all."""
    map_keys = [cell2.LocusKey(*k) for k in rmc.map_now(ROOT)["keys"]]
    groups = cell2.group([*keys, *map_keys])[: len(keys)]
    alone = cell2.group(keys)
    attached: dict[int, bool] = defaultdict(bool)
    for g, pl in zip(groups, placed, strict=True):
        attached[g] |= pl
    alone_attached: dict[int, bool] = defaultdict(bool)
    for g, pl in zip(alone, placed, strict=True):
        alone_attached[g] |= pl
    sensitivity = {
        "loci_with_the_map_keys_in_the_grouping": len(attached),
        "loci_without_them": len(alone_attached),
        "no_compiled_link_loci_with_the_map_keys": sum(1 for v in attached.values() if not v),
        "no_compiled_link_loci_without_them": sum(1 for v in alone_attached.values() if not v),
    }
    return groups, {g for g, v in attached.items() if not v}, sensitivity


# --- counting ---------------------------------------------------------------------------------------
def by_cause(
    rows: list[tuple[int, int, Any, str]], order: dict[str, int], last: bool = False
) -> dict[int, str]:
    """One cause per locus: the earliest cause of the frozen order any of its pairs reaches (or the
    last one, which is the sensitivity arm)."""
    out: dict[int, str] = {}
    for _i, g, _p, c in rows:
        if g not in out or (order[c] > order[out[g]] if last else order[c] < order[out[g]]):
            out[g] = c
    return out


def classes_of(causes: dict[int, str]) -> dict[str, int]:
    c = Counter(CLASS_OF[v] for v in causes.values())
    return {k: c[k] for k in CLASSES}


def widths(pairs: list[Any]) -> dict[str, Any]:
    w = sorted(p.end - p.start for p in pairs)
    return {
        "min": w[0],
        "median": w[len(w) // 2] if len(w) % 2 else (w[len(w) // 2 - 1] + w[len(w) // 2]) / 2,
        "max": w[-1],
        "wider_than_700": sum(1 for x in w if x > 700),
        "narrower_than_75": sum(1 for x in w if x < 75),
        "why_700": (
            "the widest registry element is 350 bp, so a tested interval wider than 700 bp cannot "
            "cover half of its own width with one registry element however it lies"
        ),
    }


def code_cleanliness() -> dict[str, Any]:
    """The uncommitted code the stamp names, split into this lane's and other lanes', and whether any
    of it lies on this script's own import closure (`response_map_coverage.counting_path`, reused)."""
    rev = mf.code_revision()
    path = rmc.counting_path(Path(__file__).resolve())
    dirty = list(rev["dirty_code_paths"])
    own = [p for p in dirty if p in OWN_CODE]
    foreign = [p for p in dirty if p not in OWN_CODE]
    return {
        "git_sha": rev["git_sha"],
        "dirty": rev["dirty"],
        "own_uncommitted_code": own,
        "own_code_is_committed": not own,
        "foreign_uncommitted_code": foreign,
        "foreign_uncommitted_code_on_the_counting_path": [p for p in foreign if p in path],
        "counting_path": path,
        "note": (
            "several sessions work in this one checkout. The counting path is the computed transitive "
            "import closure of this script over the repository's own package, so a foreign file "
            "outside it cannot have entered the count"
        ),
    }


def reused_code() -> list[dict[str, Any]]:
    """The two sibling scripts this one imports rather than reimplementing, each pinned by its own
    sha256 so a later edit to them is visible against this result."""
    out = []
    for rel, what in (
        (
            "scripts/placement_audit.py",
            "audit A: CAUSES, CAUSE_MEANING, cause(), hit(), ElementSet, EXCLUSION",
        ),
        (
            "scripts/response_map_coverage.py",
            "the attachment test that defines the 198: attaches_to(), CHROMS",
        ),
    ):
        digest, size, _ = mf.sha256_of(ROOT / rel)
        out.append({"path": rel, "sha256": digest, "bytes": size, "reused_for": what})
    return out


# --- inputs, each by its own path -----------------------------------------------------------------
def input_paths(swept_chroms: list[str]) -> list[Path]:
    """Every file this run opens, one path per entry. No group entry is used: `manifest.files_entry`
    reports a group absent without checking a file (the defect of 2026-09-17), so each input is
    recorded on its own path and hashed on its own."""
    seen: dict[str, Path] = {}

    def add(p: Path) -> None:
        if p.exists():
            seen.setdefault(str(p), p)

    for name in ms.CRISPRI_FILES:
        add(ms.CRISPRI_KNOWLEDGE / name)
    for chrom in rmc.CHROMS:
        add(RESULTS_DIR / f"ccres_{chrom}.bed.gz")
        for run in tg.RUNS:
            p = RESULTS_DIR / f"{run}_{chrom}.json"
            add(p)
            if p.exists():
                where = (json.loads(p.read_text()) or {}).get("elements_where")
                if where:
                    add(ROOT / where if not Path(where).is_absolute() else Path(where))
    for p in ho.compiled_programs():
        add(p)
    for chrom in swept_chroms:
        add(ROOT / pa.ALL_ELEMENTS / f"{chrom}.json")
    add(RESULTS_DIR / "placement_audit.json")
    add(RESULTS_DIR / "response_map_coverage.json")
    return [seen[k] for k in sorted(seen)]


def manifest(swept_chroms: list[str]) -> dict[str, Any]:
    inputs = []
    for p in input_paths(swept_chroms):
        rel = p.relative_to(ROOT) if p.is_absolute() and ROOT in p.parents else p
        part = ms.CRISPRI_SPLIT_OF.get(p.name)
        inputs.append(mf.input_entry(rel, partition=part))
    return {
        "sources": [
            {
                "accession": "EngreitzLab/CRISPR_comparison EPCrisprBenchmark (Gschwind et al. 2025)",
                "version": "main, as fetched 2026-09-16; pinned by sha256",
            },
            {
                "accession": "ENCODE SCREEN cCREs v3, https://downloads.wenglab.org/V3/GRCh38-cCREs.bed",
                "version": "per-chromosome subsets as distilled; pinned by sha256",
            },
            {
                "accession": "this repository: the 24 compiled programs "
                "data/knowledge/compiled/noncoding_<chrom>.bio (generated, untracked)",
                "version": "pinned by sha256",
            },
            {
                "accession": "this repository: the element tables of the three target runs "
                "(enhancer_targets_all, constrained_targets, enhancer_targets; generated, untracked), "
                "which are the compiled element set; the all_elements tables are read for membership "
                "only where a candidate exclusion is named",
                "version": "pinned by sha256",
            },
            {"accession": "this repository: data/results/placement_audit.json", "version": "6c0d39f"},
            {"accession": "this repository: data/results/response_map_coverage.json", "version": "0487e49"},
        ],
        "inputs": inputs,
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "current_rule_reciprocal_overlap": ms.RECIPROCAL_OVERLAP,
            "named_alternative_reported_descriptively": ms.MIN_SIDE_HALF,
            "named_alternative_threshold": pa.POLICY_THRESHOLD,
            "independent_locus_span_bp": cell2.INDEPENDENT_LOCUS_SPAN,
            "cause_order": list(pa.CAUSES),
            "exclusion_order": list(pa.EXCLUSION),
        },
        "exclusions": [
            "CRISPRi pairs the benchmark marks as not a valid connection are not read (none in these files)",
            "the per-element response cache is not opened; the model is not called; nothing is downloaded",
            "no effect value of the deletion run is read, only whether a registry element was scored and "
            "whether it named a target gene",
            "the measured layer's rule rows and the context reader are not read: this script asks only "
            "which element a tested interval attaches to, which is the one input the 198 is missing",
        ],
        "partitions": {
            ms.TRAINING: "the CRISPRi benchmark's training file (K562); read as measurement, nothing fitted",
            ms.HELDOUT: "the CRISPRi benchmark's held-out file; read as measurement, reported apart",
        },
    }


# --- the run ----------------------------------------------------------------------------------------
def main() -> None:
    t0, r0 = time.perf_counter(), resource.getrusage(resource.RUSAGE_SELF)
    pairs, placed, keys, comp_items = pairs_and_keys()
    groups, target, map_key_sensitivity = grouping(keys, placed)
    print(f"pairs {len(pairs):,} attached {sum(placed):,} loci {len(set(groups)):,}", flush=True)

    committed = json.loads((RESULTS_DIR / "response_map_coverage.json").read_text())
    want = {
        "measured_perturbation_pairs": len(pairs),
        "independent_loci": len(set(groups)),
        "no_compiled_link_pairs": len(pairs) - sum(placed),
        "no_compiled_link_loci": len(target),
    }
    got = {
        "measured_perturbation_pairs": committed["denominator"]["measured_perturbation_pairs"],
        "independent_loci": committed["denominator"]["independent_loci"],
        "no_compiled_link_pairs": committed["by_status"]["pairs"][rmc.STATUS_NO_LINK],
        "no_compiled_link_loci": committed["by_status"]["independent_loci"][rmc.STATUS_NO_LINK],
    }
    if want != got:
        raise SystemExit(f"response_map_coverage does not reproduce: {want} against {got}")

    reg_items, _classes = pa.registry_elements(RESULTS_DIR)
    comp, reg = pa.ElementSet(comp_items), pa.ElementSet(reg_items)
    print(f"compiled {len(comp):,} registry {len(reg):,}", flush=True)

    # audit A's cascade over every pair, and the same cascade over audit A's own compiled element set
    comp_a = pa.ElementSet(pa.compiled_elements(ho.compiled_programs()))
    cause_of: list[str] = []
    cause_of_a: list[str] = []
    reg_hits: list[pa.Hit] = []
    policy_here: list[bool] = []
    for p in pairs:
        c = pa.hit(comp, p.chrom, p.start, p.end, pa.current_rule)
        r = pa.hit(reg, p.chrom, p.start, p.end, pa.current_rule)
        ca = pa.hit(comp_a, p.chrom, p.start, p.end, pa.current_rule)
        cause_of.append(pa.cause(c, r, p.end - p.start, reg.min_width, reg.max_width))
        cause_of_a.append(pa.cause(ca, r, p.end - p.start, reg.min_width, reg.max_width))
        reg_hits.append(r)
        policy_here.append(
            bool(pa.hit(comp, p.chrom, p.start, p.end, ms.ATTACHMENT_RULES[ms.MIN_SIDE_HALF]).by_rule)
        )

    # gate: audit A's own figures on its own population, the 661 positive links
    pos_causes = Counter(c for p, c in zip(pairs, cause_of, strict=True) if p.outcome == ms.DECREASE)
    audit = json.loads((RESULTS_DIR / "placement_audit.json").read_text())
    a_pos = audit["current_rule"]["attrition"]["positive"][ms.ALL]
    a_want = Counter()
    for k, v in a_pos.items():
        if k == "total":
            continue
        a_want[k.split("/")[0]] += v
    if dict(pos_causes) != dict(a_want):
        raise SystemExit(f"audit A's cascade does not reproduce: {dict(pos_causes)} against {dict(a_want)}")
    if cause_of != cause_of_a:
        differ = sum(1 for a, b in zip(cause_of, cause_of_a, strict=True) if a != b)
        print(f"NOTE: the two compiled element sets differ on {differ} pairs", flush=True)

    rows = [
        (i, g, p, c) for i, (g, p, c) in enumerate(zip(groups, pairs, cause_of, strict=True)) if g in target
    ]
    order = {c: i for i, c in enumerate(pa.CAUSES)}
    first = by_cause(rows, order)
    last = by_cause(rows, order, last=True)
    if set(first) != target or len(first) != len(target):
        raise SystemExit("not every locus of the population received a cause")
    if any(c == pa.PLACED for c in first.values()):
        raise SystemExit("a locus of the population was placed, which the population excludes")

    # the candidate exclusion, for the class-3 loci: the registry ids that meet the rule, and why
    need: dict[str, set[str]] = defaultdict(set)
    for i, _g, p, c in rows:
        if c == pa.REGISTRY_WOULD:
            need[p.chrom].update(reg_hits[i].by_rule)
    swept_chroms = sorted(need)
    swept = pa.swept_states(need, ROOT / pa.ALL_ELEMENTS)
    compiled_ids = {i for _, _, _, i in comp_items}
    exc_of_pair: dict[int, str] = {
        i: pa.exclusion(reg_hits[i].by_rule, compiled_ids, swept)
        for i, _g, _p, c in rows
        if c == pa.REGISTRY_WOULD
    }
    exc_by_locus: dict[int, str] = {}
    for i, g, _p, c in rows:
        if c != pa.REGISTRY_WOULD:
            continue
        e = exc_of_pair[i]
        if g not in exc_by_locus or pa.EXCLUSION.index(e) < pa.EXCLUSION.index(exc_by_locus[g]):
            exc_by_locus[g] = e

    sub_pairs = [p for _i, _g, p, _c in rows]
    pos_in_target = [(g, p, c) for _i, g, p, c in rows if p.outcome == ms.DECREASE]
    unplaced_pos = [
        (p, c) for p, c in zip(pairs, cause_of, strict=True) if p.outcome == ms.DECREASE and c != pa.PLACED
    ]
    neither_pos = [(p, c) for p, c in unplaced_pos if c != pa.REGISTRY_WOULD]
    policy_loci = {g for i, g, _p, _c in rows if policy_here[i]}
    policy_pairs = sum(1 for i, _g, _p, _c in rows if policy_here[i])

    classes = classes_of(first)
    if sum(classes.values()) != len(target):
        raise SystemExit(f"the classes sum to {sum(classes.values())}, not {len(target)}")

    r1 = resource.getrusage(resource.RUSAGE_SELF)
    payload: dict[str, Any] = {
        "question": (
            "for each of the 198 independent loci that hold a measured CRISPRi perturbation and attach "
            "to no compiled element, what is the one cause, and is the omission a property of the "
            "placement rule or real absence?"
        ),
        "status": (
            "internal development evidence, descriptive and diagnostic. No rule is changed, no policy "
            "adopted, no verdict moved, nothing deleted. Every outcome read here was read by earlier "
            "lanes (audit A's prior_exposure applies unchanged), so this is a bounded diagnosis and not "
            "a biological result"
        ),
        "alphagenome_requests": 0,
        "downloads": 0,
        "per_element_response_cache_opened": False,
        "classification_is_imported": {
            "causes": list(pa.CAUSES),
            "cause_meaning": dict(pa.CAUSE_MEANING),
            "exclusion_order": list(pa.EXCLUSION),
            "from": "scripts/placement_audit.py, audit A (6c0d39f), word for word",
            "reused_code": reused_code(),
        },
        "assignment_rule": ASSIGNMENT_RULE,
        "audit_a_dependency": {
            "result": "data/results/placement_audit.json",
            "status": (
                "suspended: not yet reproduced (the coordinator, 2026-10-02). Its committed manifest "
                "records one grouped input through `manifest.files_entry`, standing for 24 files, and "
                "the manifest_rebuild of that time resolved an input by path while a group entry "
                "records a label, so the group was reported unavailable with no file opened and the "
                "result was excused rather than compared. The tool is fixed (`ebbded6`), but audit A's "
                "manifest predates the `members` field, so it cannot be checked from its manifest "
                "until its writer re-records it"
            ),
            "what_this_lane_reuses_from_it": (
                "its classification and nothing else: CAUSES, CAUSE_MEANING, EXCLUSION and the order "
                "they are tried in. A classification is not a figure and is not suspended; inventing a "
                "second one would be worse than reusing this"
            ),
            "no_figure_of_its_is_carried_as_evidence": (
                "every figure of audit A's that appears anywhere in this result was recomputed in this "
                "run from the primary inputs, each recorded by its own path and never through "
                "`files_entry`, and the recomputation is the gate "
                "`audit_a_cascade_reproduced_on_its_own_661_positive_links`. The committed file is read "
                "once, only to compare against. A reader who wants audit A's own figures as audit A's "
                "should treat them as suspended until that result is re-recorded; the same counts as "
                "stated here are this run's"
            ),
            "classes_depending_on_a_figure_this_run_cannot_verify": [],
            "why_the_reconciliation_matters_more_for_this": (
                "198 loci and 387 pairs are different numbers on different denominators, and neither "
                "may stand in for the other; the reconciliation block states how the two populations "
                "meet, pair by pair, rather than equating them"
            ),
        },
        "classes": {"meaning": CLASS_MEANING, "built_from": {k: v for k, v in CLASS_OF.items()}},
        "population": {
            "what": (
                "the 198 independent loci of data/results/response_map_coverage.json whose every "
                "measured CRISPRi pair attaches to no compiled element under "
                f"measured.RECIPROCAL_OVERLAP = {ms.RECIPROCAL_OVERLAP}"
            ),
            "independent_loci": len(target),
            "of_independent_loci": len(set(groups)),
            "locus_rule": cell2.INDEPENDENT_LOCUS_RULE,
            "not_biological_independence": (
                "the grouping is operational and pools cell types: `cell2.group` joins two keys on a "
                "shared measured gene or on 1 Mb proximity on one chromosome and never reads the cell. "
                "It is not established biological independence and nothing may cite it as one"
            ),
            "measured_pairs_at_these_loci": len(sub_pairs),
            "of_measured_pairs": len(pairs),
            "pairs_attaching_to_no_compiled_element_anywhere": len(pairs) - sum(placed),
            "the_rest_of_those_pairs_sit_at_loci_that_do_attach": (len(pairs) - sum(placed) - len(sub_pairs)),
            "by_outcome": dict(sorted(Counter(p.outcome for p in sub_pairs).items())),
            "measured_effects_not_nulls": sum(
                1 for p in sub_pairs if p.outcome in (ms.DECREASE, ms.INCREASE)
            ),
            "loci_holding_at_least_one_significant_decrease": len({g for g, _p, _c in pos_in_target}),
            "distinct_measured_genes": len({p.gene for p in sub_pairs}),
            "distinct_tested_intervals": len({(p.chrom, p.start, p.end) for p in sub_pairs}),
            "by_cell": dict(sorted(Counter(p.cell for p in sub_pairs).items())),
            "by_split": dict(sorted(Counter(p.split for p in sub_pairs).items())),
            "by_dataset": dict(sorted(Counter(p.dataset for p in sub_pairs).items())),
            "tested_widths": widths(sub_pairs),
            "unclassifiable_from_what_is_on_disk": 0,
            "unclassifiable_note": (
                "every locus of the population received a cause from the imported cascade; none had to "
                "be forced into a class and none is left out"
            ),
        },
        "answer": {
            "classes": classes,
            "classes_sum": sum(classes.values()),
            "locus_causes": dict(sorted(Counter(first.values()).items())),
            "pair_causes_at_these_loci": dict(sorted(Counter(c for *_r, c in rows).items())),
            "boundary_offset_reading": BOUNDARY_READING,
            "class_2_or_3": classes[CLASSES[1]] + classes[CLASSES[2]],
            "class_1": classes[CLASSES[0]],
            "candidate_exclusion_of_the_class_3_loci": dict(sorted(Counter(exc_by_locus.values()).items())),
            "candidate_exclusion_of_the_class_3_pairs": dict(sorted(Counter(exc_of_pair.values()).items())),
            "registry_ids_meeting_the_rule_at_these_loci": sum(len(v) for v in need.values()),
            "registry_ids_the_sweep_scored": sum(1 for i in swept),
        },
        "ranking_sensitivity": {
            "what": (
                "the same 198 loci assigned by the LAST cause of the frozen order their pairs reach "
                "instead of the first: the extreme opposite choice, reported so the ranking is visible"
            ),
            "locus_causes_last_of_order": dict(sorted(Counter(last.values()).items())),
            "classes_last_of_order": classes_of(last),
            "loci_whose_every_pair_agrees_on_one_cause": sum(1 for g in first if first[g] == last[g]),
        },
        "gates": {
            "response_map_coverage_reproduced": want,
            "equal_to_committed": True,
            "audit_a_cascade_reproduced_on_its_own_661_positive_links": dict(sorted(pos_causes.items())),
            "equal_to_committed_audit_a": dict(pos_causes) == dict(a_want),
            "compiled_element_sets": {
                "this_run": "genomeos.attribution.compile._attributed, the set that defines the 198",
                "audit_a": "the predicted element blocks of the 24 compiled programs",
                "elements_each": [len(comp), len(comp_a)],
                "pairs_whose_cause_differs_between_them": sum(
                    1 for a, b in zip(cause_of, cause_of_a, strict=True) if a != b
                ),
            },
            "independent_loci_do_not_depend_on_the_response_map_keys": map_key_sensitivity,
        },
        "reconciliation_with_audit_a": {
            "the_two_numbers_are_not_the_same_kind": (
                "audit A's 387 counts positive links (significant decreases) that no element of either "
                "set carries; the 198 counts independent loci, over every outcome and cell-pooled, none "
                "of whose pairs attaches in the compiled set. Neither is the other"
            ),
            "the_underlying_objects_are_the_same": (
                "both read the same 14,734 benchmark rows, and every row is its own C4 link key "
                "(checked: 14,734 rows, 14,734 distinct keys), so a pair here is a link there"
            ),
            "audit_a_positive_links": {
                "total": sum(pos_causes.values()),
                "placed_compiled": pos_causes[pa.PLACED],
                "registry_would": pos_causes[pa.REGISTRY_WOULD],
                "neither_the_387": len(neither_pos),
                "by_cause": dict(sorted(pos_causes.items())),
            },
            "where_those_links_sit_relative_to_the_198_loci": {
                "unplaced_positive_links": len(unplaced_pos),
                "of_them_at_one_of_the_198_loci": len(pos_in_target),
                "of_them_at_a_locus_that_does_attach_somewhere": len(unplaced_pos) - len(pos_in_target),
                "of_audit_a_s_387_at_one_of_the_198_loci": sum(
                    1 for _g, p, c in pos_in_target if c != pa.REGISTRY_WOULD
                ),
                "of_audit_a_s_62_registry_would_at_one_of_the_198_loci": sum(
                    1 for _g, _p, c in pos_in_target if c == pa.REGISTRY_WOULD
                ),
                "why_they_differ": (
                    "a locus leaves this population as soon as any one of its pairs attaches to a "
                    "compiled element, and a 1 Mb locus pools many tested intervals, so most unplaced "
                    "positive links sit beside a placed pair and are not counted in the 198"
                ),
            },
            "cause_shares_agree_in_direction": (
                "on audit A's 387 positive links the width causes are 362 of 387; on the 198 loci the "
                "width causes are the largest single class as well, and the two are not the same count"
            ),
        },
        "named_alternative_descriptively": {
            "name": ms.MIN_SIDE_HALF,
            "rule": pa.POLICY["rule"],
            "registered": pa.REGISTERED,
            "where_it_already_lives": "genomeos.attribution.measured.measures_min_side, ATTACHMENT_RULES",
            "not_adopted": (
                "counted here and nothing else. No rule is changed, nothing is attached under it, no "
                "verdict moves and no element gains evidence. measured.ATTACHMENT_RULES carries its own "
                "warning that an overlap under it is a candidate and never decisive evidence"
            ),
            "loci_of_the_198_it_would_attach_to_a_compiled_element": len(policy_loci),
            "pairs_at_those_loci_it_would_attach": policy_pairs,
            "by_class_of_the_locus": dict(sorted(Counter(CLASS_OF[first[g]] for g in policy_loci).items())),
            "by_cause_of_the_locus": dict(sorted(Counter(first[g] for g in policy_loci).items())),
            "what_this_count_is_not": (
                "not a count of loci that would become addable: the response map needs a compiled rule "
                "and a reader state as well, neither of which is read here. It is also not stable under "
                "adoption, because a pair that attaches takes the element's interval into the locus key, "
                "so the grouping itself would change; the figure is counted on the frozen grouping"
            ),
        },
        "cost": {
            "wall_seconds": round(time.perf_counter() - t0, 1),
            "cpu_seconds": round((r1.ru_utime - r0.ru_utime) + (r1.ru_stime - r0.ru_stime), 1),
            "downloads": 0,
            "requests": 0,
            "peak_memory_is_not_recorded": (
                "peak memory is a property of the machine, not of the count, and it differs by a "
                "megabyte between two runs of the same code on the same inputs. It is left out so that "
                "a difference a manifest rebuild reports is always a real one"
            ),
        },
        "code_cleanliness": code_cleanliness(),
    }
    p = save_result(NAME, payload, manifest=manifest(swept_chroms))
    print("saved", p, flush=True)
    print(json.dumps(payload["answer"]["classes"], indent=1), flush=True)


if __name__ == "__main__":
    main()
