# SPDX-License-Identifier: AGPL-3.0-or-later
"""Item 12 S4, applied once: the unchanged compiled labels judged claim by claim against every held-out
source, with target accuracy, role accuracy and coverage reported apart, per axis, with denominators.

    uv run python scripts/s4_correctness_run.py              # v1, S4's result as committed
    uv run python scripts/s4_correctness_run.py --rule v2    # the same claims and sources under v2
    uv run python scripts/s4_correctness_run.py --rule v3    # the same claims and sources under v3

Registered 2026-09-29 before it was run (genomeos/attribution/correctness.py, docs/ATTRIBUTION.md "S4
registered"). An internal development benchmark only. No model request; the per-element response cache
is never opened. Writes data/results/attribution_correctness.json.

The direction rule is versioned since 2026-09-29 (correctness.DIRECTION_RULE, docs/ATTRIBUTION.md "The
direction rule, versioned"). Without `--rule` this script pins v1, the rule S4's committed result was
judged by, and writes exactly what it wrote. `--rule v2` judges the same claims with the same sources
under v2 and writes a new result, data/results/attribution_correctness_v2.json, beside the v1 file,
which it reads and never writes: it reruns v1 in the same process, refuses to write if v1 does not
reproduce the committed file, and reports every count under both rules side by side, with each claim
whose verdict moved.

The target rule is versioned too (correctness.TARGET_RULE, docs/ATTRIBUTION.md "The target rule,
versioned"). `--rule v3` judges the same claims with the same sources under v1, v2 and v3, refuses to
write unless v1 reproduces attribution_correctness.json key by key and v2 reproduces the judged body
and side-by-side counts of attribution_correctness_v2.json, and writes a new result,
data/results/attribution_correctness_v3.json, beside both, which it reads and never writes. It reports
every count under the three rules side by side, the target and direction counts in the owner's terms
(correctness.IN_CONTEXT_TERMS) with the caution beside them, and each claim whose verdict moved
between v2 and v3.

`--compiled DIR` and `--chroms` exist so the pipeline can be exercised on a synthetic program; the
registered run uses neither.
"""

from __future__ import annotations

import argparse
import json
import time
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

from genomeos import manifest as mf
from genomeos.attribution import correctness as co
from genomeos.attribution import holdout as ho
from genomeos.attribution import measured as ms
from genomeos.results import RESULTS_DIR, save_result

NAME = "attribution_correctness"
NAMES = {co.RULE_V1: NAME, co.RULE_V2: f"{NAME}_v2"}
NAME_V3 = f"{NAME}_v3"  # the v3 re-judge, beside the v1 and v2 results
REPRESSION_TRACE = RESULTS_DIR / "repression_trace.json"
#: what records the run rather than the result, so a rerun of v1 may differ there
RUN_KEYS = ("result", "date", "seconds", "result_manifest")
ID1 = ("EH38E3426791", "ID1")  # the element and gene of the claim Albert named (repression trace)


@mf.depends_on_models("alphagenome")  # the unchanged labels are the compiled predicted layer
def manifest(programs: list[Path], parameters: dict[str, Any], also: tuple[Path, ...] = ()) -> dict[str, Any]:
    inputs = [mf.input_entry(p, partition=None) for p in programs]
    for name in ms.CRISPRI_FILES:
        p = ms.CRISPRI_KNOWLEDGE / name
        if p.exists():
            inputs.append(mf.input_entry(p, partition=ms.CRISPRI_SPLIT_OF[name]))
    for acc in ho.mpra.FILES.values():
        p = ho.mpra.KNOWLEDGE / f"{acc}.bed.gz"
        if p.exists():
            inputs.append(mf.input_entry(p, partition=None))
    for p in (ho.vista.locus_path(ho.vista.KNOWLEDGE), ms.satmut_knowledge.DATA_PATH):
        if p.exists():
            inputs.append(mf.input_entry(p, partition=None))
    inputs += [mf.input_entry(p, partition=None) for p in sorted(ho.GTEX_KNOWLEDGE.glob("hits_*.tsv"))]
    inputs += [mf.input_entry(p, partition=None) for p in ho.gencode_paths()]
    inputs += [mf.input_entry(p, partition=None) for p in also]  # v2: the committed v1 result, the trace
    return {
        "sources": [
            {
                "accession": "EngreitzLab/CRISPR_comparison EPCrisprBenchmark (Gschwind et al. 2025)",
                "version": "main, as fetched; pinned by sha256",
            },
            {
                "accession": f"ENCODE {ho.mpra.LIBRARY} lentiMPRA ({', '.join(ho.mpra.FILES.values())})",
                "version": "as fetched 2026-09-12; pinned by sha256",
            },
            {"accession": "VISTA Enhancer Browser locus table", "version": "pinned by sha256"},
            {"accession": "GEO GSE126550 (Kircher et al. 2019)", "version": "pinned by sha256"},
            {
                "accession": "GTEx v8 single-tissue cis-eQTL, distilled by genomeos/attribution/eqtl.py",
                "version": "pinned by sha256",
            },
            {"accession": "GENCODE v50 (per chromosome)", "version": "pinned by sha256"},
            {
                "accession": "this repository: the 24 compiled programs "
                "data/knowledge/compiled/noncoding_<chrom>.bio (generated, untracked)",
                "version": "pinned by sha256",
            },
        ],
        "inputs": inputs,
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": parameters,
        "exclusions": [
            "the measured twins (`*_measured` blocks and experimental rules) are the evidence and are not "
            "read as labels",
            "`unknown`, `unassigned` and unchosen `|` alternatives are counted as not claims, never judged",
            "CRISPRi pairs the benchmark marks as not a valid connection are not read",
            "observation kinds holdout does not load (contact, conservation, allelic readouts, the model, "
            "annotation and registry) are in the table and judge nothing in this run",
            "the AlphaGenome per-element response cache is not opened; the model is not called",
        ],
        "partitions": {
            "held_out_source": "every holdout source judges the claims; the unchanged labels read none of "
            "them, and holdout.check_provenance is called for each",
            ms.TRAINING: "the CRISPRi benchmark's training file: judges like any source",
            ms.HELDOUT: "the CRISPRi benchmark's held-out file: judges (evaluation), never evidence for a "
            "labelling; verdicts decided with one of its pairs are counted per axis",
        },
    }


def judged_under(rule: str, a: argparse.Namespace, units: dict, refs: dict) -> co.Report:
    nc: dict[str, Counter] = defaultdict(Counter)
    claims = co.compiled_claims(a.compiled, a.chroms, nc)
    return co.judge(claims, co.unchanged_labels(), units=units, references=refs, not_claims=nc, rule=rule)


def show(body: dict[str, Any]) -> None:
    for q, block in body["quantities"].items():
        for axis, s in block.items():
            print(q, axis, f"{s['numerator']:,} / {s['denominator']:,}", s["beside"], flush=True)


def v1_params(a: argparse.Namespace) -> dict[str, Any]:
    return {
        "registration": co.registration(),
        "holdout": ho.registration(),
        "compiled": str(a.compiled),
        "chroms": a.chroms or "all",
    }


def v1_payload(report: co.Report, units: dict, t0: float) -> dict[str, Any]:
    """S4's payload exactly as committed: the v1 run writes it, and the v2 run compares it."""
    body = report.to_dict()
    return {
        "question": "of the claims the unchanged compiled labels state, on each of five axes, how many can "
        "any held-out observation judge, and of those how many does it establish or refute",
        "status": co.STATUS,
        "reuse": co.REUSE,
        "may_be_called": co.MAY_BE_CALLED,
        "may_not_be_called": list(co.MAY_NOT_BE_CALLED),
        "expected_before_the_run": co.EXPECTED,
        "registration": co.registration(),
        "units_per_source": {s: len(units[s]) for s in ho.sources(units)},
        **body,
        "alphagenome_requests": 0,
        "per_element_response_cache_opened": False,
        "seconds": round(time.time() - t0, 1),
    }


def programs_of(a: argparse.Namespace) -> list[Path]:
    programs = ho.compiled_programs(a.compiled)
    if a.chroms:
        programs = [p for p in programs if p.stem.removeprefix("noncoding_") in set(a.chroms)]
    return programs


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--compiled", type=Path, default=ho.COMPILED_DIR)
    ap.add_argument("--chroms", nargs="*", default=None)
    ap.add_argument("--dry", action="store_true", help="print the quantities, write nothing")
    ap.add_argument(
        "--rule",
        choices=co.RULES,
        default=co.RULE_V1,
        help="the direction rule (default v1: S4 as committed)",
    )
    ap.add_argument("--results", type=Path, default=RESULTS_DIR, help=argparse.SUPPRESS)
    a = ap.parse_args(argv)
    if a.rule == co.RULE_V3:
        return main_v3(a)
    if a.rule == co.RULE_V2:
        return main_v2(a)
    t0 = time.time()
    units = ho.all_units()
    refs = ho.crispri_references()
    report = judged_under(co.RULE_V1, a, units, refs)  # S4's rule, pinned (correctness.RULE_RECORD)
    show(report.to_dict())
    if a.dry:
        return
    payload = v1_payload(report, units, t0)
    p = save_result(NAME, payload, a.results, manifest=manifest(programs_of(a), v1_params(a)))
    print("saved", p, round(time.time() - t0, 1), "s", flush=True)


# --- v2 ----------------------------------------------------------------------------------------------
def _text(d: dict[str, Any], drop: set[str]) -> str:
    """The result as its file holds it, keys in order, without what records the run."""
    return json.dumps({k: v for k, v in d.items() if k not in drop}, indent=2)


def _flat(tally: dict[str, Any]) -> dict[str, int]:
    """A tally's counts as `block.key` paths (the counters inside it, and its integer fields)."""
    out: dict[str, int] = {}
    for k, v in tally.items():
        if isinstance(v, int):
            out[k] = v
        elif isinstance(v, dict) and k != "cross_cell_findings":
            out.update({f"{k}.{kk}": vv for kk, vv in v.items() if isinstance(vv, int)})
    return out


def _decided_in(v: co.Verdict) -> str:
    cells = {co.norm_cell(o.cell) == co.norm_cell(v.claim.cell) for o in v.deciding}
    if not v.deciding:
        return "not decided"
    if cells == {True}:
        return "the stated cell"
    return "another cell only" if cells == {False} else "the stated cell and another"


def _claim_dict(c: co.Claim) -> dict[str, Any]:
    d = co.Verdict(c, co.NOT_JUDGED).to_dict()
    return {k: d[k] for k in ("element", "locus", "block", "axis", "value", "gene", "cell")}


def _pair(v: co.Verdict) -> tuple[Any, ...]:
    return (v.verdict, v.reason, v.detail, v.explanations, v.deciding)


def moved_claims(rep1: co.Report, rep2: co.Report, units: dict) -> list[tuple[co.Verdict, co.Verdict]]:
    """Every claim either rule judged, or v2 did not assess, judged again under both rules with the
    observations `judge` gave it; a claim moves when its verdict, reason, detail or deciding observations
    differ. A claim that moved between two other not-judged reasons is in no list: `unexplained` finds
    it in the counts."""
    ev = co.Evidence(units, rep2.sources)
    claims = {v.claim for v in (*rep1.judged, *rep2.judged, *rep2.not_assessed)}
    out = []
    for c in sorted(
        claims, key=lambda c: (c.chrom, c.start, c.end, c.element, c.axis, c.value, c.gene, c.cell)
    ):
        obs = ev.at(c.chrom, c.start, c.end)
        a, b = co.verdict_of(c, obs, co.RULE_V1), co.verdict_of(c, obs, co.RULE_V2)
        if _pair(a) != _pair(b):
            out.append((a, b))
    return out


def unexplained(rep1: co.Report, rep2: co.Report, moved: list[tuple[co.Verdict, co.Verdict]]) -> list[str]:
    """Count changes between the two rules that the moved claims do not account for. The v1 tally minus
    the moved claims' v1 verdicts plus their v2 verdicts must equal the v2 tally, count by count, apart
    from activity `refutable`, which v1 never computed (REFUTABLE_RULE) and is checked against the
    direction verdicts decided in the stated cell."""
    out = []
    for axis in co.AXES:
        m1, m2 = co.AxisTally(rule=co.RULE_V1), co.AxisTally(rule=co.RULE_V2)
        for a, b in moved:
            if a.claim.axis == axis:
                m1.add(a, set())
                m2.add(b, set())
        t1, t2 = _flat(rep1.axes[axis].to_dict()), _flat(rep2.axes[axis].to_dict())
        f1, f2 = _flat(m1.to_dict()), _flat(m2.to_dict())
        for k in sorted(set(t1) | set(t2)):
            if axis == co.ACTIVITY and k == "refutable":
                continue
            want = t1.get(k, 0) - f1.get(k, 0) + f2.get(k, 0)
            if t2.get(k, 0) != want:
                out.append(
                    f"{axis}.{k}: v1 {t1.get(k, 0)}, v2 {t2.get(k, 0)}, moved claims account for {want}"
                )
    in_cell = sum(
        1
        for v in rep2.judged
        if v.claim.axis == co.ACTIVITY
        and v.claim.value in co.DIRECTION
        and _decided_in(v) == "the stated cell"
    )
    if rep2.axes[co.ACTIVITY].refutable != in_cell:
        out.append(
            f"activity.refutable: v2 {rep2.axes[co.ACTIVITY].refutable}, decided in the stated cell {in_cell}"
        )
    return out


def side_by_side(rep1: co.Report, rep2: co.Report) -> dict[str, Any]:
    out = {}
    for axis in co.AXES:
        t1, t2 = rep1.axes[axis], rep2.axes[axis]
        row = {}
        for rule, t in ((co.RULE_V1, t1), (co.RULE_V2, t2)):
            row[rule] = {
                "claims": t.claims,
                "judged": t.judged(),
                **{v: t.verdicts[v] for v in co.VERDICTS},
                "not_judged_by_reason": {r: t.not_judged_by_reason[r] for r in co.REASONS_V2},
                "refutable": t.refutable if (axis != co.ACTIVITY or rule != co.RULE_V1) else None,
            }
        out[axis] = row
    return out


def separate_finding_id1(rep2: co.Report) -> dict[str, Any]:
    """ID1 (Albert, 2026-09-29): the v2 verdict of its direction claim, and the repression trace's
    reading of the model's own K562 value, carried word for word from the trace's result."""
    trace = json.loads(REPRESSION_TRACE.read_text())
    case = next(c for c in trace["four"] if (c["element"], c["gene"]) == ID1)
    v2 = [v for v in (*rep2.judged, *rep2.not_assessed) if (v.claim.element, v.claim.gene) == ID1]
    return {
        "claim": {"element": case["element"], "locus": case["locus"], "gene": case["gene"], **case["claim"]},
        "v1_verdict": case["committed_s4_verdict"],
        "v2_verdicts": [v.to_dict() for v in v2 if v.claim.axis == co.ACTIVITY],
        "repression_trace_class_c": case["why"]["c"],
        "reading": "a separate finding, kept apart from the verdict: the model's cached K562 prediction for "
        "ID1 disagrees with the K562 experiment. That does not directly refute the compiled whole-blood "
        "claim, which v2 does not assess in this context",
    }


def main_v2(a: argparse.Namespace) -> None:
    t0 = time.time()
    units = ho.all_units()
    refs = ho.crispri_references()
    rep2 = judged_under(co.RULE_V2, a, units, refs)
    body = rep2.to_dict()
    show(body)
    t1 = time.time()
    rep1 = judged_under(co.RULE_V1, a, units, refs)
    registered = a.compiled == ho.COMPILED_DIR and not a.chroms
    committed_path = a.results / f"{NAME}.json"
    reproduces = None
    if registered:
        committed = json.loads(committed_path.read_text())
        rerun = v1_payload(rep1, units, t1)
        drop = set(RUN_KEYS)
        reproduces = _text(committed, drop) == _text(json.loads(json.dumps(rerun, default=str)), drop)
        if not reproduces:
            raise SystemExit(f"v1 does not reproduce {committed_path}: nothing written")
    moved = moved_claims(rep1, rep2, units)
    loose = unexplained(rep1, rep2, moved)
    print("moved", len(moved), "unexplained", loose, flush=True)
    if a.dry:
        return
    params = {**v1_params(a), "judge_rule": co.RULE_V2, "rule_registration": co.rule_registration()}
    payload = {
        "question": "the same claims and sources as S4, judged under the direction rule v2: which counts "
        "move from S4's v1 result, and are the moves exactly the verdicts v1 decided only in another cell",
        "rule": co.RULE_V2,
        "status": co.STATUS,
        "reuse": co.REUSE,
        "may_be_called": co.MAY_BE_CALLED,
        "may_not_be_called": list(co.MAY_NOT_BE_CALLED),
        "rule_registration": co.rule_registration(),
        "expected_before_the_run": co.EXPECTED_V2,
        "seen_before_registration": co.SEEN_BEFORE_V2,
        "historical_result": {
            "path": str(committed_path),
            "rule": co.RULE_V1,
            "kept": "unchanged, beside this result; S4's original reading stands as the v1 result",
            "v1_rerun_reproduces_it": reproduces,
            "compared_without": list(RUN_KEYS),
        },
        "side_by_side": side_by_side(rep1, rep2),
        "moved": {
            "claims": len(moved),
            "by_axis": dict(Counter(a_.claim.axis for a_, _ in moved)),
            "v1_verdicts": dict(Counter(a_.verdict for a_, _ in moved)),
            "v1_decided_in": dict(Counter(_decided_in(a_) for a_, _ in moved)),
            "v2_verdicts": dict(
                Counter(f"{b.verdict}:{b.reason}" if b.reason else b.verdict for _, b in moved)
            ),
            "v2_cross_cell": dict(
                Counter("/".join(sorted({f for _, f in b.cross_cell or ()})) or "none" for _, b in moved)
            ),
            "unexplained_count_changes": loose,
            "each": [
                {
                    "claim": _claim_dict(a_.claim),
                    co.RULE_V1: {**a_.to_dict(), "decided_in": _decided_in(a_)},
                    co.RULE_V2: b.to_dict(),
                }
                for a_, b in moved
            ],
        },
        "separate_findings": {"ID1": separate_finding_id1(rep2)} if registered else {},
        "registration": co.registration(),
        "units_per_source": {s: len(units[s]) for s in ho.sources(units)},
        **body,
        "alphagenome_requests": 0,
        "per_element_response_cache_opened": False,
        "seconds": round(time.time() - t0, 1),
    }
    also = (committed_path, REPRESSION_TRACE) if registered else ()
    p = save_result(NAMES[co.RULE_V2], payload, a.results, manifest=manifest(programs_of(a), params, also))
    print("saved", p, round(time.time() - t0, 1), "s", flush=True)


# --- v3 ----------------------------------------------------------------------------------------------
#: the keys of a report that record what was judged; the v3 run compares them with the committed v2 file
JUDGED_KEYS = ("rule", "quantities", "also_reported", "axes", "not_claims", "judged", "not_assessed")
V3_MAY_BE_CALLED = (
    "an internal development benchmark reading of the same S4 claims under v3, in the owner's terms "
    "(correctness.IN_CONTEXT_TERMS): per axis, the claims supported and refuted in the stated context, the "
    "claims unassessed in context, and the claims supported somewhere, each a count of claims with "
    "correctness.IN_CONTEXT_CAUTION beside it, beside the v1 and v2 readings"
)
V3_MAY_NOT_BE_CALLED = (
    *co.MAY_NOT_BE_CALLED,
    "general accuracy: every N of N here, 39 of 39 on target and on direction among them, describes the "
    "small subset assessable in the stated context",
    "evidence that the targets unassessed in context are right or wrong in the cell they name",
    "a correction of S4's committed result or of the v2 re-judge, which stand as the v1 and v2 results",
)


def moved_between(ra: co.Report, rb: co.Report, units: dict) -> list[tuple[co.Verdict, co.Verdict]]:
    """`moved_claims` between the rules of any two reports: every claim either judged or did not assess,
    judged again under both rules with the observations `judge` gave it; a claim moves when its verdict,
    reason, detail or deciding observations differ."""
    ev = co.Evidence(units, rb.sources)
    claims = {v.claim for v in (*ra.judged, *ra.not_assessed, *rb.judged, *rb.not_assessed)}
    out = []
    for c in sorted(
        claims, key=lambda c: (c.chrom, c.start, c.end, c.element, c.axis, c.value, c.gene, c.cell)
    ):
        obs = ev.at(c.chrom, c.start, c.end)
        a, b = co.verdict_of(c, obs, ra.rule), co.verdict_of(c, obs, rb.rule)
        if _pair(a) != _pair(b):
            out.append((a, b))
    return out


def unexplained_between(
    ra: co.Report, rb: co.Report, moved: list[tuple[co.Verdict, co.Verdict]]
) -> list[str]:
    """`unexplained` between the rules of any two reports: the first tally minus the moved claims' first
    verdicts plus their second must equal the second tally, count by count. Activity `refutable` is
    skipped only against v1, which never computed it (REFUTABLE_RULE)."""
    out = []
    for axis in co.AXES:
        ma, mb = co.AxisTally(rule=ra.rule), co.AxisTally(rule=rb.rule)
        for a, b in moved:
            if a.claim.axis == axis:
                ma.add(a, set())
                mb.add(b, set())
        ta, tb = _flat(ra.axes[axis].to_dict()), _flat(rb.axes[axis].to_dict())
        fa, fb = _flat(ma.to_dict()), _flat(mb.to_dict())
        for k in sorted(set(ta) | set(tb)):
            if axis == co.ACTIVITY and k == "refutable" and co.RULE_V1 in (ra.rule, rb.rule):
                continue
            want = ta.get(k, 0) - fa.get(k, 0) + fb.get(k, 0)
            if tb.get(k, 0) != want:
                out.append(
                    f"{axis}.{k}: {ra.rule} {ta.get(k, 0)}, {rb.rule} {tb.get(k, 0)}, "
                    f"moved claims account for {want}"
                )
    return out


def side_by_side_v3(reps: dict[str, co.Report]) -> dict[str, Any]:
    """Every axis's counts under each rule, as `side_by_side` writes them for two."""
    out: dict[str, Any] = {}
    for axis in co.AXES:
        out[axis] = {}
        for rule, rep in reps.items():
            t = rep.axes[axis]
            out[axis][rule] = {
                "claims": t.claims,
                "judged": t.judged(),
                **{v: t.verdicts[v] for v in co.VERDICTS},
                "not_judged_by_reason": {r: t.not_judged_by_reason[r] for r in co.REASONS_V2},
                "refutable": t.refutable if (axis != co.ACTIVITY or rule != co.RULE_V1) else None,
            }
    return out


def in_context_side_by_side(reps: dict[str, co.Report]) -> dict[str, Any]:
    """The target and direction counts in the owner's terms under each rule, from the verdicts
    (correctness.in_context_counts); under v1 a verdict decided with another cell is counted apart."""
    keys = (*co.IN_CONTEXT_TERMS, co.DECIDED_WITH_ANOTHER_CELL)
    out: dict[str, Any] = {"scope": co.IN_CONTEXT_SCOPE, "caution": co.IN_CONTEXT_CAUTION}
    for rule, rep in reps.items():
        rows = co.in_context_counts((*rep.judged, *rep.not_assessed))
        out[rule] = {axis: {k: rows[axis][k] for k in keys} for axis in (co.TARGET, co.ACTIVITY)}
    return out


def _moved_block(moved: list[tuple[co.Verdict, co.Verdict]], ra: str, rb: str, loose: list[str]) -> dict:
    return {
        "from": ra,
        "to": rb,
        "claims": len(moved),
        "by_axis": dict(Counter(a_.claim.axis for a_, _ in moved)),
        f"{ra}_verdicts": dict(Counter(a_.verdict for a_, _ in moved)),
        f"{ra}_decided_in": dict(Counter(_decided_in(a_) for a_, _ in moved)),
        f"{rb}_verdicts": dict(
            Counter(f"{b.verdict}:{b.reason}" if b.reason else b.verdict for _, b in moved)
        ),
        f"{rb}_cross_cell": dict(
            Counter("/".join(sorted({f for _, f in b.cross_cell or ()})) or "none" for _, b in moved)
        ),
        "unexplained_count_changes": loose,
    }


def main_v3(a: argparse.Namespace) -> None:
    t0 = time.time()
    units = ho.all_units()
    refs = ho.crispri_references()
    rep3 = judged_under(co.RULE_V3, a, units, refs)
    body = rep3.to_dict()
    show(body)
    t1 = time.time()
    rep1 = judged_under(co.RULE_V1, a, units, refs)
    rep2 = judged_under(co.RULE_V2, a, units, refs)
    registered = a.compiled == ho.COMPILED_DIR and not a.chroms
    paths = {co.RULE_V1: a.results / f"{NAME}.json", co.RULE_V2: a.results / f"{NAMES[co.RULE_V2]}.json"}
    reproduces: dict[str, bool | None] = {co.RULE_V1: None, co.RULE_V2: None}
    if registered:
        drop = set(RUN_KEYS)
        committed1 = json.loads(paths[co.RULE_V1].read_text())
        rerun1 = json.loads(json.dumps(v1_payload(rep1, units, t1), default=str))
        reproduces[co.RULE_V1] = _text(committed1, drop) == _text(rerun1, drop)
        committed2 = json.loads(paths[co.RULE_V2].read_text())
        body2 = json.loads(json.dumps(rep2.to_dict(), default=str))
        reproduces[co.RULE_V2] = (
            all(body2[k] == committed2[k] for k in JUDGED_KEYS)
            and json.loads(json.dumps(side_by_side(rep1, rep2))) == committed2["side_by_side"]
        )
        if not all(reproduces.values()):
            raise SystemExit(
                f"an earlier rule does not reproduce its committed file {reproduces}: nothing written"
            )
    moved = moved_between(rep2, rep3, units)
    loose = unexplained_between(rep2, rep3, moved)
    moved13 = moved_between(rep1, rep3, units)
    loose13 = unexplained_between(rep1, rep3, moved13)
    print("moved v2->v3", len(moved), "unexplained", loose, "| v1->v3", len(moved13), loose13, flush=True)
    if a.dry:
        return
    params = {**v1_params(a), "judge_rule": co.RULE_V3, "rule_registration": co.rule_registration_v3()}
    reps = {co.RULE_V1: rep1, co.RULE_V2: rep2, co.RULE_V3: rep3}
    payload = {
        "question": "the same claims and sources as S4, judged under v3, where the target claim is judged "
        "only in its stated cell as the direction is under v2: which counts move from the v2 result, and are "
        "the moves exactly the targets v1 and v2 established only from another cell",
        "rule": co.RULE_V3,
        "status": co.STATUS,
        "reuse": co.REUSE,
        "may_be_called": V3_MAY_BE_CALLED,
        "may_not_be_called": list(V3_MAY_NOT_BE_CALLED),
        "in_context_caution": co.IN_CONTEXT_CAUTION,
        "rule_registration": co.rule_registration_v3(),
        "expected_before_the_run": co.EXPECTED_V3,
        "seen_before_registration": co.SEEN_BEFORE_V3,
        "earlier_results": {
            co.RULE_V1: {
                "path": str(paths[co.RULE_V1]),
                "rule": co.RULE_V1,
                "kept": "unchanged, beside this result; S4's original reading stands as the v1 result",
                "rerun_reproduces_it": reproduces[co.RULE_V1],
                "compared": "the whole file, key by key, in order",
                "compared_without": list(RUN_KEYS),
            },
            co.RULE_V2: {
                "path": str(paths[co.RULE_V2]),
                "rule": co.RULE_V2,
                "kept": "unchanged, beside this result; the v2 re-judge stands as the v2 result",
                "rerun_reproduces_it": reproduces[co.RULE_V2],
                "compared": [*JUDGED_KEYS, "side_by_side"],
            },
        },
        "side_by_side": side_by_side_v3(reps),
        "in_context_side_by_side": in_context_side_by_side(reps),
        "moved": {
            **_moved_block(moved, co.RULE_V2, co.RULE_V3, loose),
            "each": [
                {
                    "claim": _claim_dict(a_.claim),
                    co.RULE_V2: {**a_.to_dict(), "decided_in": _decided_in(a_)},
                    co.RULE_V3: b.to_dict(),
                }
                for a_, b in moved
            ],
        },
        "moved_from_v1": _moved_block(moved13, co.RULE_V1, co.RULE_V3, loose13),
        "registration": co.registration(),
        "units_per_source": {s: len(units[s]) for s in ho.sources(units)},
        **body,
        "alphagenome_requests": 0,
        "per_element_response_cache_opened": False,
        "seconds": round(time.time() - t0, 1),
    }
    also = tuple(paths.values()) if registered else ()
    p = save_result(NAME_V3, payload, a.results, manifest=manifest(programs_of(a), params, also))
    print("saved", p, round(time.time() - t0, 1), "s", flush=True)


if __name__ == "__main__":
    main()
