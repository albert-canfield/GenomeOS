# SPDX-License-Identifier: AGPL-3.0-or-later
"""Item 12 S4, applied once: the unchanged compiled labels judged claim by claim against every held-out
source, with target accuracy, role accuracy and coverage reported apart, per axis, with denominators.

    uv run python scripts/s4_correctness_run.py              # v1, S4's result as committed
    uv run python scripts/s4_correctness_run.py --rule v2    # the same claims and sources under v2

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


if __name__ == "__main__":
    main()
