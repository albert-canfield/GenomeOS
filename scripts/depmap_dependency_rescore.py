# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""Re-score the therapeutic benchmark with the deep-deletion dependency wiring off and on.

The measurement the registration in docs/THERAPEUTICS.md ("The dependency a deep
deletion creates", 2026-10-10) asked for, offline and from the local caches:

1. the ten benchmark cases (scripts/therapeutic_benchmark.py, `run_case`), each run
   with `dependency.assess` replaced by an empty answer and then as built, every row
   field compared between the two and against the committed
   data/results/therapeutic_benchmark_v2.json, and the dependency list of each run;
2. the one fixture that carries a deep deletion, the CDKN2A glioma in
   data/demo/alterations, off and on, with no lineage stated and with CNS/Brain;
3. the demo tumour, off and on.

For every run with the wiring off and on, each target candidate's gene, overall
score, evidence tier, mechanism reach, class and position are compared.

It writes nothing under data/: the comparison goes to stdout and, with --out, to a
JSON file the caller names. The committed benchmark result is not rewritten.

Run: `uv run python scripts/depmap_dependency_rescore.py [--out FILE]`
"""

from __future__ import annotations

import importlib.util
import json
import sys
import tempfile
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos.therapeutics import analyse_vcf  # noqa: E402
from genomeos.therapeutics import dependency as dep  # noqa: E402

BENCHMARK = Path(__file__).with_name("therapeutic_benchmark.py")
COMMITTED = Path("data/results/therapeutic_benchmark_v2.json")
FIXTURE = Path("data/demo/alterations")
DEMO = Path("data/demo/cancer_tumour.vcf")


def _bench():
    spec = importlib.util.spec_from_file_location("therapeutic_benchmark", BENCHMARK)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


class _Wiring:
    """`dependency.assess` as built, or replaced by the empty answer."""

    def __init__(self) -> None:
        self.built = dep.assess

    def __call__(self, on: bool) -> None:
        dep.assess = self.built if on else (lambda *a, **k: dep.empty())


def targets(analysis: dict[str, Any]) -> list[dict[str, Any]]:
    return [
        {
            "rank": i,
            "gene": c.gene,
            "overall": c.scores.overall,
            "evidence_tier": c.evidence_tier,
            "mechanism_reach": c.mechanism_reach,
            "target_class": c.target_class,
            "best_mechanism": c.best_mechanism.mechanism if c.best_mechanism else None,
        }
        for i, c in enumerate(analysis["candidates"], 1)
    ]


def dependency_rows(analysis: dict[str, Any]) -> list[dict[str, Any]]:
    out = []
    for c in analysis.get("dependency_candidates") or []:
        t = c["lineage"]
        out.append(
            {
                "gene": c["gene"],
                "lost_gene": c["lost_gene"],
                "verdict": t["verdict"],
                "lineage": t.get("name"),
                "reason": t.get("reason"),
                **{
                    k: t.get(k)
                    for k in (
                        "lines_lost",
                        "lines_intact",
                        "median_effect_lost",
                        "median_effect_intact",
                        "cliffs_delta",
                        "p",
                        "q",
                        "tests_in_this_analysis",
                    )
                    if k in t
                },
                "sweep_q": c["sweep"]["q"],
                "reading": c["reading"],
            }
        )
    return out


def run_benchmark(bench, wiring: _Wiring, vep_dir: Path) -> dict[str, Any]:
    captured: list[dict[str, Any]] = []
    built = bench.analyse_vcf

    def capture(*a, **k):
        analysis = built(*a, **k)
        captured.append(analysis)
        return analysis

    bench.analyse_vcf = capture
    committed = {r["case"]: r for r in json.loads(COMMITTED.read_text())["rows"]}
    cases = []
    try:
        for case in bench.CASES:
            runs = {}
            for on in (False, True):
                wiring(on)
                captured.clear()
                row = bench.run_case(case, False, None, vep_dir)
                row.pop("seconds", None)
                runs[on] = (row, captured[-1])
            (off_row, off_a), (on_row, on_a) = runs[False], runs[True]
            ref = {k: v for k, v in committed[case["case"]].items() if k != "seconds"}
            cases.append(
                {
                    "case": case["case"],
                    "gene": case["gene"],
                    "rank": on_row["rank"],
                    "score": on_row.get("score"),
                    "verdict": on_row.get("verdict"),
                    "row_fields_moved_off_to_on": sorted(k for k in on_row if on_row[k] != off_row.get(k)),
                    "row_fields_differing_from_committed": sorted(
                        k for k in set(on_row) | set(ref) if on_row.get(k) != ref.get(k)
                    ),
                    "targets_identical_off_on": targets(off_a) == targets(on_a),
                    "deep_deletions": [r["gene"] for r in on_a.get("dependencies") or []],
                    "dependency_candidates": dependency_rows(on_a),
                }
            )
    finally:
        bench.analyse_vcf = built
        wiring(True)
    return {"cases": cases}


def run_tumour(
    name: str, wiring: _Wiring, vep_dir: Path, lineages: tuple[str, ...], **inputs
) -> dict[str, Any]:
    out: dict[str, Any] = {"input": name, "runs": []}
    for lineage in lineages:
        runs = {}
        for on in (False, True):
            wiring(on)
            runs[on] = analyse_vcf(
                net=False, top_genes=8, indirect=True, vep_append_to=vep_dir, lineage=lineage, **inputs
            )
        wiring(True)
        out["runs"].append(
            {
                "lineage": lineage or None,
                "targets_identical_off_on": targets(runs[False]) == targets(runs[True]),
                "targets": targets(runs[True]),
                "dependencies": [
                    {k: r[k] for k in ("gene", "event", "status", "partners", "same_event")}
                    for r in runs[True].get("dependencies") or []
                ],
                "dependency_candidates": dependency_rows(runs[True]),
                "dependency_candidates_with_wiring_off": len(runs[False].get("dependency_candidates") or []),
            }
        )
    return out


def main() -> int:
    out_path = sys.argv[sys.argv.index("--out") + 1] if "--out" in sys.argv else None
    bench = _bench()
    wiring = _Wiring()
    with tempfile.TemporaryDirectory(prefix="genomeos-vep-appends-") as tmp:
        vep_dir = Path(tmp)
        result = {
            "registration": "docs/THERAPEUTICS.md, 'The dependency a deep deletion creates, wired into the "
            "evidence, pre-registered 2026-10-10' (c2f1fa9)",
            "network": False,
            "benchmark": run_benchmark(bench, wiring, vep_dir),
            "cdkn2a_glioma_fixture": run_tumour(
                "data/demo/alterations/cdkn2a_deleted.vcf + .cnv",
                wiring,
                vep_dir,
                ("", "CNS/Brain"),
                tumour_vcf=str(FIXTURE / "cdkn2a_deleted.vcf"),
                cnv=str(FIXTURE / "cdkn2a_deleted.cnv"),
            ),
            "demo_tumour": run_tumour(str(DEMO), wiring, vep_dir, ("",), tumour_vcf=str(DEMO)),
        }
    for c in result["benchmark"]["cases"]:
        print(
            f"{c['gene']:7} rank {c['rank']} score {c['score']} | moved off->on: "
            f"{c['row_fields_moved_off_to_on'] or 'none'} | vs committed: "
            f"{c['row_fields_differing_from_committed'] or 'none'} | targets identical "
            f"{c['targets_identical_off_on']} | deep deletions {c['deep_deletions'] or 'none'} | "
            f"dependency candidates {len(c['dependency_candidates'])}"
        )
    for key in ("cdkn2a_glioma_fixture", "demo_tumour"):
        for r in result[key]["runs"]:
            print(
                f"{key} lineage={r['lineage']}: targets identical {r['targets_identical_off_on']}; "
                f"targets {[t['gene'] for t in r['targets']]}; off-run dependency candidates "
                f"{r['dependency_candidates_with_wiring_off']}"
            )
            for d in r["dependency_candidates"]:
                stats = {k: d.get(k) for k in ("lines_lost", "lines_intact", "cliffs_delta", "p", "q")}
                print(f"  {d['gene']} from {d['lost_gene']}: {d['verdict']} {stats}")
                print(f"    {d['reading']}")
    if out_path:
        Path(out_path).write_text(json.dumps(result, indent=1, default=str))
        print(f"written: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
