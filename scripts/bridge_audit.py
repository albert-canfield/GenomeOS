# SPDX-License-Identifier: AGPL-3.0-or-later
"""The R3 double-counting audit over the compiled programs, answering the question registered in 06b7557.

    uv run python scripts/bridge_audit.py [--chroms chr21] [--no-save]

For every compiled program in data/knowledge/compiled and every cell context its rules are gated on,
it runs `attribution.bridge.parameterize` in classify-only mode and counts, per context, what a run in
that cell would integrate: mechanisms that carry both a predicted and a measured rule (which a plain
runtime run would integrate side by side), genes with more than one mechanism (not identifiable from
single-removal observations under the mean rule), experimental links whose strength is a maximum over
several pairs, and every other unresolved reason. No gene rates are supplied, so every mechanism that
survives the other checks is reported as `gene_parameter_missing`: the compiled genes are stubs.
Reads only the compiled text; no AlphaGenome request, no assay cache.
"""

from __future__ import annotations

import argparse
import time
from collections import Counter
from pathlib import Path

from genomeos import manifest as mf
from genomeos.attribution import bridge
from genomeos.evidence import COMPILED_DIR
from genomeos.lang.parser import parse
from genomeos.results import save_result

CHROMS = [f"chr{c}" for c in [*range(1, 23), "X", "Y"]]


def audit_module(m) -> dict:
    """Counts for one parsed program over every cell context its regulatory rules name."""
    contexts = sorted({r.when.get("cell_type") for r in m.rules if r.when.get("cell_type")})
    reasons: Counter = Counter()
    out = Counter()
    for cell in contexts:
        ctx = {"cell_type": cell}
        b = bridge.parameterize(m, ctx, build=False)
        reasons.update(i.reason for i in b.unresolved)
        # superseded by the bridge: only on genes it could consider, i.e. with one mechanism
        out["superseded_rules"] += len(b.superseded)
        by_key: dict = {}
        for r in m.active_rules(ctx):
            if r.action in (bridge.Action.ACTIVATE, bridge.Action.INHIBIT):
                by_key.setdefault(bridge.mechanism_key(r), []).append(r)
        out["active_regulatory_rules"] += sum(len(v) for v in by_key.values())
        # the audit's question, asked of the rules themselves: both forms of one relation in one run
        out["mechanisms_with_predicted_and_measured"] += sum(
            1 for v in by_key.values() if {"predicted", "experimental"} <= {bridge._kind(r) for r in v}
        )
        out["mechanisms"] += len(by_key)
        genes = Counter(k[1] for k in by_key)
        out["gene_cell_pairs"] += len(genes)
        out["gene_cell_pairs_with_several_mechanisms"] += sum(1 for n in genes.values() if n > 1)
        out["mechanisms_on_those_genes"] += sum(n for n in genes.values() if n > 1)
    out["contexts"] = len(contexts)
    return {"counts": dict(out), "unresolved": dict(reasons)}


def manifest(chroms: list[str]) -> dict:
    """The provenance contract (R9): the compiled programs read, pinned by sha256."""
    paths = [Path(COMPILED_DIR) / f"noncoding_{c}.bio" for c in chroms]
    return {
        "sources": [
            {
                "accession": "the compiled non-coding programs, scripts/compile_genome_programs.py",
                "version": "as compiled from this revision's attribution/compile.py, pinned by sha256",
            }
        ],
        "inputs": [mf.input_entry(p, partition=None) for p in paths if p.exists()],
        "assembly": "GRCh38",
        "coordinates": "n/a: counts of rules and mechanisms, no interval is reported",
        "parameters": {
            "precedence": list(bridge.PRECEDENCE),
            "citation_tolerance": bridge.CITATION_TOLERANCE,
            "experimental_summary": bridge.EXPERIMENTAL_SUMMARY,
            "gene_rates_supplied": "none: the compiled genes are stubs",
        },
        "exclusions": [],
        "partitions": "n/a: the audit counts rules as compiled; held-out rules are counted with the rest",
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--chroms", default="")
    ap.add_argument("--no-save", action="store_true")
    args = ap.parse_args(argv)
    chroms = args.chroms.split(",") if args.chroms else CHROMS
    t0 = time.time()
    per: dict[str, dict] = {}
    for chrom in chroms:
        path = Path(COMPILED_DIR) / f"noncoding_{chrom}.bio"
        if not path.exists():
            per[chrom] = {"missing": str(path)}
            continue
        per[chrom] = audit_module(parse(path.read_text(), path.stem))
        print(chrom, per[chrom], flush=True)
    pooled: Counter = Counter()
    reasons: Counter = Counter()
    for v in per.values():
        pooled.update(v.get("counts", {}))
        reasons.update(v.get("unresolved", {}))
    out = {
        "result": "bridge_audit",
        "review_item": "R3",
        "question": bridge.AUDIT_QUESTION,
        "registered_in": "06b7557 (amended b986238)",
        "pooled": dict(pooled),
        "unresolved_by_reason": dict(reasons),
        "per_chromosome": per,
        "reading": (
            "a context is a cell_type a compiled rule is gated on; one rule can count once per program"
            " because each rule names one context. gene_parameter_missing counts mechanisms that pass"
            " every other check and fail only because the compiled genes are stubs"
        ),
        "seconds": round(time.time() - t0, 1),
    }
    if not args.no_save:
        print("saved", save_result("bridge_audit", out, manifest=manifest(chroms)))
    print({"pooled": out["pooled"], "unresolved": out["unresolved_by_reason"]})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
