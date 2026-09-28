# SPDX-License-Identifier: AGPL-3.0-or-later
"""The R3 double-counting audit over the compiled programs, answering the question registered in 06b7557.

    uv run python scripts/bridge_audit.py [--chroms chr21] [--no-save]

For every compiled program in data/knowledge/compiled and every cell context its rules are gated on,
it runs `attribution.bridge.parameterize` in classify-only mode and counts, per context, what a run in
that cell would integrate: mechanisms that carry both a predicted and a measured rule (which a plain
runtime run would integrate side by side), genes with more than one mechanism (not identifiable from
single-removal observations under the mean rule), experimental links whose strength is a maximum over
several pairs, and every other unresolved reason.

Since 2026-09-28 it runs each chromosome under three rate tiers and reports all three side by side
(the registration is in ff2062a and docs/DESIGN-MINIMAL-CELL.md):

  none      no gene rates at all, the state before this work: every mechanism that survives the other
            checks fails on `gene_parameter_missing`, because the compiled genes are stubs
  measured  the gene's own measured transcription rate from data/results/gene_rates.json, carried from
            MOUSE NIH 3T3 fibroblasts to a human symbol. This is the headline count, and every gene in
            it is borrowing a mouse constant
  borrowed  the genome median of those rates for a gene with no measurement of its own. It fixes no
            absolute number: such a pair is simulable in RELATIVE units only and is never added to the
            measured count

Reads only the compiled text and the rate table; no AlphaGenome request, no assay cache.
"""

from __future__ import annotations

import argparse
import json
import time
from collections import Counter
from pathlib import Path

from genomeos import manifest as mf
from genomeos.attribution import bridge
from genomeos.evidence import COMPILED_DIR
from genomeos.lang.parser import parse
from genomeos.results import RESULTS_DIR, save_result

CHROMS = [f"chr{c}" for c in [*range(1, 23), "X", "Y"]]
RATE_TABLE = RESULTS_DIR / "gene_rates.json"
TIERS = ("none", "measured", "borrowed")


def load_rates() -> tuple[dict[str, float], float, dict[str, str]]:
    """gene symbol -> measured T, the median a borrowing gene is lent, and each gene's species+cell.

    The third return value is the provenance registered in `bridge.RATE_ROW_PROVENANCE`: it lets the
    audit count simulable pairs by the SPECIES of the rate they used, so a pair standing on a mouse
    constant is never reported as a human result. With one source it is a constant; it is read per row
    so that a second source cannot be pooled into it silently.
    """
    d = json.loads(RATE_TABLE.read_text())
    rates = {g: float(r["transcription_rate"]) for g, r in d["rates"].items()}
    species = {g: str(r.get("species_cell", "unstated")) for g, r in d["rates"].items()}
    return rates, float(d["borrowed_median_transcription_rate"]), species


def audit_module(
    m,
    rates: dict[str, float] | None = None,
    borrowed: float | None = None,
    species: dict[str, str] | None = None,
) -> dict:
    """Counts for one parsed program over every cell context its regulatory rules name.

    The shape counts (rules, mechanisms, gene-cell pairs) do not depend on the rates and are reported
    once; the resolution counts are reported per tier.
    """
    contexts = sorted({r.when.get("cell_type") for r in m.rules if r.when.get("cell_type")})
    out = Counter()
    per_tier: dict[str, Counter] = {t: Counter() for t in TIERS}
    reasons_by_tier: dict[str, Counter] = {t: Counter() for t in TIERS}
    for cell in contexts:
        ctx = {"cell_type": cell}
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
        for tier in TIERS:
            kw = {}
            if tier == "measured":
                kw = {"rates": rates}
            elif tier == "borrowed":
                kw = {"rates": rates, "borrowed_rate": borrowed}
            b = bridge.parameterize(m, ctx, build=False, **kw)
            reasons_by_tier[tier].update(i.reason for i in b.unresolved)
            c = per_tier[tier]
            # superseded by the bridge: only on genes it could consider, i.e. with one mechanism
            c["superseded_rules"] += len(b.superseded)
            c["simulable_mechanisms"] += len(b.strengths)
            # one mechanism per gene reaches the fit, so a simulable mechanism IS a simulable pair
            c["simulable_gene_cell_pairs"] += len({k[1] for k in b.strengths})
            c["simulable_compiled_rules"] += sum(len(by_key.get(k, ())) for k in b.strengths)
            c["genes_on_their_own_measured_rate"] += sum(1 for t in b.tier.values() if t == "measured")
            c["genes_on_a_borrowed_median"] += sum(1 for t in b.tier.values() if t == "borrowed_median")
            # the species of the rate each simulable pair stands on. A borrowed median is not counted
            # here at all: it fixes no absolute number, so it has no species to report
            sp = species or {}
            c["pairs_on_a_human_measured_rate"] += sum(
                1 for g, t in b.tier.items() if t == "measured" and sp.get(g, "").startswith("Homo")
            )
            c["pairs_on_a_borrowed_species_rate"] += sum(
                1 for g, t in b.tier.items() if t == "measured" and not sp.get(g, "").startswith("Homo")
            )
            # the denominator of the registered inhibitor bound: only an inhibitory mechanism can fail
            # it (an activator's split always leaves a positive basal rate), so the ratio it is judged
            # on is rate_split_infeasible over the inhibitory mechanisms that got as far as the split
            c["inhibitory_mechanisms_at_the_split"] += sum(
                1 for ob in b.observations.values() if ob.rho > 1.0
            )
    out["contexts"] = len(contexts)
    return {
        "counts": dict(out),
        "by_tier": {t: dict(per_tier[t]) for t in TIERS},
        "unresolved_by_tier": {t: dict(reasons_by_tier[t]) for t in TIERS},
        "unresolved": dict(reasons_by_tier["none"]),  # the historical field: the no-rates tier
    }


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
        "inputs": [mf.input_entry(p, partition=None) for p in paths if p.exists()]
        + ([mf.input_entry(RATE_TABLE, partition=None)] if RATE_TABLE.exists() else []),
        "assembly": "GRCh38",
        "coordinates": "n/a: counts of rules and mechanisms, no interval is reported",
        "parameters": {
            "precedence": list(bridge.PRECEDENCE),
            "citation_tolerance": bridge.CITATION_TOLERANCE,
            "experimental_summary": bridge.EXPERIMENTAL_SUMMARY,
            "gene_rates_supplied": (
                "three tiers side by side: none (the compiled genes are stubs), measured"
                f" (data/results/gene_rates.json), borrowed (its median). {bridge.RATE_SPECIES}"
            ),
            "rate_tiers": bridge.RATE_TIERS,
            "rate_row_provenance": bridge.RATE_ROW_PROVENANCE,
            "human_rate_survey": bridge.HUMAN_RATE_SURVEY_VERDICT,
            "rate_units": bridge.MEASURED_RATE_UNITS,
            "declared_strength": bridge.DECLARED_STRENGTH,
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
    rates, borrowed, species = load_rates()
    per: dict[str, dict] = {}
    for chrom in chroms:
        path = Path(COMPILED_DIR) / f"noncoding_{chrom}.bio"
        if not path.exists():
            per[chrom] = {"missing": str(path)}
            continue
        per[chrom] = audit_module(
            parse(path.read_text(), path.stem), rates=rates, borrowed=borrowed, species=species
        )
        print(chrom, per[chrom]["by_tier"], flush=True)
    pooled: Counter = Counter()
    reasons: Counter = Counter()
    tiers: dict[str, Counter] = {t: Counter() for t in TIERS}
    tier_reasons: dict[str, Counter] = {t: Counter() for t in TIERS}
    for v in per.values():
        pooled.update(v.get("counts", {}))
        reasons.update(v.get("unresolved", {}))
        for t in TIERS:
            tiers[t].update(v.get("by_tier", {}).get(t, {}))
            tier_reasons[t].update(v.get("unresolved_by_tier", {}).get(t, {}))
    out = {
        "result": "bridge_audit",
        "review_item": "R3",
        "question": bridge.AUDIT_QUESTION,
        "registered_in": "06b7557 (amended b986238)",
        "pooled": dict(pooled),
        "unresolved_by_reason": dict(reasons),
        "rate_tiers": bridge.RATE_TIERS,
        "rate_species": bridge.RATE_SPECIES,
        "rate_expectation": bridge.RATE_EXPECTED,
        "borrowed_median_transcription_rate": borrowed,
        "genes_with_a_measured_rate": len(rates),
        "genes_with_a_human_measured_rate": sum(1 for v in species.values() if v.startswith("Homo")),
        "rate_row_provenance": bridge.RATE_ROW_PROVENANCE,
        "human_rate_survey_verdict": bridge.HUMAN_RATE_SURVEY_VERDICT,
        "human_rate_survey_falsifier": bridge.HUMAN_RATE_SURVEY_FALSIFIER,
        "pooled_by_tier": {t: dict(tiers[t]) for t in TIERS},
        "unresolved_by_tier": {t: dict(tier_reasons[t]) for t in TIERS},
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
    print(
        {"pooled": out["pooled"], "by_tier": out["pooled_by_tier"], "unresolved": out["unresolved_by_tier"]}
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
