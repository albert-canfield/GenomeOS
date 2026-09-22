# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""Marker logic: what a binder must see, and what it must not, as a Boolean gate.

One antigen is rarely selective. CD20 sits on healthy spleen at 247 nTPM and
`--scan` turns it down at 2.2x for exactly that, the target of the most
successful antibody in oncology failing a selectivity screen. That refusal is
correct about the single marker and wrong about the target space, because a cell
can be addressed by a *combination*: two
antigens that must both be present, or one that must be present while another,
carried by the healthy tissue at risk, must be absent. Bispecific binders,
logic-gated CARs and dual-antigen engagers are all built on exactly that, and a
pipeline that only ever scores one gene at a time cannot propose any of them.

So this module scores gates rather than genes:

```
A            one marker, the baseline every combination has to beat
A AND B      both markers on the same cell
A AND NOT C  the marker, with a healthy-tissue marker that vetoes the cell
```

Two numbers are kept apart and never merged into one measurement:

* **tumour coverage** - the fraction of tumour cells the gate would address,
  measured in DepMap cancer cell lines of one lineage, each marker called
  present above a declared log2(TPM+1) threshold;
* **healthy burden** - the highest level the gate reaches in the essential
  healthy tissues of the packaged Human Protein Atlas atlas (the eight
  `safety.ESSENTIAL` organs), computed per tissue: an AND gate falls towards
  the *lower* of its markers but never below a quarter of the higher one,
  because bulk tissue RNA cannot say the two markers are on different cells,
  and a NOT gate is vetoed outright in a tissue that carries its exclusion
  marker, because a measured presence is evidence in a way an absence in 20
  bulk tissues is not.

The score is `coverage / (1 + burden / 10 nTPM)`: coverage paid for by what the
gate would hit in tissue that cannot be replaced. It is a sort order, not a
verdict, exactly as elsewhere in this module.

Cell lines DepMap itself calls non-cancerous are dropped before anything is
scored. A first pass kept them, and 42 fibroblast lines then supplied "tumour
coverage" for gates that had nothing to do with cancer; that pass was discarded
without its benchmark verdict being read.

`PREREGISTERED_CLAIM` is the falsification test, written before the first
ranking was read: combinations containing an approved or clinical surface
target must rank above random combinations drawn from the same pool. If the
approved targets of oncology come out indistinguishable from random pairs of
membrane proteins, the score is measuring nothing.

**No design follows.** The output names marker logic and nothing else: no
binder, no construct, no format, no cargo.
"""

from __future__ import annotations

import random
from pathlib import Path
from typing import Any

from . import atlas, depmap
from .evidence import Evidence, prediction
from .expression import CONCERN_NTPM
from .providers import CRITICAL_TISSUES, hpa_column
from .safety import ESSENTIAL
from .synthetic_lethality import median, rank_sum_p, ranks

#: A marker counts as present on a line at or above this log2(TPM+1), i.e. about
#: 7 TPM. RNA is not surface protein and this threshold is a stated policy.
PRESENT_LOG2TPM = 3.0
#: A NOT marker vetoes a healthy tissue at or above this nTPM in the atlas.
PROTECT_NTPM = 10.0
#: How much of a tissue's worst marker an AND gate is still assumed to reach.
#: Bulk tissue RNA cannot say whether two markers sit on the *same* healthy
#: cell, so requiring a second marker is credited with at most a four-fold
#: reduction in exposure and never with abolishing it. Without this floor, any
#: gate paired with a transcript the atlas reports as absent scored a burden of
#: zero, and thousands of gates tied at the top of the ranking on an absence
#: that 20 bulk tissues cannot establish.
CO_EXPRESSION_FLOOR = 0.25
#: Cell lines DepMap itself calls non-cancerous are not tumour coverage.
EXCLUDE_DISEASE = "Non-Cancerous"
#: Lineages below this many screened lines are not scored: a coverage over six
#: cell lines is not a coverage.
MIN_LINES = 15
#: How many single markers enter the combination pool, plus the benchmark genes.
POOL_SIZE = 250
#: Placebo markers, forced into the pool exactly as the benchmark genes are, so
#: that the null is matched on the route into the pool and not only on set size.
PLACEBO_POOL_SIZE = 64
#: Thresholds the post-hoc sensitivity analysis re-runs the claim at. Declared
#: as post-hoc: the pre-registered verdict is the one at PRESENT_LOG2TPM.
SENSITIVITY_THRESHOLDS = (4.0, 5.0, 6.0)
SENSITIVITY_POOL = 150
#: The random comparison set for the pre-registered claim, and the smallest
#: effect that counts as the benchmark set ranking above it.
RANDOM_COMBINATIONS = 2000
MIN_BENCHMARK_DELTA = 0.10
SEED = 20260916

#: Approved or clinical surface targets, with the therapy that made them one.
#: This is the benchmark set: it is fixed here, before any ranking was looked
#: at, and every one of them is a target a regulator or a trial has already
#: accepted as reachable.
BENCHMARK: dict[str, str] = {
    "ERBB2": "trastuzumab, pertuzumab, trastuzumab deruxtecan (approved)",
    "EGFR": "cetuximab, panitumumab (approved)",
    "CD19": "blinatumomab, tisagenlecleucel (approved)",
    "MS4A1": "rituximab, obinutuzumab (approved)",
    "CD22": "inotuzumab ozogamicin, moxetumomab (approved)",
    "TNFRSF17": "belantamab, teclistamab, idecabtagene (approved; BCMA)",
    "CD38": "daratumumab, isatuximab (approved)",
    "CD33": "gemtuzumab ozogamicin (approved)",
    "TACSTD2": "sacituzumab govitecan (approved; TROP2)",
    "FOLH1": "lutetium-177 PSMA-617 (approved; PSMA)",
    "MSLN": "mesothelin ADCs and CAR T (clinical)",
    "CEACAM5": "tusamitamab ravtansine (clinical)",
    "DLL3": "tarlatamab (approved; small-cell lung)",
    "CLDN18": "zolbetuximab (approved; claudin-18.2)",
    "GPC3": "glypican-3 CAR T (clinical)",
    "CD70": "CD70 CAR T and ADCs (clinical)",
}

PREREGISTERED_CLAIM = (
    "Marker combinations containing at least one approved or clinical surface target (the BENCHMARK "
    f"set of {len(BENCHMARK)} genes) score above random combinations drawn from the same pool of "
    "surface proteins. Declared in advance: supported only if the one-sided rank-test p is at or "
    f"below 0.05, Cliff's delta is at least {MIN_BENCHMARK_DELTA}, and that delta is above the "
    "95th percentile of a placebo null in which gene sets of the same size are drawn at random "
    "from the same pool. A comparison set is not neutral because it was called random: the null "
    "says what the comparison reads when the labels carry no information, and the result is read "
    "against that number rather than against zero."
)
#: How many placebo gene sets the null draws.
NULL_GENE_SETS = 200

LIMITS: tuple[str, ...] = (
    "tumour coverage is measured in cancer cell lines, which have no stroma, no immune "
    "infiltrate and no tissue architecture, and are a growth-selected subset of the disease",
    "both sides are RNA: a transcript is not a protein and a protein is not a protein on the "
    "outward-facing surface",
    "the two sides are different measurements on different scales - log2(TPM+1) in cell lines "
    "against consensus nTPM in healthy tissue - so coverage and burden are reported separately "
    "and never subtracted from one another",
    "an AND gate's healthy burden is floored at a quarter of its highest marker in each tissue: "
    "bulk tissue RNA cannot say whether two markers sit on the same healthy cell, so requiring a "
    "second marker is credited with reducing exposure and never with abolishing it",
    "the atlas covers 20 tissues and 8 of them are weighed as essential here: a tissue not "
    "queried is not a tissue where the marker is absent",
)


def essential_columns() -> dict[str, str]:
    """Essential healthy tissue -> the atlas column that carries its nTPM."""
    return {name: hpa_column(field) for field, (name, _w) in CRITICAL_TISSUES.items() if name in ESSENTIAL}


def healthy_levels(genes: set[str]) -> dict[str, dict[str, float]]:
    """Gene -> essential tissue -> nTPM, from the packaged atlas."""
    columns = essential_columns()
    out: dict[str, dict[str, float]] = {}
    for gene in genes:
        row = atlas.row(gene)
        if not row:
            continue
        levels = {t: float(row[c]) for t, c in columns.items() if isinstance(row.get(c), int | float)}
        if levels:
            out[gene] = levels
    return out


def surface_universe() -> list[str]:
    """Every gene in the packaged membrane and CD-marker atlas."""
    return atlas.universe()


class Gates:
    """Presence bitmaps per gene and per lineage, so a gate is two integer ops."""

    def __init__(
        self,
        expression: dict[str, dict[str, float]],
        lineage_of: dict[str, str],
        threshold: float = PRESENT_LOG2TPM,
        min_lines: int = MIN_LINES,
    ) -> None:
        lines = sorted(m for m in expression if lineage_of.get(m))
        self.lines = lines
        self.threshold = threshold
        self.bit = {m: i for i, m in enumerate(lines)}
        width = (len(lines) + 7) // 8
        bits: dict[str, bytearray] = {}
        for i, model in enumerate(lines):
            for gene, value in expression[model].items():
                if value >= threshold:
                    buffer = bits.get(gene)
                    if buffer is None:
                        buffer = bits[gene] = bytearray(width)
                    buffer[i >> 3] |= 1 << (i & 7)
        self.present: dict[str, int] = {
            gene: int.from_bytes(bytes(buffer), "little") for gene, buffer in bits.items()
        }
        counts: dict[str, int] = {}
        masks: dict[str, int] = {}
        for i, model in enumerate(lines):
            lineage = lineage_of[model]
            counts[lineage] = counts.get(lineage, 0) + 1
            masks[lineage] = masks.get(lineage, 0) | (1 << i)
        self.lineages = {k: masks[k] for k, n in counts.items() if n >= min_lines}
        self.lineage_size = {k: counts[k] for k in self.lineages}

    def mask(self, required: tuple[str, ...], excluded: tuple[str, ...] = ()) -> int:
        out = (1 << len(self.lines)) - 1
        for gene in required:
            out &= self.present.get(gene, 0)
        for gene in excluded:
            out &= ~self.present.get(gene, 0)
        return out

    def coverage(self, mask: int, lineage: str) -> float:
        hit = (mask & self.lineages[lineage]).bit_count()
        return hit / self.lineage_size[lineage]


def burden(
    required: tuple[str, ...],
    excluded: tuple[str, ...],
    levels: dict[str, dict[str, float]],
    protect_ntpm: float = PROTECT_NTPM,
    co_expression_floor: float = CO_EXPRESSION_FLOOR,
) -> tuple[float, str | None]:
    """The worst essential healthy tissue the gate still reaches, and which one.

    Per tissue, an AND gate can only hit a cell carrying every required marker,
    so its exposure falls towards the lowest of them - but bulk tissue RNA
    cannot say whether the two markers sit on the same cell, so the fall is
    floored at a quarter of the highest marker rather than taken all the way to
    the lowest. A NOT marker the tissue carries above `protect_ntpm` vetoes that
    tissue outright, because a measured presence is evidence in a way an absence
    in 20 bulk tissues is not. The maximum over tissues is what a binder would
    have to live with.
    """
    tissues = set()
    for gene in (*required, *excluded):
        tissues |= set(levels.get(gene, {}))
    worst = 0.0
    where: str | None = None
    for tissue in tissues:
        values = [levels.get(g, {}).get(tissue) for g in required]
        if not values or any(v is None for v in values):
            continue
        if any((levels.get(g, {}).get(tissue) or 0.0) >= protect_ntpm for g in excluded):
            continue
        present = [float(v) for v in values if v is not None]
        exposure = max(min(present), co_expression_floor * max(present))
        if exposure > worst:
            worst, where = exposure, tissue
    return worst, where


def score(coverage: float, healthy: float, concern: float = CONCERN_NTPM) -> float:
    """Coverage paid for by the essential-tissue burden. A sort order, not a verdict."""
    return coverage / (1.0 + healthy / concern)


Gate = tuple[tuple[str, ...], tuple[str, ...], str]


def single_gates(genes: list[str]) -> list[Gate]:
    """One gate per marker: the baseline every combination has to beat."""
    return [((g,), (), "A") for g in genes]


def combinations(pool: list[str]) -> list[Gate]:
    """The gate space: every single marker, every AND pair, every AND-NOT pair."""
    out: list[Gate] = single_gates(pool)
    for i, a in enumerate(pool):
        for b in pool[i + 1 :]:
            out.append(((a, b), (), "A and B"))
    for a in pool:
        for c in pool:
            if a != c:
                out.append(((a,), (c,), "A and not C"))
    return out


def rank_combinations(
    gates: Gates,
    levels: dict[str, dict[str, float]],
    gate_space: list[Gate],
    min_coverage: float = 0.0,
) -> list[dict[str, Any]]:
    """Score every gate in every lineage and keep each gate's best lineage."""
    rows: list[dict[str, Any]] = []
    for required, excluded, kind in gate_space:
        mask = gates.mask(required, excluded)
        if mask == 0:
            continue
        healthy, where = burden(required, excluded, levels)
        best: dict[str, Any] | None = None
        for lineage in gates.lineages:
            coverage = gates.coverage(mask, lineage)
            if coverage <= min_coverage:
                continue
            value = score(coverage, healthy)
            if best is None or value > best["score"]:
                best = {
                    "required": list(required),
                    "excluded": list(excluded),
                    "logic": kind,
                    "lineage": lineage,
                    "lines": gates.lineage_size[lineage],
                    "coverage": round(coverage, 3),
                    "healthy_burden_ntpm": round(healthy, 1),
                    "worst_essential_tissue": where,
                    "score": round(value, 4),
                }
        if best is not None:
            rows.append(best)
    rows.sort(key=lambda r: -r["score"])
    return rows


def evidence_for(row: dict[str, Any]) -> Evidence:
    """A ranked gate is a prediction about marker logic, capped accordingly.

    RNA in cell lines against RNA in healthy tissue cannot establish surface
    protein on either side, so confidence stops at 0.5 however good the score.
    """
    gate = " and not ".join([" and ".join(row["required"]), *row["excluded"]])
    confidence = round(min(0.5, 0.15 + 0.7 * row["score"]), 3)
    return prediction(
        "GenomeOS marker-selectivity screen (DepMap cell-line RNA against the packaged HPA atlas)",
        f"{gate} addresses {row['coverage']:.0%} of {row['lineage']} lines while the worst essential "
        f"healthy tissue carries {row['healthy_burden_ntpm']:g} nTPM"
        + (f" ({row['worst_essential_tissue']})" if row["worst_essential_tissue"] else ""),
        confidence,
    )


def _advantage(
    rows: list[dict[str, Any]],
    genes: set[str],
    rng: random.Random,
    sample: int,
) -> tuple[float, float, list[int], list[int]]:
    """One rank test: gates containing `genes` against a random draw of the rest."""
    marked = [i for i, r in enumerate(rows) if set(r["required"]) & genes]
    others = [i for i, r in enumerate(rows) if not (set(r["required"]) & genes)]
    if not marked or not others:
        return 1.0, 0.0, marked, []
    drawn = sorted(rng.sample(others, min(sample, len(others))))
    values = {f"m{i}": -rows[i]["score"] for i in marked} | {f"r{i}": -rows[i]["score"] for i in drawn}
    rank_of, ties = ranks(values)
    p, delta, _u = rank_sum_p(rank_of, ties, [f"m{i}" for i in marked])
    return p, -delta, marked, drawn


def placebo_null(
    rows: list[dict[str, Any]],
    candidates: list[str],
    size: int,
    draws: int = NULL_GENE_SETS,
    sample: int = RANDOM_COMBINATIONS,
    seed: int = SEED,
) -> dict[str, Any]:
    """The same test with gene sets drawn at random: what the comparison reads on noise.

    A comparison set is not neutral because it was called random. Gates
    containing *any* fixed set of genes are not a random sample of gates - they
    share a pool, a lineage structure and a coverage distribution - so the
    advantage of the benchmark set has to be read against the advantage a
    meaningless set of the same size gets for free. The candidates are the
    placebo markers, which entered the pool by the same forced route as the
    benchmark genes rather than by scoring their way in, so the null is matched
    on the one thing that otherwise separates the two sets.
    """
    rng = random.Random(seed + 1)
    deltas: list[float] = []
    for _ in range(draws):
        if len(candidates) < size:
            break
        genes = set(rng.sample(candidates, size))
        _p, delta, _m, _d = _advantage(rows, genes, rng, sample)
        deltas.append(delta)
    ordered = sorted(deltas)
    return {
        "gene_sets_drawn": len(deltas),
        "set_size": size,
        "deltas": [round(d, 4) for d in ordered],
        "mean_delta": round(sum(deltas) / len(deltas), 4) if deltas else None,
        "median_delta": None if not deltas else round(median(deltas) or 0.0, 4),
        "p95_delta": ordered[int(0.95 * (len(ordered) - 1))] if ordered else None,
        "max_delta": ordered[-1] if ordered else None,
        "flat": bool(deltas and abs(sum(deltas) / len(deltas)) < 0.05),
        "what_it_controls": "the comparison itself: how much advantage a gene set of the same size "
        "gets from the structure of the gate space alone",
    }


def benchmark_test(
    rows: list[dict[str, Any]],
    benchmark: dict[str, str] = BENCHMARK,
    sample: int = RANDOM_COMBINATIONS,
    seed: int = SEED,
    pool: list[str] | None = None,
    placebo: list[str] | None = None,
    null_draws: int = NULL_GENE_SETS,
    levels: dict[str, dict[str, float]] | None = None,
) -> dict[str, Any]:
    """The pre-registered claim: do approved targets outrank random combinations?"""
    rng = random.Random(seed)
    p, delta, marked, drawn = _advantage(rows, set(benchmark), rng, sample)
    n = len(rows)
    percentile = [100.0 * (n - i) / n for i in range(n)]
    median_marked = median([percentile[i] for i in marked])
    median_random = median([percentile[i] for i in drawn])
    if pool is None:
        pool = sorted({g for r in rows for g in r["required"]})
    candidates = placebo if placebo else [g for g in pool if g not in benchmark]
    null = placebo_null(
        rows,
        candidates,
        size=len(set(benchmark) & set(pool)) or len(benchmark),
        draws=null_draws,
        sample=sample,
        seed=seed,
    )
    above_null = bool(null["p95_delta"] is not None and delta > null["p95_delta"])
    supported = bool(p <= 0.05 and delta >= MIN_BENCHMARK_DELTA and above_null)
    return {
        "claim": PREREGISTERED_CLAIM,
        "best_gate_per_gene": best_gate_comparison(rows, sorted(set(benchmark) & set(pool)), candidates),
        "delta_against_the_null": None
        if null["mean_delta"] is None
        else round(delta - null["mean_delta"], 3),
        "delta_percentile_in_null": _percentile_in(delta, null["deltas"]),
        "combinations_with_a_benchmark_target": len(marked),
        "random_comparison_combinations": len(drawn),
        "seed": seed,
        "median_percentile_benchmark": None if median_marked is None else round(median_marked, 1),
        "median_percentile_random": None if median_random is None else round(median_random, 1),
        "median_score_benchmark": median([rows[i]["score"] for i in marked]),
        "median_score_random": median([rows[i]["score"] for i in drawn]),
        "p_one_sided": p,
        "cliffs_delta": round(delta, 3),
        "minimum_delta_declared": MIN_BENCHMARK_DELTA,
        "test": "Mann-Whitney rank sum, one-sided, on the negated score so that the alternative is "
        "benchmark combinations scoring higher",
        "placebo_null": null,
        "placebo_markers": len(candidates),
        "delta_above_null_p95": above_null,
        "supported": supported,
        "best_per_benchmark_gene": {
            gene: _best_gate(rows, gene, percentile, levels)
            for gene in sorted(benchmark)
            if gene in set(pool)
        },
        "benchmark_genes_not_in_the_pool": sorted(set(benchmark) - set(pool)),
    }


def _percentile_in(value: float, distribution: list[float]) -> float | None:
    """Where a value falls in a sorted distribution, as a percentile."""
    if not distribution:
        return None
    below = sum(1 for d in distribution if d < value)
    return round(100.0 * below / len(distribution), 1)


def best_gate_comparison(
    rows: list[dict[str, Any]],
    benchmark_genes: list[str],
    placebo_genes: list[str],
) -> dict[str, Any]:
    """Post-hoc: each gene's *best* gate, benchmark set against the matched placebo set.

    Not the pre-registered statistic. The pre-registration compared every gate
    containing a benchmark gene against random gates, and that statistic averages
    over the hundreds of pairings a target programme would never look at - a
    target is chosen with its best partner, not with its mean partner. This is
    the same comparison at one gate per gene, with the placebo markers as the
    control, and it is reported as post-hoc because it was written after the
    pre-registered verdict was read.
    """
    n = len(rows)
    best: dict[str, float] = {}
    for i, r in enumerate(rows):
        for gene in r["required"]:
            if gene not in best:
                best[gene] = 100.0 * (n - i) / n
    marked = [g for g in benchmark_genes if g in best]
    control = [g for g in placebo_genes if g in best]
    if not marked or not control:
        return {"post_hoc": True, "benchmark_genes": len(marked), "placebo_genes": len(control)}
    values = {f"m{g}": -best[g] for g in marked} | {f"p{g}": -best[g] for g in control}
    rank_of, ties = ranks(values)
    p, delta, _u = rank_sum_p(rank_of, ties, [f"m{g}" for g in marked])
    return {
        "post_hoc": True,
        "what": "one gate per gene, the best one, benchmark markers against the placebo markers that "
        "entered the pool the same way",
        "benchmark_genes": len(marked),
        "placebo_genes": len(control),
        "median_best_percentile_benchmark": round(median([best[g] for g in marked]) or 0.0, 1),
        "median_best_percentile_placebo": round(median([best[g] for g in control]) or 0.0, 1),
        "p_one_sided": p,
        "cliffs_delta": round(-delta, 3),
    }


def _best_gate(
    rows: list[dict[str, Any]],
    gene: str,
    percentile: list[float],
    levels: dict[str, dict[str, float]] | None,
) -> dict[str, Any] | None:
    """The gate a benchmark gene reaches, with why it sits where it does."""
    for i, r in enumerate(rows):
        if gene in r["required"]:
            out = {
                k: r.get(k)
                for k in (
                    "required",
                    "excluded",
                    "logic",
                    "lineage",
                    "coverage",
                    "healthy_burden_ntpm",
                    "worst_essential_tissue",
                    "score",
                )
            }
            out["rank"] = i + 1
            out["percentile"] = round(percentile[i], 1)
            if levels and gene in levels:
                worst = max(levels[gene].items(), key=lambda kv: kv[1])
                out["own_worst_essential_tissue"] = {"tissue": worst[0], "nTPM": worst[1]}
            return out
    return None


def build_pool(
    singles: list[dict[str, Any]],
    gates: Gates,
    levels: dict[str, dict[str, float]],
    pool_size: int = POOL_SIZE,
    placebo_size: int = PLACEBO_POOL_SIZE,
    seed: int = SEED,
) -> tuple[list[str], list[str]]:
    """The markers gates are built from: the best singles, the benchmark, and a placebo set.

    The benchmark genes are forced in whether or not they score their way in,
    which is the honest thing to do and also an asymmetry: gates built on a
    marker that earned its place in the pool are not comparable to gates built
    on one that was carried in. So an equally arbitrary placebo set is carried
    in the same way, and it is the placebo set the null draws from.
    """
    ranked = [r["required"][0] for r in singles]
    top = ranked[:pool_size]
    forced = sorted((set(BENCHMARK) & set(gates.present) & set(levels)) - set(top))
    rest = [g for g in ranked[pool_size:] if g not in BENCHMARK]
    rng = random.Random(seed + 2)
    placebo = sorted(rng.sample(rest, min(placebo_size, len(rest))))
    return sorted(set(top) | set(forced) | set(placebo)), placebo


def threshold_sensitivity(
    expression: dict[str, dict[str, float]],
    lineage_of: dict[str, str],
    levels: dict[str, dict[str, float]],
    thresholds: tuple[float, ...] = SENSITIVITY_THRESHOLDS,
    pool_size: int = SENSITIVITY_POOL,
    log=None,
) -> dict[str, Any]:
    """Post-hoc: the same claim at stricter presence thresholds.

    Not part of the pre-registration. It exists because a marker called present
    at about 7 TPM is a different object from a marker at 60 TPM, and an
    approved target is always the second kind, so the threshold is the first
    thing a failed claim should be read against.
    """
    out: list[dict[str, Any]] = []
    for threshold in thresholds:
        gates = Gates(expression, lineage_of, threshold=threshold)
        if not gates.lineages:
            continue
        scored = sorted(set(gates.present) & set(levels))
        singles = rank_combinations(gates, levels, single_gates(scored))
        pool, placebo = build_pool(singles, gates, levels, pool_size)
        rows = rank_combinations(gates, levels, combinations(pool))
        bench = benchmark_test(rows, pool=pool, placebo=placebo, null_draws=50)
        if log:
            print(
                f"  sensitivity at log2(TPM+1) >= {threshold}: delta {bench['cliffs_delta']}, "
                f"p {bench['p_one_sided']:.3g}, null mean {bench['placebo_null']['mean_delta']}, "
                f"{'supported' if bench['supported'] else 'not supported'}",
                file=log,
                flush=True,
            )
        out.append(
            {
                "present_log2_tpm_plus_1": threshold,
                "markers_expressed_anywhere": len(scored),
                "combinations_scored": len(rows),
                "cliffs_delta": bench["cliffs_delta"],
                "p_one_sided": bench["p_one_sided"],
                "median_percentile_benchmark": bench["median_percentile_benchmark"],
                "placebo_null_mean_delta": bench["placebo_null"]["mean_delta"],
                "placebo_null_p95_delta": bench["placebo_null"]["p95_delta"],
                "supported": bench["supported"],
                "top_gate": None
                if not rows
                else {
                    k: rows[0][k]
                    for k in ("required", "excluded", "lineage", "coverage", "healthy_burden_ntpm", "score")
                },
            }
        )
    return {
        "post_hoc": True,
        "why": "the pre-registered verdict is the one at the declared threshold; this says how far "
        "that verdict depends on where 'present' was drawn",
        "rows": out,
    }


def run(
    net: bool = True,
    log=None,
    cache_dir: Path = depmap.CACHE,
    pool_size: int = POOL_SIZE,
) -> dict[str, Any]:
    """The whole screen: universe, gates, ranking, and the pre-registered benchmark."""
    universe = surface_universe()
    levels = healthy_levels(set(universe))
    if log:
        print(
            f"  {len(universe)} surface genes in the packaged atlas, {len(levels)} with essential-tissue "
            f"values",
            file=log,
            flush=True,
        )
    table = depmap.matrix("expression", set(universe), net=net, cache_dir=cache_dir, decimals=2, log=log)
    model_rows = depmap.models(net=net, cache_dir=cache_dir, log=log)
    lineage_of = {
        m: (row.get("OncotreeLineage") or "")
        for m, row in model_rows.items()
        if row.get("OncotreeLineage") and row.get("OncotreePrimaryDisease") != EXCLUDE_DISEASE
    }
    excluded_lines = sum(
        1 for row in model_rows.values() if row.get("OncotreePrimaryDisease") == EXCLUDE_DISEASE
    )
    gates = Gates(table["values"], lineage_of)
    if log:
        print(
            f"  {len(gates.lines)} cell lines in {len(gates.lineages)} lineages of at least "
            f"{MIN_LINES} lines",
            file=log,
            flush=True,
        )
    scored = sorted(set(table["genes"]) & set(levels))
    singles = rank_combinations(gates, levels, single_gates(scored))
    pool, placebo = build_pool(singles, gates, levels, pool_size)
    if log:
        print(
            f"  pool: {len(pool)} markers ({len(singles)} single markers scored, {len(placebo)} placebo)",
            file=log,
            flush=True,
        )
    rows = rank_combinations(gates, levels, combinations(pool))
    benchmark = benchmark_test(rows, pool=pool, placebo=placebo, levels=levels)
    sensitivity = threshold_sensitivity(table["values"], lineage_of, levels, log=log)
    ranked = [
        {**r, "evidence": evidence_for(r).to_dict(), "confidence": evidence_for(r).confidence}
        for r in rows[:40]
    ]
    single_by_gene = {r["required"][0]: r for r in singles}
    gains = []
    for row in rows[:200]:
        base = max(
            (single_by_gene[g]["score"] for g in row["required"] if g in single_by_gene),
            default=None,
        )
        if base is not None and row["logic"] != "A" and row["score"] > base:
            gains.append(
                {
                    **{k: row[k] for k in ("required", "excluded", "logic", "lineage", "coverage")},
                    "score": row["score"],
                    "best_single_marker_score": base,
                    "gain": round(row["score"] - base, 4),
                }
            )
    return {
        "claim": PREREGISTERED_CLAIM,
        "thresholds": {
            "present_log2_tpm_plus_1": PRESENT_LOG2TPM,
            "not_marker_protects_at_ntpm": PROTECT_NTPM,
            "co_expression_floor": CO_EXPRESSION_FLOOR,
            "min_lines_per_lineage": MIN_LINES,
            "excluded_lines": f"{excluded_lines} models DepMap calls {EXCLUDE_DISEASE}",
            "pool_size": pool_size,
            "score": "coverage / (1 + healthy_burden_nTPM / 10), where the burden is the worst "
            "essential healthy tissue the gate still reaches, floored at a quarter of its "
            "highest required marker",
            "essential_tissues": list(essential_columns()),
        },
        "benchmark_set": BENCHMARK,
        "data": {
            "provenance": depmap.provenance(),
            "healthy_atlas": {
                "source": "Human Protein Atlas consensus tissue RNA, CC BY-SA 4.0, packaged in "
                "genomeos/therapeutics/data/normal_tissue.json.gz",
                "genes": len(universe),
                "genes_with_essential_tissue_values": len(levels),
            },
            "cell_lines": len(gates.lines),
            "lineages": {k: gates.lineage_size[k] for k in sorted(gates.lineages)},
            "genes_in_both": len(set(table["genes"]) & set(levels)),
            "expression_genes_missing_from_depmap": len(table["missing"]),
        },
        "single_markers_scored": len(singles),
        "combinations_scored": len(rows),
        "top_single_markers": [
            {k: r[k] for k in ("required", "lineage", "coverage", "healthy_burden_ntpm", "score")}
            for r in singles[:15]
        ],
        "top_combinations": ranked,
        "combinations_that_beat_their_best_single_marker": len(gains),
        "gains_by_logic": {
            kind: sum(1 for g in gains if g["logic"] == kind) for kind in sorted({g["logic"] for g in gains})
        },
        "biggest_gains_over_a_single_marker": gains[:40],
        "benchmark": benchmark,
        "threshold_sensitivity": sensitivity,
        "limits": list(LIMITS),
        "not_a_design": (
            "The output names marker logic only. No binder, format, construct, cargo or protocol "
            "follows from it here."
        ),
    }


__all__ = [
    "BENCHMARK",
    "PREREGISTERED_CLAIM",
    "Gates",
    "benchmark_test",
    "burden",
    "combinations",
    "essential_columns",
    "evidence_for",
    "healthy_levels",
    "rank_combinations",
    "run",
    "score",
    "surface_universe",
]
