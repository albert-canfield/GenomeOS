# SPDX-License-Identifier: AGPL-3.0-or-later
"""One table over every chromosome's all-elements summary, as the chain lands them.

    uv run python scripts/enhancer_targets_all_genome_wide.py [--no-control]

Reads data/results/enhancer_targets_all_chr*.json (summaries only; the element tables stay local) and
writes enhancer_targets_all_genome_wide.json: per chromosome the counts and rates, the totals over the
chromosomes complete so far with the requests spent, the spread of each rate across chromosomes rather
than one pooled number, the controls each rate has to be read against, and whether the rates drift with
chromosome size (chr21, chr22 and chrY are acrocentric or nearly gene-free, and a reading taken on them
alone has already misled one lane). Every earlier reading is kept in `history`, one entry per distinct
set of complete chromosomes, so the pooled rates can be watched moving as the sweep widens.

The node-containment rate carries its own control, measured on this element set rather than borrowed
from the archive: the same elements and the same most-moved coding gene, scored against as many
boundaries placed uniformly at random (20 draws, seed 7), using the same `infer_domains` caller the
scorer used. It needs the local element tables and takes about half a minute; `--no-control` skips it.
"""

from __future__ import annotations

import bisect
import json
import random
import sys
from pathlib import Path
from statistics import median

from genomeos.molecules.proteome import CHROM_LENGTHS
from genomeos.results import save_result

ARCHIVE = Path("data/knowledge/alphagenome/all_elements")
SHUFFLES = 20

# Chromosomes whose architecture is not the genome's: two acrocentric arms and the male-specific
# chromosome, all three small and gene-poor. Rates are reported with and without them.
UNUSUAL = ("chr21", "chr22", "chrY")

# What each rate has to be read against. Both come from measurements this project already made, and
# both qualify a rate that reads as a result on its own.
CONTROLS = {
    "fraction_with_target": {
        "control": 0.87,
        "what": "a derived layer names some target at 52 of 60 matched random windows (87%), "
        "against 11 of 12 published loci",
        "source": "docs/LOCI-BENCHMARK.md, data/results/loci_benchmark.json (commits 303d0b0, ef80075)",
        "reading": "a naming rate is not evidence by itself; what separates a real element from a "
        "matched window is a direction and a variable position on constrained sequence",
    },
    "coding_target_inside_domain": {
        "control": 0.791,
        "measured": 0.817,
        "excess": 0.026,
        "what": "over the whole deletion archive of 113,399 elements the same caller keeps a coding "
        "target inside the element's own CTCF node for 81.7% of elements, against 79.1% for the same "
        "number of boundaries placed at random: 2.6 points, not nine times out of ten",
        "source": "data/results/domains_oriented_comparison.json (commit dbad593)",
        "reading": "the node's advantage over a random partition at the same resolution is small, and "
        "it is measured on a model's reading rather than on observed enhancer-gene pairs",
    },
    "coding_target_agrees_with_nearest": {
        "control": None,
        "what": "no matched-window control has been measured for nearest-TSS agreement; the nearest "
        "TSS heuristic names the published target at 8 of 12 loci, the same count the derived layers "
        "reach on a different 8",
        "source": "docs/LOCI-BENCHMARK.md (commit 303d0b0)",
        "reading": "quote it as agreement between two readings of the same element, not as accuracy",
    },
}


def _spearman(xs: list[float], ys: list[float]) -> float | None:
    """Rank correlation, no ties expected among chromosome lengths or rates."""
    n = len(xs)
    if n < 4:
        return None

    def ranks(vs: list[float]) -> list[float]:
        order = sorted(range(n), key=lambda i: vs[i])
        r = [0.0] * n
        i = 0
        while i < n:
            j = i
            while j + 1 < n and vs[order[j + 1]] == vs[order[i]]:
                j += 1
            avg = (i + j) / 2 + 1
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r

    rx, ry = ranks(xs), ranks(ys)
    mx, my = sum(rx) / n, sum(ry) / n
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry, strict=True))
    den = (sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry)) ** 0.5
    return round(num / den, 3) if den else None


def _rates(counts: dict) -> dict:
    """The three headline rates from one chromosome's or one pool's counts."""
    agree = counts["verdicts"].get("agrees with nearest TSS in domain", 0)
    inside = agree + counts["verdicts"].get("another gene in the same domain", 0)
    coding = counts["coding_named"]
    return {
        "fraction_with_target": round(counts["named"] / counts["scored"], 3) if counts["scored"] else None,
        "coding_target_agrees_with_nearest": round(agree / coding, 3) if coding else None,
        "coding_target_inside_domain": round(inside / coding, 3) if coding else None,
    }


def _pool(rows: list[dict]) -> dict:
    tot = {"elements": 0, "scored": 0, "named": 0, "coding_named": 0, "strong": 0, "silencer_like": 0}
    verdicts: dict[str, int] = {}
    for r in rows:
        for k in ("elements", "scored", "named", "coding_named", "strong", "silencer_like"):
            tot[k] += r["counts"][k]
        for k, n in r["counts"]["verdicts"].items():
            verdicts[k] = verdicts.get(k, 0) + n
    tot["verdicts"] = verdicts
    return tot


def _spread(rows: list[dict], key: str) -> dict | None:
    vals = [(r["rates"][key], r["chrom"]) for r in rows if r["rates"][key] is not None]
    if not vals:
        return None
    vals.sort()
    return {
        "n": len(vals),
        "min": vals[0][0],
        "lowest": vals[0][1],
        "median": round(median(v for v, _ in vals), 3),
        "max": vals[-1][0],
        "highest": vals[-1][1],
        "range": round(vals[-1][0] - vals[0][0], 3),
        "per_chromosome": {c: v for v, c in vals},
    }


def aggregate(results_dir: Path = Path("data/results")) -> dict:
    rows: dict[str, dict] = {}
    complete: list[dict] = []
    requests = 0
    for f in sorted(results_dir.glob("enhancer_targets_all_chr*.json")):
        d = json.loads(f.read_text())
        s = d.get("summary") or {}
        vc = s.get("verdicts_coding") or {}
        coding_named = sum(n for k, n in vc.items() if k != "no predicted effect")
        counts = {
            "elements": d.get("elements_total") or 0,
            "scored": d.get("scored") or 0,
            "named": s.get("with_predicted_target") or 0,
            "coding_named": coding_named,
            "strong": s.get("strong") or 0,
            "silencer_like": s.get("silencer_like") or 0,
            "verdicts": dict(vc),
        }
        rates = _rates(counts)
        rows[d["chrom"]] = {
            "elements": d.get("elements_total"),
            "scored": d.get("scored", 0),
            "complete": bool(d.get("complete")),
            "length": CHROM_LENGTHS.get(d["chrom"]),
            "fraction_with_target": s.get("fraction_with_target"),
            "coding_target_agrees_with_nearest": s.get("coding_target_agrees_with_nearest"),
            "coding_target_inside_domain": rates["coding_target_inside_domain"],
            "fraction_inside_domain": s.get("fraction_inside_domain"),
            "strong": s.get("strong"),
            "silencer_like": s.get("silencer_like"),
            "requests_this_run": d.get("requests_this_run"),
        }
        requests += d.get("requests_this_run") or 0
        if d.get("complete"):
            complete.append({"chrom": d["chrom"], "counts": counts, "rates": rates})

    keys = ("fraction_with_target", "coding_target_agrees_with_nearest", "coding_target_inside_domain")
    pooled = _pool(complete)
    genome = {
        **{k: v for k, v in pooled.items() if k != "verdicts"},
        **_rates(pooled),
        "verdicts_coding": dict(sorted(pooled["verdicts"].items(), key=lambda kv: -kv[1])),
    }
    usual = [r for r in complete if r["chrom"] not in UNUSUAL]
    unusual = [r for r in complete if r["chrom"] in UNUSUAL]

    sized = [r for r in complete if r["chrom"] in CHROM_LENGTHS]
    lengths = [CHROM_LENGTHS[r["chrom"]] for r in sized]
    density = [r["counts"]["elements"] / CHROM_LENGTHS[r["chrom"]] * 1e6 for r in sized]
    trend = {
        "axis": "hg38 chromosome length against the chromosome's own rate, over the complete chromosomes",
        "chromosomes": len(lengths),
        **{k: _spearman(lengths, [r["rates"][k] for r in sized]) for k in keys},
        "against_element_density": {
            "axis": "registry elements per Mb, the stand-in for gene density this run already holds",
            **{k: _spearman(density, [r["rates"][k] for r in sized]) for k in keys},
            "per_chromosome": dict(
                sorted(
                    ((r["chrom"], round(d, 1)) for r, d in zip(sized, density, strict=True)),
                    key=lambda p: -p[1],
                )
            ),
        },
        "excluding_acrocentric_and_chrY": {
            "chromosomes": len(usual),
            "pooled": _rates(_pool(usual)) if usual else None,
        },
        "acrocentric_and_chrY_only": {
            "chromosomes": len(unusual),
            "chroms": [r["chrom"] for r in unusual],
            "pooled": _rates(_pool(unusual)) if unusual else None,
        },
        "reading": "a rate that moves with chromosome size is a rate the small chromosomes bought; "
        "chr21 and chr22 are acrocentric and chrY is nearly gene-free, and the compression lane's "
        "per-chromosome layers already looked positive there and went negative genome-wide",
    }

    return {
        "chromosomes": rows,
        "complete_chromosomes": len(complete),
        "chromosomes_complete": sorted(r["chrom"] for r in complete),
        "genome": genome,
        "spread": {k: _spread(complete, k) for k in keys},
        "controls": CONTROLS,
        "size_trend": trend,
        "requests_total": requests,
        "cells": ["K562", "HepG2", "GM12878", "IMR-90"],
        "elements_where": "data/knowledge/alphagenome/all_elements/<chrom>.json (local)",
        "evidence": "predicted: AlphaGenome deletion effect per element, per gene and per cell line; "
        "inferred: CTCF-only nodes",
    }


def _inside(starts: list[int], pairs: list[tuple[int, int]]) -> int:
    """How many (element midpoint, target TSS) pairs fall in one node of a partition at `starts`."""
    return sum(1 for m, t in pairs if bisect.bisect_right(starts, m) == bisect.bisect_right(starts, t))


BIN = 10_000_000  # the scale at which boundary placement and element placement are compared

# What the sign of the node's excess over random might be a property of. Each is measurable from what
# the control already loads, and each is a different account of why real boundaries beat scattered ones
# on some chromosomes and lose on others.
SHAPE = (
    "boundaries_per_mb",
    "median_node_kb",
    "node_length_cv",
    "median_pair_span_kb",
    "coding_genes_per_mb",
    "elements_per_mb",
    "boundaries_follow_the_elements",
    "inside_random",
)


def _shape(chrom: str, length: int, starts: list[int], pairs: list[tuple[int, int]], genes: int) -> dict:
    """The chromosome's partition and its element set, described so the sign of the excess can be
    regressed on something. `boundaries_follow_the_elements` is the mechanism worth testing: real CTCF
    boundaries cluster where genes are, and so do the elements and their targets, so a chromosome whose
    boundaries crowd into the element-dense stretches spends its cuts where the pairs are and can lose
    to boundaries scattered uniformly, most of which land in empty sequence and cut nothing."""
    mb = length / 1e6
    nodes = [b - a for a, b in zip(starts, [*starts[1:], length], strict=True)]
    mean = sum(nodes) / len(nodes)
    var = sum((x - mean) ** 2 for x in nodes) / len(nodes)
    bins = length // BIN + 1
    per_bin_boundaries = [0] * bins
    per_bin_pairs = [0] * bins
    for e in starts[1:]:
        per_bin_boundaries[e // BIN] += 1
    for m, _ in pairs:
        per_bin_pairs[m // BIN] += 1
    return {
        "boundaries_per_mb": round(len(starts[1:]) / mb, 2),
        "median_node_kb": round(median(nodes) / 1e3, 1),
        "node_length_cv": round(var**0.5 / mean, 3) if mean else None,
        "median_pair_span_kb": round(median(abs(m - t) for m, t in pairs) / 1e3, 1),
        "coding_genes_per_mb": round(genes / mb, 2),
        "elements_per_mb": round(len(pairs) / mb, 1),
        "boundaries_follow_the_elements": _spearman(per_bin_boundaries, per_bin_pairs),
    }


def node_control(chroms: list[str]) -> dict:
    """The node-containment rate against as many boundaries placed at random, on this element set.

    dbad593 measured 81.7% against 79.1% over the archive as it stood, a 2.6-point excess. That control
    was placed against the archive's own elements; this one is placed against the elements the fold is
    reporting, chromosome by chromosome, so the excess can be read where it is earned and where it is
    not. The containment figure it recomputes is the scorer's own verdict rebuilt from the element
    tables, which is also a check that the two agree.
    """
    from genomeos.genome import Annotation, default_gencode
    from genomeos.genome.domains import infer_domains
    from genomeos.genome.regulatory import load_ccres

    rows: dict[str, dict] = {}
    skipped: dict[str, str] = {}
    for chrom in chroms:
        arch = ARCHIVE / f"{chrom}.json"
        length = CHROM_LENGTHS.get(chrom)
        if not arch.exists() or not length:
            skipped[chrom] = "no local element table" if length else "no length"
            continue
        gff = default_gencode({chrom})
        ccres = load_ccres(chrom)
        if not gff or not ccres:
            skipped[chrom] = "no annotation or cCREs"
            continue
        ann = Annotation.from_gff3(gff, {chrom})
        coding = {
            g.symbol: g for g in ann.genes.values() if g.locus.chrom == chrom and g.type == "protein_coding"
        }
        pairs = []
        for e in json.loads(arch.read_text()):
            gene = (e.get("predicted_coding") or {}).get("gene")
            g = coding.get(gene)
            if g is not None:
                tss = g.locus.end - 1 if g.locus.strand.value == "-" else g.locus.start
                pairs.append(((e["start"] + e["end"]) // 2, tss))
        if not pairs:
            skipped[chrom] = "no element names a coding gene"
            continue
        starts = [d.start for d in infer_domains(chrom, length, ccres, ann)]
        edges = starts[1:]
        rng = random.Random(7)
        null = [
            _inside([0, *sorted(rng.randint(1, length - 1) for _ in edges)], pairs) for _ in range(SHUFFLES)
        ]
        ins, tot = _inside(starts, pairs), len(pairs)
        rnd = sum(null) / len(null) / tot
        rows[chrom] = {
            "coding_named": tot,
            "boundaries": len(edges),
            "inside": round(ins / tot, 4),
            "inside_random": round(rnd, 4),
            "excess": round(ins / tot - rnd, 4),
            **_shape(chrom, length, starts, pairs, len(coding)),
        }
    tot = sum(r["coding_named"] for r in rows.values())
    ins = sum(r["inside"] * r["coding_named"] for r in rows.values())
    rnd = sum(r["inside_random"] * r["coding_named"] for r in rows.values())
    won = sum(1 for r in rows.values() if r["excess"] > 0)
    return {
        "what": "the same elements and the same most-moved coding gene, against as many boundaries "
        f"placed uniformly at random ({SHUFFLES} draws, seed 7), on the CTCF-only caller the scorer used",
        "chromosomes": len(rows),
        "coding_named": tot,
        "inside": round(ins / tot, 4) if tot else None,
        "inside_random": round(rnd / tot, 4) if tot else None,
        "excess": round((ins - rnd) / tot, 4) if tot else None,
        "chromosomes_where_the_node_beats_random": won,
        "what_predicts_the_sign": _predicts_the_sign(rows),
        "per_chromosome": dict(sorted(rows.items(), key=lambda kv: -kv[1]["excess"])),
        "skipped": skipped,
        "reading": "an excess that is positive on some chromosomes and negative on others is not a "
        "property of CTCF nodes; it is a property of where the boundaries happen to fall",
    }


# Fixed on 2026-09-15 from the twelve chromosomes complete at the time (chr12 to chr22 and chrY), where
# boundary density is the one shape measure that tracks the sign of the node's excess over random
# (Spearman +0.82, +0.76 with chrY dropped) and separates it at 11 of 12. Every chromosome that lands
# afterwards is a held-out test of it, and the fold scores each one as it arrives. The rule is fitted,
# not derived, and it is recorded here so that it cannot be quietly refitted when it fails.
SIGN_RULE = {
    "fixed_on": [
        "chr12",
        "chr13",
        "chr14",
        "chr15",
        "chr16",
        "chr17",
        "chr18",
        "chr19",
        "chr20",
        "chr21",
        "chr22",
        "chrY",
    ],  # noqa: E501
    "fixed_at": "2026-09-15",
    "rule": "the node beats as many random boundaries where the caller cuts the chromosome more finely "
    "than 5.8 boundaries per Mb, and loses below it; within 0.3 per Mb of the line the excess is near "
    "zero and the call is refused",
    "threshold_per_mb": 5.8,
    "too_close_per_mb": 0.3,
    "fitted_errors": ["chr19, 5.77 per Mb, predicted to lose and it wins by 0.7 points"],
    "why": "an excess that is a function of how finely the caller cuts is a statement about resolution "
    "rather than about CTCF: where CTCF-only elements are sparse the nodes carry nothing a partition of "
    "the same count placed at random does not",
}


# One word per chromosome, so a reader who has not followed the night can see which chromosomes the
# rule was fitted on and which were tests of it without reconstructing either set.
STATUSES = (
    "held out: right",
    "held out: WRONG",
    "held out: refused",
    "fitted on, not a test",
    "still to land",
)


def sign_call(per_mb: float | None) -> str | None:
    if per_mb is None:
        return None
    if abs(per_mb - SIGN_RULE["threshold_per_mb"]) < SIGN_RULE["too_close_per_mb"]:
        return "too close to call"
    return "positive" if per_mb > SIGN_RULE["threshold_per_mb"] else "negative"


def boundary_density(chroms: list[str], known: dict | None = None) -> dict[str, float]:
    """Boundaries per Mb for each chromosome, which needs no element table and so can be computed for
    chromosomes the sweep has not reached. Cached from the previous fold: it does not change."""
    from genomeos.genome import Annotation, default_gencode
    from genomeos.genome.domains import infer_domains
    from genomeos.genome.regulatory import load_ccres

    out = dict(known or {})
    for chrom in chroms:
        length = CHROM_LENGTHS.get(chrom)
        if chrom in out or not length:
            continue
        gff, ccres = default_gencode({chrom}), load_ccres(chrom)
        if not gff or not ccres:
            continue
        doms = infer_domains(chrom, length, ccres, Annotation.from_gff3(gff, {chrom}))
        out[chrom] = round(len(doms[1:]) / (length / 1e6), 2)
    return out


def sign_prediction(control: dict, density: dict[str, float]) -> dict:
    """The rule's standing predictions, and its score on every chromosome that has landed since it was
    fixed. A chromosome in `fixed_on` is in sample and is not counted towards the held-out score."""
    measured = control.get("per_chromosome") or {}
    rows: dict[str, dict] = {}
    by_status: dict[str, list[str]] = {k: [] for k in STATUSES}
    for chrom, per_mb in sorted(density.items(), key=lambda kv: -kv[1]):
        call = sign_call(per_mb)
        row = {"boundaries_per_mb": per_mb, "predicted": call}
        if chrom not in measured:
            status = "still to land"
        elif chrom in SIGN_RULE["fixed_on"]:
            status = "fitted on, not a test"
        elif call == "too close to call":
            status = "held out: refused"
        else:
            right = (call == "positive") == (measured[chrom]["excess"] > 0)
            status = "held out: right" if right else "held out: WRONG"
        if chrom in measured:
            row["excess"] = measured[chrom]["excess"]
            row["actual"] = "positive" if measured[chrom]["excess"] > 0 else "negative"
        row["status"] = status
        rows[chrom] = row
        by_status[status].append(chrom)
    right, wrong = len(by_status["held out: right"]), len(by_status["held out: WRONG"])
    refused, landed = len(by_status["held out: refused"]), len(by_status["still to land"])
    return {
        **SIGN_RULE,
        "scoreboard": f"held out: {right} right, {wrong} wrong, {refused} refused; "
        f"{len(by_status['fitted on, not a test'])} fitted on and not tests; {landed} still to land",
        "held_out_chromosomes_scored": right + wrong,
        "held_out_right": right,
        "held_out_wrong": wrong,
        "held_out_refused": refused,
        "verdict": "FAILED on " + ", ".join(by_status["held out: WRONG"])
        if wrong
        else ("holding" if right else "not yet tested"),
        "chromosomes_by_status": by_status,
        "standing_predictions": {c: rows[c]["predicted"] for c in by_status["still to land"]},
        "per_chromosome": rows,
    }


def _predicts_the_sign(rows: dict[str, dict]) -> dict:
    """Rank the candidate accounts of why the excess changes sign. chrY carries 125 elements against
    tens of thousands elsewhere and is the extreme of every shape measure, so every correlation is
    given with it and without it: an account that only works because of chrY is not an account."""

    def run(rs: dict[str, dict]) -> dict:
        ex = [r["excess"] for r in rs.values()]
        out = {}
        for k in SHAPE:
            xs = [r[k] for r in rs.values()]
            out[k] = None if any(x is None for x in xs) else _spearman(xs, ex)
        return {"chromosomes": len(rs), **out}

    return {
        "axis": "Spearman of each chromosome's shape against its excess over random boundaries",
        "all": run(rows),
        "without_chrY": run({k: v for k, v in rows.items() if k != "chrY"}),
    }


def with_history(out: dict, previous: dict | None) -> dict:
    """Keep every earlier reading beside the current one, so the pooled rates can be watched moving as
    chromosomes land instead of being recomputed and forgotten. One entry per distinct chromosome set."""
    history = list((previous or {}).get("history") or [])
    entry = {
        "chromosomes": out["complete_chromosomes"],
        "chromosomes_complete": out["chromosomes_complete"],
        "scored": out["genome"]["scored"],
        **{k: out["genome"][k] for k in out["spread"]},
    }
    if not history or history[-1]["chromosomes_complete"] != entry["chromosomes_complete"]:
        history.append(entry)
    else:
        history[-1] = entry
    return {**out, "history": history}


def main() -> int:
    out = aggregate()
    path = Path("data/results/enhancer_targets_all_genome_wide.json")
    previous = json.loads(path.read_text()) if path.exists() else None
    if "--no-control" not in sys.argv:
        control = node_control(out["chromosomes_complete"])
        cached = ((previous or {}).get("node_excess_prediction") or {}).get("per_chromosome") or {}
        density = boundary_density(
            [c for c in CHROM_LENGTHS if c != "chrM"],
            {c: r["boundaries_per_mb"] for c, r in cached.items() if r.get("boundaries_per_mb")},
        )
        out["controls"] = {
            **out["controls"],
            "coding_target_inside_domain": {
                **out["controls"]["coding_target_inside_domain"],
                "measured_on_this_set": control,
            },
        }
        out["node_excess_prediction"] = sign_prediction(control, density)
    out = with_history(out, previous)
    save_result("enhancer_targets_all_genome_wide", out)
    g, sp = out["genome"], out["spread"]
    print(
        f"{out['complete_chromosomes']} chromosomes complete: {g['scored']:,} elements, "
        f"{g['fraction_with_target']} name a gene (87% at matched random windows), "
        f"nearest TSS {g['coding_target_agrees_with_nearest']}, "
        f"inside the node {g['coding_target_inside_domain']} (random boundaries 0.791, +2.6 points); "
        f"{out['requests_total']:,} requests so far"
    )
    nc = out["controls"]["coding_target_inside_domain"].get("measured_on_this_set")
    if nc and nc["chromosomes"]:
        print(
            f"  node against random boundaries on this set: {nc['inside']} against {nc['inside_random']}, "
            f"{nc['excess']:+} over {nc['coding_named']:,} elements, the node ahead on "
            f"{nc['chromosomes_where_the_node_beats_random']} of {nc['chromosomes']} chromosomes"
        )
    for k, s in sp.items():
        if s:
            print(f"  {k}: {s['min']} ({s['lowest']}) to {s['max']} ({s['highest']}), median {s['median']}")
    t = out["size_trend"]
    dens = t["against_element_density"]
    if len(out["history"]) > 1:
        a, b = out["history"][-2], out["history"][-1]
        moved = ", ".join(f"{k} {a[k]} -> {b[k]}" for k in sp)
        print(f"  since {a['chromosomes']} chromosomes: {moved}")
    print(f"  over {t['chromosomes']} chromosomes, rho with length / with elements per Mb:")
    for k in sp:
        print(f"    {k}: {t[k]} / {dens[k]}")
    p = out.get("node_excess_prediction")
    if p:
        print(
            f"  the excess-sign rule ({p['threshold_per_mb']} boundaries per Mb) is {p['verdict']}"
            f" -- {p['scoreboard']}"
        )
        for chrom in (
            p["chromosomes_by_status"]["held out: WRONG"]
            + p["chromosomes_by_status"]["held out: right"]
            + p["chromosomes_by_status"]["held out: refused"]
        ):
            r = p["per_chromosome"][chrom]
            print(
                f"    {chrom}: {r['boundaries_per_mb']}/Mb predicted {r['predicted']}, "
                f"got {r['excess']:+} -- {r['status']}"
            )
        standing = ", ".join(f"{c} {v}" for c, v in list(p["standing_predictions"].items())[:6])
        print(f"    still to land: {standing}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
