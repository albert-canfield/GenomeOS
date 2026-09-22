# SPDX-License-Identifier: AGPL-3.0-or-later
"""Direction, held against measured perturbations: does deleting an element raise or lower its target?

`crispri.py` asks whether the deletion layer helps *rank* enhancer-gene pairs. It never asks the
prior question: when the layer names a target, does it get the **sign** right. Every enhancer-gene
call GenomeOS makes carries a direction — `action: "activates"` or `"represses"`, from the sign of a
predicted log2 fold change on deleting the element — and that direction has never been held against a
measurement at scale. It was registered undecidable at n=2 on the third locus set
(docs/LOCI-BENCHMARK.md section 21), and the fourth frame then showed what an unchecked layer costs.

The measurement exists and is already on disk. The ENCODE CRISPRi benchmark carries `EffectSize`, the
**signed** measured effect of silencing the element on the gene, beside the `Significant` flag. The
`Regulated` column that `crispri.py` uses is `Significant AND EffectSize < 0`, so the upward half of
the measured signal is discarded before that module sees it. Here `Significant` is read two-sided and
the sign is kept.

The model's side is the finished genome-wide deletion sweep (`data/knowledge/alphagenome/all_elements`),
which stores a signed `predicted_by_cell` / `predicted_coding_by_cell` log2 fold change per cell line.
Reading a sign out of it costs nothing: no AlphaGenome request is made anywhere in this module.

**The trap this module is built around.** A sign agreement looks like a coin flip, so 0.5 looks like
its chance level. It is not, and reporting it against 0.5 would manufacture a success. Two selections
run the other way:

- `Regulated` and the screens' own power make the measured side lopsided: downward effects dominate.
- The sweep keeps **one** target per element, and an upward (typically indirect) effect's gene is
  almost never that one target, so upward pairs are filtered out of the answerable set a second time.

The result is that the subset the question can be asked of may carry only one measured sign, and then
agreement is not a comparison at all — it is the model's own rate of emitting that sign, and a caller
that says "down" and nothing else scores 1.000. So the binding comparison here is not 0.5 but
`marginal_down_rate`: the share of **all** swept elements whose predicted change in the same cell line
is negative. `verdict` is decided against both, and `UNDECIDABLE` is a first-class outcome, so a high
agreement that a constant-sign caller would equal cannot be read as direction skill.
"""

from __future__ import annotations

import bisect
import csv
import gzip
import json
import math
from collections import defaultdict
from pathlib import Path
from typing import Any

from genomeos.attribution.crispri import HELDOUT, KNOWLEDGE, REACH, TRAINING

ELEMENTS = Path("data/knowledge/alphagenome/all_elements")
HALF_WINDOW = 524_288  # MODEL_WINDOW // 2, the scorer's reach each way (benchmark/loci.py read_reach)
CHANCE = 0.5
CHANCE_MARGIN = 0.15  # agreement must clear 0.5 by this much
MARGINAL_MARGIN = 0.10  # and clear the constant-sign caller by this much
EVIDENCE = (
    "experimental: signed EffectSize of CRISPRi enhancer-gene screens, ENCODE benchmark "
    "(EngreitzLab/CRISPR_comparison, Gschwind et al. 2025); predicted: signed AlphaGenome deletion "
    "log2 fold change per cell line, from the finished genome-wide sweep. No request"
)
PREREGISTERED = (
    "on the held-out arm only, in K562 and GM12878 (the only held-out cell types with an AlphaGenome "
    "line in the sweep), sign agreement between the measured EffectSize and the sweep's predicted "
    "log2 fold change in the screen's own cell line clears 0.5 by at least 0.15 AND clears the "
    "model's marginal down-rate over the whole sweep by at least 0.10, with the Wilson lower bound "
    "above that marginal rate. Clearing 0.5 but not the marginal rate is registered in advance as "
    "UNDECIDABLE, not as a pass"
)


def wilson(k: int, n: int, z: float = 1.96) -> dict[str, Any]:
    """Wilson 95% interval for a proportion; the normal interval is wrong at k == n."""
    if not n:
        return {"rate": None, "ci95": None, "k": k, "n": n}
    p = k / n
    d = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / d
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return {
        "rate": round(p, 4),
        "ci95": [round(max(0.0, centre - half), 4), round(min(1.0, centre + half), 4)],
        "k": k,
        "n": n,
    }


def rows(name: str, knowledge: Path = KNOWLEDGE) -> list[dict[str, str]]:
    with gzip.open(knowledge / name, "rt") as f:
        return [r for r in csv.DictReader(f, delimiter="\t") if r["ValidConnection"] == "TRUE"]


def signed(table: list[dict[str, str]]) -> list[dict[str, str]]:
    """Pairs with a significant, non-zero, two-sided measured effect.

    `Significant` and not `Regulated`: `Regulated` is `Significant AND EffectSize < 0`, which throws
    away every upward measurement and would leave the direction question with one answer by
    construction.
    """
    return [r for r in table if r["Significant"] == "TRUE" and float(r["EffectSize"]) != 0]


def in_reach(row: dict[str, str]) -> bool:
    """Could the scorer have seen this gene at all: the same arithmetic as loci.read_reach."""
    return abs(float(row["distanceToTSS"])) <= HALF_WINDOW


def overlapping(els: list[dict[str, Any]], starts: list[int], start: int, end: int) -> list[dict[str, Any]]:
    out = []
    j = bisect.bisect_left(starts, end) - 1
    while j >= 0 and els[j]["start"] > start - REACH:
        if els[j]["end"] > start:
            out.append(els[j])
        j -= 1
    return out


def model_value(els: list[dict[str, Any]], gene: str, cell: str) -> dict[str, Any] | None:
    """The sweep's signed predicted change for this gene in this cell, by the registered rule.

    Among every overlapping element whose top target is `gene` and that carries a value for `cell`,
    take the largest absolute predicted change — the model's strongest statement — breaking ties
    toward `predicted_coding`, then by element id. Every match is kept in `matches` so that a pair
    whose matches disagree in sign can be counted rather than hidden by the choice.
    """
    matches = []
    for e in els:
        for rank, (key, by_cell) in enumerate(
            (("predicted_coding", "predicted_coding_by_cell"), ("predicted", "predicted_by_cell"))
        ):
            p = e.get(key)
            if not p or p["gene"] != gene:
                continue
            v = (e.get(by_cell) or {}).get(cell)
            if v is not None:
                matches.append({"key": key, "value": float(v), "element": e["id"], "rank": rank})
    if not matches:
        return None
    best = min(matches, key=lambda m: (-abs(m["value"]), m["rank"], m["element"]))
    signs = {1 if m["value"] > 0 else (-1 if m["value"] < 0 else 0) for m in matches}
    return {**best, "matches": len(matches), "matches_disagree": len(signs) > 1}


def collect(table: list[dict[str, str]], cells: tuple[str, ...], elements: Path = ELEMENTS) -> dict[str, Any]:
    """Join the signed pairs of one arm to the sweep, one chromosome at a time.

    Returns the answerable pairs and, beside them, every denominator the coverage confound needs:
    an element the sweep already scored is not a random element, so the headline is meaningless
    without the strata it was drawn from.
    """
    pairs = signed(table)
    by_chrom: dict[str, list[dict[str, str]]] = defaultdict(list)
    for r in pairs:
        by_chrom[r["chrom"]].append(r)
    answered: list[dict[str, Any]] = []
    other: list[dict[str, Any]] = []
    strata: dict[str, int] = defaultdict(int)
    for chrom in sorted(by_chrom):
        p = elements / f"{chrom}.json"
        els = sorted(json.loads(p.read_text()), key=lambda e: e["start"]) if p.exists() else []
        starts = [e["start"] for e in els]
        for r in by_chrom[chrom]:
            cell = r["CellType"]
            strata[f"signed::{cell}"] += 1
            if cell not in cells:
                strata[f"refused_no_model_line::{cell}"] += 1
                continue
            if not in_reach(r):
                strata[f"out_of_reach::{cell}"] += 1
                continue
            ov = overlapping(els, starts, int(r["chromStart"]), int(r["chromEnd"]))
            strata[f"input_presence::{cell}"] += bool(ov)
            if not ov:
                strata[f"element_not_in_the_sweep::{cell}"] += 1
                continue
            m = model_value(ov, r["measuredGeneSymbol"], cell)
            effect = float(r["EffectSize"])
            row = {
                "cell": cell,
                "chrom": chrom,
                "element": [int(r["chromStart"]), int(r["chromEnd"])],
                "gene": r["measuredGeneSymbol"],
                "measured": round(effect, 4),
                "measured_sign": "down" if effect < 0 else "up",
                "distance_to_tss": int(float(r["distanceToTSS"])),
                "dataset": r["Dataset"],
            }
            if m is None:
                strata[f"covered_but_gene_not_the_top_target::{cell}"] += 1
                named = sorted(
                    {
                        e[k]["gene"]
                        for e in ov
                        for k in ("predicted_coding", "predicted")
                        if e.get(k) and e[k].get("gene")
                    }
                )
                other.append({**row, "model_named_instead": named[:6]})
                continue
            strata[f"answerable::{cell}"] += 1
            answered.append(
                {
                    **row,
                    "predicted": round(m["value"], 4),
                    **{
                        "predicted_sign": "down" if m["value"] < 0 else ("up" if m["value"] > 0 else "zero"),
                        "via": m["key"],
                        "matches": m["matches"],
                        "matches_disagree": m["matches_disagree"],
                    },
                }
            )
        del els, starts
    return {"answered": answered, "not_top_target": other, "strata": dict(strata), "signed": len(pairs)}


def marginal_down_rate(cells: tuple[str, ...], elements: Path = ELEMENTS) -> dict[str, Any]:
    """The constant-sign caller's score: what share of the whole sweep predicts a decrease.

    The comparison the headline is judged against. Taken over every element of every chromosome, by
    the same key precedence the join uses, so it is the rate a caller that ignored the measurement
    entirely and echoed the model's habit would achieve on a measured-down subset.
    """
    out: dict[str, Any] = {}
    counts: dict[str, dict[str, int]] = {c: defaultdict(int) for c in cells}
    total = 0
    for p in sorted(elements.glob("*.json")):
        els = json.loads(p.read_text())
        total += len(els)
        for e in els:
            for cell in cells:
                v = (e.get("predicted_coding_by_cell") or {}).get(cell)
                key = "predicted_coding"
                if v is None:
                    v = (e.get("predicted_by_cell") or {}).get(cell)
                    key = "predicted"
                if v is None:
                    counts[cell]["no_value"] += 1
                    continue
                counts[cell]["with_value"] += 1
                counts[cell][f"via_{key}"] += 1
                counts[cell]["down" if v < 0 else ("up" if v > 0 else "zero")] += 1
        del els
    for cell in cells:
        c = counts[cell]
        n = c["down"] + c["up"]
        out[cell] = {
            **wilson(c["down"], n),
            "zero": c["zero"],
            "no_value": c["no_value"],
            "counts": dict(c),
        }
    out["elements_in_the_sweep"] = total
    return out


def agreement(answered: list[dict[str, Any]]) -> dict[str, Any]:
    """Sign agreement, with the zero-prediction pairs held out of the denominator and counted."""
    usable = [r for r in answered if r["predicted_sign"] != "zero"]
    k = sum(r["predicted_sign"] == r["measured_sign"] for r in usable)
    return {
        **wilson(k, len(usable)),
        "predicted_zero_excluded": len(answered) - len(usable),
        "measured_down": sum(r["measured_sign"] == "down" for r in usable),
        "measured_up": sum(r["measured_sign"] == "up" for r in usable),
        "predicted_down": sum(r["predicted_sign"] == "down" for r in usable),
        "predicted_up": sum(r["predicted_sign"] == "up" for r in usable),
        "both_measured_signs_present": bool(
            any(r["measured_sign"] == "down" for r in usable)
            and any(r["measured_sign"] == "up" for r in usable)
        ),
    }


def judge(head: dict[str, Any], marginal: float | None) -> dict[str, Any]:
    """The registered rule, applied without discretion.

    PASS needs both margins. At or below 0.5 + 0.15 is FAIL outright. Above it but not clear of the
    constant-sign caller is UNDECIDABLE, which is registered as its own outcome precisely so that a
    high agreement on a one-sign subset cannot be written up as a success.
    """
    rate, ci = head.get("rate"), head.get("ci95")
    if rate is None or not ci:
        return {"verdict": "refused", "why": "no answerable pair with a non-zero predicted sign"}
    if not head["both_measured_signs_present"]:
        base = (
            "the answerable subset carries one measured sign only, so agreement is the model's own "
            "rate of emitting that sign"
        )
    else:
        base = "the answerable subset carries both measured signs"
    if rate <= CHANCE + CHANCE_MARGIN:
        return {
            "verdict": "failed",
            "why": f"agreement {rate} does not clear {CHANCE + CHANCE_MARGIN}; {base}",
        }
    if marginal is None:
        return {"verdict": "undecidable", "why": f"no marginal rate to compare against; {base}"}
    if rate - marginal >= MARGINAL_MARGIN and ci[0] > marginal:
        return {
            "verdict": "passed",
            "why": f"agreement {rate} clears both 0.65 and the marginal {round(marginal, 4)}; {base}",
        }
    return {
        "verdict": "undecidable",
        "why": (
            f"agreement {rate} clears 0.65 but not the constant-sign caller's {round(marginal, 4)} by "
            f"{MARGINAL_MARGIN} with the interval clear of it, so the direction claim is unsupported "
            f"rather than refuted; {base}"
        ),
    }


def direction(cells: tuple[str, ...] = ("K562", "GM12878"), elements: Path = ELEMENTS) -> dict[str, Any]:
    """The whole registered measurement. Held-out arm decides; the training arm is fitted-on company."""
    held = collect(rows(HELDOUT), cells, elements)
    train = collect(rows(TRAINING), cells, elements)
    marginal = marginal_down_rate(cells, elements)
    head = agreement(held["answered"])
    per_cell = {
        c: agreement([r for r in held["answered"] if r["cell"] == c])
        for c in cells
        if any(r["cell"] == c for r in held["answered"])
    }
    pooled_marginal = None
    usable = [r for r in held["answered"] if r["predicted_sign"] != "zero"]
    if usable:
        # the constant-sign caller weighted by the cells the answerable set actually draws from
        w = defaultdict(int)
        for r in usable:
            w[r["cell"]] += 1
        num = sum(w[c] * (marginal[c]["rate"] or 0) for c in w)
        pooled_marginal = num / sum(w.values())
    return {
        "result": "crispri_direction",
        "evidence": EVIDENCE,
        "preregistered": PREREGISTERED,
        "requests": 0,
        "cells": list(cells),
        "headline": head,
        "marginal_down_rate": marginal,
        "pooled_marginal_down_rate": None if pooled_marginal is None else round(pooled_marginal, 4),
        **judge(head, pooled_marginal),
        "per_cell": per_cell,
        "training_arm_fitted_on_not_a_test": {
            "agreement": agreement(train["answered"]),
            "strata": train["strata"],
            "note": (
                "the 2026-09-16 lane fitted logistic weights on this arm "
                "(data/results/crispri_benchmark.json), so it is reported for shape only and never "
                "decides the verdict"
            ),
        },
        "coverage": {
            "heldout": held["strata"],
            "heldout_signed_pairs": held["signed"],
            "note": (
                "an element the sweep already scored is not a random element (LOCI-BENCHMARK sections 19, "
                "21, 22): input_presence and every stratum between the signed pairs and the answerable "
                "subset are reported beside the headline, not after it"
            ),
        },
        "what_it_named_instead": {
            "pairs": len(held["not_top_target"]),
            "rows": held["not_top_target"],
        },
        "answered": held["answered"],
    }
