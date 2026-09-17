# SPDX-License-Identifier: AGPL-3.0-or-later
"""Does §17's one positive transfer, or is it a property of the episomal reporter?

Section 17 of GRAMMAR-BY-COMPARISON.md tested motif *arrangement* against motif *counts* on 51,376
lentiMPRA elements and failed the arrangement claim in all three cell lines. The one large positive it
left behind was unclaimed beyond that assay: strict, family-collapsed JASPAR site counts (relative
score 0.95, one site per TFClass family unit per position) lift Spearman with measured activity by
+0.234 in K562 and +0.235 in HepG2 over GC and CpG composition alone. A lentiMPRA element is 200 bp of
DNA on a reporter integrated out of its own chromosome; nothing in §17 says whether the counts read
enhancers or read that construct. This module asks the transfer question twice.

**Reading 1, VISTA: a different assay.** The VISTA Enhancer Browser has tested about 2,400 human
sequences in e11.5 mouse embryos (attribution/vista.py): each is positive in named tissues or negative,
in a living animal, in its own context. Six nested models are fitted on the training chromosomes and
scored once on the held-out ones (VISTA_HELD_OUT, fixed before any element was scored):

  length | + composition | + conservation | + counts | conservation + counts | + shuffled counts.

VISTA negatives are **not matched** to the positives -- they are sequences that failed to drive
expression, not a background drawn to match length, GC or constraint -- so a positive-against-negative
gap here is not an effect size, and the covariates of both classes are reported side by side so the
reader sees exactly what is unmatched. The central question is therefore not the gap but whether counts
add anything **over conservation**, since the elements were chosen for conservation in the first place.
Conservation is Zoonomia phyloP over 241 placental mammals (attribution/constraint.py), read from the
cached VISTA rows. The shuffled-counts model is the null: counts recomputed on a dinucleotide shuffle of
the same element, which keeps GC, CpG and every dinucleotide count and destroys only where the sites
are, so it measures how much of any count signal is composition wearing a family name.

**Reading 2, the real unknown: a different question, decided on measurement.** The 882
constrained-unknown blocks that are not copies (the organiser's "real unknown", area I's sharpest
attribution target) were measured on 2026-09-16 with an AlphaGenome deletion sweep and came out BELOW
the neutral tier per element, 0.248 against 0.293, z -5.7. Scoring motif counts against that sweep's own
call would be two readings of one model agreeing with each other, so the primary reading here is
**measured**: 161 of the 882 blocks hold at least one ENCODE4 lentiMPRA element, 227 elements in all,
and those elements carry a measured log2(RNA/DNA). Against the 387 elements inside neutral-tier blocks
that is a measured tier comparison, and §17's count model -- refitted with every element inside a
real-unknown or neutral block held out **by location** -- can be scored against measured activity inside
the real unknown. The predicted level over all 882 blocks, window by window, against 882 length-matched
neutral blocks, is reported as a prediction. Agreement with the deletion sweep's call is reported last,
labelled as agreement between two model readings and never as validation.

Evidence: VISTA status and tissue are `experimental`; phyloP is `experimental`; lentiMPRA activity is
`experimental`; a site is `predicted` (a matrix score, not a footprint); the sweep's lead is `predicted`
(an AlphaGenome deletion, and the same model names a target at 87% of matched random windows); every
comparison is `inferred`.
"""

from __future__ import annotations

import bisect
import json
import math
import random
from array import array
from collections import defaultdict
from pathlib import Path
from typing import Any

from genomeos.attribution.motif_grammar import (
    CELLS,
    INNER_FOLDS,
    LAMBDAS,
    SITE_THRESHOLD,
    accumulate,
    add,
    bootstrap_models,
    composition,
    fit,
    fold_assignment,
    interval,
    load_sites,
    predict,
    scan_sites,
    sites_path,
    spearman,
)
from genomeos.results import RESULTS_DIR, load_result, save_result

PREREGISTERED_VISTA = (
    "Strict family-collapsed JASPAR site counts per kb, the feature set §17 found predictive on "
    "lentiMPRA, are tested on VISTA's in-vivo transgenic assay. Decided on the held-out chromosomes "
    "(VISTA_HELD_OUT), by AUROC with 95% bootstrap intervals over elements: "
    "(1) AUROC(counts) above 0.5 with its interval above 0.5; "
    "(2) AUROC(counts) - AUROC(composition) above 0 with its interval above 0; "
    "(3) AUROC(conservation + counts) - AUROC(conservation) above 0 with its interval above 0 -- this is "
    "the transfer claim, that counts add over conservation on an in-vivo assay, and it is the central "
    "one because VISTA elements were chosen for conservation; "
    "(4) the null: AUROC(counts) - AUROC(shuffled counts) above 0 with its interval above 0, without "
    "which any count signal is composition under another name. The claim holds only if (1), (3) and (4) "
    "all hold. VISTA negatives are not matched to the positives on length, GC or constraint, so (1) is "
    "not an effect size and the covariates of both classes are reported beside it. Tissue group, "
    "reported separately and not part of the claim: among held-out positives a one-vs-rest count model "
    "ranks the elements of at least one tissue group above the other positives, interval above 0.5, "
    "against a label-permutation null."
)
PREREGISTERED_BLOCKS = (
    "Reading 2 is decided on MEASURED activity, not on the deletion sweep's own call, because that call "
    "is an AlphaGenome output and scoring a motif model against it would be two readings of one model. "
    "Coverage is reported first and nothing is invented for the blocks that have none: of the 882 "
    "real-unknown blocks, how many hold at least one ENCODE4 lentiMPRA element, and how many elements. "
    "(1) PRIMARY, measured level: lentiMPRA log2(RNA/DNA) of the elements inside real-unknown blocks "
    "against the elements inside neutral-tier blocks, per cell line, as an AUROC with a rank-sum z. The "
    "deletion sweep puts the real unknown BELOW the neutral tier per element (0.248 against 0.293, z "
    "-5.7); if the measured reporter puts it above, the two readings contradict each other and the "
    "contradiction is the finding, not split and not tuned away. "
    "(2) PRIMARY, measured transfer: Spearman of §17's (b) count model with measured activity on those "
    "same in-block elements against the (a) composition model, with every element inside a real-unknown "
    "or neutral-tier block held out of the fit BY LOCATION. The claim is that rho(b) - rho(a) is above 0 "
    "by its 95% bootstrap interval inside the real unknown, as it is +0.234 over the whole library. "
    "(3) SECONDARY, predicted level: the same count model's mean prediction over every 200 bp window of "
    "all 882 blocks against 882 length-matched neutral-tier blocks, each beside a dinucleotide-shuffled "
    "null of its own windows. This is a prediction and is reported as one. "
    "(4) SECONDARY and NOT VALIDATION: the AUROC of that predicted score against whether the sweep's "
    "AlphaGenome deletion moved a gene in the block. Both sides are model outputs, so agreement here is "
    "agreement between two model readings and is never read as either confirming the other."
)

VISTA_KNOWLEDGE = Path("data/knowledge/vista")
MPRA_KNOWLEDGE = Path("data/knowledge/mpra")
ELEMENT_TABLES = Path("data/knowledge/alphagenome/all_elements")
REFERENCE = Path("data/reference")
CACHE = Path("data/knowledge/motif_transfer")
CHROMS = tuple(f"chr{c}" for c in list(range(1, 23)) + ["X", "Y"])

# Held out before anything was scored: §17's four chromosomes (chr8, chr9, chr21, chr22) plus five
# more, because VISTA holds 2,223 elements against lentiMPRA's 51,376 and four chromosomes leave 233.
VISTA_HELD_OUT = ("chr4", "chr5", "chr8", "chr9", "chr13", "chr18", "chr21", "chr22", "chrX")
CHUNK = 1900  # the k-mer index packs positions in 11 bits, so no scanned piece may reach 2,048 bases
CHUNK_OVERLAP = 60  # wider than any JASPAR profile, so no site is lost at a chunk boundary
WINDOW = 200  # a lentiMPRA element, and therefore the unit the transferred model can read
MIN_ACGT = 0.9  # a window with more N than this is not scored; it is counted as dropped
SHUFFLE_WINDOWS = 20  # windows per block carried into the dinucleotide-shuffled null
COUNTS_PER_KB = True  # VISTA elements vary 180-fold in length; counts are a density, not a total
BOOTSTRAP = 1000
PERMUTATIONS = 1000
SEED = 20260917
MIN_LOG2 = 0.1  # the deletion sweep's own threshold for "this moved a gene"
PRIMARY_CELL = "K562"  # §17's largest count gain; all three cell lines are reported
GROUPS = ("neural", "heart", "limb and mesenchyme", "liver", "pancreas")
MIN_GROUP_TRAIN = 20  # positives of a group in the training set before it gets a model
MIN_GROUP_TEST = 8  # positives of a group among the held-out before its AUROC is reported

EVIDENCE = {
    "vista": "experimental: VISTA Enhancer Browser, transgenic mouse e11.5 (LBNL); hg38 coordinates",
    "conservation": (
        "experimental: Zoonomia phyloP over 241 placental mammals, per element mean and fraction of "
        "bases at or above 2.27"
    ),
    "activity": "experimental: ENCODE4 lentiMPRA (ENCSR106SZM), log2(RNA/DNA) per 200 bp element",
    "sites": (
        f"predicted: JASPAR 2026 CORE vertebrates at relative score >= {SITE_THRESHOLD}, collapsed to "
        "one site per TFClass family unit per position"
    ),
    "sweep_lead": (
        "predicted: an AlphaGenome deletion moved a gene by at least 0.1 log2; the same model names a "
        "target at 87% of matched random windows, so a lead is a lead and not a measurement"
    ),
    "comparison": (
        "inferred: ridge fitted on training chromosomes or on locations held out, AUROC and Spearman on "
        "the held-out ones, bootstrap over elements or blocks"
    ),
}

FEATURE_GROUPS = ("intercept", "length", "conservation", "composition", "counts", "counts_shuffled")
MODELS: dict[str, tuple[str, ...]] = {
    "length": ("intercept", "length"),
    "composition": ("intercept", "length", "composition"),
    "conservation": ("intercept", "length", "conservation"),
    "counts": ("intercept", "length", "composition", "counts"),
    "both": ("intercept", "length", "conservation", "composition", "counts"),
    "counts_shuffled": ("intercept", "length", "composition", "counts_shuffled"),
}
COMPARISONS = (
    ("composition", "length"),
    ("conservation", "length"),
    ("counts", "composition"),
    ("counts", "counts_shuffled"),
    ("both", "conservation"),
    ("both", "counts"),
)


# ---------------------------------------------------------------- statistics


def auroc(scores: list[float], labels: list[int]) -> float | None:
    """Area under the ROC curve by the rank sum, ties averaged: the chance level is 0.5."""
    n1 = sum(1 for y in labels if y)
    n0 = len(labels) - n1
    if not n1 or not n0:
        return None
    order = sorted(range(len(scores)), key=lambda i: scores[i])
    rank = [0.0] * len(scores)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and scores[order[j + 1]] == scores[order[i]]:
            j += 1
        r = (i + j) / 2 + 1
        for t in range(i, j + 1):
            rank[order[t]] = r
        i = j + 1
    s = sum(rank[i] for i in range(len(labels)) if labels[i])
    return (s - n1 * (n1 + 1) / 2) / (n1 * n0)


def auprc(scores: list[float], labels: list[int]) -> float | None:
    """Average precision, each tie group scored at the precision reached at the end of the group."""
    n1 = sum(1 for y in labels if y)
    if not n1 or n1 == len(labels):
        return None
    order = sorted(range(len(scores)), key=lambda i: -scores[i])
    total = 0.0
    tp = fp = 0
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and scores[order[j + 1]] == scores[order[i]]:
            j += 1
        hits = sum(1 for t in range(i, j + 1) if labels[order[t]])
        tp += hits
        fp += (j - i + 1) - hits
        total += hits * tp / (tp + fp)
        i = j + 1
    return total / n1


def _point(name: str, obs: dict[str, float | None]) -> float | None:
    if name in obs:
        return obs[name]
    a, b = name.split(" - ")
    if obs.get(a) is None or obs.get(b) is None:
        return None
    return obs[a] - obs[b]


def bootstrap_auc(
    labels: list[int],
    preds: dict[str, list[float]],
    comparisons=(),
    draws: int = BOOTSTRAP,
    seed: int = SEED,
) -> dict[str, Any]:
    """AUROC and AUPRC per model and the AUROC differences between models, resampling units."""
    rng = random.Random(seed)
    n = len(labels)
    names = list(preds)
    obs_roc = {m: auroc(preds[m], labels) for m in names}
    obs_prc = {m: auprc(preds[m], labels) for m in names}
    roc: dict[str, list[float]] = defaultdict(list)
    prc: dict[str, list[float]] = defaultdict(list)
    for _ in range(draws):
        idx = [rng.randrange(n) for _ in range(n)]
        y = [labels[i] for i in idx]
        got: dict[str, float | None] = {}
        for m in names:
            p = preds[m]
            resampled = [p[i] for i in idx]
            got[m] = auroc(resampled, y)
            a = auprc(resampled, y)
            if a is not None:
                prc[m].append(a)
        for m in names:
            if got[m] is not None:
                roc[m].append(got[m])
        for a_name, b_name in comparisons:
            if got.get(a_name) is not None and got.get(b_name) is not None:
                roc[f"{a_name} - {b_name}"].append(got[a_name] - got[b_name])
    out: dict[str, Any] = {}
    for key, vals in roc.items():
        if not vals:
            continue
        lo, hi = interval(vals)
        reference = 0.0 if " - " in key else 0.5
        out[key] = {
            "auroc": _round(_point(key, obs_roc)),
            "ci95": [round(lo, 4), round(hi, 4)],
            "above_chance": round(sum(1 for v in vals if v > reference) / len(vals), 3),
        }
        if key in prc and prc[key]:
            plo, phi = interval(prc[key])
            out[key]["auprc"] = _round(obs_prc[key])
            out[key]["auprc_ci95"] = [round(plo, 4), round(phi, 4)]
    return out


def permutation_null(
    scores: list[float], labels: list[int], draws: int = PERMUTATIONS, seed: int = SEED
) -> dict[str, Any]:
    """AUROC with the labels permuted: what the same score gets on a label carrying no information."""
    rng = random.Random(seed + 1)
    perm = list(labels)
    vals = []
    observed = auroc(scores, labels)
    for _ in range(draws):
        rng.shuffle(perm)
        a = auroc(scores, perm)
        if a is not None:
            vals.append(a)
    if not vals:
        return {"draws": 0}
    lo, hi = interval(vals)
    return {
        "draws": len(vals),
        "mean": round(sum(vals) / len(vals), 4),
        "ci95": [round(lo, 4), round(hi, 4)],
        "share_at_or_above_observed": (
            round(sum(1 for v in vals if v >= observed) / len(vals), 4) if observed is not None else None
        ),
    }


def _round(v: float | None, places: int = 4) -> float | None:
    return round(v, places) if v is not None else None


def _mean(xs: list[float | None]) -> float | None:
    vals = [x for x in xs if x is not None]
    return sum(vals) / len(vals) if vals else None


def _median(xs: list[float | None]) -> float | None:
    s = sorted(x for x in xs if x is not None)
    return s[len(s) // 2] if s else None


def mann_whitney(a: list[float], b: list[float]) -> dict[str, Any]:
    """A level comparison between two groups as an AUROC with a normal-approximation z (0.5 is equal)."""
    scores = list(a) + list(b)
    labels = [1] * len(a) + [0] * len(b)
    u = auroc(scores, labels)
    if u is None:
        return {"auroc": None, "n": [len(a), len(b)]}
    n1, n0 = len(a), len(b)
    sd = math.sqrt(n1 * n0 * (n1 + n0 + 1) / 12) / (n1 * n0)
    return {
        "auroc": round(u, 4),
        "z": round((u - 0.5) / sd, 2) if sd else None,
        "n": [n1, n0],
        "median": [_round(_median(a)), _round(_median(b))],
        "mean": [_round(_mean(a)), _round(_mean(b))],
    }


# ---------------------------------------------------------------- scanning long sequences


def chunk_sequence(seq: str, size: int = CHUNK, overlap: int = CHUNK_OVERLAP) -> list[tuple[int, str]]:
    """A long sequence as overlapping pieces under the index's position limit, with their offsets."""
    if len(seq) <= size:
        return [(0, seq)]
    step = size - overlap
    out = []
    start = 0
    while start < len(seq):
        out.append((start, seq[start : start + size]))
        if start + size >= len(seq):
            break
        start += step
    return out


def merge_sites(sites: list[tuple[str, int, int, int]]) -> list[tuple[str, int, int, int]]:
    """Sites from overlapping chunks in one frame: one site per family unit per position.

    `collapse` already did this inside each chunk with the matrix scores; across a chunk boundary the
    same site is seen twice, so the filter here is by coordinate alone, keeping the earlier of any two
    overlapping sites of one unit.
    """
    kept: dict[str, list[tuple[int, int]]] = defaultdict(list)
    out = []
    for u, s, e, d in sorted(set(sites), key=lambda x: (x[0], x[1], x[2], x[3])):
        if any(s < ke and ks < e for ks, ke in kept[u]):
            continue
        kept[u].append((s, e))
        out.append((u, s, e, d))
    return sorted(out, key=lambda x: (x[1], x[2], x[0]))


def unit_counts(sites: list[tuple[str, int, int, int]]) -> dict[str, int]:
    out: dict[str, int] = defaultdict(int)
    for s in sites:
        out[s[0]] += 1
    return dict(out)


def scan_counts(
    motifs: list, seqs: list[str], families: dict[str, dict[str, str]], batch_bases: int = 200_000
) -> list[dict[str, int]]:
    """Strict site counts per family unit for sequences of any length, chunked and scanned in batches."""
    pieces = [(i, off, sub) for i, s in enumerate(seqs) for off, sub in chunk_sequence(s)]
    found: list[list[tuple[str, int, int, int]]] = [[] for _ in seqs]
    batch: list[tuple[int, int, str]] = []
    bases = 0

    def flush() -> None:
        nonlocal batch, bases
        if not batch:
            return
        scanned = scan_sites(motifs, [p[2] for p in batch], families)
        for (i, off, _), st in zip(batch, scanned, strict=True):
            found[i].extend((u, s + off, e + off, d) for u, s, e, d in st)
        batch, bases = [], 0

    for p in pieces:
        batch.append(p)
        bases += len(p[2])
        if bases >= batch_bases:
            flush()
    flush()
    return [unit_counts(merge_sites(f)) for f in found]


def acgt_fraction(seq: str) -> float:
    s = seq.upper()
    return sum(s.count(b) for b in "ACGT") / len(s) if s else 0.0


# ---------------------------------------------------------------- §17's count model, refitted


def count_row(gc: list[float], counts: dict[str, int], units: list[str]) -> dict[int, float]:
    """§17's (b) row: [1, GC, GC squared, CpG o/e, log1p counts per unit, total, distinct families].

    The first four columns are §17's (a), so a fit on the leading four columns is the composition model.
    """
    row: dict[int, float] = {0: 1.0}
    for i, v in enumerate(gc):
        if v:
            row[1 + i] = v
    for i, u in enumerate(units):
        n = counts.get(u, 0)
        if n:
            row[4 + i] = math.log1p(n)
    total = sum(counts.values())
    if total:
        row[4 + len(units)] = math.log1p(total)
        row[5 + len(units)] = math.log1p(len(counts))
    return row


def count_columns(units: list[str]) -> int:
    return 6 + len(units)


def _sum_blocks(blocks, chroms: list[str], k: int, cells) -> tuple[array, dict[str, array]]:
    xtx = array("d", bytes(8 * k * k))
    xty = {cell: array("d", bytes(8 * k)) for cell in cells}
    for c in chroms:
        add(xtx, blocks[c][0])
        for cell in cells:
            add(xty[cell], blocks[c][1][cell])
    return xtx, xty


def mpra_tables(
    units: list[str], knowledge: Path = MPRA_KNOWLEDGE
) -> tuple[dict[str, list[dict[int, float]]], dict[str, dict[str, list[float]]], dict[str, list[dict]]]:
    """Every lentiMPRA element with a cached site scan and activity in all three lines, as (b) rows."""
    from genomeos.attribution.mpra import load_rows

    rows: dict[str, list[dict[int, float]]] = {}
    ys: dict[str, dict[str, list[float]]] = {}
    meta: dict[str, list[dict]] = {}
    for c in CHROMS:
        if not sites_path(c, knowledge).exists():
            continue
        by_key = {r["key"]: r for r in load_rows(c, knowledge)}
        els = [
            e
            for e in load_sites(c, knowledge)
            if e["key"] in by_key and all(x in by_key[e["key"]]["activity"] for x in CELLS)
        ]
        rows[c] = [count_row(e["gc"], unit_counts([tuple(s) for s in e["sites"]]), units) for e in els]
        ys[c] = {cell: [by_key[e["key"]]["activity"][cell] for e in els] for cell in CELLS}
        meta[c] = [
            {
                "key": e["key"],
                "chrom": c,
                "start": by_key[e["key"]]["start"],
                "end": by_key[e["key"]]["end"],
                "gc": e["gc"][0],
                "sites": sum(unit_counts([tuple(s) for s in e["sites"]]).values()),
            }
            for e in els
        ]
    return rows, ys, meta


def choose_mpra_lambda(blocks, rows, ys, folds: dict[str, int], k: int, cols: int) -> dict[str, Any]:
    """The penalty per cell line by mean fold-by-chromosome Spearman inside the training set."""
    train = list(folds)
    total_x, total_y = _sum_blocks(blocks, train, k, CELLS)
    scores: dict[str, dict[float, list[float]]] = {cell: defaultdict(list) for cell in CELLS}
    for f in range(max(folds.values()) + 1):
        te = [c for c in train if folds[c] == f]
        xtx = array("d", total_x)
        xty = {cell: array("d", total_y[cell]) for cell in CELLS}
        for c in te:
            add(xtx, blocks[c][0], -1.0)
            for cell in CELLS:
                add(xty[cell], blocks[c][1][cell], -1.0)
        xs = [r for c in te for r in rows[c]]
        for lam in LAMBDAS:
            for cell in CELLS:
                w = fit(xtx, xty[cell], k, cols, lam)
                if w is None:
                    continue
                rho = spearman(predict(w, xs), [v for c in te for v in ys[c][cell]])
                if rho is not None:
                    scores[cell][lam].append(rho)
    out = {}
    for cell, rec in scores.items():
        means = {lam: sum(v) / len(v) for lam, v in rec.items()}
        if not means:
            continue
        lam = max(means, key=lambda x: means[x])
        out[cell] = {
            "lambda": lam,
            "inner_spearman": round(means[lam], 4),
            "by_lambda": {str(x): round(v, 4) for x, v in means.items()},
        }
    return out


def mpra_model(
    units: list[str],
    knowledge: Path = MPRA_KNOWLEDGE,
    exclude: frozenset[str] = frozenset(),
    progress=None,
) -> dict[str, Any]:
    """§17's (a) and (b) models refitted with a set of element keys held out by location.

    The vocabulary is §17's own 40 family units, so nothing is selected here; the penalty is chosen by
    the same fold-by-chromosome cross-validation inside whatever remains after the exclusion.
    """
    say = progress or (lambda _m: None)
    rows, ys, meta = mpra_tables(units, knowledge)
    k = count_columns(units)
    keep: dict[str, list[int]] = {
        c: [i for i, m in enumerate(meta[c]) if m["key"] not in exclude] for c in rows
    }
    train_rows = {c: [rows[c][i] for i in keep[c]] for c in rows if keep[c]}
    train_ys = {c: {cell: [ys[c][cell][i] for i in keep[c]] for cell in CELLS} for c in rows if keep[c]}
    blocks = {c: accumulate(train_rows[c], train_ys[c], k) for c in train_rows}
    folds = fold_assignment({c: len(train_rows[c]) for c in train_rows}, INNER_FOLDS)
    chosen = {
        name: choose_mpra_lambda(blocks, train_rows, train_ys, folds, k, cols)
        for name, cols in (("a", 4), ("b", k))
    }
    xtx, xty = _sum_blocks(blocks, list(train_rows), k, CELLS)
    weights: dict[str, dict[str, list[float] | None]] = {}
    for name, cols in (("a", 4), ("b", k)):
        weights[name] = {}
        for cell in CELLS:
            lam = chosen[name].get(cell, {}).get("lambda", LAMBDAS[0])
            w = fit(xtx, xty[cell], k, cols, lam)
            weights[name][cell] = [round(v, 8) for v in w] if w else None
        say(
            f"  ({name}) inner spearman "
            + ", ".join(f"{c} {v['inner_spearman']}" for c, v in chosen[name].items())
        )
    return {
        "family_units": units,
        "columns": {"a": 4, "b": k},
        "elements_available": sum(len(m) for m in meta.values()),
        "elements_excluded_by_location": sum(len(meta[c]) - len(keep[c]) for c in rows),
        "training_elements": sum(len(r) for r in train_rows.values()),
        "training_chromosomes": sorted(train_rows, key=lambda c: (len(c), c)),
        "inner_cross_validation": chosen,
        "weights": weights,
        "rows": rows,
        "activity": ys,
        "meta": meta,
    }


# ---------------------------------------------------------------- reading 1: VISTA


def vista_rows(chrom: str, knowledge: Path = VISTA_KNOWLEDGE) -> list[dict[str, Any]]:
    """The cached VISTA annotation of one chromosome (attribution/vista.py wrote it; no network here)."""
    p = knowledge / f"rows_{chrom}.json"
    return json.loads(p.read_text()) if p.exists() else []


def vista_counts_path(chrom: str, cache: Path = CACHE) -> Path:
    return cache / f"vista_counts_{chrom}.json"


def build_vista_counts(
    chrom: str,
    knowledge: Path = VISTA_KNOWLEDGE,
    reference: Path = REFERENCE,
    cache: Path = CACHE,
    seed: int = SEED,
) -> int:
    """Strict site counts per family unit for every VISTA element of one chromosome, real and shuffled."""
    from genomeos.coords import Locus
    from genomeos.genome.genome import Genome
    from genomeos.genome.motifs import dinucleotide_shuffle, load_families, load_motifs

    rows = vista_rows(chrom, knowledge)
    fa = reference / f"{chrom}.fa.gz"
    if not rows or not fa.exists():
        return 0
    genome = Genome.from_fasta(fa)
    seqs = [str(genome.fetch(Locus(chrom, r["start"], r["end"]))).upper() for r in rows]
    rng = random.Random(seed)
    shuf = [dinucleotide_shuffle(s, rng) for s in seqs]
    motifs, families = load_motifs(relative=SITE_THRESHOLD), load_families()
    real = scan_counts(motifs, seqs, families)
    fake = scan_counts(motifs, shuf, families)
    payload = {
        "chrom": chrom,
        "threshold": SITE_THRESHOLD,
        "elements": [
            {
                "id": r["id"],
                "length": r["length"],
                "acgt": round(acgt_fraction(s), 4),
                "gc": composition(s),
                "gc_shuffled": composition(t),
                "counts": a,
                "counts_shuffled": b,
            }
            for r, s, t, a, b in zip(rows, seqs, shuf, real, fake, strict=True)
        ],
    }
    p = vista_counts_path(chrom, cache)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(payload))
    return len(rows)


def load_vista_counts(chrom: str, cache: Path = CACHE) -> dict[str, dict[str, Any]]:
    p = vista_counts_path(chrom, cache)
    return {e["id"]: e for e in json.loads(p.read_text())["elements"]} if p.exists() else {}


def _scaled_counts(counts: dict[str, int], units: list[str], per_kb: float) -> list[float]:
    out = [math.log1p(counts.get(u, 0) * per_kb) for u in units]
    out.append(math.log1p(sum(counts.values()) * per_kb))
    out.append(math.log1p(len(counts)))
    return out


def vista_feature_groups(
    row: dict[str, Any], element: dict[str, Any], units: list[str], fallback: dict[str, float]
) -> dict[str, list[float]]:
    """The six named feature blocks of one VISTA element; conservation falls back to the training mean."""
    length = max(element["length"], 1)
    per_kb = 1000.0 / length if COUNTS_PER_KB else 1.0
    phylop = row.get("phylop_mean")
    frac = row.get("constrained_fraction")
    return {
        "intercept": [1.0],
        "length": [math.log(length)],
        "conservation": [
            phylop if phylop is not None else fallback["phylop_mean"],
            frac if frac is not None else fallback["constrained_fraction"],
        ],
        "composition": list(element["gc"]),
        "counts": _scaled_counts(element["counts"], units, per_kb),
        "counts_shuffled": _scaled_counts(element["counts_shuffled"], units, per_kb),
    }


def model_row(groups: dict[str, list[float]], model: tuple[str, ...]) -> dict[int, float]:
    row: dict[int, float] = {}
    i = 0
    for name in model:
        for v in groups[name]:
            if v:
                row[i] = v
            i += 1
    return row


def model_columns(groups: dict[str, list[float]], model: tuple[str, ...]) -> int:
    return sum(len(groups[name]) for name in model)


def choose_penalty(
    blocks: dict[str, tuple[array, dict[str, array]]],
    rows: dict[str, list[dict[int, float]]],
    labels: dict[str, list[float]],
    folds: dict[str, int],
    k: int,
) -> dict[str, Any]:
    """The penalty with the best mean fold-by-chromosome AUROC inside the training chromosomes."""
    train = list(folds)
    total_x, total_y = _sum_blocks(blocks, train, k, ("y",))
    means: dict[float, list[float]] = defaultdict(list)
    for f in range(max(folds.values()) + 1):
        te = [c for c in train if folds[c] == f]
        xtx = array("d", total_x)
        xty = array("d", total_y["y"])
        for c in te:
            add(xtx, blocks[c][0], -1.0)
            add(xty, blocks[c][1]["y"], -1.0)
        xs = [r for c in te for r in rows[c]]
        ys = [int(v) for c in te for v in labels[c]]
        for lam in LAMBDAS:
            w = fit(xtx, xty, k, k, lam)
            if w is None:
                continue
            a = auroc(predict(w, xs), ys)
            if a is not None:
                means[lam].append(a)
    scored = {lam: sum(v) / len(v) for lam, v in means.items() if v}
    if not scored:
        return {"lambda": LAMBDAS[0], "inner_auroc": None, "by_lambda": {}}
    lam = max(scored, key=lambda x: scored[x])
    return {
        "lambda": lam,
        "inner_auroc": round(scored[lam], 4),
        "by_lambda": {str(x): round(v, 4) for x, v in scored.items()},
    }


def fit_and_score(
    rows: dict[str, list[dict[int, float]]],
    labels: dict[str, list[float]],
    train: list[str],
    test: list[str],
    k: int,
) -> tuple[dict[str, Any], list[float]]:
    """One feature set: penalty inside the training chromosomes, predictions on the held-out ones."""
    blocks = {c: accumulate(rows[c], {"y": labels[c]}, k) for c in train + test}
    folds = fold_assignment({c: len(rows[c]) for c in train}, min(INNER_FOLDS, len(train)))
    chosen = choose_penalty(blocks, rows, labels, folds, k)
    xtx, xty = _sum_blocks(blocks, train, k, ("y",))
    w = fit(xtx, xty["y"], k, k, chosen["lambda"])
    preds = predict(w, [r for c in test for r in rows[c]]) if w else []
    return chosen, preds


def vista_dataset(
    units: list[str],
    knowledge: Path = VISTA_KNOWLEDGE,
    cache: Path = CACHE,
    held_out=VISTA_HELD_OUT,
) -> dict[str, Any]:
    """Every VISTA element with a cached site count, its label, its feature blocks and the coverage."""
    chroms = [c for c in CHROMS if vista_counts_path(c, cache).exists()]
    coverage = {
        "chromosomes_with_a_vista_table": sum(1 for c in CHROMS if vista_rows(c, knowledge)),
        "chromosomes_scanned": len(chroms),
        "elements_in_table": 0,
        "without_sequence_or_sites": 0,
        "without_conservation": 0,
        "status_not_positive_or_negative": 0,
        "used": 0,
    }
    raw: dict[str, list[tuple[dict, dict]]] = {}
    phylop: list[float] = []
    frac: list[float] = []
    for c in chroms:
        counts = load_vista_counts(c, cache)
        got = []
        for r in vista_rows(c, knowledge):
            coverage["elements_in_table"] += 1
            if r["status"] not in ("positive", "negative"):
                coverage["status_not_positive_or_negative"] += 1
                continue
            e = counts.get(r["id"])
            if e is None:
                coverage["without_sequence_or_sites"] += 1
                continue
            if r.get("phylop_mean") is None or r.get("constrained_fraction") is None:
                coverage["without_conservation"] += 1
            elif c not in held_out:
                phylop.append(r["phylop_mean"])
                frac.append(r["constrained_fraction"])
            got.append((r, e))
        raw[c] = got
    fallback = {"phylop_mean": _mean(phylop) or 0.0, "constrained_fraction": _mean(frac) or 0.0}
    groups: dict[str, list[dict[str, list[float]]]] = {}
    labels: dict[str, list[float]] = {}
    meta: dict[str, list[dict[str, Any]]] = {}
    for c, got in raw.items():
        groups[c] = [vista_feature_groups(r, e, units, fallback) for r, e in got]
        labels[c] = [1.0 if r["status"] == "positive" else 0.0 for r, _ in got]
        meta[c] = [
            {
                "id": r["id"],
                "length": r["length"],
                "status": r["status"],
                "tissue_groups": sorted(set(r.get("groups") or ())),
                "acgt": e["acgt"],
                "gc": e["gc"][0],
                "cpg_oe": e["gc"][2],
                "sites_per_kb": round(sum(e["counts"].values()) * 1000 / max(r["length"], 1), 3),
                "phylop_mean": r.get("phylop_mean"),
                "constrained_fraction": r.get("constrained_fraction"),
            }
            for r, e in got
        ]
        coverage["used"] += len(got)
    return {
        "chromosomes": chroms,
        "groups": groups,
        "labels": labels,
        "meta": meta,
        "coverage": coverage,
        "conservation_fallback": {k: round(v, 4) for k, v in fallback.items()},
    }


def covariates(meta: list[dict[str, Any]]) -> dict[str, Any]:
    """What is unmatched between VISTA positives and negatives, stated rather than assumed away."""
    out: dict[str, Any] = {}
    for status in ("positive", "negative"):
        rs = [m for m in meta if m["status"] == status]
        out[status] = {
            "elements": len(rs),
            "median_length": _median([float(m["length"]) for m in rs]),
            "median_gc": _round(_median([m["gc"] for m in rs])),
            "median_cpg_oe": _round(_median([m["cpg_oe"] for m in rs])),
            "median_phylop_mean": _round(_median([m["phylop_mean"] for m in rs])),
            "median_constrained_fraction": _round(_median([m["constrained_fraction"] for m in rs])),
            "median_sites_per_kb": _round(_median([m["sites_per_kb"] for m in rs])),
        }
    labels = [1 if m["status"] == "positive" else 0 for m in meta]
    out["positives_above_negatives_auroc"] = {
        key: _round(auroc([float(m[field]) for m in meta], labels))
        for key, field in (
            ("length", "length"),
            ("gc", "gc"),
            ("cpg_oe", "cpg_oe"),
            ("phylop_mean", "phylop_mean"),
            ("constrained_fraction", "constrained_fraction"),
            ("sites_per_kb", "sites_per_kb"),
        )
        if all(m.get(field) is not None for m in meta)
    }
    out["reading"] = (
        "VISTA negatives are sequences that failed to drive expression in the assay, not a background "
        "matched to the positives on any of these, so a separation between the classes is not an effect "
        "size; the column that matters is whether counts add over conservation"
    )
    return out


TISSUE_MODELS = ("length", "composition", "counts", "both")
TISSUE_COMPARISONS = (("counts", "composition"), ("both", "conservation"))


def vista_tissue_reading(
    data: dict[str, Any], train: list[str], test: list[str], models=TISSUE_MODELS
) -> dict[str, Any]:
    """Among positives only: does a one-vs-rest count model rank a tissue group above other positives?

    The same nested sets as the status reading, so a group AUROC can be read against what length and
    composition alone already get: a group whose elements are simply longer or GC-richer than the rest
    would be predicted by the baselines too.
    """
    groups, labels, meta = data["groups"], data["labels"], data["meta"]
    pos: dict[str, list[int]] = {
        c: [i for i, y in enumerate(labels[c]) if y] for c in groups if c in train + test
    }
    pos = {c: v for c, v in pos.items() if v}
    if not pos:
        return {}
    sample = next(iter(groups.values()))[0]
    rows = {
        name: {c: [model_row(groups[c][i], MODELS[name]) for i in pos[c]] for c in pos} for name in models
    }
    out: dict[str, Any] = {
        "feature_sets": {name: list(MODELS[name]) for name in models},
        "positives_train": sum(len(pos[c]) for c in train if c in pos),
        "positives_held_out": sum(len(pos[c]) for c in test if c in pos),
        "positives_without_a_mapped_group": sum(
            1 for c in pos for i in pos[c] if not meta[c][i]["tissue_groups"]
        ),
        "minimum_positives": {"train": MIN_GROUP_TRAIN, "held_out": MIN_GROUP_TEST},
        "by_group": {},
    }
    for g in GROUPS:
        y = {c: [1.0 if g in meta[c][i]["tissue_groups"] else 0.0 for i in pos[c]] for c in pos}
        n_train = int(sum(sum(y[c]) for c in train if c in y))
        n_test = int(sum(sum(y[c]) for c in test if c in y))
        rec: dict[str, Any] = {"positives_train": n_train, "positives_held_out": n_test}
        if n_train < MIN_GROUP_TRAIN or n_test < MIN_GROUP_TEST:
            rec["skipped"] = "too few positives of this group"
            out["by_group"][g] = rec
            continue
        tr = [c for c in train if c in y]
        te = [c for c in test if c in y]
        preds: dict[str, list[float]] = {}
        for name in models:
            chosen, p = fit_and_score(rows[name], y, tr, te, model_columns(sample, MODELS[name]))
            preds[name] = p
            rec[f"lambda_{name}"] = chosen["lambda"]
            rec[f"inner_auroc_{name}"] = chosen["inner_auroc"]
        truth = [int(v) for c in te for v in y[c]]
        scored = bootstrap_auc(truth, preds, (("counts", "composition"),))
        rec["by_feature_set"] = scored
        rec.update({k: v for k, v in scored["counts"].items()})
        rec["permutation_null"] = permutation_null(preds["counts"], truth)
        out["by_group"][g] = rec
    return out


def run_vista(
    knowledge: Path = VISTA_KNOWLEDGE,
    cache: Path = CACHE,
    results_dir: Path = RESULTS_DIR,
    held_out=VISTA_HELD_OUT,
    score_held_out: bool = True,
    progress=None,
) -> dict[str, Any]:
    """Reading 1: the six nested models on VISTA, fitted on the training chromosomes, scored once."""
    say = progress or (lambda _m: None)
    grammar = load_result("motif_grammar", results_dir)
    if not grammar:
        raise FileNotFoundError("no motif_grammar result: §17's vocabulary is what transfers")
    units = list(grammar["family_units"])
    data = vista_dataset(units, knowledge, cache, held_out)
    chroms = data["chromosomes"]
    train = [c for c in chroms if c not in held_out and data["labels"][c]]
    test = [c for c in chroms if c in held_out and data["labels"][c]]
    n_train = sum(len(data["labels"][c]) for c in train)
    n_test = sum(len(data["labels"][c]) for c in test)
    say(f"  VISTA: {n_train} training elements, {n_test} held out over {len(test)} chromosomes")
    sample = next(iter(data["groups"].values()))[0]
    out: dict[str, Any] = {
        "preregistered": PREREGISTERED_VISTA,
        "site_threshold": SITE_THRESHOLD,
        "family_units": units,
        "counts_per_kb": COUNTS_PER_KB,
        "training_chromosomes": train,
        "held_out_chromosomes": test,
        "training_elements": n_train,
        "held_out_elements": n_test,
        "coverage": data["coverage"],
        "conservation_fallback": data["conservation_fallback"],
        "models": {m: list(g) for m, g in MODELS.items()},
        "columns": {m: model_columns(sample, g) for m, g in MODELS.items()},
        "positives": {
            "train": int(sum(sum(data["labels"][c]) for c in train)),
            "held_out": int(sum(sum(data["labels"][c]) for c in test)),
        },
        "covariates_all": covariates([m for c in chroms for m in data["meta"][c]]),
        "covariates_held_out": covariates([m for c in test for m in data["meta"][c]]),
        "evidence": EVIDENCE,
    }
    inner: dict[str, Any] = {}
    preds: dict[str, list[float]] = {}
    for name, model in MODELS.items():
        k = model_columns(sample, model)
        rows = {c: [model_row(g, model) for g in data["groups"][c]] for c in chroms}
        chosen, p = fit_and_score(rows, data["labels"], train, test if score_held_out else [], k)
        inner[name] = chosen
        if score_held_out:
            preds[name] = p
        say(f"  ({name}) inner AUROC {chosen['inner_auroc']}, lambda {chosen['lambda']}")
    out["inner_cross_validation"] = inner
    if not score_held_out or not test:
        out["held_out"] = None
        return out
    truth = [int(v) for c in test for v in data["labels"][c]]
    out["held_out"] = bootstrap_auc(truth, preds, COMPARISONS)
    out["permutation_null"] = {m: permutation_null(preds[m], truth) for m in ("counts", "both")}
    out["tissue"] = vista_tissue_reading(data, train, test)
    out["verdict"] = vista_verdict(out)
    for name in MODELS:
        say(f"  held out ({name}) AUROC {out['held_out'][name]['auroc']} {out['held_out'][name]['ci95']}")
    return out


def vista_verdict(out: dict[str, Any]) -> dict[str, Any]:
    """The pre-registered VISTA claim, decided by the four intervals it named."""
    h = out["held_out"]

    def above(key: str, reference: float) -> bool:
        return key in h and h[key]["ci95"][0] > reference

    separates = above("counts", 0.5)
    over_composition = above("counts - composition", 0.0)
    over_conservation = above("both - conservation", 0.0)
    over_shuffle = above("counts - counts_shuffled", 0.0)
    holds = separates and over_conservation and over_shuffle
    tissue = out.get("tissue") or {}
    scored = [r for r in (tissue.get("by_group") or {}).values() if "ci95" in r]
    any_group = any((r.get("ci95") or [0.0])[0] > 0.5 for r in scored)
    counts_add_to_group = any(
        ((r.get("by_feature_set") or {}).get("counts - composition") or {}).get("ci95", [0.0])[0] > 0
        for r in scored
    )
    return {
        "counts_separate_positives_from_negatives": separates,
        "counts_add_over_composition": over_composition,
        "counts_add_over_conservation": over_conservation,
        "counts_beat_their_shuffled_control": over_shuffle,
        "preregistered_claim_holds": holds,
        "a_tissue_group_is_predicted": any_group,
        "counts_add_over_composition_for_the_tissue_group": counts_add_to_group,
        "tissue_statement": (
            "the tissue group is predictable from sequence and counts add to it over composition"
            if any_group and counts_add_to_group
            else "the tissue group is predictable from sequence, but length and composition alone get "
            "there: motif counts add nothing detectable to which tissue a positive drives"
            if any_group
            else "no tissue group is predicted above chance"
        ),
        "statement": (
            "strict family-collapsed motif counts separate VISTA positives from negatives in vivo and "
            "add over conservation: §17's positive transfers to a different assay"
            if holds
            else "counts separate the classes but do not add over conservation, or do not beat their own "
            "dinucleotide-shuffled control: §17's positive does not transfer as a claim about enhancers "
            "in vivo"
            if separates
            else "strict family-collapsed motif counts do not separate VISTA positives from negatives: "
            "§17's positive does not transfer to an in-vivo assay"
        ),
    }


# ---------------------------------------------------------------- reading 2: the real unknown


def moves(element: dict[str, Any], min_log2: float = MIN_LOG2) -> bool:
    """Did deleting this element move a gene, by the deletion sweep's own threshold?"""
    pred = element.get("predicted") or {}
    log2 = pred.get("log2_fold_change")
    return bool(pred.get("gene")) and log2 is not None and abs(log2) >= min_log2


def scored_elements(chrom: str, tables: Path = ELEMENT_TABLES) -> list[dict[str, Any]]:
    p = tables / f"{chrom}.json"
    return sorted(json.loads(p.read_text()), key=lambda e: e["start"]) if p.exists() else []


def mpra_positions(knowledge: Path = MPRA_KNOWLEDGE) -> dict[str, list[tuple[int, str]]]:
    """Where every lentiMPRA element sits, by midpoint: enough to place one inside a block."""
    from genomeos.attribution.mpra import load_rows

    out: dict[str, list[tuple[int, str]]] = {}
    for c in CHROMS:
        if not sites_path(c, knowledge).exists():
            continue
        out[c] = sorted(((r["start"] + r["end"]) // 2, r["key"]) for r in load_rows(c, knowledge))
    return out


def tier_blocks(
    chroms=CHROMS,
    results_dir: Path = RESULTS_DIR,
    tables: Path = ELEMENT_TABLES,
    positions: dict[str, list[tuple[int, str]]] | None = None,
) -> dict[str, list[dict[str, Any]]]:
    """The real-unknown blocks and the neutral tier, with the sweep's call and lentiMPRA coverage attached.

    The real unknown is the organiser's constrained-unknown tier with the copies taken out -- the same
    882 blocks the peer measured. The neutral tier is left as the peer left it, copies included, because
    it is their control that this reading is held against.
    """
    from genomeos.attribution import organise

    out: dict[str, list[dict[str, Any]]] = {"real_unknown": [], "neutral": []}
    for chrom in chroms:
        try:
            blocks = organise.blocks(chrom, results_dir)
        except FileNotFoundError:
            continue
        els = scored_elements(chrom, tables)
        mpra = (positions or {}).get(chrom, [])
        mids = [m[0] for m in mpra]
        for b in blocks:
            tier = b.get("tier")
            if tier == "constrained_unknown" and not b.get("copy"):
                key = "real_unknown"
            elif tier == "neutral":
                key = "neutral"
            else:
                continue
            inside = [e for e in els if b["start"] <= (e["start"] + e["end"]) // 2 < b["end"]]
            i = bisect.bisect_left(mids, b["start"])
            j = bisect.bisect_left(mids, b["end"])
            out[key].append(
                {
                    "chrom": chrom,
                    "block": f"{chrom}:{b['start']}-{b['end']}",
                    "start": b["start"],
                    "end": b["end"],
                    "length": b["length"],
                    "case": b.get("case"),
                    "tested_elements": len(inside),
                    "moving_elements": sum(1 for e in inside if moves(e)),
                    "mpra_keys": [k for _, k in mpra[i:j]],
                }
            )
    return out


def block_coverage(tiers: dict[str, list[dict[str, Any]]]) -> dict[str, Any]:
    """How much measured and how much predicted coverage each tier has, before anything is scored."""
    return {
        tier: {
            "blocks": len(blocks),
            "blocks_with_a_measured_element": sum(1 for b in blocks if b["mpra_keys"]),
            "blocks_without_any_measured_element": sum(1 for b in blocks if not b["mpra_keys"]),
            "measured_elements": sum(len(b["mpra_keys"]) for b in blocks),
            "blocks_with_an_element_the_sweep_tested": sum(1 for b in blocks if b["tested_elements"]),
            "blocks_with_a_sweep_lead": sum(1 for b in blocks if b["moving_elements"]),
            "median_length_kb": _round((_median([float(b["length"]) for b in blocks]) or 0) / 1000, 1),
        }
        for tier, blocks in tiers.items()
    }


def measured_reading(
    model: dict[str, Any],
    tiers: dict[str, list[dict[str, Any]]],
    cells=CELLS,
) -> dict[str, Any]:
    """The primary, measured reading: lentiMPRA activity inside the blocks, and what the model says of it.

    Every element inside a real-unknown or a neutral block was held out of the fit by location, so the
    predictions here are out of sample; the measured activity itself needs no model at all.
    """
    meta, activity = model["meta"], model["activity"]
    index: dict[str, tuple[str, int]] = {}
    for c, ms in meta.items():
        for i, m in enumerate(ms):
            index[m["key"]] = (c, i)
    picked: dict[str, list[tuple[str, int]]] = {}
    for tier, blocks in tiers.items():
        seen: list[tuple[str, int]] = []
        for b in blocks:
            seen.extend(index[k] for k in b["mpra_keys"] if k in index)
        picked[tier] = seen
    out: dict[str, Any] = {
        "coverage": block_coverage(tiers),
        "measured_elements_with_activity_in_all_three_lines": {t: len(v) for t, v in picked.items()},
        "elements_excluded_by_location": model["elements_excluded_by_location"],
        "training_elements": model["training_elements"],
        "gc_of_measured_elements": {
            tier: _round(_median([meta[c][i]["gc"] for c, i in picked[tier]])) for tier in picked
        },
        "sites_of_measured_elements": {
            tier: _median([float(meta[c][i]["sites"]) for c, i in picked[tier]]) for tier in picked
        },
    }
    if not picked["real_unknown"] or not picked["neutral"]:
        out["level"] = None
        out["transfer"] = None
        return out
    out["level"] = {
        cell: mann_whitney(
            [activity[c][cell][i] for c, i in picked["real_unknown"]],
            [activity[c][cell][i] for c, i in picked["neutral"]],
        )
        for cell in cells
    }
    out["transfer"] = _transfer(model, picked["real_unknown"], cells)
    out["transfer_neutral_tier"] = _transfer(model, picked["neutral"], cells)
    return out


def _transfer(model: dict[str, Any], picked: list[tuple[str, int]], cells) -> dict[str, Any]:
    """Spearman of (a) and (b) with measured activity on a set of elements, and their difference."""
    rows, activity = model["rows"], model["activity"]
    out: dict[str, Any] = {}
    for cell in cells:
        preds: dict[str, list[float]] = {}
        for name in ("a", "b"):
            w = model["weights"][name].get(cell)
            if w is None:
                continue
            p = predict(w, [rows[c][i] for c, i in picked])
            if len(set(p)) < 2:  # a constant prediction has no rank correlation to report
                continue
            preds[name] = p
        truth = [activity[c][cell][i] for c, i in picked]
        if len(preds) < 2 or len(truth) < 10:
            out[cell] = {"elements": len(truth), "skipped": "too few elements or a constant prediction"}
            continue
        out[cell] = bootstrap_models(truth, preds)
        out[cell]["elements"] = len(truth)
    return out


def match_by_length(target: list[dict[str, Any]], pool: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """One control block per target block, nearest in log length, each used once (largest target first).

    A block's length drives how many windows it holds and how many elements the sweep tested in it, so
    the control has to carry the same length distribution or it is not a control.
    """
    free = sorted(pool, key=lambda b: b["length"])
    lengths = [math.log(max(b["length"], 1)) for b in free]
    used = [False] * len(free)
    out = []
    for t in sorted(target, key=lambda b: -b["length"]):
        want = math.log(max(t["length"], 1))
        best, gap = None, None
        for i, ln in enumerate(lengths):
            if used[i]:
                continue
            d = abs(ln - want)
            if gap is None or d < gap:
                best, gap = i, d
        if best is None:
            break
        used[best] = True
        out.append(free[best])
    return out


def windows_of(seq: str, size: int = WINDOW, min_acgt: float = MIN_ACGT) -> tuple[list[str], int]:
    """Non-overlapping windows of one block, and how many were dropped for N bases."""
    kept, dropped = [], 0
    for i in range(0, len(seq) - size + 1, size):
        w = seq[i : i + size].upper()
        if acgt_fraction(w) < min_acgt:
            dropped += 1
            continue
        kept.append(w)
    return kept, dropped


def score_windows(
    windows: list[str],
    motifs: list,
    families: dict[str, dict[str, str]],
    units: list[str],
    weights: dict[str, list[float] | None],
) -> dict[str, list[float]]:
    """The transferred count model's prediction per window, per cell line."""
    sites = scan_sites(motifs, windows, families)
    rows = [
        count_row(composition(w), unit_counts([tuple(x) for x in st]), units)
        for w, st in zip(windows, sites, strict=True)
    ]
    return {cell: predict(w, rows) for cell, w in weights.items() if w}


def score_blocks(
    blocks: list[dict[str, Any]],
    units: list[str],
    weights: dict[str, list[float] | None],
    reference: Path = REFERENCE,
    shuffle_windows: int = SHUFFLE_WINDOWS,
    seed: int = SEED,
    progress=None,
) -> list[dict[str, Any]]:
    """Every block scored window by window, with a dinucleotide-shuffled null on a sample of windows."""
    from genomeos.coords import Locus
    from genomeos.genome.genome import Genome
    from genomeos.genome.motifs import dinucleotide_shuffle, load_families, load_motifs

    say = progress or (lambda _m: None)
    motifs, families = load_motifs(relative=SITE_THRESHOLD), load_families()
    rng = random.Random(seed)
    by_chrom: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for b in blocks:
        by_chrom[b["chrom"]].append(b)
    out = []
    for chrom in sorted(by_chrom, key=lambda c: (len(c), c)):
        fa = reference / f"{chrom}.fa.gz"
        if not fa.exists():
            for b in by_chrom[chrom]:
                out.append({**b, "windows": 0, "dropped_windows": 0, "missing": "no reference sequence"})
            continue
        genome = Genome.from_fasta(fa)
        for b in by_chrom[chrom]:
            seq = str(genome.fetch(Locus(chrom, b["start"], b["end"])))
            kept, dropped = windows_of(seq)
            rec: dict[str, Any] = {**b, "windows": len(kept), "dropped_windows": dropped}
            if not kept:
                rec["missing"] = "no window with enough ACGT"
                out.append(rec)
                continue
            for cell, values in score_windows(kept, motifs, families, units, weights).items():
                rec[f"mean_{cell}"] = round(sum(values) / len(values), 5)
                rec[f"max_{cell}"] = round(max(values), 5)
            pick = kept if len(kept) <= shuffle_windows else rng.sample(kept, shuffle_windows)
            fake = [dinucleotide_shuffle(w, rng) for w in pick]
            for cell, values in score_windows(fake, motifs, families, units, weights).items():
                rec[f"shuffled_mean_{cell}"] = round(sum(values) / len(values), 5)
            out.append(rec)
        say(f"  {chrom}: {len(by_chrom[chrom])} blocks scored")
    return out


def blocks_cache_path(cache: Path = CACHE) -> Path:
    return cache / "block_scores.json"


def predicted_level(real: list[dict[str, Any]], control: list[dict[str, Any]], cells=CELLS) -> dict[str, Any]:
    """Secondary: the level the count model *predicts*, against a length-matched neutral tier and a null."""

    def col(rows: list[dict[str, Any]], key: str) -> list[float]:
        return [r[key] for r in rows if r.get(key) is not None]

    out: dict[str, Any] = {
        "median_length_kb": {
            "real_unknown": _round((_median([float(r["length"]) for r in real]) or 0) / 1000, 1),
            "matched_neutral": _round((_median([float(r["length"]) for r in control]) or 0) / 1000, 1),
        },
        "reading": "a prediction of the count model, not a measurement of activity",
    }
    for cell in cells:
        if not col(real, f"mean_{cell}"):
            continue
        out[cell] = {
            "real_unknown_vs_matched_neutral": mann_whitney(
                col(real, f"mean_{cell}"), col(control, f"mean_{cell}")
            ),
            "real_unknown_vs_its_shuffle": mann_whitney(
                col(real, f"mean_{cell}"), col(real, f"shuffled_mean_{cell}")
            ),
            "matched_neutral_vs_its_shuffle": mann_whitney(
                col(control, f"mean_{cell}"), col(control, f"shuffled_mean_{cell}")
            ),
        }
    return out


def sweep_agreement(rows: list[dict[str, Any]], cell: str = PRIMARY_CELL) -> dict[str, Any]:
    """Secondary, and not validation: the predicted score against the sweep's own AlphaGenome call."""
    usable = [r for r in rows if r.get(f"mean_{cell}") is not None]
    tested = [r for r in usable if r["tested_elements"]]
    out: dict[str, Any] = {
        "cell_line": cell,
        "reading": (
            "both sides are model outputs -- a JASPAR count model and an AlphaGenome deletion -- so this "
            "is agreement between two model readings, never one validating the other"
        ),
        "coverage": {
            "blocks": len(rows),
            "without_a_scorable_window": sum(1 for r in rows if r.get(f"mean_{cell}") is None),
            "with_a_tested_element": len(tested),
            "without_a_tested_element": len(usable) - len(tested),
            "with_a_lead": sum(1 for r in usable if r["moving_elements"]),
        },
    }
    for name, subset in (("tested_blocks_only", tested), ("all_blocks", usable)):
        if not subset:
            continue
        labels = [1 if r["moving_elements"] else 0 for r in subset]
        preds = {
            "count_model_mean": [r[f"mean_{cell}"] for r in subset],
            "count_model_max": [r[f"max_{cell}"] for r in subset],
            "block_length": [float(r["length"]) for r in subset],
            "tested_elements": [float(r["tested_elements"]) for r in subset],
        }
        rec = bootstrap_auc(labels, preds, (("count_model_mean", "block_length"),))
        rec["blocks"] = len(subset)
        rec["with_a_lead"] = sum(labels)
        rec["permutation_null"] = permutation_null(preds["count_model_mean"], labels)
        out[name] = rec
    return out


def blocks_verdict(measured: dict[str, Any], level: dict[str, Any], agreement: dict[str, Any]) -> dict:
    """The pre-registered block claim, decided on the measured readings, with the contradiction named."""
    lvl = (measured.get("level") or {}).get(PRIMARY_CELL) or {}
    z = lvl.get("z")
    above = z is not None and z > 2
    below = z is not None and z < -2
    transfer = (measured.get("transfer") or {}).get(PRIMARY_CELL) or {}
    diff = transfer.get("b-a") or {}
    transfers = bool(diff.get("ci95") and diff["ci95"][0] > 0)
    pred = (level.get(PRIMARY_CELL) or {}).get("real_unknown_vs_matched_neutral") or {}
    agree = (agreement.get("tested_blocks_only") or {}).get("count_model_mean") or {}
    return {
        "cell_line": PRIMARY_CELL,
        "measured_elements": (measured.get("coverage") or {})
        .get("real_unknown", {})
        .get("measured_elements"),
        "measured_level_auroc": lvl.get("auroc"),
        "measured_level_z": z,
        "measured_transfer_b_minus_a": diff.get("point"),
        "measured_transfer_ci95": diff.get("ci95"),
        "count_model_transfers_inside_the_real_unknown": transfers,
        "predicted_level_auroc": pred.get("auroc"),
        "predicted_level_z": pred.get("z"),
        "sweep_agreement_auroc": agree.get("auroc"),
        "sweep_agreement_ci95": agree.get("ci95"),
        "measured_statement": (
            "the measured reporter reads the real unknown as MORE active than the neutral tier, which "
            "contradicts the deletion sweep's per-element reading (0.248 against 0.293, z -5.7); the "
            "contradiction is the finding and is not split"
            if above
            else "the measured reporter reads the real unknown as LESS active than the neutral tier, "
            "agreeing in direction with the deletion sweep's per-element reading"
            if below
            else "the measured reporter reads the real unknown and the neutral tier at the same level: "
            "it neither confirms nor contradicts the deletion sweep's per-element reading"
        ),
        "transfer_statement": (
            "inside the real unknown, strict motif counts still predict measured activity over "
            "composition alone: §17's positive holds on this sequence"
            if transfers
            else "inside the real unknown, strict motif counts do not measurably beat composition alone; "
            "with this few measured elements that is as likely to be the sample size as the sequence"
        ),
    }


def run_blocks(
    model: dict[str, Any],
    tiers: dict[str, list[dict[str, Any]]],
    reference: Path = REFERENCE,
    cache: Path = CACHE,
    reuse: bool = True,
    progress=None,
) -> dict[str, Any]:
    """Reading 2: measured first, then the predicted level and last the agreement with the sweep."""
    say = progress or (lambda _m: None)
    units = model["family_units"]
    say(
        f"  {len(tiers['real_unknown'])} real-unknown blocks, {len(tiers['neutral'])} neutral; "
        f"{sum(1 for b in tiers['real_unknown'] if b['mpra_keys'])} hold a measured element"
    )
    measured = measured_reading(model, tiers)
    p = blocks_cache_path(cache)
    if reuse and p.exists():
        cached = json.loads(p.read_text())
        real, control = cached["real_unknown"], cached["matched_neutral"]
        say(f"  reusing {p}: {len(real)} real-unknown and {len(control)} matched neutral blocks scored")
    else:
        control_blocks = match_by_length(tiers["real_unknown"], tiers["neutral"])
        real = score_blocks(tiers["real_unknown"], units, model["weights"]["b"], reference, progress=say)
        control = score_blocks(control_blocks, units, model["weights"]["b"], reference, progress=say)
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(json.dumps({"real_unknown": real, "matched_neutral": control}))
    level = predicted_level(real, control)
    agreement = sweep_agreement(real)
    return {
        "preregistered": PREREGISTERED_BLOCKS,
        "min_log2": MIN_LOG2,
        "window": WINDOW,
        "shuffle_windows": SHUFFLE_WINDOWS,
        "blocks": {"real_unknown": len(real), "matched_neutral": len(control)},
        "windows": {
            "real_unknown": sum(r.get("windows", 0) for r in real),
            "matched_neutral": sum(r.get("windows", 0) for r in control),
            "dropped_for_n_bases": sum(r.get("dropped_windows", 0) for r in real + control),
            "blocks_without_a_window": sum(1 for r in real + control if not r.get("windows")),
        },
        "by_case": {
            case: {
                "blocks": len(rows),
                "with_a_measured_element": sum(1 for r in rows if r["mpra_keys"]),
                "with_a_tested_element": sum(1 for r in rows if r["tested_elements"]),
                "with_a_sweep_lead": sum(1 for r in rows if r["moving_elements"]),
                f"predicted_mean_{PRIMARY_CELL}": _round(
                    _mean([r.get(f"mean_{PRIMARY_CELL}") for r in rows])
                ),
            }
            for case in sorted({r["case"] for r in real if r["case"]})
            for rows in ([r for r in real if r["case"] == case],)
        },
        "measured": measured,
        "predicted_level": level,
        "agreement_with_the_sweep": agreement,
        "verdict": blocks_verdict(measured, level, agreement),
        "evidence": EVIDENCE,
    }


# ---------------------------------------------------------------- the run


def run(
    knowledge: Path = VISTA_KNOWLEDGE,
    mpra_knowledge: Path = MPRA_KNOWLEDGE,
    reference: Path = REFERENCE,
    cache: Path = CACHE,
    results_dir: Path = RESULTS_DIR,
    readings: tuple[str, ...] = ("vista", "blocks"),
    score_held_out: bool = True,
    reuse: bool = True,
    progress=None,
) -> dict[str, Any]:
    """Both readings: the transfer of §17's count positive to VISTA and to the real unknown."""
    say = progress or (lambda _m: None)
    out: dict[str, Any] = {
        "preregistered": {"vista": PREREGISTERED_VISTA, "real_unknown": PREREGISTERED_BLOCKS},
        "site_threshold": SITE_THRESHOLD,
        "evidence": EVIDENCE,
    }
    if "vista" in readings:
        say("reading 1: VISTA positives against negatives, in vivo")
        out["vista"] = run_vista(knowledge, cache, results_dir, VISTA_HELD_OUT, score_held_out, say)
    if "blocks" in readings:
        grammar = load_result("motif_grammar", results_dir)
        if not grammar:
            raise FileNotFoundError("no motif_grammar result: §17's model is what transfers")
        units = list(grammar["family_units"])
        say("placing the lentiMPRA library inside the tiers, to hold those elements out by location")
        tiers = tier_blocks(CHROMS, results_dir, ELEMENT_TABLES, mpra_positions(mpra_knowledge))
        exclude = frozenset(k for rows in tiers.values() for b in rows for k in b["mpra_keys"])
        say(f"  {len(exclude)} lentiMPRA elements lie inside a real-unknown or neutral block")
        if not score_held_out:
            out["real_unknown"] = {"preregistered": PREREGISTERED_BLOCKS, "coverage": block_coverage(tiers)}
            return out
        model = mpra_model(units, mpra_knowledge, exclude, say)
        say("reading 2: the 882 blocks of the real unknown, measured first")
        out["mpra_model"] = {
            k: v for k, v in model.items() if k not in ("rows", "activity", "meta", "weights")
        }
        out["real_unknown"] = run_blocks(model, tiers, reference, cache, reuse, say)
    return out


def run_and_save(
    knowledge: Path = VISTA_KNOWLEDGE,
    mpra_knowledge: Path = MPRA_KNOWLEDGE,
    reference: Path = REFERENCE,
    cache: Path = CACHE,
    results_dir: Path = RESULTS_DIR,
    readings: tuple[str, ...] = ("vista", "blocks"),
    score_held_out: bool = True,
    reuse: bool = True,
    name: str | None = None,
    progress=None,
) -> dict[str, Any]:
    out = run(
        knowledge, mpra_knowledge, reference, cache, results_dir, readings, score_held_out, reuse, progress
    )
    save_result(
        name or ("motif_transfer" if score_held_out else "motif_transfer_preregistration"),
        out,
        results_dir,
    )
    return out


__all__ = [
    "PREREGISTERED_BLOCKS",
    "PREREGISTERED_VISTA",
    "VISTA_HELD_OUT",
    "auprc",
    "auroc",
    "bootstrap_auc",
    "build_vista_counts",
    "chunk_sequence",
    "count_row",
    "covariates",
    "match_by_length",
    "measured_reading",
    "merge_sites",
    "model_row",
    "mpra_model",
    "permutation_null",
    "predicted_level",
    "run",
    "run_and_save",
    "run_blocks",
    "run_vista",
    "scan_counts",
    "sweep_agreement",
    "unit_counts",
    "windows_of",
]
