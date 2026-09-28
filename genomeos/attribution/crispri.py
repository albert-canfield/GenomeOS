# SPDX-License-Identifier: AGPL-3.0-or-later
"""The enhancer-to-gene attribution held against measured perturbations: CRISPRi screens.

Every enhancer-to-gene target GenomeOS names came from a predicted deletion (AlphaGenome) or a
heuristic (the nearest TSS inside the CTCF node). Neither had been held against a measured
perturbation. The ENCODE enhancer-gene benchmark (EngreitzLab/CRISPR_comparison, Gschwind et al.
2025) collects CRISPRi screens in which an element was silenced and every nearby gene's expression
measured: 10,356 valid pairs in K562 (Nasser 2021 from Ulirsch 2016, Gasperini 2019, Schraivogel
2020) for training, and 4,378 held-out pairs in five cell types. A pair is positive when the
element's silencing lowers the gene significantly.

The two tables are streamed once into data/knowledge/crispri (1.4 MB, local, never committed).
Each pair is joined to the all-enhancer deletion table (data/knowledge/alphagenome/all_elements)
by overlap, and three predictors are compared on the same pairs:

- distance: the inverse distance from the element to the gene's TSS;
- activity over distance: sqrt(DNase x H3K27ac) measured in the screen's own cell, divided by
  distance (the activity-by-contact idea with a power-law contact of exponent 1; no Hi-C);
- the deletion: the predicted expression drop for the pair's own gene in the screen's cell line.

Until 2026-09-22 that drop was read from the compact all_elements table, which keeps one gene per
element, so a pair whose measured gene is not that gene scored zero — a zero that said "the table
had nothing to say" and was consumed as "the model predicts no effect". The same sweep wrote a
per-element response cache (`ElementCache`, data/knowledge/alphagenome/elements) with threshold 0.0,
carrying every gene in the scorer's 1 Mb window with a signed log2 fold change on each of K562,
HepG2, GM12878 and IMR-90's own track, so the drop is now read from there at no model request.
`top_target` still means what it always meant and stays a separate feature; `deletion_census`
counts the pairs the compact table could not answer and how many of them the cache answers.

The pre-registered test (fixed in PREREGISTERED before the held-out pairs were scored): a logistic
model fitted on the K562 training pairs with activity, distance and the deletion beats the same
model without the deletion in AUPRC on the held-out K562 pairs and on the held-out GM12878 pairs.
The other held-out cell types have no AlphaGenome line in the table and are refused, with that
reason. The screens are `experimental` evidence; the deletion stays `predicted`.

The second question of the module is the contact term itself. 1/distance is a power law standing in
for a measurement that exists: `score_contact` replaces it with the measured contact between the
element's bin and the TSS's bin, read from a released 4D Nucleome in-situ Hi-C matrix of the screen's
own cell line at 5 kb (genome/hic_contact.py, by range request, nothing downloaded). Its
pre-registration is PREREGISTERED_CONTACT, fixed before the held-out pairs were scored: activity x
measured contact beats activity over distance on held-out K562. Both contact features are
`experimental`; the deletion stays `predicted`.
"""

from __future__ import annotations

import bisect
import csv
import gzip
import json
import math
import random
import urllib.request
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from genomeos.genome.epigenome import solve

BASE_URL = "https://raw.githubusercontent.com/EngreitzLab/CRISPR_comparison/main/resources/crispr_data/"
TRAINING = "EPCrisprBenchmark_combined_data.training_K562.GRCh38.tsv.gz"
HELDOUT = "EPCrisprBenchmark_combined_data.heldout_5_cell_types.GRCh38.tsv.gz"
KNOWLEDGE = Path("data/knowledge/crispri")
ELEMENTS = Path("data/knowledge/alphagenome/all_elements")
ELEMENT_CACHE = Path("data/knowledge/alphagenome/elements")
EVIDENCE = (
    "experimental: CRISPRi enhancer-gene screens, ENCODE benchmark (EngreitzLab/CRISPR_comparison, "
    "Gschwind et al. 2025); predicted: AlphaGenome deletion per element"
)
MODEL_CELLS = ("K562", "HepG2", "GM12878", "IMR-90")  # the lines the deletion table was scored in
MIN_DISTANCE = 1_000  # a TSS closer than this counts as 1 kb, so contact stays finite
REACH = 10_000  # the longest element in the deletion table, for the overlap search
BOOTSTRAPS = 200
FEATURES = {
    "distance": ("log_distance",),
    "activity + distance": ("log_distance", "log_activity", "activity_over_distance"),
    "activity + distance + deletion": (
        "log_distance",
        "log_activity",
        "activity_over_distance",
        "top_target",
        "deletion_drop",
    ),
}
PREREGISTERED = (
    "fitted on the K562 training pairs, 'activity + distance + deletion' has a higher AUPRC than "
    "'activity + distance' on the held-out K562 pairs and on the held-out GM12878 pairs"
)

# --- measured contact instead of 1/distance (the same test, a different contact term) --------------
CONTACT_EVIDENCE = (
    "experimental: 4D Nucleome in-situ Hi-C contact matrix at 5 kb, the file's own balancing vector, "
    "read by range request; observed over expected from the file's expected-value vector"
)
CONTACT_CELLS = ("K562", "GM12878")  # the cell lines with a matrix; the rest are refused, with reason
CONTACT_FEATURES = {
    "distance": ("log_distance",),
    "contact": ("log_contact",),
    "activity + distance": ("log_distance", "log_activity", "activity_over_distance"),
    "activity x contact": ("log_activity", "log_contact", "activity_x_contact"),
    "activity + distance + deletion": FEATURES["activity + distance + deletion"],
    "activity + contact + deletion": (
        "log_activity",
        "log_contact",
        "activity_x_contact",
        "top_target",
        "deletion_drop",
    ),
    "activity + distance + contact": (
        "log_distance",
        "log_activity",
        "activity_over_distance",
        "log_contact",
        "log_oe",
    ),
}
PREREGISTERED_CONTACT = (
    "fitted on the K562 training pairs that have a measured contact, 'activity x contact' (measured "
    "4DN Hi-C contact at 5 kb in the screen's own cell line, balanced) has a higher AUPRC than "
    "'activity + distance' (the same activity over 1/distance) on the held-out K562 pairs"
)

# --- the published baseline on the same pairs, and the coverage-matched arm (registered 2026-09-27) --
#: Figures copied from the primary source, never from memory: Gschwind et al., "An encyclopedia of
#: human enhancer-gene regulatory interactions", Nature 2026, doi:10.1038/s41586-026-10781-4
#: (PMC13471189), Supplementary Table 3, sheets "Combined K562 data (training)" and "Held-out
#: benchmarks". The training pair set is this module's TRAINING file (10,356 pairs, 471 positives,
#: unweighted AUPRC, hold-one-chromosome-out for supervised models, a pair with no prediction given
#: the minimum score). The held-out pair set is this module's HELDOUT file (4,378 pairs, 190
#: positives, five cell types pooled), and its AUPRC is weighted: every pair by its column
#: `direct_vs_indirect_negative` (the 190 positives sum to 157.39, the paper's "157.39 weighted").
#: The AUPRC is the benchmark's own estimator (yardstick pr_curve, first and last points dropped,
#: trapezoid; EngreitzLab/CRISPR_comparison workflow/scripts/crisprComparisonBootstrapFunctions.R).
PUBLISHED_SOURCE = (
    "Gschwind et al. 2026, Nature, doi:10.1038/s41586-026-10781-4, Supplementary Table 3 "
    "(https://pmc.ncbi.nlm.nih.gov/articles/PMC13471189/)"
)
PUBLISHED = {
    "training": {  # AUPRC [95% bootstrap interval], unweighted, 10,356 pairs / 471 positives
        "ENCODE-rE2G_Extended": (0.7373, 0.6931, 0.7750),
        "ENCODE-rE2G": (0.6622, 0.6156, 0.7059),
        "ABC_A=DNase x H3K27ac, C=ENCODE Hi-C": (0.6128, 0.5580, 0.6570),
        "ABC_A=DNase, C=Average ENCODE Hi-C": (0.5649, 0.5109, 0.6103),
        "In element (DHS) & DNase RPM x distance (norm.)": (0.5272, 0.4786, 0.5745),
        "Distance to TSS": (0.4359, 0.3873, 0.4824),
    },
    "heldout": {  # weighted AUPRC [95% bootstrap interval], 4,378 pairs / 190 positives, pooled
        "ENCODE-rE2G": (0.5562, 0.4679, 0.6312),
        "ABC_A=DNase, C=Average ENCODE Hi-C": (0.4654, 0.3782, 0.5414),
        "EPIraction": (0.3827, 0.3025, 0.4558),
        "Distance to TSS": (0.3631, 0.2812, 0.4379),
    },
}
#: The estimator was checked against the source BEFORE this registration, on the one predictor that
#: needs no model: raw distance to TSS gives 0.4359 unweighted on TRAINING (published 0.4359) and
#: 0.3631 weighted on HELDOUT (published 0.3631). No model of this module had been scored with it.
PREREGISTERED_PUBLISHED: dict[str, str] = {
    "written": (
        "2026-09-27, before any model of this module was scored on the full pair sets, weighted, or "
        "on the coverage arms below. The model set is FEATURES, frozen; its weights are the ones "
        "score() fits on the covered K562 training pairs; nothing is refitted"
    ),
    "second_cell_type": (
        "AlphaGenome carries its own RNA-seq tracks for HCT116 (EFO:0002824, 3 tracks), Jurkat "
        "(CLO:0007045, 1) and WTC11 (EFO:0009747, 3), read from the model's output_metadata on "
        "2026-09-27; no proxy is needed. The refusal is a fact about the sweep's cache, which kept "
        "four cell lines' values per gene. Only HCT116 can carry a result: 396 held-out pairs, 34 "
        "regulated, 363 on a deleted element (Jurkat 6 regulated covered, WTC11 14, GM12878 14). "
        "Scoring it with the frozen feature needs one deletion request per registry element "
        "overlapping a covered HCT116 pair: 705 (363 if the CRISPR element itself were deleted, "
        "which changes the feature and is not the frozen model). The lane's budget is 20, so HCT116 "
        "is registered here and NOT run. When bought: the same code path with HCT116 added to the "
        "cells kept per gene, weights unchanged; passes if the gain on the 363 covered HCT116 pairs "
        "is above zero; stated as replicated only if its chromosome-bootstrap interval excludes zero. "
        "A gain at or below zero means the result is a K562 result and must be described as one"
    ),
    "published_training": (
        "hold-one-chromosome-out on all 10,356 TRAINING pairs (fitted on the covered pairs of the "
        "other chromosomes, the frozen recipe; an uncovered pair carries deletion features of zero, "
        "the benchmark's minimum-score rule), unweighted AUPRC by the benchmark's estimator, set "
        "beside PUBLISHED['training']. Comparable pair set and split; NOT comparable in one respect "
        "that is stated with the number: this module's deletion features were chosen after reading "
        "the single predictors on these same pairs (ENCODE-rE2G's features were also chosen by "
        "sequential selection on this set, its Fig. 4c, so both carry a selection on the data)"
    ),
    "published_heldout": (
        "the frozen weights applied to all 4,378 HELDOUT pairs, pooled, weighted AUPRC by the "
        "benchmark's estimator, set beside PUBLISHED['heldout']. Expected handicap stated before "
        "the number: 2,392 of the 4,378 pairs (HCT116, Jurkat, WTC11) carry no deletion value, so "
        "the deletion can only act on the K562 and GM12878 pairs. Bands against ENCODE-rE2G's "
        "published interval [0.468, 0.631]: above 0.631 'higher than the published ENCODE-rE2G "
        "figure on the same pairs'; inside it 'in the range of ENCODE-rE2G'; below 0.468 'below "
        "ENCODE-rE2G: the gain is over a baseline weaker than the published state of the art'. "
        "'activity + distance' is placed against ABC's [0.378, 0.541] by the same bands, which says "
        "whether the project's own baseline is weak. The K562-only weighted figure is reported as "
        "this module's alone: the source gives no per-cell-type held-out figure"
    ),
    "coverage_matched": (
        "held-out K562. Arm 1, all 1,918 pairs, uncovered ones with deletion features of zero: the "
        "gain with its chromosome-bootstrap interval. Arm 2, the covered pairs with the regulated "
        "arm thinned at random to the non-regulated arm's coverage rate, 1,000 draws: the median "
        "gain and the share of draws above zero. Arm 3, the control: 'activity + distance' plus a "
        "covered-or-not indicator in place of the deletion features, fitted on all training pairs "
        "and scored on all 1,918. The gain is 'invariant to coverage' only if arm 1's interval "
        "excludes zero, arm 2's share is at least 0.95, and arm 3's gain is below arm 1's; otherwise "
        "the headline must say 'on the pairs the sweep covered'"
    ),
}


@dataclass
class Pair:
    chrom: str
    start: int
    end: int
    gene: str
    cell: str
    dataset: str
    distance: float
    dhs: float
    h3k27ac: float
    regulated: bool
    tss: int | None = None
    features: dict[str, float] = field(default_factory=dict)
    covered: bool = False
    contact: dict[str, Any] | None = None
    weight: float = 1.0  # the benchmark's probability that the pair's effect is direct
    category: str = ""  # the benchmark's elementChromatinCategory of the tested element
    # the benchmark's five power columns (PowerAtEffectSize10..50), as measured.CrispriPair carries
    # them: None where a table lacks the column or leaves it empty. Carried, read by no scorer.
    power_at_effect_size_10: float | None = None
    power_at_effect_size_15: float | None = None
    power_at_effect_size_20: float | None = None
    power_at_effect_size_25: float | None = None
    power_at_effect_size_50: float | None = None

    @property
    def element(self) -> tuple[str, str, int, int]:
        return (self.cell, self.chrom, self.start, self.end)

    @property
    def midpoint(self) -> int:
        return (self.start + self.end) // 2


def fetch(name: str, knowledge: Path = KNOWLEDGE) -> Path:
    """Stream one benchmark table into the cache unless it is there already."""
    knowledge.mkdir(parents=True, exist_ok=True)
    p = knowledge / name
    if not p.exists():
        with urllib.request.urlopen(BASE_URL + name, timeout=120) as r:
            p.write_bytes(r.read())
    return p


def _weight(value: str | None) -> float:
    """The benchmark's direct-effect weight of a pair; 1 where the column is absent or empty."""
    try:
        return float(value) if value not in (None, "", "NA") else 1.0
    except ValueError:
        return 1.0


def _power(value: str | None) -> float | None:
    """A power column's value; None for an absent, empty or NA cell, never a made-up zero."""
    if value in (None, "", "NA"):
        return None
    try:
        return float(value)
    except ValueError:
        return None


def parse(lines: Any) -> list[Pair]:
    """Valid pairs of one benchmark table; pairs flagged as promoter or exon overlaps are dropped."""
    out = []
    for r in csv.DictReader(lines, delimiter="\t"):
        if r["ValidConnection"] != "TRUE":
            continue
        out.append(
            Pair(
                chrom=r["chrom"],
                start=int(r["chromStart"]),
                end=int(r["chromEnd"]),
                gene=r["measuredGeneSymbol"],
                cell=r["CellType"],
                dataset=r["Dataset"],
                distance=float(r["distanceToTSS"]),
                dhs=float(r["DHS.RPM"] or 0),
                h3k27ac=float(r["H3K27ac.RPM"] or 0),
                regulated=r["Regulated"] == "TRUE",
                tss=int(r["startTSS"]) if r.get("startTSS") else None,
                weight=_weight(r.get("direct_vs_indirect_negative")),
                category=r.get("elementChromatinCategory") or "",
                power_at_effect_size_10=_power(r.get("PowerAtEffectSize10")),
                power_at_effect_size_15=_power(r.get("PowerAtEffectSize15")),
                power_at_effect_size_20=_power(r.get("PowerAtEffectSize20")),
                power_at_effect_size_25=_power(r.get("PowerAtEffectSize25")),
                power_at_effect_size_50=_power(r.get("PowerAtEffectSize50")),
            )
        )
    return out


def load(name: str, knowledge: Path = KNOWLEDGE) -> list[Pair]:
    with gzip.open(knowledge / name, "rt") as f:
        return parse(f)


class DeletionTable:
    """The all-enhancer deletion table, one chromosome at a time, searched by overlap."""

    def __init__(self, root: Path = ELEMENTS) -> None:
        self.root = root
        self._by_chrom: dict[str, tuple[list[dict[str, Any]], list[int]]] = {}

    def overlapping(self, chrom: str, start: int, end: int) -> list[dict[str, Any]]:
        if chrom not in self._by_chrom:
            p = self.root / f"{chrom}.json"
            els = sorted(json.loads(p.read_text()), key=lambda e: e["start"]) if p.exists() else []
            self._by_chrom[chrom] = (els, [e["start"] for e in els])
        els, starts = self._by_chrom[chrom]
        out = []
        j = bisect.bisect_left(starts, end) - 1
        while j >= 0 and els[j]["start"] > start - REACH:
            if els[j]["end"] > start:
                out.append(els[j])
            j -= 1
        return out


class ElementCache:
    """The sweep's per-element response cache: every gene it scored, not only the top target.

    `all_elements` is a compact derived table that keeps one gene per element, so a pair whose
    measured gene is not that gene scored a drop of zero — a zero that said "the table had nothing
    to say", not "no effect is predicted". The same run wrote this cache with `threshold=0.0`
    (scripts/enhancer_targets_all.py, `worker_scorer`), so every gene in the scorer's 1 Mb window
    carries a signed log2 fold change on each of K562, HepG2, GM12878 and IMR-90's own track,
    uncensored. One chromosome's archive is held at a time, so `annotate` walks chromosome by
    chromosome; elements the archive does not carry are looked up as loose per-element files.
    """

    def __init__(self, root: Path = ELEMENT_CACHE) -> None:
        self.root = root
        self._chrom: str | None = None
        self._archive: dict[str, Any] = {}

    def _load(self, chrom: str) -> None:
        if self._chrom == chrom:
            return
        p = self.root / f"{chrom}.json.gz"
        archive: dict[str, Any] = {}
        if p.exists():
            with gzip.open(p, "rt") as fh:
                archive = json.load(fh)
        self._chrom, self._archive = chrom, archive

    def value(self, chrom: str, element_id: str, gene: str, cell: str) -> float | None:
        """This gene's signed predicted log2 fold change on this cell's own track, or None."""
        self._load(chrom)
        hit = self._archive.get(element_id)
        if hit is None:
            p = self.root / chrom / f"{element_id}.json"
            if p.exists():
                hit = json.loads(p.read_text())
        for g in (hit or {}).get("genes") or []:
            if g["gene"] == gene:
                v = (g.get("by_cell") or {}).get(cell)
                if v is not None:
                    return float(v)
        return None


def deletion_values(
    elements: list[dict[str, Any]],
    gene: str,
    cell: str,
    cache: ElementCache | None = None,
    chrom: str = "",
) -> tuple[float, list[float]]:
    """(1 if the gene is an overlapping element's top predicted target, every signed change found).

    `top` is what the compact table records: the gene IS the element's single top predicted target.
    The values are the sweep's own per-element responses for this gene on this cell's own track when
    a `cache` is given, and the compact table's single entry otherwise — so with a cache a pair whose
    gene is not the top target carries the change the sweep predicted for it, and without one the
    values are exactly the old ones. An empty list means nothing was scored for this gene in this
    cell, which is a different statement from a scored change of zero.
    """
    top, values = 0.0, []
    for e in elements:
        compact: list[float] = []
        for key, by_cell in (
            ("predicted_coding", "predicted_coding_by_cell"),
            ("predicted", "predicted_by_cell"),
        ):
            p = e.get(key)
            if not p or p["gene"] != gene:
                continue
            top = 1.0
            v = (e.get(by_cell) or {}).get(cell)
            if v is not None:
                compact.append(float(v))
        cached = cache.value(chrom, e["id"], gene, cell) if cache is not None else None
        values.extend([cached] if cached is not None else compact)
    return top, values


def deletion_drop(values: list[float]) -> float:
    """The predicted drop: the largest fall, floored at zero.

    The deletion's log2 fold change is negative for an activating element, so the drop is its
    negation, floored at zero: a predicted rise is not evidence for the activation the screens call.
    """
    return max([0.0, *(-v for v in values)])


def deletion_for(
    elements: list[dict[str, Any]],
    gene: str,
    cell: str,
    cache: ElementCache | None = None,
    chrom: str = "",
) -> tuple[float, float]:
    """(1 if the gene is an overlapping element's top predicted target, the predicted drop in `cell`)."""
    top, values = deletion_values(elements, gene, cell, cache, chrom)
    return top, deletion_drop(values)


def annotate(pairs: list[Pair], table: DeletionTable, cache: ElementCache | None = None) -> None:
    """Fill each pair's features; a pair counts as covered when a deleted element overlaps it.

    The pairs are walked chromosome by chromosome so that the per-element cache, which holds one
    chromosome's archive at a time, is read once per chromosome rather than once per pair.
    """
    by_chrom: dict[str, list[Pair]] = defaultdict(list)
    for p in pairs:
        by_chrom[p.chrom].append(p)
    for chrom in sorted(by_chrom):
        for p in by_chrom[chrom]:
            _annotate_one(p, table, cache)


def _annotate_one(p: Pair, table: DeletionTable, cache: ElementCache | None) -> None:
    els = table.overlapping(p.chrom, p.start, p.end)
    p.covered = bool(els)
    d = max(MIN_DISTANCE, p.distance)
    activity = math.sqrt(max(p.dhs, 0.0) * max(p.h3k27ac, 0.0))
    top, values = deletion_values(els, p.gene, p.cell, cache, p.chrom) if p.cell in MODEL_CELLS else (0.0, [])
    p.features = {
        "node_nearest": float(any((e.get("inferred") or {}).get("gene") == p.gene for e in els)),
        "log_distance": math.log(d),
        "log_activity": math.log1p(activity),
        "activity_over_distance": math.log1p(activity) - math.log(d),
        "top_target": top,
        "deletion_drop": deletion_drop(values),
        # not a model column: whether the sweep scored this gene in this cell at all, so that a
        # drop of zero can be told apart from a table that had nothing to say about the pair
        "deletion_answered": float(bool(values)),
    }


def deletion_census(pairs: list[Pair]) -> dict[str, Any]:
    """How many covered pairs the compact table could not answer, and how many carry a value now.

    The named defect this counts: `deletion_drop` was 0 for every pair whose measured gene is not
    the element's single top predicted target, whatever the sweep predicted for that gene. The zero
    meant "the table had nothing to say" and was consumed as "the model predicts no effect".
    `structural_zero` is that population; `answered` is how many of them the per-element cache
    scored, split into the ones that carry a fall (a non-zero `deletion_drop` now) and the ones the
    sweep predicted would rise, which still read zero but now mean it.
    """
    out: dict[str, Any] = {}
    for cell in sorted({p.cell for p in pairs if p.covered}):
        rows = [p for p in pairs if p.covered and p.cell == cell]
        zeros = [p for p in rows if not p.features["top_target"]]
        answered = [p for p in zeros if p.features.get("deletion_answered")]
        fell = [p for p in answered if p.features["deletion_drop"] > 0]
        out[cell] = {
            "covered": len(rows),
            "regulated": sum(p.regulated for p in rows),
            "gene_is_the_top_target": len(rows) - len(zeros),
            "structural_zero": len(zeros),
            "structural_zero_regulated": sum(p.regulated for p in zeros),
            "structural_zero_answered": len(answered),
            "structural_zero_answered_regulated": sum(p.regulated for p in answered),
            "structural_zero_answered_as_a_fall": len(fell),
            "structural_zero_answered_as_a_rise_or_flat": len(answered) - len(fell),
            "share_answered": round(len(answered) / len(zeros), 4) if zeros else None,
        }
    return out


def average_precision(scores: list[float], labels: list[bool]) -> float | None:
    """Area under the precision-recall curve as the mean precision at each positive (ties by order)."""
    npos = sum(labels)
    if not npos:
        return None
    hits, total = 0, 0.0
    for k, (_, y) in enumerate(sorted(zip(scores, labels, strict=True), key=lambda t: -t[0]), 1):
        if y:
            hits += 1
            total += hits / k
    return total / npos


def auroc(scores: list[float], labels: list[bool]) -> float | None:
    """Mann-Whitney area under the ROC curve, ties counted half."""
    order = sorted(range(len(scores)), key=lambda i: scores[i])
    ranks = [0.0] * len(scores)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and scores[order[j + 1]] == scores[order[i]]:
            j += 1
        for k in range(i, j + 1):
            ranks[order[k]] = (i + j) / 2 + 1
        i = j + 1
    npos = sum(labels)
    nneg = len(labels) - npos
    if not npos or not nneg:
        return None
    return (sum(r for r, y in zip(ranks, labels, strict=True) if y) - npos * (npos + 1) / 2) / (npos * nneg)


def logistic_fit(x: list[list[float]], y: list[bool], lam: float = 1e-3, rounds: int = 25) -> list[float]:
    """Logistic regression by Newton steps, intercept first, a small ridge on the other weights."""
    k = len(x[0]) + 1
    rows = [[1.0, *r] for r in x]
    w = [0.0] * k
    for _ in range(rounds):
        grad = [0.0] * k
        hess = [[0.0] * k for _ in range(k)]
        for r, yi in zip(rows, y, strict=True):
            z = sum(a * b for a, b in zip(w, r, strict=True))
            p = 1 / (1 + math.exp(-max(-30.0, min(30.0, z))))
            g, h = (1.0 if yi else 0.0) - p, p * (1 - p)
            for i in range(k):
                grad[i] += g * r[i]
                hi = hess[i]
                for j in range(i, k):
                    hi[j] += h * r[i] * r[j]
        for i in range(k):
            for j in range(i):
                hess[i][j] = hess[j][i]
            if i:
                hess[i][i] += lam
                grad[i] -= lam * w[i]
        step = solve(hess, grad)
        step = damped_step(rows, y, w, step, lam)  # an undamped step leaves for 1e12 on a near-separable fit
        w = [a + b for a, b in zip(w, step, strict=True)]
        if max(abs(s) for s in step) < 1e-6:
            break
    return w


def penalised_objective(rows: list[list[float]], y: list[bool], w: list[float], lam: float) -> float:
    """The log-likelihood the fit maximises, with the same ridge the Hessian carries; the intercept is spared.

    `rows` already carry the leading 1.0, so `w[0]` is the intercept and is not penalised. The linear
    predictor is clamped the way the sigmoid in `logistic_fit` clamps it, so the two agree on a
    saturated row instead of disagreeing about a point neither can represent.
    """
    total = 0.0
    for r, yi in zip(rows, y, strict=True):
        z = max(-30.0, min(30.0, sum(a * b for a, b in zip(w, r, strict=True))))
        total += (z if yi else 0.0) - math.log1p(math.exp(z))
    return total - 0.5 * lam * sum(v * v for v in w[1:])


def damped_step(
    rows: list[list[float]], y: list[bool], w: list[float], step: list[float], lam: float
) -> list[float]:
    """The Newton step halved until it stops lowering the objective, or zeroed if it never does.

    Newton's method is not globally convergent for logistic regression. Where the classes nearly
    separate, the full step can saturate every linear predictor; a saturated row whose label is on
    the wrong side then offers a curvature of about 9.4e-14 against a gradient of order one, and the
    next solve divides one by the other and sends the weights to 1e12. The iterate does not return,
    the step never falls under the convergence threshold, and the fit ends at a point whose
    predictions are all 0 while still ordering the rows correctly, which is why nothing downstream
    that reads a ranking could ever see it.

    The full step is accepted whenever it does not lower the objective, which is every iteration of a
    run that converges, so this reproduces the undamped fit exactly where the undamped fit was right
    rather than merely closely. A zero step ends the iteration through the caller's own check.
    """
    base = penalised_objective(rows, y, w, lam)
    scale = 1.0
    for _ in range(20):
        trial = [a + scale * b for a, b in zip(w, step, strict=True)]
        if penalised_objective(rows, y, trial, lam) >= base:
            return [scale * s for s in step]
        scale /= 2
    return [0.0] * len(step)


def logistic_score(w: list[float], x: list[list[float]]) -> list[float]:
    return [w[0] + sum(a * b for a, b in zip(w[1:], r, strict=True)) for r in x]


def matrix(pairs: list[Pair], names: tuple[str, ...]) -> list[list[float]]:
    return [[p.features[n] for n in names] for p in pairs]


def metrics(scores: list[float], labels: list[bool]) -> dict[str, Any]:
    ap, au = average_precision(scores, labels), auroc(scores, labels)
    return {
        "pairs": len(labels),
        "positives": sum(labels),
        "baseline_precision": round(sum(labels) / len(labels), 4) if labels else None,
        "auprc": None if ap is None else round(ap, 4),
        "auroc": None if au is None else round(au, 4),
    }


def gain_interval(
    a: list[float], b: list[float], pairs: list[Pair], seed: int = 0, n: int = BOOTSTRAPS
) -> dict[str, Any]:
    """AUPRC of `a` minus `b`, with a 95% interval from resampling whole chromosomes."""
    by_chrom: dict[str, list[int]] = defaultdict(list)
    for i, p in enumerate(pairs):
        by_chrom[p.chrom].append(i)
    chroms = sorted(by_chrom)
    labels = [p.regulated for p in pairs]
    point = (average_precision(a, labels) or 0) - (average_precision(b, labels) or 0)
    rng = random.Random(seed)
    diffs = []
    for _ in range(n):
        idx = [i for c in (rng.choice(chroms) for _ in chroms) for i in by_chrom[c]]
        lab = [labels[i] for i in idx]
        if not any(lab):
            continue
        diffs.append(
            (average_precision([a[i] for i in idx], lab) or 0)
            - (average_precision([b[i] for i in idx], lab) or 0)
        )
    diffs.sort()
    if not diffs:
        return {"gain": round(point, 4), "ci95": None, "resamples": 0}
    lo, hi = diffs[int(0.025 * len(diffs))], diffs[min(len(diffs) - 1, int(0.975 * len(diffs)))]
    return {"gain": round(point, 4), "ci95": [round(lo, 4), round(hi, 4)], "resamples": len(diffs)}


def one_call_per_element(pairs: list[Pair], threshold: float) -> dict[str, Any]:
    """Each method names one tested gene per element; how often is that gene the regulated one.

    The deletion calls only when its top target is among the tested genes and the predicted drop in
    the screen's cell exceeds `threshold`; distance and activity over distance always call. The
    comparison on the deletion's own elements says whether it chooses the gene better or only
    chooses easier elements.
    """
    groups: dict[tuple[str, str, int, int], list[Pair]] = defaultdict(list)
    for p in pairs:
        groups[p.element].append(p)

    def closest(g: list[Pair]) -> Pair:
        return min(g, key=lambda p: p.distance)

    def by_activity(g: list[Pair]) -> Pair:
        return max(g, key=lambda p: p.features["activity_over_distance"])

    def by_deletion(g: list[Pair]) -> Pair | None:
        best = max(g, key=lambda p: (p.features["top_target"], p.features["deletion_drop"]))
        ok = best.features["top_target"] and best.features["deletion_drop"] > threshold
        return best if ok else None

    def tally(picker: Any, keys: list[tuple[str, str, int, int]]) -> dict[str, Any]:
        calls = [c for c in (picker(groups[k]) for k in keys) if c is not None]
        right = sum(c.regulated for c in calls)
        return {
            "calls": len(calls),
            "right": right,
            "precision": round(right / len(calls), 4) if calls else None,
        }

    def by_node(g: list[Pair]) -> Pair | None:
        return next((p for p in g if p.features["node_nearest"]), None)

    keys = sorted(groups)
    called = [k for k in keys if by_deletion(groups[k]) is not None]
    return {
        "threshold": threshold,
        "elements": len(keys),
        "elements_with_a_regulated_gene": sum(any(p.regulated for p in groups[k]) for k in keys),
        "all_elements": {
            "closest gene": tally(closest, keys),
            "activity over distance": tally(by_activity, keys),
            "nearest TSS in node": tally(by_node, keys),
            "deletion": tally(by_deletion, keys),
        },
        "deletion_elements": {
            "closest gene": tally(closest, called),
            "activity over distance": tally(by_activity, called),
            "deletion": tally(by_deletion, called),
        },
    }


def element_level(pairs: list[Pair]) -> dict[str, Any]:
    """Does the element regulate any tested gene: measured activity, distance and the deletion per element."""
    groups: dict[tuple[str, str, int, int], list[Pair]] = defaultdict(list)
    for p in pairs:
        groups[p.element].append(p)
    labels, cols = [], defaultdict(list)
    for g in groups.values():
        labels.append(any(p.regulated for p in g))
        f = g[0].features
        cols["activity"].append(f["log_activity"])
        cols["closest distance"].append(-min(p.features["log_distance"] for p in g))
        cols["deletion drop"].append(max(p.features["deletion_drop"] for p in g))
    return {name: metrics(v, labels) for name, v in cols.items()}


def coverage_by_arm(pairs: list[Pair]) -> dict[str, Any]:
    """Coverage of each arm, and the two baselines on the whole set beside the covered subset.

    The sweep deleted registry elements lying inside a CTCF node, so coverage is a property of the
    element and selects on it: every pair of a covered element is covered (partly_covered is the
    proof). An arm covered better than the other would make any rate read on the covered subset a
    statement about the selection instead of about the predictors, so both arms and both baselines
    are reported. The comparison between predictors stays fair whatever the coverage, because they
    are scored on the same pairs; what coverage governs is what the subset is a sample of.
    """
    groups: dict[tuple[str, str, int, int], list[Pair]] = defaultdict(list)
    for p in pairs:
        groups[p.element].append(p)
    arms = {}
    for arm, flag in (("regulated", True), ("not regulated", False)):
        rows = [p for p in pairs if p.regulated is flag]
        covered = sum(p.covered for p in rows)
        arms[arm] = {
            "pairs": len(rows),
            "covered": covered,
            "fraction": round(covered / len(rows), 4) if rows else None,
        }
    covered_pairs = [p for p in pairs if p.covered]
    baselines = {}
    for name, rows in (("all valid pairs", pairs), ("on a deleted element", covered_pairs)):
        lab = [p.regulated for p in rows]
        baselines[name] = {
            "distance": metrics([-p.features["log_distance"] for p in rows], lab),
            "activity over distance": metrics([p.features["activity_over_distance"] for p in rows], lab),
        }
    return {
        "arms": arms,
        "elements": len(groups),
        "elements_fully_covered": sum(all(p.covered for p in g) for g in groups.values()),
        "elements_partly_covered": sum(
            any(p.covered for p in g) and not all(p.covered for p in g) for g in groups.values()
        ),
        "baselines": baselines,
    }


def score(
    training: list[Pair],
    heldout: list[Pair],
    table: DeletionTable,
    cache: ElementCache | None = None,
) -> dict[str, Any]:
    """The whole comparison: coverage, leave-chromosome-out on training, the pre-registered held-out test."""
    annotate(training + heldout, table, cache)  # both arms in one chromosome-major pass over the cache
    train = [p for p in training if p.covered]
    labels = [p.regulated for p in train]

    loco: dict[str, list[float]] = {}
    for name, cols in FEATURES.items():
        s = [0.0] * len(train)
        for c in sorted({p.chrom for p in train}):
            fit_rows = [p for p in train if p.chrom != c]
            w = logistic_fit(matrix(fit_rows, cols), [p.regulated for p in fit_rows])
            idx = [i for i, p in enumerate(train) if p.chrom == c]
            for i, v in zip(idx, logistic_score(w, matrix([train[i] for i in idx], cols)), strict=True):
                s[i] = v
        loco[name] = s

    weights = {name: logistic_fit(matrix(train, cols), labels) for name, cols in FEATURES.items()}
    held: dict[str, Any] = {}
    for cell in sorted({p.cell for p in heldout}):
        rows = [p for p in heldout if p.cell == cell and p.covered]
        if cell not in MODEL_CELLS:
            held[cell] = {"refused": "no AlphaGenome line for this cell type in the deletion table"}
            continue
        if not any(p.regulated for p in rows):
            held[cell] = {"refused": "no regulated pair on a deleted element"}
            continue
        lab = [p.regulated for p in rows]
        s = {name: logistic_score(weights[name], matrix(rows, cols)) for name, cols in FEATURES.items()}
        with_deletion, without = s["activity + distance + deletion"], s["activity + distance"]
        held[cell] = {
            "models": {name: metrics(v, lab) for name, v in s.items()},
            "deletion_gain": gain_interval(with_deletion, without, rows),
            "passes": (average_precision(with_deletion, lab) or 0) > (average_precision(without, lab) or 0),
        }
    seen = defaultdict(list)
    for p in training:
        seen[p.chrom].append((p.start, p.end))
    disjoint: dict[str, Any] = {}
    for cell, h in held.items():
        if "passes" not in h:
            continue
        rows = [
            p
            for p in heldout
            if p.cell == cell and p.covered and not any(p.start < e and p.end > s for s, e in seen[p.chrom])
        ]
        lab = [p.regulated for p in rows]
        s = {name: logistic_score(weights[name], matrix(rows, cols)) for name, cols in FEATURES.items()}
        disjoint[cell] = {
            "models": {name: metrics(v, lab) for name, v in s.items()},
            "deletion_gain": gain_interval(
                s["activity + distance + deletion"], s["activity + distance"], rows
            ),
        }
    judged = [v for v in held.values() if "passes" in v]
    return {
        "evidence": EVIDENCE,
        "preregistered": PREREGISTERED,
        "verdict": "passed" if judged and all(v["passes"] for v in judged) else "failed",
        "deletion_reader": (
            "the sweep's per-element response cache (data/knowledge/alphagenome/elements), every "
            "gene in the scorer's 1 Mb window on the cell's own track"
            if cache is not None
            else "the compact all_elements table, one gene per element"
        ),
        "deletion_census": {
            "training": deletion_census(training),
            "heldout": deletion_census(heldout),
        },
        "coverage": {
            "training_pairs": len(training),
            "training_pairs_on_a_deleted_element": len(train),
            "training_positives_on_a_deleted_element": sum(labels),
            "heldout_pairs": len(heldout),
            "heldout_pairs_on_a_deleted_element": sum(p.covered for p in heldout),
        },
        "coverage_by_arm": {
            "training": coverage_by_arm(training),
            "heldout": coverage_by_arm(heldout),
        },
        "training_single_predictors": {
            "distance": metrics([-p.features["log_distance"] for p in train], labels),
            "activity over distance": metrics([p.features["activity_over_distance"] for p in train], labels),
            "deletion (top target x drop)": metrics(
                [p.features["top_target"] * (1 + p.features["deletion_drop"]) for p in train], labels
            ),
        },
        "training_leave_chromosome_out": {
            "models": {name: metrics(s, labels) for name, s in loco.items()},
            "deletion_gain": gain_interval(
                loco["activity + distance + deletion"], loco["activity + distance"], train
            ),
        },
        "weights": {
            name: dict(zip(("intercept", *cols), (round(v, 4) for v in weights[name]), strict=True))
            for name, cols in FEATURES.items()
        },
        "heldout": held,
        "heldout_elements_not_in_training": disjoint,
        "one_call_per_element_k562": {str(t): one_call_per_element(train, t) for t in (0.0, 0.1, 0.2)},
        "element_level": {
            "K562 training": element_level(train),
            **{
                f"{c} held out": element_level([p for p in heldout if p.cell == c and p.covered])
                for c in sorted({p.cell for p in heldout})
                if any(p.regulated for p in heldout if p.cell == c and p.covered)
            },
        },
    }


def annotate_contact(
    pairs: list[Pair], source: Any, progress: Any = None, every: int = 1_000
) -> dict[str, Any]:
    """Fill each pair's measured-contact features from a Hi-C matrix source; report what was found.

    The pairs are visited in genomic order per cell and chromosome, because the source reads whole
    compressed blocks of the matrix and one block covers megabases: in that order each block is read
    once. A pair whose cell line has no matrix, whose gene has no TSS column, or whose bin the
    balancing leaves undefined gets zeros and `contact` stays None, and is counted here rather than
    silently scored as a measured contact of zero.
    """
    order = sorted(range(len(pairs)), key=lambda i: (pairs[i].cell, pairs[i].chrom, pairs[i].midpoint))
    counts: dict[str, int] = defaultdict(int)
    for n, i in enumerate(order, 1):
        p = pairs[i]
        c, why = None, ""
        if p.cell not in getattr(source, "matrices", {}):
            why = "no_matrix_for_this_cell_line"
        elif p.tss is None:
            why = "no_tss_column_for_this_pair"
        else:
            c = source.contact(p.cell, p.chrom, p.midpoint, p.tss)
            why = "" if c is not None else "no_balanced_bin_in_the_matrix"
        p.contact = c
        activity = math.log1p(math.sqrt(max(p.dhs, 0.0) * max(p.h3k27ac, 0.0)))
        if c is None:
            counts["pairs_without_contact"] += 1
            counts[f"without_contact_{why}"] += 1
            p.features.update({"log_contact": 0.0, "log_oe": 0.0, "activity_x_contact": 0.0, "same_bin": 0.0})
        else:
            counts["pairs_with_contact"] += 1
            counts[f"balanced_by_{c['norm']}"] += 1
            counts["same_bin"] += int(bool(c["same_bin"]))
            counts["zero_contact"] += int(c["observed"] <= 0)
            p.features.update(
                {
                    "log_contact": math.log1p(max(c["observed"], 0.0)),
                    "log_oe": math.log1p(max(c["oe"] or 0.0, 0.0)),
                    "activity_x_contact": activity + math.log1p(max(c["observed"], 0.0)),
                    "same_bin": float(bool(c["same_bin"])),
                }
            )
        if progress and n % every == 0:
            progress(f"contacts: {n}/{len(order)} pairs read ({counts['pairs_with_contact']} measured)")
    return dict(counts)


def contact_coverage_by_arm(pairs: list[Pair]) -> dict[str, Any]:
    """Which arm the measured contact reaches, and whether the subset that has one is the easier one.

    The same denominator discipline as `coverage_by_arm`, for the contact feature. A Hi-C bin can be
    unmappable, blacklisted or simply one the balancing did not converge on, and if that happened
    more often on one arm than the other, a rate read on the pairs that have a contact would be a
    statement about the matrix's coverage rather than about the predictors. Both arms are reported,
    and both distance baselines are read twice: on every pair of a deleted element and on the subset
    of those that also have a contact. If the baselines move between the two, the subset is not the
    same problem, and any lift on it has to be read with that in mind.
    """
    covered = [p for p in pairs if p.covered]
    arms = {}
    for arm, flag in (("regulated", True), ("not regulated", False)):
        rows = [p for p in covered if p.regulated is flag]
        with_contact = sum(p.contact is not None for p in rows)
        arms[arm] = {
            "pairs_on_a_deleted_element": len(rows),
            "with_a_measured_contact": with_contact,
            "fraction": round(with_contact / len(rows), 4) if rows else None,
        }
    baselines = {}
    for name, rows in (
        ("on a deleted element", covered),
        ("and with a measured contact", [p for p in covered if p.contact is not None]),
    ):
        lab = [p.regulated for p in rows]
        baselines[name] = {
            "distance": metrics([-p.features["log_distance"] for p in rows], lab),
            "activity over distance": metrics([p.features["activity_over_distance"] for p in rows], lab),
        }
    same_bin = [p for p in covered if p.contact is not None and p.features.get("same_bin")]
    return {
        "arms": arms,
        "baselines": baselines,
        "same_bin_pairs": len(same_bin),
        "same_bin_regulated": sum(p.regulated for p in same_bin),
    }


def contact_reading(judged: dict[str, Any], training: dict[str, Any]) -> str:
    """What the intervals say, beside the pass or fail of the pre-registered rule.

    A rule can pass on a point estimate whose interval straddles zero while the same comparison on
    the larger training sample is clearly negative, which is this measurement's own outcome. The
    verdict field answers the pre-registration; this field answers the question.
    """
    held = (judged.get("contact_gain") or {}).get("ci95")
    train = (training.get("contact_gain") or {}).get("ci95")
    parts = []
    if held and held[0] <= 0 <= held[1]:
        parts.append("the held-out interval straddles zero, so the held-out comparison settles nothing")
    elif held:
        parts.append("the held-out interval is clear of zero")
    if train and train[1] < 0:
        parts.append("the same substitution on the larger training sample is clearly negative")
    elif train and train[0] > 0:
        parts.append("the training sample agrees")
    return "; ".join(parts) if parts else "no interval to read"


def score_contact(
    training: list[Pair], heldout: list[Pair], table: DeletionTable, source: Any, progress: Any = None
) -> dict[str, Any]:
    """Measured contact in place of 1/distance, on the pairs the deletion comparison already used.

    Every model here is fitted and scored on the same pairs: on a deleted element, so the deletion
    features exist, and with a measured contact, so the contact features exist. The pre-registered
    claim is PREREGISTERED_CONTACT, fixed in the code before the held-out pairs were scored.
    """
    annotate(training, table)
    annotate(heldout, table)
    found = {
        "training": annotate_contact(training, source, progress),
        "heldout": annotate_contact(heldout, source, progress),
    }
    train = [p for p in training if p.covered and p.contact is not None]
    labels = [p.regulated for p in train]

    loco: dict[str, list[float]] = {}
    for name, cols in CONTACT_FEATURES.items():
        s = [0.0] * len(train)
        for c in sorted({p.chrom for p in train}):
            fit_rows = [p for p in train if p.chrom != c]
            w = logistic_fit(matrix(fit_rows, cols), [p.regulated for p in fit_rows])
            idx = [i for i, p in enumerate(train) if p.chrom == c]
            for i, v in zip(idx, logistic_score(w, matrix([train[i] for i in idx], cols)), strict=True):
                s[i] = v
        loco[name] = s

    weights = {name: logistic_fit(matrix(train, cols), labels) for name, cols in CONTACT_FEATURES.items()}
    held: dict[str, Any] = {}
    for cell in sorted({p.cell for p in heldout}):
        rows = [p for p in heldout if p.cell == cell and p.covered and p.contact is not None]
        if cell not in CONTACT_CELLS:
            held[cell] = {"refused": "no 4DN Hi-C matrix chosen for this cell type"}
            continue
        if cell not in MODEL_CELLS:
            held[cell] = {"refused": "no AlphaGenome line for this cell type in the deletion table"}
            continue
        if not any(p.regulated for p in rows):
            held[cell] = {"refused": "no regulated pair with both a deletion and a measured contact"}
            continue
        lab = [p.regulated for p in rows]
        s = {n: logistic_score(weights[n], matrix(rows, cols)) for n, cols in CONTACT_FEATURES.items()}
        contact_ap = average_precision(s["activity x contact"], lab) or 0
        distance_ap = average_precision(s["activity + distance"], lab) or 0
        held[cell] = {
            "models": {n: metrics(v, lab) for n, v in s.items()},
            "contact_gain": gain_interval(s["activity x contact"], s["activity + distance"], rows),
            "contact_gain_with_deletion": gain_interval(
                s["activity + contact + deletion"], s["activity + distance + deletion"], rows
            ),
            "contact_added_to_distance": gain_interval(
                s["activity + distance + contact"], s["activity + distance"], rows
            ),
            "passes": contact_ap > distance_ap,
        }
    loco_summary = {
        "models": {name: metrics(v, labels) for name, v in loco.items()},
        "contact_gain": gain_interval(loco["activity x contact"], loco["activity + distance"], train),
        "contact_gain_with_deletion": gain_interval(
            loco["activity + contact + deletion"], loco["activity + distance + deletion"], train
        ),
        "contact_added_to_distance": gain_interval(
            loco["activity + distance + contact"], loco["activity + distance"], train
        ),
    }
    judged = held.get("K562", {})  # the pre-registered claim names held-out K562
    return {
        "evidence": f"{EVIDENCE}; {CONTACT_EVIDENCE}",
        "contact_source": source.provenance() if hasattr(source, "provenance") else {},
        "preregistered": PREREGISTERED_CONTACT,
        "verdict": (
            "passed"
            if judged.get("passes")
            else ("failed" if "passes" in judged else "refused: no held-out K562 pairs to judge")
        ),
        "reading": contact_reading(judged, loco_summary),
        "coverage": {
            "training_pairs": len(training),
            "training_pairs_scored": len(train),
            "training_positives_scored": sum(labels),
            "heldout_pairs": len(heldout),
            "heldout_pairs_scored": sum(bool(p.covered and p.contact is not None) for p in heldout),
            "contacts_found": found,
        },
        "contact_coverage_by_arm": {
            "K562 training": contact_coverage_by_arm(training),
            **{
                f"{cell} held out": contact_coverage_by_arm([p for p in heldout if p.cell == cell])
                for cell in CONTACT_CELLS
                if any(p.cell == cell for p in heldout)
            },
        },
        "training_single_predictors": {
            "distance": metrics([-p.features["log_distance"] for p in train], labels),
            "activity over distance": metrics([p.features["activity_over_distance"] for p in train], labels),
            "measured contact": metrics([p.features["log_contact"] for p in train], labels),
            "measured contact over expected": metrics([p.features["log_oe"] for p in train], labels),
            "activity x measured contact": metrics([p.features["activity_x_contact"] for p in train], labels),
        },
        "training_leave_chromosome_out": loco_summary,
        "weights": {
            name: dict(zip(("intercept", *cols), (round(v, 4) for v in weights[name]), strict=True))
            for name, cols in CONTACT_FEATURES.items()
        },
        "heldout": held,
        "same_bin_pairs": {
            "training": sum(1 for p in train if p.features.get("same_bin")),
            "note": (
                "an element and its TSS inside one 5 kb bin share the diagonal cell of the matrix, which "
                "is the measured self-contact of that bin, not a contact between two places"
            ),
        },
    }


# --- the published baseline on the same pairs, and the coverage-matched arm (PREREGISTERED_PUBLISHED) --
SECOND_CELL_COST = {
    "HCT116": {"tracks": "RNA-seq x3, EFO:0002824", "pairs": 396, "regulated": 34, "covered": 363},
    "WTC11": {"tracks": "RNA-seq x3, EFO:0009747", "pairs": 1921, "regulated": 15, "covered": 1616},
    "Jurkat": {"tracks": "RNA-seq x1, CLO:0007045", "pairs": 75, "regulated": 7, "covered": 64},
}
REQUEST_BUDGET = 20
MATCH_DRAWS = 1_000


def benchmark_auprc(
    scores: list[float], labels: list[bool], weights: list[float] | None = None
) -> float | None:
    """AUPRC by the ENCODE benchmark's own estimator, optionally with case weights.

    yardstick's pr_curve, one point per distinct score from the highest down, with the first row
    (recall 0) dropped as the pipeline does, then the last point dropped (the pipeline's recall-100%
    endpoint), then the trapezoid over recall. Checked against the published figure for distance to
    TSS on both files before any model was scored with it (see PUBLISHED).
    """
    w = weights or [1.0] * len(labels)
    total = sum(x for x, y in zip(w, labels, strict=True) if y)
    if not total or all(labels):
        return None
    order = sorted(range(len(scores)), key=lambda i: -scores[i])
    tp = fp = 0.0
    points: list[tuple[float, float]] = []
    i = 0
    while i < len(order):
        s = scores[order[i]]
        while i < len(order) and scores[order[i]] == s:
            j = order[i]
            if labels[j]:
                tp += w[j]
            else:
                fp += w[j]
            i += 1
        points.append((tp / total, tp / (tp + fp) if tp + fp else 1.0))
    points = points[:-1]
    return sum((b[0] - a[0]) * (a[1] + b[1]) / 2 for a, b in zip(points, points[1:], strict=False))


def weighted_gain(
    a: list[float], b: list[float], pairs: list[Pair], weighted: bool, seed: int = 0, n: int = BOOTSTRAPS
) -> dict[str, Any]:
    """Benchmark-estimator AUPRC of `a` minus `b`, 95% interval from resampling whole chromosomes."""
    by_chrom: dict[str, list[int]] = defaultdict(list)
    for i, p in enumerate(pairs):
        by_chrom[p.chrom].append(i)
    chroms = sorted(by_chrom)
    labels = [p.regulated for p in pairs]
    w = [p.weight if weighted else 1.0 for p in pairs]

    def diff(idx: list[int]) -> float:
        lab, ww = [labels[i] for i in idx], [w[i] for i in idx]
        return (benchmark_auprc([a[i] for i in idx], lab, ww) or 0) - (
            benchmark_auprc([b[i] for i in idx], lab, ww) or 0
        )

    point = diff(list(range(len(pairs))))
    rng = random.Random(seed)
    diffs = []
    for _ in range(n):
        idx = [i for c in (rng.choice(chroms) for _ in chroms) for i in by_chrom[c]]
        if any(labels[i] for i in idx) and not all(labels[i] for i in idx):
            diffs.append(diff(idx))
    diffs.sort()
    ci = [
        round(diffs[int(0.025 * len(diffs))], 4),
        round(diffs[min(len(diffs) - 1, int(0.975 * len(diffs)))], 4),
    ]
    return {"gain": round(point, 4), "ci95": ci if diffs else None, "resamples": len(diffs)}


def band(value: float, published: tuple[float, float, float]) -> str:
    """Where a figure sits against a published figure's 95% interval."""
    _, lo, hi = published
    return (
        "above the published interval"
        if value > hi
        else ("below the published interval" if value < lo else "inside the published interval")
    )


def bench_metrics(scores: list[float], pairs: list[Pair], weighted: bool) -> dict[str, Any]:
    labels = [p.regulated for p in pairs]
    w = [p.weight for p in pairs] if weighted else None
    v = benchmark_auprc(scores, labels, w)
    return {
        "pairs": len(pairs),
        "positives": sum(labels),
        "weighted_positives": round(sum(p.weight for p in pairs if p.regulated), 2) if weighted else None,
        "auprc": None if v is None else round(v, 4),
    }


def coverage_arms(
    heldout_k562: list[Pair], training: list[Pair], weights: dict[str, list[float]], seed: int = 0
) -> dict[str, Any]:
    """The three registered arms for the held-out K562 gain; the headline's estimator (average precision)."""
    with_d, without = FEATURES["activity + distance + deletion"], FEATURES["activity + distance"]

    def gain_on(rows: list[Pair]) -> float:
        lab = [p.regulated for p in rows]
        a = logistic_score(weights["activity + distance + deletion"], matrix(rows, with_d))
        b = logistic_score(weights["activity + distance"], matrix(rows, without))
        return (average_precision(a, lab) or 0) - (average_precision(b, lab) or 0)

    rows = heldout_k562
    s_with = logistic_score(weights["activity + distance + deletion"], matrix(rows, with_d))
    s_without = logistic_score(weights["activity + distance"], matrix(rows, without))
    arm1 = {
        "pairs": len(rows),
        "positives": sum(p.regulated for p in rows),
        "uncovered_pairs": sum(not p.covered for p in rows),
        **gain_interval(s_with, s_without, rows),
    }

    covered = [p for p in rows if p.covered]
    reg_all = [p for p in rows if p.regulated]
    neg_all = [p for p in rows if not p.regulated]
    neg_rate = sum(p.covered for p in neg_all) / len(neg_all)
    reg_cov = [p for p in covered if p.regulated]
    neg_cov = [p for p in covered if not p.regulated]
    keep = round(neg_rate * len(reg_all))
    rng = random.Random(seed)
    gains = sorted(gain_on(neg_cov + rng.sample(reg_cov, keep)) for _ in range(MATCH_DRAWS))
    arm2 = {
        "regulated_coverage_before": round(len(reg_cov) / len(reg_all), 4),
        "non_regulated_coverage": round(neg_rate, 4),
        "regulated_kept": keep,
        "regulated_covered": len(reg_cov),
        "draws": MATCH_DRAWS,
        "median_gain": round(gains[len(gains) // 2], 4),
        "range": [round(gains[0], 4), round(gains[-1], 4)],
        "share_above_zero": round(sum(g > 0 for g in gains) / len(gains), 4),
    }

    control_cols = ("log_distance", "log_activity", "activity_over_distance", "covered")
    w_ctrl = logistic_fit(matrix(training, control_cols), [p.regulated for p in training])
    w_base = logistic_fit(matrix(training, without), [p.regulated for p in training])
    s_ctrl = logistic_score(w_ctrl, matrix(rows, control_cols))
    s_base = logistic_score(w_base, matrix(rows, without))
    arm3 = {
        "model": "activity + distance + covered (fitted on all training pairs)",
        "against": "activity + distance (fitted on all training pairs)",
        **gain_interval(s_ctrl, s_base, rows),
    }
    ci = arm1["ci95"]
    invariant = bool(ci and ci[0] > 0) and arm2["share_above_zero"] >= 0.95 and arm3["gain"] < arm1["gain"]
    return {
        "arm1_all_pairs": arm1,
        "arm2_coverage_matched": arm2,
        "arm3_coverage_indicator_control": arm3,
        "invariant_to_coverage": invariant,
    }


def score_published(
    training: list[Pair],
    heldout: list[Pair],
    table: DeletionTable,
    cache: ElementCache | None = None,
) -> dict[str, Any]:
    """PREREGISTERED_PUBLISHED: the frozen model on the published pair sets, and the coverage arms."""
    annotate(training + heldout, table, cache)
    for p in training + heldout:
        p.features["covered"] = float(p.covered)
    train = [p for p in training if p.covered]
    weights = {
        name: logistic_fit(matrix(train, cols), [p.regulated for p in train])
        for name, cols in FEATURES.items()
    }

    raw_distance = {
        "training_unweighted": round(
            benchmark_auprc([-abs(p.distance) for p in training], [p.regulated for p in training]) or 0, 4
        ),
        "heldout_weighted": round(
            benchmark_auprc(
                [-abs(p.distance) for p in heldout],
                [p.regulated for p in heldout],
                [p.weight for p in heldout],
            )
            or 0,
            4,
        ),
        "published": {
            "training": PUBLISHED["training"]["Distance to TSS"][0],
            "heldout": PUBLISHED["heldout"]["Distance to TSS"][0],
        },
    }

    # 1. training, hold-one-chromosome-out over all 10,356 pairs, fitted on the covered pairs elsewhere
    loco: dict[str, list[float]] = {name: [0.0] * len(training) for name in FEATURES}
    for c in sorted({p.chrom for p in training}):
        fit_rows = [p for p in train if p.chrom != c]
        idx = [i for i, p in enumerate(training) if p.chrom == c]
        for name, cols in FEATURES.items():
            w = logistic_fit(matrix(fit_rows, cols), [p.regulated for p in fit_rows])
            for i, v in zip(idx, logistic_score(w, matrix([training[i] for i in idx], cols)), strict=True):
                loco[name][i] = v
    train_models = {name: bench_metrics(s, training, weighted=False) for name, s in loco.items()}
    ours_train = train_models["activity + distance + deletion"]["auprc"] or 0.0
    training_block = {
        "models": train_models,
        "deletion_gain": weighted_gain(
            loco["activity + distance + deletion"], loco["activity + distance"], training, weighted=False
        ),
        "published": PUBLISHED["training"],
        "against_encode_re2g": band(ours_train, PUBLISHED["training"]["ENCODE-rE2G"]),
        "against_encode_re2g_extended": band(ours_train, PUBLISHED["training"]["ENCODE-rE2G_Extended"]),
        "baseline_against_abc": band(
            train_models["activity + distance"]["auprc"] or 0.0,
            PUBLISHED["training"]["ABC_A=DNase, C=Average ENCODE Hi-C"],
        ),
    }

    # 2. held-out, all 4,378 pairs pooled, weighted, frozen weights
    s = {name: logistic_score(weights[name], matrix(heldout, cols)) for name, cols in FEATURES.items()}
    pooled = {name: bench_metrics(v, heldout, weighted=True) for name, v in s.items()}
    ours_held = pooled["activity + distance + deletion"]["auprc"] or 0.0
    per_cell = {}
    for cell in sorted({p.cell for p in heldout}):
        idx = [i for i, p in enumerate(heldout) if p.cell == cell]
        rows = [heldout[i] for i in idx]
        if not any(p.regulated for p in rows):
            continue
        per_cell[cell] = {
            "models": {
                name: bench_metrics([v[i] for i in idx], rows, weighted=True) for name, v in s.items()
            },
            "deletion_gain": weighted_gain(
                [s["activity + distance + deletion"][i] for i in idx],
                [s["activity + distance"][i] for i in idx],
                rows,
                weighted=True,
            ),
            "deletion_available": cell in MODEL_CELLS,
        }
    heldout_block = {
        "models": pooled,
        "deletion_gain": weighted_gain(
            s["activity + distance + deletion"], s["activity + distance"], heldout, weighted=True
        ),
        "pairs_without_a_deletion_value": sum(p.cell not in MODEL_CELLS for p in heldout),
        "published": PUBLISHED["heldout"],
        "against_encode_re2g": band(ours_held, PUBLISHED["heldout"]["ENCODE-rE2G"]),
        "baseline_against_abc": band(
            pooled["activity + distance"]["auprc"] or 0.0,
            PUBLISHED["heldout"]["ABC_A=DNase, C=Average ENCODE Hi-C"],
        ),
        "per_cell_type_weighted": per_cell,
    }

    # 3. the coverage arms on held-out K562
    arms = coverage_arms([p for p in heldout if p.cell == "K562"], training, weights)

    # 5, NOT registered: found after the registered numbers were read (see POSITIVE_FILTER)
    diagnostic = dnase_only_diagnostic(training, heldout)

    # 4. the second cell type, costed from the join itself, not bought
    cost: dict[str, Any] = {}
    for cell in SECOND_CELL_COST:
        ids = {e["id"] for p in heldout if p.cell == cell for e in table.overlapping(p.chrom, p.start, p.end)}
        cost[cell] = {**SECOND_CELL_COST[cell], "requests_frozen_feature": len(ids)}
    return {
        "evidence": EVIDENCE,
        "preregistered": PREREGISTERED_PUBLISHED,
        "published_source": PUBLISHED_SOURCE,
        "estimator_check_raw_distance": raw_distance,
        "training_published_split": training_block,
        "heldout_published_pairs": heldout_block,
        "coverage_arms_k562_heldout": arms,
        "second_cell_type": {
            "cost": cost,
            "budget": REQUEST_BUDGET,
            "run": False,
            "why": (
                "HCT116 needs more deletion requests than the lane's budget; costed and handed to "
                "the coordinator"
            ),
        },
        "post_hoc_positive_filter": {
            "registered": False,
            "finding": POSITIVE_FILTER,
            "chromatin_category_by_label": chromatin_by_label(heldout),
            "dnase_only_diagnostic": diagnostic,
        },
        "alphagenome_requests": 0,
    }


#: Found on 2026-09-27 AFTER the registered held-out numbers were read, and therefore not part of the
#: registration: the held-out file's positives were filtered on chromatin at the tested element
#: (Gschwind et al. 2026, Methods: "filters based on effect size and chromatin state"), so every one
#: of its 190 positives lies in an H3K27ac element while 1,438 of its 4,188 negatives do not. This
#: module's activity term is sqrt(DNase x H3K27ac) read from the same file, so on the held-out set
#: it is partly reading the rule that made a pair positive. ENCODE-rE2G's held-out figure is its
#: DNase-only model, which does not read H3K27ac. The registered held-out band is therefore not a
#: fair comparison in this module's favour, and the diagnostic below re-fits both models with DNase
#: alone, as the published model has it. It is labelled post hoc everywhere it is quoted.
POSITIVE_FILTER = (
    "held-out positives were selected on H3K27ac at the element; this module's activity reads H3K27ac; "
    "ENCODE-rE2G's held-out model reads DNase only, so the registered held-out band flatters this module"
)
DNASE_FEATURES = {
    "dnase + distance": ("log_distance", "log_dnase", "dnase_over_distance"),
    "dnase + distance + deletion": (
        "log_distance",
        "log_dnase",
        "dnase_over_distance",
        "top_target",
        "deletion_drop",
    ),
}


def chromatin_by_label(pairs: list[Pair]) -> dict[str, Any]:
    """The benchmark's chromatin category of the tested element, counted per label."""
    out: dict[str, Any] = {}
    for label, flag in (("regulated", True), ("not regulated", False)):
        counts: dict[str, int] = defaultdict(int)
        for p in pairs:
            if p.regulated is flag:
                counts[p.category or "unknown"] += 1
        out[label] = dict(sorted(counts.items()))
    return out


def dnase_only_diagnostic(training: list[Pair], heldout: list[Pair]) -> dict[str, Any]:
    """Post hoc: the same two models with DNase in place of sqrt(DNase x H3K27ac). Pairs already annotated."""
    for p in training + heldout:
        d = max(MIN_DISTANCE, p.distance)
        p.features["log_dnase"] = math.log1p(max(p.dhs, 0.0))
        p.features["dnase_over_distance"] = p.features["log_dnase"] - math.log(d)
    train = [p for p in training if p.covered]
    lab = [p.regulated for p in train]
    weights = {n: logistic_fit(matrix(train, c), lab) for n, c in DNASE_FEATURES.items()}
    loco: dict[str, list[float]] = {n: [0.0] * len(training) for n in DNASE_FEATURES}
    for c in sorted({p.chrom for p in training}):
        fit_rows = [p for p in train if p.chrom != c]
        idx = [i for i, p in enumerate(training) if p.chrom == c]
        for n, cols in DNASE_FEATURES.items():
            w = logistic_fit(matrix(fit_rows, cols), [p.regulated for p in fit_rows])
            for i, v in zip(idx, logistic_score(w, matrix([training[i] for i in idx], cols)), strict=True):
                loco[n][i] = v
    s = {n: logistic_score(weights[n], matrix(heldout, c)) for n, c in DNASE_FEATURES.items()}
    k562 = [p for p in heldout if p.cell == "K562" and p.covered]
    sk = {n: logistic_score(weights[n], matrix(k562, c)) for n, c in DNASE_FEATURES.items()}
    with_d, without = "dnase + distance + deletion", "dnase + distance"
    return {
        "training_published_split": {
            "models": {n: bench_metrics(v, training, weighted=False) for n, v in loco.items()},
            "deletion_gain": weighted_gain(loco[with_d], loco[without], training, weighted=False),
        },
        "heldout_pooled_weighted": {
            "models": {n: bench_metrics(v, heldout, weighted=True) for n, v in s.items()},
            "deletion_gain": weighted_gain(s[with_d], s[without], heldout, weighted=True),
            "against_encode_re2g": band(
                benchmark_auprc(s[with_d], [p.regulated for p in heldout], [p.weight for p in heldout])
                or 0.0,
                PUBLISHED["heldout"]["ENCODE-rE2G"],
            ),
            "baseline_against_abc": band(
                benchmark_auprc(s[without], [p.regulated for p in heldout], [p.weight for p in heldout])
                or 0.0,
                PUBLISHED["heldout"]["ABC_A=DNase, C=Average ENCODE Hi-C"],
            ),
        },
        "heldout_k562_covered_headline_estimator": {
            "models": {n: metrics(v, [p.regulated for p in k562]) for n, v in sk.items()},
            "deletion_gain": gain_interval(sk[with_d], sk[without], k562),
        },
    }
