# SPDX-License-Identifier: AGPL-3.0-or-later
"""The experimental layer of the compiled genome: what an assay measured over a compiled element.

The compiled per-chromosome programs state 940,803 facts and not one of them is `experimental`:
a region's role comes from its budget tier (`inferred`) and an element's target from a deletion in
AlphaGenome (`predicted`). Meanwhile the project holds real measurements over some of that same
sequence. This module attaches them, element by element, and counts honestly how thin the layer is.

Three assays, each kept under its own name:

- **CRISPRi** (ENCODE enhancer-gene benchmark, Gschwind et al. 2025): an element was silenced in its
  own chromosome and nearby genes were measured. A pair is `regulated` or not, and **a pair measured
  as not regulated is evidence, not an absence of evidence**: it is carried through to the program as
  an experimental fact, never dropped and never turned back into UNKNOWN.
- **lentiMPRA** (ENCODE4 joint library ENCSR106SZM, K562 / HepG2 / WTC11): how much a 200 bp sequence
  drives transcription from a reporter. Episomal, so it measures the sequence and not the locus.
- **VISTA** (LBNL, transgenic mouse e11.5): whether a sequence is an enhancer in a living embryo,
  positive or negative.

**The overlap rule, stated as a constant.** A measurement is *of* a compiled element only when the
tested interval and the element are largely the same piece of DNA: `RECIPROCAL_OVERLAP = 0.5`, that
is, the overlap covers at least half of the element **and** at least half of the tested interval.
Nothing is upgraded on proximity, similarity or a single shared base. `sensitivity()` reports the
matches the rule makes at 0.25, 0.5 and 0.75 so the choice can be seen rather than trusted, and it
also reports `containment` — measured intervals that merely swallow the element — as a number that is
deliberately **not** raised to experimental.

**Two names, not one.** A measured element becomes a second block, `<id>_measured`, beside the
predicted `<id>`; the prediction is never overwritten. In these rows the prediction lives under
`predicted_*` and the measurement under `measured`, and the two are compared in a third field,
`agreement`, which is a reading and not a fact.

**The agreement rate is not computed over a subset the prediction chose.** For CRISPRi the
denominator is every element the screen touched, split three ways: the predicted gene was tested and
regulated (agrees), tested and not regulated (disagrees), or not tested at all. The narrower rate
over "the predicted gene was tested" is reported too, under its own name, because it is the number a
careless reader would quote.
"""

from __future__ import annotations

import bisect
import csv
import gzip
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from genomeos.attribution import mpra, vista
from genomeos.attribution.crispri import KNOWLEDGE as CRISPRI_KNOWLEDGE
from genomeos.results import RESULTS_DIR

# --- the rule -------------------------------------------------------------------------------------
RECIPROCAL_OVERLAP = 0.5  # the measurement is of this element only at or above this, both ways
OVERLAP_SENSITIVITY = (0.25, 0.5, 0.75)  # the rule reported at three values, never at one
# the weaker bar that decides what could have been found at all, before asking what was
ELIGIBILITY_RULE = "an assay's tested interval shares at least one base with the element"

CRISPRI_FILES = (
    "EPCrisprBenchmark_combined_data.training_K562.GRCh38.tsv.gz",
    "EPCrisprBenchmark_combined_data.heldout_5_cell_types.GRCh38.tsv.gz",
)
MPRA_ACTIVE = mpra.ACTIVE  # log2(RNA/DNA) at or above which a reporter element counts as active
REACH = 200_000  # the longest measured interval any assay holds, for the bounded overlap scan

# an endogenous perturbation or an in-vivo assay; a reporter measures the sequence, not the locus
PERTURBATION_CONFIDENCE = 0.9
REPORTER_CONFIDENCE = 0.75

SOURCES = {
    "crispri": (
        "CRISPRi enhancer-gene screens, ENCODE benchmark (EngreitzLab/CRISPR_comparison, "
        "Gschwind et al. 2025)"
    ),
    "lentimpra": "ENCODE4 lentiMPRA, joint library ENCSR106SZM (Ahituv lab), log2(RNA/DNA) per 200 bp",
    "vista": "VISTA Enhancer Browser, transgenic mouse e11.5 (LBNL), hg38 coordinates",
}
ASSAYS = tuple(SOURCES)


def reciprocal_overlap(a_start: int, a_end: int, b_start: int, b_end: int) -> float:
    """The smaller of the two overlap fractions: 1.0 for identical intervals, 0.0 for disjoint ones."""
    ov = min(a_end, b_end) - max(a_start, b_start)
    if ov <= 0 or a_end <= a_start or b_end <= b_start:
        return 0.0
    return min(ov / (a_end - a_start), ov / (b_end - b_start))


def measures(
    a_start: int, a_end: int, b_start: int, b_end: int, fraction: float = RECIPROCAL_OVERLAP
) -> bool:
    """Whether a measurement of (b_start, b_end) is a measurement of the element (a_start, a_end)."""
    return reciprocal_overlap(a_start, a_end, b_start, b_end) >= fraction


def contains(a_start: int, a_end: int, b_start: int, b_end: int) -> bool:
    """Whether the tested interval swallows the element. Counted, and never upgraded on."""
    return b_start <= a_start and b_end >= a_end


# --- the assays -----------------------------------------------------------------------------------
@dataclass(frozen=True, slots=True)
class CrispriPair:
    """One element-gene pair a screen measured. `regulated` False is a measurement, not a gap."""

    chrom: str
    start: int
    end: int
    gene: str
    cell: str
    dataset: str
    reference: str
    regulated: bool
    significant: bool
    effect_size: float
    p_adjusted: float


def parse_crispri(lines: Any, chrom: str | None = None) -> tuple[list[CrispriPair], int]:
    """Valid pairs and the count the benchmark itself marks invalid (promoter or exon overlaps).

    The invalid ones are returned as a number rather than silently skipped: they are not measured
    negatives, they are pairs the benchmark says are not a test of an enhancer-to-gene link.
    """
    out: list[CrispriPair] = []
    invalid = 0
    for r in csv.DictReader(lines, delimiter="\t"):
        if chrom and r["chrom"] != chrom:
            continue
        if r.get("ValidConnection") != "TRUE":
            invalid += 1
            continue
        out.append(
            CrispriPair(
                chrom=r["chrom"],
                start=int(r["chromStart"]),
                end=int(r["chromEnd"]),
                gene=r["measuredGeneSymbol"],
                cell=r["CellType"],
                dataset=r.get("Dataset", ""),
                reference=r.get("Reference", ""),
                regulated=r.get("Regulated") == "TRUE",
                significant=r.get("Significant") == "TRUE",
                effect_size=float(r["EffectSize"] or 0.0),
                p_adjusted=float(r["pValueAdjusted"] or 1.0),
            )
        )
    return out, invalid


def load_crispri(
    chrom: str | None = None, knowledge: Path = CRISPRI_KNOWLEDGE
) -> tuple[list[CrispriPair], int]:
    """Both benchmark tables from the local cache. Nothing is fetched; a missing table is no pairs."""
    pairs: list[CrispriPair] = []
    invalid = 0
    for name in CRISPRI_FILES:
        p = knowledge / name
        if not p.exists():
            continue
        with gzip.open(p, "rt") as fh:
            got, bad = parse_crispri(fh, chrom)
        pairs.extend(got)
        invalid += bad
    pairs.sort(key=lambda p: (p.start, p.end, p.gene, p.cell))
    return pairs, invalid


def load_mpra(chrom: str, knowledge: Path = mpra.KNOWLEDGE) -> list[mpra.Element]:
    """The lentiMPRA elements of one chromosome from the local cache, without fetching."""
    into: dict[str, mpra.Element] = {}
    for cell, acc in mpra.FILES.items():
        p = knowledge / f"{acc}.bed.gz"
        if not p.exists():
            continue
        with gzip.open(p, "rt") as fh:
            mpra.parse(fh, cell, chrom, into)
    out = [e for e in into.values() if e.activity]
    out.sort(key=lambda e: e.start)
    return out


def load_vista(chrom: str, knowledge: Path = vista.KNOWLEDGE) -> list[vista.VistaElement]:
    """The VISTA elements of one chromosome from the local cache, without fetching."""
    p = vista.locus_path(knowledge)
    if not p.exists():
        return []
    with gzip.open(p, "rt") as fh:
        return vista.parse_loci(fh, chrom)


# --- the layer ------------------------------------------------------------------------------------
@dataclass
class Layer:
    """Every measurement over one chromosome, each assay under its own name and searchable by overlap."""

    chrom: str
    crispri: list[CrispriPair] = field(default_factory=list)
    lentimpra: list[mpra.Element] = field(default_factory=list)
    vista: list[vista.VistaElement] = field(default_factory=list)
    crispri_invalid: int = 0
    _starts: dict[str, list[int]] = field(default_factory=dict, repr=False)

    def __post_init__(self) -> None:
        self.crispri = sorted(self.crispri, key=lambda p: p.start)
        self.lentimpra = sorted(self.lentimpra, key=lambda e: e.start)
        self.vista = sorted(self.vista, key=lambda e: e.start)
        self._starts = {name: [x.start for x in getattr(self, name)] for name in ASSAYS}

    @classmethod
    def load(cls, chrom: str) -> Layer:
        pairs, invalid = load_crispri(chrom)
        return cls(
            chrom=chrom,
            crispri=pairs,
            lentimpra=load_mpra(chrom),
            vista=load_vista(chrom),
            crispri_invalid=invalid,
        )

    @property
    def empty(self) -> bool:
        return not (self.crispri or self.lentimpra or self.vista)

    def counts(self) -> dict[str, int]:
        return {name: len(getattr(self, name)) for name in ASSAYS}

    def touched_by(self, start: int, end: int, reach: int = REACH) -> list[str]:
        """The assays whose tested intervals share even one base with this element.

        This is the **eligibility** denominator, and it is deliberately weaker than the overlap rule.
        An element no assay ever covered cannot be raised to experimental, and it must not sit in the
        same bucket as one that was covered and disagreed: the first is a statement about where the
        screens were pointed, the second about the prediction. `ELIGIBILITY_RULE` names the bar.
        """
        return [a for a in ASSAYS if self.near(a, start, end, reach)]

    def near(self, name: str, start: int, end: int, reach: int = REACH) -> list:
        """Everything of one assay whose interval can touch (start, end), by a bounded scan."""
        items, starts = getattr(self, name), self._starts[name]
        lo = bisect.bisect_left(starts, start - reach)
        hi = bisect.bisect_right(starts, end)
        return [x for x in items[lo:hi] if x.end > start and x.start < end]

    def for_element(
        self, start: int, end: int, fraction: float = RECIPROCAL_OVERLAP, reach: int = REACH
    ) -> dict[str, Any]:
        """The measurements *of* this interval, by assay. An assay with no match is absent from the
        result; an assay that measured no effect is present, with the genes it found no effect on."""
        out: dict[str, Any] = {}

        near = self.near("crispri", start, end, reach)
        pairs = [(p, reciprocal_overlap(start, end, p.start, p.end)) for p in near]
        pairs = [(p, f) for p, f in pairs if f >= fraction]
        if pairs:
            regulated = sorted({p.gene for p, _ in pairs if p.regulated})
            tested = sorted({p.gene for p, _ in pairs})
            out["crispri"] = {
                "genes_tested": tested,
                "genes_regulated": regulated,
                "genes_not_regulated": [g for g in tested if g not in set(regulated)],
                "cells": sorted({p.cell for p, _ in pairs}),
                "overlap": round(max(f for _, f in pairs), 3),
                "pairs": [
                    {
                        "gene": p.gene,
                        "cell": p.cell,
                        "dataset": p.dataset,
                        "regulated": p.regulated,
                        "effect_size": round(p.effect_size, 4),
                        "p_adjusted": p.p_adjusted,
                        "overlap": round(f, 3),
                    }
                    for p, f in sorted(pairs, key=lambda x: (x[0].gene, x[0].cell))
                ],
            }

        hits = [
            (e, reciprocal_overlap(start, end, e.start, e.end))
            for e in self.near("lentimpra", start, end, reach)
        ]
        hits = [(e, f) for e, f in hits if f >= fraction]
        if hits:
            best: dict[str, float] = {}
            for e, _ in hits:
                for cell, value in e.activity.items():
                    best[cell] = max(best.get(cell, value), value)
            out["lentimpra"] = {
                "activity": {c: round(v, 4) for c, v in sorted(best.items())},
                "cells_active": sorted(c for c, v in best.items() if v >= MPRA_ACTIVE),
                "cells_silent": sorted(c for c, v in best.items() if v < MPRA_ACTIVE),
                "max_activity": round(max(best.values()), 4) if best else None,
                "overlap": round(max(f for _, f in hits), 3),
                "elements": sorted({e.name for e, _ in hits}),
            }

        vs = [
            (e, reciprocal_overlap(start, end, e.start, e.end)) for e in self.near("vista", start, end, reach)
        ]
        vs = [(e, f) for e, f in vs if f >= fraction]
        if vs:
            out["vista"] = {
                "positive": sorted(e.id for e, _ in vs if e.status == "positive"),
                "negative": sorted(e.id for e, _ in vs if e.status != "positive"),
                "tissues": sorted({t for e, _ in vs for t in e.tissues}),
                "overlap": round(max(f for _, f in vs), 3),
            }
        return out


# --- prediction against measurement ---------------------------------------------------------------
AGREES, DISAGREES, NOT_TESTED = "agrees", "disagrees", "predicted_gene_not_tested"


def agreement(predicted_gene: str, measured: dict[str, Any]) -> dict[str, Any]:
    """What each assay says about the compiled claim, assay by assay and then pooled.

    The compiled claim of an element block is that deleting it moves `predicted_gene`. CRISPRi asks
    exactly that question and can answer it three ways. lentiMPRA and VISTA ask the weaker question
    of whether the sequence acts at all, so their verdicts are recorded under their own names and the
    pooled verdict says how many assays fell each way rather than hiding the mixture.
    """
    out: dict[str, Any] = {}
    c = measured.get("crispri")
    if c:
        if predicted_gene in c["genes_regulated"]:
            out["crispri"] = AGREES
        elif predicted_gene in c["genes_not_regulated"]:
            out["crispri"] = DISAGREES
        else:
            out["crispri"] = NOT_TESTED
        out["crispri_detail"] = (
            "another gene regulated" if c["genes_regulated"] else "no gene regulated in the screen"
        )
    m = measured.get("lentimpra")
    if m:
        out["lentimpra"] = AGREES if m["cells_active"] else DISAGREES
    v = measured.get("vista")
    if v:
        out["vista"] = AGREES if v["positive"] else DISAGREES
    verdicts = [out[a] for a in ASSAYS if out.get(a) in (AGREES, DISAGREES)]
    out["assays_agreeing"] = sum(1 for x in verdicts if x == AGREES)
    out["assays_disagreeing"] = sum(1 for x in verdicts if x == DISAGREES)
    out["verdict"] = (
        "mixed"
        if out["assays_agreeing"] and out["assays_disagreeing"]
        else AGREES
        if out["assays_agreeing"]
        else DISAGREES
        if out["assays_disagreeing"]
        else NOT_TESTED
    )
    return out


def confidence_of(measured: dict[str, Any]) -> float:
    """The strongest assay present decides: a perturbation or an embryo outranks a reporter."""
    if measured.get("crispri") or measured.get("vista"):
        return PERTURBATION_CONFIDENCE
    return REPORTER_CONFIDENCE


def sources_of(measured: dict[str, Any]) -> str:
    return "; ".join(SOURCES[a] for a in ASSAYS if a in measured)


def rows(
    chrom: str,
    elements: list[dict[str, Any]] | None = None,
    layer: Layer | None = None,
    fraction: float = RECIPROCAL_OVERLAP,
    results_dir: Path = RESULTS_DIR,
) -> list[dict[str, Any]]:
    """One row per compiled element that a measurement of that same element exists for."""
    if elements is None:
        from genomeos.attribution.targets import attributed

        elements = attributed(chrom, results_dir)
    layer = layer if layer is not None else Layer.load(chrom)
    out: list[dict[str, Any]] = []
    for e in elements:
        m = layer.for_element(e["start"], e["end"], fraction)
        if not m:
            continue
        pc = e.get("predicted_coding") or {}
        out.append(
            {
                "id": e["id"],
                "chrom": chrom,
                "start": e["start"],
                "end": e["end"],
                "length": e["end"] - e["start"],
                "domain": e.get("domain", ""),
                "predicted_gene": pc.get("gene", ""),
                "predicted_log2_fold_change": pc.get("log2_fold_change"),
                "predicted_action": pc.get("action"),
                "measured": m,
                "agreement": agreement(pc.get("gene", ""), m),
                "assays": sorted(m),
                "confidence": confidence_of(m),
            }
        )
    out.sort(key=lambda r: r["start"])
    return out


def eligibility(elements: list[dict[str, Any]], layer: Layer) -> dict[str, Any]:
    """What could have been found at all, before asking what was: the second of three denominators.

    A small census invites the reader to conclude something about the genome. It is a statement about
    where the assays were pointed, and the way to make that unmistakable is to say how many compiled
    elements lie inside any assay's footprint under `ELIGIBILITY_RULE` - a single shared base - before
    saying how many met the overlap rule. An element nothing ever covered cannot be raised, and it
    belongs in its own bucket, not with the ones that were covered and disagreed.
    """
    per = {a: 0 for a in ASSAYS}
    eligible = 0
    for e in elements:
        touched = layer.touched_by(e["start"], e["end"])
        for a in touched:
            per[a] += 1
        eligible += bool(touched)
    n = len(elements)
    return {
        "rule": ELIGIBILITY_RULE,
        "compiled_elements": n,
        "elements_in_an_assay_footprint": eligible,
        "elements_no_assay_ever_covered": n - eligible,
        "share_eligible": round(eligible / n, 5) if n else 0.0,
        "by_assay": per,
    }


# --- the census -----------------------------------------------------------------------------------
def census(
    chrom: str, measured_rows: list[dict[str, Any]], elements: list[dict[str, Any]], layer: Layer
) -> dict[str, Any]:
    """How thin the experimental layer is on one chromosome, in counts and shares.

    Three denominators, in order, because the census is a statement about the search and not about the
    genome: how many compiled elements there are; how many were **eligible** for any measurement at
    all (`eligibility`); and only then how many were raised, agreed and disagreed.

    Two names for two questions that are easy to confuse: `facts_with_a_measured_counterpart` counts
    the compiled *predicted* facts (an element block and its rule) that a measurement of the same
    element now sits beside; `experimental_facts_added` counts the blocks the compiler writes as a
    consequence. Nothing is overwritten, so no fact is "converted" and neither number is that.
    """
    n_elements = len(elements)
    eligible = eligibility(elements, layer)
    by_assay = {a: sum(1 for r in measured_rows if a in r["measured"]) for a in ASSAYS}
    agree = {a: sum(1 for r in measured_rows if r["agreement"].get(a) == AGREES) for a in ASSAYS}
    disagree = {a: sum(1 for r in measured_rows if r["agreement"].get(a) == DISAGREES) for a in ASSAYS}
    not_tested = sum(1 for r in measured_rows if r["agreement"].get("crispri") == NOT_TESTED)
    regulated_pairs = sum(
        1 for r in measured_rows for p in r["measured"].get("crispri", {}).get("pairs", []) if p["regulated"]
    )
    negative_pairs = sum(
        1
        for r in measured_rows
        for p in r["measured"].get("crispri", {}).get("pairs", [])
        if not p["regulated"]
    )
    measured_genes = {
        g for r in measured_rows for g in r["measured"].get("crispri", {}).get("genes_regulated", [])
    }
    crispri_matched = by_assay["crispri"]
    tested = agree["crispri"] + disagree["crispri"]
    n_eligible = eligible["elements_in_an_assay_footprint"]
    return {
        "chrom": chrom,
        # 1. what there is, 2. what could have been found, 3. what was
        "compiled_elements": n_elements,
        "eligibility": eligible,
        "elements_in_an_assay_footprint": n_eligible,
        "elements_no_assay_ever_covered": eligible["elements_no_assay_ever_covered"],
        "elements_with_any_measurement": len(measured_rows),
        "coverage_elements_measured": round(len(measured_rows) / n_elements, 5) if n_elements else 0.0,
        "raised_share_of_the_eligible": (round(len(measured_rows) / n_eligible, 5) if n_eligible else None),
        "elements_by_assay": by_assay,
        "measurements_available": layer.counts(),
        "crispri_pairs_the_benchmark_calls_invalid": layer.crispri_invalid,
        # the facts, under two names because they answer two questions
        "facts_with_a_measured_counterpart": 2 * len(measured_rows),
        "experimental_facts_added": len(measured_rows) + regulated_pairs_blocks(measured_rows),
        "experimental_element_blocks": len(measured_rows),
        "experimental_rule_blocks": regulated_pairs_blocks(measured_rows),
        "genes_named_by_a_measurement": sorted(measured_genes),
        # agreement, with the denominator that is not chosen by the prediction
        "agrees": agree,
        "disagrees": disagree,
        "crispri_predicted_gene_not_tested": not_tested,
        "agreement_rate_over_all_matched_elements": (
            round(agree["crispri"] / crispri_matched, 4) if crispri_matched else None
        ),
        "agreement_denominator_all_matched_elements": crispri_matched,
        "agreement_rate_where_the_predicted_gene_was_tested": (
            round(agree["crispri"] / tested, 4) if tested else None
        ),
        "agreement_denominator_predicted_gene_tested": tested,
        "crispri_pairs_regulated": regulated_pairs,
        "crispri_pairs_measured_as_not_regulated": negative_pairs,
    }


def regulated_pairs_blocks(measured_rows: list[dict[str, Any]]) -> int:
    """One experimental rule per (element, gene) a screen measured as regulated, deduplicated."""
    return sum(len(r["measured"].get("crispri", {}).get("genes_regulated", [])) for r in measured_rows)


def pool(censuses: list[dict[str, Any]]) -> dict[str, Any]:
    """The same census over several chromosomes, summed where summing is meaningful."""
    sums = defaultdict(int)
    by_assay: dict[str, int] = defaultdict(int)
    agree: dict[str, int] = defaultdict(int)
    disagree: dict[str, int] = defaultdict(int)
    available: dict[str, int] = defaultdict(int)
    eligible_by_assay: dict[str, int] = defaultdict(int)
    keys = (
        "compiled_elements",
        "elements_in_an_assay_footprint",
        "elements_no_assay_ever_covered",
        "elements_with_any_measurement",
        "facts_with_a_measured_counterpart",
        "experimental_facts_added",
        "experimental_element_blocks",
        "experimental_rule_blocks",
        "crispri_predicted_gene_not_tested",
        "crispri_pairs_regulated",
        "crispri_pairs_measured_as_not_regulated",
        "crispri_pairs_the_benchmark_calls_invalid",
    )
    for c in censuses:
        for k in keys:
            sums[k] += c.get(k) or 0
        for a in ASSAYS:
            by_assay[a] += c["elements_by_assay"].get(a, 0)
            agree[a] += c["agrees"].get(a, 0)
            disagree[a] += c["disagrees"].get(a, 0)
            available[a] += c["measurements_available"].get(a, 0)
            eligible_by_assay[a] += (c.get("eligibility") or {}).get("by_assay", {}).get(a, 0)
    n, m = sums["compiled_elements"], sums["elements_with_any_measurement"]
    e = sums["elements_in_an_assay_footprint"]
    tested = agree["crispri"] + disagree["crispri"]
    return {
        "chromosomes": len(censuses),
        **{k: sums[k] for k in keys},
        "eligibility": {
            "rule": ELIGIBILITY_RULE,
            "compiled_elements": n,
            "elements_in_an_assay_footprint": e,
            "elements_no_assay_ever_covered": sums["elements_no_assay_ever_covered"],
            "share_eligible": round(e / n, 5) if n else 0.0,
            "by_assay": dict(eligible_by_assay),
        },
        "coverage_elements_measured": round(m / n, 5) if n else 0.0,
        "raised_share_of_the_eligible": round(m / e, 5) if e else None,
        "elements_by_assay": dict(by_assay),
        "measurements_available": dict(available),
        "agrees": dict(agree),
        "disagrees": dict(disagree),
        "agreement_rate_over_all_matched_elements": (
            round(agree["crispri"] / by_assay["crispri"], 4) if by_assay["crispri"] else None
        ),
        "agreement_denominator_all_matched_elements": by_assay["crispri"],
        "agreement_rate_where_the_predicted_gene_was_tested": (
            round(agree["crispri"] / tested, 4) if tested else None
        ),
        "agreement_denominator_predicted_gene_tested": tested,
    }


def sensitivity(
    elements: list[dict[str, Any]], layer: Layer, values: tuple[float, ...] = OVERLAP_SENSITIVITY
) -> dict[str, Any]:
    """The matches the overlap rule makes at each value, and the containments it refuses.

    `containment_not_upgraded` is the number of elements some measured interval swallows without
    meeting the rule: a 2 kb VISTA sequence around a 300 bp element measures the sequence, not the
    element, and this lane does not raise a fact on it. The number is reported so the refusal is
    visible.
    """
    out: dict[str, Any] = {"rule": "reciprocal overlap", "used": RECIPROCAL_OVERLAP, "at": {}}
    for f in values:
        per = {a: 0 for a in ASSAYS}
        any_assay = 0
        for e in elements:
            m = layer.for_element(e["start"], e["end"], f)
            for a in m:
                per[a] += 1
            any_assay += bool(m)
        out["at"][f"{f}"] = {"any": any_assay, **per}
    swallowed = {a: 0 for a in ASSAYS}
    for e in elements:
        s, t = e["start"], e["end"]
        for a in ASSAYS:
            items = layer.near(a, s, t)
            if any(contains(s, t, x.start, x.end) for x in items) and not any(
                measures(s, t, x.start, x.end, RECIPROCAL_OVERLAP) for x in items
            ):
                swallowed[a] += 1
    out["containment_not_upgraded"] = swallowed
    return out


# --- the program ----------------------------------------------------------------------------------
def basis_text(row: dict[str, Any]) -> str:
    """What the assays measured, in words, with every measured negative named."""
    parts = []
    c = row["measured"].get("crispri")
    if c:
        bits = [f"CRISPRi in {', '.join(c['cells'])} tested {len(c['genes_tested'])} genes"]
        if c["genes_regulated"]:
            bits.append("regulated " + ", ".join(c["genes_regulated"]))
        if c["genes_not_regulated"]:
            bits.append("measured no effect on " + ", ".join(c["genes_not_regulated"]))
        parts.append(", ".join(bits))
    m = row["measured"].get("lentimpra")
    if m:
        act = ", ".join(f"{k} {v:+.2f}" for k, v in m["activity"].items())
        parts.append(
            f"lentiMPRA log2(RNA/DNA) {act}, active in {len(m['cells_active'])} of "
            f"{len(m['activity'])} cells at {MPRA_ACTIVE}"
        )
    v = row["measured"].get("vista")
    if v:
        parts.append(
            f"VISTA {len(v['positive'])} positive, {len(v['negative'])} negative"
            + (f" in {', '.join(v['tissues'])}" if v["tissues"] else "")
        )
    a = row["agreement"]
    n_a, n_d = a["assays_agreeing"], a["assays_disagreeing"]
    parts.append(
        f"against the predicted target {row['predicted_gene'] or 'none'}, "
        f"{n_a} assay{'' if n_a == 1 else 's'} agrees and {n_d} disagrees ({a['verdict']})"
    )
    parts.append(f"reciprocal overlap rule {RECIPROCAL_OVERLAP}")
    return "; ".join(parts)


def rule_lines(row: dict[str, Any]) -> list[tuple[str, str, float, str]]:
    """(gene, action, strength, cell) for every element-gene link a screen measured as regulated.

    A measured negative gets no rule: a rule would state a relation the assay says is not there. It
    is kept as an experimental fact on the measured element block instead, named gene by gene.
    """
    c = row["measured"].get("crispri")
    if not c:
        return []
    out = []
    for gene in c["genes_regulated"]:
        hit = [p for p in c["pairs"] if p["gene"] == gene and p["regulated"]]
        strongest = max(hit, key=lambda p: abs(p["effect_size"]))
        # the screen silences the element: a gene that falls was being activated by it
        action = "activates" if strongest["effect_size"] < 0 else "inhibits"
        out.append((gene, action, round(min(1.0, abs(strongest["effect_size"])), 3), strongest["cell"]))
    return out
