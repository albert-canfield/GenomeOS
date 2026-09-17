# SPDX-License-Identifier: AGPL-3.0-or-later
"""The experimental layer of the compiled genome: what an assay measured over a compiled element.

The compiled per-chromosome programs state 940,803 facts and not one of them is `experimental`:
a region's role comes from its budget tier (`inferred`) and an element's target from a deletion in
AlphaGenome (`predicted`). Meanwhile the project holds real measurements over some of that same
sequence. This module attaches them, element by element, and counts honestly how thin the layer is.

Four assays, each kept under its own name:

- **CRISPRi** (ENCODE enhancer-gene benchmark, Gschwind et al. 2025): an element was silenced in its
  own chromosome and nearby genes were measured. A pair is `regulated` or not, and **a pair measured
  as not regulated is evidence, not an absence of evidence**: it is carried through to the program as
  an experimental fact, never dropped and never turned back into UNKNOWN.
- **lentiMPRA** (ENCODE4 joint library ENCSR106SZM, K562 / HepG2 / WTC11): how much a 200 bp sequence
  drives transcription from a reporter. Episomal, so it measures the sequence and not the locus.
- **VISTA** (LBNL, transgenic mouse e11.5): whether a sequence is an enhancer in a living embryo,
  positive or negative.
- **saturation mutagenesis** (Kircher et al. 2019, GSE126550, read through `knowledge/satmut.py`):
  nearly every single-base substitution of 21 regulatory elements, in an MPRA. It is different in
  kind from the other three, which ask whether an ELEMENT does something: this one asks which BASES
  inside it matter, and the two questions do not have the same answers.

**A base-level assay against an element-level prediction, and what may be concluded from it.** The
compiled claim is that deleting the element moves a gene. Saturation mutagenesis never deletes the
element and never measures that gene, so its verdicts are asymmetric on purpose:

- *agrees* — at least one measured base inside the element is **functional** (some substitution at it
  is significant at the portal's own defaults). The element demonstrably contains bases whose identity
  changes activity, which is measured support of the same weak kind lentiMPRA's "active" is: the
  sequence does something.
- *it cannot disagree*, and `disagrees["satmut"]` is therefore **always zero, by construction rather
  than by result** (`SATMUT_CANNOT_DISAGREE` says so in the census, so the zero cannot be read as
  "never contradicted"). An element every one of whose measured bases is inert is recorded under its
  own name, `bases_measured_none_functional` — **an element whose bases mostly do not matter is not an
  element that does nothing**. Single-base substitution cannot see a function carried redundantly
  across a site, and a null is depth-dependent besides (the functional share of the same 21 elements
  runs from 1% to 77% with barcode depth), so a base-level null is evidence about THOSE BASES and
  about nothing larger. It keeps `experimental` as its evidence kind and is stated in the measured
  block by name, exactly as a CRISPRi negative is.
- `bases_not_measured` is the fourth outcome and the never-looked one: the experiment's interval met
  the overlap rule but covers none of the element's own bases.

No `rule` comes from it either: the experiment is named for a gene by its authors, and that name is
curation, not a measurement of a target, so it is carried as provenance and never as a target.

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
from functools import lru_cache
from pathlib import Path
from typing import Any

from genomeos.attribution import mpra, vista
from genomeos.attribution.crispri import KNOWLEDGE as CRISPRI_KNOWLEDGE
from genomeos.knowledge import satmut as satmut_knowledge
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
    "satmut": (
        "saturation mutagenesis MPRA, Kircher et al. 2019 (GSE126550), single-base substitutions of "
        "21 regulatory elements, GRCh38"
    ),
}
ASSAYS = tuple(SOURCES)
BASE_LEVEL = ("satmut",)  # assays that measure bases inside an element rather than the element
SATMUT_CANNOT_DISAGREE = (
    "zero by construction, not by result: saturation mutagenesis substitutes one base at a time in a "
    "reporter and never deletes the element or measures the predicted gene, so it can support the "
    "compiled claim and cannot contradict it. An element whose measured bases are all inert is "
    "counted under bases_measured_none_functional, which is evidence about those bases only"
)


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


@dataclass(frozen=True, slots=True)
class SatmutElement:
    """One saturation-mutagenesis experiment: nearly every substitution of one regulatory element.

    `bases` is the position table of `knowledge/satmut.py`, 1-based genome position to whether some
    substitution there was significant (`functional`), whether it also moved activity by the strong
    bar, and the largest effect. The table is kept per base rather than summarised per element because
    the whole point of this assay is that it answers a question about bases: a compiled element that
    meets the overlap rule usually covers only part of the tested interval, and what is said about it
    must be said over the bases inside it and not over the experiment as a whole.
    """

    chrom: str
    start: int  # 0-based half-open, from the measured positions themselves
    end: int
    experiment: str  # the primary experiment of the locus, as satmut.loci chooses it
    repeats: int  # experiments of the same locus, this one included; the others are not counted twice
    bases: dict[int, dict]


@lru_cache(maxsize=2)
def _satmut_primaries(path: Path) -> tuple[SatmutElement, ...]:
    """Every locus's primary experiment, parsed once. A missing cache is no measurement, not a fetch.

    `knowledge.satmut.load` downloads when its file is absent; this checks first, because the layer
    must be readable on a machine that holds no satmut cache and must then say it measured nothing.
    Several experiments of one locus (TERT has four) are the same DNA measured again, so only the
    primary is carried and `repeats` records the rest, exactly as satmut's own pooled statistics do.
    """
    if not path.exists():
        return ()
    experiments = satmut_knowledge.by_element(satmut_knowledge.load(path))
    out: list[SatmutElement] = []
    for lead, group in satmut_knowledge.loci(experiments).items():
        table = satmut_knowledge.base_table(experiments[lead])
        if not table:
            continue
        out.append(
            SatmutElement(
                chrom=experiments[lead][0]["chrom"],
                start=min(table) - 1,
                end=max(table),
                experiment=lead,
                repeats=len(group),
                bases=table,
            )
        )
    return tuple(sorted(out, key=lambda e: (e.chrom, e.start)))


def load_satmut(chrom: str, path: Path = satmut_knowledge.DATA_PATH) -> list[SatmutElement]:
    """The saturation-mutagenesis experiments of one chromosome from the local cache, without fetching."""
    return [e for e in _satmut_primaries(path) if e.chrom == chrom]


# --- the layer ------------------------------------------------------------------------------------
@dataclass
class Layer:
    """Every measurement over one chromosome, each assay under its own name and searchable by overlap."""

    chrom: str
    crispri: list[CrispriPair] = field(default_factory=list)
    lentimpra: list[mpra.Element] = field(default_factory=list)
    vista: list[vista.VistaElement] = field(default_factory=list)
    satmut: list[SatmutElement] = field(default_factory=list)
    crispri_invalid: int = 0
    _starts: dict[str, list[int]] = field(default_factory=dict, repr=False)

    def __post_init__(self) -> None:
        self.crispri = sorted(self.crispri, key=lambda p: p.start)
        self.lentimpra = sorted(self.lentimpra, key=lambda e: e.start)
        self.vista = sorted(self.vista, key=lambda e: e.start)
        self.satmut = sorted(self.satmut, key=lambda e: e.start)
        self._starts = {name: [x.start for x in getattr(self, name)] for name in ASSAYS}

    @classmethod
    def load(cls, chrom: str) -> Layer:
        pairs, invalid = load_crispri(chrom)
        return cls(
            chrom=chrom,
            crispri=pairs,
            lentimpra=load_mpra(chrom),
            vista=load_vista(chrom),
            satmut=load_satmut(chrom),
            crispri_invalid=invalid,
        )

    @property
    def empty(self) -> bool:
        return not (self.crispri or self.lentimpra or self.vista or self.satmut)

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

        sm = [
            (e, reciprocal_overlap(start, end, e.start, e.end))
            for e in self.near("satmut", start, end, reach)
        ]
        sm = [(e, f) for e, f in sm if f >= fraction]
        if sm:
            out["satmut"] = self._satmut_of(start, end, sm)
        return out

    @staticmethod
    def _satmut_of(start: int, end: int, hits: list[tuple[SatmutElement, float]]) -> dict[str, Any]:
        """What a base-level assay says about one element: its own bases, never the whole experiment.

        Only the positions inside (start, end) are counted. `bases_measured` and `bases_functional` are
        two different zeros and are named apart for that reason: the first can be zero because the
        overlap fell on the unmeasured flank of the experiment, the second because every base that was
        looked at turned out not to matter.
        """
        inside: dict[int, dict] = {}
        for e, _ in hits:
            for pos, b in e.bases.items():
                if start <= pos - 1 < end:
                    inside.setdefault(pos, b)
        functional = [b for b in inside.values() if b["functional"]]
        strong = [b for b in inside.values() if b["strong"]]
        n = len(inside)
        return {
            "experiments": sorted(e.experiment for e, _ in hits),
            "repeat_experiments": sum(e.repeats - 1 for e, _ in hits),
            "bases_in_element": end - start,
            "bases_measured": n,
            "bases_functional": len(functional),
            "bases_strong": len(strong),
            "share_of_the_element_measured": round(n / (end - start), 4) if end > start else None,
            "share_functional_of_measured": round(len(functional) / n, 4) if n else None,
            "share_strong_of_measured": round(len(strong) / n, 4) if n else None,
            "strongest_effect": (
                round(max((abs(b["effect"]) for b in functional), default=0.0), 4) if functional else None
            ),
            "overlap": round(max(f for _, f in hits), 3),
        }


# --- prediction against measurement ---------------------------------------------------------------
AGREES, DISAGREES, NOT_TESTED = "agrees", "disagrees", "predicted_gene_not_tested"
# the two outcomes a base-level assay has that an element-level one does not, neither of them a verdict
# on the element: the first is measured-and-inert, the second is never-looked
BASES_INERT = "bases_measured_none_functional"
BASES_NOT_MEASURED = "bases_not_measured"


def agreement(predicted_gene: str, measured: dict[str, Any]) -> dict[str, Any]:
    """What each assay says about the compiled claim, assay by assay and then pooled.

    The compiled claim of an element block is that deleting it moves `predicted_gene`. CRISPRi asks
    exactly that question and can answer it three ways. lentiMPRA and VISTA ask the weaker question
    of whether the sequence acts at all, so their verdicts are recorded under their own names and the
    pooled verdict says how many assays fell each way rather than hiding the mixture.

    Saturation mutagenesis asks a fourth question - which bases inside the element matter - and only
    one of its answers is a verdict on the element. A functional base is support; an element whose
    measured bases are all inert is `BASES_INERT`, which is a measurement of those bases and not a
    contradiction of the element, and it is therefore kept out of `assays_disagreeing` rather than
    folded into it. See `SATMUT_CANNOT_DISAGREE`.
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
    s = measured.get("satmut")
    if s:
        if not s["bases_measured"]:
            out["satmut"] = BASES_NOT_MEASURED
            out["satmut_detail"] = (
                "the experiment met the overlap rule but measured none of this element's own bases"
            )
        elif s["bases_functional"]:
            out["satmut"] = AGREES
            out["satmut_detail"] = (
                f"{s['bases_functional']} of {s['bases_measured']} measured bases are functional "
                f"({s['bases_strong']} strongly): the element contains bases whose identity changes "
                "activity in a reporter"
            )
        else:
            out["satmut"] = BASES_INERT
            out["satmut_detail"] = (
                f"all {s['bases_measured']} measured bases are inert, which is a measurement of those "
                "bases; it is not a measurement of the element and not a disagreement"
            )
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
    """The strongest assay present decides: a perturbation or an embryo outranks a reporter.

    Saturation mutagenesis is a reporter assay - the element is read out episomally, one substitution
    at a time - so it lands with lentiMPRA and not with CRISPRi, however fine its resolution.
    """
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
        # the base-level assay, under its own name: its zero in `disagrees` is a property of the assay
        "satmut": satmut_counts(measured_rows, eligible),
    }


def satmut_counts(measured_rows: list[dict[str, Any]], eligible: dict[str, Any]) -> dict[str, Any]:
    """The base-level assay's own counters, kept apart from the element-level ones.

    Four numbers that must not be added together: how many elements a base-level measurement was made
    over; how many of those hold a functional base; how many were measured base by base and found
    inert (measured-and-absent); and how many were inside the footprint but never raised at all
    (never-looked). The last is the difference between the eligibility count and the raised one, and it
    is named here so that neither zero can be mistaken for the other.
    """
    verdicts = [r["agreement"].get("satmut") for r in measured_rows if "satmut" in r["measured"]]
    sat = [r["measured"]["satmut"] for r in measured_rows if "satmut" in r["measured"]]
    footprint = eligible["by_assay"].get("satmut", 0)
    return {
        "elements_measured_base_by_base": len(sat),
        "elements_with_a_functional_base": sum(1 for v in verdicts if v == AGREES),
        "elements_measured_and_every_base_inert": sum(1 for v in verdicts if v == BASES_INERT),
        "elements_matched_but_no_base_of_theirs_measured": sum(
            1 for v in verdicts if v == BASES_NOT_MEASURED
        ),
        "elements_in_the_footprint_never_raised": footprint - len(sat),
        "bases_measured": sum(s["bases_measured"] for s in sat),
        "bases_functional": sum(s["bases_functional"] for s in sat),
        "bases_strong": sum(s["bases_strong"] for s in sat),
        "experiments_matched": sorted({x for s in sat for x in s["experiments"]}),
        "why_it_can_never_disagree": SATMUT_CANNOT_DISAGREE,
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
    satmut_keys = (
        "elements_measured_base_by_base",
        "elements_with_a_functional_base",
        "elements_measured_and_every_base_inert",
        "elements_matched_but_no_base_of_theirs_measured",
        "elements_in_the_footprint_never_raised",
        "bases_measured",
        "bases_functional",
        "bases_strong",
    )
    satmut: dict[str, Any] = dict.fromkeys(satmut_keys, 0)
    experiments: set[str] = set()
    for c in censuses:
        for k in keys:
            sums[k] += c.get(k) or 0
        s = c.get("satmut") or {}
        for k in satmut_keys:
            satmut[k] += s.get(k) or 0
        experiments |= set(s.get("experiments_matched") or ())
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
        "satmut": {
            **satmut,
            "experiments_matched": sorted(experiments),
            "why_it_can_never_disagree": SATMUT_CANNOT_DISAGREE,
        },
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
    s = row["measured"].get("satmut")
    if s:
        # the negative case is named as loudly as the positive one, and neither is stated as a verdict
        # on the element: this assay measures bases
        parts.append(
            f"saturation mutagenesis ({', '.join(s['experiments'])}) measured "
            f"{s['bases_measured']} of the element's {s['bases_in_element']} bases, "
            f"{s['bases_functional']} of them functional and {s['bases_strong']} strongly so"
            + (
                ", so no base measured here matters and that is evidence about those bases only"
                if s["bases_measured"] and not s["bases_functional"]
                else ""
            )
            + (
                " - the experiment overlaps this element but measured none of its bases, which is a "
                "gap in the assay and not a finding about the element"
                if not s["bases_measured"]
                else ""
            )
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

    Saturation mutagenesis makes no rule either, and for a different reason: it measures a reporter's
    activity, not a gene's, so it names no relation at all. Its experiments carry gene names given by
    their authors (SORT1, IRF4), which is curation and not a measured target.
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
