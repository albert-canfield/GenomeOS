# SPDX-License-Identifier: AGPL-3.0-or-later
"""The `not_open_in_reader` rules, typed and located: what kind of candidate defect each one is.

`genomeos.attribution.context_evidence` gave every compiled rule one of three states and the census
counted them: of 440,589 compiled rules, 81,635 could be read at all and **55,084 of those assert
their element in a cell where the reader does not detect it open** (`not_open_in_reader`). That is one
number over a very mixed population. This module takes the same rules, in the same order, and carries
four properties beside each one so the number becomes a list a person can work through:

    cell            the rule's own `when: cell_type` label, with the reader biosample it maps to and
                    whether the base-rate comparison found that assignment informative at all
    element_class   `class:` as the compiler already writes it (`compile.derived_class`), not a new
                    one, with the compiler's own `activity:` and `origin:` axis values beside it,
                    because `class:` collapses to one value over the attributed elements
    effect_band     `strength` as the AlphaGenome deletion prediction already carries it
                    (`predict.enhancer_target`: strong at or above STRONG_EFFECT, weak below)
    distance_band   the element-to-gene distance in the bands `attribution.executor` already uses

**Nothing here is recomputed and no cut-off is introduced.** The state comes from
`context_evidence.state_for` unchanged, the mapping table from `context_evidence.mapping`, the effect
band from the cached prediction's own `strength` field, the distance bands from `executor._band`, and
the informativeness of a cell's assignment from `context_evidence_baserate.json`'s own `difference`
against the tolerance that result fixed before any share was computed. A test holds this module's rule
enumeration to `context_evidence.rule_loci` locus for locus, so the two cannot answer differently
about which rules are being described.

**What this module is.** Descriptive and diagnostic. It locates and types candidate defects; it
establishes about no rule that the rule is wrong. The limitations it inherits are carried into every
result it writes, unchanged in force:

* `not_open_in_reader` means *not detected open at the reader's registered call*, never *closed*
  (`context_evidence.NOT_CLOSED`). A narrowPeak set is a call set and absence is absence of a call.
* reader v1's openness comes from ENCODE peaks and AlphaGenome was trained on ENCODE, so this whole
  axis is a **consistency check between two readings of the same chromatin**, never independent
  evidence and never validation (`context_evidence.NOT_VALIDATION`). The 1.81x base-rate figure of
  `context_evidence_baserate.json` is descriptive and not independent for the same reason, and is
  quoted here with that label or not at all.

**Move no verdict, delete no rule, change no compiled label.** This module writes a list beside the
compiled programs. It imports `compile` and `context_evidence` read-only and changes neither.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

from genomeos.attribution import context_evidence as ce
from genomeos.attribution import executor as ex
from genomeos.attribution.target_calibration import tss_distance
from genomeos.predict.enhancer_target import STRONG_EFFECT
from genomeos.results import RESULTS_DIR, load_result

# ---- where each rule comes from ----------------------------------------------------------------

#: A rule emitted beside an attributed element: the AlphaGenome deletion prediction is the rule's own
#: evidence, so the rule has a predicted effect of its own.
SOURCE_PREDICTED = "predicted_alphagenome_deletion"
#: A rule emitted in the experimental layer: a CRISPRi screen measured it, and its strength is the
#: screen's effect size, not AlphaGenome's. These rules carry **no** AlphaGenome effect of their own
#: and are counted as their own population in the effect-size breakdown rather than given a band.
SOURCE_MEASURED = "measured_crispri"
SOURCES = (SOURCE_PREDICTED, SOURCE_MEASURED)

# ---- the effect-size band, as the prediction already carries it ---------------------------------

STRONG, WEAK = "strong", "weak"
EFFECT_BANDS = (STRONG, WEAK)
EFFECT_BAND_CALL = (
    f"`strength` as genomeos.predict.enhancer_target already writes it on the cached prediction: "
    f"strong at or above STRONG_EFFECT = {STRONG_EFFECT} in |log2 fold change| of the predicted "
    f"expression move on deleting the element, weak below. It is an effect-size band and never a "
    f"certainty; this module reads the field and does not recompute it."
)
NO_ALPHAGENOME_EFFECT = (
    "a rule of the experimental layer: a CRISPRi screen measured it, so it carries the screen's "
    "effect size and no AlphaGenome deletion effect of its own. It is counted here and given no band."
)

# ---- the distance bands, as the project already uses them ---------------------------------------

#: `genomeos.attribution.executor` bands an element-to-gene distance in three, and its own
#: WIDE_MAX_DISTANCE is where the project stops banding ("inside the 1 Mb prediction window, as in E2
#: and E3"). Both are imported rather than restated, so this module cannot band differently.
INNER_BAND_OF = ex._band
OUTER_LIMIT = ex.WIDE_MAX_DISTANCE
OUTER_BAND = f"beyond {OUTER_LIMIT // 1000} kb"
DISTANCE_BANDS = (
    INNER_BAND_OF(0),
    INNER_BAND_OF(OUTER_LIMIT // 50),
    INNER_BAND_OF(OUTER_LIMIT // 5),
    OUTER_BAND,
)
DISTANCE_CALL = (
    "the distance from the element's midpoint to the target gene's GENCODE v50 TSS, by "
    "genomeos.attribution.target_calibration.tss_distance (floored at its MIN_DISTANCE), banded by "
    "genomeos.attribution.executor._band unchanged. That function bands up to executor."
    f"WIDE_MAX_DISTANCE = {OUTER_LIMIT}, which is where the project's own banding stops, so "
    f"everything above it is reported as one residual row, {OUTER_BAND!r}, and not as a new band."
)
#: Half the 1 Mb window the deletion scorer reads. A symbol whose TSS resolves further away than this
#: cannot be the gene the scorer moved, so the distance is a symbol-resolution artefact rather than a
#: real element-to-gene distance. Counted and named; no rule is dropped for it.
SCORER_HALF_WINDOW = 500_000
BEYOND_THE_WINDOW = (
    "the distance from the element to the TSS this gene symbol resolves to exceeds half the 1 Mb "
    "window the AlphaGenome deletion scorer reads, so the scorer cannot have read that TSS: the "
    "symbol resolved to another locus of the same symbol. This is a fault in the distance, not in "
    "the rule, and the rule is counted everywhere else unchanged."
)
NO_DISTANCE = "the target gene symbol carries no GENCODE v50 TSS on this chromosome, so it has no distance"

# ---- whether the assigned cell's assignment carries information at all --------------------------

BASERATE_RESULT = "context_evidence_baserate"
#: The three assigned cells lane-context2 singled out in docs/ATTRIBUTION.md: SK-N-SH, where the
#: difference is inside the registered tolerance and the assignment carries no information about
#: openness, and astrocyte and ovary, "both only just outside it". Named here as a citation; every
#: number about them is read from the base-rate result, never restated.
CELLS_LANE_CONTEXT2_SINGLED_OUT = ("SK-N-SH", "astrocyte", "ovary")
INFORMATIVE_CALL = (
    "a cell's assignment clears the registered tolerance when, over the rules assigned to it, the "
    "share of them open in that cell exceeds the mean share open in the other twelve reader "
    "biosamples by more than the tolerance context_evidence_baserate.json fixed before any share was "
    "computed. Both the difference and the tolerance are read from that result. This says the "
    "assignment is informative about where the element is open; it says nothing about whether the "
    "rule is right, and the base-rate comparison is itself descriptive and not independent."
)


@dataclass(frozen=True)
class Informativeness:
    """Per assigned cell, what the base-rate comparison found, read from its result and not recomputed."""

    tolerance: float
    difference: dict[str, float]
    clears: frozenset[str]
    source: str

    def of(self, biosample: str | None) -> float | None:
        return None if biosample is None else self.difference.get(biosample)

    def informative(self, biosample: str | None) -> bool:
        return biosample is not None and biosample in self.clears


def informativeness(results_dir: Path = RESULTS_DIR) -> Informativeness:
    """The assigned-cell differences and the tolerance, from `context_evidence_baserate.json`."""
    r = load_result(BASERATE_RESULT, results_dir)
    if not r:
        raise FileNotFoundError(
            f"no {BASERATE_RESULT} result under {results_dir}: the informativeness of an assigned "
            "cell is read from it and is never chosen here"
        )
    tol = float(r["tolerance"])
    diff = {cell: float(v["difference"]) for cell, v in r["per_assigned_cell"].items()}
    return Informativeness(
        tolerance=tol,
        difference=diff,
        clears=frozenset(c for c, d in diff.items() if d > tol),
        source=f"data/results/{BASERATE_RESULT}.json",
    )


# ---- one rule, with the four properties ---------------------------------------------------------


@dataclass(frozen=True)
class Rule:
    """One compiled rule with the properties this profile types it by. Nothing here is a verdict."""

    source: str
    element: str
    chrom: str
    start: int
    end: int
    cell: str
    gene: str
    element_class: str
    activity_axis: str
    origin_axis: str
    effect: float | None
    effect_band: str | None
    distance: int | None

    @property
    def locus(self) -> str:
        return f"{self.chrom}:{self.start}-{self.end}"

    @property
    def distance_band(self) -> str:
        return distance_band(self.distance)

    @property
    def beyond_the_scorer_window(self) -> bool:
        return self.distance is not None and self.distance > SCORER_HALF_WINDOW

    def row(self) -> dict[str, Any]:
        """The rule as a result records it: located, typed, and with no state of its own."""
        return {
            "chromosome": self.chrom,
            "element": self.element,
            "locus": self.locus,
            "gene": self.gene,
            "cell": self.cell,
            "element_class": self.element_class,
            "activity_axis": self.activity_axis,
            "origin_axis": self.origin_axis,
            "source": self.source,
            "alphagenome_log2_fold_change": self.effect,
            "effect_band": self.effect_band,
            "distance_to_tss_bp": self.distance,
            "distance_band": self.distance_band,
        }


def axis_value(axes: dict[str, list[list[str]]], name: str) -> str:
    """One compiled axis as the compiler writes it: `,` where every group holds, `|` where a group is
    unresolved alternatives. The same join `compile.axis_lines` puts after the axis name, so this reads
    the element's own axis text and does not restate it."""
    groups = axes.get(name) or []
    return ", ".join("|".join(g) for g in groups) or "unknown"


def distance_band(distance: int | None) -> str:
    """The project's own band for an element-to-gene distance, with one residual row above its top."""
    if distance is None:
        return NO_DISTANCE
    return OUTER_BAND if distance >= OUTER_LIMIT else INNER_BAND_OF(distance)


def rules(chrom: str, results_dir: Path = RESULTS_DIR, layer: Any = None) -> list[Rule]:
    """Every rule `compile.compile_chromosome` emits for this chromosome, in the order it emits them.

    This mirrors the compiler block for block, exactly as `context_evidence.rule_loci` does, and adds
    the gene, the class the compiler writes, the prediction's own effect and the distance to the TSS.
    `tests/test_not_open_profile.py` holds the loci of this list to `rule_loci`'s, element for
    element, so the two cannot describe different rules.
    """
    from genomeos.attribution import compile as cp
    from genomeos.attribution.measured import rule_links
    from genomeos.attribution.pilot_bio import gene_tss

    elements = cp._attributed(chrom, results_dir)
    ccre = cp._ccres(chrom, results_dir)
    rep = cp._Interspersed(chrom, results_dir)
    tss = gene_tss(chrom)
    out: list[Rule] = []
    sequence: dict[str, dict] = {}
    for e in elements:
        axes = cp.element_axes(e, ccre, rep)
        sequence[e["id"]] = axes
        pc = e["predicted_coding"]
        effect = float(pc["log2_fold_change"])
        out.append(
            Rule(
                source=SOURCE_PREDICTED,
                element=e["id"],
                chrom=chrom,
                start=e["start"],
                end=e["end"],
                cell=cp.context(pc.get("tissue")),
                gene=pc["gene"],
                element_class=cp.derived_class(axes),
                activity_axis=axis_value(axes, "activity"),
                origin_axis=axis_value(axes, "origin"),
                effect=effect,
                # the band the cached prediction already carries; `or WEAK` is the compiler's own
                # reading of a missing field (compile.py: `pc.get('strength') or 'weak'`)
                effect_band=pc.get("strength") or WEAK,
                distance=_distance(e["start"], e["end"], tss.get(pc["gene"])),
            )
        )
    _layer, measured_rows = cp._measured_rows(chrom, elements, results_dir, layer)
    for row in measured_rows:
        axes = cp.measured_axes(row, sequence.get(row["id"]))
        cls = cp.derived_class(axes)
        activity = axis_value(axes, "activity")
        origin = axis_value(axes, "origin")
        for gene, _action, _strength, cell, _split in rule_links(row):
            out.append(
                Rule(
                    source=SOURCE_MEASURED,
                    element=f"{row['id']}_measured",
                    chrom=chrom,
                    start=row["start"],
                    end=row["end"],
                    cell=cp.context(cell),
                    gene=gene,
                    element_class=cls,
                    activity_axis=activity,
                    origin_axis=origin,
                    effect=None,
                    effect_band=None,
                    distance=_distance(row["start"], row["end"], tss.get(gene)),
                )
            )
    return out


def _distance(start: int, end: int, tss: int | None) -> int | None:
    d = tss_distance((start + end) // 2, tss)
    return None if d is None else int(d)


def limitations() -> dict[str, str]:
    """The limitations this profile inherits, imported from the reading it describes, not restated."""
    return {
        "not_closed": ce.NOT_CLOSED,
        "not_validation": ce.NOT_VALIDATION,
        "base_rate_is_descriptive_and_not_independent": (
            "the 1.81x figure of context_evidence_baserate.json - the assigned cell open on 0.3252 of "
            "the 81,635 assessable rules against 0.1796 of the 979,455 rule-and-biosample pairs in "
            "the other twelve - is descriptive and not independent, for the same reason: both "
            "readings are of the same ENCODE chromatin the model was trained on. It is quoted with "
            "that label or not at all."
        ),
        "descriptive_only": (
            "this profile locates and types candidate defects. It establishes about no rule that the "
            "rule is wrong, no verdict is moved, no rule is deleted and no compiled label is changed."
        ),
        "informative_is_not_right": (
            "a cell whose assignment clears the base-rate tolerance is informative about where the "
            "element is open, not correct about the rule. A not_open_in_reader rule in such a cell is "
            "a sharper candidate than one in a cell whose assignment carries no information, and "
            "nothing more than that."
        ),
    }


def registration() -> dict[str, Any]:
    """Every call this profile makes, from the code itself, so a result cannot be read without them."""
    return {
        "sources_of_a_rule": {SOURCE_PREDICTED: "", SOURCE_MEASURED: NO_ALPHAGENOME_EFFECT},
        "effect_bands": list(EFFECT_BANDS),
        "effect_band_call": EFFECT_BAND_CALL,
        "strong_effect": STRONG_EFFECT,
        "distance_bands": list(DISTANCE_BANDS),
        "distance_call": DISTANCE_CALL,
        "beyond_the_scorer_window": BEYOND_THE_WINDOW,
        "scorer_half_window": SCORER_HALF_WINDOW,
        "element_class_call": (
            "`class:` as genomeos.attribution.compile.derived_class writes it on the compiled element, "
            "from the same axes the compiler builds: the first ENCODE registry role the element holds "
            "for certain, else unknown. No classification is invented here. Because every attributed "
            "element is a distal enhancer-like registry cCRE, that field takes one value over this "
            "population and separates nothing, so the compiler's own `activity:` and `origin:` axis "
            "values are reported beside it - also read off the compiled element, also not invented."
        ),
        "informative_call": INFORMATIVE_CALL,
        "cells_lane_context2_singled_out": list(CELLS_LANE_CONTEXT2_SINGLED_OUT),
        "state_call": (
            "genomeos.attribution.context_evidence.state_for, imported unchanged, with its own "
            "mapping table and reader v1's own openness call. No threshold, mapping or registered "
            "definition is changed, and the chromosome counts are checked against the census result."
        ),
        "limitations": limitations(),
        "alphagenome_requests": 0,
    }
