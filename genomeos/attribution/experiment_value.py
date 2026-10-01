# SPDX-License-Identifier: AGPL-3.0-or-later
"""The next experiment's information value: the questions a design must answer before a volume.

Review item 12's S8, adopted 2026-09-28: *a representative arm (how common a function is) and a
disagreement-selected arm (which model is right) kept apart, selection probabilities recorded;
independent loci, assay feasibility and distinguishable hypotheses before volume; the
second-cell-type arm stays inconclusive, and re-reading held-out data is not fresh validation.*

Three lanes in the week to 2026-09-29 each stopped on a question nobody had asked first: one fell
below a coverage floor after its scoring code was written, one found no matched alternative that
moved enough links, and one found its two readouts not comparable. In each case the answer was
available from counts alone, before any prediction was joined and before any volume was named. This
module is that order, made mechanical: gates first, each refusable with its reason, and a volume
only when every gate is met.

**What this module does not do.** It selects no experiment, reads no dataset, makes no model request
and names no price. It takes a design written down as records and either refuses it with a reason or
returns the two arms' selections and their volumes. Every dataset, assay, cell line and context
reaches it as a record carrying its own *registered reading*; no name and no number of any dataset
is written here, so a reading can never be strengthened by this code.

**The gates, in order.** A design is refused at the first gate it fails, and a refusal is never
answered with a sample size:

1. `eligibility` -- blinded. The unit of analysis named, the units counted, the label prevalence,
   and the number of independent loci by the convention below. No prediction is joined: the gate
   refuses outright if the units it is handed carry predictions, scores or hypotheses. This is the
   count five lanes lacked, and a design must be refusable here before any volume is reported.
2. `feasibility` -- per candidate, whether the assay is available at all; and, for the design as a
   whole, whether every candidate states a detection limit and declares whether its readout is
   comparable with the readout it will be set against. An unavailable assay drops one candidate; an
   undeclared detection limit or an undeclared comparability is a hole in the design, so it refuses
   the design rather than the candidate.
3. `distinguishability` -- per candidate, each hypothesis's predicted outcome in the assay's own
   units, the assay's detection limit, and whether the result would separate the hypotheses.
   Two hypotheses predicting the same outcome are indistinguishable whatever the assay can see; a
   separation below the detection limit is refused because the assay cannot report it.
4. `freshness` -- prior exposure travels. A candidate already read is development evidence and can
   never be reported as fresh validation, so a design claiming fresh validation is refused when any
   surviving candidate has been read. Each candidate's context carries its registered replication
   reading word for word, and only a context registered as `replicated` is counted as a
   replication: a context registered `inconclusive` stays inconclusive wherever it appears.

**Then, and only then, the two arms.** They answer different questions and are kept apart because
their selection probabilities come from different designs; `pool` refuses to count them together.
The representative arm draws loci uniformly, so its counts estimate how common a function is. The
disagreement arm takes the candidates where the hypotheses separate most, so it says which model is
right and is representative of nothing. Every selected candidate records the probability with which
it was selected.

**Volume comes last**, and is an allocation of the registered target over the eligible independent
loci. It is arithmetic, not a power calculation.

Every number returned carries the population it was counted over (`Count`), and every floor is the
design's own, registered with it: this module sets no floor value, so no threshold here can be moved
after an outcome is seen.
"""

from __future__ import annotations

import random
from collections.abc import Iterable, Sequence
from dataclasses import dataclass, field
from typing import Any

#: The seed for the representative arm's draw. Fixed so a selection can be reproduced and audited.
SEED = 20261001

# ---------------------------------------------------------------------------
# The registered convention for an independent locus
# ---------------------------------------------------------------------------

#: Two records on one chromosome no farther apart than this are one locus. The span is the order of a
#: mammalian topologically associating domain, which is why two elements that close are not two
#: independent draws; it is a round, pre-set bound, not a figure fitted to an outcome.
#:
#: The same rule is registered for the second-cell-type count in `genomeos/attribution/cell2.py`
#: (`INDEPENDENT_LOCUS_SPAN`, `INDEPENDENT_LOCUS_RULE`, 2026-10-01). The project states one
#: convention, and these two must not diverge: that module's count is stratified by cell type and
#: reads measured positives, this one groups whatever records a design offers, but the grouping rule
#: and the span are the same. Whether the two should share a single constant is the coordinator's to
#: decide; no import is made here, so that a generic design check does not depend on a
#: dataset-specific module.
LOCUS_WINDOW = 1_000_000

LOCUS_CONVENTION = (
    "Independent locus, the convention registered with this module and used by every count it "
    "returns: two records are the same locus if they name the same target gene, or if they lie on "
    "the same chromosome within "
    f"{LOCUS_WINDOW:,} bp of each other. The grouping is transitive, so a chain of records each "
    "within the window of the next is one locus however far the chain reaches. The reason for the "
    "rule is that both ways of being close make two observations carry the same information: two "
    "elements read out on one gene share that gene's measurement, and two elements inside one "
    "window share the local chromatin, the local variants and, in most assays, the same reference "
    "and mapping behaviour, so counting them as two would count one piece of evidence twice. "
    "This is an operational grouping for counting, chosen to be conservative. It is NOT established "
    "biological independence: elements in one window can act on different genes through different "
    "mechanisms, elements megabases apart or on different chromosomes can share a regulator, a "
    "compartment or a trans factor, and nothing here measures whether two loci are independent in "
    "any biological sense. A count of independent loci under this convention is an upper bound on "
    "how much independent evidence the design holds, never a demonstration of it."
)

# ---------------------------------------------------------------------------
# Names used in the records and in the result
# ---------------------------------------------------------------------------

#: The representative arm: how common a function is. Loci drawn uniformly, so counts generalise to
#: the eligible population and to nothing wider.
REPRESENTATIVE = "representative"

#: The disagreement-selected arm: which model is right. Candidates ranked by how far the hypotheses
#: separate relative to what the assay can see. Representative of no population.
DISAGREEMENT = "disagreement"

ARMS: tuple[str, ...] = (REPRESENTATIVE, DISAGREEMENT)

ARM_QUESTIONS: dict[str, str] = {
    REPRESENTATIVE: "how common a function is",
    DISAGREEMENT: "which model is right",
}

MERGE_REFUSED = (
    "The representative and disagreement-selected arms cannot be merged. Their selection "
    "probabilities come from different designs: the representative arm draws loci uniformly, so its "
    "counts estimate a prevalence in the eligible population, while the disagreement arm takes the "
    "largest separations, so its candidates are chosen by the very quantity under test. Pooling "
    "them would report a prevalence over a sample selected on the outcome. Count each arm "
    "separately and say which question it answers."
)

#: What a design claims the experiment will be. Carried into every refusal and every label.
DEVELOPMENT_EVIDENCE = "development evidence"
FRESH_VALIDATION = "fresh validation"
REPLICATION = "replication"
CLAIMS: tuple[str, ...] = (DEVELOPMENT_EVIDENCE, FRESH_VALIDATION, REPLICATION)

#: A context's registered replication reading. Only `REPLICATED` is ever counted as a replication.
REPLICATED = "replicated"
INCONCLUSIVE = "inconclusive"
NOT_ATTEMPTED = "not_attempted"
REPLICATION_STATUSES: tuple[str, ...] = (REPLICATED, INCONCLUSIVE, NOT_ATTEMPTED)

#: Gate names, in the order they run.
GATE_ELIGIBILITY = "eligibility"
GATE_FEASIBILITY = "feasibility"
GATE_DISTINGUISHABILITY = "distinguishability"
GATE_FRESHNESS = "freshness"
GATES: tuple[str, ...] = (
    GATE_ELIGIBILITY,
    GATE_FEASIBILITY,
    GATE_DISTINGUISHABILITY,
    GATE_FRESHNESS,
)

#: Attributes that must not reach the blinded eligibility count. A unit carrying any of them has had
#: a prediction or a model score joined to it, which is the failure the gate exists to prevent. The
#: list holds only fields that carry what a model or a hypothesis said; a field describing the
#: measurement, such as the assay or its detection limit, is not unblinding and is not listed, so the
#: refusal's reason stays true to what was found.
UNBLINDING_FIELDS: tuple[str, ...] = (
    "hypotheses",
    "hypothesis",
    "predictions",
    "prediction",
    "predicted_outcome",
    "predicted",
    "score",
    "scores",
    "expected_effect",
)

VOLUME_BASIS = (
    "The registered target number of loci allocated over the eligible independent loci, times the "
    "registered units per locus. It is arithmetic over the counts above, not a power calculation: "
    "it states how much the design asks for, not what effect the experiment could detect. What the "
    "assay could resolve is the per-candidate separation and detection limit recorded beside it."
)


# ---------------------------------------------------------------------------
# Records
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Count:
    """A number and the population it was counted over. There is no other way to report a number."""

    value: float
    population: str

    def __post_init__(self) -> None:
        if not self.population or not self.population.strip():
            raise ValueError("a count must name the population it was counted over")

    def as_dict(self) -> dict[str, Any]:
        return {"value": self.value, "population": self.population}

    def __str__(self) -> str:
        return f"{self.value} ({self.population})"


@dataclass(frozen=True)
class EligibilityUnit:
    """One unit of the blinded eligibility count.

    It has no field for a prediction, a score or a hypothesis, and that is the point: the gate 1
    count is taken over these records, so it cannot be computed from anything a model said.
    `label` is the design's own outcome label, used for prevalence only.
    """

    unit_id: str
    chromosome: str
    position: int
    target_gene: str
    label: bool

    @property
    def ident(self) -> str:
        return self.unit_id


@dataclass(frozen=True)
class Assay:
    """What the measurement can see, as the design registered it. No assay is named in this module.

    `detection_limit` is in `units` and is the smallest difference the assay reports as a
    difference. `readout_comparable` must be declared: `None` means the design has not said whether
    this readout can be set against the readout it will be compared with, which is a hole in the
    design and not a property of the candidate.
    """

    name: str
    units: str
    detection_limit: float | None
    available: bool
    availability_reading: str
    readout_comparable: bool | None = None
    comparability_reading: str = ""


@dataclass(frozen=True)
class Hypothesis:
    """One hypothesis and the outcome it predicts, in the assay's own units."""

    name: str
    predicted_outcome: float
    units: str


@dataclass(frozen=True)
class Context:
    """The cell type, tissue or setting a candidate sits in, with its registered reading.

    `registered_reading` is carried into the result word for word and is never rewritten here.
    `replication_status` is the registered reading's own verdict; `counts_as_replication` is derived
    from it, so a caller cannot promote an inconclusive context by asserting that it replicates.
    """

    name: str
    replication_status: str
    registered_reading: str

    def __post_init__(self) -> None:
        if self.replication_status not in REPLICATION_STATUSES:
            raise ValueError(
                f"a context's replication status must be one of {REPLICATION_STATUSES}, "
                f"not {self.replication_status!r}"
            )
        if not self.registered_reading.strip():
            raise ValueError("a context must carry its registered reading")

    @property
    def counts_as_replication(self) -> bool:
        return self.replication_status == REPLICATED


@dataclass(frozen=True)
class Candidate:
    """One thing the experiment could measure, with its hypotheses, its assay and its history.

    `prior_reads` holds the registered reading of every earlier read of this candidate. A non-empty
    `prior_reads` makes the candidate development evidence for good; no later use can restore it to
    fresh validation.
    """

    candidate_id: str
    chromosome: str
    position: int
    target_gene: str
    assay: Assay
    hypotheses: tuple[Hypothesis, ...]
    context: Context
    prior_reads: tuple[str, ...] = ()

    @property
    def ident(self) -> str:
        return self.candidate_id

    @property
    def was_read(self) -> bool:
        return bool(self.prior_reads)

    @property
    def evidence_role(self) -> str:
        """Development evidence if it has ever been read. Never fresh validation in that case."""
        return DEVELOPMENT_EVIDENCE if self.was_read else FRESH_VALIDATION


@dataclass(frozen=True)
class Floors:
    """The design's own registered floors. This module ships no values for them.

    Every floor is the design's, registered with it and echoed into the result, so that a floor
    moved after an outcome was seen is visible in the diff of the design rather than hidden in this
    module. `registered_at` is whatever the design uses to say when the floors were fixed.
    """

    min_units: int
    min_positive_units: int
    min_label_prevalence: float
    min_independent_loci: int
    min_distinguishing_candidates: int
    target_loci: int
    units_per_locus: int
    registered_at: str

    def as_dict(self) -> dict[str, Any]:
        return {
            "min_units": self.min_units,
            "min_positive_units": self.min_positive_units,
            "min_label_prevalence": self.min_label_prevalence,
            "min_independent_loci": self.min_independent_loci,
            "min_distinguishing_candidates": self.min_distinguishing_candidates,
            "target_loci": self.target_loci,
            "units_per_locus": self.units_per_locus,
            "registered_at": self.registered_at,
        }


@dataclass(frozen=True)
class Design:
    """A proposed experiment, as records. Nothing here is read from a dataset by this module."""

    name: str
    unit_of_analysis: str
    claim: str
    floors: Floors
    eligibility: Sequence[Any] = ()
    candidates: Sequence[Candidate] = ()
    source_reading: str = ""

    def __post_init__(self) -> None:
        if self.claim not in CLAIMS:
            raise ValueError(f"a design's claim must be one of {CLAIMS}, not {self.claim!r}")


@dataclass(frozen=True)
class GateResult:
    """One gate's verdict. A refusal always carries its reason; a pass states what it read."""

    gate: str
    passed: bool
    reason: str
    counts: tuple[Count, ...] = ()
    detail: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not self.reason.strip():
            raise ValueError(f"gate {self.gate!r} must state its reason")

    def as_dict(self) -> dict[str, Any]:
        return {
            "gate": self.gate,
            "passed": self.passed,
            "reason": self.reason,
            "counts": [c.as_dict() for c in self.counts],
            "detail": self.detail,
        }


@dataclass(frozen=True)
class Selection:
    """One selected candidate and the probability with which it was selected.

    The probability is required. A selection cannot be recorded without it, because a selected
    candidate whose inclusion probability is unknown cannot be weighted, and a count over such
    candidates estimates nothing.
    """

    candidate_id: str
    arm: str
    selection_probability: float
    scheme: str
    representative_of_population: bool

    def __post_init__(self) -> None:
        if self.arm not in ARMS:
            raise ValueError(f"an arm must be one of {ARMS}, not {self.arm!r}")
        p = self.selection_probability
        if p is None or not isinstance(p, (int, float)) or isinstance(p, bool):
            raise ValueError("a selection cannot be recorded without its selection probability")
        if not (0.0 < float(p) <= 1.0):
            raise ValueError(
                "a selection probability must be greater than 0 and at most 1, "
                f"not {self.selection_probability!r}"
            )
        if not self.scheme.strip():
            raise ValueError("a selection must name the scheme that selected it")

    def as_dict(self) -> dict[str, Any]:
        return {
            "candidate_id": self.candidate_id,
            "arm": self.arm,
            "selection_probability": self.selection_probability,
            "scheme": self.scheme,
            "representative_of_population": self.representative_of_population,
        }


@dataclass(frozen=True)
class InformationValue:
    """What one candidate would tell the design, and whether the assay could tell it."""

    candidate_id: str
    units: str
    predicted: tuple[tuple[str, float], ...]
    detection_limit: float | None
    separation: float | None
    separates: bool
    refusal: str
    evidence_role: str
    context_name: str
    context_reading: str
    counts_as_replication: bool
    prior_reads: tuple[str, ...] = ()

    @property
    def may_be_reported_as_fresh_validation(self) -> bool:
        return self.evidence_role == FRESH_VALIDATION and not self.prior_reads

    def as_dict(self) -> dict[str, Any]:
        return {
            "candidate_id": self.candidate_id,
            "units": self.units,
            "predicted": [{"hypothesis": n, "predicted_outcome": v} for n, v in self.predicted],
            "detection_limit": self.detection_limit,
            "separation": self.separation,
            "separates": self.separates,
            "refusal": self.refusal,
            "evidence_role": self.evidence_role,
            "context": self.context_name,
            "context_registered_reading": self.context_reading,
            "counts_as_replication": self.counts_as_replication,
            "prior_reads": list(self.prior_reads),
            "may_be_reported_as_fresh_validation": self.may_be_reported_as_fresh_validation,
        }


# ---------------------------------------------------------------------------
# The registered locus convention, in code
# ---------------------------------------------------------------------------


def _ident(item: Any, fallback: int = 0) -> str:
    """The identifier of a record, whatever kind of record it is."""
    for name in ("ident", "unit_id", "candidate_id"):
        value = getattr(item, name, None)
        if value:
            return str(value)
    return str(fallback)


def independent_loci(items: Sequence[Any]) -> list[tuple[str, ...]]:
    """Group records into loci by `LOCUS_CONVENTION`; return each locus's identifiers, sorted.

    Same target gene, or same chromosome within `LOCUS_WINDOW`, is one locus, and the grouping is
    transitive. The number of groups is the number of independent loci under the convention, which
    is an operational grouping and not established biological independence.
    """
    idents = [_ident(it, i) for i, it in enumerate(items)]
    parent = list(range(len(items)))

    def find(i: int) -> int:
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    def union(i: int, j: int) -> None:
        a, b = find(i), find(j)
        if a != b:
            parent[max(a, b)] = min(a, b)

    # Same target gene is one locus, wherever the records lie.
    by_gene: dict[str, int] = {}
    for i, it in enumerate(items):
        gene = str(getattr(it, "target_gene", "") or "")
        if not gene:
            continue
        if gene in by_gene:
            union(by_gene[gene], i)
        else:
            by_gene[gene] = i

    # Same chromosome within the window is one locus. Sorting by position and linking neighbours
    # gives exactly the transitive closure of the within-window relation.
    order = sorted(
        range(len(items)),
        key=lambda i: (
            str(getattr(items[i], "chromosome", "")),
            int(getattr(items[i], "position", 0)),
        ),
    )
    for a, b in zip(order, order[1:], strict=False):
        ia, ib = items[a], items[b]
        if str(getattr(ia, "chromosome", "")) != str(getattr(ib, "chromosome", "")):
            continue
        if abs(int(getattr(ib, "position", 0)) - int(getattr(ia, "position", 0))) <= LOCUS_WINDOW:
            union(a, b)

    groups: dict[int, list[str]] = {}
    for i in range(len(items)):
        groups.setdefault(find(i), []).append(idents[i])
    return [tuple(sorted(g)) for g in sorted(groups.values(), key=lambda g: sorted(g))]


def locus_of(items: Sequence[Any]) -> dict[str, int]:
    """Map each identifier to its locus index under `LOCUS_CONVENTION`."""
    return {ident: n for n, group in enumerate(independent_loci(items)) for ident in group}


# ---------------------------------------------------------------------------
# Gate 1: the blinded eligibility count
# ---------------------------------------------------------------------------


def _unblinding_fields(unit: Any) -> tuple[str, ...]:
    return tuple(f for f in UNBLINDING_FIELDS if getattr(unit, f, None) not in (None, (), [], {}, ""))


def gate_eligibility(design: Design) -> GateResult:
    """Count units, label prevalence and independent loci with no prediction joined.

    Refused if the unit of analysis is not named, if any unit carries a prediction, a score or a
    hypothesis, or if any registered floor is unmet. This runs before anything in
    `design.candidates` is looked at, so a design can be refused here before a volume exists.
    """
    population = f"{design.unit_of_analysis or 'unnamed units'} offered by design {design.name!r}"
    if not design.unit_of_analysis or not design.unit_of_analysis.strip():
        return GateResult(
            GATE_ELIGIBILITY,
            False,
            "the unit of analysis is not named, so no count can state the population it covers",
        )

    units = list(design.eligibility)
    joined: dict[str, tuple[str, ...]] = {}
    for n, unit in enumerate(units):
        fields = _unblinding_fields(unit)
        if fields:
            joined[_ident(unit, n)] = fields
    if joined:
        named = ", ".join(f"{k} ({'/'.join(v)})" for k, v in sorted(joined.items())[:5])
        return GateResult(
            GATE_ELIGIBILITY,
            False,
            "a prediction was joined to the eligibility count: "
            f"{len(joined)} of {len(units)} units carry a model-derived field ({named}). The "
            "eligibility count is blinded by construction, so it must be taken over records that "
            "have none.",
            detail={"units_with_a_prediction": sorted(joined)},
        )

    n_units = len(units)
    n_positive = sum(1 for u in units if bool(getattr(u, "label", False)))
    prevalence = (n_positive / n_units) if n_units else 0.0
    loci = independent_loci(units)
    counts = (
        Count(n_units, population),
        Count(n_positive, f"{population}, with the label present"),
        Count(round(prevalence, 6), f"share of {population} with the label present"),
        Count(len(loci), f"independent loci among {population}, by the registered convention"),
    )
    detail = {
        "unit_of_analysis": design.unit_of_analysis,
        "locus_convention": LOCUS_CONVENTION,
        "loci": [list(g) for g in loci],
    }

    floors = design.floors
    for value, floor, what in (
        (n_units, floors.min_units, "units"),
        (n_positive, floors.min_positive_units, "units with the label present"),
        (prevalence, floors.min_label_prevalence, "label prevalence"),
        (len(loci), floors.min_independent_loci, "independent loci"),
    ):
        if value < floor:
            return GateResult(
                GATE_ELIGIBILITY,
                False,
                f"below the registered floor for {what}: {value} among {population}, floor "
                f"{floor} (registered {floors.registered_at}). No volume is reported for a design "
                "that is not eligible.",
                counts=counts,
                detail=detail,
            )

    return GateResult(
        GATE_ELIGIBILITY,
        True,
        f"blinded eligibility met: {n_units} units, {n_positive} with the label present "
        f"(prevalence {prevalence:.4f}), {len(loci)} independent loci among {population}. No "
        "prediction was joined.",
        counts=counts,
        detail=detail,
    )


# ---------------------------------------------------------------------------
# Gate 2: assay feasibility
# ---------------------------------------------------------------------------


def gate_feasibility(design: Design) -> GateResult:
    """Can the assay be run, and has the design said what it can see and what it compares with?

    An unavailable assay drops that candidate, with its registered availability reading. An
    unstated detection limit or an undeclared readout comparability refuses the whole design: both
    are omissions in the design, not properties of a candidate, and a comparison whose
    comparability was never declared is the failure this gate exists to catch.
    """
    population = f"candidates offered by design {design.name!r}"
    candidates = list(design.candidates)
    if not candidates:
        return GateResult(GATE_FEASIBILITY, False, f"no candidate is offered among {population}")

    missing_limit = [
        c.candidate_id for c in candidates if c.assay.detection_limit is None or c.assay.detection_limit <= 0
    ]
    if missing_limit:
        return GateResult(
            GATE_FEASIBILITY,
            False,
            f"{len(missing_limit)} of {len(candidates)} candidates among {population} state no "
            f"positive detection limit ({', '.join(sorted(missing_limit)[:5])}). Without it there "
            "is no way to say whether a predicted difference could be seen.",
            detail={"without_a_detection_limit": sorted(missing_limit)},
        )

    undeclared = [c.candidate_id for c in candidates if c.assay.readout_comparable is None]
    if undeclared:
        return GateResult(
            GATE_FEASIBILITY,
            False,
            f"{len(undeclared)} of {len(candidates)} candidates among {population} do not declare "
            f"whether their readout is comparable with the readout they will be set against "
            f"({', '.join(sorted(undeclared)[:5])}). Readout comparability is declared before the "
            "experiment or the comparison cannot be read afterwards.",
            detail={"comparability_undeclared": sorted(undeclared)},
        )

    not_comparable = [c.candidate_id for c in candidates if c.assay.readout_comparable is False]
    unavailable = {c.candidate_id: c.assay.availability_reading for c in candidates if not c.assay.available}
    feasible = [c for c in candidates if c.assay.available and c.assay.readout_comparable is True]
    counts = (
        Count(len(candidates), population),
        Count(
            len(feasible), f"{population} whose assay is available and whose readout is declared comparable"
        ),
        Count(len(unavailable), f"{population} whose assay is unavailable"),
        Count(len(not_comparable), f"{population} whose readout is declared not comparable"),
    )
    detail = {
        "feasible": [c.candidate_id for c in feasible],
        "unavailable": unavailable,
        "readout_not_comparable": sorted(not_comparable),
    }
    if not feasible:
        return GateResult(
            GATE_FEASIBILITY,
            False,
            f"no candidate among {population} has an available assay with a readout declared "
            "comparable, so there is nothing the design could measure.",
            counts=counts,
            detail=detail,
        )
    return GateResult(
        GATE_FEASIBILITY,
        True,
        f"{len(feasible)} of {len(candidates)} candidates among {population} are feasible: assay "
        f"available, detection limit stated, readout declared comparable. Dropped: "
        f"{len(unavailable)} with no available assay, {len(not_comparable)} whose readout is "
        "declared not comparable.",
        counts=counts,
        detail=detail,
    )


# ---------------------------------------------------------------------------
# Gate 3: distinguishable hypotheses
# ---------------------------------------------------------------------------


def information_value(candidate: Candidate) -> InformationValue:
    """Each hypothesis's predicted outcome, the detection limit, and whether they separate.

    Refused, in order: fewer than two hypotheses; predictions not in the assay's units; identical
    predictions (indistinguishable whatever the assay can see); a separation below the detection
    limit (the assay cannot report it).
    """
    assay = candidate.assay
    predicted = tuple((h.name, h.predicted_outcome) for h in candidate.hypotheses)
    base = {
        "candidate_id": candidate.candidate_id,
        "units": assay.units,
        "predicted": predicted,
        "detection_limit": assay.detection_limit,
        "evidence_role": candidate.evidence_role,
        "context_name": candidate.context.name,
        "context_reading": candidate.context.registered_reading,
        "counts_as_replication": candidate.context.counts_as_replication,
        "prior_reads": tuple(candidate.prior_reads),
    }

    if len(candidate.hypotheses) < 2:
        return InformationValue(
            separation=None,
            separates=False,
            refusal=(
                f"{len(candidate.hypotheses)} hypothesis on candidate {candidate.candidate_id}: "
                "a single hypothesis has nothing to be separated from, so the result carries no "
                "information about which model is right."
            ),
            **base,
        )

    wrong_units = sorted({h.units for h in candidate.hypotheses if h.units != assay.units})
    if wrong_units:
        return InformationValue(
            separation=None,
            separates=False,
            refusal=(
                f"candidate {candidate.candidate_id} predicts in {wrong_units} but the assay reads "
                f"in {assay.units!r}: the predictions and the detection limit are not on one scale, "
                "so no separation can be compared with what the assay can see."
            ),
            **base,
        )

    values = [h.predicted_outcome for h in candidate.hypotheses]
    separation = max(values) - min(values)
    limit = float(assay.detection_limit) if assay.detection_limit is not None else None

    if separation == 0:
        pair = ", ".join(h.name for h in candidate.hypotheses)
        return InformationValue(
            separation=separation,
            separates=False,
            refusal=(
                f"the hypotheses on candidate {candidate.candidate_id} ({pair}) predict the same "
                f"outcome ({values[0]} {assay.units}), so no result can tell them apart. This is "
                "not a limit of the assay."
            ),
            **base,
        )
    if limit is None or separation < limit:
        return InformationValue(
            separation=separation,
            separates=False,
            refusal=(
                f"the separation on candidate {candidate.candidate_id} is {separation} "
                f"{assay.units}, below the assay's detection limit of {limit} {assay.units}: the "
                "assay could not report the difference, so the candidate is refused."
            ),
            **base,
        )
    return InformationValue(
        separation=separation,
        separates=True,
        refusal="",
        **base,
    )


def gate_distinguishability(
    design: Design, feasible: Sequence[Candidate]
) -> tuple[GateResult, tuple[InformationValue, ...]]:
    """Keep only candidates whose hypotheses separate by more than the assay's detection limit."""
    population = f"feasible candidates of design {design.name!r}"
    values = tuple(information_value(c) for c in feasible)
    separating = [v for v in values if v.separates]
    identical = [v.candidate_id for v in values if v.separation == 0]
    below = [v.candidate_id for v in values if v.separation not in (None, 0) and not v.separates]
    counts = (
        Count(len(values), population),
        Count(len(separating), f"{population} whose hypotheses separate above the detection limit"),
        Count(len(identical), f"{population} whose hypotheses predict the same outcome"),
        Count(len(below), f"{population} whose separation is below the detection limit"),
    )
    detail = {
        "separating": [v.candidate_id for v in separating],
        "identical_predictions": sorted(identical),
        "below_the_detection_limit": sorted(below),
        "refusals": {v.candidate_id: v.refusal for v in values if v.refusal},
    }
    floor = design.floors.min_distinguishing_candidates
    if len(separating) < floor:
        return (
            GateResult(
                GATE_DISTINGUISHABILITY,
                False,
                f"{len(separating)} of {len(values)} {population} would separate the hypotheses, "
                f"below the registered floor of {floor} (registered "
                f"{design.floors.registered_at}). Dropped: {len(identical)} with identical "
                f"predictions, {len(below)} with a separation below the detection limit.",
                counts=counts,
                detail=detail,
            ),
            values,
        )
    return (
        GateResult(
            GATE_DISTINGUISHABILITY,
            True,
            f"{len(separating)} of {len(values)} {population} would separate the hypotheses by more "
            f"than the assay's detection limit. Dropped: {len(identical)} with identical "
            f"predictions, {len(below)} below the detection limit.",
            counts=counts,
            detail=detail,
        ),
        values,
    )


# ---------------------------------------------------------------------------
# Gate 4: prior exposure and what may be called a replication
# ---------------------------------------------------------------------------


def gate_freshness(design: Design, values: Sequence[InformationValue]) -> GateResult:
    """Prior exposure travels, and only a context registered as replicated counts as a replication.

    A design claiming fresh validation is refused if any surviving candidate has been read: a
    re-read of data already seen is development evidence, however it is later described. A design
    claiming a replication is refused if no surviving candidate's context carries a registered
    `replicated` reading; a context registered `inconclusive` stays inconclusive and is never
    counted, and its registered reading is carried into the result word for word.
    """
    surviving = [v for v in values if v.separates]
    population = f"distinguishing candidates of design {design.name!r}"
    read_before = [v for v in surviving if v.prior_reads]
    replications = [v for v in surviving if v.counts_as_replication]
    contexts: dict[str, dict[str, Any]] = {}
    for v in surviving:
        contexts.setdefault(
            v.context_name,
            {
                "registered_reading": v.context_reading,
                "counts_as_replication": v.counts_as_replication,
            },
        )
    counts = (
        Count(len(surviving), population),
        Count(len(read_before), f"{population} already read, which are development evidence"),
        Count(
            len(surviving) - len(read_before),
            f"{population} never read, which could be fresh validation",
        ),
        Count(
            len(replications),
            f"{population} in a context whose registered reading is {REPLICATED!r}",
        ),
    )
    detail = {
        "claim": design.claim,
        "development_evidence": {v.candidate_id: list(v.prior_reads) for v in read_before},
        "contexts": contexts,
        "counted_as_replication": [v.candidate_id for v in replications],
    }

    if design.claim == FRESH_VALIDATION and read_before:
        named = ", ".join(sorted(v.candidate_id for v in read_before)[:5])
        return GateResult(
            GATE_FRESHNESS,
            False,
            f"the design claims {FRESH_VALIDATION!r}, but {len(read_before)} of {len(surviving)} "
            f"{population} have been read before ({named}). A candidate already read is "
            f"{DEVELOPMENT_EVIDENCE}, and re-reading it is not fresh validation. Either drop those "
            f"candidates or register the design as {DEVELOPMENT_EVIDENCE!r}.",
            counts=counts,
            detail=detail,
        )
    if design.claim == REPLICATION and not replications:
        readings = "; ".join(f"{name}: {c['registered_reading']}" for name, c in sorted(contexts.items()))
        return GateResult(
            GATE_FRESHNESS,
            False,
            f"the design claims {REPLICATION!r}, but no context among {population} carries a "
            f"registered {REPLICATED!r} reading. The registered readings are -- {readings}. A "
            "context registered as inconclusive stays inconclusive and is not counted as a "
            "replication.",
            counts=counts,
            detail=detail,
        )
    return GateResult(
        GATE_FRESHNESS,
        True,
        f"prior exposure recorded: {len(read_before)} of {len(surviving)} {population} are "
        f"{DEVELOPMENT_EVIDENCE} and can never be reported as fresh validation; "
        f"{len(replications)} sit in a context whose registered reading is {REPLICATED!r}, and "
        "every other context keeps its own registered reading.",
        counts=counts,
        detail=detail,
    )


# ---------------------------------------------------------------------------
# The two arms
# ---------------------------------------------------------------------------


def select_representative(
    design: Design, values: Sequence[InformationValue], seed: int = SEED
) -> tuple[Selection, ...]:
    """Draw loci uniformly, then one candidate uniformly inside each drawn locus.

    The inclusion probability of a candidate is (loci requested / loci available) times
    (1 / candidates in its locus), and it is recorded on the selection. Because the draw does not
    look at the separation, counts over this arm estimate how common a function is in the eligible
    population -- and in no wider population.
    """
    separating = {v.candidate_id: v for v in values if v.separates}
    pool_candidates = [c for c in design.candidates if c.candidate_id in separating]
    groups = independent_loci(pool_candidates)
    available = len(groups)
    if available == 0:
        return ()
    requested = min(design.floors.target_loci, available)
    rng = random.Random(seed)
    drawn = rng.sample(range(available), requested)
    scheme = (
        f"{requested} of {available} independent loci drawn uniformly without replacement "
        f"(seed {seed}), then one candidate drawn uniformly inside each locus"
    )
    out: list[Selection] = []
    for index in sorted(drawn):
        members = sorted(groups[index])
        chosen = members[rng.randrange(len(members))]
        out.append(
            Selection(
                candidate_id=chosen,
                arm=REPRESENTATIVE,
                selection_probability=(requested / available) * (1.0 / len(members)),
                scheme=scheme,
                representative_of_population=True,
            )
        )
    return tuple(out)


def select_disagreement(design: Design, values: Sequence[InformationValue]) -> tuple[Selection, ...]:
    """Take the candidates where the hypotheses separate most, relative to the detection limit.

    One candidate per locus, ranked, so the arm spends its volume on distinct loci. Selection is
    deterministic, so the probability recorded is 1.0 -- and the arm is representative of no
    population, because its candidates are chosen by the quantity under test.
    """
    separating = [v for v in values if v.separates]
    if not separating:
        return ()
    by_id = {c.candidate_id: c for c in design.candidates}
    pool_candidates = [by_id[v.candidate_id] for v in separating if v.candidate_id in by_id]
    locus = locus_of(pool_candidates)

    def margin(v: InformationValue) -> float:
        limit = float(v.detection_limit) if v.detection_limit else 1.0
        return float(v.separation or 0.0) / limit

    ranked = sorted(separating, key=lambda v: (-margin(v), v.candidate_id))
    scheme = (
        "deterministic rank on the separation divided by the assay's detection limit, one candidate "
        "per independent locus; selected with probability 1, representative of no population"
    )
    taken: set[int] = set()
    out: list[Selection] = []
    for v in ranked:
        if len(out) >= design.floors.target_loci:
            break
        index = locus.get(v.candidate_id)
        if index is None or index in taken:
            continue
        taken.add(index)
        out.append(
            Selection(
                candidate_id=v.candidate_id,
                arm=DISAGREEMENT,
                selection_probability=1.0,
                scheme=scheme,
                representative_of_population=False,
            )
        )
    return tuple(out)


def pool(selections: Iterable[Selection]) -> int:
    """Count selections from one arm. Refuses a pool that mixes the arms."""
    items = list(selections)
    arms = {s.arm for s in items}
    if len(arms) > 1:
        raise ValueError(MERGE_REFUSED)
    return len(items)


# ---------------------------------------------------------------------------
# Volume, last
# ---------------------------------------------------------------------------


def volume(design: Design, arm: str, selections: Sequence[Selection], loci_available: int) -> dict[str, Any]:
    """How much the arm asks for. Only ever called when every gate has passed."""
    if arm not in ARMS:
        raise ValueError(f"an arm must be one of {ARMS}, not {arm!r}")
    population = f"independent loci eligible for the {arm} arm of design {design.name!r}"
    requested = len(selections)
    per_locus = design.floors.units_per_locus
    return {
        "arm": arm,
        "question": ARM_QUESTIONS[arm],
        "unit_of_analysis": design.unit_of_analysis,
        "loci_available": Count(loci_available, population).as_dict(),
        "loci_requested": Count(requested, f"{population}, selected").as_dict(),
        "units_per_locus": per_locus,
        "units_requested": Count(
            requested * per_locus,
            f"{design.unit_of_analysis} requested by the {arm} arm of design {design.name!r}",
        ).as_dict(),
        "basis": VOLUME_BASIS,
    }


# ---------------------------------------------------------------------------
# The whole order
# ---------------------------------------------------------------------------


def plan(design: Design, seed: int = SEED) -> dict[str, Any]:
    """Run the gates in order; refuse with a reason, or return the two arms and their volumes.

    A refused design has `refused` set to the failing gate's reason and carries no `arms` key and no
    volume anywhere in the result. Only a design that passes all four gates is answered with a
    volume, and the two arms are reported separately.
    """
    out: dict[str, Any] = {
        "design": design.name,
        "unit_of_analysis": design.unit_of_analysis,
        "claim": design.claim,
        "source_reading": design.source_reading,
        "floors": design.floors.as_dict(),
        "locus_convention": LOCUS_CONVENTION,
        "locus_window": LOCUS_WINDOW,
        "merge_refused": MERGE_REFUSED,
        "gates": [],
        "refused": None,
    }

    eligibility = gate_eligibility(design)
    out["gates"].append(eligibility.as_dict())
    if not eligibility.passed:
        out["refused"] = eligibility.reason
        out["refused_at"] = GATE_ELIGIBILITY
        return out

    feasibility = gate_feasibility(design)
    out["gates"].append(feasibility.as_dict())
    if not feasibility.passed:
        out["refused"] = feasibility.reason
        out["refused_at"] = GATE_FEASIBILITY
        return out

    feasible_ids = set(feasibility.detail.get("feasible", ()))
    feasible = [c for c in design.candidates if c.candidate_id in feasible_ids]
    distinguishability, values = gate_distinguishability(design, feasible)
    out["gates"].append(distinguishability.as_dict())
    out["information_value"] = [v.as_dict() for v in values]
    if not distinguishability.passed:
        out["refused"] = distinguishability.reason
        out["refused_at"] = GATE_DISTINGUISHABILITY
        return out

    freshness = gate_freshness(design, values)
    out["gates"].append(freshness.as_dict())
    if not freshness.passed:
        out["refused"] = freshness.reason
        out["refused_at"] = GATE_FRESHNESS
        return out

    separating_ids = {v.candidate_id for v in values if v.separates}
    pool_candidates = [c for c in design.candidates if c.candidate_id in separating_ids]
    loci_available = len(independent_loci(pool_candidates))
    representative = select_representative(design, values, seed=seed)
    disagreement = select_disagreement(design, values)
    out["arms"] = {
        REPRESENTATIVE: {
            "question": ARM_QUESTIONS[REPRESENTATIVE],
            "selections": [s.as_dict() for s in representative],
            "volume": volume(design, REPRESENTATIVE, representative, loci_available),
        },
        DISAGREEMENT: {
            "question": ARM_QUESTIONS[DISAGREEMENT],
            "selections": [s.as_dict() for s in disagreement],
            "volume": volume(design, DISAGREEMENT, disagreement, loci_available),
        },
    }
    return out
