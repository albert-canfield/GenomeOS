"""Two antigens can be more selective than either one.

There is rarely a perfectly tumour-specific antigen. There is often a pair
whose *combination* is specific:

    healthy tissue 1   A+ B-
    healthy tissue 2   A- B+
    tumour             A+ B+

`A AND B` then addresses the tumour and neither healthy tissue. This module
builds and assesses such combinations against the healthy-tissue profiles it
has, and reports the selectivity gain with the caveat that decides whether it
is real: consensus tissue RNA is bulk, so two genes expressed in one tissue
are not necessarily expressed in the same *cell*. Establishing that needs
single-cell data, and the requirement is stated rather than assumed away.

Output is molecular logic and evidence only. Multispecific binders,
logic-gated cell therapies and conditional binders are downstream design
problems that GenomeOS does not attempt.
"""

from __future__ import annotations

from itertools import combinations
from typing import Any

from .evidence import derived
from .expression import CONCERN_NTPM
from .model import TargetLogic, TherapeuticTargetCandidate

OPERATORS = ("SINGLE", "AND", "OR", "AND_NOT")


def _levels(c: TherapeuticTargetCandidate) -> dict[str, float]:
    return {t.tissue: (t.value or 0.0) for t in c.normal_tissue.tissues if t.value is not None}


def _weights(c: TherapeuticTargetCandidate) -> dict[str, float]:
    return {t.tissue: t.weight for t in c.normal_tissue.tissues}


def _load(levels: dict[str, float], weights: dict[str, float]) -> tuple[str, float]:
    if not levels:
        return "", 0.0
    tissue = max(levels, key=lambda t: levels[t] * weights.get(t, 0.5))
    return tissue, levels[tissue] * weights.get(tissue, 0.5)


def pairs(candidates: list[TherapeuticTargetCandidate], limit: int = 3) -> list[TargetLogic]:
    """AND combinations of surface candidates, ranked by selectivity gain."""
    usable = [
        c
        for c in candidates
        if c.target_class in ("direct_surface", "pathway_induced_surface") and c.normal_tissue.known
    ]
    out: list[TargetLogic] = []
    for a, b in combinations(usable, 2):
        la, lb = _levels(a), _levels(b)
        shared = sorted(set(la) & set(lb))
        if not shared:
            continue
        wa = _weights(a)
        both = {t: min(la[t], lb[t]) for t in shared}
        worst_a, load_a = _load({t: la[t] for t in shared}, wa)
        worst_b, load_b = _load({t: lb[t] for t in shared}, wa)
        worst_both, load_both = _load(both, wa)
        best_single = max(load_a, load_b)
        gain = None
        if load_both > 0:
            gain = round(best_single / load_both, 2)
        elif best_single > 0:
            gain = None  # a zero denominator is not a gain, it is a different statement
        excluded = [t for t in shared if max(la[t], lb[t]) >= CONCERN_NTPM and both[t] < CONCERN_NTPM]
        still = [t for t in shared if both[t] >= CONCERN_NTPM]
        logic = TargetLogic(
            operator="AND",
            targets=[a.gene, b.gene],
            rationale=(
                f"requiring both {a.gene} and {b.gene} lowers the worst weighted healthy-tissue load "
                f"from {best_single:.1f} ({worst_a if load_a >= load_b else worst_b}) to "
                f"{load_both:.1f} ({worst_both or 'none'})"
            ),
            normal_tissues_excluded=excluded,
            normal_tissues_still_at_risk=[f"{t} ({both[t]:g} nTPM)" for t in still[:6]],
            selectivity_gain=gain,
        )
        logic.evidence.append(
            derived(
                "GenomeOS combination logic",
                f"{a.gene} AND {b.gene}: {len(excluded)} healthy tissues drop below the "
                f"{CONCERN_NTPM:g} nTPM concern threshold when both antigens are required",
                0.4,
            )
        )
        logic.evidence.append(
            derived(
                "GenomeOS combination logic",
                "bulk tissue RNA cannot show that two genes are expressed in the same cell; "
                "single-cell or co-staining data is required before this combination is treated as real",
                0.0,
            )
        )
        out.append(logic)
    out.sort(key=lambda x: -(x.selectivity_gain or 0))
    return out[:limit]


def single(c: TherapeuticTargetCandidate) -> TargetLogic:
    """The default: one antigen, with what it costs stated."""
    return TargetLogic(
        operator="SINGLE",
        targets=[c.gene],
        rationale="one antigen; selectivity rests entirely on the tumour-versus-normal difference of "
        "this target alone",
        normal_tissues_still_at_risk=list(c.normal_tissue.tissues_at_risk),
    )


def and_not(positive: TherapeuticTargetCandidate, negative: TherapeuticTargetCandidate) -> TargetLogic | None:
    """`A AND NOT B`: a healthy-tissue marker used as a veto.

    Only meaningful when the tumour is known not to carry B, which needs
    tumour expression data. Without it the combination is proposed with the
    requirement attached rather than asserted.
    """
    la, lb = _levels(positive), _levels(negative)
    shared = sorted(set(la) & set(lb))
    if not shared:
        return None
    vetoed = [t for t in shared if la[t] >= CONCERN_NTPM and lb[t] >= CONCERN_NTPM]
    if not vetoed:
        return None
    logic = TargetLogic(
        operator="AND_NOT",
        targets=[positive.gene, negative.gene],
        rationale=(
            f"healthy tissues expressing both {positive.gene} and {negative.gene} could be spared by a "
            f"binder that is vetoed by {negative.gene}"
        ),
        normal_tissues_excluded=vetoed,
    )
    logic.evidence.append(
        derived(
            "GenomeOS combination logic",
            f"requires evidence that this tumour does not express {negative.gene}; without tumour "
            "RNA-seq or protein data that cannot be established, so the veto is a hypothesis",
            0.0,
        )
    )
    return logic


def describe(logic: TargetLogic | None) -> str:
    if logic is None:
        return "none"
    if logic.operator == "SINGLE":
        return logic.targets[0]
    if logic.operator == "AND_NOT":
        return f"{logic.targets[0]} AND NOT {logic.targets[1]}"
    return f" {logic.operator} ".join(logic.targets)


def summary(all_logic: list[TargetLogic]) -> dict[str, Any]:
    return {
        "supported_operators": list(OPERATORS),
        "combinations": [x.to_dict() for x in all_logic],
        "limitation": (
            "combinations are assessed against bulk healthy-tissue RNA; co-expression in the same cell "
            "is not established by that data and requires single-cell or co-staining evidence"
        ),
    }
