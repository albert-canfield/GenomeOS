"""Targets the tumour never mutated: the expression-driven route.

Everything upstream of this module starts from an alteration — a variant, a
copy-number event, a rearrangement — and asks what can reach it. That is not
how the two best-known cell therapies were found. CD19 and BCMA are not
mutated, not amplified and not rearranged; they are lineage antigens that the
malignant cell carries in quantity and most healthy tissue does not. A pipeline
that can only reach a gene through an alteration cannot propose either of them,
which is a hole in the middle of the target space rather than a missing
refinement.

The scan closes it with the two things now available: the patient's own tumour
RNA, and the packaged healthy-tissue atlas over the membrane and CD-marker
universe. A gene qualifies when the tumour carries a lot of its transcript and
healthy tissue carries far less, and both halves of that sentence are
measurements.

Three things it refuses. It runs only on patient RNA: a cohort describes a
cancer type and cannot say what this tumour displays, so without patient RNA
the answer is that no expression-driven target can be proposed, not a weaker
guess. It never calls a transcript a protein, still less a protein on the
surface. And it says out loud what the CD19 story actually teaches: an
unaltered antigen is shared with the healthy lineage that normally carries it,
so the price of hitting it is that lineage — B-cell aplasia is not a side
effect of the approach, it is the approach working.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from . import atlas
from .evidence import Evidence, derived
from .evidence import patient as patient_evidence
from .providers import CRITICAL_TISSUES, hpa_column

#: Transcript level below which a gene is not worth calling a target, whatever
#: its ratio: a fivefold rise over nothing is still nothing.
FLOOR = 10.0

#: Tumour-to-healthy ratio a gene must clear. Threefold is deliberately low for
#: a screen; what separates the hits is reported, not hidden in the cutoff.
MIN_RATIO = 3.0


@dataclass(slots=True)
class ExpressionHit:
    """One gene the tumour carries far more of than healthy tissue does."""

    gene: str
    tumour: float
    unit: str
    normal_max: float
    normal_tissue: str
    ratio: float
    weighted_worst: tuple[str, float] | None = None
    rank: float | None = None
    specificity: str | None = None
    cd_marker: bool = False
    passes: bool = True
    rejected_because: str = ""
    evidence: list[Evidence] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return {
            "gene": self.gene,
            "tumour": self.tumour,
            "unit": self.unit,
            "highest_healthy_tissue": self.normal_tissue,
            "highest_healthy_nTPM": self.normal_max,
            "ratio": round(self.ratio, 2),
            "most_weighted_healthy_tissue": (
                {"tissue": self.weighted_worst[0], "nTPM": self.weighted_worst[1]}
                if self.weighted_worst
                else None
            ),
            "rank_in_sample": self.rank,
            "population_specificity": self.specificity,
            "cd_marker": self.cd_marker,
            "passes": self.passes,
            "rejected_because": self.rejected_because,
            "evidence": [e.to_dict() for e in self.evidence],
            "level": (
                "the tumour value is this patient's RNA; the healthy values are population-level "
                "consensus nTPM. Both are per-million normalisations and the ratio is an "
                "approximation, not a like-for-like measurement"
            ),
        }


def _healthy(row: dict[str, Any]) -> tuple[float, str, tuple[str, float] | None]:
    """The highest healthy tissue, and the highest once the weights are applied."""
    best, best_tissue = 0.0, "none above background"
    worst_weighted, worst_load = None, -1.0
    for field_id, (tissue, weight) in CRITICAL_TISSUES.items():
        value = row.get(hpa_column(field_id))
        if not isinstance(value, int | float):
            continue
        if value > best:
            best, best_tissue = float(value), tissue
        if value * weight > worst_load:
            worst_load, worst_weighted = value * weight, (tissue, float(value))
    return best, best_tissue, worst_weighted


def scan(
    patient_rna: Any,
    limit: int = 8,
    floor: float = FLOOR,
    min_ratio: float = MIN_RATIO,
    universe: list[str] | None = None,
    rejected: int = 6,
) -> list[ExpressionHit]:
    """Genes this tumour carries far more of than healthy tissue does.

    Every gene above the floor is returned, the ones that cleared the ratio
    first and up to `rejected` of the ones that did not, each carrying why.
    The rejections are the more instructive half: in a B-cell malignancy the
    genes that fail here are CD20, CD79A, CD79B and CD22 — the targets of the
    most successful antibodies in oncology — because a lineage antigen has a
    healthy home full of the same lineage, and the screen measures selectivity
    against healthy tissue rather than usefulness. A threshold here is a sort
    order, never a verdict.

    Returns an empty list, not a guess, when there is no patient RNA or no
    packaged atlas to compare it against.
    """
    if not getattr(patient_rna, "available", False):
        return []
    genes = universe if universe is not None else atlas.universe()
    if not genes:
        return []
    hits: list[ExpressionHit] = []
    for gene in genes:
        value = patient_rna.values.get(gene)
        if value is None or value < floor:
            continue
        row = atlas.row(gene)
        if row is None:
            continue
        normal_max, tissue, weighted = _healthy(row)
        ratio = (value + 1.0) / (normal_max + 1.0)
        hit = ExpressionHit(
            gene=gene,
            tumour=float(value),
            unit=getattr(patient_rna, "unit", "TPM"),
            normal_max=normal_max,
            normal_tissue=tissue,
            ratio=ratio,
            weighted_worst=weighted,
            rank=patient_rna.rank(gene) if hasattr(patient_rna, "rank") else None,
            specificity=row.get("RNA tissue specificity"),
            cd_marker=bool(row.get("cd_marker")),
            passes=ratio >= min_ratio,
        )
        if not hit.passes:
            hit.rejected_because = (
                f"{ratio:.1f}x is below the {min_ratio:g}x the screen asks for: healthy {tissue} "
                f"already carries {normal_max:g} nTPM of it. That is a statement about selectivity "
                "and not about usefulness — an antigen shared with a healthy lineage can still be "
                "the right target if losing that lineage is survivable"
            )
        hit.evidence.append(
            patient_evidence(
                f"tumour RNA-seq ({getattr(patient_rna, 'path', '') or 'supplied table'})",
                f"{gene} measured at {value:g} {hit.unit} in this tumour",
                0.9,
            )
        )
        hit.evidence.append(
            derived(
                "GenomeOS expression scan",
                f"{gene} is {ratio:.1f}x the highest of {len(CRITICAL_TISSUES)} queried healthy tissues "
                f"({tissue} at {normal_max:g} nTPM); it is raised, not altered, and raised transcript "
                "is not protein on the surface",
                0.5,
            )
        )
        hits.append(hit)
    hits.sort(key=lambda h: -h.ratio)
    kept = [h for h in hits if h.passes][:limit]
    turned_down = [h for h in hits if not h.passes][:rejected]
    return kept + turned_down


#: What an unaltered, over-expressed antigen costs. Carried on every candidate
#: the scan produces, because it is the lesson of the one that worked.
LINEAGE_CAVEAT = (
    "this gene is not altered in this tumour: it is proposed because the tumour's RNA carries far "
    "more of it than healthy tissue does. An unaltered antigen is shared with the healthy lineage "
    "that normally carries it, and hitting it hits that lineage too — CD19 therapy works and B-cell "
    "aplasia is not a side effect of it but the same event seen from the other side"
)

TRANSCRIPT_CAVEAT = (
    "the measurement is transcript in the tumour against transcript in healthy tissue. It does not "
    "show that the protein is made, that it reaches the surface, or how much of it is there; that "
    "needs proteomics, and surface proteomics for the last of the three"
)

UNIT_CAVEAT = (
    "the patient's value and the healthy values come from different experiments and different "
    "normalisations (TPM against consensus nTPM), so the ratio is an approximation and a twofold "
    "difference in it means little"
)


def limitations(hit: ExpressionHit) -> list[str]:
    out = [LINEAGE_CAVEAT, TRANSCRIPT_CAVEAT, UNIT_CAVEAT]
    if hit.weighted_worst and hit.weighted_worst[1] >= 10.0:
        out.append(
            f"healthy {hit.weighted_worst[0]} carries {hit.weighted_worst[1]:g} nTPM of it, and that "
            "tissue is among the ones this pipeline weighs hardest"
        )
    if hit.rank is None:
        out.append(
            "the supplied RNA table is too small for a within-sample rank, so the only comparison "
            "available is the cross-experiment one"
        )
    return out
