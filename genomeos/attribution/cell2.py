# SPDX-License-Identifier: AGPL-3.0-or-later
"""Eligibility of a second cell type for the CRISPRi deletion result: how many independent loci?

The CRISPRi deletion result is measured in one cell type, K562. Whether a second held-out cell type
can carry a stratified test is decided by the number of *independent* loci among its measured
positives, not by its pair count: four checks before this one (N1, the wiring diagnostic,
lane-contrast and lane-entex) each counted pairs or positives and none counted loci, so each could
have reported a stratum whose positives sit on a handful of places in the genome.

This module registers the counting convention before any count is taken, and takes the count from
labels only. It reads a pair's cell type, chromosome, element interval, measured gene and measured
`regulated` label, and nothing else: `locus_key` is the only path from a `crispri.Pair` into the
grouping, and `group` sees nothing but those keys. No prediction, no deletion value, no model score
and no performance metric is read, computed or written here.

The grouping is operational, not established biological independence. Two positives that this module
calls one locus may still be two separate regulatory events, and two it calls separate may share a
domain, a TAD or a trans factor. It is a conservative bookkeeping device for deciding whether a
stratum has enough distinct places in the genome to be worth measuring at all; it is not a claim
about biology, and nothing downstream may cite a locus count as an independence guarantee.
"""

from __future__ import annotations

from collections.abc import Iterable
from typing import Any, NamedTuple

#: The genomic span inside which two positives of the same chromosome count as one locus. 1 Mb is the
#: order of a mammalian topologically associating domain, so two elements this close are routinely
#: regulated together and their labels cannot be treated as two independent draws. The figure is a
#: round, pre-set bound chosen before any count was taken, not fitted to an outcome; it is recorded
#: here so that no later run can move it after seeing a number.
INDEPENDENT_LOCUS_SPAN = 1_000_000

#: The convention in words, carried verbatim into every result and summary that cites a locus count.
INDEPENDENT_LOCUS_RULE = (
    "two measured positives are the same locus when they share the measured gene, or when their "
    "elements lie within 1 Mb of one another on the same chromosome; the relation is chained within a "
    "chromosome, so a run of elements each within 1 Mb of the next is one locus however far its ends "
    "are apart. This is an operational grouping for deciding whether a stratum has enough distinct "
    "places in the genome to measure, not established biological independence"
)

#: What the count is allowed to read. Blinding condition of the lane: a predicted value, a deletion
#: value, a model score or an AUPRC may not enter this module.
BLINDED_FIELDS = ("cell", "chrom", "start", "end", "gene", "regulated")

#: The pre-set floor. Pooled independent loci among the candidate second cell types' positives below
#: this is a no-go: the lane records it and stops, and the floor is never moved after a count.
POOLED_LOCUS_FLOOR = 20

#: The held-out cell type the result already covers, so not a candidate for a second one.
PRIMARY_CELL = "K562"


class LocusKey(NamedTuple):
    """Everything the grouping is permitted to see about one measured positive."""

    cell: str
    chrom: str
    start: int
    end: int
    gene: str


def locus_key(pair: Any) -> LocusKey:
    """The blinded key of one pair: the only path from a `crispri.Pair` into the grouping."""
    return LocusKey(pair.cell, pair.chrom, int(pair.start), int(pair.end), pair.gene)


def gap(a: LocusKey, b: LocusKey) -> int:
    """Base pairs between two element intervals on one chromosome; 0 where they overlap or touch."""
    return max(0, max(a.start, b.start) - min(a.end, b.end))


def group(keys: Iterable[LocusKey], span: int = INDEPENDENT_LOCUS_SPAN) -> list[int]:
    """The locus group of each key, by INDEPENDENT_LOCUS_RULE. Group ids are the index of the first
    key of each group, so the result is stable under the input order it was given."""
    ks = list(keys)
    parent = list(range(len(ks)))

    def find(i: int) -> int:
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    def union(i: int, j: int) -> None:
        a, b = find(i), find(j)
        if a != b:
            parent[max(a, b)] = min(a, b)

    by_gene: dict[str, int] = {}
    by_chrom: dict[str, list[int]] = {}
    for i, k in enumerate(ks):
        if k.gene in by_gene:
            union(by_gene[k.gene], i)
        else:
            by_gene[k.gene] = i
        by_chrom.setdefault(k.chrom, []).append(i)

    for members in by_chrom.values():
        # sorted by start, the running maximum end is the nearest edge of everything already in the
        # cluster, so a key farther than `span` from it is farther than `span` from all of them.
        order = sorted(members, key=lambda i: (ks[i].start, ks[i].end))
        reach = ks[order[0]].end
        head = order[0]
        for i in order[1:]:
            if ks[i].start - reach <= span:
                union(head, i)
            else:
                head = i
            reach = max(reach, ks[i].end)

    return [find(i) for i in range(len(ks))]


def count_loci(keys: Iterable[LocusKey], span: int = INDEPENDENT_LOCUS_SPAN) -> int:
    """How many independent loci a set of blinded keys holds."""
    return len(set(group(keys, span)))


def _quantiles(values: list[float]) -> dict[str, Any]:
    """Minimum, median and maximum of a list; nulls for an empty one. No interpolation: the median of
    an even-length list is the mean of the two middle values."""
    if not values:
        return {"min": None, "median": None, "max": None}
    s = sorted(values)
    mid = len(s) // 2
    median = s[mid] if len(s) % 2 else (s[mid - 1] + s[mid]) / 2
    return {"min": s[0], "median": median, "max": s[-1]}


def locus_shapes(keys: Iterable[LocusKey], span: int = INDEPENDENT_LOCUS_SPAN) -> dict[str, Any]:
    """How wide each locus is and how many positives it holds, so a reader can see whether one long
    chain is carrying the count. Added on the reviewer's instruction of 2026-10-01; it reports beside
    INDEPENDENT_LOCUS_RULE and changes neither the rule nor POOLED_LOCUS_FLOOR."""
    ks = list(keys)
    if not ks:
        return {
            "loci": 0,
            "span_mb": _quantiles([]),
            "positives_per_locus": _quantiles([]),
            "largest_locus_positives": 0,
            "largest_locus_share_of_positives": None,
            "largest_locus_span_mb": None,
            "positives_per_locus_counts": {},
        }
    ids = group(ks, span)
    members: dict[int, list[LocusKey]] = {}
    for gid, k in zip(ids, ks, strict=True):
        members.setdefault(gid, []).append(k)
    spans_mb, sizes = [], []
    widest_of_largest = None
    largest = 0
    for group_keys in members.values():
        width = max(k.end for k in group_keys) - min(k.start for k in group_keys)
        spans_mb.append(round(width / 1_000_000, 4))
        sizes.append(len(group_keys))
        if len(group_keys) > largest:
            largest = len(group_keys)
            widest_of_largest = round(width / 1_000_000, 4)
    counts: dict[str, int] = {}
    for n in sorted(sizes):
        counts[str(n)] = counts.get(str(n), 0) + 1
    return {
        "loci": len(members),
        "span_mb": _quantiles(spans_mb),
        "span_mb_note": (
            "the span of a locus is the distance from the leftmost start to the rightmost end of its "
            "elements, so a chained locus can be wider than the 1 Mb joining span"
        ),
        "positives_per_locus": _quantiles([float(n) for n in sizes]),
        "positives_per_locus_counts": counts,
        "positives_per_locus_counts_note": "locus size in positives -> how many loci have that size",
        "largest_locus_positives": largest,
        "largest_locus_share_of_positives": round(largest / len(ks), 4),
        "largest_locus_span_mb": widest_of_largest,
    }


def deletion_available(cell: str, model_cells: Iterable[str]) -> bool:
    """Whether a deletion value exists for this cell type at all. A cell outside the lines the
    deletion table was scored in cannot carry a deletion gain: see crispri.gain_where_available and
    crispri.UNAVAILABLE_GAIN, which refuse a number there rather than report an artefact."""
    return cell in tuple(model_cells)


def eligibility(
    pairs: Iterable[Any],
    model_cells: Iterable[str],
    primary: str = PRIMARY_CELL,
    span: int = INDEPENDENT_LOCUS_SPAN,
) -> dict[str, Any]:
    """Pairs, measured positives and independent loci per held-out cell type, and pooled over the
    candidate second cell types. Labels, coordinates, genes and cell types only."""
    ps = list(pairs)
    cells = sorted({p.cell for p in ps})
    per_cell: dict[str, Any] = {}
    for cell in cells:
        of_cell = [p for p in ps if p.cell == cell]
        positives = [locus_key(p) for p in of_cell if p.regulated]
        per_cell[cell] = {
            "population": f"the {len(of_cell)} valid held-out pairs of {cell}",
            "pairs": len(of_cell),
            "positives": len(positives),
            "independent_loci": count_loci(positives, span),
            "distinct_measured_genes_among_positives": len({k.gene for k in positives}),
            "distinct_elements_among_positives": len({(k.chrom, k.start, k.end) for k in positives}),
            "chromosomes_among_positives": len({k.chrom for k in positives}),
            "deletion_value_available": deletion_available(cell, model_cells),
            "candidate_second_cell_type": cell != primary,
            "locus_shapes": locus_shapes(positives, span),
        }

    candidates = sorted(c for c in cells if c != primary)
    pooled_keys = [locus_key(p) for p in ps if p.regulated and p.cell != primary]
    # Grouping across cell types as well, by genome position and gene, ignoring which cell measured
    # the positive: two cell types that measured the same place did not sample two places, so this is
    # the conservative count and the one the floor is read against.
    genome_wide = [k._replace(cell="") for k in pooled_keys]
    # The subset the question is actually about. A candidate without a deletion value cannot carry a
    # deletion gain at all, so a pooled count that leans on such cells names loci no gain could be
    # measured at. This is reported beside the floor verdict and does not change it: the floor is read
    # against the pooled count of all candidates, as registered.
    with_feature = [c for c in candidates if per_cell[c]["deletion_value_available"]]
    feature_keys = [k._replace(cell="") for k in pooled_keys if deletion_available(k.cell, model_cells)]
    return {
        "convention": INDEPENDENT_LOCUS_RULE,
        "convention_is_operational": (
            "an operational grouping for deciding whether a stratum has enough distinct places in the "
            "genome to measure, not established biological independence"
        ),
        "span_bp": span,
        "blinded_fields": list(BLINDED_FIELDS),
        "primary_cell_type": primary,
        "candidate_second_cell_types": candidates,
        "per_cell_type": per_cell,
        "pooled_over_candidates": {
            "population": (
                "the measured positives of the held-out cell types other than "
                f"{primary}: {', '.join(candidates)}"
            ),
            "pairs": sum(per_cell[c]["pairs"] for c in candidates),
            "positives": len(pooled_keys),
            "independent_loci_genome_wide": count_loci(genome_wide, span),
            "independent_loci_summed_per_cell_type": sum(per_cell[c]["independent_loci"] for c in candidates),
            "locus_shapes_genome_wide": locus_shapes(genome_wide, span),
            "which_is_read_against_the_floor": "independent_loci_genome_wide",
            "why": (
                "a locus measured in two cell types is one place in the genome, not two, so the "
                "genome-wide grouping is the conservative count; the summed count is reported beside "
                "it and is never the one the floor is read against"
            ),
        },
        "pooled_over_candidates_with_a_deletion_value": {
            "population": (
                "the measured positives of the candidate held-out cell types that have a deletion value "
                "at all, that is those in the lines the deletion table was scored in: "
                f"{', '.join(with_feature) if with_feature else 'none'}"
            ),
            "cell_types": with_feature,
            "pairs": sum(per_cell[c]["pairs"] for c in with_feature),
            "positives": len(feature_keys),
            "independent_loci_genome_wide": count_loci(feature_keys, span),
            "locus_shapes_genome_wide": locus_shapes(feature_keys, span),
            "read_against_the_floor": False,
            "why_it_is_reported": (
                "a candidate without a deletion value cannot carry a deletion gain, so loci contributed "
                "by such cells are places where the quantity in question cannot be measured. This "
                "figure is reported beside the registered floor verdict and does not change it"
            ),
        },
        "floor": POOLED_LOCUS_FLOOR,
    }


def verdict(counted: dict[str, Any], floor: int = POOLED_LOCUS_FLOOR) -> dict[str, Any]:
    """The pre-set floor applied to the pooled count. Below the floor is a no-go: the lane records it
    and stops, computing no gain and making no model request.

    The verdict carries the population it is about and the shape of the loci that carry it, because a
    bare count says neither which cell types are pooled nor whether one chained locus holds most of
    the positives. Neither figure changes the verdict: the floor is read against the pooled count
    exactly as registered."""
    pool = counted["pooled_over_candidates"]
    pooled = pool["independent_loci_genome_wide"]
    shapes = pool.get("locus_shapes_genome_wide", {})
    primary = counted.get("primary_cell_type", PRIMARY_CELL)
    candidates = counted.get("candidate_second_cell_types", [])
    return {
        "floor": floor,
        "floor_set_before_the_count": True,
        "pooled_independent_loci_genome_wide": pooled,
        "population": pool.get("population"),
        "cell_types_pooled": list(candidates),
        "primary_cell_type_excluded": primary,
        "primary_cell_type_is_in_the_pool": False,
        "population_note": (
            f"the pooled count is over the candidate second cell types only; {primary}, which the "
            "committed result already covers, is excluded from it, so the figure is about a second "
            "cell type and not about the benchmark as a whole"
        ),
        "positives_pooled": pool.get("positives"),
        "what_carries_the_count": {
            "largest_locus_positives": shapes.get("largest_locus_positives"),
            "largest_locus_share_of_positives": shapes.get("largest_locus_share_of_positives"),
            "largest_locus_span_mb": shapes.get("largest_locus_span_mb"),
            "positives_per_locus": shapes.get("positives_per_locus"),
            "positives_per_locus_counts": shapes.get("positives_per_locus_counts"),
            "span_mb": shapes.get("span_mb"),
            "note": (
                "reported beside the verdict so a reader can see whether one chained locus carries the "
                "margin over the floor; it does not move the floor or the convention"
            ),
        },
        "meets_floor": pooled >= floor,
        "decision": "eligible to draft a registration" if pooled >= floor else "no-go, stop",
    }


#: Cell types of the candidate pool whose deletion reading this project has already seen, with where
#: the reading sits. Disclosed because the decision to pool these cell types was taken after those
#: readings were seen, which is prior exposure: a stratum that has already been scored cannot be
#: offered as a fresh test of the same question. The readings themselves are cited, never used in any
#: count here, and this lane re-derived none of them.
EXPOSED_CELL_TYPES = {
    "GM12878": (
        "already scored: a held-out stratum of data/results/crispri_published.json "
        "(heldout_published_pairs), reported gain +0.0146"
    ),
    "HCT116": (
        "already scored: data/results/crispri_published.json second_cell_type_hct116, the lane-hct116 "
        "run, reported gain +0.022, verdict 'passes, not replicated'"
    ),
}

#: Which side a locus falls on when one merged locus holds positives from an exposed and an unexposed
#: cell type. It is counted as exposed. Once any part of a locus has been scored, a later reading at
#: that locus is not independent of the reading already seen, so counting it as fresh would overstate
#: the never-exposed stratum. The rule is stated before the split is computed.
MERGED_LOCUS_SIDE = (
    "a locus holding positives from both an exposed and an unexposed cell type is counted as exposed, "
    "because once any part of a locus has been scored a later reading there is not independent of the "
    "reading already seen; counting it as fresh would overstate the never-exposed stratum"
)


def exposure(
    pairs: Iterable[Any],
    exposed: dict[str, str] | None = None,
    primary: str = PRIMARY_CELL,
    span: int = INDEPENDENT_LOCUS_SPAN,
) -> dict[str, Any]:
    """The pooled loci split into those whose cell types have already been scored and those never
    scored, computed under the pooled convention itself and never by adding per-cell-type counts.

    Per-cell-type counts do not sum to the pooled count: loci that appear in more than one cell type
    merge under the convention, so a split made by addition would double-count them."""
    exposed = EXPOSED_CELL_TYPES if exposed is None else exposed
    ps = [p for p in pairs if p.regulated and p.cell != primary]
    keys = [locus_key(p) for p in ps]
    ids = group([k._replace(cell="") for k in keys], span)
    cells_of: dict[int, set[str]] = {}
    size_of: dict[int, int] = {}
    for gid, k in zip(ids, keys, strict=True):
        cells_of.setdefault(gid, set()).add(k.cell)
        size_of[gid] = size_of.get(gid, 0) + 1
    exposed_ids = [g for g, cs in cells_of.items() if cs & set(exposed)]
    fresh_ids = [g for g in cells_of if g not in set(exposed_ids)]
    mixed = [g for g in exposed_ids if cells_of[g] - set(exposed)]
    multi = [g for g, cs in cells_of.items() if len(cs) > 1]
    fresh_cells = sorted({c for g in fresh_ids for c in cells_of[g]})
    # Why the per-cell-type counts sum to more than the pooled count: each pooled locus absorbs the
    # per-cell-type loci that fall inside it, and every absorption after the first is one the summed
    # figure counts twice. Computed, so the arithmetic can be checked rather than asserted.
    pieces: set[tuple[int, str, int]] = set()
    for cell in sorted({k.cell for k in keys}):
        of_cell = [(i, k) for i, k in enumerate(keys) if k.cell == cell]
        cell_ids = group([k for _, k in of_cell], span)
        for (i, _), piece in zip(of_cell, cell_ids, strict=True):
            pieces.add((ids[i], cell, piece))
    absorbed: dict[int, int] = {}
    for pooled_id, _, _ in pieces:
        absorbed[pooled_id] = absorbed.get(pooled_id, 0) + 1
    reconcile = {
        "independent_loci_summed_per_cell_type": len(pieces),
        "independent_loci_pooled_genome_wide": len(cells_of),
        "difference": len(pieces) - len(cells_of),
        "pooled_loci_that_absorb_more_than_one_per_cell_type_locus": sum(
            1 for n in absorbed.values() if n > 1
        ),
        "note": (
            "the summed figure counts a locus once for every separate per-cell-type group inside it; "
            "the pooled figure counts it once. The difference is the number of absorptions, which is "
            "why a split between exposed and unexposed loci must never be made by addition"
        ),
    }
    fresh_keys = [k._replace(cell="") for k in keys if k.cell not in exposed]
    return {
        "computed_under": (
            "the pooled convention itself: the pooled positives are grouped once, genome-wide, and "
            "each locus is then assigned a side. No per-cell-type count is added to another"
        ),
        "why_addition_would_be_wrong": (
            "per-cell-type locus counts do not sum to the pooled count, because a locus present in "
            "more than one cell type merges under the convention; adding the per-cell figures would "
            "count such a locus once per cell type"
        ),
        "exposed_cell_types": dict(exposed),
        "merged_locus_side": MERGED_LOCUS_SIDE,
        "pooled_loci": len(cells_of),
        "pooled_positives": len(keys),
        "already_exposed_loci": len(exposed_ids),
        "already_exposed_positives": sum(size_of[g] for g in exposed_ids),
        "never_exposed_loci": len(fresh_ids),
        "never_exposed_positives": sum(size_of[g] for g in fresh_ids),
        "never_exposed_cell_types": fresh_cells,
        "loci_in_more_than_one_cell_type": len(multi),
        "reconciliation_of_the_summed_and_pooled_counts": reconcile,
        "loci_counted_as_exposed_only_because_they_merge_with_an_exposed_cell_type": len(mixed),
        "never_scored_stratum_on_its_own": {
            "population": (
                "the measured positives of the candidate cell types never scored by this project: "
                f"{', '.join(c for c in sorted(set(k.cell for k in keys)) if c not in exposed) or 'none'}"
            ),
            "positives": len(fresh_keys),
            "independent_loci_genome_wide": count_loci(fresh_keys, span),
            "locus_shapes_genome_wide": locus_shapes(fresh_keys, span),
            "status": (
                "descriptive only: this stratum is reported because the pool's exposure makes it the "
                "only unexposed part, and it is not judged against the pooled floor"
            ),
        },
    }


# --- the power gate, registered before it is computed (reviewer, 2026-10-02) -------------------------

#: The share of subsamples whose interval must exclude zero for a run to be worth paying for. Below
#: this line a run is more likely inconclusive than not, and the recommendation is not to spend. The
#: line is 0.5 and is written down here before the share is computed, so it cannot be moved afterwards.
POWER_GATE_SHARE_FLOOR = 0.5

#: What the gate is, in words, fixed before any share is computed.
POWER_GATE_DESIGN = (
    "subsample the primary cell type's held-out pairs to the second-cell structure: the same number of "
    "independent loci, the same number of measured positives, and the same positives-per-locus "
    "profile, under the locus convention registered here. At the primary cell type's observed effect, "
    "report the share of subsamples whose interval excludes zero. A share below "
    "POWER_GATE_SHARE_FLOOR means the run is more likely inconclusive than not, and the recommendation "
    "is not to spend"
)

#: Reading fixed in advance, so neither outcome can be restated once it is seen.
READINGS = {
    "lower_bound_above_zero": "replicated outside K562 on these cell types",
    "interval_covering_zero": "no difference detected",
    "never_say": (
        "no effect: an interval covering zero is a failure to detect a difference at this number of "
        "loci, not evidence that the difference is absent"
    ),
}


def profile_match(
    target: dict[str, int], source: Iterable[LocusKey], span: int = INDEPENDENT_LOCUS_SPAN
) -> dict[str, Any]:
    """Can the primary cell type's loci supply the second cell structure's positives-per-locus profile
    exactly, one source locus per required locus? Labels only: it compares locus sizes and nothing else.

    `target` is a positives-per-locus count map as `locus_shapes` reports it. The match is greedy from
    the largest required size, which is exact for this question: a required size can only be met by a
    source locus of that size, so the sizes are matched independently."""
    available: dict[int, int] = {}
    for n in [s for s in _sizes(source, span)]:
        available[n] = available.get(n, 0) + 1
    shortfall: dict[str, dict[str, int]] = {}
    for size_text, needed in sorted(target.items(), key=lambda kv: -int(kv[0])):
        size = int(size_text)
        have = available.get(size, 0)
        if have < needed:
            shortfall[size_text] = {"required": needed, "available": have, "missing": needed - have}
        available[size] = max(0, have - needed)
    return {
        "required_positives_per_locus_counts": dict(target),
        "available_positives_per_locus_counts": {
            str(n): c for n, c in sorted(_count_sizes(source, span).items())
        },
        "exact_profile_available": not shortfall,
        "shortfall_by_locus_size": shortfall,
        "rule": (
            "a required locus of n positives can only be met by a source locus of exactly n positives, "
            "so the sizes are matched independently and the largest first; any relaxation of the "
            "profile is a change to the gate and must be registered before it is used, not chosen "
            "after seeing this"
        ),
    }


def _sizes(keys: Iterable[LocusKey], span: int = INDEPENDENT_LOCUS_SPAN) -> list[int]:
    """The number of positives in each locus of `keys`."""
    ks = list(keys)
    counts: dict[int, int] = {}
    for gid in group(ks, span):
        counts[gid] = counts.get(gid, 0) + 1
    return sorted(counts.values(), reverse=True)


def _count_sizes(keys: Iterable[LocusKey], span: int = INDEPENDENT_LOCUS_SPAN) -> dict[int, int]:
    """How many loci of `keys` hold each number of positives."""
    out: dict[int, int] = {}
    for n in _sizes(keys, span):
        out[n] = out.get(n, 0) + 1
    return out
